#if defined(_MSC_VER)
#pragma warning ( disable : 4786 )
#endif


#include "itkImageFileReader.h"
#include "itkImageFileWriter.h"
#include "itkRescaleIntensityImageFilter.h"

#include "itkImageConstIterator.h"
#include "itkConstantBoundaryCondition.h"
#include "itkConstNeighborhoodIterator.h"
#include "itkImageRegionIterator.h"

#include "itkGradientRecursiveGaussianImageFilter.h"

#include "math.h"
#include "string.h"

int main( int argc, char ** argv )
{


  //%%%%%%%%%%%%%%%%%%% ARGUMENTS PARSING %%%%%%%%%%%%%%%%%%%

  if ( argc < 8 )
    {
    std::cerr << "Missing parameters. " << std::endl;
    std::cerr << "Usage: " << std::endl;
    std::cerr << argv[0]
              << " inputphi(image) inputEdge(image)"
              << " nu rho deltat curavtureFactor NbIter"
              << " outputImageFile(image)"
              << std::endl;
    return -1;
    }
  // arguments reading and parsing
  double nu = atof(argv[3]);
  double rho = atof(argv[4]);
  double deltat = atof(argv[5]);
  double curvatureFactor = atof(argv[6]);
  int nbiter = atoi(argv[7]);
  std::stringstream outputDir;
  outputDir << argv[8];

  //%%%%%%%%%%%%%%%%%%% TYPEDEFS %%%%%%%%%%%%%%%%%%%

  // dimension of input image (and phi)
  const   unsigned int   Dimension = 2;

  // scalar field type
  typedef double                              PixelType;
  typedef itk::Image< PixelType, Dimension > ImageType;

  //vector field type
  typedef itk::CovariantVector< PixelType, Dimension  >   VectorPixelType;
  typedef itk::Image< VectorPixelType, Dimension  >   VectorImageType;
  typedef itk::GradientRecursiveGaussianImageFilter< ImageType,VectorImageType>
          GradientFilterType;

  // neighborhood iterators
  typedef itk::ConstantBoundaryCondition< ImageType >  BoundaryConditionType;
  typedef itk::ConstNeighborhoodIterator< ImageType, BoundaryConditionType >
          NeighborhoodIteratorType;
  typedef itk::ImageRegionIterator< ImageType>        IteratorType;
  // simple iterators
  typedef itk::ImageRegionConstIterator< ImageType > ConstIteratorType;
  typedef itk::ImageRegionConstIterator< VectorImageType >
          ConstVectorIteratorType;

  // reader
  typedef itk::ImageFileReader< ImageType >  ReaderType;
  // inputs rescaler
  typedef itk::RescaleIntensityImageFilter<
               ImageType, ImageType > RescaleFilterType;

  // output image type
  typedef unsigned char                            WritePixelType;
  typedef itk::Image< WritePixelType, Dimension >  WriteImageType;
  // writer
  typedef itk::ImageFileWriter< WriteImageType >   WriterType;
  // output rescaler
  typedef itk::RescaleIntensityImageFilter<
               ImageType, WriteImageType > RescaleOutputFilterType;



  //%%%%%%%%%%%%%%%%%%% INPUTS READING AND RESCALING %%%%%%%%%%%%%%%%%%%

  // reading phi0, store the rescaled version in inputphi
  ReaderType::Pointer reader = ReaderType::New();
  reader->SetFileName( argv[1] );
  try
    {
    reader->Update();
    }
  catch ( itk::ExceptionObject &err)
    {
    std::cout << "ExceptionObject caught !" << std::endl;
    std::cout << err << std::endl;
    return -1;
    }  

  // rescale phi0
  RescaleFilterType::Pointer rescalerInputs = RescaleFilterType::New();
  rescalerInputs->SetOutputMinimum( 0.01671 );
  rescalerInputs->SetOutputMaximum( 2.0241 ); // value taken from matlab
  rescalerInputs->SetInput(reader->GetOutput());
  rescalerInputs->Update();

  // pointer to phi
  ImageType::Pointer inputphi = rescalerInputs->GetOutput();
  // stop propagating the update process
  inputphi->DisconnectPipeline();



  // reading edge detector image, store the rescaled version in inputG
    // we use the same reader
  reader->SetFileName( argv[2] );
  try
    {
    reader->Update();
    }
  catch ( itk::ExceptionObject &err)
    {
    std::cout << "ExceptionObject caught !" << std::endl;
    std::cout << err << std::endl;
    return -1;
    }

  //rescale g
  rescalerInputs->SetOutputMinimum( 0. );
  rescalerInputs->SetOutputMaximum( 1. );
  rescalerInputs->SetInput(reader->GetOutput());
  rescalerInputs->Update();
  // pointer to edge indicator g
  ImageType::Pointer inputG = rescalerInputs->GetOutput();
  // stop propagating the update process
  inputG->DisconnectPipeline();
  // iterator on g (scalar)
  ConstIteratorType edgit( inputG, inputG->GetRequestedRegion() );



  //%%%%%%%%%%%%%%%%%%% COMPUTATIONS %%%%%%%%%%%%%%%%%%%

  // computation of grad(g) (DimensionxDimension)
  GradientFilterType::Pointer gradient = GradientFilterType::New();
  gradient->SetInput( inputG );
  gradient->Update();
  // iterator on grad(g) (vector)
  ConstVectorIteratorType edgderivit( gradient->GetOutput(),
                                      gradient->GetOutput()
                                              ->GetLargestPossibleRegion() );
    // gradient vector at a given point
  VectorPixelType edgederivVector;


  // ITERATIVE LOOP :

  // allocation of phi(t+1)
  ImageType::Pointer outputphi = ImageType::New();
  outputphi->SetRegions(inputphi->GetRequestedRegion());
  outputphi->SetSpacing(inputphi->GetSpacing());
  outputphi->Allocate();

  // offsets definition for neighborhood computations :
  NeighborhoodIteratorType::OffsetType minusx = {{-1,0}};
  NeighborhoodIteratorType::OffsetType plusx = {{1,0}};
  NeighborhoodIteratorType::OffsetType minusy = {{0,-1}};
  NeighborhoodIteratorType::OffsetType plusy = {{0,1}};
  NeighborhoodIteratorType::OffsetType plusxplusy = {{1,1}};
  NeighborhoodIteratorType::OffsetType minusxminusy = {{-1,-1}};
  NeighborhoodIteratorType::OffsetType plusxminusy = {{1,-1}};
  NeighborhoodIteratorType::OffsetType minusxplusy = {{-1,1}};

  // get the spacing in x and y for derivatives computations
  ImageType::SpacingType inputSpacing = reader->GetOutput()->GetSpacing();
  double spacingx = inputSpacing[0];
  double spacingy = inputSpacing[1];

  // temporary variables for the iterations
    //derivatives
  PixelType dplusx, dminusx, dplusy, dminusy,
            dx, dy,
            d2x, d2y,
            dxy,
            dxsquare, dysquare;
    //equation members
  PixelType Hg, Upwind, nextphi;

#ifdef _DEBUG
  // make sure smart pointers are smart
  std::cout << "inputphi reference count : "
            << inputphi->GetReferenceCount()
            << std::endl;
  std::cout << "output reference count : "
            << outputphi->GetReferenceCount()
            << std::endl;
  // see when we start iterating
  std::cout << "begin of iteration loop" << std::endl;
#endif  

  // Writing images for the video

  std::stringstream OutName;

  // rescaling output (phi(t=tmax))
  std::cout << "begin of iteration loop" << std::endl;

  int numImg = 0;
  for (int iter = 1; iter<=nbiter; ++iter)
    {



    // input iterator (on phi)
    NeighborhoodIteratorType::RadiusType radius;
    radius.Fill(1);
    NeighborhoodIteratorType it( radius, inputphi,
                                 inputphi->GetRequestedRegion() );

    // output iterator ( on phi(t+1) )
    IteratorType out(outputphi, inputphi->GetRequestedRegion());

    // we change the evolution of the surface from graph to level set
    if (iter == nbiter/2)
    curvatureFactor = curvatureFactor/1000000;


    // one iteration :
    for ( it.GoToBegin(), edgit.GoToBegin(), edgderivit.GoToBegin(),
            out.GoToBegin();
          !it.IsAtEnd();
          ++it, ++out, ++edgderivit, ++edgit )
      {

      // derivatives computation :
      dminusx = (it.GetCenterPixel()-it.GetPixel(minusx))/spacingx;
      dminusy = (it.GetCenterPixel()-it.GetPixel(minusy))/spacingy;

      dplusx = (it.GetPixel(plusx)-it.GetCenterPixel())/spacingx;
      dplusy = (it.GetPixel(plusy)-it.GetCenterPixel())/spacingy;

      dx = (it.GetPixel(plusx) - it.GetPixel(minusx))/(2*spacingx);
      dy = (it.GetPixel(plusy) - it.GetPixel(minusy))/(2*spacingy);

      d2x = (dplusx-dminusx)/spacingx;
      d2y = (dplusy-dminusy)/spacingy;

      dxy = (it.GetPixel(plusxplusy)
             + it.GetPixel(minusxminusy)
             - it.GetPixel(plusxminusy)
             - it.GetPixel(minusxplusy))
            / (4*spacingx*spacingy);

      dxsquare = dx*dx;
      dysquare = dy*dy;


      // evolution equation
        // Hg = g*H (with H mean curvature of graph of phi)
      Hg = edgit.Get() * (   (curvatureFactor + dxsquare)*d2y
                           + (curvatureFactor + dysquare)*d2x
                           - 2 * (dx*dy*dxy))
           / (curvatureFactor + dxsquare + dysquare ) ;

        // grad(g)
      edgederivVector = edgderivit.Get();

        // upwind differntiation term : grad(I).*grad(g)
      Upwind = std::min(-edgederivVector[0],(double)0.)*dplusx
             + std::max(-edgederivVector[0],(double)0.)*dminusx
             + std::min(-edgederivVector[1],(double)0.)*dplusy
             + std::max(-edgederivVector[1],(double)0.)*dminusy;


      // value of phi(t+1)
      nextphi = it.GetCenterPixel()+deltat*(nu*Hg-rho*Upwind);

      out.Set(nextphi);

      }

    if ( (iter % 100) == 0 || (iter==1) || (iter==nbiter) )
      {

      RescaleOutputFilterType::Pointer rescaler = RescaleOutputFilterType::New();
      WriterType::Pointer writer = WriterType::New();

      rescaler->SetOutputMinimum(   0 );
      rescaler->SetOutputMaximum( 255 );
      rescaler->SetInput(inputphi);
      if (iter==nbiter)
        rescaler->SetInput(outputphi);
      rescaler->Update();

      OutName.str("");
      OutName << outputDir.str()
              << "phi"
              << numImg
              << ".png"
              << std::ends;
      std::cout << OutName.str() << std::endl ;
      // write output in png
      writer->SetFileName( OutName.str() );
      writer->SetInput(rescaler->GetOutput());
      try
        {
        writer->Update();
        }
      catch ( itk::ExceptionObject &err)
        {
        std::cout << "ExceptionObject caught !" << std::endl;
        std::cout << err << std::endl;
        return -1;
        }
      ++numImg;
      }



    // input is now output
      //(smart pointer dereference automatically
      // the memory previously pointed by inputphi,
      // and tempImagePoint gets out of scope at the end of the iteration)
    ImageType::Pointer tempImagePoint = inputphi;
    inputphi = outputphi;
    outputphi = tempImagePoint;


#ifdef _DEBUG

    // make sure smart pointers are smart
    std::cout << "inputphi reference count : "
              << inputphi->GetReferenceCount()
              << std::endl;
    std::cout << "temp reference count : "
              << tempImagePoint->GetReferenceCount()
              << std::endl;
    std::cout << "output reference count : "
              << outputphi->GetReferenceCount()
              << std::endl;
    std::cout << "inputphi reference count : "
              << inputphi->GetReferenceCount()
              << std::endl;
    std::cout << "temp reference count : "
              << tempImagePoint->GetReferenceCount()
              << std::endl;
#endif

    }


  // see when we are done iterating
  std::cout << "end of iteration loop" << std::endl;
  // make sure smart pointers are smart
  std::cout << "inputphi reference count : "
            << inputphi->GetReferenceCount()
            << std::endl;
  std::cout << "output reference count : "
            << outputphi->GetReferenceCount()
            << std::endl;

  return 0;
}
