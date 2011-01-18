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


#include "itkMembraneEdgeDetectorFilter.h"


int main( int argc, char ** argv )
{


  //%%%%%%%%%%%%%%%%%%% ARGUMENTS PARSING %%%%%%%%%%%%%%%%%%%

  if ( argc < 8 )
    {
    std::cerr << "Missing parameters. " << std::endl;
    std::cerr << "Usage: " << std::endl;
    std::cerr << argv[0]
              << " inputphi(image)"
              << " Beta Sigma Pow dT NbIter"
              << " outputImageFile(image)"
              << std::endl;
    return -1;
    }
  // arguments reading and parsing
  double Beta =   atof(argv[2]);
  double Sigma =  atof(argv[3]);
  double Pow =    atof(argv[4]);
  double deltat = atof(argv[5]);
  int nbiter =    atoi(argv[6]);






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

  typedef itk::MembraneEdgeDetectorFilter<
      ImageType, ImageType > EdgeDetectorType;


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

  RescaleFilterType::Pointer inputRescale = RescaleFilterType::New();
  inputRescale->SetInput(reader->GetOutput());
  inputRescale->SetOutputMinimum( 0.);
  inputRescale->SetOutputMaximum( 1.);
  inputRescale->Update();

  // pointer to inputImage
  ImageType::Pointer inputImage = reader->GetOutput();
  // stop propagating the update process
  inputImage->DisconnectPipeline();


  //EDGE DETECTION AND GRADIENT
  EdgeDetectorType::Pointer edgeDetector = EdgeDetectorType::New();
  GradientFilterType::Pointer gradient = GradientFilterType::New();
  ImageType::Pointer edgeImage;
  edgeDetector->SetSigma(Sigma);
  edgeDetector->SetPow(Pow);
  edgeDetector->SetBeta(Beta);



  //%%%%%%%%%%%%%%%%%%% COMPUTATIONS %%%%%%%%%%%%%%%%%%%



  // ITERATIVE LOOP :

  // allocation of phi(t+1)
  ImageType::Pointer outputImage = ImageType::New();
  outputImage->SetRegions(inputImage->GetRequestedRegion());
  outputImage->SetSpacing(inputImage->GetSpacing());
  outputImage->Allocate();

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
  PixelType Hg, Upwind, nextValue;


  for (int iter = 0; iter<nbiter; ++iter)
    {
    // EDGE DETECTOR (g) :
    edgeDetector->SetInput(inputImage);
    edgeDetector->Update();
    // g :
    edgeImage = edgeDetector->GetOutput();
    // iterator on g (scalar)
    ConstIteratorType edgit( edgeImage,
                             edgeImage->GetLargestPossibleRegion() );

    // computation of grad(g) (DimensionxDimension)
    gradient->SetInput( edgeImage );
    gradient->Update();
    // gradient vector at a given point
    VectorPixelType edgederivVector;
    // iterator on grad(g) (vector)
    ConstVectorIteratorType edgderivit( gradient->GetOutput(),
                                        gradient->GetOutput()
                                                ->GetLargestPossibleRegion());


    // input iterator (on I)
    NeighborhoodIteratorType::RadiusType radius;
    radius.Fill(1);
    NeighborhoodIteratorType it( radius, inputImage,
                                 inputImage->GetRequestedRegion() );

    // output iterator ( on I(t+1) )
    IteratorType out(outputImage, inputImage->GetRequestedRegion());


    // ONE ITERATION :
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
      Hg = edgit.Get() * (   (1. + dxsquare)*d2y
                           + (1. + dysquare)*d2x
                           - 2 * (dx*dy*dxy))
           / (1. + dxsquare + dysquare ) ;

        // grad(g)
      edgederivVector = edgderivit.Get();

        // upwind differntiation term : grad(I).*grad(g)
      Upwind = std::min(-edgederivVector[0],(double)0.)*dplusx
             + std::max(-edgederivVector[0],(double)0.)*dminusx
             + std::min(-edgederivVector[1],(double)0.)*dplusy
             + std::max(-edgederivVector[1],(double)0.)*dminusy;


      // value of phi(t+1)
      nextValue = it.GetCenterPixel()+deltat*(Hg-Upwind);

      out.Set(nextValue);
      }
    // input is now output
      //(smart pointer dereference automatically
      // the memory previously pointed by inputphi,
      // and tempImagePoint gets out of scope at the end of the iteration)
    ImageType::Pointer tempImagePoint = inputImage;
    inputImage = outputImage;
    outputImage = tempImagePoint;
    }


  // rescaling output (phi(t=tmax))
  RescaleOutputFilterType::Pointer rescaler = RescaleOutputFilterType::New();
  rescaler->SetOutputMinimum(   0 );
  rescaler->SetOutputMaximum( 255 );
  rescaler->SetInput(outputImage);

  // write output in png
  WriterType::Pointer writer = WriterType::New();
  writer->SetFileName( argv[7] );
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

  return 0;
}
