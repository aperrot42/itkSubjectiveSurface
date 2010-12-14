
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

int main( int argc, char ** argv )
{

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


  //%%%%%%%%%%%%%%%%%%% TYPEDEFS %%%%%%%%%%%%%%%%%%%

  // dimension of input image (and phi)
  const   unsigned int   Dimension = 3;

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
  RescaleFilterType::Pointer rescalerInputPhy = RescaleFilterType::New();
  rescalerInputPhy->SetOutputMinimum( 0.01671 );
  rescalerInputPhy->SetOutputMaximum( 2.0241 ); // value taken from matlab
  rescalerInputPhy->SetInput(reader->GetOutput());
  rescalerInputPhy->Update();

  // pointer to phi
  ImageType::Pointer inputphi = rescalerInputPhy->GetOutput();
  // stop propagating the update process
  inputphi->DisconnectPipeline();



  // reading g, store the rescaled version in inputg
  ReaderType::Pointer readerEdge = ReaderType::New();
  readerEdge->SetFileName( argv[2] );
  try
    {
    readerEdge->Update();
    }
  catch ( itk::ExceptionObject &err)
    {
    std::cout << "ExceptionObject caught !" << std::endl;
    std::cout << err << std::endl;
    return -1;
    }

  // rescale inputG
  RescaleFilterType::Pointer rescalerInputG = RescaleFilterType::New();
  rescalerInputG->SetOutputMinimum( 0.01671 );
  rescalerInputG->SetOutputMaximum( 2.0241 ); // value taken from matlab
  rescalerInputG->SetInput(reader->GetOutput());
  rescalerInputG->Update();

  // pointer to phi
  ImageType::Pointer inputG = rescalerInputG->GetOutput();
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


  // offsets definition :
  NeighborhoodIteratorType::OffsetType minusx = {{-1,0,0}};
  NeighborhoodIteratorType::OffsetType plusx = {{1,0,0}};

  NeighborhoodIteratorType::OffsetType minusy = {{0,-1,0}};
  NeighborhoodIteratorType::OffsetType plusy = {{0,1,0}};

  NeighborhoodIteratorType::OffsetType minusz = {{0,0,-1}};
  NeighborhoodIteratorType::OffsetType plusz = {{0,0,1}};


  NeighborhoodIteratorType::OffsetType plusxplusy = {{1,1,0}};
  NeighborhoodIteratorType::OffsetType minusxminusy = {{-1,-1,0}};
  NeighborhoodIteratorType::OffsetType plusxminusy = {{1,-1,0}};
  NeighborhoodIteratorType::OffsetType minusxplusy = {{-1,1,0}};

  NeighborhoodIteratorType::OffsetType plusxplusz = {{1,0,1}};
  NeighborhoodIteratorType::OffsetType minusxminusz = {{-1,0,-1}};
  NeighborhoodIteratorType::OffsetType plusxminusz = {{1,0,-1}};
  NeighborhoodIteratorType::OffsetType minusxplusz = {{-1,0,1}};

  NeighborhoodIteratorType::OffsetType plusyplusz = {{0,1,1}};
  NeighborhoodIteratorType::OffsetType minusyminusz = {{0,-1,-1}};
  NeighborhoodIteratorType::OffsetType plusyminusz = {{0,1,-1}};
  NeighborhoodIteratorType::OffsetType minusyplusz = {{0,-1,1}};


  // get the spacing in x and y for derivatives computations
  ImageType::SpacingType inputSpacing = reader->GetOutput()->GetSpacing();
  float spacingx = inputSpacing[0];
  float spacingy = inputSpacing[1];
  float spacingz = inputSpacing[2];

  // temporary variables for the iterations
  PixelType dplusx, dminusx, dplusy, dminusy, dplusz, dminusz,
            dx, dy, dz,
            d2x, d2y, d2z,
            dxy, dxz, dyz,
            dxsquare, dysquare, dzsquare;
  //equation members
  PixelType Hg, Upwind, nextphi;


  std::cout << "begin of iteration loop" << std::endl;

  for (int iter = 0; iter<nbiter; ++iter)
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
  for (it.GoToBegin(),edgit.GoToBegin(),edgderivit.GoToBegin(), out.GoToBegin();
       !it.IsAtEnd();
       ++it, ++out, ++edgderivit, ++edgit)
    {

    // derivatives computation :
    dminusx = (it.GetCenterPixel()-it.GetPixel(minusx))/spacingx;
    dminusy = (it.GetCenterPixel()-it.GetPixel(minusy))/spacingy;
    dminusz = (it.GetCenterPixel()-it.GetPixel(minusz))/spacingz;


    dplusx = (it.GetPixel(plusx)-it.GetCenterPixel())/spacingx;
    dplusy = (it.GetPixel(plusy)-it.GetCenterPixel())/spacingy;
    dplusz = (it.GetPixel(plusz)-it.GetCenterPixel())/spacingz;


    dx = (it.GetPixel(plusx) - it.GetPixel(minusx))/(2*spacingx);
    dy = (it.GetPixel(plusy) - it.GetPixel(minusy))/(2*spacingy);
    dz = (it.GetPixel(plusz) - it.GetPixel(minusz))/(2*spacingz);


    d2x = dplusx-dminusx;
    d2y = dplusy-dminusy;
    d2z = dplusz-dminusz;


    dxy = (it.GetPixel(plusxplusy)
           + it.GetPixel(minusxminusy)
           - it.GetPixel(plusxminusy)
           - it.GetPixel(minusxplusy))
          / (4*spacingx*spacingy);
    dxz = (it.GetPixel(plusxplusz)
           + it.GetPixel(minusxminusz)
           - it.GetPixel(plusxminusz)
           - it.GetPixel(minusxplusz))
          / (4*spacingx*spacingz);
    dyz = (it.GetPixel(plusyplusz)
          + it.GetPixel(minusyminusz)
          - it.GetPixel(plusyminusz)
          - it.GetPixel(minusyplusz))
          / (4*spacingy*spacingz);


    dxsquare = dx*dx;
    dysquare = dy*dy;
    dzsquare = dz*dz;


    // evolution equation
    Hg = edgit.Get() * (  (curvatureFactor + dxsquare + dysquare)*d2z
                         +(curvatureFactor + dxsquare + dzsquare)*d2y
                         +(curvatureFactor + dzsquare + dysquare)*d2x
                         -2 * (dx*dz*dxz+dx*dy*dxy+dy*dz*dyz))
            / (curvatureFactor + dxsquare + dysquare + dzsquare) ;

    edgederivVector = edgderivit.Get();

    Upwind =  std::min(-edgederivVector[0],0.)*dplusx
            + std::max(-edgederivVector[0],0.)*dminusx
            + std::min(-edgederivVector[1],0.)*dplusy
            + std::max(-edgederivVector[1],0.)*dminusy
            + std::min(-edgederivVector[2],0.)*dplusz
            + std::max(-edgederivVector[2],0.)*dminusz;

    nextphi = it.GetCenterPixel()+deltat*(nu*Hg-rho*Upwind);

    out.Set(nextphi);
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
#endif
  }

#ifdef _DEBUG
// see when we are done iterating
std::cout << "end of iteration loop" << std::endl;
// make sure smart pointers are smart
std::cout << "inputphi reference count : "
          << inputphi->GetReferenceCount()
          << std::endl;
std::cout << "output reference count : "
          << outputphi->GetReferenceCount()
          << std::endl;
#endif


// rescaling output (phi(t=tmax))
RescaleOutputFilterType::Pointer rescaler = RescaleOutputFilterType::New();
rescaler->SetOutputMinimum(   0 );
rescaler->SetOutputMaximum( 255 );
rescaler->SetInput(outputphi);

// write output in png
WriterType::Pointer writer = WriterType::New();
writer->SetFileName( argv[8] );
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
