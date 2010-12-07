#if defined(_MSC_VER)
#pragma warning ( disable : 4786 )
#endif


#include "itkImageFileReader.h"
#include "itkImageFileWriter.h"
#include "itkRescaleIntensityImageFilter.h"

#include "itkImageConstIterator.h"
#include "itkConstNeighborhoodIterator.h"
#include "itkImageRegionIterator.h"

#include "itkGradientRecursiveGaussianImageFilter.h"

#include "math.h"


int main( int argc, char ** argv )
{
  if ( argc < 7 )
    {
    std::cerr << "Missing parameters. " << std::endl;
    std::cerr << "Usage: " << std::endl;
    std::cerr << argv[0]
              << " inputPhy(image) inputEdge(image) nu rho deltat curavtureFactor outputImageFile(image)"
              << std::endl;
    return -1;
    }

  double nu = atof(argv[3]);
  double rho = atof(argv[4]);
  double deltat = atof(argv[5]);
  double curvatureFactor = atof(argv[6]);



  const   unsigned int   Dimension = 2;


  // scalar field type
  typedef double                              PixelType;
  typedef itk::Image< PixelType, Dimension > ImageType;
  typedef itk::ImageFileReader< ImageType >  ReaderType;
  //vector field type
  typedef itk::CovariantVector< PixelType, Dimension  >   VectorPixelType;
  typedef itk::Image< VectorPixelType, Dimension  >   VectorImageType;
  typedef itk::GradientRecursiveGaussianImageFilter< ImageType,VectorImageType>
          GradientFilterType;


  // neighborhood iterators
  typedef itk::ConstNeighborhoodIterator< ImageType > NeighborhoodIteratorType;
  typedef itk::ImageRegionIterator< ImageType>        IteratorType;
  // simple iterators
  typedef itk::ImageRegionConstIterator< ImageType > ConstIteratorType;
  typedef itk::ImageRegionConstIterator< VectorImageType > ConstVectorIteratorType;


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



  // reading edge detector image (Dimension)
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

  // rescaling inputs
  typedef itk::RescaleIntensityImageFilter<
               ImageType, ImageType > RescaleFilterType;

  RescaleFilterType::Pointer rescalerInputs = RescaleFilterType::New();



  // rescale phy
  rescalerInputs->SetOutputMinimum( 0. );
  rescalerInputs->SetOutputMaximum( 2. ); // value taken from matlab
  rescalerInputs->SetInput(reader->GetOutput());
  rescalerInputs->Update();
  // pointer to phy
  ImageType::Pointer inputPhy = rescalerInputs->GetOutput();
  // stop propagating the update process
  inputPhy->DisconnectPipeline();



  //rescale g
  rescalerInputs->SetOutputMinimum( 0. );
  rescalerInputs->SetOutputMaximum( 1. );
  rescalerInputs->SetInput(readerEdge->GetOutput());
  rescalerInputs->Update();
  // pointer to edge indicator g
  ImageType::Pointer inputG = rescalerInputs->GetOutput();
  // stop propagating the update process
  inputG->DisconnectPipeline();
  // iterator on g (scalar)
  ConstIteratorType edgit( inputG, inputG->GetRequestedRegion() );


  // computation of grad(g) (DimensionxDimension)
  GradientFilterType::Pointer gradient = GradientFilterType::New();
  gradient->SetInput( inputG );
  gradient->Update();
  // iterator on grad(g) (vector)
  ConstVectorIteratorType edgderivit( gradient->GetOutput(),gradient->GetOutput()->GetLargestPossibleRegion() );
    // gradient vector at a given point
  VectorPixelType edgederivVector;





  // allocation of phy(t+1)
  ImageType::Pointer outputPhy = ImageType::New();
  outputPhy->SetRegions(inputPhy->GetRequestedRegion());
  outputPhy->SetSpacing(inputPhy->GetSpacing());
  outputPhy->Allocate();


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
  PixelType Hg, Upwind, nextPhy;



#ifdef _DEBUG
  std::cout << "inputPhy reference count : " << inputPhy->GetReferenceCount() << std::endl;
  std::cout << "output reference count : " << outputPhy->GetReferenceCount() << std::endl;

  std::cout << "begin of iteration loop" << std::endl;
#endif

  for (int iter = 0; iter<10000; ++iter)
    {
    // input iterator (on phy)
    NeighborhoodIteratorType::RadiusType radius;
    radius.Fill(1);
    NeighborhoodIteratorType it( radius, inputPhy,
                                 inputPhy->GetRequestedRegion() );

    // output iterator ( on phy(t+1) )
    IteratorType out(outputPhy, inputPhy->GetRequestedRegion());

    if (iter == 5000)
    curvatureFactor = curvatureFactor/1000000; //to be deleted


    // one iteration :
    for (it.GoToBegin(),edgit.GoToBegin(),edgderivit.GoToBegin(), out.GoToBegin();
       !it.IsAtEnd();
       ++it, ++out, ++edgderivit, ++edgit)
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

      Hg = edgit.Get() * (   (curvatureFactor + dxsquare)*d2y
                           + (curvatureFactor + dysquare)*d2x
                           - 2 * (dx*dy*dxy))
           / (curvatureFactor + dxsquare + dysquare ) ;//ok

      edgederivVector = edgderivit.Get();

        // upwind differntiation term : grad(I).*grad(g)
      Upwind =   std::min(-edgederivVector[0],(double)0.)*dplusx + std::max(-edgederivVector[0],(double)0.)*dminusx
               + std::min(-edgederivVector[1],(double)0.)*dplusy + std::max(-edgederivVector[1],(double)0.)*dminusy;

      // value of phy(t+1)
      nextPhy = it.GetCenterPixel()+deltat*(nu*Hg-rho*Upwind);

      out.Set(nextPhy);
      }
    // input is now output (smart pointer dereference automatically the memory previously pointed by inputPhy)
    ImageType::Pointer tempImagePoint = inputPhy;
    inputPhy = outputPhy;
    outputPhy = tempImagePoint;

    if (iter < 9999)
    outputPhy->FillBuffer(0);


#ifdef _DEBUG
    std::cout << "inputPhy reference count : " << inputPhy->GetReferenceCount() << std::endl;
    std::cout << "temp reference count : " << tempImagePoint->GetReferenceCount() << std::endl;
    std::cout << "output reference count : " << outputPhy->GetReferenceCount() << std::endl;
#endif
    }

#ifdef _DEBUG
  std::cout << "end of iteration loop" << std::endl;
  std::cout << "inputPhy reference count : " << inputPhy->GetReferenceCount() << std::endl;
  std::cout << "output reference count : " << outputPhy->GetReferenceCount() << std::endl;
#endif



  // write output image
  typedef unsigned char                            WritePixelType;
  typedef itk::Image< WritePixelType, Dimension >  WriteImageType;
  typedef itk::ImageFileWriter< WriteImageType >   WriterType;

  typedef itk::RescaleIntensityImageFilter<
               ImageType, WriteImageType > RescaleOutputFilterType;

  RescaleOutputFilterType::Pointer rescaler = RescaleOutputFilterType::New();

  rescaler->SetOutputMinimum(   0 );
  rescaler->SetOutputMaximum( 255 );
  rescaler->SetInput(outputPhy);

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
