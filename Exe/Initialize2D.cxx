#if defined(_MSC_VER)
#pragma warning ( disable : 4786 )
#endif


#include "itkImageFileReader.h"
#include "itkImageFileWriter.h"
#include "itkRescaleIntensityImageFilter.h"

#include "itkVector.h"
#include "itkArray.h"
#include "itkEuclideanDistanceMetric.h"



#include "itkImageRegionIteratorWithIndex.h"
#include "itkEuclideanDistanceMetric.h"

#include "math.h"


int main( int argc, char ** argv )
{


  //%%%%%%%%%%%%%%%%%%% ARGUMENTS PARSING %%%%%%%%%%%%%%%%%%%

  if ( argc < 5 )
    {
    std::cerr << "Missing parameters. " << std::endl;
    std::cerr << "Usage: " << std::endl;
    std::cerr << argv[0]
              << " inputimage"
              << " x y"
              << " 0-1(peak-distance)"
              << " outputImageFile(image)"
              << std::endl;
    return -1;
    }
  // arguments reading and parsing
  double centerX = atof(argv[2]);
  double centerY = atof(argv[3]);
  int func = atoi(argv[4]);


  //%%%%%%%%%%%%%%%%%%% TYPEDEFS %%%%%%%%%%%%%%%%%%%

  // dimension of input image (and phi)
  const   unsigned int   Dimension = 2;

  // scalar field type
  typedef double                              PixelType;
  typedef itk::Image< PixelType, Dimension > ImageType;

  // reader
  typedef itk::ImageFileReader< ImageType >  ReaderType;

  // writer
  typedef itk::ImageFileWriter< ImageType >   WriterType;

  // simple iterator
  typedef itk::ImageRegionIteratorWithIndex< ImageType > IteratorType;




  //%%%%%%%%%%%%%%%%%%% INPUT READING %%%%%%%%%%%%%%%%%%%

  // reading filtered input
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

  // create empty image
  ImageType::Pointer initialImage = ImageType::New();
  initialImage->SetRegions(reader->GetOutput()->GetLargestPossibleRegion());
  initialImage->SetSpacing(reader->GetOutput()->GetSpacing());
  initialImage->Allocate();
  initialImage->FillBuffer(static_cast<PixelType>(0.));


  //n distance metric initialization
  typedef itk::Array< PixelType > MeasurementVectorType;
  typedef itk::Statistics::EuclideanDistanceMetric< MeasurementVectorType >
  DistanceMetricType;

  DistanceMetricType::Pointer distanceMetric = DistanceMetricType::New();

  DistanceMetricType::OriginType originPoint( 2 );
  MeasurementVectorType queryPoint( 2 );
  originPoint[0] = centerX;
  originPoint[1] = centerY;

  distanceMetric->SetOrigin( originPoint );

  // iterator on initial image
  IteratorType it( initialImage, initialImage->GetRequestedRegion() );

  // compute distance at each point
  for ( it.GoToBegin(); !it.IsAtEnd(); ++it )
    {
    ImageType::IndexType idx = it.GetIndex();
    queryPoint[0] = idx[0];
    queryPoint[1] = idx[1];

    it.Set( distanceMetric->Evaluate( queryPoint ));
    }


  // write output
  WriterType::Pointer writer = WriterType::New();
  writer->SetFileName( argv[5] );
  writer->SetInput(initialImage);
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
