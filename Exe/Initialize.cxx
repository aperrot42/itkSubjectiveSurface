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

  if ( argc < 6 )
    {
    std::cerr << "Missing parameters. " << std::endl;
    std::cerr << "Usage: " << std::endl;
    std::cerr << argv[0]
              << " inputimage"
              << " x y z"
              << " 0-1(peak-distance)"
              << " outputImageFile(image)"
              << std::endl;
    return -1;
    }
  // arguments reading and parsing
  double centerX = atof(argv[2]);
  double centerY = atof(argv[3]);
  double centerZ = atof(argv[4]);
  int func = atoi(argv[5]);


  //%%%%%%%%%%%%%%%%%%% TYPEDEFS %%%%%%%%%%%%%%%%%%%

  // dimension of input image (and phi)
  const   unsigned int   Dimension = 3;

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

  DistanceMetricType::OriginType originPoint( Dimension );
  MeasurementVectorType queryPoint( Dimension );
  originPoint[0] = centerX;
  originPoint[1] = centerY;
  originPoint[2] = centerZ;

  distanceMetric->SetOrigin( originPoint );

  // iterator on initial image
  IteratorType it( initialImage, initialImage->GetRequestedRegion() );



  if (func == 0) // compute peak function
    {
    for ( it.GoToBegin(); !it.IsAtEnd(); ++it )
      {
      ImageType::IndexType idx = it.GetIndex();
      for (unsigned int i = 0; i<Dimension;++i)
        {
        queryPoint[i] = idx[i];
        }
      it.Set( static_cast<PixelType>
              (1./(distanceMetric->Evaluate( queryPoint ) + 0.25 )) );
      }
    }
  else // compute revert distance map
    {
    PixelType max = 0;
    PixelType distance;
    for ( it.GoToBegin(); !it.IsAtEnd(); ++it )
      {
      // compute distance from point
      ImageType::IndexType idx = it.GetIndex();
      for (unsigned int i = 0; i<Dimension;++i)
        {
        queryPoint[i] = idx[i];
        }
      distance = distanceMetric->Evaluate( queryPoint );
      // compute max distance (we could use the fact that it should
      // be located at the farthest angle of the volume but we may have
      // non parallepipedic volumes)
      if (distance > max)
          max = distance;
      // negative distance map
      it.Set( static_cast<PixelType>
              (-(distance)) );
      }
    // get a positive distance map by adding the max distance
    for ( it.GoToBegin(); !it.IsAtEnd(); ++it )
      {
      it.Set( it.Get() + max );
      }
    }


  // write output
  WriterType::Pointer writer = WriterType::New();
  writer->SetFileName( argv[6] );
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
