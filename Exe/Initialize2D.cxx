#if defined(_MSC_VER)
#pragma warning ( disable : 4786 )
#endif


#include "itkImageFileReader.h"
#include "itkImageFileWriter.h"
#include "itkRescaleIntensityImageFilter.h"
#include "itkInitialSubjectiveSurfaceSource.h"


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
  typedef itk::InitialSubjectiveSurfaceSource< ImageType > SurfaceSourceType;




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


  SurfaceSourceType::Pointer surfaceSource = SurfaceSourceType::New();
  surfaceSource->SetSize(reader->GetOutput()->GetLargestPossibleRegion().GetSize());
  surfaceSource->SetSpacing(reader->GetOutput()->GetSpacing());
  surfaceSource->SetOrigin(reader->GetOutput()->GetOrigin());

  ImageType::PointType focusPoint;
  focusPoint[0]= centerX;
  focusPoint[1]= centerY;
  //FocusPoint[2]= 0;


  surfaceSource->SetFunction( func ); //peak function
  surfaceSource->SetFocusPoint(focusPoint);
  // Run the pipeline
  surfaceSource->Update();

  // write output
  WriterType::Pointer writer = WriterType::New();
  writer->SetFileName( argv[5] );
  writer->SetInput(surfaceSource->GetOutput());
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
