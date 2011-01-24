#if defined(_MSC_VER)
#pragma warning ( disable : 4786 )
#endif


#include "itkImageFileReader.h"
#include "itkImageFileWriter.h"
#include "itkRescaleIntensityImageFilter.h"

#include "math.h"

#include "itkNucleusEdgeDetectorFilter.h"
#include "itkMembraneEdgeDetectorFilter.h"

int main( int argc, char ** argv )
{


  //%%%%%%%%%%%%%%%%%%% ARGUMENTS PARSING %%%%%%%%%%%%%%%%%%%

  if ( argc < 6 )
    {
    std::cerr << "Missing parameters. " << std::endl;
    std::cerr << "Usage: " << std::endl;
    std::cerr << argv[0]
              << " inputFilteredMembrane(image)"
              << " sigma b n"
              << " detector type(0:nucleus 1:membrane)"
              << " outputImageFile(image)"
              << std::endl;
    return -1;
    }
  // arguments reading and parsing
  double sigma = atof(argv[2]);
  double b = atof(argv[3]);
  double n = atof(argv[4]);
  int detectortype = atoi(argv[5]);

  // for the membrane, the formula of the edge detector is :
  // g(x,y)  =  1 / (1+(|gauss(I(x,y),sigma)|/b)^n)



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

  typedef itk::NucleusEdgeDetectorFilter<ImageType,ImageType>
      NucleusEdgeFilterType;
  typedef itk::MembraneEdgeDetectorFilter<ImageType,ImageType>
      MembraneEdgeFilterType;


  //%%%%%%%%%%%%%%%%%%% INPUT READING %%%%%%%%%%%%%%%%%%%

  // reading input
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


  // edge filters
  NucleusEdgeFilterType::Pointer nucEdgeFilter = NucleusEdgeFilterType::New();
  MembraneEdgeFilterType::Pointer meEdgeFilter = MembraneEdgeFilterType::New();
  // writer
  WriterType::Pointer writer = WriterType::New();
  if (detectortype == 0)
    {
    nucEdgeFilter->SetInput(reader->GetOutput());
    nucEdgeFilter->SetBeta((double)b);
    nucEdgeFilter->SetSigma((double)sigma);
    nucEdgeFilter->SetPow((double)n);
    nucEdgeFilter->Update();
    writer->SetInput(nucEdgeFilter->GetOutput());
    }
  else
    {
    meEdgeFilter->SetInput(reader->GetOutput());
    meEdgeFilter->SetBeta((double)b);
    meEdgeFilter->SetSigma((double)sigma);
    meEdgeFilter->SetPow((double)n);
    meEdgeFilter->Update();
    writer->SetInput(meEdgeFilter->GetOutput());
    }

  // write output
  writer->SetFileName( argv[6] );
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
