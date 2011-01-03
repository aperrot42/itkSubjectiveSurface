#if defined(_MSC_VER)
#pragma warning ( disable : 4786 )
#endif


#include "itkImageFileReader.h"
#include "itkImageFileWriter.h"
#include "itkRescaleIntensityImageFilter.h"


#include "itkRecursiveGaussianImageFilter.h"
#include "itkAbsImageFilter.h"
#include "itkPowImageFilter.h"
#include "itkMultiplyByConstantImageFilter.h"
#include "itkAddConstantToImageFilter.h"

#include "math.h"


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
              << " outputImageFile(image)"
              << std::endl;
    return -1;
    }
  // arguments reading and parsing
  double sigma = atof(argv[2]);
  double b = atof(argv[3]);
  double n = atof(argv[4]);

  // for the membrane, the formula of the edge detector is :
  // g(x,y)  =  1 / (1+(|g(I(x,y),sigma)|/b)^n)



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

  // take the gradient of the input image
  typedef itk::RecursiveGaussianImageFilter< ImageType,ImageType>
          GaussianFilterType;

  // take abs of the gradient
  typedef itk::AbsImageFilter< ImageType,ImageType>
      AbsFilterType;

  // multiply by b
  typedef itk::MultiplyByConstantImageFilter< ImageType, PixelType, ImageType >
      MultiplyFilterType;

  // take the image at power n
  typedef itk::PowImageFilter< ImageType, ImageType >
      PowerFilterType;

  // 1 + image
  typedef itk::AddConstantToImageFilter< ImageType, PixelType, ImageType >
      AddConstantFilterType;



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


  //%%%%%%%%%%%%%%%%%%% EDGE DETECTION %%%%%%%%%%%%%%%%%%%
  GaussianFilterType::Pointer gaussianX = GaussianFilterType::New();
  gaussianX->SetInput(reader->GetOutput());
  gaussianX->SetSigma(sigma);
  gaussianX->SetDirection(0);
  gaussianX->SetOrder( GaussianFilterType::ZeroOrder );

  GaussianFilterType::Pointer gaussianY = GaussianFilterType::New();
  gaussianY->SetInput(gaussianX->GetOutput());
  gaussianY->SetSigma(sigma);
  gaussianY->SetDirection(1);
  gaussianY->SetOrder( GaussianFilterType::ZeroOrder );

  AbsFilterType::Pointer absolute = AbsFilterType::New();
  absolute->SetInput(gaussianY->GetOutput());

  MultiplyFilterType::Pointer multiply =  MultiplyFilterType::New();
  multiply->SetInput(absolute->GetOutput());
  multiply->SetConstant(1/b);

  PowerFilterType::Pointer power = PowerFilterType::New();
  power->SetInput(multiply->GetOutput());
  power->SetPower(n);

  AddConstantFilterType::Pointer addConstant = AddConstantFilterType::New();
  addConstant->SetInput(power->GetOutput());
  addConstant->SetConstant(1);

  addConstant->Update();

  PowerFilterType::Pointer inverse = PowerFilterType::New();
  inverse->SetInput(addConstant->GetOutput());
  inverse->SetPower(-1);


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

  // write output
  WriterType::Pointer writer = WriterType::New();
  writer->SetFileName( argv[5] );
  writer->SetInput(inverse->GetOutput());
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
