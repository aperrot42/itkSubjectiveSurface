
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


#include "itkNucleusEdgeDetectorFilter.h"
#include "itkSubjectiveSurfaceEvolutionFilter.h"

#include "math.h"

int main( int argc, char ** argv )
{

  if ( argc < 8 )
    {
    std::cerr << "Missing parameters. " << std::endl;
    std::cerr << "Usage: " << std::endl;
    std::cerr << argv[0]
              << " inputphi(image) inputImage(image)"
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




  // subjective surfaces and edge detection
  typedef itk::SubjectiveSurfaceEvolutionFilter<ImageType,ImageType>
      SubjectiveSurfaceEvolutionFilterType;
  typedef itk::NucleusEdgeDetectorFilter<ImageType,ImageType>
      EdgeFilterType;
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


  ReaderType::Pointer readerImage = ReaderType::New();
  readerImage->SetFileName( argv[2] );
  try
    {
    readerImage->Update();
    }
  catch ( itk::ExceptionObject &err)
    {
    std::cout << "ExceptionObject caught !" << std::endl;
    std::cout << err << std::endl;
    return -1;
    }




  // membrane edge filter

  EdgeFilterType::Pointer edgeFilter = EdgeFilterType::New();
  edgeFilter->SetInput(readerImage->GetOutput());
  edgeFilter->SetBeta((double)255);
  edgeFilter->SetSigma((double)0.2);
  edgeFilter->SetPow((double)1);
  edgeFilter->Update();


  // subjective surfaces
  SubjectiveSurfaceEvolutionFilterType::Pointer SubjSurfFilter =
      SubjectiveSurfaceEvolutionFilterType::New();
  SubjSurfFilter->SetInput(inputphi);
  SubjSurfFilter->SetEdgeMap(edgeFilter->GetOutput());

  SubjSurfFilter->SetNumberIteration( nbiter );
  SubjSurfFilter->SetDeltaT( deltat );
  SubjSurfFilter->Seta( curvatureFactor );
  SubjSurfFilter->SetNu( nu );
  SubjSurfFilter->SetRho( rho );


  SubjSurfFilter->Update();


// rescaling output (phi(t=tmax))
RescaleOutputFilterType::Pointer rescaler = RescaleOutputFilterType::New();
rescaler->SetOutputMinimum(   0 );
rescaler->SetOutputMaximum( 255 );
rescaler->SetInput(edgeFilter->GetOutput());

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
