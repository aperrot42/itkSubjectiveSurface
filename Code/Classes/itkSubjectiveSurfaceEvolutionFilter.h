#ifndef __itkSubjectiveSurfaceEvolutionFilter_h
#define __itkSubjectiveSurfaceEvolutionFilter_h


#include "itkImageToImageFilter.h"


#include "itkConceptChecking.h"


#include "itkConstantBoundaryCondition.h"

#include "itkGradientRecursiveGaussianImageFilter.h"
#include "itkImageDuplicator.h"

namespace itk
{
/** \class SubjectiveSurfaceEvolutionFilter
 *
 * \brief Simple Edge detector for membrane confocal images
 *
 * This filter produces an output image whose pixels
 * take values according to the function below :
 * g(x,y)  =  1 / (1+( | I(x,y)*Gauss(Sigma) | / b )^Pow)
 *
 * \ingroup IntensityImageFilters
 */


template < class InputImageType, class OutputImageType >
class ITK_EXPORT SubjectiveSurfaceEvolutionFilter :
  public ImageToImageFilter< InputImageType, OutputImageType >
{
public:
  /** Standard class typedefs. */
  typedef SubjectiveSurfaceEvolutionFilter Self;
  typedef ImageToImageFilter< InputImageType, OutputImageType > Superclass;
  typedef SmartPointer< Self >       Pointer;
  typedef SmartPointer< const Self > ConstPointer;

  /** Method for creation through the object factory. */
  itkNewMacro(Self);

  /** Run-time type information (and related methods). */
  itkTypeMacro(SubjectiveSurfaceEvolutionFilter, ImageToImageFilter);


  /** Pixel types. */
  typedef typename InputImageType::PixelType  InputPixelType;
  typedef typename OutputImageType::PixelType OutputPixelType;
  typedef itk::CovariantVector< OutputPixelType, InputImageType::ImageDimension >
      VectorPixelType;

  /** Image types. */
  typedef itk::Image< VectorPixelType, InputImageType::ImageDimension  >
      VectorImageType;


  /** Iterators types */
  // neighborhood iterators
  typedef itk::ConstantBoundaryCondition< InputImageType >  BoundaryConditionType;
  typedef itk::ConstNeighborhoodIterator< InputImageType, BoundaryConditionType >
          NeighborhoodIteratorType;
  // simple iterators
  typedef itk::ImageRegionIterator< InputImageType>        IteratorType;
  typedef itk::ImageRegionConstIterator< InputImageType > ConstIteratorType;
  typedef itk::ImageRegionConstIterator< VectorImageType >
          ConstVectorIteratorType;



  void SetInput( const InputImageType * image1 );

  void SetEdgeMap( const InputImageType * image1 );


  /** Get Number of Iterations */
  itkGetConstReferenceMacro(NumberIteration, unsigned int);
  /** Set Number of Iterations */
  itkSetMacro(NumberIteration, unsigned int);

  /** Get time interval between two iterations (DeltaT) */
  itkGetConstReferenceMacro(DeltaT, unsigned int);
  /** Set time interval between two iterations (DeltaT) */
  itkSetMacro(DeltaT, unsigned int);

  /** Get coefficient of diffusive term (nu) */
  itkGetConstReferenceMacro(Nu, InputPixelType);
  /** Set coefficient of diffusive term (nu) */
  itkSetMacro(Nu, InputPixelType);

  /** Get coefficient of Level Set term (Rho) */
  itkGetConstReferenceMacro(Rho, InputPixelType);
  /** Set coefficient of Level Set term (Rho) */
  itkSetMacro(Rho, InputPixelType);


#ifdef ITK_USE_CONCEPT_CHECKING
  /** Begin concept checking */

  itkConceptMacro( InputOStreamWritableCheck,
                   ( Concept::OStreamWritable< InputPixelType > ) );
  itkConceptMacro( OutputOStreamWritableCheck,
                   ( Concept::OStreamWritable< OutputPixelType > ) );
  /** End concept checking */
#endif


protected:
  SubjectiveSurfaceEvolutionFilter();
  virtual ~SubjectiveSurfaceEvolutionFilter() {}
  void PrintSelf(std::ostream & os, Indent indent) const;


  /** Filters typedef */

  /**  Type of the Gradient Filter. */
  typedef itk::GradientRecursiveGaussianImageFilter< InputImageType, VectorImageType >
      GradientFilterType;
  /**  Type of the Duplicator Filter. */
  typedef itk::ImageDuplicator< InputImageType >
      ImageDuplicatorType;

  typename GradientFilterType::Pointer        m_ImageDuplicator;
  typename ImageDuplicatorType::Pointer       m_Gradient;
  typename InputImageType::Pointer            m_InputImage;
  typename InputImageType::Pointer            m_OutputImage;
  typename InputImageType::Pointer            m_TempImage;
  typename InputImageType::Pointer            m_EdgeDetector;

  int                    m_NumIter;
  int                    m_DeltaT;
  InputPixelType    m_Nu;
  InputPixelType    m_Rho;

  void GenerateData();


private:
  SubjectiveSurfaceEvolutionFilter(const Self &); //purposely not implemented
  void operator=(const Self &);             //purposely not implemented


};
} // end namespace itk

#ifndef ITK_MANUAL_INSTANTIATION
#include "itkSubjectiveSurfaceEvolutionFilter.txx"
#endif

#endif










