#ifndef __itkMembraneEdgeDetectorFilter_h
#define __itkMembraneEdgeDetectorFilter_h


#include "itkImageToImageFilter.h"


#include "itkConceptChecking.h"


#include "itkRecursiveGaussianImageFilter.h"
#include "itkAbsImageFilter.h"
#include "itkPowImageFilter.h"
#include "itkMultiplyByConstantImageFilter.h"
#include "itkAddConstantToImageFilter.h"


namespace itk
{
/** \class MembraneEdgeDetectorFilter
 *
 * \brief Simple Edge detector for membrane confocal images
 *
 * This filter produces an output image whose pixels
 * take values according to the function below :
 * g(x,y)  =  1 / (1+( | I(x,y)*Gauss(Sigma) | / b )^Pow)
 *
 * \ingroup IntensityImageFilters
 */


template < class TInputImage, class TOutputImage >
class ITK_EXPORT MembraneEdgeDetectorFilter :
  public ImageToImageFilter< TInputImage, TOutputImage >
{
public:
  /** Standard class typedefs. */
  typedef MembraneEdgeDetectorFilter Self;
  typedef ImageToImageFilter< TInputImage, TOutputImage > Superclass;
  typedef SmartPointer< Self >       Pointer;
  typedef SmartPointer< const Self > ConstPointer;

  /** Method for creation through the object factory. */
  itkNewMacro(Self);

  /** Run-time type information (and related methods). */
  itkTypeMacro(MembraneEdgeDetectorFilter, ImageToImageFilter);

  /** Pixel types. */
  typedef typename TInputImage::PixelType  InputPixelType;
  typedef typename TOutputImage::PixelType OutputPixelType;



  /** Set the "outside" pixel value. The default value
   * NumericTraits<OutputPixelType>::Zero. */
  itkSetMacro(OutsideValue, OutputPixelType);

  /** Get Sigma */
  itkGetConstReferenceMacro(Sigma, OutputPixelType);
  /** Set Sigma */
  itkSetMacro(Sigma, OutputPixelType);

  /** Get Beta */
  itkGetConstReferenceMacro(Beta, OutputPixelType);
  /** Set Beta */
  itkSetMacro(Beta, OutputPixelType);

  /** Get Pow */
  itkGetConstReferenceMacro(Pow, int);
  /** Set Pow */
  itkSetMacro(Pow, int);


#ifdef ITK_USE_CONCEPT_CHECKING
  /** Begin concept checking */
  itkConceptMacro( InputOStreamWritableCheck,
                   ( Concept::OStreamWritable< InputPixelType > ) );
  itkConceptMacro( OutputOStreamWritableCheck,
                   ( Concept::OStreamWritable< OutputPixelType > ) );
  /** End concept checking */
#endif
protected:
  MembraneEdgeDetectorFilter();
  virtual ~MembraneEdgeDetectorFilter() {}
  void PrintSelf(std::ostream & os, Indent indent) const;


  // take the gradient of the input image
  typedef itk::RecursiveGaussianImageFilter< TInputImage,TOutputImage>
      GaussianFilterType;
  // take abs of the gradient
  typedef itk::AbsImageFilter< TOutputImage,TOutputImage>
      AbsFilterType;
  // multiply by Beta
  typedef itk::MultiplyByConstantImageFilter< TOutputImage, OutputPixelType, TOutputImage>
      MultiplyFilterType;
  // take the image at power Pow
  typedef itk::PowImageFilter< TOutputImage,TOutputImage>
      PowerFilterType;
  // 1 + image
  typedef itk::AddConstantToImageFilter< TOutputImage, OutputPixelType, TOutputImage>
      AddConstantFilterType;


  typename GaussianFilterType::Pointer    m_GaussianFilterX,
                                          m_GaussianFilterY,
                                          m_GaussianFilterZ;
  typename AbsFilterType::Pointer         m_AbsFilter;
  typename MultiplyFilterType::Pointer    m_MultiplyFilter;
  typename PowerFilterType::Pointer       m_PowerFilter,
                                          m_InvertFilter;
  typename AddConstantFilterType::Pointer m_AddCstFilter;
  OutputPixelType m_Beta;
  OutputPixelType m_Sigma;
  int m_Pow;


  void GenerateData();

private:
  MembraneEdgeDetectorFilter(const Self &); //purposely not implemented
  void operator=(const Self &);             //purposely not implemented


};
} // end namespace itk

#ifndef ITK_MANUAL_INSTANTIATION
#include "itkMembraneEdgeDetectorFilter.txx"
#endif

#endif
