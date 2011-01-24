#ifndef __itkNucleusEdgeDetectorFilter_h
#define __itkNucleusEdgeDetectorFilter_h


#include "itkImageToImageFilter.h"


#include "itkConceptChecking.h"


#include "itkGradientMagnitudeRecursiveGaussianImageFilter.h"
#include "itkAbsImageFilter.h"
#include "itkPowImageFilter.h"
#include "itkMultiplyByConstantImageFilter.h"
#include "itkAddConstantToImageFilter.h"
#include "itkRescaleIntensityImageFilter.h"

namespace itk
{
/** \class NucleusEdgeDetectorFilter
 *
 * \brief Simple Edge detector for nuclei confocal images
 *
 * This filter produces an output image whose pixels
 * take values according to the function below :
 * g(x,y)  =  1 / (1+( | I(x,y)*derivative(Gauss(Sigma)) | / b )^Pow)
 *
 * \ingroup IntensityImageFilters
 */


template < class TInputImage, class TOutputImage >
class ITK_EXPORT NucleusEdgeDetectorFilter :
  public ImageToImageFilter< TInputImage, TOutputImage >
{
public:
  /** Standard class typedefs. */
  typedef NucleusEdgeDetectorFilter Self;
  typedef ImageToImageFilter< TInputImage, TOutputImage > Superclass;
  typedef SmartPointer< Self >       Pointer;
  typedef SmartPointer< const Self > ConstPointer;

  /** Method for creation through the object factory. */
  itkNewMacro(Self);

  /** Run-time type information (and related methods). */
  itkTypeMacro(NucleusEdgeDetectorFilter, ImageToImageFilter);

  /** Pixel types. */
  typedef typename TInputImage::PixelType  InputPixelType;
  typedef typename TOutputImage::PixelType OutputPixelType;


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
  NucleusEdgeDetectorFilter();
  virtual ~NucleusEdgeDetectorFilter() {}
  void PrintSelf(std::ostream & os, Indent indent) const;


  // take the gradient of the input image
  typedef GradientMagnitudeRecursiveGaussianImageFilter< TInputImage,TOutputImage>
      GradientGaussianFilterType;
  // take abs of the gradient
  typedef AbsImageFilter< TOutputImage,TOutputImage>
      AbsFilterType;
  // multiply by Beta
  typedef MultiplyByConstantImageFilter< TOutputImage, OutputPixelType, TOutputImage>
      MultiplyFilterType;
  // take the image at power Pow
  typedef PowImageFilter< TOutputImage,TOutputImage>
      PowerFilterType;
  // 1 + image
  typedef AddConstantToImageFilter< TOutputImage, OutputPixelType, TOutputImage>
      AddConstantFilterType;
  typedef itk::RescaleIntensityImageFilter< TInputImage,TOutputImage>
      RescaleFilterType;


  typename GradientGaussianFilterType::Pointer    m_GradientGaussianFilter;
  typename AbsFilterType::Pointer         m_AbsFilter;
  typename MultiplyFilterType::Pointer    m_MultiplyFilter;
  typename PowerFilterType::Pointer       m_PowerFilter,
                                          m_InvertFilter;
  typename AddConstantFilterType::Pointer m_AddCstFilter;
  typename RescaleFilterType::Pointer     m_RescaleFilter;
  OutputPixelType m_Beta;
  OutputPixelType m_Sigma;
  int m_Pow;


  void GenerateData();

private:
  NucleusEdgeDetectorFilter(const Self &); //purposely not implemented
  void operator=(const Self &);             //purposely not implemented


};
} // end namespace itk

#ifndef ITK_MANUAL_INSTANTIATION
#include "itkNucleusEdgeDetectorFilter.txx"
#endif

#endif
