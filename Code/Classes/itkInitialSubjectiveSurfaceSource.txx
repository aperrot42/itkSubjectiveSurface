#ifndef __itkInitialSubjectiveSurfaceSource_txx
#define __itkInitialSubjectiveSurfaceSource_txx

#include "itkInitialSubjectiveSurfaceSource.h"
#include "itkEuclideanDistanceMetric.h"
#include "itkImageRegionIterator.h"
#include "itkProgressReporter.h"
#include "itkObjectFactory.h"



namespace itk
{
/**
 *
 */
template< class TOutputImage >
InitialSubjectiveSurfaceSource<TOutputImage>::
InitialSubjectiveSurfaceSource()
{
  m_Size.Fill(64);
  m_Spacing.Fill(1.0);
  m_Origin.Fill(0.0);
  m_Direction.SetIdentity();

  m_Function = 0; // defaut function is peak
  m_FocusPoint.Fill(0);
}



template< class TOutputImage >
InitialSubjectiveSurfaceSource<TOutputImage>::
~InitialSubjectiveSurfaceSource()
{
}



template< class TOutputImage >
void
InitialSubjectiveSurfaceSource<TOutputImage>::
GenerateData()
{
  unsigned int i;
  typename TOutputImage::Pointer outputPtr = this->GetOutput();

  // allocate the output buffer
  outputPtr->SetBufferedRegion( outputPtr->GetRequestedRegion() );
  outputPtr->Allocate();


  //distance metric

  // Create and initialize a new distance function
  typedef PointType MeasurementVectorType;
  typedef typename itk::Statistics::EuclideanDistanceMetric<
                                        MeasurementVectorType >
                                          DistanceMetricType;
  typename DistanceMetricType::Pointer m_DistanceMetric =
      DistanceMetricType::New();
  typename DistanceMetricType::OriginType m_DistanceOrigin(
      NDimensions);
  for (i=0; i<NDimensions; ++i)
    {
    m_DistanceOrigin[i] = m_FocusPoint[i];
    }

  m_DistanceMetric->SetOrigin(m_DistanceOrigin);


  // Create an iterator that will walk the output region
  typedef ImageRegionIterator<TOutputImage> OutputIterator;
  OutputIterator outIt = OutputIterator(outputPtr,
                                        outputPtr->GetRequestedRegion());

  // The value produced by the spatial function
  double value;

  // The position at which the function is evaluated
  Point<double, TOutputImage::ImageDimension> evalPoint;



  if (m_Function == 0) // compute peak function
    {

    ProgressReporter progress(this, 0,
                              outputPtr->GetRequestedRegion()
                                           .GetNumberOfPixels());

    // Walk the output image, evaluating the spatial function at each pixel
    for (outIt.GoToBegin(); !outIt.IsAtEnd(); ++outIt)
      {
      typename TOutputImage::IndexType index = outIt.GetIndex();
      outputPtr->TransformIndexToPhysicalPoint(index, evalPoint );
      value = m_DistanceMetric->Evaluate(evalPoint);
      outIt.Set( static_cast<OutputImagePixelType>(1./(value + 0.25) ));
      progress.CompletedPixel();
      }
  }
  else // compute revert distance map
    {
    OutputImagePixelType max = 0;
    ProgressReporter progress(this, 0,
                              (outputPtr->GetRequestedRegion()
                                           .GetNumberOfPixels())*2);
    for (outIt.GoToBegin(); !outIt.IsAtEnd(); ++outIt)
      {

      typename TOutputImage::IndexType index = outIt.GetIndex();
      outputPtr->TransformIndexToPhysicalPoint(index, evalPoint );
      value = m_DistanceMetric->Evaluate(evalPoint);

      // compute max distance (we could use the fact that it should
      // be located at the farthest angle of the volume but we may have
      // non parallepipedic volumes)
      if (value > max)
          max = value;
      // negative distance map
      outIt.Set( static_cast<OutputImagePixelType>(-(value)) );
      progress.CompletedPixel();
      }
    // get a positive distance map by adding the max distance
    for (outIt.GoToBegin(); !outIt.IsAtEnd(); ++outIt)
      {
      outIt.Set( outIt.Get() + max );
      progress.CompletedPixel();
      }
    }
}




template< class TOutputImage >
void
InitialSubjectiveSurfaceSource<TOutputImage>::
GenerateOutputInformation()
{
  TOutputImage *output;

  typename TOutputImage::IndexType index = { { 0 } };

  output = this->GetOutput(0);

  typename TOutputImage::RegionType largestPossibleRegion;
  largestPossibleRegion.SetSize(m_Size);
  largestPossibleRegion.SetIndex(index);
  output->SetLargestPossibleRegion(largestPossibleRegion);

  output->SetSpacing(m_Spacing);
  output->SetOrigin(m_Origin);
  output->SetDirection(m_Direction);
}




/**
 *
 */
template< class TOutputImage >
void
InitialSubjectiveSurfaceSource<TOutputImage>::
PrintSelf(std::ostream & os, Indent indent) const
{
  Superclass::PrintSelf(os, indent);

    unsigned int i;

    os << indent << "Size: [";
    for (i=0; i < NDimensions - 1; i++)
      {
      os << m_Size[i] << ", ";
      }
    os << "]" << std::endl;

    os << indent << "Origin: [";
    for (i=0; i < NDimensions - 1; i++)
      {
      os << m_Origin[i] << ", ";
      }
    os << "]" << std::endl;

    os << indent << "Spacing: " << m_Spacing << std::endl;

    os << indent << "Direction:" << std::endl;
    os << m_Direction << std::endl;

}

}

#endif
