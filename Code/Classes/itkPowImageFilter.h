/*=========================================================================
 *
 *  Copyright Insight Software Consortium
 *
 *  Licensed under the Apache License, Version 2.0 (the "License");
 *  you may not use this file except in compliance with the License.
 *  You may obtain a copy of the License at
 *
 *         http://www.apache.org/licenses/LICENSE-2.0.txt
 *
 *  Unless required by applicable law or agreed to in writing, software
 *  distributed under the License is distributed on an "AS IS" BASIS,
 *  WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 *  See the License for the specific language governing permissions and
 *  limitations under the License.
 *
 *=========================================================================*/
#ifndef __itkPowImageFilter_h
#define __itkPowImageFilter_h

#include "itkUnaryFunctorImageFilter.h"
#include "vnl/vnl_math.h"

namespace itk
{
/** \class PowImageFilter
 * \brief Computes the function pow(x,n) pixel-wise
 *
 * Every output pixel is equal to vcl_pow(x,n). where x is the intensity of the
 * homologous input pixel, and n is a user-provided constant.
 *
 * \ingroup IntensityImageFilters  Multithreaded
 *
 */

namespace Function
{
template< class TInput, class TOutput >
class PowerFunction
{
public:
  PowerFunction() { m_Power = 2.; }
  ~PowerFunction() {}

  bool operator!=(const PowerFunction & other) const
  {
    if ( m_Power != other.m_Power )
      {
      return true;
      }
    return false;
  }

  bool operator==(const PowerFunction & other) const
  {
    return !( *this != other );
  }

  inline TOutput operator()(const TInput & A) const
  {
    return static_cast< TOutput >( vcl_pow( static_cast< double >( A ), m_Power ) );
  }

  void SetPower(double power)
  {
    m_Power = power;
  }

  double GetPower() const
  {
    return m_Power;
  }

private:
  double m_Power;
};
}


template< class TInputImage, class TOutputImage >
class ITK_EXPORT PowImageFilter:
  public
  UnaryFunctorImageFilter< TInputImage, TOutputImage,
                           Function::PowerFunction<
                             typename TInputImage::PixelType,
                             typename TOutputImage::PixelType >   >
{
public:
  /** Standard class typedefs. */
  typedef PowImageFilter Self;
  typedef UnaryFunctorImageFilter<
    TInputImage, TOutputImage,
    Function::PowerFunction< typename TInputImage::PixelType,
                           typename TOutputImage::PixelType > >  Superclass;

  typedef SmartPointer< Self >       Pointer;
  typedef SmartPointer< const Self > ConstPointer;

  /** Method for creation through the object factory. */
  itkNewMacro(Self);

  /** Runtime information support. */
  itkTypeMacro(PowImageFilter,
               UnaryFunctorImageFilter);

  void SetPower(double power)
  {
    if ( power == this->GetFunctor().GetPower() )
      {
      return;
      }
    this->GetFunctor().SetPower(power);
    this->Modified();
  }

#ifdef ITK_USE_CONCEPT_CHECKING
  /** Begin concept checking */
  itkConceptMacro( InputConvertibleToDoubleCheck,
                   ( Concept::Convertible< typename TInputImage::PixelType, double > ) );
  itkConceptMacro( DoubleConvertibleToOutputCheck,
                   ( Concept::Convertible< double, typename TOutputImage::PixelType > ) );
  /** End concept checking */
#endif
protected:
  PowImageFilter() {}
  virtual ~PowImageFilter() {}
private:
  PowImageFilter(const Self &); //purposely not implemented
  void operator=(const Self &);         //purposely not implemented
};
} // end namespace itk

#endif
