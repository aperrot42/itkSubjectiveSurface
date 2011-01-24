#ifndef __itkNucleusEdgeDetectorFilter_txx
#define __itkNucleusEdgeDetectorFilter_txx

#include "itkNucleusEdgeDetectorFilter.h"

namespace itk
{
/**
 *
 */
template< class TInputImage, class TOutputImage >
NucleusEdgeDetectorFilter< TInputImage, TOutputImage >
::NucleusEdgeDetectorFilter()
{
  m_Beta   = 1;
  m_Sigma  = 1;
  m_Pow    = 1;

  m_GradientGaussianFilter = GradientGaussianFilterType::New();

  m_AbsFilter = AbsFilterType::New();

  m_MultiplyFilter = MultiplyFilterType::New();

  m_InvertFilter = PowerFilterType::New();
  m_InvertFilter->SetPower(-1);

  m_RescaleFilter = RescaleFilterType::New();

  m_AddCstFilter = AddConstantFilterType::New();
}



template< class TInputImage, class TOutputImage >
void
NucleusEdgeDetectorFilter<TInputImage, TOutputImage>::
GenerateData()
{


  m_GradientGaussianFilter->SetInput( this->GetInput() );
  m_GradientGaussianFilter->SetSigma(m_Sigma);

  m_AbsFilter->SetInput(m_GradientGaussianFilter->GetOutput());

  m_RescaleFilter->SetInput(m_AbsFilter->GetOutput());
  m_RescaleFilter->SetOutputMinimum( 0. );
  m_RescaleFilter->SetOutputMaximum( 1. );

  m_MultiplyFilter->SetInput( m_RescaleFilter->GetOutput() );
  m_MultiplyFilter->SetConstant( -1. );

  m_AddCstFilter->SetInput(m_MultiplyFilter->GetOutput());
  m_AddCstFilter->SetConstant( 1. );
  m_AddCstFilter->Update();

  this->GraftOutput( m_AddCstFilter->GetOutput() );

}



/**
 *
 */
template< class TInputImage, class TOutputImage >
void
NucleusEdgeDetectorFilter<TInputImage, TOutputImage>::
PrintSelf(std::ostream & os, Indent indent) const
{
  Superclass::PrintSelf(os, indent);

  os << indent << "Beta: "
     << static_cast< typename NumericTraits< OutputPixelType >::PrintType >( m_Beta ) << std::endl;
  os << indent << "Sigma: "
     << static_cast< typename NumericTraits< OutputPixelType >::PrintType >( m_Sigma ) << std::endl;
  os << indent << "Pow: "
     << static_cast< int >( this->GetPow() ) << std::endl;
}

}
#endif
