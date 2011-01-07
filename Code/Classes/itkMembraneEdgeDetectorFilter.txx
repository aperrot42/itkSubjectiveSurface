#ifndef __itkMembraneEdgeDetectorFilter_txx
#define __itkMembraneEdgeDetectorFilter_txx

#include "itkMembraneEdgeDetectorFilter.h"

namespace itk
{
/**
 *
 */
template< class TInputImage, class TOutputImage >
MembraneEdgeDetectorFilter< TInputImage, TOutputImage >
::MembraneEdgeDetectorFilter()
{
  m_Beta   = 1;
  m_Sigma  = 1;
  m_Pow    = 1;

  m_GaussianFilterX = GaussianFilterType::New();
  m_GaussianFilterY = GaussianFilterType::New();
  m_GaussianFilterZ = GaussianFilterType::New();

  m_AbsFilter = AbsFilterType::New();

  m_MultiplyFilter = MultiplyFilterType::New();

  m_PowerFilter = PowerFilterType::New();

  m_AddCstFilter = AddConstantFilterType::New();
  m_AddCstFilter->SetConstant(1);

  m_InvertFilter = PowerFilterType::New();
  m_InvertFilter->SetPower(-1);
}



template< class TInputImage, class TOutputImage >
void
MembraneEdgeDetectorFilter<TInputImage, TOutputImage>::
GenerateData()
{

if (TInputImage::Dimension >= 1)
  {
  m_GaussianFilterX->SetInput( this->GetInput() );
  m_GaussianFilterX->SetSigma(m_Sigma);
  m_GaussianFilterX->SetDirection(0);
  m_GaussianFilterX->SetOrder( GaussianFilterType::ZeroOrder );

  m_AbsFilter->SetInput(m_GaussianFilterX->GetOutput());
  }

if (TInputImage::Dimension >= 2)
  {
  m_GaussianFilterY->SetInput(m_GaussianFilterX->GetOutput());
  m_GaussianFilterY->SetSigma(m_Sigma);
  m_GaussianFilterY->SetDirection(1);
  m_GaussianFilterY->SetOrder( GaussianFilterType::ZeroOrder );

  m_AbsFilter->SetInput(m_GaussianFilterY->GetOutput());
  }


if (TInputImage::Dimension >= 3)
  {
  m_GaussianFilterZ->SetInput( m_GaussianFilterY->GetOutput() );
  m_GaussianFilterZ->SetSigma( m_Sigma );
  m_GaussianFilterZ->SetDirection(2);
  m_GaussianFilterZ->SetOrder( GaussianFilterType::ZeroOrder );

  m_AbsFilter->SetInput( m_GaussianFilterZ->GetOutput() );
  }


  m_MultiplyFilter->SetInput( m_AbsFilter->GetOutput() );
  m_MultiplyFilter->SetConstant( 1/m_Beta );

  m_PowerFilter->SetInput( m_MultiplyFilter->GetOutput() );
  m_PowerFilter->SetPower( m_Pow);

  m_AddCstFilter->SetInput( m_PowerFilter->GetOutput() );
  m_AddCstFilter->Update();

  m_InvertFilter->SetInput( m_AddCstFilter->GetOutput() );
  m_InvertFilter->Update();


  this->GraftOutput( m_InvertFilter->GetOutput() );

}



/**
 *
 */
template< class TInputImage, class TOutputImage >
void
MembraneEdgeDetectorFilter<TInputImage, TOutputImage>::
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
