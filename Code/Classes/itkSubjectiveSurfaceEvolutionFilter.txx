#ifndef __itkSubjectiveSurfaceEvolutionFilter_txx
#define __itkSubjectiveSurfaceEvolutionFilter_txx

#include "itkSubjectiveSurfaceEvolutionFilter.h"

namespace itk
{
/**
 *
 */
template< class InputImageType, class OutputImageType >
SubjectiveSurfaceEvolutionFilter< InputImageType, OutputImageType >
::SubjectiveSurfaceEvolutionFilter()
{
  // we need the input phi and the feature map
  this->SetNumberOfRequiredInputs( 2 );


  m_ImageDuplicator = ImageDuplicatorType::New();
  // set phi as input of the duplicator
  m_ImageDuplicator->SetInput(this->GetInput(0));

  m_Gradient = GradientFilterType::New();
}

template< class InputImageType, class OutputImageType >
void
SubjectiveSurfaceEvolutionFilter<InputImageType, OutputImageType>::
SetInput( const InputImageType * image1 )
{
// Process object is not const-correct
// so the const casting is required.
SetNthInput(0, const_cast
<InputImageType *>( image1 ));
}

template< class InputImageType, class OutputImageType >
void
SubjectiveSurfaceEvolutionFilter<InputImageType, OutputImageType>::
SetEdgeMap( const InputImageType * image1 )
{
// Process object is not const-correct
// so the const casting is required.
SetNthInput(1, const_cast
<InputImageType *>( image1 ));
}


template< class InputImageType, class OutputImageType >
void
SubjectiveSurfaceEvolutionFilter<InputImageType, OutputImageType>::
GenerateData()
{
//if (InputImageType::ImageDimension >= 1)
//  {
//  m_GaussianFilterX->SetInput( this->GetInput() );
//  m_GaussianFilterX->SetSigma(m_Sigma);
//  m_GaussianFilterX->SetDirection(0);
//  m_GaussianFilterX->SetOrder( GaussianFilterType::ZeroOrder );

//  m_AbsFilter->SetInput(m_GaussianFilterX->GetOutput());
//  }

//if (InputImageType::ImageDimension >= 2)
//  {
//  m_GaussianFilterY->SetInput(m_GaussianFilterX->GetOutput());
//  m_GaussianFilterY->SetSigma(m_Sigma);
//  m_GaussianFilterY->SetDirection(1);
//  m_GaussianFilterY->SetOrder( GaussianFilterType::ZeroOrder );

//  m_AbsFilter->SetInput(m_GaussianFilterY->GetOutput());
//  }


if (InputImageType::ImageDimension >= 3)
  {
  // duplicate phi for const correctness
  m_ImageDuplicator->Update();
  m_InputImage = m_ImageDuplicator->GetOutput();

  // allocation of phi(t+1)
  m_OutputImage = ImageType::New();
  m_OutputImage->SetRegions(m_InputImage->GetRequestedRegion());
  m_OutputImage->SetSpacing(m_InputImage->GetSpacing());
  m_OutputImage->Allocate();



  // pointer to edgemap (g)
  m_EdgeDetector = this->GetInput(1);
  // iterator on g (scalar)
  ConstIteratorType edgit( m_EdgeDetector, m_EdgeDetector->GetRequestedRegion() );


  // computation of grad(g) (DimensionxDimension)
  m_Gradient->SetInput( m_EdgeDetector );
  m_Gradient->Update();
  // iterator on grad(g) (vector)
  ConstVectorIteratorType edgderivit( m_Gradient->GetOutput(),
                                      m_Gradient->GetOutput()
                                              ->GetLargestPossibleRegion() );
  // gradient vector at a given point
  VectorPixelType edgederivVector;


  // offsets definition :
  NeighborhoodIteratorType::OffsetType minusx = {{-1,0,0}};
  NeighborhoodIteratorType::OffsetType plusx = {{1,0,0}};

  NeighborhoodIteratorType::OffsetType minusy = {{0,-1,0}};
  NeighborhoodIteratorType::OffsetType plusy = {{0,1,0}};

  NeighborhoodIteratorType::OffsetType minusz = {{0,0,-1}};
  NeighborhoodIteratorType::OffsetType plusz = {{0,0,1}};


  NeighborhoodIteratorType::OffsetType plusxplusy = {{1,1,0}};
  NeighborhoodIteratorType::OffsetType minusxminusy = {{-1,-1,0}};
  NeighborhoodIteratorType::OffsetType plusxminusy = {{1,-1,0}};
  NeighborhoodIteratorType::OffsetType minusxplusy = {{-1,1,0}};

  NeighborhoodIteratorType::OffsetType plusxplusz = {{1,0,1}};
  NeighborhoodIteratorType::OffsetType minusxminusz = {{-1,0,-1}};
  NeighborhoodIteratorType::OffsetType plusxminusz = {{1,0,-1}};
  NeighborhoodIteratorType::OffsetType minusxplusz = {{-1,0,1}};

  NeighborhoodIteratorType::OffsetType plusyplusz = {{0,1,1}};
  NeighborhoodIteratorType::OffsetType minusyminusz = {{0,-1,-1}};
  NeighborhoodIteratorType::OffsetType plusyminusz = {{0,1,-1}};
  NeighborhoodIteratorType::OffsetType minusyplusz = {{0,-1,1}};


  // get the spacing in x and y for derivatives computations
  ImageType::SpacingType inputSpacing = reader->GetOutput()->GetSpacing();
  float spacingx = inputSpacing[0];
  float spacingy = inputSpacing[1];
  float spacingz = inputSpacing[2];

  // temporary variables for the iterations
  PixelType dplusx, dminusx, dplusy, dminusy, dplusz, dminusz,
            dx, dy, dz,
            d2x, d2y, d2z,
            dxy, dxz, dyz,
            dxsquare, dysquare, dzsquare;
  //equation members
  PixelType Hg, Upwind, nextphi;


  std::cout << "begin of iteration loop" << std::endl;

  for (int iter = 0; iter< m_NumIter; ++iter)
    {
    // input iterator (on phi)
    NeighborhoodIteratorType::RadiusType radius;
    radius.Fill(1);
    NeighborhoodIteratorType it( radius, m_InputImage,
                                 m_InputImage->GetRequestedRegion() );

    // output iterator ( on phi(t+1) )
    IteratorType out(m_OutputImage, m_InputImage->GetRequestedRegion());

    // we change the evolution of the surface from graph to level set
    if (iter == nbiter/2)
    curvatureFactor = curvatureFactor/1000000;



  // one iteration :
  for (it.GoToBegin(),edgit.GoToBegin(),edgderivit.GoToBegin(), out.GoToBegin();
       !it.IsAtEnd();
       ++it, ++out, ++edgderivit, ++edgit)
    {

    // derivatives computation :
    dminusx = (it.GetCenterPixel()-it.GetPixel(minusx))/spacingx;
    dminusy = (it.GetCenterPixel()-it.GetPixel(minusy))/spacingy;
    dminusz = (it.GetCenterPixel()-it.GetPixel(minusz))/spacingz;


    dplusx = (it.GetPixel(plusx)-it.GetCenterPixel())/spacingx;
    dplusy = (it.GetPixel(plusy)-it.GetCenterPixel())/spacingy;
    dplusz = (it.GetPixel(plusz)-it.GetCenterPixel())/spacingz;


    dx = (it.GetPixel(plusx) - it.GetPixel(minusx))/(2*spacingx);
    dy = (it.GetPixel(plusy) - it.GetPixel(minusy))/(2*spacingy);
    dz = (it.GetPixel(plusz) - it.GetPixel(minusz))/(2*spacingz);


    d2x = dplusx-dminusx;
    d2y = dplusy-dminusy;
    d2z = dplusz-dminusz;


    dxy = (it.GetPixel(plusxplusy)
           + it.GetPixel(minusxminusy)
           - it.GetPixel(plusxminusy)
           - it.GetPixel(minusxplusy))
          / (4*spacingx*spacingy);
    dxz = (it.GetPixel(plusxplusz)
           + it.GetPixel(minusxminusz)
           - it.GetPixel(plusxminusz)
           - it.GetPixel(minusxplusz))
          / (4*spacingx*spacingz);
    dyz = (it.GetPixel(plusyplusz)
          + it.GetPixel(minusyminusz)
          - it.GetPixel(plusyminusz)
          - it.GetPixel(minusyplusz))
          / (4*spacingy*spacingz);


    dxsquare = dx*dx;
    dysquare = dy*dy;
    dzsquare = dz*dz;


    // evolution equation
    Hg = edgit.Get() * (  (curvatureFactor + dxsquare + dysquare)*d2z
                         +(curvatureFactor + dxsquare + dzsquare)*d2y
                         +(curvatureFactor + dzsquare + dysquare)*d2x
                         -2 * (dx*dz*dxz+dx*dy*dxy+dy*dz*dyz))
            / (curvatureFactor + dxsquare + dysquare + dzsquare) ;

    edgederivVector = edgderivit.Get();

    Upwind =  std::min(-edgederivVector[0],0.)*dplusx
            + std::max(-edgederivVector[0],0.)*dminusx
            + std::min(-edgederivVector[1],0.)*dplusy
            + std::max(-edgederivVector[1],0.)*dminusy
            + std::min(-edgederivVector[2],0.)*dplusz
            + std::max(-edgederivVector[2],0.)*dminusz;

    nextphi = it.GetCenterPixel()+m_DeltaT*(nu*Hg-rho*Upwind);

    out.Set(nextphi);
    }
  // input is now output
    //(smart pointer dereference automatically
    // the memory previously pointed by m_InputImage,
    // and m_TempImage gets out of scope at the end of the iteration)
  m_TempImage = m_InputImage;
  m_InputImage = m_OutputImage;
  m_OutputImage = m_TempImage;
  }
  this->GraftOutput( outputImage  );
  }
}




/**
 *
 */
template< class InputImageType, class OutputImageType >
void
SubjectiveSurfaceEvolutionFilter<InputImageType, OutputImageType>::
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
