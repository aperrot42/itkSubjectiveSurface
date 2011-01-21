#ifndef __itkSubjectiveSurfaceEvolutionFilter_txx
#define __itkSubjectiveSurfaceEvolutionFilter_txx

#include "itkSubjectiveSurfaceEvolutionFilter.h"


#include "itkImage.h"


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


  // set phi as input of the duplicator
  m_ImageDuplicator->SetInputImage(this->GetInput(0));



  // depending on dimension of ImageType, only one version of
  // GenerateData will actually be instantiated
  // The DimType<dim> * parameter never used, but it's a low-cost way
  // to force compile-time type discrimination
  this->GenerateData(static_cast< DimensionType<InputImageType::ImageDimension>*>(0));

  this->GraftOutput( m_OutputImage  );

}


  template< class InputImageType, class OutputImageType >
  void
  SubjectiveSurfaceEvolutionFilter<InputImageType, OutputImageType>::
  GenerateData(DimensionType<2> *)
  {

  if (InputImageType::ImageDimension == 2)
    {
    std::cout << "processing2Dimage"<<std::endl;
    // duplicate phi for const correctness
    m_ImageDuplicator->Update();
    m_InputImage = m_ImageDuplicator->GetOutput();

    // allocation of phi(t+1)
    m_OutputImage = InputImageType::New();
    m_OutputImage->SetRegions(m_InputImage->GetRequestedRegion());
    m_OutputImage->SetSpacing(m_InputImage->GetSpacing());
    m_OutputImage->Allocate();

    // pointer to edgemap (g)
    m_EdgeDetector = const_cast
                     <InputImageType *>(this->GetInput(1));
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



  // offsets definition for neighborhood computations :
  NeighborhoodOffsetType minusx = {{-1,0}};
  NeighborhoodOffsetType plusx = {{1,0}};
  NeighborhoodOffsetType minusy = {{0,-1}};
  NeighborhoodOffsetType plusy = {{0,1}};
  NeighborhoodOffsetType plusxplusy = {{1,1}};
  NeighborhoodOffsetType minusxminusy = {{-1,-1}};
  NeighborhoodOffsetType plusxminusy = {{1,-1}};
  NeighborhoodOffsetType minusxplusy = {{-1,1}};

  // get the spacing in x and y for derivatives computations
  typename InputImageType::SpacingType inputSpacing = m_InputImage->GetSpacing();
  float spacingx = inputSpacing[0];
  float spacingy = inputSpacing[1];


  // temporary variables for the iterations
    //derivatives
  InputPixelType dplusx, dminusx, dplusy, dminusy,
            dx, dy,
            d2x, d2y,
            dxy,
            dxsquare, dysquare;
    //equation members
  InputPixelType Hg, Upwind, nextphi;


  for (unsigned int iter = 0; iter< m_NumberIteration; ++iter)
    {
    // input iterator (on phi)
    typename NeighborhoodIteratorType::RadiusType radius;
    radius.Fill(1);
    NeighborhoodIteratorType it( radius, m_InputImage,
                                 m_InputImage->GetRequestedRegion() );

    // output iterator ( on phi(t+1) )
    IteratorType out(m_OutputImage, m_InputImage->GetRequestedRegion());

//    // we change the evolution of the surface from graph to level set
//    if (iter == m_NumberIteration/2)
//    m_a = m_a/1000000;



  // one iteration :
  for (it.GoToBegin(),edgit.GoToBegin(),edgderivit.GoToBegin(), out.GoToBegin();
       !it.IsAtEnd();
       ++it, ++out, ++edgderivit, ++edgit)
      {

      // derivatives computation :
      dminusx = (it.GetCenterPixel()-it.GetPixel(minusx))/spacingx;
      dminusy = (it.GetCenterPixel()-it.GetPixel(minusy))/spacingy;

      dplusx = (it.GetPixel(plusx)-it.GetCenterPixel())/spacingx;
      dplusy = (it.GetPixel(plusy)-it.GetCenterPixel())/spacingy;

      dx = (it.GetPixel(plusx) - it.GetPixel(minusx))/(2*spacingx);
      dy = (it.GetPixel(plusy) - it.GetPixel(minusy))/(2*spacingy);

      d2x = (dplusx-dminusx)/spacingx;
      d2y = (dplusy-dminusy)/spacingy;

      dxy = (it.GetPixel(plusxplusy)
             + it.GetPixel(minusxminusy)
             - it.GetPixel(plusxminusy)
             - it.GetPixel(minusxplusy))
            / (4*spacingx*spacingy);

      dxsquare = dx*dx;
      dysquare = dy*dy;


      // evolution equation
        // Hg = g*H (with H mean curvature of graph of phi)
      Hg = edgit.Get() * (   (m_a  + dxsquare)*d2y
                           + (m_a  + dysquare)*d2x
                           - 2 * (dx*dy*dxy))
           / (m_a  + dxsquare + dysquare ) ;

        // grad(g)
      edgederivVector = edgderivit.Get();

        // upwind differntiation term : grad(I).*grad(g)
      Upwind = std::min(-edgederivVector[0],(double)0.)*dplusx
             + std::max(-edgederivVector[0],(double)0.)*dminusx
             + std::min(-edgederivVector[1],(double)0.)*dplusy
             + std::max(-edgederivVector[1],(double)0.)*dminusy;


      // value of phi(t+1)
      nextphi = it.GetCenterPixel()+m_DeltaT*(m_Nu*Hg-m_Rho*Upwind);

      out.Set(nextphi);
      }
    // input is now output
      //(smart pointer dereference automatically
      // the memory previously pointed by inputphi,
      // and tempImagePoint gets out of scope at the end of the iteration)
  m_TempImage = m_InputImage;
  m_InputImage = m_OutputImage;
  m_OutputImage = m_TempImage;
    }
  }


}








  template< class InputImageType, class OutputImageType >
  void
  SubjectiveSurfaceEvolutionFilter<InputImageType, OutputImageType>::
  GenerateData(DimensionType<3> *)
  {

if (InputImageType::ImageDimension >= 3)
  {
  // duplicate phi for const correctness
  m_ImageDuplicator->Update();
  m_InputImage = m_ImageDuplicator->GetOutput();

  // allocation of phi(t+1)
  m_OutputImage = InputImageType::New();
  m_OutputImage->SetRegions(m_InputImage->GetRequestedRegion());
  m_OutputImage->SetSpacing(m_InputImage->GetSpacing());
  m_OutputImage->Allocate();



  // pointer to edgemap (g)
  m_EdgeDetector = const_cast
                   <InputImageType *>(this->GetInput(1));
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
  NeighborhoodOffsetType minusx = {{-1,0,0}};
  NeighborhoodOffsetType plusx = {{1,0,0}};

  NeighborhoodOffsetType minusy = {{0,-1,0}};
  NeighborhoodOffsetType plusy = {{0,1,0}};

  NeighborhoodOffsetType minusz = {{0,0,-1}};
  NeighborhoodOffsetType plusz = {{0,0,1}};


  NeighborhoodOffsetType plusxplusy = {{1,1,0}};
  NeighborhoodOffsetType minusxminusy = {{-1,-1,0}};
  NeighborhoodOffsetType plusxminusy = {{1,-1,0}};
  NeighborhoodOffsetType minusxplusy = {{-1,1,0}};

  NeighborhoodOffsetType plusxplusz = {{1,0,1}};
  NeighborhoodOffsetType minusxminusz = {{-1,0,-1}};
  NeighborhoodOffsetType plusxminusz = {{1,0,-1}};
  NeighborhoodOffsetType minusxplusz = {{-1,0,1}};

  NeighborhoodOffsetType plusyplusz = {{0,1,1}};
  NeighborhoodOffsetType minusyminusz = {{0,-1,-1}};
  NeighborhoodOffsetType plusyminusz = {{0,1,-1}};
  NeighborhoodOffsetType minusyplusz = {{0,-1,1}};


  // get the spacing in x and y for derivatives computations
  typename InputImageType::SpacingType inputSpacing = m_InputImage->GetSpacing();
  float spacingx = inputSpacing[0];
  float spacingy = inputSpacing[1];
  float spacingz = inputSpacing[2];

  // temporary variables for the iterations
  InputPixelType dplusx, dminusx, dplusy, dminusy, dplusz, dminusz,
            dx, dy, dz,
            d2x, d2y, d2z,
            dxy, dxz, dyz,
            dxsquare, dysquare, dzsquare;
  //equation members
  InputPixelType Hg, Upwind, nextphi;


  for (unsigned int iter = 0; iter< m_NumberIteration; ++iter)
    {
    // input iterator (on phi)
    typename NeighborhoodIteratorType::RadiusType radius;
    radius.Fill(1);
    NeighborhoodIteratorType it( radius, m_InputImage,
                                 m_InputImage->GetRequestedRegion() );

    // output iterator ( on phi(t+1) )
    IteratorType out(m_OutputImage, m_InputImage->GetRequestedRegion());

//    // we change the evolution of the surface from graph to level set
//    if (iter == m_NumberIteration/2)
//      m_a = m_a/1000000;



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
      Hg = edgit.Get() * (  (m_a + dxsquare + dysquare)*d2z
                           +(m_a + dxsquare + dzsquare)*d2y
                           +(m_a + dzsquare + dysquare)*d2x
                           -2 * (dx*dz*dxz+dx*dy*dxy+dy*dz*dyz))
              / (m_a + dxsquare + dysquare + dzsquare) ;

      edgederivVector = edgderivit.Get();

      Upwind =  std::min(-edgederivVector[0],0.)*dplusx
              + std::max(-edgederivVector[0],0.)*dminusx
              + std::min(-edgederivVector[1],0.)*dplusy
              + std::max(-edgederivVector[1],0.)*dminusy
              + std::min(-edgederivVector[2],0.)*dplusz
              + std::max(-edgederivVector[2],0.)*dminusz;

      nextphi = it.GetCenterPixel()+m_DeltaT*(m_Nu*Hg-m_Rho*Upwind);

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
  }

}

  // General Case ND
  template< class InputImageType, class OutputImageType >
  template <unsigned dimension>
  void
  SubjectiveSurfaceEvolutionFilter<InputImageType, OutputImageType>::
  GenerateData(DimensionType <dimension> *)
  {
    std::cout << "this filter does not process images of dimension : "
              << InputImageType::ImageDimension
              << std::endl;
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

  os << indent << "Nu "
     << static_cast< typename NumericTraits< OutputPixelType >::PrintType >( m_Nu ) << std::endl;
  os << indent << "Rho: "
     << static_cast< typename NumericTraits< OutputPixelType >::PrintType >( m_Rho ) << std::endl;
  os << indent << "a: "
     << static_cast< typename NumericTraits< OutputPixelType >::PrintType >( m_a ) << std::endl;
  os << indent << "Number of iterations: "
     << static_cast< int >( m_NumberIteration ) << std::endl;
}

}
#endif
