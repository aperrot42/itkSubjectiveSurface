#ifndef __itkInitialSubjectiveSurfaceSource_h
#define __itkInitialSubjectiveSurfaceSource_h



#include "itkImageSource.h"
#include "itkFixedArray.h"
#include "itkSize.h"


namespace itk
{
/** \class InitialSubjectiveSurfaceSource
 *
 * \brief Class for generating a peak function
 * or an inverted distance fonction from a
 * focus point.
 */

template< class TOutputImage >
class ITK_EXPORT InitialSubjectiveSurfaceSource :
      public ImageSource< TOutputImage >
 {

 public:
   /** Standard class typedefs. */
   typedef InitialSubjectiveSurfaceSource        Self;
   typedef ImageSource<TOutputImage>  Superclass;
   typedef SmartPointer<Self>         Pointer;
   typedef SmartPointer<const Self>   ConstPointer;

   /** Typedef for the output image PixelType. */
   typedef typename TOutputImage::PixelType OutputImagePixelType;

   /** Typedef to describe the output image region type. */
   typedef typename TOutputImage::RegionType OutputImageRegionType;

   /** Spacing typedef support.  Spacing holds the size of a pixel.  The
    * spacing is the geometric distance between image samples. */
   typedef typename TOutputImage::SpacingType SpacingType;

   /** Origin typedef support.  The origin is the geometric coordinates
    * of the index (0,0). */
   typedef typename TOutputImage::PointType PointType;

   /** Direction typedef support.  The direction is the direction
    * cosines of the image. */
   typedef typename TOutputImage::DirectionType DirectionType;

   /** Dimensionality of the output image */
   itkStaticConstMacro(NDimensions, unsigned int, TOutputImage::ImageDimension);

   /** Type used to store gaussian parameters. */
   typedef FixedArray<double, itkGetStaticConstMacro(NDimensions)> ArrayType;

   /** Size type matches that used for images */
   typedef typename TOutputImage::SizeType         SizeType;
   typedef typename TOutputImage::SizeValueType    SizeValueType;

   /** Run-time type information (and related methods). */
   itkTypeMacro(GaussianImageSource,ImageSource);

   /** Method for creation through the object factory. */
   itkNewMacro(Self);


   /** Specify the size of the output image. */
   itkSetMacro(Size, SizeType);
   itkSetVectorMacro(Size, unsigned long, NDimensions);

   /** Get the size of the output image. */
   itkGetConstReferenceMacro(Size, SizeType);

   /** Specify the spacing of the output image. */
   itkSetMacro(Spacing, SpacingType);
   itkSetVectorMacro(Spacing, const float, NDimensions);

   /** Get the spacing of the output image. */
   itkGetConstReferenceMacro(Spacing, SpacingType);

   /** Specify the origin of the output image. */
   itkSetMacro(Origin, PointType);
   itkSetVectorMacro(Origin, const float, NDimensions);

   /** Get the origin of the output image. */
   itkGetConstReferenceMacro(Origin, PointType);

   /** Specify the direction of the output image. */
   itkSetMacro(Direction, DirectionType);
   itkGetConstReferenceMacro(Direction, DirectionType);

   /** Specify the point of focus (max of the initial subjective surface)*/
   itkSetMacro(FocusPoint, PointType);
   itkGetConstReferenceMacro(FocusPoint, PointType);

   /** Specify the function to be generated :
    *  0: peak at focuspoint
    *  1: inverted distance from focuspoint */
   itkSetMacro(Function, unsigned char);
   itkGetMacro(Function, unsigned char);


protected:
  InitialSubjectiveSurfaceSource();
  ~InitialSubjectiveSurfaceSource();
  void PrintSelf(std::ostream & os, Indent indent) const;
  void GenerateData();
  virtual void GenerateOutputInformation();



 private:


  SizeType      m_Size;              //size of the output image
  SpacingType   m_Spacing;           //spacing
  PointType     m_Origin;            //origin
  DirectionType m_Direction;         // direction


  PointType      m_FocusPoint; // focus point (peak of the function)
  unsigned char  m_Function; // variable for funtion selection




  InitialSubjectiveSurfaceSource(const Self &); //purposely not implemented
  void operator=(const Self &);             //purposely not implemented


};
} // end namespace itk

#ifndef ITK_MANUAL_INSTANTIATION
#include "itkInitialSubjectiveSurfaceSource.txx"
#endif

#endif
