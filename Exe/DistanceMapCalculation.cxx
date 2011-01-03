/*=========================================================================

  Program:   Insight Segmentation & Registration Toolkit
  Module:    $RCSfile: itkMultiScaleHessianSmoothed3DToMembranenessMeasureImageFilterTest.cxx,v $
  Language:  C++
  Date:      $Date: 2007/04/01 21:19:46 $
  Version:   $Revision: 1.5 $

  Copyright (c) Insight Software Consortium. All rights reserved.
  See ITKCopyright.txt or http://www.itk.org/HTML/Copyright.htm for details.

     This software is distributed WITHOUT ANY WARRANTY; without even
     the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
     PURPOSE.  See the above copyright notices for more information.

=========================================================================*/
#include "itkImageFileReader.h"
#include "itkImageFileWriter.h"
#include "itkRescaleIntensityImageFilter.h"
#include "itkSignedMaurerDistanceMapImageFilter.h"
#include "itkNumericTraits.h"

int main(int argc, char* argv [] )
{
  if ( argc < 3 )
    {
    std::cerr << "Missing Parameters: "
              << argv[0] << std::endl
              << "InputImage(itkimage) OutputImage(itkimage)" << std::endl;
    return EXIT_FAILURE;
    }
  // Define the dimension of the images
  const int Dimension = 3;

  // Declare the types of the images
  typedef float       InputPixelType;
  typedef itk::Image< InputPixelType, Dimension>  InputImageType;
  typedef float       OutputPixelType;
  typedef itk::Image< OutputPixelType, Dimension>   OutputImageType;

  // input reader
  typedef itk::ImageFileReader< InputImageType  > ReaderType;

  // Output writer
  typedef itk::ImageFileWriter< OutputImageType > WriterType;

  // distance function filter
  typedef itk::SignedMaurerDistanceMapImageFilter< InputImageType, InputImageType >
    DistanceFilterType;


  std::cout << "reading input image" << std::endl;
  ReaderType::Pointer reader = ReaderType::New();
  reader->SetFileName ( argv[1] );
  reader->Update();


  std::cout << "distance map computing" << std::endl;
  DistanceFilterType::Pointer distanceFilter
    = DistanceFilterType::New();
  distanceFilter->SetInput(reader->GetOutput());
  distanceFilter->SetUseImageSpacing(true);
  distanceFilter->SetInsideIsPositive(false);
  distanceFilter->SetSquaredDistance (false);
  distanceFilter->SetBackgroundValue(itk::NumericTraits< OutputPixelType >::Zero );
  distanceFilter->Update();

  WriterType::Pointer writer = WriterType::New();
  writer->SetFileName( argv[2] );
  writer->SetInput ( distanceFilter->GetOutput() );

  std::cout << "writing output image" << std::endl;
  try
    {
    writer->Update();
    }
  catch( itk::ExceptionObject & err )
    {
    std::cerr << "Exception caught: " << err << std::endl;
    return EXIT_FAILURE;
    }


  return EXIT_SUCCESS;
}

