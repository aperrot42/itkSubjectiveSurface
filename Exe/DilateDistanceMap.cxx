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
#include "itkGrayscaleDilateImageFilter.h"
#include "itkNumericTraits.h"
#include "itkBinaryBallStructuringElement.h"
#include "itkImage.h"

int main(int argc, char* argv [] )
{
  if ( argc < 4 )
    {
    std::cerr << "Missing Parameters: "
              << argv[0] << std::endl
              << "InputImage(itkimage) "
              << "OutputImage(itkimage) "
              << "smallestNucleiRadius(double)" << std::endl;
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

  // binary ball
  typedef itk::BinaryBallStructuringElement< InputPixelType,
  Dimension >	StructuringElementType;
  // dilate filter
  typedef itk::GrayscaleDilateImageFilter< InputImageType,  InputImageType, StructuringElementType >
    DilateFilterType;

  std::cout << "reading input image" << std::endl;
  ReaderType::Pointer reader = ReaderType::New();
  reader->SetFileName ( argv[1] );
  reader->Update();

  float   m_smallestRadius = atof(argv[3]);

  StructuringElementType::RadiusType structuringRadius;
  StructuringElementType structuringElement;
  
  InputImageType::Pointer inputImage = reader->GetOutput();

  for (int i=0; i< Dimension;++i)
    {
    
    structuringRadius[i] = static_cast<unsigned int>(
                             m_smallestRadius
                             /(inputImage->GetSpacing()[i]) );
    }
  std::cout << structuringRadius << std::endl;

  structuringElement.SetRadius(structuringRadius); // 3x3 structuring element

  structuringElement.CreateStructuringElement();

  DilateFilterType::Pointer grayscaleDilateFilter = DilateFilterType::New();
  grayscaleDilateFilter->SetKernel( structuringElement );

  std::cout << "dilate input image with a sphere of radius : " 
            << structuringRadius
            << std::endl;
  grayscaleDilateFilter->SetNumberOfThreads(6);
  grayscaleDilateFilter->SetInput(inputImage);
  grayscaleDilateFilter->Update();

  WriterType::Pointer writer = WriterType::New();
  writer->SetFileName( argv[2] );
  writer->SetInput ( grayscaleDilateFilter->GetOutput() );

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
