/*=========================================================================

  Program:   Insight Segmentation & Registration Toolkit
  Module:    $RCSfile: itkDiffeomorphicBSplineTransformTest.cxx,v $
  Language:  C++
  Date:      $Date: 2008-06-20 15:25:58 $
  Version:   $Revision: 1.2 $

  Copyright (c) Insight Software Consortium. All rights reserved.
  See ITKCopyright.txt or http://www.itk.org/HTML/Copyright.htm for details.

     This software is distributed WITHOUT ANY WARRANTY; without even 
     the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR 
     PURPOSE.  See the above copyright notices for more information.

=========================================================================*/
#if defined(_MSC_VER)
#pragma warning ( disable : 4786 )
#endif

#include <iostream>

#include "itkDiffeomorphicBSplineTransform.h"
#include "itkImage.h"
#include "itkImageFileReader.h"


int main( int argc, char * argv[] )
{
  
  const    unsigned int    ImageDimension = 3;
  typedef  signed short    PixelType;
  
  typedef itk::Image< PixelType, ImageDimension >  FixedImageType;
  typedef itk::Image< PixelType, ImageDimension >  MovingImageType;
  
  typedef itk::DiffeomorphicBSplineTransform<double,ImageDimension, 3>  TransformType;
  TransformType::Pointer  transform = TransformType::New();
  
  
  typedef itk::ImageFileReader< FixedImageType  > FixedImageReaderType;
  typedef itk::ImageFileReader< MovingImageType > MovingImageReaderType;
  
  FixedImageReaderType::Pointer  fixedImageReader  = FixedImageReaderType::New();
  MovingImageReaderType::Pointer movingImageReader = MovingImageReaderType::New();
  
  fixedImageReader->SetFileName(  argv[1] );
  movingImageReader->SetFileName( argv[2] );
  
  FixedImageType::ConstPointer fixedImage = fixedImageReader->GetOutput();
  
  try
    {
    fixedImageReader->Update();
    movingImageReader->Update();
    }
  catch( itk::ExceptionObject & err )
    {
    std::cerr << "ExceptionObject caught !" << std::endl;
    std::cerr << err << std::endl;
    return -1;
    }
  
  FixedImageType::RegionType fixedRegion = fixedImage->GetBufferedRegion();
  
  typedef TransformType::RegionType RegionType;
  RegionType bsplineRegion;
  RegionType::SizeType   gridSizeOnImage;
  RegionType::SizeType   gridBorderSize;
  RegionType::SizeType   totalGridSize;
  
  gridSizeOnImage.Fill( 12 );
  gridBorderSize.Fill( 3 );    // Border for spline order = 3 ( 1 lower, 2 upper )
  totalGridSize = gridSizeOnImage + gridBorderSize;

  bsplineRegion.SetSize( totalGridSize );

  typedef TransformType::SpacingType SpacingType;
  SpacingType spacing = fixedImage->GetSpacing();

  typedef TransformType::OriginType OriginType;
  OriginType origin = fixedImage->GetOrigin();;

  FixedImageType::SizeType fixedImageSize = fixedRegion.GetSize();

  for(unsigned int r=0; r<ImageDimension; r++)
    {
    spacing[r] *= floor( static_cast<double>(fixedImageSize[r] - 1)  /
                  static_cast<double>(gridSizeOnImage[r] - 1) );
    origin[r]  -=  spacing[r];
    }
  
  transform->DebugOn();
  transform->SetGridSpacing( spacing );
  transform->SetGridOrigin( origin );
  transform->SetGridRegion( bsplineRegion );
  transform->SetNumberOfTimeSteps(10);
  
  TransformType::ParametersType parameters(transform->GetNumberOfParameters());
  parameters.Fill(0.);
  transform->SetParameters(parameters);

  TransformType::InputPointType input_point;
  input_point.Fill(0.);
  input_point = transform->TransformPoint(input_point);

  std::cout<< input_point << std::endl;


  const TransformType::JacobianType & jaco = 
    transform->GetJacobian(input_point);

  transform->InsertTransform(0);

    
  return EXIT_SUCCESS;
}
