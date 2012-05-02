/*=========================================================================
 
 Program:   Insight Segmentation & Registration Toolkit
 Module:    $RCSfile: regTemporalImageMetric.txx,v $
 Language:  C++
 Date:      $Date: 2007-11-12 20:00:37 $
 Version:   $Revision: 1.30 $
 
 Copyright (c) Insight Software Consortium. All rights reserved.
 See ITKCopyright.txt or http://www.itk.org/HTML/Copyright.htm for details.
 
 This software is distributed WITHOUT ANY WARRANTY; without even 
 the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR 
 PURPOSE.  See the above copyright notices for more information.
 
 =========================================================================*/
#ifndef __regTemporalImageMetric_txx
#define __regTemporalImageMetric_txx

// First, make sure that we include the configuration file.
// This line may be removed once the ThreadSafeTransform gets
// integrated into ITK.
#include "itkConfigure.h"

// Second, redirect to the optimized version if necessary

#include "regTemporalImageMetric.h"

namespace itk
{
  
  /**
   * Constructor
   */
  template <class TImage> 
  TemporalImageMetric<TImage>
  ::TemporalImageMetric()
  {
    
    m_Image		  = 0; // has to be provided by the user.
    m_Transform     = 0; // has to be provided by the user.
    m_Interpolator  = 0; // has to be provided by the user.
    m_GradientImage = 0; // will receive the output of the filter;
    m_ComputeGradient = true; // metric computes gradient by default
    m_GradientImage = NULL; // computed at initialization
    m_NumberOfSamples = 500;
    m_SamplesContainer.resize(0);
    m_NumberOfValidSamples=0;
    m_Initialized = false;
    m_ReferenceTime = 0.;
    m_InterpolatorIsBSpline = false;
    m_MovingMaskTimeSpecified=false;
    
    
  }
  
  /**
   * Destructor
   */
  template <class TImage> 
  TemporalImageMetric<TImage>
  ::~TemporalImageMetric()
  {
    
  }
  
  ///**
  // * Set the parameters that define a unique transform
  // */
  template <class TImage> 
  void
  TemporalImageMetric<TImage>
  ::SetTransformParameters( const ParametersType & parameters ) const
  {
    if( !m_Transform )
    {
      itkExceptionMacro(<<"Transform has not been assigned");
    }
    m_Transform->SetParameters( parameters );
  }
  
  
  ///**
  // * Initialize
  // */
  template <class TImage> 
  void
  TemporalImageMetric<TImage>
  ::Initialize(void) throw ( ExceptionObject )
  {
    
	const unsigned int SpaceDimension = ImageDimension - 1;
    
    if( !m_Transform )
    {
      itkExceptionMacro(<<"Transform is not present");
    }
    
    if( !m_Interpolator )
    {
      itkExceptionMacro(<<"Interpolator is not present");
    }
    
    if( !m_Image )
    {
      itkExceptionMacro(<<"Input image is not present");
    }
    
    if( m_ImageRegion.GetNumberOfPixels() == 0 )
    {
      itkExceptionMacro(<<"FixedImageRegion is empty");
    }
    
    // If the image is provided by a source, update the source.
    if( m_Image->GetSource() )
    {
      m_Image->GetSource()->Update();
    }
    
    // Make sure the ImageRegion is within the Image buffered region
    if ( !m_ImageRegion.Crop( m_Image->GetBufferedRegion() ) )
    {
      itkExceptionMacro(
                        <<"ImageRegion does not overlap the image buffered region" );
    }
    
    m_Interpolator->SetInputImage( m_Image );
    
    if ( m_ComputeGradient )
    {
      this->ComputeGradient();
    }
    
    if (m_ReferenceTime == 0.)
    {
      m_ReferenceTime = m_Image->GetOrigin()[SpaceDimension];
    }



    // If there are any observers on the metric, call them to give the
    // user code a chance to set parameters on the metric
    m_Initialized = true;
    this->InvokeEvent( InitializeEvent() );
  }
  
  
  ///*
  // * Compute the gradient image and assign it to m_GradientImage.
  // */
  template <class TImage> 
  void
  TemporalImageMetric<TImage>
  ::ComputeGradient() 
  {
    GradientImageFilterPointer gradientFilter
    = GradientImageFilterType::New();
    
    gradientFilter->SetInput( m_Image );
    
    const typename ImageType::SpacingType&
    spacing = m_Image->GetSpacing();
    double maximumSpacing=0.0;
    for(unsigned int i=0; i<ImageDimension; i++)
    {
      if( spacing[i] > maximumSpacing )
      {
        maximumSpacing = spacing[i];
      }
    }
    gradientFilter->SetSigma( maximumSpacing );
    gradientFilter->SetNormalizeAcrossScale( true );
    
#ifdef ITK_USE_ORIENTED_IMAGE_DIRECTION
    gradientFilter->SetUseImageDirection( true );
#endif
    
    gradientFilter->Update();
    
    m_GradientImage = gradientFilter->GetOutput();
  }
  
  ///**
  // * PrintSelf
  // */
  template <class TImage> 
  void
  TemporalImageMetric<TImage>
  ::PrintSelf(std::ostream& os, Indent indent) const
  {
    Superclass::PrintSelf( os, indent );
    os << indent << "ComputeGradient: "
    << static_cast<typename NumericTraits<bool>::PrintType>(m_ComputeGradient)
    << std::endl;
    os << indent << "Image: " << m_Image.GetPointer()   << std::endl;
    os << indent << "Gradient Image: " << m_GradientImage.GetPointer() 
    << std::endl;
    os << indent << "Transform:    " << m_Transform.GetPointer()    << std::endl;
    os << indent << "Interpolator: " << m_Interpolator.GetPointer() << std::endl;
    os << indent << "ImageRegion: " << m_ImageRegion << std::endl;
    os << indent << "Space Image Mask: " << m_SpaceImageMask.GetPointer() 
    << std::endl;
    os << indent << "Space Moving Image Mask: " << m_SpaceMovingImageMask.GetPointer() 
    << std::endl;
    os << indent << "Number of valid Samples: " << m_NumberOfValidSamples 
    << std::endl;
    os << indent << "Number of Samples: " << m_NumberOfSamples 
    << std::endl;
    
  }
   
  ///**
  // * Sample first image domain
  // */
  template <class TImage> 
  void
  TemporalImageMetric<TImage>
  ::SampleReferenceImageDomain(bool denseSampling, unsigned int timeMargin) const
  {
    
    // Iterator type.
    typedef ImageRandomConstIteratorWithIndex<ImageType> IteratorType;
    
    ImageConstPointer image = this->m_Image;
    
    // Defining iterator region
    typename ImageType::RegionType iteratorRegion = this->GetImageRegion();
    typename ImageType::RegionType::SizeType iteratorRegionSize = iteratorRegion.GetSize();
    typename ImageType::RegionType::IndexType iteratorRegionIndex = iteratorRegion.GetIndex();    
    
    if (!denseSampling)
    {
      // In case of a fixed reference, the size of the region in time is set to 1
      //   and the frame closest to the reference time is taken as index 
      
      iteratorRegionSize[ImageDimension-1] = 1;
      double contTime = ( m_ReferenceTime - m_Image->GetOrigin()[ImageDimension-1] )
                            / m_Image->GetSpacing()[ImageDimension-1];
      iteratorRegionIndex[ImageDimension-1] = (long int)(floor(contTime));
      
    }
	
    else 
    {
      // We are in the case of the 4D sampling
      if (iteratorRegionSize[ImageDimension-1] > timeMargin) 
      {
        iteratorRegionSize[ImageDimension-1] -= timeMargin;
      }
      else 
      {
        itkExceptionMacro(<< "At least four images are required in the input image sequence." << std::endl);
      }
    }
	

    iteratorRegion.SetSize( iteratorRegionSize );
    iteratorRegion.SetIndex( iteratorRegionIndex);
   
    std::cout << "iteratorRegion" << iteratorRegion << std::endl;
    
    // Creating iterator
    IteratorType ti( image, iteratorRegion );
    
    // Number of valid samples counted inside the mask
    m_NumberOfValidSamples = 0;
    // Maximum number of times we will call the random iterator
    unsigned long maxCount = 50 * m_NumberOfSamples;
    ti.SetNumberOfSamples( maxCount );
    
    // Resize vector of spatial samples and create iterator
    m_SamplesContainer.resize(m_NumberOfSamples);
    
    typename std::vector<DiffeoSpatialSample>::iterator itVectorSamples
    = m_SamplesContainer.begin();
    
    ImageIndexType index;
    InputPointType point;
    
    ti.GoToBegin();
    
    // TMP Saving sample positions to a mask image for debugging purposes
    typedef Image<unsigned char, SpaceImageDimension> SamplesImageType;
    typename SamplesImageType::Pointer imageSamples = SamplesImageType::New();
    typename SamplesImageType::SizeType sSize;
    for (unsigned int d=0; d<SpaceImageDimension; d++) {
      sSize[d] = image->GetBufferedRegion().GetSize()[d];
    }
    typename SamplesImageType::RegionType sRegion;
    sRegion.SetSize(sSize);
    imageSamples->SetRegions(sRegion);
    typename SamplesImageType::SpacingType sSpacing;
    for (unsigned int d=0; d<SpaceImageDimension; d++) {
      sSpacing[d] = image->GetSpacing()[d];
    }
    imageSamples->SetSpacing(sSpacing);
    typename SamplesImageType::PointType sOrigin;
    for (unsigned int d=0; d<SpaceImageDimension; d++) {
      sOrigin[d] = image->GetOrigin()[d];
    }
    imageSamples->SetOrigin(sOrigin);
    imageSamples->Allocate();
    imageSamples->FillBuffer(0);
    // END TMP
    
    unsigned long count = 1;
    while(itVectorSamples != m_SamplesContainer.end())
    {
      if (count >= maxCount)
      {
        // In case we are over the maximum number of times we were allowed to call
        // the iterator, we resize the vector to the current number of samples
        // and exit
        m_SamplesContainer.resize(this->m_NumberOfValidSamples);
        break;
      }
      
      index = ti.GetIndex();
      image->TransformIndexToPhysicalPoint(index, point);
      
      SpacePointType point2;
      
      for (unsigned int d=0; d<SpaceImageDimension; d++)
      {
        point2[d] = point[d];
      }
      
      if( this->m_SpaceImageMask.IsNotNull() 
            && !this->m_SpaceImageMask->IsInside( point2 ) )
      {
        ++ti;
        count++;
        continue;
      }
      
      const RealType fixedValue = ti.Value();
      
      itVectorSamples->m_FixedValue = fixedValue;
      itVectorSamples->m_FixedPoint = point;
      
      this->m_NumberOfValidSamples++;
      
      // TMP 
      typename SamplesImageType::IndexType sIndex;
      for (unsigned int d=0; d<SpaceImageDimension; d++)
      {
        sIndex[d] = ti.GetIndex()[d];
      }
      imageSamples->SetPixel(sIndex, 1);
      // END TMP
      
	  

      ++ti;
      ++itVectorSamples;
      count++;
    }
  
    std::cout << "Number of samples is : "<<m_SamplesContainer.size()<<std::endl;
    
    // TMP
    typedef ImageFileWriter<SamplesImageType> SamplesWriterType;
    typename SamplesWriterType::Pointer sWriter = SamplesWriterType::New();
    sWriter->SetFileName("imageSamples.mhd");
    sWriter->SetInput(imageSamples);
    sWriter->Update();
    // END TMP
    
  }
  
  
  

   ///**
    
  template <class TImage> 
  void
  TemporalImageMetric<TImage>
  ::SetInterpolator(InterpolatorType* interpolator) 
  {

    InterpolatorTypeBS* bsInterp = dynamic_cast<InterpolatorTypeBS*>(interpolator);
    if (bsInterp != 0)
    {
      m_InterpolatorBS = bsInterp;
      m_Interpolator = bsInterp;
      m_InterpolatorIsBSpline = true;
      m_ComputeGradient = false;
    }
    else 
    {
      m_Interpolator = interpolator;
    }

  
  }
    
    
    
  
  
} // end namespace itk

#endif
