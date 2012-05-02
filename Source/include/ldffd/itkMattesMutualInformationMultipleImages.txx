/*=========================================================================

  Program:   Insight Segmentation & Registration Toolkit
  Module:    $RCSfile: itkMattesMutualInformationMultipleImages.txx,v $
  Language:  C++
  Date:      $Date: 2008-07-01 10:39:17 $
  Version:   $Revision: 1.3 $

  Copyright (c) Insight Software Consortium. All rights reserved.
  See ITKCopyright.txt or http://www.itk.org/HTML/Copyright.htm for details.

     This software is distributed WITHOUT ANY WARRANTY; without even 
     the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR 
     PURPOSE.  See the above copyright notices for more information.

=========================================================================*/
#ifndef _itkMattesMutualInformationMultipleImages_txx
#define _itkMattesMutualInformationMultipleImages_txx

#include "itkMattesMutualInformationMultipleImages.h"
#include "itkBSplineInterpolateImageFunction.h"
#include "itkCovariantVector.h"
#include "itkImageRandomConstIteratorWithIndex.h"
#include "itkImageRegionConstIterator.h"
#include "itkImageRegionIterator.h"
#include "itkImageIterator.h"
#include "vnl/vnl_math.h"
#include "itkBSplineDeformableTransform.h"

namespace itk
{

/**
 * Constructor
 */
template < class TFixedImage, class TMovingImage >
  MattesMutualInformationMultipleImages<TFixedImage,TMovingImage>
::MattesMutualInformationMultipleImages()
{

  m_NumberOfSpatialSamples = 500;
  m_NumberOfHistogramBins = 50;

  this->SetComputeGradient(false); // don't use the default gradient for now

  m_InterpolatorIsBSpline = false;
  m_TransformIsBSpline    = false;

  // Initialize PDFs to NULL
  m_JointPDF = NULL;
  m_JointPDFDerivatives = NULL;

  typename BSplineTransformType::Pointer transformer =
    BSplineTransformType::New();
  this->SetTransform (transformer);

  /* FIXME: This should be removed one we do not herit
   *          anymore from ImageToImageMetric. */
  typename BSplineInterpolatorType::Pointer interpolator =
    BSplineInterpolatorType::New();
  this->SetInterpolator (interpolator);

  // Intializing vectors of images and interpolators
  m_MovingImages.resize(0);
  m_BSplineInterpolators.resize(0);
  m_DerivativeCalculators.resize(0);
  m_Interpolators.resize(0);

  // Initialize memory
  m_MovingImageNormalizedMin = 0.0;
  m_FixedImageNormalizedMin = 0.0;
  m_MovingImageTrueMin = 0.0;
  m_MovingImageTrueMax = 0.0;
  m_FixedImageBinSize = 0.0;
  m_MovingImageBinSize = 0.0;
  m_CubicBSplineDerivativeKernel = NULL;
  m_BSplineInterpolator = NULL;
  m_DerivativeCalculator = NULL;
  m_NumParametersPerDim = 0;
  m_NumBSplineWeights = 0;
  m_BSplineTransform = NULL;
  m_NumberOfParameters = 0;
  m_UseAllPixels = false;
  m_ReseedIterator = false;
  m_RandomSeed = -1;
}


/**
 * Print out internal information about this class
 */
template < class TFixedImage, class TMovingImage  >
void
MattesMutualInformationMultipleImages<TFixedImage,TMovingImage>
::PrintSelf(std::ostream& os, Indent indent) const
{

  Superclass::PrintSelf(os, indent);

  os << indent << "NumberOfSpatialSamples: ";
  os << m_NumberOfSpatialSamples << std::endl;
  os << indent << "NumberOfHistogramBins: ";
  os << m_NumberOfHistogramBins << std::endl;
  os << indent << "UseAllPixels: ";
  os << m_UseAllPixels << std::endl;

  // Debugging information
  os << indent << "NumberOfParameters: ";
  os << m_NumberOfParameters << std::endl;
  os << indent << "FixedImageNormalizedMin: ";
  os << m_FixedImageNormalizedMin << std::endl;
  os << indent << "MovingImageNormalizedMin: ";
  os << m_MovingImageNormalizedMin << std::endl;
  os << indent << "MovingImageTrueMin: ";
  os << m_MovingImageTrueMin << std::endl;
  os << indent << "MovingImageTrueMax: ";
  os << m_MovingImageTrueMax << std::endl;
  os << indent << "FixedImageBinSize: "; 
  os << m_FixedImageBinSize << std::endl;
  os << indent << "MovingImageBinSize: ";
  os << m_MovingImageBinSize << std::endl;
  os << indent << "InterpolatorIsBSpline: ";
  os << m_InterpolatorIsBSpline << std::endl;
  os << indent << "TransformIsBSpline: ";
  os << m_TransformIsBSpline << std::endl;

}


/**
 * Initialize
 */
template <class TFixedImage, class TMovingImage> 
void
  MattesMutualInformationMultipleImages<TFixedImage,TMovingImage>
::Initialize(void) throw ( ExceptionObject )
{

  this->Superclass::Initialize();

  // Cache the number of transformation parameters
  m_NumberOfParameters = this->m_Transform->GetNumberOfParameters();

  /**
   * Compute the minimum and maximum for the FixedImage over
   * the FixedImageRegion.
   *
   * NB: We can't use StatisticsImageFilter to do this because
   * the filter computes the min/max for the largest possible region
   */
  double fixedImageMin = NumericTraits<double>::max();
  double fixedImageMax = NumericTraits<double>::NonpositiveMin();

  typedef ImageRegionConstIterator<FixedImageType> FixedIteratorType;
  FixedIteratorType fixedImageIterator( 
    this->m_FixedImage, this->GetFixedImageRegion() );

  for ( fixedImageIterator.GoToBegin(); 
    !fixedImageIterator.IsAtEnd(); ++fixedImageIterator )
    {

    double sample = static_cast<double>( fixedImageIterator.Get() );

    if ( sample < fixedImageMin )
      {
      fixedImageMin = sample;
      }

    if ( sample > fixedImageMax )
      {
      fixedImageMax = sample;
      }
    }


  /**
   * Compute the minimum and maximum for all the entire moving images
   */
  double movingImageMin = NumericTraits<double>::max();
  double movingImageMax = NumericTraits<double>::NonpositiveMin();

  for (unsigned int movImage = 0; movImage< m_MovingImages.size(); movImage++)
    {

    typedef ImageRegionConstIterator<MovingImageType> MovingIteratorType;
    MovingIteratorType movingImageIterator(
      this->m_MovingImages[movImage]->image, this->m_MovingImages[movImage]->image->GetBufferedRegion() );

    for ( movingImageIterator.GoToBegin(); 
      !movingImageIterator.IsAtEnd(); ++movingImageIterator)
      {
      double sample = static_cast<double>( movingImageIterator.Get() );

      if ( sample < movingImageMin )
        {
        movingImageMin = sample;
        }

      if ( sample > movingImageMax )
        {
        movingImageMax = sample;
        }
      }

    }

  m_MovingImageTrueMin = movingImageMin;
  m_MovingImageTrueMax = movingImageMax;

  itkDebugMacro( " FixedImageMin: " << fixedImageMin << 
    " FixedImageMax: " << fixedImageMax << std::endl );
  itkDebugMacro( " MovingImageMin: " << movingImageMin << 
    " MovingImageMax: " << movingImageMax << std::endl );


  /**
   * Compute binsize for the histograms.
   *
   * The binsize for the image intensities needs to be adjusted so that 
   * we can avoid dealing with boundary conditions using the cubic 
   * spline as the Parzen window.  We do this by increasing the size
   * of the bins so that the joint histogram becomes "padded" at the 
   * borders. Because we are changing the binsize, 
   * we also need to shift the minimum by the padded amount in order to 
   * avoid minimum values filling in our padded region.
   *
   * Note that there can still be non-zero bin values in the padded region,
   * it's just that these bins will never be a central bin for the Parzen
   * window.
   *
   */
  const int padding = 2;  // this will pad by 2 bins

  m_FixedImageBinSize = ( fixedImageMax - fixedImageMin ) /
    static_cast<double>( m_NumberOfHistogramBins - 2 * padding );
  m_FixedImageNormalizedMin = fixedImageMin / m_FixedImageBinSize - 
    static_cast<double>( padding );

  m_MovingImageBinSize = ( movingImageMax - movingImageMin ) /
    static_cast<double>( m_NumberOfHistogramBins - 2 * padding );
  m_MovingImageNormalizedMin = movingImageMin / m_MovingImageBinSize -
    static_cast<double>( padding );


  itkDebugMacro( "FixedImageNormalizedMin: " << m_FixedImageNormalizedMin );
  itkDebugMacro( "MovingImageNormalizedMin: " << m_MovingImageNormalizedMin );
  itkDebugMacro( "FixedImageBinSize: " << m_FixedImageBinSize );
  itkDebugMacro( "MovingImageBinSize; " << m_MovingImageBinSize );


  if( m_UseAllPixels )
    {
    m_NumberOfSpatialSamples = 
      this->GetFixedImageRegion().GetNumberOfPixels();
    }

  /**
   * Allocate memory for the fixed image sample container.
   */
  m_FixedImageSamples.resize( m_NumberOfSpatialSamples );

  /**
   * Allocate memory for the marginal PDF and initialize values
   * to zero. The marginal PDFs are stored as std::vector.
   */
  m_FixedImageMarginalPDF.resize( m_NumberOfHistogramBins, 0.0 );
  m_MovingImageMarginalPDF.resize( m_NumberOfHistogramBins, 0.0 );

  /**
   * Allocate memory for the joint PDF and joint PDF derivatives.
   * The joint PDF and joint PDF derivatives are store as itk::Image.
   */
  m_JointPDF = JointPDFType::New();
  m_JointPDFDerivatives = JointPDFDerivativesType::New();

  // Instantiate a region, index, size
  JointPDFRegionType            jointPDFRegion;
  JointPDFIndexType              jointPDFIndex;
  JointPDFSizeType              jointPDFSize;

  JointPDFDerivativesRegionType  jointPDFDerivativesRegion;
  JointPDFDerivativesIndexType  jointPDFDerivativesIndex;
  JointPDFDerivativesSizeType    jointPDFDerivativesSize;

  // For the joint PDF define a region starting from {0,0} 
  // with size {m_NumberOfHistogramBins, m_NumberOfHistogramBins}.
  // The dimension represents fixed image parzen window index
  // and moving image parzen window index, respectively.
  jointPDFIndex.Fill( 0 ); 
  jointPDFSize.Fill( m_NumberOfHistogramBins ); 

  jointPDFRegion.SetIndex( jointPDFIndex );
  jointPDFRegion.SetSize( jointPDFSize );

  // Set the regions and allocate
  m_JointPDF->SetRegions( jointPDFRegion );
  m_JointPDF->Allocate();

  // For the derivatives of the joint PDF define a region starting from {0,0,0} 
  // with size {m_NumberOfParameters,m_NumberOfHistogramBins, 
  // m_NumberOfHistogramBins}. The dimension represents transform parameters,
  // fixed image parzen window index and moving image parzen window index,
  // respectively. 
  jointPDFDerivativesIndex.Fill( 0 ); 
  jointPDFDerivativesSize[0] = m_NumberOfParameters;
  jointPDFDerivativesSize[1] = m_NumberOfHistogramBins;
  jointPDFDerivativesSize[2] = m_NumberOfHistogramBins;

  jointPDFDerivativesRegion.SetIndex( jointPDFDerivativesIndex );
  jointPDFDerivativesRegion.SetSize( jointPDFDerivativesSize );

  // Set the regions and allocate
  m_JointPDFDerivatives->SetRegions( jointPDFDerivativesRegion );
  m_JointPDFDerivatives->Allocate();


  /**
   * Setup the kernels used for the Parzen windows.
   */
  m_CubicBSplineKernel = CubicBSplineFunctionType::New();
  m_CubicBSplineDerivativeKernel = CubicBSplineDerivativeFunctionType::New();    


  if( m_UseAllPixels )
    {
    /** 
     * Take all the pixels within the fixed image region)
     * to create the sample points list.
     */
    this->SampleFullFixedImageDomain( m_FixedImageSamples );
    }
  else
    {
    /** 
     * Uniformly sample the fixed image (within the fixed image region)
     * to create the sample points list.
     */
    this->SampleFixedImageDomain( m_FixedImageSamples );
    }

  /**
   * Pre-compute the fixed image parzen window index for 
   * each point of the fixed image sample points list.
   */
  this->ComputeFixedImageParzenWindowIndices( m_FixedImageSamples );


  // If the vector of interpolators is empty, fill it with BSpline interpolators
  if  ( m_Interpolators.size() == 0 )
    {
    m_Interpolators.resize( m_MovingImages.size() );
    for (unsigned int movImage=0; movImage < m_MovingImages.size(); movImage++)
      {
      this->m_Interpolators[movImage] = 
        BSplineInterpolatorType::New();
      this->m_Interpolators[movImage]->SetInputImage( this->m_MovingImages[movImage]->image );
      }
    }

  // Check if all the interpolators are BSpline
  m_InterpolatorIsBSpline = true;
  for (unsigned int movImage=0; movImage < m_MovingImages.size(); movImage++)
    {
    BSplineInterpolatorType * testPtr = dynamic_cast<BSplineInterpolatorType *>(
      this->m_Interpolators[movImage].GetPointer() );
    if (!testPtr)
      {
      m_InterpolatorIsBSpline = false;
      break;
      }
    }

  // Checking sizes
  if ( m_DerivativeCalculators.size() != m_MovingImages.size() )
    {
    m_DerivativeCalculators.resize( m_MovingImages.size() );
    } 
  if ( m_BSplineInterpolators.size() != m_MovingImages.size() )
    {
    m_BSplineInterpolators.resize( m_MovingImages.size() );
    } 

  // Derivative calculator 
  for (unsigned int movImage=0; movImage < m_MovingImages.size(); movImage++)
    {
    if (!m_InterpolatorIsBSpline)
      {
      m_DerivativeCalculators[movImage] = DerivativeFunctionType::New();
      m_DerivativeCalculators[movImage]->SetInputImage( this->m_MovingImages[movImage]->image );
      m_BSplineInterpolators[movImage] = 0;
      }
    else
      {
      m_BSplineInterpolators[movImage] = dynamic_cast<BSplineInterpolatorType *>(
        this->m_Interpolators[movImage].GetPointer() );
      m_DerivativeCalculators[movImage] = 0;
      }
    }


}


/**
 * Uniformly sample the fixed image domain using a random walk
 */
template < class TFixedImage, class TMovingImage >
void
  MattesMutualInformationMultipleImages<TFixedImage,TMovingImage>
::SampleFixedImageDomain( FixedImageSpatialSampleContainer& samples )
{

  // Set up a random interator within the user specified fixed image region.
  typedef ImageRandomConstIteratorWithIndex<FixedImageType> RandomIterator;
  RandomIterator randIter( this->m_FixedImage, this->GetFixedImageRegion() );

  randIter.SetNumberOfSamples( m_NumberOfSpatialSamples );
  randIter.GoToBegin();

  typename FixedImageSpatialSampleContainer::iterator iter;
  typename FixedImageSpatialSampleContainer::const_iterator end=samples.end();

  if( this->m_FixedImageMask )
    {

    typename Superclass::InputPointType inputPoint;

    iter=samples.begin();
    int count = 0;
    int samples_found = 0;
    int maxcount = m_NumberOfSpatialSamples * 10;
    while( iter != end )
      {

      if ( count > maxcount )
        {
#if 0
        itkExceptionMacro( 
          "Drew too many samples from the mask (is it too small?): "
          << maxcount << std::endl );
#else
        samples.resize(samples_found);
        //        this->SetNumberOfSpatialSamples(sample_found);
        break;
#endif
        }
      count++;

      // Get sampled index
      FixedImageIndexType index = randIter.GetIndex();
      // Check if the Index is inside the mask, translate index to point
      this->m_FixedImage->TransformIndexToPhysicalPoint( index, inputPoint );

      // If not inside the mask, ignore the point
      if( !this->m_FixedImageMask->IsInside( inputPoint ) )
        {
        ++randIter; // jump to another random position
        continue;
        }

      // Get sampled fixed image value
      (*iter).FixedImageValue = randIter.Get();

      // Translate index to point
      (*iter).FixedImagePointValue = inputPoint;
      samples_found++;
      // Jump to random position
      ++randIter;
      ++iter;
      }
    }
  else
    {
    for( iter=samples.begin(); iter != end; ++iter )
      {
      // Get sampled index
      FixedImageIndexType index = randIter.GetIndex();
      // Get sampled fixed image value
      (*iter).FixedImageValue = randIter.Get();
      // Translate index to point
      this->m_FixedImage->TransformIndexToPhysicalPoint( index,
        (*iter).FixedImagePointValue );
      // Jump to random position
      ++randIter;

      }
    }
}

/**
 * Sample the fixed image domain using all pixels in the Fixed image region
 */
template < class TFixedImage, class TMovingImage >
void
  MattesMutualInformationMultipleImages<TFixedImage,TMovingImage>
::SampleFullFixedImageDomain( FixedImageSpatialSampleContainer& samples )
{

  // Set up a region interator within the user specified fixed image region.
  typedef ImageRegionConstIteratorWithIndex<FixedImageType> RegionIterator;
  RegionIterator regionIter( this->m_FixedImage, this->GetFixedImageRegion() );

  regionIter.GoToBegin();

  typename FixedImageSpatialSampleContainer::iterator iter;
  typename FixedImageSpatialSampleContainer::const_iterator end=samples.end();

  if( this->m_FixedImageMask )
    {
    typename Superclass::InputPointType inputPoint;

    iter=samples.begin();
    unsigned long nSamplesPicked = 0;

    while( iter != end && !regionIter.IsAtEnd() )
      {
      // Get sampled index
      FixedImageIndexType index = regionIter.GetIndex();
      // Check if the Index is inside the mask, translate index to point
      this->m_FixedImage->TransformIndexToPhysicalPoint( index, inputPoint );

      // If not inside the mask, ignore the point
      if( !this->m_FixedImageMask->IsInside( inputPoint ) )
        {
        ++regionIter; // jump to next pixel
        continue;
        }

      // Get sampled fixed image value
      (*iter).FixedImageValue = regionIter.Get();
      // Translate index to point
      (*iter).FixedImagePointValue = inputPoint;

      ++regionIter;
      ++iter;
      ++nSamplesPicked;
      }

    // If we picked fewer samples than the desired number, 
    // resize the container
    if (nSamplesPicked != this->m_NumberOfSpatialSamples)
      {
      this->m_NumberOfSpatialSamples = nSamplesPicked;
      samples.resize(this->m_NumberOfSpatialSamples);
      }
    }
  else // not restricting sample throwing to a mask
    {

    // cannot sample more than the number of pixels in the image region
    if (  this->m_NumberOfSpatialSamples 
      > this->GetFixedImageRegion().GetNumberOfPixels())
      {
      this->m_NumberOfSpatialSamples 
        = this->GetFixedImageRegion().GetNumberOfPixels();
      samples.resize(this->m_NumberOfSpatialSamples);
      }

    for( iter=samples.begin(); iter != end; ++iter )
      {
      // Get sampled index
      FixedImageIndexType index = regionIter.GetIndex();
      // Get sampled fixed image value
      (*iter).FixedImageValue = regionIter.Get();
      // Translate index to point
      this->m_FixedImage->TransformIndexToPhysicalPoint( index,
        (*iter).FixedImagePointValue );
      ++regionIter;
      }
    }
}

/**
 * 
 */
template < class TFixedImage, class TMovingImage >
void
MattesMutualInformationMultipleImages<TFixedImage,TMovingImage>
::ComputeFixedImageParzenWindowIndices( 
  FixedImageSpatialSampleContainer& samples )
{

  typename FixedImageSpatialSampleContainer::iterator iter;
  typename FixedImageSpatialSampleContainer::const_iterator end=samples.end();

  for( iter=samples.begin(); iter != end; ++iter )
    {

    // Determine parzen window arguments (see eqn 6 of Mattes paper [2]).  
    double windowTerm =
      static_cast<double>( (*iter).FixedImageValue ) / m_FixedImageBinSize -
      m_FixedImageNormalizedMin;
    unsigned int pindex = static_cast<unsigned int>( vcl_floor(windowTerm ) );

    // Make sure the extreme values are in valid bins
    if ( pindex < 2 )
      {
      pindex = 2;
      }
    else if ( pindex > (m_NumberOfHistogramBins - 3) )
      {
      pindex = m_NumberOfHistogramBins - 3;
      }

    (*iter).FixedImageParzenWindowIndex = pindex;

    }

}

/**
 * Get the match Measure
 */
template < class TFixedImage, class TMovingImage  >
typename MattesMutualInformationMultipleImages<TFixedImage,TMovingImage>
::MeasureType
MattesMutualInformationMultipleImages<TFixedImage,TMovingImage>
::GetValue( const ParametersType& parameters ) const
{

  // Reset marginal pdf to all zeros.
  // Assumed the size has already been set to NumberOfHistogramBins
  // in Initialize().
  for ( unsigned int j = 0; j < m_NumberOfHistogramBins; j++ )
    {
    m_FixedImageMarginalPDF[j]  = 0.0;
    m_MovingImageMarginalPDF[j] = 0.0;
    }

  // Reset the joint pdfs to zero
  m_JointPDF->FillBuffer( 0.0 );

  // Set up the parameters in the transform
  this->m_Transform->SetParameters( parameters );

  // Declare iterators for iteration over the sample container
  typename FixedImageSpatialSampleContainer::const_iterator fiter;
  typename FixedImageSpatialSampleContainer::const_iterator fend = 
    m_FixedImageSamples.end();

  unsigned long nSamples=0;
  unsigned long nFixedImageSamples=0;

  // mapped coordinates and image values
  MovingImageValuesArrayType movingImageValues( this->m_MovingImages.size() );
  MovingImagePointArrayType mappedPoints;

  for ( fiter = m_FixedImageSamples.begin(); fiter != fend; ++fiter )
    {
    bool sampleOk;

    // Compute set of mapped coordinates and
    // get moving image value
    this->TransformPoint (nFixedImageSamples, parameters, mappedPoints, 
      sampleOk, movingImageValues );

    ++nFixedImageSamples;

    if( sampleOk )
      {

      ++nSamples;
      // Since a zero-order BSpline (box car) kernel is used for
      // the fixed image marginal pdf, we need only increment the
      // fixedImageParzenWindowIndex by value of 1.0.
      m_FixedImageMarginalPDF[(*fiter).FixedImageParzenWindowIndex] += 
        static_cast<PDFValueType>( 1 );

      // For each moving image, we need to compute the contribution to the joint 
      // histogram
      for (unsigned int movImage=0; movImage<m_MovingImages.size(); movImage++)
        {

        double movingImageValue = movingImageValues[movImage];

        /**
         * Compute this sample's contribution to the marginal and 
         * joint distributions.
         */

        // Determine parzen window arguments (see eqn 6 of Mattes paper [2]).    
        double movingImageParzenWindowTerm =
          movingImageValue / m_MovingImageBinSize - m_MovingImageNormalizedMin;
        unsigned int movingImageParzenWindowIndex = 
          static_cast<unsigned int>( vcl_floor(movingImageParzenWindowTerm ) );

        // Make sure the extreme values are in valid bins
        if ( movingImageParzenWindowIndex < 2 )
          {
          movingImageParzenWindowIndex = 2;
          }
        else if ( movingImageParzenWindowIndex > (m_NumberOfHistogramBins - 3) )
          {
          movingImageParzenWindowIndex = m_NumberOfHistogramBins - 3;
          }

        /**
         * The region of support of the parzen window determines which bins
         * of the joint PDF are effected by the pair of image values.
         * Since we are using a cubic spline for the moving image parzen
         * window, four bins are affected.  The fixed image parzen window is
         * a zero-order spline (box car) and thus effects only one bin.
         *
         *  The PDF is arranged so that moving image bins corresponds to the 
         * zero-th (column) dimension and the fixed image bins corresponds
         * to the first (row) dimension.
         *
         */

        // Pointer to affected bin to be updated
        JointPDFValueType *pdfPtr = m_JointPDF->GetBufferPointer() +
          ( (*fiter).FixedImageParzenWindowIndex 
            * m_JointPDF->GetOffsetTable()[1] );

        // Move the pointer to the first affected bin
        int pdfMovingIndex = static_cast<int>( movingImageParzenWindowIndex ) - 1;
        pdfPtr += pdfMovingIndex;

        for (; pdfMovingIndex <= static_cast<int>( movingImageParzenWindowIndex )
          + 2;
          pdfMovingIndex++, pdfPtr++ )
          {

          // Update PDF for the current intensity pair
          double movingImageParzenWindowArg = 
            static_cast<double>( pdfMovingIndex ) - 
            static_cast<double>( movingImageParzenWindowTerm );

          *(pdfPtr) += static_cast<PDFValueType>( 
            m_CubicBSplineKernel->Evaluate( movingImageParzenWindowArg ) );

          }  //end parzen windowing for loop

        } //end if-block check sampleOk//end iterating over moving image
      } //end if-block check sampleOk
    }// end iterating over fixed image spatial sample container for loop


  itkDebugMacro( "Ratio of voxels mapping into moving image buffer: " 
    << nSamples << " / " << m_NumberOfSpatialSamples 
    << std::endl );

  if( nSamples < m_NumberOfSpatialSamples / 4 )
    {
    itkExceptionMacro( "Too many samples map outside moving image buffer: "
      << nSamples << " / " << m_NumberOfSpatialSamples 
      << std::endl );
    }


  /**
   * Normalize the PDFs, compute moving image marginal PDF
   *
   */
  typedef ImageRegionIterator<JointPDFType> JointPDFIteratorType;
  JointPDFIteratorType jointPDFIterator ( m_JointPDF,
    m_JointPDF->GetBufferedRegion() );

  jointPDFIterator.GoToBegin();

  // Compute joint PDF normalization factor 
  // (to ensure joint PDF sum adds to 1.0)
  double jointPDFSum = 0.0;

  while( !jointPDFIterator.IsAtEnd() )
    {
    jointPDFSum += jointPDFIterator.Get();
    ++jointPDFIterator;
    }

  if ( jointPDFSum == 0.0 )
    {
    itkExceptionMacro( "Joint PDF summed to zero" );
    }


  // Normalize the PDF bins
  jointPDFIterator.GoToEnd();
  while( !jointPDFIterator.IsAtBegin() )
    {
    --jointPDFIterator;
    jointPDFIterator.Value() /= static_cast<PDFValueType>( jointPDFSum );
    }


  // Normalize the fixed image marginal PDF
  double fixedPDFSum = 0.0;
  for( unsigned int bin = 0; bin < m_NumberOfHistogramBins; bin++ )
    {
    fixedPDFSum += m_FixedImageMarginalPDF[bin];
    }

  if ( fixedPDFSum == 0.0 )
    {
    itkExceptionMacro( "Fixed image marginal PDF summed to zero" );
    }

  for( unsigned int bin=0; bin < m_NumberOfHistogramBins; bin++ )
    {
    m_FixedImageMarginalPDF[bin] /= static_cast<PDFValueType>( fixedPDFSum );
    }


  // Compute moving image marginal PDF by summing over fixed image bins.
  typedef ImageLinearIteratorWithIndex<JointPDFType> JointPDFLinearIterator;
  JointPDFLinearIterator linearIter( m_JointPDF,
    m_JointPDF->GetBufferedRegion() );

  linearIter.SetDirection( 1 );
  linearIter.GoToBegin();
  unsigned int movingIndex = 0;

  while( !linearIter.IsAtEnd() )
    {

    double sum = 0.0;

    while( !linearIter.IsAtEndOfLine() )
      {
      sum += linearIter.Get();
      ++linearIter;
      }

    m_MovingImageMarginalPDF[movingIndex] = static_cast<PDFValueType>(sum);

    linearIter.NextLine();
    ++movingIndex;

    }

  /**
   * Compute the metric by double summation over histogram.
   */

  // Setup pointer to point to the first bin
  JointPDFValueType * jointPDFPtr = m_JointPDF->GetBufferPointer();

  // Initialze sum to zero
  double sum = 0.0;

  for( unsigned int fixedIndex = 0;
    fixedIndex < m_NumberOfHistogramBins;
    ++fixedIndex )
    {

    double fixedImagePDFValue = m_FixedImageMarginalPDF[fixedIndex];

    for( unsigned int movingIndex = 0;
      movingIndex < m_NumberOfHistogramBins; 
      ++movingIndex, jointPDFPtr++ )      
      {

      double movingImagePDFValue = m_MovingImageMarginalPDF[movingIndex];
      double jointPDFValue = *(jointPDFPtr);

      // check for non-zero bin contribution
      if( jointPDFValue > 1e-16 &&  movingImagePDFValue > 1e-16 )
        {

        double pRatio = vcl_log(jointPDFValue / movingImagePDFValue );
        if( fixedImagePDFValue > 1e-16)
          {
          sum += jointPDFValue * ( pRatio - vcl_log(fixedImagePDFValue ) );
          }

        }  // end if-block to check non-zero bin contribution
      }  // end for-loop over moving index
    }  // end for-loop over fixed index

  return static_cast<MeasureType>( -1.0 * sum );

}


/**
 * Get the both Value and Derivative Measure
 */
template < class TFixedImage, class TMovingImage  >
void
MattesMutualInformationMultipleImages<TFixedImage,TMovingImage>
::GetValueAndDerivative(
  const ParametersType& parameters,
  MeasureType& value,
  DerivativeType& derivative) const
{

  std::cout << "Inside get value and derivative" << std::endl;

  // This is really awfull because too much specific
  // I need to write an abstract class that will
  // encapsulate all necessary methods to interact with this class
  typedef itk::DiffeomorphicBSplineTransform< CoordinateRepresentationType,
  itkGetStaticConstMacro(FixedImageDimension), 3> MyTransformType;
  const MyTransformType* diffTransfo = dynamic_cast<const MyTransformType*>(this->GetTransform());

  // Set output values to zero
  value = NumericTraits< MeasureType >::Zero;
  derivative = DerivativeType( this->GetNumberOfParameters() );
  derivative.Fill( NumericTraits< MeasureType >::Zero );

  // Reset marginal pdf to all zeros.
  // Assumed the size has already been set to NumberOfHistogramBins
  // in Initialize().
  for ( unsigned int j = 0; j < m_NumberOfHistogramBins; j++ )
    {
    m_FixedImageMarginalPDF[j]  = 0.0;
    m_MovingImageMarginalPDF[j] = 0.0;
    }

  // Reset the joint pdfs to zero
  m_JointPDF->FillBuffer( 0.0 );
  m_JointPDFDerivatives->FillBuffer( 0.0 );

  // Set up the parameters in the transform
  this->m_Transform->SetParameters( parameters );

  // Declare iterators for iteration over the sample container
  typename FixedImageSpatialSampleContainer::const_iterator fiter;
  typename FixedImageSpatialSampleContainer::const_iterator fend = 
    m_FixedImageSamples.end();

  unsigned long nSamples=0;
  unsigned long nFixedImageSamples=0;

  // mapped coordinates and image values
  MovingImageValuesArrayType movingImageValues( this->m_MovingImages.size() );
  MovingImagePointArrayType mappedPoints;

  for ( fiter = m_FixedImageSamples.begin(); fiter != fend; ++fiter )
    {
    bool sampleOk;

    // Compute set of mapped coordinates and
    // get moving image value
    this->TransformPoint (nFixedImageSamples, parameters, mappedPoints, 
      sampleOk, movingImageValues );

    if( sampleOk )
      {

      ++nSamples;

      // Since a zero-order BSpline (box car) kernel is used for
      // the fixed image marginal pdf, we need only increment the
      // fixedImageParzenWindowIndex by value of 1.0.
      m_FixedImageMarginalPDF[(*fiter).FixedImageParzenWindowIndex] +=
        static_cast<PDFValueType>( 1 );

      // For each moving image, we need to compute the contribution to the joint 
      // histogram
      unsigned int oldTimePoint = 0;//m_MovingImages[0]->timePoint;
      // We are currently assuming that the fixed image is associated
      // to time 0. If we keep m_MovingImages[0]->timePoint as initial time point
      // we end up accumulating jacobians in an empty interval and the gradient
      // with respect to parameters of the first transformation is allways 0.
      for (unsigned int movImage=0; movImage<m_MovingImages.size(); movImage++)
        {
        double movingImageValue = movingImageValues[movImage];
        MovingImagePointType mappedPoint = mappedPoints[movImage];

        // Get moving image derivative at the mapped position
        ImageDerivativesType movingImageGradientValue;

        this->ComputeImageDerivatives( movImage, mappedPoint, movingImageGradientValue );
        /**
         * Compute this sample's contribution to the marginal 
         *   and joint distributions.
         *
         */
        // Determine parzen window arguments (see eqn 6 of Mattes paper [2]).    
        double movingImageParzenWindowTerm =
          movingImageValue / m_MovingImageBinSize - m_MovingImageNormalizedMin;
        unsigned int movingImageParzenWindowIndex = 
          static_cast<unsigned int>( vcl_floor(movingImageParzenWindowTerm ) );

        // Make sure the extreme values are in valid bins     
        if ( movingImageParzenWindowIndex < 2 )
          {
          movingImageParzenWindowIndex = 2;
          }
        else if ( movingImageParzenWindowIndex > (m_NumberOfHistogramBins - 3) )
          {
          movingImageParzenWindowIndex = m_NumberOfHistogramBins - 3;
          }

        /**
         * The region of support of the parzen window determines which bins
         * of the joint PDF are effected by the pair of image values.
         * Since we are using a cubic spline for the moving image parzen
         * window, four bins are effected.  The fixed image parzen window is
         * a zero-order spline (box car) and thus effects only one bin.
         *
         *  The PDF is arranged so that moving image bins corresponds to the 
         * zero-th (column) dimension and the fixed image bins corresponds
         * to the first (row) dimension.
         *
         */

        // Pointer to affected bin to be updated
        JointPDFValueType *pdfPtr = m_JointPDF->GetBufferPointer() +
          ( (*fiter).FixedImageParzenWindowIndex * m_NumberOfHistogramBins );

        // Move the pointer to the fist affected bin
        int pdfMovingIndex = static_cast<int>( movingImageParzenWindowIndex ) - 1;
        pdfPtr += pdfMovingIndex;

        // Physical jacobian accumulation
        diffTransfo->AccumulatePhysicalJacobiansInTimeInterval(mappedPoint, oldTimePoint, m_MovingImages[movImage]->timePoint);

        // Update joint PDF and joint PDF derivatives
        for (; pdfMovingIndex <= static_cast<int>( movingImageParzenWindowIndex )
          + 2;
          pdfMovingIndex++, pdfPtr++ )
          {

          // Update PDF for the current intensity pair
          double movingImageParzenWindowArg = 
            static_cast<double>( pdfMovingIndex ) - 
            static_cast<double>(movingImageParzenWindowTerm);

          *(pdfPtr) += static_cast<PDFValueType>( 
            m_CubicBSplineKernel->Evaluate( movingImageParzenWindowArg ) );

          // Compute the cubicBSplineDerivative for later repeated use.
          double cubicBSplineDerivativeValue = 
            m_CubicBSplineDerivativeKernel->Evaluate( 
              movingImageParzenWindowArg );

          // Compute PDF derivative contribution.
          this->ComputePDFDerivatives( movImage,
            nFixedImageSamples,
            pdfMovingIndex, 
            movingImageGradientValue, 
            cubicBSplineDerivativeValue );
          }  //end parzen windowing for loop

        oldTimePoint = m_MovingImages[movImage]->timePoint;
        } // end moving image loop
      } //end if-block check sampleOk
    ++nFixedImageSamples;
    } // end iterating over fixed image spatial sample container for loop

  itkDebugMacro( "Ratio of voxels mapping into moving image buffer: " 
    << nSamples << " / " << m_NumberOfSpatialSamples 
    << std::endl );

  if( nSamples < m_NumberOfSpatialSamples / 4 )
    {
    itkExceptionMacro( "Too many samples map outside moving image buffer: "
      << nSamples << " / " << m_NumberOfSpatialSamples 
      << std::endl );
    }

  this->m_NumberOfPixelsCounted = nSamples;

  /**
   * Normalize the PDFs, compute moving image marginal PDF
   *
   */
  typedef ImageRegionIterator<JointPDFType> JointPDFIteratorType;
  JointPDFIteratorType jointPDFIterator ( m_JointPDF,
    m_JointPDF->GetBufferedRegion() );

  jointPDFIterator.GoToBegin();


  // Compute joint PDF normalization factor 
  // (to ensure joint PDF sum adds to 1.0)
  double jointPDFSum = 0.0;

  while( !jointPDFIterator.IsAtEnd() )
    {
    jointPDFSum += jointPDFIterator.Get();
    ++jointPDFIterator;
    }

  if ( jointPDFSum == 0.0 )
    {
    itkExceptionMacro( "Joint PDF summed to zero" );
    }


  // Normalize the PDF bins
  jointPDFIterator.GoToEnd();
  while( !jointPDFIterator.IsAtBegin() )
    {
    --jointPDFIterator;
    jointPDFIterator.Value() /= static_cast<PDFValueType>( jointPDFSum );
    }


  // Normalize the fixed image marginal PDF
  double fixedPDFSum = 0.0;
  for( unsigned int bin = 0; bin < m_NumberOfHistogramBins; bin++ )
    {
    fixedPDFSum += m_FixedImageMarginalPDF[bin];
    }

  if ( fixedPDFSum == 0.0 )
    {
    itkExceptionMacro( "Fixed image marginal PDF summed to zero" );
    }

  for( unsigned int bin=0; bin < m_NumberOfHistogramBins; bin++ )
    {
    m_FixedImageMarginalPDF[bin] /= static_cast<PDFValueType>( fixedPDFSum );
    }


  // Compute moving image marginal PDF by summing over fixed image bins.
  typedef ImageLinearIteratorWithIndex<JointPDFType> JointPDFLinearIterator;
  JointPDFLinearIterator linearIter( 
    m_JointPDF, m_JointPDF->GetBufferedRegion() );

  linearIter.SetDirection( 1 );
  linearIter.GoToBegin();
  unsigned int movingIndex = 0;

  while( !linearIter.IsAtEnd() )
    {
    double sum = 0.0;

    while( !linearIter.IsAtEndOfLine() )
      {
      sum += linearIter.Get();
      ++linearIter;
      }

    m_MovingImageMarginalPDF[movingIndex] = static_cast<PDFValueType>(sum);

    linearIter.NextLine();
    ++movingIndex;
    }

  // Normalize the joint PDF derivatives by the test image binsize and nSamples
  typedef ImageRegionIterator<JointPDFDerivativesType> 
    JointPDFDerivativesIteratorType;
  JointPDFDerivativesIteratorType jointPDFDerivativesIterator (
    m_JointPDFDerivatives, 
    m_JointPDFDerivatives->GetBufferedRegion() 
  );

  jointPDFDerivativesIterator.GoToBegin();

  double nFactor = 1.0 / ( m_MovingImageBinSize 
    * static_cast<double>( nSamples ) );

  while( !jointPDFDerivativesIterator.IsAtEnd() )
    {
    jointPDFDerivativesIterator.Value() *= nFactor;
    ++jointPDFDerivativesIterator;
    }

  /**
   * Compute the metric by double summation over histogram.
   */

  // Setup pointer to point to the first bin
  JointPDFValueType * jointPDFPtr = m_JointPDF->GetBufferPointer();

  // Initialize sum to zero
  double sum = 0.0;

  for( unsigned int fixedIndex = 0;
    fixedIndex < m_NumberOfHistogramBins;
    ++fixedIndex )
    {
    double fixedImagePDFValue = m_FixedImageMarginalPDF[fixedIndex];

    for( unsigned int movingIndex = 0; movingIndex < m_NumberOfHistogramBins; 
      ++movingIndex, jointPDFPtr++ )      
      {
      double movingImagePDFValue = m_MovingImageMarginalPDF[movingIndex];
      double jointPDFValue = *(jointPDFPtr);

      // check for non-zero bin contribution
      if( jointPDFValue > 1e-16 &&  movingImagePDFValue > 1e-16 )
        {

        double pRatio = vcl_log(jointPDFValue / movingImagePDFValue );

        if( fixedImagePDFValue > 1e-16)
          {
          sum += jointPDFValue * ( pRatio - vcl_log(fixedImagePDFValue ) );
          }

        // move joint pdf derivative pointer to the right position
        JointPDFValueType * derivPtr = 
          m_JointPDFDerivatives->GetBufferPointer() 
          + ( fixedIndex 
            * m_JointPDFDerivatives->GetOffsetTable()[2] 
          ) 
          + ( movingIndex 
            * m_JointPDFDerivatives->GetOffsetTable()[1] 
          );

        double debug_maxderivabs=0.;
        double debug_maxderiv=0.;

        for( unsigned int parameter=0;
          parameter < m_NumberOfParameters; 
          ++parameter, derivPtr++ )
          {

          // Ref: eqn 23 of Thevenaz & Unser paper [3]
          derivative[parameter] -= (*derivPtr) * pRatio;

          if (fabs(*derivPtr) > debug_maxderivabs)
            {
            debug_maxderivabs=fabs(*derivPtr);
            debug_maxderiv=(*derivPtr);
            }

          }  // end for-loop over parameters

        }  // end if-block to check non-zero bin contribution
      }  // end for-loop over moving index
    }  // end for-loop over fixed index

  value = static_cast<MeasureType>( -1.0 * sum );

  double debug_maxderivabs=0.;
  double debug_maxderiv=0.;
  for( unsigned int parameter=0;
    parameter < m_NumberOfParameters;
    ++parameter )
    {
    if (fabs(derivative[parameter]) > debug_maxderivabs)
      {
      debug_maxderivabs=fabs(derivative[parameter]);
      debug_maxderiv=derivative[parameter];
      }
    } 

}

/**
 * Get the match measure derivative
 */
template < class TFixedImage, class TMovingImage  >
void
MattesMutualInformationMultipleImages<TFixedImage,TMovingImage>
::GetDerivative( const ParametersType& parameters,
  DerivativeType & derivative ) const
{
  MeasureType value;
  // call the combined version
  this->GetValueAndDerivative( parameters, value, derivative );
}

/**
 * Compute image derivatives using a central difference function
 * if we are not using a BSplineInterpolator, which includes
 * derivatives.
 */
template < class TFixedImage, class TMovingImage >
void
MattesMutualInformationMultipleImages<TFixedImage,TMovingImage>
::ComputeImageDerivatives( const unsigned int imageIndex,
  const MovingImagePointType& mappedPoint, 
  ImageDerivativesType& gradient ) const
{

  if( m_InterpolatorIsBSpline )
    {
    // Computed moving image gradient using derivative BSpline kernel.
    gradient = m_BSplineInterpolators[imageIndex]->EvaluateDerivative( mappedPoint );
    }
  else
    {
    // For all generic interpolator use central differencing.
    gradient = m_DerivativeCalculators[imageIndex]->Evaluate( mappedPoint );
    }

}

/**
 * Transform a point from FixedImage domain to the sequence
 * of moving images
 */ 
template < class TFixedImage, class TMovingImage >
void
MattesMutualInformationMultipleImages<TFixedImage,TMovingImage>
::TransformPoint( 
  unsigned int sampleNumber, 
  const ParametersType& parameters,
  MovingImagePointArrayType& mappedPoints,
  bool& sampleOk,
  MovingImageValuesArrayType& movingImageValues ) const
{

  // This is really awfull because too much specific
  // I need to write an abstract class that will
  // encapsulate all necessary methods to interact with this class
  typedef itk::DiffeomorphicBSplineTransform< CoordinateRepresentationType,
  itkGetStaticConstMacro(FixedImageDimension), 3> MyTransformType;

  const MyTransformType* diffTransfo = dynamic_cast<const MyTransformType*>(this->GetTransform());

  // check if sizes are correct
  if ( mappedPoints.size() != m_MovingImages.size() )
    mappedPoints.resize( m_MovingImages.size() );

  if ( movingImageValues.size() != m_MovingImages.size() )
    movingImageValues.resize(  m_MovingImages.size() );

  if (!diffTransfo)
    {
    itkExceptionMacro( << "This transform type is not supported." << std::endl);
    return;
    }

  // From now, we can assume we are dealing with a diffemorphic BSpline
  // transform ... 
  const MovingImagePointArrayType & completeSetOfMappedPoints = 
    diffTransfo->TransformPoint(
      m_FixedImageSamples[sampleNumber].FixedImagePointValue, 
      sampleOk);

  if (!sampleOk)
    {
    return;
    }

  // NB: sampleOk will be valid if the sample falls inside the valid region
  // for all BSpline transformations. Another condition is that
  // all mapped coordinates are inside moving image mask + that all mapped
  // samples fit inside the corresponding images
  //

  for (unsigned int movImage=0; movImage<m_MovingImages.size(); movImage++)
    {
    // Get time index corresponding to the current image
    unsigned int timeIndex= m_MovingImages[movImage]->timePoint - 1;//FIXME Check this...
    if (timeIndex < completeSetOfMappedPoints.size())
      {
      mappedPoints[movImage] = completeSetOfMappedPoints[timeIndex];
      }
    else
      {
      itkExceptionMacro( << "Invalid index in array of mapped points." << std::endl);
      }

    sampleOk = sampleOk && this->m_Interpolators[movImage]->IsInsideBuffer( mappedPoints[movImage] );

    if ( this->m_MovingImageMask )
      {
      // Check if mapped point is within the support region of the moving image
      // mask
      sampleOk = sampleOk && this->m_MovingImageMask->IsInside( mappedPoints[movImage] );
      }


    if ( !sampleOk )
      {
      break;
      }

    // Computing intensity at mapped point
    movingImageValues[movImage] = this->m_Interpolators[movImage]->Evaluate( mappedPoints[movImage] );


    // Make sure the interpolated value is within bounds. We might add an epsilon tolerance for
    // admitting an incorrect intensity value and reject the value if the difference is too big. 
    if ( movingImageValues[movImage] < m_MovingImageTrueMin )
      {
      movingImageValues[movImage] = m_MovingImageTrueMin;
      }

    if ( movingImageValues[movImage] > m_MovingImageTrueMax )
      {
      movingImageValues[movImage] = m_MovingImageTrueMax; 
      }
    }
}

/**
 * Compute PDF derivatives contribution for each parameter
 */
template < class TFixedImage, class TMovingImage >
void
MattesMutualInformationMultipleImages<TFixedImage,TMovingImage>
::ComputePDFDerivatives( unsigned int indexImage,
  unsigned int sampleNumber, 
  int pdfMovingIndex,
  const ImageDerivativesType& movingImageGradientValue,
  double cubicBSplineDerivativeValue ) const
{
  // This is really awfull because too much specific
  // I need to write an abstract class that will
  // encapsulate all necessary methods to interact with this class
  typedef itk::DiffeomorphicBSplineTransform< CoordinateRepresentationType,
  itkGetStaticConstMacro(FixedImageDimension), 3> MyTransformType;

  const MyTransformType* diffTransfo = dynamic_cast<const MyTransformType*>(this->GetTransform());

  // Update bins in the PDF derivatives for the current intensity pair
  JointPDFValueType * derivPtr = m_JointPDFDerivatives->GetBufferPointer() +
    ( m_FixedImageSamples[sampleNumber].FixedImageParzenWindowIndex
      * m_JointPDFDerivatives->GetOffsetTable()[2] ) +
    ( pdfMovingIndex * m_JointPDFDerivatives->GetOffsetTable()[1] );



  // Compute the transform Jacobian.
  typedef typename TransformType::JacobianType JacobianType;
  std::vector<long unsigned int> jacoColumns;

  const JacobianType& jacobian = 
    diffTransfo->GetJacobian( m_FixedImageSamples[sampleNumber].FixedImagePointValue, 
      this->m_MovingImages[indexImage]->timePoint, false, jacoColumns);

  for (unsigned int index=0; index<jacoColumns.size(); index++)
    {
    unsigned int fullJacoIndex = jacoColumns[index];
    double innerProduct = 0.0;
    for ( unsigned int dim = 0; dim < FixedImageDimension; dim++ )
      {
      innerProduct += jacobian[dim][fullJacoIndex] * 
        movingImageGradientValue[dim];
      }

    // Incrementing the corresponding bin in the PDF derivative
    JointPDFValueType * ptr = derivPtr + fullJacoIndex;
    *(ptr) -= innerProduct * cubicBSplineDerivativeValue;
    }

}


// Method to reinitialize the seed of the random number generator
template < class TFixedImage, class TMovingImage  > void
  MattesMutualInformationMultipleImages<TFixedImage,TMovingImage>
::ReinitializeSeed()
{
  Statistics::MersenneTwisterRandomVariateGenerator::GetInstance()->SetSeed();
}

// Method to reinitialize the seed of the random number generator
template < class TFixedImage, class TMovingImage  > void
  MattesMutualInformationMultipleImages<TFixedImage,TMovingImage>
::ReinitializeSeed(int seed)
{
  Statistics::MersenneTwisterRandomVariateGenerator::GetInstance()->SetSeed(seed);
}

} // end namespace itk

#endif
