/*=========================================================================

  Program:   Insight Segmentation & Registration Toolkit
  Module:    $RCSfile: regTemporalImageMetric.h,v $
  Language:  C++
  Date:      $Date: 2007-11-12 20:00:37 $
  Version:   $Revision: 1.25 $

  Copyright (c) Insight Software Consortium. All rights reserved.
  See ITKCopyright.txt or http://www.itk.org/HTML/Copyright.htm for details.

     This software is distributed WITHOUT ANY WARRANTY; without even 
     the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR 
     PURPOSE.  See the above copyright notices for more information.

=========================================================================*/
#ifndef __regTemporalImageMetric_h
#define __regTemporalImageMetric_h

// First make sure that the configuration is available.
// This line can be removed once the optimized versions
// gets integrated into the main directories.
#include "itkConfigure.h"

#include "itkPoint.h"
#include "itkIndex.h"
#include "itkImageBase.h"
#include "regTimeDiffeomorphicTransform.h"
#include "itkInterpolateImageFunction.h"
#include "itkSingleValuedCostFunction.h"
#include "itkExceptionObject.h"
#include "itkGradientRecursiveGaussianImageFilter.h"
#include "itkSpatialObject.h"

#include "itkImageRandomConstIteratorWithIndex.h"
#include "itkImageFileWriter.h"
#include "itkBSplineInterpolateImageFunction.h"

#include "vtkPolyData.h"

#include "regDiffeomorphicContinuousBSplineTransform.h"

namespace itk
{
  
/** \class TemporalImageMetric
 * \brief Computes similarity between regions of two images.
 *
 * This Class is templated over the type of images.
 * It expects a Transform and an Interpolator to be plugged in.
 * This particular class is the base class for a hierarchy of 
 * similarity metrics.
 *
 * This class computes a value that measures the similarity 
 * between the image a t=0 and the transformed  image.
 * The Interpolator is used to compute intensity values on 
 * non-grid positions resulting from mapping points through 
 * the Transform.
 * 
 *
 * \ingroup RegistrationMetrics
 *
 */

template <class TImage> 
class ITK_EXPORT TemporalImageMetric : public SingleValuedCostFunction 
{
public:
  /** Standard class typedefs. */
  typedef TemporalImageMetric           Self;
  typedef SingleValuedCostFunction     Superclass;
  typedef SmartPointer<Self>           Pointer;
  typedef SmartPointer<const Self>     ConstPointer;

  /** Type used for representing point components  */
  typedef typename Superclass::ParametersValueType CoordinateRepresentationType;

  /** Run-time type information (and related methods). */
  itkTypeMacro(TemporalImageMetric, SingleValuedCostFunction);

  /** New macro for creation of through the object factory.*/
  itkNewMacro( Self );

  virtual MeasureType GetValue( const ParametersType & parameters ) const {return 0.;};
  virtual void GetDerivative( const ParametersType & parameters,
                              DerivativeType & derivative ) const {};
  
  /**  Type of the Image. */
  typedef TImage                                ImageType;
  typedef typename ImageType::ConstPointer      ImageConstPointer;
  typedef typename ImageType::RegionType        ImageRegionType;
  typedef typename ImageType::IndexType			ImageIndexType;
  typedef typename TImage::PixelType            ImagePixelType;

  /** Constants for the image dimensions */
  itkStaticConstMacro(ImageDimension,
                      unsigned int,
                      TImage::ImageDimension);
  itkStaticConstMacro(SpaceImageDimension,
                      unsigned int,
                      ImageDimension-1);
  
  /**  Type of the Transform Base class */
  typedef TimeDiffeomorphicTransform<CoordinateRepresentationType,
    itkGetStaticConstMacro(ImageDimension),
    itkGetStaticConstMacro(ImageDimension)>
    TransformType;

  typedef typename TransformType::Pointer            TransformPointer;
  typedef typename TransformType::InputPointType     InputPointType;
  typedef typename TransformType::OutputPointType    OutputPointType;
  typedef typename TransformType::ParametersType     TransformParametersType;
  typedef typename TransformType::JacobianType       TransformJacobianType;

  /**  Type of the Interpolator Base class */
  typedef InterpolateImageFunction<
    ImageType,
    CoordinateRepresentationType > InterpolatorType;
  typedef BSplineInterpolateImageFunction<
  ImageType,
  CoordinateRepresentationType > InterpolatorTypeBS;

  /** Gaussian filter to compute the gradient of the Moving Image */
  typedef typename NumericTraits<ImagePixelType>::RealType
    RealType;
  typedef CovariantVector<RealType,
    itkGetStaticConstMacro(ImageDimension)>
    GradientPixelType;
  typedef Image<GradientPixelType,
    itkGetStaticConstMacro(ImageDimension)>
    GradientImageType;
  typedef SmartPointer<GradientImageType>            GradientImagePointer;
  typedef GradientRecursiveGaussianImageFilter< ImageType,
    GradientImageType >
    GradientImageFilterType;
  typedef typename GradientImageFilterType::Pointer  GradientImageFilterPointer;

  /** Interpolator types. */
  typedef typename InterpolatorType::Pointer         InterpolatorPointer;
  typedef typename InterpolatorTypeBS::Pointer       InterpolatorBSPointer;

  /**  Type for the mask of the fixed image. Only pixels that are "inside"
       this mask will be considered for the computation of the metric */
  typedef SpatialObject< itkGetStaticConstMacro(SpaceImageDimension) >   SpaceImageMaskType;
  typedef typename  SpaceImageMaskType::Pointer							 SpaceImageMaskPointer;

  typedef typename  SpaceImageMaskType::PointType						 SpacePointType;

  /**  Type of the measure. */
  typedef typename Superclass::MeasureType                    MeasureType;

  /**  Type of the derivative. */
  typedef typename Superclass::DerivativeType                 DerivativeType;

  /**  Type of the parameters. */
  typedef typename Superclass::ParametersType                 ParametersType;

  /** Connect the Image.  */
  itkSetConstObjectMacro( Image, ImageType );

  /** Get the Image. */
  itkGetConstObjectMacro( Image, ImageType );

  /** Connect the Transform. */
  itkSetObjectMacro( Transform, TransformType );

  /** Get a pointer to the Transform.  */
  itkGetConstObjectMacro( Transform, TransformType );
 
  /** Connect the Interpolator. */
  //itkSetObjectMacro( Interpolator, InterpolatorType );
  void SetInterpolator(InterpolatorType* interpolator);
  
  /** Get a pointer to the Interpolator.  */
  itkGetConstObjectMacro( Interpolator, InterpolatorType );

  /** Get the number of pixels considered in the computation. */
  itkGetConstReferenceMacro( NumberOfValidSamples, unsigned long );

  /** Set the region over which the metric will be computed */
  itkSetMacro( ImageRegion, ImageRegionType );

  /** Get the region over which the metric will be computed */
  itkGetConstReferenceMacro( ImageRegion, ImageRegionType );
 
  /** Set/Get the fixed image mask. */
  itkSetObjectMacro( SpaceImageMask, SpaceImageMaskType );
  itkGetConstObjectMacro( SpaceImageMask, SpaceImageMaskType );

  /** Set/Get the moving image mask. */
  itkSetObjectMacro( SpaceMovingImageMask, SpaceImageMaskType );
  itkGetConstObjectMacro( SpaceMovingImageMask, SpaceImageMaskType );


  /** Set/Get gradient computation. */
  itkSetMacro( ComputeGradient, bool);
  itkGetConstReferenceMacro( ComputeGradient, bool);
  itkBooleanMacro(ComputeGradient);

  /** Computes the gradient image and assigns it to m_GradientImage */
  virtual void ComputeGradient();

  /** Get Gradient Image. */
  itkGetConstObjectMacro( GradientImage, GradientImageType );

  /** Set the parameters defining the Transform. */
  virtual void SetTransformParameters( const ParametersType & parameters ) const;

  /** Return the number of parameters required by the Transform */
  virtual unsigned int GetNumberOfParameters(void) const
  { return m_Transform->GetNumberOfParameters(); }

  /** Initialize the Metric by making sure that all the components
   *  are present and plugged together correctly     */
  virtual void Initialize(void) throw ( ExceptionObject );

  /** Sample first image domain */
  virtual void SampleReferenceImageDomain(bool denseSampling, unsigned int timeMargin=1) const;

  /** Set/Get the number of samples.*/
  itkGetMacro(NumberOfSamples, unsigned int);
  itkSetMacro(NumberOfSamples, unsigned int);

  /** Set/Get reference time. */
  itkSetMacro(ReferenceTime, double);
  itkGetMacro(ReferenceTime, double);
 
protected: 
  TemporalImageMetric();
  virtual ~TemporalImageMetric(); 
  void PrintSelf(std::ostream& os, Indent indent) const;

  bool m_Initialized;

  mutable unsigned long       m_NumberOfValidSamples;

  ImageConstPointer           m_Image;
  ImageRegionType             m_ImageRegion; 
  
  mutable TransformPointer    m_Transform;
  InterpolatorPointer         m_Interpolator;
  InterpolatorBSPointer       m_InterpolatorBS;
  

  bool                        m_ComputeGradient;
  GradientImagePointer        m_GradientImage;

  mutable SpaceImageMaskPointer   m_SpaceImageMask;
  mutable SpaceImageMaskPointer   m_SpaceMovingImageMask;

  unsigned long m_NumberOfSamples;  

  class DiffeoSpatialSample
	{
	public:
		InputPointType m_FixedPoint;
		RealType       m_FixedValue;
	};

  mutable std::vector<DiffeoSpatialSample> m_SamplesContainer;
  
   // Time of the reference frame
  double m_ReferenceTime;
  
  // Flag for detecting bspline interpolator
  bool m_InterpolatorIsBSpline;
 
  void SetMovingMaskTimeInterval(double timeInterval[2])
  {
	  m_MovingMaskTimeSpecified=true;
	  m_MovingMaskTimeInterval[0]=timeInterval[0];
	  m_MovingMaskTimeInterval[1]=timeInterval[1];
  }

  bool IsWithinMovingMask(OutputPointType& point) const
  {
	  bool skip_sample = false;
	  SpacePointType point2;
	  // FIXME We need to find a better strategy here
	  // for avoiding such conversion
	  for (unsigned int d=0; d<this->SpaceImageDimension; d++)
	  {
		  point2[d] = point[d];
	  }

	  double time=point[SpaceImageDimension];
	
	  if (m_MovingMaskTimeSpecified &&
		  ((time<m_MovingMaskTimeInterval[0])||(time>m_MovingMaskTimeInterval[1])) )
	  {
		  // Out of the time interval so don't skip this sample
		  skip_sample = false;
	  }
	  else if( this->m_SpaceMovingImageMask && !this->m_SpaceMovingImageMask->IsInside( point2 ) )
	  {
		  skip_sample = true;
	  }

	  return skip_sample;

  }

private:
  TemporalImageMetric(const Self&); //purposely not implemented
  void operator=(const Self&); //purposely not implemented

  double m_MovingMaskTimeInterval[2];
  bool m_MovingMaskTimeSpecified;
   

};

} // end namespace itk

#ifndef ITK_MANUAL_INSTANTIATION
#include "regTemporalImageMetric.txx"
#endif

#endif
