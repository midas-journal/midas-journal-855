/*=========================================================================

  Program:   Insight Segmentation & Registration Toolkit
  Module:    $RCSfile: regHybridSeqNonSeqMetric.h,v $
  Language:  C++
  Date:      $Date: 2008-02-03 04:05:28 $
  Version:   $Revision: 1.31 $

  Copyright (c) Insight Software Consortium. All rights reserved.
  See ITKCopyright.txt or http://www.itk.org/HTML/Copyright.htm for details.

     This software is distributed WITHOUT ANY WARRANTY; without even 
     the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR 
     PURPOSE.  See the above copyright notices for more information.

=========================================================================*/
#ifndef __regHybridSeqNonSeqMetric_h
#define __regHybridSeqNonSeqMetric_h

// First make sure that the configuration is available.
// This line can be removed once the optimized versions
// gets integrated into the main directories.
#include "itkConfigure.h"

#include "regTemporalImageMetric.h"
#include "itkImageRandomConstIteratorWithIndex.h"

namespace itk
{
/** \class HybridSeqNonSeqMetric
 * \brief Computes similarity  ..... To be completed
 *
 * \ingroup RegistrationMetrics
 */
template < class TImage> 
class ITK_EXPORT HybridSeqNonSeqMetric : 
    public TemporalImageMetric< TImage>
{
public:

  /** Standard class typedefs. */
  typedef HybridSeqNonSeqMetric    Self;
  typedef SingleValuedCostFunction  Superclass;

  typedef SmartPointer<Self>         Pointer;
  typedef SmartPointer<const Self>   ConstPointer;

  typedef TemporalImageMetric<TImage> MetricType;
  typedef SmartPointer<MetricType>    MetricPointer;

    /** Types transferred from the base class */
  typedef typename MetricType::InterpolatorPointer		InterpolatorPointer;
  typedef typename MetricType::TransformType            TransformType;
  typedef typename MetricType::TransformPointer         TransformPointer;
  typedef typename MetricType::TransformParametersType  TransformParametersType;
  typedef typename MetricType::TransformJacobianType    TransformJacobianType;

  typedef typename MetricType::MeasureType              MeasureType;
  typedef typename MetricType::DerivativeType           DerivativeType;
  typedef typename MetricType::ImageType                ImageType;
  typedef typename MetricType::ImageRegionType			ImageRegionType;
  typedef typename MetricType::ImageConstPointer        ImageConstPointer;
  typedef typename MetricType::SpaceImageMaskPointer    SpaceImageMaskPointer;

  /** Method for creation through the object factory. */
  itkNewMacro(Self);
 
  /** Run-time type information (and related methods). */
  itkTypeMacro(HybridSeqNonSeqMetric, SingleValuedCostFunction);

  /** Initialize the Metric by making sure that all the components
   *  are present and plugged together correctly     */
  virtual void Initialize(void) throw ( ExceptionObject )
  {
	  // Set image
	  m_MetricNonSequential->SetImage( this->m_Image );
	  m_MetricSequential->SetImage( this->m_Image );

	  m_MetricNonSequential->SetImageRegion( this->m_ImageRegion );
	  m_MetricSequential->SetImageRegion( this->m_ImageRegion );

	  m_MetricNonSequential->SetInterpolator( this->m_Interpolator );
	  m_MetricSequential->SetInterpolator( this->m_Interpolator );

	  //const unsigned int numberOfTimePoints = this->m_ImageRegion.GetSize()[this->SpaceImageDimension];
	  //unsigned long nSamplesSeq = this->m_NumberOfSamples * (numberOfTimePoints-1);

	  m_MetricNonSequential->SetNumberOfSamples( this->m_NumberOfSamples );
	  m_MetricSequential->SetNumberOfSamples( this->m_NumberOfSamples );

	  m_MetricNonSequential->SetSpaceImageMask( this->m_SpaceImageMask );
	  m_MetricSequential->SetSpaceImageMask( this->m_SpaceImageMask );

	  m_MetricNonSequential->SetSpaceMovingImageMask( this->m_SpaceMovingImageMask );
	  m_MetricSequential->SetSpaceMovingImageMask( this->m_SpaceMovingImageMask );

	  m_MetricNonSequential->SetReferenceTime( this->m_ReferenceTime );
	  m_MetricSequential->SetReferenceTime( this->m_ReferenceTime );

	  m_MetricNonSequential->SetTransform( this->m_Transform );
	  m_MetricSequential->SetTransform( this->m_Transform );

	  m_MetricSequential->Initialize();
	  m_MetricNonSequential->Initialize();
  }
	
  /** Get the derivatives of the match measure. */
  void GetDerivative( const TransformParametersType & parameters,
                      DerivativeType & derivative ) const
  {
	  double tempValue;
	  this->GetValueAndDerivative(parameters, tempValue, derivative);

  }

  /**  Get the value for single valued optimizers. */
  MeasureType GetValue( const TransformParametersType & parameters ) const
  {
	  double tempValue;

	  // computing derivative of first term
	  tempValue = m_Weight * m_MetricSequential->GetValue(parameters);

	  // computing derivative of second term
	  tempValue += m_MetricNonSequential->GetValue(parameters);

	  return tempValue;

  }

  /**  Get value and derivatives for multiple valued optimizers. */
  void GetValueAndDerivative( const TransformParametersType & parameters,
                              MeasureType& value, DerivativeType& derivative ) const
  {
	  double tempValue;

	  // Adjusting size of internal derivative if required
	  if (m_InternalDerivative.Size()!=parameters.Size())
	  {
		  m_InternalDerivative.SetSize(parameters.Size());
	  }

	  // computing derivative of first term
	  m_MetricSequential->GetValueAndDerivative(parameters, tempValue, m_InternalDerivative);

	  // computing derivative of second term
	  m_MetricNonSequential->GetValueAndDerivative(parameters, value, derivative);

	  // summing the two
	  for (unsigned int p=0; p<parameters.Size();p++)
	  {
		  derivative[p] += m_InternalDerivative[p] * m_Weight;
	  }
	  value += tempValue * m_Weight;

  }

  /**  Get value and derivatives for multiple valued optimizers. */
  void ComputeWeight( const TransformParametersType & parameters )
  {
	  double seqValue=m_MetricSequential->GetValue(parameters);
	  double non_seqValue=m_MetricNonSequential->GetValue(parameters);

	  m_Weight = non_seqValue / seqValue;

	  std::cout << "Seq/Non-seq weight is " << m_Weight << std::endl;

  }

  itkSetMacro(Weight, double);
  itkGetMacro(Weight, double);

  itkSetObjectMacro(MetricNonSequential, MetricType);
  itkGetObjectMacro(MetricNonSequential, MetricType);

  itkSetObjectMacro(MetricSequential, MetricType);
  itkGetObjectMacro(MetricSequential, MetricType);

protected:
  HybridSeqNonSeqMetric()
  {
	  m_MetricSequential=0;
	  m_MetricNonSequential=0;
	  m_Weight=4.;

  }
  virtual ~HybridSeqNonSeqMetric() {};

  MetricPointer m_MetricSequential;
  MetricPointer m_MetricNonSequential;
  
  double m_Weight;
  mutable DerivativeType m_InternalDerivative;

private:
  HybridSeqNonSeqMetric(const Self&); //purposely not implemented
  void operator=(const Self&); //purposely not implemented

};

} // end namespace itk


#endif


