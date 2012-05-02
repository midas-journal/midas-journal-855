/*=========================================================================

  Program:   Insight Segmentation & Registration Toolkit
  Module:    $RCSfile: regMeanSquaresTemporalImageMetricSequential2.h,v $
  Language:  C++
  Date:      $Date: 2008-02-03 04:05:28 $
  Version:   $Revision: 1.31 $

  Copyright (c) Insight Software Consortium. All rights reserved.
  See ITKCopyright.txt or http://www.itk.org/HTML/Copyright.htm for details.

     This software is distributed WITHOUT ANY WARRANTY; without even 
     the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR 
     PURPOSE.  See the above copyright notices for more information.

=========================================================================*/
#ifndef __regMeanSquaresTemporalImageMetricSequential2_h
#define __regMeanSquaresTemporalImageMetricSequential2_h

// First make sure that the configuration is available.
// This line can be removed once the optimized versions
// gets integrated into the main directories.
#include "itkConfigure.h"

#include "itkPoint.h"
#include "itkCovariantVector.h"
#include "regTemporalImageMetric.h"
#include "itkImageRandomConstIteratorWithIndex.h"
#include "regMeanSquaresTemporalImageMetric.h"

namespace itk
{
/** \class MeanSquaresTemporalImageMetricSequential2
 * \brief Computes similarity ..... To be completed
 *
 * \ingroup RegistrationMetrics
 */
template < class TImage> 
class ITK_EXPORT MeanSquaresTemporalImageMetricSequential2 : 
    public MeanSquaresTemporalImageMetricBase< TImage>
{
public:

  /** Standard class typedefs. */
  typedef MeanSquaresTemporalImageMetricSequential2    Self;
  typedef MeanSquaresTemporalImageMetricBase<TImage>  Superclass;

  typedef SmartPointer<Self>         Pointer;
  typedef SmartPointer<const Self>   ConstPointer;

    /** Types transferred from the base class */
  typedef typename Superclass::RealType                 RealType;
  typedef typename Superclass::TransformType            TransformType;
  typedef typename Superclass::TransformPointer         TransformPointer;
  typedef typename Superclass::TransformParametersType  TransformParametersType;
  typedef typename Superclass::TransformJacobianType    TransformJacobianType;
  typedef typename Superclass::GradientPixelType        GradientPixelType;
  typedef typename Superclass::GradientImageType        GradientImageType;
  typedef typename Superclass::InputPointType           InputPointType;
  typedef typename Superclass::OutputPointType          OutputPointType;
  typedef typename Superclass::SpacePointType           SpacePointType;

  typedef typename Superclass::MeasureType              MeasureType;
  typedef typename Superclass::DerivativeType           DerivativeType;
  typedef typename Superclass::ImageType                ImageType;
  typedef typename Superclass::SpaceImageMaskType       SpaceImageMaskType;
  typedef typename Superclass::ImageConstPointer        ImageConstPointer;


  /** Method for creation through the object factory. */
  itkNewMacro(Self);
 
  /** Run-time type information (and related methods). */
  itkTypeMacro(MeanSquaresTemporalImageMetricSequential2, MeanSquaresTemporalImageMetricBase);

  /** Initialize the Metric by making sure that all the components
   *  are present and plugged together correctly     */
  virtual void Initialize(void) throw ( ExceptionObject );
  
  
  /** Types transferred from the transform class */
  typedef typename TransformType::SparseJacobianType SparseJacobianType;


  /** Get the derivatives of the match measure. */
  void GetDerivative( const TransformParametersType & parameters,
                      DerivativeType & derivative ) const;

  /**  Get the value for single valued optimizers. */
  MeasureType GetValue( const TransformParametersType & parameters ) const;

  /**  Get value and derivatives for multiple valued optimizers. */
  void GetValueAndDerivative( const TransformParametersType & parameters,
                              MeasureType& Value, DerivativeType& Derivative ) const;

protected:
  MeanSquaresTemporalImageMetricSequential2();
  virtual ~MeanSquaresTemporalImageMetricSequential2() {};

  // Time increment for sequential metrics
  double m_sequentialTimeIncrement;

private:
  MeanSquaresTemporalImageMetricSequential2(const Self&); //purposely not implemented
  void operator=(const Self&); //purposely not implemented

};

} // end namespace itk

#ifndef ITK_MANUAL_INSTANTIATION
#include "regMeanSquaresTemporalImageMetricSequential2.txx"
#endif

#endif


