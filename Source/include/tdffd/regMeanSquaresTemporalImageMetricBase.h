/*=========================================================================

  Program:   Insight Segmentation & Registration Toolkit
  Module:    $RCSfile: regMeanSquaresTemporalImageMetricBase.h,v $
  Language:  C++
  Date:      $Date: 2008-02-03 04:05:28 $
  Version:   $Revision: 1.31 $

  Copyright (c) Insight Software Consortium. All rights reserved.
  See ITKCopyright.txt or http://www.itk.org/HTML/Copyright.htm for details.

     This software is distributed WITHOUT ANY WARRANTY; without even 
     the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR 
     PURPOSE.  See the above copyright notices for more information.

=========================================================================*/
#ifndef __regMeanSquaresTemporalImageMetricBase_h
#define __regMeanSquaresTemporalImageMetricBase_h

// First make sure that the configuration is available.
// This line can be removed once the optimized versions
// gets integrated into the main directories.
#include "itkConfigure.h"

#include "itkPoint.h"
#include "itkCovariantVector.h"
#include "regTemporalImageMetric.h"
#include "itkImageRandomConstIteratorWithIndex.h"

namespace itk
{
/** \class MeanSquaresTemporalImageMetricBase
 * \brief Computes similarity  ..... To be completed
 *
 * \ingroup RegistrationMetrics
 */
template < class TImage> 
class ITK_EXPORT MeanSquaresTemporalImageMetricBase : 
    public TemporalImageMetric< TImage>
{
public:

  /** Standard class typedefs. */
  typedef MeanSquaresTemporalImageMetricBase    Self;
  typedef TemporalImageMetric<TImage>  Superclass;

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

  /** Run-time type information (and related methods). */
  itkTypeMacro(MeanSquaresTemporalImageMetricBase, TemporalImageMetric);
  
  /** Types transferred from the transform class */
  typedef typename Superclass::TransformType::SparseJacobianType SparseJacobianType;

  /** Initialize the Metric by making sure that all the components
   *  are present and plugged together correctly     */
  virtual void Initialize(void) throw ( ExceptionObject );
  
  /** Get the derivatives of the match measure. */
  virtual void GetDerivative( const TransformParametersType & parameters,
                      DerivativeType & derivative ) const = 0;

  /**  Get the value for single valued optimizers. */
  virtual MeasureType GetValue( const TransformParametersType & parameters ) const = 0;

  /**  Get value and derivatives for multiple valued optimizers. */
  virtual void GetValueAndDerivative( const TransformParametersType & parameters,
                              MeasureType& Value, DerivativeType& Derivative ) const = 0;

  void GetSampleContributionOnTimeInterval(InputPointType inputPoint, double finalTime, RealType fixedValue, 
                                           MeasureType & value, DerivativeType  & derivative, SparseJacobianType& jacobian,
                                           unsigned long & numberOfValidPropagatedPoints) const;

protected:
  MeanSquaresTemporalImageMetricBase(){};
  virtual ~MeanSquaresTemporalImageMetricBase() {};
  
  
private:
  
  MeanSquaresTemporalImageMetricBase(const Self&); //purposely not implemented
  void operator=(const Self&); //purposely not implemented

};

} // end namespace itk

#ifndef ITK_MANUAL_INSTANTIATION
#include "regMeanSquaresTemporalImageMetricBase.txx"
#endif

#endif


