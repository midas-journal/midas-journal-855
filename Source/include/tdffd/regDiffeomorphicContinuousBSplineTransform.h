/*=========================================================================

  Program:   Insight Segmentation & Registration Toolkit
  Module:    $RCSfile: regDiffeomorphicContinuousBSplineTransform.h,v $
  Language:  C++
  Date:      $Date: 2008-04-11 16:28:11 $
  Version:   $Revision: 1.38 $

  Copyright (c) Insight Software Consortium. All rights reserved.
  See ITKCopyright.txt or http://www.itk.org/HTML/Copyright.htm for details.

     This software is distributed WITHOUT ANY WARRANTY; without even 
     the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR 
     PURPOSE.  See the above copyright notices for more information.

=========================================================================*/

#ifndef __regDiffeomorphicContinuousBSplineTransform_h
#define __regDiffeomorphicContinuousBSplineTransform_h

#include <iostream>
#include "itkTransform.h"
#include "itkTimeProbe.h"
#include "regTimeDiffeomorphicTransform.h"
#include "regBSplineField.h"


namespace itk
{

/** \class regDiffeomorphicContinuousBSplineTransform
   * \brief Continous Deformable transform using a BSpline representation with optimization
   *
   * This class encapsulates a deformable transform of points from one 
   * N+1-dimensional one space to another (N+1)-dimensional space.
   * The deformation field is modeled using B-splines.
   *
   * T(x,t)= x(t) + \int_{0}^t v(x,\tau) d\tau
   *
   * where x is the location at time 0, v is the velocity field and T(x,t) denotes the corresponding 
   * location at time t
   *  
   * If x is the location at time t_0:
   *  
   *  T(x,t_0, t_1)= (x + \int_{t_0}^{t_1} v(x,\tau) d\tau,t1) 
   * 
   *  is the corresponding location at time t_1
   *
   * \ingroup Transforms
   */
template <
  class TScalarType = double,          // Data type for scalars
  unsigned int NDimensions = 4,        // Number of space dimensions 
  unsigned int VSplineOrder = 3 >      // Spline order
class ITK_EXPORT DiffeomorphicContinuousBSplineTransform : 
    public TimeDiffeomorphicTransform< TScalarType, NDimensions, NDimensions >
{
public:


  /** Standard class typedefs. */
  typedef DiffeomorphicContinuousBSplineTransform                         Self;
  typedef TimeDiffeomorphicTransform< TScalarType, NDimensions, NDimensions > Superclass;
  typedef SmartPointer<Self>                                 Pointer;
  typedef SmartPointer<const Self>                           ConstPointer;

  /** New macro for creation of through the object factory.*/
  itkNewMacro( Self );

  /** Run-time type information (and related methods). */
  itkTypeMacro( regDiffeomorphicContinuousBSplineTransform, TimeDiffeomorphicTransform );

  /** Name of the transform. */
  virtual std::string GetTransformTypeAsString() const;

  /** Dimension of the space domain space. */
  itkStaticConstMacro(SpaceDimension, unsigned int, NDimensions-1);

  /** Dimension of the spatio-temporal space. */
  itkStaticConstMacro(SpacePlusTimeDimension, unsigned int, NDimensions);

  /** The BSpline order. */
  itkStaticConstMacro(SplineOrder, unsigned int, VSplineOrder);

  /** Standard scalar type for this class. */
  typedef typename Superclass::ScalarType ScalarType;

  /** Standard parameters container. */
  typedef typename Superclass::ParametersType ParametersType;

  /** Standard and Sparse Jacobian containers. */
  typedef typename Superclass::JacobianType JacobianType;
  typedef typename Superclass::SparseJacobianType SparseJacobianType;
  typedef typename Superclass::SparseJacobianColumnType SparseJacobianColumnType;

  /** Spatio-temporal point type. */
  typedef Point<TScalarType, itkGetStaticConstMacro(SpaceDimension)> SpacePointType;

  /** Standard coordinate point type for this class. */
  typedef Point<TScalarType,
    itkGetStaticConstMacro(SpacePlusTimeDimension)> InputPointType;
  typedef Point<TScalarType,
    itkGetStaticConstMacro(SpacePlusTimeDimension)> OutputPointType;

  typedef TScalarType TimePointType;

  void SetParameters(const ParametersType & parameters);

  void SetFixedParameters(const ParametersType & parameters);

  void SetParametersByValue(const ParametersType & parameters);

  void SetIdentity();

  /** Get the Transformation Parameters. */
  virtual const ParametersType& GetParameters(void) const;

  /** Get the Transformation Fixed Parameters. */
  virtual const ParametersType& GetFixedParameters(void) const;


  /** Set/Get Velocity computation object. */
  typedef BSplineField<ScalarType, SpacePlusTimeDimension, SplineOrder>  VelocityType;
  typedef typename VelocityType::Pointer VelocityPointerType;
  itkGetObjectMacro( Velocity, VelocityType );
  itkSetObjectMacro( Velocity, VelocityType );

  /** Velocity-related types. */
  typedef typename VelocityType::RegionType RegionType;
  typedef typename VelocityType::SpacingType SpacingType;
  typedef typename VelocityType::OriginType OriginType;
  typedef typename VelocityType::DirectionType DirectionType;
  typedef typename VelocityType::ImageType ImageType;
  typedef typename VelocityType::ImagePointer ImagePointer;
  typedef typename VelocityType::SizeType SizeType;
  typedef typename VelocityType::PixelType PixelType;
  typedef typename VelocityType::ClassicalJacobianType ClassicalJacobianType; 

  /** Velocity-related functions. */
  void SetGridRegion( const RegionType& region );
  RegionType GetGridRegion(); 
  void SetGridSpacing( const SpacingType& spacing );
  SpacingType GetGridSpacing();
  void SetGridDirection( const DirectionType & direction );
  void SetGridOrigin( const OriginType& origin );
  OriginType GetGridOrigin();
  const ImagePointer * GetCoefficientImage() const;
  
  /** Typedef of the bulk transform. */
  typedef Transform<ScalarType,itkGetStaticConstMacro(SpaceDimension),
    itkGetStaticConstMacro(SpaceDimension)> BulkTransformType;
  typedef typename BulkTransformType::ConstPointer BulkTransformPointer;

  /** This method specifies the bulk transform to be applied. 
   * The default is the identity transform.
   */
  itkSetConstObjectMacro( BulkTransform, BulkTransformType );
  itkGetConstObjectMacro( BulkTransform, BulkTransformType );

  /** Transform points by a BSpline deformable transformation. */
  OutputPointType  TransformPoint(const InputPointType  &point, const TimePointType &endTimePoint ) const;

  /** Compute the Jacobian Matrix (w.r.t. parameters) of the transformation at one point */
  const JacobianType& GetJacobian(const InputPointType  &point,
				  const TimePointType &endTimePoint ) const;
  
  // This method updates the current jacobian at the time of the 
  //  input point to the final time, without resetting the jacobian
  void GetIncrementalSparseJacobian( const InputPointType & point,
             OutputPointType & outPoint,
				     const TimePointType & endTimePoint, 
				     SparseJacobianType & jaco, 
				     bool & valid ) const;

  /** Return the number of parameters that completely define the Transfom */
  virtual unsigned int GetNumberOfParameters(void) const;

  /** Return the number of parameters per dimension */
  unsigned int GetNumberOfParametersPerDimension(void) const;

  /** Set/Get time step size */
  itkSetMacro(TimeStep, double);
  itkGetConstMacro(TimeStep, double);

  /** Set/Get minimum time step size */
  itkSetMacro(MinimumTimeStep, double);
  itkGetConstMacro(MinimumTimeStep, double);
  
  /** This class supports incremental jacobian. */
  bool SupportIncrementalJacobian()
  {
    return true;
  }
  
  /** Transform point and compute volume change at the same time. */
    void TransformPointAndGetPhysicalJacobian( const InputPointType & point, const TimePointType & endTimePoint, ClassicalJacobianType & output, OutputPointType& outputPoint) const;
    
protected:
  /** Print contents of an BSplineDeformableTransform. */
  void PrintSelf(std::ostream &os, Indent indent) const;

  DiffeomorphicContinuousBSplineTransform();
  virtual ~DiffeomorphicContinuousBSplineTransform();

  /** Velocity object. */
  VelocityPointerType m_Velocity;

  /** Internal parameters buffer. */
  ParametersType          m_InternalParametersBuffer;


private:
  DiffeomorphicContinuousBSplineTransform(const Self&); //purposely not implemented
  void operator=(const Self&); //purposely not implemented

  /** The bulk transform. */
  BulkTransformPointer  m_BulkTransform;

  mutable double m_TimeStep;
  double m_MinimumTimeStep;

  /** Total jacobian. */
  mutable JacobianType m_TotalJacobian; 
  
  /** Cache variables for advoiding frequent reallocation. */
  mutable SparseJacobianColumnType incrJaco_columnOut;
  mutable typename VelocityType::WeightsType incrJaco_weights;
  mutable typename VelocityType::ParameterIndexArrayType incrJaco_indices;
  mutable typename VelocityType::ClassicalJacobianType incrJaco_physJacobian;
  mutable SpacePointType incrJaco_velocity;

}; //class DiffeomorphicContinuousBSplineTransform


}  // namespace itk


#if ITK_TEMPLATE_TXX
# include "regDiffeomorphicContinuousBSplineTransform.txx"
#endif

#endif /* __regDiffeomorphicContinuousBSplineTransform_h */
