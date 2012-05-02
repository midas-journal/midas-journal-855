/*=========================================================================

  Program:   Insight Segmentation & Registration Toolkit
  Module:    $RCSfile: regContinuousBSplineTransform.h,v $
  Language:  C++
  Date:      $Date: 2008-04-11 16:28:11 $
  Version:   $Revision: 1.38 $

  Copyright (c) Insight Software Consortium. All rights reserved.
  See ITKCopyright.txt or http://www.itk.org/HTML/Copyright.htm for details.

     This software is distributed WITHOUT ANY WARRANTY; without even 
     the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR 
     PURPOSE.  See the above copyright notices for more information.

=========================================================================*/

#ifndef __regContinuousBSplineTransform_h
#define __regContinuousBSplineTransform_h

#include <iostream>
#include "itkTransform.h"
#include "itkTimeProbe.h"
#include "regTimeDiffeomorphicTransform.h"
#include "regBSplineField.h"


namespace itk
{

/** \class regContinuousBSplineTransform
  */
template <
  class TScalarType = double,          // Data type for scalars
  unsigned int NDimensions = 4,        // Number of space dimensions 
  unsigned int VSplineOrder = 3 >      // Spline order
class ITK_EXPORT ContinuousBSplineTransform : 
    public TimeDiffeomorphicTransform< TScalarType, NDimensions, NDimensions >
{
public:


  /** Standard class typedefs. */
  typedef ContinuousBSplineTransform                         Self;
  typedef TimeDiffeomorphicTransform< TScalarType, NDimensions, NDimensions > Superclass;
  typedef SmartPointer<Self>                                 Pointer;
  typedef SmartPointer<const Self>                           ConstPointer;

  /** New macro for creation of through the object factory.*/
  itkNewMacro( Self );

  /** Run-time type information (and related methods). */
  itkTypeMacro( regContinuousBSplineTransform, TimeDiffeomorphicTransform );

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
	
  /** Set/Get BSplineField computation object. */
  typedef BSplineField<ScalarType, SpacePlusTimeDimension, SplineOrder>  BSplineFieldType;
  typedef typename BSplineFieldType::Pointer BSplineFieldPointerType;
  itkGetObjectMacro( BSplineField, BSplineFieldType );
  itkSetObjectMacro( BSplineField, BSplineFieldType );

  /** BSplineField-related types. */
  typedef typename BSplineFieldType::RegionType RegionType;
  typedef typename BSplineFieldType::SpacingType SpacingType;
  typedef typename BSplineFieldType::OriginType OriginType;
  typedef typename BSplineFieldType::DirectionType DirectionType;
  typedef typename BSplineFieldType::ImageType ImageType;
  typedef typename BSplineFieldType::ImagePointer ImagePointer;
  typedef typename BSplineFieldType::SizeType SizeType;
  typedef typename BSplineFieldType::PixelType PixelType;

  /** BSplineField-related functions. */
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

  /** Return the number of parameters that completely define the Transfom */
  virtual unsigned int GetNumberOfParameters(void) const;

  /** Return the number of parameters per dimension */
  unsigned int GetNumberOfParametersPerDimension(void) const;

  void GetSparseJacobian( const InputPointType & point, const TimePointType & endTimePoint, SparseJacobianType& jaco ) const;

  /** This class does not support incremental jacobian. */
  bool SupportIncrementalJacobian()
  {
    return false;
  }
  
protected:
  /** Print contents of an BSplineDeformableTransform. */
  void PrintSelf(std::ostream &os, Indent indent) const;

  ContinuousBSplineTransform();
  virtual ~ContinuousBSplineTransform();

  /** BSplineField object. */
  BSplineFieldPointerType m_BSplineField;

  /** Internal parameters buffer. */
  ParametersType          m_InternalParametersBuffer;
  
private:
  ContinuousBSplineTransform(const Self&); //purposely not implemented
  void operator=(const Self&); //purposely not implemented

  /** The bulk transform. */
  BulkTransformPointer  m_BulkTransform;
  
  /** Cache variables for advoiding frequent reallocation. */
  mutable SparseJacobianColumnType incrJaco_columnOut;
  mutable typename BSplineFieldType::WeightsType incrJaco_weights;
  mutable typename BSplineFieldType::ParameterIndexArrayType incrJaco_indices;
  mutable SpacePointType incrJaco_velocity;
 
}; //class ContinuousBSplineTransform


}  // namespace itk


#if ITK_TEMPLATE_TXX
# include "regContinuousBSplineTransform.txx"
#endif

#endif /* __regContinuousBSplineTransform_h */
