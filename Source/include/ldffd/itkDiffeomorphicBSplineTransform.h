/*=========================================================================

  Program:   Insight Segmentation & Registration Toolkit
  Module:    $RCSfile: itkDiffeomorphicBSplineTransform.h,v $
  Language:  C++
  Date:      $Date: 2008-06-20 15:25:58 $
  Version:   $Revision: 1.2 $

  Copyright (c) Insight Software Consortium. All rights reserved.
  See ITKCopyright.txt or http://www.itk.org/HTML/Copyright.htm for details.

     This software is distributed WITHOUT ANY WARRANTY; without even 
     the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR 
     PURPOSE.  See the above copyright notices for more information.

=========================================================================*/

#ifndef __itkDiffeomorphicBSplineTransform_h
#define __itkDiffeomorphicBSplineTransform_h

#include <iostream>
#include "itkTransform.h"
#include "itkExceptionObject.h"
#include "itkMatrix.h"
#include "itkBSplineDeformableTransform.h"
#include "BSplineDeformableTransformOpt.h"

namespace itk
{

/** \brief Translation transformation of a vector space (e.g. space coordinates)
 *
 * The same functionality could be obtained by using the Affine tranform,
 * but with a large difference in performace.
 *
 * \ingroup Transforms
 */
template <
  class TScalarType=double,          // Data type for scalars (float or double)
  unsigned int NDimensions=3,        // Number of dimensions
  unsigned int BSplineOrder=3>         // BSpline Order
class ITK_EXPORT DiffeomorphicBSplineTransform : 
          public Transform< TScalarType, NDimensions, NDimensions >
{
public:
  /** Standard class typedefs. */
  typedef DiffeomorphicBSplineTransform Self;
  typedef Transform< TScalarType, NDimensions, NDimensions > Superclass;
  typedef SmartPointer<Self>        Pointer;
  typedef SmartPointer<const Self>  ConstPointer;
  
  /** New macro for creation of through the object factory.*/
  itkNewMacro( Self );
  
  /** Run-time type information (and related methods). */
  itkTypeMacro( DiffeomorphicBSplineTransform, Transform );
  
  /** Dimension of the domain space. */
  itkStaticConstMacro(SpaceDimension, unsigned int, NDimensions);
  
  /** Standard scalar type for this class. */
  typedef typename Superclass::ScalarType ScalarType;

  /** Standard parameters container. */
  typedef typename Superclass::ParametersType ParametersType;

  /** Standard Jacobian container. */
  typedef typename Superclass::JacobianType JacobianType;
  
  /** Standard coordinate point type for this class. */
  typedef Point<TScalarType, itkGetStaticConstMacro(SpaceDimension)> InputPointType;
  typedef Point<TScalarType, itkGetStaticConstMacro(SpaceDimension)> OutputPointType;

  /** Internal BSpline type. */
  typedef BSplineDeformableTransformOpt<TScalarType, SpaceDimension, BSplineOrder>
    InternalTransformType;
  typedef typename InternalTransformType::Pointer
    InternalTransformPointerType;

  /** Internal jacobian type. */
  typedef vnl_matrix<double> InternalJacobianType;
  
  /** Compute the mapped coordinate of an input point. */
  OutputPointType TransformPoint(const InputPointType  &point ) const;
  OutputPointType TransformPoint(const InputPointType  &point, unsigned int timeStep ) const;
  
  bool TransformPoint(const InputPointType  &point, std::vector<OutputPointType>& outputPoints) const;

  const std::vector<OutputPointType> & TransformPoint(const InputPointType  &point, bool& valid) const;
  
  /** This methods sets the parameters for the transform
   * value specified by the user (by value or not). */
  void SetParameters(const ParametersType & parameters);
  void SetParametersByValue(const ParametersType & parameters);
  void WrapParameters();

  /** Get the Transformation Parameters. */
  virtual const ParametersType& GetParameters(void) const;
  
  /** Compute the Jacobian Matrix of the transformation at one point */
  virtual const JacobianType & GetJacobian(const InputPointType  &point ) const;
  
  /** Set the parameters to the IdentityTransform. */
  void SetIdentity(void);
  
  /** Return the number of parameters that completely define the Transfom  */
  virtual unsigned int GetNumberOfParameters(void) const 
  { return m_GridRegion.GetNumberOfPixels() * m_NumberOfTimeSteps * SpaceDimension; }
  
  /** Get number of parameters by transform. */
  virtual unsigned int GetNumberOfParametersByTransform(void) const 
  { return m_GridRegion.GetNumberOfPixels() * SpaceDimension; }
  
  /** Set/Get time step. */
  itkSetMacro(NumberOfTimeSteps, unsigned int);
  itkGetConstMacro(NumberOfTimeSteps, unsigned int);
  
  /** Typedefs from internal transforms. */
  typedef typename InternalTransformType::SpacingType SpacingType;
  typedef typename InternalTransformType::OriginType OriginType;
  typedef typename InternalTransformType::RegionType RegionType;
  typedef typename InternalTransformType::WeightsType WeightsType;
  typedef typename InternalTransformType::ParameterIndexArrayType 
    ParameterIndexArrayType;

  /** Set/Get for configuring internal transformations. */
  itkSetMacro(GridSpacing, SpacingType);
  itkSetMacro(GridOrigin, OriginType);
  itkSetMacro(GridRegion, RegionType);
  itkGetConstMacro(GridSpacing, SpacingType);
  itkGetConstMacro(GridOrigin, OriginType);
  itkGetConstMacro(GridRegion, RegionType);

  const ParametersType & GetFixedParameters( void ) const;
  void SetFixedParameters( const ParametersType & parameters );

  // Reset all internal jacobians
  void ResetAllPhysicalJacobians();

  /** Access the vector of internal transforms. */
  virtual const std::vector<InternalTransformPointerType> & GetInternalTransforms()
  {return m_InternalTransforms;}
  
  const InternalTransformType* GetInternalTransform(unsigned int index)
  {return m_InternalTransforms[index].GetPointer();}

  void InsertTransform(unsigned int pos, bool reshapeParameters=true);
  
  const std::vector<InternalJacobianType> &
    AccumulatePhysicalJacobiansInTimeInterval(const InputPointType & point_time_t1, unsigned int time1, unsigned int time2) const;

  const JacobianType & 
    GetJacobian( const InputPointType & point, 
		 unsigned int time, bool accumulate, std::vector<long unsigned int>&) const;
  
protected:
  DiffeomorphicBSplineTransform();
  ~DiffeomorphicBSplineTransform();
  
  /** Initialize internal transformations. */
  virtual void InitializeInternalTransforms();
  
  /** Print contents of an DiffemorphicBSplineTransform. */
  void PrintSelf(std::ostream &os, Indent indent) const;
  
  std::vector<InternalTransformPointerType> m_InternalTransforms;
  mutable std::vector<InternalJacobianType> m_InternalJacobians;
  mutable std::vector<WeightsType> m_InternalWeights;
  mutable std::vector<ParameterIndexArrayType> m_InternalIndexes;
  mutable std::vector<OutputPointType> m_InternalPoints;

  // array storing parameters passed to intenal transforms
  std::vector<typename InternalTransformType::ParametersType> m_InternalParameters;
  
  SpacingType m_GridSpacing;
  OriginType m_GridOrigin;
  RegionType m_GridRegion;
  
  unsigned int m_NumberOfTimeSteps;
  bool m_InternalTransformsInitialized;
  
  /** Keep a pointer to the input parameters. */
  const ParametersType*  m_InputParametersPointer;

  /** Internal parameters buffer. 
  * Stores a copy of the input parameters if the parameters
  *  are passed by value. */
  ParametersType          m_InternalParametersBuffer;
  
  private:
  DiffeomorphicBSplineTransform(const Self&); //purposely not implemented
  void operator=(const Self&); //purposely not implemented
  
  }; //class DiffeomorphicBSplineTransform
  
}  // namespace itk


#if ITK_TEMPLATE_TXX
# include "itkDiffeomorphicBSplineTransform.txx"
#endif

#endif /* __itkDiffeomorphicBSplineTransform_h */
