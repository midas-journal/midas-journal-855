/*=========================================================================

  Program:   Insight Segmentation & Registration Toolkit
  Module:    $RCSfile: itkTimeDiffeomorphicTransform.h,v $
  Language:  C++
  Date:      $Date: 2008-06-29 12:58:58 $
  Version:   $Revision: 1.64 $

  Copyright (c) Insight Software Consortium. All rights reserved.
  See ITKCopyright.txt or http://www.itk.org/HTML/Copyright.htm for details.

     This software is distributed WITHOUT ANY WARRANTY; without even 
     the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR 
     PURPOSE.  See the above copyright notices for more information.

=========================================================================*/
#ifndef __regTimeDiffeomorphicTransform_h
#define __regTimeDiffeomorphicTransform_h

#include "itkTransformBase.h"
#include "itkPoint.h"
#include "itkVector.h"
#include "itkCovariantVector.h"
#include "vnl/vnl_vector_fixed.h"
#include "itkArray.h"
#include "itkArray2D.h"
#include "vnl/vnl_sparse_matrix.h"
#include "itkTimeProbe.h"

#include "itkObjectFactory.h"


namespace itk
{

// Useful types

class ltCol
{
public:
  bool operator()(unsigned int col1, unsigned int col2) const
  {
    return col1<col2;
  }
  ltCol(){}
  ~ltCol(){}
};

/** \class TimeDiffeomorphicTransform
 * \brief TimeDiffeomorphicTransform points and vector from an input space to an output space.
 *
 * This abstract class define the generic interface for a geometrical 
 * transformation from one space to another. The class provides methods
 * for mapping points, vectors and covariant vectors from the input space 
 * to the output space. 
 *
 * Given that transformation are not necesarily invertible, this basic
 * class does not provide the methods for back transfromation. Back transform
 * methods are implemented in derived classes where appropriate.
 * 
 * \par Registration Framework Support
 * Typically a TimeDiffeomorphicTransform class has several methods for setting its 
 * parameters. For use in the registration framework, the parameters must
 * also be represented by an array of doubles to allow communication
 * with generic optimizers. The Array of transformation parameters is set using
 * the SetParameters() method.
 *
 * Another requirement of the registration framework is the computation
 * of the transform Jacobian. In general, a ImageToImageMetric requires
 * the knowledge of the Jacobian in order to compute the metric derivatives.
 * The Jacobian is a matrix whose element are the partial derivatives
 * of the output point with respect to the array of parameters that defines
 * the transform.
 *
 * \ingroup TimeDiffeomorphicTransforms
 *
 */
template <class TScalarType,
          unsigned int NInputDimensions=4, 
          unsigned int NOutputDimensions=4>
class ITK_EXPORT  TimeDiffeomorphicTransform  : public TransformBase
	/** Transform<TScalarType, NInputDimensions, NOutputDimensions>*/
{
public:
  /** Standard class typedefs. */
  typedef TimeDiffeomorphicTransform  Self;
  typedef TransformBase Superclass;
  typedef SmartPointer< Self >   Pointer;
  typedef SmartPointer< const Self >  ConstPointer;
  
  /** New method for creating an object using a factory. */
  itkNewMacro(Self);

  /** Run-time type information (and related methods). */
  itkTypeMacro( TimeDiffeomorphicTransform, TransformBase );

  /** Dimension of the domain space. */
  itkStaticConstMacro(InputSpaceDimension, unsigned int, NInputDimensions);
  itkStaticConstMacro(OutputSpaceDimension, unsigned int, NOutputDimensions);

  /** Get the size of the input space */
  unsigned int GetInputSpaceDimension(void) const {return NInputDimensions;}

  /** Get the size of the output space */
  unsigned int GetOutputSpaceDimension(void) const {return NOutputDimensions;}

  /** Type of the scalar representing coordinate and vector elements. */
  typedef  TScalarType     ScalarType;

  /** Type of the input parameters. */
  typedef  typename Superclass::ParametersType         ParametersType;

  /** Type of the Jacobian matrix. */
  typedef  Array2D< double >                           JacobianType;

  /** Type of the Sparse Jacobian matrix. */
  typedef  vnl_vector<double> SparseJacobianColumnType;
  typedef  std::map< unsigned int, SparseJacobianColumnType, ltCol > SparseJacobianType;

  /** Standard vector type for this class. */
  typedef Vector<TScalarType, NInputDimensions>  InputVectorType;
  typedef Vector<TScalarType, NOutputDimensions> OutputVectorType;
  
  /** Standard covariant vector type for this class */
  typedef CovariantVector<TScalarType, NInputDimensions>  InputCovariantVectorType;
  typedef CovariantVector<TScalarType, NOutputDimensions> OutputCovariantVectorType;
  
  /** Standard vnl_vector type for this class. */
  typedef vnl_vector_fixed<TScalarType, NInputDimensions>  InputVnlVectorType;
  typedef vnl_vector_fixed<TScalarType, NOutputDimensions> OutputVnlVectorType;
  
  /** Standard coordinate point type for this class */
  typedef Point<TScalarType, NInputDimensions> InputPointType;
  typedef Point<TScalarType, NOutputDimensions> OutputPointType;
  
  /**  Method to transform a point. */
  virtual OutputPointType TransformPoint(const InputPointType  &, const TScalarType &) const
  { return OutputPointType(); } 

  /**  Method to transform a vector. */
  virtual OutputVectorType    TransformVector(const InputVectorType &) const
  { return OutputVectorType(); }

  /**  Method to transform a vnl_vector. */
  virtual OutputVnlVectorType TransformVector(const InputVnlVectorType &) const
  { return OutputVnlVectorType(); }

  /**  Method to transform a CovariantVector. */
  virtual OutputCovariantVectorType TransformCovariantVector(
    const InputCovariantVectorType &) const
  { return OutputCovariantVectorType(); } 

  /** Set the transformation parameters and update internal transformation.
   * SetParameters gives the transform the option to set it's
   * parameters by keeping a reference to the parameters, or by
   * copying.  To force the transform to copy it's parameters call
   * SetParametersByValue.
   * \sa SetParametersByValue
   */
  virtual void SetParameters( const ParametersType & ) 
  { itkExceptionMacro( << "Subclasses should override this method" ) }

  /** Set the transformation parameters and update internal transformation. 
   * This method forces the transform to copy the parameters.  The
   * default implementation is to call SetParameters.  This call must
   * be overridden if the transform normally implements SetParameters
   * by keeping a reference to the parameters.
   * \sa SetParameters
   */
  virtual void SetParametersByValue ( const ParametersType & p ) 
  { this->SetParameters ( p ); }

  /** Get the Transformation Parameters. */
  virtual const ParametersType& GetParameters(void) const
  { 
    itkExceptionMacro( << "Subclasses should override this method" );
    // Next line is needed to avoid errors due to: 
    // "function must return a value".
    return this->m_Parameters; 
  }

  /** Set the fixed parameters and update internal transformation. */
  virtual void SetFixedParameters( const ParametersType & ) 
  { itkExceptionMacro( << "Subclasses should override this method" ) }

  /** Get the Fixed Parameters. */
  virtual const ParametersType& GetFixedParameters(void) const
  {
    itkExceptionMacro( << "Subclasses should override this method" );
    // Next line is needed to avoid errors due to: 
    // "function must return a value".
    return this->m_FixedParameters;
  }

  /** Compute the Jacobian of the transformation
   *
   * This method computes the Jacobian matrix of the transformation
   * at a given input point. The rank of the Jacobian will also indicate 
   * if the transform is invertible at this point.
   *
   * The Jacobian is be expressed as a matrix of partial derivatives of the
   * output point components with respect to the parameters that defined
   * the transform:
   *
   * \f[
   *
      J=\left[ \begin{array}{cccc}
      \frac{\partial x_{1}}{\partial p_{1}} & 
      \frac{\partial x_{1}}{\partial p_{2}} & 
      \cdots  & \frac{\partial x_{1}}{\partial p_{m}}\\
      \frac{\partial x_{2}}{\partial p_{1}} & 
      \frac{\partial x_{2}}{\partial p_{2}} & 
      \cdots  & \frac{\partial x_{2}}{\partial p_{m}}\\
      \vdots  & \vdots  & \ddots  & \vdots \\
      \frac{\partial x_{n}}{\partial p_{1}} & 
      \frac{\partial x_{n}}{\partial p_{2}} & 
      \cdots  & \frac{\partial x_{n}}{\partial p_{m}}
      \end{array}\right] 
   *
   * \f]
   * **/
  virtual const JacobianType & GetJacobian(const InputPointType  &, 
					   const TScalarType &) const
  {
    itkExceptionMacro( << "Subclass should override this method" );
    // Next line is needed to avoid errors due to: 
    // "function must return a value".
    return this->m_Jacobian;
  } 

  virtual void GetSparseJacobian(const InputPointType  &, 
						       const TScalarType &, SparseJacobianType & jaco) const
  {
    itkExceptionMacro( << "Subclass should override this method" );
  } 
  
  // This method returns true if the transform can handle sparse jacobian
  // false otherwise
  virtual bool SupportIncrementalJacobian()
  {
    itkExceptionMacro( << "Subclass should override this method" );
    // Next line is needed to avoid errors due to: 
    // "function must return a value".
    return false;
  }

  // This method updates the current jacobian at the time of the input
  // point to the final time, without resetting the jacobian
  virtual void GetIncrementalSparseJacobian(const InputPointType  & point,
              OutputPointType & outPoint,                              
					    const TScalarType & endTimePoint, 
					    SparseJacobianType & jaco, 
					    bool  & valid) const 
  {
    itkExceptionMacro( << "Subclass should override this method" );
    return; 
  } 
  

  /** Return the number of parameters that completely define the Transfom  */
  virtual unsigned int GetNumberOfParameters(void) const 
  { return this->m_Parameters.Size(); }

  /** Returns a boolean indicating whether it is possible or not to compute the
   * inverse of this current Transform. If it is possible, then the inverse of
   * the transform is returned in the inverseTransform variable passed by the
   * user.  The inverse is recomputed if this current transform has been modified.
   * This method is intended to be overriden by derived classes.
   * 
   */
  bool GetInverse(Self * inverseTransform) const {return false;}

  /** Generate a platform independant name */
  virtual std::string GetTransformTypeAsString() const;

  /** Indicates if this transform is linear. A transform is defined to be
   * linear if the transform of a linear combination of points is equal to the
   * linear combination (with the same coefficients) of the individual
   * transforms of each point. The transform T will be linear if given two
   * points P and Q, and scalar coefficients a and b, then
   *
   *           T( a*P + b*Q ) = a * T(P) + b * T(Q)
   * 
   * By default, we assume this to NOT be the case for most transforms.
   * However, transforms for which this is true will overload and reimplement
   * this method accordingly.
   * 
   **/
  virtual bool IsLinear() const { return false; }

protected:
  TimeDiffeomorphicTransform(); 
  TimeDiffeomorphicTransform(unsigned int Dimension, unsigned int NumberOfParameters);
  virtual ~TimeDiffeomorphicTransform() {}

  mutable ParametersType     m_Parameters;
  mutable ParametersType     m_FixedParameters;
  mutable JacobianType       m_Jacobian;
  
private:
  TimeDiffeomorphicTransform(const Self&); //purposely not implemented
  void operator=(const Self&); //purposely not implemented

public:
  mutable TimeProbe probeAll;
  mutable TimeProbe probeVel;
  mutable TimeProbe probePhysJ;
  mutable TimeProbe probeMult1;
  mutable TimeProbe probeMult2;
  
};
  
} // end namespace itk

#if ITK_TEMPLATE_EXPLICIT
# include "Templates/regTimeDiffeomorphicTransform+-.h"
#endif

#if ITK_TEMPLATE_TXX
# include "regTimeDiffeomorphicTransform.txx"
#endif

#endif
