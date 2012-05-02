/*=========================================================================

  Program:   Insight Segmentation & Registration Toolkit
  Module:    $RCSfile: itkBSplineInterpolationWeightFunctionDerivative.h,v $
  Language:  C++
  Date:      $Date: 2008-06-20 15:25:58 $
  Version:   $Revision: 1.2 $

  Copyright (c) Insight Software Consortium. All rights reserved.
  See ITKCopyright.txt or http://www.itk.org/HTML/Copyright.htm for details.

     This software is distributed WITHOUT ANY WARRANTY; without even 
     the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR 
     PURPOSE.  See the above copyright notices for more information.

=========================================================================*/
#ifndef __itkBSplineInterpolationWeightFunctionDerivative_h
#define __itkBSplineInterpolationWeightFunctionDerivative_h

#include "itkFunctionBase.h"
#include "itkContinuousIndex.h"
#include "itkBSplineKernelFunction.h"
#include "itkArray.h"
#include "itkIndex.h"
#include "itkArray2D.h"
#include "itkBSplineDerivativeKernelFunction.h"
#include "itkBSplineInterpolationWeightFunction.h"
#include "vnl/vnl_matrix.h"

namespace itk
{

/** \class BSplineInterpolationWeightFunctionDerivative
 * \brief 
 *
 * \ingroup Functions ImageInterpolators
 */
template <
class TCoordRep = float, 
unsigned int VSpaceDimension = 2,
unsigned int VSplineOrder = 3 
>
class ITK_EXPORT BSplineInterpolationWeightFunctionDerivative : 
  public FunctionBase< ContinuousIndex<TCoordRep,VSpaceDimension>, 
                       Array<double> > 
{
public:
  /** Standard class typedefs. */
  typedef BSplineInterpolationWeightFunctionDerivative Self;
  typedef FunctionBase< ContinuousIndex<TCoordRep,VSpaceDimension>,
            Array<double> >                  Superclass;
  typedef SmartPointer<Self>                 Pointer;
  typedef SmartPointer<const Self>           ConstPointer;

  /** New macro for creation of through the object factory.*/
  itkNewMacro( Self );
  
  /** Run-time type information (and related methods). */
  itkTypeMacro(BSplineInterpolationWeightFunctionDerivative, FunctionBase);

  /** Space dimension. */
  itkStaticConstMacro(SpaceDimension, unsigned int, VSpaceDimension);

  /** Spline order. */
  itkStaticConstMacro(SplineOrder, unsigned int, VSplineOrder);

  /** OutputType typedef support. */
  typedef Array<double> WeightsType;

  /** Index and size typedef support. */
  typedef Index<VSpaceDimension> IndexType;
  typedef Size<VSpaceDimension>  SizeType;

  /** ContinuousIndex typedef support. */
  typedef ContinuousIndex<TCoordRep,VSpaceDimension> ContinuousIndexType;

  /** Evaluate the weights at specified ContinousIndex position.
   * Subclasses must provide this method. */
  virtual WeightsType Evaluate( const ContinuousIndexType & index ) const;
  
  /** Function to compute derivative of the weights
   * in one direction for jacobian computation. 
   */
  virtual void EvaluateDerivative( const ContinuousIndexType & index,
				   WeightsType & weights, IndexType & startIndex) const;

  /** Get support region size. */
  itkGetMacro( SupportSize, SizeType );
  
  /** Get number of weights. */
  itkGetMacro( NumberOfWeights, unsigned long );
  
protected:
  BSplineInterpolationWeightFunctionDerivative();
  ~BSplineInterpolationWeightFunctionDerivative() {}
  void PrintSelf(std::ostream& os, Indent indent) const;
  
private:
  BSplineInterpolationWeightFunctionDerivative(const Self&); //purposely not implemented
  void operator=(const Self&); //purposely not implemented

  /** Number of weights. */
  unsigned long  m_NumberOfWeights; 

  /** Size of support region. */
  SizeType m_SupportSize;

  /** Lookup table type. */
  typedef Array2D<unsigned long> TableType;

  /** Table mapping linear offset to indices. */
  TableType m_OffsetToIndexTable;

  /** Interpolation kernel type. */
  typedef BSplineKernelFunction<itkGetStaticConstMacro(SplineOrder)> KernelType;
      
  /** Interpolation kernel. */
  typename KernelType::Pointer m_Kernel;

  /** Kernel derivatives. */
  typedef BSplineDerivativeKernelFunction<itkGetStaticConstMacro(SplineOrder)> KernelDerivativeType;
  typename KernelDerivativeType::Pointer m_KernelDerivative;
  
};

} // end namespace itk

#if ITK_TEMPLATE_TXX
# include "itkBSplineInterpolationWeightFunctionDerivative.txx"
#endif


#endif
