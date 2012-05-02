/*=========================================================================

  Program:   Insight Segmentation & Registration Toolkit
  Module:    $RCSfile: itkBSplineSecondDerivativeKernelFunction.h,v $
  Language:  C++
  Date:      $Date: 2008-06-25 11:00:19 $
  Version:   $Revision: 1.7 $

  Copyright (c) Insight Software Consortium. All rights reserved.
  See ITKCopyright.txt or http://www.itk.org/HTML/Copyright.htm for details.

     This software is distributed WITHOUT ANY WARRANTY; without even 
     the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR 
     PURPOSE.  See the above copyright notices for more information.

=========================================================================*/
#ifndef __itkBSplineSecondDerivativeKernelFunction_h
#define __itkBSplineSecondDerivativeKernelFunction_h

#include "itkKernelFunction.h"
#include "itkBSplineKernelFunction.h"

namespace itk
{

/** \class BSplineSecondDerivativeKernelFunction
 * \brief Derivative of a BSpline kernel used for density estimation and 
 *  nonparameteric regression.
 *
 * This class encapsulates the derivative of a BSpline kernel for
 * density estimation or nonparameteric regression.
 * See documentation for KernelFunction for more details.
 *
 * This class is templated over the spline order.
 * \warning Evaluate is only implemented for spline order 1 to 4
 *
 * \sa KernelFunction
 *
 * \ingroup Functions
 */
template <unsigned int VSplineOrder = 3>
class ITK_EXPORT BSplineSecondDerivativeKernelFunction : public KernelFunction
{
public:
  /** Standard class typedefs. */
  typedef BSplineSecondDerivativeKernelFunction Self;
  typedef KernelFunction                  Superclass;
  typedef SmartPointer<Self>              Pointer;

  /** Method for creation through the object factory. */
  itkNewMacro(Self); 

  /** Run-time type information (and related methods). */
  itkTypeMacro(BSplineSecondDerivativeKernelFunction, KernelFunction); 

  /** Enum of for spline order. */
  itkStaticConstMacro(SplineOrder, unsigned int, VSplineOrder);

  /** Evaluate the function. */
  inline double Evaluate( const double & u ) const
    {
    return ( m_KernelFunction->Evaluate( u + 1. ) + 
      m_KernelFunction->Evaluate( u - 1. ) );
    }

protected:

  typedef BSplineKernelFunction<itkGetStaticConstMacro(SplineOrder) - 2>
      KernelType;

  BSplineSecondDerivativeKernelFunction()
    {
    m_KernelFunction = KernelType::New();
    }

  ~BSplineSecondDerivativeKernelFunction(){};
  void PrintSelf(std::ostream& os, Indent indent) const
    { 
    Superclass::PrintSelf( os, indent ); 
    os << indent  << "Spline Order: " << SplineOrder << std::endl;
    }  

private:
  BSplineSecondDerivativeKernelFunction(const Self&); //purposely not implemented
  void operator=(const Self&); //purposely not implemented

  typename KernelType::Pointer  m_KernelFunction;

};

} // end namespace itk

#endif
