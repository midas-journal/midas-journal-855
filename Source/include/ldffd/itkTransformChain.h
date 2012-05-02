/*=========================================================================

  Program:   Insight Segmentation & Registration Toolkit
  Module:    $RCSfile: itkTransformChain.h,v $
  Language:  C++
  Date:      $Date: 2008-06-20 15:25:58 $
  Version:   $Revision: 1.2 $

  Copyright (c) Insight Software Consortium. All rights reserved.
  See ITKCopyright.txt or http://www.itk.org/HTML/Copyright.htm for details.

     This software is distributed WITHOUT ANY WARRANTY; without even
     the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
     PURPOSE.  See the above copyright notices for more information.

=========================================================================*/

#ifndef __itkTransformChain_h
#define __itkTransformChain_h

#include <iostream>
#include "itkTransform.h"
#include "itkExceptionObject.h"
#include "itkMatrix.h"

namespace itk
{

template <
    class TScalarType=double,          // Data type for scalars (float or double)
    unsigned int NDimensions=3>        // Number of dimensions
class ITK_EXPORT TransformChain :
          public Transform< TScalarType, NDimensions, NDimensions >
{
public:
  /** Standard class typedefs. */
  typedef TransformChain Self;
  typedef Transform< TScalarType, NDimensions, NDimensions > Superclass;
  typedef SmartPointer<Self>        Pointer;
  typedef SmartPointer<const Self>  ConstPointer;

  /** New macro for creation of through the object factory.*/
  itkNewMacro( Self );

  /** Run-time type information (and related methods). */
  itkTypeMacro( TransformChain, Transform );

  /** Dimension of the domain space. */
  itkStaticConstMacro(SpaceDimension, unsigned int, NDimensions);
  //itkStaticConstMacro(ParametersDimension, unsigned int, NDimensions);

  /** Standard scalar type for this class. */
  typedef typename Superclass::ScalarType ScalarType;

  /** Standard parameters container. */
  typedef typename Superclass::ParametersType ParametersType;

  /** Standard coordinate point type for this class. */
  typedef Point<TScalarType, itkGetStaticConstMacro(SpaceDimension)> InputPointType;
  typedef Point<TScalarType, itkGetStaticConstMacro(SpaceDimension)> OutputPointType;

   OutputPointType     TransformPoint(const InputPointType  &point ) const;
  /** Set the parameters to the IdentityTransform */
  void SetIdentity(void);

  void PushBackTransformation(Superclass* transform)
  {
    m_Transforms.push_back(transform);
  }

  unsigned int GetNumberOfTransforms(void) const
  {
    return m_Transforms.size();
  }

  Superclass* GetTransform(unsigned int pos)
  {
    return m_Transforms.at(pos);
  }

  void Erase();

protected:
  TransformChain();
  ~TransformChain();
  /** Print contents of an TransformChain. */
  void PrintSelf(std::ostream &os, Indent indent) const;

private:
  TransformChain(const Self&); //purposely not implemented
  void operator=(const Self&); //purposely not implemented

  std::vector<typename Superclass::Pointer> m_Transforms;

}; //class TransformChain

}  // namespace itk



#if ITK_TEMPLATE_TXX
# include "itkTransformChain.txx"
#endif

#endif /* __itkTransformChain_h */
