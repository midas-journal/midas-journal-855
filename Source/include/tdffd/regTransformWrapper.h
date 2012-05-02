/*=========================================================================

  Program   Insight Segmentation & Registration Toolkit
  Module    $RCSfile itkTransformWrapper.h,v $
  Language  C++
  Date      $Date 2008-06-20 152558 $
  Version   $Revision 1.2 $

  Copyright (c) Insight Software Consortium. All rights reserved.
  See ITKCopyright.txt or httpwww.itk.orgHTMLCopyright.htm for details.

     This software is distributed WITHOUT ANY WARRANTY; without even
     the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
     PURPOSE.  See the above copyright notices for more information.

=========================================================================*/

#ifndef __regTransformWrapper_h
#define __regTransformWrapper_h

#include <iostream>
#include <itkTransform.h>
#include <itkExceptionObject.h>

namespace itk
{

template <class TScalarType=double,          // Data type for scalars (float or double)
    unsigned int NDimensions=3>        // Number of dimensions
class ITK_EXPORT TransformWrapper
	:  public Transform <TScalarType, NDimensions, NDimensions> 
{
public:
  // Standard class typedefs. 
  typedef TransformWrapper Self;
  typedef Transform <TScalarType, NDimensions, NDimensions>  Superclass;
  typedef SmartPointer<Self>        Pointer;
  typedef SmartPointer<const Self > ConstPointer;

  // New macro for creation of through the object factory.
  itkNewMacro( Self );

  // Run-time type information (and related methods). 
  itkTypeMacro( TransformWrapper, Transform );
 
  // Dimension of the domain space. 
  itkStaticConstMacro(SpaceDimension, unsigned int, NDimensions);
  itkStaticConstMacro(SpaceTimeDimension, unsigned int, NDimensions+1);
  itkStaticConstMacro(ParametersDimension, unsigned int, NDimensions);

  // Time diff transform type
  typedef TimeDiffeomorphicTransform<double, SpaceTimeDimension, SpaceTimeDimension> InternalTranformType;
  typedef typename InternalTranformType::Pointer InternalTranformPointerType;

  // Standard scalar type for this class. 
  typedef typename Superclass::ScalarType ScalarType;

  // Standard parameters container. 
  typedef typename Superclass::ParametersType ParametersType;

  // Standard coordinate point type for this class. 
  typedef Point<TScalarType, itkGetStaticConstMacro(SpaceDimension)> PointType;
  typedef Point<TScalarType, itkGetStaticConstMacro(SpaceTimeDimension)> PointTypeExt;
 
  PointType     TransformPoint(const PointType  &point ) const
  {
	 PointTypeExt point2, outPoint2;
	 PointType outPoint;
	 for (unsigned int d=0; d<SpaceDimension; d++)
	 {
	   point2[d] = point[d];
	 }
	 point2[SpaceDimension] = m_InitTime;
	 outPoint2 =  m_Transform->TransformPoint(point2, m_FinalTime);
	 for (unsigned int d=0; d<SpaceDimension; d++)
	 {
	   outPoint[d] = outPoint2[d];
	 }
	 return outPoint;
  }

  void SetFinalTime (double time) 
  {
	  m_FinalTime = time;
  }
  
  void SetInitTime (double time) 
  {
	  m_InitTime = time;
  }

  void SetTransform (InternalTranformType* tr)
  {
	  m_Transform = tr;
  }

protected:
	TransformWrapper():Superclass(NDimensions, 0)
	{
		m_Transform = 0;
		m_FinalTime = 0.;
		m_InitTime = 0.;
	}
	~TransformWrapper(){};
  void PrintSelf(std::ostream &os, Indent indent) const
	{
		std::cout<<"Coucou gamin"<<std::endl;
	}

  InternalTranformPointerType m_Transform;
  double m_FinalTime;
  double m_InitTime;

private:
  TransformWrapper(const Self&); //purposely not implemented
  void operator=(const Self&); //purposely not implemented

}; //class TransformWrapper

} // namespace itk

#endif // __itkTransformWrapper_h 
