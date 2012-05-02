/*=========================================================================

  Program:   Insight Segmentation & Registration Toolkit
  Module:    $RCSfile: itkTransformChain.txx,v $
  Language:  C++
  Date:      $Date: 2008-06-20 15:25:58 $
  Version:   $Revision: 1.2 $

  Copyright (c) Insight Software Consortium. All rights reserved.
  See ITKCopyright.txt or http://www.itk.org/HTML/Copyright.htm for details.

     This software is distributed WITHOUT ANY WARRANTY; without even
     the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
     PURPOSE.  See the above copyright notices for more information.

=========================================================================*/
#ifndef _itkTransformChain_txx
#define _itkTransformChain_txx

#include "itkTransformChain.h"


namespace itk
{

// Constructor with default arguments
template<class TScalarType, unsigned int NDimensions>
TransformChain<TScalarType, NDimensions>::
TransformChain():Superclass(SpaceDimension,0)
{
m_Transforms.resize(0);
}


// Destructor
template<class TScalarType, unsigned int NDimensions>
TransformChain<TScalarType, NDimensions>::
~TransformChain()
{
  return;
}

// Print self
template<class TScalarType, unsigned int NDimensions>
void
TransformChain<TScalarType, NDimensions>::
PrintSelf(std::ostream &os, Indent indent) const
{
  Superclass::PrintSelf(os,indent);
}

// Transform a point
template<class TScalarType, unsigned int NDimensions>
typename TransformChain<TScalarType, NDimensions>::OutputPointType
TransformChain<TScalarType, NDimensions>::
TransformPoint(const InputPointType &point) const
{
  OutputPointType currentPoint = point;
  for (unsigned int index=0; index<m_Transforms.size(); index++)
    {
    currentPoint = m_Transforms[index]->TransformPoint(currentPoint);
    }
  return currentPoint;
}

// Set the parameters for an Identity transform of this class
template<class TScalarType, unsigned int NDimensions>
void
TransformChain<TScalarType, NDimensions>::
SetIdentity()
{
  for (unsigned int index=0; index<m_Transforms.size(); index++)
    {
    m_Transforms[index]->SetIdentity();
    }
}

// Set the parameters for an Identity transform of this class
template<class TScalarType, unsigned int NDimensions>
void
TransformChain<TScalarType, NDimensions>::
Erase()
{
    m_Transforms.resize(0);
}


} // namespace

#endif
