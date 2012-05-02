/*=========================================================================

  Program:   Insight Segmentation & Registration Toolkit
  Module:    $RCSfile: itkBSplineInterpolationWeightFunctionDerivative.txx,v $
  Language:  C++
  Date:      $Date: 2008-06-20 15:25:58 $
  Version:   $Revision: 1.2 $

  Copyright (c) Insight Software Consortium. All rights reserved.
  See ITKCopyright.txt or http://www.itk.org/HTML/Copyright.htm for details.

     This software is distributed WITHOUT ANY WARRANTY; without even 
     the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR 
     PURPOSE.  See the above copyright notices for more information.

=========================================================================*/
#ifndef __itkBSplineInterpolationWeightFunctionDerivative_txx
#define __itkBSplineInterpolationWeightFunctionDerivative_txx

#include "itkBSplineInterpolationWeightFunctionDerivative.h"
#include "itkImage.h"
#include "itkMatrix.h"
#include "itkImageRegionConstIteratorWithIndex.h"



namespace itk
{

/** Constructor */
template <class TCoordRep, unsigned int VSpaceDimension, 
                                                  unsigned int VSplineOrder>
BSplineInterpolationWeightFunctionDerivative<TCoordRep, VSpaceDimension, VSplineOrder>
::BSplineInterpolationWeightFunctionDerivative()
{

  // Initialize the number of weights;
  m_NumberOfWeights = 
    static_cast<unsigned long>( vcl_pow(static_cast<double>(SplineOrder + 1),
                                     static_cast<double>(SpaceDimension) ) );

  // Initialize support region is a hypercube of length SplineOrder + 1
  m_SupportSize.Fill( SplineOrder + 1 );

  // Initialize offset to index lookup table
  m_OffsetToIndexTable.set_size( m_NumberOfWeights, SpaceDimension );

  typedef Image<char,SpaceDimension> CharImageType;
  typename CharImageType::Pointer tempImage = CharImageType::New();
  tempImage->SetRegions( m_SupportSize );
  tempImage->Allocate();
  tempImage->FillBuffer( 0 );



  typedef ImageRegionConstIteratorWithIndex<CharImageType> IteratorType;
  IteratorType iterator( tempImage, tempImage->GetBufferedRegion() );
  unsigned long counter = 0;

  while ( !iterator.IsAtEnd() )
    {
    for(unsigned  int j = 0; j < SpaceDimension; j++ )
      {
      m_OffsetToIndexTable[counter][j] = iterator.GetIndex()[j];
      }
    ++counter;
    ++iterator;
    }  


  // Initialize the interpolation kernel
  m_Kernel = KernelType::New();
  m_KernelDerivative = KernelDerivativeType::New();

}


/**
 * Standard "PrintSelf" method
 */
template <class TCoordRep, unsigned int VSpaceDimension, 
                                                  unsigned int VSplineOrder>
void
BSplineInterpolationWeightFunctionDerivative<TCoordRep, VSpaceDimension, VSplineOrder>
::PrintSelf(
  std::ostream& os, 
  Indent indent) const
{
  Superclass::PrintSelf( os, indent );
  
  os << indent << "NumberOfWeights: " << m_NumberOfWeights << std::endl;
  os << indent << "SupportSize: " << m_SupportSize << std::endl;
}


/** Compute weights for interpolation at continous index position */
template <class TCoordRep, unsigned int VSpaceDimension, 
                                                  unsigned int VSplineOrder>
typename BSplineInterpolationWeightFunctionDerivative<TCoordRep, VSpaceDimension, 
                                                               VSplineOrder>
::WeightsType
BSplineInterpolationWeightFunctionDerivative<TCoordRep, VSpaceDimension, VSplineOrder>
::Evaluate(
  const ContinuousIndexType& index ) const
{

  WeightsType weights( SpaceDimension * m_NumberOfWeights );
  IndexType startIndex;

  this->EvaluateDerivative( index, weights, startIndex);

  return weights;
}


/** Compute weights DERIVATIVES for interpolation at continous index position */
template <class TCoordRep, unsigned int VSpaceDimension, 
                                                 unsigned int VSplineOrder>
void BSplineInterpolationWeightFunctionDerivative<TCoordRep, VSpaceDimension, 
                                                              VSplineOrder>
::EvaluateDerivative(
  const ContinuousIndexType & index,
  WeightsType & weights, 
  IndexType & startIndex) const
{

  unsigned int j, k;
  
  // Find the starting index of the support region
  for ( j = 0; j < SpaceDimension; j++ )
    {
    startIndex[j] = static_cast<typename IndexType::IndexValueType>(
	  //NOTE: bSplineFloor was removed in itk-3.16. They use instead math::floor. We did the same. MBianculli
      Math::Floor( index[j] - static_cast<double>( SplineOrder + 1 ) / 2.  ) + 1 );
    }
  
  // Compute the weights
  Matrix<double,SpaceDimension,SplineOrder + 1> weights1D;
  Matrix<double,SpaceDimension,SplineOrder + 1> weights1Dderivative;
  for ( j = 0; j < SpaceDimension; j++ )
    {
    double x = index[j] - static_cast<double>( startIndex[j] );
    
    for( k = 0; k <= SplineOrder; k++ )
      {
      double kernelValue = m_Kernel->Evaluate( x );
      weights1D[j][k] = kernelValue;
      weights1Dderivative[j][k] = m_KernelDerivative->Evaluate(x);
      
      x -= 1.0;
      }
    }
  
  // Fill in the output weights vector
  // This vector contains de BSplines derivatives : m_NumberOfWeights
  // for dimension 0, m_NumberOfWeights for dimension 1, ...
  weights.Fill(1.);
  for (unsigned int direction=0; direction < SpaceDimension; direction++)
    {
    for ( k = 0; k < m_NumberOfWeights; k++ )
      {
      
      for ( j = 0; j < SpaceDimension; j++ )
	{
	if ( j == direction)
	  {
	  weights[k + m_NumberOfWeights * direction] 
	    *= weights1Dderivative[j][ m_OffsetToIndexTable[k][j] ];
	  }
	else
	  {
	  weights[k + m_NumberOfWeights * direction] 
	    *= weights1D[j][ m_OffsetToIndexTable[k][j] ];
	  }
	}
      }
    }
}

} // end namespace itk

#endif
