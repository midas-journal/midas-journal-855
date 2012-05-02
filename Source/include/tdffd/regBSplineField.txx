#ifndef __regBSplineField_txx
#define __regBSplineField_txx

#include "regBSplineField.h"
#include "itkImageRegionIterator.h"

namespace itk
{

/**
 * Constructor
 */
template<class TScalarType, unsigned int NDimensions, unsigned int VSplineOrder>
BSplineField<TScalarType, NDimensions, VSplineOrder>
::BSplineField()
{

  // Instantiate a weights function
  m_WeightsFunction = WeightsFunctionType::New();
  m_SupportSize = m_WeightsFunction->GetSupportSize();
  
  // Default grid size is zero
  typename RegionType::SizeType size;
  typename RegionType::IndexType index;
  size.Fill( 0 );
  index.Fill( 0 );
  m_GridRegion.SetSize( size );
  m_GridRegion.SetIndex( index );
  m_GridOrigin.Fill( 0.0 );  // default origin is all zeros
  m_GridSpacing.Fill( 1.0 ); // default spacing is all ones
  m_GridDirection.SetIdentity(); // default spacing is all ones
  
    // Initialize coeffient images
  for ( unsigned int j = 0; j < SpaceDimension; j++ )
    {
    m_WrappedImage[j] = ImageType::New();
    m_WrappedImage[j]->SetRegions( m_GridRegion );
    m_WrappedImage[j]->SetOrigin( m_GridOrigin.GetDataPointer() );
    m_WrappedImage[j]->SetSpacing( m_GridSpacing.GetDataPointer() );
    m_WrappedImage[j]->SetDirection( m_GridDirection );
    m_CoefficientImage[j] = NULL;
    }

  // Setup variables for computing interpolation
  m_Offset = SplineOrder / 2;
  if ( SplineOrder % 2 ) 
    {
    m_SplineOrderOdd = true;
    }
  else
    {
    m_SplineOrderOdd = false;
    }
  m_ValidRegion = m_GridRegion;
  
  m_weightsFunctionDerivative = WeightFunctionDerivativeType::New();
  m_weightsFunctionSecondDerivative = WeightFunctionSecondDerivativeType::New();
  
  m_InputParametersPointer = NULL;
  
  DirectionType scale;
  for( unsigned int i=0; i<SpacePlusTimeDimension; i++)
    {
    scale[i][i] = this->GetGridSpacing()[i];
    }
  
  m_IndexToPoint = m_GridDirection * scale;
  m_PointToIndex = m_IndexToPoint.GetInverse();
  
  
}

/**
 * Set parameters
 */
template<class TScalarType, unsigned int NDimensions, unsigned int VSplineOrder>
void
BSplineField<TScalarType, NDimensions, VSplineOrder>
::SetParameters( const ParametersType & params ) 
{
  m_InputParametersPointer = & params;
}


/**
 * Get parameters
 */
template<class TScalarType, unsigned int NDimensions, unsigned int VSplineOrder>
const typename BSplineField<TScalarType, NDimensions, VSplineOrder>::ParametersType&  
BSplineField<TScalarType, NDimensions, VSplineOrder>
::GetParameters(void) const
{
	
	if (m_InputParametersPointer == NULL)
	{
		itkExceptionMacro( << "Parameters were incorrectly set." << std::endl );
	}
	
	return (*m_InputParametersPointer);

}

/**
 * Get Number of parameters
 */
template<class TScalarType, unsigned int NDimensions, unsigned int VSplineOrder>
unsigned int 
BSplineField<TScalarType, NDimensions, VSplineOrder>
::GetNumberOfParameters(void) const
{
  // The number of parameters equal SpaceDimension * number of
  // of pixels in the grid region.
  return ( static_cast<unsigned int>( SpaceDimension ) *
	   static_cast<unsigned int>( m_GridRegion.GetNumberOfPixels() ) );	   
}


// Get the number of parameters per dimension
template<class TScalarType, unsigned int NDimensions, unsigned int VSplineOrder>
unsigned int 
BSplineField<TScalarType, NDimensions, VSplineOrder>
::GetNumberOfParametersPerDimension(void) const
{
  // The number of parameters per dimension equal number of
  // of pixels in the grid region.
  return ( static_cast<unsigned int>( m_GridRegion.GetNumberOfPixels() ) );

}
	   
// Set the grid region
template<class TScalarType, unsigned int NDimensions, unsigned int VSplineOrder>
void
BSplineField<TScalarType, NDimensions, VSplineOrder>
::SetGridRegion( const RegionType& region )
{

  if ( m_GridRegion != region )
    {

    m_GridRegion = region;

    // set regions for each coefficient and jacobian image
    for ( unsigned int j = 0; j < SpaceDimension; j++ )
      {
      m_WrappedImage[j]->SetRegions( m_GridRegion );
      }

    // Set the valid region
    // If the grid spans the interval [start, last].
    // The valid interval for evaluation is [start+offset, last-offset]
    // when spline order is even.
    // The valid interval for evaluation is [start+offset, last-offset)
    // when spline order is odd.
    // Where offset = vcl_floor(spline / 2 ).
    // Note that the last pixel is not included in the valid region
    // with odd spline orders.
    typename RegionType::SizeType size = m_GridRegion.GetSize();
    typename RegionType::IndexType index = m_GridRegion.GetIndex();

    for ( unsigned int j = 0; j < SpacePlusTimeDimension; j++ )
      {
      index[j] += 
	static_cast< typename RegionType::IndexValueType >( m_Offset );
      size[j] -= 
	static_cast< typename RegionType::SizeValueType> ( 2 * m_Offset );
      m_ValidRegionFirst[j] = index[j];
      m_ValidRegionLast[j] = index[j] +
	static_cast< typename RegionType::IndexValueType >( size[j] ) - 1;
      }
    m_ValidRegion.SetSize( size );
    m_ValidRegion.SetIndex( index );
    }
}


// Set the grid spacing
template<class TScalarType, unsigned int NDimensions, unsigned int VSplineOrder>
void
BSplineField<TScalarType, NDimensions, VSplineOrder>
::SetGridSpacing( const SpacingType& spacing )
{
  if ( m_GridSpacing != spacing )
    {
    m_GridSpacing = spacing;

    // set spacing for each coefficient and jacobian image
    for ( unsigned int j = 0; j < SpaceDimension; j++ )
      {
      m_WrappedImage[j]->SetSpacing( m_GridSpacing.GetDataPointer() );
      }

    DirectionType scale;
    for( unsigned int i=0; i<SpacePlusTimeDimension; i++)
      {
      scale[i][i] = m_GridSpacing[i];
      }

    m_IndexToPoint = m_GridDirection * scale;
    m_PointToIndex = m_IndexToPoint.GetInverse();

    this->Modified();
    }

}

// Set the grid direction
template<class TScalarType, unsigned int NDimensions, unsigned int VSplineOrder>
void
BSplineField<TScalarType, NDimensions, VSplineOrder>
::SetGridDirection( const DirectionType & direction )
{
  if ( m_GridDirection != direction )
    {
    m_GridDirection = direction;

    // set direction for each coefficient and jacobian image
    for ( unsigned int j = 0; j < SpaceDimension; j++ )
      {
      m_WrappedImage[j]->SetDirection( m_GridDirection );
      }

    DirectionType scale;
    for( unsigned int i=0; i<SpacePlusTimeDimension; i++)
      {
      scale[i][i] = m_GridSpacing[i];
      }

    m_IndexToPoint = m_GridDirection * scale;
    m_PointToIndex = m_IndexToPoint.GetInverse();

    this->Modified();
    }

}


// Set the grid origin
template<class TScalarType, unsigned int NDimensions, unsigned int VSplineOrder>
void
BSplineField<TScalarType, NDimensions, VSplineOrder>
::SetGridOrigin( const OriginType& origin )
{
  if ( m_GridOrigin != origin )
    {
    m_GridOrigin = origin;

    // set spacing for each coefficient and jacobianimage
    for ( unsigned int j = 0; j < SpaceDimension; j++ )
      {
      m_WrappedImage[j]->SetOrigin( m_GridOrigin.GetDataPointer() );
      }

    this->Modified();
    }

}


// Wrap flat parameters as images
template<class TScalarType, unsigned int NDimensions, unsigned int VSplineOrder>
void
BSplineField<TScalarType, NDimensions, VSplineOrder>
::WrapAsImages()
{

  /**
       * Wrap flat parameters array into SpaceDimension number of ITK images
       * NOTE: For efficiency, parameters are not copied locally. The parameters
       * are assumed to be maintained by the caller.
       */
  PixelType * dataPointer =
    const_cast<PixelType *>(( m_InputParametersPointer->data_block() ));
  unsigned int numberOfPixels = m_GridRegion.GetNumberOfPixels();

  for ( unsigned int j = 0; j < SpaceDimension; j++ )
    {
    m_WrappedImage[j]->GetPixelContainer()->
      SetImportPointer( dataPointer, numberOfPixels );
    dataPointer += numberOfPixels;
    m_CoefficientImage[j] = m_WrappedImage[j];
    }
    
}


// Set the B-Spline coefficients using input images
template<class TScalarType, unsigned int NDimensions, unsigned int VSplineOrder>
void 
BSplineField<TScalarType, NDimensions, VSplineOrder>
::SetCoefficientImage( ImagePointer images[] )
{
  if ( images[0] )
    {
    this->SetGridRegion( images[0]->GetBufferedRegion() );
    this->SetGridSpacing( images[0]->GetSpacing() );
    this->SetGridDirection( images[0]->GetDirection() );
    this->SetGridOrigin( images[0]->GetOrigin() );

    for( unsigned int j = 0; j < SpaceDimension; j++ )
      {
      m_CoefficientImage[j] = images[j];
      }

    // Clean up buffered parameters
    m_InputParametersPointer  = NULL;
    }
}  


// 
template<class TScalarType, unsigned int NDimensions, unsigned int VSplineOrder>
bool 
BSplineField<TScalarType, NDimensions, VSplineOrder>
::InsideValidRegion( 
  const ContinuousIndexType& index ) const
{
  bool inside = true;

#ifndef ITK_USE_CENTERED_PIXEL_COORDINATES_CONSISTENTLY
  if ( !m_ValidRegion.IsInside( index ) )
    {
    inside = false;
    }
#endif

  if ( inside && m_SplineOrderOdd )
    {
    typedef typename ContinuousIndexType::ValueType ValueType;
    for( unsigned int j = 0; j < SpacePlusTimeDimension; j++ )
      {
      if ( index[j] >= static_cast<ValueType>( m_ValidRegionLast[j] ) )
	{ 
	inside = false;
	break;
	}
#ifdef ITK_USE_CENTERED_PIXEL_COORDINATES_CONSISTENTLY
       if ( index[j] < static_cast<ValueType>( m_ValidRegionFirst[j] ) )
        {
        inside = false;
        break;
        }
 #endif
      }
    }

  return inside;
}

// Displace a point
template<class TScalarType, unsigned int NDimensions, unsigned int VSplineOrder>
void
BSplineField<TScalarType, NDimensions, VSplineOrder>
::GetPointVelocity( 
  const PointType & point, 
  SpacePointType & outputVelocity, 
  WeightsType & weights, 
  ParameterIndexArrayType & indices,
  bool& inside ) const
{

  //It is assumed that  m_CoefficientImage[0] is filled 

  unsigned int j;
  IndexType supportIndex;

  outputVelocity.Fill( NumericTraits<TScalarType>::Zero );

  inside = true;
  ContinuousIndexType index;

  this->TransformPointToContinuousIndex( point, index );

  // NOTE: if the support region does not lie totally within the grid
  // we assume zero displacement 
  inside = this->InsideValidRegion( index );
  if ( !inside )
    {
    return;
    }

  // Compute interpolation weights
  m_WeightsFunction->Evaluate( index, weights, supportIndex );
  

  // TO FIX UP
/*
  for ( j = 0; j < SpacePlusTimeDimension; j++ )
    {
    if ( supportIndex[j] < m_CoefficientImage[0]->GetBufferedRegion().GetIndex()[j] )
	{
	supportIndex[j] = m_CoefficientImage[0]->GetBufferedRegion().GetIndex()[j];
	}
    }
*/
  // !!! Debugging !!!

  for ( j = 0; j < SpacePlusTimeDimension; j++ )
    {
    if ( supportIndex[j] < m_CoefficientImage[0]->GetBufferedRegion().GetIndex()[j] )
	{
	std::cout << "THIS SHOULD NOT HAPPEN : index is " << index << " support index is " << supportIndex << " and buffered region starts at " << m_CoefficientImage[0]->GetBufferedRegion().GetIndex() << std::endl; 
	}
    }


  // For each dimension, correlate coefficient with weights
  RegionType supportRegion;
  supportRegion.SetSize( m_SupportSize );
  supportRegion.SetIndex( supportIndex );

  typedef ImageRegionConstIterator<ImageType> IteratorType;
  IteratorType m_Iterator[ SpaceDimension ];
  unsigned long counter = 0;
  const PixelType * basePointer = m_CoefficientImage[0]->GetBufferPointer();

  for ( j = 0; j < SpaceDimension; j++ )
    {
    m_Iterator[j] = IteratorType( m_CoefficientImage[j], supportRegion );
    m_Iterator[j].GoToBegin();
    }

  while ( ! m_Iterator[0].IsAtEnd() && counter < indices.Size() )
    {

    // multiply weigth with coefficient
    for ( j = 0; j < SpaceDimension; j++ )
      {
      outputVelocity[j] += static_cast<TScalarType>( 
	weights[counter] * m_Iterator[j].Get());
      }

    // populate the indices array
    indices[counter] = &(m_Iterator[0].Value()) - basePointer;


    // Increment iterator
    for ( j = 0; j < SpaceDimension; j++ )
      {
      ++( m_Iterator[j] );
      }

    // go to next coefficient in the support region
    ++ counter;


    }//END WHILE


}

// Compute the Jacobian in one position 
template<class TScalarType, unsigned int NDimensions, unsigned int VSplineOrder>
void
BSplineField<TScalarType, NDimensions,VSplineOrder>
::GetVelocityJacobian( const PointType & point, WeightsType& weights, ParameterIndexArrayType& indexes) const
{

  RegionType supportRegion;
  supportRegion.SetSize( m_SupportSize );
  const PixelType * basePointer = m_CoefficientImage[0]->GetBufferPointer();

  ContinuousIndexType index;

  this->TransformPointToContinuousIndex( point, index ); 

  // NOTE: if the support region does not lie totally within the grid
  // we assume zero displacement and return the input point
  if ( !this->InsideValidRegion( index ) )
    {
    weights.Fill(0.0);
    indexes.Fill(0);
    return;
    }

  // Compute interpolation weights
  IndexType supportIndex;

  m_WeightsFunction->Evaluate( index, weights, supportIndex );

  // For each dimension, copy the weight to the support region
  supportRegion.SetIndex( supportIndex );
  unsigned long counter = 0;

  typedef ImageRegionIterator<ImageType> IteratorType;
  IteratorType m_Iterator = IteratorType( m_CoefficientImage[0], supportRegion );


  while ( ! m_Iterator.IsAtEnd() )
    {


    indexes[counter] = &(m_Iterator.Value()) - basePointer;

    // go to next coefficient in the support region
    ++ counter;
    ++m_Iterator;

    }

}


template<class TScalarType, unsigned int NDimensions, unsigned int VSplineOrder>
void
BSplineField<TScalarType, NDimensions,VSplineOrder>
::TransformPointToContinuousIndex( const PointType & point, ContinuousIndexType & index ) const
{
  unsigned int j;

  Vector<double, SpacePlusTimeDimension> tvector;

  for ( j = 0; j < SpacePlusTimeDimension; j++ )
    {
    tvector[j] = point[j] - this->m_GridOrigin[j];
    }

  Vector<double, SpacePlusTimeDimension> cvector;

  cvector = m_PointToIndex * tvector;

  for ( j = 0; j < SpacePlusTimeDimension; j++ )
    {
    index[j] = static_cast< typename ContinuousIndexType::CoordRepType >( cvector[j] );
    }
}


template<class TScalarType, unsigned int NDimensions, unsigned int VSplineOrder>
unsigned int 
BSplineField<TScalarType, NDimensions,VSplineOrder>
::GetNumberOfAffectedWeights() const
{
  return m_WeightsFunction->GetNumberOfWeights();
}



// Compute classical jacobian
template<class TScalarType, unsigned int NDimensions, unsigned int VSplineOrder>
bool
BSplineField<TScalarType, NDimensions,VSplineOrder>
::GetClassicalJacobian
(const PointType & point, ClassicalJacobianType& jacobianMatrix, double timeStep) const
{

  // Clean jacobian matrix
  jacobianMatrix.fill(0.);

  RegionType supportRegion;
  supportRegion.SetSize( m_SupportSize );
  const PixelType * basePointer = m_CoefficientImage[0]->GetBufferPointer();

  ContinuousIndexType index;
  this->TransformPointToContinuousIndex( point, index );

  // Derivative weights and indexes
  const unsigned int numberOfWeights 
    =  m_weightsFunctionDerivative->GetNumberOfWeights();
  WeightsType weights
    ( numberOfWeights* SpacePlusTimeDimension); 
  weights.Fill(0.0);
  ParameterIndexArrayType 
    indexes( numberOfWeights );

  // NOTE: if the support region does not lie totally within the grid
  // we assume zero displacement and return the input point
  if ( !this->InsideValidRegion( index ) )
    {
    jacobianMatrix.set_identity();
    return false;
    }

  // Compute interpolation weights
  IndexType supportIndex;

  m_weightsFunctionDerivative->EvaluateDerivative( index, weights, supportIndex);

  // For each dimension, copy the weight to the support region
  supportRegion.SetIndex( supportIndex );
  unsigned long counter = 0;

  typedef ImageRegionIterator<ImageType> IteratorType;

  IteratorType m_Iterator = IteratorType( m_CoefficientImage[0], supportRegion );

  while ( ! m_Iterator.IsAtEnd() )
    {
    indexes[counter] = &(m_Iterator.Value()) - basePointer;

    // go to next coefficient in the support region
    ++ counter;
    ++m_Iterator;
    }

  // Number of grid points
  const unsigned int numberOfGridPoints 
    = this->GetNumberOfParametersPerDimension();

  for (unsigned int dir=0; dir<SpaceDimension; dir++)
    {
    for (unsigned int dim=0; dim<SpaceDimension; dim++) //
      {
      for (unsigned int counter=0; counter<indexes.Size(); counter++) 
    {

	unsigned int param_index = indexes[counter] + dim * numberOfGridPoints;

	if (param_index < this->GetNumberOfParameters())
	  {
    // The lines in the jacobian matrix are related to the component of the velocity field and the 
    // columns to the direction regarding which we compute the derivative	 
	  jacobianMatrix[dim][dir] += 1. / ( m_GridSpacing[dir] ) *  
	    weights[counter + dir * numberOfWeights] * this->GetParameters()[param_index];
	  }
	else
	  {
	  itkExceptionMacro( << "Attempting to access unvalid index in vector or parameters : " 
			     << "Index : " << param_index << ", vector size : "
			     << this->GetNumberOfParameters()
			     << std::endl );
	  }
	}//ROF counter

      // The transformation being defined as T(x) = v(x)*m_TimeStep + x, we need
      // to multiply by timeStep and add 1. on the diagonal of the velocity jacobian to get
      // the transformation jacobian.

      jacobianMatrix[dim][dir] *= timeStep;
      if (dir==dim)
        jacobianMatrix[dir][dim] += 1.;

      }// ROF dim
    }// ROF dir

  return true;
}


// Compute divergence of the field
template<class TScalarType, unsigned int NDimensions, unsigned int VSplineOrder>
bool
BSplineField<TScalarType, NDimensions,VSplineOrder>
::GetVelocityDivergence ( const PointType & inputPoint,
							double & div, 
							WeightsType & weights,
    					    ParameterIndexArrayType & indexes ) const
{

  RegionType supportRegion;
  supportRegion.SetSize( m_SupportSize );
  const PixelType * basePointer = m_CoefficientImage[0]->GetBufferPointer();

  ContinuousIndexType index;
  this->TransformPointToContinuousIndex( inputPoint, index );

  // Derivative weights and indexes
  const unsigned int numberOfWeights 
    =  m_weightsFunctionDerivative->GetNumberOfWeights();
    
  if (weights.Size() != numberOfWeights* SpacePlusTimeDimension)
    {
    weights.SetSize(numberOfWeights* SpacePlusTimeDimension);
    }
     
  weights.Fill(0.0);
  if (indexes.Size() != numberOfWeights)
    { 
    indexes.SetSize( numberOfWeights );
    }

  // NOTE: if the support region does not lie totally within the grid
  // we return a null divergence 
  if ( !this->InsideValidRegion( index ) )
    {
    div = 0.;
    return false;
    }

  // Compute interpolation weights
  IndexType supportIndex;

  m_weightsFunctionDerivative->EvaluateDerivative( index, weights, supportIndex);

  // Record the index values at this point
  supportRegion.SetIndex( supportIndex );
  unsigned long counter = 0;

  typedef ImageRegionIterator<ImageType> IteratorType;

  IteratorType m_Iterator = IteratorType( m_CoefficientImage[0], supportRegion );

  while ( ! m_Iterator.IsAtEnd() )
    {
    indexes[counter] = &(m_Iterator.Value()) - basePointer;

    // go to next coefficient in the support region
    ++ counter;
    ++m_Iterator;
    }

  // Number of grid points
  const unsigned int numberOfGridPoints 
    = this->GetNumberOfParametersPerDimension();

  // Computing divergence
  div = 0.;  
  for (unsigned int dir=0; dir<SpaceDimension; dir++)
    {

	double gridSpacing = this->GetGridSpacing()[dir];
    
    for (unsigned int index=0; index<indexes.Size(); index++)
		{
		unsigned int param_index = indexes[index] + dir * numberOfGridPoints;
		
		div += weights[index + dir * numberOfWeights] * this->GetParameters()[param_index] / gridSpacing;
		}
		
    }// ROF dir

  return true;
}

// Compute divergence of the field
template<class TScalarType, unsigned int NDimensions, unsigned int VSplineOrder>
bool
BSplineField<TScalarType, NDimensions,VSplineOrder>
::GetVelocityDivergenceHighOrder ( const PointType & inputPoint,
  WeightsType & weights) const
{
  ContinuousIndexType index;
  this->TransformPointToContinuousIndex( inputPoint, index );

  // Derivative weights and indexes
  const unsigned int numberOfWeights 
    =  m_weightsFunctionSecondDerivative->GetNumberOfWeights();

  if (weights.Size() != numberOfWeights * SpacePlusTimeDimension * SpacePlusTimeDimension)
    {
    weights.SetSize(numberOfWeights * SpacePlusTimeDimension * SpacePlusTimeDimension);
    }

  weights.Fill(0.0);
  // NOTE: if the support region does not lie totally within the grid
  // we return a null divergence 
  if ( !this->InsideValidRegion( index ) )
    {
    return false;
    }

  // Compute interpolation weights
  IndexType supportIndex;

  m_weightsFunctionSecondDerivative->EvaluateSecondDerivative( index, weights, supportIndex);

  return true;
}


} // end namespace itk

#endif


