/*=========================================================================

  Program:   Insight Segmentation & Registration Toolkit
  Module:    $RCSfile: itkDiffeomorphicBSplineTransform.txx,v $
  Language:  C++
  Date:      $Date: 2008-06-20 15:25:58 $
  Version:   $Revision: 1.2 $

  Copyright (c) Insight Software Consortium. All rights reserved.
  See ITKCopyright.txt or http://www.itk.org/HTML/Copyright.htm for details.

     This software is distributed WITHOUT ANY WARRANTY; without even 
     the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR 
     PURPOSE.  See the above copyright notices for more information.

=========================================================================*/
#ifndef _itkDiffeomorphicBSplineTransform_txx
#define _itkDiffeomorphicBSplineTransform_txx

#include "itkDiffeomorphicBSplineTransform.h"

namespace itk
{

// Constructor with default arguments
template<class TScalarType, unsigned int NDimensions, unsigned int BSplineOrder>
DiffeomorphicBSplineTransform<TScalarType, NDimensions, BSplineOrder>::
DiffeomorphicBSplineTransform():Superclass(SpaceDimension, 0)
{
  m_InternalTransforms.resize(0);
  m_GridSpacing.Fill(0.);
  m_GridOrigin.Fill(0.);

  typename RegionType::SizeType dummySize;
  dummySize.Fill(0);
  m_GridRegion.SetSize(dummySize);
  typename RegionType::IndexType dummyIndex;
  dummyIndex.Fill(0);  	 
  m_GridRegion.SetIndex(dummyIndex);

  m_InternalTransformsInitialized = false;

  m_NumberOfTimeSteps=0;

  m_InternalJacobians.resize(0);
  m_InternalWeights.resize(0);
  m_InternalIndexes.resize(0);

  m_InternalParameters.resize(0);


  this->m_FixedParameters.SetSize ( 1 + NDimensions * 3 );
  this->m_FixedParameters.Fill ( 0.0 );


}

// Destructor
template<class TScalarType, unsigned int NDimensions, unsigned int BSplineOrder>
DiffeomorphicBSplineTransform<TScalarType, NDimensions, BSplineOrder>::
~DiffeomorphicBSplineTransform()
{
  return;
}

// Set the parameters
template <class TScalarType, unsigned int NDimensions, unsigned int BSplineOrder>
void
DiffeomorphicBSplineTransform<TScalarType, NDimensions, BSplineOrder>
::InitializeInternalTransforms()
{
  if (m_NumberOfTimeSteps == 0)
    itkExceptionMacro(<< "Number of time steps has not been set correcly." );
  
  if (!m_InternalTransformsInitialized)
    {
    if (m_InternalTransforms.size() != m_NumberOfTimeSteps)
      m_InternalTransforms.resize(m_NumberOfTimeSteps, 0);
    
    for (unsigned int index=0; index<m_InternalTransforms.size(); index++)
      {
      if (m_InternalTransforms.at(index).IsNull())
	{
	m_InternalTransforms.at(index) = InternalTransformType::New();
	m_InternalTransforms.at(index)->SetGridRegion(m_GridRegion);
	m_InternalTransforms.at(index)->SetGridSpacing(m_GridSpacing);
	m_InternalTransforms.at(index)->SetGridOrigin(m_GridOrigin);
	}
      }
    
    const unsigned int numberOfWeights 
      = m_InternalTransforms[0]->GetNumberOfWeights();
    
    // m_InternalJacobians stores intermediate results for computing
    // the whole transform jacobian. Each matrix in the vector has the
    // (SpaceDimension, numberOfWeights) dimensions
    if (m_InternalJacobians.size() != m_NumberOfTimeSteps)
      {
      m_InternalJacobians.resize( m_NumberOfTimeSteps, 
				  InternalJacobianType( SpaceDimension,
							SpaceDimension, 0.) );
      m_InternalWeights.resize( m_NumberOfTimeSteps,
				WeightsType(numberOfWeights));
      
      m_InternalIndexes.resize( m_NumberOfTimeSteps,
				ParameterIndexArrayType(numberOfWeights));
      
      m_InternalPoints.resize( m_NumberOfTimeSteps );
      }
    
    if (this->m_Jacobian.cols() != this->GetNumberOfParameters())
      {
      this->m_Jacobian.SetSize(SpaceDimension,this->GetNumberOfParameters());
      }
  
    }
  
  m_InternalTransformsInitialized = true;
  this->Modified();
  
}

// Insert an additional transform before pos
template <class TScalarType, unsigned int NDimensions, unsigned int BSplineOrder>
void
DiffeomorphicBSplineTransform<TScalarType, NDimensions, BSplineOrder>
::InsertTransform(unsigned int pos, bool reshapeParameters)
{
  
  // Increase number of time steps
  m_NumberOfTimeSteps++;

  if (pos>= m_InternalTransforms.size())
    {
    itkExceptionMacro( << "Uncorrect position for transformation insertion."
		       << std::endl);
    }
  
  // Allocate internal stuff
  m_InternalTransformsInitialized = false;
  InitializeInternalTransforms();
  
  if (reshapeParameters)
    {
    // Reshape set of parameters
    const ParametersType* currentParameters = m_InputParametersPointer;
    ParametersType* newParameters = new ParametersType(this->GetNumberOfParameters());
    newParameters->Fill(0.);
    
    // Copy of old parameters in the new array
    const typename ParametersType::ValueType* data_block_current = currentParameters->data_block();
    typename ParametersType::ValueType* data_block_new = newParameters->data_block();
    
    unsigned int numberOfParametersByTransform = m_InternalTransforms[0]->GetNumberOfParameters();
    
    // Copy first part of parameters (before pos)
    unsigned int index=0;
    while ( index<numberOfParametersByTransform*pos )
      {
      *data_block_new = *data_block_current;
      
      data_block_current++;
      data_block_new++;
      index++;
      }
    
    // Skip in new parameters the parameters of the added transform
    data_block_new += numberOfParametersByTransform;
    
    // Copy segond part (after pos)
    while (index<currentParameters->Size())
      {
      *data_block_new = *data_block_current;
      
      data_block_current++;
      data_block_new++;
      index++;
      }
    
    m_InputParametersPointer = newParameters;
    WrapParameters();
    
    }
}

// Set the parameters by value
template <class TScalarType, unsigned int NDimensions, unsigned int BSplineOrder>
void
DiffeomorphicBSplineTransform<TScalarType, NDimensions, BSplineOrder>
::SetParametersByValue( const ParametersType & parameters )
{
  //  Copy parameters instead of keeping a reference
  m_InternalParametersBuffer = parameters;
  m_InputParametersPointer = &m_InternalParametersBuffer;

  WrapParameters();
}


// Set the parameters
template <class TScalarType, unsigned int NDimensions, unsigned int BSplineOrder>
void
DiffeomorphicBSplineTransform<TScalarType, NDimensions, BSplineOrder>
::SetParameters( const ParametersType & parameters )
{
  // Keep a reference to the input parameters
  m_InputParametersPointer = &parameters;
  WrapParameters();
}

// Wrap parameters 
template <class TScalarType, unsigned int NDimensions, unsigned int BSplineOrder>
void
DiffeomorphicBSplineTransform<TScalarType, NDimensions, BSplineOrder>
::WrapParameters( )
{
  // Check if internal parameters has been allocated or not
  if (m_InternalParameters.size() != m_NumberOfTimeSteps)
    m_InternalParameters.resize(m_NumberOfTimeSteps);

  // Check if internal transforms have been initialized
  if (!m_InternalTransformsInitialized)
    InitializeInternalTransforms();
    
  const double* ParametersRawBuffer = m_InputParametersPointer->data_block();
  for (unsigned int index=0; index<m_NumberOfTimeSteps; index++)
    {

    typename InternalTransformType::ParametersType 
      params(ParametersRawBuffer,
        m_InternalTransforms.at(index)->GetNumberOfParameters());

    m_InternalParameters[index].SetData(params.data_block(), 
                                m_InternalTransforms.at(index)->GetNumberOfParameters());

    // Acutal setting of parameters
    m_InternalTransforms.at(index)->SetParameters(m_InternalParameters[index]);

    //Moving to the next parameters block
    ParametersRawBuffer += 
      m_InternalTransforms.at(index)->GetNumberOfParameters();
    }

  this->Modified();   
}

// Get the parameters
template <class TScalarType, unsigned int NDimensions, unsigned int BSplineOrder>
const typename DiffeomorphicBSplineTransform<TScalarType, NDimensions, BSplineOrder>::ParametersType &
DiffeomorphicBSplineTransform<TScalarType, NDimensions, BSplineOrder>
::GetParameters( void ) const
{  
  return (*m_InputParametersPointer);
}

// Print self
template<class TScalarType, unsigned int NDimensions, unsigned int BSplineOrder>
void
DiffeomorphicBSplineTransform<TScalarType, NDimensions, BSplineOrder>::
PrintSelf(std::ostream &os, Indent indent) const
{
  Superclass::PrintSelf(os,indent);
}

// Transform a point
template<class TScalarType, unsigned int NDimensions, unsigned int BSplineOrder>
typename DiffeomorphicBSplineTransform<TScalarType, NDimensions, BSplineOrder>::OutputPointType
DiffeomorphicBSplineTransform<TScalarType, NDimensions, BSplineOrder>::
TransformPoint(const InputPointType &point) const 
{

  // Check if internal transforms have been initialized
  if (!m_InternalTransformsInitialized)
    itkExceptionMacro(<< "Parameters have not been set correctly." << std::endl
		      << "Please use SetParameters() or SetIdentity prior to "
		      << "calling this method" << std::endl);
  
  // Chaining internal transforms
  typename InternalTransformType::OutputPointType currentPoint = point;
  
  for (unsigned int index=0;index<m_NumberOfTimeSteps; index++)
    {
    currentPoint = m_InternalTransforms[index]->TransformPoint(currentPoint);
    }
  
  return currentPoint;
}


// Transform a point
template<class TScalarType, unsigned int NDimensions, unsigned int BSplineOrder>
const std::vector<typename DiffeomorphicBSplineTransform<TScalarType, NDimensions, BSplineOrder>::OutputPointType> &
DiffeomorphicBSplineTransform<TScalarType, NDimensions, BSplineOrder>::
TransformPoint(const InputPointType  &point, bool& valid) const
{
  
  // Check if internal transforms have been initialized
  if (!m_InternalTransformsInitialized)
    itkExceptionMacro(<< "Parameters have not been set correctly." << std::endl
		      << "Please use SetParameters() or SetIdentity prior to "
		      << "calling this method" << std::endl);
  
  // Chaining internal transforms
  const typename InternalTransformType::OutputPointType* currentPoint = & point;
  
  for (unsigned int index=0; index < m_NumberOfTimeSteps; index++)
    {
    m_InternalPoints[index] = m_InternalTransforms[index]->TransformPoint(*currentPoint);
    currentPoint =  & m_InternalPoints[index];
    }
  
  valid = true;

  return m_InternalPoints;

}

// Transform a point
template<class TScalarType, unsigned int NDimensions, unsigned int BSplineOrder>
typename DiffeomorphicBSplineTransform<TScalarType, NDimensions, BSplineOrder>::OutputPointType
DiffeomorphicBSplineTransform<TScalarType, NDimensions, BSplineOrder>::
TransformPoint(const InputPointType &point, unsigned int timeStep) const 
{

  // Check if internal transforms have been initialized
  if (!m_InternalTransformsInitialized)
    itkExceptionMacro(<< "Parameters have not been set correctly." << std::endl
		      << "Please use SetParameters() or SetIdentity prior to "
		      << "calling this method" << std::endl);
  
  if (timeStep > m_NumberOfTimeSteps)
    itkExceptionMacro(<< "Incorrect time step." << std::endl);
  
  // Chaining internal transforms
  const typename InternalTransformType::OutputPointType* currentPoint = & point;
  for (unsigned int index=0;index<timeStep; index++)
    {
    m_InternalPoints[index] = m_InternalTransforms[index]->TransformPoint(*currentPoint);
    currentPoint=&m_InternalPoints[index];
    }
  
  return * currentPoint;
  
}

// Reset internal jacobian to the identity
template<class TScalarType, unsigned int NDimensions, unsigned int BSplineOrder>
void
DiffeomorphicBSplineTransform< TScalarType, NDimensions, BSplineOrder >::
ResetAllPhysicalJacobians()
{
  
  for (unsigned int transfo=0;transfo<m_InternalTransforms.size(); transfo++)
    {
    m_InternalJacobians[transfo].set_identity();
    }
}

// Accumulate physical jacobians from at time point up to a given time
// point
template<class TScalarType, unsigned int NDimensions, unsigned int BSplineOrder>
const std::vector<typename DiffeomorphicBSplineTransform<TScalarType, NDimensions, BSplineOrder>
		  ::InternalJacobianType> &
DiffeomorphicBSplineTransform< TScalarType, NDimensions, BSplineOrder >::
AccumulatePhysicalJacobiansInTimeInterval(const InputPointType & point_time_t1, 
					  unsigned int time1, unsigned int time2) const
{
  
  // Check if internal transforms have been initialized
  if (!m_InternalTransformsInitialized)
    itkExceptionMacro(<< "Parameters have not been set correctly." << std::endl
		      << "Please use SetParameters() or SetIdentity prior to "
		      << "calling this method" << std::endl);
 
 
  // Check that the time interval is valid
  if ( ( time1 <= time2 ) && ( time1 < m_InternalTransforms.size())
    && ( time2 <= m_InternalTransforms.size()) )
    {
    // Chaining internal transforms and pre-computing information
    // required for jacobian computation.
    typename InternalTransformType::OutputPointType currentPoint = point_time_t1;

    for (unsigned int transfo=time1;transfo<time2; transfo++)
      {
      // Get the ITK jacobian of this transformation
 m_InternalTransforms[transfo]
   ->GetJacobian( currentPoint, m_InternalWeights[transfo],m_InternalIndexes[transfo]);


      // Get the classical jacobian of this transformation
      bool valid = false;
      valid = m_InternalTransforms[transfo]
        ->GetClassicalJacobian( currentPoint, m_InternalJacobians[transfo] );

      // Compute accumulation of classical jacobians for all transformations
      for (unsigned int transfo2=time1+1; transfo2<transfo; transfo2++)
        {
        m_InternalJacobians[transfo2] 
          = m_InternalJacobians[transfo] * m_InternalJacobians[transfo2];
        }

      // Apply the transformation
      currentPoint = m_InternalTransforms[transfo]->TransformPoint(currentPoint);
      }
    }
  else
    {
    itkExceptionMacro( << "Time interval is incorrect" << std::endl);
    }

  return m_InternalJacobians;

}

// Compute the Jacobian in one position 
template<class TScalarType, unsigned int NDimensions, unsigned int BSplineOrder>
const typename DiffeomorphicBSplineTransform<TScalarType, NDimensions, BSplineOrder>::JacobianType & 
DiffeomorphicBSplineTransform< TScalarType, NDimensions, BSplineOrder >::
GetJacobian( const InputPointType & point, unsigned int time, 
	     bool accumulate, std::vector<long unsigned int>& jacoColumnIndexes) const
{
  // Check if internal transforms have been initialized
  if (!m_InternalTransformsInitialized)
    itkExceptionMacro(<< "Parameters have not been set correctly." << std::endl
		      << "Please use SetParameters() or SetIdentity prior to "
		      << "calling this method" << std::endl);

  this->m_Jacobian.Fill( 0.0 );
  jacoColumnIndexes.resize( 0 );

  // Calling phys. jacobians chaining function for all time points
  // prior to time.
  if (accumulate)
    AccumulatePhysicalJacobiansInTimeInterval( point, 0, time );


  // We know have to fill the jacobian matrix
  // The jacobian for a specific transform in the chain is the product
  // of all classical jacobians upstream and the ITK jacobian of this
  // transform. Since for BSplineDeformable transforms, the jacobian
  // has the form of a sparse matrix where non-null columns have only
  // one non-null element equal to the BSpline weight, we can simplify
  // this product.
  for ( unsigned int transfo=0;
	transfo <time;
	transfo++ )
    {
    // For each transform, the m_InternalIndexes stores the indexes of
    // non-null columns in the jacobian
    for (unsigned int index=0; index<m_InternalIndexes[transfo].Size(); index++)
      {
      unsigned int controlPointIndex = m_InternalIndexes[transfo][index];
      typename WeightsType::ValueType weight 
        = m_InternalWeights[transfo][index];
    
      /*
      if ( (transfo==0)&&(weight>1.e-8) )
          std::cout<<"Weight for transform 0: "<<weight<<std::endl;
      */
      
      // Copy info in the final jacobian
      for (unsigned int dim=0; dim<SpaceDimension; dim++)
        {
        // Position of the current transformation parameters in the
        // concatenation of all transformation parameters
        unsigned int fullJacoColumn 
          = transfo * m_InternalTransforms[0]->GetNumberOfParameters();
	
        // Each control point acts on the three dimensions. The offset
        // between each dimension is the number of parameters per
        // dimension for one transform in the chain
        fullJacoColumn += controlPointIndex 
          + dim * m_InternalTransforms[0]->GetNumberOfParametersPerDimension();
	
	      jacoColumnIndexes.push_back(fullJacoColumn);
	
        if (transfo == time - 1 )
          {
          // Treat the case of the last transformation 
          // which does not need any mutplication by 
          // downstream jacobians
          this->m_Jacobian[dim][fullJacoColumn]
            = weight;
          }
        else
          {
          for (unsigned int row=0; row<SpaceDimension; row++)
            {
	          this->m_Jacobian[row][fullJacoColumn]
              = m_InternalJacobians[transfo + 1][row][dim] * weight;
            }
          }
        }//ROF dim
      }//ROF control point index
    }// ROF transform index
  
  return this->m_Jacobian;
  
}


// Compute the Jacobian in one position 
template<class TScalarType, unsigned int NDimensions, unsigned int BSplineOrder>
const typename DiffeomorphicBSplineTransform<TScalarType, NDimensions, BSplineOrder>::JacobianType & 
DiffeomorphicBSplineTransform< TScalarType, NDimensions, BSplineOrder >::
GetJacobian( const InputPointType & point) const
{
  
  // Check if internal transforms have been initialized
  if (!m_InternalTransformsInitialized)
    itkExceptionMacro(<< "Parameters have not been set correctly." << std::endl
      << "Please use SetParameters() or SetIdentity prior to "
      << "calling this method" << std::endl);
  
  
  this->m_Jacobian.Fill( 0.0 );
  
  // Calling phys. jacobians chaining function for the entire time
  // range. This will update as well m_InternalIndexes for all
  // transformation between time steps time1 and time2
  AccumulatePhysicalJacobiansInTimeInterval( point, 0, this->m_InternalTransforms.size() );

  // We know have to fill the jacobian matrix
  // The jacobian for a specific transform in the chain is the product
  // of all classical jacobians upstream and the ITK jacobian of this
  // transform. Since for BSplineDeformable transforms, the jacobian
  // has the form of a sparse matrix where non-null columns have only
  // one non-null element equal to the BSpline weight, we can simplify
  // this product.
  for (unsigned int transfo=0;
    transfo < m_InternalTransforms.size(); 
    transfo++)
    {
    // For each transform, the m_InternalIndexes stores the indexes of
    // non-null columns in the jacobian
    for (unsigned int index=0; index<m_InternalIndexes[transfo].Size(); index++)
      {
      unsigned int controlPointIndex = m_InternalIndexes[transfo][index];
      typename WeightsType::ValueType weight 
        = m_InternalWeights[transfo][index];

      // Copy info in the final jacobian
      for (unsigned int dim=0; dim<SpaceDimension; dim++)
        {
        // Position of the current transformation parameters in the
        // concatenation of all transformation parameters
        unsigned int fullJacoColumn 
          = transfo * m_InternalTransforms[0]->GetNumberOfParameters();
        // Each control point acts on the three dimensions. The offset
        // between each dimension is the number of parameters per
        // dimension for one transform in the chain
        fullJacoColumn += controlPointIndex 
          + dim * m_InternalTransforms[0]->GetNumberOfParametersPerDimension();

        if (transfo == m_InternalTransforms.size() - 1 )
          {
          // Treat the case of the last transformation 
          // which does not need any mutplication by 
          // downstream jacobians
          this->m_Jacobian[dim][fullJacoColumn]
            = weight;
          }
        else
          {
          for (unsigned int row=0; row<SpaceDimension; row++)
            {
           this->m_Jacobian[row][fullJacoColumn]
              = m_InternalJacobians[transfo + 1][row][dim] * weight;
            }
          }
        }//ROF dim
      }//ROF control point index
    }// ROF transform index

  return this->m_Jacobian;
}

// Set the parameters for an Identity transform of this class
template<class TScalarType, unsigned int NDimensions, unsigned int BSplineOrder>
void
DiffeomorphicBSplineTransform<TScalarType, NDimensions, BSplineOrder>::
SetIdentity()
{
  
  // Check if internal transforms have been initialized
  if (!m_InternalTransformsInitialized)
    InitializeInternalTransforms();
  
  for (unsigned int index=0;index<m_InternalTransforms.size(); index++)
    {
    m_InternalTransforms[index]->SetIdentity();
    m_InputParametersPointer->Fill(0.);
    }
}

// Set the Fixed Parameters
template<class TScalarType, unsigned int NDimensions, unsigned int BSplineOrder>
void
DiffeomorphicBSplineTransform< TScalarType, NDimensions, BSplineOrder >
::SetFixedParameters( const ParametersType & parameters )
{
 
  // check if the number of parameters match the
  // expected number of parameters
  if ( parameters.Size() != 1+NDimensions*3 )
    {
    itkExceptionMacro(<< "Mismatched between parameters size "
                      << parameters.size() 
                      << " and number of fixed parameters "
                      << 1 + NDimensions*3 );
    }

  /********************************************************** 
    Fixed Parameters store the following information:
        Grid Size
        Grid Origin
        Grid Spacing
     The size of these is equal to the  NInputDimensions
  **********************************************************/
  
  // Number of time steps is the first fixed parameter recorded
  m_NumberOfTimeSteps=static_cast<unsigned int>(parameters[0]);

  /*** Set the Grid Parameters ***/
  typename InternalTransformType::SizeType   gridSize;
  for (unsigned int i=0;i<NDimensions;i++)
    {
    gridSize[i] = static_cast<int> (parameters[1+i]);
    }
  RegionType bsplineRegion;
  bsplineRegion.SetSize( gridSize );
  
  /*** Set the Origin Parameters ***/
  OriginType origin;
  for (unsigned int i=0;i<NDimensions;i++)
    {
    origin[i] = parameters[1+NDimensions+i];
    }
  
  /*** Set the Spacing Parameters ***/
  SpacingType spacing;
  for (unsigned int i=0;i<NDimensions;i++)
    {
    spacing[i] = parameters[1+2*NDimensions+i];
    }

  
  this->SetGridSpacing( spacing );
  this->SetGridOrigin( origin );
  this->SetGridRegion( bsplineRegion );

  this->Modified();
}


// Get the parameters
template<class TScalarType, unsigned int NDimensions, unsigned int BSplineOrder>
const 
typename DiffeomorphicBSplineTransform< TScalarType, NDimensions, BSplineOrder >
::ParametersType &
DiffeomorphicBSplineTransform< TScalarType, NDimensions, BSplineOrder >
::GetFixedParameters( void ) const
{
  RegionType resRegion = this->GetGridRegion(  );
  
  this->m_FixedParameters[0] = this->GetNumberOfTimeSteps();//m_NumberOfTimeSteps;

  for (unsigned int i=0;i<NDimensions;i++)
    {
    this->m_FixedParameters[1+i] = (resRegion.GetSize())[i];
    }
  for (unsigned int i=0;i<NDimensions;i++)
    {
    this->m_FixedParameters[1+NDimensions+i] = (this->GetGridOrigin())[i];
    } 
  for (unsigned int i=0;i<NDimensions;i++)
    {
    this->m_FixedParameters[1+2*NDimensions+i] =  (this->GetGridSpacing())[i];
    }
  
  return (this->m_FixedParameters);
}



}// namespace

#endif
