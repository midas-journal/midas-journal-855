/*=========================================================================
 
 Program:   Insight Segmentation & Registration Toolkit
 Module:    $RCSfile: regDiffeomorphicContinuousBSplineTransform.txx,v $
 Language:  C++
 Date:      $Date: 2008-05-08 23:22:35 $
 Version:   $Revision: 1.41 $
 
 Copyright (c) Insight Software Consortium. All rights reserved.
 See ITKCopyright.txt or http://www.itk.org/HTML/Copyright.htm for details.
 
 This software is distributed WITHOUT ANY WARRANTY; without even 
 the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR 
 PURPOSE.  See the above copyright notices for more information.
 
 =========================================================================*/
#ifndef __regDiffeomorphicContinuousBSplineTransform_txx
#define __regDiffeomorphicContinuousBSplineTransform_txx

#include "regDiffeomorphicContinuousBSplineTransform.h"
#include "itkContinuousIndex.h"
#include "itkImageRegionIteratorWithIndex.h"
#include "itkImageRegionIterator.h"
#include "itkImageRegionConstIterator.h"
#include "itkImageRegionConstIteratorWithIndex.h"
#include "itkIdentityTransform.h"

namespace itk
{
  
  // Constructor with default arguments
  template<class TScalarType, unsigned int NDimension, unsigned int VSplineOrder>
  DiffeomorphicContinuousBSplineTransform<TScalarType, NDimension, VSplineOrder>
  ::DiffeomorphicContinuousBSplineTransform(void):Superclass(SpacePlusTimeDimension,0)
  {
    
    // Velocity object
    m_Velocity = VelocityType::New();
    
    // Instantiate an identity transform
    typedef IdentityTransform<ScalarType, SpaceDimension> IdentityTransformType;
    typename IdentityTransformType::Pointer id = IdentityTransformType::New();
    m_BulkTransform = id;
    
    m_TimeStep = 0.25;//1.0;
    m_MinimumTimeStep=1.0;
    
    m_InternalParametersBuffer = ParametersType(0);
    // Make sure the parameters pointer is not NULL after construction.
    m_Velocity->SetParameters( m_InternalParametersBuffer );
    
    /** Fixed Parameters store the following information:
     *     Grid Size
     *     Grid Origin
     *     Grid Spacing
     *     Grid Direction
     *  The size of these is equal to the  NInputDimensions
     */
    
    this->m_FixedParameters.SetSize ( SpacePlusTimeDimension * (SpacePlusTimeDimension + 3) );  
    this->m_FixedParameters.Fill ( 0.0 );
    
    for (unsigned int i=0; i<SpacePlusTimeDimension; i++)
    {
      this->m_FixedParameters[2*SpacePlusTimeDimension+i] = m_Velocity->GetGridSpacing()[i];
    }
    
    for (unsigned int di=0; di<SpacePlusTimeDimension; di++)
    {
      for (unsigned int dj=0; dj<SpacePlusTimeDimension; dj++)
      {
        this->m_FixedParameters[3*SpacePlusTimeDimension+(di*SpacePlusTimeDimension+dj)] = m_Velocity->GetGridDirection()[di][dj];
      }
    }
    
    /** Cache variables for advoiding frequent reallocation. */
    incrJaco_columnOut = SparseJacobianColumnType(SpaceDimension, 0.);
    incrJaco_weights = typename VelocityType::WeightsType( m_Velocity->GetWeightsFunction()->GetNumberOfWeights() );
    incrJaco_physJacobian = typename VelocityType::ClassicalJacobianType(SpaceDimension,SpaceDimension);
    incrJaco_indices = typename VelocityType::ParameterIndexArrayType( m_Velocity->GetWeightsFunction()->GetNumberOfWeights() );
    
  }
  
  // Destructor
  template<class TScalarType, unsigned int NDimensions, unsigned int VSplineOrder>
  DiffeomorphicContinuousBSplineTransform<TScalarType, NDimensions, VSplineOrder>
  ::~DiffeomorphicContinuousBSplineTransform()
  {
  }
  
  
  // Name of the transform.
  template<class TScalarType, unsigned int NDimensions, unsigned int VSplineOrder>
  std::string 
  DiffeomorphicContinuousBSplineTransform<TScalarType, NDimensions, VSplineOrder>
  ::GetTransformTypeAsString() const
  {
    OStringStream n;
    n << GetNameOfClass();
    n << "_";
    if ( typeid ( TScalarType ) == typeid ( float ) )
    {
      n << "float";
    }
    else if ( typeid ( TScalarType ) == typeid ( double ) )
    {
      n << "double";
    }
    else
    {
      n << "other";
    }
    n << "_" << SpacePlusTimeDimension << "_" << SplineOrder;
    return n.str();
  }
  
  
  // Get the number of parameters
  template<class TScalarType, unsigned int NDimensions, unsigned int VSplineOrder>
  unsigned int
  DiffeomorphicContinuousBSplineTransform<TScalarType, NDimensions, VSplineOrder>
  ::GetNumberOfParameters(void) const
  {
    return m_Velocity->GetNumberOfParameters();
  }
  
  
  // Get the number of parameters per dimension
  template<class TScalarType, unsigned int NDimensions, unsigned int VSplineOrder>
  unsigned int
  DiffeomorphicContinuousBSplineTransform<TScalarType, NDimensions, VSplineOrder>
  ::GetNumberOfParametersPerDimension(void) const
  {
    return m_Velocity->GetNumberOfParametersPerDimension();
  }
  
  
  // Set the grid region
  template<class TScalarType, unsigned int NDimensions, unsigned int VSplineOrder>
  void
  DiffeomorphicContinuousBSplineTransform<TScalarType, NDimensions, VSplineOrder>
  ::SetGridRegion( const RegionType& region )
  {
    m_Velocity->SetGridRegion(region);
    
    // Input parameters point to internal buffer => using default parameters.
    if ( & m_Velocity->GetParameters() == & m_InternalParametersBuffer)
    {
      // Check if we need to resize the default parameter buffer.
      if ( m_InternalParametersBuffer.GetSize() != this->GetNumberOfParameters() )
      {
        m_InternalParametersBuffer.SetSize( this->GetNumberOfParameters() );
        // Fill with zeros for identity.
        m_InternalParametersBuffer.Fill( 0 );
      }
    }
    
    /**
     * Allocate memory for Jacobian and wrap into SpaceDimension number
     * of ITK images
     */
    this->m_Jacobian.set_size( SpaceDimension, this->GetNumberOfParameters() );
    this->m_Jacobian.Fill( NumericTraits<PixelType>::Zero );
    
    // Allocating total jacobian  
    m_TotalJacobian.set_size( SpaceDimension, this->GetNumberOfParameters() );
    this->m_TotalJacobian.Fill( NumericTraits<PixelType>::Zero );
    
    this->Modified();
  }
  
  
  // Set the grid spacing
  template<class TScalarType, unsigned int NDimensions, unsigned int VSplineOrder>
  void
  DiffeomorphicContinuousBSplineTransform<TScalarType, NDimensions, VSplineOrder>
  ::SetGridSpacing( const SpacingType& spacing )
  {
    m_Velocity->SetGridSpacing(spacing);
  }
  
  // Set the grid direction
  template<class TScalarType, unsigned int NDimensions, unsigned int VSplineOrder>
  void
  DiffeomorphicContinuousBSplineTransform<TScalarType, NDimensions, VSplineOrder>
  ::SetGridDirection( const DirectionType & direction )
  {
    m_Velocity->SetGridDirection( direction );
    
  }
  
  
  // Set the grid origin
  template<class TScalarType, unsigned int NDimensions, unsigned int VSplineOrder>
  void
  DiffeomorphicContinuousBSplineTransform<TScalarType, NDimensions, VSplineOrder>
  ::SetGridOrigin( const OriginType& origin )
  {
    m_Velocity->SetGridOrigin( origin );
  }
  
  // Get the grid origin
  template<class TScalarType, unsigned int NDimensions, unsigned int VSplineOrder>
  typename DiffeomorphicContinuousBSplineTransform<TScalarType, NDimensions, VSplineOrder>::OriginType
  DiffeomorphicContinuousBSplineTransform<TScalarType, NDimensions, VSplineOrder>
  ::GetGridOrigin()
  {
    return m_Velocity->GetGridOrigin();
  }
  
  // Get the grid spacing
  template<class TScalarType, unsigned int NDimensions, unsigned int VSplineOrder>
  typename DiffeomorphicContinuousBSplineTransform<TScalarType, NDimensions, VSplineOrder>::SpacingType
  DiffeomorphicContinuousBSplineTransform<TScalarType, NDimensions, VSplineOrder>
  ::GetGridSpacing()
  {
    return m_Velocity->GetGridSpacing();
  }
  
  // Get the grid region
  template<class TScalarType, unsigned int NDimensions, unsigned int VSplineOrder>
  typename DiffeomorphicContinuousBSplineTransform<TScalarType, NDimensions, VSplineOrder>::RegionType
  DiffeomorphicContinuousBSplineTransform<TScalarType, NDimensions, VSplineOrder>
  ::GetGridRegion()
  {
    return m_Velocity->GetGridRegion();
  }
  
  // Get Coefficient Image
  template<class TScalarType, unsigned int NDimensions, unsigned int VSplineOrder>
  const 
  typename DiffeomorphicContinuousBSplineTransform<TScalarType, NDimensions, VSplineOrder>::ImagePointer * 
  DiffeomorphicContinuousBSplineTransform<TScalarType, NDimensions, VSplineOrder>
  ::GetCoefficientImage() const
  {
    return m_Velocity->GetCoefficientImage();
  }
  
  
  
  // Reset transform to identity
  template<class TScalarType, unsigned int NDimensions, unsigned int VSplineOrder>
  void
  DiffeomorphicContinuousBSplineTransform<TScalarType, NDimensions, VSplineOrder>
  ::SetIdentity()
  {
    // Check if we need to resize the default parameter buffer.
    if ( m_InternalParametersBuffer.GetSize() != this->GetNumberOfParameters() )
    {
      m_InternalParametersBuffer.SetSize( this->GetNumberOfParameters() );
    }
    
    // Fill with zeros for identity.
    m_InternalParametersBuffer.Fill( 0 );
    
    // Pass parameters to velocity
    m_Velocity->SetParameters( &m_InternalParametersBuffer );
    
  }
  
  // Set the parameters
  template<class TScalarType, unsigned int NDimensions, unsigned int VSplineOrder>
  void
  DiffeomorphicContinuousBSplineTransform<TScalarType, NDimensions, VSplineOrder>
  ::SetParameters( const ParametersType & parameters )
  {
    
    // check if the number of parameters match the
    // expected number of parameters
    if ( parameters.Size() != this->GetNumberOfParameters() )
    {
      itkExceptionMacro(<<"Mismatched between parameters size "
                        << parameters.size() 
                        << " and number of parameters "
                        <<  this->GetNumberOfParameters() );
    }
    
    // Clean up buffered parameters
    m_InternalParametersBuffer = ParametersType( 0 );
    
    // Keep a reference to the input parameters
    m_Velocity->SetParameters( parameters );
    
    // Wrap flat array as images of coefficients
    m_Velocity->WrapAsImages();
    
    // Modified is always called since we just have a pointer to the
    // parameters and cannot know if the parameters have changed.
    this->Modified();
  }
  
  
  // Set the parameters
  template<class TScalarType, unsigned int NDimensions, unsigned int VSplineOrder>
  void
  DiffeomorphicContinuousBSplineTransform<TScalarType, NDimensions, VSplineOrder>
  ::SetParametersByValue( const ParametersType & parameters )
  {
    
    // check if the number of parameters match the
    // expected number of parameters
    if ( parameters.Size() != this->GetNumberOfParameters() )
    {
      itkExceptionMacro(<<"Mismatched between parameters size "
                        << parameters.size() 
                        << " and number of parameters "
                        <<  this->GetNumberOfParameters() );
    }
    
    m_InternalParametersBuffer = parameters;
    m_Velocity->SetParameters( m_InternalParametersBuffer );
    
    // wrap flat array as images of coefficients
    m_Velocity->WrapAsImages();
    
    // Modified is always called since we just have a pointer to the
    // parameters and cannot know if the parameters have changed.
    this->Modified();
    
  }
  
  
  // Set the Fixed Parameters
  template<class TScalarType, unsigned int NDimensions, unsigned int VSplineOrder>
  void
  DiffeomorphicContinuousBSplineTransform<TScalarType, NDimensions, VSplineOrder>
  ::SetFixedParameters( const ParametersType & passedParameters )
  {
    
    ParametersType parameters( SpacePlusTimeDimension * (3 + SpacePlusTimeDimension) );
    
    // check if the number of parameters match the
    // expected number of parameters
    if ( passedParameters.Size() == SpacePlusTimeDimension * 3 )
    {
      parameters.Fill( 0.0 );
      for(unsigned int i=0; i<3 * SpacePlusTimeDimension; i++)
      {
        parameters[i] = passedParameters[i];
      }
      for (unsigned int di=0; di<SpacePlusTimeDimension; di++)
      {
        parameters[3*SpacePlusTimeDimension+(di*SpacePlusTimeDimension+di)] = 1;
      }
    }
    else if ( passedParameters.Size() != SpacePlusTimeDimension * (3 + SpacePlusTimeDimension) )
    {
      itkExceptionMacro(<< "Mismatched between parameters size "
                        << passedParameters.size() 
                        << " and number of fixed parameters "
                        << SpacePlusTimeDimension * (3 + SpacePlusTimeDimension) );
    }
    else
    {
      for(unsigned int i=0; i<SpacePlusTimeDimension * (3 + SpacePlusTimeDimension); i++)
      {
        parameters[i] = passedParameters[i];
      }
    }
    
    /********************************************************* 
     Fixed Parameters store the following information:
     Grid Size
     Grid Origin
     Grid Spacing
     Grid Direction
     The size of these is equal to the  NInputDimensions
     *********************************************************/
    
    /** Set the Grid Parameters */
    SizeType   gridSize;
    for (unsigned int i=0; i<SpacePlusTimeDimension; i++)
    {
      gridSize[i] = static_cast<int> (parameters[i]);
    }
    RegionType bsplineRegion;
    bsplineRegion.SetSize( gridSize );
    
    /** Set the Origin Parameters */
    OriginType origin;
    for (unsigned int i=0; i<SpacePlusTimeDimension; i++)
    {
      origin[i] = parameters[SpacePlusTimeDimension+i];
    }
    
    /** Set the Spacing Parameters */
    SpacingType spacing;
    for (unsigned int i=0; i<SpacePlusTimeDimension; i++)
    {
      spacing[i] = parameters[2*SpacePlusTimeDimension+i];
    }
    
    /** Set the Direction Parameters */
    DirectionType direction;
    for (unsigned int di=0; di<SpacePlusTimeDimension; di++)
    {
      for (unsigned int dj=0; dj<SpacePlusTimeDimension; dj++)
      {
        direction[di][dj] = parameters[3*SpacePlusTimeDimension+(di*SpacePlusTimeDimension+dj)];
      }
    }
    
    
    m_Velocity->SetGridSpacing( spacing );
    m_Velocity->SetGridDirection( direction );
    m_Velocity->SetGridOrigin( origin );
    m_Velocity->SetGridRegion( bsplineRegion );
    
    this->Modified();
  }
  
  
  // Get the parameters
  template<class TScalarType, unsigned int NDimensions, unsigned int VSplineOrder>
  const 
  typename DiffeomorphicContinuousBSplineTransform<TScalarType, NDimensions, VSplineOrder>
  ::ParametersType &
  DiffeomorphicContinuousBSplineTransform<TScalarType, NDimensions, VSplineOrder>
  ::GetParameters( void ) const
  {
    return m_Velocity->GetParameters();
  }
  
  
  // Get the fixed parameters
  template<class TScalarType, unsigned int NDimensions, unsigned int VSplineOrder>
  const 
  typename DiffeomorphicContinuousBSplineTransform<TScalarType, NDimensions, VSplineOrder>
  ::ParametersType &
  DiffeomorphicContinuousBSplineTransform<TScalarType, NDimensions, VSplineOrder>
  ::GetFixedParameters( void ) const
  {
    RegionType resRegion = m_Velocity->GetGridRegion(  );
    
    for (unsigned int i=0; i<SpacePlusTimeDimension; i++)
    {
      this->m_FixedParameters[i] = (resRegion.GetSize())[i];
    }
    for (unsigned int i=0; i<SpacePlusTimeDimension; i++)
    {
      this->m_FixedParameters[SpacePlusTimeDimension+i] = (m_Velocity->GetGridOrigin())[i];
    } 
    for (unsigned int i=0; i<SpacePlusTimeDimension; i++)
    {
      this->m_FixedParameters[2*SpacePlusTimeDimension+i] =  (m_Velocity->GetGridSpacing())[i];
    }
    for (unsigned int di=0; di<SpacePlusTimeDimension; di++)
    {
      for (unsigned int dj=0; dj<SpacePlusTimeDimension; dj++)
      {
        this->m_FixedParameters[3*SpacePlusTimeDimension+(di*SpacePlusTimeDimension+dj)] = (m_Velocity->GetGridDirection())[di][dj];
      }
    }
    
    return (this->m_FixedParameters);
  }
  
  
  // Print self
  template<class TScalarType, unsigned int NDimensions, unsigned int VSplineOrder>
  void
  DiffeomorphicContinuousBSplineTransform<TScalarType, NDimensions, VSplineOrder>
  ::PrintSelf(std::ostream &os, Indent indent) const
  {
    this->Superclass::PrintSelf(os, indent);
  }
  
  // Transform a point
  template<class TScalarType, unsigned int NDimensions, unsigned int VSplineOrder>
  typename DiffeomorphicContinuousBSplineTransform<TScalarType, NDimensions, VSplineOrder>
  ::OutputPointType
  DiffeomorphicContinuousBSplineTransform<TScalarType, NDimensions, VSplineOrder>
  ::TransformPoint(const InputPointType &point, const TimePointType &endTimePoint) const 
  {
    typename VelocityType::WeightsType weights( m_Velocity->GetWeightsFunction()->GetNumberOfWeights() );
    typename VelocityType::ParameterIndexArrayType indices( m_Velocity->GetWeightsFunction()->GetNumberOfWeights() );
    OutputPointType outputPoint;
    bool inside;
    unsigned int j, numberOfTimeSteps;
    InputPointType spaceTimePoint;
    SpacePointType bulkPoint,velocity; 
    double firstTimeStep,lastTimeStep;
    
    
    // Check for bulk transform
    if ( m_BulkTransform )
    {
      for ( j = 0; j < SpaceDimension; j++ )
      {    
        bulkPoint[j] = point[j];    
      }
      bulkPoint = m_BulkTransform->TransformPoint( bulkPoint );
    }
    
    spaceTimePoint=point;
    firstTimeStep=point[SpaceDimension];
    
    numberOfTimeSteps = (unsigned int) floor( vnl_math_abs((endTimePoint-firstTimeStep)/m_TimeStep));
    double timeStepWithSign = m_TimeStep;
    if (endTimePoint < firstTimeStep)
    {
      timeStepWithSign *= -1.;
    }
    
    if ( this->GetCoefficientImage()[0] )
    {
      for (unsigned int i= 0; i<numberOfTimeSteps; i++)
      {
        m_Velocity->GetPointVelocity( spaceTimePoint, velocity, weights, indices, inside);
        
        for ( j = 0; j < SpaceDimension; j++ )
        {
          spaceTimePoint[j] += velocity[j] * timeStepWithSign;
        }
        
        // Update the time point
        spaceTimePoint[SpaceDimension] +=timeStepWithSign;
      }
      
      // We need to cover the time interval between the current step and endTimePoint
      if (  ( (timeStepWithSign>=0)&&(spaceTimePoint[SpaceDimension] < endTimePoint) )
          ||( (timeStepWithSign<0)&&(spaceTimePoint[SpaceDimension] > endTimePoint) ) )
      {
        lastTimeStep= endTimePoint-spaceTimePoint[SpaceDimension];
        
        m_Velocity->GetPointVelocity( spaceTimePoint, velocity, weights, indices, inside);
        
        for ( j = 0; j < SpaceDimension; j++ )
        {
          spaceTimePoint[j] += velocity[j]*lastTimeStep;
        }
        
        // The last time point is the endTimePoint
        spaceTimePoint[SpaceDimension]=endTimePoint;
      }
      
      
      outputPoint = spaceTimePoint;
      
      // Bulk transform point  
      if ( m_BulkTransform )
      {
        for (unsigned int d=0; d<SpaceDimension; d++)
        {
          outputPoint[d] = spaceTimePoint[d] - point[d] + bulkPoint[d];
        }
      }
    }
    else
    {
      
      itkWarningMacro( << "B-spline coefficients have not been set" );
      outputPoint=point;
    }
    
    return outputPoint;
  }
  
  
  template<class TScalarType, unsigned int NDimensions, unsigned int VSplineOrder>
  void
  DiffeomorphicContinuousBSplineTransform<TScalarType, NDimensions, VSplineOrder>
  ::GetIncrementalSparseJacobian(
                                 const InputPointType & point, 
                                 OutputPointType & outPoint,
                                 const TimePointType & endTimePoint, 
                                 SparseJacobianType & jaco, 
                                 bool & valid ) const
  {
    InputPointType spaceTimePoint;
    
    bool inside, endTimePointReached;
    unsigned int i,j,k;
    unsigned long pos;
    double determinant, currentTimeStep;
        
    typename SparseJacobianType::iterator it;
    double* blockWeights;
    long unsigned int* blockIndices;
    
    // Initial time and number of time steps
    spaceTimePoint=point;
    double firstTimeStep=(double)point[SpaceDimension];
    
    valid = true;
    endTimePointReached = false;
    currentTimeStep = m_TimeStep;
    if (endTimePoint < firstTimeStep)
    {
      currentTimeStep *= -1.;
    }
    
    while ( ! endTimePointReached )
    {
      
      
      // Set the current time step: In case we are at the last integration step (i.e. the next 
      // time point would exceed the final time), we reduce the time step accordingly and set
      // endTimePointReached to true 
      if ( ( (currentTimeStep>0)&&(spaceTimePoint[SpaceDimension] + currentTimeStep >= endTimePoint) )
          || ( (currentTimeStep<=0)&&(spaceTimePoint[SpaceDimension] + currentTimeStep <= endTimePoint) ) )
      {
        currentTimeStep = endTimePoint - spaceTimePoint[SpaceDimension];
        endTimePointReached = true;
      }
      
      m_Velocity->GetPointVelocity(spaceTimePoint, incrJaco_velocity, incrJaco_weights, incrJaco_indices, inside);

	  // If the point is not inside the valid region, the transform is the identity and 
      // jaco = jaco
      if (inside)
      {
        // Compute physicaljacobian(transform)=dT(x)/dx = Id+physicaljacobian(velocity)*timeStep
        inside = m_Velocity->GetClassicalJacobian(spaceTimePoint, incrJaco_physJacobian, currentTimeStep);
        
        determinant = static_cast<PixelType>(vnl_determinant<double>(incrJaco_physJacobian));
        
        if ( determinant <= 0. )
        {
          //valid = false;
          outPoint=point;
          return;
        }
        
        // We actually want to compute this :
        //    jaco = physJacobian*jaco + m_TimeStep * partialJacobian;
        // The first term is obtained by taking all non-zero column of the jacobian, 
        //  and replacing them by a linear combination of the columns of the phys. jacobian,
        //  the coefficients of the linear combination appearing in the currently considered
        //  jacobian column.
        
        it = jaco.begin();
        
        while ( it != jaco.end() )
        {
          
          for (j=0; j<SpaceDimension; j++)
          {
            // reset the result
            incrJaco_columnOut[j] = 0.;
            // compute the result         
            for (k=0; k<SpaceDimension; k++)
            {
              incrJaco_columnOut[j] += it->second[k] * incrJaco_physJacobian[j][k];
            }
          }
           
          it->second = incrJaco_columnOut;
          
          ++it;
        }
        
        
        // The second term is obtained by checking if the column contains non-zero elements
        // If yes, add the bspline interpolation weight on the correct line
        // If not, create a new column.
        
        blockWeights = incrJaco_weights.data_block();
        blockIndices = incrJaco_indices.data_block();
        
        for (i=0; i < incrJaco_weights.Size(); i++, ++blockIndices, ++blockWeights)
        {
          for (j=0; j<SpaceDimension; j++)
          {
            pos = (*blockIndices) + j * this->GetNumberOfParametersPerDimension();
            it = jaco.find(pos);
            
            if ( it != jaco.end() )
            {
              // This column was already existing in the map
              it->second[j] += currentTimeStep * (*blockWeights);//tmpWeight;
            }
            else
            {
              // The element was not existing in the map so we create it and
              // insert it in the map.
              incrJaco_columnOut[j] = currentTimeStep * (*blockWeights);//tmpWeight;
              jaco[pos] = incrJaco_columnOut;
            }
          } 
        }
        
        
        //Update point
        for ( j = 0; j < SpaceDimension; j++ )
        {
          spaceTimePoint[j]+=incrJaco_velocity[j]*currentTimeStep;
        }
      } // end of inside
      
      // Update time step
      spaceTimePoint[SpaceDimension] += currentTimeStep;
      
      
    } // while ( ! endTimePointReached )
    
    outPoint = spaceTimePoint;
    
    return;
    
  }
  
  template<class TScalarType, unsigned int NDimensions, unsigned int VSplineOrder>
  const 
  typename DiffeomorphicContinuousBSplineTransform<TScalarType, NDimensions, VSplineOrder>
  ::JacobianType & 
  DiffeomorphicContinuousBSplineTransform<TScalarType, NDimensions, VSplineOrder>
  ::GetJacobian( const InputPointType & point, const TimePointType & endTimePoint ) const
  {
    
    // In case the user still wants to access the jacobian in the traditional way, we reset
    // the full jacobian matrix and copy all info from the sparse jacobian to the full jacobian.
    
    bool valid = true; 
    m_TotalJacobian.Fill( NumericTraits<PixelType>::Zero );
    
    SparseJacobianType sparseJaco;
    OutputPointType outPoint;
    
    this->GetIncrementalSparseJacobian(point, outPoint, endTimePoint, sparseJaco, valid );
    
    while (!valid && m_TimeStep >= m_MinimumTimeStep)
    {
      m_TimeStep *=1/2;
      this->GetIncrementalSparseJacobian(point, outPoint, endTimePoint, sparseJaco, valid);
      itkDebugMacro ( << "Dividing time step by 2. Current time step is " << m_TimeStep << std::endl );
    }
    
    if (!valid)
      itkExceptionMacro( << "Determinant of physical jacobian is negative! ");
    
    // Copy back sparse jacobian into full jacobian
    typename SparseJacobianType::iterator it = sparseJaco.begin();
    
    while ( it != sparseJaco.end() )
    {
      unsigned int column = (*it).first;
      for (unsigned int d=0; d<SpaceDimension; d++)
      {
        m_TotalJacobian[d][column] = (*it).second[d];
      }
      ++it;
    }
    
    return m_TotalJacobian;
  }
    
  template<class TScalarType, unsigned int NDimensions, unsigned int VSplineOrder>
  void 
  DiffeomorphicContinuousBSplineTransform<TScalarType, NDimensions, VSplineOrder>
  ::TransformPointAndGetPhysicalJacobian( const InputPointType & point, const TimePointType & endTimePoint, ClassicalJacobianType & output, OutputPointType& outputPoint) const
  {
      // >>>>>>>
      ClassicalJacobianType tmpOutput(SpaceDimension,SpaceDimension);
      output.set_size(SpaceDimension,SpaceDimension);
      output.set_identity();
      // <<<<<<<
      
      // FIXME: Code refactoring for avoiding duplications
      typename VelocityType::WeightsType weights( m_Velocity->GetWeightsFunction()->GetNumberOfWeights() );
      typename VelocityType::ParameterIndexArrayType indices( m_Velocity->GetWeightsFunction()->GetNumberOfWeights() );
      bool inside;
      unsigned int j, numberOfTimeSteps;
      InputPointType spaceTimePoint;
      SpacePointType bulkPoint,velocity; 
      double firstTimeStep,lastTimeStep;
      
      // Check for bulk transform
      if ( m_BulkTransform )
      {
          for ( j = 0; j < SpaceDimension; j++ )
          {    
              bulkPoint[j] = point[j];    
          }
          bulkPoint = m_BulkTransform->TransformPoint( bulkPoint );
      }
      
      spaceTimePoint=point;
      firstTimeStep=point[SpaceDimension];
      
      numberOfTimeSteps = (unsigned int) floor( vnl_math_abs((endTimePoint-firstTimeStep)/m_TimeStep));
      double timeStepWithSign = m_TimeStep;
      if (endTimePoint < firstTimeStep)
      {
          timeStepWithSign *= -1.;
      }
      
      if ( this->GetCoefficientImage()[0] )
      {
          for (unsigned int i= 0; i<numberOfTimeSteps; i++)
          {
              m_Velocity->GetPointVelocity( spaceTimePoint, velocity, weights, indices, inside);
              // >>>>>>>
              m_Velocity->GetClassicalJacobian (spaceTimePoint, tmpOutput, m_TimeStep);
              output = tmpOutput*output;
              // <<<<<<<
              
              for ( j = 0; j < SpaceDimension; j++ )
              {
                  spaceTimePoint[j] += velocity[j] * timeStepWithSign;
              }
              
              // Update the time point
              spaceTimePoint[SpaceDimension] +=timeStepWithSign;
          }
          
          // We need to cover the time interval between the current step and endTimePoint
          if (  ( (timeStepWithSign>=0)&&(spaceTimePoint[SpaceDimension] < endTimePoint) )
              ||( (timeStepWithSign<0)&&(spaceTimePoint[SpaceDimension] > endTimePoint) ) )
          {
              lastTimeStep= endTimePoint-spaceTimePoint[SpaceDimension];
              
              m_Velocity->GetPointVelocity( spaceTimePoint, velocity, weights, indices, inside);
              // >>>>>>>
              m_Velocity->GetClassicalJacobian (spaceTimePoint, tmpOutput, m_TimeStep);
              output = tmpOutput*output;
              // <<<<<<<
              
              for ( j = 0; j < SpaceDimension; j++ )
              {
                  spaceTimePoint[j] += velocity[j]*lastTimeStep;
              }
              
              // The last time point is the endTimePoint
              spaceTimePoint[SpaceDimension]=endTimePoint;
          }
          
          outputPoint = spaceTimePoint;
          
          // Bulk transform point  
          if ( m_BulkTransform )
          {
              for (unsigned int d=0; d<SpaceDimension; d++)
              {
                  outputPoint[d] = spaceTimePoint[d] - point[d] + bulkPoint[d];
              }
          }
      }
      else
      {
          
          itkWarningMacro( << "B-spline coefficients have not been set" );
          outputPoint=point;
      }
      
  }    

    
  
} // namespace

#endif
