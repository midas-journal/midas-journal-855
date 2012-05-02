/*=========================================================================
 
 Program:   Insight Segmentation & Registration Toolkit
 Module:    $RCSfile: regContinuousBSplineTransform.txx,v $
 Language:  C++
 Date:      $Date: 2008-05-08 23:22:35 $
 Version:   $Revision: 1.41 $
 
 Copyright (c) Insight Software Consortium. All rights reserved.
 See ITKCopyright.txt or http://www.itk.org/HTML/Copyright.htm for details.
 
 This software is distributed WITHOUT ANY WARRANTY; without even 
 the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR 
 PURPOSE.  See the above copyright notices for more information.
 
 =========================================================================*/
#ifndef __regContinuousBSplineTransform_txx
#define __regContinuousBSplineTransform_txx

#include "regContinuousBSplineTransform.h"
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
  ContinuousBSplineTransform<TScalarType, NDimension, VSplineOrder>
  ::ContinuousBSplineTransform(void):Superclass(SpacePlusTimeDimension,0)
  {
    
    // BSplineField object
    m_BSplineField = BSplineFieldType::New();
    
    // Instantiate an identity transform
    typedef IdentityTransform<ScalarType, SpaceDimension> IdentityTransformType;
    typename IdentityTransformType::Pointer id = IdentityTransformType::New();
    m_BulkTransform = id;
    
    m_InternalParametersBuffer = ParametersType(0);
    // Make sure the parameters pointer is not NULL after construction.
    m_BSplineField->SetParameters( m_InternalParametersBuffer );
    
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
      this->m_FixedParameters[2*SpacePlusTimeDimension+i] = m_BSplineField->GetGridSpacing()[i];
    }
    
    for (unsigned int di=0; di<SpacePlusTimeDimension; di++)
    {
      for (unsigned int dj=0; dj<SpacePlusTimeDimension; dj++)
      {
        this->m_FixedParameters[3*SpacePlusTimeDimension+(di*SpacePlusTimeDimension+dj)] = m_BSplineField->GetGridDirection()[di][dj];
      }
    }
    
    // Cache variables allocation
    incrJaco_columnOut = SparseJacobianColumnType(SpaceDimension, 0.);
    incrJaco_weights = typename BSplineFieldType::WeightsType( m_BSplineField->GetWeightsFunction()->GetNumberOfWeights() );
    incrJaco_indices = typename BSplineFieldType::ParameterIndexArrayType( m_BSplineField->GetWeightsFunction()->GetNumberOfWeights() );
    
  }
  
  // Destructor
  template<class TScalarType, unsigned int NDimensions, unsigned int VSplineOrder>
  ContinuousBSplineTransform<TScalarType, NDimensions, VSplineOrder>
  ::~ContinuousBSplineTransform()
  {
  }
  
  
  // Name of the transform.
  template<class TScalarType, unsigned int NDimensions, unsigned int VSplineOrder>
  std::string 
  ContinuousBSplineTransform<TScalarType, NDimensions, VSplineOrder>
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
  ContinuousBSplineTransform<TScalarType, NDimensions, VSplineOrder>
  ::GetNumberOfParameters(void) const
  {
    return m_BSplineField->GetNumberOfParameters();
  }
  
  
  // Get the number of parameters per dimension
  template<class TScalarType, unsigned int NDimensions, unsigned int VSplineOrder>
  unsigned int
  ContinuousBSplineTransform<TScalarType, NDimensions, VSplineOrder>
  ::GetNumberOfParametersPerDimension(void) const
  {
    return m_BSplineField->GetNumberOfParametersPerDimension();
  }
  
  
  // Set the grid region
  template<class TScalarType, unsigned int NDimensions, unsigned int VSplineOrder>
  void
  ContinuousBSplineTransform<TScalarType, NDimensions, VSplineOrder>
  ::SetGridRegion( const RegionType& region )
  {
    m_BSplineField->SetGridRegion(region);
    
    // Input parameters point to internal buffer => using default parameters.
    if ( & m_BSplineField->GetParameters() == & m_InternalParametersBuffer)
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
        
    this->Modified();
  }
  
  
  // Set the grid spacing
  template<class TScalarType, unsigned int NDimensions, unsigned int VSplineOrder>
  void
  ContinuousBSplineTransform<TScalarType, NDimensions, VSplineOrder>
  ::SetGridSpacing( const SpacingType& spacing )
  {
    m_BSplineField->SetGridSpacing(spacing);
  }
  
  // Set the grid direction
  template<class TScalarType, unsigned int NDimensions, unsigned int VSplineOrder>
  void
  ContinuousBSplineTransform<TScalarType, NDimensions, VSplineOrder>
  ::SetGridDirection( const DirectionType & direction )
  {
    m_BSplineField->SetGridDirection( direction );
  }
  
  
  // Set the grid origin
  template<class TScalarType, unsigned int NDimensions, unsigned int VSplineOrder>
  void
  ContinuousBSplineTransform<TScalarType, NDimensions, VSplineOrder>
  ::SetGridOrigin( const OriginType& origin )
  {
    m_BSplineField->SetGridOrigin( origin );
  }
  
  // Get the grid origin
  template<class TScalarType, unsigned int NDimensions, unsigned int VSplineOrder>
  typename ContinuousBSplineTransform<TScalarType, NDimensions, VSplineOrder>::OriginType
  ContinuousBSplineTransform<TScalarType, NDimensions, VSplineOrder>
  ::GetGridOrigin()
  {
    return m_BSplineField->GetGridOrigin();
  }
  
  // Get the grid spacing
  template<class TScalarType, unsigned int NDimensions, unsigned int VSplineOrder>
  typename ContinuousBSplineTransform<TScalarType, NDimensions, VSplineOrder>::SpacingType
  ContinuousBSplineTransform<TScalarType, NDimensions, VSplineOrder>
  ::GetGridSpacing()
  {
    return m_BSplineField->GetGridSpacing();
  }
  
  // Get the grid region
  template<class TScalarType, unsigned int NDimensions, unsigned int VSplineOrder>
  typename ContinuousBSplineTransform<TScalarType, NDimensions, VSplineOrder>::RegionType
  ContinuousBSplineTransform<TScalarType, NDimensions, VSplineOrder>
  ::GetGridRegion()
  {
    return m_BSplineField->GetGridRegion();
  }
  
  // Get Coefficient Image
  template<class TScalarType, unsigned int NDimensions, unsigned int VSplineOrder>
  const 
  typename ContinuousBSplineTransform<TScalarType, NDimensions, VSplineOrder>::ImagePointer * 
  ContinuousBSplineTransform<TScalarType, NDimensions, VSplineOrder>
  ::GetCoefficientImage() const
  {
    return m_BSplineField->GetCoefficientImage();
  }
  
  
  
  // Reset transform to identity
  template<class TScalarType, unsigned int NDimensions, unsigned int VSplineOrder>
  void
  ContinuousBSplineTransform<TScalarType, NDimensions, VSplineOrder>
  ::SetIdentity()
  {
    // Check if we need to resize the default parameter buffer.
    if ( m_InternalParametersBuffer.GetSize() != this->GetNumberOfParameters() )
    {
      m_InternalParametersBuffer.SetSize( this->GetNumberOfParameters() );
    }
    
    // Fill with zeros for identity.
    m_InternalParametersBuffer.Fill( 0 );
    
    // Pass parameters to BSplineField
    m_BSplineField->SetParameters( &m_InternalParametersBuffer );
    
  }
  
  // Set the parameters
  template<class TScalarType, unsigned int NDimensions, unsigned int VSplineOrder>
  void
  ContinuousBSplineTransform<TScalarType, NDimensions, VSplineOrder>
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
    m_BSplineField->SetParameters( parameters );
    
    // Wrap flat array as images of coefficients
    m_BSplineField->WrapAsImages();
    
    // Modified is always called since we just have a pointer to the
    // parameters and cannot know if the parameters have changed.
    this->Modified();
  }
  
  
  // Set the parameters
  template<class TScalarType, unsigned int NDimensions, unsigned int VSplineOrder>
  void
  ContinuousBSplineTransform<TScalarType, NDimensions, VSplineOrder>
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
    m_BSplineField->SetParameters( m_InternalParametersBuffer );
    
    // wrap flat array as images of coefficients
    m_BSplineField->WrapAsImages();
    
    // Modified is always called since we just have a pointer to the
    // parameters and cannot know if the parameters have changed.
    this->Modified();
    
  }
  
  
  // Set the Fixed Parameters
  template<class TScalarType, unsigned int NDimensions, unsigned int VSplineOrder>
  void
  ContinuousBSplineTransform<TScalarType, NDimensions, VSplineOrder>
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
    
    
    m_BSplineField->SetGridSpacing( spacing );
    m_BSplineField->SetGridDirection( direction );
    m_BSplineField->SetGridOrigin( origin );
    m_BSplineField->SetGridRegion( bsplineRegion );
    
    this->Modified();
  }
  
  
  // Get the parameters
  template<class TScalarType, unsigned int NDimensions, unsigned int VSplineOrder>
  const 
  typename ContinuousBSplineTransform<TScalarType, NDimensions, VSplineOrder>
  ::ParametersType &
  ContinuousBSplineTransform<TScalarType, NDimensions, VSplineOrder>
  ::GetParameters( void ) const
  {
    return m_BSplineField->GetParameters();
  }
  
  
  // Get the fixed parameters
  template<class TScalarType, unsigned int NDimensions, unsigned int VSplineOrder>
  const 
  typename ContinuousBSplineTransform<TScalarType, NDimensions, VSplineOrder>
  ::ParametersType &
  ContinuousBSplineTransform<TScalarType, NDimensions, VSplineOrder>
  ::GetFixedParameters( void ) const
  {
    RegionType resRegion = m_BSplineField->GetGridRegion(  );
    
    for (unsigned int i=0; i<SpacePlusTimeDimension; i++)
    {
      this->m_FixedParameters[i] = (resRegion.GetSize())[i];
    }
    for (unsigned int i=0; i<SpacePlusTimeDimension; i++)
    {
      this->m_FixedParameters[SpacePlusTimeDimension+i] = (m_BSplineField->GetGridOrigin())[i];
    } 
    for (unsigned int i=0; i<SpacePlusTimeDimension; i++)
    {
      this->m_FixedParameters[2*SpacePlusTimeDimension+i] =  (m_BSplineField->GetGridSpacing())[i];
    }
    for (unsigned int di=0; di<SpacePlusTimeDimension; di++)
    {
      for (unsigned int dj=0; dj<SpacePlusTimeDimension; dj++)
      {
        this->m_FixedParameters[3*SpacePlusTimeDimension+(di*SpacePlusTimeDimension+dj)] = (m_BSplineField->GetGridDirection())[di][dj];
      }
    }
    
    return (this->m_FixedParameters);
  }
  
  
  // Print self
  template<class TScalarType, unsigned int NDimensions, unsigned int VSplineOrder>
  void
  ContinuousBSplineTransform<TScalarType, NDimensions, VSplineOrder>
  ::PrintSelf(std::ostream &os, Indent indent) const
  {
    this->Superclass::PrintSelf(os, indent);
  }
  
  // Transform a point
  template<class TScalarType, unsigned int NDimensions, unsigned int VSplineOrder>
  typename ContinuousBSplineTransform<TScalarType, NDimensions, VSplineOrder>
  ::OutputPointType
  ContinuousBSplineTransform<TScalarType, NDimensions, VSplineOrder>
  ::TransformPoint(const InputPointType &point, const TimePointType &endTimePoint) const 
  {
    typename BSplineFieldType::WeightsType weights( m_BSplineField->GetWeightsFunction()->GetNumberOfWeights() );
    typename BSplineFieldType::ParameterIndexArrayType indices( m_BSplineField->GetWeightsFunction()->GetNumberOfWeights() );
    OutputPointType outputPoint;
    bool inside;
    unsigned int j;
    
    SpacePointType bulkPoint,outdisp; 
    InputPointType point2 = point;
    point2[SpaceDimension] = endTimePoint;
     
    // Check for bulk transform
    if ( m_BulkTransform )
    {
      for ( j = 0; j < SpaceDimension; j++ )
      {    
        bulkPoint[j] = point[j];    
      }
      bulkPoint = m_BulkTransform->TransformPoint( bulkPoint );
    }
    
    m_BSplineField->GetPointVelocity( point2, outdisp, weights, indices, inside);
    
    OutputPointType output;
    for (unsigned int d=0; d<SpaceDimension; d++) {
      output[d] = point[d]+outdisp[d];
    }
    output[SpaceDimension] = endTimePoint;
    
    return output;
  }
  
 
  template<class TScalarType, unsigned int NDimensions, unsigned int VSplineOrder>
  const 
  typename ContinuousBSplineTransform<TScalarType, NDimensions, VSplineOrder>
  ::JacobianType & 
  ContinuousBSplineTransform<TScalarType, NDimensions, VSplineOrder>
  ::GetJacobian( const InputPointType & point, const TimePointType & endTimePoint ) const
  {
    
    // In case the user still wants to access the jacobian in the traditional way, we reset
    // the full jacobian matrix and copy all info from the sparse jacobian to the full jacobian.
    this->m_Jacobian.Fill( NumericTraits<PixelType>::Zero );
    
    SparseJacobianType sparseJaco;
    this->GetSparseJacobian(point, endTimePoint, sparseJaco);
    
    // Copy back sparse jacobian into full jacobian
    typename SparseJacobianType::const_iterator it = sparseJaco.begin();
    
    while ( it != sparseJaco.end() )
    {
      unsigned int column = (*it).first;
      for (unsigned int d=0; d<SpaceDimension; d++)
      {
        this->m_Jacobian[d][column] = (*it).second[d];
      }
      ++it;
    }
    
    return this->m_Jacobian;
  }
  
  
  template<class TScalarType, unsigned int NDimensions, unsigned int VSplineOrder>
  void
  ContinuousBSplineTransform<TScalarType, NDimensions, VSplineOrder>
  ::GetSparseJacobian( const InputPointType & point, const TimePointType & endTimePoint, SparseJacobianType& jacobian ) const
  {
    
    // In case the user still wants to access the jacobian in the traditional way, we reset
    // the full jacobian matrix and copy all info from the sparse jacobian to the full jacobian.
    
    bool valid = true; 
    jacobian.clear();
    
    InputPointType point2 = point;
    point2[SpaceDimension]=endTimePoint;
    
    m_BSplineField->GetPointVelocity(point2, incrJaco_velocity, incrJaco_weights, incrJaco_indices, valid);
    
    double* blockWeights = incrJaco_weights.data_block();
    long unsigned int* blockIndices = incrJaco_indices.data_block();
    unsigned int i,j;
    long unsigned int pos;    
    typename SparseJacobianType::iterator it;
    
    for (i=0; i < incrJaco_weights.Size(); i++, ++blockIndices, ++blockWeights)
    {
      for (j=0; j<SpaceDimension; j++)
      {
        
          pos = (*blockIndices) + j * this->GetNumberOfParametersPerDimension();
          it = jacobian.find(pos);
          
          if ( it != jacobian.end() )
          {
            // This column was already existing in the map
            it->second[j] += (*blockWeights);
          }
          else
          {
            // The element was not existing in the map so we create it and
            // insert it in the map.
            incrJaco_columnOut[j] = (*blockWeights);
            jacobian[pos] = incrJaco_columnOut;
          }
        
      } 
      
    }
    
  }

  
} // namespace

#endif
