#ifndef _regBSplineFieldRefine_txx
#define _regBSplineFieldRefine_txx

// First make sure that the configuration is available.
// This line can be removed once the optimized versions
// gets integrated into the main directories.
#include "itkConfigure.h"
#include "itkTimeProbe.h"
#include "regBSplineFieldRefine.h"
#include "itkIdentityTransform.h"
#include "itkImageRegionIterator.h"

namespace itk
{

/*
 * Constructor
 */
template < class TScalarType, unsigned int NDimensions, unsigned int VSplineOrder>
BSplineFieldRefine< TScalarType, NDimensions, VSplineOrder>
::BSplineFieldRefine()
{
  m_BSplineOutputField = BSplineFieldType::New();
  m_BSplineInputField = 0; 
  m_OutputParameters = ParametersType(0);
}

/*
 * Update
 */
template < class TScalarType, unsigned int NDimensions, unsigned int VSplineOrder>
void
BSplineFieldRefine< TScalarType, NDimensions, VSplineOrder>
::Update(void)
{ 
  // Check that input has been passed
  if ( m_BSplineInputField.IsNull() )
  {
    itkExceptionMacro ( << "No input has been given" << std::endl);
  }	
  
  typedef itk::IdentityTransform<double, itk::GetImageDimension<ImageType>::ImageDimension> IdentityTransformType;
  typename IdentityTransformType::Pointer identity = IdentityTransformType::New();
  
  typename InterpolatorType::Pointer interpol = InterpolatorType::New();
    
  const unsigned int spaceDimension = BSplineFieldType::SpaceDimension;
    
  for (unsigned int d=0; d<spaceDimension; d++)
  {
    typename ImageType::Pointer coeff = m_BSplineInputField->GetCoefficientImage()[d];
    
    ResamplerPointerType upsampler = ResamplerType::New();
    
    upsampler->SetInterpolator( interpol );
    upsampler->SetInput( coeff );
    upsampler->SetTransform( identity );
    upsampler->SetSize( m_BSplineOutputField->GetGridRegion().GetSize() );
    upsampler->SetOutputSpacing( m_BSplineOutputField->GetGridSpacing() );
    upsampler->SetOutputOrigin( m_BSplineOutputField->GetGridOrigin() );

    DecompositionPointerType decomposition = DecompositionType::New();

    decomposition->SetSplineOrder( VSplineOrder );
    decomposition->SetInput( upsampler->GetOutput() );
    decomposition->Update();
    
    typedef ImageRegionIterator<ImageType> IteratorType;
    RegionType regionIn = decomposition->GetOutput()->GetBufferedRegion();
    RegionType regionOut = m_BSplineOutputField->GetCoefficientImage()[d]->GetBufferedRegion();
    
    if ( regionIn.GetNumberOfPixels() != regionOut.GetNumberOfPixels())
    {
      itkExceptionMacro( << "Input and output regions do not match. " << std::endl);
    }
    
    IteratorType it( decomposition->GetOutput(), regionIn );
	IteratorType itOut ( m_BSplineOutputField->GetCoefficientImage()[d], regionOut);
	
	for (it.GoToBegin(),itOut.GoToBegin(); !(it.IsAtEnd()||itOut.IsAtEnd()); ++it,++itOut) 
	{
	  itOut.Set(it.Get());
	}
  }
}


/*
 * Configure Output
 */
template < class TScalarType, unsigned int NDimensions, unsigned int VSplineOrder>
void
BSplineFieldRefine< TScalarType, NDimensions, VSplineOrder>
::ConfigureOutput(const SpacingType& spacing, const OriginType& origin, const RegionType& region)
{
   m_BSplineOutputField->SetGridSpacing(spacing);
   m_BSplineOutputField->SetGridOrigin(origin);
   m_BSplineOutputField->SetGridRegion(region);
   const unsigned int outputNumberOfParameters = m_BSplineOutputField->GetNumberOfParameters();
   m_OutputParameters.SetSize(outputNumberOfParameters);
   m_BSplineOutputField->SetParameters(m_OutputParameters);
   m_BSplineOutputField->WrapAsImages();
}

/*
 * SetInput
 */
template < class TScalarType, unsigned int NDimensions, unsigned int VSplineOrder>
void
BSplineFieldRefine< TScalarType, NDimensions, VSplineOrder>
::SetInput(BSplineFieldType* field)
{
   m_BSplineInputField = field;
}

/*
 * GetOutput
 */
template < class TScalarType, unsigned int NDimensions, unsigned int VSplineOrder>
typename BSplineFieldRefine< TScalarType, NDimensions, VSplineOrder>::BSplineFieldType*
BSplineFieldRefine< TScalarType, NDimensions, VSplineOrder>
::GetOutput(void)
{
   return m_BSplineOutputField;
}

} // typename

#endif

