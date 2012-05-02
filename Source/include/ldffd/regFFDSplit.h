/*
** Splitting BSpline transformation
**
*/

#ifndef   	REGFFDSPLIT_H_
# define   	REGFFDSPLIT_H_

#include "itkImageRegion.h"
#include "itkImage.h"
#include "itkBSplineDeformableTransform.h"
#include "itkBSplineResampleImageFunction.h"


template <class TransformType, class ImageType>
  typename TransformType::Pointer
  RefineFFDTransformation( TransformType* transfo, unsigned int SplineOrder,
							typename TransformType::RegionType bsplineRegion, 
							typename TransformType::OriginType originHigh,
							typename TransformType::SpacingType spacingHigh)
{
  const unsigned int SpaceDimension = transfo->GetInputSpaceDimension();


  typename TransformType::Pointer  transformHigh = TransformType::New();
  
  transformHigh->SetGridSpacing( spacingHigh );
  transformHigh->SetGridOrigin( originHigh );
  transformHigh->SetGridRegion( bsplineRegion );

  typedef typename TransformType::ParametersType ParametersType;
  ParametersType parametersHigh( transformHigh->GetNumberOfParameters() );
  parametersHigh.Fill( 0.0 );

  unsigned int counter = 0;

  for ( unsigned int k = 0; k < SpaceDimension; k++ )
    {
    typedef typename TransformType::ImageType ParametersImageType;
    typedef itk::ResampleImageFilter<ParametersImageType,ParametersImageType> ResamplerType;
    typename ResamplerType::Pointer upsampler = ResamplerType::New();

    typedef itk::BSplineResampleImageFunction<ParametersImageType,double> FunctionType;
    typename FunctionType::Pointer function = FunctionType::New();

    typedef itk::IdentityTransform<double, itk::GetImageDimension<ImageType>::ImageDimension> IdentityTransformType;
    typename IdentityTransformType::Pointer identity = IdentityTransformType::New();

    upsampler->SetInput( transfo->GetCoefficientImage()[k] );
    upsampler->SetInterpolator( function );
    upsampler->SetTransform( identity );
    upsampler->SetSize( transformHigh->GetGridRegion().GetSize() );
    upsampler->SetOutputSpacing( transformHigh->GetGridSpacing() );
    upsampler->SetOutputOrigin( transformHigh->GetGridOrigin() );

    typedef itk::BSplineDecompositionImageFilter<ParametersImageType,ParametersImageType>
      DecompositionType;
    typename DecompositionType::Pointer decomposition = DecompositionType::New();

    decomposition->SetSplineOrder( SplineOrder );
    decomposition->SetInput( upsampler->GetOutput() );
    decomposition->Update();

    typename ParametersImageType::Pointer newCoefficients = decomposition->GetOutput();

    // copy the coefficients into the parameter array
    typedef itk::ImageRegionIterator<ParametersImageType> Iterator;
    Iterator it( newCoefficients, transformHigh->GetGridRegion() );
    while ( !it.IsAtEnd() )
      {
      parametersHigh[ counter++ ] = it.Get();
      ++it;
      }

    }

  transformHigh->SetParametersByValue( parametersHigh );

  return transformHigh;

}



// FIXME: This version is here for backward compatibility and should be removed !!
template <class TransformType, class ImageType>
  typename TransformType::Pointer
RefineFFDTransformation( TransformType* transfo, unsigned int SplineOrder)
{
  const unsigned int SpaceDimension = transfo->GetInputSpaceDimension();

  typename TransformType::Pointer  transformHigh = TransformType::New();
  typedef typename TransformType::RegionType RegionType;
  typedef typename TransformType::ParametersType ParametersType;
  typedef typename TransformType::SpacingType SpacingType;
  typedef typename TransformType::OriginType OriginType;

  typename RegionType::SizeType   gridHighSize = transfo->GetGridRegion().GetSize();
  SpacingType spacingHigh = transfo->GetGridSpacing();
  OriginType  originHigh  = transfo->GetGridOrigin();

  for (unsigned int d=0;d<SpaceDimension; d++)
    {
		gridHighSize[d] = (gridHighSize[d]-1)*2 + 1;
		spacingHigh[d] /= 2.;
    }

  RegionType bsplineRegion;
  bsplineRegion.SetSize( gridHighSize );

  transformHigh->SetGridSpacing( spacingHigh );
  transformHigh->SetGridOrigin( originHigh );
  transformHigh->SetGridRegion( bsplineRegion );

  ParametersType parametersHigh( transformHigh->GetNumberOfParameters() );
  parametersHigh.Fill( 0.0 );

  unsigned int counter = 0;

  for ( unsigned int k = 0; k < SpaceDimension; k++ )
    {
    typedef typename TransformType::ImageType ParametersImageType;
    typedef itk::ResampleImageFilter<ParametersImageType,ParametersImageType> ResamplerType;
    typename ResamplerType::Pointer upsampler = ResamplerType::New();

    typedef itk::BSplineResampleImageFunction<ParametersImageType,double> FunctionType;
    typename FunctionType::Pointer function = FunctionType::New();

    typedef itk::IdentityTransform<double, itk::GetImageDimension<ImageType>::ImageDimension> IdentityTransformType;
    typename IdentityTransformType::Pointer identity = IdentityTransformType::New();

    upsampler->SetInput( transfo->GetCoefficientImage()[k] );
    upsampler->SetInterpolator( function );
    upsampler->SetTransform( identity );
    upsampler->SetSize( transformHigh->GetGridRegion().GetSize() );
    upsampler->SetOutputSpacing( transformHigh->GetGridSpacing() );
    upsampler->SetOutputOrigin( transformHigh->GetGridOrigin() );

    typedef itk::BSplineDecompositionImageFilter<ParametersImageType,ParametersImageType>
      DecompositionType;
    typename DecompositionType::Pointer decomposition = DecompositionType::New();

    decomposition->SetSplineOrder( SplineOrder );
    decomposition->SetInput( upsampler->GetOutput() );
    decomposition->Update();

    typename ParametersImageType::Pointer newCoefficients = decomposition->GetOutput();

    // copy the coefficients into the parameter array
    typedef itk::ImageRegionIterator<ParametersImageType> Iterator;
    Iterator it( newCoefficients, transformHigh->GetGridRegion() );
    while ( !it.IsAtEnd() )
      {
      parametersHigh[ counter++ ] = it.Get();
      ++it;
      }

    }

  transformHigh->SetParametersByValue( parametersHigh );

  return transformHigh;

}

#endif /* !REGFFDSPLIT_H_ */
