/*=========================================================================

Program:   Insight Segmentation & Registration Toolkit
Module:    $RCSfile: regMeanSquaresTemporalImageMetric.txx,v $
Language:  C++
Date:      $Date: 2008-02-03 19:00:36 $
Version:   $Revision: 1.51 $

Copyright (c) Insight Software Consortium. All rights reserved.
See ITKCopyright.txt or http://www.itk.org/HTML/Copyright.htm for details.

This software is distributed WITHOUT ANY WARRANTY; without even 
the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR 
PURPOSE.  See the above copyright notices for more information.

=========================================================================*/
#ifndef _regMeanSquaresTemporalImageMetric_txx
#define _regMeanSquaresTemporalImageMetric_txx

// First make sure that the configuration is available.
// This line can be removed once the optimized versions
// gets integrated into the main directories.
#include "itkConfigure.h"
#include "itkTimeProbe.h"

#include "regMeanSquaresTemporalImageMetric.h"

namespace itk
{

	/*
	* Constructor
	*/
	template <class TImage> 
	MeanSquaresTemporalImageMetric<TImage>
		::MeanSquaresTemporalImageMetric()
	{

	}

	/*
	* Get the match Measure
	*/
	template <class TImage> 
	typename MeanSquaresTemporalImageMetric<TImage>::MeasureType
		MeanSquaresTemporalImageMetric<TImage>
		::GetValue( const TransformParametersType & parameters ) const
	{

		MeasureType measure = NumericTraits< MeasureType >::Zero;
		DerivativeType derivative( this->GetNumberOfParameters());
		derivative.Fill(NumericTraits< MeasureType >::Zero);
		GetValueAndDerivative(parameters, measure, derivative);
		return measure;
	}




	/*
	* Get the Derivative Measure
	*/
	template < class TImage> 
	void
		MeanSquaresTemporalImageMetric<TImage>
		::GetDerivative( const TransformParametersType & parameters,
		DerivativeType & derivative  ) const
	{

		MeasureType measure = NumericTraits< MeasureType >::Zero;
		GetValueAndDerivative(parameters, measure, derivative);
	}



	/*
	* Get both the match Measure and the Derivative Measure 
	*/
	template <class TImage> 
	void
		MeanSquaresTemporalImageMetric<TImage>
		::GetValueAndDerivative(const TransformParametersType & parameters, 
		MeasureType & value, DerivativeType  & derivative) const
	{

		if( !this->m_Initialized )
		{
			itkExceptionMacro(<<"You forgot to call Initialize()");
		}

		// Debug
		//
		std::cout << "    min p=" << parameters.min_value() << "    max p=" << parameters.max_value() << std::endl;
		// Debug


		ImageConstPointer image = this->m_Image;

		const unsigned int Dimension = ImageType::ImageDimension;
		const unsigned int SpaceDimension = Dimension - 1;

		const unsigned int ParametersDimension = this->GetNumberOfParameters();

		// Computing SSD and its derivatives
		MeasureType measure = NumericTraits< MeasureType >::Zero;

		derivative = DerivativeType( ParametersDimension );
		derivative.Fill( NumericTraits<ITK_TYPENAME DerivativeType::ValueType>::Zero );

		unsigned long numberOfValidPropagatedPoints = 0;

		double timeOrigin = image->GetOrigin()[SpaceDimension];
		double timeEnd = timeOrigin + 
			(double)image->GetSpacing()[SpaceDimension] * 
			(double)(this->GetImageRegion().GetSize()[SpaceDimension]-1);

		this->SetTransformParameters( parameters );

		typedef typename TransformType::SparseJacobianType SparseJacobianType;
		SparseJacobianType jacobian;

		typedef typename Superclass::DiffeoSpatialSample SampleType;
		typename std::vector<SampleType>::const_iterator ti = this->m_SamplesContainer.begin();

		while ( ti != this->m_SamplesContainer.end())
		{
			InputPointType inputPoint = ti->m_FixedPoint;

			// Going backward
			GetSampleContributionOnTimeInterval(inputPoint, timeOrigin,  ti->m_FixedValue,
				measure, derivative, jacobian,
				numberOfValidPropagatedPoints);

			// Going forward
			GetSampleContributionOnTimeInterval(inputPoint, timeEnd,  ti->m_FixedValue,
				measure, derivative, jacobian,
				numberOfValidPropagatedPoints);

			++ti;

		}

		if( !numberOfValidPropagatedPoints )
		{
			itkExceptionMacro(<<"All the points mapped to outside of the moving image");
		}
		else
		{
			for(unsigned int i=0; i<ParametersDimension; i++)
			{
				derivative[i] /= static_cast<double>( numberOfValidPropagatedPoints );
			}
			measure /= static_cast<double>( numberOfValidPropagatedPoints );
		}

		std::cout << "SSD : " << measure << std::endl;

		// By default, only the SSD value is returned

		value = measure;
	}


	/*
	* Initialize the metric 
	*/
	template <class TImage> 
	void
		MeanSquaresTemporalImageMetric<TImage>
		::Initialize(void) throw ( ExceptionObject )
	{
		Superclass::Initialize();
		// Sparse fixed image domain
		this->SampleReferenceImageDomain(false, 0);
	}





} // end namespace itk


#endif

