/*=========================================================================

Program:   Insight Segmentation & Registration Toolkit
Module:    $RCSfile: regMeanSquaresTemporalImageMetricBase.txx,v $
Language:  C++
Date:      $Date: 2008-02-03 19:00:36 $
Version:   $Revision: 1.51 $

Copyright (c) Insight Software Consortium. All rights reserved.
See ITKCopyright.txt or http://www.itk.org/HTML/Copyright.htm for details.

This software is distributed WITHOUT ANY WARRANTY; without even 
the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR 
PURPOSE.  See the above copyright notices for more information.

=========================================================================*/
#ifndef _regMeanSquaresTemporalImageMetricBase_txx
#define _regMeanSquaresTemporalImageMetricBase_txx

// First make sure that the configuration is available.
// This line can be removed once the optimized versions
// gets integrated into the main directories.
#include "itkConfigure.h"
#include "itkTimeProbe.h"

#include "regMeanSquaresTemporalImageMetricBase.h"

namespace itk
{


	/*
	* Initialize the metric 
	*/
	template <class TImage> 
	void
		MeanSquaresTemporalImageMetricBase<TImage>
		::Initialize(void) throw ( ExceptionObject )
	{
		Superclass::Initialize();
	}


	/*
	*
	*
	*/	 
	template <class TImage> 
	void
		MeanSquaresTemporalImageMetricBase<TImage>
		::GetSampleContributionOnTimeInterval(InputPointType inputPoint, double finalTime,  RealType fixedValue,
		MeasureType & value, DerivativeType  & derivative, SparseJacobianType& jacobian,
		unsigned long & numberOfValidPropagatedPoints) const
	{

		ImageConstPointer image = this->m_Image;

		const unsigned int Dimension = ImageType::ImageDimension;
		const unsigned int SpaceDimension = Dimension - 1;

		// Reset the jacobian
		jacobian.clear();

		bool lastTimeStepReached = false;
		double currentTime = inputPoint[SpaceDimension];
		double timeStep = image->GetSpacing()[SpaceDimension];
		bool valid = true; 

		if (finalTime<currentTime){
			timeStep *= -1;
		}else if (finalTime==currentTime) {
			lastTimeStepReached=true;
		}

		OutputPointType transformedPoint;
		SpacePointType transformedPoint2;

		while ( !lastTimeStepReached ) 
		{
			currentTime += timeStep;

			// Checking if we are at the last time step
			if ( ( (timeStep>0)&&(currentTime + timeStep > finalTime) )
				|| ( (timeStep<=0)&&(currentTime + timeStep < finalTime) ) )
			{
				lastTimeStepReached=true;
			}

			this->m_Transform->GetIncrementalSparseJacobian( inputPoint, transformedPoint, currentTime, jacobian, valid);

			// FIXME: Write adaptive time step strategy
			if (!valid)
			{
				itkExceptionMacro( << "Determinant of physical jacobian is negative! ");
			}

			if( this->m_Interpolator->IsInsideBuffer( transformedPoint ) )
			{

				bool skip_sample = this->IsWithinMovingMask(transformedPoint);

				if (!skip_sample) {

					const RealType movingValue  = this->m_Interpolator->Evaluate( transformedPoint );	        

					numberOfValidPropagatedPoints++;

					const RealType diff = movingValue - fixedValue; 

					value += diff * diff;

					GradientPixelType gradient;

					if (this->m_ComputeGradient) 
					{
						// Get the gradient by NearestNeighboorInterpolation: 
						// which is equivalent to round up the point components.
						typedef typename OutputPointType::CoordRepType CoordRepType;
						typedef ContinuousIndex<CoordRepType,ImageType::ImageDimension> 
							ImageContinuousIndexType;

						ImageContinuousIndexType tempIndex;
						this->m_Image->TransformPhysicalPointToContinuousIndex( transformedPoint, tempIndex );

						typename ImageType::IndexType mappedIndex; 
						mappedIndex.CopyWithRound( tempIndex );

						gradient = this->GetGradientImage()->GetPixel( mappedIndex );
					}
					else if (this->m_InterpolatorIsBSpline) 
					{
						gradient = this->m_InterpolatorBS->EvaluateDerivative(transformedPoint);
					}
					else
					{
						itkExceptionMacro(<< "Unsupported case for the computation of image derivatives. " <<std::endl);
					}

					// Project the jacobian on the image gradient
					typename SparseJacobianType::iterator it = jacobian.begin();

					while (it != jacobian.end())
					{

						double sum = 0;
						unsigned long col = it->first;

						for (unsigned int dim=0; dim<SpaceDimension; dim++)
						{
							sum += it->second[dim] * gradient[dim];
						}
						derivative[ col ] += 2.0 * diff * sum;

						++it; 	  
					}

				} // FI skip_sample
			}// FI Inside buffer

			inputPoint=transformedPoint;


		}// ELIHW lastTimeStepReached

	}



} // end namespace itk


#endif

