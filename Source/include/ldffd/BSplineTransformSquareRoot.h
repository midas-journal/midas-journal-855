/*
** BSplineTransformSquareRoot.h
** 
** Made by Mathieu De Craene
** Login   <mdecraene@cilab-mathieu>
** 
*/

#ifndef   	BSPLINETRANSFORMSQUAREROOT_H_
# define   	BSPLINETRANSFORMSQUAREROOT_H_

#include "itkImageRegion.h"
#include "itkImage.h"
#include "BSplineDeformableTransformOpt.h"
#include "itkImageRandomConstIteratorWithIndex.h"

template <class ImageType, class TransformType>
typename TransformType::ParametersType 
ComputeSquareRootTransform(const TransformType * transform, 
			   const ImageType * inputImage,  
			   const typename ImageType::RegionType inputRegion,
			   unsigned int numberOfIterations,
			   double learningRate)
{
  
  
  // Initializing output transformation
  typename TransformType::Pointer outputTransform = TransformType::New();
  
  // This should copy the information about the grid region/spacing/origin
  outputTransform->SetFixedParameters(transform->GetFixedParameters());
  typename TransformType::ParametersType currentParameters = transform->GetParameters();
  
  for (unsigned int index=0; index<transform->GetNumberOfParameters(); index++)
    {
    currentParameters[index] /= 2.;
    }

  return currentParameters;

/*
  
  // Starting to iterate along gradient direction
  typedef typename TransformType::OutputPointType PointType;
  typedef typename TransformType::ParametersType ParametersType;
  
  ParametersType paramBuffer = ParametersType(outputTransform->GetNumberOfParameters());
  
  // Random iterator over image region
  typedef typename itk::ImageRandomConstIteratorWithIndex< ImageType > RandomIteratorType;
  RandomIteratorType iter(inputImage, inputRegion);

  iter.SetNumberOfSamples( inputRegion.GetNumberOfPixels() / 50 );

  // Classical jacobian matrix
  typename TransformType::ClassicalJacobianType jacobian(ImageDimension, ImageDimension, 0.);
  PointType transformDiff; 
  for (unsigned int it=0; it<numberOfIterations; it++)
    {
    paramBuffer.Fill(0.);
    outputTransform->SetParameters(currentParameters);

    iter.GoToBegin();
    PointType inCoord;
  
    double diffMagnitude=0.;
  
    for (iter.GoToBegin(); !iter.IsAtEnd(); ++iter)
      {
      // input coordinate corresponding to this sample
      inputImage->TransformIndexToPhysicalPoint(iter.GetIndex(), inCoord);

      // Computing the 2 output coordinates
      PointType outCoord1 = outputTransform->TransformPoint(inCoord);
      PointType outCoord2 = outputTransform->TransformPoint(outCoord1);
      // ... and the output coordinate obtained in 1 step only
      PointType outCoordOneStep = transform->TransformPoint(inCoord);

      // Compute difference between the transformed coordinate using
      // the full step and the transformed coordinate doing the two half steps.
      for  (unsigned int dim=0; dim<ImageDimension; dim++)
        {
        transformDiff[dim] = outCoord2[dim] - outCoordOneStep[dim];
        diffMagnitude += transformDiff[dim] * transformDiff[dim]; 
        }

      // Call to classical jacobian
      bool valid = outputTransform->GetClassicalJacobian(outCoord1, jacobian);

      // call to ITK jacobian of the two transformations
      //
      typename TransformType::WeightsType 
        weights1(outputTransform->GetNumberOfWeights()), 
        weights2(outputTransform->GetNumberOfWeights());
      typename TransformType::ParameterIndexArrayType
        indexes1(outputTransform->GetNumberOfWeights()), 
        indexes2(outputTransform->GetNumberOfWeights());

      outputTransform->GetJacobian(inCoord, indexes1, weights1);
      outputTransform->GetJacobian(outCoord1, indexes2, weights2);

      // Now we need to calculate the product of 
      // CLASSICAL_JAC * ITK_JAC1 and sum it with ITK_JAC2. We can try
      // to do this by exploiting the sparse nature of the jacobians
      // matrices.
      //
      for (unsigned int index=0; index<indexes1.Size(); index++)
        {
        double weight1=weights1[index];  
        double weight2=weights2[index];
        // Each index acts on ImageDimension parameters. 
        // For each dimension(dim), we compute the index of
        // the corresponding parameter (totalIndex). 
        // This index corresponds 
        // to a non-zero column in CLASSICAL_JAC * ITK_JAC1 
        // ... And This column is the dim_th column of CLASSICAL_JAC
        // multiplied by weight.
        for  (unsigned int dim=0; dim<ImageDimension; dim++)
          {
          // full index corresponding to index1
          unsigned int totalIndex1 = indexes1[index] 
            + outputTransform->GetNumberOfParametersPerDimension() * dim;
          // Extract the dim_th column of CLASSICAL_JAC and multiply it by weight
          double sum = 0.;
          for  (unsigned int dim2=0; dim2<ImageDimension; dim2++)
            {
            sum+= jacobian[dim2][dim] * weight1 *  transformDiff[dim2];
            }
          paramBuffer[totalIndex1] += sum; 
          // full index corresponding to index2
          unsigned int totalIndex2 = indexes2[index]
            + outputTransform->GetNumberOfParametersPerDimension() * dim;
          // We now have the contribution of ITK_JAC2 (easier) 
          paramBuffer[totalIndex2] += transformDiff[dim]*weight2;
          }
        }

      }//ELIHW iter is not at end

      diffMagnitude /=  (double)iter.GetNumberOfSamples();

      std::cout<<it<<" "<<diffMagnitude<<std::endl;

      // Updating parameters
      double meanChange = 0.;
      for (unsigned int index=0;index<currentParameters.Size(); index++)
        {
        double update = -1. * learningRate
          * paramBuffer[index] / (double)iter.GetNumberOfSamples();
        currentParameters[index] += update;
        meanChange += update*update;
        }
      meanChange /= double(currentParameters.Size());
      meanChange = sqrt(meanChange);
      
      std::cout << "Mean parameter change: " << meanChange << std::endl;

      // Appliying new set of parameters to output transformation
      outputTransform->SetParameters(currentParameters);

    }//ROF iterations


  return outputTransform->GetParameters();
*/
}

template <class ImageType, class TransformType>
typename TransformType::ParametersType 
AddOneStepInDiffemorphicTransform(TransformType * transform,
				  const unsigned int pos,
				  const ImageType * sqrt_inputImage,  
				  const typename ImageType::RegionType sqrt_inputRegion,
				  unsigned int sqrt_numberOfIterations,
				  double sqrt_learningRate)
{

  // Compute square root of the pos_th transform
  const typename TransformType::InternalTransformType*
    internalTr = transform->GetInternalTransform(pos);

  typename TransformType::ParametersType sqrtParams = 
    ComputeSquareRootTransform
    <ImageType, typename TransformType::InternalTransformType>(internalTr, 
						      sqrt_inputImage,  
						      sqrt_inputRegion,
						      sqrt_numberOfIterations, 
						      sqrt_learningRate); 
  

  // Defining the right parameters with the result 
  // of the square root computation
  const typename TransformType::ParametersType& oldParameters = transform->GetParameters();

  
  // Insert transform at the right position
  transform->InsertTransform(pos, false);
  
  // New set of parameters
  typename TransformType::ParametersType params(transform->GetNumberOfParameters());
  
  const typename 
    TransformType::ParametersType::ValueType* data_block_current = oldParameters.data_block();
  typename TransformType::ParametersType::ValueType* data_block_new = params.data_block();
  unsigned int numberOfParametersByTransform = internalTr->GetNumberOfParameters();
  
  // Copy first part of parameters (before pos)
  unsigned int index=0;
  while ( index<numberOfParametersByTransform*pos )
    {
    *data_block_new = *data_block_current;
    
    data_block_current++;
    data_block_new++;
    index++;
    }
  
  // std::cout << "First part of parameters copied" << std::endl;

  // Copy parameters of square root transform (cover a lenght of 2
  // times the number of parameters of one BSpline transform)
  unsigned int counter = 0;
  while ( index < numberOfParametersByTransform * (pos+2) )
    {
    unsigned int indexSqrt = 
      counter%sqrtParams.Size();
    *data_block_new 
      = sqrtParams[indexSqrt];
    data_block_new++;
    index++;
    counter++;
    }
  
  // std::cout << "Second part of parameters copied" << std::endl;
  
  // Copy segond part (after pos)
  while (index<params.Size())
    {
    *data_block_new = *data_block_current;
    
    data_block_current++;
    data_block_new++;
    index++;
    }
  
  // std::cout << "Third part of parameters copied" << std::endl;

  return params;

}

#endif 	    /* !BSPLINETRANSFORMSQUAREROOT_H_ */
