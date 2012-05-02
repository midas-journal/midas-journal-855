#ifndef __regBSplineFieldRefine_h
#define __regBSplineFieldRefine_h

#include <iostream>
#include "itkObject.h"
#include "itkObjectFactory.h"

#include "regBSplineField.h"
#include "itkBSplineDecompositionImageFilter.h"
#include "itkBSplineResampleImageFunction.h"
#include "itkResampleImageFilter.h"

namespace itk
{
  
/** \class BSplineFieldRefine
 *
 *
 */
template <
  class TScalarType = double,          // Data type for scalars
  unsigned int NDimensions = 4,        // Number of space dimensions 
  unsigned int VSplineOrder = 3 >      // Spline order
class ITK_EXPORT  BSplineFieldRefine  : public Object
{

public:

  /** Standard class typedefs. */
  typedef BSplineFieldRefine  Self;
  typedef Object Superclass;
  typedef SmartPointer< Self >   Pointer;
  typedef SmartPointer< const Self >  ConstPointer;

  /** New macro for creation of through the object factory.*/
  itkNewMacro( Self );

  /** BSpline object field and related types. */
  typedef BSplineField<TScalarType, NDimensions, VSplineOrder> BSplineFieldType;
  typedef typename BSplineFieldType::Pointer BSplineFieldPointerType;
  typedef typename BSplineFieldType::SpacingType SpacingType;
  typedef typename BSplineFieldType::OriginType OriginType;
  typedef typename BSplineFieldType::RegionType RegionType;
  typedef typename BSplineFieldType::ImageType ImageType;
  typedef typename BSplineFieldType::SizeType SizeType;
  typedef typename BSplineFieldType::ParametersType ParametersType;

  typedef BSplineResampleImageFunction<ImageType,double> InterpolatorType;
  typedef ResampleImageFilter<ImageType,ImageType> ResamplerType;
  typedef typename ResamplerType::Pointer ResamplerPointerType;
  typedef BSplineDecompositionImageFilter<ImageType,ImageType> DecompositionType;
  typedef typename DecompositionType::Pointer DecompositionPointerType;

  void SetInput(BSplineFieldType*);
  void ConfigureOutput(const SpacingType&, const OriginType&, const RegionType&);
  void Update(void);
  BSplineFieldType* GetOutput(void);

protected:

  BSplineFieldPointerType m_BSplineInputField;
  BSplineFieldPointerType m_BSplineOutputField;
  ParametersType m_OutputParameters;

  BSplineFieldRefine(); 
  virtual ~BSplineFieldRefine() {};

private:
  BSplineFieldRefine(const Self&); //purposely not implemented
  void operator=(const Self&); //purposely not implemented

};

}// end namespace itk

#if ITK_TEMPLATE_TXX
#include "regBSplineFieldRefine.txx"
#endif

#endif
