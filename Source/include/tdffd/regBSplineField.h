#ifndef __regBSplineField_h
#define __regBSplineField_h

#include <iostream>
#include "itkObject.h"
#include "itkObjectFactory.h"

#include "itkPoint.h"
#include "itkVector.h"
#include "itkCovariantVector.h"
#include "vnl/vnl_vector_fixed.h"
#include "itkArray.h"
#include "itkArray2D.h"

#include "itkImage.h"
#include "itkImageRegion.h"
#include "itkBSplineInterpolationWeightFunction.h"
#include "itkBSplineInterpolationWeightFunctionDerivative.h"
#include "itkBSplineInterpolationWeightFunctionSecondDerivative.h"
#include "itkImageRegionConstIterator.h"
#include "itkImageRegionIterator.h"

namespace itk
{
  
/** \class BSplineField
 *
 * This class is ...
 *
 */
template <
  class TScalarType = double,          // Data type for scalars
  unsigned int NDimensions = 4,        // Number of space dimensions 
  unsigned int VSplineOrder = 3 >      // Spline order
class ITK_EXPORT  BSplineField  : public Object
{

public:

  /** Standard class typedefs. */
  typedef BSplineField  Self;
  typedef Object Superclass;
  typedef SmartPointer< Self >   Pointer;
  typedef SmartPointer< const Self >  ConstPointer;

  /** New macro for creation of through the object factory.*/
  itkNewMacro( Self );

  /** Type of the input parameters. */
  typedef  Array< double >           ParametersType;

  /** Parameter index array type. */
  typedef  Array<unsigned long>      ParameterIndexArrayType;

  /** Dimension of the space domain space. */
  itkStaticConstMacro(SpaceDimension, unsigned int, NDimensions-1);

  /** Dimension of the spatio-temporal space. */
  itkStaticConstMacro(SpacePlusTimeDimension, unsigned int, NDimensions);

  /** BSpline order. */
  itkStaticConstMacro(SplineOrder, unsigned int, VSplineOrder);

  /** Point type. */
  typedef Point<TScalarType, SpacePlusTimeDimension> PointType;
  typedef Point<TScalarType, SpaceDimension> SpacePointType;

  /* Classical jacobian type. By classical, we mean the jacobian defined
   * as d(T(x))/dx
   */
  typedef vnl_matrix<double> ClassicalJacobianType;

  /** Run-time type information (and related methods). */
  itkTypeMacro( BSplineField, Object );

  /** Return the number of parameters that completely define the Transfom  */
  virtual unsigned int GetNumberOfParameters(void) const; 

  /** Return the number of parameters per dimension */
  unsigned int GetNumberOfParametersPerDimension(void) const;

  /** Get the Transformation Parameters. */
  virtual const ParametersType& GetParameters(void) const;

  /** Set the transformation parameters and update internal transformation. */
  virtual void SetParameters( const ParametersType & );

  /** Parameters as SpaceDimension number of images. */
  typedef typename ParametersType::ValueType PixelType;
  typedef Image<PixelType,itkGetStaticConstMacro(SpacePlusTimeDimension)> ImageType;

  typedef typename ImageType::Pointer ImagePointer;

  /** Get the array of coefficient images. */
  virtual ImagePointer * GetCoefficientImage()
  { return m_CoefficientImage; }
  virtual const ImagePointer * GetCoefficientImage() const
  { return m_CoefficientImage; }

  virtual void SetCoefficientImage( ImagePointer images[] );  

  /** Typedefs for specifying the extend to the grid. */
  typedef ImageRegion<itkGetStaticConstMacro(SpacePlusTimeDimension)>    RegionType;

  typedef typename RegionType::IndexType  IndexType;
  typedef typename RegionType::SizeType   SizeType;
  typedef typename ImageType::SpacingType SpacingType;
  typedef typename ImageType::DirectionType DirectionType;
  typedef typename ImageType::PointType   OriginType;

  /** This method specifies the region over which the grid resides. */
  virtual void SetGridRegion( const RegionType& region );
  itkGetMacro( GridRegion, RegionType );
  itkGetConstMacro( GridRegion, RegionType );

  /** This method specifies the grid spacing or resolution. */
  virtual void SetGridSpacing( const SpacingType& spacing );
  itkGetMacro( GridSpacing, SpacingType );
  itkGetConstMacro( GridSpacing, SpacingType );

  /** This method specifies the grid directions . */
  virtual void SetGridDirection( const DirectionType & spacing );
  itkGetMacro( GridDirection, DirectionType );
  itkGetConstMacro( GridDirection, DirectionType );

  /** This method specifies the grid origin. */
  virtual void SetGridOrigin( const OriginType& origin );
  itkGetMacro( GridOrigin, OriginType );
  itkGetConstMacro( GridOrigin, OriginType );

   /** Interpolation weights function type. */
  typedef BSplineInterpolationWeightFunction<TScalarType,
    itkGetStaticConstMacro(SpacePlusTimeDimension),
    itkGetStaticConstMacro(SplineOrder)>                      WeightsFunctionType;
  typedef typename WeightsFunctionType::WeightsType           WeightsType;
  typedef typename WeightsFunctionType::ContinuousIndexType   ContinuousIndexType;

  typedef BSplineInterpolationWeightFunctionDerivative<TScalarType, 
    itkGetStaticConstMacro(SpacePlusTimeDimension),
    itkGetStaticConstMacro(SplineOrder)> WeightFunctionDerivativeType;
  typedef typename WeightFunctionDerivativeType::Pointer 
    WeightFunctionDerivativePointerType;
 
  typedef BSplineInterpolationWeightFunctionSecondDerivative<TScalarType, 
  itkGetStaticConstMacro(SpacePlusTimeDimension),
  itkGetStaticConstMacro(SplineOrder)> WeightFunctionSecondDerivativeType;
  typedef typename WeightFunctionSecondDerivativeType::Pointer 
  WeightFunctionSecondDerivativePointerType;

  /** Set Weight function. */
  itkSetObjectMacro( WeightsFunction, WeightsFunctionType );
  itkGetObjectMacro( WeightsFunction, WeightsFunctionType );

  /** Get point velocity. */ 
  virtual void GetPointVelocity(const PointType  & point, SpacePointType & outputVelocity, 
									WeightsType & weights, 
									ParameterIndexArrayType & indices,
									bool& inside )  const;

   /** This is the derivative of the veloctiy field with respect to 
    ** the parameters, neglecting the dependency of x(p)	*/

  virtual void GetVelocityJacobian( const PointType & inputPoint,
										WeightsType & weights,
									    ParameterIndexArrayType & indices ) const;

  virtual bool GetVelocityDivergence ( const PointType & inputPoint,
										double & div, 
										WeightsType & weights,
									    ParameterIndexArrayType & indices ) const;
 
  virtual bool GetVelocityDivergenceHighOrder ( const PointType & inputPoint,
                                               WeightsType & weights) const;

  bool GetClassicalJacobian (const PointType & point, ClassicalJacobianType& jacobianMatrix, double timeStep) const;

  /** Wrap flat array into images of coefficients. */
  void WrapAsImages();

  /** Return the region of the grid wholly within the support region */
  itkGetConstReferenceMacro( ValidRegion, RegionType );

  void TransformPointToContinuousIndex( const PointType & point, ContinuousIndexType & index ) const;
    unsigned int GetNumberOfAffectedWeights() const;

protected:

 /** Variables defining the coefficient grid extend. */
  RegionType    m_GridRegion;
  SpacingType   m_GridSpacing;
  DirectionType m_GridDirection;
  OriginType    m_GridOrigin;
  RegionType    m_ValidRegion;

  DirectionType m_PointToIndex;
  DirectionType m_IndexToPoint;

  /** Variables defining the interpolation support region. */
  unsigned long m_Offset;
  bool          m_SplineOrderOdd;
  SizeType      m_SupportSize;
  IndexType     m_ValidRegionLast;
  IndexType     m_ValidRegionFirst;

  /** Array holding images wrapped from the flat parameters. */
  ImagePointer   m_WrappedImage[SpaceDimension];

  /** Array of images representing the B-spline coefficients 
   *  in each spatial dimension. */
  ImagePointer   m_CoefficientImage[SpaceDimension];

  /** Keep a pointer to the input parameters. */
  const ParametersType *  m_InputParametersPointer;

  /** Pointer to function used to compute Bspline interpolation weights. */
  typename WeightsFunctionType::Pointer  m_WeightsFunction;

  /** Pointer to function used to compute Bspline interpolation weights. */
  WeightFunctionDerivativePointerType m_weightsFunctionDerivative;
 
  /** Pointer to function used to compute second order Bspline interpolation weights. */
  WeightFunctionSecondDerivativePointerType m_weightsFunctionSecondDerivative;

  /** Check if a continuous index is inside the valid region. */
  bool InsideValidRegion( const ContinuousIndexType& index ) const;

  BSplineField(); 
  virtual ~BSplineField() {};

private:
  BSplineField(const Self&); //purposely not implemented
  void operator=(const Self&); //purposely not implemented
  

};

}// end namespace itk

#if ITK_TEMPLATE_TXX
#include "regBSplineField.txx"
#endif

#endif
