#ifndef _blStrainComputationFilter_H
#define _blStrainComputationFilter_H

#include "itkProcessObject.h"
#include "itkDataObject.h"

#include "vtkSmartPointer.h"
#include "vtkPolyData.h"
#include "vtkFloatArray.h"
#include "vtkPointData.h"
#include "vtkPoints.h"
#include "vtkPolyDataNormals.h"
#include "vtkCell.h"
#include "vtkCellData.h"
#include "vtkCellType.h"
#include "DataWrappers.h"


#include <vnl/vnl_vector.h>	

class ITK_EXPORT regStrainComputationFilter : public itk::ProcessObject
{
protected:
	//! Image dimension is currently set to 3 by default. 
	const static unsigned int m_ImageDimension = 3;

	//! Struct with all strain values
	struct StrainValues{

		double currentStrainValRad;
		double currentStrainValCirc;
		double currentStrainValLong;

		StrainValues( )
		{
			currentStrainValRad = 0.;
			currentStrainValCirc = 0.;
			currentStrainValLong = 0.;
		}
	};
public:
	//! Type of strains
	enum DirectionStrainProjection {
		Radial,
		Longitudinal,
		Circumferential,
	};

	typedef regStrainComputationFilter Self;
	typedef itk::ProcessObject ParentType;
	typedef itk::SmartPointer<Self> Pointer;
	typedef itk::SmartPointer<const Self>  ConstPointer;
	typedef vtkPolyData MeshType;
    typedef MeshContainer MeshCollectionType;
	
	///** Run-time type information (and related methods). */
	itkTypeMacro(Self,ParentType);
	itkFactorylessNewMacro(Self);

public:
	//! Constructor
	regStrainComputationFilter();

	//! This filter takes a mesh as input
	virtual void SetInput(
		const MeshCollectionType* inputMesh );

	//! Access the input of the filter (wrapper object around vtkPolyData)
	MeshCollectionType* GetInput();

	//! Update filter
	virtual void Update();

    // Set/Get axis computation method
    itkSetMacro(AxisCompMethod, unsigned int);
    itkGetMacro(AxisCompMethod, unsigned int);
    
    void SetLongAxis (const vnl_vector<double>& la)
	{
		m_LongAxis = la;
        this->Modified();
	}


protected:
	//!
	virtual void GenerateData();
	
	//! This method is called by Generate Data 
	virtual void ComputeStrain();

	//! Initial preprocess for each time step
	void PreProcessTimeStep( );
	
	virtual void ComputeStrainOnMesh( 
		vtkPolyData* inputMesh,
		unsigned int cellId,
		vtkCell* cell,
		vnl_vector<double>  &radialDirection,
		vnl_vector<double>  &circumDirection,
		vnl_vector<double>  &longDirection,
		StrainValues &result );


protected:
	//! Filter for computing normal directions
	vtkSmartPointer<vtkPolyDataNormals> m_NormalsFilter;
    
	unsigned int m_CurrentTimeStep;
	unsigned int m_AxisCompMethod;
    
    vnl_vector<double> m_LongAxis;
};



#endif //_blStrainComputationFilter_H
