#ifndef regDeformPolyDataFilterTDFFD_H
#define regDeformPolyDataFilterTDFFD_H

#include "itkProcessObject.h"
#include "itkDataObject.h"
#include "DataWrappers.h"

// Image types and other standard things
#include "itkImage.h"
#include "itkImageFileReader.h"
#include "itkTransformFileReader.h"
#include "itkImageMaskSpatialObject.h"
#include "vtkPolyData.h"
#include "vtkPolyDataReader.h"
#include "vtkPolyDataWriter.h"
#include "vtkPolyDataNormals.h"
#include "vtkStructuredPointsReader.h"
#include "vtkUnstructuredGrid.h"
#include "vtkUnstructuredGridWriter.h"
#include "vtkIdList.h"
#include <vtkImageData.h>
#include <vtkStructuredPoints.h>
#include <vtkImageDilateErode3D.h>
#include <vtkSmartPointer.h>
#include <vtkDoubleArray.h>
#include <vtkFloatArray.h>
#include <vtkPointData.h>
#include <vnl/vnl_cross.h>

// Our temporal diff matching stuff
#include "regDiffeomorphicContinuousBSplineTransform.h"
#include "regTransformWrapper.h"
#include <vtkClipPolyData.h>
#include <vtkImplicitVolume.h>

class ITK_EXPORT regDeformPolyDataFilterTDFFD : public itk::ProcessObject
{
private:
	//! Image dimension is currently set to 3 by default. FIXME: This should be a template parameter.
	const static unsigned int m_SpaceDimension = 3;
	const static int m_BSplineOrder = 3;
	const static unsigned int m_Dimension = m_SpaceDimension + 1;
public:
	typedef regDeformPolyDataFilterTDFFD Self;
	typedef itk::ProcessObject ParentType;
	typedef itk::SmartPointer<Self> Pointer;
	typedef itk::SmartPointer<const Self>  ConstPointer;
	typedef vtkPolyData MeshType;
	typedef vtkUnstructuredGrid VolMeshType;
	typedef itk::ImageMaskSpatialObject< m_SpaceDimension > MaskSpatialObjectType;
	typedef std::vector<vtkSmartPointer<MeshType> > MeshCollectionType;
	typedef std::vector<vtkSmartPointer<VolMeshType> > UgCollectionType;
	
	// Image types
	typedef  unsigned char    PixelType;
	typedef itk::Image< PixelType, m_SpaceDimension >  ImageType;
	// TDFFD types
	typedef itk::DiffeomorphicContinuousBSplineTransform<double,m_Dimension, m_BSplineOrder>  TransformType;
    typedef TransformType::Pointer TransformPointerType;
	
	///** Run-time type information (and related methods). */
	itkTypeMacro(regDeformPolyDataFilterTDFFD,ParentType);
	itkFactorylessNewMacro(regDeformPolyDataFilterTDFFD);
	
	//! Constructor
	regDeformPolyDataFilterTDFFD();
	
    // Set/Gets time origin/spacing and number of time steps
	itkSetMacro(TimeSpacing, double);
    itkGetMacro(TimeSpacing, double);
    itkSetMacro(TimeOrigin, double);
    itkGetMacro(TimeOrigin, double);
    itkSetMacro(NumberOfTimeSteps, unsigned int);
    itkGetMacro(NumberOfTimeSteps, unsigned int);
    
    // Access output
	MeshCollectionType GetOutput();
	UgCollectionType GetOutputUg();

    // Set/Get axis computation method
    itkSetMacro(AxisCompMethod, unsigned int);
    itkGetMacro(AxisCompMethod, unsigned int);

	void SetLongAxis (const vnl_vector<double>& la)
	{
		m_LongAxis = la;
        this->Modified();
	}
    
    // Set/Get the transform
    itkSetObjectMacro(InputTransform, TransformType);
    itkGetObjectMacro(InputTransform, TransformType);
    
    // Set/Get input mesh
    void SetInputMesh (MeshType* mesh)
    {
        m_InputMesh = mesh;
        
        MeshContainer::Pointer container = MeshContainer::New();
        container->AddInput(mesh);
        this->ProcessObject::SetNthInput(0, static_cast< MeshContainer* >(container) );
        this->Modified();
    }
    
    // Set/Get vol input mesh
    void SetInputVolMesh (VolMeshType* mesh)
    {
        // Process object is not const-correct so the const_cast is required here
        m_InputVolMesh = mesh;
        this->Modified();
    }
    
    // Get volumetric mode
    itkGetMacro(VolumetricMode, bool);

    
    void Update();

private:

    void GenerateData();
	MaskSpatialObjectType::Pointer m_Mask;
	unsigned int m_NumberOfTimeSteps;
	double m_TimeSpacing;
    double m_TimeOrigin;
	unsigned int m_AxisCompMethod;
	std::vector< vtkSmartPointer<MeshType> > m_SmartOutput;
	std::vector< vtkSmartPointer<VolMeshType> > m_SmartOutput2;

	bool m_VolumetricMode;

	vnl_vector<double> m_LongAxis;
    
	// Keep track of created strain arrays
	std::vector< vtkSmartPointer<vtkFloatArray> >	m_LongitudinalStrains;
	std::vector< vtkSmartPointer<vtkFloatArray> >	m_CircumferentialStrains;
	std::vector< vtkSmartPointer<vtkFloatArray> >	m_RadialStrains;
	std::vector< vtkSmartPointer<vtkFloatArray> >	m_FullStrains;
	std::vector< vtkSmartPointer<vtkFloatArray> >	m_InsideArrays;
	std::vector< vtkSmartPointer<vtkPoints> >		m_OutputPoints;
    
    // Inputs
    vtkSmartPointer<MeshType>       m_InputMesh;
    vtkSmartPointer<VolMeshType>    m_InputVolMesh;
    TransformPointerType            m_InputTransform;
    
    
};

#endif //regDeformPolyDataFilterTDFFD_H
