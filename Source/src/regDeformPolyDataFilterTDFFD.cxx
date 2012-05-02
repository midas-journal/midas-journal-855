#include "regDeformPolyDataFilterTDFFD.h"
#include <vnl/vnl_cross.h>


/**
 */
regDeformPolyDataFilterTDFFD::regDeformPolyDataFilterTDFFD()
{
	this->SetNumberOfRequiredInputs(1);
    this->SetNumberOfRequiredOutputs( 0 );
	this->m_Mask= 0;
    
	this->m_NumberOfTimeSteps=0;
	this->m_TimeSpacing=1.;
	this->m_TimeOrigin=0.;
    
	this->m_AxisCompMethod=2;
    
	this->m_SmartOutput.resize(0);
    
	this->m_VolumetricMode=false;
    
	m_LongAxis.set_size(m_SpaceDimension);
	m_LongAxis.fill(0.);
    
    m_InputVolMesh = 0;
    m_InputMesh = 0;
    
	return;
}

void regDeformPolyDataFilterTDFFD::Update()
{
    
	// Notify start event observers
	this->InvokeEvent( itk::StartEvent() );
    
	// Actually do something
	this->GenerateData();
    
	// Notify end event observers
	this->InvokeEvent( itk::EndEvent() );
    
}


/**
 This method effectively produces the ouput of the filter
 - set of deformed meshes using the input 3D or 4D transformation
 - strain computed on each node of the deformed meshes
 */
void regDeformPolyDataFilterTDFFD::GenerateData()
{
    
	m_LongitudinalStrains.clear();
	m_CircumferentialStrains.clear();
	m_RadialStrains.clear();
	m_FullStrains.clear();
	m_InsideArrays.clear();
	m_OutputPoints.clear();
    
	m_SmartOutput2.clear();
	m_SmartOutput.clear();
    
	MeshType* inputMesh = m_InputMesh;
	VolMeshType* inputVolMesh = m_InputVolMesh;
    TransformPointerType transfo = m_InputTransform;
    
    // Sanity checks
    if (inputMesh == 0)
        itkExceptionMacro(<< "An input mesh has to be set.");
    
    if (transfo.IsNull())
        itkExceptionMacro(<< "An input transform has to be set.");
    
    if (m_LongAxis.one_norm() == 0.)
        itkExceptionMacro(<< "Long axis has to be set in volumetric mode");
    
    if (m_InputVolMesh != 0)
        this->m_VolumetricMode = true;
    
	// Longitudinal, radial and circumferential directions for strain computation
	vnl_vector<double> longDirection(m_SpaceDimension);
	vnl_vector<double> radialDirection(m_SpaceDimension);
	vnl_vector<double> circDirection(m_SpaceDimension);
    
    // Compute the normals from the input surface mesh
    vtkSmartPointer<vtkPolyDataNormals> normalsFilter = vtkSmartPointer<vtkPolyDataNormals>::New();
    normalsFilter->SetInput(inputMesh);
    normalsFilter->ComputePointNormalsOn();
    normalsFilter->Update();
    
    vtkFloatArray* normals = 
    vtkFloatArray::SafeDownCast(normalsFilter->GetOutput()->GetPointData()->GetNormals());
    
	// Accessing input points
	vtkPoints* points = 0;
	if (!m_VolumetricMode)
	{
		points = inputMesh->GetPoints();
	}
	else
	{
		points = inputVolMesh->GetPoints();
	}
    
	// Formatting output: strain components + strain tensor + output points
	const unsigned int timeLength = m_NumberOfTimeSteps + 1;
	for (unsigned int index=0; index<timeLength; index++)
	{
		vtkSmartPointer<vtkFloatArray> s_array = vtkSmartPointer<vtkFloatArray>::New();
		s_array->SetName("longStrain");
		
		vtkSmartPointer<vtkFloatArray> rs_array = vtkSmartPointer<vtkFloatArray>::New();
		rs_array->SetName("radStrain");
        
		vtkSmartPointer<vtkFloatArray> cs_array = vtkSmartPointer<vtkFloatArray>::New();
		cs_array->SetName("circStrain");
        
		m_LongitudinalStrains.push_back(s_array);
		m_CircumferentialStrains.push_back(cs_array);
		m_RadialStrains.push_back(rs_array);
        
		vtkSmartPointer<vtkFloatArray> tensor_array = vtkSmartPointer<vtkFloatArray>::New();
		tensor_array->SetName("tensors");
		tensor_array->SetNumberOfComponents(9);
        
		m_FullStrains.push_back(tensor_array);
        
		vtkSmartPointer<vtkPoints> p = vtkSmartPointer<vtkPoints>::New();
		p->SetNumberOfPoints(points->GetNumberOfPoints());
		m_OutputPoints.push_back(p);
        
        
		if (m_Mask.IsNotNull()){
			vtkSmartPointer<vtkFloatArray> inside_array = vtkSmartPointer<vtkFloatArray>::New();
			inside_array->SetName("isinside");
			m_InsideArrays.push_back(inside_array);
		}
        
	}
    
	// Compute transformed coordinates and strain for each point in the input dataset
	for (int pid=0; pid < points->GetNumberOfPoints(); pid++)
	{
        
		double point[3];
        points->GetPoint(pid, point);
        
        
		TransformType::InputPointType inputPoint; 
		inputPoint.Fill(0.);
		TransformType::OutputPointType outputPoint; 
		outputPoint.Fill(0.);
        
		double normalToSurface[m_SpaceDimension];
        vtkIdType pid2;
        if (this->m_VolumetricMode)
        {
            // find the closest id in the surface mesh to this point
            pid2 = this->m_InputMesh->FindPoint(point);
        }
        else
        {
            pid2=pid;
        }
        
        if ( pid2<normals->GetMaxId() )
        {
            normals->GetTuple(pid2,normalToSurface);
            for (unsigned int d=0; d<m_SpaceDimension; d++)
            {
                radialDirection[d]=normalToSurface[d];
            }
        }
        else
        {
            for (unsigned int d=0; d<m_SpaceDimension; d++)
            {
                radialDirection[d]=0.;
            }
        }
            
		// Computation of local direction for strain 
		radialDirection.normalize();
        
		switch (this->m_AxisCompMethod)
		{
            case 0:
                
                circDirection[0] = -radialDirection[1];
                circDirection[1] = radialDirection[0];
                circDirection[2] = 0;
                circDirection.normalize();
                longDirection = vnl_cross_3d(circDirection,radialDirection);
                longDirection.normalize();
                break;
                
            case 1:
                circDirection = vnl_cross_3d(radialDirection, m_LongAxis); 
                circDirection.normalize();
                
                longDirection = vnl_cross_3d(circDirection,radialDirection); 
                longDirection.normalize();
                break;
                
            case 2:
                longDirection=m_LongAxis;
                radialDirection = radialDirection - dot_product(radialDirection,longDirection)*longDirection;
                radialDirection.normalize();
                circDirection = vnl_cross_3d(radialDirection, m_LongAxis); 
                circDirection.normalize();
                break;
		}
        
		vnl_matrix<double> longDirectionCopy(m_SpaceDimension,1);
		vnl_matrix<double> circDirectionCopy(m_SpaceDimension,1);
		vnl_matrix<double> radDirectionCopy(m_SpaceDimension,1);
        
		for (unsigned int d=0; d<m_SpaceDimension; d++)
		{
			longDirectionCopy(d,0)=longDirection[d];
			radDirectionCopy(d,0)=radialDirection[d];
			circDirectionCopy(d,0)=circDirection[d];
		}    
        
		// Applying transform to current point
		for (unsigned int d=0; d<m_SpaceDimension; d++)
		{
			inputPoint[d]=point[d];
		}
		inputPoint[m_SpaceDimension] = 0;
        
		// Physical Jacobian
		TransformType::VelocityType::ClassicalJacobianType physJacobian(m_SpaceDimension,m_SpaceDimension);
        
		TransformType::VelocityType::ClassicalJacobianType physJacobianAcc(m_SpaceDimension,m_SpaceDimension);
		physJacobianAcc.set_identity();
        
		TransformType::VelocityType::ClassicalJacobianType Ident(m_SpaceDimension,m_SpaceDimension);
		Ident.set_identity();
        
		double timeIncr = m_TimeSpacing;
		for (double time = m_TimeOrigin; time<= (double)(m_NumberOfTimeSteps-1); time+=timeIncr)
		{
			unsigned int timeIndex = (unsigned int)time;
            
			outputPoint = transfo->TransformPoint(inputPoint, time);
			transfo->GetVelocity()->GetClassicalJacobian(inputPoint, physJacobian,  m_TimeOrigin+timeIncr*m_TimeSpacing);
			physJacobianAcc =  physJacobian*physJacobianAcc;
            
			// Checking if point is inside the mask
			float isInside = 1.;
			if (m_Mask.IsNotNull())
			{
				MaskSpatialObjectType::PointType mPoint;
				mPoint[0] = outputPoint[0];
				mPoint[1] = outputPoint[1];
				mPoint[2] = outputPoint[2];
                
				if (! m_Mask->IsInside(mPoint))
					isInside = 0.;
                
				m_InsideArrays[timeIndex]->InsertTuple1(pid, isInside);
                
			}
            
			// Set transformed point 
			m_OutputPoints[timeIndex]->SetPoint(pid, outputPoint[0],
                                                outputPoint[1],
                                                outputPoint[2] );
            
			// Projecting strain on longitudinal direction
			vnl_matrix<double> strain(m_SpaceDimension, m_SpaceDimension);
			strain = (physJacobianAcc.transpose()*physJacobianAcc - Ident)/2.;
            
			// Strain for tensor
			vnl_matrix<double> strainForTensor(m_SpaceDimension, m_SpaceDimension);
			strainForTensor = physJacobianAcc.transpose()*physJacobianAcc;
            
			// Copy strain in ug
			m_FullStrains[timeIndex]->InsertTuple9(pid, strainForTensor[0][0], strainForTensor[0][1], strainForTensor[0][2],
                                                   strainForTensor[1][0], strainForTensor[1][1], strainForTensor[1][2],
                                                   strainForTensor[2][0], strainForTensor[2][1], strainForTensor[2][2]);
			
			// Longitudinal strain
			vnl_matrix<double> strainProj = longDirectionCopy.transpose() * strain  * longDirectionCopy;
			double strainValue = strainProj(0,0);
			m_LongitudinalStrains[timeIndex]->InsertTuple1(pid, strainValue);
            
			// Radial strain
			strainProj = radDirectionCopy.transpose() * strain  * radDirectionCopy;
			strainValue = strainProj(0,0);
			m_RadialStrains[timeIndex]->InsertTuple1(pid, strainValue);
			
			// Circ strain
			strainProj = circDirectionCopy.transpose() * strain  * circDirectionCopy;
			strainValue = strainProj(0,0);
			m_CircumferentialStrains[timeIndex]->InsertTuple1(pid, strainValue);
            
			// Getting ready for next iteration
			inputPoint = outputPoint;
            
		}// for (unsigned int time = 1; ... ; time++)
        
	}// for (int pid=0; ... ; pid++)
    
    
	// Configure output
	const double timeIncr = 1;
	for (double time = 0.; time<= (double)(m_NumberOfTimeSteps-1); time+=timeIncr)
	{
        
		unsigned int timeIndex = (unsigned int)time;
		if (!this->m_VolumetricMode)
        {
            // Generate the deformed polydata
            vtkSmartPointer<vtkPolyData> deformedPD = vtkSmartPointer<vtkPolyData>::New();
            deformedPD->DeepCopy( inputMesh );
            deformedPD->SetPoints( m_OutputPoints[timeIndex] );
            deformedPD->GetPointData()->AddArray( m_RadialStrains[timeIndex] );
            deformedPD->GetPointData()->AddArray( m_CircumferentialStrains[timeIndex] );
            deformedPD->GetPointData()->AddArray( m_LongitudinalStrains[timeIndex] );
            
            m_SmartOutput.push_back(deformedPD);
        }
		else
        {
            // Same with unstructured grid
            vtkSmartPointer<vtkUnstructuredGrid> ug = vtkSmartPointer<vtkUnstructuredGrid>::New();
            ug->DeepCopy(inputVolMesh);
            ug->SetPoints(m_OutputPoints[timeIndex]);
            ug->GetPointData()->AddArray( m_RadialStrains[timeIndex] );
            ug->GetPointData()->AddArray( m_CircumferentialStrains[timeIndex] );
            ug->GetPointData()->AddArray( m_LongitudinalStrains[timeIndex] );
            
            m_SmartOutput2.push_back(ug);
            
        }
        
	}
    
	return;
}


/**
 */
regDeformPolyDataFilterTDFFD::MeshCollectionType regDeformPolyDataFilterTDFFD::GetOutput()
{
	return m_SmartOutput;
}

/**
 */

regDeformPolyDataFilterTDFFD::UgCollectionType regDeformPolyDataFilterTDFFD::GetOutputUg()
{
	return m_SmartOutput2;
}

