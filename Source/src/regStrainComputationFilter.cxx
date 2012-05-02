// Copyright 2006 P	ompeu Fabra University (Computational Imaging Laboratory), Barcelona, Spain. Web: www.cilab.upf.edu.
// This software is distributed WITHOUT ANY WARRANTY;
// without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.

#include "regStrainComputationFilter.h"
#include <vnl/vnl_cross.h>



/**
*/
regStrainComputationFilter::regStrainComputationFilter()
{
	this->SetNumberOfRequiredInputs(1);
	SetNumberOfRequiredOutputs( 0 );
	m_NormalsFilter = vtkSmartPointer<vtkPolyDataNormals>::New();

	this->m_AxisCompMethod = 2;


	return;
}

/**
*/
void regStrainComputationFilter::SetInput( const MeshCollectionType* inputMesh )
{
	// Process object is not const-correct so the const_cast is required here
	this->ProcessObject::SetNthInput(0,  const_cast< MeshCollectionType* >(inputMesh) );

	return;
}
/**
This method effectively produces the ouput of the filter
- set of deformed meshes using the input 3D or 4D transformation
- strain computed on each node of the deformed meshes
*/
void regStrainComputationFilter::GenerateData()
{
	// In case we want to compute strain, delegate the computation to the ComputeStrain
	// function
	ComputeStrain();

	return;
}

void regStrainComputationFilter::ComputeStrain()
{
	MeshCollectionType* inputMeshWrapper = dynamic_cast<MeshCollectionType*>(this->ProcessObject::GetInput(0));
	vtkPolyData* inputMesh = 0;
	inputMesh = inputMeshWrapper->pds[ 0 ];

	// filter to compute normal directions
	// Computing CELL normals for input poly data
	m_NormalsFilter->SetInput(inputMesh);
	m_NormalsFilter->ComputeCellNormalsOn();
	// Avoid splitting cells
	m_NormalsFilter->SplittingOff( );
	m_NormalsFilter->Update();
	vtkCellData* normals2 = m_NormalsFilter->GetOutput()->GetCellData();

	// Specific preprocess
    m_CurrentTimeStep=0;
	if ((inputMesh->GetCellData()->GetArray("radialStrain") == NULL) || 
		(inputMesh->GetCellData()->GetArray("circStrain") == NULL) ||
		(inputMesh->GetCellData()->GetArray("longStrain") == NULL) ||
		(inputMesh->GetCellData()->GetVectors() == NULL) )
	{
	  PreProcessTimeStep( );
	}

	// From 1 to the end
	for (m_CurrentTimeStep=1; m_CurrentTimeStep < inputMeshWrapper->pds.size(); m_CurrentTimeStep++)
	{

		if ( inputMeshWrapper->pds[m_CurrentTimeStep]->GetNumberOfCells() != 
			 inputMesh->GetNumberOfCells() )
			break; // This should not be happening

		// Specific preprocess
		PreProcessTimeStep( );

		unsigned int normalsCounter = 0;
		for ( int cellId=0; cellId < m_NormalsFilter->GetOutput()->GetNumberOfCells(); cellId++)

		{
			vtkCell* cell = inputMeshWrapper->pds.at(m_CurrentTimeStep)->GetCell(cellId);

			// If the cell is triangular, process it. (It does not make sense to compute the normal
			// on non polygonal cells)
			if ( cell->GetCellType() == VTK_TRIANGLE ) 
			{
				// Getting normal direction of this cell
				double normalToSurface[m_ImageDimension];

				normals2->GetNormals()->GetTuple(normalsCounter, normalToSurface);

				vnl_vector<double>  radialDirection(m_ImageDimension,0.0);
				vnl_vector<double>  circDirection(m_ImageDimension,0.0);
				vnl_vector<double>  longDirection(m_ImageDimension,0.0);


                switch (this->m_AxisCompMethod)
                {
                    case 0:
                        
                        circDirection[0] = -normalToSurface[1];
                        circDirection[1] = normalToSurface[0];
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
        
				// Comptue strain
				StrainValues result;
				

				inputMeshWrapper->pds[m_CurrentTimeStep]->GetCellData()->GetVectors()
					->InsertTuple3(normalsCounter,	radialDirection[0],
					radialDirection[1],
					radialDirection[2]);

				ComputeStrainOnMesh( 
					inputMesh, 
					cellId, 
					cell, 
					radialDirection, 
					circDirection, 
					longDirection,
					result );

				vtkDataArray* currentArray = inputMeshWrapper->pds[m_CurrentTimeStep]->GetCellData()->GetArray("radialStrain");
				if (currentArray != NULL)
				{
				currentArray->InsertTuple1(normalsCounter, result.currentStrainValRad);
				}
				currentArray = inputMeshWrapper->pds[m_CurrentTimeStep]->GetCellData()->GetArray("circStrain");
				if (currentArray != NULL)
				{
				currentArray->InsertTuple1(normalsCounter, result.currentStrainValCirc);
				}
				currentArray = inputMeshWrapper->pds[m_CurrentTimeStep]->GetCellData()->GetArray("longStrain");
				if (currentArray != NULL)
				{
				currentArray->InsertTuple1(normalsCounter, result.currentStrainValLong);
				}
				

				// increment the counter of visited triangular cells
				normalsCounter++;

			} // END if ( cell->GetCellType() == VTK_TRIANGLE ) 			
		} //END for ( int cellId=0; cellId < m_NormalsFilter->GetOutput()->GetNumberOfCells(); cellId++)

	}
	// Normal exit
	return;
}

void regStrainComputationFilter::Update()
{
	MeshCollectionType* inputMeshWrapper = dynamic_cast<MeshCollectionType*>(this->ProcessObject::GetInput(0));

	// Make sure input is available
	if ( inputMeshWrapper == 0 )
	{
		itkExceptionMacro(<< "No input mesh");
	}

	// Notify start event observers
	this->InvokeEvent( itk::StartEvent() );

	// Actually do something
	this->GenerateData();

	// Notify end event observers
	this->InvokeEvent( itk::EndEvent() );

}

regStrainComputationFilter::MeshCollectionType* 
regStrainComputationFilter::GetInput()
{
	MeshCollectionType* inputMeshWrapper = 
		dynamic_cast<MeshCollectionType*>(this->ProcessObject::GetInput(0));
	return inputMeshWrapper;
}

void regStrainComputationFilter::ComputeStrainOnMesh(
	vtkPolyData* inputMesh,
	unsigned int cellId,
	vtkCell* cell,
	vnl_vector<double>  &radialDirection,
	vnl_vector<double>  &circumDirection,
	vnl_vector<double>  &longDirection,
	StrainValues &result )
{
	
	double point0[m_ImageDimension], point1[m_ImageDimension], point2[m_ImageDimension];
	double point0init[m_ImageDimension], point1init[m_ImageDimension], point2init[m_ImageDimension];
	double* displacements = new double[3*m_ImageDimension]; //Number of points in a triangle * ImageDimension
	double pcoord[m_ImageDimension];
	// Going through deformed meshes
	double* deriv = new double[m_ImageDimension*m_ImageDimension];

	cell->GetParametricCenter(pcoord);
	cell->GetPoints()->GetPoint(0, point0);
	cell->GetPoints()->GetPoint(1, point1);
	cell->GetPoints()->GetPoint(2, point2);

	unsigned int numberOfCells = inputMesh->GetNumberOfCells();

	if (cellId < numberOfCells)
	{
		inputMesh->GetCell(cellId)->GetPoints()->GetPoint(0, point0init);
		inputMesh->GetCell(cellId)->GetPoints()->GetPoint(1, point1init);
		inputMesh->GetCell(cellId)->GetPoints()->GetPoint(2, point2init);
	}

	for (unsigned int d=0; d<m_ImageDimension;d++)
		displacements[d] = point0[d] - point0init[d];
	for (unsigned int d=0; d<m_ImageDimension;d++)
		displacements[d+m_ImageDimension] = point1[d] - point1init[d];
	for (unsigned int d=0; d<m_ImageDimension;d++)
		displacements[d+2*m_ImageDimension] = point2[d] - point2init[d];

	cell->Derivatives(0, pcoord, displacements, 3, deriv);

	vnl_matrix<double> strain(deriv, m_ImageDimension, m_ImageDimension);
	vnl_matrix<double> identity(m_ImageDimension, m_ImageDimension);

	identity.set_identity();
    
	strain = (strain.transpose()*strain + strain + strain.transpose())/2.;

	vnl_matrix<double> radDirectionCopy(m_ImageDimension, 1);
	vnl_matrix<double> circDirectionCopy(m_ImageDimension, 1);
	vnl_matrix<double> longDirectionCopy(m_ImageDimension, 1);
	for (unsigned int d=0; d<m_ImageDimension; d++)
	{
		radDirectionCopy(d,0) = radialDirection[d];
		circDirectionCopy(d,0) = circumDirection[d];
		longDirectionCopy(d,0) = longDirection[d];
	}

	vnl_matrix<double> currentStrainRad 
		= radDirectionCopy.transpose() * strain * radDirectionCopy;
	vnl_matrix<double> currentStrainCirc 
		= circDirectionCopy.transpose() * strain * circDirectionCopy;
	vnl_matrix<double> currentStrainLong
		= longDirectionCopy.transpose() * strain * longDirectionCopy;

	result.currentStrainValRad = currentStrainRad(0,0);
	result.currentStrainValCirc = currentStrainCirc(0,0);
	result.currentStrainValLong = currentStrainLong(0,0);

	delete [] displacements;
	delete [] deriv;
}

void regStrainComputationFilter::PreProcessTimeStep()
{
	MeshCollectionType* inputMeshWrapper = dynamic_cast<MeshCollectionType*>(this->ProcessObject::GetInput(0));


	//Add the scalar strain values to the initial mesh. All values are 0
	vtkSmartPointer<vtkFloatArray> rs_array = vtkSmartPointer<vtkFloatArray>::New();
	rs_array->SetName("radialStrain");
	vtkSmartPointer<vtkFloatArray> cs_array = vtkSmartPointer<vtkFloatArray>::New();
	cs_array->SetName("circStrain");
	vtkSmartPointer<vtkFloatArray> ls_array = vtkSmartPointer<vtkFloatArray>::New();
	ls_array->SetName("longStrain");
	vtkSmartPointer<vtkFloatArray> vectors = vtkSmartPointer<vtkFloatArray>::New();
	vectors->SetNumberOfComponents(3);
	vectors->SetName("vectors");
	if (inputMeshWrapper->pds[ m_CurrentTimeStep ]->GetCellData()->GetArray("radialStrain") == NULL)
	  inputMeshWrapper->pds[ m_CurrentTimeStep ]->GetCellData()->AddArray(rs_array);
    if (inputMeshWrapper->pds[ m_CurrentTimeStep ]->GetCellData()->GetArray("longStrain") == NULL)
	  inputMeshWrapper->pds[ m_CurrentTimeStep ]->GetCellData()->AddArray(ls_array);
	if (inputMeshWrapper->pds[ m_CurrentTimeStep ]->GetCellData()->GetArray("circStrain") == NULL)
	  inputMeshWrapper->pds[ m_CurrentTimeStep ]->GetCellData()->AddArray(cs_array);
	if (inputMeshWrapper->pds[ m_CurrentTimeStep ]->GetCellData()->GetVectors() == NULL)
	  inputMeshWrapper->pds[ m_CurrentTimeStep ]->GetCellData()->SetVectors(vectors);

	for(int i = 0; i<inputMeshWrapper->pds[ m_CurrentTimeStep ]->GetNumberOfCells(); i++ )
	{	
		vtkDataArray* currentArray = inputMeshWrapper->pds[ m_CurrentTimeStep ]->GetCellData()->GetArray("radialStrain");
		if (currentArray != NULL)
		{
			currentArray->InsertTuple1(i,0.);
		}
		currentArray = inputMeshWrapper->pds[ m_CurrentTimeStep ]->GetCellData()->GetArray("circStrain");
		if (currentArray != NULL)
		{
			currentArray->InsertTuple1(i, 0.);
		}
		currentArray = inputMeshWrapper->pds[ m_CurrentTimeStep ]->GetCellData()->GetArray("longStrain");
		if (currentArray != NULL)
		{
			currentArray->InsertTuple1(i, 0.);
		}
		
		currentArray = inputMeshWrapper->pds[ m_CurrentTimeStep ]->GetCellData()->GetVectors();
		if (currentArray != NULL)
		{
			currentArray->InsertTuple3(i, 0.,0.,0.);
		}
	}
}

