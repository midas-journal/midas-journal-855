#include <itkTransformFactory.h>
#include <itkTransformFileWriter.h>
#include <itkTransformFileReader.h>
#include <vtkPolyData.h>
#include <vtkPolyDataReader.h>
#include <vtkUnstructuredGrid.h>
#include <vtkUnstructuredGridReader.h>
#include <vtkPolyDataWriter.h>
#include <regDeformPolyDataFilterTDFFD.h>

// Our temporal diff matching stuff
#include "regDiffeomorphicContinuousBSplineTransform.h"


int main( int argc, char * argv[] )
{
    
    const unsigned int FN_MAX_LENGTH = 1024;
    
    // Parsing arguments
    if (argc < 5)
    {
        std::cout << "Usage: " << argv[0] 
        << " inputPolyData.vtk inputTransform.dof numberOfTimeSteps outputMeshPrefix [-la 0 0 1]"
        << std::endl;
        
        return EXIT_FAILURE;
    }
    
    const    unsigned int SpaceDimension = 3;
    const    unsigned int Dimension = SpaceDimension+1;
    const    unsigned int BSplineOrder = 3;
    typedef itk::DiffeomorphicContinuousBSplineTransform<double,Dimension, BSplineOrder>  TransformType;
    
    const char* inputPDFileName = argv[1];
    const char* inputTransformFileName = argv[2];
    const unsigned int numberOfTimeSteps = atoi(argv[3]);
    const char* outputMeshPrefix = argv[4];
    
    vnl_vector<double> longAxis(SpaceDimension);
    longAxis[0]=0.; longAxis[1]=0.; longAxis[2]=1.; 
    char* ugFileName = 0;
    
    int arg=5;
    while (arg<argc) {
        
        if (strcmp(argv[arg],"-la")==0)
        {
            longAxis[0]=atof(argv[arg+1]);
            longAxis[1]=atof(argv[arg+2]);
            longAxis[2]=atof(argv[arg+3]);
            arg += 4;
        }
        else if (strcmp(argv[arg],"-ug")==0)
        {
            ugFileName = argv[arg+1];
            arg+=2;
        }
        else
        {
            std::cout << "ERROR PARSING ARGUMENTS" << std::endl;
            return EXIT_FAILURE;
        }
    }

    
    // Load the input polydata
    vtkSmartPointer<vtkPolyDataReader> reader = vtkSmartPointer<vtkPolyDataReader>::New();
    reader->SetFileName(inputPDFileName);
    reader->Update();
    vtkPolyData* inputPD = reader->GetOutput();
    
    // If provided load input ug
    vtkSmartPointer<vtkUnstructuredGridReader> readerug = vtkSmartPointer<vtkUnstructuredGridReader>::New();
    vtkUnstructuredGrid* inputUG = 0;
    if (ugFileName != 0)
    {
        readerug->SetFileName(ugFileName);
        readerug->Update();
        inputUG = readerug->GetOutput();
    }
    
    
    // Load the input transformation
    TransformType::Pointer transform=0;
    itk::TransformFactory<TransformType>::RegisterTransform();
    itk::TransformFileReader::Pointer tr_reader;
    tr_reader = itk::TransformFileReader::New();
    tr_reader->SetFileName( inputTransformFileName );
    
    // Actual reading
    try
    {
        tr_reader->Update();
    }
    catch( itk::ExceptionObject & excp )
    {
        std::cerr << "Error while reading the transform file" << std::endl;
        std::cerr << excp << std::endl;
        return EXIT_FAILURE;
    }
    
    // Accessing the output
    typedef itk::TransformFileReader::TransformListType * TransformListType;
    TransformListType transforms = tr_reader->GetTransformList();
    std::cout << "Number of transforms = " << transforms->size() << std::endl;
    itk::TransformFileReader::TransformListType::const_iterator it 
    = transforms->begin();
    transform = 
    static_cast<TransformType*>((*it).GetPointer());
    
    // strain filter
    typedef regDeformPolyDataFilterTDFFD FilterType;
    FilterType::Pointer filter=FilterType::New();
    filter->SetInputTransform(transform);
    filter->SetInputMesh(inputPD);
    filter->SetNumberOfTimeSteps(numberOfTimeSteps);
    filter->SetLongAxis(longAxis);
    if (inputUG != 0) {
        filter->SetInputVolMesh(inputUG);
    }
    filter->Update();//GenerateData();
    
    // Writing output meshes
    char outputFileName [FN_MAX_LENGTH];
    if (!filter->GetVolumetricMode())
    {
        vtkSmartPointer<vtkPolyDataWriter> writer = vtkSmartPointer<vtkPolyDataWriter>::New();
        
        for (unsigned int ts=0; ts<numberOfTimeSteps; ts++)
        {
            sprintf(outputFileName, "%s%03d.vtk",outputMeshPrefix,ts);
            writer->SetFileName(outputFileName);
            writer->SetInput(filter->GetOutput()[ts]);
            writer->Update();
        }
    }
    else
    {
        vtkSmartPointer<vtkUnstructuredGridWriter> writer = vtkSmartPointer<vtkUnstructuredGridWriter>::New();
        
        for (unsigned int ts=0; ts<numberOfTimeSteps; ts++)
        {
            sprintf(outputFileName, "%s%03d.vtk",outputMeshPrefix,ts);
            writer->SetFileName(outputFileName);
            writer->SetInput(filter->GetOutputUg()[ts]);
            writer->Update();
        }
    }
    
    return EXIT_SUCCESS;    
}
