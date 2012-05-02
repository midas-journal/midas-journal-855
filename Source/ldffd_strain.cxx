#include <vtkPolyData.h>
#include <vtkPolyDataReader.h>
#include <vtkPolyDataWriter.h>
#include <regStrainComputationFilter.h>
#include "itkTransformFileReader.h"
#include "itkTransformFactory.h"
#include "itkDiffeomorphicBSplineTransform.h"

// Our temporal diff matching stuff
#include "regDiffeomorphicContinuousBSplineTransform.h"

const    unsigned int SpaceDimension = 3; 
typedef itk::DiffeomorphicBSplineTransform<double, SpaceDimension, 3>  TransformType;


// Function that deforms a polydata according to a LDFFD transform
std::vector<vtkSmartPointer<vtkPolyData> > deformPD(vtkPolyData* input, TransformType* tr);

int main( int argc, char * argv[] )
{
    
    const unsigned int FN_MAX_LENGTH = 1024;
    
    // Parsing arguments
    if (argc < 4)
    {
        std::cout << "Usage: " << argv[0] 
        << " inputPolyData.vtk inputTransform.dof outputMeshPrefix [-la 0 0 1]"
        << std::endl;
        
        return EXIT_FAILURE;
    }
    
    
    const char* inputPDFileName = argv[1];
    const char* inputTransformFileName = argv[2];   
    const char* outputMeshPrefix = argv[3];
    
    
    vnl_vector<double> longAxis(SpaceDimension);
    longAxis[0]=0.; longAxis[1]=0.; longAxis[2]=1.;
    
    int arg=3;
    while (arg<argc) {
        
        if (strcmp(argv[arg],"-la")==0)
        {
            longAxis[0]=atof(argv[arg+1]);
            longAxis[1]=atof(argv[arg+2]);
            longAxis[2]=atof(argv[arg+3]);
            arg += 4;
        }
        else
        {
            std::cout << "ERROR PARSING ARGUMENTS" << std::endl;
            return EXIT_FAILURE;
        }
    }
    
    // Reading input transform
    // Read transform from file
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
    TransformType::Pointer transform =
    static_cast<TransformType*>((*it).GetPointer());

    
    // Reading input pd
    // Load the input polydata
    vtkSmartPointer<vtkPolyDataReader> reader = vtkSmartPointer<vtkPolyDataReader>::New();
    reader->SetFileName(inputPDFileName);
    reader->Update();
    vtkPolyData* inputPD = reader->GetOutput();
    
    // Generate a vector of pds
    std::vector<vtkSmartPointer<vtkPolyData> > pds = deformPD(inputPD, transform);
    
    // input meshes for strain filter 
    MeshContainer::Pointer cont = MeshContainer::New();
    for (int index=0; index<pds.size(); index++)
    {
        cont->AddInput(pds[index]);
    }
    
    // strain filter
    typedef regStrainComputationFilter FilterType;
    FilterType::Pointer filter=FilterType::New();
    filter->SetLongAxis(longAxis);
    filter->SetInput(cont);
    filter->Update();
    
    
    // Writing output meshes
    char outputFileName [FN_MAX_LENGTH];
   
        vtkSmartPointer<vtkPolyDataWriter> writer = vtkSmartPointer<vtkPolyDataWriter>::New();
        
        for (unsigned int ts=0; ts<cont->pds.size(); ts++)
        {
            sprintf(outputFileName, "%s%03d.vtk",outputMeshPrefix,ts);
            writer->SetFileName(outputFileName);
            writer->SetInput(cont->pds[ts]);
            writer->Update();
        }
   
    return EXIT_SUCCESS;    
}


std::vector<vtkSmartPointer<vtkPolyData> > deformPD(vtkPolyData* input, TransformType* tr)
{
    std::vector<vtkSmartPointer<vtkPolyData> > output(0);
    double point[SpaceDimension];
    TransformType::InputPointType point2;
    TransformType::OutputPointType point3;
    
    for (unsigned int ts=0; ts<tr->GetNumberOfTimeSteps(); tr++) {
        
        
        vtkSmartPointer<vtkPolyData> pd = vtkSmartPointer<vtkPolyData>::New();
        pd->DeepCopy(input);
        
        for (int pid=0; pid<input->GetNumberOfPoints(); pid++) {
            
            input->GetPoints()->GetPoint(pid, point);
            point2=point;
            point3=tr->TransformPoint(point2, ts);
            pd->GetPoints()->SetPoint(pid, point3.GetVnlVector().data_block());
        }
        
        output.push_back(pd);
    }
    
    
    return output;

}


