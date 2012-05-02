// Image types and other standard things
#include <itkImage.h>
#include <itkImageFileReader.h>
#include <itkImageFileWriter.h>
#include <itkImageRegionIterator.h>
#include <itkLinearInterpolateImageFunction.h>
#include <itkBSplineInterpolateImageFunction.h>
#include <itkWarpImageFilter.h>
#include <itkTransformFactory.h>
#include <itkTransformFileWriter.h>
#include <itkTransformFileReader.h>
#include <itkImageMaskSpatialObject.h>
#include <itkTransformFactory.h>
#include <vtkPolyDataReader.h>

// Our temporal diff matching stuff
#include "regDiffeomorphicContinuousBSplineTransform.h"
#include "regTemporalImageMetric.h"
#include "regMeanSquaresTemporalImageMetricSequential2.h"
#include "regMeanSquaresTemporalImageMetric.h"
#include "regMeanSquaresTemporalImageMetricDense.h"
#include "regHybridSeqNonSeqMetric.h" 
#include "regTransformWrapper.h"
#include "regBSplineFieldRefine.h"
#include "vnl/vnl_math.h"

// Optimizer
#include <itkLBFGSBOptimizer.h>

// Iterators
#include <itkImageRegionConstIteratorWithIndex.h>

// Generate displacement field and writing the transformation
#include "itkTransformToDeformationFieldFilter.h"
#include <itkTransformFileWriter.h>
#include "itkTimeProbe.h"


#include <itkDanielssonDistanceMapImageFilter.h>

#define __UseTimeProbes 

class MyGlobalTypesAndParameters
{
public:
    static const unsigned int    SpaceDimension = 3;
    static const unsigned int    Dimension = SpaceDimension + 1;
    static const unsigned int    BSplineOrder = 3;
    
    typedef itk::DiffeomorphicContinuousBSplineTransform<double,Dimension, BSplineOrder>  TransformType;	
    typedef itk::LBFGSBOptimizer OptimizerType;
    typedef itk::ImageRegion<Dimension> RegionType;
    typedef itk::ImageMaskSpatialObject< SpaceDimension > MaskSpatialObjectType;
    // Main image type
    typedef  short  PixelType;
    typedef itk::Image< PixelType, Dimension >  ImageType;
    
    // Global parameters
    char* inputImage;
    char* outputImageFilePrefix;
    char* outputTransformFile;
    char* outputDispFileName;
    char* inputMaskFileName;
    char* inputMovingMaskFileName;
    
    char* inputTransformationFileName;
    char* incompressibilityPointSetFileName;
    double incompressibilityWeight;
    
    unsigned int gridSize[Dimension];
    bool splitTransformation;
    float maximumVelocityInsideMaskX;
    float maximumVelocityInsideMaskY;
    float maximumVelocityInsideMaskZ;
    bool metricTypeSequential;
    bool metricTypeDense;
    double referenceTime;
    double minTimeStep;
    double timeStep;
    RegionType imageRegion;
    unsigned int numberOfSamples;
    unsigned int padding[Dimension];
    bool timeMultiRes;
    float seqWeight;
    
    MyGlobalTypesAndParameters(){
        // Default parameters
        inputImage = 0;
        outputImageFilePrefix = 0;
        outputTransformFile = 0;
        outputDispFileName = 0;
        inputMaskFileName = 0;
        inputMovingMaskFileName = 0;
        inputTransformationFileName = 0;
        incompressibilityPointSetFileName = 0;
        maximumVelocityInsideMaskX=10.f;
        maximumVelocityInsideMaskY=10.f;
        maximumVelocityInsideMaskZ=10.f;
        splitTransformation=false;
        metricTypeSequential=false;
        metricTypeDense=false;
        referenceTime=0.;
        timeStep=2.0;
        minTimeStep=10.0;
        numberOfSamples=2000;
        for (unsigned int d=0; d<Dimension; d++)
        {
            padding[d] = 0;
        }
        incompressibilityWeight=0.;
        seqWeight=0.;
        
        for (unsigned int d=0; d<Dimension; d++)
        {
            gridSize[d] = 0;
        }
        
        RegionType::SizeType size;
        size.Fill(0);
        imageRegion.SetSize(size);
        RegionType::IndexType index;
        index.Fill(0);
        imageRegion.SetIndex(index);
        
        timeMultiRes=false;
        
    };
    
    ~MyGlobalTypesAndParameters(){};
    
};

// Parse command line
bool ParseCommandLine(int argc, char * argv[], MyGlobalTypesAndParameters& params);

void usage(){
    
    std::cout << "Usage : inputImage outputImageFilePrefix outputTransformFile" << std::endl;
    std::cout << "Optional arguments : " << std::endl;
    std::cout << " -inputTransformation transformFile" << std::endl;
    std::cout << " -splitTransformation" << std::endl;
    std::cout << " -timeMultiRes" << std::endl;
    std::cout << " -regionSize SizeX SizeY SizeZ SizeT" << std::endl;
    std::cout << " -inputMask MaskFileName" << std::endl;
    std::cout << " -inputMovingMask MaskMovingFileName" << std::endl;
    std::cout << " -imageRegion indexX indexY indexZ indexT sizeX sizeY sizeZ sizeT" << std::endl;
    std::cout << " -MetricSequential" << std::endl;
    std::cout << " -MetricDense" << std::endl;
    std::cout << " -weightSequential weight" << std::endl;
    std::cout << " -referenceTime refTime" << std::endl;
    std::cout << " -numSamples numberOfSamples" << std::endl;
    std::cout << " -maxv maxVelocityValueInsideMask" << std::endl;
    std::cout << " -maxvx maxVelocityValueInsideMaskX" << std::endl;
    std::cout << " -maxvy maxVelocityValueInsideMaskY" << std::endl;
    std::cout << " -maxvz maxVelocityValueInsideMaskZ" << std::endl;
    std::cout << " -timeStep timeStep" << std::endl;
    std::cout << " -minTimeStep minTimeStep" << std::endl;
    std::cout << " -padding numberOfVoxels" << std::endl;
    std::cout << " -paddingx numberOfVoxelsInX" << std::endl;
    std::cout << " -paddingy numberOfVoxelsInY" << std::endl;
    std::cout << " -paddingz numberOfVoxelsInZ" << std::endl;
    std::cout << " -dispFileName outputDispFilePrefix" << std::endl;
    
}

int main( int argc, char * argv[] )
{
    
    MyGlobalTypesAndParameters params;
    
    if ( (argc < 4)||(!ParseCommandLine(argc, argv, params)) )
    {
        std::cerr << "Problem encountered while parsing arguments. " << std::endl;
        usage();   
        exit(EXIT_FAILURE);
    }
    
    const    unsigned int SpaceDimension = MyGlobalTypesAndParameters::SpaceDimension;
    const    unsigned int Dimension = MyGlobalTypesAndParameters::Dimension;
    const    unsigned int BSplineOrder = MyGlobalTypesAndParameters::BSplineOrder;
    typedef  MyGlobalTypesAndParameters::MaskSpatialObjectType MaskSpatialObjectType;
    
    typedef MyGlobalTypesAndParameters::PixelType PixelType;
    typedef MyGlobalTypesAndParameters::ImageType ImageType;
    typedef MyGlobalTypesAndParameters::TransformType TransformType;
    
    typedef MaskSpatialObjectType::ImageType MaskImageType;
    
    // Reading input image
    typedef itk::ImageFileReader<ImageType> ReaderType;
    ReaderType::Pointer reader =   ReaderType::New();
    reader->SetFileName(params.inputImage);
    
    try
    {
        reader->Update();
    }
    catch( itk::ExceptionObject & excp )
    {
        std::cerr << "Error while reading the transform file" << std::endl;
        std::cerr << excp << std::endl;
        return EXIT_FAILURE;
    }
    
    ImageType::Pointer image = reader->GetOutput();
    
    // Reading input mask
    MaskSpatialObjectType::Pointer mask = 0;
    MaskSpatialObjectType::Pointer movingMask = 0;
    
    if (params.inputMaskFileName)
    {
        typedef itk::ImageFileReader<MaskImageType> MaskReaderType;
        MaskReaderType::Pointer maskReader =   MaskReaderType::New();
        maskReader->SetFileName(params.inputMaskFileName);
        
        try
        {
            maskReader->Update();
        }
        catch( itk::ExceptionObject & excp )
        {
            std::cerr << "Error while reading the mask file" << std::endl;
            std::cerr << excp << std::endl;
            return EXIT_FAILURE;
        }
        
        MaskImageType::Pointer maskImage = maskReader->GetOutput();
        mask = MaskSpatialObjectType::New();
        mask->SetImage(maskImage);
        mask->Initialize();
    }
    
    if (params.inputMovingMaskFileName)
    {
        typedef itk::ImageFileReader<MaskImageType> MaskReaderType;
        MaskReaderType::Pointer maskReader =   MaskReaderType::New();
        maskReader->SetFileName(params.inputMovingMaskFileName);
        std::cout << "Reading Moving image mask" << std::endl;
        try
        {
            maskReader->Update();
        }
        catch( itk::ExceptionObject & excp )
        {
            std::cerr << "Error while reading the mask file" << std::endl;
            std::cerr << excp << std::endl;
            return EXIT_FAILURE;
        }
        
        MaskImageType::Pointer maskImage = maskReader->GetOutput();
        
        movingMask = MaskSpatialObjectType::New();
        movingMask->SetImage(maskImage);
        movingMask->Initialize();
    }
    
    
    // Metric
    typedef itk::TemporalImageMetric<ImageType> GenericMetricType;
    typedef itk::MeanSquaresTemporalImageMetricSequential2<ImageType> MetricTypeSeq;
    typedef itk::MeanSquaresTemporalImageMetricDense<ImageType> MetricTypeDense;
    typedef itk::MeanSquaresTemporalImageMetric<ImageType> MetricTypeFirst;
    typedef itk::HybridSeqNonSeqMetric<ImageType> MetricTypeHybrid;
    
    typedef itk::BSplineInterpolateImageFunction< ImageType,
    double > InterpolatorType;
    InterpolatorType::Pointer interpolator = InterpolatorType::New();
    
    GenericMetricType::Pointer metric = 0;
    MetricTypeSeq::Pointer metricSeq = 0;
    MetricTypeFirst::Pointer metricFirst = 0;
    MetricTypeDense::Pointer metricDense = 0;
    MetricTypeHybrid::Pointer metricMutante = 0;
    
    if (params.metricTypeSequential==true) 
    {
        std::cout << "Hybrid ToFirst + Sequential metric" << std::endl;
        metricFirst = MetricTypeFirst::New();
        metricSeq = MetricTypeSeq::New();
        metricMutante = MetricTypeHybrid::New();
        metricMutante->SetMetricSequential(metricSeq);
        metricMutante->SetMetricNonSequential(metricFirst);
        metric = metricMutante;
    }
    else if (params.metricTypeDense==true)
    {
        std::cout << "Hybrid Dense + Sequential Metric" << std::endl;
        metricDense = MetricTypeDense::New();
        metricSeq = MetricTypeSeq::New();
        metricMutante = MetricTypeHybrid::New();
        metricMutante->SetMetricSequential(metricSeq);
        metricMutante->SetMetricNonSequential(metricDense);
        
        metric = metricMutante;
    }
    else
    {
        std::cout << "Classical Non sequential metric" << std::endl;
        metricFirst = MetricTypeFirst::New();
        metric = metricFirst;
    }
    
    metric->SetImage( image );
    
    
    if (params.imageRegion.GetSize()[0] == 0){
        
        
        if (mask.IsNotNull()){
            MaskSpatialObjectType::ImageType::RegionType tmpRegion =
            mask->GetAxisAlignedBoundingBoxRegion();
            
            std::cout << tmpRegion << std::endl;
            
            ImageType::RegionType::SizeType regionSize;
            for (unsigned int d=0; d<SpaceDimension; d++) {
                regionSize[d] = tmpRegion.GetSize()[d];
            }
            regionSize[SpaceDimension] = image->GetBufferedRegion().GetSize()[SpaceDimension];
            
            ImageType::RegionType::IndexType index;
            for (unsigned int d=0; d<SpaceDimension; d++) {
                index[d] = tmpRegion.GetIndex()[d];
            }
            index[SpaceDimension]=0;
            
            ImageType::RegionType region;
            region.SetIndex(index);
            region.SetSize(regionSize);
            
            std::cout<<"Setting image region to : "<<region<<std::endl;
            
            metric->SetImageRegion(region);
        }else {
            metric->SetImageRegion( image->GetBufferedRegion() );
        }
        
    }
    else
    {
        std::cout << "Setting image region to " << params.imageRegion << std::endl;
        metric->SetImageRegion( params.imageRegion );
    }
    metric->SetInterpolator( interpolator );
    std::cout << "Setting number of samples to " << params.numberOfSamples << std::endl;
    metric->SetNumberOfSamples(params.numberOfSamples);
    
    if (mask.IsNotNull())
    {
        metric->SetSpaceImageMask(mask);
    }
    
    if (movingMask.IsNotNull())
    {
        metric->SetSpaceMovingImageMask(movingMask);
    }
    
    if (params.referenceTime != 0.)
    {
        double referenceTime = params.referenceTime;
        metric->SetReferenceTime(referenceTime);
        std::cout << "Setting refererence time to : " << referenceTime << std::endl;
    }
    
    ///////// Configure Transform ////////
    TransformType::Pointer transform;
    
    if (params.inputTransformationFileName)
    {
		// Load the input transformation
		// Read transform from file
		itk::TransformFactory<TransformType>::RegisterTransform();
		itk::TransformFileReader::Pointer tr_reader;
		tr_reader = itk::TransformFileReader::New();
		tr_reader->SetFileName( params.inputTransformationFileName );
        
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
        
		// Check if we need to refine input transformation by a factor 2
		if (params.splitTransformation)
		{
			typedef itk::BSplineFieldRefine<double, Dimension, BSplineOrder> RefineType;
			RefineType::Pointer refine = RefineType::New();
			refine->SetInput(transform->GetVelocity());
			RefineType::SpacingType spacingBefore = transform->GetGridSpacing();
			RefineType::OriginType originBefore = transform->GetGridOrigin();
			RefineType::SizeType sizeBefore = transform->GetGridRegion().GetSize();
            
			// Double size after removing 3 control points taking into account the borders
			RefineType::SizeType sizeOnImageBefore;
			
            unsigned int maxDim;
            
            if (params.timeMultiRes)
            {
                maxDim=Dimension;
            }
            else
            {
                maxDim=SpaceDimension;
            }
            
            for (unsigned int d=0; d<maxDim; d++)
			{
				if (sizeBefore[d] >= 3)
				{
					sizeOnImageBefore[d] = sizeBefore[d] - 3;
				}
				else
				{
					sizeOnImageBefore[d] = sizeBefore[d];
				}
			}
			RefineType::SizeType sizeAfter = sizeBefore;
			RefineType::OriginType originAfter = originBefore;
			RefineType::SpacingType spacingAfter = spacingBefore;
			
            
            for (unsigned int d=0; d<maxDim; d++)
			{
				sizeAfter[d] = 2 * sizeOnImageBefore[d] + 3;
				spacingAfter[d] = spacingBefore[d] / 2.f;
				originAfter[d] = originBefore[d] + spacingBefore[d] - spacingAfter[d];
			}
			RefineType::RegionType regionAfter;
			regionAfter.SetSize(sizeAfter);
			
			refine->ConfigureOutput(spacingAfter, originAfter, regionAfter);
			refine->Update();
			transform->SetVelocity(refine->GetOutput());
			transform->SetParametersByValue(refine->GetOutput()->GetParameters());
            
		}// Split
        
        /////////// Debug: Save initial transformation
        itk::TransformFileWriter::Pointer trWriterTemp = itk::TransformFileWriter::New();
        trWriterTemp->SetInput(transform);
        trWriterTemp->SetFileName("initialTransform.dof");
        try
        {  
            trWriterTemp->Update();
        }
        catch( itk::ExceptionObject & excp )
        {
            std::cerr << "Error when saving intial transformation" << std::endl;
            std::cerr << excp << std::endl;
            return EXIT_FAILURE;
        }
        
        
    }// Input transform specified
    else
    {
        
        // Specifying transformation region, origin and spacing
        typedef TransformType::RegionType RegionType;
        RegionType bsplineRegion;
        TransformType::OriginType originBSpline;
        TransformType::SpacingType spacingBSpline;
        
        RegionType::SizeType   gridSizeOnImage;
        RegionType::SizeType   gridBorderSize;
        RegionType::SizeType   totalGridSize;
        
        if (params.gridSize[0] == 0)
        {
            gridSizeOnImage.Fill(5); // Default value
        }
        else
        {
            for (unsigned int d=0; d<Dimension; d++)
            {
                gridSizeOnImage[d] = params.gridSize[d];
            }
        }
        
        gridBorderSize.Fill( 3 );    // Border for spline order = 3 ( 1 lower, 2 upper )
        totalGridSize = gridSizeOnImage + gridBorderSize;
        bsplineRegion.SetSize( totalGridSize );
        
        for(unsigned int r=0; r<Dimension; r++)
        {
            
            spacingBSpline[r] = image->GetSpacing()[r]
            * static_cast<double>(metric->GetImageRegion().GetSize()[r] - 1 + 2 * params.padding[r]) 
            / static_cast<double>(gridSizeOnImage[r] - 1) ;
            
            originBSpline[r]  =  image->GetOrigin()[r] + image->GetSpacing()[r] * static_cast<double> (metric->GetImageRegion().GetIndex()[r]) 
            -  image->GetSpacing()[r] * static_cast<double> (params.padding[r])
            - spacingBSpline[r];
        }
        
        
        // Actually creating the transformation
        transform = TransformType::New();
        transform->SetGridRegion(bsplineRegion);
        transform->SetGridSpacing(spacingBSpline);
        transform->SetGridOrigin(originBSpline);
        // parameters
        TransformType::ParametersType parameters(transform->GetNumberOfParameters());
        parameters.Fill(0.);
        transform->SetParametersByValue( parameters );
        
    }
    
    transform->SetTimeStep(image->GetSpacing()[SpaceDimension] / params.timeStep);
    transform->SetMinimumTimeStep(image->GetSpacing()[SpaceDimension] / params.minTimeStep);
    std::cout <<"Time step: "<< 1/params.timeStep<<std::endl;
    
    ///////// End: Configure Transform /////////
    
    // Metric initialization
    //
    metric->SetTransform(transform);
    metric->Initialize();
    
    // Optimizer
    typedef MyGlobalTypesAndParameters::OptimizerType OptimizerType;
    OptimizerType::Pointer optim = OptimizerType::New();
    optim->SetInitialPosition(transform->GetParameters());
    
    // Bound selection 
    OptimizerType::BoundSelectionType selection(transform->GetNumberOfParameters());
    selection.Fill(2); // This means we are going to specify lower and upper bounds on all parameters
    
    // Get the indexes of parameters that are on the borders 
    OptimizerType::BoundValueType jamesBoundsLow(transform->GetNumberOfParameters());
    OptimizerType::BoundValueType jamesBoundsHigh(transform->GetNumberOfParameters());
    
    // Setting bound
    for (unsigned long index=0; index<transform->GetNumberOfParameters(); index++)
    {
        unsigned int dim=(unsigned int)(index%transform->GetNumberOfParametersPerDimension());
        if (dim==0)
        {
            jamesBoundsHigh[index]=params.maximumVelocityInsideMaskX;
            jamesBoundsLow[index]=-jamesBoundsHigh[index];
        }else if (dim==1) {
            jamesBoundsHigh[index]=params.maximumVelocityInsideMaskY;
            jamesBoundsLow[index]=-jamesBoundsHigh[index];
        }else {
            jamesBoundsHigh[index]=params.maximumVelocityInsideMaskZ;
            jamesBoundsLow[index]=-jamesBoundsHigh[index];
        }
    }
    
    optim->SetBoundSelection(selection);
    optim->SetMaximumNumberOfIterations(100);
    optim->MinimizeOn();
    optim->SetCostFunction(metric);
    
    for (unsigned int i=0; i<transform->GetParameters().GetSize(); i++) {
        jamesBoundsLow[i]  = jamesBoundsLow[i] + transform->GetParameters()[i];
        jamesBoundsHigh[i] = jamesBoundsHigh[i] + transform->GetParameters()[i];
        
    }
    
    optim->SetLowerBound(jamesBoundsLow);
    optim->SetUpperBound(jamesBoundsHigh);
    
    // In case we use an hybrid metric, adjust the weight automatically
    
    const TransformType::ParametersType parametersCopy = transform->GetParameters();
    if (metricMutante.IsNotNull())
    {
        // In case we use an hybrid metric, we need to adjust the weight between
        // sequential and non-sequential terms
        if (params.seqWeight==0.)
        {
            // FIXME Don't understand why I cannot pass directly transform->GetParameters()
            metricMutante->ComputeWeight(parametersCopy);
        }
        else
        {
            metricMutante->SetWeight(params.seqWeight);
        }
    }
    
    try
    {
        optim->StartOptimization();
    }
    catch( itk::ExceptionObject & excp )
    {
        std::cerr << "Error during optimization" << std::endl;
        std::cerr << excp << std::endl;
        return EXIT_FAILURE;
    }
    
    // Recover last transformation parameters
    transform->SetParameters(optim->GetCurrentPosition());
    
    // Writing what we got
    itk::TransformFileWriter::Pointer trWriter = itk::TransformFileWriter::New();
    trWriter->SetInput(transform);
    trWriter->SetFileName(params.outputTransformFile);
    try
    {  
        trWriter->Update();
    }
    catch( itk::ExceptionObject & excp )
    {
        std::cerr << "Error during optimization" << std::endl;
        std::cerr << excp << std::endl;
        return EXIT_FAILURE;
    }
    
    // If by miracle we reach this line ... We do a resampling of the input image
    // using our tricky transform wrapper to get a sequence of 3D images
    typedef itk::TransformWrapper<double, SpaceDimension> TrWrapperType;
    TrWrapperType::Pointer trWrapper = TrWrapperType::New();
    trWrapper->SetTransform(transform);
    trWrapper->SetInitTime(image->GetOrigin()[SpaceDimension]);
    
    typedef itk::Image<PixelType, SpaceDimension> OutputImageType;
    typedef itk::ImageFileWriter<OutputImageType> WriterType;
    
    PixelType* pixelContainer = image->GetPixelContainer()->GetBufferPointer();
    unsigned long numberOfVoxelsPerTimePoint = image->GetBufferedRegion().GetNumberOfPixels() 
    / image->GetBufferedRegion().GetSize()[SpaceDimension];
    
    // Moving pixel container ( the first image to deform is the second one)
    pixelContainer += numberOfVoxelsPerTimePoint;
    
    for (unsigned int time = 1; time<image->GetBufferedRegion().GetSize()[SpaceDimension]; time++)
    {
        double contTime = image->GetOrigin()[SpaceDimension] + 
        image->GetSpacing()[SpaceDimension] * (double)time;
        trWrapper->SetFinalTime(contTime);
        
        // This will point to a 3D image containing the 4D image data at the current time
        OutputImageType::Pointer fakeImage = OutputImageType::New();
        OutputImageType::RegionType fakeRegion;
        OutputImageType::SizeType fakeSize;
        OutputImageType::IndexType fakeIndex;
        OutputImageType::PointType  fakeOrigin;
        OutputImageType::SpacingType  fakeSpacing;
        for (unsigned int d=0; d<SpaceDimension; d++)
        {
            fakeOrigin[d] = image->GetOrigin()[d];
            fakeIndex[d] = image->GetBufferedRegion().GetIndex()[d];
            fakeSize[d] = image->GetBufferedRegion().GetSize()[d];
            fakeSpacing[d] = image->GetSpacing()[d];
        }
        
        fakeRegion.SetIndex(fakeIndex);
        fakeRegion.SetSize(fakeSize);
        fakeImage->SetRegions(fakeRegion);
        fakeImage->SetOrigin(fakeOrigin);
        fakeImage->SetSpacing(fakeSpacing);
        fakeImage->GetPixelContainer()->SetImportPointer( pixelContainer, numberOfVoxelsPerTimePoint );
        
        // Generating displacement field
        typedef itk::Vector<float, SpaceDimension> DisplacementPixelType;
        typedef itk::Image<DisplacementPixelType, SpaceDimension> OutputDispType;
        typedef itk::TransformToDeformationFieldFilter< OutputDispType, double > DeformationFieldGeneratorType;
        
        DeformationFieldGeneratorType::Pointer genDef = DeformationFieldGeneratorType::New();
        genDef->SetTransform(trWrapper);
        genDef->SetOutputRegion(fakeImage->GetBufferedRegion());
        genDef->SetOutputSpacing(fakeImage->GetSpacing());
        genDef->SetOutputOrigin(fakeImage->GetOrigin());
        genDef->Update();
        
        typedef itk::ImageFileWriter<DeformationFieldGeneratorType::OutputImageType> DispWriterType;
        DispWriterType::Pointer dispWriter = DispWriterType::New();
        dispWriter->SetInput(genDef->GetOutput());
        
        if (params.outputDispFileName!=0)
        {
            char outputFileName [2048];
            sprintf(outputFileName, "%s%04d.vtk", params.outputDispFileName, time);
            dispWriter->SetFileName(outputFileName);
            dispWriter->Update();
        }
        
        typedef itk::WarpImageFilter<OutputImageType, OutputImageType, OutputDispType> WarpFilterType;
        WarpFilterType::Pointer warp = WarpFilterType::New();
        warp->SetInput(fakeImage);
        warp->SetOutputSpacing(fakeImage->GetSpacing());
        warp->SetOutputOrigin(fakeImage->GetOrigin());
        warp->SetDeformationField(genDef->GetOutput());
        warp->Update();
        
        // Writing this magnificent result
        char outputFileName [2048];
        WriterType::Pointer writer = WriterType::New();
        writer->SetInput(warp->GetOutput());
        sprintf(outputFileName, "%s%03d.mhd", params.outputImageFilePrefix, time);
        
        writer->SetFileName(outputFileName);
        try
        {
            writer->Update();
        }
        catch( itk::ExceptionObject & excp )
        {
            std::cerr << "Error while applying the resampling" << std::endl;
            std::cerr << excp << std::endl;
            return EXIT_FAILURE;
        }
        
        // Moving pixel container 
        pixelContainer += numberOfVoxelsPerTimePoint;
        
    }
    
} 


bool ParseCommandLine(int inargc, char * inargv[], MyGlobalTypesAndParameters& params)
{
    
	char** argv = inargv;
	int argc = inargc;
    
	params.inputImage = argv[1];				
	argc--; 
	argv++;
	params.outputImageFilePrefix = argv[1];		
	argc--; 
	argv++;
	params.outputTransformFile = argv[1];		
	argc--; 
	argv++;
    
	if (argc <= 2){
		std::cerr << std::endl << "Unspecified parameter " << argv[1] << std::endl << std::endl;
		usage();
	}
    
	bool ok;
    
	while (argc >= 2){
        
		ok = false;
        
		if ((ok == false) && (strcmp(argv[1], "-regionSize") == 0)){
			if (argc < 5)
			{
				std::cerr << std::endl << "Missing parameters " << argv[1] << std::endl << std::endl;
				usage();
				return false;
			}
			else
			{
				argc--; argv++;
				for (unsigned int d=0; d<MyGlobalTypesAndParameters::Dimension; d++)
				{
                    params.gridSize[d] = atoi (argv[1]); argc--; argv++;
				}
                
				ok = true;
			}
		}
        
		if ((ok == false) && (strcmp(argv[1], "-inputTransformation") == 0))
		{
			argc--; argv++;
			params.inputTransformationFileName = argv[1]; argc--; argv++;
			ok = true;
		}
        
		if ((ok == false) && (strcmp(argv[1], "-splitTransformation") == 0))
		{
			argc--; argv++;
			params.splitTransformation = true;
		    ok = true;
		}
        
		if ((ok == false) && (strcmp(argv[1], "-inputMask") == 0))
		{
			argc--; argv++;
			params.inputMaskFileName = argv[1]; argc--; argv++;
			ok = true;
		}
		if ((ok == false) && (strcmp(argv[1], "-inputMovingMask") == 0))
		{
			argc--; argv++;
			params.inputMovingMaskFileName = argv[1]; argc--; argv++;
			ok = true;
		}
        
		if ((ok == false) && (strcmp(argv[1], "-inputMeshInvert") == 0))
		{
			argc--; argv++;
			params.incompressibilityPointSetFileName = argv[1]; argc--; argv++;
			ok = true;
		}
        
		if ((ok == false) && (strcmp(argv[1], "-imageRegion") == 0))
		{
			argc--; argv++;
			MyGlobalTypesAndParameters::RegionType::IndexType index;
			MyGlobalTypesAndParameters::RegionType::SizeType size;
			for (unsigned int d=0; d<MyGlobalTypesAndParameters::Dimension; d++)
            {
				index[d] = atoi (argv[1]); argc--; argv++;
            }
			for (unsigned int d=0; d<MyGlobalTypesAndParameters::Dimension; d++)
            {
				size[d] = atoi (argv[1]); argc--; argv++;
            }
			params.imageRegion.SetSize(size);
			params.imageRegion.SetIndex(index);
			ok = true;
		}
        
		if ((ok == false) && (strcmp(argv[1], "-maxv") == 0))
		{
			argc--; argv++;
			params.maximumVelocityInsideMaskX = atof(argv[1]); argc--; argv++;
            params.maximumVelocityInsideMaskY=params.maximumVelocityInsideMaskX;
            params.maximumVelocityInsideMaskZ=params.maximumVelocityInsideMaskX;
			
			ok = true;
		}
        
        if ((ok == false) && (strcmp(argv[1], "-maxvx") == 0))
		{
			argc--; argv++;
			params.maximumVelocityInsideMaskX = atof(argv[1]); argc--; argv++;
			
			ok = true;
		}
        
        if ((ok == false) && (strcmp(argv[1], "-maxvy") == 0))
		{
			argc--; argv++;
			params.maximumVelocityInsideMaskY = atof(argv[1]); argc--; argv++;
			
			ok = true;
		}
        
        if ((ok == false) && (strcmp(argv[1], "-maxvz") == 0))
		{
			argc--; argv++;
			params.maximumVelocityInsideMaskZ = atof(argv[1]); argc--; argv++;
			
			ok = true;
		}
        
        if ((ok == false) && (strcmp(argv[1], "-MetricSequential") == 0))
		{
			argc--; argv++;
			params.metricTypeSequential = true;
			ok = true;
		}
        
        if ((ok == false) && (strcmp(argv[1], "-MetricDense") == 0))
		{
			argc--; argv++;
			params.metricTypeDense = true;
			ok = true;
		}
        
        if ((ok == false) && (strcmp(argv[1], "-weightSequential") == 0))
		{
			argc--; argv++;
			params.seqWeight = atof(argv[1]); argc--; argv++;
			ok = true;
		}
        
		if ((ok == false) && (strcmp(argv[1], "-referenceTime") == 0))
		{
			argc--; argv++;
			params.referenceTime = atof(argv[1]); argc--; argv++;
			
			ok = true;
		}
		if ((ok == false) && (strcmp(argv[1], "-timeStep") == 0))
        {
            argc--; argv++;
            params.timeStep = atof(argv[1]); argc--; argv++;
            
            ok = true;
        }
		if ((ok == false) && (strcmp(argv[1], "-minTimeStep") == 0))
        {
            argc--; argv++;
            params.minTimeStep = atof(argv[1]); argc--; argv++;
            
            ok = true;
        }
        
        if ((ok == false) && (strcmp(argv[1], "-numSamples") == 0))
        {
            argc--; argv++;
            params.numberOfSamples = atoi(argv[1]); argc--; argv++;
            ok = true;
        }
        
        if ((ok == false) && (strcmp(argv[1], "-dispFileName") == 0))
        {
            argc--; argv++;
            params.outputDispFileName= argv[1]; argc--; argv++;
            ok = true;
        }
        
        if ((ok == false) && (strcmp(argv[1], "-padding") == 0))
        {
            argc--; argv++;
            for (unsigned int d=0; d<MyGlobalTypesAndParameters::SpaceDimension; d++)
            {
                params.padding[d] = atoi(argv[1]);
            }
            params.padding[MyGlobalTypesAndParameters::SpaceDimension] = 0; // We don't padd the temporal dimension
            argc--; argv++;
            ok = true;
        }
        
        if ((ok == false) && (strcmp(argv[1], "-paddingx") == 0))
        {
            argc--; argv++;
            params.padding[0] = atoi(argv[1]); 
            argc--; argv++;
            ok = true;
        }
        
        if ((ok == false) && (strcmp(argv[1], "-paddingy") == 0))
        {
            argc--; argv++;
            params.padding[1] = atoi(argv[1]); 
            argc--; argv++;
            ok = true;
        }
        
        if ((ok == false) && (strcmp(argv[1], "-paddingz") == 0))
        {
            argc--; argv++;
            params.padding[2] = atoi(argv[1]); 
            argc--; argv++;
            ok = true;
        }
        
        
        if ((ok == false) && (strcmp(argv[1], "-incompressibilityWeight") == 0))
        {
            argc--; argv++;
            params.incompressibilityWeight= atoi(argv[1]); argc--; argv++;
            ok = true;
        }    
        
        if ((ok == false) && (strcmp(argv[1], "-timeMultiRes") == 0))
        {
            argc--; argv++;
            params.timeMultiRes=true; 
            ok = true;
        } 
		
        if (ok == false)
		{
			std::cerr << std::endl << "Cannot parse argument " << argv[1] << std::endl << std::endl;
			usage();
			return false;
		}
        
        
        
        
	}
    
	return true;
    
    
}

