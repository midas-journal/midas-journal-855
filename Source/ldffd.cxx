/*=========================================================================

Program:   Insight Segmentation & Registration Toolkit
Module:    $$
Language:  C++
Date:      $$
Version:   $$

Copyright (c) Insight Software Consortium. All rights reserved.
See ITKCopyright.txt or http://www.itk.org/HTML/Copyright.htm for details.

This software is distributed WITHOUT ANY WARRANTY; without even
the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
PURPOSE.  See the above copyright notices for more information.

Contact: mathieu.de_craene@philips.com

=========================================================================*/
#if defined(_MSC_VER)
#pragma warning ( disable : 4786 )
#endif

// Just for test
#include <itkImageRandomConstIteratorWithIndex.h>


#include <iostream>
#include <iomanip>
#include "itkSmoothingRecursiveGaussianImageFilter.h"
#include "itkNumericTraits.h"
#include "itkTransformFileReader.h"
#include "itkTransformFactory.h"

#include "itkDiffeomorphicBSplineTransform.h"
#include "itkImage.h"
#include "itkImageFileReader.h"
#include "itkImageFileWriter.h"
#include "itkMattesMutualInformationMultipleImages.h"

#include "itkLinearInterpolateImageFunction.h"
#include "itkImage.h"
#include "itkTimeProbesCollectorBase.h"

#include "itkLBFGSBOptimizer.h"
#include <itkRegularStepGradientDescentOptimizer.h>
#include "itkResampleImageFilter.h"
#include "itkTransformFileWriter.h"
#include "itkVector.h"
#include "itkImageRegionIterator.h"
#include "itkImageRegionConstIteratorWithIndex.h"
#include "vxl/vcl/vcl_cmath.h"
#include "vxl/vcl/vcl_algorithm.h"
#include "BSplineTransformSquareRoot.h"
#include "itkImageMaskSpatialObject.h"

#include "itkHistogramMatchingImageFilter.h"
#include "regFFDSplit.h"

//  The following section of code implements a Command observer
//  used to monitor the evolution of the registration process.
//
#include "itkCommand.h"
class CommandIterationUpdate : public itk::Command
{
public:
  typedef  CommandIterationUpdate   Self;
  typedef  itk::Command             Superclass;
  typedef  itk::SmartPointer<Self>  Pointer;
  itkNewMacro( Self );

  typedef itk::LBFGSBOptimizer     OptimizerType;
  typedef   const OptimizerType   *    OptimizerPointer;
  typedef itk::DiffeomorphicBSplineTransform<double, 3, 3>  TransformType;

  // Transformation object
  TransformType* m_Transform;
  TransformType::Pointer m_InternalTransform;
  TransformType::ParametersType m_Parameters;
  std::string m_FileNamePrefix;
  unsigned long m_IterationCounter;

  CommandIterationUpdate()
    {
    m_InternalTransform = TransformType::New();
    m_FileNamePrefix="";
    m_IterationCounter=0;
    }

  void Execute(itk::Object *caller, const itk::EventObject & event)
    {
    Execute( (const itk::Object *)caller, event);
    }

  void Execute(const itk::Object * object, const itk::EventObject & event)
    {
    OptimizerPointer optimizer =
      dynamic_cast< OptimizerPointer >( object );

    if (optimizer!=0)
      {
      if( !(itk::IterationEvent().CheckEvent( &event )) )
        {
        return;
        }

      std::cout << optimizer->GetCurrentIteration() << "   ";
      std::cout << optimizer->GetCachedValue() << std::endl;

      m_Parameters = optimizer->GetCachedCurrentPosition();
      if (m_Transform)
        {
        m_InternalTransform->SetFixedParameters(m_Transform->GetFixedParameters());
        }
      m_InternalTransform->SetParameters(m_Parameters);

      if (  m_FileNamePrefix != "" )
        {

        // Writing output transformation
        itk::TransformFileWriter::Pointer trwriter = itk::TransformFileWriter::New();
        char iterationString[1024];
        sprintf(iterationString, "%03lu", m_IterationCounter);
        std::string outputFileName = m_FileNamePrefix + iterationString + ".txt";
        trwriter->SetFileName( outputFileName.c_str() );
        trwriter->SetInput(m_InternalTransform);
        try
          {
          trwriter->Update();
          }
        catch ( itk::ExceptionObject & err )
          {
          std::cerr << "An error occured while writing image file" << std::endl;
          std::cerr << err << std::endl;
          }
        }

      }//if (optimizer!=0)


      m_IterationCounter++;

    }
};


// Config object containing all the variables defined by the user
typedef struct configStruct
{

  // input and output file names
  std::string fixedImageFileName;
  std::vector< std::string > movingImageFileNames;
  std::string outputImageFileNamePrefix;
  std::string outputTransformFileName;
  std::string initTransformFileName;
  std::string fixedImageMaskFileName;

  // fixed image region stuff
  itk::ImageRegion<3> fixedImageRegion;
  bool fixedImageRegionSpecified;

  // number of time steps
  unsigned int initialNumberOfTimeSteps;
  unsigned int maxNumberOfTimeSteps;

  // grid size on image
  float gridSpacingOnImage;

  // maximum displacement for each transform
  // (as a fraction of the smalled grid spacing
  // along image dimensions)
  double maxRelativeDisplacement;

  // options for BSpline Sqrt computation
  double sqrtLearningRate;
  unsigned int sqrtNumberOfIterations;

  // whether we want to split the tranfsormation
  bool splitTransform;

  // whether we want to unbound the tranfsormation parameters
  bool unleashed;

  // wether we want to equalize the histogram of all moving images
  // to the fixed imag
  bool equalize;
    
  // Smoothing of input images
  double sigma;


} *configStructPointer;

// Utility function to parse all parameters passed by the user
bool parseConfigFile(char* filename, configStructPointer config);


// Function optimizing the diffeomorphic transform for a given number of
// time steps
template <class TMetric>
typename TMetric::ParametersType
PerformOptimization(
  TMetric* metric,
  typename TMetric::TransformType* transform,
  typename TMetric::InterpolatorType* interpolator,
  const typename TMetric::FixedImageType* fixedImage,
  std::vector<typename TMetric::MovingImageType::Pointer> movingImages,
  typename TMetric::FixedImageType::RegionType fixedRegion,
  typename TMetric::FixedImageMaskType* fixedMask,
  configStructPointer config,
  CommandIterationUpdate* observer,
  double minimalSpacing,
  typename TMetric::ParametersType onTheBorder
)
{

  typedef typename TMetric::ParametersType
    ParametersType;
  typedef itk::LBFGSBOptimizer       OptimizerType;
  typename OptimizerType::Pointer      optimizer     = OptimizerType::New();

  // Setting up registration components
  metric->SetMovingImage( movingImages[0] );
  metric->SetFixedImage( fixedImage );

  unsigned int time = 0;
  for (unsigned int movImage=0; movImage<movingImages.size(); movImage++)
    {
    time++;
    metric->SetMovingImage(movImage, movingImages[movImage], time);
    }

  metric->SetInterpolator( interpolator );
  metric->SetTransform( transform );
  metric->SetFixedImageRegion( fixedRegion );
  if ( fixedMask )
    {
    metric->SetFixedImageMask(fixedMask);
    std::cout << "mask passed to metric" << std::endl;
    }
  else
    {
    std::cout << "No mask has been set" << std::endl;
    }

  metric->Initialize();

  // Setup the optimizer
  optimizer->SetCostFunction( metric );
  optimizer->SetInitialPosition( transform->GetParameters() );

  // Configuring optimizer
  //
  
  
  // Bounds
  if (!config->unleashed)
    {
    OptimizerType::BoundSelectionType boundSelect( transform->GetNumberOfParameters() );
    OptimizerType::BoundValueType upperBound( transform->GetNumberOfParameters() );
    OptimizerType::BoundValueType lowerBound( transform->GetNumberOfParameters() );

    std::cout << "Minimal Spacing: " << minimalSpacing << std::endl;

    boundSelect.Fill( 2 );
    upperBound.Fill( minimalSpacing * 0.4 );
    lowerBound.Fill( -1. * minimalSpacing * 0.4 );

    for (unsigned int index=0; index<transform->GetNumberOfParameters();index++)
      if (onTheBorder[index] != 0.)
        {
        lowerBound[index] *= 0.5;
        upperBound[index] *= 0.5;
        }

    optimizer->SetBoundSelection( boundSelect );
    optimizer->SetUpperBound( upperBound );
    optimizer->SetLowerBound( lowerBound );
    }
  else
    {
    OptimizerType::BoundSelectionType boundSelect( transform->GetNumberOfParameters() );
    OptimizerType::BoundValueType upperBound( transform->GetNumberOfParameters() );
    OptimizerType::BoundValueType lowerBound( transform->GetNumberOfParameters() );

    boundSelect.Fill(0);
    upperBound.Fill(0.);
    lowerBound.Fill(0.);
    
    optimizer->SetBoundSelection( boundSelect );
    optimizer->SetUpperBound( upperBound );
    optimizer->SetLowerBound( lowerBound );
    }
  
  optimizer->SetCostFunctionConvergenceFactor( 1e+7 );
  optimizer->SetProjectedGradientTolerance( 1e-4 );
  optimizer->SetMaximumNumberOfIterations( 75 );
  optimizer->SetMaximumNumberOfEvaluations( 250 );
  optimizer->SetMaximumNumberOfCorrections( 12 );

  // Register the observer with the optimizer.
  optimizer->AddObserver( itk::IterationEvent(), observer );

  // Add a time probe
  itk::TimeProbesCollectorBase collector;
  std::cout << std::endl << "Starting Registration" << std::endl;

  try
    {
    collector.Start( "Registration" );
    optimizer->StartOptimization();
    collector.Stop( "Registration" );
    }
  catch( itk::ExceptionObject & err )
    {
    std::cerr << "ExceptionObject caught !" << std::endl;
    std::cerr << err << std::endl;
    throw err;
    }

  // Report the time taken by the registration
  collector.Report();

  return optimizer->GetCurrentPosition();

}



int main( int argc, char * argv[] )
{

  configStructPointer config = new configStruct;

  if (argc != 2)
    {
    std::cout << "incorrect number of arguments (" << argc << ")" <<std::endl;
    std::cout << "usage: " << argv[0]
      << " configFile";
    std::cout << std::endl;
    return -1;
    }

  if ( ! parseConfigFile(argv[1], config) )
  {
    std::cerr << "Error while reading configuration file" << std::endl;
    return -1;
  }

  const    unsigned int    ImageDimension = 3;
  const    unsigned int    BSplineOrder = 3;
  typedef  float    PixelType;

  typedef itk::Image< PixelType, ImageDimension >  ImageType;

  typedef itk::LBFGSBOptimizer       OptimizerType;

  typedef itk::MattesMutualInformationMultipleImages< ImageType,
          ImageType > MetricType;

  typedef itk::LinearInterpolateImageFunction< ImageType,
          double > InterpolatorType;

  typedef itk::DiffeomorphicBSplineTransform<double,ImageDimension, 3>  TransformType;

  typedef itk::ImageFileReader< ImageType  > FixedImageReaderType;
  typedef itk::ImageFileReader< ImageType > MovingImageReaderType;
  typedef itk::SmoothingRecursiveGaussianImageFilter<ImageType, ImageType >  SmootherType;

  typedef itk::ResampleImageFilter< ImageType,
          ImageType >    ResampleFilterType;


  typedef itk::ImageFileWriter<ImageType> WriterType;

  MetricType::Pointer         metric        = MetricType::New();
  OptimizerType::Pointer      optimizer     = OptimizerType::New();
  InterpolatorType::Pointer   interpolator  = InterpolatorType::New();

  // Fixed image reading
  FixedImageReaderType::Pointer fixedImageReader = FixedImageReaderType::New();
  fixedImageReader->SetFileName(  config->fixedImageFileName );
  

  try
    {
    fixedImageReader->Update();
    }
  catch( itk::ExceptionObject & err )
    {
    std::cerr << "ExceptionObject caught !" << std::endl;
    std::cerr << err << std::endl;
    return -1;
    }

  ImageType::ConstPointer fixedImage;
  
  if (config->sigma == 0.)
  {
    fixedImage = fixedImageReader->GetOutput();
  }else
  {
      SmootherType::Pointer smoother = SmootherType::New();
      smoother->SetInput(fixedImageReader->GetOutput());
      smoother->SetSigma(config->sigma);
      smoother->Update();
      fixedImage=smoother->GetOutput();
  }
    
  // Reading all the moving images
  std::vector<ImageType::Pointer> movingImages;

  // Typedef for possible equalization
  typedef itk::HistogramMatchingImageFilter< ImageType, ImageType > HistoMatchingFilterType;
  
  for (unsigned int movImage=0; movImage < config->movingImageFileNames.size(); movImage++)
  {
  ImageType::Pointer movingImageTmp = 0; 
  MovingImageReaderType::Pointer movingImageReader = MovingImageReaderType::New();
  movingImageReader->SetFileName(  config->movingImageFileNames[movImage] );
  movingImageReader->Update();
  movingImageTmp =  movingImageReader->GetOutput();
  
  if (config->equalize)
    {
      std::cout << "Equalizing moving image " << movImage << "." << std::endl;  
      HistoMatchingFilterType::Pointer matcher = HistoMatchingFilterType::New();
      matcher->SetInput( movingImageTmp );
      matcher->SetReferenceImage( fixedImage );
      matcher->SetNumberOfHistogramLevels(100);
      matcher->SetNumberOfMatchPoints(3);
      matcher->ThresholdAtMeanIntensityOff();
      matcher->Update();
      movingImageTmp = matcher->GetOutput();
    }
    
    if (config->sigma == 0.)
    {
          movingImages.push_back( movingImageTmp );
    }else
    {
          SmootherType::Pointer smoother = SmootherType::New();
          smoother->SetInput(movingImageTmp);
          smoother->SetSigma(config->sigma);
          smoother->Update();
          movingImages.push_back( smoother->GetOutput() );
    }
 
  
  
  }

  // Dealing with fixed image region ...
  ImageType::RegionType fixedRegion = fixedImage->GetBufferedRegion();

  // Checking if a region has been specified by user
  if (config->fixedImageRegionSpecified)
    {
    // If yes, check if it is inside the fixed image buffered region
    if( fixedRegion.IsInside(config->fixedImageRegion) )
      {
      fixedRegion = config->fixedImageRegion;
      }
    else
      {
      std::cerr << "Input region is not contained in fixed image buffered region." << std::endl;
      std::cerr << "Continuting with fixed image buffered region as fixed region." << std::endl;
      std::cerr << "Input region : " << config->fixedImageRegion << std::endl;
      std::cerr << "Fixed region : " << fixedRegion << std::endl;

      }
    }

  // Fixed image mask
  typedef itk::Image< unsigned char, ImageDimension > ImageMaskType;
  typedef itk::ImageMaskSpatialObject< ImageDimension >  MaskType;
  MaskType::Pointer spatialObjectMask = 0;

  if ( config->fixedImageMaskFileName != "" )
    {
    spatialObjectMask = MaskType::New();

    typedef itk::ImageFileReader<ImageMaskType> ImageMaskReaderType;
    ImageMaskReaderType::Pointer maskReader = ImageMaskReaderType::New();
    maskReader->SetFileName(config->fixedImageMaskFileName.c_str());

    try
      {
      maskReader->Update();
      std::cout << "Mask is read" << std::endl;
      }
    catch( itk::ExceptionObject & err )
      {
      std::cerr << "ExceptionObject caught !" << std::endl;
      std::cerr << err << std::endl;
      return -1;
      }

    spatialObjectMask->SetImage( maskReader->GetOutput() );
    }

  // Specifying transformation region, origin and spacing
  typedef TransformType::RegionType RegionType;
  RegionType bsplineRegion;
  RegionType::SizeType   gridSizeOnImage;
  RegionType::SizeType   gridBorderSize;
  RegionType::SizeType   totalGridSize;

  typedef TransformType::SpacingType SpacingType;
  SpacingType spacing = fixedImage->GetSpacing();
  ImageType::SizeType fixedImageSize = fixedRegion.GetSize();

  // Computing dimension of the BSpline grid in each dimension
  // from the grid spacing specified in the configuration file
  for (unsigned int d=0; d<ImageDimension; d++)
    {
    gridSizeOnImage[d] = (unsigned int)
      floor( static_cast<double>(fixedImageSize[d] - 1)  *
        static_cast<double>(spacing[d]) /
        static_cast<double>(config->gridSpacingOnImage) );
    }

  std::cout << "Grid size on image is " << gridSizeOnImage << std::endl;

  gridBorderSize.Fill( 3 );    // Border for spline order = 3 ( 1 lower, 2 upper )
  totalGridSize = gridSizeOnImage + gridBorderSize;

  bsplineRegion.SetSize( totalGridSize );

  typedef TransformType::OriginType OriginType;

  OriginType origin;
  fixedImage->TransformIndexToPhysicalPoint(fixedRegion.GetIndex(), origin);

  for(unsigned int r=0; r<ImageDimension; r++)
    {
    spacing[r] *= floor( static_cast<double>(fixedImageSize[r] - 1)  /
      static_cast<double>(gridSizeOnImage[r] - 1) );
    origin[r]  -=  spacing[r];
    }

  // Creating transform or reading it from file
  TransformType::Pointer  transform = 0;
  TransformType::ParametersType parameters;

  // Loading initial transform if specified
  if (config->initTransformFileName == "")
    {
    // Create transform from scratch
    transform = TransformType::New();
    transform->SetGridSpacing( spacing );
    transform->SetGridOrigin( origin );
    transform->SetGridRegion( bsplineRegion );
    transform->SetNumberOfTimeSteps(movingImages.size());

    parameters = TransformType::ParametersType(transform->GetNumberOfParameters());
    parameters.Fill(0.);
    transform->SetParameters(parameters);

    }
  else
    {
    // Read transform from file
    itk::TransformFactory<TransformType>::RegisterTransform();
    itk::TransformFileReader::Pointer tr_reader;
    tr_reader = itk::TransformFileReader::New();
    tr_reader->SetFileName( config->initTransformFileName );
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
    }


    // Splitting transforms if needed
    TransformType::Pointer transformOld = 0;
    if (config->splitTransform)
    {
        const std::vector<TransformType::InternalTransformPointerType>&
          internalTransforms = transform->GetInternalTransforms();

        typedef TransformType::InternalTransformPointerType InternalTransformPointerType;
        typedef TransformType::InternalTransformType InternalTransformType;

        std::vector<InternalTransformPointerType> outputTransforms;
        outputTransforms.resize( internalTransforms.size() );

        // increasing the spatial resolution of each FFD transformation
        // by a factor 2
        for (unsigned int tr=0; tr<outputTransforms.size(); tr++)
        {
            outputTransforms[tr]
               = RefineFFDTransformation<InternalTransformType,ImageType>(internalTransforms[tr].GetPointer(), BSplineOrder);

        }

        TransformType::Pointer transform_temp = TransformType::New();
        transform_temp->SetGridSpacing( outputTransforms[0]->GetGridSpacing() );
        transform_temp->SetGridOrigin( outputTransforms[0]->GetGridOrigin() );
        transform_temp->SetGridRegion( outputTransforms[0]->GetGridRegion() );
        transform_temp->SetNumberOfTimeSteps( transform->GetNumberOfTimeSteps() );

        parameters = TransformType::ParametersType(transform_temp->GetNumberOfParameters());
        unsigned int init = 0;
        for (unsigned int tr=0; tr<outputTransforms.size(); tr++)
        {
            for (unsigned int index=0; index<outputTransforms[tr]->GetNumberOfParameters(); index++)
            {
                parameters[init+index] = outputTransforms[tr]->GetParameters()[index];
            }
            init += outputTransforms[tr]->GetNumberOfParameters();
        }

        // Backuping transform in transformOld and setting transform to the new transformation
        transform_temp->SetParametersByValue(parameters);
        transformOld = transform;
        transform = transform_temp;
    }

  // Computing the smallest grid spacing of the transformation
  // among all dimensions
  double smallestGridSpacing=transform->GetGridSpacing()[0];
  for (unsigned int d=1; d<ImageDimension; d++)
    {
    if (transform->GetGridSpacing()[d] < smallestGridSpacing)
      smallestGridSpacing=transform->GetGridSpacing()[d];
    }

    
  // Setting up metric
  metric->SetNumberOfHistogramBins( 50 );
  const unsigned int numberOfSamples = vcl_min(10000,
					       (int) (fixedRegion.GetNumberOfPixels() / 10));

  std::cout << "Metric configuration: number of samples being used: " << numberOfSamples << std::endl;

  metric->SetNumberOfSpatialSamples( numberOfSamples );
  metric->SetNumberOfMovingImages( movingImages.size() );
  //  Given that the Mattes Mutual Information metric uses a random iterator in
  //  order to collect the samples from the images, it is usually convenient to
  //  initialize the seed of the random number generator.
  metric->ReinitializeSeed( 76926294 );

  // Create observer
  CommandIterationUpdate::Pointer observer = CommandIterationUpdate::New();
  observer->m_FileNamePrefix = "myIntermediateTransformResult";
  observer->m_Transform = transform;

  // Finding out which control points are on the border
  TransformType::InternalTransformType::ImagePointer coeff 
    = transform->GetInternalTransforms()[0]->GetCoefficientImage()[0];
  typedef itk::ImageRegionConstIteratorWithIndex<TransformType::InternalTransformType::ImageType> IteratorType;

  IteratorType it (coeff, coeff->GetBufferedRegion());
  TransformType::InternalTransformType::ImageType::PixelType*  first = coeff->GetBufferPointer();

  TransformType::ParametersType onTheBorder(transform->GetNumberOfParameters());
  onTheBorder.Fill(0.);
  
  unsigned int numberOfParametersPerTimePointAndDimension =
    transform->GetInternalTransforms()[0]->GetCoefficientImage()[0]->GetBufferedRegion().GetNumberOfPixels();
  
  for ( it.GoToBegin(); !it.IsAtEnd();++it)
    {
    unsigned long index = &( it.Value() ) - first;
    TransformType::InternalTransformType::ImageType::IndexType itIndex = it.GetIndex();

    // Checking if it is on the border or not
    bool isOnBorder = false;
    for (unsigned int d=0; d<coeff->GetImageDimension();d++)
      {
      if (itIndex[d] == 0)
        isOnBorder=true;
      }

    // If yes, we write all parametric positions to 1 for this control point
    // in the onTheBorder vector.
    if (isOnBorder)
      for (unsigned int d=0; d<coeff->GetImageDimension();d++)
        {
        const unsigned long habemusIndex = index + d*numberOfParametersPerTimePointAndDimension;
        onTheBorder[ habemusIndex ] = 1.;
        }
    }

  // We still have to constrain the same positions at all next time points
  unsigned long counterBorder = 0;
  for (unsigned int p=0; p<transform->GetNumberOfParametersByTransform(); p++)
    {
    if (onTheBorder[p] != 0.)
      {
      counterBorder++;
      for (unsigned int t=1; t<transform->GetNumberOfTimeSteps(); t++)
        {
        const unsigned long habemusIndex = p+t*transform->GetNumberOfParametersByTransform();
        onTheBorder[habemusIndex] = 1.;
        counterBorder++;
        }
      }
    }
  std::cout << "Total Number Of Parameters : " << transform->GetNumberOfParameters() << std::endl;
  std::cout << "Number of 'on the border' parameters  : " <<counterBorder << std::endl;
    
  // Call registration
  TransformType::ParametersType currentParameters;

  currentParameters =
        PerformOptimization
        <MetricType>
        ( metric, transform, interpolator, fixedImage, movingImages,
          fixedRegion, spatialObjectMask, config, observer, smallestGridSpacing, onTheBorder);


  // Passing final parameters
  transform->SetParameters( currentParameters );
     
    
     

  // Writing output transformation
  itk::TransformFileWriter::Pointer trwriter = itk::TransformFileWriter::New();
  trwriter->SetFileName( config->outputTransformFileName );
  trwriter->SetInput(transform);

  try
    {
    trwriter->Update();
    }
  catch ( itk::ExceptionObject & err )
    {
    std::cerr << "An error occured while writing image file" << std::endl;
    std::cerr << err << std::endl;
    return EXIT_FAILURE;
    }

  // Final resampling
  unsigned int ts=1;

  for (unsigned int movImage=0; movImage < movingImages.size(); movImage++)
    {

    transform->SetNumberOfTimeSteps(ts);

    ResampleFilterType::Pointer resample = ResampleFilterType::New();
    resample->SetTransform( transform );
    resample->SetInput( movingImages[movImage] );
    resample->SetOutputParametersFromImage(fixedImage);
    resample->SetDefaultPixelValue( 100 );

    std::ostringstream outputFileName;
    outputFileName << config->outputImageFileNamePrefix
      << std::dec << std::setw(3) << movImage << ".mhd";

    WriterType::Pointer      writer =  WriterType::New();

    writer->SetFileName( outputFileName.str().c_str() );
    writer->SetInput( resample->GetOutput() );
    try
      {
      writer->Update();
      }
    catch( itk::ExceptionObject & err )
      {
      std::cerr << "ExceptionObject caught !" << std::endl;
      std::cerr << err << std::endl;
      return EXIT_FAILURE;
      }

    ts++;

    }

  return EXIT_SUCCESS;

}

bool parseConfigFile(char* filename, configStructPointer config)
{

  // Initialization of the config object
  config->fixedImageFileName="";
  config->movingImageFileNames.resize(0);
  config->outputImageFileNamePrefix="";
  config->outputTransformFileName="";
  config->initTransformFileName="";
  config->fixedImageMaskFileName="";

  config->fixedImageRegionSpecified=false;
  config->initialNumberOfTimeSteps=4;
  config->maxNumberOfTimeSteps=8;
  config->gridSpacingOnImage=10;
  config->maxRelativeDisplacement=0.4;
  config->sqrtLearningRate=0.01;
  config->sqrtNumberOfIterations=10;
  config->splitTransform=false;
  config->unleashed=false;
  config->equalize = false;
  config->sigma = 0.;
  
  // Parsing input file
  std::ifstream initFile(filename);
  if( initFile.fail() )
    {
    std::cout << "could not open file: " << filename << std::endl;
    return false;
    }

  while( !initFile.eof() )
    {

    if ( (std::ws(initFile).peek())=='#' )
      {
      initFile.ignore(itk::NumericTraits<int>::max(), '\n');
      }
    else
      {

      std::string dummy;
      initFile >> dummy;

      // Fixed image file name
      if(dummy == "-fixedImageFileName")
        {
        initFile >> dummy;
        config->fixedImageFileName = dummy;
        std::cout << "PARAM: fixed image file name is "<<config->fixedImageFileName<<std::endl;
        }

      // Moving image file name
      else if(dummy == "-movingImageFileName")
        {
        initFile >> dummy;
        unsigned int tmpSize = config->movingImageFileNames.size();
        config->movingImageFileNames.push_back( dummy );
        std::cout << "PARAM: adding moving image file name "<<config->movingImageFileNames.at(tmpSize)<<std::endl;
        }

      // Output image file name
      else if(dummy == "-outputImageFileNamePrefix")
        {
        initFile >> dummy;
        config->outputImageFileNamePrefix = dummy;
        std::cout << "PARAM: output image file prefix is "<<config->outputImageFileNamePrefix<<std::endl;
        }

      // Fixed image mask
      else if(dummy == "-fixedImageMaskFileName")
        {
        initFile >> dummy;
        config->fixedImageMaskFileName = dummy;
        std::cout << "PARAM: fixed image mask file name is "<<config->fixedImageMaskFileName<<std::endl;
        }

      // Fixed image region
      else if(dummy == "-fixedImageRegion")
        {
        // 3 Indexes defining where the region is starting
        itk::Index<3> startIndex;
        for (unsigned int dim=0; dim<3; dim++)
          {
          initFile >> dummy;
          startIndex[dim]= atoi(dummy.c_str());
          }

        config->fixedImageRegion.SetIndex(startIndex);
        // 3 Indexes defining the size of the region
        itk::Size<3> size;
        for (unsigned int dim=0; dim<3; dim++)
          {
          initFile >> dummy;
          size[dim]= atoi(dummy.c_str());
          }

        config->fixedImageRegion.SetSize(size);
        config->fixedImageRegionSpecified=true;

        std::cout << "PARAM: fixed image region is "<<config->fixedImageRegion<<std::endl;
        }

      // Initial number of time steps
      else if(dummy == "-initialNumberOfTimeSteps")
        {
        initFile >> dummy;
        config->initialNumberOfTimeSteps = atoi(dummy.c_str());
        std::cout << "PARAM: intial number of time steps is "<<config->initialNumberOfTimeSteps<<std::endl;
        }

      // Max number of time steps
      else if(dummy == "-maxNumberOfTimeSteps")
        {
        initFile >> dummy;
        config->maxNumberOfTimeSteps = atoi(dummy.c_str());
        std::cout << "PARAM: max number of time steps is "<<config->maxNumberOfTimeSteps<<std::endl;
        }

      // Output transformation filename
      else if(dummy == "-outputTransformFileName")
        {
        initFile >> dummy;
        config->outputTransformFileName = dummy.c_str();
        std::cout << "PARAM: output transform file name is  "<<config->outputTransformFileName<<std::endl;
        }

      // Initial transformation filename
      else if(dummy == "-inputTransformFileName")
        {
        initFile >> dummy;
        config->initTransformFileName = dummy.c_str();
        std::cout << "PARAM: initial transform file name is  "<<config->initTransformFileName<<std::endl;
        }

      // Initial transformation filename
      else if(dummy == "-inputTransformFileName")
        {
        initFile >> dummy;
        config->initTransformFileName = dummy.c_str();
        std::cout << "PARAM: initial transform file name is  "<<config->initTransformFileName<<std::endl;
        }

      // BSpline grid resolution
      else if(dummy == "-gridSpacingOnImage")
        {
        initFile >> dummy;
        config->gridSpacingOnImage = atof(dummy.c_str());
        std::cout << "PARAM: bspline grid spacing is  "<<config->gridSpacingOnImage<<std::endl;
        }

      // Maximum relative displacement
      else if(dummy == "-maxRelativeDisplacement")
        {
        initFile >> dummy;
        config->maxRelativeDisplacement = (double)(atof(dummy.c_str()));
        std::cout << "PARAM: maximum relative displacement is  "<<config->maxRelativeDisplacement<<std::endl;
        }

      // Sqrt learning rate
      else if(dummy == "-sqrtLearningRate")
        {
        initFile >> dummy;
        config->sqrtLearningRate = (double)(atof(dummy.c_str()));
        std::cout << "PARAM: sqrt learning rate is  "<<config->sqrtLearningRate<<std::endl;
        }

      // Sqrt number of iterations
      else if(dummy == "-sqrtNumberOfIterations")
        {
        initFile >> dummy;
        config->sqrtNumberOfIterations = atoi(dummy.c_str());
        std::cout << "PARAM: sqrt number of iterations is  "<<config->sqrtNumberOfIterations<<std::endl;
        }

      // splitTransform or not
      else if (dummy == "-splitTransform")
        {
        config->splitTransform=true;
        }

      // splitTransform or not
      else if (dummy == "-unleashed")
        {
        config->unleashed=true;
        }

      // equalize 
      else if (dummy == "-equalize")
        {
      config->equalize=true;
        }
          
      // smoothing 
      else if (dummy == "-sigma")
        {
          initFile >> dummy;
          config->sigma = (double)(atof(dummy.c_str()));
        }          

      }
    }

  // Checking if mandatory arguments have been specified.
  if (config->fixedImageFileName=="")
      return false;
  if (config->movingImageFileNames.size() == 0)
      return false;
  if (config->outputImageFileNamePrefix =="")
      return false;
  if (config->outputTransformFileName=="")
      return false;

  initFile.close();

  return true;

}


