/****************************************************************/
/*               DO NOT MODIFY THIS HEADER                      */
/* MOOSE - Multiphysics Object Oriented Simulation Environment  */
/*                                                              */
/*           (c) 2010 Battelle Energy Alliance, LLC             */
/*                   ALL RIGHTS RESERVED                        */
/*                                                              */
/*          Prepared by Battelle Energy Alliance, LLC           */
/*            Under Contract No. DE-AC07-05ID14517              */
/*            With the U. S. Department of Energy               */
/*                                                              */
/*            See COPYRIGHT for full restrictions               */
/****************************************************************/

// MOOSE includes
#include "ImageSamplerItk.h"
#include "MooseApp.h"
#include "ImageMeshItk.h"

#include "itkGrayscaleFillholeImageFilter.h"
#include "itkGradientAnisotropicDiffusionImageFilter.h"

#include "itkCannyEdgeDetectionImageFilter.h"  
#include "itkImageFileReader.h"
#include "itkImageFileWriter.h"
#include "itkCastImageFilter.h"
#include "itkRescaleIntensityImageFilter.h"


#include <iostream>
#include "itkCannyEdgeDetectionImageFilter.h"
#include "itkImageFileReader.h"
#include "itkImageFileWriter.h"
#include "itkSimpleFilterWatcher.h"
#include "itkMath.h"
#include "itkTestingMacros.h"

#include "itkThresholdImageFilter.h"
#include "itkConnectedThresholdImageFilter.h"
#include "itkCurvatureFlowImageFilter.h"
#include "itkConnectedThresholdImageFilter.h"

#include "itkImage.h"
#include "itkImageFileWriter.h"
#include "itkConnectedThresholdImageFilter.h"
#include "itkCastImageFilter.h"

#include "itkThresholdImageFilter.h"
#include "itkImage.h"
#include "itkImageFileReader.h"
#include "itkImageFileWriter.h"
#include "itkCastImageFilter.h"





typedef   float           FloatPixelType;
  typedef itk::Image< FloatPixelType, 3 > FloatImageType;


template <>
InputParameters
validParams<ImageSamplerItk>()
{
  // Define the general parameters
  InputParameters params = emptyInputParameters();
  params += validParams<FileDicomChoose>();
  params.addParam<Point>("origin", "Origin of the image (defaults to mesh origin)");
  params.addParam<unsigned int>(
      "component",
      "The image RGB-component to return, leaving this blank will result in a greyscale value "
      "for the image to be created. The component number is zero based, i.e. 0 returns the first "
      "(RED) component of the image.");
  // Threshold parameters
  params.addParam<double>("threshold", "The threshold value");
  params.addParam<double>(
      "upper_value", 1, "The value to set for data greater than the threshold value");
  params.addParam<double>(
      "lower_value", 0, "The value to set for data less than the threshold value");
  params.addParamNamesToGroup("threshold upper_value lower_value", "Threshold");
  return params;
}

ImageSamplerItk::ImageSamplerItk(const InputParameters & parameters)
: FileDicomChoose(parameters),
  _is_pars(parameters),
  _is_console((parameters.getCheckedPointerParam<MooseApp *>("_moose_app"))->getOutputWarehouse())
  {
  #ifndef LIBMESH_HAVE_VTK
  // This should be impossible to reach, the registration of ImageSampler is also guarded with
  // LIBMESH_HAVE_VTK
  mooseError("libMesh must be configured with ITK enabled to utilize ImageSamplerItk");
  #endif
  }

void
ImageSamplerItk::setupImageSampler(MooseMesh & mesh)
{
  // Don't warn that mesh or _is_pars are unused when ITK is not enabled.
  libmesh_ignore(mesh);
  libmesh_ignore(_is_pars);

  // Get access to the Mesh object
  MeshTools::BoundingBox bbox = MeshTools::bounding_box(mesh.getMesh());

  // Set the dimensions from the Mesh if not set by the User
  if (_is_pars.isParamValid("dimensions"))
    _physical_dims = _is_pars.get<Point>("dimensions");

  else
  {
    _physical_dims(0) = bbox.max()(0) - bbox.min()(0);
  #if LIBMESH_DIM > 1
    _physical_dims(1) = bbox.max()(1) - bbox.min()(1);
  #endif
  #if LIBMESH_DIM > 2
    _physical_dims(2) = bbox.max()(2) - bbox.min()(2);
  #endif
  }

  // Set the origin from the Mesh if not set in the input file
  if (_is_pars.isParamValid("origin"))
    _origin = _is_pars.get<Point>("origin");
  else
  {
    _origin(0) = bbox.min()(0);
  #if LIBMESH_DIM > 1
    _origin(1) = bbox.min()(1);
  #endif
  #if LIBMESH_DIM > 2
    _origin(2) = bbox.min()(2);
  #endif
  }


  if (_status != 0)
  {
    // We don't have parameters, so see if we can get them from ImageMesh
    ImageMeshItk * image_mesh = dynamic_cast<ImageMeshItk *>(&mesh);
    if (!image_mesh)
      mooseError("Not able to cast image_mesh.");
  }

  try
  {
    _is_console << "Read...     " << std::endl<<std::endl;
    imageSize =reader->GetOutput()->GetLargestPossibleRegion().GetSize();
    _is_console  <<"DICOM Serie Dimension:    " << imageSize<<std::endl;
    const ImageType::SpacingType& inputSpacing=reader->GetOutput()->GetSpacing();
    _is_console <<"DICOM Serie Spacing:      " << inputSpacing << std::endl;
    const ImageType::PointType & origin = reader->GetOutput()->GetOrigin();
    _is_console  <<"DICOM Serie Origin:       "<<  origin << std::endl;

    _is_console << "                            ...image read!" << std::endl;
    _is_console << "Applying the filters...    " << std::endl<<std::endl;

 
    scaledImage=(reader->GetOutput());
    scaledImage->SetSpacing( inputSpacing );
    scaledImage->GetLargestPossibleRegion();
    scaledImage->Update();

    for (unsigned int i = 0; i < 3; ++i)
    {  
     _voxel.push_back(_physical_dims(i) / (imageSize[i]));
    }

   ImageType::SpacingType spacing2;
    
    spacing2[0] = _voxel[0]; // spacing along X
    spacing2[1] = _voxel[1];
    spacing2[2] = _voxel[2];
    
    scaledImage->SetSpacing(spacing2);

    ImageType::PointType newOrigin;
    
    newOrigin[0] = _origin(0); // spacing along X
    newOrigin[1] = _origin(1);
    newOrigin[2] = _origin(2);

    scaledImage->SetOrigin( newOrigin);       
    _is_console  <<"DICOM Serie New Origin:   " << newOrigin<<std::endl<<std::endl;





typedef itk::GrayscaleFillholeImageFilter<ImageType,ImageType>  FillholeFilterType;
FillholeFilterType::Pointer  fillhole = FillholeFilterType::New();

typedef itk::CastImageFilter< ImageType,OutputImageType > CastFilterType;
CastFilterType::Pointer caster = CastFilterType::New();
CastFilterType::Pointer caster2 = CastFilterType::New();

typedef itk::CurvatureFlowImageFilter< ImageType, ImageType > CurvatureFlowImageFilterType;
CurvatureFlowImageFilterType::Pointer smoothing = CurvatureFlowImageFilterType::New();

typedef itk::ConnectedThresholdImageFilter< ImageType, ImageType > ConnectedFilterType;
ConnectedFilterType::Pointer connectedThreshold = ConnectedFilterType::New();



// smoothing->SetNumberOfIterations( 0 ); smoothing->SetTimeStep( 0.125 );

// smoothing->SetInput( reader->GetOutput() );

  fillhole->SetInput( scaledImage );
  fillhole->Update();

connectedThreshold->SetInput(scaledImage); 
connectedThreshold->SetLower(0 ); 
connectedThreshold->SetUpper( 250 );
connectedThreshold->SetReplaceValue( 255 );
connectedThreshold->SetInput( scaledImage ); 



ImageType::IndexType  index;
  index[0]=101;
  index[1]=177;
  //index[2]=4;

std::cout << "Seed index = " << index << std::endl;
connectedThreshold->SetSeed(index);
connectedThreshold->Update();


caster->SetInput(connectedThreshold->GetOutput());
caster->Update();

imageSize =connectedThreshold->GetOutput()->GetLargestPossibleRegion().GetSize();   
std::cout<<imageSize<<": that one was ImageSize"<<std::endl;

unsigned char  pixel_value; 
pixel_value= scaledImage->GetPixel( index ); 
std::cout<<"pixel_value : "<<static_cast<unsigned>(pixel_value) <<std::endl;
caster2->SetInput( scaledImage);


//caster->SetInput( connectedThreshold->GetOutput() );
// writer->SetInput( caster->GetOutput() );

 //typedef itk::ThresholdImageFilter< ImageType >  FilterType;
 //FilterType::Pointer filter = FilterType::New();
  // filter->SetInput( fillhole->GetOutput() );
  // filter->SetOutsideValue( 0 );
  // filter->ThresholdBelow( 180 );
  // filter->ThresholdAbove( 180 );
  // filter->Update();



  // typedef itk::ThresholdImageFilter<InternalImageType> ThresholdFilterType;
  // typename ThresholdFilterType::Pointer  thresholdfilter = ThresholdFilterType::New();


  // thresholdfilter->SetInput(fillhole->GetOutput());
  // thresholdfilter->SetOutsideValue(255);

  // thresholdfilter->ThresholdAbove(100);
  // thresholdfilter->Update();






typedef itk::GradientMagnitudeImageFilter<ImageType, ImageType >  GradientMagnitudeImageFilterType;
 
//   GradientMagnitudeImageFilterType::Pointer gradientMagnitudeImageFilter2 = GradientMagnitudeImageFilterType::New();
//   gradientMagnitudeImageFilter2->SetInput(thresholdfilter->GetOutput()   );
//   gradientMagnitudeImageFilter2->Update();


// connectedThreshold->SetLower( 100.0 ); 
// connectedThreshold->SetUpper( 230.0 );
// connectedThreshold->SetReplaceValue( 255 );
// connectedThreshold->SetInput( gradientMagnitudeImageFilter2->GetOutput() ); 

// InternalImageType::IndexType  index;
//   index[0]=60;
//   index[0]=116;
// connectedThreshold->SetSeed( index);



 // rescaler = RescaleFilterType::New();
 //  rescaler->SetOutputMinimum(   0 );
 // rescaler->SetOutputMaximum( 255 );
 // rescaler->SetInput( connectedThreshold->GetOutput() );
 //  rescaler->Update();


  itk::TIFFImageIOFactory::RegisterOneFactory();
  typedef itk::ImageFileWriter< ImageType >  Writer2Type2;
  Writer2Type2::Pointer writer2 = Writer2Type2::New();

 writer2->SetFileName("outputFilenameFloat2.tiff" );
     writer2->SetInput( connectedThreshold->GetOutput() );
 
      try
      {
      writer2->Update();
      }
   catch (itk::ExceptionObject & e)
      {
      std::cerr << e << std::endl;
      mooseError("Exception in file writer");
      }

  typedef itk::ImageFileWriter< OutputImageType >  WriterType;
  WriterType::Pointer writer3 = WriterType::New();

    writer3->SetFileName("outputFilename.tiff" );
     writer3->SetInput( caster2->GetOutput() );
 
      try
      {
      writer3->Update();
      }
   catch (itk::ExceptionObject & e)
      {
      std::cerr << e << std::endl;
      mooseError("Exception in file writer");
      }

  // Parse arguments
  std::string strThreshold = ".009";
  float threshold = 0.0;
  std::stringstream ssThreshold;
  ssThreshold << strThreshold;
  ssThreshold >> threshold;
 
  std::string strLevel ="0.9";
  float level = 0.0;
  std::stringstream ssLevel;
  ssLevel << strLevel;
  ssLevel >> level;
 
  // Output arguments
  std::cout << "Running with:" << std::endl
            << "Threshold: " << threshold << std::endl
            << "Level: " << level << std::endl;
 



// typedef itk::GradientAnisotropicDiffusionImageFilter<ImageType, ImageType > DiffusionFilterType;

// DiffusionFilterType::Pointer diffusion = DiffusionFilterType::New();
// diffusion->SetInput(gradienMagnitudeFilter->GetOutput() );
//  diffusion->SetNumberOfIterations(4 ); 
// diffusion->SetConductanceParameter( 1.0);
//  diffusion->SetTimeStep(0.125);

 // typedef itk::GradientMagnitudeImageFilter<UnsignedCharImageType, FloatImageType >  GradientMagnitudeImageFilterType;
 //typedef itk::GradientMagnitudeImageFilter<ImageType, FloatImageType >  GradientMagnitudeImageFilterType;
 
 //  GradientMagnitudeImageFilterType::Pointer gradientMagnitudeImageFilter = GradientMagnitudeImageFilterType::New();
 //  gradientMagnitudeImageFilter->SetInput(thresholdfilter->GetOutput()   );
 //  gradientMagnitudeImageFilter->Update();
 // // gradientMagnitudeImageFilter->SetUsePrincipleComponents(on);





 //  PerformSegmentation(gradientMagnitudeImageFilter->GetOutput(), threshold, level);
 
 //  // Fixed parameters

 //  PerformSegmentation(gradientMagnitudeImageFilter->GetOutput(), .00125, .125);
 //  PerformSegmentation(gradientMagnitudeImageFilter->GetOutput(), .0025, .25);
 //  PerformSegmentation(gradientMagnitudeImageFilter->GetOutput(), .005, .5);
 //  PerformSegmentation(gradientMagnitudeImageFilter->GetOutput(), .0075, .75);
 //  PerformSegmentation(gradientMagnitudeImageFilter->GetOutput(), .009, .9);

//

    _is_console << "                        ...filtering finished" << std::endl;

    _bounding_box.min() = _origin;
    _bounding_box.max() = _origin + _physical_dims;

    if (_is_pars.isParamValid("component"))
    {
      unsigned int n =  scaledImage->GetNumberOfComponentsPerPixel();
      std::cout<<"Number Of Components Per Pixel: "<< n <<std::endl;
     _component = _is_pars.get<unsigned int>("component");

      if (_component >= n)
      mooseError("'component' parameter must be empty or have a value of 0 to ", n - 1);
  
    }
    
    else
    _component = 0;
  }

  catch (itk::ExceptionObject &ex)
  { 
    mooseError("exception in reading dicom series");
  }

}

Real
ImageSamplerItk::sample(const Point & p)
{
  // Do nothing if the point is outside of the image domain
  if (!_bounding_box.contains_point(p))
    return 0.0;
  // Determine pixel coordinates
  std::vector<int> x(3, 0);
  for (int i = 0; i < LIBMESH_DIM; ++i)
  {
    // Compute position, only if voxel size is greater than zero
    if (_voxel[i] == 0)
     { x[i] = 0;   }
    else
    {
      x[i] = std::floor((p(i) - _origin(i)) / _voxel[i]);
       // If the point falls on the mesh extents the index needs to be decreased by one
      if (x[i] == imageSize[i])
        x[i]--;
    }
  }

  pixelIndex[0]=x[0];
  pixelIndex[1]=x[1];
  pixelIndex[2]=x[2];

  pixelValue = scaledImage->GetPixel(pixelIndex);
  
  return pixelValue;

}

void ImageSamplerItk::PerformSegmentation(ImageType::Pointer image, const float threshold, const float level)
{
  itk::TIFFImageIOFactory::RegisterOneFactory();

//  typedef itk::MorphologicalWatershedImageFilter<FloatImageType, LabeledImageType> MorphologicalWatershedFilterType;
//   MorphologicalWatershedFilterType::Pointer watershedFilter = MorphologicalWatershedFilterType::New();
//  // watershedFilter->SetThreshold(threshold);
//   watershedFilter->SetLevel(level);
//   watershedFilter->SetInput(image);  
//   watershedFilter->Update();


//   //  typedef itk:: WatershedImageFilter<FloatImageType> WatershedFilterType;
//   // WatershedFilterType::Pointer watershedFilter = WatershedFilterType::New();
//   // watershedFilter->SetThreshold(threshold);
//   // watershedFilter->SetLevel(level);
//   // watershedFilter->SetInput(image);  
//   // watershedFilter->Update();

// std::cout<<"sono arrivato H"<<std::endl;


//   typedef itk::ScalarToRGBColormapImageFilter<LabeledImageType, RGBImageType> RGBFilterType;
//   RGBFilterType::Pointer colormapImageFilter = RGBFilterType::New();
//   colormapImageFilter->SetInput(watershedFilter->GetOutput());
//   colormapImageFilter->SetColormap( RGBFilterType::Jet );
//   colormapImageFilter->Update();


// typedef itk::Functor::ScalarToRGBPixelFunctor<unsigned long> ColorMapFunctorType;
//   typedef itk::UnaryFunctorImageFilter<LabeledImageType, RGBImageType, ColorMapFunctorType> ColorMapFilterType;

// ColorMapFilterType::Pointer colormapper = ColorMapFilterType::New();
// colormapper->SetInput(watershedFilter->GetOutput());
//   colormapper->Update();



// std::cout<<"sono arrivato i"<<std::endl;

//  std::stringstream ss;
//   ss << "output_" << threshold << "_" << level << ".tiff";
 
//   typedef itk::ImageFileWriter<RGBImageType> FileWriterType2;
//   FileWriterType2::Pointer writer = FileWriterType2::New();
//   writer->SetFileName(ss.str());
//   writer->SetInput(colormapImageFilter->GetOutput());
//   writer->Update();

// std::cout<<"sono arrivato l"<<std::endl;
 
}
