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
#include "itkSliceBySliceImageFilter.h"


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
{itk::TIFFImageIOFactory::RegisterOneFactory();


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
    const InternalImageType::SpacingType& inputSpacing=reader->GetOutput()->GetSpacing();
    _is_console <<"DICOM Serie Spacing:      " << inputSpacing << std::endl;
    const InternalImageType::PointType & origin = reader->GetOutput()->GetOrigin();
    _is_console  <<"DICOM Serie Origin:       "<<  origin << std::endl;

    _is_console << "                            ...image read!" << std::endl;
    _is_console << "Applying the filters...    " << std::endl<<std::endl;




 typedef itk::CastImageFilter< ShortImageType, InternalImageType > CastFilterType;

   CastFilterType::Pointer castFilter = CastFilterType::New();
   castFilter->SetInput( reader->GetOutput() );


 typedef itk::RescaleIntensityImageFilter<   InternalImageType, InternalImageType > RescaleFilterType;
  RescaleFilterType::Pointer rescaler = RescaleFilterType::New();
  rescaler->SetOutputMinimum(   0 );
  rescaler->SetOutputMaximum( 255 ); //255 means white
  rescaler->SetInput(castFilter->GetOutput() );



    scaledImage=rescaler->GetOutput();
    scaledImage->SetSpacing( inputSpacing );
    scaledImage->GetLargestPossibleRegion();
    scaledImage->Update(); 

    for (unsigned int i = 0; i < 3; ++i)
    {  
     _voxel.push_back(_physical_dims(i) / (imageSize[i]));
    }

   InternalImageType::SpacingType spacing2;
    
    spacing2[0] = _voxel[0]; // spacing along X
    spacing2[1] = _voxel[1];
    spacing2[2] = _voxel[2];
    
    scaledImage->SetSpacing(spacing2);

    InternalImageType::PointType newOrigin;
    
    newOrigin[0] = _origin(0); // spacing along X
    newOrigin[1] = _origin(1);
    newOrigin[2] = _origin(2);

    scaledImage->SetOrigin( newOrigin);       
    _is_console  <<"DICOM Serie New Origin:   " << newOrigin<<std::endl<<std::endl;



     


typedef itk::GrayscaleFillholeImageFilter<InternalImageType,InternalImageType>  FillholeFilterType;
FillholeFilterType::Pointer  fillhole = FillholeFilterType::New();



typedef itk::CurvatureFlowImageFilter< InternalImageType, InternalImageType > CurvatureFlowImageFilterType;
CurvatureFlowImageFilterType::Pointer smoothing = CurvatureFlowImageFilterType::New();

typedef itk::ConnectedThresholdImageFilter< InternalImageType, InternalImageType > ConnectedFilterType;
ConnectedFilterType::Pointer connectedThreshold = ConnectedFilterType::New();

 typedef itk::CastImageFilter< ShortImageType, OutputImageType > CastFilterType2;
CastFilterType2::Pointer caster3 = CastFilterType2::New();
CastFilterType2::Pointer caster2 = CastFilterType2::New();

 typedef itk::CastImageFilter< InternalImageType, OutputImageType > CastFilterType3;
CastFilterType3::Pointer caster4 = CastFilterType3::New();



// smoothing->SetNumberOfIterations( 0 ); smoothing->SetTimeStep( 0.125 );

// smoothing->SetInput( reader->GetOutput() );



   InternalImageType::SizeType size2 =castFilter->GetOutput()->GetLargestPossibleRegion().GetSize();

  InternalImageType::IndexType index2 =castFilter->GetOutput()->GetLargestPossibleRegion().GetIndex();

  InternalImageType::PixelType maxValue;
  InternalImageType::PixelType minValue;

             for(int i = index2[0]; i<size2[0]; i++)

             {

                    for(int j = index2[1]; j<size2[1]; j++)

                    {

                           InternalImageType::IndexType currentIndex;

                           currentIndex[0] = i;

                           currentIndex[1] = j;

                           InternalImageType::PixelType currentValue =castFilter->GetOutput()->GetPixel(currentIndex);

                           if(currentValue>maxValue)

                           {


//std::cout<<"i: "<<i<<"j: "<<j<<std::endl;
                                  maxValue = currentValue;

                           }

                           if(currentValue<minValue)

                           {

                                  minValue = currentValue;

                           }

                    }

             }


             std::cout<<"maxValue "<<maxValue<<std::endl;
             std::cout<<"minValue "<<minValue<<std::endl;

// connectedThreshold->SetInput(scaledImage); 
// connectedThreshold->SetLower(80 ); 
// connectedThreshold->SetUpper( 130 );
// connectedThreshold->SetReplaceValue( 255 ); //255=white means selected
// connectedThreshold->SetInput( scaledImage ); 


connectedThreshold->SetInput(scaledImage); 
connectedThreshold->SetLower(90 ); 
connectedThreshold->SetUpper( 112 );
connectedThreshold->SetReplaceValue( 255 ); //255=white means selected
connectedThreshold->SetInput( scaledImage ); 



InternalImageType::IndexType  index;


// //shoulder 
//    index[0]=142;
//    index[1]=138;
//    index[2]=4;


// //shoulder from above
//    index[0]=278;
//    index[1]=230;
//    index[2]=4;

//brain 1 
  // index[0]=101;
  // index[1]=177;
  // index[2]=4;


//brain 2 
  index[0]=60;
  index[1]=116;
  index[2]=8;

std::cout << "Seed index = " << index << std::endl;
connectedThreshold->SetSeed(index);
connectedThreshold->Update();


caster4->SetInput(connectedThreshold->GetOutput());
caster4->Update();

imageSize =connectedThreshold->GetOutput()->GetLargestPossibleRegion().GetSize();   
std::cout<<imageSize<<": that one was ImageSize"<<std::endl;

unsigned char  pixel_value; 
pixel_value= scaledImage->GetPixel( index ); 

std::cout<<"pixel_value : "<<static_cast<unsigned>(pixel_value) <<std::endl;

pixel_value= connectedThreshold->GetOutput()->GetPixel( index ); 

std::cout<<"pixel_value : "<<static_cast<unsigned>(pixel_value) <<std::endl;

caster4->SetInput(connectedThreshold->GetOutput());
caster4->Update();






// typedef itk::GradientMagnitudeImageFilter<InternalImageType, InternalImageType >  GradientMagnitudeImageFilterType;


//   itk::TIFFImageIOFactory::RegisterOneFactory();
//   typedef itk::ImageFileWriter< InternalImageType >  Writer2Type2;
//   Writer2Type2::Pointer writer2 = Writer2Type2::New();

//  writer2->SetFileName("outputFilenameFloat2.tiff" );
//      writer2->SetInput( connectedThreshold->GetOutput() );
 
//       try
//       {
//       writer2->Update();
//       }
//    catch (itk::ExceptionObject & e)
//       {
//       std::cerr << e << std::endl;
//       mooseError("Exception in file writer");
//       }

  typedef itk::ImageFileWriter< OutputImageType >  WriterType;
  WriterType::Pointer writer3 = WriterType::New();

    writer3->SetFileName("outputFilename.tiff" );
     writer3->SetInput( caster4->GetOutput() );
 
      try
      {
      writer3->Update();
      }
   catch (itk::ExceptionObject & e)
      {
      std::cerr << e << std::endl;
      mooseError("Exception in file writer");
      }

  
    _is_console << "                        ...filtering finished" << std::endl;

    _bounding_box.min() = _origin;
    _bounding_box.max() = _origin + _physical_dims;



    _is_console << "                        ...filtering finished 2" << std::endl;

    if (_is_pars.isParamValid("component"))
    {

      unsigned int n =  scaledImage->GetNumberOfComponentsPerPixel();
      std::cout<<"Number Of Components Per Pixel: "<< n <<std::endl;
     _component = _is_pars.get<unsigned int>("component");

   _is_console << "                        ...filtering finished" << std::endl;
      if (_component >= n)
      mooseError("'component' parameter must be empty or have a value of 0 to ", n - 1);
  
    }
    
    else{
    _component = 0;}
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

