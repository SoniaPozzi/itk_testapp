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
#include <iostream>

// ITK headers
#include "itkGDCMImageIO.h"
#include "itkGDCMSeriesFileNames.h"
#include "itkTIFFImageIO.h"


#include <itkPoint.h>



template <>
InputParameters
validParams<ImageSamplerItk>()
{
  // Define the general parameters
  InputParameters params = emptyInputParameters();
  params += validParams<FileDicomChoose>();

  params.addParam<Point>("origin", "Origin of the image (defaults to mesh origin)");
  params.addParam<Point>("dimensions",
                         "x,y,z dimensions of the image (defaults to mesh dimensions)");
  params.addParam<unsigned int>(
      "component",
      "The image RGB-component to return, leaving this blank will result in a greyscale value "
      "for the image to be created. The component number is zero based, i.e. 0 returns the first "
      "(RED) component of the image.");

  // Shift and Scale (application of these occurs prior to threshold)
  params.addParam<double>("shift", 0, "Value to add to all pixels; occurs prior to scaling");
  params.addParam<double>(
      "scale", 1, "Multiplier to apply to all pixel values; occurs after shifting");
  params.addParamNamesToGroup("shift scale", "Rescale");

  // Threshold parameters
  params.addParam<double>("threshold", "The threshold value");
  params.addParam<double>(
      "upper_value", 1, "The value to set for data greater than the threshold value");
  params.addParam<double>(
      "lower_value", 0, "The value to set for data less than the threshold value");
  params.addParamNamesToGroup("threshold upper_value lower_value", "Threshold");

  // Flip image
  params.addParam<bool>("flip_x", false, "Flip the image along the x-axis");
  params.addParam<bool>("flip_y", false, "Flip the image along the y-axis");
  params.addParam<bool>("flip_z", false, "Flip the image along the z-axis");
  params.addParamNamesToGroup("flip_x flip_y flip_z", "Flip");

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

  itk::TIFFImageIOFactory::RegisterOneFactory();

  std::cout<<"Reading DICOM using ITK:"<<std::endl;

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
      mooseError("No file range parameters were provided and the Mesh is not an ImageMesh.");

    // Get the ImageMesh's parameters.  This should work, otherwise
    // errors would already have been thrown...
    //filenames = image_mesh->filenames();
    //file_suffix = image_mesh->fileSuffix();
  }
  else
  {
    // Use our own parameters (using 'this' b/c of conflicts with filenames the local variable)
   // filenames = this->filenames();
    //file_suffix = fileSuffix();
  }


  //  reader->SetImageIO( dicomIO );
  //  nameGenerator->SetUseSeriesDetails( true );
  //  nameGenerator->SetDirectory("/Users/sonia/Projects/internshipINL/itk_testapp/tests/image_function/new_stack/");

   try{

  //   std::cout <<"The directory: " << std::endl;
  //   std::cout << "/Users/sonia/Projects/internshipINL/itk_testapp/tests/image_function/new_stack/"  << std::endl;
  //   std::cout << "Contains the following DICOM Series: "<< std::endl;

  //  typedef std::vector< std::string >    SeriesIdContainer;
  //   const SeriesIdContainer & seriesUID = nameGenerator->GetSeriesUIDs();
  //   SeriesIdContainer::const_iterator seriesItr = seriesUID.begin();
  //   SeriesIdContainer::const_iterator seriesEnd = seriesUID.end();
  //   while( seriesItr != seriesEnd )
  //     {
  //     std::cout << seriesItr->c_str() << std::endl;
  //     ++seriesItr;
  //     }


  //  std::string seriesIdentifier;
  // if( 0) // If no optional series identifier
  //   {
  //   seriesIdentifier = "se definiamo da fuori la serie da leggere";
  //   }
  // else
  //   {
  //   seriesIdentifier = seriesUID.begin()->c_str();
  //   }

  // std::cout << "Now reading series: " << std::endl << std::endl;
  // std::cout << seriesIdentifier << std::endl << std::endl;
 
  // fileNames = nameGenerator->GetFileNames( seriesIdentifier );
  // reader->SetFileNames( fileNames );
std::cout<<"sono quo"<<std::endl;
  // try
  //   {
  //   reader->Update();
  //   reader->GetOutput();
  //   }
  // catch(itk::ExceptionObject &ex)
  //   {
  //   std::cout << ex << std::endl;
  //   }

  imageSize =reader->GetOutput()->GetLargestPossibleRegion().GetSize();
  std::cout<<"DICOM Serie Dimension:    " <<imageSize<<std::endl;
  const ImageType::SpacingType& inputSpacing=reader->GetOutput()->GetSpacing();
  std::cout << "DICOM Serie Spacing:    " << inputSpacing << std::endl;
  const ImageType::PointType & origin = reader->GetOutput()->GetOrigin();
  std::cout << "DICOM Serie Origin:    "<<  origin << std::endl;

  _is_console << "          ...image read finished" << std::endl;
  _is_console << "Applying the filters...    " << std::endl;

  rescaler = RescaleFilterType::New();
  rescaler->SetOutputMinimum(   0 );
  rescaler->SetOutputMaximum( 255 );
  rescaler->SetInput( reader->GetOutput() );
  rescaler->Update();



  std::cout << "Scaling Voxels Ratio    "<<  std::endl; //_voxel <<

  scaledImage=(rescaler->GetOutput());
  scaledImage->SetSpacing( inputSpacing );
  scaledImage->GetLargestPossibleRegion();


std::cout << "Scaling Voxels Ratio    "<<  std::endl; //_voxel <<


std::cout<<"phy prima"<<_physical_dims(2)<<std::endl;
//_physical_dims(2)=_physical_dims(2)*imageSize[2]*inputSpacing[2]/(imageSize[1]*inputSpacing[1]);
std::cout<<"phy dopo"<<_physical_dims(2)<<std::endl;
 for (unsigned int i = 0; i < 3; ++i)
  {   _voxel.push_back(_physical_dims(i) / (imageSize[i]));
   }

  std::cout << "Scaling Voxels Ratio    "<<  std::endl; //_voxel << 

  WriteImageType::SpacingType spacing2;
  spacing2[0] = _voxel[0]; // spacing along X
  spacing2[1] = _voxel[1];
  spacing2[2] = _voxel[2];
  scaledImage->SetSpacing(spacing2);


  WriteImageType::PointType newOrigin;
  newOrigin[0] = _origin(0); // spacing along X
  newOrigin[1] = _origin(1);
  newOrigin[2] = _origin(2);

  scaledImage->SetOrigin( newOrigin);
  std::cout << "New Origin = "<<  newOrigin << std::endl;

  _is_console << "          ...filtering finished" << std::endl;

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

catch (itk::ExceptionObject &ex){ 
  mooseError("exception in reading dicom series");
}

}

Real
ImageSamplerItk::sample(const Point & p)
{
    // writer->SetFileName("outputFilenameFiltred.tiff" );
    // writer->SetInput( scaledImage );

    //  try
    //  {
    //  writer->Update();
    //  }
    // catch (itk::ExceptionObject & e)
    //  {
    //  std::cerr << e << std::endl;
    //  mooseError("Exception in file writer");
    //  }

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
 
  typedef itk::Point< double, WriteImageType::ImageDimension > PointType;
  PointType coordinate;

  pixelIndex[0]=x[0];
  pixelIndex[1]=x[1];
  pixelIndex[2]=x[2];

   // std::cout<<"pixel index"<<pixelIndex<<std::endl;
    pixelValue = scaledImage->GetPixel(pixelIndex);
   //std::cout<<"pixel Values"<<static_cast<unsigned int>(pixelValue)<<std::endl;
  return pixelValue;

}

