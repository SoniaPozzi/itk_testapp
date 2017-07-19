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

// Software Guide : BeginCodeSnippet
#include "itkImageMaskSpatialObject.h"
// Software Guide : EndCodeSnippet
#include "itkImageRegionIterator.h"
#include "itkAffineTransform.h"

#include "itkImageRegionIterator.h"
#include "itkBinaryImageToShapeLabelMapFilter.h"
#include "itkBinaryImageToLabelMapFilter.h"
#include "itkObjectByObjectLabelMapFilter.h"
#include "itkLabelMapToLabelImageFilter.h"
#include "itkLabelImageToLabelMapFilter.h"
#include "itkScalarToRGBColormapImageFilter.h"
#include "itkConnectedComponentImageFilter.h"
#include "itkImageRegionIterator.h"
#include "itkBinaryImageToLabelMapFilter.h"

#include "itkImageFileWriter.h"
#include "itkImageRegionIterator.h"
#include "itkBinaryImageToShapeLabelMapFilter.h"

// Software Guide : BeginCodeSnippet
#include "itkImageFileReader.h"
#include "itkImageFileWriter.h"
// Software Guide : EndCodeSnippet
//  Software Guide : BeginLatex
//
//  The filter used to extract a region from an image is the
//  \doxygen{ExtractImageFilter}. Its header is included below.  This filter
//  is capable of extracting a slice from the input image.
//
//  \index{itk::ExtractImageFilter!header}
//
//  Software Guide : EndLatex
// Software Guide : BeginCodeSnippet
#include "itkExtractImageFilter.h"
// Software Guide : EndCodeSnippet
//  Software Guide : BeginLatex
//
//  The filter used to place the processed image in a region of the output
//  image is the \doxygen{PasteImageFilter}. Its header is included below.
//  This filter is capable of inserting the processed image into the
//  destination image.
//
//  \index{itk::PasteImageFilter!header}
//
//  Software Guide : EndLatex
// Software Guide : BeginCodeSnippet
#include "itkPasteImageFilter.h"
// Software Guide : EndCodeSnippet
// Software Guide : BeginCodeSnippet
#include "itkMedianImageFilter.h"

#include "itkImageFileReader.h"
#include "itkImageFileWriter.h"
#include "itkSimpleFilterWatcher.h"

#include "itkLabelImageToLabelMapFilter.h"
#include "itkLabelSelectionLabelMapFilter.h"
#include "itkAutoCropLabelMapFilter.h"
#include "itkLabelMapToLabelImageFilter.h"
#include "itkImageFileReader.h"
#include "itkImageFileWriter.h"
#include "itkSimpleFilterWatcher.h"

#include "itkLabelImageToLabelMapFilter.h"
#include "itkAutoCropLabelMapFilter.h"
#include "itkLabelMapToLabelImageFilter.h"

#include "itkImageFileWriter.h"
#include "itkImageRegionIterator.h"
#include "itkBinaryImageToLabelMapFilter.h"
#include "itkLabelMapToLabelImageFilter.h"
#include "itkLabelStatisticsImageFilter.h"



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
  itk::TIFFImageIOFactory::RegisterOneFactory();

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


//////////////////////


for (unsigned int i = 0; i < 3; ++i)
{  
 _voxel.push_back(outputImageSpacing[i]);
}

filteredImage->Update();

OutputImageType::SpacingType spacing;

spacing[0] = _voxel[0]; // spacing along X
spacing[1] = _voxel[1];
spacing[2] = _voxel[2];

//filteredImage->SetSpacing(spacing);
filteredImage->Update();

OutputImageType::PointType newOrigin;

newOrigin[0] = _origin(0); 
newOrigin[1] = _origin(1);
newOrigin[2] = _origin(2);

filteredImage->SetOrigin(newOrigin);
filteredImage->Update();


const OutputImageType::PointType & croppedorigin = filteredImage->GetOrigin();
_is_console  <<"   Filterd DICOM Serie Origin moved to:    "<<  croppedorigin << std::endl;
_is_console << "                        ...filtering finished" << std::endl;

_bounding_box.min() = _origin;
_bounding_box.max() = _origin + _physical_dims;

    if (_is_pars.isParamValid("component"))
    {

      unsigned int n =  filteredImage->GetNumberOfComponentsPerPixel();
      std::cout<<"Number Of Components Per Pixel: "<< n <<std::endl;
     _component = _is_pars.get<unsigned int>("component");

   _is_console << "                        ...filtering finished" << std::endl;
      if (_component >= n)
      mooseError("'component' parameter must be empty or have a value of 0 to ", n - 1);
  
    }
    
    else{
    _component = 0;}


}

Real
ImageSamplerItk::sample(const Point & p)
{


 //   _is_console <<"   Filterd DICOM Serie Spacing:             " << croppedSpacing << std::endl;
  if (!_bounding_box.contains_point(p))
    return 0.0;
  // Determine pixel coordinates
  std::vector<int> x(3, 0);
//  std::vector<int> res(3, 0);
//  res[0]=55;
//  res[1]=16;
//  res[2]=0;
  for (int i = 0; i < LIBMESH_DIM; ++i)
  {
    // Compute position, only if voxel size is greater than zero
    if (_voxel[i] == 0)
     { x[i] = 0;   }
    else
    {
      x[i] = std::floor((p(i) - _origin(i)) / _voxel[i]);
       // If the point falls on the mesh extents the index needs to be decreased by one
      if (x[i] >= outputImageSize[i])
        x[i]--;
    }
  }

 OutputImageType::IndexType pixelIndex;
 OutputPixelType pixelValue;

  pixelIndex[0]=x[0];
  pixelIndex[1]=x[1];
  pixelIndex[2]=x[2];

  pixelValue = filteredImage->GetPixel(pixelIndex);

    return pixelValue;

}

