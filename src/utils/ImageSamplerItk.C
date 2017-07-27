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
#include "itkImage.h"


template <>
InputParameters
validParams<ImageSamplerItk>()
{
  // Define the general parameters
  InputParameters params = emptyInputParameters();
  params += validParams<FileDicomChoose>();
  params.addParam<Point>("origin", "Origin of the image (defaults to mesh origin)");
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

    for (unsigned int i = 0; i < 3; ++i)
    {  
     _voxel.push_back(  _physical_dims(i) / ( outputImageSize[i]) );
    }

    filteredImage -> Update();

    OutputImageType::SpacingType spacing;

    spacing[0] = _voxel[0]; // spacing along X
    spacing[1] = _voxel[1];
    spacing[2] = _voxel[2];

    filteredImage -> Update();

    OutputImageType::PointType newOrigin;

    newOrigin[0] = _origin(0); 
    newOrigin[1] = _origin(1);
    newOrigin[2] = _origin(2);

    filteredImage -> SetOrigin( newOrigin );
    filteredImage -> Update();

    _is_console << "Filtered and Cropped DICOM Serie Origin moved to: " << filteredImage -> GetOrigin() << std::endl;

    _bounding_box.min() = _origin;
    _bounding_box.max() = _origin + _physical_dims;
    
}

Real
ImageSamplerItk::sample(const Point & p)
{
      if (!_bounding_box.contains_point(p))
      return 0.0;


      std::cout<<"============"<<std::endl;
      typedef itk::Point< double, 3 > PointType;

    //  OutputImageType::IndexType pixel;
     itk::ContinuousIndex<double, 3> pixel;
      OutputPixelType pixelValue;

      std::vector<double> x(3, 0);
      for (int i = 0; i < LIBMESH_DIM; ++i)
      {
          if (_voxel[i] == 0)
          { pixel[i] = 0;   }
          else
          {
          x[i] = ((p(i) - _origin(i))/(  _voxel[i] ));
          if (x[i] == outputImageSize[i])
          x[i]--;

          }
      }

      pixel[0]=x[0];
      pixel[1]=x[1];
      pixel[2]=x[2];


      itk::LinearInterpolateImageFunction<OutputImageType, double>::Pointer linearInterpolator =
      itk::LinearInterpolateImageFunction<OutputImageType, double>::New();
      linearInterpolator->SetInputImage(smoothedImage);

      //std::cout << "Value at 1.3: " << linearInterpolator->EvaluateAtContinuousIndex(pixel) << std::endl;

      pixelValue=  linearInterpolator->EvaluateAtContinuousIndex(pixel) ;

    //pixelValue=  smoothedImage->GetPixel(pixel) ;

      std::cout << pixel<<"Value at: " << static_cast<unsigned> (pixelValue) << std::endl;




    return pixelValue;

  }

