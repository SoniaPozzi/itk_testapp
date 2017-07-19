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

#include <cstdlib> 
#include "ImageMeshItk.h"
#include "MooseApp.h"

// ITK headers
#include <itkPoint.h>

// libMesh includes
#include "libmesh/mesh_generation.h"



template <>
InputParameters
validParams<ImageMeshItk>()
{
  InputParameters params = validParams<GeneratedMesh>();
  params += validParams<FileDicomChoose>();
  params.addRequiredParam<std::string>("dicomDirectory", "Name of single image file to extract mesh parameters from.  If provided, a 2D mesh is created.");
  params.addRequiredParam<std::string>("file_base", "A filename in the DICOM series to use.");
  params.addClassDescription("Generated mesh with the aspect ratio of a given image stack.");
  // Add ImageMesh-specific params
  params.addParam<bool>(
      "scale_to_one", true, "Whether or not to scale the image so its max dimension is 1");
  params.addRequiredParam<std::vector<Real>>("cells_per_pixel_vector","Mulitplier of image dimension to find nx, ny nz.");

  return params;
}

ImageMeshItk::ImageMeshItk(const InputParameters & params)
  : GeneratedMesh(params),
  FileDicomChoose(params),
    _scale_to_one(getParam<bool>("scale_to_one")),
    _cells_per_pixel_vector(getParam<std::vector<Real>>("cells_per_pixel_vector"))
{
    _dicomDirectory = params.get<std::string>("dicomDirectory");
    _file_base = params.get<std::string>("file_base");
}
  // Failure is not an option

ImageMeshItk::ImageMeshItk(const ImageMeshItk & other_mesh)
  : GeneratedMesh(other_mesh),
  FileDicomChoose(other_mesh.parameters()),
    _scale_to_one(getParam<bool>("scale_to_one")),
    _cells_per_pixel_vector(getParam<std::vector<Real>>("cells_per_pixel_vector"))
{
}

MooseMesh &
ImageMeshItk::clone() const
{
  return *(new ImageMeshItk(*this));
}

void
ImageMeshItk::buildMesh()
{
    buildMesh3D();
}

void
ImageMeshItk::buildMesh3D()
{


  int xpixels = 0, ypixels = 0, zpixels = 0;


  xpixels = outputImageSize[0], ypixels = outputImageSize[1], zpixels =outputImageSize[2];

  _xmax = outputImageSize[0]*outputImageSpacing[0];
  _ymax = outputImageSize[1]*outputImageSpacing[1];
  _zmax = outputImageSize[2]*outputImageSpacing[2];

  if (_scale_to_one)
  {
    Real max = std::max(std::max(_xmax, _ymax), _zmax);
    _xmax /= max;
    _ymax /= max;
    _zmax /= max;
   
  }

  std::cout<<"xmax"<<_xmax<<_ymax<<_zmax<<std::endl;
  // Compute the number of cells in the x and y direction based on
  // the user's cells_per_pixel parameter.  Note: we use ints here
  // because the GeneratedMesh params object uses ints for these...
  _nx = static_cast<int>(_cells_per_pixel_vector[0] * xpixels);
  _ny = static_cast<int>(_cells_per_pixel_vector[1] * ypixels);
  _nz = static_cast<int>(_cells_per_pixel_vector[2]* zpixels);

  std::cout<<"Mesh dimensions: nx "<<_nx<<", ny "<<_ny<<", nz "<<_nz<<std::endl;

  // Actually build the Mesh
  MeshTools::Generation::build_cube(dynamic_cast<UnstructuredMesh &>(getMesh()),
                                    _nx,
                                    _ny,
                                    _nz,
                                    /*xmin=*/0.,
                                    /*xmax=*/_xmax,
                                    /*ymin=*/0.,
                                    /*ymax=*/_ymax,
                                    /*zmin=*/0.,
                                    /*zmax=*/_zmax,
                                    HEX8);
   
}




