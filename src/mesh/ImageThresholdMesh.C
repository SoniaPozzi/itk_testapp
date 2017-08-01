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

// STL includes
#include <cmath> // provides round, not std::round (see http://www.cplusplus.com/reference/cmath/round/)

// MOOSE includes
#include "ImageThresholdMesh.h"
#include "MooseMesh.h"


template <>
InputParameters
validParams<ImageThresholdMesh>()
{
  InputParameters params = validParams<ElementDeleterBase>();
  params += validParams<ImageSamplerItk>();
  params.addRequiredParam<unsigned>("lower_threshold", "Image lower threshold to create the mesh");
  params.addRequiredParam<unsigned>("upper_threshold", "Image upper threshold to create the mesh");
  return params;
}

ImageThresholdMesh::ImageThresholdMesh(const InputParameters & parameters)
  : ElementDeleterBase(parameters), ImageSamplerItk(parameters),
   _lower_threshold(getParam<unsigned>("lower_threshold")),
   _upper_threshold(getParam<unsigned>("upper_threshold"))
{
}


void
ImageThresholdMesh::modify()
{
  // Check that we have access to the mesh
  if (!_mesh_ptr)
    mooseError("_mesh_ptr must be initialized before calling SubdomainBoundingBox::modify()");

  // Initialize the ImageSampler
  setupImageSampler(*_mesh_ptr);
  ElementDeleterBase::modify();

}


bool
ImageThresholdMesh::shouldDelete(const Elem * elem)
{ 
  bool outside=0;
  if (round(ImageSamplerItk::sample((*elem).centroid()))< _lower_threshold || round(ImageSamplerItk::sample((*elem).centroid()))> _upper_threshold  ) 
    outside=1;

  return outside  == 1;
}



