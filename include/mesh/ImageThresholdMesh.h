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

#ifndef IMAGETHRESHOLDMESH_H
#define IMAGETHRESHOLDMESH_H

// MOOSE includes
#include "ElementDeleterBase.h"
#include "ImageSamplerItk.h"

// Forward declarations
class ImageThresholdMesh;

template <>
InputParameters validParams<ImageThresholdMesh>();

/**
 * MeshModifier for defining a Subdomains based on Image data.
 */
class ImageThresholdMesh : public ElementDeleterBase, public ImageSamplerItk
{
public:
  ImageThresholdMesh(const InputParameters & parameters);


  unsigned int _lower_threshold;
  unsigned int _upper_threshold;

protected:
	
  virtual void modify() override;
  virtual bool shouldDelete(const Elem * elem) override;
};

#endif // IMAGETHRESHOLDMESH_H
