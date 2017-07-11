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

#include "ImageFunctionItk.h"
#include "MooseMesh.h"

template <>
InputParameters
validParams<ImageFunctionItk>()
{
  // Define the general parameters
  InputParameters params = validParams<Function>();
  params += validParams<ImageSamplerItk>();
  params.addClassDescription("Function with values sampled from a given image stack");
  return params;
}     

ImageFunctionItk::ImageFunctionItk(const InputParameters & parameters)
  : ImageSamplerItk(parameters), Function(parameters)
{
}


ImageFunctionItk::~ImageFunctionItk() {}  

void
ImageFunctionItk::initialSetup()
{

  std::cout<<"ImageFunctionItk::initialSetup"<<std::endl;
  FEProblemBase * fe_problem = this->getParam<FEProblemBase *>("_fe_problem_base");
  MooseMesh & mesh = fe_problem->mesh();
  setupImageSampler(mesh);
}

Real
ImageFunctionItk::value(Real /*t*/, const Point & p)
{

  return sample(p);
}
