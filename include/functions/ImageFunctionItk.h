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

#ifndef IMAGEFUNCTIONITK_H
#define IMAGEFUNCTIONITK_H

// MOOSE includes
#include "Function.h"
#include "ImageSamplerItk.h"

// Forward declarations
class ImageFunctionItk;

template <>
InputParameters validParams<ImageFunctionItk>();

/**
 * A function for extracting data from an image or stack of images
 */
class ImageFunctionItk : public ImageSamplerItk, public Function
{
public:
  /**
   * Class constructor
   * @param parameters The parameters object holding data for the class to use.
   */
  ImageFunctionItk(const InputParameters & parameters);

  /**
   * Class destructor
   */
  virtual ~ImageFunctionItk();

  /**
   * Initialize the ImageSampler
   */
  virtual void initialSetup() override;

  /**
   * Return the pixel value for the given point
   * @param t Time (unused)
   * @param p The point at which to extract pixel data
   */
  virtual Real value(Real t, const Point & p) override;
};

#endif // IMAGEFUNCTIONITK_H
