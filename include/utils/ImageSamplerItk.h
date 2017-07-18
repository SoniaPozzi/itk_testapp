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

#ifndef IMAGESAMPLERITK_H
#define IMAGESAMPLERITK_H

// MOOSE includes
#include "FileDicomChoose.h"
#include "ConsoleStream.h"  

// libmesh includes
#include "libmesh/mesh_tools.h"

// VTK includes
#ifdef LIBMESH_HAVE_VTK

// Some VTK header files have extra semi-colons in them, and clang
// loves to warn about it...
#include "libmesh/ignore_warnings.h"

#endif


#include "itkScalarToRGBPixelFunctor.h"
#include "itkUnaryFunctorImageFilter.h"
#include "itkVectorCastImageFilter.h"
#include "itkVectorGradientAnisotropicDiffusionImageFilter.h"
#include "itkWatershedImageFilter.h"
#include "itkRescaleIntensityImageFilter.h"
#include "itkScalarToRGBColormapImageFilter.h"
#include "itkGradientMagnitudeImageFilter.h"
#include "itkGradientAnisotropicDiffusionImageFilter.h"
#include "itkVectorMagnitudeImageFilter.h"
#include "itkMorphologicalWatershedImageFilter.h"
 #include "itkTIFFImageIOFactory.h"

#include "itkGradientMagnitudeRecursiveGaussianImageFilter.h"
#include <itkThresholdImageFilter.h>



// Forward declarations
class ImageSamplerItk;
class MooseMesh;

template <>
InputParameters validParams<ImageSamplerItk>();

/**
 * A helper class for reading and sampling images using VTK.
 */
class ImageSamplerItk : public FileDicomChoose
{
public:
  /**
   * Constructor.
   *
   * Use this object as an interface, being sure to also add the parameters to the
   * child class.
   *
   * @see ImageFunction
   */
  ImageSamplerItk(const InputParameters & parameters);

  /**
   * Return the pixel value for the given point
   * @param p The point at which to extract pixel data
   */
  virtual Real sample(const Point & p);
 
  /**
   * Perform initialization of image data
   */
  virtual void setupImageSampler(MooseMesh & mesh);

protected:
  /**
   * Apply image re-scaling using the vtkImageShiftAndRescale object
   */
  void  ItkImageSampler(MooseMesh & mesh);




private:
  /// Origin of image
  Point _origin;

  /// Pixel dimension of image
  std::vector<int> _dims;

  /// Physical dimensions of image
  Point _physical_dims;

  /// Physical pixel size
  std::vector<double> _voxel;
  std::vector<double> _ratios;

  /// Component to extract

  unsigned int _component;

  /// Bounding box for testing points
  MeshTools::BoundingBox _bounding_box;

  /// Parameters for interface
  const InputParameters & _is_pars;

  /// Create a console stream object for this helper class
  ConsoleStream _is_console;

protected:





};




#endif // IMAGESAMPLERDICOM_H
