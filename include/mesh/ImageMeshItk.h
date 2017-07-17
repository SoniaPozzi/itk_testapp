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

#ifndef IMAGEMESHITK_H
#define IMAGEMESHITK_H

#include "GeneratedMesh.h"
#include "FileDicomChoose.h"

#include "itkImageFileWriter.h"
#include "itkImageSeriesWriter.h"
#include "itkRescaleIntensityImageFilter.h"

#include "itkImageSeriesWriter.h"

#include "itkImageFileReader.h"
#include "itkImageFileWriter.h"
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





#include "itkGrayscaleFillholeImageFilter.h"

#include "itkCannyEdgeDetectionImageFilter.h"  
#include "itkCastImageFilter.h"
#include "itkRescaleIntensityImageFilter.h"

#include "itkThresholdImageFilter.h"
#include "itkConnectedThresholdImageFilter.h"
#include "itkCurvatureFlowImageFilter.h"


#include "itkLabelImageToLabelMapFilter.h"
#include "itkLabelSelectionLabelMapFilter.h"
#include "itkAutoCropLabelMapFilter.h"
#include "itkLabelMapToLabelImageFilter.h"
#include "itkSimpleFilterWatcher.h"

class ImageMeshItk;

template <>
InputParameters validParams<ImageMeshItk>();

/**
 * A 2D GeneratedMesh where xmin, xmax, etc. are determined from an input image file.
 */
class ImageMeshItk : public GeneratedMesh,  public FileDicomChoose
{
  public:
  ImageMeshItk(const InputParameters & parameters);
  ImageMeshItk(const ImageMeshItk & other_mesh);

  virtual MooseMesh & clone() const override;
  virtual void buildMesh() override;



  protected:

  /**
   * buildMesh() calls this helper function to build 3D ImageMeshes from stacks of images.
   */
  void buildMesh3D();

  const bool _scale_to_one;
  bool _has_chooseSeries;


std::vector<Real>  _cells_per_pixel_vector;

typedef unsigned char                      OutputPixelType;
typedef itk::Image<OutputPixelType, 3 >    OutputImageType;
OutputImageType::Pointer scaledImage2=OutputImageType::New();
typedef itk::Image<itk::IdentifierType, 3> LabeledImageType;




};

#endif /* IMAGEMESHITK_H */