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
#include "FileRangeBuilder.h"
#include "ConsoleStream.h"  

// libmesh includes
#include "libmesh/mesh_tools.h"

// VTK includes
#ifdef LIBMESH_HAVE_VTK

// Some VTK header files have extra semi-colons in them, and clang
// loves to warn about it...
#include "libmesh/ignore_warnings.h"

#include "vtkSmartPointer.h"
#include "vtkDICOMImageReader.h"
#include "vtkImageData.h"
#include "vtkStringArray.h"
#include "vtkImageThreshold.h"
#include "vtkImageNormalize.h"
#include "vtkImageCast.h"
#include "vtkImageShiftScale.h"
#include "vtkImageMagnitude.h"
#include "vtkImageFlip.h"


#include "vtkRenderWindow.h"
#include "vtkRenderWindowInteractor.h"
#include "vtkRenderer.h"
#include "vtkImageViewer2.h"


#include "libmesh/restore_warnings.h"


#include "itkImage.h"
#include "itkGDCMImageIO.h"
#include "itkGDCMSeriesFileNames.h"
#include "itkImageSeriesReader.h"
#include "itkImageFileWriter.h"
#include "itkImageSeriesWriter.h"
#include "itkImageIOBase.h"

#include "itkRescaleIntensityImageFilter.h"



// Software Guide : EndCodeSnippet

#endif

// Forward declarations
class ImageSamplerItk;
class MooseMesh;

template <>
InputParameters validParams<ImageSamplerItk>();

/**
 * A helper class for reading and sampling images using VTK.
 */
class ImageSamplerItk : public FileRangeBuilder
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
  void vtkMagnitude();


private:
#ifdef LIBMESH_HAVE_VTK

  /// List of file names to extract data
  vtkSmartPointer<vtkStringArray> _files;

  /// Complete image data
  vtkImageData * _data;

  /// VTK-6 seems to work better in terms of "algorithm outputs" rather than vtkImageData pointers...
  vtkAlgorithmOutput * _algorithm;

  /// Complete image data
  vtkSmartPointer<vtkDICOMImageReader> _image;

  //vtkSmartPointer<vtkDICOMFileSorter> _sorter;

 vtkSmartPointer<vtkImageViewer2>  _viewer;

  /// Pointer to thresholding filter
  vtkSmartPointer<vtkImageThreshold> _image_threshold;

  /// Pointer to the shift and scaling filter
  vtkSmartPointer<vtkImageShiftScale> _shift_scale_filter;

  /// Pointer to the magnitude filter
  vtkSmartPointer<vtkImageMagnitude> _magnitude_filter;

  /// Pointers to image flipping filter.  May be used for x, y, or z.
  vtkSmartPointer<vtkImageFlip> _flip_filter;
#endif

  

/**
 * Helper method for flipping image
 * @param axis Flag for determing the flip axis: "x=0", "y=1", "z=2"
 * @return A smart pointer the flipping filter
 */
#ifdef LIBMESH_HAVE_VTK
  vtkSmartPointer<vtkImageFlip> imageFlip(const int & axis);
#endif

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
#ifdef LIBMESH_HAVE_VTK
  unsigned int _component;
#endif



  /// Bounding box for testing points
  MeshTools::BoundingBox _bounding_box;

  /// Parameters for interface
  const InputParameters & _is_pars;

  /// Create a console stream object for this helper class
  ConsoleStream _is_console;


protected:
 

  typedef short    PixelType;
  typedef unsigned char WritePixelType;
  const unsigned int      Dimension = 3;
  typedef std::vector< std::string >   FileNamesContainer;

  typedef itk::Image<PixelType, 3>     ImageType;
  typedef itk::ImageSeriesReader<ImageType >        ReaderType;

  typedef itk::Image< WritePixelType, 3 > WriteImageType;
  typedef itk::ImageFileWriter< WriteImageType >  Writer2Type;
  typedef itk::RescaleIntensityImageFilter< ImageType, WriteImageType > RescaleFilterType;
  typedef itk::GDCMImageIO     ImageIOType;
  typedef itk::GDCMSeriesFileNames NamesGeneratorType;

  ReaderType::Pointer reader = ReaderType::New();
  ImageType::SizeType imageSize;
  FileNamesContainer fileNames;
  RescaleFilterType::Pointer rescaler;
  Writer2Type::Pointer writer = Writer2Type::New();
  WriteImageType::Pointer scaledImage=WriteImageType::New();
  ImageIOType::Pointer dicomIO = ImageIOType::New();
  NamesGeneratorType::Pointer nameGenerator = NamesGeneratorType::New();

 WriteImageType::IndexType pixelIndex;
 WriteImageType::PixelType pixelValue;

};

#endif // IMAGESAMPLERDICOM_H
