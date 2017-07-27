/****************************************************************/
/* MOOSE - Multiphysics Object Oriented Simulation Environment  */
/*                                                              */
/*          All contents are licensed under LGPL V2.1           */
/*             See LICENSE for full restrictions                */
/****************************************************************/

#ifndef FILEDICOMCHOOSE_H
#define FILEDICOMCHOOSE_H

// MOOSE includes
#include "Moose.h"
#include "InputParameters.h"
#include <iostream>

// itk includes

#include "itkGDCMImageIO.h"
#include "itkGDCMSeriesFileNames.h"
#include "itkImageSeriesReader.h"
#include "itkImageIOBase.h"


#include "itkDICOMSeriesFileNames.h"
#include "itkCastImageFilter.h"
#include "itkJoinSeriesImageFilter.h"
#include "itkBSplineResampleImageFunction.h"

// Forward declarations
class FileDicomChoose;
class InputParameters;

template <typename T>
InputParameters validParams();

/**
 * To be called in the validParams functions of classes that need to
 * operate on ranges of files.  Adds several non-required parameters
 * that are parsed in the parseFileRange function.
 */
template <>
InputParameters validParams<FileDicomChoose>();

/**
 * Augments an InputParameters object with file range information.
 * Creates and adds a vector<string> with the list of filenames to the
 * params object for use by the calling object.  The params object
 * passed in must contain suitable information for building the list
 * of filenames in the range.  Returns a non-zero error code if there
 * is an error while parsing.
 */
class FileDicomChoose
{
  public:

  FileDicomChoose(const InputParameters & params);
  virtual ~FileDicomChoose() = default;
  std::string & dicomDirectory() { return _dicomDirectory; }

  protected:
  
  int status(){ return _status; }
  void errorCheck();

  int _status;

  std::vector<Real>  _filtering_params;
  std::vector<unsigned int>  _seed_index;
  std::vector<Real>  _lower_upper_threshold_values;
  
  std::string _dicomDirectory;
  std::string _file_base;

  bool _has_chooseSeries;


  const unsigned int Dimension = 3;

  typedef  signed short    ShortPixelType;
  typedef  float           FloatPixelType;
  typedef unsigned char    OutputPixelType;

  typedef itk::Image< ShortPixelType ,  3 >      ShortImageType;
  typedef itk::Image< FloatPixelType ,  3 >      InternalImageType;
  typedef itk::Image< OutputPixelType , 3 >      OutputImageType;
  typedef itk::Image< OutputPixelType , 2 >      OutputImageType2D;   
  typedef itk::Image< ShortPixelType ,  2 >      ShortImageType2D;

  typedef itk::ImageSeriesReader< ShortImageType > ReaderType;
  typedef itk::ImageSeriesReader< ShortImageType2D > ReaderType2;
  typedef itk::CastImageFilter<ShortImageType,ShortImageType2D>  ImageTypecast;
  typedef itk::JoinSeriesImageFilter<ShortImageType2D, ShortImageType> JoinSeriesImageFilterType;
  typedef itk::GDCMImageIO ImageIOType;
  typedef itk::GDCMSeriesFileNames NamesGeneratorType;
  typedef std::vector< std::string > FileNamesContainer;
  typedef itk::BSplineResampleImageFunction<OutputImageType,  double>    InterpolatorType;


  ImageIOType::Pointer dicomIO = ImageIOType::New();
  ReaderType::Pointer reader = ReaderType::New();
  NamesGeneratorType::Pointer nameGenerator = NamesGeneratorType::New();
  FileNamesContainer fileNames;
  std::string finalSeriesIdentifier;


  InternalImageType::Pointer scaledImage = InternalImageType::New();  
  OutputImageType::Pointer filteredImage =OutputImageType::New();
  OutputImageType::SizeType outputImageSize;
  OutputImageType::SpacingType outputImageSpacing;
  InterpolatorType::Pointer   interpolator  = InterpolatorType::New();

};

#endif // FILEDICOMCHOOSE_H
