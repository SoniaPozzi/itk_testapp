/****************************************************************/
/* MOOSE - Multiphysics Object Oriented Simulation Environment  */
/*                                                              */
/*          All contents are licensed under LGPL V2.1           */
/*             See LICENSE for full restrictions                */
/****************************************************************/

// MOOSE includes
#include "FileDicomChoose.h"


// itk includes
#include "itkCastImageFilter.h"
#include "itkRescaleIntensityImageFilter.h"
#include "itkGrayscaleFillholeImageFilter.h"
#include "itkThresholdImageFilter.h"
#include "itkConnectedThresholdImageFilter.h"
#include "itkCurvatureFlowImageFilter.h"
#include "itkImageMaskSpatialObject.h"
#include "itkRegionOfInterestImageFilter.h"
#include "itkTIFFImageIOFactory.h"
#include "itkImageSeriesWriter.h"
#include <itkFileTools.h>

template <>
InputParameters
validParams<FileDicomChoose>()
{
  InputParameters params = emptyInputParameters();
  params.addRequiredParam<std::string>("dicomDirectory", "Name of single image file to extract mesh parameters from.  If provided, a 2D mesh is created.");
  params.addRequiredParam<std::string>("file_base", "A filename in the DICOM series to use.");
  params.addRequiredParam<std::vector<Real>>("filtering_params","Number of iterations (suggested 5) and time step (suggested  0.125) for itk CurvatureFlowImageFilter ");
  params.addRequiredParam<std::vector<unsigned int>>("seed_index","Seed pixel index  (center of threshold) for itk ConnectedFilterType ");
  params.addRequiredParam<std::vector<Real>>("lower_upper_threshold_values","Lower and upper threshold values for itk ConnectedFilterType ");
 
  return params;
}

FileDicomChoose::FileDicomChoose(const InputParameters & params) : _status(0)
{

  itk::TIFFImageIOFactory::RegisterOneFactory();
  bool has_directory = params.isParamValid( "dicomDirectory" ),
       has_file_base = params.isParamValid( "file_base" );

  if (has_directory)
      _dicomDirectory = params.get<std::string>( "dicomDirectory" );
  
  if (has_file_base)
    _file_base = params.get<std::string>( "file_base" );
  
  // 1.) one of the two is not provided
  if ( !has_directory || !has_file_base )
  {
    _status = 1;
    return;
  }
  

  _filtering_params = params.get<std::vector<Real>>("filtering_params");
  _seed_index = params.get<std::vector<unsigned int>>("seed_index");
  _lower_upper_threshold_values = params.get<std::vector<Real>>("lower_upper_threshold_values");

if (_lower_upper_threshold_values[0]>_lower_upper_threshold_values[1])
   mooseError("Lower threshold value bigger than upper threhold value");


  try{  
    
    /// Reading all the series in the folder to check if there is the one related to _file_base
    reader -> SetImageIO( dicomIO );
    nameGenerator -> SetUseSeriesDetails( true );
    nameGenerator -> SetDirectory( _dicomDirectory );

    typedef std::vector< std::string >    SeriesIdContainer;
    const SeriesIdContainer & seriesUID = nameGenerator->GetSeriesUIDs();
    SeriesIdContainer::const_iterator seriesItr = seriesUID.begin();
    SeriesIdContainer::const_iterator seriesEnd = seriesUID.end();

    bool ret_value = false;

    while( seriesItr != seriesEnd && ret_value==false )
    {
      const ReaderType::FileNamesContainer & filenames = nameGenerator->GetFileNames( seriesItr -> c_str() );
      unsigned int numberOfFilenames =  filenames.size();
      
      for(unsigned int fni = 0; fni<numberOfFilenames; fni++)
      {
        std::string seriesName=filenames[fni];
        if (seriesName.find(_file_base) != std::string::npos) 
        {

          std::cout << "The directory: " << _dicomDirectory << std::endl;  
          std::cout << "contains the DICOM Series related to '" <<  _file_base << "'!" << std::endl;
          std::cout << "'The DICOM Series is: " << seriesItr -> c_str() << std::endl;
          std::cout <<  "There are "<< numberOfFilenames << " files in Serie: " << std::endl;
      
          reader -> SetFileNames(filenames);
          finalSeriesIdentifier = seriesItr -> c_str();
          ret_value = true;
    
          for(unsigned int fni = 0; fni<numberOfFilenames; fni++)
         {

          std::string seriesName = filenames[fni];
          std::cout << "          # " << fni << " Filename = " << seriesName << std::endl;
          
          }
  
          break;

        }
      }

      ++seriesItr;
    }

    /// In case the related series is not found, elencate the series available and give an error

    if(!ret_value)
    {  
      seriesItr = seriesUID.begin();
      mooseWarning("");
      std::cout << "The directory: " << _dicomDirectory << std::endl;  
      std::cout << "does not contain the DICOM Series related to '"<< _file_base<<"'! It has the following: "<<std::endl;
      
      while( seriesItr != seriesEnd )
     {     

        const ReaderType::FileNamesContainer & filenames = nameGenerator->GetFileNames(seriesItr->c_str());
        unsigned int numberOfFilenames =  filenames.size();
        std::cout <<std::endl;
        std::cout <<  "DICOM Series: " << seriesItr->c_str() <<";   Files in Serie   " << numberOfFilenames << std::endl;
                    
        for(unsigned int fni = 0; fni<numberOfFilenames; fni++)
          { 
            std::string seriesName=filenames[fni];
            std::cout << "          # " << fni << " Filename = " << seriesName << std::endl; 
          }   

          ++seriesItr;

          }

        mooseError("Set file_base to a filename in the DICOM series to use.");
      }

      std::cout << std::endl;  

      reader->Update();  
      std::cout  << "---------------------------------------------" << std::endl<<std::endl;
      std::cout << " Reading the DICOM Series of interest  " << std::endl << std::endl;
      std::cout << "   DICOM Serie Dimension:    " << reader -> GetOutput() -> GetLargestPossibleRegion().GetSize() << std::endl;
      std::cout << "   DICOM Serie Spacing:      " << reader -> GetOutput() -> GetSpacing() << std::endl;
      std::cout << "   DICOM Serie Origin:       " << reader -> GetOutput() -> GetOrigin() << std::endl << std::endl;
 
      /// Cast filter to convert from short data type to float data type

      typedef itk::CastImageFilter< ShortImageType, InternalImageType > CastFilterType;
      typedef itk::RescaleIntensityImageFilter<   InternalImageType, InternalImageType > RescaleFilterType;
      CastFilterType::Pointer castFilter = CastFilterType::New();
      RescaleFilterType::Pointer rescaler = RescaleFilterType::New();
      InternalImageType::Pointer scaledImage = InternalImageType::New();
      
      /// Castering from short to float data typpe requires rescale

      castFilter -> SetInput( reader -> GetOutput() );

      rescaler -> SetOutputMinimum(   0 );
      rescaler -> SetOutputMaximum( 255 ); //255 means white
      rescaler -> SetInput( castFilter -> GetOutput() );

      scaledImage = rescaler -> GetOutput();
      scaledImage -> GetLargestPossibleRegion();
      scaledImage -> Update(); 


      typedef itk::CastImageFilter< InternalImageType, OutputImageType > CastFilterOutType;
      typedef itk::ImageFileWriter< OutputImageType >  WriterType;
      CastFilterOutType::Pointer casterOut = CastFilterOutType::New();
      CastFilterOutType::Pointer caster = CastFilterOutType::New();
       
      // Images  are written in a sequence of tiff images to help locating Seed index and pixel value
       
      caster -> SetInput(scaledImage);
      caster -> Update();

     itk::FileTools::CreateDirectory("Output");
      const ReaderType::FileNamesContainer & filenamesfin = nameGenerator -> GetFileNames( finalSeriesIdentifier );
      unsigned int finalNumberOfFilenames =  filenamesfin.size();
      std::string format =  std::string( "Output/" ) + std::string( _file_base ) + std::string( "-rescaled-%d.tiff" );
      itk::NumericSeriesFileNames::Pointer fnames = itk::NumericSeriesFileNames::New();
      fnames->SetStartIndex( 0 );
      fnames->SetEndIndex(   finalNumberOfFilenames-1 );
      fnames->SetIncrementIndex( 1 );
      fnames->SetSeriesFormat( format.c_str() );

      typedef itk::Image< OutputPixelType, 2 > OutputImageType2d;
      typedef itk::ImageSeriesWriter< OutputImageType, OutputImageType2d >    WriterType2d;
      WriterType2d::Pointer outputWriter = WriterType2d::New();
      outputWriter->SetInput( caster -> GetOutput() );
      outputWriter->SetFileNames( fnames -> GetFileNames() );

      try
      {
        outputWriter -> Update();
      }
     catch (itk::ExceptionObject & e)
      {
        mooseError("Exception in file writer");
      }


      std::cout << " Open  Output/" << _file_base << "-rescaled-#.tiff to see the images in rescaled gray scale" << std::endl << std::endl;

      std::cout  << "---------------------------------------------" << std::endl<<std::endl;
      std::cout  << "Applying ITK filters  " << std::endl<<std::endl;
      std::cout  << " Applying connectedThreshold ITK filter:  " << std::endl;

      typedef itk::CurvatureFlowImageFilter< InternalImageType, InternalImageType > CurvatureFlowImageFilterType;
      typedef itk::ConnectedThresholdImageFilter< InternalImageType, InternalImageType > ConnectedFilterType;
      CurvatureFlowImageFilterType::Pointer smoothing = CurvatureFlowImageFilterType::New();
      ConnectedFilterType::Pointer connectedThreshold = ConnectedFilterType::New();

      smoothing -> SetNumberOfIterations( _filtering_params[0] ); 
      smoothing -> SetTimeStep( _filtering_params[1]  );
      smoothing -> SetInput( scaledImage );

      connectedThreshold -> SetInput( smoothing -> GetOutput() ); 
      connectedThreshold -> SetLower( _lower_upper_threshold_values[0] ); 
      connectedThreshold -> SetUpper( _lower_upper_threshold_values[1] );
      connectedThreshold -> SetReplaceValue( 255 ); //Considered threshold set to 255 (white)
      connectedThreshold -> SetInput( scaledImage ); 

 

   

      if ( _seed_index.size() < 1)
      {

       mooseError("Add at least a seed index for ConnectedFilterType.");

      }

      InternalImageType::IndexType  seedIndex;

      seedIndex[0]=_seed_index[0];
      seedIndex[1]=_seed_index[1];
      seedIndex[2]=_seed_index[2];

      connectedThreshold -> SetSeed( seedIndex );

      InternalImageType::PixelType pixel_value;
      pixel_value= scaledImage->GetPixel( seedIndex ); 

      std::cout   << "  selected Seed index:  " << seedIndex <<  " with Pixel Value: " << pixel_value << std::endl;


      if( _seed_index.size() % 3 != 0)
      {
       
       mooseError("A pixel index is missing or redundant in seed_index");

      }
      else{

      int numberOfSeeds = _seed_index.size()/3;

 
      InternalImageType::IndexType  addSeedIndex;

      for (unsigned i(1); i < numberOfSeeds; ++i)
     {
      
      addSeedIndex[0]=_seed_index[3*i];
      addSeedIndex[1]=_seed_index[3*i+1];
      addSeedIndex[2]=_seed_index[3*i+2];

      connectedThreshold -> AddSeed( addSeedIndex );

      pixel_value= scaledImage->GetPixel( addSeedIndex ); 

      std::cout   << "  selected Additional Seed index:  " << addSeedIndex <<  " with Pixel Value: " << pixel_value << std::endl;

}
      }

      std::cout<<std::endl;
      
      connectedThreshold -> Update();


      /// Cast filter to convert from float to unsigned char

      casterOut -> SetInput( connectedThreshold -> GetOutput() );
      casterOut -> UpdateLargestPossibleRegion(); 
      casterOut -> Update();


      std::cout  <<" Cropping of the Region of Interest (ITK filter):  " << std::endl;

      typedef itk::ImageMaskSpatialObject< 3 >        ImageMaskSpatialObjectType;
      typedef ImageMaskSpatialObjectType::ImageType   OutputImageType;
      typedef itk::RegionOfInterestImageFilter< OutputImageType, OutputImageType > FilterType;

      ImageMaskSpatialObjectType::Pointer       ImageMaskSpatialObject  = ImageMaskSpatialObjectType::New();
      ImageMaskSpatialObject->SetImage ( casterOut -> GetOutput() );
      OutputImageType::RegionType boundingBoxRegion = ImageMaskSpatialObject -> GetAxisAlignedBoundingBoxRegion();     
      std::cout   << "  Bounding Box Region Informations: " << boundingBoxRegion << std::endl;

      FilterType::Pointer filterRegion = FilterType::New();
      filterRegion -> SetRegionOfInterest( boundingBoxRegion );
      filterRegion -> SetInput( casterOut -> GetOutput() );
      filterRegion -> Update();

      filteredImage = filterRegion -> GetOutput();

      outputImageSize = filteredImage -> GetLargestPossibleRegion().GetSize(); 
      outputImageSpacing = filteredImage -> GetSpacing();



       std::cout  << " Applying FillholeFilterType ITK filter  " << std::endl<<std::endl;

      typedef itk::GrayscaleFillholeImageFilter<OutputImageType,OutputImageType>  FillholeFilterType;
      FillholeFilterType::Pointer  fillholes = FillholeFilterType::New();

      fillholes -> SetInput( filteredImage );
      fillholes -> Update();
    
      std::cout << " Information about the filtered and cropped DICOM Series  " << std::endl << std::endl;
      std::cout << "  DICOM Serie Dimension:    " << filteredImage -> GetLargestPossibleRegion().GetSize() << std::endl;
      std::cout << "  DICOM Serie Spacing:      " << filteredImage -> GetSpacing() << std::endl;
      std::cout << "  DICOM Serie Origin:       " << filteredImage -> GetOrigin() << std::endl << std::endl; 
  
      /// Visualize the results of the filtering and cropping writing a tiff image
      
      WriterType::Pointer writer = WriterType::New();
      writer -> SetInput( filteredImage ); //_file_base
      std::string name = std::string( "Output/" ) + std::string( _file_base ) + std::string( "-filtered-and-cropped.tiff" );
      writer -> SetFileName( name );
      
      try
      {
        writer -> Update();
      }
     catch (itk::ExceptionObject & e)
      {
        mooseError("Exception in file writer");
      }
    
     std::cout << " Open  Output/" << _file_base <<"-filtered-and-cropped.tiff to see the filtered and cropped DICOM Series" << std::endl << std::endl;
      std::cout  <<"---------------------------------------------" << std::endl;
      std::cout << "---------------------------------------------" << std::endl;
    }
  
    catch (itk::ExceptionObject &ex)
    {

      mooseError("exception in reading dicom series");
    
    }
  }

void
FileDicomChoose::errorCheck()
{
  switch (_status)
  {
    case 0:
   
      return;
   
    case 1:
   
      mooseError("not provided");
   
      break;
  }
}
