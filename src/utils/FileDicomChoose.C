/****************************************************************/
/* MOOSE - Multiphysics Object Oriented Simulation Environment  */
/*                                                              */
/*          All contents are licensed under LGPL V2.1           */
/*             See LICENSE for full restrictions                */
/****************************************************************/

// MOOSE includes
#include "FileDicomChoose.h"
#include "itkFileTools.h"

#include "itkImageFileReader.h"
#include "itkImageFileWriter.h"

// itk includes
#include "itkRescaleIntensityImageFilter.h"
#include "itkGrayscaleFillholeImageFilter.h"
#include "itkThresholdImageFilter.h"
#include "itkResampleImageFilter.h"
#include "itkNumericSeriesFileNames.h"

#include "itkConnectedThresholdImageFilter.h"
#include "itkCurvatureFlowImageFilter.h"
#include "itkImageMaskSpatialObject.h"
#include "itkRegionOfInterestImageFilter.h"
#include "itkTIFFImageIOFactory.h"
#include "itkImageSeriesWriter.h"
#include "itkNeighborhoodConnectedImageFilter.h"
#include "itkConfidenceConnectedImageFilter.h"
#include "itkVotingBinaryIterativeHoleFillingImageFilter.h"
#include "itkBinaryThresholdImageFilter.h"
#include "itkExtractImageFilter.h"
#include "itkAffineTransform.h"
#include "itkBSplineResampleImageFunction.h"


#include "itkMeanImageFilter.h"
#include "itkAddImageFilter.h"

#include "itkSiemensVisionImageIOFactory.h"
#include "itkMetaDataDictionary.h"
#include "itkMetaDataObject.h"
#include "DICOMAppHelper.h"

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

  if(_lower_upper_threshold_values[0]>_lower_upper_threshold_values[1])
     
    mooseError("Lower threshold value bigger than upper threhold value");

  try{  
    
    /// Reading all the series in the folder to check if there is the one related to _file_base
    reader -> SetImageIO( dicomIO );

    nameGenerator -> SetUseSeriesDetails( true );
    nameGenerator -> AddSeriesRestriction("0008|0021" );
    nameGenerator -> SetDirectory( _dicomDirectory );

    typedef std::vector< std::string >    SeriesIdContainer;
    const SeriesIdContainer & seriesUID = nameGenerator->GetSeriesUIDs();
    SeriesIdContainer::const_iterator seriesItr = seriesUID.begin();
    SeriesIdContainer::const_iterator seriesEnd = seriesUID.end();

    if ( 1 )
    { 
      //starting from the series name: working no problem
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
          
            typedef itk::GDCMImageIO ImageIOType;
            ImageIOType::Pointer dicomIO = ImageIOType::New(); reader->SetImageIO( dicomIO );
            reader->Update();

            finalSeriesIdentifier = seriesItr -> c_str();
            ret_value = true;
      
            for(unsigned int fni = 0; fni<numberOfFilenames; fni++)
            {

            std::string seriesName = filenames[fni];
            std::cout << "          # " << fni << " Filename = " << seriesName << std::endl;
            
            ////print tag info
            // typedef itk::MetaDataDictionary DictionaryType;
            // const DictionaryType & dictionary = dicomIO->GetMetaDataDictionary();
            // typedef itk::MetaDataObject< std::string > MetaDataStringType;
            // DictionaryType::ConstIterator itr = dictionary.Begin(); 
            // DictionaryType::ConstIterator end = dictionary.End();
            // while( itr != end )
            // {
            //   itk::MetaDataObjectBase::Pointer entry = itr->second; MetaDataStringType::Pointer entryvalue = dynamic_cast<MetaDataStringType *>( entry.GetPointer() );
            //   if( entryvalue )
            //   {
            //     std::string tagkey = itr->first;
            //     std::string tagvalue = entryvalue->GetMetaDataObjectValue(); std::cout << tagkey << " = " << tagvalue << std::endl;
            //   }
              
            //   ++itr; 
            // }
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
      }
    else
    {   
      //trying to consider different temporal Siemens series: to be fixed 
      std::string seriesName;

      typedef itk::ExtractImageFilter< ShortImageType, ShortImageType2D > FilterType;
      FilterType::Pointer extract = FilterType::New();
      extract -> InPlaceOn();
      extract -> SetDirectionCollapseToSubmatrix();

      unsigned int inputImageNumber = 0;
      typedef itk::NumericSeriesFileNames    NameGeneratorType;

      NameGeneratorType::Pointer nameGenerator = NameGeneratorType::New();

      nameGenerator -> SetSeriesFormat( _dicomDirectory + "MR000000_%d.dcm" );
      nameGenerator -> SetStartIndex( 4 );
      nameGenerator -> SetEndIndex(5);
      nameGenerator -> SetIncrementIndex( 1 );
      std::vector<std::string> names = nameGenerator->GetFileNames();

      // List the files

      std::vector<std::string>::iterator nit;
      for (nit = names.begin();
      nit != names.end();
      nit++)
      {
      std::cout << "File: " << (*nit).c_str() << std::endl;
      }

      JoinSeriesImageFilterType::Pointer joinFilter = JoinSeriesImageFilterType::New();
      joinFilter -> SetOrigin( 0.0 );
      joinFilter -> SetSpacing( 1.0 );
      joinFilter -> SetOrigin( reader -> GetOutput() -> GetOrigin()[2] );
      joinFilter -> SetSpacing( reader -> GetOutput() -> GetSpacing()[2] );
      reader -> SetImageIO( dicomIO );

      std::vector<ShortImageType::RegionType> region_types(names.size());
      std::vector<ShortImageType::SizeType>    size_types(names.size());
      std::vector<ShortImageType::IndexType> start_types(names.size());
      std::vector<ShortImageType::RegionType> desiredRegion_types(names.size());
      std::vector<ReaderType::Pointer> reader_types;
      std::vector<FilterType::Pointer> extract_types;
      std::vector<ShortImageType2D::Pointer> inputImageTile_types;

      const unsigned int sliceNumber = 0;

      for(unsigned int fni = 0; fni<names.size(); fni++)
      {  

        reader_types.push_back(ReaderType::New());
        reader_types[fni] -> SetImageIO( dicomIO );
        reader_types[fni] -> SetFileName( names[fni] );
        reader_types[fni] -> UpdateLargestPossibleRegion();
        reader_types[fni] -> Update();

        region_types[fni] = reader_types[fni]->GetOutput() -> GetLargestPossibleRegion();

        std::cout<<"Dicom size"<<reader_types[fni] -> GetOutput() -> GetLargestPossibleRegion() << std::endl;
        extract_types.push_back(FilterType::New());

        extract_types[fni] -> InPlaceOn();
        extract_types[fni] -> SetDirectionCollapseToSubmatrix();

        size_types[fni] = region_types[fni].GetSize();
        start_types[fni] = region_types[fni].GetIndex();
        size_types[fni][2] = 0;
        start_types[fni][2] = sliceNumber;

        desiredRegion_types[fni].SetSize( size_types[fni] ); 
        desiredRegion_types[fni].SetIndex( start_types[fni] );

        extract_types[fni] -> SetInput(   reader_types[fni]->GetOutput());
        extract_types[fni] -> SetExtractionRegion( desiredRegion_types[fni] );
        extract_types[fni] -> Update();

        inputImageTile_types.push_back(ShortImageType2D::New());

        inputImageTile_types[fni] = extract_types[fni] -> GetOutput();
        inputImageTile_types[fni] -> Allocate();

        std::cout << "Slice Dicom size" << inputImageTile_types[fni] -> GetLargestPossibleRegion() << std::endl;
        inputImageNumber++;

        joinFilter -> SetInput( inputImageNumber , inputImageTile_types[fni] );

      }

      joinFilter -> Update();
      ShortImageType::Pointer newDicom = ShortImageType::New();
      newDicom = joinFilter -> GetOutput();

    }

    std::cout  << "---------------------------------------------" << std::endl<<std::endl;
    std::cout << " Reading the DICOM Series of interest  " << std::endl << std::endl;
    std::cout << "   DICOM Serie Dimension:    " << reader->GetOutput() -> GetLargestPossibleRegion().GetSize() << std::endl;
    std::cout << "   DICOM Serie Spacing:      " << reader->GetOutput() -> GetSpacing() << std::endl;
    std::cout << "   DICOM Serie Origin:       " << reader->GetOutput() -> GetOrigin() << std::endl << std::endl;

    /// Cast filter to convert from short data type to float data type. Castering from short to float data type requires rescale

    typedef itk::CastImageFilter< ShortImageType, InternalImageType > CastFilterType;
    typedef itk::RescaleIntensityImageFilter<   InternalImageType, InternalImageType > RescaleFilterType;
    CastFilterType::Pointer castFilter = CastFilterType::New();
    RescaleFilterType::Pointer rescaler = RescaleFilterType::New();

    castFilter -> SetInput( reader->GetOutput() );

    rescaler -> SetOutputMinimum(   0 );
    rescaler -> SetOutputMaximum( 255 ); //255 means white
    rescaler -> SetInput( castFilter -> GetOutput() );
   
    scaledImage = rescaler -> GetOutput();
    scaledImage -> GetLargestPossibleRegion();
    scaledImage -> Update(); 

    // Images  are written in a sequence of tiff images to help locating Seed index and pixel value
    
    typedef itk::CastImageFilter< InternalImageType, OutputImageType > CastFilterOutType;
    typedef itk::ImageFileWriter< OutputImageType >  WriterType; ///qui
    CastFilterOutType::Pointer casterOut = CastFilterOutType::New();
    CastFilterOutType::Pointer caster = CastFilterOutType::New();
     
    caster -> SetInput(scaledImage);
    caster -> Update();

    itk::FileTools::CreateDirectory("Output");
   
    const ReaderType::FileNamesContainer & names = nameGenerator -> GetFileNames( finalSeriesIdentifier );
    unsigned int finalNumberOfFilenames =  names.size();
     
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
    outputWriter->SetFileNames( fnames->GetFileNames() );

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
    typedef itk::NeighborhoodConnectedImageFilter<InternalImageType, InternalImageType > NeighborhoodConnectedFilterType;

    typedef itk::ConfidenceConnectedImageFilter<InternalImageType, InternalImageType > ConfidenceonnectedFilterType;

    CurvatureFlowImageFilterType::Pointer smoothing = CurvatureFlowImageFilterType::New();
    ConnectedFilterType::Pointer connectedThreshold = ConnectedFilterType::New();
    NeighborhoodConnectedFilterType::Pointer neighborhoodConnectedThreshold = NeighborhoodConnectedFilterType::New();
    ConfidenceonnectedFilterType::Pointer confidenceConnectedThreshold = ConfidenceonnectedFilterType::New();

    smoothing -> SetNumberOfIterations( _filtering_params[0] ); 
    smoothing -> SetTimeStep( _filtering_params[1]  );
    smoothing -> SetInput( scaledImage );

    connectedThreshold -> SetInput( smoothing -> GetOutput() ); 
    connectedThreshold -> SetLower( _lower_upper_threshold_values[0] ); 
    connectedThreshold -> SetUpper( _lower_upper_threshold_values[1] );
    connectedThreshold -> SetReplaceValue( 255 ); //Considered threshold set to 255 (white)
    connectedThreshold -> SetInput( scaledImage ); 

    if ( _seed_index.size() < 1)

       mooseError("Add at least a seed index for ConnectedFilterType.");

    InternalImageType::IndexType  seedIndex;

    seedIndex[0] = _seed_index[0];
    seedIndex[1] = _seed_index[1];
    seedIndex[2] = _seed_index[2];

    connectedThreshold -> SetSeed( seedIndex );

    InternalImageType::PixelType pixel_value;
    pixel_value = scaledImage -> GetPixel( seedIndex ); 

    std::cout   << "  selected Seed index:  " << seedIndex <<  " with Pixel Value: " << pixel_value << std::endl;

    if( _seed_index.size() % 3 != 0)
      
      mooseError("A pixel index is missing or redundant in seed_index");

    else
    {

      int numberOfSeeds = _seed_index.size()/3;

      for (unsigned i(1); i < numberOfSeeds; ++i)
      {
      
        seedIndex[0] = _seed_index[3*i];
        seedIndex[1] = _seed_index[3*i+1];
        seedIndex[2] = _seed_index[3*i+2];

        connectedThreshold -> AddSeed( seedIndex );
        pixel_value = scaledImage->GetPixel( seedIndex ); 

        std::cout   << "  selected Additional Seed index:  " << seedIndex <<  " with Pixel Value: " << pixel_value << std::endl;

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

    std::cout << " Information about the cropped DICOM Series  " << std::endl << std::endl;
    std::cout << "  DICOM Serie Dimension:    " << filteredImage -> GetLargestPossibleRegion().GetSize() << std::endl;
    std::cout << "  DICOM Serie Spacing:      " << filteredImage -> GetSpacing() << std::endl;
    std::cout << "  DICOM Serie Origin:       " << filteredImage -> GetOrigin() << std::endl << std::endl; 

    std::cout  <<" Resampling the cropped series by spline to fill the holes in series direction (ITK filter)... " << std::endl;
    
    //1. creating a image interpolating in series direction (to fill holes)

    typedef itk::ResampleImageFilter<  OutputImageType, OutputImageType >  ResampleFilterType;
    typedef itk::IdentityTransform<double, 3>  TransformType;
    typedef itk::BinaryThresholdImageFilter<  OutputImageType, OutputImageType > BinaryType;

    
    ResampleFilterType::Pointer resampler = ResampleFilterType::New();
    TransformType::Pointer   transform  = TransformType::New();
    BinaryType::Pointer binaryThreshold = BinaryType::New();

    transform -> SetIdentity();
    interpolator -> SetSplineOrder(5);

    resampler -> SetTransform( transform );
    resampler -> SetInterpolator( interpolator );
    resampler -> SetOutputOrigin( filteredImage->GetOrigin() );
    resampler -> SetOutputSpacing( filteredImage->GetSpacing() );
    resampler -> SetOutputDirection( filteredImage->GetDirection() );
    resampler -> SetSize( filteredImage->GetLargestPossibleRegion().GetSize());
    resampler -> SetInput( filteredImage );
    resampler -> UpdateLargestPossibleRegion();
    resampler -> Update();

    // the obtained result is converted in binary by thresholding  
    binaryThreshold -> SetInput(  resampler->GetOutput()  );
    binaryThreshold -> SetOutsideValue( 0 );
    binaryThreshold -> SetInsideValue(  255  );
    binaryThreshold -> SetLowerThreshold( 1 );
    binaryThreshold -> SetUpperThreshold( 255 );
    binaryThreshold -> Update();

    OutputImageType::SizeType radius;
    radius[0] = 2; // radius along x
    radius[1] = 2; // radius along y
    radius[2] = 0; // radius along z 
  
    std::cout << " Filling the holes in the other directions by a voting strategy (ITK filter)... " << std::endl;

    //2. creating a image interpolating in other directions (to fill holes)

    typedef itk::VotingBinaryIterativeHoleFillingImageFilter< OutputImageType > VotingType;
    VotingType::Pointer voting = VotingType::New();

    voting -> SetInput( binaryThreshold->GetOutput() );
    voting -> SetRadius( radius );
    voting -> SetMajorityThreshold( 2 );
    voting -> SetBackgroundValue( 0 );
    voting -> SetForegroundValue( 255);
    voting -> SetMaximumNumberOfIterations( 120 );
    voting -> Update();
        

    std::cout   << "  The filter used " << voting -> GetCurrentNumberOfIterations() << " iterations " << std::endl;
    std::cout   << "  and changed a total of " << voting -> GetNumberOfPixelsChanged() << " pixels" << std::endl;

    // 3. Sum of the two images (1. and 2.): the initial interpolated in series direction one, resampler -> GetOutput() and the one with filled holes, i.e. voting -> GetOutput. 
    // The image voting -> GetOutput is rescaled (to have value 100 as holes filler). The two images are summed and the sum is rescaled.

    typedef itk::AddImageFilter<OutputImageType,  OutputImageType,  OutputImageType > AddFilterType;
    typedef itk::RescaleIntensityImageFilter< OutputImageType, OutputImageType >   RescalerOutType;
    AddFilterType::Pointer addFilter = AddFilterType::New();
    RescalerOutType::Pointer rescalerOut = RescalerOutType::New();
    RescalerOutType::Pointer rescalerOutSum = RescalerOutType::New();

    rescalerOut -> SetOutputMinimum( 0 );
    rescalerOut -> SetOutputMaximum( 100 );
    rescalerOut -> SetInput( voting -> GetOutput() );

    addFilter -> SetInput1( rescalerOut -> GetOutput() );
    addFilter -> SetInput2( resampler -> GetOutput() );
    addFilter -> Update();

    rescalerOutSum -> SetOutputMinimum(  0 );
    rescalerOutSum -> SetOutputMaximum( 255 );
    rescalerOutSum -> SetInput( addFilter -> GetOutput() );

    std::cout  << " Computing the mean over the pixels in the neighborhood to smooth the image (ITK filter)...  " << std::endl << std::endl;

    typedef itk::MeanImageFilter<  OutputImageType, OutputImageType >  MeanType;
    MeanType::Pointer mean = MeanType::New();

    OutputImageType::SizeType indexRadius;
    indexRadius[0] = 2; // radius along x
    indexRadius[1] = 2; // radius along y
    indexRadius[2] = 0; // radius along y
    mean -> SetInput( rescalerOutSum->GetOutput()  );
    mean -> SetRadius( indexRadius );

    //Finally save the obtained image.

    filteredImage = mean -> GetOutput();

    WriterType::Pointer writer = WriterType::New();
    writer -> SetInput( filteredImage  ); //_file_base
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

    std::cout << "---------------------------------------------" << std::endl<<std::endl; 
    std::cout << " Open  Output/" << _file_base <<"-filtered-and-cropped.tiff to see the filtered and cropped DICOM Series" << std::endl<<std::endl;
    std::cout  <<"---------------------------------------------" << std::endl<<std::endl;
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
