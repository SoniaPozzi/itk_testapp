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
#include "itkDicomImageIOFactory.h"
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

/////
#include "itkImage.h"
#include "itkGDCMImageIO.h"
#include "itkGDCMSeriesFileNames.h"
#include "itkImageSeriesReader.h"
#include "itkImageFileWriter.h"
#include "itkNiftiImageIO.h"
#include "itkPermuteAxesImageFilter.h"



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

    if ( 0 )
    { 

    nameGenerator -> AddSeriesRestriction("0008|0031");
    nameGenerator -> SetDirectory( _dicomDirectory );

    typedef std::vector< std::string >    SeriesIdContainer;
    const SeriesIdContainer & seriesUID = nameGenerator->GetSeriesUIDs();
    SeriesIdContainer::const_iterator seriesItr = seriesUID.begin();
    SeriesIdContainer::const_iterator seriesEnd = seriesUID.end();

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

          nameGenerator -> SetDirectory( _dicomDirectory );
          ReaderType::Pointer reader = ReaderType::New();

          // 3D and 4D Image Types  
          ShortImageType4D::Pointer image4D = ShortImageType4D::New();
          ShortImageType::Pointer image3D = ShortImageType::New();

          // set read I/O
          reader->SetImageIO(dicomIO);

          // restriction on DICOM  (acquisition time)
          const std::string entry_id = "0018|1090 ";

          nameGenerator->AddSeriesRestriction(entry_id);

          // Series ID

          typedef std::vector< std::string >    SeriesIdContainer;
          const SeriesIdContainer & seriesUID = nameGenerator->GetSeriesUIDs();
          SeriesIdContainer::const_iterator seriesItr = seriesUID.begin();
          SeriesIdContainer::const_iterator seriesEnd = seriesUID.end();

          // container with file names 
          typedef std::vector< std::string >   FileNamesContainer;

          // figure out how many time points there are based on number of series 
          unsigned int timepoints = 0;
          while (seriesItr != seriesEnd){
          timepoints++;
          seriesItr++;
          }
          std::cout << "Number of time points  : " << timepoints << std::endl;

          // read first image
          seriesItr = seriesUID.begin();
          std::cout << "Reading first image : " << std::endl;
          std::cout << seriesItr->c_str() << std::endl;

          std::string seriesIdentifier;
          seriesIdentifier = seriesItr->c_str();

          // generate file names 
          FileNamesContainer fileNames;
          fileNames = nameGenerator->GetFileNames(seriesIdentifier);
          FileNamesContainer::iterator fileiterator = fileNames.begin();


          std::cout<<"sono qui"<<seriesIdentifier<<std::endl;
          // get Image info
          //itk::ImageIOBase::Pointer imageIO = getImageIO(fileiterator->c_str());

          // read file 
          reader->SetFileNames(fileNames);

          try
          {
          reader->Update();
          }
          catch (itk::ExceptionObject &ex)
          {
          std::cout << ex << std::endl;
          mooseError("Exception in file reader");
          }


          // spacing 
          const ShortImageType::SpacingType& spacing3D = reader->GetOutput()->GetSpacing();
          std::cout << "Spacing 3D = " << spacing3D[0] << ", " << spacing3D[1] << ", " << spacing3D[2] << std::endl;
          // origin info 
          const ShortImageType::PointType origin3D = reader->GetOutput()->GetOrigin();
          std::cout << "Origin  3D = " << origin3D[0] << ", " << origin3D[1] << ", " << origin3D[2] << std::endl;
          // get 3D size 
          const ShortImageType::SizeType size3D = reader->GetOutput()->GetBufferedRegion().GetSize();

          // decleare 4D parameters 
          ShortImageType4D::SpacingType spacing4D;
          ShortImageType4D::PointType origin4D;
          ShortImageType4D::SizeType size4D;

          // copy parameters from 3D
          for (int i = 0; i < 3; ++i){
          spacing4D[i] = spacing3D[i];
          origin4D[i] = origin3D[i];
          size4D[i] = size3D[i];
          }

          // add 4-th dimension 
          spacing4D[3] = 1;
          origin4D[3] = 0;
          size4D[3] = timepoints;

          // start 
          ShortImageType4D::IndexType start4D;
          start4D.Fill(0);

          // set  region 
          ShortImageType4D::RegionType region4D(start4D, size4D);

          std::cout << "Spacing 4D = " << spacing4D[0] << ", " << spacing4D[1] << ", " << spacing4D[2] << ", " << spacing4D[3] << std::endl;
          std::cout << "Size 4D = " << size4D[0] << ", " << size4D[1] << ", " << size4D[2] << ", " << size4D[3] << std::endl;

          // allocate 4D image  
          image4D->SetRegions(region4D);
          image4D->SetSpacing(spacing4D);
          image4D->SetOrigin(origin4D);
          image4D->Allocate();

          // point again to first series ID   
          seriesItr = seriesUID.begin();

          // iterators to loop into images 
          typedef itk::ImageRegionConstIterator< ShortImageType >  Iterator3D;
          typedef itk::ImageRegionIterator< ShortImageType4D >  Iterator4D;

          // define pointer to 4d image
          Iterator4D it4(image4D, image4D->GetBufferedRegion());
          it4.GoToBegin();

          // Loop to read the Dicom volumes one by one   
          unsigned short int idx = 0;
          while (seriesItr != seriesEnd){

          seriesIdentifier = seriesItr->c_str();
          std::cout << "Reading series " << std::endl;
          std::cout << seriesItr->c_str() << std::endl;

          // generate file names 
          fileNames = nameGenerator->GetFileNames(seriesIdentifier);

          for(unsigned int fni = 0; fni<fileNames.size(); fni++)
          {  std::cout<<"FN "<<fileNames[fni]<<std::endl;

            }

          // read volume files 
          reader->SetFileNames(fileNames);

          // set image3D 
          image3D = reader->GetOutput();
          image3D->SetRegions(reader->GetOutput()->GetRequestedRegion());
          image3D->CopyInformation(reader->GetOutput());
          image3D->Allocate();

          std::cout << "reading image volume " << idx << std::endl << std::endl;
          try
          {
          reader->Update();
          }
          catch (itk::ExceptionObject &ex)
          {
          std::cout << ex << std::endl;
               mooseError("Exception in file reader");
          }

          // point to the current image     
          Iterator3D it3(image3D, image3D->GetBufferedRegion());
          it3.GoToBegin();

          while (!it3.IsAtEnd())
          {
          it4.Set(it3.Get());
          ++it3;
          ++it4;
          }

          // increment iterator
          seriesItr++; 
          idx++;

          std::string format2 =  std::string( "Output/Rescaled-%d.dcm" );
          itk::NumericSeriesFileNames::Pointer fnames2 = itk::NumericSeriesFileNames::New();
          fnames2->SetStartIndex( 0 );
          fnames2->SetEndIndex(  0);
          fnames2->SetIncrementIndex( 1 );
          fnames2->SetSeriesFormat( format2.c_str() );

          typedef itk::ImageSeriesWriter<  ShortImageType, ShortImageType>  SeriesWriterType;
          SeriesWriterType::Pointer seriesWriter = SeriesWriterType::New();

          seriesWriter->SetInput( image3D );
          seriesWriter->SetImageIO( dicomIO );

          // Software Guide : BeginCodeSnippet
          nameGenerator->SetOutputDirectory( _dicomDirectory );

          seriesWriter->SetFileNames(fnames2->GetFileNames());

          seriesWriter->SetMetaDataDictionaryArray(
                           reader->GetMetaDataDictionaryArray() );  
          try{ seriesWriter->Update();}
          catch (itk::ExceptionObject &ex)
          {
          std::cout << ex << std::endl;
          mooseError("Exception in file writer");
          }



          }


          typedef itk::PermuteAxesImageFilter< ShortImageType4D > SwapFilterType;
          SwapFilterType::Pointer swap = SwapFilterType::New();
          swap->SetInput( image4D );

          SwapFilterType::PermuteOrderArrayType order;
          order[0] = 0;
          order[1] = 1;
          order[2] = 3;
          order[3] = 2;
          swap->SetOrder( order );


          const ShortImageType4D::SizeType size4Dr = swap->GetOutput()->GetBufferedRegion().GetSize();
          std::cout<<"new size"<<size4Dr <<std::endl;

          std::string format3 =  std::string( "Output/-rescaled-%d.nii.gz" );
          itk::NumericSeriesFileNames::Pointer fnames3 = itk::NumericSeriesFileNames::New();
          fnames3->SetStartIndex( 0 );
          fnames3->SetEndIndex(   25-1 );
          fnames3->SetIncrementIndex( 1 );
          fnames3->SetSeriesFormat( format3.c_str() );

          std::cout << " creating nifti file ...   " << std::endl << std::endl;
          typedef itk::Image<short, 4>  DWI;
          itk::NiftiImageIO::Pointer nifti_io = itk::NiftiImageIO::New();

          itk::ImageFileWriter<DWI>::Pointer dwi_writer = itk::ImageFileWriter<DWI>::New();
          dwi_writer->SetFileName("test.nii.gz");
          dwi_writer->SetInput(swap->GetOutput());

          dwi_writer->SetImageIO(nifti_io);
          dwi_writer->Update();


return ;


      //we try to read the cardiac series:


      // //trying to consider different temporal Siemens series: to be fixed 
      // std::string seriesName;

      // typedef itk::ExtractImageFilter< ShortImageType, ShortImageType2D > FilterType;
      // FilterType::Pointer extract = FilterType::New();
      // extract -> InPlaceOn();
      // extract -> SetDirectionCollapseToSubmatrix();

      // unsigned int inputImageNumber = 0;
      // typedef itk::NumericSeriesFileNames    NameGeneratorType;

      // NameGeneratorType::Pointer nameGenerator = NameGeneratorType::New();

      // nameGenerator -> SetSeriesFormat( _dicomDirectory + "MR000000_%d.dcm" );
      // nameGenerator -> SetStartIndex( 4 );
      // nameGenerator -> SetEndIndex(5);
      // nameGenerator -> SetIncrementIndex( 1 );
      // std::vector<std::string> names = nameGenerator->GetFileNames();

      // // List the files

      // std::vector<std::string>::iterator nit;
      // for (nit = names.begin();
      // nit != names.end();
      // nit++)
      // {
      // std::cout << "File: " << (*nit).c_str() << std::endl;
      // }

      // JoinSeriesImageFilterType::Pointer joinFilter = JoinSeriesImageFilterType::New();
      // joinFilter -> SetOrigin( 0.0 );
      // joinFilter -> SetSpacing( 1.0 );
      // joinFilter -> SetOrigin( reader -> GetOutput() -> GetOrigin()[2] );
      // joinFilter -> SetSpacing( reader -> GetOutput() -> GetSpacing()[2] );
      // reader -> SetImageIO( dicomIO );

      // std::vector<ShortImageType::RegionType> region_types(names.size());
      // std::vector<ShortImageType::SizeType>    size_types(names.size());
      // std::vector<ShortImageType::IndexType> start_types(names.size());
      // std::vector<ShortImageType::RegionType> desiredRegion_types(names.size());
      // std::vector<ReaderType::Pointer> reader_types;
      // std::vector<FilterType::Pointer> extract_types;
      // std::vector<ShortImageType2D::Pointer> inputImageTile_types;

      // const unsigned int sliceNumber = 0;

      // for(unsigned int fni = 0; fni<names.size(); fni++)
      // {  

      //   reader_types.push_back(ReaderType::New());
      //   reader_types[fni] -> SetImageIO( dicomIO );
      //   reader_types[fni] -> SetFileName( names[fni] );
      //   reader_types[fni] -> UpdateLargestPossibleRegion();
      //   reader_types[fni] -> Update();

      //   region_types[fni] = reader_types[fni]->GetOutput() -> GetLargestPossibleRegion();

      //   std::cout<<"Dicom size"<<reader_types[fni] -> GetOutput() -> GetLargestPossibleRegion() << std::endl;
      //   extract_types.push_back(FilterType::New());

      //   extract_types[fni] -> InPlaceOn();
      //   extract_types[fni] -> SetDirectionCollapseToSubmatrix();

      //   size_types[fni] = region_types[fni].GetSize();
      //   start_types[fni] = region_types[fni].GetIndex();
      //   size_types[fni][2] = 0;
      //   start_types[fni][2] = sliceNumber;

      //   desiredRegion_types[fni].SetSize( size_types[fni] ); 
      //   desiredRegion_types[fni].SetIndex( start_types[fni] );

      //   extract_types[fni] -> SetInput(   reader_types[fni]->GetOutput());
      //   extract_types[fni] -> SetExtractionRegion( desiredRegion_types[fni] );
      //   extract_types[fni] -> Update();

      //   inputImageTile_types.push_back(ShortImageType2D::New());

      //   inputImageTile_types[fni] = extract_types[fni] -> GetOutput();
      //   inputImageTile_types[fni] -> Allocate();

      //   std::cout << "Slice Dicom size" << inputImageTile_types[fni] -> GetLargestPossibleRegion() << std::endl;
      //   inputImageNumber++;

      //   joinFilter -> SetInput( inputImageNumber , inputImageTile_types[fni] );

      // }

      // joinFilter -> Update();
      // ShortImageType::Pointer newDicom = ShortImageType::New();
      // newDicom = joinFilter -> GetOutput();

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


inline itk::ImageIOBase::IOPixelType FileDicomChoose::pixel_type(itk::ImageIOBase::Pointer imageIO){
  return imageIO->GetPixelType();
}


inline itk::ImageIOBase::Pointer FileDicomChoose::getImageIO(std::string input){
  itk::ImageIOBase::Pointer imageIO = itk::ImageIOFactory::CreateImageIO(input.c_str(), itk::ImageIOFactory::ReadMode);


std::cout<<"qui"<<std::endl;
  imageIO->SetFileName(input.c_str());
  imageIO->ReadImageInformation();
std::cout<<"quo"<<std::endl;
  return imageIO;
}


