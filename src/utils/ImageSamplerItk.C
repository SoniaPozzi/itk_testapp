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

// MOOSE includes
#include "ImageSamplerItk.h"
#include "MooseApp.h"
#include "ImageMesh.h"
#include <iostream>

// ITK headers
#include "itkGDCMImageIO.h"
#include "itkGDCMSeriesFileNames.h"
#include "itkTIFFImageIO.h"
#include "itkTIFFImageIOFactory.h"
#include "itkRescaleIntensityImageFilter.h"

#include "itkPointSet.h"
#include <itkPoint.h>



template <>
InputParameters
validParams<ImageSamplerItk>()
{
  // Define the general parameters
  InputParameters params = emptyInputParameters();
  params += validParams<FileRangeBuilder>();

  params.addParam<Point>("origin", "Origin of the image (defaults to mesh origin)");
  params.addParam<Point>("dimensions",
                         "x,y,z dimensions of the image (defaults to mesh dimensions)");
  params.addParam<unsigned int>(
      "component",
      "The image RGB-component to return, leaving this blank will result in a greyscale value "
      "for the image to be created. The component number is zero based, i.e. 0 returns the first "
      "(RED) component of the image.");

  // Shift and Scale (application of these occurs prior to threshold)
  params.addParam<double>("shift", 0, "Value to add to all pixels; occurs prior to scaling");
  params.addParam<double>(
      "scale", 1, "Multiplier to apply to all pixel values; occurs after shifting");
  params.addParamNamesToGroup("shift scale", "Rescale");

  // Threshold parameters
  params.addParam<double>("threshold", "The threshold value");
  params.addParam<double>(
      "upper_value", 1, "The value to set for data greater than the threshold value");
  params.addParam<double>(
      "lower_value", 0, "The value to set for data less than the threshold value");
  params.addParamNamesToGroup("threshold upper_value lower_value", "Threshold");

  // Flip image
  params.addParam<bool>("flip_x", false, "Flip the image along the x-axis");
  params.addParam<bool>("flip_y", false, "Flip the image along the y-axis");
  params.addParam<bool>("flip_z", false, "Flip the image along the z-axis");
  params.addParamNamesToGroup("flip_x flip_y flip_z", "Flip");

  return params;
}

ImageSamplerItk::ImageSamplerItk(const InputParameters & parameters)
  : FileRangeBuilder(parameters),
#ifdef LIBMESH_HAVE_VTK
    _data(NULL),
    _algorithm(NULL),
#endif
    _is_pars(parameters),
    _is_console((parameters.getCheckedPointerParam<MooseApp *>("_moose_app"))->getOutputWarehouse())

{
#ifndef LIBMESH_HAVE_VTK
  // This should be impossible to reach, the registration of ImageSampler is also guarded with
  // LIBMESH_HAVE_VTK
  mooseError("libMesh must be configured with VTK enabled to utilize ImageSampler");
#endif
}

void
ImageSamplerItk::setupImageSampler(MooseMesh & mesh)
{

  std::cout<<"USING ITK...................."<<std::endl;

  ItkImageSampler(mesh);
  // Don't warn that mesh or _is_pars are unused when VTK is not enabled.
  std::cout<<"Back to VTK..................."<<std::endl;
  libmesh_ignore(mesh);
  libmesh_ignore(_is_pars);

#ifdef LIBMESH_HAVE_VTK
  // Get access to the Mesh object
  MeshTools::BoundingBox bbox = MeshTools::bounding_box(mesh.getMesh());

  // Set the dimensions from the Mesh if not set by the User
  if (_is_pars.isParamValid("dimensions"))
    _physical_dims = _is_pars.get<Point>("dimensions");

  else
  {
    _physical_dims(0) = bbox.max()(0) - bbox.min()(0);
#if LIBMESH_DIM > 1
    _physical_dims(1) = bbox.max()(1) - bbox.min()(1);
#endif
#if LIBMESH_DIM > 2
    _physical_dims(2) = bbox.max()(2) - bbox.min()(2);
#endif
  }

  // Set the origin from the Mesh if not set in the input file
  if (_is_pars.isParamValid("origin"))
    _origin = _is_pars.get<Point>("origin");
  else
  {
    _origin(0) = bbox.min()(0);
#if LIBMESH_DIM > 1
    _origin(1) = bbox.min()(1);
#endif
#if LIBMESH_DIM > 2
    _origin(2) = bbox.min()(2);
#endif
  }

  // An array of filenames, to be filled in
  std::vector<std::string> filenames;

  // The file suffix, to be determined
  std::string file_suffix;

  // Try to parse our own file range parameters.  If that fails, then
  // see if the associated Mesh is an ImageMesh and use its.  If that
  // also fails, then we have to throw an error...
  //
  // The parseFileRange method sets parameters, thus a writable reference to the InputParameters
  // object must be obtained from the warehouse. Generally, this should be avoided, but
  // this is a special case.
  if (_status != 0)
  {
    // We don't have parameters, so see if we can get them from ImageMesh
    ImageMesh * image_mesh = dynamic_cast<ImageMesh *>(&mesh);
    if (!image_mesh)
      mooseError("No file range parameters were provided and the Mesh is not an ImageMesh.");

    // Get the ImageMesh's parameters.  This should work, otherwise
    // errors would already have been thrown...
    filenames = image_mesh->filenames();
    file_suffix = image_mesh->fileSuffix();
  }
  else
  {
    // Use our own parameters (using 'this' b/c of conflicts with filenames the local variable)
    filenames = this->filenames();
    file_suffix = fileSuffix();
  }

  // Storage for the file names
  _files = vtkSmartPointer<vtkStringArray>::New();

  for (const auto & filename : filenames){
    _files->InsertNextValue(filename);
    std::cout<<filename<<std::endl;}

  // Error if no files where located
  if (_files->GetNumberOfValues() == 0)
    mooseError("No image file(s) located");

  // Read the image stack.  Hurray for VTK not using polymorphism in a
  // smart way... we actually have to explicitly create the type of
  // reader based on the file extension, using an if-statement...
  if (file_suffix == "dicom" || file_suffix == "dcm")
    _image = vtkSmartPointer<vtkDICOMImageReader>::New();

  else
    mooseError("Un-supported file type '", file_suffix, "'");

  // Now that _image is set up, actually read the images
  // Indicate that data read has started
  _is_console << "Reading image(s)..." << std::endl;

  _image->SetDirectoryName("/Users/sonia/Projects/internshipINL/itk_testapp/tests/image_function/new_stack/");
  
/////// e' come se fossi arrivata qui

  _image->Update();
  _data = _image->GetOutput();
  _algorithm = _image->GetOutputPort();


  // Set the image dimensions and voxel size member variable
  int * dims = _data->GetDimensions();

  for (unsigned int i = 0; i < 3; ++i)
  {
    _dims.push_back(dims[i]);
   int prova= (_dims[i]);
   std::cout<<prova<<std::endl;
    _voxel.push_back(_physical_dims(i) / _dims[i]);
   std::cout<<"voxel dopo"<<std::endl;
  std::cout<<_voxel[i]<<std::endl;
  }

  // Set the dimensions of the image and bounding box
  _data->SetSpacing(_voxel[0], _voxel[1], _voxel[2]);
  _data->SetOrigin(_origin(0), _origin(1), _origin(2));
  _bounding_box.min() = _origin;
  _bounding_box.max() = _origin + _physical_dims;

  // Indicate data read is completed
  _is_console << "          ...image read finished" << std::endl;

  // Set the component parameter
  // If the parameter is not set then vtkMagnitude() will applied
  if (_is_pars.isParamValid("component"))
  {
    unsigned int n = _data->GetNumberOfScalarComponents();


std::cout<<"Number Of Components Per Pixel"<<n <<std::endl;
    _component = _is_pars.get<unsigned int>("component");
    if (_component >= n)
      mooseError("'component' parameter must be empty or have a value of 0 to ", n - 1);
  }
  else
    _component = 0;

  // Apply filters, the toggling on and off of each filter is handled internally
  vtkMagnitude();
#endif
}

Real
ImageSamplerItk::sample(const Point & p)
{


#ifdef LIBMESH_HAVE_VTK

  // Do nothing if the point is outside of the image domain
  if (!_bounding_box.contains_point(p))
    return 0.0;

  // Determine pixel coordinates
  std::vector<int> x(3, 0);
  for (int i = 0; i < LIBMESH_DIM; ++i)
  {
    // Compute position, only if voxel size is greater than zero
    if (_voxel[i] == 0)
      x[i] = 0;

    else
    {
      x[i] = std::floor((p(i) - _origin(i)) / _voxel[i]);

      // If the point falls on the mesh extents the index needs to be decreased by one
      if (x[i] == _dims[i])
        x[i]--;
    }
  }

if (p(0)==1.0 &&p(1)==1.0 &&p(2)==1.0 ){
  std::cout<<"VTK SAMPLE"<<std::endl;
  std::cout<<p<<std::endl;
  std::cout<< _data->GetScalarComponentAsDouble(x[0], x[1], x[2], _component)<<std::endl;

}

  // Return the image data at the given point
  return _data->GetScalarComponentAsDouble(x[0], x[1], x[2], _component);




#else
  libmesh_ignore(p); // avoid un-used parameter warnings
  return 0.0;
#endif
}

void
ImageSamplerItk::vtkMagnitude()
{
#ifdef LIBMESH_HAVE_VTK
  // Do nothing if 'component' is set
  if (_is_pars.isParamValid("component"))
    return;

  // Apply the greyscale filtering
  _magnitude_filter = vtkSmartPointer<vtkImageMagnitude>::New();
  _magnitude_filter->SetInputConnection(_algorithm);
  _magnitude_filter->Update();

  // Update the pointers
  _data = _magnitude_filter->GetOutput();
  _algorithm = _magnitude_filter->GetOutputPort();
#endif
}


void
ImageSamplerItk::ItkImageSampler(MooseMesh & mesh)
{


  itk::TIFFImageIOFactory::RegisterOneFactory();

  // see https://itk.org/Doxygen46/html/IO_2DicomSeriesReadImageWrite2_8cxx-example.html
  typedef signed short    PixelType;
  const unsigned int      Dimension = 3;
  typedef itk::Image<PixelType, Dimension>     ImageType;
 // typedef itk::Image<PixelType, Dimension>     OutputImageType;

  typedef itk::ImageSeriesReader<ImageType >        ReaderType;
  ReaderType::Pointer reader = ReaderType::New();

 // A GDCMImageIO object is created and connected to the reader. This object is
 // the one that is aware of the internal intricacies of the DICOM format.
   typedef itk::GDCMImageIO     ImageIOType;
   ImageIOType::Pointer dicomIO = ImageIOType::New();
   reader->SetImageIO( dicomIO );
   
   typedef itk::GDCMSeriesFileNames NamesGeneratorType;
   NamesGeneratorType::Pointer nameGenerator = NamesGeneratorType::New();

   nameGenerator->SetUseSeriesDetails( true );
   //nameGenerator->AddSeriesRestriction("0008|0021" );
   nameGenerator->SetDirectory("/Users/sonia/Projects/internshipINL/itk_testapp/tests/image_function/new_stack/");

try{

      std::cout<<"DICOM" <<std::endl;
    std::cout<<"/////////////////////////////" <<std::endl;
    std::cout << std::endl << "The directory: " << std::endl;
    std::cout << std::endl << "/Users/sonia/Projects/internshipINL/itk_testapp/tests/image_function/new_stack/"  << std::endl;
    std::cout << "Contains the following DICOM Series: ";
    std::cout << std::endl << std::endl;

   typedef std::vector< std::string >    SeriesIdContainer;
    const SeriesIdContainer & seriesUID = nameGenerator->GetSeriesUIDs();
    SeriesIdContainer::const_iterator seriesItr = seriesUID.begin();
    SeriesIdContainer::const_iterator seriesEnd = seriesUID.end();
    while( seriesItr != seriesEnd )
      {
      std::cout << seriesItr->c_str() << std::endl;
      ++seriesItr;
      }


     std::string seriesIdentifier;
    if( 0) // If no optional series identifier
      {
      seriesIdentifier = "se definiamo da fuori la serie da leggere";
      }
    else
      {
      seriesIdentifier = seriesUID.begin()->c_str();
      }

    std::cout << std::endl << std::endl;
    std::cout << "Now reading series: " << std::endl << std::endl;
    std::cout << seriesIdentifier << std::endl << std::endl;

     ////////////////////////////////////////////////////////

    typedef std::vector< std::string >   FileNamesContainer;
    FileNamesContainer fileNames = nameGenerator->GetFileNames( seriesIdentifier );
    reader->SetFileNames( fileNames );

    try
      {
      reader->Update();
      reader->GetOutput();
      }
    catch(itk::ExceptionObject &ex)
      {
      std::cout << ex << std::endl;
      }
  
  ///////// info from initial image /////
    std::cout<<"DICOM Serie info" <<std::endl;
    std::cout<<"/////////////////////////////" <<std::endl;
    ImageType::SizeType imageSize =reader->GetOutput()->GetLargestPossibleRegion().GetSize();
    std::cout<<"Dimensions =" <<imageSize<<std::endl;
    const ImageType::SpacingType& inputSpacing =reader->GetOutput()->GetSpacing();
    std::cout << "Spacing = " << inputSpacing << std::endl;
    const ImageType::PointType & origin = reader->GetOutput()->GetOrigin();
    std::cout << "Origin = "<<  origin << std::endl;
    std::cout<<"/////////////////////////////" <<std::endl;

  ////// apply the filter 

  typedef itk::RescaleIntensityImageFilter< ImageType, WriteImageType > RescaleFilterType;
  RescaleFilterType::Pointer rescaler = RescaleFilterType::New();
  rescaler->SetOutputMinimum(   0 );
  rescaler->SetOutputMaximum( 255 );


   /////////write filtred image





  rescaler->SetInput( reader->GetOutput() );
  scaledImage=(rescaler->GetOutput());
  scaledImage->SetSpacing( inputSpacing );
   scaledImage->GetLargestPossibleRegion();



      for (unsigned int i = 0; i < 3; ++i){
     // _voxel2.push_back(_physical_dims(i) / imageSize[i]); //uncomment
      _voxel2.push_back(5.0);  //comment
       std::cout<<"voxel dopo"<<std::endl;
        std::cout<<_voxel2[i]<<std::endl;
      }

       ImageType::SpacingType spacing2;
       spacing2[0] = _voxel2[0]; // spacing along X
       spacing2[1] = _voxel2[1];
       spacing2[2] = _voxel2[2];
       scaledImage->SetSpacing(spacing2);


      ImageType::PointType newOrigin;
      newOrigin[0] = _origin(0); // spacing along X
       newOrigin[1] = _origin(1);
       newOrigin[2] = _origin(2);

      scaledImage->SetOrigin( newOrigin);
      std::cout << "New Origin = "<<  newOrigin << std::endl;


  writer->SetFileName("outputFilenameFiltred.tiff" );
 writer->SetInput( scaledImage );

     try
     {
     writer->Update();
     }
    catch (itk::ExceptionObject & e)
     {
     std::cerr << e << std::endl;
     mooseError("Exception in file writer");
     }



std::cout<<"IMMAGINE PRODOTTA"<<std::endl;




    unsigned int n =  scaledImage->GetNumberOfComponentsPerPixel();
    std::cout<<"Number Of Components Per Pixel: "<<n <<std::endl;
 

 }

catch (itk::ExceptionObject &ex){
    mooseError("exception in reading dicom series");

}


//////////////////////////////

}



Real
ImageSamplerItk::itksample(const Point & p)
{


  writer2->SetFileName("outputFilenameFiltred2.tiff" );
 writer2->SetInput( scaledImage );

     try
     {
     writer2->Update();
     }
    catch (itk::ExceptionObject & e)
     {
     std::cerr << e << std::endl;
     mooseError("Exception in file writer");
     }



//std::cout<<"IMMAGINE PRODOTTA"<<std::endl;




  // Do nothing if the point is outside of the image domain
  if (!_bounding_box.contains_point(p))
    return 0.0;

  // Determine pixel coordinates
  std::vector<int> x(3, 0);
  for (int i = 0; i < LIBMESH_DIM; ++i)
  {
    // Compute position, only if voxel size is greater than zero
    if (_voxel[i] == 0){
      x[i] = 0;
      // std::cout<<"printo x nel primp caso"<<  x[i]<<std::endl;
     }

    else
    {
      x[i] = std::floor((p(i) - _origin(i)) / _voxel[i]);

      // std::cout<<"printo p nel secondo caso: "<<  p(i)<<std::endl;

      // If the point falls on the mesh extents the index needs to be decreased by one
      if (x[i] == _dims[i])
        x[i]--;
    }
  }

  
// if (p(0)==0.8125 &&p(1)==0.9375 &&p(2)==1.0 ){

//   std::cout<<"ITK SAMPLE"<<std::endl;
//   std::cout<<"punto"<<p<<std::endl;
//    std::cout<<"coordinate: "<<x[0]<<", "<<x[1]<<", "<<x[2]<<std::endl;}


 WriteImageType::IndexType pixelIndex;

  typedef itk::Point< double, WriteImageType::ImageDimension > PointType;
  PointType coordinate;
  coordinate[0]=x[0];
  coordinate[1]=x[1];
  coordinate[2]=x[2];

 
const bool isInside2 =scaledImage->TransformPhysicalPointToIndex(coordinate , pixelIndex );
 if ( isInside2 )
    {
         std::cout<<"Coordinates"<<coordinate<<std::endl;
    WriteImageType::PixelType pixelValue = scaledImage->GetPixel( pixelIndex);

     std::cout<<"pixel"<<pixelValue<<std::endl;
  }
  
  //std::cout<< _data->GetScalarComponentAsDouble(x[0], x[1], x[2], _component)<<std::endl;



  // Return the image data at the given point
  return _data->GetScalarComponentAsDouble(x[0], x[1], x[2], _component);


}
