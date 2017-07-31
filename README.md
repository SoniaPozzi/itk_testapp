ItkTestApp
=====

"Fork ItkTestApp" to create a new MOOSE-based application.

For more information see: [http://mooseframework.org/create-an-app/](http://mooseframework.org/create-an-app/)


Inside this application are created the following:

FileDicomChoose:
	-Mainly based on itk, FileDicomChoose open and filters the chosen itk series.


A Mesh Constructor (ImageMeshItk):
	-Based on FileDicomChoose, ImageMeshItk construct a mesh starting from the cropped and filtered series' dimension.
	
	

A Function (ImageFunctionItk):
	 - Based on ImageSampler,  ImageFunctionItk fills the mesh with the pixels values from the cropped and filtered series . Based on FileDicomChoose,  ImageSampler fills the mesh with the pixels values from the cropped and filtered series .
	
 
 
 A Mesh Modifier  (ImageThresholdMesh):
 	- Based on ImageFunctionItk, ImageThresholdMesh  deletes from the mesh all the elements with image function  value outside some thresholded region.
	
	

	
		