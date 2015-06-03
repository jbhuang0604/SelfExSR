Mean Shift Image Processor Class ver1.0 README
----------------------------------------------

Class Overview:
===============

The mean shift image processor class is designed to offer the following functionality:

	(1) Perform image segmentation and edge-preserving filtering using the mean shift algorithm.
	(2) Perform (1) using a general kernel and/or an arbitrary input data space.

Table of Contents:
-------------------
(A) Image Segmentation and Filtering
(B) Synergistic Image Segmentation
(C) Using a General Kernel
(D) Using an Arbitrary Input Data Space
(E) The Class Error Handler
(F) Current Version Information
(G) References
(H) Contact Information

================================================================================================

(A) Image Segmentation and Filtering

-------------------------------------------------------------------------------------------------

Mean shift based image segmentation and filtering is performed using use the following methods:

	msImageProcess::Filter		- filters the image
	msImageProcessor::Segment	- segments the image

The input image processed by these methods is defined via the method,

msImageProcessor::DefineImage	- this uploads the RGB data into the msImageProcessor class for processing

To obtain the output call:

	msImageProcessor::GetResults	- returns filtered or segmented image in RGB space
	msImageProcessor::GetBoundaries	- returns the boundaries of regions resulting from filtering
					  or segmentation
	msImageProcessor::GetRegions	- returns the classification structure that maps each
					  data point in the image to a given mode, and also
					  the number of points in the image correlating to each mode.

NOTE:
-----

The modes returned by GetRegions are not in the RGB space. If DefineImage was used, they are in the LUV space. The modes may be converted from LUV to RGB (and visa versa) using the space conversion methods of the msImageProcessor class:

	msImageProcessor::RGBtoLUV	- converts data points from the RGB data space to LUV
	msImageProcessor::LUVtoRGB	- converts data points from the LUV data space to RGB

Alternatively, mean shift may be applyed to data that lies in a space other than LUV. This may be accomplished through the use of the method MeanShift::DefineLInput (see section D).

================================================================================================

(B) Synergistic Image Segmentation

-------------------------------------------------------------------------------------------------

A weight map may be provided to the mean shift image processor class, used to perform synergistic image segmentation as described in the paper [3]. One may specify a weight map by calling either of the following methods:

	MeanShift::SetWeightMap			- defines the weight map used to specify a weighted kernel during
						  mean shift; the weight map may only be used for data that lies
						  on a lattice (e.g. an image)
	msImageProcessor::SetWeightMap		- specifies a weight map to be used for performing synergistic image
						  segmentation

Each of the above methods accept a floating point array of size L elements containing the weight map. When using the mean shift base class L is the number of data points in the specified data set; when using the image processor class L = height x width, where height and width are the dimensions of the image. The method msImageProcessor::SetWeightMap accepts an additional parameter, namely t_e, a threshold value used during the transitive closure step of the image segmentation algorithm. See the paper [3] for details.

================================================================================================

(C) Using a General Kernel

-------------------------------------------------------------------------------------------------

A general kernel can be used to perform mean shift filtering and segmentation by calling the inherited method:

	MeanShift::DefineKernel	- defines an N-dimensional kernel having kp subspaces, in which each subspace
				  can be of one of three types: Uniform, Gaussian, or UserDefined.

DefineImage, used to define the input image when performing image segmentation or filtering, defines a Uniform kernel having two subspaces (one spatial (x,y) and one range (L,U,V)) each subspace having bandwidths sigmaS and sigmaR respectively. By skimming the method definition one may get an idea of how to define a general kernel.

NOTE:
----

For data that is defined on a lattice, it is always assumed that the spatial domain is treated as a single subspace. Also, DefineKernel() must be called *after* DefineImage() when these methods are used together.

================================================================================================

(D) Using an Arbitrary Input Data Space

-------------------------------------------------------------------------------------------------

Mean shift filtering and segmentation may be performed on an arbitary image data space. Such data is defined through calling the inherited method:

	MeanShift::DefineLInput	- specifies input defined on a lattice

DefineImage() calls this method using the LUV data it generates. Through the use of the above methods, mean shift may be applied to an arbitrary input data space using a general kernel. In doing so, one must ensure that the dimension of the input data space and kernel are the same (N). If their dimensions do not agree an error will be flagged.

================================================================================================

(F) The Class Error Handler

-------------------------------------------------------------------------------------------------

The mean shift image processor class uses an error message string and error-level flag to perform error handling. These two variables, MeanShift::ErrorMessage and MeanShift::ErrorLevel, are public data members of the class.

Upon the occurance of an error,

	* An error message is copied into the error message string.
	* The error level of the class is set to EL_ERROR.

The following example demonstrates the use of the error handling mechanism described above.

msImageProcessor iProc;

...

iProc.Segment(sigmaS, sigmaR, minRegion, SPEEDUP);
if(iProc.ErrorLevel == EL_ERROR)
{
	fprintf(stderr, iProc.ErrorMessage);	
	exit(1);
}

...

================================================================================================

(G) Current Version Information

-------------------------------------------------------------------------------------------------

The current version of the code was tested under both UNIX and Windows environments.

================================================================================================

(H) References

-------------------------------------------------------------------------------------------------

[1] D. Comanicu, P. Meer: "Mean shift: A robust approach toward feature space analysis".
    IEEE Trans. Pattern Anal. Machine Intell., May 2002.

[2] P. Meer, B. Georgescu: "Edge detection with embedded confidence". IEEE Trans. Pattern Anal.
    Machine Intell., 28, 2001.

[3] C. Christoudias, B. Georgescu, P. Meer: "Synergism in low level vision". 16th International
    Conference of Pattern Recognition, Track 1 - Computer Vision and Robotics, Quebec City,
    Canada, August 2001.

================================================================================================

(I) Contact Information

-------------------------------------------------------------------------------------------------

Personal Contact Information
----------------------------

Email:

	cmch@caip.rutgers.edu		(Chris Christoudias)
	georgesc@caip.rutgers.edu	(Bogdan Georgescu)

Laboratory Contact Information
------------------------------

Laboratory Website:

	www.caip.rutgers.edu/riul/

================================================================================================