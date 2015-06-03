/*******************************************************

                 Mean Shift Analysis Library
	=============================================


	The mean shift library is a collection of routines
	that use the mean shift algorithm. Using this algorithm,
	the necessary output will be generated needed
	to analyze a given input set of data.

  Mean Shift Image Processor Class:
  ================================

	The following class inherits from the mean shift library
	in order to perform the specialized tasks of image
	segmentation and filtering.
	
	The prototype of the Mean Shift	Image Processor Class
	is provided below. Its definition is provided in
	'msImageProcessor.cc'.

The theory is described in the papers:

  D. Comaniciu, P. Meer: Mean Shift: A robust approach toward feature
									 space analysis.

  C. Christoudias, B. Georgescu, P. Meer: Synergism in low level vision.

and they are is available at:
  http://www.caip.rutgers.edu/riul/research/papers/

Implemented by Chris M. Christoudias, Bogdan Georgescu
********************************************************/

#ifndef msImageProcessor_H
#define msImageProcessor_H

//include mean shift library
#include	"ms.h"

//include prototypes of additional strucuters
//used for image segmentation...

//include region list used to store boundary pixel
//indeces for each region
#include	"rlist.h"

//include region adjacency list class used for
//region pruning and transitive closure
#include	"RAList.h"

//define constants

	//image pruning
#define	TOTAL_ITERATIONS	14
#define BIG_NUM				0xffffffff	//BIG_NUM = 2^32-1
#define NODE_MULTIPLE		10

	//data space conversion...
const double Xn			= 0.95050;
const double Yn			= 1.00000;
const double Zn			= 1.08870;
//const double Un_prime	= 0.19780;
//const double Vn_prime	= 0.46830;
const double Un_prime	= 0.19784977571475;
const double Vn_prime	= 0.46834507665248;
const double Lt			= 0.008856;

	//RGB to LUV conversion
const double XYZ[3][3] = {	{  0.4125,  0.3576,  0.1804 },
							{  0.2125,  0.7154,  0.0721 },
							{  0.0193,  0.1192,  0.9502 }	};

	//LUV to RGB conversion
const double RGB[3][3] = {	{  3.2405, -1.5371, -0.4985 },
							{ -0.9693,  1.8760,  0.0416 },
							{  0.0556, -0.2040,  1.0573 }	};

//define data types
typedef unsigned char byte;

//define enumerations
enum imageType {GRAYSCALE, COLOR};

//define prototype
class msImageProcessor: public MeanShift {

public:

  /*/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\*/
  /* Class Constructor and Destructor */
  /*\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/*/

  msImageProcessor( void );        //Default Constructor
 ~msImageProcessor( void );        //Class Destructor

 /*/\/\/\/\/\/\/\/\/\/\/\/\/\*/
 /* Input Image Declaration  */
 /*\/\/\/\/\/\/\/\/\/\/\/\/\/*/

  //--\\||//--\\||//--\\||//--\\||//--\\||//--\\||//--\\||//
  //<--------------------------------------------------->|//
  //|                                                    |//
  //|	Method Name:								     |//
  //|   ============								     |//
  //|				  * Define Image *                   |//
  //|                                                    |//
  //<--------------------------------------------------->|//
  //|                                                    |//
  //|	Description:								     |//
  //|	============								     |//
  //|                                                    |//
  //|   Uploads an image to be segmented by the image    |//
  //|   segmenter class.                                 |//
  //|                                                    |//
  //|   An image is defined by specifying the folloing:  |//
  //|                                                    |//
  //|   <* data *>                                       |//
  //|   A one dimensional unsigned char array of RGB     |//
  //|   vectors.                                         |//
  //|                                                    |//
  //|   <* type *>                                       |//
  //|   Specifies the image type: COLOR or GREYSCALE.    |//
  //|                                                    |//
  //|   <* height *>                                     |//
  //|   The image height.                                |//
  //|                                                    |//
  //|   <* width *>                                      |//
  //|   The image width.                                 |//
  //|                                                    |//
  //|   This method uploads the image and converts its   |//
  //|   data into the LUV space. If another conversion   |//
  //|   is desired data may be uploaded into this class  |//
  //|   via the procedure MeanShift::UploadInput().      |//
  //|                                                    |//
  //<--------------------------------------------------->|//
  //|                                                    |//
  //|	Usage:      								     |//
  //|   ======      								     |//
  //|		DefineImage(data, type, height, width)       |//
  //|                                                    |//
  //<--------------------------------------------------->|//
  //--\\||//--\\||//--\\||//--\\||//--\\||//--\\||//--\\||//

  void DefineImage(byte*,imageType, int, int);
  void DefineBgImage(byte*, imageType , int , int );


 /*/\/\/\/\/\/\*/
 /* Weight Map */
 /*\/\/\/\/\/\/*/

  //--\\||//--\\||//--\\||//--\\||//--\\||//--\\||//--\\||//
  //<--------------------------------------------------->|//
  //|                                                    |//
  //|	Method Name:								     |//
  //|   ============								     |//
  //|			     * Set Weight Map *                  |//
  //|                                                    |//
  //<--------------------------------------------------->|//
  //|                                                    |//
  //|	Description:								     |//
  //|	============								     |//
  //|                                                    |//
  //|   Uploads weight map specifying for each pixel     |//
  //|   in the image a value between 0 and 1 - 1 indica- |//
  //|   ting the presence of an edge and 0 the absense   |//
  //|   of an edge.                                      |//
  //|                                                    |//
  //|   The arguments to this method are:                |//
  //|                                                    |//
  //|   <* weightMap *>                                  |//
  //|   A floating point array of size (height x width)  |//
  //|   specifying at location (i,j) the edge strength   |//
  //|   of pixel (i,j). (e.g. pixel (i,j) has an edge    |//
  //|   strength of weightMap[j*width+i]).               |//
  //|                                                    |//
  //|   <* epsilon *>                                    |//
  //|   A floating point number specifying the threshold |//
  //|   used to fuse regions during transitive closure.  |//
  //|                                                    |//
  //|   Note: DefineImage must be called prior to call-  |//
  //|         ing this method. DefineImage is used to    |//
  //|         define the dimensions of the image.        |//
  //|                                                    |//
  //<--------------------------------------------------->|//
  //|                                                    |//
  //|	Usage:      								     |//
  //|   ======      								     |//
  //|		SetWeightMap(weightMap, epsilon)             |//
  //|                                                    |//
  //<--------------------------------------------------->|//
  //--\\||//--\\||//--\\||//--\\||//--\\||//--\\||//--\\||//

  void SetWeightMap(float*, float);

  //--\\||//--\\||//--\\||//--\\||//--\\||//--\\||//--\\||//
  //<--------------------------------------------------->|//
  //|                                                    |//
  //|	Method Name:								     |//
  //|   ============								     |//
  //|			   * Remove Weight Map *                 |//
  //|                                                    |//
  //<--------------------------------------------------->|//
  //|                                                    |//
  //|	Description:								     |//
  //|	============								     |//
  //|                                                    |//
  //|   Removes weight map. An error is NOT flagged      |//
  //|   if a weight map was not defined prior to calling |//
  //|   this method.                                     |//
  //|                                                    |//
  //<--------------------------------------------------->|//
  //|                                                    |//
  //|	Usage:      								     |//
  //|   ======      								     |//
  //|		RemoveWeightMap(void)                        |//
  //|                                                    |//
  //<--------------------------------------------------->|//
  //--\\||//--\\||//--\\||//--\\||//--\\||//--\\||//--\\||//

  void RemoveWeightMap(void);

 /*/\/\/\/\/\/\/\/\/\*/
 /* Image Filtering  */
 /*\/\/\/\/\/\/\/\/\/*/

  //--\\||//--\\||//--\\||//--\\||//--\\||//--\\||//--\\||//
  //<--------------------------------------------------->|//
  //|                                                    |//
  //|	Method Name:								     |//
  //|   ============								     |//
  //|                   *  Filter  *                     |//
  //|                                                    |//
  //<--------------------------------------------------->|//
  //|                                                    |//
  //|	Description:								     |//
  //|	============								     |//
  //|                                                    |//
  //|   Apply mean shift filter to the defined image,    |//
  //|   defined either via MeanShift::DefineLInput or    |//
  //|   msImageProcessor::DefineImage. The resulting     |//
  //|   segmented image is stored in the private data    |//
  //|   members of the image segmenter class which can   |//
  //|   be obtained by calling image msImageProcessor::  |//
  //|   GetResults().                                    |//
  //|                                                    |//
  //|   The arguments to this method are:                |//
  //|                                                    |//
  //|   <* sigmaS *>                                     |//
  //|   The spatial radius of the mean shift window.     |//
  //|                                                    |//
  //|   <* sigmaR *>                                     |//
  //|   The range radius of the mean shift window.       |//
  //|                                                    |//
  //|   <* speedUpLevel *>                               |//
  //|   Determines if a speed up optimization should be  |//
  //|   used to perform image filtering. A value of      |//
  //|   NO_SPEEDUP turns this optimization off and a     |//
  //|   value of SPEEDUP turns this optimization on.     |//
  //|                                                    |//
  //<--------------------------------------------------->|//
  //|                                                    |//
  //|	Usage:      								     |//
  //|   ======      								     |//
  //|		Filter(sigmaS, sigmaR, speedUpLevel)         |//
  //|                                                    |//
  //<--------------------------------------------------->|//
  //--\\||//--\\||//--\\||//--\\||//--\\||//--\\||//--\\||//

  void Filter(int, float, SpeedUpLevel);

 /*/\/\/\/\/\/\/\/\/\/\/\*/
 /* Image Region Fusing  */
 /*\/\/\/\/\/\/\/\/\/\/\/*/

  //--\\||//--\\||//--\\||//--\\||//--\\||//--\\||//--\\||//
  //<--------------------------------------------------->|//
  //|                                                    |//
  //|	Method Name:								     |//
  //|   ============								     |//
  //|				  *  Fuse Regions  *                 |//
  //|                                                    |//
  //<--------------------------------------------------->|//
  //|                                                    |//
  //|	Description:								     |//
  //|	============								     |//
  //|                                                    |//
  //|   Fuses the regions of a filtered image,           |//
  //|   defined either via MeanShift::DefineLInput or    |//
  //|   msImageProcessor::DefineImage. The resulting     |//
  //|   segmented image is stored in the private data    |//
  //|   members of the image segmenter class which can   |//
  //|   be obtained by calling image msImageProcessor::  |//
  //|   GetResults().                                    |//
  //|                                                    |//
  //|   The arguments to this method are:                |//
  //|                                                    |//
  //|   <* sigmaR *>                                     |//
  //|   The range radius that defines similar color      |//
  //|   amongst image regions.                           |//
  //|                                                    |//
  //|   <* minRegion *>                                  |//
  //|   The minimum density a region may have in the     |//
  //|   resulting segmented image. All regions have      |//
  //|   point density < minRegion are pruned from the    |//
  //|   image.                                           |//
  //|                                                    |//
  //<--------------------------------------------------->|//
  //|                                                    |//
  //|	Usage:      								     |//
  //|   ======      								     |//
  //|		FuseRegions(sigmaR, minRegion)               |//
  //|                                                    |//
  //<--------------------------------------------------->|//
  //--\\||//--\\||//--\\||//--\\||//--\\||//--\\||//--\\||//

  void FuseRegions(float, int);

 /*/\/\/\/\/\/\/\/\/\/\*/
 /* Image Segmentation */
 /*\/\/\/\/\/\/\/\/\/\/*/

  //--\\||//--\\||//--\\||//--\\||//--\\||//--\\||//--\\||//
  //<--------------------------------------------------->|//
  //|                                                    |//
  //|	Method Name:								     |//
  //|   ============								     |//
  //|				     *  Segment  *                   |//
  //|                                                    |//
  //<--------------------------------------------------->|//
  //|                                                    |//
  //|	Description:								     |//
  //|	============								     |//
  //|                                                    |//
  //|   Segments the defined image, defined either via   |//
  //|   MeanShift::DefineLInput or msImageProcessor::De- |//
  //|   fineImage. The resulting segmented image is      |//
  //|   stored in the private data members of the image  |//
  //|   processor class which can be obtained by calling |//
  //|   ImageSegmenter::GetResults().                    |//
  //|                                                    |//
  //|   The arguments to this method are:                |//
  //|                                                    |//
  //|   <* sigmaS *>                                     |//
  //|   The spatial radius of the mean shift window.     |//
  //|                                                    |//
  //|   <* sigmaR *>                                     |//
  //|   The range radius of the mean shift window.       |//
  //|                                                    |//
  //|   <* minRegion *>                                  |//
  //|   The minimum density a region may have in the     |//
  //|   resulting segmented image. All regions have      |//
  //|   point density < minRegion are pruned from the    |//
  //|   image.                                           |//
  //|                                                    |//
  //|   <* speedUpLevel *>                               |//
  //|   Determines if a speed up optimization should be  |//
  //|   used to perform image filtering. A value of      |//
  //|   NO_SPEEDUP turns this optimization off and a     |//
  //|   value of SPEEDUP turns this optimization on.     |//
  //|                                                    |//
  //<--------------------------------------------------->|//
  //|                                                    |//
  //|	Usage:      								     |//
  //|   ======      								     |//
  //|		Segment(sigmaS, sigmaR, minRegion,           |//
  //|                       speedUpLevel)                |//
  //|                                                    |//
  //<--------------------------------------------------->|//
  //--\\||//--\\||//--\\||//--\\||//--\\||//--\\||//--\\||//
 
  void Segment(int, float, int, SpeedUpLevel);

  /*/\/\/\/\/\/\/\/\/\/\/\/\*/
  /* Data Space Conversion  */
  /*\/\/\/\/\/\/\/\/\/\/\/\/*/

  //--\\||//--\\||//--\\||//--\\||//--\\||//--\\||//--\\||//
  //<--------------------------------------------------->|//
  //|                                                    |//
  //|	Method Name:								     |//
  //|   ============								     |//
  //|				 *  RGB To LUV  *                    |//
  //|                                                    |//
  //<--------------------------------------------------->|//
  //|                                                    |//
  //|	Description:								     |//
  //|	============								     |//
  //|                                                    |//
  //|   Converts an RGB vector to LUV.                   |//
  //|                                                    |//
  //|   The arguments to this method are:                |//
  //|                                                    |//
  //|   <* rgbVal *>                                     |//
  //|   An unsigned char array containing the RGB        |//
  //|   vector.                                          |//
  //|                                                    |//
  //|   <* luvVal *>                                     |//
  //|   A floating point array containing the LUV        |//
  //|   vector.                                          |//
  //|                                                    |//
  //<--------------------------------------------------->|//
  //|                                                    |//
  //|	Usage:      								     |//
  //|   ======      								     |//
  //|		RGBtoLUV(rgbVal, luvVal)                     |//
  //|                                                    |//
  //<--------------------------------------------------->|//
  //--\\||//--\\||//--\\||//--\\||//--\\||//--\\||//--\\||//

  void RGBtoLUV(byte*, float*);

  //--\\||//--\\||//--\\||//--\\||//--\\||//--\\||//--\\||//
  //<--------------------------------------------------->|//
  //|                                                    |//
  //|	Method Name:								     |//
  //|   ============								     |//
  //|				 *  LUV To RGB  *                    |//
  //|                                                    |//
  //<--------------------------------------------------->|//
  //|                                                    |//
  //|	Description:								     |//
  //|	============								     |//
  //|                                                    |//
  //|   Converts an LUV vector to RGB.                   |//
  //|                                                    |//
  //|   The arguments to this method are:                |//
  //|                                                    |//
  //|   <* luvVal *>                                     |//
  //|   A floating point array containing the LUV        |//
  //|   vector.                                          |//
  //|                                                    |//
  //|   <* rgbVal *>                                     |//
  //|   An unsigned char array containing the RGB        |//
  //|   vector.                                          |//
  //|                                                    |//
  //<--------------------------------------------------->|//
  //|                                                    |//
  //|	Usage:      								     |//
  //|   ======      								     |//
  //|		LUVtoRGB(luvVal, rgbVal)                     |//
  //|                                                    |//
  //<--------------------------------------------------->|//
  //--\\||//--\\||//--\\||//--\\||//--\\||//--\\||//--\\||//

  void LUVtoRGB(float*, byte*);

  /*/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\*/
  /*  Filtered and Segmented Image Output */
  /*\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/*/

  //--\\||//--\\||//--\\||//--\\||//--\\||//--\\||//--\\||//
  //<--------------------------------------------------->|//
  //|                                                    |//
  //|	Method Name:								     |//
  //|   ============								     |//
  //|			      *  Get Raw Data  *                 |//
  //|                                                    |//
  //<--------------------------------------------------->|//
  //|                                                    |//
  //|	Description:								     |//
  //|	============								     |//
  //|                                                    |//
  //|   Returns the resulting filtered or segmented im-  |//
  //|   age data after calling Filter() or Segment().    |//
  //|                                                    |//
  //|   The arguments to this method are:                |//
  //|                                                    |//
  //|   <* outputImageData *>                            |//
  //|   A floating point array containing the vector     |//
  //|   data of the filtered or segmented image.         |//
  //|                                                    |//
  //|   NOTE: If DefineImage was used to specify the     |//
  //|         the input to this class, outputImageData   |//
  //|         is in the LUV data space.                  |//
  //|                                                    |//
  //<--------------------------------------------------->|//
  //|                                                    |//
  //|	Usage:      								     |//
  //|   ======      								     |//
  //|		GetResults(outputImageData)                  |//
  //|                                                    |//
  //<--------------------------------------------------->|//
  //--\\||//--\\||//--\\||//--\\||//--\\||//--\\||//--\\||//

  void GetRawData(float*);

  //--\\||//--\\||//--\\||//--\\||//--\\||//--\\||//--\\||//
  //<--------------------------------------------------->|//
  //|                                                    |//
  //|	Method Name:								     |//
  //|   ============								     |//
  //|				 *  Get Results  *                   |//
  //|                                                    |//
  //<--------------------------------------------------->|//
  //|                                                    |//
  //|	Description:								     |//
  //|	============								     |//
  //|                                                    |//
  //|   Returns the resulting filtered or segmented im-  |//
  //|   age after calling Filter() or Segment().         |//
  //|                                                    |//
  //|   The arguments to this method are:                |//
  //|                                                    |//
  //|   <* outputImage *>                                |//
  //|   An unsigned char array containing the RGB        |//
  //|   vector data of the output image.                 |//
  //|                                                    |//
  //|   To obtain the un-converted (LUV) data space      |//
  //|   output one may use                               |//
  //|   msImageProcessor::GetRawData().                  |//
  //<--------------------------------------------------->|//
  //|                                                    |//
  //|	Usage:      								     |//
  //|   ======      								     |//
  //|		GetResults(outputImage)                      |//
  //|                                                    |//
  //<--------------------------------------------------->|//
  //--\\||//--\\||//--\\||//--\\||//--\\||//--\\||//--\\||//

  void GetResults(byte*);

  //--\\||//--\\||//--\\||//--\\||//--\\||//--\\||//--\\||//
  //<--------------------------------------------------->|//
  //|                                                    |//
  //|	Method Name:								     |//
  //|   ============								     |//
  //|				 *  Get Boundaries  *                |//
  //|                                                    |//
  //<--------------------------------------------------->|//
  //|                                                    |//
  //|	Description:								     |//
  //|	============								     |//
  //|                                                    |//
  //|   Returns the boundaries of each region of the     |//
  //|   segmented image using a region list object,      |//
  //|   available after filtering or segmenting the      |//
  //|   defined image.                                   |//
  //|                                                    |//
  //<--------------------------------------------------->|//
  //|                                                    |//
  //|	Usage:      								     |//
  //|   ======      								     |//
  //|		regionList = GetBoundaries()                 |//
  //|                                                    |//
  //<--------------------------------------------------->|//
  //--\\||//--\\||//--\\||//--\\||//--\\||//--\\||//--\\||//

  RegionList *GetBoundaries( void );

  //--\\||//--\\||//--\\||//--\\||//--\\||//--\\||//--\\||//
  //<--------------------------------------------------->|//
  //|                                                    |//
  //|	Method Name:								     |//
  //|   ============								     |//
  //|			        * Get Regions *                  |//
  //|                                                    |//
  //<--------------------------------------------------->|//
  //|                                                    |//
  //|	Description:								     |//
  //|	============								     |//
  //|                                                    |//
  //|   Returns the regions of the processed image.      |//
  //|   Each region in the image is uniquely character-  |//
  //|   ized by its location and color (e.g. RGB).       |//
  //|   GetRegions() therefore returns the following     |//
  //|   information about the regions of the processed   |//
  //|   image:                                           |//
  //|                                                    |//
  //|   <* regionCount *>                                |//
  //|   An integer that specifies the number of regions  |//
  //|   contained in the processed image.                |//
  //|                                                    |//
  //|   <* modes *>                                      |//
  //|   A floating point array of length regionCount*N   |//
  //|   containing the feature space component of each   |//
  //|   region (e.g. LUV), and indexed by region label.  |//
  //|                                                    |//
  //|   <* labels *>                                     |//
  //|   An integer array of length (height*width) which  |//
  //|   contains at every pixel location (x,y) a label   |//
  //|   relating that pixel to a region whose mode is    |//
  //|   specified by modes and whose area is specified   |//
  //|   by modePointCounts.                              |//
  //|                                                    |//
  //|   <* modePointCounts *>                            |//
  //|   An integer array of length regionCount and ind-  |//
  //|   exed by region label, that specifies the region  |//
  //|   area (in pixels) for each segmented image reg-   |//
  //|   ion. (e.g. Area of region having label specif-   |//
  //|   ified by l, has area modePointCounts[l] (pix-    |//
  //|   els).)                                           |//
  //|                                                    |//
  //|   NOTE: Memory for the above integer and floating  |//
  //|         point arrays is allocated inside this      |//
  //|         method.                                    |//
  //|                                                    |//
  //|         Also modes stored by the modes array are   |//
  //|         not in the RGB space. Instead if the       |//
  //|         method DefineImage was used, these data    |//
  //|         points are in the LUV space, and if the    |//
  //|         method DefineLInput was used these data    |//
  //|         points are in whatever space you specified |//
  //|         them to be in when calling DefineLInput.   |//
  //|                                                    |//
  //<--------------------------------------------------->|//
  //|                                                    |//
  //|	Usage:      								     |//
  //|   ======      								     |//
  //|		regionCount = GetRegions(labels, modes       |//
  //|                                modePointCounts)    |//
  //|                                                    |//
  //<--------------------------------------------------->|//
  //--\\||//--\\||//--\\||//--\\||//--\\||//--\\||//--\\||//

  int GetRegions(int**, float**, int**);
  int *GetLabels() { return labels; };


  void SetSpeedThreshold(float);
private:

  //========================
  // *** Private Methods ***
  //========================

	/*/\/\/\/\/\/\/\/\/\*/
	/*  Image Filtering */
	/*\/\/\/\/\/\/\/\/\/*/

	void NonOptimizedFilter(float, float);	// filters the image applying mean shift to each point
											// Advantage	: most accurate
											// Disadvantage	: time expensive
   void NewNonOptimizedFilter(float, float);

	void OptimizedFilter1(float, float);	// filters the image using previous mode information
											// to avoid re-applying mean shift to some data points
											// Advantage	: maintains high level of accuracy,
											//				  large speed up compared to non-optimized
											//				  version
											// Disadvantage	: POSSIBLY not as accurate as non-optimized
											//				  version
   void NewOptimizedFilter1(float, float);


	void OptimizedFilter2(float, float);	//filter the image using previous mode information
											//and window traversals to avoid re-applying mean shift to
											//some data points
											// Advantage	: huge speed up - maintains accuracy good enough
											//				  for segmentation
											// Disadvantage	: not as accurate as previous filters
   void NewOptimizedFilter2(float, float);

	
	/*/\/\/\/\/\/\/\/\/\/\/\*/
	/* Image Classification */
	/*\/\/\/\/\/\/\/\/\/\/\/*/

	void Connect( void );					// classifies mean shift filtered image regions using
											// private classification structure of this class

	void Fill(int, int);					// used by Connect to perform label each region in the
											// mean shift filtered image using an eight-connected
											// fill

	/*/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\*/
	/* Transitive Closure and Image Pruning */
	/*\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/*/

	void BuildRAM( void );					// build a region adjacency matrix using the region list
											// object

	void DestroyRAM( void );				// destroy the region adjacency matrix: de-allocate its memory
											// initialize it for re-use

	void TransitiveClosure( void );			// use the RAM to apply transitive closure to the image modes

	void ComputeEdgeStrengths( void );		// computes the weights of the weighted graph using the weight
											// map

	//Usage: Prune(minRegion)
	void Prune(int);						// use the RAM to prune the image of spurious regions (regions
											// whose area is less than minRegion pixels, where minRegion is
											// an argument of this method)

	/*/\/\/\/\/\/\/\/\/\/\/\/\/\/\*/
	/*  Region Boundary Detection */
	/*\/\/\/\/\/\/\/\/\/\/\/\/\/\/*/

	void DefineBoundaries( void );			// defines the boundaries of each region using the classified segmented
											// image storing the resulting boundary locations for each region using
											// a region list object

	/*/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\*/
	/*  Image Data Searching/Distance Calculation */
	/*\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/*/

	//Usage: InWindow(modeIndex1, modeIndex2)
	bool InWindow(int, int);				//returns true if the range data of the specified data points
											//are within the defined search window (defined by kernel
											//bandwidth h[1])

	float SqDistance(int, int);				// computes the normalized square distance between two modes 

	/*/\/\/\/\/\/\/\/\/\/\*/
	/* Memory Management  */
	/*\/\/\/\/\/\/\/\/\/\/*/

	void InitializeOutput( void );			//Allocates memory needed by this class to perform image
											//filtering and segmentation

	void DestroyOutput( void );				//De-allocates memory needed by this class to perform image
											//filtering and segmentation

  //=============================
  // *** Private Data Members ***
  //=============================

   //##########################################
   //#######    IMAGE CLASSIFICATION   ########
   //##########################################

	/////////Image Boundaries/////////
	RegionList		*regionList;			// stores the boundary locations for each region

	/////////Image Regions////////
	int				regionCount;			// stores the number of connected regions contained by the
											// image

	/////////8 Connected Neighbors/////////
	int				neigh[8];

	/////////Index Table/////////////////
	int				*indexTable;			//used during fill algorithm

	/////////LUV_data/////////////////
   //int            *LUV_data;           //stores modes in integer format on lattice
	float				*LUV_data;				//stores modes in float format on lattice
   float          LUV_treshold;        //in float mode this determines what "close" means between modes


   //##########################################
   //#######   OUTPUT DATA STORAGE     ########
   //##########################################

	////////Raw Data (1 to 1 correlation with input)////////
	float			*msRawData;				// Raw data output of mean shift algorithm
											// to the location of the data point on the lattice

	////////Data Modes////////
	int				*labels;				// assigns a label to each data point associating it to
											// a mode in modes (e.g. a data point having label l has
											// mode modes[l])

	float			*modes;					// stores the mode data of the input data set, indexed by labels

	int				*modePointCounts;		// stores for each mode the number of point mapped to that mode,
											// indexed by labels

   //##########################################
   //#######  REGION ADJACENCY MATRIX  ########
   //##########################################

	//////////Region Adjacency List/////////
	RAList			*raList;				// an array of RAList objects containing an entry for each
											// region of the image

	//////////RAMatrix Free List///////////
	RAList			*freeRAList;			// a pointer to the head of a region adjacency list object
											// free list

	RAList			*raPool;				// a pool of RAList objects used in the construction of the
											// RAM

   //##############################################
   //#######  COMPUTATION OF EDGE STRENGTHS #######
   //##############################################

	//////////Epsilon///////////
	float			epsilon;				//Epsilon used for transitive closure

	//////Visit Table//////
	unsigned char	*visitTable;			//Table used to keep track of which pixels have been
											//already visited upon computing the boundary edge strengths

   //##########################################
   //#######       IMAGE PRUNING       ########
   //##########################################

	////////Transitive Closure/////////
	float			rR2;					//defines square range radius used when clustering pixels
											//together, thus defining image regions

   float speedThreshold; // the % of window radius used in new optimized filter 2.
};

#endif
