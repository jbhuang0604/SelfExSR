/*******************************************************

                 Mean Shift Analysis Library
	=============================================


	The mean shift library is a collection of routines
	that use the mean shift algorithm. Using this algorithm,
	the necessary output will be generated needed
	to analyze a given input set of data.

  Region List Class:
  =================

	During segmentation, data regions are defined. The 
	RegionList class provides a mechanism for doing so, as
	well as defines some basic operations, such as region
	growing or small region pruning, on the defined regions.
	It is defined below. Its prototype is given in "region.h".

The theory is described in the papers:

  D. Comaniciu, P. Meer: Mean Shift: A robust approach toward feature
									 space analysis.

  C. Christoudias, B. Georgescu, P. Meer: Synergism in low level vision.

and they are is available at:
  http://www.caip.rutgers.edu/riul/research/papers/

Implemented by Chris M. Christoudias, Bogdan Georgescu
********************************************************/


#include	"rlist.h"
#include	<stdio.h>
#include	<stdlib.h>

/*@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@*/
/*@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@*/
/*@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@      PUBLIC METHODS     @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@*/
/*@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@*/
/*@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@*/

  /*/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\*/
  /*** Class Constructor and Destructor ***/
  /*\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/*/

/*******************************************************/
/*Constructor                                          */
/*******************************************************/
/*Constructor                                          */
/*******************************************************/
/*Pre:                                                 */
/*      - modesPtr is a pointer to an array of modes   */
/*      - maxRegions_ is the maximum number of regions */
/*        that can be defined                          */
/*      - L_ is the number of data points being class- */
/*        ified by the region list class               */
/*      - N is the dimension of the data set being cl- */
/*        assified by the region list class            */
/*Post:                                                */
/*      - a region list object has been properly init- */
/*        ialized.                                     */
/*******************************************************/

RegionList::RegionList(int maxRegions_, int L_, int N_)
{

	//Obtain maximum number of regions that can be
	//defined by user
	if((maxRegions = maxRegions_) <= 0)
		ErrorHandler("RegionList", "Maximum number of regions is zero or negative.", FATAL);

	//Obtain dimension of data set being classified by
	//region list class
	if((N = N_) <= 0)
		ErrorHandler("RegionList", "Dimension is zero or negative.", FATAL);

	//Obtain length of input data set...
	if((L = L_) <= 0)
		ErrorHandler("RegionList", "Length of data set is zero or negative.", FATAL);

	//Allocate memory for index table
	if(!(indexTable = new int [L]))
		ErrorHandler("RegionList", "Not enough memory.", FATAL);

	//Allocate memory for region list array
	if(!(regionList = new REGION [maxRegions]))
		ErrorHandler("RegionList", "Not enough memory.", FATAL);

	//Initialize region list...
	numRegions		= freeRegion = 0;

	//Initialize indexTable
	freeBlockLoc	= 0;

	//done.
	return;

}

/*******************************************************/
/*Destructor                                           */
/*******************************************************/
/*Destroys region list object.                         */
/*******************************************************/
/*Post:                                                */
/*      - region list object has been properly dest-   */
/*        oyed.                                        */
/*******************************************************/

RegionList::~RegionList( void )
{
	//de-allocate memory...
	delete [] regionList;
	delete [] indexTable;

	//done.
	return;
}

  /*/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\*/
  /***  Region List Manipulation  ***/
  /*\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/*/

/*******************************************************/
/*Add Region                                           */
/*******************************************************/
/*Adds a region to the region list.                    */
/*******************************************************/
/*Pre:                                                 */
/*      - label is a positive integer used to uniquely */
/*        identify a region                            */
/*      - pointCount is the number of N-dimensional    */
/*        data points that exist in the region being   */
/*        classified.                                  */
/*      - indeces is a set of indeces specifying the   */
/*        data points contained within this region     */
/*      - pointCount must be > 0                       */
/*Post:                                                */
/*      - a new region labeled using label and contai- */
/*        ning pointCount number of points has been    */
/*        added to the region list.                    */
/*******************************************************/

void RegionList::AddRegion(int label, int pointCount, int *indeces)
{

	//make sure that there is enough room for this new region 
	//in the region list array...
	if(numRegions >= maxRegions)
		ErrorHandler("AddRegion", "Not enough memory allocated.", FATAL);

	//make sure that label is positive and point Count > 0...
	if((label < 0)||(pointCount <= 0))
		ErrorHandler("AddRegion", "Label is negative or number of points in region is invalid.", FATAL);

	//make sure that there is enough memory in the indexTable
	//for this region...
	if((freeBlockLoc + pointCount) > L)
		ErrorHandler("AddRegion", "Adding more points than what is contained in data set.", FATAL);

	//place new region into region list array using
	//freeRegion index
	regionList[freeRegion].label		= label;
	regionList[freeRegion].pointCount	= pointCount;
	regionList[freeRegion].region		= freeBlockLoc;

	//copy indeces into indexTable using freeBlock...
	int i;
	for(i = 0; i < pointCount; i++)
		indexTable[freeBlockLoc+i] = indeces[i];

	//increment freeBlock to point to the next free
	//block
	freeBlockLoc	+= pointCount;

	//increment freeRegion to point to the next free region
	//also, increment numRegions to indicate that another
	//region has been added to the region list
	freeRegion++;
	numRegions++;

	//done.
	return;

}

/*******************************************************/
/*Reset                                                */
/*******************************************************/
/*Resets the region list.                              */
/*******************************************************/
/*Post:                                                */
/*      - the region list has been reset.              */
/*******************************************************/

void RegionList::Reset( void )
{

	//reset region list
	freeRegion = numRegions = freeBlockLoc = 0;

	//done.
	return;

}

  /*/\/\/\/\/\/\/\/\/\/\*/
  /*  Query Region List */
  /*\/\/\/\/\/\/\/\/\/\/*/

/*******************************************************/
/*Get Number Regions                                   */
/*******************************************************/
/*Returns the number of regions stored by region list. */
/*******************************************************/
/*Post:                                                */
/*      - the number of regions stored by the region   */
/*        list is returned.                            */
/*******************************************************/

int	RegionList::GetNumRegions( void )
{
	// return region count
	return numRegions;
}

/*******************************************************/
/*Get Label                                            */
/*******************************************************/
/*Returns the label of a specified region.             */
/*******************************************************/
/*Pre:                                                 */
/*      - regionNum is an index into the region list   */
/*        array.                                       */
/*Post:                                                */
/*      - the label of the region having region index  */
/*        specified by regionNum has been returned.    */
/*******************************************************/

int	RegionList::GetLabel(int regionNum)
{
	//return the label of a specified region
	return regionList[regionNum].label;
}

/*******************************************************/
/*Get Region Count                                     */
/*******************************************************/
/*Returns the point count of a specified region.       */
/*******************************************************/
/*Pre:                                                 */
/*      - regionNum is an index into the region list   */
/*        array.                                       */
/*Post:                                                */
/*      - the number of points that classify the       */
/*        region whose index is specified by regionNum */
/*        is returned.                                 */
/*******************************************************/

int RegionList::GetRegionCount(int regionNum)
{
	//return the region count of a specified region
	return regionList[regionNum].pointCount;
}

/*******************************************************/
/*Get Region Indeces                                   */
/*******************************************************/
/*Returns the point indeces specifying a region.       */
/*******************************************************/
/*Pre:                                                 */
/*      - regionNum is an index into the region list   */
/*        array.                                       */
/*Post:                                                */
/*      - the region indeces specifying the points     */
/*        contained by the region specified by region- */
/*        Num are returned.                            */
/*******************************************************/

int *RegionList::GetRegionIndeces(int regionNum)
{
	//return point indeces using regionNum
	return &indexTable[regionList[regionNum].region];
}

/*@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@*/
/*@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@*/
/*@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@     PRIVATE METHODS     @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@*/
/*@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@*/
/*@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@*/

  /*/\/\/\/\/\/\/\/\/\/\/\*/
  /*  Class Error Handler */
  /*\/\/\/\/\/\/\/\/\/\/\/*/

/*******************************************************/
/*Error Handler                                        */
/*******************************************************/
/*Class error handler.                                 */
/*******************************************************/
/*Pre:                                                 */
/*      - functName is the name of the function that   */
/*        caused an error                              */
/*      - errmsg is the error message given by the     */
/*        calling function                             */
/*      - status is the error status: FATAL or NON-    */
/*        FATAL                                        */
/*Post:                                                */
/*      - the error message errmsg is flagged on beh-  */
/*        ave of function functName.                   */
/*      - if the error status is FATAL then the program*/
/*        is halted, otherwise execution is continued, */
/*        error recovery is assumed to be handled by   */
/*        the calling function.                        */
/*******************************************************/

void RegionList::ErrorHandler(char *functName, char* errmsg, ErrorType status)
{

	//flag error message on behalf of calling function, error format
	//specified by the error status...
	if(status == NONFATAL)
		fprintf(stderr, "\n%s Error: %s\n", functName, errmsg);
	else
	{
		fprintf(stderr, "\n%s Fatal Error: %s\n\nAborting Program.\n\n", functName, errmsg);
		exit(1);
	}

}

/*@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@*/
/*@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@*/
/*@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@ END OF CLASS DEFINITION @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@*/
/*@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@*/
/*@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@*/
