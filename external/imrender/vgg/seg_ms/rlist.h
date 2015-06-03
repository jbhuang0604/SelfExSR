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
	The prototype for the RegionList class is provided below. It
	is defined in "region.cc".

The theory is described in the papers:

  D. Comaniciu, P. Meer: Mean Shift: A robust approach toward feature
									 space analysis.

  C. Christoudias, B. Georgescu, P. Meer: Synergism in low level vision.

and they are is available at:
  http://www.caip.rutgers.edu/riul/research/papers/

Implemented by Chris M. Christoudias, Bogdan Georgescu
********************************************************/

#ifndef RLIST_H
#define RLIST_H

//include global type definitions
#include	"tdef.h"

//define region structure
struct REGION {
	int			label;
	int			pointCount;
	int			region;

};

//region class prototype...
class RegionList {

public:

  /*/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\*/
  /* Class Constructor and Destructor */
  /*\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/*/

  //--\\||//--\\||//--\\||//--\\||//--\\||//--\\||//--\\||//
  //<--------------------------------------------------->|//
  //|                                                    |//
  //|	Method Name:								     |//
  //|   ============								     |//
  //|               *  Class Constructor  *              |//
  //|                                                    |//
  //<--------------------------------------------------->|//
  //|                                                    |//
  //|	Description:								     |//
  //|	============								     |//
  //|                                                    |//
  //|   Constructs a region list object.                 |//
  //|                                                    |//
  //|   Its arguments are:                               |//
  //|                                                    |//
  //|   <* maxRegions *>                                 |//
  //|   The maximum amount of regions that can be class- |//
  //|   ified by the region list.                        |//
  //|                                                    |//
  //|   <* L *>                                          |//
  //|   The length of the input data set being class-    |//
  //|   ified by the region list object.                 |//
  //|                                                    |//
  //|   <* N *>                                          |//
  //|   The dimension of the input data set being class- |//
  //|   ified by the region list object.                 |//
  //|                                                    |//
  //<--------------------------------------------------->|//
  //|                                                    |//
  //|	Usage:      								     |//
  //|   ======      								     |//
  //|     RegionList(maxRegions, L, N)                   |//
  //|                                                    |//
  //<--------------------------------------------------->|//
  //--\\||//--\\||//--\\||//--\\||//--\\||//--\\||//--\\||//

	RegionList(int, int, int);

	// Class Destructor
	~RegionList( void );

  /*/\/\/\/\/\/\/\/\/\/\/\/\/\/\*/
  /*  Region List Manipulation  */
  /*\/\/\/\/\/\/\/\/\/\/\/\/\/\/*/

  //--\\||//--\\||//--\\||//--\\||//--\\||//--\\||//--\\||//
  //<--------------------------------------------------->|//
  //|                                                    |//
  //|	Method Name:								     |//
  //|   ============								     |//
  //|                *  Add Region  *                    |//
  //|                                                    |//
  //<--------------------------------------------------->|//
  //|                                                    |//
  //|	Description:								     |//
  //|	============								     |//
  //|                                                    |//
  //|   Adds a region to the region list.                |//
  //|                                                    |//
  //|   Its arguments are:                               |//
  //|                                                    |//
  //|   <* label *>                                      |//
  //|                                                    |//
  //|   A positive integer used to uniquely identify     |//
  //|   a region.                                        |//
  //|                                                    |//
  //|   <* pointCount *>                                 |//
  //|   A positive integer that specifies the number of  |//
  //|   N-dimensional data points that exist in the re-  |//
  //|   gion being classified.                           |//
  //|                                                    |//
  //|   <* indeces *>                                    |//
  //|   An integer array that specifies the set of ind-  |//
  //|   eces of the data points that are contianed with- |//
  //|   in this region.                                  |//
  //|                                                    |//
  //<--------------------------------------------------->|//
  //|                                                    |//
  //|	Usage:      								     |//
  //|   ======      								     |//
  //|     AddRegion(label, pointCount, indeces)          |//
  //|                                                    |//
  //<--------------------------------------------------->|//
  //--\\||//--\\||//--\\||//--\\||//--\\||//--\\||//--\\||//

	void AddRegion(int, int, int*);

  //--\\||//--\\||//--\\||//--\\||//--\\||//--\\||//--\\||//
  //<--------------------------------------------------->|//
  //|                                                    |//
  //|	Method Name:								     |//
  //|   ============								     |//
  //|                    *  Reset  *                     |//
  //|                                                    |//
  //<--------------------------------------------------->|//
  //|                                                    |//
  //|	Description:								     |//
  //|	============								     |//
  //|                                                    |//
  //|   Resets the region list for re-use (for new       |//
  //|   classification).                                 |//
  //|                                                    |//
  //<--------------------------------------------------->|//
  //--\\||//--\\||//--\\||//--\\||//--\\||//--\\||//--\\||//

	void Reset( void );	

  /*/\/\/\/\/\/\/\/\/\/\*/
  /*  Query Region List */
  /*\/\/\/\/\/\/\/\/\/\/*/

  //--\\||//--\\||//--\\||//--\\||//--\\||//--\\||//--\\||//
  //<--------------------------------------------------->|//
  //|                                                    |//
  //|	Method Name:								     |//
  //|   ============								     |//
  //|          *  Get Number of Regions  *               |//
  //|                                                    |//
  //<--------------------------------------------------->|//
  //|                                                    |//
  //|	Description:								     |//
  //|	============								     |//
  //|                                                    |//
  //|   Returns the number of regions stored by the      |//
  //|   region list.                                     |//
  //|                                                    |//
  //<--------------------------------------------------->|//
  //--\\||//--\\||//--\\||//--\\||//--\\||//--\\||//--\\||//

	int	GetNumRegions ( void );

  //--\\||//--\\||//--\\||//--\\||//--\\||//--\\||//--\\||//
  //<--------------------------------------------------->|//
  //|                                                    |//
  //|	Method Name:								     |//
  //|   ============								     |//
  //|                  *  Get Label  *                   |//
  //|                                                    |//
  //<--------------------------------------------------->|//
  //|                                                    |//
  //|	Description:								     |//
  //|	============								     |//
  //|                                                    |//
  //|   Returns the label of a specified region.         |//
  //|                                                    |//
  //|   Its arguments are:                               |//
  //|                                                    |//
  //|   <* regionNumber *>                               |//
  //|   The index of the region in the region list       |//
  //|   array.                                           |//
  //|                                                    |//
  //<--------------------------------------------------->|//
  //|                                                    |//
  //|	Usage:      								     |//
  //|   ======      								     |//
  //|     label = GetLabel(regionNumber)                 |//
  //|                                                    |//
  //<--------------------------------------------------->|//
  //--\\||//--\\||//--\\||//--\\||//--\\||//--\\||//--\\||//

	int	GetLabel(int);

  //--\\||//--\\||//--\\||//--\\||//--\\||//--\\||//--\\||//
  //<--------------------------------------------------->|//
  //|                                                    |//
  //|	Method Name:								     |//
  //|   ============								     |//
  //|                *  Get Region Count  *              |//
  //|                                                    |//
  //<--------------------------------------------------->|//
  //|                                                    |//
  //|	Description:								     |//
  //|	============								     |//
  //|                                                    |//
  //|   Returns number of data points contained by a sp- |//
  //|   ecified region.                                  |//
  //|                                                    |//
  //|   Its arguments are:                               |//
  //|                                                    |//
  //|   <* regionNumber *>                               |//
  //|   The index of the region in the region list       |//
  //|   array.                                           |//
  //|                                                    |//
  //<--------------------------------------------------->|//
  //|                                                    |//
  //|	Usage:      								     |//
  //|   ======      								     |//
  //|     pointCount = GetRegionCount(regionNumber)      |//
  //|                                                    |//
  //<--------------------------------------------------->|//
  //--\\||//--\\||//--\\||//--\\||//--\\||//--\\||//--\\||//

	int GetRegionCount(int);

  //--\\||//--\\||//--\\||//--\\||//--\\||//--\\||//--\\||//
  //<--------------------------------------------------->|//
  //|                                                    |//
  //|	Method Name:								     |//
  //|   ============								     |//
  //|               *  Get Region Indeces  *             |//
  //|                                                    |//
  //<--------------------------------------------------->|//
  //|                                                    |//
  //|	Description:								     |//
  //|	============								     |//
  //|                                                    |//
  //|   Returns a pointer to a set of grid location ind- |//
  //|   eces specifying the data points belonging to a   |//
  //|   specified region.                                |//
  //|                                                    |//
  //|   Its arguments are:                               |//
  //|                                                    |//
  //|   <* regionNumber *>                               |//
  //|   The index of the region in the region list       |//
  //|   array.                                           |//
  //|                                                    |//
  //<--------------------------------------------------->|//
  //|                                                    |//
  //|	Usage:      								     |//
  //|   ======      								     |//
  //|     indeces = GetRegionIndeces(regionNumber)       |//
  //|                                                    |//
  //<--------------------------------------------------->|//
  //--\\||//--\\||//--\\||//--\\||//--\\||//--\\||//--\\||//

	int*GetRegionIndeces(int);

private:

  /*/\/\/\/\/\/\/\/\/\/\/\*/
  /*  Class Error Handler */
  /*\/\/\/\/\/\/\/\/\/\/\/*/

	void ErrorHandler(char*, char*, ErrorType);

  //=============================
  // *** Private Data Members ***
  //=============================

	//#####################################
	//### REGION LIST PARTITIONED ARRAY ###
	//#####################################

	REGION		*regionList;			//array of maxRegions regions
	int			minRegion;

	int			maxRegions;				//defines the number maximum number of regions
										//allowed (determined by user during class construction)
	int			numRegions;				//the number of regions currently stored by the
										//region list
	int			freeRegion;				//an index into the regionList pointing to the next
										//available region in the regionList

	//#####################################
	//###         INDEX TABLE           ###
	//#####################################

	int			*indexTable;			//an array of indexes that point into an external structure
										//specifying which points belong to a region
	int			freeBlockLoc;			//points to the next free block of memory in the indexTable

	//#####################################
	//###     INPUT DATA PARAMETERS     ###
	//#####################################

	//Dimension of data set
	int			N;						//dimension of data set being classified by region list
										//class

	//Length of the data set
	int			L;						//number of points contained by the data set being classified by
										//region list class

};

#endif



