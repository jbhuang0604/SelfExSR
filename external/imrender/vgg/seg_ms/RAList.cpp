/*******************************************************

                 Mean Shift Analysis Library
	=============================================

	The mean shift library is a collection of routines
	that use the mean shift algorithm. Using this algorithm,
	the necessary output will be generated needed
	to analyze a given input set of data.

  Region Adjacency List:
  =====================

	The Region Adjacency List class is used by the Image 
	Processor class in the construction of a Region Adjacency
	Matrix, used by	this class to applying transitive closure
	and to prune spurious regions during image segmentation.

	The definition of the RAList class is provided below. Its
	prototype is provided in "RAList.h".

The theory is described in the papers:

  D. Comaniciu, P. Meer: Mean Shift: A robust approach toward feature
									 space analysis.

  C. Christoudias, B. Georgescu, P. Meer: Synergism in low level vision.

and they are is available at:
  http://www.caip.rutgers.edu/riul/research/papers/

Implemented by Chris M. Christoudias, Bogdan Georgescu
********************************************************/
//include Region Adjacency List class prototype
#include	"RAList.h"

//include needed libraries
#include	<stdio.h>
#include	<assert.h>
#include	<stdlib.h>

/*@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@*/
/*@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@*/
/*@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@      PUBLIC METHODS     @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@*/
/*@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@*/
/*@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@*/

	/*/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\*/
	/* Class Constructor and Destructor */
	/*\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/*/

/*******************************************************/
/*Class Constructor                                    */
/*******************************************************/
/*Constructs a RAList object.                          */
/*******************************************************/
/*Post:                                                */
/*      - a RAlist object has been properly constru-   */
/*        cted.                                        */
/*******************************************************/

RAList::RAList( void )
{
	//initialize label and link
	label			= -1;
	next			= NULL;

	//initialize edge strenght weight and count
	edgeStrength	= 0;
	edgePixelCount	= 0;
}

/*******************************************************/
/*Class Destructor                                     */
/*******************************************************/
/*Destructrs a RAList object.                          */
/*******************************************************/
/*Post:                                                */
/*      - the RAList object has been properly dest-    */
/*        ructed.                                      */
/*******************************************************/

RAList::~RAList( void )
{
	//do nothing
}

/*******************************************************/
/*Insert                                               */
/*******************************************************/
/*Insert a region node into the region adjacency list. */
/*******************************************************/
/*Pre:                                                 */
/*      - entry is a node representing a connected re- */
/*        gion                                         */
/*Post:                                                */
/*      - entry has been inserted into the region adj- */
/*        acency list if it does not already exist     */
/*        there.                                       */
/*      - if the entry already exists in the region    */
/*        adjacency list 1 is returned otherwise 0 is  */
/*        returned.                                    */
/*******************************************************/

int RAList::Insert(RAList *entry)
{

	//if the list contains only one element
	//then insert this element into next
	if(!next)
	{
		//insert entry
		next		= entry;
		entry->next = NULL;

		//done
		return 0;
	}

	//traverse the list until either:

	//(a) entry's label already exists - do nothing
	//(b) the list ends or the current label is
	//    greater than entry's label, thus insert the entry
	//    at this location

	//check first entry
	if(next->label > entry->label)
	{
		//insert entry into the list at this location
		entry->next	= next;
		next		= entry;

		//done
		return 0;
	}

	//check the rest of the list...
	exists	= 0;
	cur		= next;
	while(cur)
	{
		if(entry->label == cur->label)
		{
			//node already exists
			exists = 1;
			break;
		}
		else if((!(cur->next))||(cur->next->label > entry->label))
		{
			//insert entry into the list at this location
			entry->next	= cur->next;
			cur->next	= entry;
			break;
		}

		//traverse the region adjacency list
		cur = cur->next;
	}

	//done. Return exists indicating whether or not a new node was
	//      actually inserted into the region adjacency list.
	return (int)(exists);

}

/*@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@*/
/*@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@*/
/*@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@ END OF CLASS DEFINITION @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@*/
/*@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@*/
/*@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@*/

