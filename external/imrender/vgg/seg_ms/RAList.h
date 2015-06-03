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

	The prototype for the RAList class is provided below. Its
	defition is provided in "RAList.cc".

The theory is described in the papers:

  D. Comaniciu, P. Meer: Mean Shift: A robust approach toward feature
									 space analysis.

  C. Christoudias, B. Georgescu, P. Meer: Synergism in low level vision.

and they are is available at:
  http://www.caip.rutgers.edu/riul/research/papers/

Implemented by Chris M. Christoudias, Bogdan Georgescu
********************************************************/

#ifndef RALIST_H
#define RALIST_H

//define Region Adjacency List class prototype
class RAList {

public:

	//============================
	// *** Public Data Members ***
	//============================

	////////////RAM Label//////////
	int		label;

	////////////RAM Weight/////////
	float	edgeStrength;
	int		edgePixelCount;

	////////////RAM Link///////////
	RAList	*next;

	//=======================
	// *** Public Methods ***
	//=======================

	/*/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\*/
	/* Class Constructor and Destructor */
	/*\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/*/

	//***Class Constrcutor***
	RAList( void );

	//***Class Destructor***
	~RAList( void );

	/*/\/\/\/\/\/\/\/\/\/\/\/\*/
	/*  RAM List Manipulation */
	/*\/\/\/\/\/\/\/\/\/\/\/\/*/

	//Usage: Insert(entry)
	int Insert(RAList*);		//Insert a region node into the region adjecency list

private:

	//=============================
	// *** Private Data Members ***
	//=============================

	///////current and previous pointer/////
	RAList	*cur, *prev;

	////////flag///////////
	unsigned char exists;

};

#endif
