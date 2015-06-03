/*******************************************************

  Mean Shift Analysis Library
  =============================================
  
	The mean shift library is a collection of routines
	that use the mean shift algorithm. Using this algorithm,
	the necessary output will be generated needed
	to analyze a given input set of data.
	
  MeanShift Base Class:
  ====================
	  
	The mean shift library of routines is realized
	via the creation of a MeanShift base class. This class
	provides a mechanism for calculating the mean shift vector
	at a specified data point, using an arbitrary N-dimensional
	data set, and a user-defined kernel.
		
	For image processing the mean shift base class also allows
	for the definition of a data set that is on a two-dimensional
	lattice. The amount of time needed to compute the mean shift
	vector using such a data set is much less than that of an
	arbitrary one. Because images usually contain many data points,
	defining the image input data points as being on a lattice
	greatly improves computation time and makes algorithms such
	as image filtering practical.
		  
	The definition of the MeanShift class is provided below. Its
	prototype is provided in 'ms.h'.
			
The theory is described in the papers:

  D. Comaniciu, P. Meer: Mean Shift: A robust approach toward feature
									 space analysis.

  C. Christoudias, B. Georgescu, P. Meer: Synergism in low level vision.

and they are is available at:
  http://www.caip.rutgers.edu/riul/research/papers/

Implemented by Chris M. Christoudias, Bogdan Georgescu
********************************************************/


//Include Needed Libraries

#include	"ms.h"
#include	<string.h>
#include	<stdlib.h>
#include	<stdio.h>
#include	<math.h>

/*@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@*/
/*@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@*/
/*@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@      PUBLIC METHODS     @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@*/
/*@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@*/
/*@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@*/

  /*/\/\/\/\/\/\/\/\/\/\/\/\/\/\*/
  /*** Constructor/Destructor ***/
  /*\/\/\/\/\/\/\/\/\/\/\/\/\/\/*/

/*******************************************************/
/*Class Constructor                                    */
/*******************************************************/
/*Post:                                                */
/*      The MeanShift class has been properly          */
/*      initialized.                                   */
/*******************************************************/

MeanShift::MeanShift( void )
{
	
	//intialize input data set parameters...
	P							= NULL;
	L							= 0;
	N							= 0;
	kp							= 0;
	
	//initialize input data set storage structures...
	data						= NULL;
	
	//initialize input data set kd-tree
	root						= NULL;
	forest						= NULL;
	range						= NULL;
	
	//intialize lattice structure...
	height						= 0;
	width						= 0;
	
	//intialize kernel strucuture...
	h							= NULL;
	kernel						= NULL;
	w							= NULL;
	offset						= NULL;
	increment					= NULL;
	uniformKernel				= false;
	
	//initialize weight function linked list...
	head						= cur	= NULL;
	
	//intialize mean shift processing data structures...
	uv							= NULL;

	//set lattice weight map to null
	weightMap					= NULL;

	//indicate that the lattice weight map is undefined
	weightMapDefined			= false;
	
	//allocate memory for error message buffer...
	ErrorMessage				= new char [256];
	
	//initialize error status to OKAY
	ErrorStatus					= EL_OKAY;
	
	//Initialize class state...
	class_state.INPUT_DEFINED	= false;
	class_state.KERNEL_DEFINED	= false;
	class_state.LATTICE_DEFINED	= false;
	class_state.OUTPUT_DEFINED	= false;
	
}

/*******************************************************/
/*Class Destructor                                     */
/*******************************************************/
/*Post:                                                */
/*      The MeanShift class has been properly          */
/*      destroyed.                                     */
/*******************************************************/

MeanShift::~MeanShift( void )
{
	delete [] ErrorMessage;
   if (weightMap)
   {
      delete [] weightMap;
   }

	//de-allocate memory used to store
	//user defined weight functions
	ClearWeightFunctions();
	
	//de-allocate memory used for kernel
	DestroyKernel();
	
	//de-allocate memory used for input
	ResetInput();
	
}

  /*/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\*/
  /*** Creation/Initialization of Mean Shift Kernel ***/
  /*\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/*/

/*******************************************************/
/*Define Kernel                                        */
/*******************************************************/
/*Creats custom user defined Kernel to be used by the  */
/*mean shift procedure.                                */
/*******************************************************/
/*Pre:                                                 */
/*      - kernel is an array of kernelTypes specifying */
/*        the type of kernel to be used on each sub-   */
/*        space of the input data set x                */
/*      - h is the set of bandwidths used to define the*/
/*        the search window                            */
/*      - P is a one dimensional array of integers of  */
/*        size kp, that specifies the dimension of each*/
/*        subspace of the input data set x             */
/*      - kp is the total number of subspaces used to  */
/*        the input data set x                         */
/*Post:                                                */
/*      - the custom kernel has been created for use   */
/*        by the mean shift procedure.                 */
/*******************************************************/

void MeanShift::DefineKernel(kernelType *kernel_, float *h_, int *P_, int kp_)
{
	
	// Declare variables
	int i, kN;
	
	//if a kernel has already been created then destroy it
	if(kp)
		DestroyKernel();
	
	//Obtain kp...
	if((kp = kp_) <= 0)
	{
		ErrorHandler("MeanShift", "CreateKernel", "Subspace count (kp) is zero or negative.");
		return;
	}
	
	//Allocate memory for h, P, kernel, offset, and increment
	if((!(P = new int [kp]))||(!(h = new float [kp]))||(!(kernel = new kernelType [kp]))||
		(!(offset = new float [kp]))||(!(increment = new double [kp])))
	{
		ErrorHandler("MeanShift", "CreateKernel", "Not enough memory available to create kernel.");
		return;
	}
	
	//Populate h, P and kernel, also use P to calculate
	//the dimension (N_) of the potential input data set x
	kN = 0;
	for(i = 0; i < kp; i++)
	{
		if((h[i] = h_[i]) <= 0)
		{
			ErrorHandler("MeanShift", "CreateKernel", "Negative or zero valued bandwidths are prohibited.");
			return;
		}
		if((P[i] = P_[i]) <= 0)
		{
			ErrorHandler("MeanShift", "CreateKernel", "Negative or zero valued subspace dimensions are prohibited.");
			return;
		}
		kernel[i] = kernel_[i];
		kN	   += P[i];
	}
	
	//Allocate memory for range vector and uv using N_
	if((!(range = new float [2*kN]))||(!(uv = new double [kN])))
	{
		ErrorHandler("MeanShift", "CreateKernel", "Not enough memory available to create kernel.");
		return;
	}
	
	// Generate weight function lookup table
	// using above information and user
	// defined weight function list
	generateLookupTable();
	
	//check for errors
	if(ErrorStatus == EL_ERROR)
		return;
	
	//indicate that the kernel has been defined
	class_state.KERNEL_DEFINED	= true;
	
	//done.
	return;
	
}

/*******************************************************/
/*Add Weight Function                                  */
/*******************************************************/
/*Adds a weight function to the Mean Shift class to be */
/*used by the mean shift procedure                     */
/*******************************************************/
/*Pre:                                                 */
/*      - g(u) is the normalized weight function with  */
/*        respect to u = (norm(x-xi))^2/h^2            */
/*      - sampleNumber is the number of samples to be  */
/*        taken of g(u) over halfWindow interval       */
/*      - halfWindow is the radius of g(u) such that   */
/*        g(u) is defined for 0 <= u <= halfWindow     */
/*      - subspace is the subspace number for which    */
/*        g(u) is to be applied during the mean shift  */
/*        procedure.                                   */
/*Post:                                                */
/*      - g(u) has been added to the Mean Shift class  */
/*        private data structure to be used by the     */
/*        mean shift procedure.                        */
/*      - if a weight function has already been spec-  */
/*        ified for the specified subspace, the weight */
/*        function for this subspace has been replaced.*/
/*******************************************************/

void MeanShift::AddWeightFunction(double g(double), float halfWindow, int sampleNumber, int subspace)
{
	
	// Declare Variables
	int   i;
	double increment;
	
	// Search to see if a weight function has already been
	// defined for specified subspace, if not then insert
	// into the head of the weight function list, otherwise
	// replace entry
	
	// Perform Search
	cur = head;
	while((cur)&&(cur->subspace != subspace))
		cur = cur->next;
	
	// Entry Exists - Replace It! 
	// Otherwise insert at the head of the the weight functon list
	if(cur)
		delete cur->w;
	else
    {
		cur       = new userWeightFunct;
		cur->next = head;
		head      = cur;
    }
	
	// Generate lookup table
	increment = halfWindow/(double)(sampleNumber);
	
	cur->w = new double [sampleNumber+1];
	for(i = 0; i <= sampleNumber; i++)
		cur->w[i] = g((double)(i*increment));
	
	// Set weight function parameters
	cur->halfWindow   = halfWindow;
	cur->sampleNumber = sampleNumber;
	cur->subspace     = subspace;
	
	//done.
	return;
	
}

/*******************************************************/
/*Clear Weight Functions                               */
/*******************************************************/
/*Clears user defined weight from the Mean Shift class */
/*private data structure.                              */
/*******************************************************/
/*Post:                                                */
/*      - all user defined weight functions ahve been  */
/*        cleared from the private data structure of   */
/*        the mean shift class.                        */
/*******************************************************/

void MeanShift::ClearWeightFunctions( void )
{
	
	while(head)
    {
		delete head->w;
		cur  = head;
		head = head->next;
		delete cur;
    }
	
}

  /*/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\*/
  /*** Input Data Set Declaration ***/
  /*\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/*/

/*******************************************************/
/*Define Input                                         */
/*******************************************************/
/*Uploads input data set x into the mean shift class.  */
/*******************************************************/
/*Pre:                                                 */
/*      - x is a one dimensional array of L N-dimen-   */
/*        ional data points.                           */
/*Post:                                                */
/*      - x has been uploaded into the mean shift      */
/*        class.                                       */
/*      - the height and width of a previous data set  */
/*        has been undefined.                          */
/*******************************************************/

void MeanShift::DefineInput(float *x, int L_, int N_)
{
	
	
	//if input data is defined de-allocate memory, and
	//re-initialize the input data structure
	if((class_state.INPUT_DEFINED)||(class_state.LATTICE_DEFINED))
		ResetInput();
	
	//make sure x is not NULL...
	if(!x)
	{
		ErrorHandler("MeanShift", "UploadInput", "Input data set is NULL.");
		return;
	}
	
	//Obtain L and N
	if(((L = L_) <= 0)||((N = N_) <= 0))
	{
		ErrorHandler("MeanShift", "UploadInput", "Input data set has negative or zero length or dimension.");
		return;
	}
	
	//Allocate memory for data
	if(!(data = new float [L*N]))
	{
		ErrorHandler("MeanShift", "UploadInput", "Not enough memory.");
		return;
	}
	
	//Allocate memory for input data set, and copy
	//x into the private data members of the mean
	//shift class
	InitializeInput(x);
	
	//check for errors
	if(ErrorStatus == EL_ERROR)
		return;
	
	// Load x into the MeanShift object using
	// using a kd-tree, resulting in better
	// range searching of the input data points
	// x - also upload window centers into
	// msRawData
	CreateBST();
	
	//indicate that the input has been recently defined
	class_state.INPUT_DEFINED	= true;
	class_state.LATTICE_DEFINED	= false;
	class_state.OUTPUT_DEFINED	= false;
	
	//done.
	return;
	
}

/*******************************************************/
/*Define Lattice                                       */
/*******************************************************/
/*Defines the height and width of the input lattice.   */
/*******************************************************/
/*Pre:                                                 */
/*      - ht is the height of the lattice              */
/*      - wt is the width of the lattice               */
/*Post:                                                */
/*      - the height and width of the lattice has been */
/*        specified.                                   */
/*      - if a data set is presently loaded into the   */
/*        mean shift class, an error is flagged if the */
/*        number of elements in that data set does not */
/*        equal the product ht*wt.                     */
/*******************************************************/

void MeanShift::DefineLInput(float *x, int ht, int wt, int N_)
{
	
	//if input data is defined de-allocate memory, and
	//re-initialize the input data structure
	if((class_state.INPUT_DEFINED)||(class_state.LATTICE_DEFINED))
		ResetInput();
	
	//Obtain lattice height and width
	if(((height	= ht) <= 0)||((width	= wt) <= 0))
	{
		ErrorHandler("MeanShift", "DefineLInput", "Lattice defined using zero or negative height and/or width.");
		return;
	}
	
	//Obtain input data dimension
	if((N = N_) <= 0)
	{
		ErrorHandler("MeanShift", "DefineInput", "Input defined using zero or negative dimension.");
		return;
	}
	
	//compute the data length, L, of input data set
	//using height and width
	L		= height*width;
	
	//Allocate memory for input data set, and copy
	//x into the private data members of the mean
	//shift class
	InitializeInput(x);
	
	//check for errors
	if(ErrorStatus == EL_ERROR)
		return;

	//allocate memory for weight map
	if(!(weightMap = new float [L]))
	{
		ErrorHandler("MeanShift", "InitializeInput", "Not enough memory.");
		return;
	}

	//initialize weightMap to an array of zeros
	memset(weightMap, 0, L*(sizeof(float)));
	
	//Indicate that a lattice input has recently been
	//defined
	class_state.LATTICE_DEFINED	= true;
	class_state.INPUT_DEFINED	= false;
	class_state.OUTPUT_DEFINED	= false;
	
	//done.
	return;
	
}

/*******************************************************/
/*Set Lattice Weight Map                               */
/*******************************************************/
/*Populates the lattice weight map with specified      */
/*weight values.                                       */
/*******************************************************/
/*Pre:                                                 */
/*      - wm is a floating point array of size L       */
/*        specifying for each data point a weight      */
/*        value                                        */
/*Post:                                                */
/*      - wm has been used to populate the lattice     */
/*        weight map.                                  */
/*******************************************************/

void MeanShift::SetLatticeWeightMap(float *wm)
{
	//make sure wm is not NULL
	if(!wm)
	{
		ErrorHandler("MeanShift", "SetWeightMap", "Specified weight map is NULL.");
		return;
	}

	//populate weightMap using wm
	int i;
	for(i = 0; i < L; i++)
		weightMap[i] = wm[i];

	//indicate that a lattice weight map has been specified
	weightMapDefined	= true;

	//done.
	return;

}


/*******************************************************/
/*Remove Lattice Weight Map                            */
/*******************************************************/
/*Removes the lattice weight map.                      */
/*******************************************************/
/*Post:                                                */
/*      - the lattice weight map has been removed.     */
/*      - if a weight map did not exist NO error is    */
/*        flagged.                                     */
/*******************************************************/

void MeanShift::RemoveLatticeWeightMap(void)
{

	//only remove weight map if it exists, otherwise
	//do nothing...
	if(weightMapDefined)
	{
		//set values of lattice weight map to zero
		memset(weightMap, 0, L*sizeof(float));

		//indicate that a lattice weight map is no longer
		//defined
		weightMapDefined	= false;
	}

	//done.
	return;

}


  /*/\/\/\/\/\/\/\/\/\/\/\/\/\/\*/
  /*** Mean Shift Operations  ***/
  /*\/\/\/\/\/\/\/\/\/\/\/\/\/\/*/

/*******************************************************/
/*Mean Shift Vector                                    */
/*******************************************************/
/*Calculates the mean shift vector at a specified data */
/*point yk.                                            */
/*******************************************************/
/*Pre:                                                 */
/*      - a kernel has been created                    */
/*      - a data set has been uploaded                 */
/*      - Mh is an N dimensional mean shift vector     */
/*      - yk is an N dimensional data point            */
/*Post:                                                */
/*      - the mean shift vector at yk has been         */
/*        calculated and stored in and returned by Mh. */
/*******************************************************/

void MeanShift::msVector(double *Mh, double *yk)
{
	
	//make sure that Mh and/or yk are not NULL...
	if((!Mh)||(!yk))
	{
		ErrorHandler("MeanShift", "msVector", "Invalid argument(s) passed to this method.");
		return;
	}
	
	//make sure that a kernel has been created, data has
	//been uploaded, and that they are consistent with one
	//another...
	classConsistencyCheck(N, false);
	
	//calculate mean shift vector at yk using created kernel
	//and uploaded data set
	MSVector(Mh, yk);
	
	//done.
	return;
	
}

/*******************************************************/
/*Lattice Mean Shift Vector                            */
/*******************************************************/
/*Calculates the mean shift vector at a specified data */
/*point yk, assuming that the data set exhists on a    */
/*height x width two dimensional lattice.              */
/*******************************************************/
/*Pre:                                                 */
/*      - a kernel has been created                    */
/*      - a data set has been uploaded                 */
/*      - the height and width of the lattice has been */
/*        specified using method DefineLattice()       */
/*      - Mh is an N dimensional mean shift vector     */
/*      - yk is an N dimensional data point            */
/*Post:                                                */
/*      - the mean shift vector at yk has been         */
/*        calculated and stored in and returned by Mh. */
/*      - Mh was calculated using the defined input    */
/*        lattice.                                     */
/*******************************************************/

void MeanShift::latticeMSVector(double *Mh, double *yk)
{
	
	//make sure that Mh and/or yk are not NULL...
	if((!Mh)||(!yk))
	{
		ErrorHandler("MeanShift", "lmsVector", "Invalid argument(s) passed to this method.");
		return;
	}
	
	//make sure that a kernel has been created, data has
	//been uploaded, and that they are consistent with one
	//another...
	classConsistencyCheck(N+2, true);
	
	//calculate mean shift vector at yk using created kernel
	//and uploaded data set
	LatticeMSVector(Mh, yk);
	
	//done.
	return;
	
}

/*******************************************************/
/*Find Mode                                            */
/*******************************************************/
/*Calculates the mode of a specified data point yk.    */
/*******************************************************/
/*Pre:                                                 */
/*      - a kernel has been created                    */
/*      - a data set has been uploaded                 */
/*      - mode is the N dimensional mode of the N-dim- */
/*        ensional data point yk                       */
/*Post:                                                */
/*      - the mode of yk has been calculated and       */
/*        stored in mode.                              */
/*******************************************************/

void MeanShift::FindMode(double *mode, double *yk)
{
	
	//make sure that mode and/or yk are not NULL...
	if((!mode)||(!yk))
	{
		ErrorHandler("MeanShift", "FindMode", "Invalid argument(s) passed to this method.");
		return;
	}
	
	//make sure that a kernel has been created, data has
	//been uploaded, and that they are consistent with one
	//another...
	classConsistencyCheck(N, false);
	
	//allocate memory for Mh
	double	*Mh	= new double [N];
	
	//copy yk into mode
	int i;
	for(i = 0; i < N; i++)
		mode[i] = yk[i];
	
	//calculate mean shift vector at yk
	MSVector(Mh, yk);
	
	//calculate mvAbs = |Mh|^2
	double mvAbs = 0;
	for(i = 0; i < N; i++)
		mvAbs	+= Mh[i]*Mh[i];
	
	//shift mode until convergence (mvAbs = 0)...
	int iterationCount = 1;
	while((mvAbs >= EPSILON)&&(iterationCount < LIMIT))
	{
		//shift mode...
		for(i = 0; i < N; i++)
			mode[i]	+= Mh[i];
		
		//re-calculate mean shift vector at new
		//window location have center defined by
		//mode
		MSVector(Mh, mode);
		
		//calculate mvAbs = |Mh|^2
		mvAbs = 0;
		for(i = 0; i < N; i++)
			mvAbs	+= Mh[i]*Mh[i];
		
		//increment interation count...
		iterationCount++;
		
	}
	
	//shift mode...
	for(i = 0; i < N; i++)
		mode[i]	+= Mh[i];
	
	//de-allocate memory
	delete [] Mh;
	
	//done.
	return;
	
}

/*******************************************************/
/*Find Lattice Mode                                    */
/*******************************************************/
/*Calculates the mode of a specified data point yk,    */
/*assuming that the data set exhists on a height x     */
/*width two dimensional lattice.                       */
/*******************************************************/
/*Pre:                                                 */
/*      - a kernel has been created                    */
/*      - a data set has been uploaded                 */
/*      - the height and width of the lattice has been */
/*        specified using method DefineLattice()       */
/*      - mode is the N dimensional mode of the N-dim- */
/*        ensional data point yk                       */
/*Post:                                                */
/*      - the mode of yk has been calculated and       */
/*        stored in mode.                              */
/*      - mode was calculated using the defined input  */
/*        lattice.                                     */
/*******************************************************/

void MeanShift::FindLMode(double *mode, double *yk)
{
	
	//make sure that mode and/or yk are not NULL...
	if((!mode)||(!yk))
	{
		ErrorHandler("MeanShift", "FindLMode", "Invalid argument(s) passed to this method.");
		return;
	}
	
	//make sure the lattice height and width have been defined...
	if(!height)
	{
		ErrorHandler("MeanShift", "FindLMode", "Lattice height and width is undefined.");
		return;
	}
	
	//make sure that a kernel has been created, data has
	//been uploaded, and that they are consistent with one
	//another...
	classConsistencyCheck(N+2, true);
	
	//define gridN
	int gridN = N+2;
	
	//allocate memory for Mh
	double	*Mh	= new double [gridN];
	
	//copy yk into mode
	int i;
	for(i = 0; i < gridN; i++)
		mode[i] = yk[i];
	
	//calculate mean shift vector at yk
	LatticeMSVector(Mh, mode);
	
	//calculate mvAbs = |Mh|^2
	double mvAbs = 0;
	for(i = 0; i < gridN; i++)
		mvAbs	+= Mh[i]*Mh[i];
	
	//shift mode until convergence (mvAbs = 0)...
	int iterationCount = 1;
	while((mvAbs >= EPSILON)&&(iterationCount < LIMIT))
	{
		//shift mode...
		for(i = 0; i < gridN; i++)
			mode[i]	+= Mh[i];
		
		//re-calculate mean shift vector at new
		//window location have center defined by
		//mode
		LatticeMSVector(Mh, mode);
		
		//calculate mvAbs = |Mh|^2
		mvAbs = 0;
		for(i = 0; i < gridN; i++)
			mvAbs	+= Mh[i]*Mh[i];
		
		//increment interation count...
		iterationCount++;
		
	}
	
	//shift mode...
	for(i = 0; i < gridN; i++)
		mode[i]	+= Mh[i];
	
	//de-allocate memory
	delete [] Mh;
	
	//done.
	return;
	
}

/*@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@*/
/*@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@*/
/*@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@    PROTECTED METHODS    @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@*/
/*@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@*/
/*@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@*/

  /*/\/\/\/\/\/\/\/\/\/\/\/\/\/\*/
  /* Mean Shift: Using kd-Tree  */
  /*\/\/\/\/\/\/\/\/\/\/\/\/\/\/*/

/*******************************************************/
/*Mean Shift Vector                                    */
/*******************************************************/
/*Computes the mean shift vector at a window location  */
/*yk using input data set x using a custom, user defin-*/
/*ed kernel.                                           */
/*******************************************************/
/*Pre:                                                 */
/*      - input data has been uploaded into the private*/
/*        data members of the MeanShift class          */
/*      - a window center yk has been defined          */
/*      - uniformKernel indicates the which type of    */
/*        kernel to be used by this procedure: uniform */
/*        or general                                   */
/*Post:                                                */
/*      - the mean shift vector calculated at yk       */
/*        using a either a custom, user defined kernel */
/*        or a uniform kernel is returned              */
/*******************************************************/

void MeanShift::MSVector(double *Mh_ptr, double *yk_ptr)
{
	
	// Declare Variables
	int i,j;
	
	// Initialize mean shift vector
	for(i = 0; i < N; i++)
		Mh_ptr[i] = 0;
	
	// Initialize wsum to zero, the sum of the weights of each
	// data point found to lie within the search window (sphere)
	wsum = 0;
	
	// Build Range Vector using h[i] and yk
	
	int s = 0;
	
	// The flag uniformKernel is used to determine which
	// kernel function is to be used in the calculation
	// of the mean shift vector
	if(uniformKernel)
    {
		for(i = 0; i < kp; i++)
		{
			for(j = 0; j < P[i]; j++)
			{
				range[2*(s+j)  ] = (float)(yk_ptr[s+j] - h[i]);
				range[2*(s+j)+1] = (float)(yk_ptr[s+j] + h[i]);
			}
			s += P[i];
		}
    }
	else
    {
		for(i = 0; i < kp; i++)
		{
			for(j = 0; j < P[i]; j++)
			{
				range[2*(s+j)  ] = (float)(yk_ptr[s+j] - h[i]*float(sqrt(offset[i])));
				range[2*(s+j)+1] = (float)(yk_ptr[s+j] + h[i]*float(sqrt(offset[i])));
			}
			s += P[i];
		}
    }
	
	// Traverse through the data set x, performing the
	// weighted sum of each point xi that lies within
	// the search window (sphere) using a general,
	// user defined kernel or uniform kernel depending
	// on the uniformKernel flag
	if(uniformKernel)
		uniformSearch(root, 0, Mh_ptr, yk_ptr);
	else
		generalSearch(root, 0, Mh_ptr, yk_ptr);
	
	// Calculate the mean shift vector using Mh and wsum
	for(i = 0; i < N; i++)
    {
		
		// Divide Sum by wsum
		Mh_ptr[i] /= wsum;
		
		// Calculate mean shift vector: Mh(yk) = y(k+1) - y(k)
		Mh_ptr[i] -= yk_ptr[i];
		
    }
	
	//done.
	return;
	
}

  /*/\/\/\/\/\/\/\/\/\/\/\/\/\/\*/
  /*  Mean Shift: Using Lattice */
  /*\/\/\/\/\/\/\/\/\/\/\/\/\/\/*/

/*******************************************************/
/*Lattice Mean Shift Vector                            */
/*******************************************************/
/*Computes the mean shift vector at a specfied window  */
/*yk using the lattice data structure.                 */
/*******************************************************/
/*Pre:                                                 */
/*      - Mh_ptr and yh_ptr are arrays of doubles con- */
/*        aining N+2 elements                          */
/*      - Mh_ptr is the mean shift vector calculated   */
/*        at window center yk_ptr                      */
/*Post:                                                */
/*      - the mean shift vector at the window center   */
/*        pointed to by yk_ptr has been calculated and */
/*        stored in the memory location pointed to by  */
/*        Mh_ptr                                       */
/*******************************************************/

void MeanShift::LatticeMSVector(double *Mh_ptr, double *yk_ptr)
{
	
	// Initialize mean shift vector
	register int i;
	for(i = 0; i < N+2; i++)
		Mh_ptr[i] = 0;
	
	// Initialize wsum
	wsum = 0;
	
	// Perform lattice search summing
	// all the points that lie within the search
	// window defined using the kernel specified
	//by uniformKernel
	if(uniformKernel)
		uniformLSearch(Mh_ptr, yk_ptr);
	else
		generalLSearch(Mh_ptr, yk_ptr);
	
	// Compute mean shift vector using sum computed
	// by lattice search, wsum, and yk_ptr:
	// Mh = Mh/wsum - yk_ptr 
	
	if (wsum > 0)
	{
		for(i = 0; i < N+2; i++)
			Mh_ptr[i] = Mh_ptr[i]/wsum - yk_ptr[i];
	}
	else
	{
		for(i = 0; i < N+2; i++)
			Mh_ptr[i] = 0;
	}
	
	// done.
	return;
	
}

/*******************************************************/
/*Optimized Lattice Mean Shift Vector                  */
/*******************************************************/
/*Computes the mean shift vector at a specfied window  */
/*yk using the lattice data structure. Also the points */
/*that lie within the window are stored into the basin */
/*of attraction structure used by the optimized mean   */
/*shift algorithms.									   */
/*******************************************************/
/*Pre:                                                 */
/*      - Mh_ptr and yh_ptr are arrays of doubles con- */
/*        aining N+2 elements                          */
/*      - Mh_ptr is the mean shift vector calculated   */
/*        at window center yk_ptr                      */
/*Post:                                                */
/*      - the mean shift vector at the window center   */
/*        pointed to by yk_ptr has been calculated and */
/*        stored in the memory location pointed to by  */
/*        Mh_ptr                                       */
/*      - the data points lying within h of of yk_ptr  */
/*        have been stored into the basin of attract-  */
/*        ion data structure.                          */
/*******************************************************/

void MeanShift::OptLatticeMSVector(double *Mh_ptr, double *yk_ptr)
{
	
	// Initialize mean shift vector
	register int i;
	for(i = 0; i < N+2; i++)
		Mh_ptr[i] = 0;
	
	// Initialize wsum
	wsum = 0;
	
	// Perform lattice search summing
	// all the points that lie within the search
	// window defined using the kernel specified
	//by uniformKernel
	if(uniformKernel)
		optUniformLSearch(Mh_ptr, yk_ptr);
	else
		optGeneralLSearch(Mh_ptr, yk_ptr);
	
	// Compute mean shift vector using sum computed
	// by lattice search, wsum, and yk_ptr:
	// Mh = Mh/wsum - yk_ptr 
	
   if (wsum > 0)
   {
	   for(i = 0; i < N+2; i++)
   		Mh_ptr[i] = Mh_ptr[i]/wsum - yk_ptr[i];
   } else
   {
      for (i=0; i< N+2; i++)
         Mh_ptr[i] = 0;
   }
	
	// done.
	return;
	
}

  /*/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\*/
  /*** Kernel-Input Data Consistency  ***/
  /*\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/*/

/*******************************************************/
/*Class Consistency Check                              */
/*******************************************************/
/*Checks the state of the class prior to the applicat- */
/*ion of mean shift.                                   */
/*******************************************************/
/*Pre:                                                 */
/*      - iN is the specified dimension of the input,  */
/*        iN = N for a general input data set, iN = N  */
/*        + 2 for a input set defined using a lattice  */
/*Post:                                                */
/*      - if the kernel has not been created, an input */
/*        has not been defined and/or the specified    */
/*        input dimension (iN) does not match that of  */
/*        the kernel a fatal error is flagged.         */
/*******************************************************/

void MeanShift::classConsistencyCheck(int iN, bool usingLattice)
{
	
	//make sure that kernel has been created...
	if(class_state.KERNEL_DEFINED == false)
	{
		ErrorHandler("MeanShift", "classConsistencyCheck", "Kernel not created.");
		return;
	}
	
	//make sure input data set has been loaded into mean shift object...
	if((class_state.INPUT_DEFINED == false)&&(!usingLattice))
	{
		ErrorHandler("MeanShift", "classConsistencyCheck", "No input data specified.");
		return;
	}
	
	//make sure that the lattice is defined if it is being used
	if((class_state.LATTICE_DEFINED == false)&&(usingLattice))
	{
		ErrorHandler("MeanShift", "classConsistencyCheck", "Latice not created.");
		return;
	}
	
	//make sure that dimension of the kernel and the input data set
	//agree
	
	//calculate dimension of kernel (kN)
	int i, kN	= 0;
	for(i = 0; i < kp; i++)
		kN	+= P[i];
	
	//perform comparison...
	if(iN != kN)
	{
		ErrorHandler("MeanShift", "classConsitencyCheck", "Kernel dimension does not match defined input data dimension.");
		return;
	}
	
	//done.
	return;
	
}

  /*/\/\/\/\/\/\/\/\/\/\/\/\/\*/
  /*** Class Error Handler  ***/
  /*\/\/\/\/\/\/\/\/\/\/\/\/\/*/

/*******************************************************/
/*Error Handler                                        */
/*******************************************************/
/*Class error handler.                                 */
/*******************************************************/
/*Pre:                                                 */
/*      - className is the name of the class that fl-  */
/*        agged an error                               */
/*      - methodName is the name of the method that    */
/*        flagged an error                             */
/*      - errmsg is the error message given by the     */
/*        calling function                             */
/*Post:                                                */
/*      - the error message errmsg is flagged on beh-  */
/*        ave of method methodName belonging to class  */
/*        className:                                   */
/*                                                     */
/*        (1) ErrorMessage has been updated with the   */
/*            appropriate error message using the arg- */
/*            ments passed to this method.             */
/*        (2) ErrorStatus is set to ERROR              */
/*            (ErrorStatus = 1)                        */
/*******************************************************/

void MeanShift::ErrorHandler(char *className, char *methodName, char* errmsg)
{
	
	//store trace into error message
	strcpy(ErrorMessage, className);
	strcat(ErrorMessage, "::");
	strcat(ErrorMessage, methodName);
	strcat(ErrorMessage, " Error: ");
	
	//store message into error message
	strcat(ErrorMessage, errmsg);
	
	//set error status to ERROR
	ErrorStatus = EL_ERROR;
	
	
}

/*@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@*/
/*@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@*/
/*@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@     PRIVATE METHODS     @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@*/
/*@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@*/
/*@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@*/

  /*/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\*/
  /*** Kernel Creation/Manipulation ***/
  /*\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/*/

/*******************************************************/
/*Generate Lookup Table                                */
/*******************************************************/
/*A weight function look up table is generated.        */
/*******************************************************/
/*Pre:                                                 */
/*      - kernel is an array of kernelTypes specifying */
/*        the type of kernel to be used on each sub-   */
/*        space of the input data set x                */
/*      - kp is the total number of subspaces used to  */
/*        the input data set x                         */
/*      - the above information has been pre-loaded    */
/*        into the MeanShift class private members     */
/*Post:                                                */
/*      - a lookup table is generated for the weight   */
/*        function of the resulting kernel             */
/*      - uniformKernel is set to true if the kernel   */
/*        to be used is uniform, false is returned     */
/*        otherwise                                    */
/*      - if a user defined weight function is requred */
/*        for a given subspace but not defined in the  */
/*        user defined weight function list, an error  */
/*        is flagged and the program is halted         */
/*******************************************************/

void MeanShift::generateLookupTable( void )
{
	
	// Declare Variables
	int i,j;
	
	// Allocate memory for lookup table w
	w = new double*[kp];
	
	// Traverse through kernel generating weight function
	// lookup table w
	
	// Assume kernel is uniform
	uniformKernel = true;
	
	for(i = 0; i < kp; i++)
    {
		switch(kernel[i])
		{
			// *Uniform Kernel* has weight funciton w(u) = 1
			// therefore, a weight funciton lookup table is
			// not needed for this kernel --> w[i] = NULL indicates
			// this
		case Uniform:
			
			w        [i] = NULL;  //weight function not needed for this kernel
			offset   [i] =    1;  //uniform kernel has u < 1.0
			increment[i] =    1;  //has no meaning
			break;
			
			// *Gaussian Kernel* has weight function w(u) = constant*exp(-u^2/[2h[i]^2])
		case Gaussian:
			
			// Set uniformKernel to false
			uniformKernel = false;
			
			// generate weight function using expression,
			// exp(-u/2), where u = norm(xi - x)^2/h^2
			
			// Allocate memory for weight table
			w[i] = new double [GAUSS_NUM_ELS+1];
			
			for(j = 0; j <= GAUSS_NUM_ELS; j++)
				w[i][j] = exp(-j*GAUSS_INCREMENT/2);
			
			// Set offset = offset^2, and set increment
			offset   [i] = (float)(GAUSS_LIMIT*GAUSS_LIMIT);
			increment[i] = GAUSS_INCREMENT;
			
			// done
			break;
			
			// *User Define Kernel* uses the weight function wf(u)
		case UserDefined:
			
			// Set uniformKernel to false
			uniformKernel = false;
			
			// Search for user defined weight function
			// defined for subspace (i+1)
			cur = head;
			while((cur)&&(cur->subspace != (i+1)))
				cur = cur->next;
			
			// If a user defined subspace has not been found
			// for this subspace, flag an error
			if(cur == NULL)
			{
				fprintf(stderr, "\ngenerateLookupTable Fatal Error: User defined kernel for subspace %d undefined.\n\nAborting Program.\n\n", i+1);
				exit(1);
			}
			
			// Otherwise, copy weight function lookup table to w[i]
			w[i] = new double [cur->sampleNumber+1];
			for(j = 0; j <= cur->sampleNumber; j++)
				w[i][j] = cur->w[j];
			
			// Set offset and increment accordingly
			offset   [i] = (float)(cur->halfWindow);
			increment[i] = cur->halfWindow/(float)(cur->sampleNumber);
			
			// done
			break;
			
		default:
			
			ErrorHandler("MeanShift", "generateLookupTable", "Unknown kernel type.");
			
		}
		
    }
}

/*******************************************************/
/*Destroy Kernel                                       */
/*******************************************************/
/*Destroys and initializes kernel.                     */
/*******************************************************/
/*Post:                                                */
/*      - memory for the kernel private data members   */
/*        have been destroyed and the kernel has been  */
/*        initialized for re-use.                      */
/*******************************************************/

void MeanShift::DestroyKernel( void )
{
	
	//de-allocate memory...
	if(kernel)	delete	[] kernel;
	if     (h)	delete	[] h;
	if     (P)	delete	[] P;
	if (range)	delete	[] range;

   if (uv) delete [] uv;
   if(increment) delete [] increment;
   if (offset) delete [] offset;
   
   if (kp>0)
   {
      if (w)
      {
         int i;
         for (i=0; i<kp; i++)
            delete [] w[i];
         delete [] w;
      }
      w = NULL;
   }
   
	//intialize kernel for re-use...
	kp		= 0;
	kernel	= NULL;
	h		= NULL;
	P		= NULL;
	range	= NULL;

   increment = NULL;
   uv = NULL;
   offset = NULL;
	
	//done.
	return;
	
}

  /*/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\*/
  /*** Input Data Initialization/Destruction  ***/
  /*\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/*/

/*******************************************************/
/*Create Binary Search Tree                            */
/*******************************************************/
/*Uploads input data set x into a kd-BST.              */
/*******************************************************/
/*Pre:                                                 */
/*      - x is a one dimensional array of L N-dimensi- */
/*        onal data points                             */
/*Post:                                                */
/*      - x has been uploaded into a balanced kd-BST   */
/*        data structure for use by the mean shift     */
/*        procedure                                    */
/*******************************************************/

void MeanShift::CreateBST( void )
{
	
	// Create BST using data....
	
	// Allocate memory for tree
	forest = new tree[L];
	
	// Populate 'forest' of tree's with
	// the values stored in x
	int i;
	for(i = 0; i < L; i++)
    {
		forest[i].x      = &data[i*N];
		forest[i].right  = NULL;
		forest[i].left   = NULL;
		forest[i].parent = NULL;
    }
	
	// Build balanced Nd-tree from the
	// forest of trees generated above
	// retaining the root of this tree
	
	root = BuildKDTree(forest, L, 0, NULL);
	
	//done.
	return;
	
}

/*******************************************************/
/*Initialize Input                                     */
/*******************************************************/
/*Allocates memory for and initializes the input data  */
/*structure.                                           */
/*******************************************************/
/*Pre:                                                 */
/*      - x is a floating point array of L, N dimens-  */
/*        ional input data points                      */
/*Post:                                                */
/*      - memory has been allocated for the input data */
/*        structure and x has been stored using into   */
/*        the mean shift class using the resulting     */
/*        structure.                                   */
/*******************************************************/

void MeanShift::InitializeInput(float *x)
{
	
	//allocate memory for input data set
	if(!(data = new float [L*N]))
	{
		ErrorHandler("MeanShift", "InitializeInput", "Not enough memory.");
		return;
	}
	
	//copy x into data
	int i;
	for(i = 0; i < L*N; i++)
		data[i]	= x[i];
	
	//done.
	return;
	
}

/*******************************************************/
/*Reset Input                                          */
/*******************************************************/
/*De-allocates memory for and re-intializes input data */
/*structure.                                           */
/*******************************************************/
/*Post:                                                */
/*      - the memory of the input data structure has   */
/*        been de-allocated and this strucuture has    */
/*        been initialized for re-use.                 */
/*******************************************************/

void MeanShift::ResetInput( void )
{
	
	//de-allocate memory of input data structure (BST)
	if(data)	delete [] data;
	if(forest)	delete [] forest;
	
	//initialize input data structure for re-use
	data	= NULL;
	forest	= NULL;
	root	= NULL;
	L		= 0;
	N		= 0;
	width	= 0;
	height	= 0;
	
	//re-set class input to indicate that
	//an input is not longer stored by
	//the private data members of this class
	class_state.INPUT_DEFINED	= class_state.LATTICE_DEFINED = false;
	
}

  /*/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\*/
  /*** k-dimensional Binary Search Tree ***/
  /*\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/*/

/*******************************************************/
/*Build KD Tree (for Tree Structure)                   */
/*******************************************************/
/*Builds a KD Tree given a forest of tree's.           */
/*******************************************************/
/*Pre:                                                 */
/*      - subset is a subset of L un-ordered tree nodes*/
/*        each containing an N-dimensional data point  */
/*      - d is the depth of the subset, used to specify*/
/*        the dimension used to construct the tree at  */
/*        the given depth                              */
/*      - parent is the parent tree of subset          */
/*Post:                                                */
/*      - a balanced KD tree has been constructed using*/
/*        the forest subset, the root of this tree has */
/*        been returned                                */
/*******************************************************/

tree *MeanShift::BuildKDTree(tree *subset, int length, int d, tree* parent)
{
	
	// If the subset is a single tree
	// then return this tree otherwise
	// partition the subset and place
	// these subsets recursively into
	// the left and right sub-trees having
	// their root specified by the median
	// of this subset in dimension d
	if(length == 1)
	{
		subset->parent = parent;
		return subset;
	}
	else if(length > 1)
    {
		
		// Sort Subset
		QuickMedian(subset, 0, length-1, d);
		
		// Get Median of Subset and Partition
		// it into two sub-trees - create
		// a tree with its root being the median
		// of the subset and its left and right
		// children being the medians of the subsets
		int median            = length/2;
		subset[median].parent = parent;
		subset[median].left   = BuildKDTree(subset           , median         , (d+1)%N, &subset[median]);
		subset[median].right  = BuildKDTree(&subset[median+1], length-median-1, (d+1)%N, &subset[median]);
		
		// Output tree structure
		return &subset[median];
		
    }
	else
		return NULL;
	
	//done.
	
}

/*******************************************************/
/*Quick Median (for Tree Structure)                    */
/*******************************************************/
/*Finds the median element in an un-ordered set, re-   */
/*structuring the set such that points less than the   */
/*median point are located to the left of the median   */
/*and points greater than the median point are located */
/*to the right.                                        */
/*******************************************************/
/*Pre:                                                 */
/*      - arr is a subset of tree nodes whose leftmost */
/*        element is specified by left and rightmost   */
/*        element is specified by left                 */
/*      - d is the dimension of the data set stored by */
/*        the tree structure that is used to find      */
/*        the median                                   */
/*Post:                                                */
/*      - the median point is found and the subset     */
/*        of trees is re-ordered such that all trees   */
/*        whose data points with d dimensional value   */
/*        less than that of the median tree node are   */
/*        located to the left of the median tree node, */
/*        otherwise they are located to the right      */
/*******************************************************/

void MeanShift::QuickMedian(tree *arr, int left, int right, int d)
{
	unsigned long k;
	unsigned long n;
	float* a;
	float* temp;
	n = right-left+1;
	k = n/2 + 1;
	unsigned long i, ir, j, l, mid;
	
	l = 1;
	ir = n;
	for (;;)
	{
		if (ir <= l+1)
		{
			if (ir == l+1 && arr[ir-1].x[d] < arr[l-1].x[d])
			{
				SWAP(arr[l-1].x, arr[ir-1].x)
			}
			return;
		} else
		{
			mid = (l+ir) >> 1;
			SWAP(arr[mid-1].x, arr[l+1-1].x)
				if (arr[l-1].x[d] > arr[ir-1].x[d])
				{
					SWAP(arr[l-1].x, arr[ir-1].x)
				}
				if (arr[l+1-1].x[d] > arr[ir-1].x[d])
				{
					SWAP(arr[l+1-1].x, arr[ir-1].x)
				}
				if (arr[l-1].x[d] > arr[l+1-1].x[d])
				{
					SWAP(arr[l-1].x, arr[l+1-1].x)
				}
				i = l+1;
				j = ir;
				a = arr[l+1-1].x;
				for (;;) {
					do i++; while (arr[i-1].x[d] < a[d]);
					do j--; while (arr[j-1].x[d] > a[d]);
					if (j<i) break;
					SWAP(arr[i-1].x, arr[j-1].x)
				}
				arr[l+1-1].x = arr[j-1].x;
				arr[j-1].x = a;
				if (j>=k) ir = j-1;
				if (j<=k) l = i;
		}
	}
}

  /*/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\*/
  /*** Mean Shift: Using kd-Tree  ***/
  /*\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/*/

/*******************************************************/
/*Uniform Search                                       */
/*******************************************************/
/*Searches the input data using a kd-tree, performs the*/
/*sum on the data within the Hypercube defined by the  */
/*tree using a uniform kernel.                         */
/*******************************************************/
/*Pre:                                                 */
/*      - gt is a possibly NULL pointer to a kd tree   */
/*      - Mh_ptr is a pointer to the mean shift vector */
/*        being calculated                             */
/*      - yk_ptr is a pointer to the current window    */
/*        center location                              */
/*      - gd is the depth of the current subtree       */
/*Post:                                                */
/*      - the mean of the points within the Hypercube  */
/*        of the kd tree is computed using a uniform   */
/*        kernel                                       */
/*******************************************************/

void MeanShift::uniformSearch(tree *gt, int gd, double *Mh_ptr, double *yk_ptr)
{
	tree* c_t;
	int c_d;
	int i;
	int actionType;
	
	c_t = gt;
	c_d = gd;
	actionType = 0;
	
	double el, diff;
	int k, j, s;
	
	while (c_t != NULL)
	{
		switch(actionType) {
		case 0: // forward
			if ((c_t->x[c_d] > range[2*c_d]) && ((c_t->left) != NULL))
			{
				c_t = c_t->left;
				c_d = (c_d+1)%N;
			} else
			{
				actionType = 1; 
			}
			break;
		case 1: // backleft
			
			for(i = 0; i < N; i++)
			{
				if((c_t->x[i] < range[2*i])||(c_t->x[i] > range[2*i+1]))
					break;
			}
			
			if(i == N)
			{
				
				// ***     Visit Tree       ***
				
				// Check if xi is in the window centered about yk_ptr
				// If so - use it to compute y(k+1)
				diff = 0;
				j = 0;
				s = 0;
				while((diff < 1.0)&&(j < kp)) // Partial Distortion Search (PDS)
				{
					
					// test each sub-dimension independently
					diff  = 0;
					for(k = 0; k < P[j]; k++)
					{
						el = (c_t->x[s+k] - yk_ptr[s+k])/h[j];
						diff += el*el;
					}
					
					s += P[j];                        // next subspace
					j++;
					
				}
				
				if(diff < 1.0)
				{
					wsum += 1;
					for(j = 0; j < N; j++)
						Mh_ptr[j] += c_t->x[j];
				}
				
			}
			if ((c_t->x[c_d] < range[2*c_d+1]) && ((c_t->right) != NULL))
			{
				c_t = c_t->right;
				c_d = (c_d+1)%N;
				actionType = 0;
			} else
			{
				actionType = 2;
			}
			break;
		case 2: // backright
			c_d = (c_d+N-1)%N;
			
			if (c_t->parent == NULL)
			{
				c_t = NULL;
				break;
			}
			
			if (c_t->parent->left == c_t)
				actionType = 1;
			else
				actionType = 2;
			c_t = c_t->parent;
			break;
		}
	}
}

/*******************************************************/
/*General Search                                       */
/*******************************************************/
/*Searches the input data using a kd tree, performs the*/
/*sum on the data within the Hypercube defined by the  */
/*tree using a general kernel.                         */
/*******************************************************/
/*Pre:                                                 */
/*      - gt is a possibly NULL pointer to a kd tree   */
/*      - Mh_ptr is a pointer to the mean shift vector */
/*        being calculated                             */
/*      - yk_ptr is a pointer to the current window    */
/*        center location                              */
/*      - gd is the depth of the current subtree       */
/*Post:                                                */
/*      - the mean of the points within the Hypercube  */
/*        of the kd tree is computed using a general   */
/*        kernel                                       */
/*******************************************************/

void MeanShift::generalSearch(tree *gt, int gd, double *Mh_ptr, double *yk_ptr)
{
	tree* c_t;
	int c_d;
	int i;
	int actionType;
	
	c_t = gt;
	c_d = gd;
	actionType = 0;
	
	double el, diff, u, tw, y0, y1;
	int k, j, s, x0, x1;
	
	while (c_t != NULL)
	{
		switch(actionType) {
		case 0: // forward
			if ((c_t->x[c_d] > range[2*c_d]) && ((c_t->left) != NULL))
			{
				c_t = c_t->left;
				c_d = (c_d+1)%N;
			} else
			{
				actionType = 1; 
			}
			break;
		case 1: // backleft
			
			for(i = 0; i < N; i++)
			{
				if((c_t->x[i] < range[2*i])||(c_t->x[i] > range[2*i+1]))
					break;
			}
			
			if(i == N)
			{
				
				// ***      Visit Tree      ***
				
				// Check if xi is in the window centered about yk_ptr
				// If so - use it to compute y(k+1)
				s = 0;
				for(j = 0; j < kp; j++)
				{
					
					// test each sub-dimension independently
					diff  = 0;
					for(k = 0; k < P[j]; k++)
					{
						el = (c_t->x[s+k] - yk_ptr[s+k])/h[j];
						diff += uv[s+k] = el*el;      // Update uv and diff
						if(diff >= offset[j])         // Partial Distortion Search (PDS)
							break;
					}
					
					if(diff >= offset[j])             // PDS
						break;
					
					s += P[j];                        // next subspace
					
				}
				
				// j == kp indicates that all subspaces passed the test:
				// the data point is within the search window
				if(j == kp) j--;
				if(diff < offset[j])
				{
					
					// Initialize total weight to 1
					tw = 1;
					
					// Calculate weight factor using weight function
					// lookup tables and uv
					s = 0;
					for(j = 0; j < kp; j++)
					{
						if(kernel[j]) // not uniform kernel
						{
							// Compute u[i]
							u = 0;
							for(k = 0; k < P[j]; k++)
								u += uv[s+k];
							
							// Accumulate tw using calculated u
							// and weight function lookup table
							
							// Linear interpolate values given by
							// lookup table
							
							// Calculate x0 and x1, the points surounding
							// u
							x0 = (int)(u/increment[j]);
							x1 = x0+1;
							
							// Get y0 and y1 from the lookup table
							y0 = w[j][x0];
							y1 = w[j][x1];
							
							// Accumulate tw using linear interpolation
							tw *= (((double)(x1)*increment[j] - u)*y0+(u - (double)(x0)*increment[j])*y1)/(double)(x1*increment[j] - x0*increment[j]);
							
						}
						s += P[j];                               // next subspace
					}
					
					// Perform weighted sum using xi
					for(j = 0; j < N; j++)
						Mh_ptr[j] += tw*c_t->x[j];
					
					// Increment wsum by tw
					wsum += tw;
					
				}
			}
			if ((c_t->x[c_d] < range[2*c_d+1]) && ((c_t->right) != NULL))
			{
				c_t = c_t->right;
				c_d = (c_d+1)%N;
				actionType = 0;
			} else
			{
				actionType = 2;
			}
			break;
		case 2: // backright
			c_d = (c_d+N-1)%N;
			
			if (c_t->parent == NULL)
			{
				c_t = NULL;
				break;
			}
			
			if (c_t->parent->left == c_t)
				actionType = 1;
			else
				actionType = 2;
			c_t = c_t->parent;
			break;
      }
   }
}

  /*/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\*/
  /***  Mean Shift: Using Lattice ***/
  /*\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/*/

/*******************************************************/
/*Uniform Lattice Search                               */
/*******************************************************/
/*Performs search on data set for all points lying     */
/*within the search window defined using a uniform     */
/*kernel. Their point-wise sum and count is computed   */
/*and returned.                                        */
/*                                                     */
/*NOTE: This method is the only method in the          */
/*      MeanShift class that uses the weight           */
/*      map asside from optUniformLSearch.             */
/*******************************************************/
/*Pre:                                                 */
/*      - Mh_ptr is a length N array of doubles        */
/*      - yk_ptr is a length N array of doubles        */
/*      - Mh_ptr is the sum of the data points found   */
/*        within search window having center yk_ptr    */
/*Post:                                                */
/*      - a search on the data set using the lattice   */
/*        has been performed, and all points found to  */
/*        lie within the search window defined using   */
/*        a uniform kernel are summed and counted.     */
/*      - their point wise sum is pointed to by Mh_ptr */
/*        and their count is stored by wsum.           */
/*******************************************************/

void MeanShift::uniformLSearch(double *Mh_ptr, double *yk_ptr)
{
	
	//Declare variables
	register int	i, j, k;
	int				s, p, dataPoint, lN;
	double			diff, el, dx, dy, tx, weight;
	
	//Define lattice data dimension...
	lN	= N + 2;
	
	//Define bounds of lattice...
	
	//the lattice is a 2dimensional subspace whose
	//search window bandwidth is specified by
	//h[0]:
	tx = yk_ptr[0] - h[0] + DELTA + 0.99;
	if (tx < 0)
		LowerBoundX = 0;
	else
		LowerBoundX = (int) tx;
	tx = yk_ptr[1] - h[0] + DELTA + 0.99;
	if (tx < 0)
		LowerBoundY = 0;
	else
		LowerBoundY = (int) tx;
	tx = yk_ptr[0] + h[0] - DELTA;
	if (tx >= width)
		UpperBoundX = width-1;
	else
		UpperBoundX = (int) tx;
	tx = yk_ptr[1] + h[0] - DELTA;
	if (tx >= height)
		UpperBoundY = height - 1;
	else
		UpperBoundY = (int) tx;
	
	//Perform search using lattice
	for(i = LowerBoundY; i <= UpperBoundY; i++)
		for(j = LowerBoundX; j <= UpperBoundX; j++)
		{
			
			//get index into data array
			dataPoint = N*(i*width+j);
			
			//Determine if inside search window
			k		= 1;
			s		= 0;
			dx      = j - yk_ptr[0];
			dy      = i - yk_ptr[1];
			diff	= (dx*dx+dy*dy)/(h[0]*h[0]);
			while((diff < 1.0)&&(k != kp)) // Partial Distortion Search
			{
				//Calculate distance squared of sub-space s	
				diff = 0;
				for(p = 0; p < P[k]; p++)
				{
					el    = (data[dataPoint+p+s]-yk_ptr[p+s+2])/h[k];               
					if((!p)&&(yk_ptr[2] > 80))
						diff += 4*el*el;
					else
						diff += el*el;
				}
				
				//next subspace
				s += P[k];
				k++;
			}
			
			//if its inside search window perform sum and count
			if(diff < 1.0)
			{
				weight = 1 - weightMap[i*width+j];
				Mh_ptr[0] += weight*j;
				Mh_ptr[1] += weight*i;
				for(k = 2; k < lN; k++)
					Mh_ptr[k] += weight*data[dataPoint+k-2];
				wsum += weight;
			}
			//done.
		}
		//done.		
		return;
		
}

/*******************************************************/
/*Optimized Uniform Latice Search                      */
/*******************************************************/
/*Performs search on data set for all points lying     */
/*within the search window defined using a uniform     */
/*kernel. Their point-wise sum and count is computed   */
/*and returned. Also the points that lie within the    */
/*window are stored into the basin of attraction stru- */
/*cture used by the optimized mean shift algorithms.   */
/*                                                     */
/*NOTE: This method is the only method in the          */
/*      MeanShift class that uses the weight           */
/*      map asside from uniformLSearch.                */
/*******************************************************/
/*Pre:                                                 */
/*      - Mh_ptr is a length N array of doubles        */
/*      - yk_ptr is a length N array of doubles        */
/*      - Mh_ptr is the sum of the data points found   */
/*        within search window having center yk_ptr    */
/*Post:                                                */
/*      - a search on the data set using the latice    */
/*        has been performed, and all points found to  */
/*        lie within the search window defined using   */
/*        a uniform kernel are summed and counted.     */
/*      - their point wise sum is pointed to by Mh_ptr */
/*        and their count is stored by wsum.           */
/*      - the data points lying within h of of yk_ptr  */
/*        have been stored into the basin of attract-  */
/*        ion data structure.                          */
/*******************************************************/

void MeanShift::optUniformLSearch(double *Mh_ptr, double *yk_ptr)
{
	
	//Declare variables
	register int	i, j, k;
	int				s, p, dataPoint, pointIndx, lN;
	double			diff, el, dx, dy, tx, weight;
	
	//Define latice data dimension...
	lN	= N + 2;
	
	//Define bounds of latice...
	
	//the latice is a 2dimensional subspace whose
	//search window bandwidth is specified by
	//h[0]:
	tx = yk_ptr[0] - h[0] + DELTA + 0.99;
	if (tx < 0)
		LowerBoundX = 0;
	else
		LowerBoundX = (int) tx;
	tx = yk_ptr[1] - h[0] + DELTA + 0.99;
	if (tx < 0)
		LowerBoundY = 0;
	else
		LowerBoundY = (int) tx;
	tx = yk_ptr[0] + h[0] - DELTA;
	if (tx >= width)
		UpperBoundX = width-1;
	else
		UpperBoundX = (int) tx;
	tx = yk_ptr[1] + h[0] - DELTA;
	if (tx >= height)
		UpperBoundY = height - 1;
	else
		UpperBoundY = (int) tx;
	
	//Perform search using latice
	for(i = LowerBoundY; i <= UpperBoundY; i++)
		for(j = LowerBoundX; j <= UpperBoundX; j++)
		{
			
			//get index into data array
			pointIndx	= i*width+j;
			dataPoint	= N*(pointIndx);
			
			//Determine if inside search window
			k		= 1;
			s		= 0;
			dx      = j - yk_ptr[0];
			dy      = i - yk_ptr[1];
			diff	= (dx*dx+dy*dy)/(h[0]*h[0]);
			while((diff < 1.0)&&(k != kp)) // Partial Distortion Search
			{
				//Calculate distance squared of sub-space s	
				diff = 0;
				for(p = 0; p < P[k]; p++)
				{
					el    = (data[dataPoint+p+s]-yk_ptr[p+s+2])/h[k];
					if((!p)&&(yk_ptr[2] > 80))
						diff += 4*el*el;
					else               
					   diff += el*el;
				}
				
				//next subspace
				s += P[k];
				k++;
			}
			
			//if its inside search window perform sum and count
			if(diff < 1.0)
			{
				weight = 1 - weightMap[i*width+j];
				Mh_ptr[0] += weight*j;
				Mh_ptr[1] += weight*i;
				for(k = 2; k < lN; k++)
					Mh_ptr[k] += weight*data[dataPoint+k-2];
				wsum += weight;
				
				//set basin of attraction mode table
            if (diff < 0.5)
            {
				   if(modeTable[pointIndx] == 0)
				   {
   					pointList[pointCount++]	= pointIndx;
					   modeTable[pointIndx]	= 2;
				   }
            }
				
			}
			
			//done.
			
		}
		
		//done.		
		return;
		
}

/*******************************************************/
/*General Lattice Search                               */
/*******************************************************/
/*Performs search on data set for all points lying     */
/*within the search window defined using a general     */
/*kernel. Their point-wise sum and count is computed   */
/*and returned.                                        */
/*******************************************************/
/*Pre:                                                 */
/*      - Mh_ptr is a length N array of doubles        */
/*      - yk_ptr is a length N array of doubles        */
/*      - Mh_ptr is the sum of the data points found   */
/*        within search window having center yk_ptr    */
/*Post:                                                */
/*      - a search on the data set using the lattice   */
/*        has been performed, and all points found to  */
/*        lie within the search window defined using   */
/*        a general kernel are summed and counted      */
/*      - their point wise sum is pointed to by Mh_ptr */
/*        and their count is stored by wsum            */
/*******************************************************/

void MeanShift::generalLSearch(double *Mh_ptr, double *yk_ptr)
{
	
	//Declare variables
	register int i, j, k;
	int			 s, p, dataPoint, lN, x0, x1;
	double		 diff, el, dx, dy, tw, u, y0, y1, tx;
	
	//Define lattice data dimension...
	lN	= N + 2;
	
	//Define bounds of lattice...
	
	//the lattice is a 2dimensional subspace whose
	//search window bandwidth is specified by
	//h[0]:
	tx = yk_ptr[0] - h[0] + DELTA + 0.99;
	if (tx < 0)
		LowerBoundX = 0;
	else
		LowerBoundX = (int) tx;
	tx = yk_ptr[1] - h[0] + DELTA + 0.99;
	if (tx < 0)
		LowerBoundY = 0;
	else
		LowerBoundY = (int) tx;
	tx = yk_ptr[0] + h[0] - DELTA;
	if (tx >= width)
		UpperBoundX = width-1;
	else
		UpperBoundX = (int) tx;
	tx = yk_ptr[1] + h[0] - DELTA;
	if (tx >= height)
		UpperBoundY = height - 1;
	else
		UpperBoundY = (int) tx;
	
	//Perform search using lattice
	for(i = LowerBoundY; i <= UpperBoundY; i++)
		for(j = LowerBoundX; j <= UpperBoundX; j++)
		{
			
			//get index into data array
			dataPoint = N*(i*width+j);
			
			//Determine if inside search window
			k		= 1;
			s		= 0;
			dx      = j - yk_ptr[0];
			dy      = i - yk_ptr[1];
			uv[0]	= (dx*dx)/(h[0]*h[0]);
			uv[1]	= (dy*dy)/(h[0]*h[0]);
			diff	= uv[0] + uv[1];
			while((diff < offset[k-1])&&(k != kp)) // Partial Distortion Search
			{
				//Calculate distance squared of sub-space s	
				diff = 0;
				for(p = 0; p < P[k]; p++)
				{
					el    = (data[dataPoint+p+s]-yk_ptr[p+s+2])/h[k];
					diff += uv[p+s+2] = el*el;
				}
				
				//next subspace
				s += P[k];
				k++;
			}
			
			//if its inside search window perform weighted sum and count
			if(diff < offset[k-1])
			{
				
				// Initialize total weight to 1
				tw = 1;
				
				// Calculate weight factor using weight function
				// lookup tables and uv
				s = 0;
				for(k = 0; k < kp; k++)
				{
					if(kernel[k]) // not uniform kernel
					{
						// Compute u[i]
						u = 0;
						for(p = 0; p < P[k]; p++)
							u += uv[s+p];
						
						// Accumulate tw using calculated u
						// and weight function lookup table
						
						// Linear interpolate values given by
						// lookup table
						
						// Calculate x0 and x1, the points surounding
						// u
						x0 = (int)(u/increment[k]);
						x1 = x0+1;
						
						// Get y0 and y1 from the lookup table
						y0 = w[k][x0];
						y1 = w[k][x1];
						
						// Accumulate tw using linear interpolation
						tw *= (((double)(x1)*increment[k] - u)*y0+(u - (double)(x0)*increment[k])*y1)/(double)(x1*increment[k] - x0*increment[k]);
						
					}
					s += P[k];                               // next subspace
				}
				
				// Perform weighted sum using xi
				Mh_ptr[0]	+= tw*j;
				Mh_ptr[1]	+= tw*i;
				for(k = 0; k < N; k++)
					Mh_ptr[k+2] += tw*data[dataPoint+k];
				
				// Increment wsum by tw
				wsum += tw;
				
			}
			
			//done.
			
		}
		
		//done.		
		return;
		
}

/*******************************************************/
/*Optimized General Lattice Search                     */
/*******************************************************/
/*Performs search on data set for all points lying     */
/*within the search window defined using a general     */
/*kernel. Their point-wise sum and count is computed   */
/*and returned. Also the points that lie within the    */
/*window are stored into the basin of attraction stru- */
/*cture used by the optimized mean shift algorithms.   */
/*******************************************************/
/*Pre:                                                 */
/*      - Mh_ptr is a length N array of doubles        */
/*      - yk_ptr is a length N array of doubles        */
/*      - Mh_ptr is the sum of the data points found   */
/*        within search window having center yk_ptr    */
/*Post:                                                */
/*      - a search on the data set using the lattice   */
/*        has been performed, and all points found to  */
/*        lie within the search window defined using   */
/*        a general kernel are summed and counted      */
/*      - their point wise sum is pointed to by Mh_ptr */
/*        and their count is stored by wsum            */
/*      - the data points lying within h*offset of     */
/*        yk_ptr have been stored into the basin of    */
/*        attraction data structure.                   */
/*******************************************************/

void MeanShift::optGeneralLSearch(double *Mh_ptr, double *yk_ptr)
{
	
	//Declare variables
	register int	i, j, k;
	int				s, p, dataPoint, pointIndx, lN, x0, x1;
	double			diff, el, dx, dy, tw, u, y0, y1, tx;
	
	//Define lattice data dimension...
	lN	= N + 2;
	
	//Define bounds of lattice...
	
	//the lattice is a 2dimensional subspace whose
	//search window bandwidth is specified by
	//h[0]:
	tx = yk_ptr[0] - h[0] + DELTA + 0.99;
	if (tx < 0)
		LowerBoundX = 0;
	else
		LowerBoundX = (int) tx;
	tx = yk_ptr[1] - h[0] + DELTA + 0.99;
	if (tx < 0)
		LowerBoundY = 0;
	else
		LowerBoundY = (int) tx;
	tx = yk_ptr[0] + h[0] - DELTA;
	if (tx >= width)
		UpperBoundX = width-1;
	else
		UpperBoundX = (int) tx;
	tx = yk_ptr[1] + h[0] - DELTA;
	if (tx >= height)
		UpperBoundY = height - 1;
	else
		UpperBoundY = (int) tx;
	
	//Perform search using lattice
	for(i = LowerBoundY; i <= UpperBoundY; i++)
		for(j = LowerBoundX; j <= UpperBoundX; j++)
		{
			
			//get index into data array
			pointIndx	= i*width+j;
			dataPoint = N*(i*width+j);
			
			//Determine if inside search window
			k		= 1;
			s		= 0;
			dx      = j - yk_ptr[0];
			dy      = i - yk_ptr[1];
			uv[0]	= (dx*dx)/(h[0]*h[0]);
			uv[1]	= (dy*dy)/(h[0]*h[0]);
			diff	= uv[0] + uv[1];
			while((diff < offset[k-1])&&(k != kp)) // Partial Distortion Search
			{
				//Calculate distance squared of sub-space s	
				diff = 0;
				for(p = 0; p < P[k]; p++)
				{
					el    = (data[dataPoint+p+s]-yk_ptr[p+s+2])/h[k];
					diff += uv[p+s+2] = el*el;
				}
				
				//next subspace
				s += P[k];
				k++;
			}
			
			//if its inside search window perform weighted sum and count
			if(diff < offset[k-1])
			{
				
				// Initialize total weight to 1
				tw = 1;
				
				// Calculate weight factor using weight function
				// lookup tables and uv
				s = 0;
				for(k = 0; k < kp; k++)
				{
					if(kernel[k]) // not uniform kernel
					{
						// Compute u[i]
						u = 0;
						for(p = 0; p < P[k]; p++)
							u += uv[s+p];
						
						// Accumulate tw using calculated u
						// and weight function lookup table
						
						// Linear interpolate values given by
						// lookup table
						
						// Calculate x0 and x1, the points surounding
						// u
						x0 = (int)(u/increment[k]);
						x1 = x0+1;
						
						// Get y0 and y1 from the lookup table
						y0 = w[k][x0];
						y1 = w[k][x1];
						
						// Accumulate tw using linear interpolation
						tw *= (((double)(x1)*increment[k] - u)*y0+(u - (double)(x0)*increment[k])*y1)/(double)(x1*increment[k] - x0*increment[k]);
						
					}
					s += P[k];                               // next subspace
				}
				
				// Perform weighted sum using xi
				Mh_ptr[0]	+= tw*j;
				Mh_ptr[1]	+= tw*i;
				for(k = 0; k < N; k++)
					Mh_ptr[k+2] += tw*data[dataPoint+k];
				
				// Increment wsum by tw
				wsum += tw;
				
				//set basin of attraction mode table
				if(modeTable[pointIndx] == 0)
				{
					pointList[pointCount++]	= pointIndx;
					modeTable[pointIndx]	= 2;
				}
				
			}
			
			//done.
			
		}
		
		//done.		
		return;
		
}

/*@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@*/
/*@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@*/
/*@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@ END OF CLASS DEFINITION @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@*/
/*@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@*/
/*@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@*/


