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

	The MeanShift class prototype is provided below. Its
	definition is provided in 'ms.cc'.

The theory is described in the papers:

  D. Comaniciu, P. Meer: Mean Shift: A robust approach toward feature
									 space analysis.

  C. Christoudias, B. Georgescu, P. Meer: Synergism in low level vision.

and they are is available at:
  http://www.caip.rutgers.edu/riul/research/papers/

Implemented by Chris M. Christoudias, Bogdan Georgescu
********************************************************/

#ifndef MS_H
#define MS_H

//Included needed libraries

//Include type definitions
#include	"tdef.h"

//include mean shift system used
//for function timing and system output
#include	"msSys.h"

//Include Debugging Constant
//#define	DEBUG

//Define Prompt - Prompts user on progress of Mean Shift algorithm
//#define	PROMPT

//Define Show Progress - Prompts user on percent complete of a given
//                       mean shift algorithm
//#define SHOW_PROGRESS

//Define Progress Rate - Indicates the number of convergences before
//						 checking progress
#define PROGRESS_RATE	100

// Define Macros
#define	SWAP(d_a, d_b) temp=(d_a);(d_a)=(d_b);(d_b)=temp;

// Define Structures 

 //k-Dimensional Binary Search Tree
struct tree {
  float *x;
  tree  *right;
  tree  *left;
  tree  *parent;
};

 // User Defined Weight Function
struct userWeightFunct {
   
  double			*w;
  double			halfWindow;
  int				sampleNumber;
  int				subspace;
  userWeightFunct	*next;

};

//Define class state structure
struct ClassStateStruct {
	bool	KERNEL_DEFINED;
	bool	INPUT_DEFINED;
	bool	LATTICE_DEFINED;
	bool	OUTPUT_DEFINED;
};

// Define Constants

 // Threshold
const double	EPSILON		= 0.01;			// define threshold (approx. Value of Mh at a peak or plateau)
const double	MU				= 0.05;		// define threshold required that window is near convergence
const double	TC_DIST_FACTOR	= 0.5;		// cluster search windows near convergence that are a distance
											// h[i]*TC_DIST_FACTOR of one another (transitive closure)
const double	SQ_TC_DFACTOR	= 0.0625;	// (TC_DIST_FACTOR)^2
const int		LIMIT           = 100;		// define max. # of iterations to find mode

 // Gaussian Lookup Table
const int		GAUSS_NUM_ELS   = 16;		// take 16 samples of exp(-u/2)
const double	GAUSS_LIMIT     = 2.9;		// GAUSS_LIMIT     = c
const double	GAUSS_INCREMENT = GAUSS_LIMIT*GAUSS_LIMIT/GAUSS_NUM_ELS;
											// GAUSS_INCREMENT = (c^2)/(# of samples)

 // Numerical Analysis
const double	DELTA           = 0.00001;	// used for floating point to integer conversion

//MeanShift Prototype
class MeanShift {

 public:

  /*/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\*/
  /* Class Constructor and Destructor */
  /*\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/*/

  MeanShift( void ); //Default Constructor
 ~MeanShift( void ); //Class Destructor

  /*/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\*/
  /* Creation/Initialization of Mean Shift Kernel */
  /*\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/*/

  //--\\||//--\\||//--\\||//--\\||//--\\||//--\\||//--\\||//
  //<--------------------------------------------------->|//
  //|                                                    |//
  //|	Method Name:								     |//
  //|   ============								     |//
  //|				  *  Define Kernel  *                |//
  //|                                                    |//
  //<--------------------------------------------------->|//
  //|                                                    |//
  //|	Description:								     |//
  //|	============								     |//
  //|                                                    |//
  //|   Uploads a custom kernel into the private data    |//
  //|   members of the mean shift class. This kernel is  |//
  //|   used by the mean shift class to perform mean     |//
  //|   shift.                                           |//
  //|                                                    |//
  //|   In order to create a valid kernel the following  |//
  //|   argumens must be provided this method:           |//
  //|                                                    |//
  //|   <* kernel *>                                     |//
  //|   A one dimensional array of type kernelType used  |//
  //|   to specify the kernel type (Uniform, Gaussian,   |//
  //|   or User Defined) of a given subspace of the input|//
  //|   set. Entry i of kernel correlates to the i-th    |//
  //|   subspace of the input data set.                  |//
  //|                                                    |//
  //|   <* h *>                                          |//
  //|   A one dimensional array of floating point numb-  |//
  //|   ers that are used to normalize the input data    |//
  //|   set, each bandwidth specifying the relative imp- |//
  //|   ortance of a subspace of the input data set.     |//
  //|                                                    |//
  //|   <* kp *>                                         |//
  //|   An integer that specifies the number of sub-     |//
  //|   contained by the input data set. Both P and h    |//
  //|   therefore consist of kp entries.                 |//
  //|                                                    |//
  //<--------------------------------------------------->|//
  //|                                                    |//
  //|	Usage:      								     |//
  //|   ======      								     |//
  //|		DefineKernel(kernel, h, P, kp)               |//
  //|                                                    |//
  //<--------------------------------------------------->|//
  //--\\||//--\\||//--\\||//--\\||//--\\||//--\\||//--\\||//

  void  DefineKernel(kernelType*, float*, int*, int);

  //--\\||//--\\||//--\\||//--\\||//--\\||//--\\||//--\\||//
  //<--------------------------------------------------->|//
  //|                                                    |//
  //|	Method Name:								     |//
  //|   ============								     |//
  //|               * Add Weight Function *              |//
  //|                                                    |//
  //<--------------------------------------------------->|//
  //|                                                    |//
  //|	Description:								     |//
  //|	============								     |//
  //|                                                    |//
  //|   Each subspace specified as User Defined is un-   |//
  //|   quely defined by a correlating weight function   |//
  //|   which is user defined.                           |//
  //|                                                    |//
  //|   A weight function w(u) exhibits the following    |//
  //|   properties:                                      |//
  //|                                                    |//
  //|   (1) w(u) = w(-u)                                 |//
  //|   (2) u = ((x_i-y_k)^2)/(h^2) (see docs)           |//
  //|   (3) w(u) = 0, for |u| >= halfWindow              |//
  //|                                                    |//
  //|   To add a weight function to the mean shift class |//
  //|   the following must be specified:                 |//
  //|                                                    |//
  //|   <* g() *>                                        |//
  //|   A pointer the weight function w(u) exhibiting    |//
  //|   the above properties.                            |//
  //|                                                    |//
  //|   <* halfWindow *>                                 |//
  //|   A floating point number specifying where w(u)    |//
  //|   exists (is non zero). [See Property 3 Above]     |//
  //|                                                    |//
  //|   <* sampleNumber *>                               |//
  //|   An integer used to specify the number of samples |//
  //|   used to describe w(u). Linear interpolation is   |//
  //|   used during the mean shift calculation using the |//
  //|   the samples of w(u) to determine the value of w  |//
  //|   at a location |u| < halfWindow.                  |//
  //|                                                    |//
  //|   <* subspace *>                                   |//
  //|   An integer specifying which kernel w(u) defines. |//
  //|                                                    |//
  //|   Weight functions are accounted for every time    |//
  //|   a new kernel is created.                         |//
  //|                                                    |//
  //|   If a weight function is added to non-existing    |//
  //|   subspace of the input data set  (example: the    |//
  //|   input data set containes 3 subspaces and this    |//
  //|   method is given subspace = 4) then the weight    |//
  //|   defintion will simply be ignored by the mean     |//
  //|   shift class.                                     |//
  //|                                                    |//
  //|   If a subspace is declared as kernel type User    |//
  //|   Defined and a weight function is not defined     |//
  //|   for that subspace a fatal error will occur.      |//
  //|                                                    |//
  //<--------------------------------------------------->|//
  //|                                                    |//
  //|	Usage:      								     |//
  //|   ======      								     |//
  //|   AddWeightFunction(g(u)        , halfWindow,      |//
  //|                     sampleNumber, subspace);       |//
  //|                                                    |//
  //<--------------------------------------------------->|//
  //--\\||//--\\||//--\\||//--\\||//--\\||//--\\||//--\\||//

  void AddWeightFunction(double g(double), float, int, int);

  //--\\||//--\\||//--\\||//--\\||//--\\||//--\\||//--\\||//
  //<--------------------------------------------------->|//
  //|                                                    |//
  //|	Method Name:								     |//
  //|   ============								     |//
  //|           *  Clear Weight Functions  *             |//
  //|                                                    |//
  //<--------------------------------------------------->|//
  //|                                                    |//
  //|	Description:								     |//
  //|	============								     |//
  //|                                                    |//
  //|   Removes all user defined weight functions added  |//
  //|   using method AddWeightFunction() from the        |//
  //|   private data members of the mean shift class.    |//
  //|                                                    |//
  //<--------------------------------------------------->|//
  //--\\||//--\\||//--\\||//--\\||//--\\||//--\\||//--\\||//

  void ClearWeightFunctions( void );

  /*/\/\/\/\/\/\/\/\/\/\/\/\/\/\*/
  /* Input Data Set Declaration */
  /*\/\/\/\/\/\/\/\/\/\/\/\/\/\/*/

  //--\\||//--\\||//--\\||//--\\||//--\\||//--\\||//--\\||//
  //<--------------------------------------------------->|//
  //|                                                    |//
  //|	Method Name:								     |//
  //|   ============								     |//
  //|               * Define Input *                     |//
  //|                                                    |//
  //<--------------------------------------------------->|//
  //|                                                    |//
  //|	Description:								     |//
  //|	============								     |//
  //|                                                    |//
  //|   Uploads a one dimensional array containing L     |//
  //|   N-dimensional data points into the mean shift    |//
  //|   class.                                           |//
  //|                                                    |//
  //|   An input data set is specified by:               |//
  //|                                                    |//
  //|   <* x *>                                          |//
  //|   A pointer to a floating point array.             |//
  //|                                                    |//
  //|   <* L *>                                          |//
  //|   The number of data points stored by x.           |//
  //|                                                    |//
  //|   <* N *>                                          |//
  //|   The dimension of the data points stored by x.    |//
  //|                                                    |//
  //|   The input x has the following format:            |//
  //|                                                    |//
  //|   x = <x11, x12,..., x1N,..., xL1, xL2,..., xLN>   |//
  //|                                                    |//
  //<--------------------------------------------------->|//
  //|                                                    |//
  //|	Usage:      								     |//
  //|   ======      								     |//
  //|       DefineInput(x, L, N)                         |//
  //|                                                    |//
  //<--------------------------------------------------->|//
  //--\\||//--\\||//--\\||//--\\||//--\\||//--\\||//--\\||//

  void	DefineInput(float*, int, int);

  //--\\||//--\\||//--\\||//--\\||//--\\||//--\\||//--\\||//
  //<--------------------------------------------------->|//
  //|                                                    |//
  //|	Method Name:								     |//
  //|   ============								     |//
  //|            * Define Lattice Input *                |//
  //|                                                    |//
  //<--------------------------------------------------->|//
  //|                                                    |//
  //|	Description:								     |//
  //|	============								     |//
  //|                                                    |//
  //|   Use this method to specify define an input data  |//
  //|   set defined on a lattice.                        |//
  //|                                                    |//
  //|   The arguments of this method are:                |//
  //|                                                    |//
  //|   <* x *>                                          |//
  //|   A pointer to a floating point array containing   |//
  //|   height*width, N-dimensional data points.         |//
  //|                                                    |//
  //|   <* height *>                                     |//
  //|   An integer specifying the height of the lattice. |//
  //|                                                    |//
  //|   <* width *>                                      |//
  //|   An integer specifying the width of the lattice.  |//
  //|                                                    |//
  //|   <* N *>                                          |//
  //|   The dimension of the data points stored by x.    |//
  //|                                                    |//
  //<--------------------------------------------------->|//
  //|                                                    |//
  //|	Usage:      								     |//
  //|   ======      								     |//
  //|       DefineLInput(x, height, width, N)            |//
  //|                                                    |//
  //<--------------------------------------------------->|//
  //--\\||//--\\||//--\\||//--\\||//--\\||//--\\||//--\\||//

  void	DefineLInput(float*, int, int, int);

 /*/\/\/\/\/\/\/\/\/\/\/\*/
 /*  Lattice Weight Map  */
 /*\/\/\/\/\/\/\/\/\/\/\/*/

  //--\\||//--\\||//--\\||//--\\||//--\\||//--\\||//--\\||//
  //<--------------------------------------------------->|//
  //|                                                    |//
  //|	Method Name:								     |//
  //|   ============								     |//
  //|			  * Set Lattice Weight Map *             |//
  //|                                                    |//
  //<--------------------------------------------------->|//
  //|                                                    |//
  //|	Description:								     |//
  //|	============								     |//
  //|                                                    |//
  //|   Uploads weight map specifying for each data      |//
  //|   point a value used to weight the uniform kernel  |//
  //|   when computing mean shift.                       |//
  //|                                                    |//
  //|   The arguments to this method are:                |//
  //|                                                    |//
  //|   <* weightMap *>                                  |//
  //|   A floating point array of size L specifying for  |//
  //|   each data point a weight.                        |//
  //|                                                    |//
  //|   Note: DefineLInput must be called prior to call- |//
  //|         ing this method. DefineLInput is used to   |//
  //|         define the dimensions of the input data    |//
  //|         set.                                       |//
  //|                                                    |//
  //|                                                    |//
  //|  The weight map is used to weight the uniform      |//
  //|  kernel used to computing meanshift on a data      |//
  //|  point situated on a lattice. Alternatively, a     |//
  //|  weight function may defined, however, if speed    |//
  //|  is an issue, the lattice may be exploited to      |//
  //|  result in a faster implementation of a weighted   |//
  //|  kernel.                                           |//
  //|                                                    |//
  //<--------------------------------------------------->|//
  //|                                                    |//
  //|	Usage:      								     |//
  //|   ======      								     |//
  //|		SetLatticeWeightMap(weightMap)               |//
  //|                                                    |//
  //<--------------------------------------------------->|//
  //--\\||//--\\||//--\\||//--\\||//--\\||//--\\||//--\\||//

  void SetLatticeWeightMap(float*);

  //--\\||//--\\||//--\\||//--\\||//--\\||//--\\||//--\\||//
  //<--------------------------------------------------->|//
  //|                                                    |//
  //|	Method Name:								     |//
  //|   ============								     |//
  //|			* Remove Lattice Weight Map *            |//
  //|                                                    |//
  //<--------------------------------------------------->|//
  //|                                                    |//
  //|	Description:								     |//
  //|	============								     |//
  //|                                                    |//
  //|   Removes lattice weight map. An error is NOT      |//
  //|   flagged if a weight map was not defined prior    |//
  //|   to calling this method.                          |//
  //|                                                    |//
  //<--------------------------------------------------->|//
  //|                                                    |//
  //|	Usage:      								     |//
  //|   ======      								     |//
  //|		RemoveLatticeWeightMap(weightMap)            |//
  //|                                                    |//
  //<--------------------------------------------------->|//
  //--\\||//--\\||//--\\||//--\\||//--\\||//--\\||//--\\||//

  void RemoveLatticeWeightMap(void);

  /*/\/\/\/\/\/\/\/\/\/\/\/\*/
  /* Mean Shift Operations  */
  /*\/\/\/\/\/\/\/\/\/\/\/\/*/

  //--\\||//--\\||//--\\||//--\\||//--\\||//--\\||//--\\||//
  //<--------------------------------------------------->|//
  //|                                                    |//
  //|	Method Name:								     |//
  //|   ============								     |//
  //|               *  Mean Shift Vector  *              |//
  //|                                                    |//
  //<--------------------------------------------------->|//
  //|                                                    |//
  //|	Description:								     |//
  //|	============								     |//
  //|                                                    |//
  //|   If a kernel is created and input is uploaded,    |//
  //|   this method calcualtes the mean shift vector,    |//
  //|   Mh, at specific data point yk.                   |//
  //|                                                    |//
  //|   The arguments of this method are:                |//
  //|                                                    |//
  //|   <* Mh *>                                         |//
  //|   An array of N doubles storing the N dimensional  |//
  //|   mean shift vector.                               |//
  //|                                                    |//
  //|   <* yk *>                                         |//
  //|   An array of N doubles storing the N dimensional  |//
  //|   data point where the mean shift vector is to be  |//
  //|   calculate.                                       |//
  //|                                                    |//
  //<--------------------------------------------------->|//
  //|                                                    |//
  //|	Usage:      								     |//
  //|   ======      								     |//
  //|       msVector(Mh, yk)                             |//
  //|                                                    |//
  //<--------------------------------------------------->|//
  //--\\||//--\\||//--\\||//--\\||//--\\||//--\\||//--\\||//

  void	msVector(double*, double*);

  //--\\||//--\\||//--\\||//--\\||//--\\||//--\\||//--\\||//
  //<--------------------------------------------------->|//
  //|                                                    |//
  //|	Method Name:								     |//
  //|   ============								     |//
  //|          *  Lattice Mean Shift Vector  *           |//
  //|                                                    |//
  //<--------------------------------------------------->|//
  //|                                                    |//
  //|	Description:								     |//
  //|	============								     |//
  //|                                                    |//
  //|   If a kernel is created and input is uploaded,    |//
  //|   this method calcualtes the mean shift vector,    |//
  //|   Mh, at specific data point yk, assuming that the |//
  //|   data set exhists on a height x width two dim-    |//
  //|   ensional lattice.                                |//
  //|                                                    |//
  //|   The arguments of this method are:                |//
  //|                                                    |//
  //|   <* Mh *>                                         |//
  //|   An array of N doubles storing the N dimensional  |//
  //|   mean shift vector.                               |//
  //|                                                    |//
  //|   <* yk *>                                         |//
  //|   An array of N doubles storing the N dimensional  |//
  //|   data point where the mean shift vector is to be  |//
  //|   calculate.                                       |//
  //|                                                    |//
  //|   The height and width of the lattice must be      |//
  //|   specified using DefineLattice() method. If this  |//
  //|   is not performed prior to calling this method a  |//
  //|   fatal error will be flagged.                     |//
  //|                                                    |//
  //<--------------------------------------------------->|//
  //|                                                    |//
  //|	Usage:      								     |//
  //|   ======      								     |//
  //|        latticeMSVector(Mh, yk)                     |//
  //|                                                    |//
  //<--------------------------------------------------->|//
  //--\\||//--\\||//--\\||//--\\||//--\\||//--\\||//--\\||//

  void	latticeMSVector(double*, double*);

  //--\\||//--\\||//--\\||//--\\||//--\\||//--\\||//--\\||//
  //<--------------------------------------------------->|//
  //|                                                    |//
  //|	Method Name:								     |//
  //|   ============								     |//
  //|                 *  Find Mode  *                    |//
  //|                                                    |//
  //<--------------------------------------------------->|//
  //|                                                    |//
  //|	Description:								     |//
  //|	============								     |//
  //|                                                    |//
  //|   If a kernel is created and input is uploaded,    |//
  //|   this method calcualtes the mode of a specified   |//
  //|   data point yk.                                   |//
  //|                                                    |//
  //|   The arguments of this method are:                |//
  //|                                                    |//
  //|   <* mode *>                                       |//
  //|   An array of N doubles storing the N dimensional  |//
  //|   mode of yk.                                      |//
  //|                                                    |//
  //|   <* yk *>                                         |//
  //|   An array of N doubles storing the N dimensional  |//
  //|   data point where the mean shift vector is to be  |//
  //|   calculate.                                       |//
  //|                                                    |//
  //<--------------------------------------------------->|//
  //|                                                    |//
  //|	Usage:      								     |//
  //|   ======      								     |//
  //|       FindMode(mode, yk)                           |//
  //|                                                    |//
  //<--------------------------------------------------->|//
  //--\\||//--\\||//--\\||//--\\||//--\\||//--\\||//--\\||//

  void FindMode(double*, double*);

  //--\\||//--\\||//--\\||//--\\||//--\\||//--\\||//--\\||//
  //<--------------------------------------------------->|//
  //|                                                    |//
  //|	Method Name:								     |//
  //|   ============								     |//
  //|              *  Find Lattice Mode  *               |//
  //|                                                    |//
  //<--------------------------------------------------->|//
  //|                                                    |//
  //|	Description:								     |//
  //|	============								     |//
  //|                                                    |//
  //|   If a kernel is created and input is uploaded,    |//
  //|   this method calcualtes the mode of a specified   |//
  //|   data point yk, assuming that the data set        |//
  //|   exhists on a height x width two dimensional      |//
  //|   lattice.                                         |//
  //|                                                    |//
  //|   The arguments of this method are:                |//
  //|                                                    |//
  //|   <* mode *>                                       |//
  //|   An array of N doubles storing the N dimensional  |//
  //|   mode of yk.                                      |//
  //|                                                    |//
  //|   <* yk *>                                         |//
  //|   An array of N doubles storing the N dimensional  |//
  //|   data point where the mean shift vector is to be  |//
  //|   calculate.                                       |//
  //|                                                    |//
  //|   The height and width of the lattice must be      |//
  //|   specified using DefineLattice() method. If this  |//
  //|   is not performed prior to calling this method a  |//
  //|   fatal error will be flagged.                     |//
  //|                                                    |//
  //<--------------------------------------------------->|//
  //|                                                    |//
  //|	Usage:      								     |//
  //|   ======      								     |//
  //|       FindLMode(mode, yk)                          |//
  //|                                                    |//
  //<--------------------------------------------------->|//
  //--\\||//--\\||//--\\||//--\\||//--\\||//--\\||//--\\||//

  void FindLMode(double*, double*);

  /*/\/\/\/\/\/\/\/\/\/\/\/\/\*/
  /*  Error Handler Mechanism */
  /*/\/\/\/\/\/\/\/\/\/\/\/\/\*/

  //--\\||//--\\||//--\\||//--\\||//--\\||//--\\||//--\\||//
  //<--------------------------------------------------->|//
  //|                                                    |//
  //|	Description:								     |//
  //|	============								     |//
  //|                                                    |//
  //|   ErrorMessage is an error message that is set by  |//
  //|   a mean shift library class when an error occurs. |//
  //|                                                    |//
  //<--------------------------------------------------->|//
  //--\\||//--\\||//--\\||//--\\||//--\\||//--\\||//--\\||//

  char			*ErrorMessage;

  //--\\||//--\\||//--\\||//--\\||//--\\||//--\\||//--\\||//
  //<--------------------------------------------------->|//
  //|                                                    |//
  //|	Description:								     |//
  //|	============								     |//
  //|                                                    |//
  //|   ErrorStatus indicates if an error has occured as |//
  //|   a result of improper use of a mean shift library |//
  //|   class method or because of insufficient resour-  |//
  //|   ces. ErrorStatus is set to EL_ERROR (ErrorStatus |//
  //|   = 1) if an error has occured. If no error occur- |//
  //|   ed when calling a particular method ErrorStatus  |//
  //|   is set to EL_OKAY (ErrorStatus = 0).             |//
  //|                                                    |//
  //<--------------------------------------------------->|//
  //--\\||//--\\||//--\\||//--\\||//--\\||//--\\||//--\\||//

  ErrorLevel	ErrorStatus;

 protected:

  //==========================
  // *** Protected Methods ***
  //==========================

     /*/\/\/\/\/\/\/\/\/\/\/\/\/\/\*/
     /* Mean Shift: Using kd-Tree  */
     /*\/\/\/\/\/\/\/\/\/\/\/\/\/\/*/

   /////////////////////////////////////////
   // <<*>> Usage: MSVector(Mh, yk) <<*>> //
   /////////////////////////////////////////

   void MSVector      (double*, double*);               // Computes the mean shift vector at a specified
                                                        // window location yk in the data set x given
				                                        // the vector yk
     /*/\/\/\/\/\/\/\/\/\/\/\/\/\/\*/
     /*  Mean Shift: Using Lattice */
     /*\/\/\/\/\/\/\/\/\/\/\/\/\/\/*/

   ////////////////////////////////////////////////
   // <<*>> Usage: LatticeMSVector(Mh, yk) <<*>> //
   ////////////////////////////////////////////////

   void	LatticeMSVector     (double*, double*);			// Uses the lattice defined by DefineLattice to compute the
														// mean shift vector at a specified window location yk

   ///////////////////////////////////////////////////
   // <<*>> Usage: OptLatticeMSVector(Mh, yk) <<*>> //
   ///////////////////////////////////////////////////

   void	OptLatticeMSVector     (double*, double*);      // Uses the lattice defined by DefineLattice to compute the
														// mean shift vector at a specified window location yk using
														// the basin of attraction optimization for better performace
														// during mean shift filtering - used by a derived class

    /*/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\*/
	/* Kernel-Input Data Consistency  */
    /*\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/*/

   /////////////////////////////////////////////////
   // <<*>> Usage: classConsistencyCheck(N) <<*>> //
   /////////////////////////////////////////////////

   void classConsistencyCheck(int, bool);				// checks to see that a kernel is created and input defined, as
														// well as the specified dimension of the data set matches that of
														// the kernel, if not an error is flagged and the program is halted

   	 /*/\/\/\/\/\/\/\/\/\/\/\*/
     /* Class Error Handler  */
	 /*\/\/\/\/\/\/\/\/\/\/\/*/

   /////////////////////////////////////////////////////
   // <<*>> Usage: ErrorHandler(                <<*>> //
   //			className, functName, errMessage)     //      
   /////////////////////////////////////////////////////                           

   void ErrorHandler(char*, char*, char*);				// flags an error and halts the system


  //===============================
  // *** Protected Data Members ***
  //===============================

   //##########################################
   //#########   MEAN SHIFT SYSTEM   ##########
   //##########################################

	msSystem		msSys;								// used for function timing and system output

   //##########################################
   //######### INPUT DATA PARAMETERS ##########
   //##########################################

	int				L, N, kp, *P;						// length, dimension, subspace number, and subspace dimensions


   //##########################################
   //######### INPUT DATA STORAGE    ##########
   //##########################################

	////////Linear Storage (used by lattice and bst)////////
	float			*data;								// memory allocated for data points stored by tree nodes
														// when used by the lattice data structure data does not store
														// the lattice information; format of data:
														// data = <x11, x12, ..., x1N,...,xL1, xL2, ..., xLN>
														// in the case of the lattice the i in data(i,j) corresponds

   //##########################################
   //######## LATTICE DATA STRUCTURE ##########
   //##########################################

	////////Lattice Data Structure////////
	int				height, width;						// Height and width of lattice

   //##########################################
   //######### KERNEL DATA STRUCTURE ##########
   //##########################################

	float			*h;									// bandwidth vector

	float			*offset;							// defines bandwidth offset caused by the use of a Gaussian kernel
                                                        // (for example)

   //##########################################
   //#########  BASIN OF ATTRACTION  ##########
   //##########################################

	unsigned char	*modeTable;							// Assigns a marking to each data point specifying whether
														// or not it has been assigned a mode. These labels are:
														// modeTable[i] = 0 - data point i is not associated with a mode
														// modeTable[i] = 1 - data point i is associated with a mode
														// modeTable[i] = 2 - data point i is associated with a mode
														//                    however its mode is yet to be determined

	int				*pointList;							// a list of data points that due to basin of attraction will
														// converge to the same mode as the mode that mean shift is
														// currently being applied to

	int				pointCount;							// the number of points stored by the point list


   //##########################################
   //#########  WEIGHT MAP USED      ##########
   //#########  WHEN COMPUTING MEAN  ##########
   //#########  SHIFT ON A LATTICE   ##########
   //##########################################

	float			*weightMap;							// weight map that may be used to weight the kernel
														// upon performing mean shift on a lattice

	bool			weightMapDefined;					// used to indicate if a lattice weight map has been
														// defined

   //##########################################
   //#######        CLASS STATE        ########
   //##########################################

	ClassStateStruct	class_state;					//specifies the state of the class(i.e if data has been loaded into 
														//the class, if a kernel has been defined, etc.)

 private:

  //========================
  // *** Private Methods ***
  //========================

	 /*/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\*/
	 /* Kernel Creation/Manipulation */
	 /*\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/*/

   void	generateLookupTable ( void );                  // Generates Weight Function Lookup Table

   void	DestroyKernel       ( void );                  // Destroys mean shift kernel, re-initializes kernel


	 /*/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\*/
	 /* Input Data Initialization/Destruction  */
	 /*\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/*/

   void	CreateBST		( void );						// Upload input into a kd-BST

   void	InitializeInput	(float*);						// Allocates memory for and initializes the input data structure

   void	ResetInput		( void );						// de-allocate memory for and re-initialize input data structure
														// and mode structure

     /*/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\*/
     /* k-dimensional Binary Search Tree */
     /*\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/*/

	////////Data Search Tree/////////
   tree  *BuildKDTree (tree*,  int, int, tree* );       // Builds a kd tree given a subset of points initialized
                                                        // at depth 0 (dimension 0) (for Tree Structure)

   void  QuickMedian (tree*, int, int, int );           // Finds the median tree in a forest of trees using
                                                        // dimension d, placing the median tree in the array of tree
                                                        // nodes at L/2, in which all trees to the left of the median tree
                                                        // have values less than that of the median tree in dimension d
                                                        // and all trees having values greater than that of the median tree
                                                        // in dimension d are placd to the right of this tree -
                                                        // This algorithm is used by BuildKDTree to construct a balanced tree

     /*/\/\/\/\/\/\/\/\/\/\/\/\/\/\*/
     /* Mean Shift: Using kd-Tree  */
     /*\/\/\/\/\/\/\/\/\/\/\/\/\/\/*/

   void uniformSearch (tree*, int, double*, double*);     // uses kdbst to perform range search on input data,
                                                        // computing the weighted sum of these points using
                                                        // a uniform kernel and storing the result into Mh
                                                        // (called by uniformMSVector)

   void generalSearch (tree*, int, double*, double*);     // uses kdbst to perform range search on input data,
                                                        // computing the weighted sum of these points using
                                                        // a general kernel and storing the result into Mh
                                                        // (called by generalMSVector)

     /*/\/\/\/\/\/\/\/\/\/\/\/\/\/\*/
     /*  Mean Shift: Using Lattice */
     /*\/\/\/\/\/\/\/\/\/\/\/\/\/\/*/

   void	uniformLSearch	 (double *, double *);			// given a center location and mean shift vector, a lattice
														// search is performed to compute the mean shift vector
														// using a uniform kernel

   void optUniformLSearch(double *, double *);			// given a center location and mean shift vector, a lattice
														// search is performed to compute the mean shift vector
														// using a uniform kernel and the basin of attraction
														// optimization for better performance

   void generalLSearch	 (double *, double *);			// given a center location and mean shift vector, a lattice
														// search is performed to compute the mean shift vector
														// using a general kernel

   void optGeneralLSearch(double *, double *);			// given a center location and mean shift vector, a lattice
														// search is performed to compute the mean shift vector
														// using a general kernel and the basin of attraction
														// optimization for better performance


  //=============================
  // *** Private Data Members ***
  //=============================

   //##########################################
   //######### KERNEL DATA STRUCTURE ##########
   //##########################################

	kernelType		*kernel;							// kernel types for each subspace S[i]

	double			**w;								// weight function lookup table

	double			*increment;							// increment used by weight hashing function

	bool			uniformKernel;						// flag used to indicate if the kernel is uniform or not
	
	userWeightFunct	*head, *cur;						// user defined weight function linked list
   

   //##########################################
   //######### INPUT DATA STORAGE    ##########
   //##########################################

	////////Range Searching on General Input Data Set////////
	tree			*root;								// root of kdBST used to store input

	tree			*forest;							// memory allocated for tree nodes

	float			*range;								// range vector used to perform range search on kd tree, indexed
														// by dimension of input - format:
                                                        // range = {Lower_Limit_1, Upper_Limit_1, ..., Lower_Limit_N, Upper_Limit_N}

   //##########################################
   //######### MEAN SHIFT PROCESSING ##########
   //######### DATA STRUCTURES       ##########
   //##########################################

	double			*uv;								// stores normalized distance vector between
														// yk and xi

	double			wsum;								// sum of weights calculated at data points within the sphere

   //##########################################
   //######### LATTICE DATA STRUCTURE #########
   //##########################################

	////////Lattice Data Structure////////
	int				LowerBoundX, UpperBoundX;			// Upper and lower bounds for lattice search window
														// in the x dimension

	int				LowerBoundY, UpperBoundY;			// Upper and lower bounds for lattice search window
														// in the y dimension

};

#endif
