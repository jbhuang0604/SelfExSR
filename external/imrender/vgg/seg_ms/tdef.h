/*******************************************************

                 Mean Shift Analysis Library
	=============================================

	The mean shift library is a collection of routines
	that use the mean shift algorithm. Using this algorithm,
	the necessary output will be generated needed
	to analyze a given input set of data.

  Type Defintions:
  ===============

	This header file contains the type defintions and
	enumerations shared among the various classes of the mean
	shift library.

The theory is described in the papers:

  D. Comaniciu, P. Meer: Mean Shift: A robust approach toward feature
									 space analysis.

  C. Christoudias, B. Georgescu, P. Meer: Synergism in low level vision.

and they are is available at:
  http://www.caip.rutgers.edu/riul/research/papers/

Implemented by Chris M. Christoudias, Bogdan Georgescu
********************************************************/

#ifndef TDEF_H
#define TDEF_H

/*/\/\/\/\/\/\/\/\/\/\/\*/
/* Define Enumerations  */
/*\/\/\/\/\/\/\/\/\/\/\/*/

//Kernel
enum kernelType		{Uniform, Gaussian, UserDefined};

// kd-Tree
enum childType		{LEFT, RIGHT};

// Speed Up Level
enum SpeedUpLevel	{NO_SPEEDUP, MED_SPEEDUP, HIGH_SPEEDUP};

// Error Handler
enum ErrorLevel		{EL_OKAY, EL_ERROR, EL_HALT};
enum ErrorType		{NONFATAL, FATAL};

#endif
