//  I = find_interactions(C, O, dist)

// $Id: find_interactions.cxx,v 1.1 2007/12/07 11:27:51 ojw Exp $

#include "mex.h"
#define MAX_MEAN_INTERACTIONS 40

// Define types
#ifdef _MSC_VER
typedef unsigned __int32 uint32_t;
#else
#include <stdint.h>
#endif

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
	// Check number of arguments
	if (nrhs != 2)
		mexErrMsgTxt("2 input arguments expected.");
	if (nlhs != 1)
		mexErrMsgTxt("1 output argument expected.");

	// Check argument types and sizes are valid
	for (int i = 0; i < 1; i++) {
        if (mxIsComplex(prhs[i]))
			mexErrMsgTxt("Inputs cannot be complex.");
	}
	if (!mxIsDouble(prhs[0]))
		mexErrMsgTxt("C must be doubles.");
	if (mxGetN(prhs[0]) != 3)
		mexErrMsgTxt("C has unexpected dimensions.");
	int length = mxGetM(prhs[0]);
	if (mxGetNumberOfElements(prhs[1]) != 1)
		mexErrMsgTxt("dist must be a scalar.");
	double dist = mxGetScalar(prhs[1]);

	// Get pointers to input arrays
	const double *X = (const double *)mxGetData(prhs[0]) - 1;
	const double *Y = X + length;
	const double *Z = Y + length;

	// Create an output buffer
	plhs[0] = mxCreateNumericMatrix(2, MAX_MEAN_INTERACTIONS*length, mxUINT32_CLASS, mxREAL);
	uint32_t *Istart = (uint32_t *)mxGetData(plhs[0]);
	uint32_t *I = Istart;

	// Find interactions
	for (int a = 1; a < length; a++) {
		double x = X[a] + dist;
		double yl = Y[a] - dist;
		double yh = Y[a] + dist;
		for (int b = a+1; b <= length; b++) {
			if (X[b] > x)
				break;
			if (Y[b] < yl || Y[b] > yh) 
				continue;
			if (I >= &Istart[MAX_MEAN_INTERACTIONS*length])
				mexErrMsgTxt("Temporary buffer not big enough!!!");
			if (Z[a] < Z[b]) {
				*I++ = a;
				*I++ = b;
			} else {
				*I++ = b;
				*I++ = a;
			}
		}
	}
	// Shrink the buffer to the correct size
	int n = I - Istart;
	mxSetData(plhs[0], mxRealloc(Istart, n*sizeof(uint32_t)));
	mxSetN(plhs[0], n/2);
	return;
}