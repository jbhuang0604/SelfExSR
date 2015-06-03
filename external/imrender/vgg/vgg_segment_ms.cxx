// S = vgg_segment_ms(A, h_s, h_r, min_sz[, W]);
// $Id: vgg_segment_ms.cxx,v 1.1 2007/12/10 10:59:31 ojw Exp $

#include "msImageProcessor.h"
#include <mex.h>

// Define types
#ifdef _MSC_VER
typedef __int8 int8_t;
typedef __int32 int32_t;
typedef unsigned __int8 uint8_t;
typedef unsigned __int32 uint32_t;
#else
#include <stdint.h>
#endif

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
	// Check inputs
	if (nrhs < 4 || nrhs > 5)
		mexErrMsgTxt("Unexpected number of input arguments.");
	if (nlhs < 1 || nlhs > 2)
		mexErrMsgTxt("Unexpected number of output arguments.");
	const int *dims = mxGetDimensions(prhs[0]);
	if (!mxIsUint8(prhs[0]) || mxGetNumberOfDimensions(prhs[0]) != 3 || dims[2] != 3)
		mexErrMsgTxt("A must be an HxWx3 uint8 array.");

	if (nrhs > 4) {
		if (!mxIsSingle(prhs[4]) || mxGetM(prhs[4]) != dims[0] || mxGetN(prhs[4]) != dims[1])
			mexErrMsgTxt("Edge Weight must be an HxW single array.");
	}

	// Read in the input variables
	int sigmaS = (int)mxGetScalar(prhs[1]);
	float sigmaR = (float)mxGetScalar(prhs[2]);
	int minRegion = (int)mxGetScalar(prhs[3]);
	
	msImageProcessor im_proc;

	// Read in the input image
	const uint8_t *A = (const uint8_t *)mxGetData(prhs[0]);
	uint8_t *temp_im = new uint8_t[dims[0]*dims[1]*3];
	uint8_t *B = temp_im;
	for (int h = 0; h < dims[0]; h++) {
		for (int w = 0; w < dims[0]*dims[1]; w += dims[0]) {
			*B++ = A[h+w];
			*B++ = A[h+w+dims[0]*dims[1]];
			*B++ = A[h+w+dims[0]*dims[1]*2];
		}
	}
	im_proc.DefineImage(temp_im, COLOR, dims[0], dims[1]);
	delete temp_im;
	if (im_proc.ErrorStatus == EL_ERROR)
		mexErrMsgTxt(im_proc.ErrorMessage);

	if (nrhs > 4) {
		// Define the synergistic weight map
		const float *D = (const float *)mxGetData(prhs[4]);
		float *temp_map = new float[dims[0]*dims[1]];
		float *M = temp_map;
		for (int h = 0; h < dims[0]; h++) {
			for (int w = 0; w < dims[0]*dims[1]; w += dims[0]) {
				*M++ = D[h+w];
			}
		}
		im_proc.SetWeightMap(temp_map, 1e-7);
		delete temp_map;
		if (im_proc.ErrorStatus == EL_ERROR)
			mexErrMsgTxt(im_proc.ErrorMessage);
	}

	// Segment the image
	im_proc.Segment(sigmaS, sigmaR, minRegion, HIGH_SPEEDUP);
	if (im_proc.ErrorStatus == EL_ERROR)
		mexErrMsgTxt(im_proc.ErrorMessage);

	// Get regions
	int *labels = im_proc.GetLabels();

	// Create the output image
	plhs[0] = mxCreateNumericMatrix(dims[0], dims[1], mxUINT32_CLASS, mxREAL);
	uint32_t *C = (uint32_t *)mxGetData(plhs[0]);
	for (int h = 0; h < dims[0]; h++) {
		for (int w = 0; w < dims[0]*dims[1]; w += dims[0])
				C[h+w] = (uint32_t)(*labels++) + 1;
	}

	return;
}
