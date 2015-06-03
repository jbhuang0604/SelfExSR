// S = vgg_segment_gb(A, sigma, k, min_sz, compress);
// $Id: vgg_segment_gb.cxx,v 1.2 2007/12/11 16:33:41 ojw Exp $


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

#include "image.h"
#include "misc.h"
#include "segment-image.h"

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
	// Check inputs
	if (nrhs < 4 || nrhs > 5)
		mexErrMsgTxt("Unexpected number of input arguments.");
	if (nlhs != 1)
		mexErrMsgTxt("Unexpected number of output arguments.");
	int ndims = mxGetNumberOfDimensions(prhs[0]);
	const int *dims = mxGetDimensions(prhs[0]);
	if (!mxIsUint8(prhs[0]) || ndims != 3 || dims[2] != 3)
		mexErrMsgTxt("A must be an HxWx3 uint8 array.");

	// Read in the input variables
	float sigma = (float)mxGetScalar(prhs[1]);
	float k = (float)mxGetScalar(prhs[2]);
	int min_size = (int)mxGetScalar(prhs[3]);

	// Read in the input image
	image<rgb> *input = new image<rgb>(dims[1], dims[0]);
	const uint8_t *A = (const uint8_t *)mxGetData(prhs[0]);
	uint8_t *B = (uint8_t *)imPtr(input, 0, 0);
	for (int h = 0; h < dims[0]; h++) {
		for (int w = 0; w < dims[0]*dims[1]; w += dims[0]) {
				*B++ = A[h+w];
				*B++ = A[h+w+dims[0]*dims[1]];
				*B++ = A[h+w+dims[0]*dims[1]*2];
		}
	}

	// Create the output image
	plhs[0] = mxCreateNumericMatrix(dims[0], dims[1], mxUINT32_CLASS, mxREAL);
    uint32_t *C = (uint32_t *)mxGetData(plhs[0]);

	// Segment the image
	int num_sets;
	segment_image(input, sigma, k, min_size, &num_sets, C);
	delete input;
	
	if (nrhs > 4 && mxGetScalar(prhs[4])) {
		// Compress the labelling
		num_sets++;
		uint32_t *map = new uint32_t[num_sets];
		int labels = 0;
		for (int w = 0; w < dims[0]*dims[1]; w += dims[0]) {
			for (int h = 0; h < dims[0]; h++) {
				int label = 0;
				while (label < labels) {
					if (map[label++] == C[h+w])
						goto skip_point;
				}
				label = labels++;
				if (labels > num_sets) {
					char str[30];
					sprintf(str, "At least %d sets found.\n", labels);
					mexErrMsgTxt(str);
				}
				map[label++] = C[h+w];
skip_point:
				C[h+w] = (uint32_t)label;
			}
		}
		delete map;
	}

	return;
}

