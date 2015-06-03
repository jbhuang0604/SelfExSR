// B = vgg_interp2(A, X, Y[, method[, oobv]])

// $Id: vgg_interp2.cxx,v 1.4 2008/03/30 19:25:28 ojw Exp $
// Rewritten by ojw 20/9/06

#include <mex.h>
#include <math.h>
#ifdef _OPENMP
#include <omp.h>
#endif

// Matlab should include this but doesn't - why waste time and destroy the cache writing values you don't need
static inline mxArray *mxCreateUninitializedArray(int ndim, const int *dims, mxClassID classid, mxComplexity ComplexFlag)
{
    mxArray *outPtr = mxCreateNumericMatrix(0, 0, classid, ComplexFlag);
    int nBytes = mxGetElementSize(outPtr);
    for (int i = 0; i < ndim; ++i)
        nBytes *= dims[i];
    mxSetDimensions(outPtr, dims, ndim);
    mxSetData(outPtr, mxMalloc(nBytes));
    return outPtr;
}

// Define types
#ifdef _MSC_VER
typedef __int16 int16_t;
typedef unsigned __int8 uint8_t;
typedef unsigned __int16 uint16_t;
#else
#include <stdint.h>
#endif

struct interp_nearest { enum {method = 0}; };
struct interp_linear { enum {method = 1}; };
struct interp_cubic { enum {method = 2}; };

template<class Method> static inline void wrapper_func(void *B, const void *A, const float *X, const float *Y, const int num_points, const int w, const int h, const int col, const double oobv, const mxClassID out_class, const mxClassID in_class);
template<class Method, class T> static inline void wrapper_func2(void *B, const T *A, const float *X, const float *Y, const int num_points, const int w, const int h, const int col, const double oobv, const mxClassID out_class);
template<class T, class U> static inline void interp_func(interp_nearest *, U *B, const T *A, const float *X, const float *Y, const int num_points, const int w, const int h, const int col, const U oobv);
template<class T, class U> static inline void interp_func(interp_linear *, U *B, const T *A, const float *X, const float *Y, const int num_points, const int w, const int h, const int col, const U oobv);
template<class T, class U> static inline void interp_func(interp_cubic *, U *B, const T *A, const float *X, const float *Y, const int num_points, const int w, const int h, const int col, const U oobv);

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
	// Check number of arguments
	if (nrhs < 3 || nrhs > 5)
		mexErrMsgTxt("Unexpected number of input arguments.");
	if (nlhs > 1)
		mexErrMsgTxt("Unexpected number of output arguments.");

	// Check argument types are valid
	mxClassID in_class = mxGetClassID(prhs[0]);
	if (mxGetClassID(prhs[1]) != mxDOUBLE_CLASS || mxGetClassID(prhs[2]) != mxDOUBLE_CLASS)
		// mexErrMsgTxt("X and Y must be doubles");
	for (int i = 0; i < nrhs; i++) {
        if (mxIsComplex(prhs[i]))
			mexErrMsgTxt("Inputs cannot be complex.");
	}

	// Get and check array dimensions
	int num_points = mxGetNumberOfElements(prhs[1]);
	if (num_points != mxGetNumberOfElements(prhs[2]))
		mexErrMsgTxt("X and Y must have the same dimensions");
	int ndims = mxGetNumberOfDimensions(prhs[0]);
	if (ndims > 20)
		mexErrMsgTxt("A has an unsupported number of dimensions");
	int out_dims[20];
	out_dims[0] = mxGetM(prhs[1]);
	out_dims[1] = mxGetN(prhs[1]);
	const int *dims = mxGetDimensions(prhs[0]);
	out_dims[2] = 1;
	int nchannels = 1;
	for (int i = 2; i < ndims; i++) {
		out_dims[i] = dims[i];
		nchannels *= dims[i];
	}
    
	// Get pointers to input arrays
	const void *A = mxGetData(prhs[0]);
	const float *X = (const float *)mxGetData(prhs[1]);
	const float *Y = (const float *)mxGetData(prhs[2]);

	// Get the interpolation method
	int method = interp_linear::method;
	if (nrhs > 3) {
		// Read in the method string
		char buffer[32];
		mxGetString(prhs[3], buffer, 32);
		int k = 0;
		// Remove '*' from the start
		if (buffer[k] == '*')
			k++;
		// Remove 'bi' from the start
		if (buffer[k] == 'b' && buffer[k+1] == 'i')
			k += 2;
		switch (buffer[k]) {
			case 'n': 
				method = interp_nearest::method;
				break;
            case 'l': 
				method = interp_linear::method;
				break;
            case 'c': 
				method = interp_cubic::method;
				break;
			default:
				mexErrMsgTxt("Unsupported interpolation method");
				break;
		}
	}
	
	// Get the out of bounds value (oobv) and set the output class to the same class as the oobv.
	double oobv;
	mxClassID out_class;
	if (nrhs > 4) {
		// Get the value for oobv
		if (mxGetNumberOfElements(prhs[4]) != 1)
			mexErrMsgTxt("oobv must be a scalar.");
		oobv = mxGetScalar(prhs[4]);
		out_class = mxGetClassID(prhs[4]);
	} else {
		// Use the default value for oobv
		oobv = mxGetNaN();
		out_class = mxDOUBLE_CLASS;
	}
	
	// Create the output array
	plhs[0] = mxCreateUninitializedArray(ndims, out_dims, out_class, mxREAL);
	void *B = mxGetData(plhs[0]);

	// Call the wrapper function, which then sets the input type
	// Note the wrapper function is given the interpolation method struct as an input
	switch (method) {
        case interp_nearest::method:
			wrapper_func<interp_nearest>(B, A, X, Y, num_points, dims[1], dims[0], nchannels, oobv, out_class, in_class);
			break;
        case interp_linear::method:
			wrapper_func<interp_linear>(B, A, X, Y, num_points, dims[1], dims[0], nchannels, oobv, out_class, in_class);
			break;
        case interp_cubic::method:
			wrapper_func<interp_cubic>(B, A, X, Y, num_points, dims[1], dims[0], nchannels, oobv, out_class, in_class);
			break;
	}
	return;
}

template<class Method> static inline void wrapper_func(void *B, const void *A, const float *X, const float *Y, const int num_points, const int w, const int h, const int col, const double oobv, const mxClassID out_class, const mxClassID in_class)
{
	// Call the second wrapper function according to the input type
	switch (in_class) {
		case mxLOGICAL_CLASS:
			wrapper_func2<Method>(B, (const mxLogical *)A, X, Y, num_points, w, h, col, oobv, out_class);
			break;
		case mxUINT8_CLASS:
			wrapper_func2<Method>(B, (const uint8_t *)A, X, Y, num_points, w, h, col, oobv, out_class);
			break;
		case mxUINT16_CLASS:
			wrapper_func2<Method>(B, (const uint16_t *)A, X, Y, num_points, w, h, col, oobv, out_class);
			break;
		case mxINT16_CLASS:
			wrapper_func2<Method>(B, (const int16_t *)A, X, Y, num_points, w, h, col, oobv, out_class);
			break;
		case mxDOUBLE_CLASS:
			wrapper_func2<Method>(B, (const double *)A, X, Y, num_points, w, h, col, oobv, out_class);
			break;
		case mxSINGLE_CLASS:
			wrapper_func2<Method>(B, (const float *)A, X, Y, num_points, w, h, col, oobv, out_class);
			break;
		default:
			mexErrMsgTxt("A is of an unsupported type");
			break;
	}
	return;
}

template<class Method, class T> static inline void wrapper_func2(void *B, const T *A, const float *X, const float *Y, const int num_points, const int w, const int h, const int col, const double oobv, const mxClassID out_class)
{
	// Call the interp function according to the input and output type
	switch (out_class) {
		case mxLOGICAL_CLASS:
			interp_func((Method *)0, (mxLogical *)B, A, X, Y, num_points, w, h, col, (const mxLogical)oobv);
			break;
		case mxUINT8_CLASS:
			interp_func((Method *)0, (uint8_t *)B, A, X, Y, num_points, w, h, col, (const uint8_t)oobv);
			break;
		case mxUINT16_CLASS:
			interp_func((Method *)0, (uint16_t *)B, A, X, Y, num_points, w, h, col, (const uint16_t)oobv);
			break;
		case mxDOUBLE_CLASS:
			interp_func((Method *)0, (double *)B, A, X, Y, num_points, w, h, col, (const double)oobv);
			break;
		case mxSINGLE_CLASS:
			interp_func((Method *)0, (float *)B, A, X, Y, num_points, w, h, col, (const float)oobv);
			break;
		default:
			mexErrMsgTxt("Unsupported output type");
			break;
	}
	return;
}

// Function for correct rounding
// Add these to use numeric_limits class
#include <limits>
using namespace std;
template<class U, class T> static inline U saturate_cast(T val)
{
	if (numeric_limits<U>::is_integer && !numeric_limits<T>::is_integer) {
		if (numeric_limits<U>::is_signed)
			return val > 0 ? (val > (T)numeric_limits<U>::max() ? numeric_limits<U>::max() : static_cast<U>(val + 0.5)) : (val < (T)numeric_limits<U>::min() ? numeric_limits<U>::min() : static_cast<U>(val - 0.5));
		else
			return val > 0 ? (val > (T)numeric_limits<U>::max() ? numeric_limits<U>::max() : static_cast<U>(val + 0.5)) : 0;
	}
	return static_cast<U>(val);
}

// Nearest neighbour
template<class T, class U> static inline void interp_func(interp_nearest *, U *B, const T *A, const float *X, const float *Y, const int num_points, const int w, const int h, const int col, const U oobv)
{
	int end = num_points * col;
	int step = h * w;

	double dw = (double)w + 0.5;
	double dh = (double)h + 0.5;

	// For each of the interpolation points
    int i, j, k;
#pragma omp parallel for if (num_points > 1000) num_threads(omp_get_num_procs()) default(shared) private(i,j,k)
	for (i = 0; i < num_points; i++) {

		if (X[i] >= 0.5 && X[i] < dw && Y[i] >= 0.5 && Y[i] < dh) {
			// Find nearest neighbour
			k = h * int(X[i]-0.5) + int(Y[i]-0.5);
			for (j = i; j < end; j += num_points, k += step)
				B[j] = saturate_cast<U, T>(A[k]);
		} else {
			// Out of bounds
			for (j = i; j < end; j += num_points)
				B[j] = oobv;
		}
	}
	return;
}

// Linear interpolation
template<class T, class U> static inline void interp_func(interp_linear *, U *B, const T *A, const float *X, const float *Y, const int num_points, const int w, const int h, const int col, const U oobv)
{
	int end = num_points * col;
	int step = h * w;

	double dw = (double)w;
	double dh = (double)h;

	// For each of the interpolation points
    int i, j, k, x, y;
    double u, v, out;
#pragma omp parallel for if (num_points > 300) num_threads(omp_get_num_procs()) default(shared) private(i,j,k,u,v,x,y,out)
	for (i = 0; i < num_points; i++) {

		if (X[i] >= 1 && Y[i] >= 1) {
            if (X[i] < dw) {
                if (Y[i] < dh) {
                    // Linearly interpolate
                    x = (int)X[i];
                    y = (int)Y[i];
                    u = X[i] - x;
                    v = Y[i] - y;
                    k = h * (x - 1) + y - 1;
                    for (j = i; j < end; j += num_points, k += step) {
                        out = A[k] + (A[k+h] - A[k]) * u;
                        out += ((A[k+1] - out) + (A[k+h+1] - A[k+1]) * u) * v;
                        B[j] = saturate_cast<U, double>(out);
                    }
                } else if (Y[i] == dh) {
                    // The Y coordinate is on the boundary
                    // Avoid reading outside the buffer to avoid crashes
                    // Linearly interpolate along X
                    x = (int)X[i];
                    u = X[i] - x;
                    k = h * x - 1;
                    for (j = i; j < end; j += num_points, k += step)
                        B[j] = saturate_cast<U, double>(A[k] + (A[k+h] - A[k]) * u);
                } else {
                    // Out of bounds
                    for (j = i; j < end; j += num_points)
                        B[j] = oobv;
                }
            } else if (X[i] == dw) {
                if (Y[i] < dh) {
                    // The X coordinate is on the boundary
                    // Avoid reading outside the buffer to avoid crashes
                    // Linearly interpolate along Y
                    y = (int)Y[i];
                    v = Y[i] - y;
                    k = h * (w - 1) + y - 1;
                    for (j = i; j < end; j += num_points, k += step)
                        B[j] = saturate_cast<U, double>(A[k] + (A[k+1] - A[k]) * v);
                } else if (Y[i] == dh) {
                    // The X and Y coordinates are on the boundary
                    // Avoid reading outside the buffer to avoid crashes
                    // Output the last value in the array
                    k = h * w - 1;
                    for (j = i; j < end; j += num_points, k += step)
                        B[j] = saturate_cast<U, double>(A[k]);
                } else {
                    // Out of bounds
                    for (j = i; j < end; j += num_points)
                        B[j] = oobv;
                }
            } else {
                // Out of bounds
                for (j = i; j < end; j += num_points)
                    B[j] = oobv;
            }
		} else {
			// Out of bounds
			for (j = i; j < end; j += num_points)
				B[j] = oobv;
		}
	}
	return;
}

// Hermite cubic spline interpolation
template<class T, class U> static inline void interp_func(interp_cubic *, U *B, const T *A, const float *X, const float *Y, const int num_points, const int w, const int h, const int col, const U oobv)
{
	int end = num_points * col;
	int step = h * w;

	double dw = (double)w - 1;
	double dh = (double)h - 1;

	// For each of the interpolation points
    int i, j, k, m, n, x, y;
    double a, b[4], c[4], u[3], v[3];
#pragma omp parallel for if (num_points > 100) num_threads(omp_get_num_procs()) default(shared) private(a,b,c,i,j,k,m,n,u,v,x,y)
	for (i = 0; i < num_points; i++) {

		if (X[i] >= 2 && X[i] < dw && Y[i] >= 2 && Y[i] < dh) {
			// Bicubicly interpolate
			x = (int)X[i];
			y = (int)Y[i];
			u[0] = X[i] - x;
			v[0] = Y[i] - y;
			u[1] = u[0] * u[0];
			v[1] = v[0] * v[0];
			u[2] = u[1] * u[0];
			v[2] = v[1] * v[0];
			k = h * (x - 2) + y - 2;
			for(j = i; j < end; j += num_points, k += step) {
				for (m = 0, n = k; m < 4; m++, n += h) {
					c[0] = (double)A[n+0];
					c[1] = (double)A[n+1];
					c[2] = (double)A[n+2];
					c[3] = (double)A[n+3];
					a = (c[3] + c[1]) - (c[2] + c[0]);
					b[m] = v[2] * a + v[1] * ((c[0] - c[1]) - a) + v[0] * (c[2] - c[0]) + c[1];
				}
				a = (b[3] + b[1]) - (b[2] + b[0]);
				B[j] = saturate_cast<U, double>(u[2] * a + u[1] * ((b[0] - b[1]) - a) + u[0] * (b[2] - b[0]) + b[1]);
			}
		} else {
			// Out of bounds
			for (j = i; j < end; j += num_points)
				B[j] = oobv;
		}
	}
	return;
}
