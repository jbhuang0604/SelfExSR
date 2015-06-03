// [U P PI T TI] = ibr_gen_cliques(IA, VA, V, Kocc, method)

// $Id: ibr_gen_cliques.cxx,v 1.2 2008/07/16 20:24:01 ojw Exp $

#include <mex.h>
#include <math.h>
#include <cmath>

// Define types
#ifdef _MSC_VER
typedef __int32 int32_t;
typedef unsigned __int32 uint32_t;
#else
#include <stdint.h>
#endif

struct ssd { enum {method = 1}; };
struct sad { enum {method = 2}; };
struct rssd { enum {method = 3}; };

template<class Q> static inline void wrapper_func(const Q *IA, const mxLogical *VA, const mxLogical *V, const int *dims, double Kocc, int method, mxClassID out_class, mxArray *plhs[]);
template<class Q, class R> static inline void wrapper_func2(const Q *IA, const mxLogical *VA, const mxLogical *V, const int *dims, R Kocc, int method, mxClassID out_class, mxArray *plhs[]);
template<class Method, class Q, class R> static inline void gen_cliques(const Q *IA, const mxLogical *VA, const mxLogical *V, const int *dims, R Kocc, Method *, mxClassID out_class, mxArray *plhs[]);

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
	// Check number of arguments
	if (nrhs != 5)
		mexErrMsgTxt("Unexpected number of input arguments.");
	if (nlhs > 5)
		mexErrMsgTxt("Unexpected number of output arguments.");

	// Check argument types are valid
	mxClassID in_class = mxGetClassID(prhs[0]);
	mxClassID out_class = mxGetClassID(prhs[3]);
	if (mxGetClassID(prhs[1]) != mxLOGICAL_CLASS || mxGetClassID(prhs[2]) != mxLOGICAL_CLASS)
		mexErrMsgTxt("VA and A must be logical arrays");
	for (int i = 0; i < nrhs; i++) {
        if (mxIsComplex(prhs[i]))
			mexErrMsgTxt("Inputs cannot be complex.");
	}

	// Get and check array dimensions
	if (mxGetNumberOfDimensions(prhs[0]) != 3)
		mexErrMsgTxt("IA has an unexpected number of dimensions");
	const int *dims = mxGetDimensions(prhs[0]);
	if ((dims[0]*dims[2]) != mxGetNumberOfElements(prhs[1]) || (dims[0]*dims[2]) != mxGetNumberOfElements(prhs[2]))
		mexErrMsgTxt("VA and/or V have an unexpected number of dimensions");

	// Get input values
	const void *IA = mxGetData(prhs[0]);
	const mxLogical *VA = (const mxLogical *)mxGetData(prhs[1]);
	const mxLogical *V = (const mxLogical *)mxGetData(prhs[2]);
	double Kocc = mxGetScalar(prhs[3]);
	int method = (int)mxGetScalar(prhs[4]);

	// Call the wrapper function according to the input type
	switch (in_class) {
		case mxDOUBLE_CLASS:
			wrapper_func((const double *)IA, VA, V, dims, Kocc, method, out_class, plhs);
			break;
		case mxSINGLE_CLASS:
			wrapper_func((const float *)IA, VA, V, dims, Kocc, method, out_class, plhs);
			break;
		default:
			mexErrMsgTxt("IA is of an unsupported type");
			break;
	}
	return;
}

template<class Q> static inline void wrapper_func(const Q *IA, const mxLogical *VA, const mxLogical *V, const int *dims, double Kocc, int method, mxClassID out_class, mxArray *plhs[])
{
	// Call the interp function according to the input and output type
	switch (out_class) {
		case mxINT32_CLASS:
			wrapper_func2(IA, VA, V, dims, (int32_t)Kocc, method, out_class, plhs);
			break;
		case mxUINT32_CLASS:
			wrapper_func2(IA, VA, V, dims, (uint32_t)Kocc, method, out_class, plhs);
			break;
		case mxDOUBLE_CLASS:
			wrapper_func2(IA, VA, V, dims, (double)Kocc, method, out_class, plhs);
			break;
		case mxSINGLE_CLASS:
			wrapper_func2(IA, VA, V, dims, (float)Kocc, method, out_class, plhs);
			break;
		default:
			mexErrMsgTxt("Unsupported output type");
			break;
	}
	return;
}

template<class Q, class R> static inline void wrapper_func2(const Q *IA, const mxLogical *VA, const mxLogical *V, const int *dims, R Kocc, int method, mxClassID out_class, mxArray *plhs[])
{
	// Call the interp function according to the input and output type
	switch (method) {
        case ssd::method:
			gen_cliques(IA, VA, V, dims, Kocc, (ssd *)0, out_class, plhs);
			break;
		case sad::method:
			gen_cliques(IA, VA, V, dims, Kocc, (sad *)0, out_class, plhs);
			break;
		case rssd::method:
			gen_cliques(IA, VA, V, dims, Kocc, (rssd *)0, out_class, plhs);
			break;
		default:
			mexErrMsgTxt("Unsupported data cost method");
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

template<class T> static inline void calc_mean(T *M, const T *IA, const int *dims, int a) 
{
	// Sum the values
	for (int c = 0; c < dims[1]; c++){
		M[c] = 0; // Clear the mean
		for (int b = 0; b < dims[2]; b++)
			M[c] += IA[a+c*dims[0]+b*dims[0]*dims[1]];
		M[c] /= dims[2]; // Divide by N
	}
}

template<class T> static inline void calc_mean(T *M, const T *IA, const int *dims, int a, mxLogical *V) 
{
	int num = 0;

	// Clear the mean
	for (int c = 0; c < dims[1]; c++)
		M[c] = 0;

	// Sum the values
	for (int b = 0; b < dims[2]; b++) {
		if (V[b]) {
			for (int c = 0; c < dims[1]; c++)
				M[c] += IA[a+c*dims[0]+b*dims[0]*dims[1]];
			num++;
		}
	}

	if (num) {
		// Divide by N
		for (int c = 0; c < dims[1]; c++)
			M[c] /= num;
	} else {
		for (int c = 0; c < dims[1]; c++)
			M[c] = -1000;
	}
}

template<class Q, class R>static inline R calc_cost(ssd *, const Q *M, const Q *IA, const int *dims, int a, int b, R Kocc) {
	Q tot = 0;
	for (int c = 0; c < dims[1]; c++) {
		Q val = M[c] - IA[a+c*dims[0]+b*dims[0]*dims[1]];
		tot += val * val;
	}
	R out = saturate_cast<R, Q>(tot);
	out = out > Kocc ? Kocc : out;
	return out;
}
template<class Q, class R>static inline R calc_cost(sad *, const Q *M, const Q *IA, const int *dims, int a, int b, R Kocc) {
	Q tot = 0;
	for (int c = 0; c < dims[1]; c++) {
		Q val = M[c] - IA[a+c*dims[0]+b*dims[0]*dims[1]];
		tot += abs(val);
	}
	R out = saturate_cast<R, Q>(tot);
	out = out > Kocc ? Kocc : out;
	return out;
}
template<class Q, class R>static inline R calc_cost(rssd *, const Q *M, const Q *IA, const int *dims, int a, int b, R Kocc) {
	Q tot = 0;
	for (int c = 0; c < dims[1]; c++) {
		Q val = M[c] - IA[a+c*dims[0]+b*dims[0]*dims[1]];
		tot += val * val;
	}
	tot = sqrt(tot);
	R out = saturate_cast<R, double>(tot);
	out = out > Kocc ? Kocc : out;
	return out;
}
template<class Q, class R>static inline R calc_cost(ssd *, const Q *IA, const int *dims, int a, R Kocc) {
	Q tot = 0;
	for (int c = 0; c < dims[1]; c++) {
		Q val = IA[a+c*dims[0]] - IA[a+c*dims[0]+dims[0]*dims[1]];
		tot += val * val;
	}
	R out = saturate_cast<R, Q>(tot);
	out = out > Kocc ? Kocc : out;
	return out;
}
template<class Q, class R>static inline R calc_cost(sad *, const Q *IA, const int *dims, int a, R Kocc) {
	Q tot = 0;
	for (int c = 0; c < dims[1]; c++) {
		Q val = IA[a+c*dims[0]] - IA[a+c*dims[0]+dims[0]*dims[1]];
		tot += abs(val);
	}
	R out = saturate_cast<R, Q>(tot);
	out = out > Kocc ? Kocc : out;
	return out;
}
template<class Q, class R>static inline R calc_cost(rssd *, const Q *IA, const int *dims, int a, R Kocc) {
	Q tot = 0;
	for (int c = 0; c < dims[1]; c++) {
		Q val = IA[a+c*dims[0]] - IA[a+c*dims[0]+dims[0]*dims[1]];
		tot += val * val;
	}
	tot = sqrt(tot);
	R out = saturate_cast<R, double>(tot);
	out = out > Kocc ? Kocc : out;
	return out;
}

// Generate cliques
template<class Method, class Q, class R> static inline void gen_cliques(const Q *IA, const mxLogical *VA, const mxLogical *V, const int *dims, R Kocc, Method *, mxClassID out_class, mxArray *plhs[])
{
	// Count up numbers of cliques of each size
	int npair = 0;
	int ntriple = 0;
	for (int c = 0; c < dims[0]; c++) {
		int num_occl = dims[2];
		for (int b = 0; b < dims[2]; b++)
			num_occl -= VA[c+b*dims[0]];
		switch (num_occl) {
			case 0:
				break;
			case 1:
				npair++;
				break;
			case 2:
				ntriple++;
				break;
			default:
				npair += num_occl;
				break;
		}
	}

	// Create the output matrices
	plhs[0] = mxCreateNumericMatrix(2, (dims[0]+2*dims[0]*dims[2])/2, out_class, mxREAL);
	plhs[1] = mxCreateNumericMatrix(4, npair, out_class, mxREAL);
	plhs[2] = mxCreateNumericMatrix(2, npair, mxUINT32_CLASS, mxREAL);
	plhs[3] = mxCreateNumericMatrix(8, ntriple, out_class, mxREAL);
	plhs[4] = mxCreateNumericMatrix(3, ntriple, mxUINT32_CLASS, mxREAL);
	R *U = (R *)mxGetData(plhs[0]);
	R *P = (R *)mxGetData(plhs[1]);
	uint32_t *PI = (uint32_t *)mxGetData(plhs[2]);
	R *T = (R *)mxGetData(plhs[3]);
	uint32_t *TI = (uint32_t *)mxGetData(plhs[4]);

	// Go through each value
	Q* M = (Q *)mxMalloc(sizeof(Q)*dims[1]);
	mxLogical *vis = (mxLogical *)mxMalloc(sizeof(mxLogical)*dims[2]);
	npair = 0;
	ntriple = 0;
	int a = 0;
	if (1 || dims[2] != 2) {
		for (int label = 0; label <= 1; label++) {
			for (int node = 0; node < dims[0]/2; a++, node++) {
				// Cache visibility
				int num_occl = dims[2];
				int v1 = -1;
				int v2 = -1;
				for (int b = 0; b < dims[2]; b++) {
					vis[b] = VA[a+b*dims[0]];
					num_occl -= (int)vis[b];
					if (!vis[b]) {
						if (v1 == -1)
							v1 = b;
						else
							v2 = b;
					}
				}

				switch (num_occl) {
					case 0:
						// Calculate the mean
						calc_mean(M, IA, dims, a);
						// For each image
						for (int b = 0; b < dims[2]; b++)
							U[node*2+label] += calc_cost((Method *)0, M, IA, dims, a, b, Kocc); // Add to unary
						break;
					case 1:
						// Set up the connections
						PI[npair*2] = node + 1;
						PI[npair*2+1] = node + 1 + (1 + label) * dims[0] / 2 + v1 * dims[0];
						// Calculate the mean
						calc_mean(M, IA, dims, a, vis);
						// For each image
						for (int b = 0; b < dims[2]; b++) {
							if (vis[b])
								P[npair*4+label*2] += calc_cost((Method *)0, M, IA, dims, a, b, Kocc);
							else
								P[npair*4+label*2] += Kocc + 1;
						}
						// Calculate the other mean
						calc_mean(M, IA, dims, a);
						// For each image
						for (int b = 0; b < dims[2]; b++)
							P[npair*4+label*2+1] += calc_cost((Method *)0, M, IA, dims, a, b, Kocc);
						npair++;
						break;
					case 2:
						// Set up the connections
						TI[ntriple*3] = node + 1;
						TI[ntriple*3+1] = node + 1 + (1 + label) * dims[0] / 2 + v1 * dims[0];
						TI[ntriple*3+2] = node + 1 + (1 + label) * dims[0] / 2 + v2 * dims[0];
						// Calculate the mean
						calc_mean(M, IA, dims, a, vis);
						// For each image
						for (int b = 0; b < dims[2]; b++) {
							if (vis[b])
								T[ntriple*8+label*4] += calc_cost((Method *)0, M, IA, dims, a, b, Kocc);
							else
								T[ntriple*8+label*4] += Kocc + 1;
						}
						// Calculate another mean
						vis[v2] = 1;
						calc_mean(M, IA, dims, a, vis);
						// For each image
						for (int b = 0; b < dims[2]; b++) {
							if (vis[b])
								T[ntriple*8+label*4+1] += calc_cost((Method *)0, M, IA, dims, a, b, Kocc);
							else
								T[ntriple*8+label*4+1] += Kocc + 1;
						}
						// Calculate another mean
						vis[v2] = 0;
						vis[v1] = 1;
						calc_mean(M, IA, dims, a, vis);
						// For each image
						for (int b = 0; b < dims[2]; b++) {
							if (vis[b])
								T[ntriple*8+label*4+2] += calc_cost((Method *)0, M, IA, dims, a, b, Kocc);
							else
								T[ntriple*8+label*4+2] += Kocc + 1;
						}
						// Calculate last mean
						calc_mean(M, IA, dims, a);
						// For each image
						for (int b = 0; b < dims[2]; b++)
							T[ntriple*8+label*4+3] += calc_cost((Method *)0, M, IA, dims, a, b, Kocc);
						ntriple++;
						break;
					default:
						// Cache approximate visibility
						for (int b = 0; b < dims[2]; b++)
							vis[b] = V[a+b*dims[0]];

						// Calculate the mean
						calc_mean(M, IA, dims, a, vis);

						// For each image
						for (int b = 0; b < dims[2]; b++) {
							R data = calc_cost((Method *)0, M, IA, dims, a, b, Kocc);
							if (VA[a+b*dims[0]]) {
								// Add to unary
								U[node*2+label] += data;
							} else {
								// Add an edge
								PI[npair*2] = node + 1;
								PI[npair*2+1] = b * dims[0] + (1 + label) * dims[0] / 2 + node + 1;
								P[npair*4+label*2] = Kocc + 1;
								P[npair*4+label*2+1] = data;
								npair++;
							}
						}
						break;
				}
			}
		}
	} else {
		// Special case - dims[2] == 2
		for (int label = 0; label <= 1; label++) {
			for (int node = 0; node < dims[0]/2; a++, node++) {
				// Cache visibility
				int v1 = !VA[a];
				int v2 = !VA[a+dims[0]];

				switch (v1+v2*2) {
					case 0:
						U[node*2+label] += calc_cost((Method *)0, IA, dims, a, Kocc); // Add to unary
						break;
					case 1:
						// Set up the connections
						PI[npair*2] = node + 1;
						PI[npair*2+1] = node + 1 + (1 + label) * dims[0] / 2;
						// Occlusion cost
						P[npair*4+label*2] = Kocc + 1;
						// Data cost
						P[npair*4+label*2+1] = calc_cost((Method *)0, IA, dims, a, Kocc);
						npair++;
						break;
					case 2:
						// Set up the connections
						PI[npair*2] = node + 1;
						PI[npair*2+1] = node + 1 + (3 + label) * dims[0] / 2;
						// Occlusion cost
						P[npair*4+label*2] = Kocc + 1;
						// Data cost
						P[npair*4+label*2+1] = calc_cost((Method *)0, IA, dims, a, Kocc);
						npair++;
						break;
					case 3:
						// Set up the connections
						TI[ntriple*3] = node + 1;
						TI[ntriple*3+1] = node + 1 + (1 + label) * dims[0] / 2;
						TI[ntriple*3+2] = node + 1 + (3 + label) * dims[0] / 2;
						// Occlusion costs
						T[ntriple*8+label*4] = 2 * (Kocc + 1);
						T[ntriple*8+label*4+1] = Kocc + 1;
						T[ntriple*8+label*4+2] = Kocc + 1;
						// Data cost
						T[ntriple*8+label*4+3] = calc_cost((Method *)0, IA, dims, a, Kocc);
						ntriple++;
						break;
				}
			}
		}
	}
	mxFree(vis);
	mxFree(M);
	return;
}
