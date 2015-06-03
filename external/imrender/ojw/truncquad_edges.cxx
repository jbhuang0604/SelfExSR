// edge_costs = truncquad_edges(I, modes, EI, thresh)
// I - CxLxN array of colour libraries
// modes - 1xN cell array of [CxM] containing modes
// EI - 2xP matrix of offsets for modes

// $Id: truncquad_edges.cxx,v 1.2 2008/11/17 11:27:35 ojw Exp $

#include <mex.h>

// Define types
#ifdef _MSC_VER
typedef unsigned __int8 uint8_t;
typedef unsigned __int32 uint32_t;
#else
#include <stdint.h>
#endif

template<class T> static inline void wrapper_func(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]);

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
	// Check number of arguments
	if (nrhs != 5)
		mexErrMsgTxt("Unexpected number of input arguments.");

	// Call the function according to the type
	switch (mxGetClassID(prhs[0])) {
		case mxDOUBLE_CLASS:
			wrapper_func<double>(nlhs, plhs, nrhs, prhs);
			break;
		case mxSINGLE_CLASS:
			wrapper_func<float>(nlhs, plhs, nrhs, prhs);
			break;
		case mxUINT8_CLASS:
			wrapper_func<uint8_t>(nlhs, plhs, nrhs, prhs);
			break;
		default:
			mexErrMsgTxt("I is of an unsupported type");
			break;
	}
	return;
}

static inline void cache_colour(const double *I, double **colour_vector, int num_chans)
{
    *colour_vector = (double *)I;
}

template<class T> static inline void cache_colour(const T *I, double **colour_vector, int num_chans)
{
    for (int i = 0; i < num_chans; i++)
        (*colour_vector)[i] = (double)(I[i]);
}

template<class T> static inline void wrapper_func(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
	// Check number of arguments
	if (nlhs != 1)
		mexErrMsgTxt("Unexpected number of output arguments.");

	// Get and check array dimensions
	if (mxGetNumberOfDimensions(prhs[0]) != 3)
		mexErrMsgTxt("I has an unexpected number of dimensions.");
	const int *Idims = (const int *)mxGetDimensions(prhs[0]);
	if (!mxIsCell(prhs[1]))
		mexErrMsgTxt("modes must be a cell array.");
	if (Idims[2] != mxGetNumberOfElements(prhs[1]))
		mexErrMsgTxt("3rd dimension of I must be the same as the numel(modes).");
	if (!mxIsUint32(prhs[2]))
		mexErrMsgTxt("EI must be a uint32 matrix.");
	if (2 != mxGetM(prhs[2]))
		mexErrMsgTxt("EI must have 2 rows.");
	int num_edges = mxGetN(prhs[2]);

	// Check all modes
	int *nmodes = (int *)mxMalloc(Idims[2]*sizeof(int));
	int max_nmodes = 0;
	for (int i = 0; i < Idims[2]; i++) {
		mxArray *mode = mxGetCell(prhs[1], i);
        if (!mxIsDouble(mode))
			mexErrMsgTxt("modes must be doubles.");
        if (mxIsComplex(mode))
			mexErrMsgTxt("modes cannot be complex.");
        if (mxGetM(mode) != Idims[0])
			mexErrMsgTxt("A mode has the wrong number of channels.");
		nmodes[i] = mxGetN(mode);
		max_nmodes = nmodes[i] > max_nmodes ? nmodes[i] : max_nmodes;
	}

	// Get pointers to input arrays
	const T *I = (T *)mxGetData(prhs[0]);
	const uint32_t *EI = (uint32_t *)mxGetData(prhs[2]);
	double thresh = mxGetScalar(prhs[3]);
	double weight = mxGetScalar(prhs[4]);
	
	// Create output array
	plhs[0] = mxCreateCellMatrix(num_edges, 1);

	// Cache for all distances of the first node's modes
	// May be faster if user has sorted the indices
	double *node1_cache = (double *)mxMalloc(((max_nmodes+1)*Idims[1]+Idims[0])*sizeof(double));
    double *min_val = &node1_cache[max_nmodes*Idims[1]];
    double *colour_vector = &node1_cache[max_nmodes*Idims[1]+Idims[1]];
	int cached_node1 = -3487; // Random to avoid a bug with first cacheing
	int pitch[2];
	const double *mode[2];
	int dims[2];

	// Go through each clique...
	for (int c = 0; c < num_edges; c++, EI += 2) {
        // Load first mode
        int ind = EI[0] - 1;
        if (ind != cached_node1) {
            if (ind >= Idims[2])
                mexErrMsgTxt("Invalid index in EI.");
            pitch[0] = ind * Idims[0] * Idims[1];
            dims[0] = nmodes[ind];
            mode[0] = (const double *)mxGetData(mxGetCell(prhs[1], ind));
            
            // Cache the distance to each mode for each colour in the library
            // For each vector in the library
            for (int v = 0; v < Idims[1]; v++) {
                // Cache the colour
                cache_colour(&I[v*Idims[0]+pitch[0]], &colour_vector, Idims[0]);
                
                // Calculate the distance to each mode
                min_val[v] = thresh + 1;
                for (int a = 0; a < dims[0]; a++) {
                    double dist = 0;
                    // Go through each channel
                    for (int i = 0; i < Idims[0]; i++) {
                        double diff = colour_vector[i] - mode[0][i+a*Idims[0]];
                        dist += diff * diff;
                    }
                    min_val[v] = min_val[v] > dist ? dist : min_val[v];
                    node1_cache[a+v*dims[0]] = dist;
				}
                min_val[v] = thresh - min_val[v];
            }
            cached_node1 = ind;
        }
        
        // Load second mode
        ind = EI[1] - 1;
        if (ind >= Idims[2])
            mexErrMsgTxt("Invalid index in EI.");
        pitch[1] = ind * Idims[0] * Idims[1];
        dims[1] = nmodes[ind];
        mode[1] = (const double *)mxGetData(mxGetCell(prhs[1], ind));
        
        // Create energy matrix
        mxArray *en = mxCreateNumericArray(2, dims, mxDOUBLE_CLASS, mxREAL); // This way doesn't set data
		double *edge_cost = mxGetPr(en);

		// Initialize edge costs
		for (int d = 0; d < dims[0]*dims[1]; d++)
			edge_cost[d] = thresh;

		// For each vector in the library
		for (int v = 0; v < Idims[1]; v++) {
            if (min_val[v] < 0)
                continue; // Skip as all node1 modes are above thresh
                
            // Cache the colour
            cache_colour(&I[v*Idims[0]+pitch[1]], &colour_vector, Idims[0]);
                    
            // Calculate the distance to each mode
            for (int b = 0; b < dims[1]; b++) {
                double dist = 0;
                // Go through each channel
                for (int i = 0; i < Idims[0]; i++) {
                    double diff = colour_vector[i] - mode[1][i+b*Idims[0]];
                    dist += diff * diff;
                }
                
                // For each mode combination...
                if (dist < min_val[v]) {
                    for (int a = 0; a < dims[0]; a++) {
						double dist2 = node1_cache[a+v*dims[0]] + dist;
						edge_cost[b*dims[0]+a] = dist2 < edge_cost[b*dims[0]+a] ? dist2 : edge_cost[b*dims[0]+a];
					}
                }   
			}
		}

		// Weight edge costs
		for (int d = 0; d < dims[0]*dims[1]; d++)
			edge_cost[d] *= weight;

		mxSetCell(plhs[0], c, en);
	}
	mxFree(node1_cache);
	mxFree(nmodes);
	return;
}