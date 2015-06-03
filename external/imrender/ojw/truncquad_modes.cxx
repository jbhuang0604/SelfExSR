// [modes depth energy inliers] = truncquad_modes(I, thresh, use_variance, search_width)

// $Id: truncquad_modes.cxx,v 1.2 2008/11/17 11:27:35 ojw Exp $

#include <mex.h>
#include <memory.h>

static inline double update_energy(const double *A, double *C, const int m, const int n, const double thresh, double *temp, mxLogical *bt, int *ic);
static inline double calculate_energy(const double *A, const double *C, const int m, const int n, const double thresh);
static inline double calculate_energy(const double *A, const double *C, const int m, const int n, const double thresh, int use_variance);

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
	// Check number of arguments
	if (nrhs < 2 || nrhs > 4)
		mexErrMsgTxt("Unexpected number of input arguments.");
	if (nlhs < 1 || nlhs > 4)
		mexErrMsgTxt("Unexpected number of output arguments.");

	// Check argument types are valid
	mxClassID in_class = mxGetClassID(prhs[0]);
	if (mxGetClassID(prhs[0]) != mxDOUBLE_CLASS || mxGetClassID(prhs[1]) != mxDOUBLE_CLASS)
		mexErrMsgTxt("I and thresh must be doubles");
	for (int i = 0; i < nrhs; i++) {
        if (mxIsComplex(prhs[i]))
			mexErrMsgTxt("Inputs cannot be complex.");
	}

	// Get and check array dimensions
	int ndims = mxGetNumberOfDimensions(prhs[0]);
	if (ndims != 3)
		mexErrMsgTxt("I should have 3 dimensions");
	const int *dims = mxGetDimensions(prhs[0]);
	int d_step = dims[0] * dims[1];

	// Get pointers to input arrays
	const double *I = (const double *)mxGetData(prhs[0]);
	double thresh = mxGetScalar(prhs[1]);
	double pair_thresh = thresh * 4;
	double normalizer = 1 / (double)dims[1];

	int use_variance = -1;
	int search_width = 2 * dims[2];
	if (nrhs > 2) {
		if (mxGetNumberOfElements(prhs[2]) != 1 || !mxIsNumeric(prhs[2]))
			mexErrMsgTxt("use_variance should be a scalar.");
		use_variance = (int)mxGetScalar(prhs[2]) - 1;
		if (nrhs > 3) {
			if (mxGetNumberOfElements(prhs[3]) != 1 || !mxIsNumeric(prhs[3]))
				mexErrMsgTxt("search_width should be a scalar.");
			search_width = (int)mxGetScalar(prhs[3]);
		}
	}

	// Create the temporary data store
	int k = (dims[1]*(dims[1]-1)*dims[2])/2;
	double *modes = (double *)mxMalloc(dims[0]*(k+1)*sizeof(double));
	double *depth = (double *)mxMalloc(k*sizeof(double));
	double *energy = (double *)mxMalloc((k+((dims[1]*(dims[1]-1))/2))*sizeof(double));
	double *e_thisdepth = &energy[k-1];
	double *centre = modes;
	double *temp_buf = &centre[dims[0]*k];
	mxLogical *inliers = (mxLogical *)mxMalloc(dims[1]*k*sizeof(mxLogical));
	mxLogical *below_thresh = inliers;
	bool *seen_before = (bool *)mxMalloc(dims[1]*(dims[1]-1)*sizeof(bool));
	k = 0;
	int d_offset = 0;

	// For each depth...
	for (int d = 0; d < dims[2]; d++) {
		int j = 0;
		// Reset if seen before
		memset(seen_before, 0, sizeof(bool)*dims[1]*(dims[1]-1));
		// For each unique pair of points in the set...
		for (int p1 = 0; p1 < dims[1]-1; p1++) {
			for (int p2 = p1+1; p2 < dims[1]; p2++) {
				// Check we haven't seen this pair in a cluster already
				if (seen_before[p1*dims[1]+p2])
					continue;

				// Are they in the same cluster?
				double dist = 0;
				for (int i = 0; i < dims[0]; i++) {
					double diff = I[d_offset+p1*dims[0]+i] - I[d_offset+p2*dims[0]+i];
					dist += diff * diff;
				}
				if (dist > pair_thresh)
					continue;

				// Calculate the starting centre of the cluster
				for (int i = 0; i < dims[0]; i++)
					centre[i] = (I[d_offset+p1*dims[0]+i] + I[d_offset+p2*dims[0]+i]) * 0.5;

				// Move to the local minimum using mean shift
				int n_inliers;
				double e_up, e_curr = -1;
				do {
					e_up = e_curr;
					e_curr = update_energy(&I[d_offset], centre, dims[0], dims[1], thresh, temp_buf, below_thresh, &n_inliers);
				} while (e_up != e_curr);

				// Ensure there are two inliers
				if (n_inliers < 2)
					continue;

				int lim; // Do here not later to keep gcc happy about goto

				// Check we haven't seen this mode before 
				// (assume against the freakish possibility of identical energies for different centres)
				for (int i = 0; i < j; i++) {
					if (e_curr == e_thisdepth[i])
						goto skip_point;
				}
				e_thisdepth[j] = e_curr;
				j++;

				// Mark the other pairs in this cluster so we save time later
				for (int p3 = p1; p3 < dims[1]-1; p3++) {
					for (int p4 = p2; p4 < dims[1]; p4++) {
						seen_before[p3*dims[1]+p4] |= below_thresh[p3] & below_thresh[p4];
					}
				}

				if (use_variance < 0) {
					// Check the other depths to see if this is a colour mode
					lim = d + search_width;
					lim = lim < dims[2] ? lim : dims[2];
					for (int d2 = d+1; d2 < lim; d2++) {
						e_up = calculate_energy(&I[d2*d_step], centre, dims[0], dims[1], thresh);
						if (e_up < e_curr)
							goto skip_point;
					}
					lim = d - search_width;
					lim = lim > 0 ? lim : 0;
					for (int d2 = d-1; d2 >= lim; d2--) {
						e_up = calculate_energy(&I[d2*d_step], centre, dims[0], dims[1], thresh);
						if (e_up < e_curr)
							goto skip_point;
					}
				} else {
					// Check the other depths to see if this is a colour mode
					double e_curr2 = (e_curr - thresh * (dims[1] - n_inliers)) / (n_inliers - use_variance);
					lim = d + search_width;
					lim = lim < dims[2] ? lim : dims[2];
					for (int d2 = d+1; d2 < lim; d2++) {
						e_up = calculate_energy(&I[d2*d_step], centre, dims[0], dims[1], thresh, use_variance);
						if (e_up < e_curr2)
							goto skip_point;
					}
					lim = d - search_width;
					lim = lim > 0 ? lim : 0;
					for (int d2 = d-1; d2 >= lim; d2--) {
						e_up = calculate_energy(&I[d2*d_step], centre, dims[0], dims[1], thresh, use_variance);
						if (e_up < e_curr2)
							goto skip_point;
					}
				}

				// Add this mode to the list
				centre += dims[0];
				below_thresh += dims[1];
				depth[k] = d + 1;
				energy[k] = e_curr * normalizer; // Normalize by number of input images
				k++;
skip_point:
				continue;
			}
		}
		d_offset += d_step;
	}
	
	// Create output arrays
    mxFree(seen_before);
    if (nlhs > 3) {
        plhs[3] = mxCreateNumericMatrix(dims[1], 0, mxLOGICAL_CLASS, mxREAL);
        mxSetN(plhs[3], k);
        mxSetData(plhs[3], mxRealloc(inliers, dims[1]*k*sizeof(mxLogical)));
    } else {
        mxFree(inliers);
    }
    if (nlhs > 2) {
        plhs[2] = mxCreateNumericMatrix(1, 0, mxDOUBLE_CLASS, mxREAL);
        mxSetN(plhs[2], k);
        mxSetData(plhs[2], mxRealloc(energy, k*sizeof(double)));
    } else {
        mxFree(energy);
    }
    if (nlhs > 1) {
        plhs[1] = mxCreateNumericMatrix(1, 0, mxDOUBLE_CLASS, mxREAL);
        mxSetN(plhs[1], k);
        mxSetData(plhs[1], mxRealloc(depth, k*sizeof(double)));
    } else {
        mxFree(depth);
    }
    plhs[0] = mxCreateNumericMatrix(dims[0], 0, mxDOUBLE_CLASS, mxREAL);
    mxSetN(plhs[0], k);
    mxSetData(plhs[0], mxRealloc(modes, dims[0]*k*sizeof(double)));
	return;
}

static inline double update_energy(const double *A, double *C, const int m, const int n, const double thresh, double *temp, mxLogical *bt, int *ic)
{
	double e_out = 0;
	int in_cluster = 0;
	// Clear the temp buffer
	for (int b = 0; b < m; b++)
		temp[b] = 0;

	// For each vector in the set...
	for (int a = 0; a < n; a++) {
		// Calculate the distance from the centre
		double dist = 0;
		for (int b = 0; b < m; b++) {
			double diff = A[b] - C[b];
			dist += diff * diff;
		}
		// Truncate if necessary
		if (dist <= thresh) {
			e_out += dist;
			bt[a] = 1;
			in_cluster++;
			// Add on to the centre coordinates
			for (int b = 0; b < m; b++)
				temp[b] += A[b];
		} else {
			e_out += thresh;
			bt[a] = 0;
		}
		A += m;
	}
	// Calculate the centre coordinates
	double recip = 1 / (double)in_cluster;
	for (int b = 0; b < m; b++)
		C[b] = temp[b] * recip;

	*ic = in_cluster;
	return e_out;
}

static inline double calculate_energy(const double *A, const double *C, const int m, const int n, const double thresh)
{
	double e_out = 0;

	// For each vector in the set...
	for (int a = 0; a < n; a++) {
		// Calculate the distance from the centre
		double dist = 0;
		for (int b = 0; b < m; b++) {
			double diff = A[b] - C[b];
			dist += diff * diff;
		}
		A += m;
		// Truncate if necessary
		e_out += (dist < thresh) ? dist : thresh;
	}

	return e_out;
}

static inline double calculate_energy(const double *A, const double *C, const int m, const int n, const double thresh, int use_variance)
{
	double e_out = 0;
	int in_cluster = 0;

	// For each vector in the set...
	for (int a = 0; a < n; a++) {
		// Calculate the distance from the centre
		double dist = 0;
		for (int b = 0; b < m; b++) {
			double diff = A[b] - C[b];
			dist += diff * diff;
		}
		A += m;
		// Truncate if necessary
		if (dist <= thresh) {
			e_out += dist;
			in_cluster++;
		}
	}
	if (in_cluster < 2)
		return 1e300;

	return e_out / (in_cluster - use_variance);
}
