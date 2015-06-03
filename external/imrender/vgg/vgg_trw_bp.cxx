//[L energy lower_bound] = vgg_trw_bp(UE, PI, PE, options);

// $Id: vgg_trw_bp.cxx,v 1.3 2009/08/31 22:05:09 ojw Exp $

#include <mex.h>
#include "MRFEnergy.h"

// Define types
#ifdef _MSC_VER
typedef unsigned __int16 uint16_t;
typedef __int32 int32_t;
typedef unsigned __int32 uint32_t;
#else
#include <stdint.h>
#endif

static void erfunc(char *err) {mexErrMsgTxt(err);}
template<class TYPE> static inline void wrapper_func(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[], int options[], int *nlabels);

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
	// Check number of inputs
	if (nrhs < 3 || nrhs > 4)
		mexErrMsgTxt("Unexpected number of input arguments.");
	if (nlhs < 1 || nlhs > 3)
		mexErrMsgTxt("Unexpected number of outputs.");
	for (int i = 0; i < nrhs; i++) {
		if (mxIsComplex(prhs[i]))
			mexErrMsgTxt("Inputs must be real.");
	}

	// Check unary terms
	if (!mxIsCell(prhs[0]))
		mexErrMsgTxt("UE must be a cell array.");
	int n_nodes = mxGetNumberOfElements(prhs[0]);
	int *nlabels = new int[n_nodes];
	int max_labels = 0;
	int min_labels = 65537;
	for (int i = 0; i < n_nodes; i++) {
		mxArray *data_array = mxGetCell(prhs[0], i);
		if (!mxIsDouble(data_array))
			mexErrMsgTxt("UE cells must be doubles.");
		if (mxIsComplex(data_array))
			mexErrMsgTxt("UE cells must be real.");
		nlabels[i] = mxGetN(data_array);
		max_labels = nlabels[i] > max_labels ? nlabels[i] : max_labels;
		min_labels = nlabels[i] < min_labels ? nlabels[i] : min_labels;
	}
	if (max_labels > 65536)
		mexErrMsgTxt("A maximum of 65536 nodes per label are supported.");

	// Check options
	int options[] = {-1, max_labels, 1, 0, 30}; // ParamsPerEdge, max_labels, UseTRW, Type, max_iters
	if (nrhs == 4) {
		if (!mxIsInt32(prhs[3]))
			mexErrMsgTxt("options should be int32s");
		const int32_t *params = (const int32_t *)mxGetData(prhs[3]);
		switch (mxGetNumberOfElements(prhs[3])) {
			default:
			case 3:
				options[4] = (int)params[2];
			case 2:
				options[3] = (int)params[1];
			case 1:
				options[2] = (int)params[0] != 0;
			case 0:
				break;
		}
	}

	if (min_labels == 2 && max_labels == 2) {
		// Binary mode
		if (options[3] != 0)
			mexErrMsgTxt("For binary graphs, only the general form of energies is supported.");
		options[0] = 4; // 4 costs per edge
		if (nlhs < 2) {
			// Energy & lower bound not required, so use fast binary method
			wrapper_func<TypeBinaryFast>(nlhs, plhs, nrhs, prhs, options, nlabels);
		} else {
			// Use slow method
			mexPrintf("With binary graphs, optimization is faster if you have only 1 output argument.\n");
			wrapper_func<TypeBinary>(nlhs, plhs, nrhs, prhs, options, nlabels);
		}
	} else {
		// Non-binary mode
		// Which type of pairwise energies are used
		switch (options[3]) {
			case 0:
				// General form
				wrapper_func<TypeGeneral>(nlhs, plhs, nrhs, prhs, options, nlabels);
				break;
			case 1:
				// Truncated quadratic
				options[0] = 2; // 2 parameters per edge: alpha & trunc thresh
				wrapper_func<TypeTruncatedQuadratic>(nlhs, plhs, nrhs, prhs, options, nlabels);
				break;
			case 2:
				// Truncated linear
				options[0] = 2; // 2 parameters per edge: alpha & trunc thresh
				wrapper_func<TypeTruncatedLinear>(nlhs, plhs, nrhs, prhs, options, nlabels);
				break;
			default:
				// Type not implemented
				mexErrMsgTxt("Energy type not yet implemented. Why not give it a go yourself.");
				break;
		}
	}
	delete nlabels;
	return;
}

// Functions for general type
static inline MRFEnergy<TypeGeneral>::NodeId add_node(MRFEnergy<TypeGeneral> *mrf, int local_modes, TypeGeneral::REAL *graph_data)
{
	return mrf->AddNode(TypeGeneral::LocalSize(local_modes), TypeGeneral::NodeData(graph_data));
}
static inline void add_edge(MRFEnergy<TypeGeneral> *mrf, MRFEnergy<TypeGeneral>::NodeId node1, MRFEnergy<TypeGeneral>::NodeId node2, TypeGeneral::REAL *graph_data)
{
	mrf->AddEdge(node1, node2, TypeGeneral::EdgeData(TypeGeneral::GENERAL, graph_data));
}

// Functions for binary type
static inline MRFEnergy<TypeBinary>::NodeId add_node(MRFEnergy<TypeBinary> *mrf, int local_modes, TypeBinary::REAL *graph_data)
{
	return mrf->AddNode(TypeBinary::LocalSize(), TypeBinary::NodeData(graph_data[0], graph_data[1]));
}
static inline void add_edge(MRFEnergy<TypeBinary> *mrf, MRFEnergy<TypeBinary>::NodeId node1, MRFEnergy<TypeBinary>::NodeId node2, TypeBinary::REAL *graph_data)
{
	mrf->AddEdge(node1, node2, TypeBinary::EdgeData(graph_data[0], graph_data[2], graph_data[1], graph_data[3]));
}

// Functions for binary fast type
static inline MRFEnergy<TypeBinaryFast>::NodeId add_node(MRFEnergy<TypeBinaryFast> *mrf, int local_modes, TypeBinaryFast::REAL *graph_data)
{
	return mrf->AddNode(TypeBinaryFast::LocalSize(), TypeBinaryFast::NodeData(graph_data[0], graph_data[1]));
}
static inline void add_edge(MRFEnergy<TypeBinaryFast> *mrf, MRFEnergy<TypeBinaryFast>::NodeId node1, MRFEnergy<TypeBinaryFast>::NodeId node2, TypeBinaryFast::REAL *graph_data)
{
	mrf->AddEdge(node1, node2, TypeBinaryFast::EdgeData(graph_data[0], graph_data[2], graph_data[1], graph_data[3]));
}

// Functions for truncated quadratic type
static inline MRFEnergy<TypeTruncatedQuadratic>::NodeId add_node(MRFEnergy<TypeTruncatedQuadratic> *mrf, int local_modes, TypeTruncatedQuadratic::REAL *graph_data)
{
	return mrf->AddNode(TypeTruncatedQuadratic::LocalSize(), TypeTruncatedQuadratic::NodeData(graph_data));
}
static inline void add_edge(MRFEnergy<TypeTruncatedQuadratic> *mrf, MRFEnergy<TypeTruncatedQuadratic>::NodeId node1, MRFEnergy<TypeTruncatedQuadratic>::NodeId node2, TypeTruncatedQuadratic::REAL *graph_data)
{
	mrf->AddEdge(node1, node2, TypeTruncatedQuadratic::EdgeData(graph_data[0], graph_data[1]));
}

// Functions for truncated linear type
static inline MRFEnergy<TypeTruncatedLinear>::NodeId add_node(MRFEnergy<TypeTruncatedLinear> *mrf, int local_modes, TypeTruncatedLinear::REAL *graph_data)
{
	return mrf->AddNode(TypeTruncatedLinear::LocalSize(), TypeTruncatedLinear::NodeData(graph_data));
}
static inline void add_edge(MRFEnergy<TypeTruncatedLinear> *mrf, MRFEnergy<TypeTruncatedLinear>::NodeId node1, MRFEnergy<TypeTruncatedLinear>::NodeId node2, TypeTruncatedLinear::REAL *graph_data)
{
	mrf->AddEdge(node1, node2, TypeTruncatedLinear::EdgeData(graph_data[0], graph_data[1]));
}

template<class TYPE> static inline void wrapper_func(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[], int options[], int *nlabels)
{
	int n_nodes = mxGetNumberOfElements(prhs[0]);

	// edges
	int Pindices = mxGetM(prhs[1]);
	if (Pindices != 2 && Pindices != 3)
		mexErrMsgTxt("Unexpected dimensions for PI");
	int nP = mxGetN(prhs[1]);
	if (!mxIsCell(prhs[2]))
		mexErrMsgTxt("PE must be a cell array.");
	int nPE = mxGetNumberOfElements(prhs[2]);
	if (Pindices == 2 && nP != nPE)
		mexErrMsgTxt("Without explicit indices, PI should have the same number of PE as pairwise hase cells");
	if (!mxIsUint32(prhs[1]))
		mexErrMsgTxt("PI should be uint32s");
	const uint32_t *PI = (const uint32_t *)mxGetData(prhs[1]);

	MRFEnergy<TYPE> *mrf = new MRFEnergy<TYPE>(typename TYPE::GlobalSize(options[1]), erfunc);
	typename MRFEnergy<TYPE>::NodeId *nodes =  new typename MRFEnergy<TYPE>::NodeId[n_nodes];
	typedef typename TYPE::REAL REAL; 
	REAL *graph_data = (REAL *)mxCalloc(options[1]*options[1], sizeof(REAL));
	
	// Add unary energies
	for (int i = 0; i < n_nodes; i++) {
		mxArray *data_array = mxGetCell(prhs[0], i);
		int pitch = mxGetM(data_array);
		const double *data = mxGetPr(data_array);
		for (int j = 0; j < nlabels[i]; j++)
			graph_data[j] = (REAL)data[j*pitch];
		nodes[i] = add_node(mrf, nlabels[i], graph_data);
	}

	// Add pairwise energies
	for (int i = 0; i < nP; i++, PI += Pindices) {
		// Calculate index of energy table
		mxArray *data_array;
		if (Pindices == 2)
			data_array = mxGetCell(prhs[2], i);
		else {
			if (PI[2] < 1 || PI[2] > nPE)
				mexErrMsgTxt("Column of PI references invalid cell of PE");
			data_array = mxGetCell(prhs[2], PI[2]-1);
		}

		if (!mxIsDouble(data_array))
			mexErrMsgTxt("pairwise cells must be doubles.");
		if (mxIsComplex(data_array))
			mexErrMsgTxt("pairwise cells must be real.");
		const double *data = mxGetPr(data_array);
		int len;
		int n1 = PI[0] - 1;
		int n2 = PI[1] - 1;
		if (options[0] < 0) {
			if (mxGetM(data_array) != nlabels[n1] || mxGetN(data_array) != nlabels[n2])
				mexErrMsgTxt("PE cell has unexpected dimensions.");
			len = nlabels[n1] * nlabels[n2];
		} else {
			if (mxGetNumberOfElements(data_array) != options[0])
				mexErrMsgTxt("PE cell has unexpected number of elements.");
			len = options[0];
		}
		for (int j = 0; j < len; j++)
			graph_data[j] = (REAL)data[j];
		add_edge(mrf, nodes[n1], nodes[n2], graph_data);
	}

	mxFree(graph_data);

	// Function below is optional - it may help if, for example, nodes are added in a random order
	mrf->SetAutomaticOrdering();
	
	typename MRFEnergy<TYPE>::Options mrf_options;
	mrf_options.m_iterMax = options[4]; // maximum number of iterations
	REAL energy, lowerBound = 0;

	if (options[2]) {
		mexPrintf("Graph loaded. Starting optimization using TRW-S.\n");
		/////////////////////// TRW-S algorithm //////////////////////
		mrf->Minimize_TRW_S(mrf_options, lowerBound, energy);
	} else {
		mexPrintf("Graph loaded. Starting optimization using BP.\n");
		//////////////////////// BP algorithm ////////////////////////
		mrf->ZeroMessages(); // in general not necessary - it may be faster to start 
		// with messages computed in previous iterations
		mrf->Minimize_BP(mrf_options, energy);
	}

	// Read solution
	plhs[0] = mxCreateNumericMatrix(mxGetM(prhs[0]), mxGetN(prhs[0]), mxUINT16_CLASS, mxREAL);
	uint16_t *L = (uint16_t *)mxGetData(plhs[0]);
	for (int i = 0; i < n_nodes; i++ )
		L[i] = (uint16_t)mrf->GetSolution(nodes[i]) + 1;

	if (nlhs > 1) {
		plhs[1] = mxCreateDoubleScalar((double)energy);
		if (nlhs > 2) {
			plhs[2] = mxCreateDoubleScalar((double)lowerBound);
		}
	}

	// Clean up
	delete nodes;
	delete mrf;
}
