// Calls Valdimir Kolmogorov's QPBO code, downloaded from:
// http://www.adastral.ucl.ac.uk/~vladkolm/software.html

// L = vgg_qpbo(UE, PI, PE[, TI, TE], [params]);
// T UE 2xM matrix of unary terms for M nodes, or scalar = M if all unary terms are zero
// uint32 PI 2or3xN matrix, each column [start_node end_node [pairwise_energy_table_index]]
// T PE 4xP pairwise energy table, each column [E00 E01 E10 E11]
// uint32 TI 3or4xQ matrix, each column [node1 node2 node3 [triple_energy_table_index]]
// T TE 8xR triple clique energy table, each column [E000 E001 E010 E011 E100 E101 E110 E111]

#include "QPBO.h"

static void erfunc(char *err) {mexErrMsgTxt(err);}
template<class INPUT, class INTERNAL> static inline void wrapper_func(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[], INTERNAL infinite_edge_cost);

// Define types
#ifdef _MSC_VER
#pragma warning(disable: 4661)
typedef __int32 int32_t;
typedef __int64 int64_t;
typedef unsigned __int32 uint32_t;
#else
#include <stdint.h>
#endif

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
	// Check number of inputs
	if (nrhs < 3 || nrhs > 6)
		mexErrMsgTxt("Unexpected number of input arguments.");
	
	// Check input types are valid
	mxClassID in_class = mxGetClassID(prhs[2]);
	if (mxGetNumberOfElements(prhs[0]) > 1 && in_class != mxGetClassID(prhs[0]))
		mexErrMsgTxt("Types of the input arguments don't agree.");
	if (nrhs > 4 && mxGetNumberOfElements(prhs[4]) > 1 && in_class != mxGetClassID(prhs[4]))
		mexErrMsgTxt("Types of the input arguments don't agree.");

	// Do we want a larger internal representation
	bool largerInternal = nrhs == 4 || nrhs == 6;
    if (largerInternal) {
        const mxArray *paramsArray = prhs[nrhs-1];
        if (mxIsCell(paramsArray)) {
            if (mxGetNumberOfElements(paramsArray) < 1 || mxGetNumberOfElements(paramsArray) > 2)
                mexErrMsgTxt("options cell array has unexpected size.");
            paramsArray = mxGetCell(paramsArray, 0);
        }
		if (!mxIsInt32(paramsArray))
			mexErrMsgTxt("options should be int32s");
        largerInternal = mxGetNumberOfElements(paramsArray) > 3 && ((const int32_t *)mxGetData(paramsArray))[3];
    }

	// Call the wrapper function according to the input type
	switch (in_class) {
		case mxINT32_CLASS:
			if (largerInternal)
				wrapper_func<int32_t, int64_t>(nlhs, plhs, nrhs, prhs, (int64_t)1<<47);
			else
				wrapper_func<int32_t, int32_t>(nlhs, plhs, nrhs, prhs, (int32_t)1<<20);
			break;
		case mxINT64_CLASS:
			wrapper_func<int64_t, int64_t>(nlhs, plhs, nrhs, prhs, (int64_t)1<<47);
			break;
		case mxDOUBLE_CLASS:
			mexPrintf("Warning: using non-integer energies removes optimality guarantees!\n");
			wrapper_func<double, double>(nlhs, plhs, nrhs, prhs, (double)1e100);
			break;
		default:
			mexErrMsgTxt("Inputs are of an unsupported type");
			break;
	}
	return;
}

template<class INPUT, class INTERNAL> static inline void wrapper_func(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[], INTERNAL infinite_edge_cost)
{
	// Check inputs and outputs
	if (nlhs < 1 || nlhs > 2)
		mexErrMsgTxt("Unexpected number of outputs.");
	for (int i = 0; i < nrhs; i++) {
		if (mxIsComplex(prhs[i]))
			mexErrMsgTxt("Inputs must be real.");
	}

	// Check types and dimensions, and get pointers
	int Pindices = mxGetM(prhs[1]);
	if (Pindices != 2 && Pindices != 3)
		mexErrMsgTxt("Unexpected dimensions for PI");
	int nP = mxGetN(prhs[1]);
	int nPE = mxGetN(prhs[2]);
	if (Pindices == 2 && nP != nPE)
		mexErrMsgTxt("Without explicit indices, PI and PE should have the same number of columns");
	if (mxGetM(prhs[2]) != 4)
		mexErrMsgTxt("Unexpected dimensions for PE");
	if (!mxIsUint32(prhs[1]))
		mexErrMsgTxt("PI should be uint32s");
	int nU;
	if (mxGetNumberOfElements(prhs[0]) == 1)
		nU = (int)mxGetScalar(prhs[0]);
	else {
		if (mxGetM(prhs[0]) != 2)
			mexErrMsgTxt("Unexpected dimensions for U");
		nU = mxGetN(prhs[0]);
	}
	int nT = 0;
	int nTE = 0;
	int Tindices = 0;
	const INPUT *U = (const INPUT *)mxGetData(prhs[0]);
	const INPUT *PE = (const INPUT *)mxGetData(prhs[2]);
	const INPUT *TE = (const INPUT *)NULL;
	const uint32_t *PI = (const uint32_t *)mxGetData(prhs[1]);
	const uint32_t *TI = (const uint32_t *)NULL;
	if (nrhs > 4) {
		Tindices = mxGetM(prhs[3]);
		if (Tindices != 3 && Tindices != 4)
			mexErrMsgTxt("Unexpected dimensions for TI");
		nT = mxGetN(prhs[3]);
		nTE = mxGetN(prhs[4]);
		if (Tindices == 3 && nT != nTE)
			mexErrMsgTxt("Without explicit indices, TI and TE should have the same number of columns");
		if (mxGetM(prhs[4]) != 8)
			mexErrMsgTxt("Unexpected dimensions for TE");
		if (!mxIsUint32(prhs[3]))
			mexErrMsgTxt("TI should be uint32s");
		TE = (const INPUT *)mxGetData(prhs[4]);
		TI = (const uint32_t *)mxGetData(prhs[3]);
	}
	// Default options
	int firstNNodes = nU; // Only return labels of first N nodes
	int improveMethod = 0; // Method to improve unknown nodes: 0 - none, 1 - improve (assume 0 best), 2 - optimal splice
	int contractIters = 0; // Number of contract cycles to try
	if (nrhs == 4 || nrhs == 6) {
        const mxArray *paramsArray = prhs[nrhs-1];
        if (mxIsCell(paramsArray))
            paramsArray = mxGetCell(paramsArray, 0);
		const int32_t *params = (const int32_t *)mxGetData(paramsArray);
		switch (mxGetNumberOfElements(paramsArray)) {
			default:
			case 3:
				contractIters = (int)params[2];
			case 2:
				improveMethod = (int)params[1];
			case 1:
				firstNNodes = (int)params[0];
				firstNNodes = nU < firstNNodes ? nU : firstNNodes;
			case 0:
				break;
		}
	}

	// Define the graph
	QPBO<INTERNAL> *graph = new QPBO<INTERNAL>(nU+nT, nP+6*nT, erfunc);

	// Create graph nodes
	int start_node = graph->AddNode(nU);
	code_assert(start_node == 0);

	// Enter unary terms
	if (mxGetNumberOfElements(prhs[0]) > 1) {
		for (int nodeInd = 0; nodeInd < nU; nodeInd++) {
			// Add unary potetial
			graph->AddUnaryTerm(nodeInd, (INTERNAL)U[0], (INTERNAL)U[1]);
			U += 2;
		}
	}

	// Enter pairwise terms
	int countIrregular = 0;
	const INPUT *INPUT_ptr = PE - 4;
	const uint32_t *UL_end = &PI[Pindices*nP];
	for (int i = 0; PI < UL_end; PI += Pindices) {
		
		// Calculate index of energy table
		if (Pindices == 2)
			INPUT_ptr += 4;
		else {
			if (PI[2] < 1 || PI[2] > nPE)
				mexErrMsgTxt("Column of PI references invalid column of PE");
			INPUT_ptr = &PE[4*(PI[2]-1)];
		}

		// Check submodularity 
		//if (PE_ptr[0] + PE_ptr[3] > PE_ptr[1] + PE_ptr[2])
		//	countIrregular++;

		graph->AddPairwiseTerm(PI[0]-1, PI[1]-1, (INTERNAL)INPUT_ptr[0], (INTERNAL)INPUT_ptr[1], (INTERNAL)INPUT_ptr[2], (INTERNAL)INPUT_ptr[3]);
	}

	// Enter triple clique terms
	int nNodes = nU;
	if (nT) {
		INPUT_ptr = TE - 8;
		UL_end = &TI[Tindices*nT];
		for (; TI < UL_end; TI += Tindices) {

			// Calculate index of energy table
			if (Tindices == 3)
				INPUT_ptr += 8;
			else {
				if (TI[3] < 1 || TI[3] > nTE)
					mexErrMsgTxt("Column of TI references invalid column of TE");
				INPUT_ptr = &TE[8*(TI[3]-1)];
			}

			int node1 = TI[0] - 1;
			int node2 = TI[1] - 1;
			int	node3 = TI[2] - 1;

			INTERNAL A = (INTERNAL)INPUT_ptr[0]; // E000
			INTERNAL B = (INTERNAL)INPUT_ptr[1]; // E001
			INTERNAL C = (INTERNAL)INPUT_ptr[2]; // E010
			INTERNAL D = (INTERNAL)INPUT_ptr[3]; // E011
			INTERNAL E = (INTERNAL)INPUT_ptr[4]; // E100
			INTERNAL F = (INTERNAL)INPUT_ptr[5]; // E101
			INTERNAL G = (INTERNAL)INPUT_ptr[6]; // E110
			INTERNAL H = (INTERNAL)INPUT_ptr[7]; // E111

			INTERNAL pi = (A + D + F + G) - (B + C + E + H);

			if (pi >= 0) {
				graph->AddPairwiseTerm(node1, node2, 0, C-A, 0, G-E);
				//if (C-A < G-E)
				//	countIrregular++;
				graph->AddPairwiseTerm(node1, node3, 0, 0, E-A, F-B);
				//if (E-A < F-B)
				//	countIrregular++;
				graph->AddPairwiseTerm(node2, node3, 0, B-A, 0, D-C);
				//if (B-A < D-C)
				//	countIrregular++;

				if (pi > 0) {
					// Add node
					int node4 = graph->AddNode();
					graph->AddUnaryTerm(node4, A, A-pi);
					graph->AddPairwiseTerm(node1, node4, 0, pi, 0, 0);
					graph->AddPairwiseTerm(node2, node4, 0, pi, 0, 0);
					graph->AddPairwiseTerm(node3, node4, 0, pi, 0, 0);
					nNodes++;
				}
			} else {
				graph->AddPairwiseTerm(node1, node2, B-D, 0, F-H, 0);
				//if (F-H < B-D)
				//	countIrregular++;
				graph->AddPairwiseTerm(node1, node3, C-G, D-H, 0, 0);
				//if (D-H < C-G)
				//	countIrregular++;
				graph->AddPairwiseTerm(node2, node3, E-F, 0, G-H, 0);
				//if (G-H < E-F)
				//	countIrregular++;

				// Add node
				int node4 = graph->AddNode();
				graph->AddUnaryTerm(node4, H+pi, H);
				graph->AddPairwiseTerm(node1, node4, 0, 0, -pi, 0);
				graph->AddPairwiseTerm(node2, node4, 0, 0, -pi, 0);
				graph->AddPairwiseTerm(node3, node4, 0, 0, -pi, 0);
				nNodes++;
			}
		}
		// Merge edges
		graph->MergeParallelEdges();
	}

	// Solve for optimimum
	graph->Solve();

	// Label a few more nodes if there are several global minima, and clump unlabelled nodes into groups
	graph->ComputeWeakPersistencies();

	// Create the output array
	plhs[0] = mxCreateNumericMatrix(firstNNodes, 1, mxINT32_CLASS, mxREAL);   
	int32_t *labelOutPtr = (int32_t *)mxGetData(plhs[0]);

	// Read out labelling
	int countUnlabel = 0;
	int *listUnlabel = new int[firstNNodes];
	for (int nodeCount = 0; nodeCount < firstNNodes; nodeCount++) {
		labelOutPtr[nodeCount] = (int32_t)graph->GetLabel(nodeCount);
		if (labelOutPtr[nodeCount] < 0)
			listUnlabel[countUnlabel++] = nodeCount;
	}

	int countRegions = 0;
	int countUnlabelAfterProbe = countUnlabel;
	if (countUnlabel) {
		// Fix up stage
		// Initialise mapping for probe
		int *mapping = new int[nNodes];
		for (int i = 0; i < nNodes; i++)
			mapping[i] = i * 2;

		if (contractIters > 0) {
			// Probe options
			typename QPBO<INTERNAL>::ProbeOptions options;
			options.C = infinite_edge_cost;
			options.dilation = 1;
			options.weak_persistencies = 1;
			options.iters = contractIters;

			// Run probe steps
		    graph->Probe(mapping, options);

		    /*
			int *new_mapping = new int[nNodes];
			for (int i = 0; i < contractIters; i++) {
			    graph->Probe(new_mapping, options);
			    graph->MergeMappings(nNodes, mapping, new_mapping);
			}
			delete new_mapping;

			// Solve for weak persistencies
			graph->ComputeWeakPersistencies();
			*/

			// Read out entire labelling again (as weak persistencies may have changed)
			countUnlabelAfterProbe = 0;
			for (int nodeCount = 0; nodeCount < firstNNodes; nodeCount++) {
				labelOutPtr[nodeCount] = (int32_t)graph->GetLabel(mapping[nodeCount]/2);
				if (labelOutPtr[nodeCount] < 0)
					listUnlabel[countUnlabelAfterProbe++] = nodeCount;
				else
					labelOutPtr[nodeCount] = (labelOutPtr[nodeCount] + mapping[nodeCount]) % 2;
			}
		}

		if (improveMethod != 1) {
			// Clump unlabelled nodes into independent regions
			int *regionMap = new int[countUnlabelAfterProbe];
			for (int nodeCount=0; nodeCount<countUnlabelAfterProbe; nodeCount++) {
				int regionId = graph->GetRegion(mapping[listUnlabel[nodeCount]] /2);
				// Compress the labelling to consecutive integers
				int region;
				for (region = 0; region < countRegions; region++) {
					if (regionMap[region] == regionId)
						goto skip_point;
				}
				regionMap[countRegions] = regionId;
				countRegions++;
skip_point:
				labelOutPtr[listUnlabel[nodeCount]] = -1 - (int32_t)region;
			}
			delete regionMap;
		}

        if (improveMethod) {
            // Improve
            int *improve_order = new int[countUnlabelAfterProbe];
            if (improveMethod == 2) {
                // Call the callback to determine the labels
                if (!mxIsCell(prhs[nrhs-1]) || mxGetNumberOfElements(prhs[nrhs-1]) != 2) {
                    delete improve_order;
                    delete mapping;
                    delete listUnlabel;
                    delete graph;
                    mexErrMsgTxt("options requires function name in second cell.");
                }
                mxArray *returnArray;
                mxArray *inputArrays[2];
                inputArrays[0] = mxGetCell(prhs[nrhs-1], 1);
                inputArrays[1] = plhs[0];
                int success = mexCallMATLAB(1, &returnArray, 2, inputArrays, "feval");
                if (success < 0 || !mxIsLogical(returnArray) || mxGetNumberOfElements(returnArray) != firstNNodes) {
                    delete improve_order;
                    delete mapping;
                    delete listUnlabel;
                    delete graph;
                    if (success == 0)
                        mxDestroyArray(returnArray);
                    mexErrMsgTxt("Callback fails or returns an array of unexpected type/size.");
                }
                const mxLogical *userDefinedLabel = (const mxLogical *)mxGetData(returnArray);
                // Set the labels to the user-defined value
                for (int nodeCount=0; nodeCount<countUnlabelAfterProbe; nodeCount++) {
                    improve_order[nodeCount] = mapping[listUnlabel[nodeCount]] / 2;
                    graph->SetLabel(improve_order[nodeCount], userDefinedLabel[improve_order[nodeCount]]);
                }
                mxDestroyArray(returnArray);
            } else {
                // Create new list of nodes to improve
                for (int nodeCount=0; nodeCount<countUnlabelAfterProbe; nodeCount++) {
                    // Set the labels to 0
                    improve_order[nodeCount] = mapping[listUnlabel[nodeCount]] / 2;
                    graph->SetLabel(improve_order[nodeCount], 0);
                }
            }
            
            // Randomize order
            for (int i = 0; i < countUnlabelAfterProbe-1; i++) {
                int j = i + (int)(((double)rand() / ((double)RAND_MAX+1)) * (countUnlabelAfterProbe - i));
                code_assert(j<countUnlabelAfterProbe);
                int k = improve_order[j];
                improve_order[j] = improve_order[i];
                improve_order[i] = k;
            }
            
            // Run QPBO-I
            graph->Improve(countUnlabelAfterProbe, improve_order);
            delete improve_order;
            
            // Read out the labels
            for (int nodeCount=0; nodeCount<countUnlabelAfterProbe; nodeCount++)
                labelOutPtr[listUnlabel[nodeCount]] = (int32_t)((graph->GetLabel(mapping[listUnlabel[nodeCount]]/2) + mapping[listUnlabel[nodeCount]]) % 2);
        }
		delete mapping;
	}
	delete listUnlabel;
	delete graph;

	if (nlhs > 1) {
		// Save statistics
		plhs[1] = mxCreateNumericMatrix(3, 1, mxINT32_CLASS, mxREAL);
		int32_t *stats = (int32_t *)mxGetData(plhs[1]);
		stats[0] = (int32_t)countUnlabel;
		stats[1] = (int32_t)countRegions;
		stats[2] = (int32_t)countUnlabelAfterProbe;
	}

	// mexPrintf("Unlabeled: %d, Non-sub %d out of %d (perc. %f) max Nodes: %d Max Edges: %d\n", countUnlabel, countIrregular, maxEdges, ((double)countIrregular/(double)maxEdges), maxNodes, maxEdges);  
	return;
}
