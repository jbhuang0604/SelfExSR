%VGG_QPBO  Binary MRF energy minimization on non-submodular graphs
%
%   [L stats] = vgg_qpbo(UE, PI, PE, [TI, TE], [options])
%
% Uses the Quadratic Pseudo-Boolean Optimization (QPBO - an extension of
% graph cuts that solves the "roof duality" problem, allowing graphs with
% submodular edges to be solved) to solve binary, pairwise and triple
% clique MRF energy minimization problems.
%
% The algorithm can be viewed as fusing two labellings in such a way as to
% minimize output labellings energy:
%   L = FUSE(L0, L1)
% In particular the solution has these properties:
%   P1: Not all nodes are labelled. Nodes which are labelled form part of
%       an optimal labelling.
%   P2: If all unlabelled nodes in L are fixed to L0, then E(L) <= E(L0).
%
% This function also implements QPBO-P and QPBO-I. See "Optimizing binary
% MRFs via extended roof duality", Rother et al., CVPR 2007, for more
% details.
%
% This function uses mexified C++ code downloaded from www page of Vladimir
% Kolmogorov. To acknowledge him/reference the correct papers, see his web
% page for details.
%
% IN:
%   UE - 2xM matrix of unary terms for M nodes, or M = UE(1) if numel(UE)
%        == 1 (assumes all unary terms are zero). UE can be any type, but
%        PE and TE must be of the same type; also, optimality is not
%        guaranteed for non-integer types.
%   PI - {2,3}xN uint32 matrix, each column containing the following
%        information on an edge: [start_node end_node
%        [pairwise_energy_table_index]]. If there are only 2 rows then
%        pairwise_energy_table_index is assumed to be the column number,
%        therefore you must have N == P.
%   PE - 4xP pairwise energy table, each column containing the energies
%        [E00 E01 E10 E11] for a given edge.
%   TI - {3,4}xQ uint32 matrix, each column containing the following
%        information on a triple clique: [node1 node2 node3
%        [triple_clique_energy_table_index]]. If there are only 3 rows then
%        triple_clique_energy_table_index is assumed to be the column
%        number, therefore you must have Q == R.
%   TE - 8xR triple clique energy table, each column containing the
%        energies [E000 E001 E010 E011 E100 E101 E110 E111] for a given
%        triple clique.
%   options - 2x1 cell array. options{1} is a 1x4 uint32 vector of optional
%             parameters: 
%      FirstNNodes - Only labels of nodes 1:FirstNNodes will be output.
%                    Also limits nodes considered in QPBO-P and QPBO-I.
%                    Default: M.
%      ImproveMethod - 0: Unlabelled nodes not improved; 1: QPBO-I
%                      employed, assuming preferred label is L0. 2: QPBO-I
%                      employed, using user-defined labelling returned by
%                      callback_func.
%      ContractIters - Number of iterations of QPBO-P to do. Default: 0.
%      LargerInternal - 0: input type also used internally; 1: code uses a
%                       larger representation (i.e. more bits) for energies
%                       than that provided, reducing the possibility of
%                       edges becoming saturated. Default: 0.
%             options{2} - callback_func, a function name or handle which,
%                          given the labelling L (the first output of this
%                          function), returns a logical array of the same
%                          size defining a labelling to be transfromed by
%                          QPBO-I.
%
% OUT:
%   L - 1xFirstNNodes int32 matrix of the energy minimizing state of each
%       node. Labels {0,1} are the optimal label (in a weak sense), while
%       labels < 0 indicate unlabelled nodes. For unlabelled nodes, if L(x)
%       == L(y) then x and y are in the same "strongly connected region".
%   stats - 3x1 int32 vector of stats:
%      stats(1) - number of unlabelled pixels
%      stats(2) - number of strongly connected regions of unlabelled pixels
%      stats(3) - number of unlabelled pixels after running QPBO-P

% $Id: vgg_qpbo.m,v 1.1 2007/12/07 11:27:56 ojw Exp $

function varargout = vgg_qpbo(varargin)
funcName = mfilename;
sd = 'qpbo/';
sourceList = {['-I' sd], [funcName '.cxx'], [sd 'QPBO.cpp'],...
              [sd 'QPBO_maxflow.cpp'], [sd 'QPBO_extra.cpp'],...
              [sd 'QPBO_postprocessing.cpp']};
vgg_mexcompile_script; % Compilation happens in this script
return
