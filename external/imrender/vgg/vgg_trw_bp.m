%VGG_TRW_BP  Multi-label MRF energy minimization using TRW-S & LBP
%
%   [L energy lower_bound] = vgg_trw_bp(UE, PI, PE, [options])
%
% Uses the message passing algorithms TRW-S or LBP to solve an MRF energy
% minimization problem with binary or multiple labels.
%
% This function uses mexified C++ code written by Vladimir Kolmogorov, and
% downloaded from Microsoft Research. See "Convergent Tree-Reweighted
% Message Passing for Energy Minimization", Kolmogorov, PAMI 2006 for
% details.
%
% IN:
%   UE - HxW cell array of unary terms for H*W nodes in the graph. Each
%        cell contains a KxL double matrix, with the first row being the L
%        energies for each of the L possible labels for that node. The
%        number of labels does not need to be constant across nodes.
%   PI - {2,3}xN uint32 matrix, each column containing the following
%        information on an edge: [start_node end_node
%        [pairwise_energy_table_index]]. If there are only 2 rows then
%        pairwise_energy_table_index is assumed to be the column number,
%        therefore you must have N == P.
%   PE - 1xP cell array of pairwise functions. Each cell gives the
%        a pairwise function of the form defined by options(1). If the form
%        is general (i.e. lookup table) then the cell conatins a L1xL2
%        matrix, where L1 and L2 are the number of labels for the start and
%        end nodes respectively.
%   options - 1x3 int32 vector of optional parameters:
%      UseTRW - 0: use LBP; otherwise: use TRW-S. Default: 1.
%      Type - Pairwise energy functional type. 0: general. 1: truncated
%             quadratic - min(a*(l1-l2)^2,b). 2: truncated linear -
%             min(a*abs(l1-l2),b). Easy to add support for others. Default:
%             0.
%      MaxIters - number of iterations of LBP or TRW to do. Default: 30.
%
% OUT:
%   L - HxW uint16 matrix of the energy minimizing state of each
%       node.
%   energy - scalar giving E(L).
%   lower_bound - scalar giving a lower bound on E(L) for any L. Only
%                 calculated by TRW-S (i.e. 0 for LBP).

% $Id: vgg_trw_bp.m,v 1.2 2008/03/10 18:45:27 ojw Exp $

function varargout = vgg_trw_bp(varargin)
funcName = mfilename;
sd = 'trw-s/';
sourceList = {['-I' sd], [funcName '.cxx'], [sd 'MRFEnergy.cpp'],...
              [sd 'minimize.cpp'], [sd 'ordering.cpp'],...
              [sd 'treeProbabilities.cpp']};
vgg_mexcompile_script; % Compilation happens in this script
return
