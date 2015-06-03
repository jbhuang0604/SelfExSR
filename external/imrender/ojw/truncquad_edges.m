%TRUNCQUAD_EDGES  Compute truncated quadratic pairwise texture costs
%
%   edge_costs = truncquad_edges(I, modes, EI, thresh)
%
% Computes the truncated quadratic cost of each combination of pairs of
% modes in each edge, based on the minimum distance to a pair of pixels in
% the corresponding library.
%
% IN:
%   I - CxLxN array of colour libraries, where C is the number of colour
%       channels, L the length of the library and N the number of
%       libraries.
%   modes - 1xN cell array of [CxM] containing modes, one set per library.
%   EI - 2xP matrix of indices for modes/libaries belonging to each edge.
%   thresh - scalar cost truncation threshold.
%
% OUT:
%   edge_costs - Px1 cell array of edge cost matrices.

% $Id: truncquad_edges.m,v 1.2 2008/11/17 11:27:35 ojw Exp $

function varargout = truncquad_edges(varargin)
funcName = mfilename;
sourceList = {[funcName '.cxx']}; % Cell array of source files
vgg_mexcompile_script; % Compilation happens in this script
return