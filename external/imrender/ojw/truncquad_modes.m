%TRUNCQUAD_MODES  Compute colour modes for a truncated quadratic kernel
%
%   [modes depth energy inliers] = truncquad_modes(I, thresh[, use_variance[, search_width]])
%
% Computes the colour modes for each set of input vectors, using the
% truncated quadratic kernel.
%
% IN:
%   I - CxLxM array of L input colour vectors with C channels at M depths.
%   thresh - scalar cost truncation threshold.
%   use_variance - 0: sum cost over all input vectors; 1: sum cost over
%                  inlying input vectors / num_inliers; 2: sum cost over
%                  inliers / (num_inliers - 1). Default: 0.
%   search_width - scalar indicating depths above and below current depth
%                  to search for modes within. Default: M.
%
% OUT:
%   modes - CxN matrix of colour modes.
%   depth - 1xN list of indices of depth for each mode.
%   energy - 1xN list of costs for each mode.
%   inliers - LxN logical array indicating inliers used to calculate the
%             colour of each mode.

% $Id: truncquad_modes.m,v 1.2 2008/11/17 11:27:35 ojw Exp $

function varargout = truncquad_modes(varargin)
funcName = mfilename;
sourceList = {[funcName '.cxx']}; % Cell array of source files
vgg_mexcompile_script; % Compilation happens in this script
return