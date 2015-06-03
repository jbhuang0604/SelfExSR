%FIND_INTERACTIONS  Given ordered 3d points, find occlusions pairs
%
%   P = find_interactions(V)
%
% Given a set of 3d image coordinates (pixel coordinates plus depth),
% ordered such that the x coordinates are monotonically increasing, as are
% the y coordinates within each block of identical x coordinates, finds any
% occluding/occluded pairs in this set.
%
%IN:
%   V - Mx3 list of 3d image points, defined as [x y Z], which has been
%       sorted using sortrows(V).
%
%OUT:
%   P - 2xN list of vectors of indices for [occluding occluded]' pixel
%       pairs.

% $Id: find_interactions.m,v 1.2 2007/12/14 17:06:23 ojw Exp $

function varargout = find_interactions(varargin)
funcName = mfilename;
sourceList = {[funcName '.cxx']}; % Cell array of source files
vgg_mexcompile_script; % Compilation happens in this script
return