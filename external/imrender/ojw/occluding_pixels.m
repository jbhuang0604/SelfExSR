function P = occluding_pixels(V)
%OCCLUDING_PIXELS  Find occluding/occluded pixel pairs
%
%   P = occluding_pixels(V)
%
% Given a set of 3d image coordinates (pixel coordinates plus depth), finds
% any occluding/occluded pairs in this set.
%
%IN:
%   V - Mx3 list of 3d image points, defined as [x y Z], for which
%       occlusion interactions are to be found.
%
%OUT:
%   P - 2xN list of vectors of indices for [occluding occluded]' pixel
%       pairs.

% $Id: occluding_pixels.m,v 1.2 2007/12/14 17:06:23 ojw Exp $

[M N] = sortrows(V);
P = find_interactions(M, 0.5); % Optimized version
P = N(P);
return