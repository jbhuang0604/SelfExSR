function A = slice_cube(Acube, D)
%SLICE_CUBE  Given depths, output a slice of a colour cube
%
%   A = slice_cube(Acube, D)
%
% Given a matrix of depth indices D, outputs the values in the array Acube
% at these depths.
%
% IN:
%   Acube - LxMxNxC or CxLxMxN cube of vectors (e.g. colour vectors)
%   D - MxN matrix of indices denoting which vector along the first
%       dimension of Acube should be output.
%
% OUT:
%   A - MxNxC output array

% $Id: slice_cube.m,v 1.1 2007/12/07 11:27:52 ojw Exp $

[d h w c] = size(Acube);
if d == 3 && c ~= 3 && c ~= 1
    % Colour is probably along the first dimension
    d = h;
    h = w;
    w = c;
    c = 3;
    offsets = uint32(d .* (repmat((0:h-1)', [1 w]) + (h .* repmat((0:w-1), [h 1]))));
    A = Acube(:,offsets + uint32(D))';
else
    offsets = uint32(d .* (repmat((0:h-1)', [1 w]) + (h .* repmat((0:w-1), [h 1]))));
    Acube = reshape(Acube, [d*h*w c]);
    A = Acube(offsets + uint32(D),:);
end
A = reshape(A, [h w c]);
return