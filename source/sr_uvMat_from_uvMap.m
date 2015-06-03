function uvMat = sr_uvMat_from_uvMap(uvMap, uvPixInd)
% SR_UVMAT_FROM_UVMAP: Convert uvMap into a matrix form
%
% Input:
%   - uvMap:    [imgH] x [imgW] x [nCh]
%   - uvPixInd: pixel index - [numUvPix] x [1]
% Output:
%   - uvMat:    the data in matrix form - [numUvPix] x [nCh]
% =========================================================================

[imgH, imgW, nCh] = size(uvMap);

% Get the offset vector
offset = uint64((0:nCh-1)*imgH*imgW);

% Adding the offset
uvPixInd = bsxfun(@plus, uvPixInd, offset);

% Retrieve data from uvMap
uvMat = uvMap(uvPixInd);

end