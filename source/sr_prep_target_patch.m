function trgPatch = sr_prep_target_patch(img, patchSize)

% SR_PREP_TARGET_PATCH
%
% Prepare target patches. Target patches are axis-aligned with sizes
% patchSize x patchSize. 
%
% Input:
%   - img:       input image
%   - patchSize: typically 5 or 7
% Output:
%   - trgPatch:  target patches [patchSize*patchSize] x [3] x [numUvPix] 

% =========================================================================

[imgH, imgW, nCh] = size(img);
numUvPix = (imgH - patchSize + 1)*(imgW - patchSize + 1);

% Initialization
trgPatch = zeros(patchSize*patchSize, 3, numUvPix, 'single');

% Get target patches using im2col
for i = 1 : nCh
    trgPatch(:,i,:) = im2col(img(:, :, i), [patchSize, patchSize], 'sliding');
end

end