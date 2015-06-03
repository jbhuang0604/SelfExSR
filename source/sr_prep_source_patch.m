function [srcPatch, srcPatchScale] = sr_prep_source_patch(imgPyr, uvTform, opt)

% SC_PREP_SOURCE_PATCH
%
% Prepare source patches according to the patch transformation uvTform
%
% Input:
%   - imgPyr:  image pyramid
%   - uvTform: homography transformation
%   - opt:     parameters
% Output:
%   - srcPatch
% =========================================================================

% Estimate the scale of the transformation
srcPatchScale = sr_scale_tform(uvTform);

% Find the closest scales in the image pyramid
uvScaleLvlInd = round(-log(srcPatchScale)/log(opt.alpha));

% clamp the scale
uvScaleLvlInd = sr_clamp(uvScaleLvlInd, 1, opt.nPyrLowLvl);

% Update the uvTform with the appropriate scales
srcPatchScaleQ = opt.alpha.^(uvScaleLvlInd);
uvTform(:, [1,4,7,2,5,8]) = bsxfun(@times, uvTform(:, [1,4,7,2,5,8]), srcPatchScaleQ);

% Prepare source patch sampling position
numUvPix = size(uvTform, 1);
srcPatch = zeros(opt.pNumPix, numUvPix, 3, 'single');

for iLvlCur = 1: opt.nPyrLowLvl
    scaleCurInd = uvScaleLvlInd == iLvlCur;
    nSrcPatchScaleCur = sum(scaleCurInd);
    if(nSrcPatchScaleCur)
        % Get low-res image with the closet scale
        img = imgPyr{opt.iLvl + iLvlCur};
        
        % Get the source patch positions
        uvTformCur = uvTform(scaleCurInd, :);
        
        % Get srcPatchPos
        c1 = reshape(uvTformCur(:,1:3)', 1, 3, nSrcPatchScaleCur);
        c2 = reshape(uvTformCur(:,4:6)', 1, 3, nSrcPatchScaleCur);
        c3 = reshape(uvTformCur(:,7:9)', 1, 3, nSrcPatchScaleCur);
        
        % Get the source patch pixel positions
        srcPatchPos = bsxfun(@times, opt.refPatchPos(:,1), c1) + bsxfun(@times, opt.refPatchPos(:,2), c2);
        srcPatchPos = bsxfun(@plus, srcPatchPos, c3);
        
        % Convert back to Eucledian coordinate
        indH = srcPatchPos(opt.pMidPix, 3, :) ~= 1;
        srcPatchPos(:, 1:2,indH) = bsxfun(@rdivide, srcPatchPos(:, 1:2, indH), srcPatchPos(:, 3, indH));
      
        % Grab the color values of source patch using bilinear interpolation
        srcPatch(:, scaleCurInd, :) = vgg_interp2(img, srcPatchPos(:,1,:), srcPatchPos(:,2,:), 'linear', 0);
    end
end

srcPatch = permute(srcPatch, [1,3,2]);

end