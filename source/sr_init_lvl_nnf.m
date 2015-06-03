function NNF = sr_init_lvl_nnf(imgSize, NNF, modelPlane, opt)

% SR_INIT_LVL_NNF
%
% Initialize the nearest neighbor field for the current level
%
% Input:
%   - imgSize   : [imgHeight, imgWeight]
%   - NNF       : nearest neighbor field from previous level
%   - modelPlane: plane parameters
%   - opt       : parameters
% Output:
%   - NNF       : initialized nearest neighbor field for the current level
% =========================================================================
if(opt.iLvl == opt.origResLvl - 1)
    % Initialize the NNF for the coarest level using random sampling
    NNF = sr_init_nnf(imgSize, modelPlane, opt);
else
    % Initialize the NNF upsampling of NNF from previous level
    NNF = sr_upsample(imgSize, NNF, modelPlane, opt);
end

end

function NNF = sr_init_nnf(imgSize, modelPlane, opt)

% SR_INIT_NNF: Initialize the nearest neighbor field
% Input:
%     - holeMask: the binary mask indicate the hole pixels
%     - modelPlane: model of detected plane structures
%     - opt: parameters for synthesis process
% Output:
%   Image height and width
%     - NNF.imgH: image height
%     - NNF.imgW: image width
%   Precomputing pixel positions for PatchMatch algorithm
%     - NNF.uvPix: the patch positions containing hole pixels
%     - NNF.uvPixN: the neighboring patch positions of 4 connected pixels
%     - NNF.validPix: the patch positions containing known pixels
%     - NNF.trgPatchInd: the indice to retrieve target patches
%   Initialize NNF components:
%     - NNF.uvTform: the transformation specifying the source patch
%         - NNF.uvTform.data:   9 x numUvPix
%         - NNF.uvTform.map:    H x W x 9
%     - NNF.uvPlaneID: the plane assignment for unknown region
%         - NNF.uvPlaneID.data: 1 x numUvPix
%         - NNF.uvPlaneID.map:  H x W
%     - NNF.uvBias: the bias between source and target patches
%         - NNF.uvBias.data:    3 x numUvPix
%         - NNF.uvBias.map:     H x W x 3
%     - NNF.uvCost: the matching cost between source and target patches
%         - NNF.uvCost.data:    3 x numUvPix
%         - NNF.uvCost.map:     H x W x 3
% =========================================================================

% =========================================================================
% Initialize uvPix and uvPixN
% =========================================================================
% Initialize uvPix
NNF.imgH = imgSize(1);
NNF.imgW = imgSize(2);
NNF.uvPix = sr_get_uvpix([NNF.imgH, NNF.imgW], opt.pRad);

% Initialize uvPixN (the 4-connected neighbors for uvPix)
NNF.uvPixN = cell(4,1);
for i = 1: 4
    NNF.uvPixN{i}.sub = bsxfun(@minus, NNF.uvPix.sub, opt.propDir(i,:));
    NNF.uvPixN{i}.ind = uint64(sub2ind(imgSize, NNF.uvPixN{i}.sub(:,2), NNF.uvPixN{i}.sub(:,1)));
    NNF.uvPixN{i}.validInd = NNF.uvPix.mask(NNF.uvPixN{i}.ind);
end

% =========================================================================
% Initialize trgPatchInd: indices for target color patches
% =========================================================================
numPixel        = NNF.imgH*NNF.imgW;
pixIndMap       = reshape(1:numPixel, NNF.imgH, NNF.imgW);
indTrgPatch     = im2col(pixIndMap, [opt.pSize, opt.pSize], 'sliding');
NNF.trgPatchInd = uint64(cat(1, indTrgPatch, indTrgPatch+numPixel, ...
    indTrgPatch + 2*numPixel));

% =========================================================================
% Initialize uvPlaneID: indices for target color patches
% =========================================================================
% Initialize uvPlaneID with the fronto-parallel plane
NNF.uvPlaneID.data = modelPlane.numPlane*ones(NNF.uvPix.numUvPix, 1, 'uint8');
NNF.uvPlaneID.map  = modelPlane.numPlane*ones(imgSize, 'uint8');
% Accumulative plane probability
NNF.uvPlaneID.planeProbAcc = sr_prep_plane_prob_acc(modelPlane.planeProb, NNF.uvPix.ind);
% Number of planes
NNF.uvPlaneID.numPlane = modelPlane.numPlane;

% =========================================================================
% Initialize uvPlaneA (affine) and uvPlaneH (homography)
% =========================================================================
% Initialize uvPlaneA.data with inplace samples and scaling 2x
NNF.uvTformA.data = zeros(NNF.uvPix.numUvPix, 4, 'single');
NNF.uvTformA.data(:,1) = 2;
NNF.uvTformA.data(:,4) = 2;

NNF.uvTformA.map  = zeros(NNF.imgH, NNF.imgW, 4, 'single');
NNF.uvTformA.map  = sr_update_uvMap(NNF.uvTformA.map, NNF.uvTformA.data, NNF.uvPix.ind);

% Homograph transformation
NNF.uvTformH.data = sr_src_domain_tform(NNF.uvPlaneID.data, modelPlane, NNF.uvTformA.data, NNF.uvPix.sub, NNF.uvPix.sub);
NNF.uvTformH.map  = zeros(NNF.imgH, NNF.imgW, 9, 'single');
NNF.uvTformH.map  = sr_update_uvMap(NNF.uvTformH.map, NNF.uvTformH.data, NNF.uvPix.ind);

% =========================================================================
% Initialize bias for patch matching
% =========================================================================
NNF.uvBias.data   = zeros(1, 3, NNF.uvPix.numUvPix, 'single');
NNF.uvBias.map    = zeros(NNF.imgH, NNF.imgW, 3, 'single');

% =========================================================================
% Initialize cost for patch matching
% =========================================================================
NNF.uvCost.data   = zeros(NNF.uvPix.numUvPix, 1, 'single');
NNF.uvCost.map    = zeros(NNF.imgH, NNF.imgW,    'single');

% %% === Initialize update ===
% NNF.update.data = false(NNF.uvPix.numUvPix, 1);
% NNF.update.map  = false(NNF.imgH, NNF.imgW);
%
% %% === Initialize uvPixUpdateSrc ===
% NNF.uvPixUpdateSrc.data = zeros(1, NNF.uvPix.numUvPix);
% NNF.uvPixUpdateSrc.map  = zeros(NNF.imgH, NNF.imgW);

% =========================================================================
% Initialize patch matching weights
% =========================================================================
wPatch      = fspecial('gauss', [opt.pSize, opt.pSize], 3);
NNF.wPatch  = wPatch(:);
NNF.wSumImg = zeros(NNF.imgH, NNF.imgW, 'single');
for i = 1 : NNF.uvPix.numUvPix
    NNF.wSumImg(indTrgPatch(:,i)) = NNF.wSumImg(indTrgPatch(:,i)) + NNF.wPatch;
end

end


function NNF_H = sr_upsample(imgSize, NNF_L, modelPlane, opt)

% SC_UPSAMPLE: upsample the nearest neighbor field

% Input:
%     - holeMask: the binary mask indicate the hole pixels
%     - NNF_L: nearest neighbor field of the low resolution image
%     - optS: parameters for synthesis process
% Output:
%   Image height and width
%     - NNF_H.imgH: image height
%     - NNF_H.imgW: image width
%   Precomputing pixel positions for PatchMatch algorithm
%     - NNF_H.uvPix: the patch positions containing hole pixels
%     - NNF_H.uvPixN: the neighboring patch positions of 4 connected pixels
%     - NNF_H.validPix: the patch positions containing known pixels
%     - NNF_H.trgPatchInd: the indice to retrieve target patches
%   Initialize NNF components:
%     - NNF_H.uvTform: the transformation specifying the source patch
%         - NNF.uvTform.data:   9 x numUvPix
%         - NNF.uvTform.map:    H x W x 9
%     - NNF_H.uvPlaneID: the plane assignment for unknown region
%         - NNF.uvPlaneID.data: 1 x numUvPix
%         - NNF.uvPlaneID.map:  H x W
%     - NNF_H.uvBias: the bias between source and target patches
%         - NNF.uvBias.data:    3 x numUvPix
%         - NNF.uvBias.map:     H x W x 3
%     - NNF_H.uvCost: the matching cost between source and target patches
%         - NNF.uvCost.data:    3 x numUvPix
%         - NNF.uvCost.map:     H x W x 3

% =========================================================================
% Initialize uvPix and uvPixN
% =========================================================================
% Initialize uvPix
NNF_H.imgH = imgSize(1);  NNF_H.imgW = imgSize(2);
% NNF_H.uvPix = sr_init_level([NNF_H.imgH, NNF_H.imgW], optS.pRad);
NNF_H.uvPix = sr_get_uvpix([NNF_H.imgH, NNF_H.imgW], opt.pRad);

% Initialize uvPixN (the 4-connected neighbors for uvPix)
NNF_H.uvPixN = cell(4,1);
for i = 1: 4
    NNF_H.uvPixN{i}.sub = bsxfun(@minus, NNF_H.uvPix.sub, opt.propDir(i,:));
    NNF_H.uvPixN{i}.ind = uint64(sub2ind([NNF_H.imgH, NNF_H.imgW], NNF_H.uvPixN{i}.sub(:,2), NNF_H.uvPixN{i}.sub(:,1)));
    NNF_H.uvPixN{i}.validInd = NNF_H.uvPix.mask(NNF_H.uvPixN{i}.ind);
end

% =========================================================================
% Initialize trgPatchInd: indices for target color patches
% =========================================================================
pixIndMap = reshape(1:NNF_H.imgH*NNF_H.imgW, NNF_H.imgH, NNF_H.imgW);
indTrgPatch = im2col(pixIndMap, [opt.pSize, opt.pSize]);

NNF_H.trgPatchInd = cat(1, indTrgPatch, ...
    indTrgPatch+NNF_H.imgH*NNF_H.imgW, indTrgPatch + 2*NNF_H.imgH*NNF_H.imgW);
NNF_H.trgPatchInd = uint64(NNF_H.trgPatchInd);

% =========================================================================
% Initialize uvPixL: the correspondence between high-res and low-res pixels
% =========================================================================
imgH_H = NNF_H.imgH;    imgW_H = NNF_H.imgW;
imgH_L = NNF_L.imgH;    imgW_L = NNF_L.imgW;

sX = imgH_L/imgH_H;     sY = imgW_L/imgW_H;
uvPixL.sub = round(NNF_H.uvPix.sub*diag([sX, sY]));
uvPixL.sub(:,1) = sr_clamp(uvPixL.sub(:,1), opt.pRad+1, imgW_L - opt.pRad);
uvPixL.sub(:,2) = sr_clamp(uvPixL.sub(:,2), opt.pRad+1, imgH_L - opt.pRad);
uvPixL.ind = uint64(sub2ind([imgH_L, imgW_L], uvPixL.sub(:,2), uvPixL.sub(:,1)));

% =========================================================================
% Initialize uvPlaneID: indices for target color patches
% =========================================================================

NNF_H.uvPlaneID.data = sr_uvMat_from_uvMap(NNF_L.uvPlaneID.map, uvPixL.ind);
NNF_H.uvPlaneID.data(NNF_H.uvPlaneID.data==0) = 1;
NNF_H.uvPlaneID.map = zeros(NNF_H.imgH, NNF_H.imgW, 'uint8');
NNF_H.uvPlaneID.map = sr_update_uvMap(NNF_H.uvPlaneID.map, NNF_H.uvPlaneID.data, NNF_H.uvPix.ind);
NNF_H.uvPlaneID.numPlane = modelPlane.numPlane;

planeProbAccData = sr_prep_plane_prob_acc(modelPlane.planeProb, NNF_H.uvPix.ind);

NNF_H.uvPlaneID.planeProbAcc = planeProbAccData;
NNF_H.uvPlaneID.mLogLikelihood = -log(planeProbAccData);

% =========================================================================
% Initialize uvPlaneA (affine) and uvPlaneH (homography)
% =========================================================================

uvTformH_L = sr_uvMat_from_uvMap(NNF_L.uvTformH.map, uvPixL.ind);
uvTformH_L(:,7:8) = uvTformH_L(:,7:8)*diag([1/sX, 1/sY]);

% Refinement
refineVec = NNF_H.uvPix.sub - uvPixL.sub*diag([1/sX, 1/sY]);
uvTformH_H = sr_trans_tform(uvTformH_L, refineVec);

% Clamp
uvTformH_H(:,7) = sr_clamp(uvTformH_H(:,7), opt.pRad+1, imgW_H - opt.pRad);
uvTformH_H(:,8) = sr_clamp(uvTformH_H(:,8), opt.pRad+1, imgH_H - opt.pRad);
% uvValid_H = sr_check_valid_uv(uvTform_H(7:8,:), NNF_H.validPix.mask);

% uvInvalidInd = ~uvValid_H;
% nInvalidUv_H = sum(uvInvalidInd);
% if(nInvalidUv_H)
%     randInd = randi(size(NNF_H.validPix.ind, 2), nInvalidUv_H, 1);
%     uvRand = NNF_H.validPix.sub(:, randInd);
%     uvTform_H(7,uvInvalidInd) = uvRand(1,:);
%     uvTform_H(8,uvInvalidInd) = uvRand(2,:);
% end

% Update uvTform.map
NNF_H.uvTformH.data = uvTformH_H;
I = reshape(eye(3), 1, 1, 9);
NNF_H.uvTformH.map = repmat(I, [imgH_H, imgW_H, 1]);
NNF_H.uvTformH.map = sr_update_uvMap(NNF_H.uvTformH.map, NNF_H.uvTformH.data, NNF_H.uvPix.ind);

% Update uvTformA
uvTformA_L = sr_uvMat_from_uvMap(NNF_L.uvTformA.map, uvPixL.ind);
NNF_H.uvTformA.data = uvTformA_L;
NNF_H.uvTformA.map  = zeros(NNF_H.imgH, NNF_H.imgW, 4, 'single');
NNF_H.uvTformA.map  = sr_update_uvMap(NNF_H.uvTformA.map, NNF_H.uvTformA.data, NNF_H.uvPix.ind);

% =========================================================================
% Initialize bias for patch matching
% =========================================================================
uvBias   = sr_uvMat_from_uvMap(NNF_L.uvBias.map, uvPixL.ind);
NNF_H.uvBias.map    = zeros(imgH_H, imgW_H, 3, 'single');
NNF_H.uvBias.map    = sr_update_uvMap(NNF_H.uvBias.map, uvBias, NNF_H.uvPix.ind);
NNF_H.uvBias.data   = reshape(uvBias, 1, 3, NNF_H.uvPix.numUvPix);

% =========================================================================
% Initialize cost for patch matching
% =========================================================================
NNF_H.uvCost.map  = zeros(imgH_H, imgW_H, 'single');
NNF_H.uvCost.data = sr_uvMat_from_uvMap(NNF_L.uvCost.map, uvPixL.ind);
NNF_H.uvCost.map  = sr_update_uvMap(NNF_H.uvCost.map, NNF_H.uvCost.data, NNF_H.uvPix.ind);

% =========================================================================
% Initialize patch weights
% =========================================================================
NNF_H.wPatch = fspecial('gauss', [opt.pSize, opt.pSize], 3);
NNF_H.wPatch = NNF_H.wPatch(:);
NNF_H.wSumImg = zeros(NNF_H.imgH, NNF_H.imgW, 'single');
indMap = reshape(1:NNF_H.imgH*NNF_H.imgW, NNF_H.imgH, NNF_H.imgW);
indPatch = im2col(indMap, [opt.pSize, opt.pSize]);

for i = 1: size(indPatch, 2)
    NNF_H.wSumImg(indPatch(:,i)) = NNF_H.wSumImg(indPatch(:,i)) + NNF_H.wPatch(:);
end


end

function uvPix = sr_get_uvpix(imgSize, prad)
%  Get uvPix: a field of patch centers for patch-based synthesis
%
%  Input:
%   - imgSize: [imgHeight, imgWidth]
%   - prad   : radius of a patch, prad = floor(patchSize/2);
%  Output:
%   - uvPix  : patch center positions
% =========================================================================

% Get uvMap
uvMap = true(imgSize);
uvMap([1:prad, end-prad+1:end], :) = 0;
uvMap(:, [1:prad, end-prad+1:end]) = 0;

[Y, X] = find(uvMap);

% Get uvPix in terms of pos, ind, and mask
uvPix.sub  = single(cat(2, X(:), Y(:)));
uvPix.ind  = uint64(sub2ind(imgSize, Y(:), X(:)));
uvPix.mask = uvMap;
uvPix.numUvPix = size(uvPix.ind, 1);

end