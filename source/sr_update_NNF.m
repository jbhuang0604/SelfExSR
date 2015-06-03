% === Nearest neighbor field update ===
function [NNF, nUpdate]= sr_update_NNF(trgPatch, imgPyr, NNF, modelPlane, opt)
%
% SR_UPDATE_NNF:
%
% One iteration of random search and spatial propagation (generalized PatchMatch)
% for nearest neighbor field estimation
%
% Input:
%   - trgPatch:   target patches
%   - imgPyr:     source image pyramid
%   - NNF:        the current nearest neighbor field
%   - modelPlane: planar structure model
%   - opt:        parameters
% Output:
%   - NNF:        updated nearest neighbor field
%   - nUpdate:    nUpdate(1) = #updates from propagation
%                 nUpdate(2) = #updates from random sampling
% =========================================================================

nUpdate = zeros(1,2);

NNF.update.data = false(NNF.uvPix.numUvPix, 1);
NNF.update.map  = false(NNF.imgH, NNF.imgW);

% Random sampling
[NNF, n] = sr_random_search(trgPatch, imgPyr, NNF, modelPlane, opt);
nUpdate(2) = nUpdate(2) + n;

% propagate along four directions
for iDirect = 1:4
    [NNF, n] = sr_propagate(trgPatch, imgPyr, NNF, modelPlane, opt, iDirect);
    nUpdate(1) = nUpdate(1) + n;
end
end

% === Random search ===
function  [NNF, nUpdateTotal] = sr_random_search(trgPatch, imgPyr, NNF, modelPlane, opt)
%
% SC_RANDOM_SEARCH:
%
% Update the nearest neighbor field using coarse-to-fine random search
%
% Input:
%   - trgPatch:   target patches
%   - imgPyr:     source image pyramid
%   - NNF:        the current nearest neighbor field
%   - modelPlane: planar structure model
%   - opt:        parameters
% Output:
%   - NNF:        updated nearest neighbor field
%   - nUpdateTotal:#updates
% =========================================================================

imgSize  = [NNF.imgH, NNF.imgW];
uvPix    = NNF.uvPix;
numUvPix = uvPix.numUvPix;

searchPosRad = max(imgSize)/4;
nUpdateTotal = 0;

iter = 1;

while(searchPosRad > 1) % Coarse to fine random sampling
    iter = iter + 1;
    % =====================================================================
    % Draw random sample: srcPos, uvTformA, uvTformH, uvPlaneID
    % =====================================================================
    % Reduce search radius by half
    searchPosRad = searchPosRad/2;
    
    % Retrieve the current source patch position and affine transformation
    uvTformHCand = sr_uvMat_from_uvMap(NNF.uvTformH.map, uvPix.ind);
    uvTformACand = sr_uvMat_from_uvMap(NNF.uvTformA.map, uvPix.ind);
    
    % Draw random samples of (position: srcPos) and affine deformation (scale, rotation : uvTformD)
    [srcPosOffset, uvTformD] = sr_draw_rand_sample(searchPosRad, numUvPix, iter, opt);
    
    % Apply the position offset
    srcPos = uvTformHCand(:,7:8) + srcPosOffset;
    % Apply the affine transformation deformation
    uvTformACand = sr_apply_affine_tform(uvTformACand, uvTformD);
    
    % Draw plane ID candidate
    if(opt.usePlaneGuide)
        uvPlaneIDCand = sr_draw_plane_id(NNF.uvPlaneID.planeProbAcc);
    else
        uvPlaneIDCand = modelPlane.numPlane*ones(1, uvPix.numUvPix);
    end
    
    % Estimate the domain transformation using the updated source positions
    % and affine transformation
    uvTformHCand = sr_src_domain_tform(uvPlaneIDCand, modelPlane, uvTformACand, srcPos, NNF.uvPix.sub);
    
    % =====================================================================
    % Reject invalid samples
    % =====================================================================
    % Check if the scale of the source patch is valid
    uvTformHScale = sr_scale_tform(uvTformHCand);
    uvValidScaleInd = (uvTformHScale >= opt.minScale) & (uvTformHScale <= opt.maxScale);
    
    % Check if the source patch position is valid
    uvValidPosInd = sr_check_valid_pos(uvTformHCand(:,7:8), imgSize, opt.pRad);

    % Check if the errors are already below a threshold
    uvValidErrInd = NNF.uvCost.data > opt.errThres;
    
    % Valid uv pixels
    uvValidInd = uvValidPosInd & uvValidScaleInd & uvValidErrInd;

    uvValidPos = find(uvValidInd);
    
    if(size(uvValidPos, 1) ~= 0)
        % =====================================================================
        % Check if the randomly drawed samples reduce the matching cost
        % =====================================================================
        % Prepare valid samples
        trgPatchCur      = trgPatch(:,:,uvValidInd);
        uvCostDataCur    = NNF.uvCost.data(uvValidInd);
        uvTformHCandCur  = uvTformHCand(uvValidInd, :);
        uvTformACandCur  = uvTformACand(uvValidInd, :);
        uvPlaneIDCandCur = uvPlaneIDCand(uvValidInd);
        
        uvPixValid.sub   = uvPix.sub(uvValidInd, :);
        uvPixValid.ind   = uvPix.ind(uvValidInd);
        
        % Prepare source patches
        [srcPatch, srcPatchScale] = sr_prep_source_patch(imgPyr, uvTformHCandCur, opt);
        % Compute patch matching cost: appearance cost
        [costPatchCand, uvBiasCand] = sr_patch_cost_app(trgPatchCur, srcPatch, opt);
        % Compute patch matching cost: scale cost
        if(opt.useScaleCost)
            costScale = opt.lambdaScale*max(0, opt.scaleThres - srcPatchScale);
            costPatchCand = costPatchCand + costScale;
        end
        % Compute patch matching cost: plane compatibility cost
        if(opt.usePlaneGuide)
            costPlane = sr_patch_cost_plane(modelPlane.mLogLPlaneProb, uvPlaneIDCandCur, uvPixValid.ind, srcPos(uvValidInd,:));
            costPatchCand = costPatchCand + opt.lambdaPlane*costPlane;
        end
        
        % Check which one to update
        updateInd = (costPatchCand < uvCostDataCur);
        
        % =====================================================================
        % Update the nearest neighbor field
        % =====================================================================
        nUpdate = sum(updateInd);
        if(nUpdate~=0)
            uvUpdatePos = uvValidPos(updateInd);
            
            % === Update NNF data ===
            NNF.uvTformH.data(uvUpdatePos, :) = uvTformHCandCur(updateInd, :);
            NNF.uvTformA.data(uvUpdatePos, :) = uvTformACandCur(updateInd, :);
            NNF.uvPlaneID.data(uvUpdatePos)   = uvPlaneIDCandCur(updateInd);
            NNF.uvCost.data(uvUpdatePos)      = costPatchCand(updateInd);
            if(opt.useBiasCorrection)
                uvBiasCand = uvBiasCand(:, :, updateInd);
                NNF.uvBias.data(:,:,uvUpdatePos) = uvBiasCand;
            end
            NNF.update.data(uvUpdatePos) = 1;
            
            % === Update NNF map ===
            uvPixValidInd = uvPixValid.ind(updateInd);
            NNF.uvTformH.map   = sr_update_uvMap(NNF.uvTformH.map,  uvTformHCandCur(updateInd,:), uvPixValidInd);
            NNF.uvTformA.map   = sr_update_uvMap(NNF.uvTformA.map,  uvTformACandCur(updateInd,:), uvPixValidInd);
            NNF.uvPlaneID.map  = sr_update_uvMap(NNF.uvPlaneID.map, uvPlaneIDCandCur(updateInd),  uvPixValidInd);
            NNF.uvCost.map     = sr_update_uvMap(NNF.uvCost.map,    costPatchCand(updateInd),     uvPixValidInd);
            if(opt.useBiasCorrection)
                uvBiasCand     = squeeze(uvBiasCand)';
                NNF.uvBias.map = sr_update_uvMap(NNF.uvBias.map, uvBiasCand, uvPixValidInd);
            end
            NNF.update.map  = sr_update_uvMap(NNF.update.map, 1, uvPixValidInd);
            nUpdateTotal = nUpdateTotal + nUpdate;
            
        end
    end
end

end


function uvTformACand = sr_apply_affine_tform(uvTformA, uvTformD)

uvTformACand = zeros(size(uvTformA), 'single');
uvTformACand(:, 1:2) = bsxfun(@times, uvTformD(:,1:2), uvTformA(:,1)) + ...
    bsxfun(@times, uvTformD(:,3:4), uvTformA(:,2));
uvTformACand(:, 3:4) = bsxfun(@times, uvTformD(:,1:2), uvTformA(:,3)) + ...
    bsxfun(@times, uvTformD(:,3:4), uvTformA(:,4));
end

function validPosInd = sr_check_valid_pos(pos, imgSize, prad)

imgH = imgSize(1);
imgW = imgSize(2);

validPosInd = (pos(:,1) <= imgW - prad) & (pos(:,1) >= prad + 1) & ...
    (pos(:,2) >= prad + 1 ) & (pos(:,2) <= imgH - prad);

end

function [srcPosOffset, uvTformD] = sr_draw_rand_sample(searchPosRad, numUvPix, iter, opt)
%
% SR_DRAW_RAND_SAMPLE
%
% Nearest neighbor field estimation at the current level using the
% generalized PatchMatch algorithm
%
% Input:
%   - searchPosRad: search radius for position search
%   - numUvPix:     number of target pixels
%   - iter:         iteration (used for narrowing down the search space)
%  -  opt:          parameters
% Output:
%   - srcPosOffset: source patch position offset
%   - uvTformD:     affine perturbation matrix
% =========================================================================

% Positional offset
srcPosOffset = 2*searchPosRad*(rand(numUvPix, 2) - 0.5);

% Affine transformation
scale = opt.scaleRadA*(rand(numUvPix,1) - 0.5)/iter;  % scale
scale = scale + 1;                                    % perturb around one
theta = opt.rotRadA*  (rand(numUvPix,1) - 0.5)/iter;  % rot
sh_x  = opt.shRadA*   (rand(numUvPix,1) - 0.5)/iter;  % sheer x
sh_y  = opt.shRadA*   (rand(numUvPix,1) - 0.5)/iter;  % sheer y

% Create a affine perturbation matrix
uvTformD = zeros(numUvPix, 4, 'single');
uvTformD(:,1) = cos(theta) - sin(theta).*sh_y;
uvTformD(:,2) = sin(theta) + cos(theta).*sh_y;
uvTformD(:,3) = cos(theta).*sh_x - sin(theta);
uvTformD(:,4) = sin(theta).*sh_x + cos(theta);

uvTformD = bsxfun(@times, uvTformD, scale);

end

% === Propagation ===

function [NNF, nUpdateTotal] = sr_propagate(trgPatch, imgPyr, NNF, modelPlane, opt, indDirection)
%
% SR_PROPAGATE: update the nearest neighbor field using propagation
%
% Update the nearest neighbor field using coarse-to-fine random search
%
% Input:
%   - trgPatch:     target patches
%   - imgPyr:       source image pyramid
%   - NNF:          the current nearest neighbor field
%   - modelPlane:   planar structure model
%   - indDirection: direction of propagation
%   - opt:          parameters
% Output:
%   - NNF:        updated nearest neighbor field
%   - nUpdateTotal:#updates
% =========================================================================

nUpdateTotal = 0;

% The positions of neighboring pixels
uvPixN = NNF.uvPixN{indDirection};
% Check if the patch cost is already low
uvValidCostInd = NNF.uvCost.data > opt.errThres;
% The valid patch indices to start with
uvValidInd =  uvPixN.validInd & NNF.update.map(uvPixN.ind) & uvValidCostInd;

numUpdatePix = sum(uvValidInd);

while(numUpdatePix ~= 0)
    numUpdatePix  = 0; % exit the loop if there are no valid candidates
    % =====================================================================
    % Propare candidate patch transformation
    % =====================================================================
    % Prepare uvPix, uvPixNCur
    uvPix.sub     = NNF.uvPix.sub(uvValidInd, :);
    uvPix.ind     = NNF.uvPix.ind(uvValidInd);
    
    uvPixNCur.sub = uvPixN.sub(uvValidInd, :);         % neighbor of uvPix
    uvPixNCur.ind = uvPixN.ind(uvValidInd);
    
    trgPatchCur   = trgPatch(:,:, uvValidInd);          % target patch
    srcPosCur     = NNF.uvTformH.data(uvValidInd, 7:8); % source patch pos
    uvCostCur     = NNF.uvCost.data(uvValidInd);        % current patch matching cost
    uvPlaneIDCur  = NNF.uvPlaneID.map(uvPixNCur.ind);   % current plane ID
    
    uvValidPos    = find(uvValidInd);                   % Valid pixel positions
    
    % Get candidate uvTform candidates
    uvTformACand = sr_uvMat_from_uvMap(NNF.uvTformA.map, uvPixNCur.ind);
    uvTformHCand = sr_uvMat_from_uvMap(NNF.uvTformH.map, uvPixNCur.ind);
    srcPos       = uvTformHCand(:,7:8);
    
    % Generate candidate transformation by propagation
    uvTformHCand = sr_trans_tform(uvTformHCand, opt.propDir(indDirection,:));
    
    % =====================================================================
    % Reject invalid samples
    % =====================================================================
    % Check if the nearest neighbors are valid source patches
    uvValidSrcInd = sr_check_valid_pos(srcPos, [NNF.imgH, NNF.imgW], opt.pRad);
    
    % Check if the nearest neighbors are already the same as the existing one
    diff = abs(srcPos - srcPosCur);
    uvValidDistInd = ((diff(:,1) > 1 ) | (diff(:,2) > 1 ));
    
    % Check if the errors are already below a threshold
    uvValidErrInd = uvCostCur > opt.errThres;

    % Valid pixel indices
    uvValidInd = uvValidSrcInd & uvValidDistInd & uvValidErrInd;
    numUvValid = sum(uvValidInd);
    
    if(numUvValid ~= 0)
        % =====================================================================
        % Check if the propagated samples reduce the matching cost
        % =====================================================================
        trgPatchCur     = trgPatchCur(:,:, uvValidInd);
        uvCostCur       = uvCostCur(uvValidInd);
        uvTformHCand    = uvTformHCand(uvValidInd, :);
        uvTformACand    = uvTformACand(uvValidInd, :);
        uvPlaneIDCand   = uvPlaneIDCur(uvValidInd);
        
        uvValidPos      = uvValidPos(uvValidInd);
        uvPixValid.sub  = uvPix.sub(uvValidInd,:);
        uvPixValid.ind  = uvPix.ind(uvValidInd);
        
        % Grab source patches
        [srcPatch, srcPatchScale] = sr_prep_source_patch(imgPyr, uvTformHCand, opt);
        % Compute patch matching cost: appearance cost
        [costPatchCand, uvBiasCand] = sr_patch_cost_app(trgPatchCur, srcPatch, opt);
        % Compute patch matching cost: scale cost
        if(opt.useScaleCost)
            costScale = opt.lambdaScale*max(0, opt.scaleThres - srcPatchScale);
            costPatchCand = costPatchCand + costScale;
        end
        % Compute patch matching cost: planar compatibiltiy cost
        if(opt.usePlaneGuide)
            costPlane = sr_patch_cost_plane(modelPlane.mLogLPlaneProb, uvPlaneIDCand, uvPixValid.ind, srcPosCur(uvValidInd, :));
            costPatchCand = costPatchCand + opt.lambdaPlane*costPlane;
        end
        
        % Check which one to update
        updateInd = costPatchCand < uvCostCur;
        
        % =====================================================================
        % Update the nearest neighbor field
        % =====================================================================
        nUpdate = sum(updateInd);
        
        if(nUpdate ~= 0)
            uvUpdatePos = uvValidPos(updateInd);
            % === Update NNF data ===
            NNF.uvTformH.data(uvUpdatePos, :) = uvTformHCand(updateInd,:);
            NNF.uvTformA.data(uvUpdatePos, :) = uvTformACand(updateInd,:);
            NNF.uvPlaneID.data(uvUpdatePos)   = uvPlaneIDCand(updateInd);
            NNF.uvCost.data(uvUpdatePos)      = costPatchCand(updateInd);
            
            % Apply bias correction
            if(opt.useBiasCorrection)
                uvBiasCand = uvBiasCand(:,:,updateInd);
                NNF.uvBias.data(:,:,uvUpdatePos)   = uvBiasCand;
            end
            NNF.update.data(uvUpdatePos)   = 1; % updateInd;
            
            % === Update NNF map ===
            uvPixValidInd = uvPixValid.ind(updateInd);
            NNF.uvTformH.map    = sr_update_uvMap(NNF.uvTformH.map,  uvTformHCand(updateInd,:), uvPixValidInd);
            NNF.uvTformA.map    = sr_update_uvMap(NNF.uvTformA.map,  uvTformACand(updateInd,:), uvPixValidInd);
            NNF.uvPlaneID.map   = sr_update_uvMap(NNF.uvPlaneID.map, uvPlaneIDCand(updateInd),  uvPixValidInd);
            NNF.uvCost.map      = sr_update_uvMap(NNF.uvCost.map,    costPatchCand(updateInd),  uvPixValidInd);
            if(opt.useBiasCorrection)
                uvBiasCand = squeeze(uvBiasCand)';
                NNF.uvBias.map  = sr_update_uvMap(NNF.uvBias.map, uvBiasCand, uvPixValidInd);
            end
            NNF.update.map  = sr_update_uvMap(NNF.update.map, 1, uvPixValidInd);
            
            % === Update uvValidInd ===
            uvPixNextSub = uvPixValid.sub(updateInd,:);
            uvPixNextSub = bsxfun(@plus, uvPixNextSub, opt.propDir(indDirection,:));
            uvPixNextInd = sub2ind([NNF.imgH, NNF.imgW], uvPixNextSub(:,2), uvPixNextSub(:,1));
            
            updateMap = NNF.uvPix.mask;
            updateMap(uvPixNextInd) = 0;
            uvValidInd = ~updateMap(NNF.uvPix.ind);
            uvValidInd = uvValidInd & uvPixN.validInd;
            
            nUpdateTotal = nUpdateTotal + nUpdate;
            
            numUpdatePix = sum(uvValidInd);
        end
    end
end

end