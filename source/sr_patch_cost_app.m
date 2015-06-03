function [uvCostApp, uvBias] = sr_patch_cost_app(trgPatch, srcPatch, opt)

% SR_PATCH_COST_APP
%
% Compute appearance-based patch cost: the weighted sum of the absolute difference 
% between between source and target patches

% Input:
%   - trgPatch: target patches - [patchSize*patchSize] x [3] x [numUvPix]
%   - srcPatch: source patches - [patchSize*patchSize] x [3] x [numUvPix]
%   - opt:      parameters
% Output:
%   - uvCostApp: weighted sum of squared difference - [numUvPix] x [1]
%   - uvBias:    color bias term - [1] x [3] x [numUvPix]
% =========================================================================

% =========================================================================
% Apply bias correction
% =========================================================================
if(opt.useBiasCorrection)
    % Mean of source and target patch
    meanTrgPatch = mean(trgPatch, 1);
    meanSrcPatch = mean(srcPatch, 1);
    
    % Compute bias and clamp it to inteval [opt.minBias, opt.maxBias]
    uvBias = meanTrgPatch - meanSrcPatch;
    uvBias = sr_clamp(uvBias, opt.minBias, opt.maxBias);
        
    % Apply the bias compensation
    srcPatch = bsxfun(@plus, srcPatch, uvBias);
else
    uvBias = [];
end
patchDist = trgPatch - srcPatch;

% =========================================================================
% Compute weighted sum of squared distance
% =========================================================================
% Sum of absolute distance
if(strcmp(opt.costType, 'L1'))
    patchDist = abs(patchDist);
elseif(strcmp(opt.costType, 'L2'))
    patchDist = patchDist.^2;
end

% Apply per pixel weight (a Gaussian falloff function)
uvCostApp = bsxfun(@times, patchDist, opt.wPatch);
uvCostApp = squeeze(sum(sum(uvCostApp, 1),2));

end