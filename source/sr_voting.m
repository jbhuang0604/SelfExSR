function imgRec = sr_voting(imgPyr, NNF, opt)

% SR_VOTING: 
%
% Image reconstruction using patch-based voting
%
% Input:
%   - imgPyr: high-freq band image pyramid
%   - NNF:    the estimated nearest neighbor field
%   - opt:    parameters
% Output:
%   - imgRec: reconstructed image 

% =========================================================================

imgH = NNF.imgH;
imgW = NNF.imgW;

% =========================================================================
% Prepare source patch
% =========================================================================
[srcPatch, ~] = sr_prep_source_patch(imgPyr, NNF.uvTformH.data, opt);

% Apply bias correction
if(opt.useBiasCorrection)
    srcPatch  = bsxfun(@plus, srcPatch, NNF.uvBias.data);
end

% Patch weighting using a Gaussian falloff function
srcPatch = bsxfun(@times, srcPatch, NNF.wPatch);

% =========================================================================
% Compute weighted average from source patches
% =========================================================================
imgAcc = zeros(imgH, imgW, 3, 'single');
imgH_  = imgH - opt.pSize;
imgW_  = imgW - opt.pSize;

for i = 1: opt.pSize*opt.pSize
    imgR = col2im(srcPatch(i,1,:), [opt.pSize,opt.pSize], [imgH, imgW], 'sliding');
    imgG = col2im(srcPatch(i,2,:), [opt.pSize,opt.pSize], [imgH, imgW], 'sliding');
    imgB = col2im(srcPatch(i,3,:), [opt.pSize,opt.pSize], [imgH, imgW], 'sliding');
    [r, c] = ind2sub([opt.pSize,opt.pSize], i);
    imgAcc(r:r+imgH_,c:c+imgW_, :) = imgAcc(r:r+imgH_,c:c+imgW_, :) + cat(3, imgR, imgG, imgB);
end
imgRec = bsxfun(@rdivide, imgAcc, NNF.wSumImg);

% Iterative back-projection method
imgRec = sr_backprojection(imgRec, imgPyr{opt.origResLvl}, opt.bpKernelSigma, opt.nIterBP);

end

