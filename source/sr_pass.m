function NNF = sr_pass(imgTrg, imgSrcPyr, NNF, modelPlane, numIterLvl, opt)
%
% SR_PASS
%
% Nearest neighbor field estimation at the current level using the
% generalized PatchMatch algorithm
%
% Input:
%   - imgTrg:     target image
%   - imgSrcPyr:  source image pyramid
%   - NNF:        the current nearest neighbor field 
%   - modelPlane: planar structure model 
%   - numIterLvl: number of iterations at the current level
%   - opt:        parameters
% Output:
%   - NNF:        updated nearest neighbor field
% =========================================================================

% =========================================================================
% Compute the initial patch cost at the current level
% =========================================================================
% Prepare target patch
trgPatch = sr_prep_target_patch(imgTrg, opt.pSize);

% Prepare source patch
[srcPatch, srcPatchScale] = sr_prep_source_patch(imgSrcPyr, NNF.uvTformH.data, opt);
% Compute patch matching cost: appearance cost
[NNF.uvCost.data, NNF.uvBias.data] = sr_patch_cost_app(trgPatch, srcPatch, opt);
% Compute patch matching cost: scale cost       
if(opt.useScaleCost)
    costScale = opt.lambdaScale*max(0, opt.scaleThres - srcPatchScale);
    NNF.uvCost.data = NNF.uvCost.data + costScale;
end
% Compute patch matching cost: plane compatibility cost
if(opt.usePlaneGuide)
    costPlane = sr_patch_cost_plane(modelPlane.mLogLPlaneProb, NNF.uvPlaneID.data, NNF.uvPix.ind, NNF.uvTformH.data(:,7:8));
    NNF.uvCost.data = NNF.uvCost.data + opt.lambdaPlane*costPlane;
end

% Update cost map
NNF.uvCost.map = sr_update_uvMap(NNF.uvCost.map, NNF.uvCost.data, NNF.uvPix.ind);

% Initialize update index map (for early termination)
% NNF.update.data = true(NNF.uvPix.numUvPix, 1);
% NNF.update.map  = true(NNF.imgH, NNF.imgW);

% Initialize visualization
if(opt.visFlag)
    uvTformMap    = zeros(NNF.imgH, NNF.imgW, 3, numIterLvl);
    uvCostMap     = zeros(NNF.imgH, NNF.imgW, numIterLvl);
    uvCostMapInit = NNF.uvCost.map;
    imgRecSet     = zeros(NNF.imgH, NNF.imgW, 3, numIterLvl);
end

% =========================================================================
% Update the nearest neighbor field using PatchMatch 
% =========================================================================

for iter = 1 : numIterLvl
    [NNF, nUpdate] = sr_update_NNF(trgPatch, imgSrcPyr, NNF, modelPlane, opt);
    avgPatchCost = mean(NNF.uvCost.data);
    
    fprintf('    %3d\t%12d\t%12d\t%14f\n', iter, nUpdate(1), nUpdate(2), avgPatchCost);
      
    if(0)
        NNFVis = sr_vis_nnf(NNF, imgSrcPyr, opt);
        
%         uvTformMap(:,:,:, iter) = NNFVis.uvTfomMapVis;
%         uvCostMap(:,:, iter) = NNF.uvCost.map;
%         imgRecSet(:,:,:, iter) = imgRec;
        
%         figure(1);
%         subplot(2,4,1), imshow(imresize(imgSrcPyr{opt.origResLvl}, [NNF.imgH, NNF.imgW], 'bicubic'));
%         title('Current Image', 'fontsize', 16);
%         subplot(2,4,2), imshow(imgRec);
%         title('Super-resolution Image', 'fontsize', 16);
%         subplot(2,4,3), imshow(NNFVis.uvTformScaleVis);
%         title('Scale map', 'fontsize', 16);
%         subplot(2,4,4), imshow(NNFVis.uvTformRotVis); colormap jet
%         title('Rotation map', 'fontsize', 16);
%         subplot(2,4,5), imshow(NNFVis.uvTfomMapVis);
%         title('Nearest neighbor field', 'fontsize', 16);
%         subplot(2,4,6), imshow(NNFVis.uvPosMap); colormap jet
%         title('Position map', 'fontsize', 16);
%         subplot(2,4,7), imshow(NNFVis.uvCostMapVis); colormap jet
%         title('Cost map', 'fontsize', 16);
%         subplot(2,4,8), imshow(NNFVis.uvPlaneIDMapVis); colormap jet
%         title('Plane ID', 'fontsize', 16);
    end
end

% Visualization
if(opt.visFlag)
    resPath = 'result/nIter';
    if(~exist(resPath, 'dir'))
        mkdir(resPath);
    end
    imgID = opt.imgID;
    maxCost = max(uvCostMapInit(:));
    uvCostMap = uvCostMap/maxCost;
    
    uvCostMapInit = uvCostMapInit/maxCost;
    h2 = figure(2);
    imshow(uvCostMapInit); colormap jet;
    uvCostMapC = getframe(h2);
    imwrite(uvCostMapC.cdata,  ...
        fullfile(resPath, ['img_', num2str(imgID, '%03d'), '_costMap_iter_00.png']));
    
    % Position map
    uvPosMap = NNFVis.uvPosMap;
    
    imwrite(uvPosMap, fullfile(resPath, ...
        ['img_', num2str(imgID, '%03d'), '_uvPosMap.png']));
    
    for iter = 1: numIterLvl
        % Save uvTformMap
        imwrite(uvTformMap(:,:,:,iter), ...
            fullfile(resPath, ['img_', num2str(imgID, '%03d'),'_tformMap_iter_', num2str(iter, '%02d'), '.png']));
        
        % Save reconstructed image
        imwrite(imgRecSet(:,:,:,iter),  ...
            fullfile(resPath, ['img_', num2str(imgID, '%03d'),'_rec_iter_', num2str(iter, '%02d'), '.png']));
        
        % Save matching cost map
        h2 = figure(2);
        imshow(uvCostMap(:,:,iter)); colormap jet;
        uvCostMapC = getframe(h2);
        imwrite(uvCostMapC.cdata,  ...
            fullfile(resPath, ['img_', num2str(imgID, '%03d'), '_costMap_iter_', num2str(iter, '%02d'), '.png']));
    end
end

end