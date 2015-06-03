function NNFVis = sr_vis_nnf(NNF, imgPyrH, opt)

mask = NNF.uvPix.mask;

[imgH, imgW] = size(mask);

bdMask = edge(mask);
bdMask = bdMask(:,:,ones(3,1));
mask = mask(:,:,ones(3,1));

%% uvTfomMapVis
[NNFVis.uvTfomMapVis, NNFVis.uvPosMap]  = vis_tform_map(NNF.uvTformH.map);

%% scale map
uvTformScale = sr_scale_tform(NNF.uvTformH.data);
NNFVis.uvTformScaleVis = zeros(NNF.imgH, NNF.imgW);
NNFVis.uvTformScaleVis = sr_update_uvMap(NNFVis.uvTformScaleVis, uvTformScale, NNF.uvPix.ind);

NNFVis.uvTformScaleVis = NNFVis.uvTformScaleVis/opt.maxScale;

% NNFVis.uvTformScaleVis = (NNFVis.uvTformScaleVis - min(NNFVis.uvTformScaleVis(:)))/...
%     (opt.maxScale - opt.minScale);
% NNFVis.uvTformScaleVis = (NNFVis.uvTformScaleVis - min(NNFVis.uvTformScaleVis(:)))/max(NNFVis.uvTformScaleVis(:));

%% rotation map

uvTformRot = NNF.uvTformA.data;
uvTformRot = atan2(uvTformRot(2,:), uvTformRot(1,:));
NNFVis.uvTformRotVis = zeros(NNF.imgH, NNF.imgW);
NNFVis.uvTformRotVis = sr_update_uvMap(NNFVis.uvTformRotVis, uvTformRot, NNF.uvPix.ind);
NNFVis.uvTformRotVis = (NNFVis.uvTformRotVis + pi/4)/(pi/2);

%% uvCostMapVis

NNFVis.uvCostMapVis = (NNF.uvCost.map - min(NNF.uvCost.map(:)))/max(NNF.uvCost.map(:));

%% uvBiasMapVis

NNFVis.uvBiasMapVis = NNF.uvBias.map + 0.5;
% NNFVis.uvBiasMapVis(bdMask) = 1;

%% uvPixUpdateSrc

if(0)
    NNFVis.uvPixUpdateSrcMap = zeros(imgH, imgW, 3);
    for ch = 1:3
        NNFVis.uvPixUpdateSrcMap(:,:,ch) = im2double(NNF.uvPixUpdateSrc.map == ch);
    end
end
% NNFVis.uvGainMapVis = NNF.uvGain.map;
% NNFVis.uvGainMapVis = NNFVis.uvGainMapVis - 0.5;
% NNFVis.uvGainMapVis(bdMask) = 1;

%% uvPlaneIDMapVis
numPlane = NNF.uvPlaneID.numPlane;
NNFVis.uvPlaneIDMapVis = NNF.uvPlaneID.map/numPlane;

if(0)
    numPlane = NNF.uvPlaneID.numPlane;
    NNFVis.uvPlaneIDMapVis = NNF.uvPlaneID.map/numPlane;
    %     NNFVis.uvPlaneIDMapVis = zeros(imgH, imgW, 3);
    %     for i = numPlane:-1:1
    %         if(i == numPlane)
    %             fpIDMap = NNF.uvPlaneID.map == i;
    %             for ch = 1: 3
    %                 NNFVis.uvPlaneIDMapVis(:,:,ch) = NNFVis.uvPlaneIDMapVis(:,:,ch) +  im2double(fpIDMap);
    %             end
    %             %         NNFVis.uvPlaneIDMapVis(:,:,1) = fpIDMap;
    %             %         NNFVis.uvPlaneIDMapVis(:,:,2) = fpIDMap;
    %             %         NNFVis.uvPlaneIDMapVis(:,:,3) = fpIDMap;
    %         else
    %             NNFVis.uvPlaneIDMapVis(:,:,i) = NNFVis.uvPlaneIDMapVis(:,:,i) + ...
    %                 im2double(NNF.uvPlaneID.map == i);
    %         end
    %     end
    
    % NNFVis.uvPlaneIDMapVis()
    
    NNFVis.uvPlaneIDMapVis = im2double(NNFVis.uvPlaneIDMapVis);
end
% NNFVis.uvPlaneIDMapVis = NNF.uvPlaneID.map/numPlane;


figure(1);
subplot(2,4,1), imshow(imresize(imgPyrH{opt.origResLvl}, [NNF.imgH, NNF.imgW], 'bicubic'));
title('Current Image', 'fontsize', 16);
subplot(2,4,2), imshow(imgPyrH{opt.iLvl});
title('Super-resolution Image', 'fontsize', 16);
subplot(2,4,3), imshow(NNFVis.uvTformScaleVis);
title('Scale map', 'fontsize', 16);
subplot(2,4,4), imshow(NNFVis.uvTformRotVis); colormap jet
title('Rotation map', 'fontsize', 16);
subplot(2,4,5), imshow(NNFVis.uvTfomMapVis);
title('Nearest neighbor field', 'fontsize', 16);
subplot(2,4,6), imshow(NNFVis.uvPosMap); colormap jet
title('Position map', 'fontsize', 16);
subplot(2,4,7), imshow(NNFVis.uvCostMapVis); colormap jet
title('Cost map', 'fontsize', 16);
subplot(2,4,8), imshow(NNFVis.uvPlaneIDMapVis); colormap jet
title('Plane ID', 'fontsize', 16);




end

function [NNFMapCurVis, NNFMapVis] = vis_tform_map(NNFMap)

[imgH, imgW, ch] = size(NNFMap);

% Initialize NNFMapVis
NNFMapVis = zeros(imgH, imgW, 3, 'single');
[X, Y] = meshgrid(1:imgW, 1:imgH);
X = X/imgW;     Y = Y/imgH;

NNFMapVis(:,:,2) = 0.5;
NNFMapVis(:,:,1) = X;
NNFMapVis(:,:,3) = Y;

% Prepare the visualization of the NNF
NNFMapCurVis = zeros(imgH, imgW, 3, 'single');
NNFMapCurVis(:,:,2) = 0.5;
NNFMapCurVis(:,:,1) = NNFMap(:,:,7)/imgW;
NNFMapCurVis(:,:,3) = NNFMap(:,:,8)/imgH;

% NNFMapVis = NNFMapVis.*(1-mask) + NNFMapCurVis.*mask;
% NNFMapVis(bdMask) = 1;
% NNFMapVis = cat(2, NNFMapCurVis, NNFMapVis);

end