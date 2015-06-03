function imgTformVis = sr_vis_planar_tform(img, modelPlane)

imgTformVis = [];

[imgH, imgW, nCh] = size(img);

% Target position
pRad = 15;
n = 7;
numUvPix = n^2;
trgPos = ones(2, numUvPix);
trgPos(1,:) = trgPos(1,:)*(imgW+1)/2;
trgPos(2,:) = trgPos(2,:)*(imgH+1)/2;

% Source position
sX = linspace(1, imgW, n+2);
sY = linspace(1, imgH, n+2);

[X Y] = meshgrid(sX(2:end-1), sY(2:end-1));
srcPos = cat(1, X(:)', Y(:)');

%
uvTformA = zeros(4, numUvPix);
uvTformA([1,4],:) = 1;
uvPlaneID = ones(1, numUvPix);
uvTformH = sr_src_domain_tform(uvPlaneID, modelPlane, uvTformA, srcPos, trgPos);

% Compute planar transform
% 1 2
% 4 3

% trgPatchPos{1} = bsxfun(@plus, trgPos, [-pRad, -pRad]');
% trgPatchPos{2} = bsxfun(@plus, trgPos, [ pRad, -pRad]');
% trgPatchPos{3} = bsxfun(@plus, trgPos, [ pRad, pRad]');
% trgPatchPos{4} = bsxfun(@plus, trgPos, [-pRad, pRad]');

patchPos = [-pRad, pRad, pRad, -pRad, -pRad;
    -pRad, -pRad, pRad, pRad, -pRad;
    1, 1, 1, 1, 1]; 
trgPatchPos = bsxfun(@plus, trgPos(:,1), patchPos(1:2,:));
 
h3 = figure(3); imshow(img); hold on;
plot(trgPatchPos(1,:), trgPatchPos(2,:), 'LineWidth', 5, 'Color', 'r');
 
for i = 1: numUvPix
    H = uvTformH(:,i);
    H = reshape(H, 3, 3);
    patchPosS = H*patchPos;
    patchPosS = bsxfun(@rdivide, patchPosS, patchPosS(3,:));

    plot(patchPosS(1,:), patchPosS(2,:), 'LineWidth', 3, 'Color', 'b');
end
imgTformVis = getframe(h3);
imgTformVis = imgTformVis.cdata;
hold off; 
 
end