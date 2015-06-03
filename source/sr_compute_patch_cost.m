function uvPatchCost = sc_compute_patch_cost(img, NNF, distMap, uvPixels, validPixels, modelPlane, modelReg, optS)


% Compute the patch matching cost for all uvPixels
tic;
bbR = im2col(img(:,:,1),[optS.pSize, optS.pSize], 'sliding');
bbG = im2col(img(:,:,2),[optS.pSize, optS.pSize], 'sliding');
bbB = im2col(img(:,:,3),[optS.pSize, optS.pSize], 'sliding');
indUvPixels = sub2ind(size(distMap), uvPixels(:,2), uvPixels(:,1));
bbUvR = bbR(:,indUvPixels);
bbUvG = bbG(:,indUvPixels);
bbUvB = bbB(:,indUvPixels);
t1 = toc;
 
numUvPixel = size(uvPixels,1);
% bbUv = zeros(optS.pSize*optS.pSize*3, )
bbUvR = zeros(optS.pSize,optS.pSize, numUvPixel);
bbUvG = zeros(optS.pSize,optS.pSize, numUvPixel);
bbUvB = zeros(optS.pSize,optS.pSize, numUvPixel);

tic;
for i = 1: numUvPixel
    bbUvR(:,:,i) = img(uvPixels(i,2)-optS.pRad: uvPixels(i,2)+optS.pRad,  ...
        uvPixels(i,1)-optS.pRad: uvPixels(i,1)+optS.pRad, 1);
    bbUvG(:,:,i) = img(uvPixels(i,2)-optS.pRad: uvPixels(i,2)+optS.pRad,  ...
        uvPixels(i,1)-optS.pRad: uvPixels(i,1)+optS.pRad, 2);
    bbUvB(:,:,i) = img(uvPixels(i,2)-optS.pRad: uvPixels(i,2)+optS.pRad,  ...
        uvPixels(i,1)-optS.pRad: uvPixels(i,1)+optS.pRad, 3);
end
t2 = toc;
 
uvPatchCost = [];

end