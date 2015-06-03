function modelPlanePyr = sr_planar_structure_pyramid(scaleImgPyr, modelPlane, topLevel)

%
% SR_PLANAR_STRUCTURE_PYRAMID:
%
% Pre-compute the models of planes and regularity according to the image
% pyramid, rescale parameters according to the image sizes.
%
% Output: 
%   - modelPlanePyr: plane model
%   - modelRegPyr: regularity model

numLevel = length(scaleImgPyr);

modelPlanePyr = cell(numLevel, 1);

for iLvl = topLevel: numLevel
    % === Update plane model ===
    scaleImgCur = scaleImgPyr{iLvl}.imgScale;
    modelPlanePyr{iLvl}.numPlane = modelPlane.numPlane;

    % Update rectification matrix   
    for iPlane = 1: modelPlane.numPlane
        H = eye(3);
        vLine = modelPlane.plane{iPlane}.vLine;
        vLine(1:2) = vLine(1:2)/scaleImgCur;
        H(3,:) = vLine;
        modelPlanePyr{iLvl}.rectMat{iPlane} = H;
        
        % resize the posterior probability map
        postProb = imresize(modelPlane.postProb, ...
            [scaleImgPyr{iLvl}.imgSize(1), scaleImgPyr{iLvl}.imgSize(2)], 'bicubic');
        postProb = bsxfun(@rdivide, postProb, sum(postProb, 3));
        modelPlanePyr{iLvl}.planeProb = postProb;
        modelPlanePyr{iLvl}.mLogLPlaneProb = -log(postProb);
    end
end

end