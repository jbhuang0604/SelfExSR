function modelPlane = sr_extract_plane(imgPath, imgName, opt)

% SR_EXTRACT_PLANE:
%
% Extract plane model from an image
%
% The code is adapted and modified from the paper
%
% Jia-Bin Huang, Sing Bing Kang, Narendra Ahuja, Johannes Kopf
% Image Completion using Planar Structure Guidance
% ACM Transactions on Graphics (Proceedings of SIGGRAPH 2014).
%
% Output: modelPlane
% The model plane has the following data structure
%
% modelPlane
% modelPlane.vpData:                        Vanishing point data
% modelPlane.numPlane                       Number of planes
% modelPlane.plane{indPlane}.vLine          The vanishing line of the plane
% modelPlane.plane{indPlane}.imgPlaneProb;  The planar location density
% modelPlane.plane{indPlane}.sourceVP       Which two VPs form the plane
% modelPlane.plane{indPlane}.rotPar(vpInd)  Rotation parameters for aligning
%                                           two sets of lines with the x-axis                                           x-axis
% modelPlane.plane{indPlane}.postProb       Posteior probability of the plane

% =========================================================================
% Vanishing point detection
% =========================================================================

vpFilePath = 'cache\vpdetection';
vpFileName  = [imgName(1:end-4), '-vanishingpoints.txt'];

recomputeFlag = 1;
if(~exist(fullfile(vpFilePath, 'text', vpFileName), 'file') || recomputeFlag)
    vpExeFile = 'source\EdgeDetectTest.exe';
    vpDetectCMD = [vpExeFile, ' -indir ', imgPath, ' -infile ', imgName, ' -outdir ', vpFilePath];
    system(vpDetectCMD);
end
% Read vanishing point data
vpData = sr_read_vpdata(fullfile(vpFilePath, 'text', vpFileName));

img = imread(fullfile(imgPath, imgName));

% =========================================================================
% Plane localization
% =========================================================================

modelPlane = sr_detect_plane_from_vp(vpData, img, opt);
% figure(1), imshow(modelPlane.postProb(:,:,1:3));

visVPFlag = 0;
if(visVPFlag)
    [vpVis, planeVis] = vis_vp(img, vpData, modelPlane.postProb(:,:,1:3));
end
end

function [vpVis, planeVis] = vis_vp(img, vpData, postProb)

img = im2double(img);

[imgH, imgW, nCh] = size(img);
vpVis = zeros(imgH, imgW, nCh);
planeVis = zeros(imgH, imgW, nCh);
whiteImg = ones(imgH, imgW, nCh);
% blackImg = zeros()
alphaW = 0.5;

imgW = img*(1 - alphaW) + whiteImg*(alphaW);

% Plot VP
h1 = figure(1);
imshow(imgW); hold on;
lineWidth = 4;
for vpInd = 1: vpData.numVP
    vpCurr = vpData.vp{vpInd};
    for lineInd = 1: vpCurr.numLines
        lineCur = vpCurr.lines(lineInd,:);
        lineCur = lineCur + 1;
        if(vpInd == 1)
            plot(lineCur([1,3]), lineCur([2,4]), 'r', 'LineWidth', lineWidth);
        elseif(vpInd == 2)
            plot(lineCur([1,3]), lineCur([2,4]), 'g', 'LineWidth', lineWidth);
        elseif(vpInd == 3)
            plot(lineCur([1,3]), lineCur([2,4]), 'b', 'LineWidth', lineWidth);
        end
    end
end
hold off;
print(h1, '-dpng', fullfile('paper\planar_struct_SR_CVPR2015\figures\plane', 'vpdetection.png'));

% Plot posterior
alphaP = 0.9;
imgP = (1 - alphaP)*img + alphaP*postProb;
figure(2); imshow(imgP);
imwrite(imgP, fullfile('paper\planar_struct_SR_CVPR2015\figures\plane', 'planedetection.png'));

end
function vpData = sr_read_vpdata(fileName)

% SC_READ_VPDATA: read the data from vanishing point detection algorithm
% Input:
%   - fileName: the txt file containing the pre-computed vanishing point
%   detection code
% Output:
%   - vpData
%   The data structure of vpData
%       - vpData.numVP: number of detected vanishing points
%       - vpData.vp{i}.pos: the vanishing point position in the homogenous coordiante
%       - vpData.vp{i}.score: the score of the vanishing point
%       - vpData.vp{i}.numLines: number of lines supporting the vanishing point
%       - vpData.vp{i}.lines{j}.p1: (x1, y1): starting position
%       - vpData.vp{i}.lines{j}.p2: (x2, y2): ending position
%       - vpData.vp{i}.lines{j}.length: length of the line segment

vpData = [];

% Read data
fid = fopen(fileName);

%% Parse VP positions
temp = fscanf(fid, '%s ', [1 5]);
numVP = 0;
readVPFlag = 1;
VP = [];
while(readVPFlag)
    numVP = numVP + 1;
    vpCurr = fscanf(fid, '%g %g %g %g %g', [5 1]);
    if(~isempty(vpCurr))
        VP(:,numVP) = vpCurr;
    else
        temp = fscanf(fid, '%s ', [1 6]);
        readVPFlag = 0;
    end
end
VP = VP';

vpData.numVP = size(VP, 1);

% Save VP position data
for i = 1: vpData.numVP
    vpData.vp{i}.pos = VP(i, 1:3);
    vpData.vp{i}.score = VP(i, 4);
    vpData.vp{i}.numLines = VP(i, 5);
end

%% Parse each set of line segments for the corresponding VP
for i = 1: vpData.numVP
    numLine = fscanf(fid, '%d ', [1 1]);
    lines = fscanf(fid, '%g %g %g %g %g', [5 numLine]);
    vpData.vp{i}.lines = lines';
end

fclose(fid);
end

function modelPlane = sr_detect_plane_from_vp(vpData, img, opt)

% SC_DETECT_PLANE_FROM_VP: simple plane detection algorithm
% Input:
%     - vpData: vanishing point data
%     - img: input image
%     - mask: hole mask
% Output:
%     - modelPlane

%%

modelPlane = [];

% === Setting up ===
[imgH, imgW, ch] = size(img);
HfilterX = fspecial('gaussian', [1, opt.filterSize], opt.filterSigma);
HfilterY = HfilterX';
% fspecial('gaussian', opt.filterSize, opt.filterSigma);

img = im2double(img);

% === Supporting lines spatial support estimation ===
shapeInserter = vision.ShapeInserter('Shape', 'Lines','BorderColor', 'White');

for i = 1: vpData.numVP
    % The support lines
    imgLines = zeros(imgH, imgW);
    imgLines = step(shapeInserter, imgLines, int16(round(vpData.vp{i}.lines(:,1:4))));
    % Spatial density estimation via blurring
    imgLinesPosMap = imgLines;
    for k = 1:opt.numFilterIter
        imgLinesPosMap = imfilter(imgLinesPosMap, HfilterX, 'conv', 'replicate');
    end
    for k = 1:opt.numFilterIter
        imgLinesPosMap = imfilter(imgLinesPosMap, HfilterY, 'conv', 'replicate');
    end
    
    % Save results
    modelPlane.vp{i}.imgLines = imgLines;
    modelPlane.vp{i}.imgLinesPosMap = imgLinesPosMap;
end


% === Estimate plane support and plane parameters ===
numPlane = (vpData.numVP)*(vpData.numVP-1)/2;
% Initialize plane data
modelPlane.plane = cell(numPlane, 1);

indPlane = 1;
% A pair of vanishing points forms a plane hypothesis
for i = 1: vpData.numVP - 1
    for j = i+1: vpData.numVP
        % Compute the vanishing line
        modelPlane.plane{indPlane}.vLine = vLineFromTwoVP(vpData.vp{i}.pos, vpData.vp{j}.pos);
        % Element-wise product of two support line density
        modelPlane.plane{indPlane}.imgPlaneProb = modelPlane.vp{i}.imgLinesPosMap.*modelPlane.vp{j}.imgLinesPosMap; % Product of two probability maps
        
        %         modelPlane.plane{indPlane}.imgPlaneProb(mask) = 1e-10;
        modelPlane.plane{indPlane}.sourceVP = [i, j];
        
        indPlane = indPlane + 1;
    end
end


% === Compute rectified rotation parameters ===

for i = 1: numPlane
    for vpInd = 1: 2
        
        linesCurr = vpData.vp{modelPlane.plane{i}.sourceVP(vpInd)}.lines;
        invalidLineInd = linesCurr(:,5) == 0;
        linesCurr = linesCurr(~invalidLineInd,:);
        numLines = size(linesCurr, 1);
        
        vLineCurr = modelPlane.plane{i}.vLine;
        
        % Rectified homography
        H = eye(3);
        H(3,:) = vLineCurr;
        
    end
end


% === Add a fronto-parallel plane ===

modelPlane.plane{indPlane}.vLine = [0 0 1];
modelPlane.plane{indPlane}.imgPlaneProb = opt.fpPlaneProb*ones(imgH, imgW);
modelPlane.plane{indPlane}.score = sum(modelPlane.plane{indPlane}.imgPlaneProb(:));

numPlane = numPlane + 1;

modelPlane.numPlane = numPlane;

% === Compute posterior probability ===

planeProb = zeros(imgH, imgW, numPlane);
for i = 1 : numPlane
    planeProb(:,:,i) = modelPlane.plane{i}.imgPlaneProb;
end
planeProbSum = sum(planeProb, 3);
planeProb = bsxfun(@rdivide, planeProb, planeProbSum);
planeProb = planeProb + 0.1; % blur the posterior map
planeProb = bsxfun(@rdivide, planeProb, sum(planeProb, 3));

modelPlane.postProb = planeProb;

end

function vLine = vLineFromTwoVP(vp1, vp2)

A = cat(1, vp1, vp2);

[U S V] = svd(A, 0);
vLine = V(:,end);
vLine = vLine/vLine(3); % [h7, h8, 1]

end