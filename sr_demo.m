function imgHiRes = sr_demo(filePath, opt)
% SR_DEMO
%
% Input:
%     - filePath: filePath contains
%       (1) dataPath: path to the input images
%       (2) imgFileName: the file name of the input low-resolution image
%     - opt: algorithm parameters, obtained by
%            opt = sr_init_opt(SRF);
% Output:
%     - imgHiRes: the super-resolved images
%
% Example usage:
%   filePath.dataPath    = 'data\Urban100\image_SRF_4';
%   filePath.imgFileName = 'img_001_SRF_4_LR.png';
%
%   imgHiRes = sr_demo(filePath, opt);
% Disclaimer:
%
%   It is provided for educational/researrch purpose only.
%   If you find the software useful, please consider cite our paper.
%
%   Single Image Super-Resolution using Transformed Self-Exemplars
%   Jia-Bin Huang, Abhishek Singh, and Narendra Ahuja
%   IEEE Conference on Computer Vision and Pattern Recognition, CVPR 2015
%
% Contact:
%   Jia-Bin Huang
%   Electrical and Computer Engineering
%   University of Illinois, Urbana-Champaign
%   www.jiabinhuang.com

% =========================================================================
% Extract planar structures
% =========================================================================
fprintf('- Extract planar structures \n');
tic;
modelPlane = sr_extract_plane(filePath.dataPath, filePath.imgFileName, opt);
tAnalysis = toc;
fprintf('Done in %6.3f seconds.\n\n', tAnalysis);

% =========================================================================
% Construct image pyramid for coarse-to-fine image super-resolution
% =========================================================================
tic;
% Load the input low-resolution image:
img = im2single(imread(fullfile(filePath.dataPath, filePath.imgFileName)));
if(ndims(img)~=3) % support only RGB images
    img = cat(3, img, img, img);
end
tImgPyramid = toc;
fprintf('Done in %6.3f seconds.\n\n', tImgPyramid);

% Create image pyramid and multi-level planar structure constraints
fprintf('- Construct plane pyramid: \n');
tic;
[imgPyrH, imgPyrL, scaleImgPyr] = sr_create_img_pyramid(img, opt);
modelPlane = sr_planar_structure_pyramid(scaleImgPyr, modelPlane, opt.topLevel);
tPlanePyramid = toc;
fprintf('Done in %6.3f seconds.\n\n', tPlanePyramid);

% =========================================================================
% Super-resolution by patch-based synthesis
% =========================================================================
fprintf('- Single Image Super-Resolution using Transformed Self-Exemplars \n');
tic;
imgPyr = sr_synthesis(imgPyrH, imgPyrL, scaleImgPyr, modelPlane, filePath, opt);
tSynthesis = toc;
fprintf('Synthesis took %6.3f seconds.\n', tSynthesis);

% Get the level corresponding to the desired high-resolution image
if(mod(opt.SRF,2) == 0 )
    lvlInd = opt.origResLvl - log(opt.SRF)/log(2)*opt.nLvlToRedRes;
else
    lvlInd = opt.origResLvl - log(opt.SRF)/log(3)*opt.nLvlToRedRes;
end

imgHiRes = imgPyr{uint8(lvlInd)};

% Visualize plane parameters
% if(0)
%     imgTformVis = sr_vis_planar_tform(imgPyr{opt.origResLvl-6}, modelPlane{opt.origResLvl-6});
%
%     imgVisResName = [filePath.imgFileName(1:end-4), '_vis.png'];
%     imwrite(imgTformVis, fullfile('result\vis', imgVisResName));
% end

end