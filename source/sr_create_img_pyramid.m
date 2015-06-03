function [imgPyrH, imgPyrL, scaleImgPyr] = sr_create_img_pyramid(img, opt)

% SR_CREATE_IMG_PYRAMID
%
% Create two image pyramids: imgPyrH, imgPyrL.
% The low-passed image pyramid imgPyrL contains the low-passed version of imgPyrH. 
% This follows the use of high-freq. band and low-freq super-resolution for self-exemplars super-resolution 
% proposed in 
% 
%  G Freedman, R Fattal, Image and video upscaling from local self-examples
%  ACM Transactions on Graphics (TOG), 2011 
%
% Input:
%   - img:         input low-resolution image 
%   - opt:         pyramid option
% Output:
%   - imgPyrH:     high-freq band image pyramid
%   - imgPyrL:     low-freq band image pyramid
%   - scaleImgPyr: image dimensions in each level
% =========================================================================

% =========================================================================
% Image pyramid sizes and scales
% =========================================================================

% Image size in the high-resolution image
[imgHeight, imgWidth, ~] = size(img);

% Determine the scales in the image pyramid 
scalePyr = opt.alpha.^linspace(-opt.nPyrLowLvl, opt.nPyrLowLvl, opt.nPyrLvl);
% Add the dummy level for creating the low-pass version of the pyramid
scalePyr = [scalePyr, scalePyr(end)*opt.alpha]; 

% Image size in each layer
imgHPyr = round(imgHeight *scalePyr);
imgWPyr = round(imgWidth  *scalePyr);

% Initialize image pyramid
imgPyrH      = cell(opt.nPyrLvl, 1);
imgPyrL      = cell(opt.nPyrLvl, 1);
scaleImgPyr  = cell(opt.nPyrLvl, 1);

% =========================================================================
% Create high-freq band pyramid from the original low-resolution image
% =========================================================================

% Create the scales for the image pyramid
for k = opt.topLevel: opt.nPyrLvl + 1
    scaleImgPyr{k}.imgScale = scalePyr(k);
    scaleImgPyr{k}.imgSize  = [imgHPyr(k), imgWPyr(k)];
end

for k = opt.origResLvl: opt.nPyrLvl + 1
    imgPyrH{k} = imresize(img, [imgHPyr(k), imgWPyr(k)], opt.resampleKernel);
end

% Original resolution level
imgPyrH{opt.origResLvl} = img;
scaleImgPyr{opt.origResLvl}.imgScale = 1;
scaleImgPyr{opt.origResLvl}.imgSize = [imgHeight, imgWidth];

% =========================================================================
% Create low-freq band pyramid from the high-freq band pyramid
% =========================================================================

for k = opt.origResLvl: opt.nPyrLvl
    imgPyrL{k} = imresize(imgPyrH{k+1}, [imgHPyr(k), imgWPyr(k)], opt.resampleKernel);
end

end