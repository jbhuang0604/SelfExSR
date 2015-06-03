function [imgHiRes] = sr_backprojection(imgHiRes, imgLowRes, sigma, nIter)

% SR_BACKPROJECTION
%
% Iterative back-projection algorithm 
% 
% =========================================================================
[imgH_L, imgW_L, ~] = size(imgLowRes);
[imgH_H, imgW_H, ~] = size(imgHiRes);

% Defining backprojection kernel
f = fspecial('gaussian', 5, sigma);
f = f.^2;
f = f./sum(f(:));

% Iterative reconstruction
for iIter = 1:nIter
    imgLowRes_S = imresize(imgHiRes, [imgH_L, imgW_L], 'bicubic');  
    imgDiff_H   = imresize(imgLowRes - imgLowRes_S, [imgH_H, imgW_H], 'bicubic');
    imgHiRes    = imgHiRes + imfilter(imgDiff_H, f, 'conv', 'same', 'replicate');

    % Clamp result to be [0, 1]
    imgHiRes = sr_clamp(imgHiRes, 0, 1);
end

