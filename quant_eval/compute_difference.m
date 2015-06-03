function [psnr, ssim, ifc] = compute_difference(im, imgt, SRF)
% COMPUTE_DIFFERENCE: compute image quality
%
% Input:
%     - im:   super-resolved image
%     - imgt: groundtruth high-resolution image
%     - SRF:  super-resolution factor
% Output:
%     - psnr: Peak signal-to-noise ratio 
%     - ssim: Structural similarity index
%     - ifc:  Information fidelity criterion

% =========================================================================
% Retrieve only luminance channel
% =========================================================================
if size(im,3)>1
    im = rgb2ycbcr(im);
    im = im(:,:,1);
end

if size(imgt,3)>1
    imgt = rgb2ycbcr(imgt);
    imgt = imgt(:,:,1);
end

% =========================================================================
% Remove border pixels as some methods (e.g., A+) do not predict border pixels
% =========================================================================
cropPix     = SRF;
im          = shave(im, [cropPix, cropPix]);
imgt        = shave(imgt, [cropPix, cropPix]); 

% Convert to double (with dynamic range 255)
im          = double(im); 
imgt        = double(imgt); 

% =========================================================================
% Compute Peak signal-to-noise ratio (PSNR)
% =========================================================================
mse = mean(mean((im - imgt).^2, 1), 2);
psnr = 10*log10(255*255/mse);

% =========================================================================
% Compute Structural similarity index (SSIM index)
% =========================================================================
[ssim, ~] = ssim_index(imgt, warpedImage);

% =========================================================================
% Compute information fidelity criterion (IFC)
% =========================================================================
ifc = ifcvec(imgt, warpedImage);

end

function I = shave(I, border)
I = I(1+border(1):end-border(1), ...
      1+border(2):end-border(2), :, :);
end