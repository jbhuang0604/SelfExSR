function uvPixValidInd = sr_get_uvpix(img, opt)

% SR_GET_UVPIX
%
%

if(ndims(img) == 3)
    img = rgb2gray(img);
end

% Gradient
[imgGx, imgGy] = gradient(img);
imgGx2 = imgGx.^2;
imgGy2 = imgGy.^2;
imgGxy = imgGx.*imgGy;

% Blur
h = fspecial('gaussian',[5 1], 1.5);
imgGx2 = conv2(h, h', imgGx2, 'same');
imgGy2 = conv2(h, h', imgGy2, 'same');

imgGradEng = imgGx2 + imgGy2;
% figure(1); imagesc(imgGradEng > opt.gradThres);
imgGradEng = imgGradEng(opt.pRad+1:end-opt.pRad ,opt.pRad+1:end-opt.pRad);

uvPixValidInd = imgGradEng(:) > opt.gradThres;

end