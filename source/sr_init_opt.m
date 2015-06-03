function opt = sr_init_opt(SRF)

% SR_INIT_OPT: parameter initialization
%
% Jia-Bin Huang
% Electrical and Computer Engineering
% University of Illinois, Urbana-Champaign
% www.jiabinhuang.com
% =========================================================================


fprintf('- Initialize parameters \n');

opt = [];

% =========================================================================
% Multi-resolution processing parameters
% =========================================================================

opt.SRF = SRF;

if(mod(SRF,2) == 0)
    opt.nLvlToRedRes = 3;                    % Number of level to half of the original resolution
    opt.alpha = (1/2).^(1/opt.nLvlToRedRes); % Scale factor
    opt.coarseLvlImgScale = 1/8;             % The coarsest image scale
else
    opt.nLvlToRedRes = 5;                    % Number of level to one-third of the original resolution
    opt.alpha = (1/3).^(1/opt.nLvlToRedRes); % Scale factor
    opt.coarseLvlImgScale = 1/9;             % The coarsest image scale
end

% Number of pyramid levels
opt.nPyrLowLvl = round(log(opt.coarseLvlImgScale)/log(opt.alpha));
opt.nPyrLvl = 2*opt.nPyrLowLvl + 1;          % Total number of pyramid level
opt.origResLvl = round(opt.nPyrLvl/2);       % The level of the original image resolution
opt.resampleKernel = 'bicubic';              % Resampling kernel used for constructing the pyramid
                                             % options: 'lanczos3', 'bicubic', 'bilinear'

if(mod(SRF, 2) == 0)
    opt.topLevel = opt.origResLvl - ...      % Which level to stop
        (log(opt.SRF)/log(2))*opt.nLvlToRedRes;
else
    opt.topLevel = opt.origResLvl - ...      % Which level to stop
        (log(opt.SRF)/log(3))*opt.nLvlToRedRes;
end
opt.topLevel = uint8(opt.topLevel);

% =========================================================================
% Method parameters
% =========================================================================
% For larger SRF (e.g., SRF > 4), it is recommended to use larger opt.scaleThres to avoid undesired visual artifacts
opt.scaleThres  = 1/opt.alpha;               % Scale cost
opt.lambdaScale = 1*1e-3;                    % Weight for scale cost
opt.lambdaPlane = 1*1e-3;                    % Plane compatibility cost

% Number of iterations per level (more iterations typically improve the reconstruction quality, but slower)
opt.numIter    = 15;                         % The initial iteration number 15
opt.numIterDec = 3;                          % Number of decrements
opt.numIterMin = 3;                          % Minimum number of iterations

opt.nIterBP = 20;                            % Number of iterations for backprojection
opt.bpKernelSigma = 1;

% Method configuration
opt.useScaleCost      = 1;                   % Use scale cost
opt.usePlaneGuide     = 1;                   % Use planar guidance
opt.useAffine         = 1;                   % Use affine transformation
opt.useBiasCorrection = 1;                   % Use bias correction

opt.visFlag           = 0;                   

% =========================================================================
% Patch matching parameters
% =========================================================================

% === Patch size ===
opt.pSize = 5;
opt.pRad  = floor(opt.pSize/2);              % Patch radius
opt.pNumPix = opt.pSize*opt.pSize;           % Number of pixels in a patch
opt.pMidPix = round(opt.pNumPix/2);          % The center of the patch

% === Affine transform parameters ===
opt.scaleRadA = 1;                           % Search radius for patch scale
opt.rotRadA   = pi/4;                        % Search radius for patch rotation
opt.shRadA    = 0.05;                        % Search radius for patch sheer

% === Parameters for domain transformation ===
% This scale range is used to reject unlikely patch transformation
opt.minScale = 1/opt.alpha;                  % Mininum patch scale variation
opt.maxScale = 8;                            % Maximum patch scale variation

% === Parameters for photometric compensation ===
opt.minBias = -0.25;                         % Mininum bias compensation
opt.maxBias =  0.25;                         % Maximum bias compensation

opt.costType = 'L2';                         % Patch matching cost
wPatch = fspecial('gaussian', [opt.pSize,opt.pSize],3);
opt.wPatch = reshape(wPatch, opt.pNumPix, 1, 1);

% === Precomputed patch position in the reference position ===
[X, Y] = meshgrid(-opt.pRad:opt.pRad, -opt.pRad:opt.pRad);
opt.refPatchPos = single(cat(2, X(:), Y(:), ones(opt.pSize*opt.pSize, 1)));

% === Propagation directions ===
opt.propDir = [1 0; 0 1; -1 0; 0 -1];
opt.errThres = 0*1e-3;                        % Random sampling threshold for early termination
% opt.numRandSample = 5;                       % Number of coarse-to-fine random sampling per iteration

opt.filterSize = 100;
opt.filterSigma = 50;
opt.numFilterIter = 5;

opt.fpPlaneProb = 1e-4;

end