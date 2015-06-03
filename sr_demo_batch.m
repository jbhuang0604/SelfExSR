% SR_DEMO_BATCH
%
% Example script for super-resolving all images in folder in a batch mode
%
% Jia-Bin Huang
% Electrical and Computer Engineering
% University of Illinois, Urbana-Champaign
% www.jiabinhuang.com

% =========================================================================
% Environment settings
% =========================================================================

startup;

% Set up dataset
datasetList = {'Urban100', 'BSD100', 'Sun-Hays80', 'Set5', 'Set14'};  % Data set list: Urban100 and BSD100
numDataset = length(datasetList);        % Number of datatset
numImgDataset = [100,100, 80, 5, 14];    % Number of image in each dataset
SRFList{1} = [2,4];                      % Super-Resolution Factor (SRF) for Urban100 dataset
SRFList{2} = [2,3,4];                    % Super-Resolution Factor (SRF) for BSD100 dataset
SRFList{3} = [8];
SRFList{4} = [2,3,4];
SRFList{5} = [2,3,4];

% Set up file patch
% filePath.dataset     = [];             % Dataset name, e.g., 'Urban100'
% filePath.dataPath    = [];             % Path to images
% filePath.imgFileName = [];             % Input low-resolution image filename

% =========================================================================
% Start super-resolving images
% =========================================================================

for iDataset = 1 : numDataset
    % Current dataset
    datasetCur = datasetList{iDataset};
    
    for SRF = SRFList{iDataset}
        % Initialize the paramters for super-resolution
        opt = sr_init_opt(SRF);
        if(strcmp(datasetCur,'Sun-Hays80'))
            opt.scaleThres = 2;
        end
        % Process all images in the dataset
        parfor imgID = 1:numImgDataset(iDataset)
            filePath = [];
            filePath.dataPath   = fullfile('data', datasetCur, ['image_SRF_',num2str(SRF)]);
            filePath.resLvlPath = fullfile('data', datasetCur, ['image_SRF_',num2str(SRF),'_lvl']);
            
            % Image file name for the low-resolution image
            imgInName  = ['img_', num2str(imgID, '%03d'), '_SRF_',  num2str(SRF), '_LR.png'];
            
            % Image file name for the super-resolved image
            imgResName = ['img_', num2str(imgID, '%03d'), '_SRF_', num2str(SRF), '_SelfExSR.png'];
            
            disp(['Super-resolution ', datasetCur, ': ', imgResName]);
            
            if(~exist(fullfile(filePath.dataPath, imgResName), 'file') && ...
                    exist(fullfile(filePath.dataPath, imgInName), 'file'))
                filePath.imgFileName = imgInName;
                imgHiRes = sr_demo(filePath, opt);
                
                % Save results
                imwrite(imgHiRes, fullfile(filePath.dataPath, imgResName));
            end
        end
        
    end
    
end
