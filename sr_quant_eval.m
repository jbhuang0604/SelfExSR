% SR_QUANT_EVAL
%
% Script for computing averaged PSNR, SSIM, and IVC. 
% The script reproduces the quantitative results on our project page on
% Github: https://github.com/jbhuang0604/SelfExSR
%
% Jia-Bin Huang
% Electrical and Computer Engineering
% University of Illinois, Urbana-Champaign
% www.jiabinhuang.com

clc;
clear;
close all;

% Add pathes for running SSIM and IVC
addpath(genpath('quant_eval'));

% Dataset list
numDataset = 4;
dataset = cell(numDataset,1);
dataset{1}.name = 'Urban100';
dataset{2}.name = 'BSD100';
dataset{3}.name = 'Set5';
dataset{4}.name = 'Set14';

% Super-resolution factor
dataset{1}.SRF  = [2, 4];
dataset{2}.SRF  = [2, 3, 4];
dataset{3}.SRF  = [2, 3, 4];
dataset{4}.SRF  = [2, 3, 4];

% Number of images in the dataset
dataset{1}.numImg = 100;
dataset{2}.numImg = 100;
dataset{3}.numImg = 5;
dataset{4}.numImg = 14;

% Method list
methodList = {'Bicubic', 'ScSR', 'Kim', 'Abhishek', 'Glasner', 'SRCNN', 'A+', 'SelfExSR'};
numMethod  = length(methodList);
indValidMethod = [1, 1, 1, 1, 1, 1, 1, 1]; % Test only

% Initialize result path

resPath = fullfile('quant_eval', 'result');
if(~exist(resPath, 'dir'))
    mkdir(resPath);
end
% =========================================================================
% Start quantitative evaluation
% =========================================================================
for indDataset = 1 : numDataset
    datasetName   = dataset{indDataset}.name;
    numImgDataset = dataset{indDataset}.numImg;
    numSRF        = size(dataset{indDataset}.SRF,2);
    
    % Run each super-resolution factor
    for indSRF = 1: numSRF
        SRF = dataset{indDataset}.SRF(indSRF);
        imgPath = fullfile('data', datasetName, ['image_SRF_',num2str(SRF)]);
        
        resName = ['result_', datasetName, '_SRF_', num2str(SRF),'.mat'];
        if(~exist(fullfile(resPath, resName), 'file'))
            % Initialize result table
            SSIM_table = zeros(numImgDataset, numMethod);
            PSNR_table = zeros(numImgDataset, numMethod);
            IFC_table  = zeros(numImgDataset, numMethod);
            reverseStr = '';
            for indImg = 1: numImgDataset
                % Load groundtruth high-resolution image
                imgName = ['img_',num2str(indImg, '%03d'), '_SRF_', num2str(SRF), '_HR.png'];
                imgGT = imread(fullfile(imgPath, imgName));
                
                % Load super-resolution results
                for indMethod= 1:numMethod
                    if(indValidMethod(indMethod))
                        imgName = ['img_',num2str(indImg, '%03d'), '_SRF_', num2str(SRF), '_', methodList{indMethod}, '.png'];
                        img = imread(fullfile(imgPath, imgName), 'png');
                        
                        % Compute image quality
                        [psnr, ssim, ifc] = compute_difference(img, imgGT, SRF);
                        
                        PSNR_table(indImg,indMethod) = psnr;
                        SSIM_table(indImg,indMethod) = ssim;
                        IFC_table(indImg, indMethod)  = ifc;
                    end
                    % Display progress
                    percentDone = 100*((indImg-1)*numMethod + indMethod)/(numMethod*numImgDataset);
                    msg = sprintf('Evaluating dataset %s on SRF %d, progress: %3.2f', datasetName, SRF, percentDone);
                    fprintf([reverseStr, msg]);
                    reverseStr = repmat(sprintf('\b'), 1, length(msg));
                end
            end
            
            % Save results
            resName = ['result_', datasetName, '_SRF_', num2str(SRF),'.mat'];
            save(fullfile(resPath, resName), 'PSNR_table', 'SSIM_table', 'IFC_table');
        else
            load(fullfile(resPath, resName));
        end
        
        % Display results
        avgPSNR = mean(PSNR_table,1);
        avgSSIM = mean(SSIM_table,1);
        avgIFC  = mean(IFC_table, 1);
        fprintf('\n\n=== Quantitative results for dataset %s on SRF %d === \n\n', datasetName, SRF);
        fprintf('Peak signal-to-noise ratio (PSNR) \n')
        fprintf('     %8s\t%8s\t%8s\t%8s\t%8s\t%8s\t%8s\t%8s\t \n', methodList{1}, methodList{2}, ...
            methodList{3}, methodList{4}, methodList{5}, methodList{6}, methodList{7}, methodList{8});
        fprintf('PSNR|%8.02f\t|%8.02f\t|%8.02f\t|%8.02f\t|%8.02f\t|%8.02f\t|%8.02f\t|%8.02f\t| \n', ...
            avgPSNR(1), avgPSNR(2), avgPSNR(3), avgPSNR(4), avgPSNR(5), avgPSNR(6), avgPSNR(7), avgPSNR(8));
        fprintf('Structural similarity index (SSIM) \n')
        fprintf('     %8s\t%8s\t%8s\t%8s\t%8s\t%8s\t%8s\t%8s\t \n', methodList{1}, methodList{2}, ...
            methodList{3}, methodList{4}, methodList{5}, methodList{6}, methodList{7}, methodList{8});
        fprintf('SSIM|%8.04f\t|%8.04f\t|%8.04f\t|%8.04f\t|%8.04f\t|%8.04f\t|%8.04f\t|%8.04f\t| \n', ...
            avgSSIM(1), avgSSIM(2), avgSSIM(3), avgSSIM(4), avgSSIM(5), avgSSIM(6), avgSSIM(7), avgSSIM(8));
        fprintf('Information fidelity criterion (IFC) \n')
        fprintf('    %8s\t%8s\t%8s\t%8s\t%8s\t%8s\t%8s\t%8s\t \n', methodList{1}, methodList{2}, ...
            methodList{3}, methodList{4}, methodList{5}, methodList{6}, methodList{7}, methodList{8});
        fprintf('IFC|%8.02f\t|%8.02f\t|%8.02f\t|%8.02f\t|%8.02f\t|%8.02f\t|%8.02f\t|%8.02f\t| \n\n', ...
            avgIFC(1), avgIFC(2), avgIFC(3), avgIFC(4), avgIFC(5), avgIFC(6), avgIFC(7), avgIFC(8));
        fprintf('\n');
    end
end