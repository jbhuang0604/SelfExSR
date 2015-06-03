function [images P disps nclosest] = ojw_setup(options, Pout)
%OJW_SETUP  Cache data required for OJW's vision algorithms
%
%   [images P disps] = ojw_setup(options[, Pout])
%
% Caches the data structures required for OJW's image-based rendering and
% stereo functions.
%
% This function assumes that the current directory contains the input
% images, input.000.png etc., as well as a file called data.mat, which
% contains the input matrices, Pi, and either the 3D points, points, from
% the structure from motion (SfM) algotrithm used to calibrate the
% sequence, or the disparity values, disps, at which to sample.
%
%IN:
%   options - structure containing at least the following fields:
%      nclosest - 1xN array of indices of the input images to use in the
%                 rendered, indexed in order of distance of camera center
%                 from output camera center.
%      imout - positive integer indexing the view to be rendered. Only
%              required if size(Pout, 3) > 1.
%      dim_out - 1x4 array of output image dimensions, defined thus:
%                [start_x-1 start_y-1 width height].
%   Pout - 3x4xN projection matrices specifying the output view. If not
%          given, the input matrices are used.
%
%OUT:
%   images - 1xN cell array of input images.
%   P - 3x4xN array of input image projection matrices, relative
%       to the output image.
%   disps - 1xM vector of descending disparities to sample at in rendering.

% $Id: ojw_setup.m,v 1.6 2008/06/04 12:19:08 ojw Exp $

% Load the input variables we need
load data.mat

% Get the output matrix we want
if isempty(Pout)
    Pout = Pi;
end
if size(Pout, 3) > 1
    Pout = Pout(:,:,options.imout);
end

if isfield(options, 'dim_out') && ~isempty(options.dim_out)
    % Configure the output P matrix according to the new top-left pixel
    % position
    Pout = [1, 0, -options.dim_out(1); 0, 1, -options.dim_out(2); 0 0 1] * Pout;
end

% Make input matrices relative to output
[P nclosest] = P_r2j(Pout, Pi, options.nclosest);

if exist('points', 'var')
    % Determine range of depth in scene from SfM feature points in space
    points(4,:) = 1; % Homogenize
    points = Pout(3,:) * points;
    points = [min(points) max(points)];
    % Extend range by 20% front and back
    points = points .* [0.8 1.2];

    % Determine number of disparity levels
    % Project pixel 0,0 into all views, and ensure that minimum spacing
    % between samples is 0.5 pixels.
    WC = zeros(4, 2);
    WC(3,:) = points;
    WC(4,:) = 1;
    disp_vals = 0;
    for a = 1:size(P, 3)
        X = P(:,:,a) * WC;
        X = ojw_bsxfun(@times, X(1:2,:), 1./X(3,:));
        X = max(abs(X(:,1) - X(:,2)));
        disp_vals = max(disp_vals, X);
    end
    disp_vals = ceil(disp_vals * 2);

    % Calculate disparities
    disps = 0:disp_vals-1;
    disps = disps .* ((1-(points(1)/points(2)))/(disp_vals-1));
    disps = (1 - disps) ./ points(1);
end

% Order disparities from foreground to background
disps = sort(disps, 'descend');

% Load the input images
images = cell(size(nclosest));
for a = 1:numel(nclosest)
    images{a} = imread(sprintf('input.%3.3d.png', nclosest(a)-1));
end
return

