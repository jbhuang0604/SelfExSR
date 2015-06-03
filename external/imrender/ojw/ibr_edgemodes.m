function [A info] = ibr_edgemodes(images, P, disps, sz, options)
%IBR_EDGEMODES  Rendering method from OJW's CVPR 2007 paper
%
%   [A info] = ibr_edgemodes(images, P, disps, sz, options)
%
% Generates an output image by finding colour modes for each output pixel,
% then choosing between these modes using texture priors.
%
% This is an implementation of the algorithm described in Woodford et al.'s
% CVPR 2007 paper:
%   "Efficient New-view Synthesis using Pairwise Dictionary Priors".
%
%IN:
%   images - 1xN cell array of input images.
%   P - 3x4xN array of input image projection matrices, relative
%       to the output image.
%   disps - 1xM vector of descending disparities to sample at.
%   sz - 1x2 vector of output image dimensions: [H W].
%   options - a structure containing the following input parameters:
%      lambda - scalar indicating weight of texture prior relative to data
%               likelihood term.
%      thresh - truncation threshold (per channel) for colour modes.
%      connect - connection system of pairwise graph (4 or 8 only).
%
%OUT:
%   A - HxWxC rendered image.
%   info - structure containing the following additional output infomation:
%      D - HxW disparities at which the selected modes are from.

% $Id: ibr_edgemodes.m,v 1.5 2008/11/20 21:35:48 ojw Exp $

% Initialize variables
colors = size(images{1}, 3);
nDisps = numel(disps);
num_in = numel(images);
options.thresh = colors * options.thresh ^ 2;
P = permute(P, [2 1 3]);
I = zeros(colors, num_in, nDisps, sz(1), 2); % Depth/intensity table
mode_cell = cell(sz(1), 2);
% Edge indices
T = reshape(uint32(1:sz(1)*2), [sz(1) 2]);
EIcon = [T(1:end-1,1) T(2:end,1); T(:,1)  T(:,2)];
if options.connect == 8
    EIcon = [EIcon; T(1:end-1,1) T(2:end,2); T(2:end,1) T(1:end-1,2)];
end
% Swap and sort to make use of cacheing in truncquad_edges
[M N] = ind2sub(sz, EIcon);
M = diff(mod(int32(M+N), 2), 1, 2) < 0;
EIcon(M,:) = EIcon(M,[2 1]);
EIcon = sortrows(EIcon)';
clear T M N
WC = ones(nDisps*sz(1), 4); % Depth sampling points in world coordinates
% Initialise homogenous coordinates
WC(:,4) = repmat(disps(:), [sz(1) 1]);
WC(:,2) = reshape(repmat(1:sz(1), [nDisps 1]), [], 1);

data = cell(sz);
if options.connect == 8
    num_edges = 4*sz(1)*sz(2)-3*(sz(2)+sz(1))+2;
else
    num_edges = 2*sz(1)*sz(2)-sz(2)-sz(1);
end
edges = cell(num_edges, 1);
EI = zeros(2, num_edges, 'uint32');
k = 1;

t_start = cputime;
% Loop over every pixel to be rendered
for x = 1:sz(2)
    % Set x coordinate of world coordinates
    WC(:,1) = x;

    % For each input image...
    for a = 1:num_in
        % Calculate the coordinates in the input image
        N = WC * P(:,:,a);
        N(:,3) = 1 ./ N(:,3);
        N(:,1) = N(:,1) .* N(:,3);
        N(:,2) = N(:,2) .* N(:,3);

        % Look-up the colours for each point
        I(:,a,:,:,2) = reshape(squeeze(vgg_interp2(images{a}, N(:,1), N(:,2), 'linear', -1000))', colors, 1, nDisps, sz(1));
    end
    
    for y = 1:sz(1)
        % Calculate the modes for each pixel
        [modes depth energy] = truncquad_modes(I(:,:,:,y,2), options.thresh, 0, 1e4);
        data{y,x} = [energy; depth; modes];
        mode_cell{y,2} = modes;
    end
        
    if options.lambda > 0 && x > 1
        % Calculate the edge costs in one go
        E = truncquad_edges(reshape(I, colors, num_in*nDisps, sz(1)*2), mode_cell, EIcon, 1e100, options.lambda);
        k2 = k + numel(E) - 1;
        edges(k:k2) = E(:);
        EI(:,k:k2) = EIcon + uint32((x-2)*sz(1));
        k = k2+1;
    end
    I = circshift(I, [0 0 0 0 -1]);
    mode_cell = circshift(mode_cell, [0 -1]);
    vgg_progressbar('Computing colour modes and edge costs...', x/sz(2));
end
if options.lambda > 0
    % Calculate the edge costs for the last column
    EIcon = uint32([1:sz(1)-1; 2:sz(1)]);
    E = truncquad_edges(reshape(I, colors, num_in*nDisps, sz(1)*2), mode_cell, EIcon, 1e100, options.lambda);
    k2 = k + numel(E) - 1;
    edges(k:k2) = E(:);
    EI(:,k:k2) = EIcon + uint32((sz(2)-1)*sz(1));
end
info.times(1) = cputime - t_start;

if options.lambda > 0
    % Optimization using TRW
    [L info.energy info.lowerBound] = vgg_trw_bp(data, EI, edges, int32([1 0 100])); drawnow; % Write out result
    [A D] = slice_cell_image(data, L);
else
    % Just select the lowest energy mode
    [A D] = slice_cell_image(data);
end
info.times(2) = cputime - t_start;
A = permute(A, [2 3 1]);
if nargout > 1
    info.D = disps(D);
end
return
