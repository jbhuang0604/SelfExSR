function options = ojw_default_options(algorithm)
%OJW_DEFAULT_OPTIONS  Initialize options for OJW's vision algorithms
%
%   options = ojw_default_options(algorithm)
%
% Outputs a data structure containing all options required for the given 
% image-based rendering or stereo function, initialized to their default
% values.
%
%IN:
%   algorithm - string specifying which rendering or stereo algorithm is to
%               be used. Options are:
%      'cvpr07' - Rendering algorithm from Woodford et al.'s CVPR 2007
%                 paper, "Efficient New-view Synthesis using Pairwise
%                 Dictionary Priors".
%      'bmvc07' - Rendering algorithm from Woodford et al.'s BMVC 2007
%                 paper, "On New View Synthesis Using Multiview Stereo".
%      'cvpr08' - Stereo algorithm from Woodford et al.'s CVPR 2008 paper,
%                 "Global Stereo Reconstruction under Second Order
%                 Smoothness Priors".
%
%OUT:
%   options - structure containing default values for all general and
%             algorithm specific options.

% $Id: ojw_default_options.m,v 1.9 2008/11/17 11:27:35 ojw Exp $

% Global defaults
options.show_output = usejava('awt'); % Display progress if possible
options.imout = 1; % The output image we want to create
options.nclosest = 1:8; % Index of input images (in order of proximity to output) to use
options.dim_out = []; % Overides the default output dimensions, unless empty
options.disp_vals = []; % Overides the default number of disparity values, unless empty

% Algorithm specific defaults
switch algorithm
    case {'cvpr07', 'edgemodes'}
        options.render_func = 'ibr_edgemodes';
        % Values from the CVPR '07 paper
        options.thresh = 50 / sqrt(3); % 50 for 3-channel images
        options.lambda = 1;
        options.connect = 8; % 8-connected is slower, but better quality, than 4
    case {'bmvc07', 'occlrender'}
        options.render_func = 'ibr_occlrender';
        % Energy parameters
        options.col_thresh = 12.5; % Colour threshold parameter for truncated quadratic kernel
        options.disp_thresh = 1.9; % Disparity threshold for the smoothness prior
        options.lambda = 0.24; % Weight of priors vs data likelihood
        options.tex_weight = 6; % Ratio of texture weight to smoothness
        options.tex_thresh = sqrt(5000) / 2; % Truncation threshold for texture prior
        % Implementation settings
        options.num_loops = 2; % Number of loops through proposals
        options.smoothness_kernel = 1; % 1: truncated linear kernel; 2: truncated quadratic kernel
        options.connect = 4; % 4 connected graph only
        options.contract = 0; % Use QPBOP with n contractions
        options.improve = 0; % Use QPBOI
        options.visibility = true; % Employ the geometrical visbility contsraint
    case {'cvpr08', 'stereo'}
        options.render_func = 'ojw_stereo';
        options.planar = true; % Use a planar prior
        % Default energy parameters
        options.disp_thresh = 0.02;
        options.smoothness_kernel = 1; % 1: truncated linear kernel; 2: truncated quadratic kernel
        options.col_thresh = 30;
        options.occl_const = 0.01;
        options.lambda_l = 9;
        options.lambda_h = 108;
        options.seg_params = [4 5 0];
        % Optimization settings
        options.connect = 4; % 4 connected (bi-directional) or 8 connected (quad-directional) graph
        options.contract = 0; % Use QPBOP with n contractions
        options.improve = 4; % 0: QPBO-F, 1: QPBOI-F, 2: QPBO-R, 3: QPBO-L, 4:QPBOI-R
        options.visibility = true; % Employ the geometrical visbility contstraint
        options.compress_graph = false; % Compression makes graph smaller, but is generally slower over all
        options.max_iters = 3000; % Maximum number of iterations
        options.converge = 0.01; % Loop until percentage decrease in energy per loop is less than options.converge (=> 101 == loop once)
        options.average_over = 20; % Number of iterations to average over when checking convergence
        options.independent = false; % Use strongly-connected, rather than independent, regions
        options.window = 2; % Half-size of window to use in window matching
        options.proposal_method = 1:3; % SameUni, SegPln, Smooth
        % Data gathering parameters
        options.save_name = []; % Filename to save frames in
    otherwise
        error('Algorithm not recognised');
end
return