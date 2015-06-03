% Load default options for the desired algorithm
algs = {'cvpr07', 'bmvc07', 'cvpr08'};
options = ojw_default_options(algs{3});

% In the case of stereo, download the desired sequence
sequence = 'cones';
download_stereo(sequence);
cd(sequence)

% Change any options here
options.dim_out = [0 0 450 375]; % Output image dimensions: [start_x-1 start_y-1 width height]
options.imout = 3; % Index of projection matrix to use for output
options.nclosest = [1 7]; % Input images to use, in terms of distance of camera centres from output view

% Generate output matrices
Pout = zeros(0, 0, 1); % Use input matrix defined by options.im_out
%Pout = ojw_genview('steady', [1 63], (0:29)/29); % 30 frame steadicam sequence 
%Pout = ojw_genview('stereo', 1:63, 'l'); % Generate left views of stereo pairs

for a = 1:size(Pout, 3)
    % Call the rendering function
    [A out] = ibr_render(options, Pout(:,:,a));

    % Save the output
    clf; sc(A);
end