function [D info] = ojw_stereo(images, P, disps, sz, options)
%OJW_STEREO  Global stereo with 2nd-order smoothness prior & occlusion model
%
%   [D info] = ojw_stereo(images, P, disps, sz, options)
%
% Generates a disparity (1/depth) map for the first input image. This
% algorithm implements a "global" stereo algorithm, with asymmetrical
% occlusion modelling, using an alpha-expansion graph cuts style approach,
% but with arbitrary disparity proposals.
%
%IN:
%   images - 1xN cell array of input images, the reference image being
%            images{1}.
%   P - 3x4xN array of projection matrices for the input images, relative
%       to the output image.
%   disps - 1xM list of disparities to sample at.
%   sz - 1x2 vector of output image dimensions: [H W].
%   options - a structure containing the following input parameters:
%       col_thresh - scalar noise parameter for data likelihood.
%       occl_const - scalar occlusion cost.
%       disp_thresh - scalar disparity threshold for smoothness prior.
%       smoothness_kernel - index denoting which smoothness kernel to use.
%                           1: truncated linear; 2: truncated quadratic.
%       lambda_l - scalar smoothness prior weight for cliques crossing
%                  segmentation boundaries.
%       lambda_h - scalar smoothness prior weight for cliques not crossing
%                  segmentation boundaries.
%       seg_params - 1x3 vector of parameters for the mean-shift
%                    over-segmentation of the reference image.
%       visibility - boolean indicating whether to employ the geometrical
%                    visbility contstraint.
%       connect - scalar neighbourhood system of the graph, 4 or 8
%                 connected.
%       max_iters - scalar number of iterations to halt after, if
%                   convergence is not achieved first.
%       converge - scalar percentage decrease in energy per iteration at
%                  which optimization stops.
%       average_over - scalar number of iterations to average over when
%                      checking convergence.
%       contract - scalar number of QPBOP iterations to do.
%       improve - scalar indicating which method to use to label unlabelled
%                 nodes. 0: QPBO-F, 1: QPBOI-F, 2: QPBO-R, 3: QPBO-L,
%                 4: QPBOI-R.
%       independent - boolean indicating whether to use independent, or
%                     merely strongly-connected, regions for improve
%                     methods 2 & 4.
%
%OUT:
%   D - HxW disparity map.
%   info - structure containing other outputs from the algorithm.

% $Id: ojw_stereo.m,v 1.5 2008/11/17 11:27:35 ojw Exp $

% Crude check for a reference image
if max(abs(P([1:6 9])-[1 0 0 0 1 0 1])) > 1e-12
    error('First image must be reference image');
end

% Initialize data arrays
R = images{1}(round(P(8))+(1:sz(1)),round(P(7))+(1:sz(2)),:);
vals.I = images(2:end);
vals.P = permute(P(:,:,2:end), [2 1 3]);
vals.sz = sz;
colors = size(R, 3);
num_in = numel(images);
Rorig = uint8(R);
if colors == 1
    Rorig = repmat(Rorig, [1 1 3]);
end
vals.R = repmat(reshape(single(R), [], colors), [2 1]);
vals.d_min = disps(end);
vals.d_step = disps(1) - vals.d_min;
vals.ndisps = numel(disps);

T = reshape(uint32(1:prod(sz)), sz);
if options.planar
    % Use 2nd order smoothness prior
    SEI = [reshape(T(1:end-2,:), 1, []) reshape(T(:,1:end-2), 1, []); ...
           reshape(T(2:end-1,:), 1, []) reshape(T(:,2:end-1), 1, []); ...
           reshape(T(3:end,:), 1, []) reshape(T(:,3:end), 1, [])];
    if options.connect == 8
        SEI = [SEI [reshape(T(1:end-2,1:end-2), 1, []) reshape(T(3:end,1:end-2), 1, []); ...
                    reshape(T(2:end-1,2:end-1), 1, []) reshape(T(2:end-1,2:end-1), 1, []); ...
                    reshape(T(3:end,3:end), 1, []) reshape(T(1:end-2,3:end), 1, [])]];
    end
else
    % Use 1st order smoothness prior
    SEI = [reshape(T(1:end-1,:), 1, []) reshape(T(:,1:end-1), 1, []); ...
           reshape(T(2:end,:), 1, []) reshape(T(:,2:end), 1, [])];
    if options.connect == 8
        SEI = [SEI [reshape(T(1:end-1,1:end-1), 1, []) reshape(T(2:end,1:end-1), 1, []); ...
                    reshape(T(2:end,2:end), 1, []) reshape(T(1:end-1,2:end), 1, [])]];
    end
end
clear T

% Initialise display
vals.show_output = options.show_output;
if vals.show_output
    vals.show_output = gcf;
    set(0, 'CurrentFigure', vals.show_output);
    subplot('Position', [0 0.5 1/3 0.5]);
    sc(R, [0 255]);
end

% Segment the image using mean shift
info.segment = vgg_segment_ms(Rorig, options.seg_params(1), options.seg_params(2), options.seg_params(3));
% Find smoothness edges which don't cross segmentation boundaries
EW = reshape(~any(diff(int32(info.segment(SEI))), 1), 1, []);
EW = EW * options.lambda_h + ~EW * options.lambda_l;
EW = EW * (num_in / ((options.connect==8) + 1));
EW = reshape(repmat(EW, [4*(1+(options.planar~=0)) 1]), [], 1);

% Set up values for ibr_fuse_depths
vals.visibility = (options.visibility ~= 0) * 1e4;
vals.improve = options.improve;
vals.contract = options.contract;
vals.independent = options.independent;
vals.compress_graph = options.compress_graph;

% Set up our robust kernels
vals.ephoto = @(F) log(2) - log(exp(sum(F .^ 2, 2)*(-1/(options.col_thresh*colors)))+1);
switch options.smoothness_kernel
    case 1
        vals.esmooth = @(F) EW .* min(abs(F), options.disp_thresh);
    case 2
        EW = EW / options.disp_thresh;
        vals.esmooth = @(F) EW .* min(F.^2, options.disp_thresh^2);
    otherwise
        error('Unknown smoothness kernel specified');
end
vals.occl_val = options.occl_const + log(2);
vals.SEI = SEI;
clear T SEI EW Rorig

if nargout > 1
    % Save parameters
    info.params.disp_thresh = options.disp_thresh;
    info.params.col_thresh = options.col_thresh;
    info.params.occl_const = options.occl_const;
    info.params.lambda_l = options.lambda_l;
    info.params.lambda_h = options.lambda_h;
end

if isnumeric(options.proposal_method) && size(options.proposal_method, 1) == 1
    % Use the proposal methods:
    for a = options.proposal_method
        switch a
            case 0
                % Ordered fronto-parallel
                [D info.samedisc_optim] = ojw_stereo_optim(vals, @(n) 3, options);
                info.samedisc_optim.D = D;
            case 1
                % SameUni (random fronto-parallel)
                [D info.sameuni_optim] = ojw_stereo_optim(vals, @(n) 1, options);
                info.sameuni_optim.D = D;
            case 2
                % SegPln (prototypical segment-based stereo proposals)
                [Dproposals info.segpln_gen] = ojw_segpln(images, P, disps, R, options);
                clear R
                Dproposals = @(n) Dproposals(:,:,mod(n-1, size(Dproposals, 3))+1);
                [D info.segpln_optim] = ojw_stereo_optim(vals, Dproposals, options);
                clear Dproposals
                info.segpln_optim.D = D;
            case 3
                % Smooth*
                Dproposals = {info.segpln_optim.D, info.sameuni_optim.D, 2, 2, 2, 2};
                Dproposals = @(n) Dproposals{mod(n-1, 6)+1};
                [D info.smooth_optim] = ojw_stereo_optim(vals, Dproposals, options);
                clear Dproposals
                info.smooth_optim.D = D;
            case 4
                % Smooth
                Dproposals = {D, 2};
                Dproposals = @(n) Dproposals{(n>1)+1};
                [D info.smooth2_optim] = ojw_stereo_optim(vals, Dproposals, options);
                clear Dproposals
                info.smooth2_optim.D = D;
        end
    end
else
    % Input fixed set of proposals
    if isnumeric(options.proposal_method)
        Dproposals = @(n) options.proposal_method(:,:,mod(n-1, size(options.proposal_method, 3))+1);
    else
        Dproposals = options.proposal_method;
    end
    [D info.udprop_optim] = ojw_stereo_optim(vals, Dproposals, options);
    clear Dproposals
    info.udprop_optim.D = D;
end
return