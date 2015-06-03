function [D info] = ojw_stereo_optim(vals, Dproposals, options)
%OJW_STEREO_OPTIM  Optimization for global stereo
%
%   [D info] = ojw_stereo_optim(vals, Dproposals, options)
%
% Given a binary energy function for a stereo problem, and a set of
% proposals (or proposal index), fuses these proposals until convergence of
% the energy.
%
%IN:
%   vals - structure containing the following data relevant to the stereo
%          problem:
%       R - (H*W*2)xC columnated reference image (twice).
%       I - 1xN cell array of input images.
%       P - 4x3xN array of transposed projection matrices.
%       sz - 1x2 vector of output image dimensions: [H W].
%       d_min - scalar indicating the minimum disparity in the scene.
%       d_step - scalar indicating the range of disparities in the scene.
%       ephoto - function handle to data cost energy function.
%       esmooth - function handle to smoothness prior energy function.
%       occl_val - scalar penalty energy cost for occluded pixels.
%       SEI - Mx(2 or 3) uint32 array of smoothness clique indices.
%       visibility - boolean indicating whether to employ the geometrical
%                    visbility contstraint.
%       contract - scalar number of QPBOP iterations to do.
%       improve - scalar indicating which method to use to label unlabelled
%                 nodes. 0: QPBO-F, 1: QPBOI-F, 2: QPBO-R, 3: QPBO-L,
%                 4: QPBOI-R.
%       independent - boolean indicating whether to use independent, or
%                     merely strongly-connected, regions for improve
%                     methods 2 & 4.
%   Dproposals - function handle: Dproposals(iter) returns a proposal
%                disparity map, or one of the following proposal methods:
%                0 - Dp = rand(H, W).
%                1 - Dp = repmat(rand(1), [H W]).
%                2 - Dp = (D(1:end-2,:) + D(3:end,:)) / 2 if iter is odd
%                       & (D(:,1:end-2) + D(:,3:end)) / 2 if iter is even
%   options - a structure containing the following input parameters:
%       max_iters - scalar number of iterations to halt after, if
%                   convergence is not achieved first.
%       converge - scalar percentage decrease in energy per iteration at
%                  which optimization stops.
%       average_over - scalar number of iterations to average over when
%                      checking convergence.
%
%OUT:
%   D - HxW disparity map.
%   info - structure containing other outputs from the algorithm.

% $Id: ojw_stereo_optim.m,v 1.2 2008/11/17 11:27:35 ojw Exp $

% Create initial arrays
info.map = zeros(vals.sz, 'uint16');
max_iters = options.max_iters + options.average_over + 1;
info.energy = zeros(max_iters, numel(vals.improve));
info.numbers = zeros(4, max_iters, numel(vals.improve), 'uint32'); % Number [updated; unlabelled; independent regions]
info.timings = zeros(4, max_iters, numel(vals.improve)); % Cumulative timings of [proposal; data; smoothness; optimization; finish]

% Initialise depth map
D = rand(vals.sz) * vals.d_step + vals.d_min;

options.converge = options.converge * 0.01 * options.average_over;
iter = options.average_over + 1;
info.energy(1:options.average_over,end) = realmax;
info.energy(iter,end) = realmax / 1e20;
while (1-(info.energy(iter,end)/info.energy(iter-options.average_over,end)) > options.converge) && iter < max_iters
    iter = iter + 1;
    % Set the new (proposal) depth map
    t_start = cputime; % Start timing
    Dnew = Dproposals(iter-(options.average_over+1));
    if isscalar(Dnew)
        switch Dnew
            case 0
                % Completely random
                Dnew = rand(vals.sz);
            case 1
                % Random fronto-parallel
                Dnew = rand(1);
            case 2
                % Interpolate disparity over rows and columns
                Dnew = (D - vals.d_min) / vals.d_step;
                if mod(iter, 2)
                    Dnew(2:end-1,:) = (Dnew(1:end-2,:) + Dnew(3:end,:)) / 2;
                else
                    Dnew(:,2:end-1) = (Dnew(:,1:end-2) + Dnew(:,3:end)) / 2;
                end
            case 3
                % Ordered fronto-parallel, front-to-back
                Dnew = 1 - mod(iter-options.average_over-2, vals.ndisps) / (vals.ndisps - 1);
            case -3
                % Ordered fronto-parallel, back-to-front
                Dnew = mod(iter-options.average_over-2, vals.ndisps) / (vals.ndisps - 1);
            otherwise
                % Quit
                if iter == (options.average_over + 2)
                    [M info info.energy(iter) V] = ibr_fuse_depths(D, D, vals);
                else
                    iter = iter - 1;
                end
                break;
        end
        Dnew = Dnew * vals.d_step + vals.d_min;
    else
        Dnew = double(Dnew);
    end
    Dnew(~(Dnew>vals.d_min)) = vals.d_min;
    if isscalar(Dnew)
        Dnew = repmat(Dnew, vals.sz);
    end
    info.timings(1,iter,:) = cputime - t_start; % Time proposal generation

    % Fuse the depths
    try
        [M stats info.energy(iter,:) V] = ibr_fuse_depths(D, Dnew, vals);
    catch
        if iter == (options.average_over + 2)
            rethrow(lasterror);
        end
        % Fusing failed - perhaps user bailed early.
        % Output current state.
        info.error = lasterr;
        iter = iter - 1;
        break;
    end
    D(M) = Dnew(M);
    info.map(M) = iter;
    info.timings(2:4,iter,:) = stats.timings + info.timings(1,iter,end);
    info.numbers(:,iter,:) = stats.numbers;
    
    % Record progress of disparity map
    save_progress(options.save_name, 'D');
end

if nargout > 1
    % Save stats
    info.energy = info.energy(options.average_over+2:iter,:);
    info.numbers = info.numbers(:,options.average_over+2:iter,:);
    info.timings = info.timings(:,options.average_over+2:iter,:);
    % Save visibility & mapping
    info.vis = reshape(V, vals.sz(1), vals.sz(2), []);
    info.map = info.map - uint16(options.average_over + 1);
end
return