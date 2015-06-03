function Pout = ojw_genview(type, ind, offset)
%OJW_GENVIEW  Generate output projection matrices
%
%   Pout = ojw_genview(type, ind, offset)
%
% Produces various types of output matrices for a scene - same as input
% (e.g. for stereo), horizontally shifted from input (e.g. for 3D display)
% and linearly interpolated between two input frames (e.g. for stabilizing
% a video).
%
% This function assumes that the current directory contains a file called
% data.mat, which contains the input matrices, Pi, and either the 3D
% points, points, from the structure from motion (SfM) algotrithm used to
% calibrate the sequence, or the disparity values, disps, at which to
% sample.
%
% IN:
%   type - string, one of the following:
%      'input' - output the input projection matrices.
%      'stereo' - output the input matrices horizontally displaced.
%      'steady' - output a linear interpolation of input matrices.
%   ind - 1xN list of indices of input matrices to use for 'input' and
%         'stereo', or 1x2 list of indices for 'steady', indicating the
%         input matrices to interpolate between.
%   offset - scalar used to define the horizontal offset for 'stereo' (0 =
%            input view, -1 = left view, 1 = right view, other values shift
%            view proportionally), or 1xN list of distances along
%            interpolated path for 'steady' (0 = ind(1) view, 1 = ind(2)
%            view, other values shift views proportionally.
%
% OUT:
%   Pout - 3x4xN output projection matrices.

% $Id: ojw_genview.m,v 1.4 2008/06/04 12:19:08 ojw Exp $

% Load the input data
load data.mat

% Produce the output projection matrix
switch type
    case 'input'
        Pout = Pi(:,:,ind);
    case 'stereo'
        Pout = Pi(:,:,ind);
        
        % Determine minimum scene depth
        if exist('disps', 'var')
            Zmin = 1 / max(disps(:));
        else
            % Determine range of depth in scene from SfM feature points in space
            points(4,:) = 1; % Homogenize
            points = Pout(3,:,1) * points;
            Zmin = min(points) * 0.8; % Extend range by 20% at front
        end
        
        % Generate stereo views
        Pout = P2stereoP(Pout, Zmin, offset);
    case 'steady'
        % Linear interpolation between ind(1) and ind(2)
        Pout = P_interp(Pi(:,:,ind(1)), Pi(:,:,ind(2)), offset); 
    otherwise
        error('Type not recognised');
end

return