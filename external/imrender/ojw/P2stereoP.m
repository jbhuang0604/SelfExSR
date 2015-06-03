function P = P2stereoP(P, Zmin, view)
%P2STEREOP  Generates horizontally shifted viewpoints
%
%   Po = P2stereoP(Pi, Zmin[, view])
%
% Generates horizontally shifted viewpoints, suitable for viewing in
% stereo.
%
%IN:
%   Pi - 3x4xN input projection matrices.
%   Zmin - scalar, minimum depth in the scene.
%   view - scalar value of horizontal offset, where -1 is optimal left view
%          and 1 is optimal right view, relative to the input matrix.
%          Default: -1.
%
%OUT:
%   Po - 3x4xN output projection matrices.

% $Id: P2stereoP.m,v 1.5 2008/06/04 10:50:55 ojw Exp $

% Determine whether to render left or right view
if nargin < 2
    view = -1;
end
if ischar(view)
    % For 'l', 'm' & 'r' (left, middle and right)
    view = sign(double(view(1)) - double('m'));
end

% Calculate the new projection matrix
P(1,4,:) = P(1,4,:) + (view * 70 * Zmin); % Shift the camera centre
P(1,:,:) = P(1,:,:) - P(3,:,:) * (view * 50); % Shift the image plane
return