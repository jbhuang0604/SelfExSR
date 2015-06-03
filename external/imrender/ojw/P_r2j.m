function [Pj closest] = P_r2j(Pr, Pj, closest)
%P_R2J  Change reference frame of projection matrices 
%
%   [Pk I] = P_r2j(Pr, Pj, closest)
%
% Tranforms the coordinate frame of a set of projection matrices such that
% a reference projection matrix becomes the identity matrix eye(3, 4). Can
% also select input matrices based on camera centre proximity to the
% reference frame.
%
%IN:
%   Pr - 3x4 projection matrix of the reference frame.
%   Pj - 3x4xN array of input projection matrices.
%   closest - 1xP indices of matrices of Pj to output, once arranged in
%             order of proximity of camera centre to Pr. Default: 1:N.
%
%OUT:
%   Pk - 3x4x(min(N,P)) array of transformed projection matrices.
%   I - 1x(min(N,P)) indices of matrices of Pj that make up Pk, such
%       that Pk = T(Pj(:,:,I)), T() being the coordinate transformation.

% $Id: P_r2j.m,v 1.2 2007/12/10 10:58:31 ojw Exp $

if nargin > 2
    % Determine which Pj matrices are required based on distance of optical
    % centres from that of Pr.
    T = zeros(3, size(Pj, 3));
    for a = 1:size(Pj, 3)
        T(:,a) = -Pj(:,1:3,a) \ Pj(:,4,a);
    end
    order = ojw_bsxfun(@minus, T, -Pr(:,1:3) \ Pr(:,4));
    [T order] = sort(sum(order .^ 2));
    closest = order(closest);
    Pj = Pj(:,:,closest);
end

% Calculate Pj relative to Pr
% Assume projection matrices make 0,0 the top left corner of the top left
% pixel. Convert to Matlab form, i.e. 0.5,0.5 is the top left corner of the
% top left pixel.
T = [1 0 0.5; 0 1 0.5; 0 0 1];
Pr = inv([T * Pr; 0 0 0 1]);
for a = 1:size(Pj, 3)
    Pj(:,:,a) = T * Pj(:,:,a) * Pr;
end
return
    