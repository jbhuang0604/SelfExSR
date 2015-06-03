function P = P_interp(first, last, frames)
%P_INTERP  Linear interpolation between 2 projection matrices
%
%   P = P_interp(first, last, frames)
%
% Creates a series off projection matrices interpolated from two projection
% matrices.
%
%IN:
%   first - 3x4 The first projection matrix in the series, and the first
%           projection matrix between which the output matrices are
%           interpolated from.
%   last - 3x4 The second projection matrix between which the output
%          matrices are interpolated from.
%   frames - Nx1 The frames to interpolate along a line, where 0 is the
%            position of first and 1 the position of last.
%
%OUT:
%   P - 3x4xN The output projection matrices.

% $Id: P_interp.m,v 1.5 2008/05/09 09:48:34 ojw Exp $

[k_first r_first t_first] = KR_from_P(first);
[k_last r_last t_last] = KR_from_P(last);

% If K matrices are identical but for sign, put the sign change into the R
% matrices
a = diag(k_first) ./ diag(k_last);
b = sign(a);
if all(abs(a-b)<1e-8) && any(b==-1)
    a = diag(sign(b+0.5));
    k_first = k_first * a;
    r_first = a * r_first;
end
clear a b

% Calculate our step sizes
t_step = t_last - t_first;
r_step = r_first' * r_last;
k_step = k_last - k_first;

P = repmat(eye(3, 4), [1 1 numel(frames)]);
% Interpolate
for a = 1:numel(frames)
    % Translation
    P(:,4,a) = -t_first - t_step * frames(a);
    % Rotation
    P(:,:,a) = r_first * real(r_step ^ frames(a)) * P(:,:,a);
    % Calibration
    P(:,:,a) = (k_first + k_step * frames(a)) * P(:,:,a);
end

% Projection matrix decomposition
function [K R t] = KR_from_P(P)
t = -P(:,1:3) \ P(:,end);

st = @(M) M(end:-1:1,end:-1:1)';
[R K] = qr(st(P(:,1:3)));
R = st(R);
K = st(K);
