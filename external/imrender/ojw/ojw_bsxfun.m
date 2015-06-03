function C = ojw_bsxfun(func, A, B)
%OJW_BSXFUN  Provides backwards compatibility for Matlab's bsxfun
%
%   C = ojw_bsxfun(func, A, B)
%
% Each dimension of A and B must be equal to each other, or equal to one.
% Whenever a dimension of one of A or B is singleton (equal to 1), the
% array is virtually replicated along that dimension to match the other
% array (or diminished if the corresponding dimension of the other array is
% 0). The two resized arrays are then passed, as vectors, to func(), which
% must be able to accept as input either two column vectors of the same
% size, or one column vector and one scalar, and return as output a column
% vector of the same size as the input(s).
%
% IN:
%    func - Function handle.
%    A - Input array.
%    B - Input array of same type as A.
%
% OUT:
%    C - max(size(A),size(B)).*(size(A)>0 & size(B)>0) sized output array.
%
% See also BSXFUN.

% $Id: ojw_bsxfun.m,v 1.1 2007/12/07 11:27:52 ojw Exp $

if exist('bsxfun', 'builtin') 
    % Do the obvious - Matlab versions below 7.3 don't have it though
    C = bsxfun(func, A, B);
else
    % Backwardly compatible method
    if isscalar(A)
        C = reshape(func(A, B(:)), size(B));
        return
    elseif isscalar(B)
        C = reshape(func(A(:), B), size(A));
        return
    end
    sA = size(A);
    sB = size(B);
    % Enlarge the smaller of the two sizes
    s = numel(sA) - numel(sB);
    if s < 0
        sA = [sA ones(1, -s)];
    elseif s > 0
        sB = [sB ones(1, s)];
    end
    % Calculate the output array size
    sMax = max(sA, sB) .* (sA > 0 & sB > 0);
    % Make sure arrays have same size of dimension or size 1
    a = sA ~= sMax;
    b = sB ~= sMax;
    if any(a & sA ~= 1) || any(b & sB ~= 1)
        error('A and B must have equal or singleton dimensions');
    end
    if ~all(sMax)
        % Some entries are zero, so array is empty
        C = zeros(sMax);
        return
    end
    % Resize the arrays to have the same size, sMax
    if any(a)
        A = repmat(A, sMax ./ sA);
    end
    if any(b)
        B = repmat(B, sMax ./ sB);
    end
    % Apply the function
    C = reshape(func(A(:), B(:)), sMax);
end