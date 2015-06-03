function [A D E N] = slice_cell_image(data, L, F)
% Given data structure and a labelling, produces an output image A
if ~iscell(data)
    error('data must be a cell array');
end
if nargin > 1 && numel(L) ~= numel(data)
    error('data and L must be the same size');
end
E = zeros(size(data));
D = zeros(size(data));
N = zeros(size(data));
A = zeros([size(data{1}, 1)-2 size(data)]);
if nargin < 2
    for a = 1:numel(data)
        N(a) = size(data{a}, 2);
        if N(a) > 0
            [B L] = min(data{a}(1,:));
            B = data{a}(:,L);
            E(a) = B(1);
            D(a) = B(2);
            A(:,a) = B(3:end);
        end
    end
else
    for a = 1:numel(data)
        B = data{a}(:,L(a));
        E(a) = B(1);
        D(a) = B(2);
        A(:,a) = B(3:end);
    end
end
if nargin > 2
    A = permute(A, [2 3 1]);
    A = ibr_integrate_responses(A, F, [2 20]);
end
return