function [A info] = ibr_occlrender(images, P, disps, sz, options)
%IBR_OCCLRENDER  Rendering method from OJW's BMVC 2007 paper
%
%   [A info] = ibr_occlrender(images, P, disps, sz, options)
%
% Implementation of the image based rendering algorithm described in
% Woodford etal.'s BMVC 2007 paper:
%   "On New View Synthesis Using Multiview Stereo".
% 
% Has explicit depth reconstruction in order to provide geometrical
% occlusion modelling. Uses both smoothness and texture (optional)
% regularization.
%
%IN:
%   images - 1xN cell array of input images.
%   P - 3x4xN array of projection matrices for the input images, relative
%       to the output image.
%   disps - 1xM list of disparities to sample at.
%   sz - 1x2 vector of output image dimensions: [H W].
%   options - a structure containing the following input parameters:
%       col_thresh - scalar colour threshold parameter for data likelihood.
%       lambda - scalar smoothness prior weight.
%       disp_thresh - scalar disparity threshold for smoothness prior, in
%                     terms of the number of disparity levels. 
%       tex_weight - scalar weighting of the texture prior. 0 means no
%                    texture prior.
%       tex_thresh - scalar threshold multiplier for texture prior.
%       visibility - boolean indicating whether to employ a geometrical
%                    visibility constraint.
%
%OUT:
%   A - HxWxC rendered image.
%   info - structure containing the following additional output infomation:
%      D - HxW disparity map for output image
%      V - HxWxN boolean array indicating the visibility of the output
%          pixels in the input images. 

% $Id: ibr_occlrender.m,v 1.4 2008/07/16 16:59:20 ojw Exp $

% Set Ephoto parameters
colors = size(images{1}, 3);
num_in = numel(images);
options.col_thresh = options.col_thresh * num_in / (num_in - 1);
Kocc = options.col_thresh ^ 2 * colors;

% Set Esmooth parameters
options.disp_thresh = options.disp_thresh * mean(abs(diff(disps)));
if options.smoothness_kernel == 2
    options.disp_thresh = options.disp_thresh ^ 2;
end
options.lambda = options.lambda * Kocc * num_in / options.disp_thresh;

% Set Etexture parameters
if isempty(options.tex_thresh)
    options.tex_thresh = options.col_thresh;
end
options.tex_thresh = options.tex_thresh ^ 2 * colors * 2;
options.tex_weight = options.tex_weight / options.tex_thresh;

% Create initial arrays
tp = prod(sz);
nDisps = numel(disps);
out.energy = [];
out.ul = [];
out.time_fuse = [];
X = ones(sz(1), 1) * (1:sz(2));
Y = (1:sz(1))' * ones(1, sz(2));
WC = ones(tp, 3);
WC(:,1) = X(:);
WC(:,2) = Y(:);
clear X Y
P = permute(P, [2 1 3]);

if options.tex_weight
    % Cache all the input samples
    I = class(images{1});
    oobv = cast(-1000, I);
    I = zeros(colors, nDisps, num_in, sz(1)*sz(2), I);

    % For each image...
    for a = 1:num_in
        % Calculate image coordinates
        T = WC * P(1:3,:,a);
        % For each disparity...
        for b = 1:nDisps
            % Vary image coordinates according to disparity
            d = disps(b) * P(4,:,a);
            Z = 1 ./ (T(:,3) + d(3));
            % Sample the colours
            I(:,b,a,:) = reshape(squeeze(vgg_interp2(images{a}, (T(:,1) + d(1)) .* Z, (T(:,2) + d(2)) .* Z, 'linear', oobv))', colors, 1, 1, []);
        end
    end
    I = reshape(I, colors, nDisps*num_in, sz(1)*sz(2));
end

% Depth sampling points in world coordinates
WC(:,4) = 0;
WC = repmat(WC, [2 1]);

% Graph constants
Kocc = int32(Kocc);
Kinf = cast(2^28, class(Kocc));
oobv = single(-1000);

% Surface smoothness and texture cliques
T = reshape(uint32(1:tp), sz);
% Use fronto-parallel prior
SEI = [reshape(T(1:end-1,:), [], 1) reshape(T(2:end,:), [], 1);...
       reshape(T(:,1:end-1), [], 1) reshape(T(:,2:end), [], 1)];
if options.connect == 8
    SEI = [SEI; reshape(T(1:end-1,1:end-1), [], 1) reshape(T(2:end,2:end), [], 1);...
                reshape(T(2:end,1:end-1), [], 1) reshape(T(1:end-1,2:end), [], 1)];
    % Account for their being twice as many edges
    options.lambda = options.lambda / 2;
end
% Swap and sort to make use of cacheing in truncquad_edges
SEI_2nd = SEI(:,2);
[M N] = ind2sub(sz, SEI);
M = diff(mod(int32(M+N), 2), 1, 2) < 0;
SEI(M,:) = SEI(M,[2 1]);
[SEI M] = sortrows(SEI);
SEI = SEI';

% Visualization stuff
if options.show_output
    SEI_2nd = SEI_2nd(M);
    options.show_output = gcf;
else
    clear SEI_2nd
end
clear M N

% Initialise depth map
D = repmat(disps(1), sz);
    
for loop = 1:options.num_loops
    Dold = D;
    for d = disps(1+(loop==1):end)
        % Set the new depth map
        Dnew = repmat(d, sz);

        % Fuse the depth maps
        tic;
        % Calculate the homogenous coordinates of our two labellings
        WC(:,4) = [D(:); Dnew(:)];

        % Initialise arrays for the image samples and occlusion terms
        IA = zeros(2*tp, colors, num_in, class(oobv));
        TEI = zeros(2, 0, 'uint32');
        TE = zeros(4, 0, class(Kocc));
        V = true(2*tp, num_in);
        VA = true(2*tp, num_in);

        % For each input image...
        for a = 1:num_in
            % Calculate the coordinates in the input image
            T = WC * P(:,:,a);
            N = 1 ./ T(:,3);
            T(:,1) = T(:,1) .* N;
            T(:,2) = T(:,2) .* N;

            % Look up image samples
            IA(:,:,a) = squeeze(vgg_interp2(images{a}, T(:,1), T(:,2), 'linear', oobv));

            % Find interactions
            T(:,3) = T(:,3) ./ WC(:,4);
            N = occluding_pixels(T);
            M = abs(diff(N));
            M = M ~= tp; % & M ~= tp + 1 & M ~= 1 & M ~= sz(1) & M ~= tp + sz(1);
            N = uint32(N(:,M)); % Remove interactions between the same node

            % Mark all depths occluded by old depths
            V(N(2,N(1,:)<=tp),a) = false;
            if options.visibility
                VA(N(2,:),a) = false;
            end

            % Add the pixel interactions to the graph
            M = N(1,:) > tp;
            TEI = [TEI [N(1,:)-uint32(tp*M); (tp*2*a-tp)+N(2,:)]];
            T = zeros(4, numel(M));
            T(2,~M) = Kinf;
            T(4,M) = Kinf;
            TE = [TE T];
        end
        clear M T

        % Calculate the means for the new labelling
        A = ojw_bsxfun(@times, IA, permute(V, [1 3 2]));
        N = sum(V, 2);
        N(N==0) = 1;
        A = ojw_bsxfun(@times, sum(A, 3), 1./N);
        if options.show_output
            set(0, 'CurrentFigure', options.show_output); subplot('Position', [0 0.5 1/3 0.5]);
            sc(reshape(A(1:tp,:), [sz colors]), [0 255]);
            drawnow
        end

        % Calculate the data terms
        [U E EI TR TRI] = ibr_gen_cliques(IA, VA, V, Kocc, 1);
        clear IA VA V

        % Add surface smoothness constraints
        SE = [D(:)'; Dnew(:)'];
        SE = SE(:,SEI);
        % Fronto-parallel prior - Finite differences 1st derivative of
        % disparity
        SE = reshape(SE, 4, []);
        SE = squeeze(diff(reshape(SE([1 3; 1 4; 2 3; 2 4]',:), 2, 4, [])));
        % Apply our smoothness weighting
        if options.smoothness_kernel == 2
            SE = SE .^ 2; % Quadratic kernel
        else
            SE = abs(SE); % Linear kernel
        end
        SE = min(SE, options.disp_thresh);
        
        % Add texture constraints
        if options.tex_weight
            % Calculate texture edge costs
            A = permute(reshape(A, tp, 2, colors), [3 2 1]);
            A = mat2cell(double(A(:,:)), colors, repmat(2, [tp 1]));
            N = truncquad_edges(I, A, SEI, options.tex_thresh, options.tex_weight);
            N = cell2mat(reshape(N, 1, 1, []));
            N = reshape(permute(N, [2 1 3]), 4, []);
            SE = (1 + N) .* SE;
        end
        clear A N
        SE = cast(SE*options.lambda, class(E));

        n = size(EI, 2);
        if options.visibility
            % Concatenate data and visibility edges
            E = [E TE];
            EI = [EI TEI];
        end
        TE = TE(2,:) ~= 0;
        % Concatenate data and clique energies
        E = [E SE];
        EI = [EI SEI];
        % Fuse the two labellings, using contract and/or improve if desired
        [M stats] = vgg_qpbo(U, EI, E, TRI, TR, int32([tp options.improve options.contract options.contract>0]));
        out.ul(end+1) = stats(1);
        L = M == 1;
        E = E(:,1:n);
        EI = EI(:,1:n);

        % Update the depth map
        D(L) = Dnew(L);

        % Generate visibility maps
        TE = L(TEI(1,:)) ~= TE';
        V = true(2*tp, num_in);
        V(TEI(2,TE)-tp) = false;
        clear TEI TE
        T = (tp * L) + (1:tp)';
        L = [L; V(:)];
        for a = 1:num_in
            V(1:tp,a) = V(T);
            T = T + 2*tp;
        end
        V(tp+1:end,:) = [];

        % Calculate the energy
        U = U((0:tp-1)'*2+L(1:tp)+1);
        E = E((0:size(EI, 2)-1)'*4+L(EI(1,:))*2+L(EI(2,:))+1);
        TR = TR((0:size(TRI, 2)-1)'*8+L(TRI(1,:))*4+L(TRI(2,:))*2+L(TRI(3,:))+1);
        SE = SE((0:size(SEI, 2)-1)'*4+L(SEI(1,:))*2+L(SEI(2,:))+1);
        out.energy(end+1) = sum(U) + sum(E) + sum(TR) + sum(SE);

        if options.show_output
            % Display the output figures
            U = double(U) + accumarray(EI(1,:)', E, [tp 1]);
            U = U + accumarray(TRI(1,:)', TR, [tp 1]);
            U = reshape(U, sz);
            if options.visibility && num_in > 2
                % Take off the occlusion costs and normalize
                EI = reshape(sum(V, 2), sz);
                U = U - (num_in - EI) * double(Kocc + 1);
                E = EI ~= 0;
                U(E) = U(E) ./ EI(E);
            end
            set(0, 'CurrentFigure', options.show_output);
            subplot('Position', [1/3 0.5 1/3 0.5]);
            sc(D, 'contrast', disps([end 1]));
            subplot('Position', [2/3 0.5 1/3 0.5]);
            sc(U, 'jet');
            subplot('Position', [0 0 1/3 0.5]);
            T = reshape(sc(reshape(M, sz), 'prism'), [], 3);
            N = M < 0;
            T(N,:) = 1 - (1 - T(N,:)) * 0.3; % Lighten unlabelled pixels set to 0
            N = M > 1;
            T(N,:) = T(N,:) .* 0.3; % Darken unlabelled pixels set to 1 by optimal splice
            T(M==0,:) = 1;
            T(M==1,:) = 0;
            sc(reshape(T, [sz 3]), [0 1]);
            subplot('Position', [1/3 0 1/3 0.5]);
            sc(reshape(sum(V, 2), sz), [0 num_in], 'contrast');
            subplot('Position', [2/3 0 1/3 0.5]);
            U = accumarray(SEI_2nd, SE, [tp 1]);
            sc(reshape(-U, sz));
            drawnow;
        end
        clear U SE E EI TR TRI
        out.time_fuse(end+1) = toc;
    end
    
    if all(D==Dold)
        % Stop early if we've mad no progress
        break;
    end
end

% Generate the final image
% Calculate the homogenous coordinates of our two labellings
WC = WC(1:tp,:);
WC(:,4) = D(:);

% Initialise arrays for the image samples and occlusion terms
IA = zeros(tp, colors, num_in, class(oobv));
V = true(tp, num_in);

% For each input image...
for a = 1:num_in
    % Calculate the coordinates in the input image
    T = WC * P(:,:,a);
    N = 1 ./ T(:,3);
    T(:,1) = T(:,1) .* N;
    T(:,2) = T(:,2) .* N;

    % Look up image samples
    IA(:,:,a) = squeeze(vgg_interp2(images{a}, T(:,1), T(:,2), 'linear', oobv));

    % Find interactions
    T(:,3) = T(:,3) ./ WC(:,4);
    N = occluding_pixels(T);

    % Mark all depths occluded by old depths
    V(N(2,:),a) = false;
end

% Calculate the means for the new labelling
A = ojw_bsxfun(@times, IA, permute(V, [1 3 2]));
N = sum(V, 2);
N(N==0) = 1;
A = ojw_bsxfun(@times, sum(A, 3), 1./N);
A = reshape(A, [sz colors]);

if nargout > 1
    info.V = reshape(V, [sz num_in]);
    info.D = D;
end
return

