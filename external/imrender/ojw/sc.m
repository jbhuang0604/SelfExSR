function I = sc(I, varargin)
%SC  Display/output truecolor images with a range of colormaps
%
% Examples:
%   sc(image)
%   sc(..., limits)
%   sc(..., colormap)
%   out = sc(...)
%   sc
%
% Generates a truecolor RGB image based on the input values in 'image' and
% any maximum and minimum limits specified, using the colormap specified.
% The image is displayed on screen if there is no output argument.
% 
% SC has these advantages over MATLAB image rendering functions:
%   - images can be displayed or output; makes combining/overlaying images
%     simple.
%   - images are rendered/output in truecolor (RGB [0,1]); no nasty
%     discretization of the input data.
%   - many different truecolor colormaps for viewing various types of data.
%   - no border and automatic, integer magnification (unless figure is
%     docked or maximized) for better display.
%   - multiple images are displayed in a 3x4 grid with easy navigation
%     through arrow keys.
%
% For a demonstration, simply call SC without any input arguments.
%
% IN:
%   image - MxNxCxP or 3xMxNxP image array. MxN are the dimensions of the
%           image(s), C is the number of channels, and P the number of
%           images. If P > 1 images are diplayed in a 3x4 grid with arrow
%           keys for navigation.
%   limits - [min max] where values in image less than min will be set to
%            min and values greater than max will be set to max.
%   colormap - string containing the name of the color map to use to create
%              the output image. Default colormap: 'none', which is RGB for
%              3-channel images, grayscale otherwise. Conversion of multi-
%              channel images to intensity for intensity-based colormaps is
%              done using the L2 norm. Most Matlab colormaps are supported.
%              All colormaps can be reversed by prefixing '-' to the
%              string. This maintains integrity of the colorbar.
%              Special, non-Matlab colormaps are:
%      'contrast' - a high contrast colormap for intensity images that
%                   maintains intensity scale when converted to grayscale,
%                   for example when printing in black & white.
%      'prob' - first channel is plotted as hue, and the other channels
%               modulate intensity. Useful for laying probabilites over
%               images.
%      'prob_jet' - first channel is plotted as jet colormap, and the other
%                   channels modulate intensity.
%      'diff' - intensity values are marked blue for > 0 and red for < 0.
%               Darker colour means larger absolute value. For multi-
%               channel images, the L2 norm of the other channels sets
%               green level. 3 channel images are converted to YUV and
%               images with more that 3 channels are projected onto the
%               principle components first.
%      'compress' - compress many channels to RGB while maximizing
%                   variance.
%      'flow' - display two channels representing a 2d Cartesian vector as
%               hue for angle and intensity for magnitude (darker colour
%               indicates a larger magnitude).
%      'phase' - first channel is intensity, second channel is phase in
%                radians. Darker colour means greater intensity, hue
%                represents phase from 0 to 2 pi.
%      'stereo' - pair of concatenated images used to generate a red/cyan
%                 anaglyph.
%      'stereo_col' - pair of concatenated RGB images used to generate a
%                     colour anaglyph.
%      'rand' - gives an index image a random colormap. Useful for viewing
%               segmentations.
%      'rgb2gray' - converts an RGB image to grayscale in the same fashion
%                   as Matlab's rgb2gray (in the image processing toolbox). 
%
% OUT:
%   out - MxNx3xP truecolour (double) RGB image array in range [0, 1]
%
% See also IMAGE, IMAGESC, IMSHOW, COLORMAP, COLORBAR.

% $Id: sc.m,v 1.8 2008/06/04 10:50:55 ojw Exp $
% Copyright: Oliver Woodford, 2007

%% Check for arguments
if nargin == 0
    % If there are no input arguments then run the demo
    if nargout > 0
        error('Output expected from no inputs!');
    end
    demo; % Run the demo
    return
end

%% Size our image(s)
I = I(:,:,:,:);
[y x c n] = size(I);

%% Don't do much if I is empty
if isempty(I)
    if nargout == 0
        % Clear the current axes if we were supposed to display the image
        cla; axis off;
    else
        % Create an empty array with the correct dimensions
        I = zeros(y, x, (c~=0)*3, n);
    end
    return
end

%% Parse the input arguments coming after I (1st input)
map = 'none';
limits = [];
for i = 1:nargin-1
    if ischar(varargin{i})
        map = varargin{i};
    elseif isnumeric(varargin{i}) && length(varargin{i}) > 1
        limits = double(varargin{i}(1:2));
    end
end

%% Check if image is given with RGB colour along the first dimension
if y == 3 && c > 3
    % Flip colour to 3rd dimension
    I = permute(I, [2 3 1 4]);
    [y x c n] = size(I);
end

%% Check for multiple images
% If we have a non-singleton 4th dimension we want to display the images in
% a 3x4 grid and use buttons to cycle through them
global SC_INDEX
if n > 1
    if nargout > 0
        % Return transformed images in an YxXx3xN array
        A = zeros([y x 3 n]);
        for a = 1:n
            A(:,:,:,a) = sc(I(:,:,:,a), map, limits);
        end
        I = A;
    else
        % Create a data structure to store the data in
        state.index = 13;
        state.map = map;
        state.n = 12 * ceil(n / 12);
        state.limits = limits;
        state.I = squeeze(num2cell(I, [1 2 3]));
        state.I = [state.I; cell(state.n-n, 1)];
        fig = gcf;
        clear I % Stop the matrix getting printed
        % Call the callback to initialise the figure
        keypress_callback(fig, struct('Character', 28), state);
        % Check if we need to be able to scroll through images
        if n > 12
            % Set the callback for image navigation, and save the image
            % data in the figure
            set(fig, 'KeyPressFcn', @keypress_callback, 'Interruptible', 'off');
        end
    end
    return
end

%% If map starts with a '-' sign, invert the colourmap
if map(1) == '-'
    map(1) = [];
    reverseMap = true;
else
    reverseMap = false;
end

%% Call the rendering function
I = reshape(double(real(I)), [], c); % Only work with real doubles

% Large switch statement for all the colourmaps
switch map
%% Prism
    case 'prism'
        % Similar to the Matlab internal prism colormap, but only works on
        % index images, assigning each index (or rounded float) to a
        % different colour
        [I limits] = index_im(I);
        % Generate prism colourmap
        map = prism(6);
        if reverseMap
            map = map(end:-1:1,:); % Reverse the map
        end
        % Lookup the colours
        I =  mod(I, 6) + 1;
        I = map(I,:);
%% Rand
    case 'rand'
        % Assigns a random colour to each index
        [I limits num_vals] = index_im(I);
        % Generate random colourmap
        map = rand(num_vals, 3);
        % Lookup the colours
        I = map(I,:);
%% Diff
    case 'diff'
        % Show positive as blue and negative as red, white is 0
        switch c
            case 1
                I(:,2:3) = 0;
            case 2
                % Second channel can only have absolute value
                I(:,3) = abs(I(:,2));
            case 3
                % Diff of RGB images - convert to YUV first
                I = rgb2yuv(I);
                I(:,3) = sqrt(sum(I(:,2:end) .^ 2, 2)) ./ sqrt(2);
            otherwise
                % Use difference along principle component, and other
                % channels to modulate second channel
                I = calc_prin_comps(I);
                I(:,3) = sqrt(sum(I(:,2:end) .^ 2, 2)) ./ sqrt(c - 1);
                I(:,4:end) = [];
        end
        % Generate limits
        if isempty(limits)
            limits = [min(I(:,1)) max(I(:,1))];
        end
        limits = max(abs(limits));
        % Scale
        if c > 1
            I(:,[1 3]) = I(:,[1 3]) / limits;
        else
            I = I / (limits * 0.5);
        end
        % Colour
        M = I(:,1) > 0;
        I(:,2) = -I(:,1) .* ~M;
        I(:,1) = I(:,1) .* M;
        if reverseMap
            % Swap first two channels
            I = I(:,[2 1 3]);
        end
        %I = 1 - I * [1 0.4 1; 0.4 1 1; 1 1 0.4]; % (Green/Red)
        I = 1 - I * [1 1 0.4; 0.4 1 1; 1 0.4 1]; % (Blue/Red)
        I = min(max(I(:), 0), 1);
        limits = [-limits limits]; % For colourbar
%% Flow
    case 'flow'
        % Calculate amplitude and phase, and use 'phase'
        if c ~= 2
            error('''flow'' requires two channels');
        end
        A = sqrt(sum(I .^ 2, 2));
        if isempty(limits)
            limits = [min(A) max(A)*2];
        else
            limits = [0 max(abs(limits)*sqrt(2))*2];
        end
        I(:,1) = atan2(I(:,2), I(:,1));
        I(:,2) = A;
        if reverseMap
            % Invert the amplitude
            I(:,2) = -I(:,2);
            limits = -limits([2 1]);
        end
        I = phase_helper(I, limits, 2); % Last parameter tunes how saturated colors can get
        % Set NaNs (unknown flow) to 0
        I(isnan(I)) = reverseMap;
        limits = []; % This colourmap doesn't have a valid colourbar
%% Phase
    case 'phase'
        % Plot amplitude as intensity and angle as hue
        if c < 2
            error('''phase'' requires two channels');
        end
        if isempty(limits)
            limits = [min(I(:,1)) max(I(:,1))];
        end
        if reverseMap
            % Invert the phase 
            I(:,2) = -I(:,2);
        end
        I = I(:,[2 1]);
        I = phase_helper(I, limits, 1.3); % Last parameter tunes how saturated colors can get
        limits = []; % This colourmap doesn't have a valid colourbar
%% RGB2Grey
    case {'rgb2grey', 'rgb2gray'}
        % Compress RGB to greyscale
        [I limits] = rgb2grey(I, limits, reverseMap);
%% RGB2YUV
    case 'rgb2yuv'
        % Convert RGB to YUV - not for displaying or saving to disk!
        [I limits] = rgb2yuv(I);
%% YUV2RGB
    case 'yuv2rgb'
        % Convert YUV to RGB - undo conversion of rgb2yuv
        if c ~= 3
            error('''yuv2rgb'' requires a 3 channel image');
        end
        I = reshape(I, [], 3);
        I = I * [1 1 1; 0, -0.39465, 2.03211; 1.13983, -0.58060  0];
        I = reshape(I, y, x, 3);
        I = sc(I, limits);
        limits = []; % This colourmap doesn't have a valid colourbar
%% Prob
    case 'prob'
        % Plot first channel as grey variation of 'bled' and modulate
        % according to other channels
        if c > 1
            A = rgb2grey(I(:,2:end), [], false);
            I = I(:,1);
        else
            A = 0.5;
        end
        [I limits] = bled(I, limits, reverseMap);
        I = normalize(A + I, [-0.1 1.3]);
%% Prob_jet
    case 'prob_jet'
        % Plot first channel as 'jet' and modulate according to other
        % channels
        if c > 1
            A = rgb2grey(I(:,2:end), [], false);
            I = I(:,1);
        else
            A = 0.5;
        end
        [I limits] = jet_helper(I, limits, reverseMap);
        I = normalize(A + I, [0.2 1.8]);
%% Compress
    case 'compress'
        % Compress to RGB, maximizing variance
        I = reshape(I, [], c);
        % Normalize
        if isempty(limits)
            limits = [min(I(:)) max(I(:))];
            I = I - limits(1);
            if limits(2) ~= limits(1)
                I = I * (1 / (limits(2) - limits(1)));
            end
        else
            I = I - limits(1);
            if limits(2) ~= limits(1)
                I = I * (1 / (limits(2) - limits(1)));
            end
            I = reshape(min(max(I(:), 0), 1), [], c);
        end
        if reverseMap
            % Invert after truncation
            I = 1 - I;
        end
        % Zero mean
        meanCol = sum(I, 1) / (x * y);
        isBsx = exist('bsxfun', 'builtin');
        if isBsx
            I = bsxfun(@minus, I, meanCol);
        else
            I = I - meanCol(ones(x*y, 1, 'uint8'),:);
        end
        % Calculate top 3 principle components
        I = calc_prin_comps(I, 3);
        % Normalize each channel independently
        if isBsx
            I = bsxfun(@minus, I, min(I, [], 1));
            I = bsxfun(@times, I, 1./max(I, [], 1));
        else
            for a = 1:3
                I(:,a) = I(:,a) - min(I(:,a));
                I(:,a) = I(:,a) / max(I(:,a));
            end
        end
        % Put components in order of human eyes' response to channels
        I = I(:,[2 1 3]);
        limits = []; % This colourmap doesn't have a valid colourbar
%% Stereo (anaglyph)     
    case 'stereo'
        % Convert 2 colour images to intensity images
        % Show first channel as red and second channel as cyan
        A = rgb2grey(I(:,1:floor(end/2)), limits, false);
        I = rgb2grey(I(:,floor(end/2)+1:end), limits, false);
        if reverseMap
            I(:,2:3) = A(:,1:2); % Make first image cyan
        else
            I(:,1) = A(:,1); % Make first image red
        end
        limits = []; % This colourmap doesn't have a valid colourbar
%% Coloured anaglyph      
    case 'stereo_col'
        if c ~= 6
            error('''stereo_col'' requires a 6 channel image');
        end
        I = normalize(I, limits);
        % Red channel from one image, green and blue from the other
        if reverseMap
            I(:,1) = I(:,4); % Make second image red
        else
            I(:,2:3) = I(:,5:6); % Make first image red
        end
        I = I(:,1:3);
        limits = []; % This colourmap doesn't have a valid colourbar
%% None            
    case 'none'
        % No colour map - just output the image
        if c ~= 3
            [I limits] = grey(I, limits, reverseMap);
        else
            I = intensity(I(:), limits, reverseMap);
            limits = [];
        end
%% Grey                
    case {'gray', 'grey'}
        % Greyscale
        [I limits] = grey(I, limits, reverseMap);
%% Jet                
    case 'jet'
        % Dark blue to dark red, through green
        [I limits] = jet_helper(I, limits, reverseMap);
%% Hot                
    case 'hot'
        % Black to white, through red
        [I limits] = intensity(I, limits, reverseMap); % Intensity map
        I = I * 3;
        I = [I, I-1, I-2];
        I = min(max(I(:), 0), 1); % Truncate
%% Contrast                
    case 'contrast'
        % A high contrast, full-colour map that goes from black to white
        % linearly when converted to greyscale
        [I limits] = intensity(I, limits, reverseMap); % Intensity map
        % Calculate red and blue channels
        I = [abs(sin(I*(1.5*pi))) I sin(I*(3.5*pi)).^2];
        % Calculate green channel
        I(:,2) = I * ([-0.299; 1; -0.114] / 0.587);
        % Update red and blue channels to avoid green being saturated
        M0 = I(:,2) < 0;
        M1 = I(:,2) > 1;
        M = M0 | M1;
        K = I(:,2);
        K = K - M1;
        J = I(M,[1 3]);
        K = K(M) ./ (J * ([0.299; 0.114] / 0.587)) + 1;
        I(M,[1 3]) = J .* K(:,[1 1]);
        % Quick truncation using cached indices
        I(M0,2) = 0;
        I(M1,2) = 1;
%% HSV                
    case 'hsv'
        % Cycle through hues
        [I limits] = intensity(I, limits, reverseMap); % Intensity map
        I = hsv_helper(I);
%% Bone                
    case 'bone'
        % Greyscale with a blue tint (a tiny bit different to Matlab's)
        [I limits] = intensity(I, limits, reverseMap); % Intensity map
        J = [I-2/3, I-1/3, I];
        J = reshape(max(min(J(:), 1/3), 0), [y*x 3]) * (2 / 5);
        I = I * (13 / 15);
        I = J + I(:,[1 1 1]);
%% Colourcube                
    case {'colorcube', 'colourcube'}
        % Psychedelic colourmap inspired by Matlab's version
        [I limits] = intensity(I, limits, reverseMap); % Intensity map
        step = 4;
        I = I * (step * (1 - eps));
        J = I * step;
        K = floor(J);
        I = cat(3, mod(K, step)/(step-1), J - floor(K), mod(floor(I), step)/(step-1));
%% Cool                
    case 'cool'
        % Cyan through to magenta
        [I limits] = intensity(I, limits, reverseMap); % Intensity map
        I = [I, 1-I, ones(y*x, 1)];
%% Spring                
    case 'spring'
        % Magenta through to yellow
        [I limits] = intensity(I, limits, reverseMap); % Intensity map
        I = [ones(y*x, 1), I, 1-I];
%% Summer                
    case 'summer'
        % Darkish green through to pale yellow
        [I limits] = intensity(I, limits, reverseMap); % Intensity map
        I = [I, 0.5+I*0.5, 0.4*ones(y*x, 1)];
%% Autumn                
    case 'autumn'
        % Red through to yellow
        [I limits] = intensity(I, limits, reverseMap); % Intensity map
        I = [ones(y*x, 1), I, zeros(y*x, 1)];
%% Winter                
    case 'winter'
        % Blue through to turquoise
        [I limits] = intensity(I, limits, reverseMap); % Intensity map
        I = [zeros(y*x, 1), I, 1-I*0.5];
%% Copper                
    case 'copper'
        % Black through to copper
        [I limits] = intensity(I, limits, reverseMap); % Intensity map
        I = [I*(1/0.8), I*0.78, I*0.5];
        I = min(max(I(:), 0), 1); % Truncate
%% Pink                
    case 'pink'
        % Greyscale with a pink tint
        [I limits] = intensity(I, limits, reverseMap); % Intensity map
        J = I * (2 / 3);
        I = [I, I-1/3, I-2/3];
        I = reshape(max(min(I(:), 1/3), 0), [y*x 3]);
        I = I + J(:,[1 1 1]);
        I = sqrt(I);
%% Bled
    case 'bled'
        % Black to red, through blue
        [I limits] = bled(I, limits, reverseMap);
%% Unknown colourmap
    otherwise
        error('Colormap ''%s'' not recognised.', map);
end % End of colourmap switch statement
I = reshape(I, [y x 3]); % Reshape to correct size

%% Only display if the output isn't used
if nargout == 0
    % Get handles to the figure and axes
    hFig = gcf; hAx = gca;
    % Clear what's there
    cla(hAx, 'reset');
    % Display the image - using image() is fast
    hIm = image(I);
    % Don't print out the matrix if we've forgotten the ";"
    clear I
    % Axes invisible and equal
    set(hFig, 'Units', 'pixels');
    bgcolor = get(hFig, 'Color');
    set(hAx, 'Visible', 'on', 'Box', 'off', 'Xtick', [], 'Ytick', [], 'Color', 'none', 'XColor', bgcolor, 'YColor', bgcolor, 'DataAspectRatio', [1 1 1], 'DrawMode', 'fast', 'Layer', 'bottom', 'LineWidth', 0.001);
    % Set data for a colourbar
    if ~isempty(limits) && limits(1) ~= limits(2)
        colBar = (0:255) * ((limits(2) - limits(1)) / 255) + limits(1);
        colBar = squeeze(sc(colBar, map, limits));
        if reverseMap
            colBar = colBar(end:-1:1,:);
        end
        set(hFig, 'Colormap', colBar);
        set(hAx, 'CLim', limits);
        set(hIm, 'CDataMapping', 'scaled');
    end
    % Only resize image if it is alone in the figure
    if numel(findobj(get(hFig, 'Children'), 'Type', 'axes')) > 1
        return
    end
    % Could still be the first subplot - do another check
    axesPos = get(hAx, 'Position');
    if isequal(axesPos, get(hFig, 'DefaultAxesPosition'))
        % Default position => not a subplot
        % Fill the window
        set(hAx, 'Units', 'normalized', 'Position', [0 0 1 1]);
        axesPos = [0 0 1 1];
    end
    if ~isequal(axesPos, [0 0 1 1]) || strcmp(get(hFig, 'WindowStyle'), 'docked')
        % Figure not alone, or docked. Either way, don't resize.
        return
    end
    % Get the size of the monitor we're on
    figPosCur = get(hFig, 'Position');
    MonSz = get(0, 'MonitorPositions');
    MonOn = size(MonSz, 1);
    if MonOn > 1
        figCenter = figPosCur(1:2) + figPosCur(3:4) / 2;
        figCenter = MonSz - repmat(figCenter, [MonOn 2]);
        MonOn = all(sign(figCenter) == repmat([-1 -1 1 1], [MonOn 1]), 2);
        MonOn(1) = MonOn(1) | ~any(MonOn);
        MonOn = find(MonOn);
        if MonOn ~= 1
            
        end
        MonSz = MonSz(MonOn,:);
    end
    MonSz(3:4) = MonSz(3:4) - MonSz(1:2) + 1;
    % Check if the window is maximized
    % This is a hack which may only work on Windows! No matter, though.
    if isequal(MonSz([1 3]), figPosCur([1 3]))
        % Leave maximized
        return
    end
    % Compute the size to set the window
    MaxSz = MonSz(3:4) - [20 120];
    RescaleFactor = min(MaxSz ./ [x y]);
    if RescaleFactor > 1
        % Integer scale for enlarging, but don't make too big
        MaxSz = min(MaxSz, [1000 680]);
        RescaleFactor = max(floor(min(MaxSz ./ [x y])), 1);
    end
    figPosNew = ceil([x y] * RescaleFactor);
    % Don't move the figure if the size needs to change
    if isequal(figPosCur(3:4), figPosNew)
        return
    end
    % Keep the centre of the figure stationary
    figPosNew = [max(1, floor(figPosCur(1:2)+(figPosCur(3:4)-figPosNew)/2)) figPosNew];
    % Ensure the figure bar is in bounds
    figPosNew(1:2) = min(figPosNew(1:2), MonSz(1:2)+MonSz(3:4)-[6 101]-figPosNew(3:4));
    set(hFig, 'Position', figPosNew);
end
return

%% Keypress callback
% For displaying multiple (> 12) images using arrow keys for navigation
function keypress_callback(fig, event_data, state)
% Get the state data, if not given
if nargin < 3
    state = get(fig, 'UserData');
end
% Get the current index
index = state.index;
% Check what key was pressed and update the image index as necessary
switch event_data.Character
    case 28 % Left
        index = index - 13;
    case 29 % Right
        index = index + 11;
    case 30 % Up
        index = index + 119;
    case 31 % Down
        index = index - 121;
    otherwise
        % Another key was pressed - ignore it
        return
end
% Generate 12 valid indices
index = mod(index:index+11, state.n) + 1;
% Plot the images
figure(fig);
for a = 0:2
    for b = 0:3
        subplot('Position', [b*0.25062 a*0.334437 0.248139 0.331126]);
        sc(state.I{index(1+b+(2-a)*4)}, state.map, state.limits);
    end
end
drawnow;
% Save the current index
state.index = index(1);
set(fig, 'UserData', state);
return

%% Grey
function [I limits] = grey(I, limits, reverseMap)
% Greyscale
[I limits] = intensity(I, limits, reverseMap);
I = I(:,[1 1 1]);
return

%% RGB2grey
function [I limits] = rgb2grey(I, limits, reverseMap)
% Compress RGB to greyscale
if size(I, 2) == 3
    I = I * [0.299; 0.587; 0.114];
end
[I limits] = grey(I, limits, reverseMap);
return
        
%% RGB2YUV
function [I limits] = rgb2yuv(I)
% Convert RGB to YUV - not for displaying or saving to disk!
if size(I, 2) ~= 3
    error('rgb2yuv requires a 3 channel image');
end
I = I * [0.299, -0.14713, 0.615; 0.587, -0.28886, -0.51498; 0.114, 0.436, -0.10001];
limits = []; % This colourmap doesn't have a valid colourbar
return

%% Phase helper
function I = phase_helper(I, limits, n)
I(:,1) = mod(I(:,1)/(2*pi), 1);
I(:,2) = I(:,2) - limits(1);
I(:,2) = I(:,2) * (n / (limits(2) - limits(1)));
I(:,3) = n - I(:,2);
I(:,[2 3]) = min(max(I(:,[2 3]), 0), 1);
I = hsv2rgb(reshape(I, [], 1, 3));
return

%% Jet helper        
function [I limits] = jet_helper(I, limits, reverseMap)
% Dark blue to dark red, through green
[I limits] = intensity(I, limits, reverseMap);
I = I * 4;
I = [I-3, I-2, I-1];
I = 1.5 - abs(I);
I = reshape(min(max(I(:), 0), 1), size(I));
return

%% HSV helper
function I = hsv_helper(I)
I = I * 6;
I = abs([I-3, I-2, I-4]);
I(:,1) = I(:,1) - 1;
I(:,2:3) = 2 - I(:,2:3);
I = reshape(min(max(I(:), 0), 1), size(I));
return

%% Bled
function [I limits] = bled(I, limits, reverseMap)
% Black to red through blue
[I limits] = intensity(I, limits, reverseMap);
J = reshape(hsv_helper(I), [], 3);
if exist('bsxfun', 'builtin') 
    I = bsxfun(@times, I, J);
else
    I = J .* I(:,[1 1 1]);
end
return

%% Normalize
function [I limits] = normalize(I, limits)
if isempty(limits)
    limits = [min(I(:)) max(I(:))];
    I = I - limits(1);
    if limits(2) ~= limits(1)
        I = I * (1 / (limits(2) - limits(1)));
    end
else
    I = I - limits(1);
    if limits(2) ~= limits(1)
        I = I * (1 / (limits(2) - limits(1)));
    end
    I = reshape(min(max(I(:), 0), 1), size(I));
end
return

%% Intensity maps
function [I limits] = intensity(I, limits, reverseMap)
% Squash to 1d using L2 norm
if size(I, 2) > 1
    I = sqrt(sum(I .^ 2, 2));
end
% Determine and scale to limits
[I limits] = normalize(I, limits);
if reverseMap
    % Invert after everything
    I = 1 - I;
end
return

%% Index images
function [J limits num_vals] = index_im(I)
% Returns an index image
if size(I, 2) ~= 1
    error('Index maps only work on single channel images');
end
J = round(I);
rescaled = any(abs(I - J) > 0.01);
if rescaled
    % Appears not to be an index image. Rescale over 256 indices
    m = min(I);
    m = m * (1 - sign(m) * eps);
    I = I - m;
    I = I * (256 / max(I(:)));
    J = ceil(I);
    num_vals = 256;
elseif nargout > 2
    % Output the number of values
    J = J - (min(J) - 1);
    num_vals = max(J);
end
% These colourmaps don't have valid colourbars
limits = [];    
return

%% Calculate principle components
function I = calc_prin_comps(I, numComps)
if nargin < 2
    numComps = size(I, 2);
end
% Do SVD
[I S] = svd(I, 0);
% Calculate projection of data onto components
S = diag(S(1:numComps,1:numComps))';
if exist('bsxfun', 'builtin')
    I = bsxfun(@times, I(:,1:numComps), S);
else
    I = I(:,1:numComps) .* S(ones(size(I, 1), 1, 'uint8'),:);
end
return

%% Demo function to show capabilities of sc
function demo
%% Demo gray & lack of border
figure; fig = gcf; Z = peaks(256); sc(Z);
display_text([...
' Lets take a standard, Matlab, real-valued function:\n\n    peaks(256)\n\n'...
' Calling:\n\n    figure\n    Z = peaks(256);\n    sc(Z)\n\n'...
' gives (see figure). SC automatically scales intensity to fill the\n'...
' truecolor range of [0 1].\n\n'...
' If your figure isn''t docked, then the image will have no border, and\n'...
' will be magnified by an integer factor (in this case, 2) so that the\n'...
' image is a reasonable size.']);

%% Demo colour image display 
figure(fig); clf;
load mandrill; mandrill = ind2rgb(X, map); sc(mandrill);
display_text([...
' That wasn''t so interesting. The default colormap is ''none'', which\n'...
' produces RGB images given a 3-channel input image, otherwise it produces\n'...
' a grayscale image. So calling:\n\n    load mandrill\n'...
'    mandrill = ind2rgb(X, map);\n    sc(mandrill)\n\n gives (see figure).\n']);

%% Demo discretization
figure(fig); clf;
subplot(121); sc(Z, 'jet'); label(Z, 'sc(Z, ''jet'')');
subplot(122); imagesc(Z); axis image off; colormap(jet(64)); % Fix the fact we change the default depth
label(Z, 'imagesc(Z); axis image off; colormap(''jet'');');
display_text([...
' However, if we want to display intensity images in color we can use any\n'...
' of the Matlab colormaps implemented (most of them) to give truecolor\n'...
' images. For example, to use ''jet'' simply call:\n\n'...
'    sc(Z, ''jet'')\n\n'...
' The Matlab alternative, shown on the right, is:\n\n'...
'    imagesc(Z)\n    axis equal off\n    colormap(jet)\n\n'...
' which generates noticeable discretization artifacts.']);

%% Demo intensity colourmaps
figure(fig); clf;
subplot(221); sc(Z, 'hsv'); label(Z, 'sc(Z, ''hsv'')');
subplot(222); sc(Z, 'colorcube'); label(Z, 'sc(Z, ''colorcube'')');
subplot(223); sc(Z, 'contrast'); label(Z, 'sc(Z, ''contrast'')');
subplot(224); sc(Z-round(Z), 'diff'); label(Z, 'sc(Z-round(Z), ''diff'')');
display_text([...
' There are several other intensity colormaps to choose from. Calling:\n\n'...
'    help sc\n\n'...
' will give you a list of them. Here are several others demonstrated.\n']);

%% Demo saturation limits & colourmap reversal
figure(fig); clf;
subplot(121); sc(Z, [0 max(Z(:))], '-hot'); label(Z, 'sc(Z, [0 max(Z(:))], ''-hot'')');
subplot(122); sc(mandrill, [-0.5 0.5]); label(mandrill, 'sc(mandrill, [-0.5 0.5])');
display_text([...
' SC can also rescale intensity, given an upper and lower bound provided\n'...
' by the user, and invert most colormaps simply by prefixing a ''-'' to the\n'...
' colormap name. For example:\n\n'...
'    sc(Z, [0 max(Z(:))], ''-hot'');\n'...
'    sc(mandrill, [-0.5 0.5]);\n\n'...
' Note that the order of the colormap and limit arguments are\n'...
' interchangable.\n']);

%% Demo prob
load gatlin;
gatlin = X;
figure(fig); clf; im = cat(3, abs(Z)', gatlin(1:256,end-255:end)); sc(im, 'prob');
label(im, 'sc(cat(3, prob, gatlin), ''prob'')');
display_text([...
' SC outputs the recolored data as a truecolor RGB image. This makes it\n'...
' easy to combine colormaps, either arithmetically, or by masking regions.\n'...
' For example, we could combine an image and a probability map\n'...
' arithmetically as follows:\n\n'...
'    load gatlin\n'...
'    gatlin = X(1:256,end-255:end);\n'...
'    prob = abs(Z)'';\n'...
'    im = sc(prob, ''hsv'') .* sc(prob, ''gray'') + sc(gatlin, ''rgb2gray'');\n'...
'    sc(im, [-0.1 1.3]);\n\n'...
' In fact, that particular colormap has already been implemented in SC.\n'...
' Simply call:\n\n'...
'    sc(cat(3, prob, gatlin), ''prob'');\n']);

%% Demo colorbar
colorbar;
display_text([...
' SC also makes possible the generation of a colorbar in the normal way, \n'...
' with all the colours and data values correct. Simply call:\n\n'...
'    colorbar\n\n'...
' The colorbar doesn''t work with all colormaps, but when it does,\n'...
' inverting the colormap (using ''-map'') maintains the integrity of the\n'...
' colorbar (i.e. it works correctly) - unlike if you invert the input data.\n']);

%% Demo combine by masking
figure(fig); clf;
im = cat(4, sc(Z, [0 max(Z(:))], '-hot'), sc(Z-round(Z), 'diff'));
mask = repmat(Z > 0, [1 1 3]);
mask = cat(4, mask, ~mask);
im = sum(im .* mask, 4);
sc(im);
display_text([...
' It''s just as easy to combine generated images by masking too. Here''s an\n'...
' example:\n\n'...
'    im = cat(4, sc(Z, [0 max(Z(:))], ''-hot''),...\n'...
'                sc(Z-round(Z), ''diff''));\n'...
'    mask = repmat(Z < 0, [1 1 3]);\n'...
'    mask = cat(4, mask, ~mask);\n'...
'    im = sum(im .* mask, 4);\n'...
'    sc(im)\n\n'...
' Of course, you can also do any combination of the two methods.\n']);

%% Demo texture map
figure(fig); clf;
surf(Z, sc(Z, 'contrast'), 'edgecolor', 'none');
display_text([...
' Other benefits of SC outputting the image as an array are that the image\n'...
' can be saved straight to disk using imwrite() (if you have the image\n'...
' processing toolbox), or can be used to texture map a surface, thus:\n\n'...
'    tex = sc(Z, ''contrast'');\n'...
'    surf(Z, tex, ''edgecolor'', ''none'');\n\n']);

%% Demo compress
load mri;
mri = D;
close(fig); % Only way to get round loss of focus (bug?)
figure(fig); clf;
sc(squeeze(mri(:,:,:,1:6)), 'compress');
display_text([...
' For images with more than 3 channels, SC can compress these images to RGB\n'...
' while maintaining the maximum amount of variance in the data. For\n'...
' example, this 6 channel image:\n\n'...
'    load mri\n    mri = D;\n    sc(squeeze(mri(:,:,:,1:6), ''compress'')\n\n']);

%% Demo multiple images
figure(fig); clf; sc(mri, 'bone');
display_text([...
' Finally, SC makes it easy for you to view multiple images in a grid when\n'...
' passed in as a 4d array. For example:\n\n'...
'    sc(mri, ''bone'')\n\n'...
' Use the arrow keys to navigate through the images.']);

clc; fprintf('End of demo.\n');
return

%% Some helper functions for the demo
function display_text(str)
clc;
fprintf([str '\n\n']);
fprintf('Press a key to go on.\n');
figure(gcf);
waitforbuttonpress;
return

function label(im, str)
text(size(im, 2)/2, size(im, 1)+12, str,...
    'Interpreter', 'none', 'HorizontalAlignment', 'center', 'VerticalAlignment', 'middle');
return
