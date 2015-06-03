function download_stereo(sequence)
%DOWNLOAD_STEREO  Downloads a stereo sequence from the Middlebury website
%
%   download_stereo(sequence)
%
% Downloads the named Middlebury stereo sequence from the Middlebury stereo
% website, into a subdirectory of the sme name, if the subdirectory doesn't
% already exist.
%
%IN:
%   sequence - string containing the sequence name, e.g. 'cones'.

% $Id: download_stereo.m,v 1.2 2008/11/17 11:27:35 ojw Exp $

% Check if the dataset exists
subdir = ['./' sequence '/'];
if isdir(subdir)
    return
end

% Load the sequence specific details
switch sequence
    case 'tsukuba'
        stereo_class = 2001; % Middlebury 2001 dataset
        inn = 5;
        disp_range = [0 15];
        disparity_factor = 16;
        imname = @(n) sprintf('scene1.row3.col%d.ppm', n);
        disp_ims = [3 4];
        dispname = @(n) sprintf('truedisp.row3.col%d.pgm', disp_ims(n));
    case 'map'
        stereo_class = 2001; % Middlebury 2001 dataset
        inn = 2;
        disp_range = [0 31];
        disparity_factor = 8;
        imname = @(n) sprintf('im%d.pgm', n-1);
        disp_ims = [1 2];
        dispname = @(n) sprintf('disp%d.pgm', disp_ims(n)-1);
    case {'sawtooth','venus','bull','poster','barn1','barn2'}
        stereo_class = 2001; % Middlebury 2001 dataset
        inn = 9;
        disp_range = [0 31];
        disparity_factor = 8;
        imname = @(n) sprintf('im%d.ppm', n-1);
        disp_ims = [3 7];
        dispname = @(n) sprintf('disp%d.pgm', disp_ims(n)-1);
    case {'teddy','cones'}
        stereo_class = 2003; % Middlebury 2003 dataset
        inn = 9;
        imname = @(n) sprintf('im%d.ppm', n-1);
        disp_range = [0 59];
        disparity_factor = 4;
        disp_ims = [3 7];
        dispname = @(n) sprintf('disp%d.pgm', disp_ims(n)-1);
    case {'dolls','art','books','moebius','laundry','reindeer','dwarves','drumsticks','computer'}
        stereo_class = 2005; % Middlebury 2005 dataset
    case {'aloe','baby1','baby2','baby3','bowling1','bowling2','cloth1','cloth2','cloth3','cloth4',...
          'flowerpots','lampshade1','lampshade2','midd1','midd2','monopoly','plastic','rocks1',...
          'rocks2','wood1','wood2'}
        stereo_class = 2006; % Middlebury 2006 dataset
    otherwise
        error('Sequence ''%s'' not recognised', sequence);
end

% Values for 2005 & 2006 datasets are the same.
if stereo_class >= 2005
    inn = 7;
    disp_range = [0 85];
    disparity_factor = 3;
    imname = @(n) sprintf('view%d.png', n-1);
    disp_ims = [2 6];
    dispname = @(n) sprintf('disp%d.png', disp_ims(n)-1);
    sequence = [upper(sequence(1)) sequence(2:end)];
end
if ~exist('im_space', 'var')
    im_space = diff(disp_ims);
end

% Create the directory
mkdir(subdir);

% Download the sequence from the Middlebury website
base = sprintf('http://vision.middlebury.edu/stereo/data/scenes%d/', stereo_class);
if stereo_class == 2001
    base = [base 'data/' sequence '/'];
elseif stereo_class == 2003
    base = [base 'newdata/test/' sequence '/input/'];
else
    base = [base 'ThirdSize/' sequence '/'];
end
try
    % Download disparity images if available
    urlread([base dispname(1)]);
    urlwrite([base dispname(1)], [subdir dispname(1)]);
    urlread([base dispname(2)]);
    urlwrite([base dispname(2)], [subdir dispname(2)]);
catch
end
switch stereo_class
    case 2005
        base = [base 'Illum1/Exp1/'];
    case 2006
        base = [base 'Illum1/Exp2/'];
end
for a = 1:inn
    % Download images
    im = imname(a);
    urlwrite([base im], [subdir im]);
    % Convert to correct form
    imwrite(imread([subdir im]), sprintf('%sinput.%3.3d.png', subdir, a-1));
    delete([subdir im]);
end

% Generate the data file
% Input projection matrices
Pi = repmat([eye(3) zeros(3, 1)], [1 1 inn]);
Pi(1,4,:) = -(0:inn-1) ./ (disparity_factor * im_space);

% Disparity values
disps = disp_range(1)*disparity_factor:disp_range(2)*disparity_factor;

% Save
save([subdir 'data.mat'], '-v6', 'Pi', 'disps');