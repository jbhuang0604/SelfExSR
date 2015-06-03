%VGG_SEGMENT_MS  Mean shift image segmentation
%
%   S = vgg_segment_ms(A, h_s, h_r, min_sz[, W])
%
% Segmentation of an image using mean shift.
% 
% Uses EDISON code by Chris Christoudias and Bogdan Georgescu, downloaded
% from http://www.caip.rutgers.edu/riul/research/code/EDISON/index.html,
% and based on the following papers:
%   [1] D. Comanicu, P. Meer: "Mean shift: A robust approach toward feature
%   space analysis". IEEE Trans. Pattern Anal. Machine Intell., May 2002.
%   [2] C. Christoudias, B. Georgescu, P. Meer: "Synergism in low level
%   vision". 16th International Conference of Pattern Recognition, Track 1
%   - Computer Vision and Robotics, Quebec City, Canada, August 2001.
%
%IN:
%   A - HxWx3 uint8 image for segmentation.
%   h_s - scalar parameter on.
%   h_r - scalar parameter on.
%   min_sz - scalar indicating the minimum number of pixels per segment.
%   W - HxW single matrix containing weights from edge detection for
%       synergism with segmentation. See [2]. Default: uniform weight.
%
%OUT:
%   S - HxW uint32 segmentation matrix, each value of which gives the index
%       of the region said pixel belongs to, from 1 to num_segments.

% $Id: vgg_segment_ms.m,v 1.1 2007/12/10 10:59:31 ojw Exp $

function varargout = vgg_segment_ms(varargin)
funcName = mfilename;
sd = 'seg_ms/';
sourceList = {['-I' sd], [funcName '.cxx'], [sd 'msImageProcessor.cpp'],...
              [sd 'ms.cpp'], [sd 'rlist.cpp'], [sd 'RAList.cpp'],...
              [sd 'msSys.cpp']};
vgg_mexcompile_script; % Compilation happens in this script
return