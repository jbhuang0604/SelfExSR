%VGG_SEGMENT_GB  Graph-based image segmentation
%
%   S = vgg_segment_gb(A, sigma, K, min_sz[, compress])
%
% Segmentation of an image using the graph-based method described in:
%   "Efficient Graph-Based Image Segmentation.", Pedro F. Felzenszwalb and
%   Daniel P. Huttenlocher. International Journal of Computer Vision,
%   Volume 59, Number 2, September 2004. 
% Code downloaded from:
%   http://people.cs.uchicago.edu/~pff/segment/
% 
%IN:
%   A - HxWx3 uint8 image for segmentation.
%   sigma - scalar parameter on smoothing kernel to use prior to
%           segmentation.
%   k - scalar parameter on prefered segment size.
%   min_sz - scalar indicating the minimum number of pixels per segment.
%   compress - scalar boolean indicating whether the user wants the segment
%              indices compressed to the range [1 num_segments].
%              Default: 0.
%
%OUT:
%   S - HxW uint32 segmentation matrix, each value of which gives the index
%       of the region said pixel belongs to.

% $Id: vgg_segment_gb.m,v 1.1 2007/12/10 10:59:30 ojw Exp $

function varargout = vgg_segment_gb(varargin)
funcName = mfilename;
sd = 'seg_gb/';
sourceList = {['-I' sd], [funcName '.cxx']};
vgg_mexcompile_script; % Compilation happens in this script
return