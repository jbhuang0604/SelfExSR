%IBR_GEN_CLIQUES  Generate graph cliques given visibilities
%
%   [U P PI T TI] = ibr_gen_cliques(IA, VA, V, Kocc, method)
%
% Generates graph cliques for the algorithm described in Woodford etal.'s
% BMVC 2007 paper:
%   "On New View Synthesis Using Multiview Stereo".
%
% IN:
%   IA - 1xN cell array of input images.
%   VA - 
%   V - 
%
% OUT:
%   A - (height)x(width)xC rendered image.
%   D - (height)x(width) disparity map for output image
%   V - (height)x(width)xN boolean array indicating the visibility of the
%       output pixels in the input images.

% $Id: ibr_gen_cliques.m,v 1.1 2007/12/07 11:27:52 ojw Exp $

function varargout = ibr_gen_cliques(varargin)
funcName = mfilename;
sourceList = {[funcName '.cxx']}; % Cell array of source files
vgg_mexcompile_script; % Compilation happens in this script
return