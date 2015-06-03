function [A info] = ibr_render(options, Pout)
%IBR_RENDER  Rendering wrapper function
%
%   [A info] = ibr_render(options, Pout)
%
% Wrapper function for the various rendering methods in OJW's IBR toolbox.
% This function caches the input data and calls the relevant rendering
% function.
%
%IN:
%   options - a structure containing input parameters required by the
%             various setup and rendering functions. See help text of these
%             functions for details.
%   Pout - 3x4 projection matrix of the desired output view.
%
%OUT:
%   A - HxWxC rendered image.
%   info - structure containing additional output infomation, as given by
%          the specific rendering function used.

% $Id: ibr_render.m,v 1.2 2008/06/04 11:12:44 ojw Exp $

if nargin < 2
    Pout = [];
end

% Cache the input data
[images P disps] = ojw_setup(options, Pout);

% Get the output image dimensions
if isfield(options, 'dim_out') && ~isempty(options.dim_out)
    sz = options.dim_out([4 3]);
else
    % Set the default ouput dimensions
    sz = [size(images{1}, 1) size(images{1}, 2)];
end

% Call the rendering function
[A info] = feval(options.render_func, images, P, disps, sz, options);

return
