%VGG_MEXCOMPILE_SCRIPT  Compilation helper script for vgg mex files
%
% Script that should be placed in the m-file of a mexed function, after the
% following lines of code, thus:
%
%function varargout = function_name_here(varargin) 
%   funcName = mfilename;
%   sourceList = {[funcName '.cxx'], other_files};
%   vgg_mexcompile_script
%
% The script will compile the source inline (i.e. compile & then run) if
% the function has been called without first being compiled.
%
% If the function is being called by vgg_mexall, such that the first
% varargin{1} == 'compile' and varargin{2} == last_compilation_time, then
% the script calls vgg_mexcompile to compile the function.

% $Id: vgg_mexcompile_script.m,v 1.2 2009/09/13 20:34:58 ojw Exp $

if nargin == 2 && isequal(varargin{1}, 'compile')
    % Standard vgg_matlab compilation behaviour
    [varargout{1:nargout}] = vgg_mexcompile(funcName, sourceList, varargin{2});
else
    % Function called without first being compiled
    warning('Missing MEX-file: %s. Will attempt to compile and run.', funcName);
    retval = vgg_mexcompile(funcName, sourceList);
    if retval > 0
        [varargout{1:nargout}] = feval(funcName, varargin{:});
    else
        error('Unable to compile %s.', funcName);
    end
end