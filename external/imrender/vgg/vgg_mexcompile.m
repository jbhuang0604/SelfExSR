function okay = vgg_mexcompile(funcName, varargin)
%VGG_MEXCOMPILE  Mex compile helper function
%
%   okay = vgg_mexcompile(funcName)
%   okay = vgg_mexcompile(..., sourceList)
%   okay = vgg_mexcompile(..., lastCompiled)
%
% Compile mexed function, given an optional list of source files. Can
% optionally check if source files have been modified since the last
% compilation, and only compile if they have.
%
%IN:
%   funcName - string containg the name of the function to compile
%   sourceList - cell array of source files to be compiled. Default:
%                {[funcName '.cxx']}.
%   lastCompiled - datenum of the current mex file. Default: 0 (i.e. force
%                  compilation).
%
%OUT:
%   okay - 1: function compiled; 0: up-to-date, no need to compile; -1:
%          compilation failed.

% $Id: vgg_mexcompile.m,v 1.3 2009/09/13 20:34:58 ojw Exp $

% Set defaults for optional inputs
sourceList = [funcName '.cxx'];
lastCompiled = 0;
% Parse inputs
for a = 1:numel(varargin)
   if iscell(varargin{a}) && ischar(varargin{a}{1})
       sourceList = varargin{a};
   elseif isnumeric(varargin{a}) && isscalar(varargin{a})
       lastCompiled = varargin{a};
   end
end

% Go to the directory containing the file
currDir = cd;
sourceDir = fileparts(which(funcName));
cd(sourceDir);

if lastCompiled
    compile = false;
    okay = 0;
    % Compile if current mex file is older than any of the source files
    for a = 1:numel(sourceList)
        dirSource = dir(sourceList{a});
        if ~isempty(dirSource) && datenum(dirSource.date) > lastCompiled
            compile = true;
            break;
        end
    end
else
    compile = true;
end

% Compile if we need to
if compile
    % Set the compiler flags
    flags = ['-O -I"' sourceDir '"'];
    switch mexext
        case 'mexsol'
            flags = [flags ' CC=gcc CFLAGS=-fPIC'];
        case {'mexglx', 'mexa64'}
            str = '"-O3 -ffast-math -funroll-loops"';
            flags = sprintf('%s CXXOPTIMFLAGS=%s LDCXXOPTIMFLAGS=%s LDOPTIMFLAGS=%s', flags, str, str, str);
        otherwise
    end

    % Call mex to compile the code
    cmd = sprintf('mex %s%s', flags, sprintf(' "%s"', sourceList{:}));
    disp(cmd);
    try
        eval(cmd);
        okay = 1;
    catch
        okay = -1;
        fprintf('ERROR while compiling %s\n', funcName);
    end
end

% Return to the original directory
cd(currDir);
return
