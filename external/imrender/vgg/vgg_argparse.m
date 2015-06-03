function [opts,rem_opts] = vgg_argparse(opts,varargin)

%VGG_ARGPARSE  Parse variable arguments into a structure
%  opts = vgg_argparse(inopts,varargin)
%    inopts: structure (cells array) listing valid members and default values
%    varargin: variable arguments of form '<name>',<value>,...
%    opts: opts modified by varargin
%
%  Example:
%    function f = foo(varargin)
%    opts = vgg_argparse(struct('maxiters',10,'verbose',0), varargin)
%    ...
%
%  An unknown option (ie, present in varargin but absent in inopts)
%  causes an error. Calling the function as 
%  [opts,rem_opts] = vgg_argparse(inopts,varargin) returns the unknown
%  option(s) in rem_opts for later use rather than causes an error.
%
%  May also use OPTS = VGG_ARGPARSE(OPTS, ASTRUCT) where ASTRUCT is a struct
%  of options.

% Author: Mark Everingham <me@robots.ox.ac.uk>
% modified by werner, Jan 03
% Date: 16 Jan 02

if iscell(opts)
  opts=struct(opts{:});
end

if length(varargin) & iscell(varargin{1})
    if isempty(varargin{1})
        inopts = struct([]);
    else
        inopts=struct(varargin{1}{:});
    end
else
    if isempty(varargin)
        inopts = struct([]);
    elseif isstruct(varargin{1})
        inopts = varargin{1};
    else
        inopts=struct(varargin{:});
    end
end

rem_opts = [];
fn = fieldnames(inopts);
for i=1:length(fn)
    if isfield(opts,fn{i})
        %opts.(fn{i})=inopts.(fn{i});
        opts = setfield(opts,fn{i},getfield(inopts,fn{i}));
    else
        if nargout < 2
            error(sprintf('bad argument: ''%s''', fn{i}));
        else
            %rem_opts.(fn{i}) = inopts.(fn{i});
            rem_opts = setfield(rem_opts,fn{i},getfield(inopts,fn{i}));
        end
    end
end