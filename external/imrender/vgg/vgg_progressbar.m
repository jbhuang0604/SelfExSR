function retval = vgg_progressbar(tag, proportion, min_update_interval)
%VGG_PROGRESSBAR  Simple progress bar implementation
%
%   retval = vgg_progressbar(tag, proportion[, min_update_interval])
%
% Starts, updates and closes a progress bar according to the proportion of
% time left.
%
% IN:
%   tag - String that appears on progress bar, specific to each function
%         calling VGG_PROGRESSBAR.
%   proportion - Proportion of time elapsed, between 0 and 1.
%   min_update_interval - Minimum time (in seconds) between updates of the
%                         progress bar. The value is kept while the bar is
%                         alive. Default: 0.5.
%
% OUT:
%   retval - 2 iff time since last report > min_update_interval,
%            1 iff progress bar initialized or reset,
%            0 otherwise.

% $Id: vgg_progressbar.m,v 1.1 2007/12/07 11:27:55 ojw Exp $
% Based on Andrew Fitzgibbon's awf_progressbar

% Check the input arguments
if nargin < 2
    error('Atleast 2 input arguments expected');
end
if ~isstr(tag)
  error('First argument should be a string');
end
if ~isscalar(proportion) | proportion < 0
  error('Second argument should be a non-negative scalar');
end
if nargin > 2
    if ~isscalar(min_update_interval) | min_update_interval < 0
        error('Third argument should be a non-negative scalar');
    end
end
retval = 0;

% Only display if we can
if ~usejava('awt')
    return
end
    
% Ensure the global data structure exists
global vgg_progressbar_data
try
  vgg_progressbar_data;
catch
  vgg_progressbar_data.xx = 0;
end

% Record the time
curr_time = clock;

% Check the tag exists
tag_title = tag;
tag = tag(isstrprop(tag, 'alphanum'));
if isfield(vgg_progressbar_data, tag)
    % Cache the data structure
    info = vgg_progressbar_data.(tag);
else
    % No tag by this name
    if proportion >= 1
        % No need to create one
        return
    end
    
    % Create a data structure for this tag
    info.bar = [];
    info.min_update = 0.5; % Default seconds between updates
    info.prop = proportion;
    info.start_prop = proportion;
    info.timer = curr_time;
    info.last_update = curr_time;
    retval = 1;
end

% Update the minimum update interval if a new one is given
if nargin > 2
    info.min_update = min_update_interval;
end

if proportion >= 1
    % Close the progress bar
    close(info.bar);
    drawnow;
    vgg_progressbar_data = rmfield(vgg_progressbar_data, tag);
    return
end

% Check to see if we haven't started again
if proportion < info.prop
    % Reset the information
    info.start_prop = proportion;
    info.timer = curr_time;
    info.last_update = info.timer;
    
    % Update the progress bar
    retval = 1;
elseif etime(curr_time, info.last_update) >= info.min_update
    % An update of the progress bar is required   
    if (proportion - info.start_prop) > 0
        retval = 2;
    else
        retval = 1;
    end
end

switch retval
    case 1
        newtitle = 'Starting...';
    case 2
        info.last_update = curr_time;
        t_elapsed = etime(curr_time, info.timer);
        t_remaining = ((1 - proportion) * t_elapsed) / (proportion - info.start_prop);
        newtitle = sprintf('Elapsed: %s', timestr(t_elapsed));
        if proportion > 0.01 | t_elapsed > 30
            if t_remaining < 600
                newtitle = sprintf('%s, Remaining: %s', newtitle, timestr(t_remaining));
            else
                newtitle = sprintf('%s, ETA: %s', newtitle, datestr(datenum(curr_time) + (t_remaining * 1.15741e-5), 0));
            end
        end
    otherwise
        return
end

% Update the waitbar
if ishandle(info.bar) 
    waitbar(proportion, info.bar, newtitle);
else
    info.bar = waitbar(proportion, newtitle, 'Name', tag_title);
end

% Update our global variable with the changes to this tag
info.prop = proportion;
vgg_progressbar_data.(tag) = info;
return

% Time string function
function str = timestr(t)
s = rem(t, 60);
m = rem(floor(t/60), 60);
h = floor(t/3600);

if h > 0
    str= sprintf('%dh%02dm%02.0fs', h, m, s);
elseif m > 0
    str = sprintf('%dm%02.0fs', m, s);
else
    str = sprintf('%2.1fs', s);
end
return
