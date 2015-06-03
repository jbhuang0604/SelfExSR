function save_progress(filename, varargin)
if isempty(filename)
    return
end
% Get the global data structure
global save_progress_data
% Generate the unique tag name
tag = filename(isstrprop(filename, 'alphanum'));
% For each variable
n = numel(varargin);
varnames = cell(n, 1);
for a = 1:n
    % Get the iteration number
    try
        iter = save_progress_data.(tag).(varargin{a}) + 1;
    catch
        iter = 1; % First one
    end
    % Update the save iteration
    save_progress_data.(tag).(varargin{a}) = iter;
    % Capture the variable being saved
    var = evalin('caller', varargin{a});
    % Store the variable locally with the sequential name
    varnames{a} = sprintf('%s%d', varargin{a}, iter);
    eval([varnames{a} ' = var;']);
end
var = exist(filename, 'file');
if var == 0 || var == 7
    % Create the file; don't compress for speed & compatibility
    save(filename, '-v6', varnames{:});
else
    % Append variables to the file
    save(filename, '-append', varnames{:});
end
return