function output_dir = cpad_analysis_output_dir(analysis_name, param_set, P, info)
%CPAD_ANALYSIS_OUTPUT_DIR Return (and create) an analysis output directory.
%
%   dir = cpad_analysis_output_dir(analysis_name)
%       <root>/<analysis_name>/
%   dir = cpad_analysis_output_dir(analysis_name, param_set)
%       <root>/<analysis_name>/<param_set>/
%   dir = cpad_analysis_output_dir(analysis_name, param_set, P, info)
%       Same, and (re)writes metadata.txt / metadata.mat in that folder
%       from the parameter struct P and an optional struct of run settings.
%
% analysis_name is the analysis script's name (e.g. 'Sweep2DCad'); param_set
% is a short label of the parameter set (e.g. 'M5_kEoffEngaged0.05'). Rerunning
% the same parameter set overwrites its files; no dates are stored.
%
% <root> is CPAD_RESULTS_ROOT when set, otherwise <project>/Results.

results_root = getenv('CPAD_RESULTS_ROOT');
if isempty(results_root)
    results_root = fullfile(fileparts(mfilename('fullpath')), 'Results');
end

output_dir = fullfile(results_root, analysis_name);
if nargin >= 2 && ~isempty(param_set)
    output_dir = fullfile(output_dir, param_set);
end
if ~exist(output_dir, 'dir')
    mkdir(output_dir);
end

if nargin >= 3 && ~isempty(P)
    if nargin < 4
        info = struct();
    end
    write_metadata(output_dir, analysis_name, param_set, P, info);
end

end

function write_metadata(output_dir, analysis_name, param_set, P, info) %#ok<INUSL>
parameters = P; %#ok<NASGU>
run_settings = info; %#ok<NASGU>
save(fullfile(output_dir, 'metadata.mat'), 'analysis_name', 'param_set', ...
    'parameters', 'run_settings');

fid = fopen(fullfile(output_dir, 'metadata.txt'), 'w');
if fid < 0
    error('cpad_analysis_output_dir:MetadataOpenFailed', ...
        'Could not write metadata in %s', output_dir);
end
cleanup = onCleanup(@() fclose(fid)); %#ok<NASGU>
fprintf(fid, 'Analysis: %s\n', analysis_name);
fprintf(fid, 'Parameter set: %s\n\n', param_set);
fprintf(fid, 'Parameters (default_parameters with this run''s overrides):\n');
write_struct(fid, P, '  ');
if ~isempty(fieldnames(info))
    fprintf(fid, '\nRun settings:\n');
    write_struct(fid, info, '  ');
end
end

function write_struct(fid, S, indent)
names = fieldnames(S);
for k = 1:numel(names)
    value = S.(names{k});
    if isstruct(value) && isscalar(value)
        fprintf(fid, '%s%s:\n', indent, names{k});
        write_struct(fid, value, [indent '  ']);
    elseif ischar(value)
        fprintf(fid, '%s%s = %s\n', indent, names{k}, value);
    elseif isnumeric(value) || islogical(value)
        fprintf(fid, '%s%s = %s\n', indent, names{k}, mat2str(value, 6));
    elseif iscellstr(value)
        fprintf(fid, '%s%s = {%s}\n', indent, names{k}, strjoin(value, ', '));
    elseif isa(value, 'function_handle')
        fprintf(fid, '%s%s = %s\n', indent, names{k}, func2str(value));
    else
        fprintf(fid, '%s%s = <%s>\n', indent, names{k}, class(value));
    end
end
end
