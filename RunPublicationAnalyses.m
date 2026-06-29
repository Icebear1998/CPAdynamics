function results_root = RunPublicationAnalyses(results_root)
%RUNPUBLICATIONANALYSES Run all paper-facing CPAdynamics analyses.
%
% Windows command from a MATLAB-enabled machine:
%   matlab -batch "cd('C:\path\to\CPAdynamics'); RunPublicationAnalyses"
%
% Optional custom output folder:
%   matlab -batch "cd('C:\path\to\CPAdynamics'); RunPublicationAnalyses('D:\CPAdynamicsPublicationResults')"
%
% The runner writes to a stable publication results folder by default:
%   PublicationResults/
%     MainResults/
%     SupplementalResults/

project_root = fileparts(mfilename('fullpath'));
if nargin < 1 || isempty(results_root)
    results_root = fullfile(project_root, 'PublicationResults');
end

if isa(results_root, 'string')
    results_root = char(results_root);
end

if ~exist(results_root, 'dir')
    mkdir(results_root);
end

main_results_root = fullfile(results_root, 'MainResults');
supplemental_results_root = fullfile(results_root, 'SupplementalResults');
if ~exist(main_results_root, 'dir'), mkdir(main_results_root); end
if ~exist(supplemental_results_root, 'dir'), mkdir(supplemental_results_root); end

old_results_root = getenv('CPAD_RESULTS_ROOT');
old_force_save = getenv('CPAD_FORCE_SAVE');
old_skip_clear = getenv('CPAD_SKIP_SCRIPT_CLEAR');
cleanup_env = onCleanup(@() restore_environment(old_results_root, old_force_save, old_skip_clear));

addpath(project_root);
setenv('CPAD_FORCE_SAVE', 'true');
setenv('CPAD_SKIP_SCRIPT_CLEAR', 'true');

log_file = fullfile(results_root, 'publication_analysis_run.log');
diary(log_file);
diary on;
cleanup_diary = onCleanup(@() diary('off'));

fprintf('=== CPAdynamics Publication Analysis Runner ===\n');
fprintf('Started: %s\n', datestr(now));
fprintf('Project root: %s\n', project_root);
fprintf('Results root: %s\n\n', results_root);

write_manifest(results_root, main_results_root, supplemental_results_root);

main_scripts = { ...
    'EBindingNumberVsCad.m', ...
    'PlotEBindingProfile.m', ...
    'GeneLengthGenerateGrid.m', ...
    'GeneLengthBuildInterpolation.m', ...
    'GeneLengthAnalyze.m', ...
    'SimulateCpaAssembly.m', ...
    'SweepParameterPasUsage.m' ...
};

supplemental_scripts = { ...
    'Sweep2DkHdkEdCad.m', ...
    'ParameterSweep1D.m' ...
};

fprintf('--- Main Results ---\n');
setenv('CPAD_RESULTS_ROOT', main_results_root);
run_script_group(project_root, main_scripts);

fprintf('\n--- Supplemental Results ---\n');
setenv('CPAD_RESULTS_ROOT', supplemental_results_root);
run_script_group(project_root, supplemental_scripts);

fprintf('\n=== Publication analyses complete ===\n');
fprintf('Finished: %s\n', datestr(now));
fprintf('Results root: %s\n', results_root);
fprintf('Run log: %s\n', log_file);

clear cleanup_diary cleanup_env

end

function run_script_group(project_root, script_names)
for idx = 1:numel(script_names)
    run_analysis_script(project_root, script_names{idx});
end
end

function run_analysis_script(project_root, script_name)
script_path = fullfile(project_root, script_name);
if ~exist(script_path, 'file')
    error('RunPublicationAnalyses:MissingScript', ...
        'Required analysis script not found: %s', script_path);
end

fprintf('\n[%s] Starting %s\n', datestr(now, 'HH:MM:SS'), script_name);
timer_id = tic;
run(script_path);
elapsed_seconds = toc(timer_id);
fprintf('[%s] Completed %s in %.1f minutes\n', ...
    datestr(now, 'HH:MM:SS'), script_name, elapsed_seconds / 60);
end

function write_manifest(results_root, main_results_root, supplemental_results_root)
manifest_file = fullfile(results_root, 'publication_analysis_manifest.txt');
fid = fopen(manifest_file, 'w');
if fid < 0
    error('RunPublicationAnalyses:ManifestOpenFailed', ...
        'Could not create manifest file: %s', manifest_file);
end
cleanup_file = onCleanup(@() fclose(fid));

fprintf(fid, 'CPAdynamics Publication Analysis Manifest\n');
fprintf(fid, 'Generated: %s\n\n', datestr(now));
fprintf(fid, 'Main results folder:\n%s\n\n', main_results_root);
fprintf(fid, 'Supplemental results folder:\n%s\n\n', supplemental_results_root);

fprintf(fid, 'Main results scripts, in order:\n');
fprintf(fid, '1. EBindingNumberVsCad.m\n');
fprintf(fid, '2. PlotEBindingProfile.m\n');
fprintf(fid, '3. GeneLengthGenerateGrid.m\n');
fprintf(fid, '4. GeneLengthBuildInterpolation.m\n');
fprintf(fid, '5. GeneLengthAnalyze.m\n');
fprintf(fid, '6. SimulateCpaAssembly.m\n');
fprintf(fid, '7. SweepParameterPasUsage.m\n\n');

fprintf(fid, 'Supplemental results scripts, in order:\n');
fprintf(fid, '1. Sweep2DkHdkEdCad.m\n');
fprintf(fid, '2. ParameterSweep1D.m\n\n');

fprintf(fid, 'Windows command:\n');
fprintf(fid, 'matlab -batch "cd(''C:\\path\\to\\CPAdynamics''); RunPublicationAnalyses"\n');
end

function restore_environment(old_results_root, old_force_save, old_skip_clear)
setenv('CPAD_RESULTS_ROOT', old_results_root);
setenv('CPAD_FORCE_SAVE', old_force_save);
setenv('CPAD_SKIP_SCRIPT_CLEAR', old_skip_clear);
end
