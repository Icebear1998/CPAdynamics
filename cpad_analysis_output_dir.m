function output_dir = cpad_analysis_output_dir(analysis_type, default_root)
%CPAD_ANALYSIS_OUTPUT_DIR Return an analysis output directory and create it.
%
% The aggregate publication runner sets CPAD_RESULTS_ROOT so existing
% analyses can write into a single publication-specific output tree.

if nargin < 2 || isempty(default_root)
    default_root = 'Results';
end

results_root = getenv('CPAD_RESULTS_ROOT');
if isempty(results_root)
    results_root = default_root;
end

output_dir = fullfile(results_root, analysis_type);
if ~exist(output_dir, 'dir')
    mkdir(output_dir);
end

end
