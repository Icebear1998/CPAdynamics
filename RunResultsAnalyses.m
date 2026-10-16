function RunResultsAnalyses(sections)
%RUNRESULTSANALYSES Raw data behind every figure and value in the Results draft.
%
%   RunResultsAnalyses                         % all sections
%   RunResultsAnalyses({'GlobalSample', 'Apa'}) % selected sections
%
% R2018b has no -batch flag. From a shell:
%   matlab -nodesktop -nosplash -wait -r "cd('<repo>'); RunResultsAnalyses; exit"
%
% Each section writes CSV tables and metadata to
% Results/RunResultsAnalyses/<section>/. GeneLength also reruns the
% three-step gene-length pipeline at the strong PAS into
% Results/GeneLengthAnalysis/M5_kHoff0.01/. The figures and the values that
% fill the draft come from python/assemble_results_paper.py.
%
% Baselines: weak PAS (default_parameters, kHoff = 1) for R2-R4 and a strong
% PAS (kHoff = 0.01) for R5-R6. Parameter samples are log-uniform over
% P.ranges. Failed solves are kept as NaN rows with their error message.

if nargin < 1 || isempty(sections)
    sections = {'Abundance', 'BaselineWeakPAS', 'BaselineStrongPAS', 'Maps2D', ...
        'GlobalSample', 'SingleTightSite', 'GeneLength', 'Apa', 'ApaBinding'};
end
if ischar(sections)
    sections = {sections};
end

project_root = fileparts(mfilename('fullpath'));
addpath(project_root);
root_dir = cpad_analysis_output_dir('RunResultsAnalyses');
diary(fullfile(root_dir, 'run_log.txt'));
diary on;
cleanup_diary = onCleanup(@() diary('off')); %#ok<NASGU>

if isempty(gcp('nocreate'))
    parpool('local');
end
unreached_id = 'calculate_full_pas_cleavage_profile:ThresholdNotReached';
warning('off', unreached_id);
pctRunOnAll(sprintf('warning(''off'', ''%s'');', unreached_id));

P = default_parameters();
fprintf('=== RunResultsAnalyses started %s ===\n', datestr(now));
for k = 1:numel(sections)
    section_timer = tic;
    fprintf('\n--- %s ---\n', sections{k});
    switch sections{k}
        case 'Abundance'
            run_abundance(P);
        case 'BaselineWeakPAS'
            run_baseline_weak(P);
        case 'BaselineStrongPAS'
            run_baseline_strong(P);
        case 'Maps2D'
            run_maps(P);
        case 'GlobalSample'
            run_global_sample(P);
        case 'SingleTightSite'
            run_single_tight_site(P);
        case 'GeneLength'
            run_gene_length(project_root);
        case 'Apa'
            run_apa(P);
        case 'ApaBinding'
            run_apa_binding(P);
        otherwise
            error('RunResultsAnalyses:UnknownSection', 'Unknown section: %s', sections{k});
    end
    close all;
    fprintf('--- %s done in %.1f min ---\n', sections{k}, toc(section_timer)/60);
end
fprintf('\n=== RunResultsAnalyses finished %s ===\n', datestr(now));
end

%% ------------------------------------------------------------------------
function kHoff = strong_kHoff()
kHoff = 0.01;  % Kd ~20 nM, middle of 3-40 nM (appendix B)
end

function P = with_kHoff(P, kHoff)
% Change PAS strength. kEoff_engaged is one value for all PAS strengths
% (appendix C), so it is not rescaled.
P.kHoff = kHoff;
end

function run_abundance(P)
% Copy numbers per HeLa cell from Appendix G, Table G1 (PaxDB calibrated
% against Zeiler et al. 2012). Pol II is the model total, P.Pol_total.
protein = {'CPSF-160'; 'WDR33'; 'CPSF-30'; 'CPSF-100'; 'CPSF-73'; ...
    'CstF-50'; 'CstF-64'; 'CstF-77'; 'Active Pol II'};
copies = [102000; 214000; 292000; 241000; 187000; 332000; 371000; 268000; P.Pol_total];
group = {'CPSF'; 'CPSF'; 'CPSF'; 'CPSF'; 'CPSF'; 'CstF'; 'CstF'; 'CstF'; 'Pol II'};
source = [repmat({'Appendix G, Table G1'}, 8, 1); {'default_parameters P.Pol_total'}];
out_dir = cpad_analysis_output_dir('RunResultsAnalyses', 'Abundance', P, ...
    struct('nucleus_volume_L', 3.74e-13));
writetable(table(protein, group, copies, source), fullfile(out_dir, 'copy_numbers.csv'));
writetable(table(P.E_total, P.Pol_total, P.kEon, P.kEoff, 3.74e-13, ...
    P.ranges.E_total(1), P.ranges.E_total(2), P.ranges.Pol_total(1), P.ranges.Pol_total(2), ...
    'VariableNames', {'E_total', 'Pol_total', 'kEon', 'kEoff', 'nucleus_volume_L', ...
    'E_total_min', 'E_total_max', 'Pol_total_min', 'Pol_total_max'}), ...
    fullfile(out_dir, 'pools.csv'));
% Sensitivity ranges, so the figures draw the current appendix bounds.
names = fieldnames(P.ranges);
bounds = cell2mat(cellfun(@(n) P.ranges.(n)(:)', names, 'UniformOutput', false));
writetable(table(names, bounds(:, 1), bounds(:, 2), P_values(P, names), ...
    'VariableNames', {'param', 'min', 'max', 'default'}), fullfile(out_dir, 'ranges.csv'));
end

function v = P_values(P, names)
% Default value of each ranged parameter (derived ratios from their parts).
v = NaN(numel(names), 1);
for k = 1:numel(names)
    switch names{k}
        case 'kHd'
            v(k) = P.kHoff/P.kHon;
        case 'kEd'
            v(k) = P.kEoff/P.kEon;
        otherwise
            v(k) = P.(names{k});
    end
end
end

%% ------------------------------------------------------------------------
function out_dir = run_cad_vs_M(P, section, M_values)
% CAD, CDF, loading and profiles for each M at the PAS strength in P.
out_dir = cpad_analysis_output_dir('RunResultsAnalyses', section, P, ...
    struct('M_values', M_values));
S = cell(numel(M_values), 1);
parfor k = 1:numel(M_values)
    S{k} = solve_case(P, M_values(k));
end
write_summary(out_dir, 'cad_vs_M.csv', 'M', M_values(:), S, P.E_total);
write_cdfs(out_dir, 'cdf_by_M.csv', 'M', M_values(:), S);
write_profiles(out_dir, 'profiles_by_M.csv', 'M', M_values(:), S, P.L_a);
end

function run_baseline_strong(P)
% Same as the weak-PAS baseline, at the strong PAS used in R5-R6.
P = with_kHoff(P, strong_kHoff());
run_cad_vs_M(P, 'BaselineStrongPAS', 1:10);
end

function run_baseline_weak(P)
% CAD vs M = 1..10 at the weak PAS; CAD vs E_total and vs Pol_total for M = 1, 5.
ratios = 2.^(-2:0.5:3);
out_dir = run_cad_vs_M(P, 'BaselineWeakPAS', 1:10);

pol_values = unique([logspace(log10(P.ranges.Pol_total(1)), log10(P.ranges.Pol_total(2)), 8), P.Pol_total]);
[M_job, pol_job] = ndgrid([1 5], pol_values);
M_job = M_job(:);
pol_job = pol_job(:);
SP = cell(numel(M_job), 1);
parfor j = 1:numel(M_job)
    P_run = P;
    P_run.Pol_total = pol_job(j);
    SP{j} = solve_case(P_run, M_job(j));
end
T = summary_table(SP, repmat(P.E_total, numel(SP), 1));
T = [table(M_job, pol_job, 'VariableNames', {'M', 'Pol_total'}), T];
writetable(T, fullfile(out_dir, 'cad_vs_Poltotal.csv'));

[M_job, ratio_job] = ndgrid([1 5], ratios);
M_job = M_job(:);
ratio_job = ratio_job(:);
SE = cell(numel(M_job), 1);
parfor j = 1:numel(M_job)
    P_run = P;
    P_run.E_total = ratio_job(j)*P.Pol_total;
    SE{j} = solve_case(P_run, M_job(j));
end
T = summary_table(SE, ratio_job*P.Pol_total);
T = [table(M_job, ratio_job, 'VariableNames', {'M', 'E_to_Pol'}), T];
writetable(T, fullfile(out_dir, 'cad_vs_Etotal.csv'));
end

%% ------------------------------------------------------------------------
function run_maps(P)
% Pairwise CAD maps over the kHd, kEd and kc ranges, as in Sweep2DCad.m.
% kEoff_engaged stays at its default along both axes.
names = {'kHd', 'kEd', 'kc'};
n_points = 10;
pairs = [1 2; 3 1; 3 2];
axis_values = cell(1, 3);
for a = 1:3
    b = P.ranges.(names{a});
    axis_values{a} = logspace(log10(b(1)), log10(b(2)), n_points);
end
out_dir = cpad_analysis_output_dir('RunResultsAnalyses', 'Maps2D', P, ...
    struct('n_points', n_points, 'M_values', [1 5]));

[M_job, pair_job, ix_job, iy_job] = ndgrid([1 5], 1:size(pairs, 1), 1:n_points, 1:n_points);
M_job = M_job(:); pair_job = pair_job(:); ix_job = ix_job(:); iy_job = iy_job(:);
J = numel(M_job);
x_job = NaN(J, 1); y_job = NaN(J, 1);
x_name = cell(J, 1); y_name = cell(J, 1);
for j = 1:J
    x_name{j} = names{pairs(pair_job(j), 1)};
    y_name{j} = names{pairs(pair_job(j), 2)};
    x_job(j) = axis_values{pairs(pair_job(j), 1)}(ix_job(j));
    y_job(j) = axis_values{pairs(pair_job(j), 2)}(iy_job(j));
end
S = cell(J, 1);
parfor j = 1:J
    P_run = set_ratio(set_ratio(P, x_name{j}, x_job(j)), y_name{j}, y_job(j));
    S{j} = solve_case(P_run, M_job(j));
end
T = summary_table(S, repmat(P.E_total, J, 1));
T = [table(M_job, x_name, y_name, x_job, y_job, ...
    'VariableNames', {'M', 'x_name', 'y_name', 'x', 'y'}), T(:, {'CAD50', 'max_exit_cdf', 'ok', 'error'})];
writetable(T, fullfile(out_dir, 'maps.csv'));
writetable(table(P.kHoff/P.kHon, P.kEoff/P.kEon, P.kc, ...
    'VariableNames', {'kHd', 'kEd', 'kc'}), fullfile(out_dir, 'baseline_point.csv'));
end

function P = set_ratio(P, name, value)
switch name
    case 'kHd'
        P.kHoff = value*P.kHon;
    case 'kEd'
        P.kEoff = value*P.kEon;
    case 'kc'
        P.kc = value;
end
end

%% ------------------------------------------------------------------------
function run_global_sample(P)
% Paired M = 1 / M = 5 solves over a log-uniform sample of P.ranges, plus the
% kEoff a single site needs to match each set's M = 5 CAD50.
names = {'kHon', 'kHoff', 'kEon', 'kEoff', 'kc', 'k_e', 'k_e2', ...
    'kPon_slope', 'kPoff', 'E_total', 'Pol_total'};
n_sets = 300;
seed = 1;
bounds = zeros(numel(names), 2);
for j = 1:numel(names)
    bounds(j, :) = P.ranges.(names{j});
end
rng(seed, 'twister');
U = rand(n_sets, numel(names));
values = exp(log(bounds(:, 1))' + U.*(log(bounds(:, 2))-log(bounds(:, 1)))');
out_dir = cpad_analysis_output_dir('RunResultsAnalyses', 'GlobalSample', P, ...
    struct('n_sets', n_sets, 'seed', seed, 'sampling', 'log-uniform over P.ranges', ...
    'sampled', {names}));

columns = {'cad_M1', 'cad_M5', 'maxcdf_M1', 'maxcdf_M5', 'freeE_M1', 'freeE_M5', ...
    'loaded_M1', 'loaded_M5', 'sites_M1', 'sites_M5', 'kEoff_match', 'tighten_fold', 'both_ok'};
out = NaN(n_sets, numel(columns));
match_status = cell(n_sets, 1);
errors = cell(n_sets, 1);
parfor i = 1:n_sets
    [out(i, :), match_status{i}, errors{i}] = sample_one(P, names, values(i, :));
end
T = array2table([(1:n_sets)', values, out], 'VariableNames', [{'set'}, names, columns]);
T.match_status = match_status;
T.error = errors;
writetable(T, fullfile(out_dir, 'sample.csv'));
end

function [row, match_status, err] = sample_one(P, names, vals)
for j = 1:numel(names)
    P.(names{j}) = vals(j);
end
s1 = solve_case(P, 1);
s5 = solve_case(P, 5);
row = [s1.CAD50, s5.CAD50, s1.max_exit_cdf, s5.max_exit_cdf, ...
    s1.E_free/P.E_total, s5.E_free/P.E_total, s1.loaded_prePAS, s5.loaded_prePAS, ...
    s1.sites_prePAS, s5.sites_prePAS, NaN, NaN, double(s1.ok && s5.ok)];
err = strjoin({s1.error, s5.error}, ' | ');
match_status = 'not_attempted';
if s1.ok && s5.ok && isfinite(s5.CAD50)
    [k_match, match_status] = match_single_site(P, s5.CAD50, []);
    row(11) = k_match;
    row(12) = P.kEoff/k_match;
end
end

%% ------------------------------------------------------------------------
function run_single_tight_site(P)
% Five weak sites vs one site whose kEoff is tightened to the same CAD50.
s5 = solve_case(P, 5);
s1 = solve_case(P, 1);
[k_match, status] = match_single_site(P, s5.CAD50, []);
P_tight = P;
P_tight.kEoff = k_match;
s_tight = solve_case(P_tight, 1);
out_dir = cpad_analysis_output_dir('RunResultsAnalyses', 'SingleTightSite', P, ...
    struct('kEoff_match', k_match, 'match_status', status));
design = {'M5'; 'M1'; 'M1_tight'};
S = {s5; s1; s_tight};
T = summary_table(S, repmat(P.E_total, 3, 1));
T = [table(design, [5; 1; 1], [P.kEoff; P.kEoff; k_match], ...
    'VariableNames', {'design', 'M', 'kEoff'}), T];
writetable(T, fullfile(out_dir, 'summary.csv'));
write_profiles(out_dir, 'profiles.csv', 'design', design, S, P.L_a);
end

%% ------------------------------------------------------------------------
function run_gene_length(project_root)
% Genome-wide shared pools at the strong PAS (existing pipeline), then CAD vs
% length at those pools for three loading designs and at range extremes.
suffix = sprintf('_kHoff%g', strong_kHoff());
P_strong = with_kHoff(default_parameters(), strong_kHoff());
run_gene_length_pipeline(project_root, struct('kHoff', P_strong.kHoff, ...
    'kEoff_engaged', P_strong.kEoff_engaged), suffix);
set_dir = cpad_analysis_output_dir('GeneLengthAnalysis', ['M5' suffix]);
loaded = load(fullfile(set_dir, 'TCD_analysis.mat'), 'analysis_results');
A = loaded.analysis_results;
base = A.parameters.base_parameters;
after_PAS = A.parameters.after_PAS_length;
pools = [A.solution.R_free, A.solution.E_free];
L_values = A.TCD.gene_lengths(:);

out_dir = cpad_analysis_output_dir('RunResultsAnalyses', 'GeneLength', base, ...
    struct('source', set_dir, 'R_free', pools(1), 'E_free', pools(2)));
writetable(table(L_values, A.TCD.values(:, 1), A.TCD.max_exit_cdf(:), ...
    'VariableNames', {'L_bp', 'CAD50', 'max_exit_cdf'}), fullfile(out_dir, 'genome_cad_vs_length.csv'));
writetable(table(pools(1), pools(2), A.solution.num_active_genes, A.solution.R_total_target, ...
    A.solution.E_total_target, A.integration.probability_mass_in_grid, ...
    'VariableNames', {'R_free', 'E_free', 'num_genes', 'Pol_total', 'E_total', 'pdf_mass'}), ...
    fullfile(out_dir, 'pools.csv'));
writetable(table(A.integration.lengths(:), A.integration.weights(:), ...
    'VariableNames', {'L_bp', 'weight'}), fullfile(out_dir, 'length_weights.csv'));

% One tight site matched to M = 5 at 20 kb, at the same shared pools.
P20 = with_length(base, 20000, after_PAS);
s5_20 = solve_case(P20, 5, pools);
[k_match, status] = match_single_site(P20, s5_20.CAD50, pools);
designs = {'M5', 'M1', 'M1_tight'};
[d_job, L_job] = ndgrid(1:3, 1:numel(L_values));
d_job = d_job(:); L_job = L_job(:);
S = cell(numel(d_job), 1);
parfor j = 1:numel(d_job)
    P_run = with_length(base, L_values(L_job(j)), after_PAS);
    M = 5;
    if d_job(j) > 1, M = 1; end
    if d_job(j) == 3, P_run.kEoff = k_match; end
    S{j} = solve_case(P_run, M, pools);
end
T = summary_table(S, repmat(pools(2), numel(S), 1));
T = [table(designs(d_job)', L_values(L_job), 'VariableNames', {'design', 'L_bp'}), ...
    T(:, {'CAD50', 'max_exit_cdf', 'ok', 'error'})];
writetable(T, fullfile(out_dir, 'designs_vs_length.csv'));
writetable(table(k_match, base.kEoff/k_match, {status}, s5_20.CAD50, ...
    'VariableNames', {'kEoff_match', 'tighten_fold', 'status', 'CAD50_M5_20kb'}), ...
    fullfile(out_dir, 'tight_match.csv'));

% Ser2P and E profiles for a short and a long gene at the shared pools.
profile_lengths = [5000; 20000];
SL = cell(numel(profile_lengths), 1);
parfor k = 1:numel(profile_lengths)
    SL{k} = solve_case(with_length(base, profile_lengths(k), after_PAS), 5, pools);
end
write_profiles(out_dir, 'profiles_short_long.csv', 'L_bp', profile_lengths, SL, base.L_a);

% M = 5 curves at the range ends of the runway parameters, same pools.
P_ranges = default_parameters();
params = {'kPon_slope', 'kPoff', 'k_e', 'kc'};
[p_job, lvl_job, L_job] = ndgrid(1:numel(params), 1:2, 1:numel(L_values));
p_job = p_job(:); lvl_job = lvl_job(:); L_job = L_job(:);
value_job = NaN(numel(p_job), 1);
for j = 1:numel(p_job)
    value_job(j) = P_ranges.ranges.(params{p_job(j)})(lvl_job(j));
end
S = cell(numel(p_job), 1);
parfor j = 1:numel(p_job)
    P_run = with_length(base, L_values(L_job(j)), after_PAS);
    P_run.(params{p_job(j)}) = value_job(j);
    S{j} = solve_case(P_run, 5, pools);
end
levels = {'min', 'max'};
T = summary_table(S, repmat(pools(2), numel(S), 1));
T = [table(params(p_job)', levels(lvl_job)', value_job, L_values(L_job), ...
    'VariableNames', {'param', 'level', 'value', 'L_bp'}), ...
    T(:, {'CAD50', 'max_exit_cdf', 'ok', 'error'})];
writetable(T, fullfile(out_dir, 'extremes_vs_length.csv'));
end

function run_gene_length_pipeline(project_root, cpad_grid_overrides, cpad_grid_set_suffix) %#ok<INUSD>
% The grid script reads the two input names from this workspace.
run(fullfile(project_root, 'GeneLengthGenerateGrid.m'));
run(fullfile(project_root, 'GeneLengthBuildInterpolation.m'));
run(fullfile(project_root, 'GeneLengthAnalyze.m'));
end

function P = with_length(P, L, after_PAS)
P.PASposition = L;
P.geneLength_bp = L+after_PAS;
end

%% ------------------------------------------------------------------------
function run_apa(P)
% Proximal PAS usage at the strong PAS: CDFs by kHoff (30 kb window), and
% usage vs E/Pol II at baseline and at each parameter range end.
P_strong = with_kHoff(P, strong_kHoff());
kHoff_values = [0.0015 0.01 0.05 0.2 1 5];
ratios = 2.^(-3:0.5:4);
out_dir = cpad_analysis_output_dir('RunResultsAnalyses', 'Apa', P_strong, ...
    struct('kHoff_values', kHoff_values, 'E_to_Pol_ratios', ratios, ...
    'cdf_gene_length_bp', 50000));

S = cell(numel(kHoff_values), 1);
parfor k = 1:numel(kHoff_values)
    P_run = with_kHoff(P_strong, kHoff_values(k));
    P_run.geneLength_bp = 50000;
    S{k} = solve_case(P_run, 5);
end
write_cdfs(out_dir, 'cdf_by_kHoff.csv', 'kHoff', kHoff_values(:), S);

params = {'kEon', 'kEoff', 'kc', 'kPon_slope', 'kPoff', 'Pol_total', 'k_e', 'k_e2', 'kHon'};
case_param = [{'baseline'}, reshape(repmat(params, 2, 1), 1, [])];
case_level = [{'baseline'}, repmat({'min', 'max'}, 1, numel(params))];
case_value = NaN(1, numel(case_param));
for c = 2:numel(case_param)
    case_value(c) = P.ranges.(case_param{c})(1+strcmp(case_level{c}, 'max'));
end
[c_job, r_job] = ndgrid(1:numel(case_param), 1:numel(ratios));
c_job = c_job(:); r_job = r_job(:);
usage = NaN(numel(c_job), 3);
errors = cell(numel(c_job), 1);
parfor j = 1:numel(c_job)
    P_run = P_strong;
    if c_job(j) > 1
        P_run.(case_param{c_job(j)}) = case_value(c_job(j));
    end
    P_run.E_total = ratios(r_job(j))*P_run.Pol_total;
    s = solve_case(P_run, 5);
    errors{j} = s.error;
    if s.ok
        usage(j, :) = [interp1(s.distances_bp, s.cdf, 300), ...
            interp1(s.distances_bp, s.cdf, 1000), s.CAD50];
    end
end
T = table(case_param(c_job)', case_level(c_job)', case_value(c_job)', ratios(r_job)', ...
    usage(:, 1), usage(:, 2), usage(:, 3), errors, ...
    'VariableNames', {'param', 'level', 'value', 'E_to_Pol', 'usage_300', 'usage_1000', 'CAD50', 'error'});
writetable(T, fullfile(out_dir, 'usage_vs_Etotal.csv'));
end

function run_apa_binding(P)
% Why E does not bind more tightly: proximal usage at 300 bp vs E/Pol II over
% a kEoff sweep (M = 1 and M = 5), PAS-strength discrimination (kHoff = 1 vs
% the strong PAS) at the measured E, and one site matched to the M = 5 CAD50.
P_strong = with_kHoff(P, strong_kHoff());
kEoff_values = unique([logspace(-3, log10(P.ranges.kEoff(2)), 13), P.kEoff]);
ratios = 2.^(-3:0.5:4);
kHoff_pair = [strong_kHoff(), 1];
s5 = solve_case(P_strong, 5);
[k_match, status] = match_single_site(P_strong, s5.CAD50, []);
out_dir = cpad_analysis_output_dir('RunResultsAnalyses', 'ApaBinding', P_strong, ...
    struct('kEoff_values', kEoff_values, 'E_to_Pol_ratios', ratios, ...
    'kHoff_pair', kHoff_pair, 'kEoff_match', k_match, 'match_status', status));

% Design rows: M = 5 and M = 1 at each kEoff, plus the matched single site.
design_M = [5*ones(1, numel(kEoff_values)), ones(1, numel(kEoff_values)), 1];
design_k = [kEoff_values, kEoff_values, k_match];
design_name = [repmat({'sweep'}, 1, 2*numel(kEoff_values)), {'M1_matched'}];
[d_job, r_job] = ndgrid(1:numel(design_M), 1:numel(ratios));
d_job = d_job(:); r_job = r_job(:);
out = NaN(numel(d_job), 4);
errors = cell(numel(d_job), 1);
parfor j = 1:numel(d_job)
    P_run = P_strong;
    P_run.kEoff = design_k(d_job(j));
    P_run.E_total = ratios(r_job(j))*P_run.Pol_total;
    s = solve_case(P_run, design_M(d_job(j)));
    errors{j} = s.error;
    if s.ok
        out(j, :) = [interp1(s.distances_bp, s.cdf, 300), s.CAD50, ...
            s.E_free/P_run.E_total, s.loaded_prePAS];
    end
end
T = table(design_name(d_job)', design_M(d_job)', design_k(d_job)', ratios(r_job)', ...
    out(:, 1), out(:, 2), out(:, 3), out(:, 4), errors, ...
    'VariableNames', {'design', 'M', 'kEoff', 'E_to_Pol', 'usage_300', 'CAD50', ...
    'free_E_fraction', 'loaded_prePAS', 'error'});
writetable(T, fullfile(out_dir, 'usage_vs_Etotal_by_kEoff.csv'));

% PAS-strength discrimination at the measured E.
[d_job, h_job] = ndgrid(1:numel(design_M), 1:numel(kHoff_pair));
d_job = d_job(:); h_job = h_job(:);
usage = NaN(numel(d_job), 2);
parfor j = 1:numel(d_job)
    P_run = with_kHoff(P_strong, kHoff_pair(h_job(j)));
    P_run.kEoff = design_k(d_job(j));
    s = solve_case(P_run, design_M(d_job(j)));
    if s.ok
        usage(j, :) = [interp1(s.distances_bp, s.cdf, 300), s.CAD50];
    end
end
T = table(design_name(d_job)', design_M(d_job)', design_k(d_job)', kHoff_pair(h_job)', ...
    usage(:, 1), usage(:, 2), ...
    'VariableNames', {'design', 'M', 'kEoff', 'kHoff', 'usage_300', 'CAD50'});
writetable(T, fullfile(out_dir, 'pas_strength_by_kEoff.csv'));
end

%% ------------------------------------------------------------------------
function s = solve_case(P, M, free_pools)
% One steady state with the readouts used by the Results figures.
s = struct('ok', false, 'CAD50', NaN, 'max_exit_cdf', NaN, 'E_free', NaN, ...
    'E_bound', NaN, 'Pol_bound', NaN, 'loaded_prePAS', NaN, 'sites_prePAS', NaN, ...
    'distances_bp', [], 'cdf', [], 'avg_E', [], 'avg_Ser2P', [], 'pas_node', NaN, ...
    'error', '');
try
    if nargin < 3 || isempty(free_pools)
        [R, RHE, P_out, d] = run_full_termination_simulation(P, M);
    else
        [R, RHE, P_out, d] = run_full_termination_simulation(P, M, 'FreePools', free_pools);
    end
    [cdf, distances, cad, dg] = calculate_full_pas_cleavage_profile(R, RHE, P_out);
    n = P_out.PAS-1;                  % last node before the PAS
    e = d.R_states(:, 2);
    p = d.R_states(:, 1);
    w = sum(d.R_micro(n, :));
    s.ok = true;
    s.CAD50 = cad;
    s.max_exit_cdf = dg.max_exit_cdf;
    s.E_free = P_out.Ef_ss;
    s.E_bound = d.E_bound;
    s.Pol_bound = d.Pol_bound;
    s.loaded_prePAS = sum(d.R_micro(n, e >= 1))/w;
    s.sites_prePAS = (d.R_micro(n, :)*p)/w;
    s.distances_bp = [0; distances(:)];
    s.cdf = [0; cdf(:)];
    s.avg_E = d.avg_E_bound;
    s.avg_Ser2P = d.avg_Ser2P;
    s.pas_node = P_out.PAS;
catch ME
    s.error = sprintf('%s: %s', ME.identifier, ME.message);
end
end

function c = cad_or_inf(P, M, free_pools)
% CAD50 for root finding: Inf when unreached in the window, NaN on failure.
s = solve_case(P, M, free_pools);
if ~s.ok
    c = NaN;
elseif isnan(s.CAD50)
    c = Inf;
else
    c = s.CAD50;
end
end

function [kEoff_match, status] = match_single_site(P, target, free_pools)
% kEoff at which M = 1 reaches the target CAD50, by bisection in log10 kEoff
% down to 1e5-fold below P.kEoff. CAD50 falls as kEoff falls.
kEoff_match = NaN;
f = @(log_k) cad_or_inf(setfield(P, 'kEoff', 10^log_k), 1, free_pools)-target;
hi = log10(P.kEoff);
lo = hi-5;
f_hi = f(hi);
if isnan(f_hi)
    status = 'error';
    return;
end
if f_hi <= 0
    kEoff_match = P.kEoff;
    status = 'already_fast';
    return;
end
f_lo = f(lo);
if isnan(f_lo)
    status = 'error';
    return;
end
if f_lo > 0
    status = 'unreachable';
    return;
end
for iteration = 1:16
    mid = (lo+hi)/2;
    f_mid = f(mid);
    if isnan(f_mid)
        status = 'error';
        return;
    end
    if f_mid > 0
        hi = mid;
    else
        lo = mid;
    end
end
kEoff_match = 10^((lo+hi)/2);
status = 'matched';
end

%% ------------------------------------------------------------------------
function T = summary_table(S, E_totals)
n = numel(S);
CAD50 = NaN(n, 1); max_exit_cdf = NaN(n, 1); E_free = NaN(n, 1);
E_bound = NaN(n, 1); free_E_fraction = NaN(n, 1); loaded_prePAS = NaN(n, 1);
sites_prePAS = NaN(n, 1); ok = false(n, 1); error = cell(n, 1);
for k = 1:n
    s = S{k};
    CAD50(k) = s.CAD50; max_exit_cdf(k) = s.max_exit_cdf; E_free(k) = s.E_free;
    E_bound(k) = s.E_bound; free_E_fraction(k) = s.E_free/E_totals(k);
    loaded_prePAS(k) = s.loaded_prePAS; sites_prePAS(k) = s.sites_prePAS;
    ok(k) = s.ok; error{k} = s.error;
end
T = table(CAD50, max_exit_cdf, E_free, E_bound, free_E_fraction, ...
    loaded_prePAS, sites_prePAS, ok, error);
end

function write_summary(out_dir, name, key_name, keys, S, E_total)
T = summary_table(S, repmat(E_total, numel(S), 1));
T = [table(keys, 'VariableNames', {key_name}), T];
writetable(T, fullfile(out_dir, name));
end

function write_cdfs(out_dir, name, key_name, keys, S)
key = []; distance_bp = []; cdf = [];
for k = 1:numel(S)
    if S{k}.ok
        key = [key; repmat(keys(k), numel(S{k}.cdf), 1)]; %#ok<AGROW>
        distance_bp = [distance_bp; S{k}.distances_bp]; %#ok<AGROW>
        cdf = [cdf; S{k}.cdf]; %#ok<AGROW>
    end
end
writetable(table(key, distance_bp, cdf, 'VariableNames', {key_name, 'distance_bp', 'cdf'}), ...
    fullfile(out_dir, name));
end

function write_profiles(out_dir, name, key_name, keys, S, L_a)
key = keys([]); node = []; avg_E = []; avg_Ser2P = []; pas_node = [];
for k = 1:numel(S)
    if S{k}.ok
        n = numel(S{k}.avg_E);
        key = [key; repmat(keys(k), n, 1)]; %#ok<AGROW>
        node = [node; (1:n)']; %#ok<AGROW>
        avg_E = [avg_E; S{k}.avg_E]; %#ok<AGROW>
        avg_Ser2P = [avg_Ser2P; S{k}.avg_Ser2P]; %#ok<AGROW>
        pas_node = [pas_node; repmat(S{k}.pas_node, n, 1)]; %#ok<AGROW>
    end
end
T = table(key, node, node*L_a, pas_node, avg_E, avg_Ser2P, ...
    'VariableNames', {key_name, 'node', 'position_bp', 'pas_node', 'avg_E', 'avg_Ser2P'});
writetable(T, fullfile(out_dir, name));
end
