% GeneLengthAnalyze.m — full finite-rate gene-length/resource competition.
% Run GeneLengthGenerateGrid and GeneLengthBuildInterpolation first.
% Each gene uses the same solved free pools; local E/Pol totals are not solved
% again. CAD is measured only within the simulated downstream window.

fprintf('=== Full-Model Gene Length TCD Analysis ===\n');
interp_dir = cpad_analysis_output_dir('GeneLengthAnalysis', 'SecondVersionResults');
interp_files = dir(fullfile(interp_dir, 'full_gene_length_interpolation_*.mat'));
if isempty(interp_files)
    error('GeneLengthAnalyze:MissingFullInterpolation', ...
        'Run GeneLengthGenerateGrid.m and GeneLengthBuildInterpolation.m to create full-model data.');
end
[~, newest_idx] = max([interp_files.datenum]);
interp_filename = fullfile(interp_dir, interp_files(newest_idx).name);
loaded = load(interp_filename, 'interpolation_results');
interpolation_results = loaded.interpolation_results;
if ~isfield(interpolation_results.metadata, 'model_variant') || ...
        ~strcmp(interpolation_results.metadata.model_variant, 'full_kinetics_rapid_EH_disassembly') || ...
        ~isfield(interpolation_results.metadata, 'pool_mode') || ...
        ~strcmp(interpolation_results.metadata.pool_mode, 'fixed_free_pools')
    error('GeneLengthAnalyze:IncompatibleInterpolation', ...
        'Regenerate the grid and interpolation with the full finite-rate model.');
end
grid_parameters = interpolation_results.original_grid;
base_params = grid_parameters.base_parameters;
R_total_target = grid_parameters.R_total_base;
E_total_target = grid_parameters.E_total_base;
num_active_genes = grid_parameters.num_active_genes;
after_PAS_length = grid_parameters.after_PAS_length;

% Normalize the length distribution on the simulated range only. This is a
% truncated gene population; its probability mass is saved with the result.
L_min_int = grid_parameters.L_range(1);
L_max_int = grid_parameters.L_range(2);
L_integration_points = 200;
L_integration = logspace(log10(L_min_int), log10(L_max_int), L_integration_points);
L_integration([1 end]) = [L_min_int L_max_int];
spacing = diff(L_integration);
dL = [spacing(1)/2, (spacing(1:end-1)+spacing(2:end))/2, spacing(end)/2];
f_L = interpolation_results.functions.gene_length_pdf(L_integration);
length_probability_mass = sum(f_L.*dL);
weights = f_L.*dL/length_probability_mass;
fprintf('Population: %d genes, TSS-to-PAS lengths %.1f–%.1f kb (PDF mass %.3f)\n', ...
    num_active_genes, L_min_int/1000, L_max_int/1000, length_probability_mass);

% Full-model occupancy is linear in free Pol II at fixed free E. Solve E
% conservation on its physical interval and recover R_free analytically.
pool_timer = tic;
pool_solution = solve_full_genome_pools(interpolation_results, L_integration, ...
    weights, R_total_target, E_total_target, num_active_genes);
solve_time = toc(pool_timer);
R_free_solution = pool_solution.R_free;
E_free_solution = pool_solution.E_free;
fprintf('Shared free pools: Pol II = %.6g, E = %.6g\n', R_free_solution, E_free_solution);
fprintf('Conservation residuals: R = %.3g, E = %.3g\n', pool_solution.conservation_residuals);

% Evaluate new full-model profiles at the shared pools and the same geometry
% convention as the lookup grid. Keep diagnostics for unreached thresholds.
L_TCD_analysis = logspace(log10(L_min_int), log10(L_max_int), 20);
L_TCD_analysis([1 end]) = [L_min_int L_max_int];
TCD_thresholds = 0.5; % Fractions; for example [0.25 0.5 0.75]
n_L_TCD = numel(L_TCD_analysis);
n_thresholds = numel(TCD_thresholds);
TCD_results = NaN(n_L_TCD, n_thresholds);
max_exit_cdf = NaN(n_L_TCD, 1);
rhs_residuals = NaN(n_L_TCD, 1);
error_messages = repmat({''}, n_L_TCD, 1);
termination_profiles = cell(n_L_TCD, 1);
if isempty(gcp('nocreate')), parpool; end
parfor i = 1:n_L_TCD
    try
        P = base_params;
        P.PASposition = L_TCD_analysis(i);
        P.geneLength_bp = P.PASposition+after_PAS_length;
        [R, RHE, P_sim, details] = run_full_termination_simulation(P, P.EBindingNumber, ...
            'FreePools', [R_free_solution E_free_solution]);
        [cdf, distances, ~, diagnostics] = calculate_full_pas_cleavage_profile( ...
            R, RHE, P_sim, 'PercentCleavage', 0);
        max_exit_cdf(i) = diagnostics.max_exit_cdf;
        rhs_residuals(i) = details.rhs_max_abs;
        termination_profiles{i} = struct('distances', distances, 'profile', cdf);
        local_TCD = NaN(1, n_thresholds);
        for j = 1:n_thresholds
            [~, ~, local_TCD(j)] = calculate_full_pas_cleavage_profile( ...
                R, RHE, P_sim, 'PercentCleavage', 100*TCD_thresholds(j));
        end
        TCD_results(i, :) = local_TCD;
    catch ME
        error_messages{i} = ME.message;
        warning('GeneLengthAnalyze:SimulationFailed', ...
            'Full-model profile failed for L = %.0f bp: %s', L_TCD_analysis(i), ME.message);
    end
end

% NaNs remain in the plotted series, so unreached thresholds leave gaps.
valid_TCD = isfinite(TCD_results);
correlations = NaN(n_thresholds, 1);
fig = figure('Position', [100 100 800 600]);
ax = axes('Parent', fig);
hold(ax, 'on');
colors = lines(n_thresholds);
for j = 1:n_thresholds
    plot(ax, L_TCD_analysis/1000, TCD_results(:, j), 'o-', ...
        'LineWidth', 2, 'MarkerSize', 4, 'Color', colors(j, :), ...
        'DisplayName', sprintf('%.0f%% threshold', 100*TCD_thresholds(j)));
    valid = valid_TCD(:, j);
    fprintf('CAD_%.0f measured within the window for %d/%d lengths.\n', ...
        100*TCD_thresholds(j), nnz(valid), n_L_TCD);
    if any(valid)
        values = TCD_results(valid, j);
        fprintf('  Range %.0f–%.0f bp; median %.0f bp\n', min(values), max(values), median(values));
        if nnz(valid) > 1
            correlation_matrix = corrcoef(log10(L_TCD_analysis(valid)), values);
            correlations(j) = correlation_matrix(1, 2);
        end
    end
end
set(ax, 'XScale', 'log');
xlabel(ax, 'TSS-to-PAS distance (kb)');
ylabel(ax, 'CAD (bp)');
title(ax, 'Full finite-rate termination distance vs gene length');
legend(ax, 'Location', 'best');
grid(ax, 'on');
hold(ax, 'off');

analysis_results = struct();
analysis_results.metadata.creation_date = datestr(now);
analysis_results.metadata.model_variant = interpolation_results.metadata.model_variant;
analysis_results.metadata.pool_mode = 'shared_genome_pools';
analysis_results.metadata.source_interpolation_file = interp_filename;
analysis_results.parameters = grid_parameters;
analysis_results.integration.lengths = L_integration;
analysis_results.integration.weights = weights;
analysis_results.integration.probability_mass_in_grid = length_probability_mass;
analysis_results.solution = pool_solution;
analysis_results.solution.solve_time_seconds = solve_time;
analysis_results.TCD.gene_lengths = L_TCD_analysis;
analysis_results.TCD.thresholds = TCD_thresholds;
analysis_results.TCD.values = TCD_results;
analysis_results.TCD.valid_indices = valid_TCD;
analysis_results.TCD.termination_profiles = termination_profiles;
analysis_results.TCD.max_exit_cdf = max_exit_cdf;
analysis_results.TCD.rhs_residuals = rhs_residuals;
analysis_results.TCD.error_messages = error_messages;
analysis_results.statistics.correlation_with_log_length = correlations;
timestamp = datestr(now, 'yyyymmdd_HHMMSS');
results_filename = fullfile(interp_dir, sprintf('full_gene_length_TCD_analysis_%s.mat', timestamp));
save(results_filename, 'analysis_results', '-v7.3');
plot_filename = fullfile(interp_dir, sprintf('Full_TCD_vs_gene_length_%s.png', timestamp));
saveas(fig, plot_filename);
fprintf('Saved full-model analysis: %s\n', results_filename);

function solution = solve_full_genome_pools(interpolation_results, L_values, weights, R_total, E_total, num_genes)
% SOLVE_FULL_GENOME_POOLS Shared free pools using full-model gene occupancies.
% weights are normalized integration weights for L_values. Grid data must
% come from the fixed-free-pool, finite-rate model and cover [0,R_total],
% [0,E_total] and every requested length. No extrapolation is permitted.
% At fixed E_free, occupancy is exactly proportional to R_free. Eliminate
% R_free analytically and solve only E conservation on [0,E_total].

validateattributes(R_total, {'numeric'}, {'scalar', 'real', 'finite', 'positive'});
validateattributes(E_total, {'numeric'}, {'scalar', 'real', 'finite', 'nonnegative'});
validateattributes(num_genes, {'numeric'}, {'scalar', 'real', 'finite', 'integer', 'positive'});
validateattributes(L_values, {'numeric'}, {'vector', 'nonempty', 'real', 'finite', 'positive'});
validateattributes(weights, {'numeric'}, {'vector', 'numel', numel(L_values), 'real', 'finite', 'nonnegative'});
L_values = L_values(:);
weights = weights(:);
if abs(sum(weights)-1) > 1e-8
    error('solve_full_genome_pools:InvalidWeights', 'Integration weights must sum to one.');
end
metadata = interpolation_results.metadata;
if ~isfield(metadata, 'model_variant') || ...
        ~strcmp(metadata.model_variant, 'full_kinetics_rapid_EH_disassembly') || ...
        ~isfield(metadata, 'pool_mode') || ~strcmp(metadata.pool_mode, 'fixed_free_pools')
    error('solve_full_genome_pools:IncompatibleModel', ...
        'Rebuild the gene-length grid and interpolation using the full finite-rate model.');
end
grid = interpolation_results.original_grid;
if grid.R_free_range(1) > 0 || grid.R_free_range(2) < R_total || ...
        grid.E_free_range(1) > 0 || grid.E_free_range(2) < E_total || ...
        any(L_values < grid.L_range(1) | L_values > grid.L_range(2))
    error('solve_full_genome_pools:OutsideGrid', ...
        'Regenerate a grid covering zero through the target totals and all integration lengths.');
end
R_interp = interpolation_results.functions.R_occupied_interp;
E_interp = interpolation_results.functions.E_occupied_interp;

if E_total == 0
    E_free = 0;
else
    options = optimset('TolX', eps, 'FunValCheck', 'on');
    [E_free, ~, exitflag] = fzero(@evaluate_pools, [0, E_total], options);
    if exitflag <= 0
        error('solve_full_genome_pools:NotConverged', 'Shared free-E conservation solve did not converge.');
    end
end
[E_residual, R_free, mean_R, mean_E] = evaluate_pools(E_free);
residuals = [R_free+num_genes*mean_R-R_total, E_residual];
if any(~isfinite(residuals)) || max(abs(residuals)) > 1e-6 || ...
        R_free < 0 || R_free > R_total || E_free < 0 || E_free > E_total
    error('solve_full_genome_pools:InvalidSolution', ...
        'Shared pool solution failed conservation: R residual %.3g, E residual %.3g.', residuals(1), residuals(2));
end
solution.R_free = R_free;
solution.E_free = E_free;
solution.R_total_target = R_total;
solution.E_total_target = E_total;
solution.num_active_genes = num_genes;
solution.mean_R_occupied = mean_R;
solution.mean_E_occupied = mean_E;
solution.conservation_residuals = residuals;

    function [residual, trial_R, mean_R_at_trial, mean_E_at_trial] = evaluate_pools(trial_E)
        query_shape = ones(size(L_values));
        R_at_reference = R_interp(R_total*query_shape, trial_E*query_shape, L_values);
        E_at_reference = E_interp(R_total*query_shape, trial_E*query_shape, L_values);
        if any(~isfinite(R_at_reference)) || any(~isfinite(E_at_reference)) || ...
                any(R_at_reference < 0) || any(E_at_reference < 0)
            error('solve_full_genome_pools:InvalidInterpolation', ...
                'Occupancy interpolation returned invalid values; rebuild the grid.');
        end
        mean_R_per_free = sum(weights.*R_at_reference(:))/R_total;
        trial_R = R_total/(1+num_genes*mean_R_per_free);
        mean_R_at_trial = trial_R*mean_R_per_free;
        mean_E_at_trial = (trial_R/R_total)*sum(weights.*E_at_reference(:));
        residual = trial_E+num_genes*mean_E_at_trial-E_total;
    end
end
