% GeneLengthBuildInterpolation.m — full finite-rate occupancy grid
% Builds interpolation functions for gene length analysis
%
% This script implements Step 2 of the gene length analysis:
% - Load and validate grid data from GeneLengthGenerateGrid.m
% - Build smooth interpolation functions R_occupied(R_free, E_free, L) and E_occupied(R_free, E_free, L)
% - Define realistic gene length distribution f(L) based on human genomic data
% - Validate interpolation quality and save functions for analysis

fprintf('=== Gene Length Interpolation Builder ===\n');

%% --- LOAD GRID DATA ---

% Find the most recent grid data file
grid_dir = cpad_analysis_output_dir('GeneLengthAnalysis', 'SecondVersionResults');
if ~exist(grid_dir, 'dir')
    error('Grid data directory not found. Please run GeneLengthGenerateGrid.m first.');
end

% Look for .mat files
mat_files = dir(fullfile(grid_dir, 'full_gene_length_grid_data_*.mat'));
if isempty(mat_files)
    error('No grid data files found. Please run GeneLengthGenerateGrid.m first.');
end

% Use the most recent file
[~, newest_idx] = max([mat_files.datenum]);
grid_filename = fullfile(grid_dir, mat_files(newest_idx).name);
fprintf('Loading: %s\n', grid_filename);

load(grid_filename, 'results');

% Never mix old equilibrium grids with full finite-rate profiles.
if ~isfield(results.metadata, 'model_variant') || ...
        ~strcmp(results.metadata.model_variant, 'full_kinetics_rapid_EH_disassembly') || ...
        ~isfield(results.metadata, 'pool_mode') || ...
        ~strcmp(results.metadata.pool_mode, 'fixed_free_pools') || ...
        ~isfield(results.metadata, 'grid_layout') || ~strcmp(results.metadata.grid_layout, 'ndgrid')
    error('GeneLengthBuildInterpolation:IncompatibleGrid', ...
        'Regenerate full-model grid data with GeneLengthGenerateGrid.m.');
end

%% --- DATA VALIDATION AND CLEANING ---

% Extract data
R_free_data = results.data.R_free_vec(:);
E_free_data = results.data.E_free_vec(:);
L_data = results.data.L_vec(:);
[R_occupied_data, E_occupied_data, grid_validation] = validate_full_gene_length_grid( ...
    results.data, results.parameters.base_parameters.EBindingNumber);

% Keep every Cartesian grid row. Only solver-level roundoff is corrected;
% the validator rejects failed, nonfinite or materially negative results.
n_total = numel(results.data.success_flag);
valid_indices = true(n_total, 1);
n_valid = n_total;

fprintf('Valid data points: %d/%d (%.1f%%)\n', n_valid, n_total, n_valid/n_total*100);

if grid_validation.corrected_rows > 0
    fprintf('Corrected numerical roundoff to zero in %d rows (Pol: %d, E: %d values).\n', ...
        grid_validation.corrected_rows, grid_validation.corrected_R_values, grid_validation.corrected_E_values);
end

% Clean data
R_free_clean = R_free_data(valid_indices);
E_free_clean = E_free_data(valid_indices);
L_clean = L_data(valid_indices);
R_occupied_clean = R_occupied_data(valid_indices);
E_occupied_clean = E_occupied_data(valid_indices);

%% --- BUILD INTERPOLATION FUNCTIONS ---

% Preserve the Cartesian grid and exact linear dependence on R_free.
% Extrapolation is disabled; downstream pool and length queries stay in range.
R_axis = results.grid.R_free_values;
E_axis = results.grid.E_free_values;
L_axis = results.grid.L_values;
grid_size = [numel(R_axis), numel(E_axis), numel(L_axis)];
R_occupied_interp = griddedInterpolant({R_axis, E_axis, L_axis}, ...
    reshape(R_occupied_data, grid_size), 'linear', 'none');
E_occupied_interp = griddedInterpolant({R_axis, E_axis, L_axis}, ...
    reshape(E_occupied_data, grid_size), 'linear', 'none');

%% --- INTERPOLATION VALIDATION ---

% Test interpolation on a subset of original data
n_test = min(1000, max(1, floor(n_valid * 0.1)));  % Test on 10% of data or 1000 points, whichever is smaller
test_indices = randperm(n_valid, n_test);

R_free_test = R_free_clean(test_indices);
E_free_test = E_free_clean(test_indices);
L_test = L_clean(test_indices);
R_occupied_true = R_occupied_clean(test_indices);
E_occupied_true = E_occupied_clean(test_indices);

% Interpolated values
R_occupied_interp_test = R_occupied_interp(R_free_test, E_free_test, L_test);
E_occupied_interp_test = E_occupied_interp(R_free_test, E_free_test, L_test);

% Calculate errors
R_error = abs(R_occupied_interp_test - R_occupied_true) ./ max(abs(R_occupied_true), eps) * 100;
E_error = abs(E_occupied_interp_test - E_occupied_true) ./ max(abs(E_occupied_true), eps) * 100;

fprintf('Interpolation reconstruction check at grid nodes (not an accuracy estimate):\n');
fprintf('  R_occupied - Mean error: %.2f%%, Max error: %.2f%%\n', mean(R_error), max(R_error));
fprintf('  E_occupied - Mean error: %.2f%%, Max error: %.2f%%\n', mean(E_error), max(E_error));

if mean(R_error) > 10 || mean(E_error) > 10
    warning('High interpolation errors detected. Consider using more grid points or different interpolation method.');
end

%% --- DEFINE GENE LENGTH DISTRIBUTION ---

% Based on the histogram provided:
% - Log-normal distribution with median ~22 kb
% - 25th percentile ~7-8 kb, 75th percentile ~64-65 kb
% - Standard deviation ~0.25-0.3 Mb in linear terms
% - σ ≈ 0.68 in log10 space

% Gene length distribution parameters (log-normal)
median_length = 22000;  % 22 kb median
log_sigma = 0.68;       % Standard deviation in log10 space

% Convert to natural log parameters for lognormal distribution
mu_ln = log(median_length);  % Mean of underlying normal distribution
sigma_ln = log_sigma * log(10);  % Standard deviation of underlying normal distribution

% Create gene length distribution function
gene_length_pdf = @(L) exp(-0.5*((log(L)-mu_ln)/sigma_ln).^2)./(L*sigma_ln*sqrt(2*pi));

% Validate distribution by computing percentiles
% Analytic log-normal percentiles; no uniform-bin approximation on a log grid.
percentile_25 = exp(mu_ln + sigma_ln*sqrt(2)*erfinv(2*0.25-1));
percentile_50 = median_length;
percentile_75 = exp(mu_ln + sigma_ln*sqrt(2)*erfinv(2*0.75-1));

%% --- SAVE INTERPOLATION RESULTS ---

% Create output structure
interpolation_results = struct();
interpolation_results.metadata.creation_date = datestr(now);
interpolation_results.metadata.model_variant = results.metadata.model_variant;
interpolation_results.metadata.pool_mode = results.metadata.pool_mode;
interpolation_results.metadata.source_grid_file = grid_filename;
interpolation_results.metadata.n_valid_points = n_valid;
interpolation_results.metadata.interpolation_method = 'griddedInterpolant, linear, no extrapolation';
interpolation_results.metadata.description = 'Interpolation functions and gene length distribution for gene length analysis';

% Validation results
interpolation_results.validation.R_occupied_mean_error = mean(R_error);
interpolation_results.validation.R_occupied_max_error = max(R_error);
interpolation_results.validation.E_occupied_mean_error = mean(E_error);
interpolation_results.validation.E_occupied_max_error = max(E_error);
interpolation_results.validation.grid_roundoff = grid_validation;

% Gene length distribution parameters
interpolation_results.gene_length_distribution.type = 'log-normal';
interpolation_results.gene_length_distribution.median_bp = median_length;
interpolation_results.gene_length_distribution.log_sigma = log_sigma;
interpolation_results.gene_length_distribution.mu_ln = mu_ln;
interpolation_results.gene_length_distribution.sigma_ln = sigma_ln;
interpolation_results.gene_length_distribution.percentile_25 = percentile_25;
interpolation_results.gene_length_distribution.percentile_50 = percentile_50;
interpolation_results.gene_length_distribution.percentile_75 = percentile_75;

% Store functions (note: these are function handles, may need special handling)
interpolation_results.functions.R_occupied_interp = R_occupied_interp;
interpolation_results.functions.E_occupied_interp = E_occupied_interp;
interpolation_results.functions.gene_length_pdf = gene_length_pdf;

% Original grid information for reference
interpolation_results.original_grid = results.parameters;

% Save results
timestamp = datestr(now, 'yyyymmdd_HHMMSS');
output_filename = fullfile(grid_dir, sprintf('full_gene_length_interpolation_%s.mat', timestamp));
save(output_filename, 'interpolation_results', '-v7.3');

fprintf('Interpolation results saved to: %s\n', output_filename);

%% --- GENERATE SUMMARY PLOTS ---

% Plot 1: Gene length vs R_occupied and E_occupied (line plots with multiple R_free, E_free combinations)
figure('Position', [100, 100, 1000, 400]);

% Create a range of gene lengths for plotting
L_plot_range = logspace(log10(min(L_clean)), log10(max(L_clean)), 100);

% Define multiple R_free and E_free values across their ranges
n_lines = 4;  % Number of different (R_free, E_free) combinations to plot
R_free_values = linspace(min(R_free_clean), max(R_free_clean), n_lines);
E_free_values = linspace(min(E_free_clean), max(E_free_clean), n_lines);

% Create color map for different lines
colors = lines(n_lines);

% Plot R_occupied
subplot(1, 2, 1);
hold on;
for line_idx = 1:n_lines
    R_free_val = R_free_values(line_idx);
    E_free_val = E_free_values(line_idx);
    
    % Calculate R_occupied for this R_free, E_free combination
    R_occupied_line = zeros(size(L_plot_range));
    for i = 1:length(L_plot_range)
        R_occupied_line(i) = R_occupied_interp(R_free_val, E_free_val, L_plot_range(i));
    end
    
    semilogx(L_plot_range/1000, R_occupied_line, '-', 'LineWidth', 2, 'Color', colors(line_idx, :), ...
        'DisplayName', sprintf('R_{free}=%.0f, E_{free}=%.0f', R_free_val, E_free_val));
end
xlabel('TSS-to-PAS Distance (kb)', 'FontSize', 12);
ylabel('R_{occupied}', 'FontSize', 12);
title('R_{occupied} vs Gene Length', 'FontSize', 12, 'FontWeight', 'bold');
legend('Location', 'best', 'FontSize', 9);
grid on;
hold off;

% Plot E_occupied
subplot(1, 2, 2);
hold on;
for line_idx = 1:n_lines
    R_free_val = R_free_values(line_idx);
    E_free_val = E_free_values(line_idx);
    
    % Calculate E_occupied for this R_free, E_free combination
    E_occupied_line = zeros(size(L_plot_range));
    for i = 1:length(L_plot_range)
        E_occupied_line(i) = E_occupied_interp(R_free_val, E_free_val, L_plot_range(i));
    end
    
    semilogx(L_plot_range/1000, E_occupied_line, '-', 'LineWidth', 2, 'Color', colors(line_idx, :), ...
        'DisplayName', sprintf('R_{free}=%.0f, E_{free}=%.0f', R_free_val, E_free_val));
end
xlabel('TSS-to-PAS Distance (kb)', 'FontSize', 12);
ylabel('E_{occupied}', 'FontSize', 12);
title('E_{occupied} vs Gene Length', 'FontSize', 12, 'FontWeight', 'bold');
legend('Location', 'best', 'FontSize', 9);
grid on;
hold off;

% Save plot
line_plot_filename = fullfile(grid_dir, sprintf('gene_length_vs_occupied_%s.png', timestamp));
saveas(gcf, line_plot_filename);
fprintf('Gene length vs occupied resources plot saved to: %s\n', line_plot_filename);

% Plot 2: Sample interpolation surface (R_occupied)
figure('Position', [200, 200, 1000, 400]);

% Plot at median gene length
L_sample = 20000;
R_free_range = linspace(min(R_free_clean), max(R_free_clean), 50);
E_free_range = linspace(min(E_free_clean), max(E_free_clean), 50);
[R_grid_plot, E_grid_plot] = meshgrid(R_free_range, E_free_range);
L_grid_plot = L_sample * ones(size(R_grid_plot));

R_occupied_surface = R_occupied_interp(R_grid_plot, E_grid_plot, L_grid_plot);

subplot(1, 2, 1);
surf(R_grid_plot/1000, E_grid_plot/1000, R_occupied_surface);
xlabel('R_{free} (thousands)', 'FontSize', 10);
ylabel('E_{free} (thousands)', 'FontSize', 10);
zlabel('R_{occupied}', 'FontSize', 10);
title(sprintf('R_{occupied} at L = %.0f kb', L_sample/1000), 'FontSize', 12);
colorbar;
shading interp;

% Plot E_occupied surface
E_occupied_surface = E_occupied_interp(R_grid_plot, E_grid_plot, L_grid_plot);

subplot(1, 2, 2);
surf(R_grid_plot/1000, E_grid_plot/1000, E_occupied_surface);
xlabel('R_{free} (thousands)', 'FontSize', 10);
ylabel('E_{free} (thousands)', 'FontSize', 10);
zlabel('E_{occupied}', 'FontSize', 10);
title(sprintf('E_{occupied} at L = %.0f kb', L_sample/1000), 'FontSize', 12);
colorbar;
shading interp;

% Save plot
surface_plot_filename = fullfile(grid_dir, sprintf('interpolation_surfaces_%s.png', timestamp));
saveas(gcf, surface_plot_filename);
fprintf('Interpolation surfaces plot saved to: %s\n', surface_plot_filename);

%% --- HELPER FUNCTIONS ---

function [R_occupied, E_occupied, diagnostics] = validate_full_gene_length_grid(data, EBindingNumber)
% VALIDATE_FULL_GENE_LENGTH_GRID Accept solver roundoff, reject invalid rows.
% Sparse steady solves can leave tiny signed populations at exact-zero
% boundaries. The solver accepts roundoff down to -1e-9; an exact >= 0
% occupancy check must not subsequently discard those successful rows.
% Only roundoff is corrected. Failed/nonfinite/materially negative rows
% still prevent construction of a complete Cartesian interpolant.

validateattributes(EBindingNumber, {'numeric'}, ...
    {'scalar', 'real', 'finite', 'integer', 'positive'});
fields = {'R_free_vec', 'E_free_vec', 'L_vec', ...
          'R_occupied_vec', 'E_occupied_vec', 'success_flag'};
n = numel(data.success_flag);
for k = 1:numel(fields)
    validateattributes(data.(fields{k}), {'numeric'}, ...
        {'vector', 'real', 'nonempty', 'numel', n}, mfilename, fields{k});
    data.(fields{k}) = data.(fields{k})(:);
end
R_occupied = data.R_occupied_vec;
E_occupied = data.E_occupied_vec;

% Absolute solver tolerance plus a small floating-point allowance for sums
% of large populations. E sums are bounded in scale by M times Pol occupancy.
R_scale = max(ones(n, 1), abs(R_occupied));
E_scale = max([ones(n, 1), EBindingNumber*abs(R_occupied), abs(E_occupied)], [], 2);
R_tolerance = 1e-9 + 64*eps(R_scale);
E_tolerance = 1e-9 + 64*eps(E_scale);
failed = data.success_flag ~= 1;
nonfinite = ~isfinite(R_occupied) | ~isfinite(E_occupied);
negative = R_occupied < -R_tolerance | E_occupied < -E_tolerance;
bad_coordinates = ~isfinite(data.R_free_vec) | data.R_free_vec < 0 ...
    | ~isfinite(data.E_free_vec) | data.E_free_vec < 0 ...
    | ~isfinite(data.L_vec) | data.L_vec <= 0;
zero_pol = data.R_free_vec == 0;
zero_e = zero_pol | data.E_free_vec == 0;
bad_boundary = (zero_pol & abs(R_occupied) > R_tolerance) ...
    | (zero_e & abs(E_occupied) > E_tolerance);
invalid = failed | nonfinite | negative | bad_coordinates | bad_boundary;
if any(invalid)
    first = find(invalid, 1);
    solver_message = '';
    if isfield(data, 'error_messages') && numel(data.error_messages) >= first
        solver_message = data.error_messages{first};
    end
    error('GeneLengthBuildInterpolation:IncompleteGrid', ...
        ['Full-model interpolation requires a complete grid. Invalid rows: %d/%d ' ...
         '(failed=%d, nonfinite=%d, materially negative=%d, bad coordinates=%d, ' ...
         'inconsistent zero-pool boundary=%d; categories may overlap).\n' ...
         'First invalid row %d: R_free=%.16g, E_free=%.16g, L=%.16g, ' ...
         'R_occupied=%.16g, E_occupied=%.16g; tolerances R=%.3g, E=%.3g.\n' ...
         'Solver message: %s\nInspect these rows in results.data and regenerate if needed.'], ...
        nnz(invalid), n, nnz(failed), nnz(nonfinite), nnz(negative), ...
        nnz(bad_coordinates), nnz(bad_boundary), first, ...
        data.R_free_vec(first), data.E_free_vec(first), data.L_vec(first), ...
        R_occupied(first), E_occupied(first), R_tolerance(first), E_tolerance(first), solver_message);
end

correct_R = R_occupied < 0 | (zero_pol & R_occupied ~= 0);
correct_E = E_occupied < 0 | (zero_e & E_occupied ~= 0);
diagnostics.corrected_rows = nnz(correct_R | correct_E);
diagnostics.corrected_R_values = nnz(correct_R);
diagnostics.corrected_E_values = nnz(correct_E);
diagnostics.minimum_raw_R_occupied = min(R_occupied);
diagnostics.minimum_raw_E_occupied = min(E_occupied);
diagnostics.max_R_tolerance = max(R_tolerance);
diagnostics.max_E_tolerance = max(E_tolerance);
diagnostics.tolerance_rule = '1e-9 + 64*eps(population scale); E scale includes M*Pol';
R_occupied(correct_R) = 0;
E_occupied(correct_E) = 0;
end
