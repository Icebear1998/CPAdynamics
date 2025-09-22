% BUILD_GENE_LENGTH_INTERPOLATION.m
% Builds interpolation functions for gene length analysis
%
% This script implements Step 2 of the gene length analysis:
% - Load and validate grid data from generate_gene_length_grid.m
% - Build smooth interpolation functions R_occupied(R_free, E_free, L) and E_occupied(R_free, E_free, L)
% - Define realistic gene length distribution f(L) based on human genomic data
% - Validate interpolation quality and save functions for analysis

fprintf('=== Gene Length Interpolation Builder ===\n');
fprintf('Building interpolation functions and gene length distribution...\n\n');

%% --- LOAD GRID DATA ---
fprintf('Loading grid data...\n');

% Find the most recent grid data file
grid_dir = 'SecondVersionResults/GeneLengthAnalysis/';
if ~exist(grid_dir, 'dir')
    error('Grid data directory not found. Please run generate_gene_length_grid.m first.');
end

% Look for .mat files
mat_files = dir(fullfile(grid_dir, 'gene_length_grid_data_*.mat'));
if isempty(mat_files)
    error('No grid data files found. Please run generate_gene_length_grid.m first.');
end

% Use the most recent file
[~, newest_idx] = max([mat_files.datenum]);
grid_filename = fullfile(grid_dir, mat_files(newest_idx).name);
fprintf('Loading: %s\n', grid_filename);

load(grid_filename, 'results');
fprintf('Grid data loaded successfully!\n');
fprintf('  Grid points: %d\n', length(results.data.R_free_vec));
fprintf('  Success rate: %.1f%%\n', results.metadata.success_rate);

%% --- DATA VALIDATION AND CLEANING ---
fprintf('\nValidating and cleaning data...\n');

% Extract data
R_free_data = results.data.R_free_vec;
E_free_data = results.data.E_free_vec;
L_data = results.data.L_vec;
R_occupied_data = results.data.R_occupied_vec;
E_occupied_data = results.data.E_occupied_vec;
success_flags = results.data.success_flag;

% Remove failed simulations
valid_indices = (success_flags == 1) & ~isnan(R_occupied_data) & ~isnan(E_occupied_data);
n_valid = sum(valid_indices);
n_total = length(success_flags);

fprintf('Valid data points: %d/%d (%.1f%%)\n', n_valid, n_total, n_valid/n_total*100);

if n_valid < 0.5 * n_total
    warning('Less than 50%% of simulations succeeded. Consider adjusting parameter ranges.');
end

% Clean data
R_free_clean = R_free_data(valid_indices);
E_free_clean = E_free_data(valid_indices);
L_clean = L_data(valid_indices);
R_occupied_clean = R_occupied_data(valid_indices);
E_occupied_clean = E_occupied_data(valid_indices);

% Data statistics
fprintf('\nData statistics (valid points):\n');
fprintf('  R_free: %.0f to %.0f\n', min(R_free_clean), max(R_free_clean));
fprintf('  E_free: %.0f to %.0f\n', min(E_free_clean), max(E_free_clean));
fprintf('  L (TSS-to-PAS): %.0f to %.0f bp\n', min(L_clean), max(L_clean));
fprintf('  R_occupied: %.2f to %.2f (mean: %.2f)\n', min(R_occupied_clean), max(R_occupied_clean), mean(R_occupied_clean));
fprintf('  E_occupied: %.2f to %.2f (mean: %.2f)\n', min(E_occupied_clean), max(E_occupied_clean), mean(E_occupied_clean));

%% --- BUILD INTERPOLATION FUNCTIONS ---
fprintf('\nBuilding interpolation functions...\n');

% Create interpolation functions using scatteredInterpolant
% This handles irregular 3D grids better than griddata
fprintf('  Creating R_occupied interpolant...\n');
R_occupied_interp = scatteredInterpolant(R_free_clean, E_free_clean, L_clean, R_occupied_clean, ...
    'linear', 'nearest');  % Linear interpolation with nearest neighbor extrapolation

fprintf('  Creating E_occupied interpolant...\n');
E_occupied_interp = scatteredInterpolant(R_free_clean, E_free_clean, L_clean, E_occupied_clean, ...
    'linear', 'nearest');  % Linear interpolation with nearest neighbor extrapolation

fprintf('Interpolation functions created successfully!\n');

%% --- INTERPOLATION VALIDATION ---
fprintf('\nValidating interpolation quality...\n');

% Test interpolation on a subset of original data
n_test = min(1000, floor(n_valid * 0.1));  % Test on 10% of data or 1000 points, whichever is smaller
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
R_error = abs(R_occupied_interp_test - R_occupied_true) ./ R_occupied_true * 100;
E_error = abs(E_occupied_interp_test - E_occupied_true) ./ E_occupied_true * 100;

fprintf('Interpolation validation results:\n');
fprintf('  R_occupied - Mean error: %.2f%%, Max error: %.2f%%\n', mean(R_error), max(R_error));
fprintf('  E_occupied - Mean error: %.2f%%, Max error: %.2f%%\n', mean(E_error), max(E_error));

if mean(R_error) > 10 || mean(E_error) > 10
    warning('High interpolation errors detected. Consider using more grid points or different interpolation method.');
end

%% --- DEFINE GENE LENGTH DISTRIBUTION ---
fprintf('\nDefining gene length distribution f(L)...\n');

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

fprintf('Gene length distribution parameters:\n');
fprintf('  Distribution: Log-normal\n');
fprintf('  Median: %.1f kb\n', median_length/1000);
fprintf('  Log10 sigma: %.2f\n', log_sigma);
fprintf('  Natural log mu: %.3f\n', mu_ln);
fprintf('  Natural log sigma: %.3f\n', sigma_ln);

% Create gene length distribution function
gene_length_pdf = @(L) lognpdf(L, mu_ln, sigma_ln);

% Validate distribution by computing percentiles
L_test_range = logspace(3, 6, 10000);  % 1 kb to 1 Mb
pdf_values = gene_length_pdf(L_test_range);
cdf_values = cumsum(pdf_values) * (L_test_range(2) - L_test_range(1));
cdf_values = cdf_values / cdf_values(end);  % Normalize

% Find percentiles
percentile_25 = interp1(cdf_values, L_test_range, 0.25);
percentile_50 = interp1(cdf_values, L_test_range, 0.50);
percentile_75 = interp1(cdf_values, L_test_range, 0.75);

fprintf('Distribution validation:\n');
fprintf('  25th percentile: %.1f kb (expected: ~7-8 kb)\n', percentile_25/1000);
fprintf('  50th percentile: %.1f kb (expected: ~22 kb)\n', percentile_50/1000);
fprintf('  75th percentile: %.1f kb (expected: ~64-65 kb)\n', percentile_75/1000);

%% --- CREATE VISUALIZATION FUNCTIONS ---
fprintf('\nCreating visualization functions...\n');

% Function to plot interpolation surfaces
plot_interpolation_surfaces = @() create_interpolation_plots(R_occupied_interp, E_occupied_interp, results);

% Function to plot gene length distribution
plot_gene_length_distribution = @() create_distribution_plot(gene_length_pdf, L_test_range);

fprintf('Visualization functions created.\n');

%% --- SAVE INTERPOLATION RESULTS ---
fprintf('\nSaving interpolation results...\n');

% Create output structure
interpolation_results = struct();
interpolation_results.metadata.creation_date = datestr(now);
interpolation_results.metadata.source_grid_file = grid_filename;
interpolation_results.metadata.n_valid_points = n_valid;
interpolation_results.metadata.interpolation_method = 'scatteredInterpolant with linear interpolation';
interpolation_results.metadata.description = 'Interpolation functions and gene length distribution for gene length analysis';

% Validation results
interpolation_results.validation.R_occupied_mean_error = mean(R_error);
interpolation_results.validation.R_occupied_max_error = max(R_error);
interpolation_results.validation.E_occupied_mean_error = mean(E_error);
interpolation_results.validation.E_occupied_max_error = max(E_error);

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
interpolation_results.functions.plot_interpolation_surfaces = plot_interpolation_surfaces;
interpolation_results.functions.plot_gene_length_distribution = plot_gene_length_distribution;

% Original grid information for reference
interpolation_results.original_grid = results.parameters;

% Save results
timestamp = datestr(now, 'yyyymmdd_HHMMSS');
output_filename = fullfile(grid_dir, sprintf('gene_length_interpolation_%s.mat', timestamp));
save(output_filename, 'interpolation_results', '-v7.3');

fprintf('Interpolation results saved to: %s\n', output_filename);

%% --- GENERATE SUMMARY PLOTS ---
fprintf('\nGenerating summary plots...\n');

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
L_sample = 3000;
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

%% --- SUMMARY ---
fprintf('\n=== INTERPOLATION BUILDING COMPLETE ===\n');
fprintf('Results summary:\n');
fprintf('  Valid data points used: %d\n', n_valid);
fprintf('  Interpolation mean errors: R=%.2f%%, E=%.2f%%\n', mean(R_error), mean(E_error));
fprintf('  Gene length distribution: Log-normal (median=%.1f kb)\n', median_length/1000);
fprintf('  Output file: %s\n', output_filename);
fprintf('  Generated plots: gene length vs occupied resources (multi-line), interpolation surfaces\n');

fprintf('\nNext steps:\n');
fprintf('1. Run analyze_gene_length_TCD.m to perform the full analysis\n');
fprintf('2. Use the interpolation functions to solve conservation equations\n');
fprintf('3. Calculate TCD relationships across gene lengths\n');

fprintf('\nInterpolation building completed successfully!\n');

%% --- HELPER FUNCTIONS ---

function create_interpolation_plots(R_interp, E_interp, grid_results)
    % Create comprehensive interpolation validation plots
    figure('Position', [100, 100, 1200, 800]);
    
    % Plot original data distribution
    subplot(2, 3, 1);
    valid_idx = grid_results.data.success_flag == 1;
    scatter3(grid_results.data.R_free_vec(valid_idx)/1000, ...
             grid_results.data.E_free_vec(valid_idx)/1000, ...
             grid_results.data.L_vec(valid_idx)/1000, 10, ...
             grid_results.data.R_occupied_vec(valid_idx), 'filled');
    xlabel('R_{free} (thousands)'); ylabel('E_{free} (thousands)'); zlabel('L (kb)');
    title('Original Data: R_{occupied}');
    colorbar;
    
    % Additional validation plots can be added here
    % This is a placeholder for more comprehensive validation
end

function create_distribution_plot(pdf_func, L_range)
    % Create detailed gene length distribution plot
    figure('Position', [300, 300, 800, 600]);
    
    pdf_vals = pdf_func(L_range);
    semilogx(L_range/1000, pdf_vals*1000, 'b-', 'LineWidth', 2);
    xlabel('Gene Length (kb)');
    ylabel('Probability Density (per kb)');
    title('Gene Length Distribution');
    grid on;
end
