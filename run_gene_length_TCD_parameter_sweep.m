% RUN_GENE_LENGTH_TCD_PARAMETER_SWEEP.m
% Full pipeline for gene length TCD analysis with parameter sweeping
%
% This script runs the complete gene length analysis pipeline:
%   1. generate_gene_length_grid.m - Generate lookup data
%   2. build_gene_length_interpolation.m - Build interpolation functions
%   3. analyze_gene_length_TCD.m - Calculate TCD vs gene length
%
% Parameter sweeping options:
%   - kPon_max: Maximum promoter binding rate
%   - SD: Saturation distance (bp)
%
% The script allows you to choose which parameter to sweep and generates
% a comprehensive plot showing TCD vs gene length for all parameter values.

fprintf('=== Gene Length TCD Analysis - Parameter Sweep Pipeline ===\n\n');

%% --- USER CONFIGURATION ---

% Choose which parameter to sweep: 'kPon_max' or 'SD'
sweep_parameter = 'kPon_max';  % Change to 'SD' to sweep saturation distance

% Parameter sweep values
kPon_max_values = [0.2, 0.5, 1, 5, 10];
SD_values = [4000, 10000, 20000, 40000, 70000, 100000];

% Fixed values when not sweeping
kPon_max_default = 1;      % Default kPon_max
SD_default = 20000;        % Default SD (bp)

% Analysis options
run_parallel = true;       % Use parallel computation
TCD_threshold = 0.5;       % Which TCD threshold to plot (0.25, 0.5, or 0.75)

%% --- SETUP ---

% Select sweep values based on parameter choice
if strcmp(sweep_parameter, 'kPon_max')
    sweep_values = kPon_max_values;
    sweep_name = 'kPon_{max}';
    fixed_param_name = sprintf('SD=%dkb', SD_default/1000);
elseif strcmp(sweep_parameter, 'SD')
    sweep_values = SD_values;
    sweep_name = 'SD (kb)';
    fixed_param_name = sprintf('kPon_{max}=%.1f', kPon_max_default);
else
    error('Invalid sweep_parameter. Choose ''kPon_max'' or ''SD''.');
end

n_sweep = length(sweep_values);

fprintf('Parameter sweep configuration:\n');
fprintf('  Sweeping: %s\n', sweep_parameter);
fprintf('  Number of values: %d\n', n_sweep);
fprintf('  Values: %s\n', mat2str(sweep_values));
fprintf('  Fixed parameter: %s\n', fixed_param_name);
fprintf('  TCD threshold for plotting: %.0f%%\n\n', TCD_threshold*100);

% Create output directory for sweep results
output_dir = sprintf('SecondVersionResults/GeneLengthAnalysis/ParameterSweep_%s_%s', ...
    sweep_parameter, datestr(now, 'yyyymmdd_HHMMSS'));
if ~exist(output_dir, 'dir')
    mkdir(output_dir);
end

fprintf('Output directory: %s\n\n', output_dir);

% Start parallel pool if requested
if run_parallel && isempty(gcp('nocreate'))
    parpool;
end

%% --- PARAMETER SWEEP LOOP ---

% Pre-allocate results storage
sweep_results = cell(n_sweep, 1);
sweep_success = false(n_sweep, 1);

fprintf('Starting parameter sweep...\n');
fprintf('========================================\n\n');

for i = 1:n_sweep
    sweep_value = sweep_values(i);
    
    % Set parameters for this sweep iteration
    if strcmp(sweep_parameter, 'kPon_max')
        kPon_max_current = sweep_value;
        SD_current = SD_default;
        param_label = sprintf('kPon_max=%.1f', kPon_max_current);
    else
        kPon_max_current = kPon_max_default;
        SD_current = sweep_value;
        param_label = sprintf('SD=%dkb', SD_current/1000);
    end
    
    fprintf('--- Sweep %d/%d: %s ---\n', i, n_sweep, param_label);
    
    try
        %% Step 1: Generate Grid Data
        fprintf('Step 1/3: Generating grid data...\n');
        tic;
        [grid_results, grid_file] = run_grid_generation(kPon_max_current, SD_current, output_dir, param_label);
        fprintf('  Grid generation completed in %.1f minutes\n', toc/60);
        
        %% Step 2: Build Interpolation
        fprintf('Step 2/3: Building interpolation functions...\n');
        tic;
        [interp_results, interp_file] = run_interpolation_build(grid_file, output_dir, param_label);
        fprintf('  Interpolation completed in %.1f seconds\n', toc);
        
        %% Step 3: Analyze TCD
        fprintf('Step 3/3: Analyzing TCD vs gene length...\n');
        tic;
        [analysis_results, analysis_file] = run_TCD_analysis(interp_file, output_dir, param_label);
        fprintf('  TCD analysis completed in %.1f seconds\n', toc);
        
        %% Store results
        sweep_results{i} = struct(...
            'sweep_value', sweep_value, ...
            'param_label', param_label, ...
            'kPon_max', kPon_max_current, ...
            'SD', SD_current, ...
            'grid_file', grid_file, ...
            'interp_file', interp_file, ...
            'analysis_file', analysis_file, ...
            'analysis_results', analysis_results);
        
        sweep_success(i) = true;
        fprintf('SUCCESS: %s completed!\n\n', param_label);
        
    catch ME
        fprintf('ERROR in %s: %s\n', param_label, ME.message);
        fprintf('Skipping to next parameter value...\n\n');
        sweep_success(i) = false;
    end
end

fprintf('========================================\n');
fprintf('Parameter sweep completed!\n');
fprintf('Successful runs: %d/%d\n\n', sum(sweep_success), n_sweep);

%% --- COMPILE AND PLOT RESULTS ---

if sum(sweep_success) == 0
    error('No successful parameter sweep runs. Cannot generate plots.');
end

fprintf('Compiling results and generating plots...\n');

% Extract TCD data for plotting
figure('Position', [100, 100, 1000, 700]);

% Define color map for different parameter values
% Use jet colormap for compatibility with older MATLAB versions
colors = jet(n_sweep);

% Plot TCD vs gene length for each successful sweep
hold on;
legend_entries = {};

for i = 1:n_sweep
    if ~sweep_success(i)
        continue;
    end
    
    % Extract TCD results
    results = sweep_results{i}.analysis_results;
    L_TCD = results.TCD.gene_lengths;
    TCD_thresholds = results.TCD.thresholds;
    TCD_values = results.TCD.values;
    valid_idx = results.TCD.valid_indices;
    
    % Find the column corresponding to the desired threshold
    [~, threshold_col] = min(abs(TCD_thresholds - TCD_threshold));
    
    % Extract valid data
    L_valid = L_TCD(valid_idx);
    TCD_valid = TCD_values(valid_idx, threshold_col);
    
    % Plot with unique color and marker
    sweep_val = sweep_results{i}.sweep_value;
    if strcmp(sweep_parameter, 'kPon_max')
        label = sprintf('k_{Pon,max} = %.1f', sweep_val);
    else
        label = sprintf('SD = %d kb', sweep_val/1000);
    end
    
    semilogx(L_valid/1000, TCD_valid, 'o-', 'LineWidth', 2, 'MarkerSize', 6, ...
        'Color', colors(i, :), 'DisplayName', label);
    legend_entries{end+1} = label;
end

hold off;

% Format plot
xlabel('TSS-to-PAS Distance (kb)', 'FontSize', 14, 'FontWeight', 'bold');
ylabel(sprintf('TCD (bp, %.0f%% threshold)', TCD_threshold*100), 'FontSize', 14, 'FontWeight', 'bold');
title(sprintf('TCD vs Gene Length: %s Sweep (%s)', sweep_name, fixed_param_name), ...
    'FontSize', 16, 'FontWeight', 'bold');
legend('Location', 'best', 'FontSize', 11);
grid on;
set(gca, 'FontSize', 12);

% Save main plot
main_plot_file = fullfile(output_dir, sprintf('TCD_vs_L_%s_sweep.png', sweep_parameter));
saveas(gcf, main_plot_file);
fprintf('Main plot saved: %s\n', main_plot_file);

%% --- ADDITIONAL ANALYSIS PLOTS ---

% Plot 2: TCD at specific gene lengths vs sweep parameter
fprintf('Generating additional analysis plots...\n');

% Select representative gene lengths for analysis
representative_L = [5000, 20000, 50000, 100000];  % 5kb, 20kb, 50kb, 100kb
n_rep_L = length(representative_L);

figure('Position', [200, 200, 1000, 700]);

% Pre-allocate for TCD at representative lengths
TCD_at_L = zeros(n_sweep, n_rep_L);
successful_indices = find(sweep_success);

for i = 1:length(successful_indices)
    idx = successful_indices(i);
    results = sweep_results{idx}.analysis_results;
    
    L_TCD = results.TCD.gene_lengths;
    TCD_thresholds = results.TCD.thresholds;
    TCD_values = results.TCD.values;
    valid_idx = results.TCD.valid_indices;
    
    [~, threshold_col] = min(abs(TCD_thresholds - TCD_threshold));
    
    L_valid = L_TCD(valid_idx);
    TCD_valid = TCD_values(valid_idx, threshold_col);
    
    % Interpolate TCD at representative gene lengths
    for j = 1:n_rep_L
        if min(L_valid) <= representative_L(j) && representative_L(j) <= max(L_valid)
            TCD_at_L(idx, j) = interp1(L_valid, TCD_valid, representative_L(j), 'linear');
        else
            TCD_at_L(idx, j) = NaN;
        end
    end
end

% Plot TCD vs sweep parameter for each representative gene length
colors_L = lines(n_rep_L);
hold on;

for j = 1:n_rep_L
    valid_sweep = sweep_success & ~isnan(TCD_at_L(:, j));
    if sum(valid_sweep) > 0
        plot_x = sweep_values(valid_sweep);
        plot_y = TCD_at_L(valid_sweep, j);
        
        % Use log scale for x-axis if sweeping SD
        if strcmp(sweep_parameter, 'SD')
            loglog(plot_x/1000, plot_y, 'o-', 'LineWidth', 2, 'MarkerSize', 8, ...
                'Color', colors_L(j, :), 'DisplayName', sprintf('L = %d kb', representative_L(j)/1000));
        else
            semilogx(plot_x, plot_y, 'o-', 'LineWidth', 2, 'MarkerSize', 8, ...
                'Color', colors_L(j, :), 'DisplayName', sprintf('L = %d kb', representative_L(j)/1000));
        end
    end
end

hold off;

% Format plot
if strcmp(sweep_parameter, 'SD')
    xlabel('Saturation Distance (kb)', 'FontSize', 14, 'FontWeight', 'bold');
else
    xlabel('k_{Pon,max}', 'FontSize', 14, 'FontWeight', 'bold');
end
ylabel(sprintf('TCD (bp, %.0f%% threshold)', TCD_threshold*100), 'FontSize', 14, 'FontWeight', 'bold');
title(sprintf('TCD vs %s at Different Gene Lengths', sweep_name), 'FontSize', 16, 'FontWeight', 'bold');
legend('Location', 'best', 'FontSize', 11);
grid on;
set(gca, 'FontSize', 12);

% Save secondary plot
secondary_plot_file = fullfile(output_dir, sprintf('TCD_vs_%s_at_fixed_L.png', sweep_parameter));
saveas(gcf, secondary_plot_file);
fprintf('Secondary plot saved: %s\n', secondary_plot_file);

%% --- SAVE COMPILED RESULTS ---

fprintf('Saving compiled sweep results...\n');

% Create comprehensive results structure
compiled_results = struct();
compiled_results.metadata.creation_date = datestr(now);
compiled_results.metadata.sweep_parameter = sweep_parameter;
compiled_results.metadata.sweep_values = sweep_values;
compiled_results.metadata.n_successful = sum(sweep_success);
compiled_results.metadata.TCD_threshold_plotted = TCD_threshold;
compiled_results.metadata.output_directory = output_dir;

if strcmp(sweep_parameter, 'kPon_max')
    compiled_results.metadata.fixed_SD = SD_default;
else
    compiled_results.metadata.fixed_kPon_max = kPon_max_default;
end

compiled_results.sweep_results = sweep_results;
compiled_results.sweep_success = sweep_success;
compiled_results.representative_L = representative_L;
compiled_results.TCD_at_representative_L = TCD_at_L;

% Save compiled results
compiled_file = fullfile(output_dir, 'compiled_sweep_results.mat');
save(compiled_file, 'compiled_results', '-v7.3');
fprintf('Compiled results saved: %s\n', compiled_file);

%% --- SUMMARY ---

fprintf('\n=== PARAMETER SWEEP PIPELINE COMPLETE ===\n');
fprintf('Summary:\n');
fprintf('  Sweep parameter: %s\n', sweep_parameter);
fprintf('  Successful runs: %d/%d\n', sum(sweep_success), n_sweep);
fprintf('  Output directory: %s\n', output_dir);
fprintf('  Main plot: %s\n', main_plot_file);
fprintf('  Secondary plot: %s\n', secondary_plot_file);
fprintf('  Compiled results: %s\n', compiled_file);

fprintf('\nKey findings:\n');
for i = 1:n_sweep
    if sweep_success(i)
        results = sweep_results{i}.analysis_results;
        fprintf('  %s: TCD range = %.0f to %.0f bp\n', ...
            sweep_results{i}.param_label, ...
            min(results.TCD.values(results.TCD.valid_indices, 2)), ...
            max(results.TCD.values(results.TCD.valid_indices, 2)));
    end
end

fprintf('\nPipeline completed successfully!\n');

%% ========================================================================
%% HELPER FUNCTIONS
%% ========================================================================

function [grid_results, grid_file] = run_grid_generation(kPon_max, SD_bp, output_dir, param_label)
    % Run grid generation with specified parameters
    
    global N PAS N_PAS;
    
    % Parameter ranges (keep same as original)
    R_total_base = 70000;
    E_total_base = 70000;
    R_free_min = 10;
    R_free_max = 1000;
    R_free_points = 4;
    E_free_min = 0.1 * E_total_base;
    E_free_max = 0.9 * E_total_base;
    E_free_points = 4;
    L_min = 2500;
    L_max = 200000;
    L_points = 8;
    after_PAS_length = 5000;
    
    % Base parameters
    P_base.L_a = 100;
    k_in_global = 2;
    num_active_genes = 10000;
    P_base.k_in = k_in_global / num_active_genes;
    P_base.kEon = 0.00025;
    P_base.kEoff = 10;
    P_base.k_e = 65/100;
    P_base.k_e2 = 30/100;
    P_base.kHon = 0.2;
    P_base.kHoff = 0.0125;
    P_base.kc = 0.05;
    P_base.kPon_min = 0.01;
    P_base.kPon_max = kPon_max;  % SWEPT PARAMETER
    P_base.kPoff_min = 0.1;
    P_base.kPoff_const = 1;
    P_base.SD_bp = SD_bp;  % SWEPT PARAMETER
    P_base.kPon_option = 1;
    P_base.EBindingNumber = 1;
    
    % Create parameter grids
    R_free_values = linspace(R_free_min, R_free_max, R_free_points);
    E_free_values = linspace(E_free_min, E_free_max, E_free_points);
    L_values = logspace(log10(L_min), log10(L_max), L_points);
    
    [R_grid, E_grid, L_grid] = meshgrid(R_free_values, E_free_values, L_values);
    R_free_vec = R_grid(:);
    E_free_vec = E_grid(:);
    L_vec = L_grid(:);
    n_points = length(R_free_vec);
    
    R_occupied_vec = zeros(n_points, 1);
    E_occupied_vec = zeros(n_points, 1);
    success_flag = zeros(n_points, 1);
    
    % Parallel computation
    parfor i = 1:n_points
        R_free_i = R_free_vec(i);
        E_free_i = E_free_vec(i);
        L_i = L_vec(i);
        
        P_i = P_base;
        P_i.Pol_free = R_free_i;
        P_i.E_free = E_free_i;
        P_i.PASposition = L_i;
        P_i.geneLength_bp = L_i + after_PAS_length;
        
        [R_sol, REH_sol, avg_E_bound] = run_single_gene_simulation(P_i);
        R_occupied_i = sum(R_sol) + sum(REH_sol);
        E_occupied_i = calculate_E_occupied(R_sol, REH_sol, avg_E_bound);
        
        R_occupied_vec(i) = R_occupied_i;
        E_occupied_vec(i) = E_occupied_i;
        success_flag(i) = 1;
    end
    
    % Organize results
    grid_results = struct();
    grid_results.metadata.creation_date = datestr(now);
    grid_results.metadata.param_label = param_label;
    grid_results.metadata.success_rate = sum(success_flag)/n_points*100;
    
    grid_results.parameters.R_free_range = [R_free_min, R_free_max];
    grid_results.parameters.E_free_range = [E_free_min, E_free_max];
    grid_results.parameters.L_range = [L_min, L_max];
    grid_results.parameters.R_free_points = R_free_points;
    grid_results.parameters.E_free_points = E_free_points;
    grid_results.parameters.L_points = L_points;
    grid_results.parameters.after_PAS_length = after_PAS_length;
    grid_results.parameters.base_parameters = P_base;
    grid_results.parameters.R_total_base = R_total_base;
    grid_results.parameters.E_total_base = E_total_base;
    grid_results.parameters.num_active_genes = num_active_genes;
    
    grid_results.grid.R_free_values = R_free_values;
    grid_results.grid.E_free_values = E_free_values;
    grid_results.grid.L_values = L_values;
    
    grid_results.data.R_free_vec = R_free_vec;
    grid_results.data.E_free_vec = E_free_vec;
    grid_results.data.L_vec = L_vec;
    grid_results.data.R_occupied_vec = R_occupied_vec;
    grid_results.data.E_occupied_vec = E_occupied_vec;
    grid_results.data.success_flag = success_flag;
    
    % Save
    grid_file = fullfile(output_dir, sprintf('grid_data_%s.mat', param_label));
    save(grid_file, 'grid_results', '-v7.3');
end

function [interp_results, interp_file] = run_interpolation_build(grid_file, output_dir, param_label)
    % Build interpolation functions from grid data
    
    load(grid_file, 'grid_results');
    
    % Extract and clean data
    R_free_data = grid_results.data.R_free_vec;
    E_free_data = grid_results.data.E_free_vec;
    L_data = grid_results.data.L_vec;
    R_occupied_data = grid_results.data.R_occupied_vec;
    E_occupied_data = grid_results.data.E_occupied_vec;
    success_flags = grid_results.data.success_flag;
    
    valid_indices = (success_flags == 1) & ~isnan(R_occupied_data) & ~isnan(E_occupied_data);
    
    R_free_clean = R_free_data(valid_indices);
    E_free_clean = E_free_data(valid_indices);
    L_clean = L_data(valid_indices);
    R_occupied_clean = R_occupied_data(valid_indices);
    E_occupied_clean = E_occupied_data(valid_indices);
    
    % Build interpolants
    R_occupied_interp = scatteredInterpolant(R_free_clean, E_free_clean, L_clean, R_occupied_clean, 'linear', 'nearest');
    E_occupied_interp = scatteredInterpolant(R_free_clean, E_free_clean, L_clean, E_occupied_clean, 'linear', 'nearest');
    
    % Gene length distribution
    median_length = 22000;
    log_sigma = 0.68;
    mu_ln = log(median_length);
    sigma_ln = log_sigma * log(10);
    gene_length_pdf = @(L) lognpdf(L, mu_ln, sigma_ln);
    
    % Create results structure
    interp_results = struct();
    interp_results.metadata.creation_date = datestr(now);
    interp_results.metadata.param_label = param_label;
    interp_results.metadata.source_grid_file = grid_file;
    interp_results.metadata.n_valid_points = sum(valid_indices);
    
    interp_results.gene_length_distribution.type = 'log-normal';
    interp_results.gene_length_distribution.median_bp = median_length;
    interp_results.gene_length_distribution.log_sigma = log_sigma;
    interp_results.gene_length_distribution.mu_ln = mu_ln;
    interp_results.gene_length_distribution.sigma_ln = sigma_ln;
    
    interp_results.functions.R_occupied_interp = R_occupied_interp;
    interp_results.functions.E_occupied_interp = E_occupied_interp;
    interp_results.functions.gene_length_pdf = gene_length_pdf;
    interp_results.original_grid = grid_results.parameters;
    
    % Save
    interp_file = fullfile(output_dir, sprintf('interpolation_%s.mat', param_label));
    save(interp_file, 'interp_results', '-v7.3');
end

function [analysis_results, analysis_file] = run_TCD_analysis(interp_file, output_dir, param_label)
    % Run TCD analysis using interpolation functions
    
    load(interp_file, 'interp_results');
    
    R_occupied_interp = interp_results.functions.R_occupied_interp;
    E_occupied_interp = interp_results.functions.E_occupied_interp;
    gene_length_pdf = interp_results.functions.gene_length_pdf;
    
    R_total_target = interp_results.original_grid.R_total_base;
    E_total_target = interp_results.original_grid.E_total_base;
    num_active_genes = interp_results.original_grid.num_active_genes;
    base_params = interp_results.original_grid.base_parameters;
    
    % Integration domain
    L_min_int = 1000;
    L_max_int = 300000;
    L_integration_points = 200;
    L_integration = logspace(log10(L_min_int), log10(L_max_int), L_integration_points);
    dL = diff([L_integration(1)/1.1, L_integration]);
    
    f_L = gene_length_pdf(L_integration);
    f_L_normalized = f_L ./ sum(f_L .* dL);
    
    % Solve for self-consistent solution
    residual_function = @(x) conservation_equations_residual(x, R_occupied_interp, E_occupied_interp, ...
        L_integration, f_L_normalized, dL, R_total_target, E_total_target, num_active_genes);
    
    initial_guess = [0.3 * R_total_target, 0.3 * E_total_target];
    options = optimoptions('fsolve', 'Display', 'off', 'FunctionTolerance', 1e-6, 'MaxIterations', 100);
    
    [solution, ~, exitflag] = fsolve(residual_function, initial_guess, options);
    
    if exitflag <= 0
        error('Failed to find self-consistent solution');
    end
    
    R_free_solution = solution(1);
    E_free_solution = solution(2);
    
    % Calculate TCD for different gene lengths
    L_TCD_analysis = logspace(log10(2000), log10(200000), 20);
    n_L_TCD = length(L_TCD_analysis);
    TCD_thresholds = [0.25, 0.5, 0.75];
    n_thresholds = length(TCD_thresholds);
    
    TCD_results = zeros(n_L_TCD, n_thresholds);
    termination_profiles = cell(n_L_TCD, 1);
    
    for i = 1:n_L_TCD
        L_current = L_TCD_analysis(i);
        try
            [termination_profile, distances_bp] = calculate_termination_profile(L_current, ...
                R_free_solution, E_free_solution, base_params);
            termination_profiles{i} = struct('distances', distances_bp, 'profile', termination_profile);
            
            for j = 1:n_thresholds
                TCD_results(i, j) = calculate_TCD_from_profile(termination_profile, distances_bp, TCD_thresholds(j));
            end
        catch
            TCD_results(i, :) = NaN;
        end
    end
    
    % Create results structure
    valid_TCD = ~isnan(TCD_results(:, 1));
    
    analysis_results = struct();
    analysis_results.metadata.creation_date = datestr(now);
    analysis_results.metadata.param_label = param_label;
    analysis_results.metadata.source_interpolation_file = interp_file;
    
    analysis_results.solution.R_free = R_free_solution;
    analysis_results.solution.E_free = E_free_solution;
    analysis_results.solution.R_total_target = R_total_target;
    analysis_results.solution.E_total_target = E_total_target;
    
    analysis_results.TCD.gene_lengths = L_TCD_analysis;
    analysis_results.TCD.thresholds = TCD_thresholds;
    analysis_results.TCD.values = TCD_results;
    analysis_results.TCD.valid_indices = valid_TCD;
    analysis_results.TCD.termination_profiles = termination_profiles;
    
    % Save
    analysis_file = fullfile(output_dir, sprintf('TCD_analysis_%s.mat', param_label));
    save(analysis_file, 'analysis_results', '-v7.3');
end

%% --- SHARED HELPER FUNCTIONS (from original scripts) ---

function [R_sol, REH_sol, avg_E_bound] = run_single_gene_simulation(P)
    global N PAS N_PAS;
    syms Ef real;
    
    L_a = P.L_a;
    N = floor(P.geneLength_bp / L_a);
    PAS = floor(P.PASposition / L_a);
    N_PAS = N - PAS + 1;
    SD = floor(P.SD_bp / L_a);
    
    [r_E_BeforePas] = compute_steady_states(P, P.EBindingNumber + 1);
    
    kPon_vals = zeros(1, PAS);
    if PAS <= SD
        kPon_vals = linspace(P.kPon_min, P.kPon_min + (P.kPon_max - P.kPon_min) * (PAS / SD), PAS);
    else
        if P.kPon_option == 1
            kPon_vals(1:SD) = linspace(P.kPon_min, P.kPon_max, SD);
            kPon_vals((SD+1):PAS) = P.kPon_max;
        else
            kPon_vals(1:SD) = linspace(P.kPon_min, P.kPon_max, SD);
            slope = (P.kPon_max - P.kPon_min) / SD;
            for i = (SD+1):PAS
                kPon_vals(i) = P.kPon_max + slope * (i - SD);
            end
        end
    end
    
    RE_vals = sym(zeros(P.EBindingNumber + 1, N));
    for e = 1:(P.EBindingNumber + 1)
        for idx = 1:length(kPon_vals)
            RE_vals(e, idx) = subs(r_E_BeforePas(e), {'kPon', 'kPoff'}, {kPon_vals(idx), P.kPoff_const});
        end
        for idx = (PAS+1):N
            kPon_val = kPon_vals(end);
            RE_vals(e, idx) = subs(r_E_BeforePas(e), {'kPon', 'kPoff'}, {kPon_val, P.kPoff_const});
        end
    end
    
    P.RE_val_bind_E = matlabFunction(simplify(sum(sym(1:P.EBindingNumber)' .* RE_vals(2:end, :), 1)), 'Vars', {Ef});
    
    X_guess = 1e-6 * ones(N + N_PAS, 1);
    options = optimoptions('fsolve', 'Display', 'off', 'FunctionTolerance', 1e-8);
    
    avg_E_bound = P.RE_val_bind_E(P.E_free);
    P.kHon = P.kHon * avg_E_bound(PAS);
    X_final = fsolve(@(xx) ode_dynamics_multipleE(xx, P), X_guess, options);
    
    R_sol = X_final(1:N);
    REH_sol = X_final((N+1):(N+N_PAS));
end

function E_occupied = calculate_E_occupied(R_sol, REH_sol, avg_E_bound)
    global PAS;
    avg_E_bound_profile = avg_E_bound;
    E_bound_R = sum(R_sol .* avg_E_bound_profile');
    E_bound_REH = sum(REH_sol .* avg_E_bound_profile(PAS:end)');
    E_occupied = E_bound_R + E_bound_REH;
end

function dxdt = ode_dynamics_multipleE(X, P)
    global N PAS;
    k_in = P.k_in;
    k_e = P.k_e;
    k_e2 = P.k_e2;
    kHoff_t = P.kHoff;
    kc_t = P.kc;
    kHon_t = P.kHon;
    
    R = X(1:N);
    REH = X(N+1:end);
    dxdt = zeros(length(X), 1);
    
    n = 1;
    dxdt(n) = P.Pol_free*k_in - k_e*R(n);
    
    for n = 2:(PAS-1)
        dxdt(n) = k_e*R(n-1) - k_e*R(n);
    end
    
    n = PAS;
    j = n - PAS + 1;
    dxdt(n) = k_e*R(n-1) - k_e*R(n) - kHon_t*R(n) + kHoff_t*REH(j);
    dxdt(N+j) = -k_e2*REH(j) + kHon_t*R(n) - kHoff_t*REH(j);
    
    for n = (PAS+1):N
        j = n - PAS + 1;
        dxdt(n) = k_e*R(n-1) - k_e*R(n) - kHon_t*R(n) + kHoff_t*REH(j);
        dxdt(N+j) = k_e2*REH(j-1) - k_e2*REH(j) + kHon_t*R(n) - kHoff_t*REH(j) - kc_t*REH(j);
    end
end

function residual = conservation_equations_residual(x, R_interp, E_interp, L_vals, f_L_vals, dL_vals, R_target, E_target, num_genes)
    R_free = x(1);
    E_free = x(2);
    [R_total_calc, E_total_calc] = conservation_equations(R_free, E_free, R_interp, E_interp, L_vals, f_L_vals, dL_vals, num_genes);
    residual = [R_total_calc - R_target, E_total_calc - E_target];
end

function [R_total_calc, E_total_calc] = conservation_equations(R_free, E_free, R_interp, E_interp, L_vals, f_L_vals, dL_vals, num_genes)
    n_L = length(L_vals);
    R_occupied_vals = zeros(1, n_L);
    E_occupied_vals = zeros(1, n_L);
    
    for i = 1:n_L
        try
            R_occupied_vals(i) = R_interp(R_free, E_free, L_vals(i));
            E_occupied_vals(i) = E_interp(R_free, E_free, L_vals(i));
        catch
            R_occupied_vals(i) = 0;
            E_occupied_vals(i) = 0;
        end
    end
    
    R_integral = sum(R_occupied_vals .* f_L_vals .* dL_vals);
    E_integral = sum(E_occupied_vals .* f_L_vals .* dL_vals);
    
    R_total_calc = R_free + R_integral * num_genes;
    E_total_calc = E_free + E_integral * num_genes;
end

function [termination_profile, distances_bp] = calculate_termination_profile(tss_to_pas_distance, R_free, E_free, base_params)
    syms Ef real;
    
    P = base_params;
    P.E_free = E_free;
    P.Pol_free = R_free;
    
    after_PAS_length = 5000;
    P.PASposition = tss_to_pas_distance;
    P.geneLength_bp = tss_to_pas_distance + after_PAS_length;
    
    [R_sol, REH_sol, ~] = run_single_gene_simulation(P);
    
    flux_cleavage_per_node = P.kc * REH_sol;
    flux_R_exit = P.k_e * R_sol(end);
    flux_REH_exit = P.k_e2 * REH_sol(end);
    total_outflux = sum(flux_cleavage_per_node) + flux_R_exit + flux_REH_exit;
    
    if total_outflux > 1e-9
        cumulative_exit_flux = cumsum(flux_cleavage_per_node);
        termination_profile = cumulative_exit_flux / total_outflux;
    else
        termination_profile = zeros(size(REH_sol));
    end
    
    nodes_post_pas = 1:length(REH_sol);
    distances_bp = nodes_post_pas * P.L_a;
end

function TCD = calculate_TCD_from_profile(termination_profile, distances_bp, threshold)
    if isempty(termination_profile) || all(isnan(termination_profile))
        TCD = NaN;
        return;
    end
    
    threshold_idx = find(termination_profile >= threshold, 1, 'first');
    
    if isempty(threshold_idx)
        TCD = distances_bp(end);
    elseif threshold_idx == 1
        TCD = distances_bp(1);
    else
        x1 = termination_profile(threshold_idx-1);
        x2 = termination_profile(threshold_idx);
        d1 = distances_bp(threshold_idx-1);
        d2 = distances_bp(threshold_idx);
        TCD = d1 + (threshold - x1) * (d2 - d1) / (x2 - x1);
    end
end

