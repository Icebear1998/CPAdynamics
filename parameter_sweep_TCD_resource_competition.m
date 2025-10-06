% PARAMETER_SWEEP_TCD_RESOURCE_COMPETITION.m
% Parameter sweep for TCD with full resource competition framework
%
% This script performs parameter sweeps using the complete three-step pipeline:
% 1. Generate grid data for each parameter value (reuses generate_gene_length_grid.m logic)
% 2. Build interpolation functions (reuses build_gene_length_interpolation.m logic)
% 3. Solve conservation equations and calculate TCD (reuses analyze_gene_length_TCD.m logic)
%
% Key difference from parameter_sweep_TCD.m:
% - Includes genome-wide resource competition effects
% - Self-consistent (R_free, E_free) solution for each parameter value
% - Captures how parameter changes affect resource availability
% - More computationally intensive but scientifically rigorous
%
% IMPORTANT: This script reuses helper functions from the existing three-file pipeline
% to avoid code duplication and ensure consistency.

fprintf('=== Parameter Sweep for TCD with Resource Competition ===\n');
fprintf('Full three-step pipeline for each parameter value\n\n');

%% --- CONFIGURATION ---
fprintf('Configuration:\n');

% Choose parameter to sweep
% Options: 'SD_bp', 'kPon_max', 'kHon', 'kHoff', 'kc', 'kEon', 'EBindingNumber'
parameter_to_sweep = 'SD_bp';  % Change this to sweep different parameters

% Sweep configuration
switch parameter_to_sweep
    case 'SD_bp'
        param_values = [5000, 10000, 20000, 30000, 50000];  % Different saturation distances
        param_label = 'Saturation Distance (kb)';
        param_display = @(x) x/1000;  % Convert to kb for display
        
    case 'kPon_max'
        param_values = [0.05, 0.1, 0.2, 0.5, 1.0];
        param_label = 'k_{Pon,max}';
        param_display = @(x) x;
        
    case 'kHon'
        param_values = [0.05, 0.1, 0.2, 0.4, 0.8];
        param_label = 'k_{Hon}';
        param_display = @(x) x;
        
    case 'kHoff'
        param_values = [0.00625, 0.0125, 0.025, 0.05];
        param_label = 'k_{Hoff}';
        param_display = @(x) x;
        
    case 'kc'
        param_values = [0.025, 0.05, 0.1, 0.2];
        param_label = 'k_c (cleavage rate)';
        param_display = @(x) x;
        
    case 'kEon'
        param_values = [0.0001, 0.00025, 0.0005, 0.001];
        param_label = 'k_{Eon}';
        param_display = @(x) x;
        
    case 'EBindingNumber'
        param_values = [1, 2, 3, 4, 5, 6];
        param_label = 'E Binding Number';
        param_display = @(x) x;
        
    otherwise
        error('Unknown parameter: %s', parameter_to_sweep);
end

n_param_values = length(param_values);
fprintf('  Parameter to sweep: %s\n', parameter_to_sweep);
fprintf('  Number of values: %d\n', n_param_values);
fprintf('  Range: ');
for i = 1:n_param_values
    fprintf('%.3g', param_display(param_values(i)));
    if i < n_param_values
        fprintf(', ');
    end
end
fprintf('\n');

% Gene lengths for TCD analysis
gene_lengths_kb = [5, 10, 20, 30, 50, 100, 150];  % Representative gene lengths
n_gene_lengths = length(gene_lengths_kb);
gene_lengths_bp = gene_lengths_kb * 1000;

fprintf('  Gene lengths to analyze: ');
fprintf('%.0f ', gene_lengths_kb);
fprintf('kb\n');

% TCD threshold
TCD_threshold = 0.5;  % 50% termination threshold
fprintf('  TCD threshold: %.0f%%\n', TCD_threshold*100);

% Grid resolution (trade-off between accuracy and speed)
% Reduced from full analysis for computational efficiency
R_free_points = 4;  % Resolution for R_free in grid
E_free_points = 4;  % Resolution for E_free in grid
L_points = 8;       % Resolution for gene length in grid
fprintf('  Grid resolution: R_free=%d, E_free=%d, L=%d points\n', R_free_points, E_free_points, L_points);

% Start parallel pool
if isempty(gcp('nocreate'))
    parpool;
end

%% --- BASE PARAMETERS ---
fprintf('\nSetting up base parameters...\n');

P_base = struct();
P_base.L_a = 100;
P_base.k_in = 2;
P_base.kEon = 0.00025;
P_base.kEoff = 10;
P_base.k_e = 65/100;
P_base.k_e2 = 30/100;
P_base.kHon = 0.2;
P_base.kHoff = 0.0125;
P_base.kc = 0.05;
P_base.kPon_min = 0.01;
P_base.kPon_max = 0.1;
P_base.kPoff_const = 1;
P_base.SD_bp = 20000;
P_base.kPon_option = 1;  % Saturate at SD
P_base.EBindingNumber = 1;
P_base.after_PAS_length = 5000;

% Resource pool targets
R_total_target = 70000;
E_total_target = 70000;

fprintf('Base parameters set.\n');

%% --- OUTPUT DIRECTORY ---
output_dir = sprintf('SecondVersionResults/TCD_ParameterSweep_ResourceCompetition/%s/', parameter_to_sweep);
if ~exist(output_dir, 'dir')
    mkdir(output_dir);
end
fprintf('Output directory: %s\n', output_dir);

%% --- PARAMETER SWEEP LOOP ---
fprintf('\n=== Starting Parameter Sweep ===\n');
fprintf('This will run the full three-step pipeline for each parameter value.\n');
fprintf('Estimated time: ~5-15 minutes per parameter value\n\n');

% Pre-allocate results storage
sweep_results = struct();
sweep_results.metadata.creation_date = datestr(now);
sweep_results.metadata.parameter_name = parameter_to_sweep;
sweep_results.metadata.parameter_values = param_values;
sweep_results.metadata.base_parameters = P_base;
sweep_results.metadata.R_total_target = R_total_target;
sweep_results.metadata.E_total_target = E_total_target;

% Storage for results from each parameter
sweep_results.data = cell(n_param_values, 1);

% Storage for summary data
R_free_solutions = zeros(n_param_values, 1);
E_free_solutions = zeros(n_param_values, 1);
TCD_matrix = zeros(n_param_values, n_gene_lengths);
success_flags = ones(n_param_values, 1);

total_sweep_time = tic;

for p_idx = 1:n_param_values
    param_value = param_values(p_idx);
    
    fprintf('\n+===========================================================+\n');
    fprintf('| Parameter %d/%d: %s = %.3g %-20s|\n', p_idx, n_param_values, ...
        parameter_to_sweep, param_display(param_value), '');
    fprintf('+===========================================================+\n\n');
    
    param_time = tic;
    
    try
        % Update parameter
        P_sweep = P_base;
        P_sweep.(parameter_to_sweep) = param_value;
        
        %% STEP 1: GENERATE GRID DATA
        fprintf('STEP 1/3: Generating grid data...\n');
        step1_time = tic;
        
        grid_data = generate_grid_data(P_sweep, R_total_target, E_total_target, ...
            R_free_points, E_free_points, L_points);
        
        fprintf('  Grid generation completed in %.1f seconds\n', toc(step1_time));
        fprintf('  Valid points: %d/%d (%.1f%%)\n', grid_data.n_valid, grid_data.n_total, ...
            grid_data.n_valid/grid_data.n_total*100);
        
        %% STEP 2: BUILD INTERPOLATION
        fprintf('STEP 2/3: Building interpolation functions...\n');
        step2_time = tic;
        
        interp_data = build_interpolation(grid_data, P_sweep);
        
        fprintf('  Interpolation built in %.1f seconds\n', toc(step2_time));
        fprintf('  Mean interpolation error: R=%.2f%%, E=%.2f%%\n', ...
            interp_data.validation.R_mean_error, interp_data.validation.E_mean_error);
        
        %% STEP 3: SOLVE CONSERVATION AND CALCULATE TCD
        fprintf('STEP 3/3: Solving conservation equations and calculating TCD...\n');
        step3_time = tic;
        
        tcd_data = analyze_tcd_with_conservation(interp_data, gene_lengths_bp, ...
            TCD_threshold, R_total_target, E_total_target);
        
        fprintf('  TCD analysis completed in %.1f seconds\n', toc(step3_time));
        fprintf('  Solution: R_free = %.0f (%.1f%%), E_free = %.0f (%.1f%%)\n', ...
            tcd_data.R_free_solution, tcd_data.R_free_solution/R_total_target*100, ...
            tcd_data.E_free_solution, tcd_data.E_free_solution/E_total_target*100);
        
        % Store results
        sweep_results.data{p_idx} = struct();
        sweep_results.data{p_idx}.parameter_value = param_value;
        sweep_results.data{p_idx}.grid_data = grid_data;
        sweep_results.data{p_idx}.interp_data = interp_data;
        sweep_results.data{p_idx}.tcd_data = tcd_data;
        
        % Store summary data
        R_free_solutions(p_idx) = tcd_data.R_free_solution;
        E_free_solutions(p_idx) = tcd_data.E_free_solution;
        TCD_matrix(p_idx, :) = tcd_data.TCD_values;
        success_flags(p_idx) = 1;
        
        fprintf('\n✓ Parameter value %.3g completed in %.1f minutes\n', ...
            param_display(param_value), toc(param_time)/60);
        
    catch ME
        fprintf('\n✗ ERROR for parameter value %.3g: %s\n', param_display(param_value), ME.message);
        fprintf('  Stack trace:\n');
        for k = 1:length(ME.stack)
            fprintf('    %s (line %d)\n', ME.stack(k).name, ME.stack(k).line);
        end
        
        success_flags(p_idx) = 0;
        R_free_solutions(p_idx) = NaN;
        E_free_solutions(p_idx) = NaN;
        TCD_matrix(p_idx, :) = NaN;
    end
end

total_time = toc(total_sweep_time);

fprintf('\n+===========================================================+\n');
fprintf('| PARAMETER SWEEP COMPLETED                                 |\n');
fprintf('+===========================================================+\n');
fprintf('Total time: %.1f minutes\n', total_time/60);
fprintf('Successful parameter values: %d/%d\n', sum(success_flags), n_param_values);

%% --- SUMMARY ANALYSIS ---
fprintf('\n=== Analysis Summary ===\n');

% Add summary data to results
sweep_results.summary.R_free_solutions = R_free_solutions;
sweep_results.summary.E_free_solutions = E_free_solutions;
sweep_results.summary.TCD_matrix = TCD_matrix;
sweep_results.summary.success_flags = success_flags;
sweep_results.summary.gene_lengths_kb = gene_lengths_kb;
sweep_results.summary.TCD_threshold = TCD_threshold;
sweep_results.summary.total_time_minutes = total_time/60;

% Analyze how resources change with parameter
valid_idx = success_flags == 1;
if sum(valid_idx) > 1
    fprintf('\nResource pool trends:\n');
    corr_R = corr(param_values(valid_idx)', R_free_solutions(valid_idx));
    corr_E = corr(param_values(valid_idx)', E_free_solutions(valid_idx));
    fprintf('  R_free vs %s: r = %.3f\n', parameter_to_sweep, corr_R);
    fprintf('  E_free vs %s: r = %.3f\n', parameter_to_sweep, corr_E);
    
    fprintf('\nTCD trends (correlation with %s):\n', parameter_to_sweep);
    for g_idx = 1:n_gene_lengths
        valid_tcd = valid_idx & ~isnan(TCD_matrix(:, g_idx));
        if sum(valid_tcd) > 1
            corr_tcd = corr(param_values(valid_tcd)', TCD_matrix(valid_tcd, g_idx));
            fprintf('  L=%.0f kb: r = %.3f', gene_lengths_kb(g_idx), corr_tcd);
            if abs(corr_tcd) > 0.7
                fprintf(' (strong)');%, corr_tcd > 0 ? 'positive' : 'negative');
            end
            fprintf('\n');
        end
    end
end

%% --- VISUALIZATION ---
fprintf('\nGenerating comprehensive visualizations...\n');

% Create comprehensive figure
fig = figure('Position', [50, 50, 1600, 1000]);

% Plot 1: Resource pools vs parameter
subplot(2, 3, 1);
yyaxis left;
plot(param_display(param_values(valid_idx)), R_free_solutions(valid_idx), 'o-', ...
    'LineWidth', 2, 'MarkerSize', 10);
ylabel('R_{free}', 'FontSize', 12);
yyaxis right;
plot(param_display(param_values(valid_idx)), E_free_solutions(valid_idx), 's-', ...
    'LineWidth', 2, 'MarkerSize', 10);
ylabel('E_{free}', 'FontSize', 12);
xlabel(param_label, 'FontSize', 12);
title('Self-Consistent Resource Pools', 'FontSize', 13, 'FontWeight', 'bold');
grid on;

% Plot 2: Resource fractions
subplot(2, 3, 2);
R_frac = R_free_solutions(valid_idx) / R_total_target * 100;
E_frac = E_free_solutions(valid_idx) / E_total_target * 100;
hold on;
plot(param_display(param_values(valid_idx)), R_frac, 'o-', 'LineWidth', 2, ...
    'MarkerSize', 10, 'DisplayName', 'R_{free} %');
plot(param_display(param_values(valid_idx)), E_frac, 's-', 'LineWidth', 2, ...
    'MarkerSize', 10, 'DisplayName', 'E_{free} %');
ylabel('Free Resource Fraction (%)', 'FontSize', 12);
xlabel(param_label, 'FontSize', 12);
title('Resource Availability', 'FontSize', 13, 'FontWeight', 'bold');
legend('Location', 'best');
grid on;
hold off;

% Plot 3: TCD vs parameter for different gene lengths
subplot(2, 3, 3);
colors = jet(n_gene_lengths);
hold on;
for g_idx = 1:n_gene_lengths
    valid_tcd = valid_idx & ~isnan(TCD_matrix(:, g_idx));
    if sum(valid_tcd) > 0
        plot(param_display(param_values(valid_tcd)), TCD_matrix(valid_tcd, g_idx), ...
            'o-', 'LineWidth', 2, 'MarkerSize', 8, 'Color', colors(g_idx, :), ...
            'DisplayName', sprintf('L = %.0f kb', gene_lengths_kb(g_idx)));
    end
end
xlabel(param_label, 'FontSize', 12);
ylabel(sprintf('TCD at %.0f%% (bp)', TCD_threshold*100), 'FontSize', 12);
title('TCD vs Parameter', 'FontSize', 13, 'FontWeight', 'bold');
legend('Location', 'best', 'FontSize', 9);
grid on;
hold off;

% Plot 4: TCD vs gene length for different parameters
subplot(2, 3, 4);
colors_param = parula(n_param_values);
hold on;
for p_idx = 1:n_param_values
    if success_flags(p_idx)
        valid_tcd = ~isnan(TCD_matrix(p_idx, :));
        if sum(valid_tcd) > 0
            semilogx(gene_lengths_kb(valid_tcd), TCD_matrix(p_idx, valid_tcd), ...
                'o-', 'LineWidth', 2, 'MarkerSize', 8, 'Color', colors_param(p_idx, :), ...
                'DisplayName', sprintf('%.3g', param_display(param_values(p_idx))));
        end
    end
end
xlabel('TSS-to-PAS Distance (kb)', 'FontSize', 12);
ylabel(sprintf('TCD at %.0f%% (bp)', TCD_threshold*100), 'FontSize', 12);
title('TCD vs Gene Length', 'FontSize', 13, 'FontWeight', 'bold');
lg = legend('Location', 'best', 'FontSize', 9);
title(lg, param_label);
grid on;
hold off;

% Plot 5: Heatmap of TCD values
subplot(2, 3, 5);
TCD_plot = TCD_matrix;
TCD_plot(~valid_idx, :) = NaN;
imagesc(1:n_gene_lengths, 1:n_param_values, TCD_plot);
colorbar;
xlabel('Gene Length Index', 'FontSize', 12);
ylabel('Parameter Value Index', 'FontSize', 12);
title('TCD Heatmap (bp)', 'FontSize', 13, 'FontWeight', 'bold');
set(gca, 'XTick', 1:n_gene_lengths, 'XTickLabel', ...
    arrayfun(@(x) sprintf('%.0f', x), gene_lengths_kb, 'UniformOutput', false));
set(gca, 'YTick', 1:n_param_values, 'YTickLabel', ...
    arrayfun(@(x) sprintf('%.3g', x), param_display(param_values), 'UniformOutput', false));

% Plot 6: TCD sensitivity (normalized)
subplot(2, 3, 6);
baseline_idx = ceil(n_param_values / 2);
if success_flags(baseline_idx)
    TCD_normalized = TCD_matrix ./ TCD_matrix(baseline_idx, :);
    TCD_normalized(isinf(TCD_normalized) | isnan(TCD_normalized)) = 1;
    
    hold on;
    for g_idx = 1:n_gene_lengths
        valid_tcd = valid_idx & ~isnan(TCD_normalized(:, g_idx));
        if sum(valid_tcd) > 0
            plot(param_display(param_values(valid_tcd)), TCD_normalized(valid_tcd, g_idx), ...
                'o-', 'LineWidth', 2, 'MarkerSize', 8, 'Color', colors(g_idx, :), ...
                'DisplayName', sprintf('L = %.0f kb', gene_lengths_kb(g_idx)));
        end
    end
    yline(1, 'k--', 'LineWidth', 1.5);
    xlabel(param_label, 'FontSize', 12);
    ylabel('Normalized TCD', 'FontSize', 12);
    title('TCD Sensitivity', 'FontSize', 13, 'FontWeight', 'bold');
    legend('Location', 'best', 'FontSize', 8);
    grid on;
    hold off;
end

sgtitle(sprintf('Parameter Sweep: %s with Resource Competition', parameter_to_sweep), ...
    'FontSize', 15, 'FontWeight', 'bold');

%% --- SAVE RESULTS ---
fprintf('\nSaving results...\n');

timestamp = datestr(now, 'yyyymmdd_HHMMSS');

% Save MATLAB file with all data
mat_filename = fullfile(output_dir, sprintf('sweep_results_%s.mat', timestamp));
save(mat_filename, 'sweep_results', '-v7.3');
fprintf('  Full results saved: %s\n', mat_filename);

% Save figure
fig_filename = fullfile(output_dir, sprintf('sweep_plots_%s.png', timestamp));
saveas(fig, fig_filename);
fprintf('  Figure saved: %s\n', fig_filename);

% Save summary text file
txt_filename = fullfile(output_dir, sprintf('sweep_summary_%s.txt', timestamp));
fid = fopen(txt_filename, 'w');

fprintf(fid, '%% Parameter Sweep with Resource Competition\n');
fprintf(fid, '%% Parameter: %s\n', parameter_to_sweep);
fprintf(fid, '%% Generated: %s\n', datestr(now));
fprintf(fid, '%% Total time: %.1f minutes\n', total_time/60);
fprintf(fid, '%% Success rate: %d/%d\n', sum(success_flags), n_param_values);
fprintf(fid, '%%\n');
fprintf(fid, '%% RESOURCE POOL SOLUTIONS\n');
fprintf(fid, '%% ParamValue\tR_free\tE_free\tR_free_pct\tE_free_pct\n');
for p_idx = 1:n_param_values
    if success_flags(p_idx)
        fprintf(fid, '%.6g\t%.2f\t%.2f\t%.2f\t%.2f\n', ...
            param_display(param_values(p_idx)), ...
            R_free_solutions(p_idx), E_free_solutions(p_idx), ...
            R_free_solutions(p_idx)/R_total_target*100, ...
            E_free_solutions(p_idx)/E_total_target*100);
    end
end

fprintf(fid, '%%\n');
fprintf(fid, '%% TCD VALUES (bp) at %.0f%% threshold\n', TCD_threshold*100);
fprintf(fid, '%% ParamValue');
for g_idx = 1:n_gene_lengths
    fprintf(fid, '\tL%.0fkb', gene_lengths_kb(g_idx));
end
fprintf(fid, '\n');

for p_idx = 1:n_param_values
    fprintf(fid, '%.6g', param_display(param_values(p_idx)));
    for g_idx = 1:n_gene_lengths
        if success_flags(p_idx) && ~isnan(TCD_matrix(p_idx, g_idx))
            fprintf(fid, '\t%.2f', TCD_matrix(p_idx, g_idx));
        else
            fprintf(fid, '\tNaN');
        end
    end
    fprintf(fid, '\n');
end

fclose(fid);
fprintf('  Summary saved: %s\n', txt_filename);

%% --- FINAL SUMMARY ---
fprintf('\n+===========================================================+\n');
fprintf('| ANALYSIS COMPLETE                                         |\n');
fprintf('+===========================================================+\n\n');

fprintf('Parameter: %s\n', parameter_to_sweep);
fprintf('Range: %.3g to %.3g\n', param_display(param_values(1)), param_display(param_values(end)));
fprintf('Success rate: %d/%d (%.1f%%)\n', sum(success_flags), n_param_values, ...
    sum(success_flags)/n_param_values*100);
fprintf('Total computation time: %.1f minutes\n', total_time/60);

if sum(valid_idx) > 1
    fprintf('\nKey insights:\n');
    fprintf('  • R_free range: %.0f to %.0f (%.1f%% to %.1f%% of total)\n', ...
        min(R_free_solutions(valid_idx)), max(R_free_solutions(valid_idx)), ...
        min(R_free_solutions(valid_idx))/R_total_target*100, ...
        max(R_free_solutions(valid_idx))/R_total_target*100);
    fprintf('  • E_free range: %.0f to %.0f (%.1f%% to %.1f%% of total)\n', ...
        min(E_free_solutions(valid_idx)), max(E_free_solutions(valid_idx)), ...
        min(E_free_solutions(valid_idx))/E_total_target*100, ...
        max(E_free_solutions(valid_idx))/E_total_target*100);
    
    valid_tcd = TCD_matrix(valid_idx, :);
    valid_tcd = valid_tcd(~isnan(valid_tcd));
    if ~isempty(valid_tcd)
        fprintf('  • TCD range: %.0f to %.0f bp\n', min(valid_tcd), max(valid_tcd));
    end
end

fprintf('\nOutput files:\n');
fprintf('  %s\n', mat_filename);
fprintf('  %s\n', fig_filename);
fprintf('  %s\n', txt_filename);

fprintf('\nParameter sweep with resource competition completed successfully!\n');

%% ========================================================================
%% HELPER FUNCTIONS
%% ========================================================================
% These functions reuse code from the existing three-file pipeline to ensure
% consistency and avoid duplication:
%
% 1. generate_grid_data() - wraps logic from generate_gene_length_grid.m
% 2. build_interpolation() - wraps logic from build_gene_length_interpolation.m
% 3. analyze_tcd_with_conservation() - wraps logic from analyze_gene_length_TCD.m
% 4. run_single_gene_simulation_local() - exact copy from generate_gene_length_grid.m
% 5. conservation functions - exact copies from analyze_gene_length_TCD.m
%
% All helper functions have "_local" suffix to avoid naming conflicts when
% this script is used alongside the original three-file pipeline.

function grid_data = generate_grid_data(P_base, R_total, E_total, R_pts, E_pts, L_pts)
    % Generate grid data (Step 1 of pipeline)
    % Reuses logic from generate_gene_length_grid.m
    
    % Define ranges
    R_free_min = 0.1 * R_total;
    R_free_max = 0.9 * R_total;
    E_free_min = 0.1 * E_total;
    E_free_max = 0.9 * E_total;
    L_min = 2500;
    L_max = 200000;
    after_PAS_length = P_base.after_PAS_length;
    
    % Create grids
    R_free_values = linspace(R_free_min, R_free_max, R_pts);
    E_free_values = linspace(E_free_min, E_free_max, E_pts);
    L_values = logspace(log10(L_min), log10(L_max), L_pts);
    
    [R_grid, E_grid, L_grid] = meshgrid(R_free_values, E_free_values, L_values);
    
    R_free_vec = R_grid(:);
    E_free_vec = E_grid(:);
    L_vec = L_grid(:);
    n_points = length(R_free_vec);
    
    R_occupied_vec = zeros(n_points, 1);
    E_occupied_vec = zeros(n_points, 1);
    success_flag = zeros(n_points, 1);
    
    % Parallel computation - reuses run_single_gene_simulation from generate_gene_length_grid.m
    parfor i = 1:n_points
        try
            P_i = P_base;
            P_i.Pol_total = R_free_vec(i);
            P_i.E_total = E_free_vec(i);
            P_i.PASposition = L_vec(i);
            P_i.geneLength_bp = L_vec(i) + after_PAS_length;
            
            % Call the existing function from generate_gene_length_grid.m
            [R_sol, REH_sol, avg_E_bound] = run_single_gene_simulation_local(P_i);
            
            % Calculate occupied amounts (same as calculate_E_occupied in generate_gene_length_grid.m)
            R_occupied_vec(i) = sum(R_sol) + sum(REH_sol);
            
            PAS = floor(P_i.PASposition / P_i.L_a);
            E_occupied_vec(i) = sum(R_sol .* avg_E_bound') + sum(REH_sol .* avg_E_bound(PAS:end)');
            success_flag(i) = 1;
        catch
            R_occupied_vec(i) = NaN;
            E_occupied_vec(i) = NaN;
            success_flag(i) = 0;
        end
    end
    
    % Package results
    grid_data.R_free_vec = R_free_vec;
    grid_data.E_free_vec = E_free_vec;
    grid_data.L_vec = L_vec;
    grid_data.R_occupied_vec = R_occupied_vec;
    grid_data.E_occupied_vec = E_occupied_vec;
    grid_data.success_flag = success_flag;
    grid_data.n_valid = sum(success_flag);
    grid_data.n_total = n_points;
    grid_data.R_free_values = R_free_values;
    grid_data.E_free_values = E_free_values;
    grid_data.L_values = L_values;
end

function interp_data = build_interpolation(grid_data, P_base)
    % Build interpolation functions (Step 2 of pipeline)
    % Reuses logic from build_gene_length_interpolation.m
    
    % Clean data (same as build_gene_length_interpolation.m lines 50-65)
    valid_idx = grid_data.success_flag == 1 & ~isnan(grid_data.R_occupied_vec) & ~isnan(grid_data.E_occupied_vec);
    
    R_free_clean = grid_data.R_free_vec(valid_idx);
    E_free_clean = grid_data.E_free_vec(valid_idx);
    L_clean = grid_data.L_vec(valid_idx);
    R_occupied_clean = grid_data.R_occupied_vec(valid_idx);
    E_occupied_clean = grid_data.E_occupied_vec(valid_idx);
    
    % Build interpolants (same as build_gene_length_interpolation.m lines 80-86)
    R_occupied_interp = scatteredInterpolant(R_free_clean, E_free_clean, L_clean, ...
        R_occupied_clean, 'linear', 'nearest');
    E_occupied_interp = scatteredInterpolant(R_free_clean, E_free_clean, L_clean, ...
        E_occupied_clean, 'linear', 'nearest');
    
    % Validate (same as build_gene_length_interpolation.m lines 93-113)
    n_test = min(100, floor(length(R_free_clean) * 0.1));
    test_indices = randperm(length(R_free_clean), n_test);
    
    R_error = abs(R_occupied_interp(R_free_clean(test_indices), E_free_clean(test_indices), L_clean(test_indices)) - R_occupied_clean(test_indices)) ./ R_occupied_clean(test_indices) * 100;
    E_error = abs(E_occupied_interp(R_free_clean(test_indices), E_free_clean(test_indices), L_clean(test_indices)) - E_occupied_clean(test_indices)) ./ E_occupied_clean(test_indices) * 100;
    
    % Gene length distribution (same as build_gene_length_interpolation.m lines 128-144)
    median_length = 22000;
    log_sigma = 0.68;
    mu_ln = log(median_length);
    sigma_ln = log_sigma * log(10);
    gene_length_pdf = @(L) lognpdf(L, mu_ln, sigma_ln);
    
    % Package results
    interp_data.R_occupied_interp = R_occupied_interp;
    interp_data.E_occupied_interp = E_occupied_interp;
    interp_data.gene_length_pdf = gene_length_pdf;
    interp_data.validation.R_mean_error = mean(R_error);
    interp_data.validation.E_mean_error = mean(E_error);
    interp_data.parameters = P_base;
end

function tcd_data = analyze_tcd_with_conservation(interp_data, gene_lengths_bp, threshold, R_target, E_target)
    % Solve conservation and calculate TCD (Step 3 of pipeline)
    % Reuses logic from analyze_gene_length_TCD.m
    
    % Set up integration domain (same as analyze_gene_length_TCD.m lines 52-66)
    L_min_int = 1000;
    L_max_int = 300000;
    L_integration_points = 100;  % Reduced for speed
    
    L_integration = logspace(log10(L_min_int), log10(L_max_int), L_integration_points);
    dL = diff([L_integration(1)/1.1, L_integration]);
    
    f_L = interp_data.gene_length_pdf(L_integration);
    f_L_normalized = f_L ./ sum(f_L .* dL);
    
    % Solve for self-consistent (R_free, E_free) (same as analyze_gene_length_TCD.m lines 111-165)
    % Uses conservation_equations_residual from analyze_gene_length_TCD.m
    residual_function = @(x) conservation_residual_local(x, interp_data.R_occupied_interp, ...
        interp_data.E_occupied_interp, L_integration, f_L_normalized, dL, R_target, E_target);
    
    R_free_guess = 0.3 * R_target;
    E_free_guess = 0.3 * E_target;
    initial_guess = [R_free_guess, E_free_guess];
    
    options = optimoptions('fsolve', 'Display', 'off', 'FunctionTolerance', 1e-6, 'MaxIterations', 100);
    
    [solution, ~, exitflag] = fsolve(residual_function, initial_guess, options);
    
    if exitflag <= 0
        error('Failed to find self-consistent solution');
    end
    
    R_free_solution = solution(1);
    E_free_solution = solution(2);
    
    % Calculate TCD for each gene length (same as analyze_gene_length_TCD.m lines 185-207)
    n_gene_lengths = length(gene_lengths_bp);
    TCD_values = zeros(1, n_gene_lengths);
    
    for i = 1:n_gene_lengths
        try
            [term_profile, distances_bp] = calc_termination_profile_local(gene_lengths_bp(i), ...
                R_free_solution, E_free_solution, interp_data.parameters);
            TCD_values(i) = calc_TCD_local(term_profile, distances_bp, threshold);
        catch
            TCD_values(i) = NaN;
        end
    end
    
    % Package results
    tcd_data.R_free_solution = R_free_solution;
    tcd_data.E_free_solution = E_free_solution;
    tcd_data.TCD_values = TCD_values;
    tcd_data.gene_lengths_bp = gene_lengths_bp;
    tcd_data.threshold = threshold;
end

% Conservation and termination calculation functions
% These are copies from analyze_gene_length_TCD.m with "_local" suffix to avoid conflicts

function residual = conservation_residual_local(x, R_interp, E_interp, L_vals, f_L_vals, dL_vals, R_target, E_target)
    % Conservation equation residual (from analyze_gene_length_TCD.m line 319)
    R_free = x(1);
    E_free = x(2);
    
    [R_total_calc, E_total_calc] = conservation_equations_local(R_free, E_free, R_interp, E_interp, L_vals, f_L_vals, dL_vals);
    residual = [R_total_calc - R_target, E_total_calc - E_target];
end

function [R_total_calc, E_total_calc] = conservation_equations_local(R_free, E_free, R_interp, E_interp, L_vals, f_L_vals, dL_vals)
    % Core conservation equation implementation (from analyze_gene_length_TCD.m line 351)
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
    
    R_total_calc = R_free + R_integral;
    E_total_calc = E_free + E_integral;
end

function [term_profile, distances_bp] = calc_termination_profile_local(tss_to_pas_distance, R_free, E_free, P_base)
    % Calculate termination profile (simplified from analyze_gene_length_TCD.m line 412)
    P = P_base;
    P.Pol_total = R_free;
    P.E_total = E_free;
    P.PASposition = tss_to_pas_distance;
    P.geneLength_bp = tss_to_pas_distance + P.after_PAS_length;
    
    [R_sol, REH_sol, P_sim] = run_single_gene_simulation_local(P);
    
    % Use calculate_pas_usage_profile for termination
    [exit_cdf, distances_bp] = calculate_pas_usage_profile(R_sol, REH_sol, P_sim);
    term_profile = exit_cdf;
end

function TCD = calc_TCD_local(term_profile, distances_bp, threshold)
    % Calculate TCD from termination profile (from analyze_gene_length_TCD.m line 465)
    if isempty(term_profile) || max(term_profile) < threshold
        TCD = NaN;
        return;
    end
    
    TCD = interp1(term_profile, distances_bp, threshold, 'linear', 'extrap');
end

% Single gene simulation function - exact copy from generate_gene_length_grid.m line 277
function [R_sol, REH_sol, avg_E_bound] = run_single_gene_simulation_local(P)
    % Run single gene simulation similar to existing scripts
    global N PAS N_PAS Ef_ss;
    syms Ef real;
    
    % Set up geometry
    L_a = P.L_a;
    N = floor(P.geneLength_bp / L_a);
    PAS = floor(P.PASposition / L_a);
    N_PAS = N - PAS + 1;
    SD = floor(P.SD_bp / L_a);  % Saturation distance in nodes
    
    % Compute steady states
    [r_E_BeforePas] = compute_steady_states(P, P.EBindingNumber + 1);
    
    % Set up kPon values based on the saturation distance concept
    kPon_vals = zeros(1, PAS);
    
    if PAS <= SD
        % Case 1: PAS is before or at saturation distance
        % Linear increase from kPon_min to value at PAS (never reaches kPon_max)
        kPon_vals = linspace(P.kPon_min, P.kPon_min + (P.kPon_max - P.kPon_min) * (PAS / SD), PAS);
    else
        % Case 2: PAS is beyond saturation distance
        if P.kPon_option == 1
            % Option 1: Saturate at SD, then constant
            kPon_vals(1:SD) = linspace(P.kPon_min, P.kPon_max, SD);
            kPon_vals((SD+1):PAS) = P.kPon_max;  % Constant at kPon_max
        else
            % Option 2: Continue linear increase after SD
            kPon_vals(1:SD) = linspace(P.kPon_min, P.kPon_max, SD);
            % Continue linear increase beyond kPon_max
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
            kPon_val = kPon_vals(end);  % Use the kPon value at PAS for after-PAS region
            RE_vals(e, idx) = subs(r_E_BeforePas(e), {'kPon', 'kPoff'}, {kPon_val, P.kPoff_const});
        end
    end
    
    P.RE_val_bind_E = matlabFunction(simplify(sum(sym(1:P.EBindingNumber)' .* RE_vals(2:end, :), 1)), 'Vars', {Ef});
    
    % Solve system
    X_guess = 1e-6 * ones(N + N_PAS, 1);
    options = optimoptions('fsolve', 'Display', 'off', 'FunctionTolerance', 1e-8);
    
    % Two-step solution
    P.FirstRun = true;
    P.is_unphysical = false;
    Ef_ss = 0;
    
    X_base = fsolve(@(xx) ode_dynamics_multipleE(xx, P), X_guess, options);
    if P.is_unphysical
        error('Unphysical result in step 1');
    end
    
    % Update kHon and resolve
    avg_E_bound = P.RE_val_bind_E(Ef_ss);
    P.FirstRun = false;
    P.kHon = P.kHon * avg_E_bound(end);
    X_final = fsolve(@(xx) ode_dynamics_multipleE(xx, P), X_base, options);
    
    R_sol = X_final(1:N);
    REH_sol = X_final((N+1):(N+N_PAS));
end

