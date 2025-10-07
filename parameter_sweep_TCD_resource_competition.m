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

% Two-parameter sweep: kPon_max and SD_bp
kPon_max_values = [0.1, 0.5, 1, 2, 5, 10];  % Range from 0.1 to 10
SD_bp_values = [4000, 10000, 20000, 40000, 70000, 100000];  % Range from 4 kb to 100 kb

n_kPon = length(kPon_max_values);
n_SD = length(SD_bp_values);
n_param_combinations = n_kPon + n_SD;  % Total number of sweeps

fprintf('  Two-parameter sweep:\n');
fprintf('  Parameter 1: kPon_max\n');
fprintf('    Values: ');
for i = 1:n_kPon
    fprintf('%.2g', kPon_max_values(i));
    if i < n_kPon
        fprintf(', ');
    end
end
fprintf('\n');
fprintf('  Parameter 2: SD_bp\n');
fprintf('    Values: ');
for i = 1:n_SD
    fprintf('%.0f', SD_bp_values(i)/1000);
    if i < n_SD
        fprintf(', ');
    end
end
fprintf(' kb\n');

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
P_base.kPon_max = 1;
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
output_dir = 'SecondVersionResults/TCD_ParameterSweep_ResourceCompetition/kPon_SD_2D/';
if ~exist(output_dir, 'dir')
    mkdir(output_dir);
end
fprintf('Output directory: %s\n', output_dir);

%% --- PARAMETER SWEEP LOOP ---
fprintf('\n=== Starting Parameter Sweep ===\n');
fprintf('This will run the full three-step pipeline for each parameter combination.\n');
fprintf('Total combinations: %d (kPon_max: %d + SD_bp: %d)\n', n_param_combinations, n_kPon, n_SD);
fprintf('Estimated time: ~5-15 minutes per parameter value\n\n');

% Pre-allocate results storage
sweep_results = struct();
sweep_results.metadata.creation_date = datestr(now);
sweep_results.metadata.sweep_type = '2D: kPon_max and SD_bp';
sweep_results.metadata.kPon_max_values = kPon_max_values;
sweep_results.metadata.SD_bp_values = SD_bp_values;
sweep_results.metadata.base_parameters = P_base;
sweep_results.metadata.R_total_target = R_total_target;
sweep_results.metadata.E_total_target = E_total_target;

% Storage for results from each parameter sweep
sweep_results.kPon_sweep = cell(n_kPon, 1);
sweep_results.SD_sweep = cell(n_SD, 1);

% Storage for summary data
% For kPon_max sweep (varying kPon_max, SD_bp fixed at base value)
kPon_R_free = zeros(n_kPon, 1);
kPon_E_free = zeros(n_kPon, 1);
kPon_TCD_matrix = zeros(n_kPon, n_gene_lengths);
kPon_success = ones(n_kPon, 1);

% For SD_bp sweep (varying SD_bp, kPon_max fixed at base value)
SD_R_free = zeros(n_SD, 1);
SD_E_free = zeros(n_SD, 1);
SD_TCD_matrix = zeros(n_SD, n_gene_lengths);
SD_success = ones(n_SD, 1);

total_sweep_time = tic;

%% SWEEP 1: kPon_max (SD_bp fixed at base value)
fprintf('\n╔══════════════════════════════════════════════════════════╗\n');
fprintf('║ SWEEP 1: kPon_max (SD_bp = %.0f bp)                    ║\n', P_base.SD_bp);
fprintf('╚══════════════════════════════════════════════════════════╝\n');

for k_idx = 1:n_kPon
    kPon_value = kPon_max_values(k_idx);
    
    fprintf('\n+===========================================================+\n');
    fprintf('| kPon_max %d/%d: kPon_max = %.2g %-25s|\n', k_idx, n_kPon, kPon_value, '');
    fprintf('+===========================================================+\n\n');
    
    param_time = tic;
    
    try
        % Update parameter
        P_sweep = P_base;
        P_sweep.kPon_max = kPon_value;
        
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
        sweep_results.kPon_sweep{k_idx} = struct();
        sweep_results.kPon_sweep{k_idx}.kPon_max = kPon_value;
        sweep_results.kPon_sweep{k_idx}.grid_data = grid_data;
        sweep_results.kPon_sweep{k_idx}.interp_data = interp_data;
        sweep_results.kPon_sweep{k_idx}.tcd_data = tcd_data;
        
        % Store summary data
        kPon_R_free(k_idx) = tcd_data.R_free_solution;
        kPon_E_free(k_idx) = tcd_data.E_free_solution;
        kPon_TCD_matrix(k_idx, :) = tcd_data.TCD_values;
        kPon_success(k_idx) = 1;
        
        fprintf('\n✓ kPon_max = %.2g completed in %.1f minutes\n', kPon_value, toc(param_time)/60);
        
    catch ME
        fprintf('\n✗ ERROR for kPon_max = %.2g: %s\n', kPon_value, ME.message);
        fprintf('  Stack trace:\n');
        for k = 1:length(ME.stack)
            fprintf('    %s (line %d)\n', ME.stack(k).name, ME.stack(k).line);
        end
        
        kPon_success(k_idx) = 0;
        kPon_R_free(k_idx) = NaN;
        kPon_E_free(k_idx) = NaN;
        kPon_TCD_matrix(k_idx, :) = NaN;
    end
end

%% SWEEP 2: SD_bp (kPon_max fixed at base value)
fprintf('\n╔══════════════════════════════════════════════════════════╗\n');
fprintf('║ SWEEP 2: SD_bp (kPon_max = %.2g)                       ║\n', P_base.kPon_max);
fprintf('╚══════════════════════════════════════════════════════════╝\n');

for s_idx = 1:n_SD
    SD_value = SD_bp_values(s_idx);
    
    fprintf('\n+===========================================================+\n');
    fprintf('| SD_bp %d/%d: SD = %.0f kb %-29s|\n', s_idx, n_SD, SD_value/1000, '');
    fprintf('+===========================================================+\n\n');
    
    param_time = tic;
    
    try
        % Update parameter
        P_sweep = P_base;
        P_sweep.SD_bp = SD_value;
        
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
        sweep_results.SD_sweep{s_idx} = struct();
        sweep_results.SD_sweep{s_idx}.SD_bp = SD_value;
        sweep_results.SD_sweep{s_idx}.grid_data = grid_data;
        sweep_results.SD_sweep{s_idx}.interp_data = interp_data;
        sweep_results.SD_sweep{s_idx}.tcd_data = tcd_data;
        
        % Store summary data
        SD_R_free(s_idx) = tcd_data.R_free_solution;
        SD_E_free(s_idx) = tcd_data.E_free_solution;
        SD_TCD_matrix(s_idx, :) = tcd_data.TCD_values;
        SD_success(s_idx) = 1;
        
        fprintf('\n✓ SD_bp = %.0f kb completed in %.1f minutes\n', SD_value/1000, toc(param_time)/60);
        
    catch ME
        fprintf('\n✗ ERROR for SD_bp = %.0f kb: %s\n', SD_value/1000, ME.message);
        fprintf('  Stack trace:\n');
        for k = 1:length(ME.stack)
            fprintf('    %s (line %d)\n', ME.stack(k).name, ME.stack(k).line);
        end
        
        SD_success(s_idx) = 0;
        SD_R_free(s_idx) = NaN;
        SD_E_free(s_idx) = NaN;
        SD_TCD_matrix(s_idx, :) = NaN;
    end
end

total_time = toc(total_sweep_time);

fprintf('\n+===========================================================+\n');
fprintf('| PARAMETER SWEEP COMPLETED                                 |\n');
fprintf('+===========================================================+\n');
fprintf('Total time: %.1f minutes\n', total_time/60);
fprintf('Successful sweeps:\n');
fprintf('  kPon_max: %d/%d\n', sum(kPon_success), n_kPon);
fprintf('  SD_bp: %d/%d\n', sum(SD_success), n_SD);

%% --- SUMMARY ANALYSIS ---
fprintf('\n=== Analysis Summary ===\n');

% Add summary data to results
sweep_results.summary.kPon_max.R_free = kPon_R_free;
sweep_results.summary.kPon_max.E_free = kPon_E_free;
sweep_results.summary.kPon_max.TCD_matrix = kPon_TCD_matrix;
sweep_results.summary.kPon_max.success = kPon_success;

sweep_results.summary.SD_bp.R_free = SD_R_free;
sweep_results.summary.SD_bp.E_free = SD_E_free;
sweep_results.summary.SD_bp.TCD_matrix = SD_TCD_matrix;
sweep_results.summary.SD_bp.success = SD_success;

sweep_results.summary.gene_lengths_kb = gene_lengths_kb;
sweep_results.summary.TCD_threshold = TCD_threshold;
sweep_results.summary.total_time_minutes = total_time/60;

% Analyze kPon_max trends
kPon_valid = kPon_success == 1;
if sum(kPon_valid) > 1
    fprintf('\nkPon_max sweep trends:\n');
    corr_R = corr(kPon_max_values(kPon_valid)', kPon_R_free(kPon_valid));
    corr_E = corr(kPon_max_values(kPon_valid)', kPon_E_free(kPon_valid));
    fprintf('  R_free vs kPon_max: r = %.3f\n', corr_R);
    fprintf('  E_free vs kPon_max: r = %.3f\n', corr_E);
end

% Analyze SD_bp trends
SD_valid = SD_success == 1;
if sum(SD_valid) > 1
    fprintf('\nSD_bp sweep trends:\n');
    corr_R = corr(SD_bp_values(SD_valid)', SD_R_free(SD_valid));
    corr_E = corr(SD_bp_values(SD_valid)', SD_E_free(SD_valid));
    fprintf('  R_free vs SD_bp: r = %.3f\n', corr_R);
    fprintf('  E_free vs SD_bp: r = %.3f\n', corr_E);
end

%% --- VISUALIZATION ---
fprintf('\nGenerating visualizations...\n');

% Create two-panel figure showing TCD vs gene length for each parameter
fig = figure('Position', [50, 50, 1400, 600]);

% Panel 1: TCD vs gene length for different kPon_max values
subplot(1, 2, 1);
colors_kPon = parula(n_kPon);
hold on;
for k_idx = 1:n_kPon
    if kPon_success(k_idx)
        valid_tcd = ~isnan(kPon_TCD_matrix(k_idx, :));
        if sum(valid_tcd) > 0
            semilogx(gene_lengths_kb(valid_tcd), kPon_TCD_matrix(k_idx, valid_tcd), ...
                'o-', 'LineWidth', 2.5, 'MarkerSize', 8, 'Color', colors_kPon(k_idx, :), ...
                'DisplayName', sprintf('k_{Pon,max} = %.2g', kPon_max_values(k_idx)));
        end
    end
end
xlabel('TSS-to-PAS Distance (kb)', 'FontSize', 13, 'FontWeight', 'bold');
ylabel(sprintf('TCD at %.0f%% (bp)', TCD_threshold*100), 'FontSize', 13, 'FontWeight', 'bold');
title(sprintf('Effect of k_{Pon,max} (SD = %.0f kb)', P_base.SD_bp/1000), ...
    'FontSize', 14, 'FontWeight', 'bold');
legend('Location', 'best', 'FontSize', 10);
grid on;
set(gca, 'FontSize', 11);
hold off;

% Panel 2: TCD vs gene length for different SD_bp values
subplot(1, 2, 2);
colors_SD = parula(n_SD);
hold on;
for s_idx = 1:n_SD
    if SD_success(s_idx)
        valid_tcd = ~isnan(SD_TCD_matrix(s_idx, :));
        if sum(valid_tcd) > 0
            semilogx(gene_lengths_kb(valid_tcd), SD_TCD_matrix(s_idx, valid_tcd), ...
                'o-', 'LineWidth', 2.5, 'MarkerSize', 8, 'Color', colors_SD(s_idx, :), ...
                'DisplayName', sprintf('SD = %.0f kb', SD_bp_values(s_idx)/1000));
        end
    end
end
xlabel('TSS-to-PAS Distance (kb)', 'FontSize', 13, 'FontWeight', 'bold');
ylabel(sprintf('TCD at %.0f%% (bp)', TCD_threshold*100), 'FontSize', 13, 'FontWeight', 'bold');
title(sprintf('Effect of SD (k_{Pon,max} = %.2g)', P_base.kPon_max), ...
    'FontSize', 14, 'FontWeight', 'bold');
legend('Location', 'best', 'FontSize', 10);
grid on;
set(gca, 'FontSize', 11);
hold off;

sgtitle('TCD vs Gene Length: Two-Parameter Sweep with Resource Competition', ...
    'FontSize', 16, 'FontWeight', 'bold');

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

fprintf(fid, '%% Two-Parameter Sweep with Resource Competition\n');
fprintf(fid, '%% Parameters: kPon_max and SD_bp\n');
fprintf(fid, '%% Generated: %s\n', datestr(now));
fprintf(fid, '%% Total time: %.1f minutes\n', total_time/60);
fprintf(fid, '%% Success rates: kPon_max=%d/%d, SD_bp=%d/%d\n', ...
    sum(kPon_success), n_kPon, sum(SD_success), n_SD);
fprintf(fid, '%%\n');

fprintf(fid, '%% ===== SWEEP 1: kPon_max (SD_bp fixed at %.0f bp) =====\n', P_base.SD_bp);
fprintf(fid, '%% RESOURCE POOL SOLUTIONS\n');
fprintf(fid, '%% kPon_max\tR_free\tE_free\tR_free_pct\tE_free_pct\n');
for k_idx = 1:n_kPon
    if kPon_success(k_idx)
        fprintf(fid, '%.6g\t%.2f\t%.2f\t%.2f\t%.2f\n', ...
            kPon_max_values(k_idx), ...
            kPon_R_free(k_idx), kPon_E_free(k_idx), ...
            kPon_R_free(k_idx)/R_total_target*100, ...
            kPon_E_free(k_idx)/E_total_target*100);
    end
end

fprintf(fid, '%%\n');
fprintf(fid, '%% TCD VALUES (bp) at %.0f%% threshold - kPon_max sweep\n', TCD_threshold*100);
fprintf(fid, '%% kPon_max');
for g_idx = 1:n_gene_lengths
    fprintf(fid, '\tL%.0fkb', gene_lengths_kb(g_idx));
end
fprintf(fid, '\n');

for k_idx = 1:n_kPon
    fprintf(fid, '%.6g', kPon_max_values(k_idx));
    for g_idx = 1:n_gene_lengths
        if kPon_success(k_idx) && ~isnan(kPon_TCD_matrix(k_idx, g_idx))
            fprintf(fid, '\t%.2f', kPon_TCD_matrix(k_idx, g_idx));
        else
            fprintf(fid, '\tNaN');
        end
    end
    fprintf(fid, '\n');
end

fprintf(fid, '%%\n');
fprintf(fid, '%% ===== SWEEP 2: SD_bp (kPon_max fixed at %.2g) =====\n', P_base.kPon_max);
fprintf(fid, '%% RESOURCE POOL SOLUTIONS\n');
fprintf(fid, '%% SD_bp\tR_free\tE_free\tR_free_pct\tE_free_pct\n');
for s_idx = 1:n_SD
    if SD_success(s_idx)
        fprintf(fid, '%.0f\t%.2f\t%.2f\t%.2f\t%.2f\n', ...
            SD_bp_values(s_idx), ...
            SD_R_free(s_idx), SD_E_free(s_idx), ...
            SD_R_free(s_idx)/R_total_target*100, ...
            SD_E_free(s_idx)/E_total_target*100);
    end
end

fprintf(fid, '%%\n');
fprintf(fid, '%% TCD VALUES (bp) at %.0f%% threshold - SD_bp sweep\n', TCD_threshold*100);
fprintf(fid, '%% SD_bp');
for g_idx = 1:n_gene_lengths
    fprintf(fid, '\tL%.0fkb', gene_lengths_kb(g_idx));
end
fprintf(fid, '\n');

for s_idx = 1:n_SD
    fprintf(fid, '%.0f', SD_bp_values(s_idx));
    for g_idx = 1:n_gene_lengths
        if SD_success(s_idx) && ~isnan(SD_TCD_matrix(s_idx, g_idx))
            fprintf(fid, '\t%.2f', SD_TCD_matrix(s_idx, g_idx));
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

fprintf('Two-parameter sweep: kPon_max and SD_bp\n');
fprintf('kPon_max range: %.2g to %.2g\n', kPon_max_values(1), kPon_max_values(end));
fprintf('SD_bp range: %.0f to %.0f bp\n', SD_bp_values(1), SD_bp_values(end));
fprintf('Success rates:\n');
fprintf('  kPon_max sweep: %d/%d (%.1f%%)\n', sum(kPon_success), n_kPon, ...
    sum(kPon_success)/n_kPon*100);
fprintf('  SD_bp sweep: %d/%d (%.1f%%)\n', sum(SD_success), n_SD, ...
    sum(SD_success)/n_SD*100);
fprintf('Total computation time: %.1f minutes\n', total_time/60);

if sum(kPon_valid) > 1
    fprintf('\nkPon_max sweep insights:\n');
    kPon_R_range = [min(kPon_R_free(kPon_valid)), max(kPon_R_free(kPon_valid))];
    kPon_E_range = [min(kPon_E_free(kPon_valid)), max(kPon_E_free(kPon_valid))];
    fprintf('  • R_free range: %.0f to %.0f (%.1f%% to %.1f%%)\n', ...
        kPon_R_range(1), kPon_R_range(2), ...
        kPon_R_range(1)/R_total_target*100, kPon_R_range(2)/R_total_target*100);
    fprintf('  • E_free range: %.0f to %.0f (%.1f%% to %.1f%%)\n', ...
        kPon_E_range(1), kPon_E_range(2), ...
        kPon_E_range(1)/E_total_target*100, kPon_E_range(2)/E_total_target*100);
    
    valid_tcd = kPon_TCD_matrix(kPon_valid, :);
    valid_tcd = valid_tcd(~isnan(valid_tcd));
    if ~isempty(valid_tcd)
        fprintf('  • TCD range: %.0f to %.0f bp\n', min(valid_tcd), max(valid_tcd));
    end
end

if sum(SD_valid) > 1
    fprintf('\nSD_bp sweep insights:\n');
    SD_R_range = [min(SD_R_free(SD_valid)), max(SD_R_free(SD_valid))];
    SD_E_range = [min(SD_E_free(SD_valid)), max(SD_E_free(SD_valid))];
    fprintf('  • R_free range: %.0f to %.0f (%.1f%% to %.1f%%)\n', ...
        SD_R_range(1), SD_R_range(2), ...
        SD_R_range(1)/R_total_target*100, SD_R_range(2)/R_total_target*100);
    fprintf('  • E_free range: %.0f to %.0f (%.1f%% to %.1f%%)\n', ...
        SD_E_range(1), SD_E_range(2), ...
        SD_E_range(1)/E_total_target*100, SD_E_range(2)/E_total_target*100);
    
    valid_tcd = SD_TCD_matrix(SD_valid, :);
    valid_tcd = valid_tcd(~isnan(valid_tcd));
    if ~isempty(valid_tcd)
        fprintf('  • TCD range: %.0f to %.0f bp\n', min(valid_tcd), max(valid_tcd));
    end
end

fprintf('\nOutput files:\n');
fprintf('  %s\n', mat_filename);
fprintf('  %s\n', fig_filename);
fprintf('  %s\n', txt_filename);

fprintf('\nTwo-parameter sweep with resource competition completed successfully!\n');

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

