% PARAMETER_SWEEP_TCD.m
% Parameter analysis for Termination Commitment Distance (TCD)
%
% This script performs parameter sweeps to analyze how different model parameters
% affect the TCD across different gene lengths. It uses the same framework as
% the gene length analysis but varies one or two parameters at a time.
%
% Key analyses:
% - How does saturation distance (SD) affect TCD?
% - How does kPon_max affect TCD?
% - How do other parameters (kHon, kHoff, kc, etc.) affect TCD?
% - TCD sensitivity analysis across gene lengths

fprintf('=== Parameter Sweep for TCD Analysis ===\n');
fprintf('Analyzing how model parameters affect Termination Commitment Distance...\n\n');

%% --- CONFIGURATION ---
fprintf('Configuration:\n');

% Choose parameter to sweep
% Options: 'SD_bp', 'kPon_max', 'kHon', 'kHoff', 'kc', 'kEon', 'EBindingNumber'
parameter_to_sweep = 'SD_bp';  % Change this to sweep different parameters

% Sweep configuration
switch parameter_to_sweep
    case 'SD_bp'
        param_values = [5000, 10000, 15000, 20000, 30000, 50000];  % Different saturation distances
        param_label = 'Saturation Distance (kb)';
        param_display = @(x) x/1000;  % Convert to kb for display
        
    case 'kPon_max'
        param_values = [0.05, 0.1, 0.2, 0.5, 1.0, 2.0];
        param_label = 'k_{Pon,max}';
        param_display = @(x) x;
        
    case 'kHon'
        param_values = [0.05, 0.1, 0.2, 0.4, 0.8, 1.6];
        param_label = 'k_{Hon}';
        param_display = @(x) x;
        
    case 'kHoff'
        param_values = [0.005, 0.00625, 0.0125, 0.025, 0.05];
        param_label = 'k_{Hoff}';
        param_display = @(x) x;
        
    case 'kc'
        param_values = [0.01, 0.025, 0.05, 0.1, 0.2, 0.4];
        param_label = 'k_c (cleavage rate)';
        param_display = @(x) x;
        
    case 'kEon'
        param_values = [0.0001, 0.00025, 0.0005, 0.001, 0.002];
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

% TCD threshold (typically 50% or 75%)
TCD_threshold = 0.5;  % 50% termination threshold
fprintf('  TCD threshold: %.0f%%\n', TCD_threshold*100);

% Parallel computation flag
use_parallel = true;
if use_parallel && isempty(gcp('nocreate'))
    parpool;
end

%% --- BASE PARAMETERS ---
fprintf('\nSetting up base parameters...\n');

% Standard parameter set (will be modified during sweep)
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
P_base.SD_bp = 20000;  % Default saturation distance
P_base.kPon_option = 1;  % Saturate at SD
P_base.EBindingNumber = 1;
P_base.after_PAS_length = 5000;  % Fixed after-PAS length

% Resource pools (typical values)
P_base.Pol_total = 70000;
P_base.E_total = 70000;

fprintf('Base parameters set.\n');
fprintf('  Default SD: %.0f kb\n', P_base.SD_bp/1000);
fprintf('  Default kPon_max: %.3f\n', P_base.kPon_max);
fprintf('  Default kHon: %.3f\n', P_base.kHon);

%% --- PARAMETER SWEEP ---
fprintf('\nStarting parameter sweep...\n');
fprintf('Total simulations: %d parameter values × %d gene lengths = %d\n', ...
    n_param_values, n_gene_lengths, n_param_values * n_gene_lengths);

% Pre-allocate results
TCD_matrix = zeros(n_param_values, n_gene_lengths);
success_matrix = ones(n_param_values, n_gene_lengths);  % Track successful simulations

tic;

for p_idx = 1:n_param_values
    param_value = param_values(p_idx);
    
    fprintf('\n--- Parameter %d/%d: %s = %.3g ---\n', p_idx, n_param_values, ...
        parameter_to_sweep, param_display(param_value));
    
    % Create parameter set for this sweep value
    P_sweep = P_base;
    P_sweep.(parameter_to_sweep) = param_value;
    
    % Parallel loop over gene lengths
    TCD_row = zeros(1, n_gene_lengths);
    success_row = ones(1, n_gene_lengths);
    
    parfor (g_idx = 1:n_gene_lengths, use_parallel * 8)  % Use 8 workers if parallel
        try
            % Set up gene-specific parameters
            P_gene = P_sweep;
            P_gene.PASposition = gene_lengths_bp(g_idx);  % TSS-to-PAS distance
            P_gene.geneLength_bp = gene_lengths_bp(g_idx) + P_gene.after_PAS_length;
            
            % Run simulation and calculate TCD
            [R_sol, REH_sol, P_sim] = run_single_gene_simulation_TCD_sweep(P_gene);
            
            % Calculate termination profile using flux-based approach
            [exit_cdf, distances_bp] = calculate_pas_usage_profile(R_sol, REH_sol, P_sim);
            
            % Find TCD at specified threshold
            if max(exit_cdf) < TCD_threshold
                TCD_row(g_idx) = -1;  % Insufficient termination
            else
                TCD_row(g_idx) = interp1(exit_cdf, distances_bp, TCD_threshold, 'linear', 'extrap');
            end
            
        catch ME
            % Handle failed simulations
            fprintf('Warning: Failed for %s=%.3g, L=%.0f kb: %s\n', ...
                parameter_to_sweep, param_display(param_value), gene_lengths_kb(g_idx), ME.message);
            TCD_row(g_idx) = NaN;
            success_row(g_idx) = 0;
        end
    end
    
    % Store results
    TCD_matrix(p_idx, :) = TCD_row;
    success_matrix(p_idx, :) = success_row;
    
    % Progress report
    successful = sum(success_row);
    fprintf('  Completed: %d/%d successful\n', successful, n_gene_lengths);
end

sweep_time = toc;
fprintf('\nParameter sweep completed in %.1f minutes\n', sweep_time/60);

%% --- ANALYSIS ---
fprintf('\nAnalyzing results...\n');

% Calculate statistics
valid_count = sum(success_matrix(:));
total_count = numel(success_matrix);
success_rate = valid_count / total_count * 100;

fprintf('Success rate: %d/%d (%.1f%%)\n', valid_count, total_count, success_rate);

% Find trends
fprintf('\nTrends:\n');
for g_idx = 1:n_gene_lengths
    valid_idx = ~isnan(TCD_matrix(:, g_idx)) & TCD_matrix(:, g_idx) > 0;
    if sum(valid_idx) > 1
        % Calculate correlation between parameter and TCD
        corr_val = corr(param_values(valid_idx)', TCD_matrix(valid_idx, g_idx));
        fprintf('  L=%.0f kb: correlation = %.3f', gene_lengths_kb(g_idx), corr_val);
        
        % Describe trend
        if abs(corr_val) > 0.7
            if corr_val > 0
                fprintf(' (strong positive trend)');
            else
                fprintf(' (strong negative trend)');
            end
        elseif abs(corr_val) > 0.4
            if corr_val > 0
                fprintf(' (moderate positive trend)');
            else
                fprintf(' (moderate negative trend)');
            end
        else
            fprintf(' (weak/no trend)');
        end
        fprintf('\n');
    end
end

%% --- VISUALIZATION ---
fprintf('\nGenerating visualizations...\n');

% Create figure with multiple subplots
figure('Position', [100, 100, 1400, 900]);

% Plot 1: TCD vs Parameter for different gene lengths
subplot(2, 2, 1);
colors = jet(n_gene_lengths);
hold on;
for g_idx = 1:n_gene_lengths
    valid_idx = ~isnan(TCD_matrix(:, g_idx)) & TCD_matrix(:, g_idx) > 0;
    if sum(valid_idx) > 0
        plot(param_display(param_values(valid_idx)), TCD_matrix(valid_idx, g_idx), ...
            'o-', 'LineWidth', 2, 'MarkerSize', 8, 'Color', colors(g_idx, :), ...
            'DisplayName', sprintf('L = %.0f kb', gene_lengths_kb(g_idx)));
    end
end
xlabel(param_label, 'FontSize', 12);
ylabel(sprintf('TCD at %.0f%% threshold (bp)', TCD_threshold*100), 'FontSize', 12);
title(sprintf('TCD vs %s', parameter_to_sweep), 'FontSize', 14, 'FontWeight', 'bold');
legend('Location', 'best', 'FontSize', 10);
grid on;
hold off;

% Plot 2: TCD vs Gene Length for different parameter values
subplot(2, 2, 2);
colors_param = parula(n_param_values);
hold on;
for p_idx = 1:n_param_values
    valid_idx = ~isnan(TCD_matrix(p_idx, :)) & TCD_matrix(p_idx, :) > 0;
    if sum(valid_idx) > 0
        semilogx(gene_lengths_kb(valid_idx), TCD_matrix(p_idx, valid_idx), ...
            'o-', 'LineWidth', 2, 'MarkerSize', 8, 'Color', colors_param(p_idx, :), ...
            'DisplayName', sprintf('%s = %.3g', parameter_to_sweep, param_display(param_values(p_idx))));
    end
end
xlabel('TSS-to-PAS Distance (kb)', 'FontSize', 12);
ylabel(sprintf('TCD at %.0f%% threshold (bp)', TCD_threshold*100), 'FontSize', 12);
title(sprintf('TCD vs Gene Length'), 'FontSize', 14, 'FontWeight', 'bold');
legend('Location', 'best', 'FontSize', 9);
grid on;
hold off;

% Plot 3: Heatmap of TCD values
subplot(2, 2, 3);
TCD_plot = TCD_matrix;
TCD_plot(TCD_plot < 0) = NaN;  % Remove insufficient termination cases
imagesc(1:n_gene_lengths, 1:n_param_values, TCD_plot);
colorbar;
xlabel('Gene Length Index', 'FontSize', 12);
ylabel('Parameter Value Index', 'FontSize', 12);
title(sprintf('TCD Heatmap (bp)'), 'FontSize', 14, 'FontWeight', 'bold');
set(gca, 'XTick', 1:n_gene_lengths, 'XTickLabel', arrayfun(@(x) sprintf('%.0f', x), gene_lengths_kb, 'UniformOutput', false));
set(gca, 'YTick', 1:n_param_values, 'YTickLabel', arrayfun(@(x) sprintf('%.3g', x), param_display(param_values), 'UniformOutput', false));
set(gca, 'FontSize', 10);

% Plot 4: Normalized TCD (relative to baseline)
subplot(2, 2, 4);
% Use middle parameter value as baseline
baseline_idx = ceil(n_param_values / 2);
TCD_normalized = TCD_matrix ./ TCD_matrix(baseline_idx, :);
TCD_normalized(isinf(TCD_normalized) | isnan(TCD_normalized)) = 1;

hold on;
for g_idx = 1:n_gene_lengths
    valid_idx = ~isnan(TCD_normalized(:, g_idx));
    if sum(valid_idx) > 0
        plot(param_display(param_values(valid_idx)), TCD_normalized(valid_idx, g_idx), ...
            'o-', 'LineWidth', 2, 'MarkerSize', 8, 'Color', colors(g_idx, :), ...
            'DisplayName', sprintf('L = %.0f kb', gene_lengths_kb(g_idx)));
    end
end
yline(1, 'k--', 'LineWidth', 1.5, 'DisplayName', 'Baseline');
xlabel(param_label, 'FontSize', 12);
ylabel('Normalized TCD (relative to baseline)', 'FontSize', 12);
title(sprintf('Normalized TCD vs %s', parameter_to_sweep), 'FontSize', 14, 'FontWeight', 'bold');
legend('Location', 'best', 'FontSize', 9);
grid on;
hold off;

%% --- SAVE RESULTS ---
fprintf('\nSaving results...\n');

% Create output directory
output_dir = 'SecondVersionResults/TCD_ParameterSweep/';
if ~exist(output_dir, 'dir')
    mkdir(output_dir);
end

% Create results structure
results = struct();
results.metadata.creation_date = datestr(now);
results.metadata.sweep_time_minutes = sweep_time/60;
results.metadata.description = sprintf('Parameter sweep analysis for TCD: varying %s', parameter_to_sweep);

% Parameter info
results.sweep.parameter_name = parameter_to_sweep;
results.sweep.parameter_values = param_values;
results.sweep.parameter_label = param_label;
results.sweep.gene_lengths_kb = gene_lengths_kb;
results.sweep.gene_lengths_bp = gene_lengths_bp;
results.sweep.TCD_threshold = TCD_threshold;

% Base parameters
results.base_parameters = P_base;

% Results data
results.data.TCD_matrix = TCD_matrix;
results.data.success_matrix = success_matrix;
results.data.success_rate = success_rate;

% Statistics
results.statistics.correlations = zeros(n_gene_lengths, 1);
for g_idx = 1:n_gene_lengths
    valid_idx = ~isnan(TCD_matrix(:, g_idx)) & TCD_matrix(:, g_idx) > 0;
    if sum(valid_idx) > 1
        results.statistics.correlations(g_idx) = corr(param_values(valid_idx)', TCD_matrix(valid_idx, g_idx));
    else
        results.statistics.correlations(g_idx) = NaN;
    end
end

% Save MATLAB file
timestamp = datestr(now, 'yyyymmdd_HHMMSS');
mat_filename = fullfile(output_dir, sprintf('TCD_sweep_%s_%s.mat', parameter_to_sweep, timestamp));
save(mat_filename, 'results', '-v7.3');

% Save figure
fig_filename = fullfile(output_dir, sprintf('TCD_sweep_%s_%s.png', parameter_to_sweep, timestamp));
saveas(gcf, fig_filename);

% Save text summary
txt_filename = fullfile(output_dir, sprintf('TCD_sweep_%s_%s.txt', parameter_to_sweep, timestamp));
fid = fopen(txt_filename, 'w');

fprintf(fid, '%% TCD Parameter Sweep Analysis\n');
fprintf(fid, '%% Generated on: %s\n', results.metadata.creation_date);
fprintf(fid, '%% Parameter swept: %s\n', parameter_to_sweep);
fprintf(fid, '%% TCD threshold: %.0f%%\n', TCD_threshold*100);
fprintf(fid, '%% \n');
fprintf(fid, '%% Base Parameters:\n');
fprintf(fid, '%% L_a = %g, k_in = %g, kEon = %g, kEoff = %g\n', P_base.L_a, P_base.k_in, P_base.kEon, P_base.kEoff);
fprintf(fid, '%% k_e = %g, k_e2 = %g, kHon = %g, kHoff = %g\n', P_base.k_e, P_base.k_e2, P_base.kHon, P_base.kHoff);
fprintf(fid, '%% kc = %g, SD_bp = %g, kPon_max = %g\n', P_base.kc, P_base.SD_bp, P_base.kPon_max);
fprintf(fid, '%% EBindingNumber = %g\n', P_base.EBindingNumber);
fprintf(fid, '%% \n');
fprintf(fid, '%% Results Summary:\n');
fprintf(fid, '%% Success rate: %.1f%%\n', success_rate);
fprintf(fid, '%% \n');

% Write data table
fprintf(fid, '%% TCD Data (rows = parameter values, columns = gene lengths)\n');
fprintf(fid, '%% Parameter values: ');
for i = 1:n_param_values
    fprintf(fid, '%.6g ', param_display(param_values(i)));
end
fprintf(fid, '\n');

fprintf(fid, '%% Gene lengths (kb): ');
for i = 1:n_gene_lengths
    fprintf(fid, '%.0f ', gene_lengths_kb(i));
end
fprintf(fid, '\n');
fprintf(fid, '%% \n');

% Header
fprintf(fid, 'ParamValue\t');
for i = 1:n_gene_lengths
    fprintf(fid, 'L%.0fkb', gene_lengths_kb(i));
    if i < n_gene_lengths
        fprintf(fid, '\t');
    end
end
fprintf(fid, '\n');

% Data rows
for p_idx = 1:n_param_values
    fprintf(fid, '%.6g\t', param_display(param_values(p_idx)));
    for g_idx = 1:n_gene_lengths
        fprintf(fid, '%.2f', TCD_matrix(p_idx, g_idx));
        if g_idx < n_gene_lengths
            fprintf(fid, '\t');
        end
    end
    fprintf(fid, '\n');
end

fclose(fid);

fprintf('Results saved:\n');
fprintf('  MATLAB file: %s\n', mat_filename);
fprintf('  Figure: %s\n', fig_filename);
fprintf('  Text file: %s\n', txt_filename);

%% --- FINAL SUMMARY ---
fprintf('\n=== PARAMETER SWEEP COMPLETE ===\n');
fprintf('Parameter: %s\n', parameter_to_sweep);
fprintf('Range: %.3g to %.3g\n', param_display(param_values(1)), param_display(param_values(end)));
fprintf('Gene lengths: %.0f to %.0f kb\n', gene_lengths_kb(1), gene_lengths_kb(end));
fprintf('Success rate: %.1f%%\n', success_rate);
fprintf('\nKey findings:\n');

% Report strongest correlations
[~, sorted_idx] = sort(abs(results.statistics.correlations), 'descend');
for i = 1:min(3, n_gene_lengths)
    g_idx = sorted_idx(i);
    if ~isnan(results.statistics.correlations(g_idx))
        fprintf('  L=%.0f kb: r=%.3f with %s\n', ...
            gene_lengths_kb(g_idx), results.statistics.correlations(g_idx), parameter_to_sweep);
    end
end

fprintf('\nParameter sweep completed successfully!\n');

%% --- HELPER FUNCTIONS ---

function [R_sol, REH_sol, P_sim] = run_single_gene_simulation_TCD_sweep(P)
    % Run single gene simulation for TCD sweep
    % Similar to generate_gene_length_grid.m but optimized for speed
    
    global N PAS N_PAS Ef_ss;
    syms Ef real;
    
    % Set up geometry
    L_a = P.L_a;
    N = floor(P.geneLength_bp / L_a);
    PAS = floor(P.PASposition / L_a);
    N_PAS = N - PAS + 1;
    SD = floor(P.SD_bp / L_a);
    
    % Compute steady states
    [r_E_BeforePas] = compute_steady_states(P, P.EBindingNumber + 1);
    
    % Set up kPon values based on saturation distance
    kPon_vals = zeros(1, PAS);
    
    if PAS <= SD
        % Linear increase, never reaches kPon_max
        kPon_vals = linspace(P.kPon_min, P.kPon_min + (P.kPon_max - P.kPon_min) * (PAS / SD), PAS);
    else
        % PAS beyond SD
        if P.kPon_option == 1
            % Option 1: Saturate at SD
            kPon_vals(1:SD) = linspace(P.kPon_min, P.kPon_max, SD);
            kPon_vals((SD+1):PAS) = P.kPon_max;
        else
            % Option 2: Continue linear increase
            kPon_vals(1:SD) = linspace(P.kPon_min, P.kPon_max, SD);
            slope = (P.kPon_max - P.kPon_min) / SD;
            for i = (SD+1):PAS
                kPon_vals(i) = P.kPon_max + slope * (i - SD);
            end
        end
    end
    
    % Set up binding expressions
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
    P_sim = P;
end

