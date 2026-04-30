% SCRIPT to plot Proximal PAS Usage vs. Inter-PAS Distance
% GENERALIZED to sweep over any chosen global factor
% Uses flux-based termination profile to calculate proximal site usage
save_results = true;

% MODIFIED: Start a parallel pool of workers if one is not already running.
if isempty(gcp('nocreate')); parpool; end

% --- CONFIGURATION: CHOOSE THE PARAMETER TO SWEEP ---
sweep_param_name = 'kHoff';

switch sweep_param_name
    case 'E_total'
        sweep_param_values = [10000, 30000, 70000, 100000, 120000, 150000];
    case 'k_e'
        sweep_param_values = [45/100, 65/100, 85/100];
    case 'kc'
        sweep_param_values = [0.01, 0.05, 0.1, 0.2];
    case 'kHoff'
        sweep_param_values = [0.1, 1, 5, 10, 50, 100];
    otherwise
        error('Selected sweep parameter is not defined in the switch-case block.');
end

% --- BASE PARAMETERS (matched to CPA_multipleE_main.m) ---
P.L_a = 100;
P.k_in    = 2;
P.k_e     = 65/P.L_a;
P.k_e2    = 30/P.L_a;
P.E_total = 100000;
P.L_total = 100000;
P.Pol_total = 70000;

P.kEon = 0.0000025;
P.kEoff = 0.1;
P.kHon = 2;
P.kHoff = 1;
P.kc = 0.1;

P.kPon_min = 0.01;
P.kPon_slope = 0.005;
P.kPoff = 1;

P.geneLength_bp = 50000;
P.PASposition = 20000;
EBindingNumber = 5;

% --- SIMULATION SETUP ---
inter_pas_distances_bp = 0:1000:30000;
proximal_usage_results = zeros(length(inter_pas_distances_bp), length(sweep_param_values));

% --- PARAMETER SWEEP LOOP ---
fprintf('Starting parallel sweep over parameter: %s\n', sweep_param_name);

% MODIFIED: Changed 'for' to 'parfor' to distribute iterations across workers.
parfor p_idx = 1:length(sweep_param_values)
    % Create a copy of the parameters for this iteration.
    P_run = P;
    P_run.(sweep_param_name) = sweep_param_values(p_idx);

    % Run simulation using the shared function (uses numerical null-space,
    % robust for any EBindingNumber including >= 5)
    [R_sol, REH_sol, P_sim] = run_termination_simulation(P_run, EBindingNumber);

    % Calculate termination profile (CDF)
    % The CDF at distance X represents the fraction of polymerases that 
    % terminated by distance X, which IS the proximal site usage
    [exit_cdf, distances_bp] = calculate_pas_cleavage_profile(R_sol, REH_sol, P_sim);
    
    % Interpolate to get proximal usage at specific inter-PAS distances
    proximal_usage_results(:, p_idx) = interp1(distances_bp, exit_cdf, inter_pas_distances_bp, 'linear', 'extrap');
end
disp('All parallel simulations complete.');

% --- PLOT THE RESULTS ---
figure('Position', [100, 100, 800, 600]);
hold on;
colors = lines(length(sweep_param_values));

for p_idx = 1:length(sweep_param_values)
    semilogx(inter_pas_distances_bp, proximal_usage_results(:, p_idx) * 100, ...
         'LineWidth', 2.5, 'Color', colors(p_idx,:), ...
         'DisplayName', sprintf('%s = %.3g', strrep(sweep_param_name, '_', '\_'), sweep_param_values(p_idx)));
end

xlabel('Inter-PAS distance (bp) (log scale)', 'FontSize', 12, 'FontWeight', 'bold');
ylabel('Proximal site usage (%)', 'FontSize', 12, 'FontWeight', 'bold');
%title(sprintf('Proximal PAS Usage vs %s (EBindingNumber=%d)', strrep(sweep_param_name, '_', ' '), EBindingNumber), ...
      %'FontSize', 14, 'FontWeight', 'bold');
grid on; 
legend('show', 'Location', 'best', 'FontSize', 11); 
set(gca, 'FontSize', 11); 
box on;
ylim([0 100]);

if save_results
    % --- SAVE RESULTS ---
    % Prepare data structure for saving
    data.results_matrix = proximal_usage_results;
    data.x_values = inter_pas_distances_bp;
    data.sweep_values = sweep_param_values;
    data.sweep_param = sweep_param_name;
    data.description = 'Proximal PAS usage (%) at different inter-PAS distances';

    % Save results using the utility function
    extra_info = sprintf('EBinding%d', EBindingNumber);
    save_analysis_results('ProximalPASUsage_ParameterSweep', data, P, 'ExtraInfo', extra_info);
end
