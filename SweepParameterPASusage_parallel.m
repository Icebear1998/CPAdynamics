% --- CONFIGURATION: CHOOSE THE PARAMETER TO SWEEP ---
sweep_param_name = 'kHoff';

switch sweep_param_name
    case 'E_total'
        sweep_param_values = [10000, 30000, 70000, 100000, 120000, 150000];
    case 'kHoff'
        sweep_param_values = [0.001, 0.01, 0.1, 1, 5, 10];
    otherwise
        error('Selected sweep parameter is not defined.');
end

% --- BASE PARAMETERS ---
P.k_in = 2; P.k_e = 65/100; P.k_e2 = 30/100; P.E_total = 70000;
P.Pol_total = 70000; P.kEon = 0.00025; P.kEoff = 10; P.kHon = 0.2;
P.kHoff = 0.0125; P.kc = 0.05; P.kPon_min = 0.01; P.kPon_max = 1;
P.kPoff_const = 1; P.kPoff_max = 2; P.kPoff_min = 0.1; P.L_a = 100;
P.geneLength_bp = 25000; P.PASposition = 20000; EBindingNumber = 3;

% --- SIMULATION SETUP ---
inter_pas_distances_bp = 0:100:2500;
proximal_usage_results_cdf = zeros(length(inter_pas_distances_bp), length(sweep_param_values));

% --- Start a parallel pool of workers ---
if isempty(gcp('nocreate')); parpool; end
p = gcp; % Get the current parallel pool object
fprintf('Parallel pool started with %d workers...\n', p.NumWorkers);


% --- PARALLEL PARAMETER SWEEP LOOP ---
fprintf('Starting parallel sweep over parameter: %s\n', sweep_param_name);
parfor p_idx = 1:length(sweep_param_values)
    P_run = P; 
    P_run.(sweep_param_name) = sweep_param_values(p_idx);
    
    [R_sol, REH_sol, P_sim] = run_termination_simulation_parallel(P_run, EBindingNumber);

    % --- Flux-based CDF calculation ---
    flux_cleavage_per_node = P_sim.kc * REH_sol;
    total_cleavage_flux    = sum(flux_cleavage_per_node);
    
    if total_cleavage_flux > 1e-9
        cumulative_cleavage_flux = cumsum(flux_cleavage_per_node);
        termination_cdf = cumulative_cleavage_flux / total_cleavage_flux;
    else
        termination_cdf = zeros(size(REH_sol));
    end

    nodes_post_pas = 1:length(REH_sol);
    bp_post_pas = nodes_post_pas * P_sim.L_a;
    bp_for_interp = [0; bp_post_pas(:)];
    cdf_for_interp = [0; termination_cdf(:)];
    
    proximal_usage_results_cdf(:, p_idx) = interp1(bp_for_interp, cdf_for_interp, inter_pas_distances_bp, 'linear', 'extrap');
end
disp('All simulations complete.');

% --- PLOT THE RESULTS ---
figure('Position', [100, 100, 800, 600]);
hold on;
colors = lines(length(sweep_param_values));

for p_idx = 1:length(sweep_param_values)
    plot(inter_pas_distances_bp, proximal_usage_results_cdf(:, p_idx) * 100, ...
         'LineWidth', 2.5, 'Color', colors(p_idx,:), ...
         'DisplayName', sprintf('%s = %.3g', strrep(sweep_param_name, '_', '\_'), sweep_param_values(p_idx)));
end

plot_title = sprintf('CDF of Proximal PAS Usage vs. %s (EBindingNumber=%d)', strrep(sweep_param_name, '_', ' '), EBindingNumber);
xlabel('Inter-PAS Distance (bp)', 'FontSize', 12);
ylabel('Cumulative Proximal Site Usage (%)', 'FontSize', 12);
title(plot_title, 'FontSize', 14, 'FontWeight', 'bold');
grid on; legend('show', 'Location', 'best'); set(gca, 'FontSize', 10); box on;
ylim([0 100]);