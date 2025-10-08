% SCRIPT to plot Proximal PAS Usage vs. Inter-PAS Distance
% GENERALIZED to sweep over any chosen global factor
% UPDATED with user-proposed flux-based CDF to ensure 100% completion
save_results = true;

% MODIFIED: Start a parallel pool of workers if one is not already running.
if isempty(gcp('nocreate')); parpool; end

% --- CONFIGURATION: CHOOSE THE PARAMETER TO SWEEP ---
sweep_param_name = 'E_total';

switch sweep_param_name
    case 'E_total'
        sweep_param_values = [10000, 30000, 70000, 100000, 120000, 150000];
    case 'k_e'
        sweep_param_values = [45/100, 65/100, 85/100];
    case 'kc'
        sweep_param_values = [0.01, 0.05, 0.1, 0.2];
    case 'kHoff'
        sweep_param_values = [0.001, 0.01, 0.1, 1, 5, 10];
    otherwise
        error('Selected sweep parameter is not defined in the switch-case block.');
end

% --- BASE PARAMETERS ---
% (Parameter definitions are unchanged)
P.k_in = 2; P.k_e = 65/100; P.k_e2 = 30/100; P.E_total = 70000;
P.Pol_total = 70000; P.kEon = 0.00025; P.kEoff = 10; P.kHon = 0.2;
P.kHoff = 0.0125; P.kc = 0.05; P.kPon_min = 0.01; P.kPon_max = 1;
P.kPoff_const = 1; P.kPoff_max = 2; P.kPoff_min = 0.1; P.L_a = 100;
P.geneLength_bp = 25000; P.PASposition = 20000; EBindingNumber = 5;

% --- SIMULATION SETUP ---
inter_pas_distances_bp = 0:100:2500;
proximal_usage_results_cdf = zeros(length(inter_pas_distances_bp), length(sweep_param_values));

% --- PARAMETER SWEEP LOOP ---
fprintf('Starting parallel sweep over parameter: %s\n', sweep_param_name);

% MODIFIED: Changed 'for' to 'parfor' to distribute iterations across workers.
parfor p_idx = 1:length(sweep_param_values)
    % Create a copy of the parameters for this iteration.
    P_run = P;
    P_run.(sweep_param_name) = sweep_param_values(p_idx);
    
    % NOTE: Removed the fprintf from inside the loop as output order is not guaranteed.

    % Run simulation to get steady-state concentrations and parameters used
    [R_sol, REH_sol, P_sim] = run_termination_simulation(P_run, EBindingNumber);

    [exit_cdf, distances_bp] = calculate_pas_usage_profile(R_sol, REH_sol, P_sim);
    
    % Interpolate the CDF and store it in the results matrix.
    % parfor handles this "sliced" assignment correctly.
    
    proximal_usage_results_cdf(:, p_idx) = interp1(distances_bp, exit_cdf, inter_pas_distances_bp, 'linear', 'extrap');
end
disp('All parallel simulations complete.');

% --- PLOT THE RESULTS ---
figure('Position', [100, 100, 800, 600]);
hold on;
colors = lines(length(sweep_param_values));

for p_idx = 1:length(sweep_param_values)
    plot(inter_pas_distances_bp, proximal_usage_results_cdf(:, p_idx) * 100, ...
         'LineWidth', 2.5, 'Color', colors(p_idx,:), ...
         'DisplayName', sprintf('%s = %.3g', strrep(sweep_param_name, '_', '\_'), sweep_param_values(p_idx)));
end

plot_title = sprintf('CDF of Total Polymerase Exit vs. %s (EBindingNumber=%d)', strrep(sweep_param_name, '_', ' '), EBindingNumber);
xlabel('Distance Downstream of PAS (bp)', 'FontSize', 12);
ylabel('Cumulative Polymerase Exit (%)', 'FontSize', 12);
title(plot_title, 'FontSize', 14, 'FontWeight', 'bold');
grid on; legend('show', 'Location', 'best'); set(gca, 'FontSize', 10); box on;
ylim([0 100]);

if save_results
    % --- SAVE RESULTS ---
    % Prepare data structure for saving
    data.results_matrix = proximal_usage_results_cdf;
    data.x_values = inter_pas_distances_bp;
    data.sweep_values = sweep_param_values;
    data.sweep_param = sweep_param_name;

    % Save results using the utility function
    extra_info = sprintf('EBinding%d', EBindingNumber);
    save_analysis_results('ParameterSweep', data, P, 'ExtraInfo', extra_info);
end

%% --- Helper function to run the simulation ---
function [R_sol, REH_sol, P] = run_termination_simulation(P, EBindingNumber)
    % MODIFIED: Now returns the P struct as well, as it contains k_e, k_e2 etc.
    global N PAS N_PAS Ef_ss;
    syms Ef real;
    
    L_a = P.L_a; N = floor(P.geneLength_bp / L_a); PAS = floor(P.PASposition / L_a); N_PAS = N - PAS + 1;
    kHon_base = P.kHon;
    
    [r_E_BeforePas] = compute_steady_states(P, EBindingNumber + 1);
    kPon_vals = linspace(P.kPon_min, P.kPon_max, PAS);
    RE_vals = sym(zeros(EBindingNumber + 1, N));
    for e = 1:EBindingNumber + 1
        for idx = 1:length(kPon_vals); RE_vals(e, idx) = subs(r_E_BeforePas(e), {'kPon', 'kPoff'}, {kPon_vals(idx), P.kPoff_const}); end
        for idx = PAS+1:N; RE_vals(e, idx) = subs(r_E_BeforePas(e), {'kPon', 'kPoff'}, {P.kPon_max, P.kPoff_const}); end
    end
    P.RE_val_bind_E = matlabFunction(simplify(sum(sym(1:EBindingNumber)' .* RE_vals(2:end, :), 1)), 'Vars', {Ef});

    X_guess = 1e-6 * ones(N + N_PAS, 1);
    options = optimoptions('fsolve', 'Display', 'off', 'FunctionTolerance', 1e-8);

    P.FirstRun = true; P.is_unphysical = false; Ef_ss = 0;
    try; X_base = fsolve(@(xx) ode_dynamics_multipleE(xx, P), X_guess, options);
    catch; error('Solver failed in Step 1.'); end
    if P.is_unphysical; error('Solver returned unphysical result in Step 1.'); end

    avg_E_bound = P.RE_val_bind_E(Ef_ss);
    P.FirstRun = false; P.kHon = kHon_base * avg_E_bound(end);
    X_final = fsolve(@(xx) ode_dynamics_multipleE(xx, P), X_base, options);
    
    R_sol = X_final(1:N);
    REH_sol = X_final(N+1 : N+N_PAS);
end