% PASUAGEVSINTERPAS DISTANCE.m
% Calculate and plot Proximal PAS Usage vs. Inter-PAS Distance
% Uses flux-based CDF approach for accurate termination probability

clear; clc;
fprintf('=== Proximal PAS Usage vs Inter-PAS Distance ===\n\n');

% --- BASE PARAMETERS ---
P.L_a = 100;
P.k_in = 2;
P.kEon = 0.00025;
P.kEoff = 10;
P.k_e = 65/P.L_a;
P.k_e2 = 30/P.L_a;
P.E_total = 70000;
P.Pol_total = 70000;
P.kHon = 0.2;
P.kHoff = 0.125;
P.kc = 0.05;
P.kPon_min = 0.01;
P.kPon_max = 4;
P.kPoff_const = 1;
P.geneLength_bp = 25000;
P.PASposition = 20000;

% --- CONFIGURATION ---
EBindingNumber = 5;
save_result = false;

% --- RUN SIMULATION ---
fprintf('Running simulation for EBindingNumber = %d...\n', EBindingNumber);

try
    [R_sol, REH_sol, P_sim] = run_termination_simulation(P, EBindingNumber);
    fprintf('Simulation completed successfully.\n');
catch ME
    error('Simulation failed: %s', ME.message);
end

% --- CALCULATE TERMINATION PROFILE (CDF) ---
fprintf('Calculating termination profile...\n');

% Use standardized function to calculate PAS usage profile
[exit_cdf, distances_bp] = calculate_pas_usage_profile(R_sol, REH_sol, P_sim);

% Prepend zero for interpolation at distance = 0
distances_for_interp = [0; distances_bp(:)];
cdf_for_interp = [0; exit_cdf(:)];

% --- CALCULATE PROXIMAL PAS USAGE ---
fprintf('Calculating proximal PAS usage...\n');

% Define inter-PAS distances
inter_pas_distances_bp = 100:50:500;

% Calculate proximal PAS usage for each inter-PAS distance
proximal_usage_prob = zeros(size(inter_pas_distances_bp));

for i = 1:length(inter_pas_distances_bp)
    dist_bp = inter_pas_distances_bp(i);
    
    % CDF at this distance tells us fraction that terminated by this point (at proximal PAS)
    % So proximal usage = CDF at that distance
    proximal_usage_prob(i) = interp1(distances_for_interp, cdf_for_interp, dist_bp, 'linear', 'extrap');
end

% --- PLOT RESULTS ---
fprintf('Generating plot...\n');

figure('Position', [100, 100, 800, 600]);
hold on;

plot(inter_pas_distances_bp, proximal_usage_prob * 100, 'o-', 'LineWidth', 2.5, ...
     'MarkerSize', 6, 'Color', [0, 0.4470, 0.7410], 'DisplayName', 'Proximal PAS Usage');

% Add vertical line at 300 bp (median tandem distance)
line([300 300], [0 100], 'Color', 'r', 'LineStyle', '--', 'LineWidth', 1.5, ...
     'DisplayName', 'Typical inter-PAS distance (300 bp)');

% Customize plot
xlabel('Inter-PAS Distance (bp)', 'FontSize', 12, 'FontWeight', 'bold');
ylabel('Proximal Site Usage (%)', 'FontSize', 12, 'FontWeight', 'bold');
title(sprintf('Proximal PAS Usage vs Inter-PAS Distance (EBindingNumber=%d)', EBindingNumber), ...
      'FontSize', 14, 'FontWeight', 'bold');
grid on;
legend('show', 'Location', 'best', 'FontSize', 11);
set(gca, 'FontSize', 11);
box on;
ylim([0 100]);

hold off;

% --- SAVE RESULTS ---
if save_result
    fprintf('Saving results...\n');
    
    % Prepare data structure
    data.EBindingNumber = EBindingNumber;
    data.inter_pas_distances_bp = inter_pas_distances_bp;
    data.proximal_usage_prob = proximal_usage_prob;
    data.exit_cdf = exit_cdf;
    data.distances_bp = distances_bp;
    data.R_sol = R_sol;
    data.REH_sol = REH_sol;
    
    % Save using utility function
    save_analysis_results('PASUsagevsInterPASDistance', data, P);
end

fprintf('\n=== Analysis Complete ===\n');
fprintf('Results at key distances:\n');
fprintf('  100 bp: %.1f%% proximal usage\n', interp1(inter_pas_distances_bp, proximal_usage_prob*100, 100));
fprintf('  300 bp: %.1f%% proximal usage\n', interp1(inter_pas_distances_bp, proximal_usage_prob*100, 300));
fprintf('  500 bp: %.1f%% proximal usage\n', interp1(inter_pas_distances_bp, proximal_usage_prob*100, 500));
fprintf('  1000 bp: %.1f%% proximal usage\n', interp1(inter_pas_distances_bp, proximal_usage_prob*100, 1000));

%% --- HELPER FUNCTION ---
function [R_sol, REH_sol, P] = run_termination_simulation(P, EBindingNumber)
    % Run termination simulation for given E binding number
    global N PAS N_PAS Ef_ss;
    syms Ef real;
    
    % Setup geometry
    L_a = P.L_a;
    N = floor(P.geneLength_bp / L_a);
    PAS = floor(P.PASposition / L_a);
    N_PAS = N - PAS + 1;
    kHon_base = P.kHon;
    
    % Compute steady states
    [r_E_BeforePas] = compute_steady_states(P, EBindingNumber + 1);
    
    % Setup kPon values
    kPon_vals = linspace(P.kPon_min, P.kPon_max, PAS);
    RE_vals = sym(zeros(EBindingNumber + 1, N));
    
    for e = 1:(EBindingNumber + 1)
        for idx = 1:length(kPon_vals)
            RE_vals(e, idx) = subs(r_E_BeforePas(e), {'kPon', 'kPoff'}, {kPon_vals(idx), P.kPoff_const});
        end
        for idx = (PAS+1):N
            RE_vals(e, idx) = subs(r_E_BeforePas(e), {'kPon', 'kPoff'}, {P.kPon_max, P.kPoff_const});
        end
    end
    
    % Create E binding function
    P.RE_val_bind_E = matlabFunction(simplify(sum(sym(1:EBindingNumber)' .* RE_vals(2:end, :), 1)), 'Vars', {Ef});
    
    % Solve system - Step 1
    X_guess = 1e-6 * ones(N + N_PAS, 1);
    options = optimoptions('fsolve', 'Display', 'off', 'FunctionTolerance', 1e-8);
    
    P.FirstRun = true;
    P.is_unphysical = false;
    Ef_ss = 0;
    
    X_base = fsolve(@(xx) ode_dynamics_multipleE(xx, P), X_guess, options);
    
    if P.is_unphysical
        error('Unphysical result in Step 1');
    end
    
    % Solve system - Step 2 with adjusted kHon
    avg_E_bound = P.RE_val_bind_E(Ef_ss);
    P.FirstRun = false;
    P.kHon = kHon_base * avg_E_bound(end);
    
    X_final = fsolve(@(xx) ode_dynamics_multipleE(xx, P), X_base, options);
    
    R_sol = X_final(1:N);
    REH_sol = X_final(N+1:N+N_PAS);
end