% SCRIPT to demonstrate how a global factor (E_total) shifts APA usage

clear all;
close all;
clc;

% --- Simulation & Plotting Setup ---
P.k_in = 2;
P.k_e = 65/100;
P.k_e2 = 30/100;
P.Pol_total = 70000;
P.kEon = 0.00025;
P.kEoff = 10;
P.kHon = 0.05;
P.kHoff = 0.0025;
P.kc = 0.8;
P.kPon_min = 0.01;
P.kPon_max = 1;
P.kPoff_const = 1;
P.kPoff_max = 2;
P.kPoff_min = 0.1;
P.L_a = 100;
P.geneLength_bp = 25000;
P.PASposition = 20000;
EBindingNumber = 3;

% --- SWEEP PARAMETER: E_total ---
% Define the different global factor concentrations to test
E_total_values = [10000, 30000, 70000, 100000, 120000];

% Define the distances to calculate usage probability for
inter_pas_distances_bp = 0:100:2500;

% Matrix to store results: each column is an E_total value
% We will store the DISTAL usage probability (the read-through ratio)
distal_usage_results = zeros(length(inter_pas_distances_bp), length(E_total_values));

% --- LOOP THROUGH E_total VALUES ---
fprintf('Starting sweep over E_total values...\n');
for e_idx = 1:length(E_total_values)
    P_run = P; % Use a temporary struct for this run
    P_run.E_total = E_total_values(e_idx);
    fprintf('Running simulation for E_total = %d...\n', P_run.E_total);

    % Run the simulation to get the read-through ratio vector
    [ratio] = run_termination_simulation(P_run, EBindingNumber);

    % Interpolate the ratio at the desired distances
    nodes_post_pas = 0:(length(ratio)-1);
    bp_post_pas = nodes_post_pas * P.L_a;
    distal_usage_results(:, e_idx) = interp1(bp_post_pas, ratio, inter_pas_distances_bp, 'linear', 'extrap');
end
disp('All simulations complete.');

% --- PLOT FIGURE A: Proximal Site Usage vs. Inter-PAS Distance ---
figure('Position', [100, 100, 800, 600]);
hold on;
colors = cool(length(E_total_values)); % Get a set of distinct colors

% Proximal usage is 1 - distal usage
proximal_usage_results = 1 - distal_usage_results;

for e_idx = 1:length(E_total_values)
    plot(inter_pas_distances_bp, proximal_usage_results(:, e_idx) * 100, ...
         'LineWidth', 2.5, 'Color', colors(e_idx,:), ...
         'DisplayName', sprintf('E_{total} = %dK', E_total_values(e_idx)/1000));
end
xlabel('Inter-PAS Distance (bp)', 'FontSize', 12);
ylabel('Proximal Site Usage (%)', 'FontSize', 12);
title('Global Factor Abundance Tunes Kinetic Competition', 'FontSize', 14, 'FontWeight', 'bold');
grid on; legend('show', 'Location', 'best'); set(gca, 'FontSize', 10); box on; ylim([0 100]);

% --- PLOT FIGURE B: Predicted Distribution of Chosen PAS Distances ---
figure('Position', [950, 100, 800, 600]);
hold on;

% 1. Model the "Available" distribution of tandem sites (f_avail)
% We use a log-normal distribution, which is good for skewed data like this.
% Parameters chosen to create a peak around 300 bp.
mu = log(280); % Corresponds to a peak near 300
sigma = 0.6;
f_avail = lognpdf(inter_pas_distances_bp, mu, sigma);
f_avail = f_avail / trapz(inter_pas_distances_bp, f_avail); % Normalize

% Plot the available distribution as a shaded area
patch([inter_pas_distances_bp, fliplr(inter_pas_distances_bp)], [f_avail, zeros(size(f_avail))], ...
      'k', 'FaceAlpha', 0.15, 'EdgeColor', 'none', 'DisplayName', 'Available Sites');

% 2. Calculate and plot the "Chosen" distribution for each E_total
for e_idx = 1:length(E_total_values)
    % The probability of choosing a distal site at this distance
    P_distal_choice = distal_usage_results(:, e_idx);
    
    % The unnormalized PDF of chosen sites
    unnormalized_pdf = f_avail' .* P_distal_choice;
    
    % Normalize to get the final predicted distribution
    f_chosen = unnormalized_pdf / trapz(inter_pas_distances_bp, unnormalized_pdf);
    
    plot(inter_pas_distances_bp, f_chosen, 'LineWidth', 2.5, 'Color', colors(e_idx,:), ...
         'DisplayName', sprintf('E_{total} = %dK', E_total_values(e_idx)/1000));
end
xlabel('Inter-PAS Distance (bp)', 'FontSize', 12);
ylabel('Probability Density', 'FontSize', 12);
title('Global Factor Levels Reshape the 3'' UTR Landscape', 'FontSize', 14, 'FontWeight', 'bold');
grid on; legend('show', 'Location', 'best'); set(gca, 'FontSize', 10); box on;


%% --- Helper function to run the simulation ---
function [ratio_vector] = run_termination_simulation(P, EBindingNumber)
    % This function uses the full symbolic pre-computation as requested
    global N PAS N_PAS Ef_ss;
    syms Ef real;
    
    L_a = P.L_a;
    N = floor(P.geneLength_bp / L_a);
    PAS = floor(P.PASposition / L_a);
    N_PAS = N - PAS + 1;
    
    kHon_base = P.kHon;
    [r_E_BeforePas] = compute_steady_states(P, EBindingNumber + 1);
    kPon_vals = linspace(P.kPon_min, P.kPon_max, PAS);
    RE_vals = sym(zeros(EBindingNumber + 1, N));
    for e = 1:EBindingNumber + 1
        for idx = 1:length(kPon_vals)
            kPon_val = kPon_vals(idx);
            RE_vals(e, idx) = subs(r_E_BeforePas(e), {'kPon', 'kPoff'}, {kPon_val, P.kPoff_const});
        end
        for idx = PAS+1:N
            RE_vals(e, idx) = subs(r_E_BeforePas(e), {'kPon', 'kPoff'}, {P.kPon_max, P.kPoff_min});
        end
    end
    P.RE_val_bind_E = matlabFunction(simplify(sum(sym(1:EBindingNumber)' .* RE_vals(2:end, :), 1)), 'Vars', {Ef});

    max_iterations = 10;
    convergence_tol = 1e-6;
    X_guess = 1e-6 * ones(N + N_PAS, 1);
    options = optimoptions('fsolve', 'Display', 'off', 'FunctionTolerance', 1e-8);

    for iter = 1:max_iterations
        P.FirstRun = true; P.is_unphysical = false; Ef_ss = 0;
        X = fsolve(@(xx) ode_dynamics_multipleE(xx, P), X_guess, options);
        if P.is_unphysical; error('Solver failed.'); end
        avg_E_bound = P.RE_val_bind_E(Ef_ss);
        P.FirstRun = false; P.kHon = kHon_base * avg_E_bound(end);
        X_adj = fsolve(@(xx) ode_dynamics_multipleE(xx, P), X, options);
        solution_change = norm(X_adj - X_guess) / norm(X_adj);
        if solution_change < convergence_tol; X_guess = X_adj; break; end
        X_guess = X_adj;
    end
    
    R_sol = X_guess(1:N);
    REH_sol = X_guess(N+1 : N+N_PAS);
    ratio_vector = (REH_sol(1:end) + R_sol(PAS:end)) / (R_sol(PAS-1) + 1e-9);
end