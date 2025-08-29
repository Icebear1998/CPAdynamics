% This code is redundent, to be deleted
% --- Simulation & Plotting Setup ---
P.k_in = 2;
P.k_e = 65/100;
P.k_e2 = 30/100;
P.E_total = 70000;
P.Pol_total = 70000;
P.kEon = 0.00025;
P.kEoff = 10;
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

% --- SWEEP PARAMETER: kHon ---
% Define the different kHon values (PAS strengths) to test
kHon_values = [0.01, 0.03, 0.05, 0.1];

% Define the distances to calculate usage probability for
inter_pas_distances_bp = 0:100:2500;

% CHANGED: Matrix to store results for proximal usage
proximal_usage_results = zeros(length(inter_pas_distances_bp), length(kHon_values));

% --- LOOP THROUGH kHon VALUES ---
fprintf('Starting sweep over kHon values...\n');
for k_idx = 1:length(kHon_values)
    % Set the kHon for the current simulation run
    P.kHon = kHon_values(k_idx);
    fprintf('Running simulation for kHon = %.3f...\n', P.kHon);

    % --- Run Simulation to get the read-through ratio ---
    [ratio] = run_termination_simulation(P, EBindingNumber);

    % --- Calculate Proximal Usage Probability for this kHon value ---
    for i = 1:length(inter_pas_distances_bp)
        dist_bp = inter_pas_distances_bp(i);
        node_idx = round(dist_bp / P.L_a);
        if node_idx < 1; node_idx = 1; end
        if node_idx > length(ratio); node_idx = length(ratio); end
        
        % CHANGED: Proximal usage is 1 minus the read-through ratio (distal usage)
        proximal_usage_results(i, k_idx) = 1 - ratio(node_idx);
    end
end
disp('All simulations complete.');

% --- PLOT THE RESULTS ---
figure('Position', [100, 100, 800, 600]);
hold on;
colors = lines(length(kHon_values)); % Get a set of distinct colors

for k_idx = 1:length(kHon_values)
    % CHANGED: Plotting the proximal usage results
    plot(inter_pas_distances_bp, proximal_usage_results(:, k_idx) * 100, ...
         'LineWidth', 2, 'Color', colors(k_idx,:), ...
         'DisplayName', sprintf('kHon = %.3f', kHon_values(k_idx)));
end

% CHANGED: Customize the plot for proximal usage
xlabel('Inter-PAS Distance (bp)', 'FontSize', 12);
ylabel('Proximal Site Usage (%)', 'FontSize', 12);
title('Predicted Proximal PAS Usage vs. PAS Strength', 'FontSize', 14, 'FontWeight', 'bold');
grid on;
legend('show', 'Location', 'best');
set(gca, 'FontSize', 10);
box on;
%ylim([0 100]);


%% --- Helper function to run the simulation ---
function [ratio_vector] = run_termination_simulation(P, EBindingNumber)
    % This function now uses the stable, direct two-step solver.
    global N PAS N_PAS Ef_ss;
    syms Ef real;
    
    L_a = P.L_a;
    N = floor(P.geneLength_bp / L_a);
    PAS = floor(P.PASposition / L_a);
    N_PAS = N - PAS + 1;
    
    kHon_base = P.kHon; % Store the base kHon for this run
    
    % --- Symbolic Pre-computation ---
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

    % --- STABLE TWO-STEP SOLVER ---
    X_guess = 1e-6 * ones(N + N_PAS, 1);
    options = optimoptions('fsolve', 'Display', 'off', 'FunctionTolerance', 1e-8);

    % Step 1: Solve for the initial steady state and self-consistent Ef_ss
    P.FirstRun = true; 
    P.is_unphysical = false; 
    Ef_ss = 0;
    try
        X_base = fsolve(@(xx) ode_dynamics_multipleE(xx, P), X_guess, options);
    catch
        error('Solver failed in Step 1.');
    end
    if P.is_unphysical; error('Solver returned unphysical result in Step 1.'); end

    % Step 2: Update kHon and solve for the final steady state
    avg_E_bound = P.RE_val_bind_E(Ef_ss);
    P.FirstRun = false; 
    P.kHon = kHon_base * avg_E_bound(end);
    
    % Use the result of Step 1 as the initial guess for Step 2
    X_final = fsolve(@(xx) ode_dynamics_multipleE(xx, P), X_base, options);
    
    % --- Calculate and return the final ratio ---
    R_sol = X_final(1:N);
    REH_sol = X_final(N+1 : N+N_PAS);
    ratio_vector = (REH_sol(1:end) + R_sol(PAS:end)) / (R_sol(PAS-1) + 1e-9);
end