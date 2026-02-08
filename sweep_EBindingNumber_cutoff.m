% SWEEP_EBINDINGNUMBER_CUTOFF.m
% Plot position at which 25% cutoff vs. E Binding Number
% Simple script to analyze how E binding number affects termination distance

clear; clc;
fprintf('=== E Binding Number Sweep Analysis ===\n\n');

% --- BASE PARAMETERS ---
P.L_a = 100;
P.k_in = 2;
P.kEon = 0.00025;
P.kEoff = 10;
P.k_e = 65/100;
P.k_e2 = 30/100;
P.E_total = 100000;
P.Pol_total = 70000;
P.kHon = 0.2;
P.kHoff = 0.0125;
P.kc = 0.05;
P.kPon_min = 0.01;
P.kPon_slope = 0.1;
P.kPoff = 1;
P.geneLength_bp = 25000;
P.PASposition = 20000;


% --- SWEEP CONFIGURATION ---
EBindingNumber_values = [1, 2, 3, 4, 5];
cutoff_threshold = 0.5;  % 50% termination threshold

% Pre-allocate results
cutoff_positions = zeros(size(EBindingNumber_values));

% --- SWEEP LOOP ---
fprintf('Sweeping E Binding Number: %s\n', mat2str(EBindingNumber_values));

for i = 1:length(EBindingNumber_values)
    EBindingNumber = EBindingNumber_values(i);
    fprintf('Running EBindingNumber = %d... ', EBindingNumber);
    
    try
        % Run simulation
        [R_sol, REH_sol, P_sim] = run_termination_simulation(P, EBindingNumber);
        
        % Calculate termination profile (CDF)
        [exit_cdf, distances_bp] = calculate_pas_usage_profile(R_sol, REH_sol, P_sim);
        
        % Find position where cutoff threshold is reached
        if max(exit_cdf) >= cutoff_threshold
            cutoff_positions(i) = interp1(exit_cdf, distances_bp, cutoff_threshold, 'linear', 'extrap');
        else
            cutoff_positions(i) = NaN;  % Threshold not reached
        end
        
        fprintf('Cutoff position = %.0f bp\n', cutoff_positions(i));
        
    catch ME
        fprintf('ERROR: %s\n', ME.message);
        cutoff_positions(i) = NaN;
    end
end

% --- PLOT RESULTS ---
figure('Position', [100, 100, 800, 600]);
plot(EBindingNumber_values, cutoff_positions, 'o', 'LineWidth', 2.5, 'MarkerSize', 10, ...
     'Color', [0, 0.4470, 0.7410], 'MarkerFaceColor', [0, 0.4470, 0.7410]);

% Add data labels
for i = 1:length(EBindingNumber_values)
    if ~isnan(cutoff_positions(i))
        text(EBindingNumber_values(i), cutoff_positions(i) + 50, ...
             sprintf('%.0f', cutoff_positions(i)), ...
             'FontSize', 11, 'HorizontalAlignment', 'center');
    end
end

xlabel('Number of maxium CPA factor binding', 'FontSize', 12, 'FontWeight', 'bold');
ylabel('CAD (bp)', 'FontSize', 12, 'FontWeight', 'bold');
xlim([0.5 7.5]);
% Set x-axis ticks to only integer values
ax = gca;
ax.XTick = unique(round(ax.XTick));
% Format x-axis tick labels to display as integers
xtickformat('%d');
set(gca, 'FontSize', 11);
grid on;
box on;

fprintf('\n=== Analysis Complete ===\n');
fprintf('Results:\n');
for i = 1:length(EBindingNumber_values)
    fprintf('  EBindingNumber = %d: %.0f bp\n', EBindingNumber_values(i), cutoff_positions(i));
end

%% --- HELPER FUNCTION ---
function [R_sol, REH_sol, P] = run_termination_simulation(P, EBindingNumber)
    % Run single simulation for given E binding number
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
    
    % Setup kPon values with linear increase
    kPon_vals = P.kPon_min + P.kPon_slope * (0:N-1);
    RE_vals = sym(zeros(EBindingNumber + 1, N));
    
    for e = 1:(EBindingNumber + 1)
        for idx = 1:N
            RE_vals(e, idx) = subs(r_E_BeforePas(e), {'kPon', 'kPoff'}, {kPon_vals(idx), P.kPoff});
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

