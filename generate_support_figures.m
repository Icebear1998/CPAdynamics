% Multi-Threaded Simulation for Support Figures
% This script runs simulations with different EBindingNumbers and generates:
% 1. Average E-binding profile comparison
% 2. Pol II state distribution heatmap (for N=5 or similar)
% 3. Fraction of "Competent" Pol II

clear;
close all;
global N PAS N_PAS Ef_ss;

% Define output directory
outputDir = 'SecondVersionResults/SupportFigures/';
if ~exist(outputDir, 'dir')
    mkdir(outputDir);
end

% --- Common Parameters (Match CPA_multipleE_main.m) ---
P.L_a = 100;
P.k_in = 2;
P.kEon = 0.00025;
P.kEoff = 1;
P.k_e = 65 / P.L_a;
P.k_e2 = 30 / P.L_a;
P.E_total = 100000;
P.L_total = 100000;
P.Pol_total = 70000;
P.kc = 0.1;

% Parameters for kHon calculation (from main: kHon approx 0.2, then renormalized)
% Note: The main script recalculates kHon based on E_bound at PAS. We will replicate this.
Initial_kHon = 0.1;
P.kHoff = 0.01;

kPon_min = 0.01;
kPon_slope = 0.005;
kPoff = 1;

geneLength_bp = 25000;
PASposition = 20000;
N = floor(geneLength_bp / P.L_a);
PAS = floor(PASposition / P.L_a);
N_PAS = N - PAS + 1;
kPon_vals = kPon_min + kPon_slope * (0 : N - 1);

% Analysis Scenarios
BindingNumbers = [9];
Colors = {'r', 'g', 'b', 'm', 'c'};
avg_E_profiles = {};

fprintf('Starting simulations for Support Figures...\n');

%% 1. Loop over Binding Numbers for Profile Comparison
figure('Name', 'Average E Binding Profile', 'Position', [100, 100, 800, 600]);
hold on;

for idx = 1 : length(BindingNumbers)
    nb = BindingNumbers(idx);
    fprintf('Simulating EBindingNumber = %d...\n', nb);

    P.EBindingNumber = nb;
    P.kHon = Initial_kHon;
    
    % Reset global Ef for each binding number
    Ef_ss = 0;

    % Step A: Pre-compute Steady States (Symbolic)
    [r_E_BeforePas] = compute_steady_states(P, nb + 1);

    % Prepare RE_vals (Symbolic array)
    RE_vals = sym(zeros(nb + 1, N));
    
    % Construct RE_vals matching main script logic
    for e = 1:nb+1
        % Substitute kPoff once (scalar)
        r_E_e_fixed = subs(r_E_BeforePas(e), 'kPoff', kPoff);
        
        % Substitute kPon vector element-wise
        for i = 1:N
            RE_vals(e, i) = subs(r_E_e_fixed, 'kPon', kPon_vals(i));
        end
    end
    
    % Create handle for total E bound (function of Ef)
    % Sums weighted average: sum( (1:nb) * Probability )
    P.RE_val_bind_E = matlabFunction(simplify(sum(sym(1:nb)' .* RE_vals(2:end, :), 1)), 'Vars', {'Ef'});
    
    % Step B: Solve ODE
    P.FirstRun = true;
    P.is_unphysical = false;
    X0 = 1e-6 * ones(N + N_PAS, 1);
    options = optimoptions('fsolve', 'Display', 'off', 'FunctionTolerance', 1e-8);
    X = fsolve(@(xx) ode_dynamics_multipleE(xx, P), X0, options);
    
    % Recalculate kHon constraint (match run_termination_simulation)
    avg_E_bound_temp = P.RE_val_bind_E(Ef_ss);
    P.kHon = P.kHon * avg_E_bound_temp(end);
    P.FirstRun = false;
    X = fsolve(@(xx) ode_dynamics_multipleE(xx, P), X, options);
    
    % Step C: Extract Results
    avg_E_bound = P.RE_val_bind_E(Ef_ss);
    
    % Store for plotting (PAS-centered coordinates)
    x_coords = ((1-PAS):(N-PAS)) * P.L_a / 1000;  % Position relative to PAS in kb
    plot(x_coords, avg_E_bound, 'Color', Colors{idx}, 'LineWidth', 2, ...
        'DisplayName', sprintf('Max Sites = %d', nb));

end

% Finalize Profile Plot
xlabel('Position relative to PAS (kb)');
ylabel('Average Number of Bound E Factors');
title('Average E-Factor Binding Profile');
legend('Location', 'northwest');
xline(0, 'k--', 'PAS');
grid on;

saveas(gcf, fullfile(outputDir, 'Average_E_Binding_Comparison.png'));

fprintf('Support figures generated in %s\n', outputDir);