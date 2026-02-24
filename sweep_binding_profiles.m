% Sweep Binding Profiles Script (Parallelized)
% This script sweeps over kEoff (or kPon_slope) and plots both:
% 1. Ser2P phosphorylation profile (dashed lines)
% 2. Average E-factor binding profile (solid lines)
% for a specific E binding number on the SAME plot

clear;
close all;

% --- CONFIGURATION ---
% Switch between parameter sweeps: 'kEoff' or 'kPon'
SweepParameter = 'kEoff';  % Change to 'kPon' to sweep kPon_slope instead

% Fixed E binding number for analysis
EBindingNumber = 5;

% Define sweep values based on parameter choice
if strcmp(SweepParameter, 'kEoff')
    % Sweep kEoff values (default is 10)
    SweepValues = [5, 10, 20];
    SweepLabel = 'k_{Eoff}';
else
    % Sweep kPon_slope values (default is 0.02)
    SweepValues = [0.01, 0.02, 0.04];
    SweepLabel = 'k_{Pon} slope';
end

% Define output directory
outputDir = 'SecondVersionResults/BindingProfileSweep/';
if ~exist(outputDir, 'dir')
    mkdir(outputDir);
end

% --- COMMON PARAMETERS ---
P.L_a = 100;
P.k_in = 2;
P.kEon = 0.00025;
P.kEoff = 10;  % Default value (will be swept if selected)
P.k_e = 65 / P.L_a;
P.k_e2 = 30 / P.L_a;
P.E_total = 100000;
P.L_total = 100000;
P.Pol_total = 70000;
P.kc = 0.1;
P.kHon = 0.1;  % Initial kHon value

Initial_kHon = 0.1;
P.kHoff = 0.01;

P.kPon_min = 0.01;
P.kPon_slope = 0.02;  % Default value (will be swept if selected)
P.kPoff = 1;

P.geneLength_bp = 25000;
P.PASposition = 20000;

% Colors for different sweep values
Colors = {'r', 'g', 'b', 'm', 'c'};

fprintf('Starting %s sweep for E Binding Number = %d ...\n', SweepParameter, EBindingNumber);

%% Pre-allocate storage for results
nSweeps = length(SweepValues);
avg_P_results = cell(nSweeps, 1);
avg_E_results = cell(nSweeps, 1);
legendLabels = cell(nSweeps, 1);

%% Loop over sweep values
for idx = 1 : nSweeps
    sweepVal = SweepValues(idx);
    
    % Create local parameter copy
    P_local = P;
    
    % Update the parameter being swept
    if strcmp(SweepParameter, 'kEoff')
        P_local.kEoff = sweepVal;
        fprintf('  Simulating kEoff = %.2f...\n', sweepVal);
        legendLabels{idx} = sprintf('%s = %.1f', SweepLabel, sweepVal);
    else
        P_local.kPon_slope = sweepVal;
        fprintf('  Simulating kPon_slope = %.3f...\n', sweepVal);
        legendLabels{idx} = sprintf('%s = %.3f', SweepLabel, sweepVal);
    end
    
    % Run simulation using existing helper function (avoids global variables)
    [R_sol, REH_sol, P_local, r_E_BeforePas] = run_termination_simulation(P_local, EBindingNumber);
    
    % Get compute steady states for Ser2P calculation
    [~, r_P] = compute_steady_states(P_local, EBindingNumber + 1);
    
    % Get necessary dimensions
    N = floor(P_local.geneLength_bp / P_local.L_a);
    PAS = floor(P_local.PASposition / P_local.L_a);
    kPon_vals = P_local.kPon_min + P_local.kPon_slope * (0 : N - 1);
    
    % Calculate Ser2P Profile
    P_vals = sym(zeros(EBindingNumber + 1, N));
    
    for e = 1:EBindingNumber+1
        r_P_e_fixed = subs(r_P(e), 'kPoff', P_local.kPoff);
        for i = 1:N
            P_vals(e, i) = subs(r_P_e_fixed, 'kPon', kPon_vals(i));
        end
    end
    
    % Get Ef_ss from global (set by run_termination_simulation)
    % Since we can't use global in parfor, we need to recalculate it
    global Ef_ss;
    Ef_ss_local = calculate_Efree_from_solution(R_sol, REH_sol, P_local);
    
    avg_P = zeros(1, N);
    for i = 1:N
        num_sum = 0; 
        den_sum = 0;
        for e = 1:EBindingNumber+1
            r_P_e = P_vals(e, i);
            val = double(R_sol(i) * double(subs(r_P_e, 'Ef', Ef_ss_local)));
            
            num_sum = num_sum + (e-1) * val;
            den_sum = den_sum + val;
        end
        if den_sum > 0
            avg_P(i) = num_sum / den_sum;
        end
    end
    
    % Get average E binding directly from the simulation result
    avg_E_bound = P_local.RE_val_bind_E(Ef_ss_local);
    
    % Store results
    avg_P_results{idx} = avg_P;
    avg_E_results{idx} = avg_E_bound;
end

%% Plot Combined Results
N_global = floor(P.geneLength_bp / P.L_a);
PAS_global = floor(P.PASposition / P.L_a);

fig = figure('Name', 'Ser2P and E-Factor Binding Profiles', 'Position', [100, 100, 1000, 700]);
hold on;

x_coords = ((1-PAS_global):(N_global-PAS_global)) * P.L_a / 1000;  % Position relative to PAS in kb

% Plot all results
for idx = 1 : nSweeps
    % Plot Ser2P (dashed line)
    plot(x_coords, avg_P_results{idx}, '--', 'Color', Colors{idx}, 'LineWidth', 2, ...
        'DisplayName', sprintf('Ser2P (%s)', legendLabels{idx}));
    
    % Plot Average E Binding (solid line)
    plot(x_coords, avg_E_results{idx}, '-', 'Color', Colors{idx}, 'LineWidth', 2, ...
        'DisplayName', sprintf('Avg E (%s)', legendLabels{idx}));
end

% Finalize plot
xlabel('Position relative to PAS (kb)', 'FontSize', 12);
ylabel('Average Number of Bound Factors', 'FontSize', 12);
title(sprintf('Ser2P (dashed) and E-Factor (solid) Binding Profiles\n%s Sweep, N=%d', SweepParameter, EBindingNumber), 'FontSize', 14);
legend('Location', 'northwest', 'FontSize', 10);
xline(0, 'k--', 'PAS', 'LineWidth', 1.5);
grid on;

saveas(fig, fullfile(outputDir, sprintf('Combined_Profile_%s_Sweep.png', SweepParameter)));

fprintf('Binding profile sweep completed. Results saved to %s\n', outputDir);

%% Helper function to calculate Ef_ss from solution
function Ef_ss = calculate_Efree_from_solution(R_sol, REH_sol, P)
    % Calculate free E based on solution and conservation
    N = floor(P.geneLength_bp / P.L_a);
    PAS = floor(P.PASposition / P.L_a);
    
    % Use the RE_val_bind_E function and conservation equation
    % E_total = E_free + E_bound
    % E_bound = sum(R.*avg_E_per_R) + sum(REH.*avg_E_per_REH)
    
    % Solve for Ef that satisfies conservation
    constraint = @(Ef_candidate) calculate_constraint(Ef_candidate, R_sol, REH_sol, P, N, PAS);
    
    options = optimoptions('fsolve', 'Display', 'off', 'FunctionTolerance', 1e-8);
    Ef_ss = fsolve(constraint, P.E_total * 0.5, options);
end

function residual = calculate_constraint(Ef_candidate, R, REH, P, N, PAS)
    % Calculate E conservation residual
    REvalbindE_candidate = P.RE_val_bind_E(Ef_candidate);
    
    E_used = sum(R(1:N)'.* REvalbindE_candidate) + ...
             sum(REH' .* REvalbindE_candidate(PAS:N));
    
    E_free_implied = P.E_total - E_used;
    
    residual = Ef_candidate - E_free_implied;
end
