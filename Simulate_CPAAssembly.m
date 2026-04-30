% SimulateFigure8_CPAAssembly.m
% Simulate CPA complex assembly kinetics as a function of distance
% downstream of the PAS — analogous to Figure 8 from Chao et al. (1999)
%
% The CDF of termination events downstream of the PAS represents the
% fraction of polymerases that have committed to cleavage/polyadenylation
% within a given distance, analogous to "% Rescue of CAT Activity"
% measured in the cis-antisense rescue assay.

clear; clc;
fprintf('=== CPA Assembly vs Distance (Figure 8 Simulation) ===\n\n');

% --- BASE PARAMETERS ---
P.L_a = 100;
P.k_in    = 2;
P.k_e     = 65/P.L_a;
P.k_e2    = 30/P.L_a;
P.E_total = 100000;
P.L_total = 100000;
P.Pol_total = 70000;

% Binding/unbinding rates
P.kEon = 0.0000025;
P.kEoff = 0.1;
P.kHon = 2;
P.kHoff = 1;
P.kc = 0.1;

% Phosphorylation profile
P.kPon_min = 0.01;
P.kPon_slope = 0.005;
P.kPoff = 1;

% Gene geometry
P.geneLength_bp = 25000;
P.PASposition = 20000;

% --- CONFIGURATION ---
EBindingNumber = 5;

% --- RUN SIMULATION ---
fprintf('Running simulation (EBindingNumber = %d)...\n', EBindingNumber);

try
    [R_sol, REH_sol, P_sim] = run_termination_simulation(P, EBindingNumber);
    fprintf('Simulation completed successfully.\n');
catch ME
    error('Simulation failed: %s', ME.message);
end

% --- CALCULATE TERMINATION CDF ---
fprintf('Calculating termination profile (CDF)...\n');
[exit_cdf, distances_bp] = calculate_pas_cleavage_profile(R_sol, REH_sol, P_sim);

% Prepend zero for interpolation at distance = 0
distances_for_interp = [0; distances_bp(:)];
cdf_for_interp = [0; exit_cdf(:)];

% --- EVALUATE AT DESIRED SEPARATIONS ---
separations_bp = 0:10:1000;
rescue_fraction = interp1(distances_for_interp, cdf_for_interp, separations_bp, 'linear', 'extrap');

% --- PLOT ---
fprintf('Generating plot...\n');

figure('Position', [100, 100, 700, 500]);
plot(separations_bp, rescue_fraction * 100, '-', 'LineWidth', 2.5, 'Color', [0.2 0.4 0.8]);

xlabel('Poly(A)-AntiPoly(A) distance (bp)', 'FontSize', 13, 'FontWeight', 'bold');
ylabel('% CPA Assembly Completion', 'FontSize', 13, 'FontWeight', 'bold');
grid on;
set(gca, 'FontSize', 12);
box on;

% --- PRINT KEY VALUES ---
fprintf('\n=== Results at key distances ===\n');
key_distances = [100, 200, 300, 400, 500, 600, 800, 1000];
for d = key_distances
    val = interp1(distances_for_interp, cdf_for_interp, d, 'linear', 'extrap') * 100;
    fprintf('  %4d bp: %5.1f%% assembly\n', d, val);
end

fprintf('\n=== Simulation Complete ===\n');
