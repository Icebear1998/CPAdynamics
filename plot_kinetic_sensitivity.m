% PLOT_KINETIC_SENSITIVITY.m
% Analyze CAD vs. Gene Length for different E-factor binding kinetics
%
% This script demonstrates that the "runway effect" is a kinetic bottleneck:
% - Faster k_on: Shorter runway needed (curve shifts left)
% - Slower k_on: Longer runway needed (curve shifts right)
% - All curves plateau at long gene lengths
%
% Expected output: Three overlaid CAD curves showing kinetic sensitivity
%
% WORKFLOW: This script uses self-consistent (R_free, E_free) from the
% gene length analysis workflow. If analysis results are not available,
% it falls back to fixed values (30% of total).

clear;
close all;
global N PAS N_PAS Ef_ss;

fprintf('=== Kinetic Sensitivity Analysis ===\n');
fprintf('Calculating CAD vs. Gene Length for different k_on values...\n\n');

% --- SETUP OUTPUT DIRECTORY ---
outputDir = 'SecondVersionResults/GeneLengthAnalysis/';
if ~exist(outputDir, 'dir')
    mkdir(outputDir);
end

% --- LOAD SELF-CONSISTENT SOLUTION (if available) ---
fprintf('Attempting to load self-consistent (R_free, E_free) solution...\n');

analysis_dir = 'SecondVersionResults/GeneLengthAnalysis/';
use_self_consistent = false;
R_free_value = 21000; % Default: 30% of 70,000
E_free_value = 21000; % Default: 30% of 70,000

if exist(analysis_dir, 'dir')
    analysis_files = dir(fullfile(analysis_dir, 'gene_length_TCD_analysis_*.mat'));
    if ~isempty(analysis_files)
        [~, newest_idx] = max([analysis_files.datenum]);
        analysis_filename = fullfile(analysis_dir, analysis_files(newest_idx).name);

        try
            load(analysis_filename, 'analysis_results');
            R_free_value = analysis_results.solution.R_free;
            E_free_value = analysis_results.solution.E_free;
            use_self_consistent = true;
            fprintf('  Loaded self-consistent solution from: %s\n', ...
                analysis_files(newest_idx).name);
            fprintf('  R_free = %.0f, E_free = %.0f\n', R_free_value, E_free_value);
        catch
            fprintf('  Could not load analysis results, using default values\n');
        end
    else
        fprintf('  No analysis results found, using default values\n');
    end
else
    fprintf('  Analysis directory not found, using default values\n');
end

if ~use_self_consistent
    fprintf('  Using default: R_free = %.0f, E_free = %.0f (30%% of total)\n', ...
        R_free_value, E_free_value);
end

% --- DEFINE PARAMETERS ---
% Base parameters
P.L_a = 100;          % Node spacing (bp)
P.kEoff = 10;         % E-factor unbinding rate
P.kHon = 0.2;         % H-factor binding rate (initial)
P.kHoff = 0.0125;     % H-factor unbinding rate
P.kc = 0.05;          % Cleavage rate
P.k_e = 65 / P.L_a;   % Pol II elongation rate
P.k_e2 = 30 / P.L_a;  % REH elongation rate
P.k_in = 2;           % Pol II initiation rate
P.EBindingNumber = 2; % Max E-binding sites
P.kPon_min = 0.01;    % Min Ser2P phosphorylation (at TSS)
P.kPon_slope = 0.02;  % Slope of linear increase
P.kPoff_const = 1;    % Ser2P dephosphorylation

% Free resource pools (from self-consistent solution or default)
P.Pol_free = R_free_value;
P.E_free = E_free_value;

% Total resource pools (needed by global ode_dynamics_multipleE.m)
P.Pol_total = 70000; % Total polymerase
P.E_total = 70000;   % Total E-factors

% Flags for global ode_dynamics_multipleE.m compatibility
P.FirstRun = false;     % Use simple E_free calculation (not self-consistent solver)
P.is_unphysical = false;

% Fixed after-PAS region
after_PAS_length = 5000; % 5 kb after PAS

% --- DEFINE KINETIC SCENARIOS ---
% Three k_on values: Slow, Baseline, Fast
kEon_scenarios = [ 0.000125, 0.00025, 0.0005 ]; % 0.5x, 1x, 2x baseline
scenario_names = {'Slow (0.5\times k_{on})', 'Baseline (1\times k_{on})', ...
                  'Fast (2\times k_{on})'};
scenario_colors = {'r', 'b', 'g'};
n_scenarios = length(kEon_scenarios);

fprintf('\nKinetic scenarios:\n');
for i = 1 : n_scenarios
    fprintf('  %s: k_on = %.6f\n', scenario_names{i}, kEon_scenarios(i));
end

% --- GENE LENGTH SWEEP ---
tss_to_pas_distances = logspace(log10(2000), log10(200000), 25); % 2 kb to 200 kb
n_lengths = length(tss_to_pas_distances);

fprintf('\nGene length sweep:\n');
fprintf('  Range: %.1f kb to %.1f kb\n', tss_to_pas_distances(1) / 1000, ...
        tss_to_pas_distances(end) / 1000);
fprintf('  Number of points: %d\n', n_lengths);

% Pre-allocate results
CAD_results = zeros(n_scenarios, n_lengths);

% --- CALCULATE CAD FOR EACH SCENARIO ---
% Start parallel pool if not already running
if isempty(gcp('nocreate'))
    parpool;
end

for s = 1:n_scenarios
    fprintf('\n--- Scenario %d: %s ---\n', s, scenario_names{s});

    % Set k_on for this scenario
    P.kEon = kEon_scenarios(s);

    % Pre-allocate for this scenario (needed for parfor slicing)
    CAD_results_s = zeros(1, n_lengths);

    parfor i = 1 : n_lengths
        tss_to_pas = tss_to_pas_distances(i);
        P_local = P; % Copy P for this worker
        
        try
            % Set up gene geometry
            P_local.PASposition = tss_to_pas;
            P_local.geneLength_bp = tss_to_pas + after_PAS_length;

            % Run simulation and calculate termination profile
            [termination_profile, distances_bp] = calculate_termination_profile(P_local);

            % Calculate CAD at 50% threshold
            CAD_results_s(i) = calculate_CAD_from_profile(termination_profile, distances_bp, 0.5);

            fprintf('  Finished L = %.1f kb\n', tss_to_pas / 1000);

        catch ME
            warning('Failed for L = %.0f bp: %s', tss_to_pas, ME.message);
            CAD_results_s(i) = NaN;
        end
    end
    
    % Store results for this scenario
    CAD_results(s, :) = CAD_results_s;
    fprintf('  Scenario completed.\n');
end

fprintf('\nCalculations complete!\n');

% --- GENERATE PLOT ---
fprintf('\nGenerating Kinetic Sensitivity plot...\n');

figure('Position', [ 100, 100, 900, 600 ]);
hold on;

% Plot CAD curves for each scenario
for s = 1:n_scenarios
    semilogx(tss_to_pas_distances/1000, CAD_results(s, :)/1000, '-o', ...
        'LineWidth', 2.5, 'MarkerSize', 6, 'Color', scenario_colors{s}, ...
        'MarkerFaceColor', scenario_colors{s}, 'DisplayName', scenario_names{s});
end

% Formatting
xlabel('Gene Length (TSS to PAS, kb)', 'FontSize', 12);
ylabel('CAD (kb)', 'FontSize', 12);
title('Kinetic Sensitivity: CAD vs. Gene Length for Different k_{on}', ...
      'FontSize', 14, 'FontWeight', 'bold');
legend('Location', 'southeast');
grid on;
xlim([ tss_to_pas_distances(1) / 1000, tss_to_pas_distances(end) / 1000 ]);

hold off;

% Save figure
saveas(gcf, fullfile(outputDir, 'CAD_Kinetic_Sensitivity.png'));
fprintf('Figure saved: CAD_Kinetic_Sensitivity.png\n');

% --- ANALYSIS SUMMARY ---
fprintf('\n=== Analysis Summary ===\n');
for s = 1 : n_scenarios
    fprintf('\n%s:\n', scenario_names{s});
    fprintf('  CAD range: %.0f to %.0f bp\n', ...
        min(CAD_results(s, :)), max(CAD_results(s, :)));
    fprintf('  CAD for short genes (< 10 kb): %.0f \pm %.0f bp\n', ...
        mean(CAD_results(s, tss_to_pas_distances < 10000)), ...
        std(CAD_results(s, tss_to_pas_distances < 10000)));
    fprintf('  CAD for long genes (> 40 kb): %.0f \pm %.0f bp\n', ...
        mean(CAD_results(s, tss_to_pas_distances >= 40000)), ...
        std(CAD_results(s, tss_to_pas_distances >= 40000)));
end

fprintf('\nBiological interpretation:\n');
fprintf('  - Faster k_on shifts curve leftward (shorter runway needed)\n');
fprintf('  - Slower k_on shifts curve rightward (longer runway needed)\n');
fprintf('  - All curves plateau at long gene lengths\n');
fprintf('  - Demonstrates kinetic nature of gene length dependency\n');

fprintf('\nAnalysis complete!\n');

%% --- HELPER FUNCTIONS ---

function [termination_profile, distances_bp] = calculate_termination_profile(P)
    % Calculate termination profile for a specific gene
    global N PAS N_PAS Ef_ss;
    syms Ef real;

    % Run simulation
    [R_sol, REH_sol, P_sim] = run_single_gene_simulation(P);

    % Calculate termination profile
    flux_cleavage_per_node = P_sim.kc * REH_sol;
    flux_R_exit = P_sim.k_e * R_sol(end);
    flux_REH_exit = P_sim.k_e2 * REH_sol(end);
    total_outflux = sum(flux_cleavage_per_node) + flux_R_exit + flux_REH_exit;

    if total_outflux > 1e-9
        cumulative_exit_flux = cumsum(flux_cleavage_per_node);
        termination_profile = cumulative_exit_flux / total_outflux;
    else
        termination_profile = zeros(size(REH_sol));
    end

    % Distance from PAS
    nodes_post_pas = 1 : length(REH_sol);
    distances_bp = nodes_post_pas * P_sim.L_a;
end

function CAD = calculate_CAD_from_profile(termination_profile, distances_bp, threshold)
    % Calculate CAD for a given threshold
    
    if isempty(termination_profile) || all(isnan(termination_profile))
        CAD = NaN;
        return;
    end

    % Find where termination reaches threshold
    threshold_idx = find(termination_profile >= threshold, 1, 'first');

    if isempty(threshold_idx)
        CAD = distances_bp(end); % Never reaches threshold
    elseif threshold_idx == 1
        CAD = distances_bp(1);    % Reaches threshold immediately
    else
        % Interpolate for more precise CAD
        x1 = termination_profile(threshold_idx-1);
        x2 = termination_profile(threshold_idx);
        d1 = distances_bp(threshold_idx - 1);
        d2 = distances_bp(threshold_idx);

        CAD = d1 + (threshold - x1) * (d2 - d1) / (x2 - x1);
    end
end

function [R_sol, REH_sol, P_sim] = run_single_gene_simulation(P)
    % Run single gene simulation
    global N PAS N_PAS Ef_ss;
    syms Ef real;

    % Set up geometry
    L_a = P.L_a;
    N = floor(P.geneLength_bp / L_a);
    PAS = floor(P.PASposition / L_a);
    N_PAS = N - PAS + 1;

    % Compute steady states
    [r_E_BeforePas] = compute_steady_states(P, P.EBindingNumber + 1);

    % Set up kPon values with linear increase (like CPA_multipleE_main.m)
    kPon_vals = P.kPon_min + P.kPon_slope * (0 : N - 1);

    % Build RE_vals
    RE_vals = sym(zeros(P.EBindingNumber + 1, N));

    for e = 1:(P.EBindingNumber + 1)
        for idx = 1:N
            RE_vals(e, idx) = subs(r_E_BeforePas(e), {'kPon', 'kPoff'}, {kPon_vals(idx), P.kPoff_const});
        end
    end
    
    % Create E-binding function
    P.RE_val_bind_E = matlabFunction(simplify(sum(sym(1:P.EBindingNumber)' .* RE_vals(2:end, :), 1)), 'Vars', {'Ef'});
    
    % Solve system
    X_guess = 1e-6 * ones(N + N_PAS, 1);
    options = optimoptions('fsolve', 'Display', 'off', 'FunctionTolerance', 1e-8);
    
    % Initialize Ef_ss to the free E-factor pool determined in the main script setup
    Ef_ss = P.E_free;
    X_temp = fsolve(@(xx) ode_dynamics_multipleE(xx, P), X_guess, options);
    
    % Recalculate kHon
    avg_E_bound = P.RE_val_bind_E(Ef_ss);
    P.kHon = P.kHon * avg_E_bound(PAS);
    
    X_final = fsolve(@(xx) ode_dynamics_multipleE(xx, P), X_temp, options);
    
    R_sol = X_final(1:N);
    REH_sol = X_final((N+1):(N+N_PAS));
    
    P_sim = P;
end