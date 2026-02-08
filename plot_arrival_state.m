% PLOT_ARRIVAL_STATE.m
% Visualize E-factor occupancy at PAS vs. Gene Length
%
% This script demonstrates the "unpreparedness" mechanism for short genes:
% - Short genes: Pol II arrives at PAS with few E-factors bound
% - Long genes: Pol II arrives fully saturated with E-factors
%
% Expected output: Linear rise for short genes, plateau at N=5 for long genes
%
% WORKFLOW: This script uses self-consistent (R_free, E_free) from the
% gene length analysis workflow. If analysis results are not available,
% it falls back to fixed values (30% of total).

clear;
close all;
global N PAS N_PAS Ef_ss;

fprintf('=== Arrival State Analysis ===\n');
fprintf('Calculating E-factor occupancy at PAS vs. Gene Length...\n\n');

% --- SETUP OUTPUT DIRECTORY ---
outputDir = 'SecondVersionResults/GeneLengthAnalysis/';
if ~exist(outputDir, 'dir')
    mkdir(outputDir);
end

% --- LOAD SELF-CONSISTENT SOLUTION (if available) ---
fprintf('Attempting to load self-consistent (R_free, E_free) solution...\n');

analysis_dir = 'SecondVersionResults/GeneLengthAnalysis/';
use_self_consistent = false;
R_free_value = 21000;  % Default: 30% of 70,000
E_free_value = 21000;  % Default: 30% of 70,000

if exist(analysis_dir, 'dir')
    % Look for most recent analysis results
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
% Base parameters (consistent with CPA_multipleE_main.m)
P.L_a = 100;         % Node spacing (bp)
P.kEon = 0.00025;    % E-factor binding rate
P.kEoff = 10;        % E-factor unbinding rate
P.kHon = 0.2;        % H-factor binding rate (initial)
P.kHoff = 0.0125;    % H-factor unbinding rate
P.kc = 0.05;         % Cleavage rate
P.k_e = 65 / P.L_a;  % Pol II elongation rate
P.k_e2 = 30 / P.L_a; % REH elongation rate
P.k_in = 2;          % Pol II initiation rate
P.EBindingNumber = 5; % Max E-binding sites
P.kPon_min = 0.01;   % Min Ser2P phosphorylation (at TSS)
P.kPon_slope = 0.02; % Slope of linear increase
P.kPoff_const = 1;   % Ser2P dephosphorylation

% Free resource pools (from self-consistent solution or default)
P.Pol_free = R_free_value;
P.E_free = E_free_value;

% Total resource pools (needed by global ode_dynamics_multipleE.m)
P.Pol_total = 70000; % Total polymerase
P.E_total = 70000;   % Total E-factors

% Flags for global ode_dynamics_multipleE.m compatibility
P.FirstRun = false; % Use simple E_free calculation (not self-consistent solver)
P.is_unphysical = false;

% Fixed after-PAS region
after_PAS_length = 5000; % 5 kb after PAS

fprintf('\nParameters:\n');
fprintf('  E-binding sites (N): %d\n', P.EBindingNumber);
fprintf('  Free Pol II: %.0f\n', P.Pol_free);
fprintf('  Free E-factors: %.0f\n', P.E_free);
fprintf('  kPon: linear increase from %.3f (slope %.3f)\n', P.kPon_min, P.kPon_slope);

% --- GENE LENGTH SWEEP ---
% TSS-to-PAS distances to analyze
tss_to_pas_distances = logspace(log10(2000), log10(200000), 5); % 2 kb to 200 kb
n_lengths = length(tss_to_pas_distances);

fprintf('\nGene length sweep:\n');
fprintf('  Range: %.1f kb to %.1f kb\n', tss_to_pas_distances(1) / 1000, ...
        tss_to_pas_distances(end) / 1000);
fprintf('  Number of points: %d\n', n_lengths);

% Pre-allocate results
E_at_PAS = zeros(1, n_lengths);

% --- CALCULATE E-FACTOR OCCUPANCY AT PAS ---
fprintf('\nCalculating E-factor occupancy at PAS...\n');

% Start parallel pool if not already running
if isempty(gcp('nocreate'))
    parpool;
end

parfor i = 1 : n_lengths
    tss_to_pas = tss_to_pas_distances(i);

    % Create local parameter structure to avoid data race
    P_local = P;

    try
        % Set up gene geometry
        P_local.PASposition = tss_to_pas;
        P_local.geneLength_bp = tss_to_pas + after_PAS_length;

        % Run simulation
        [~, ~, avg_E_bound, ~] = run_single_gene_simulation(P_local);

        % Extract E-factor occupancy at PAS
        PAS_node = floor(P_local.PASposition / P_local.L_a);
        E_at_PAS(i) = avg_E_bound(PAS_node);

        fprintf('  Finished L = %.1f kb\n', tss_to_pas / 1000);

    catch ME
        warning('Failed for L = %.0f bp: %s', tss_to_pas, ME.message);
        E_at_PAS(i) = NaN;
    end
end

fprintf('Calculations complete!\n');

% --- GENERATE PLOT ---
fprintf('\nGenerating Arrival State plot...\n');

figure('Position', [ 100, 100, 800, 600 ]);
hold on;

% Plot E-factor occupancy at PAS
semilogx(tss_to_pas_distances / 1000, E_at_PAS, 'b-o', 'LineWidth', 2.5, ...
         'MarkerSize', 6, 'MarkerFaceColor', 'b', 'DisplayName', 'E-factors at PAS');

% Add horizontal line at maximum capacity
yline(P.EBindingNumber, 'r--', 'LineWidth', 2, 'DisplayName', ...
      sprintf('Maximum Capacity (N=%d)', P.EBindingNumber)); % Fixed DisplayName to use P.EBindingNumber

% Formatting
xlabel('Gene Length (TSS to PAS, kb)', 'FontSize', 12);
ylabel('Average E-Factors Bound at PAS', 'FontSize', 12);
title('Arrival State: E-Factor Occupancy at PAS vs. Gene Length', 'FontSize', ...
      14, 'FontWeight', 'bold');
legend('Location', 'southeast');
grid on;
xlim([ tss_to_pas_distances(1) / 1000, tss_to_pas_distances(end) / 1000 ]);
ylim([ 0, P.EBindingNumber + 0.5 ]);

hold off;

% Save figure
saveas(gcf, fullfile(outputDir, 'Arrival_State_vs_Gene_Length.png'));
fprintf('Figure saved: Arrival_State_vs_Gene_Length.png\n');

% --- ANALYSIS SUMMARY ---
fprintf('\n=== Analysis Summary ===\n');
fprintf('E-factor occupancy at PAS:\n');
fprintf('  Short genes (< 10 kb): %.2f \pm %.2f\n', ...
        mean(E_at_PAS(tss_to_pas_distances < 10000)), ...
        std(E_at_PAS(tss_to_pas_distances < 10000)));
fprintf('  Medium genes (10-40 kb): %.2f \pm %.2f\n', ...
        mean(E_at_PAS(tss_to_pas_distances >= 10000 & tss_to_pas_distances < 40000)), ...
        std(E_at_PAS(tss_to_pas_distances >= 10000 & tss_to_pas_distances < 40000)));
fprintf('  Long genes (> 40 kb): %.2f \pm %.2f\n', ...
        mean(E_at_PAS(tss_to_pas_distances >= 40000)), ...
        std(E_at_PAS(tss_to_pas_distances >= 40000)));

fprintf('\nBiological interpretation:\n');
fprintf('  - Short genes: Pol II arrives "unprepared" with insufficient E-factors\n');
fprintf('  - Long genes: Pol II arrives fully saturated (plateau at N=%d)\n', ...
        P.EBindingNumber);
fprintf('  - This explains the gene length dependency of CAD\n');

fprintf('\nAnalysis complete!\n');

% --- HELPER FUNCTION ---
function [R_sol, REH_sol, avg_E_bound, P_sim] = run_single_gene_simulation(P)
    % Run single gene simulation and return E-factor binding profile
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
    
    % Solve system (two-step like main script)
    X_guess = 1e-6 * ones(N + N_PAS, 1);
    options = optimoptions('fsolve', 'Display', 'off', 'FunctionTolerance', 1e-8);
    
    % First solve (using fixed E_free)
    Ef_ss = P.E_free; % Use the global E_free determined in the main script setup
    X_temp = fsolve(@(xx) ode_dynamics_multipleE(xx, P), X_guess, options);
    
    % Recalculate kHon (renormalization)
    avg_E_bound = P.RE_val_bind_E(Ef_ss);
    P.kHon = P.kHon * avg_E_bound(PAS);
    
    % Final solve
    X_final = fsolve(@(xx) ode_dynamics_multipleE(xx, P), X_temp, options);
    
    R_sol = X_final(1:N);
    REH_sol = X_final((N+1):(N+N_PAS));
    
    % Recalculate avg_E_bound with final Ef_ss
    avg_E_bound = P.RE_val_bind_E(Ef_ss);
    
    P_sim = P;
end