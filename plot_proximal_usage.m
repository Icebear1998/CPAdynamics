% PLOT_PROXIMAL_USAGE.m
% Visualize Proximal PAS Usage vs. Inter-PAS Distance
%
% This script simulates the competition between a proximal PAS and a distal PAS
% by calculating the probability of termination at the proximal site as a
% function of the distance to the next (distal) site.
%
% Mechanism:
% - The "Termination Profile" (CDF) represents the cumulative probability of
% termination by a certain distance.
% - If a proximal site is located at distance 'd' from the TSS (or relative
% to a previous site), the flux that terminates *before* the next site
% represents the usage of that proximal site.
%
% Expected output: Usage curve increasing with inter-PAS distance
%
% WORKFLOW: This script uses self-consistent (R_free, E_free) from the
% gene length analysis workflow. If analysis results are not available,
% it falls back to fixed values (30% of total).

clear;
close all;
global N PAS N_PAS Ef_ss;

fprintf('=== Proximal PAS Usage Analysis ===\n');
fprintf('Calculating Proximal Usage vs. Inter-PAS Distance...\n\n');

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
P.kEon = 0.00025;     % E-factor binding rate
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

% Simulation geometry
% We simulate a single "long enough" gene to capture the full termination profile
P.PASposition = 20000;          % PAS at 20 kb (doesn't matter much for post-PAS profile)
after_PAS_length = 5000;        % 5 kb after PAS to capture full termination
P.geneLength_bp = P.PASposition + after_PAS_length;

fprintf('\nParameters:\n');
fprintf('  Simulating gene length: %.1f kb\n', P.geneLength_bp/1000);
fprintf('  Analyzing post-PAS usage...\n');

%% --- RUN SIMULATION ---
% Start parallel pool if not already running (good practice even for single run if we expand later)
if isempty(gcp('nocreate'))
    parpool;
end

fprintf('\nRunning simulation...\n');
try
    [termination_profile, distances_bp] = calculate_termination_profile(P);
    fprintf('Simulation complete!\n');
catch ME
    error('Simulation failed: %s', ME.message);
end

%% --- CALCULATE PROXIMAL USAGE PROBABILITY ---
% We'll interpolate the termination profile at various "inter-PAS distances"
% If the distal site is at distance 'd' from the proximal site, then
% Prob(Proximal Usage) ~ CDF(d)
% (This assumes the distal site is a hard stop or we are just asking 
% "what fraction terminated BEFORE reaching d?")

inter_pas_distances = 0:5:3000; % 0 to 3 kb scan
proximal_usage = interp1([0; distances_bp(:)], [0; termination_profile(:)], inter_pas_distances, 'linear', 'extrap');

% Cap at 1.0 (100%)
proximal_usage = min(proximal_usage, 1.0);

%% --- GENERATE PLOT ---
fprintf('\nGenerating Proximal Usage plot...\n');

figure('Position', [100, 100, 800, 600]);
hold on;

% Plot usage curve
plot(inter_pas_distances, proximal_usage * 100, 'r-', 'LineWidth', 3, 'DisplayName', 'Proximal PAS Usage');

% Highlight "Typical" Distance (e.g., 300 bp median)
typical_dist = 300;
typical_usage = interp1(inter_pas_distances, proximal_usage * 100, typical_dist);

plot([typical_dist, typical_dist], [0, typical_usage], 'k--', 'LineWidth', 1.5, 'DisplayName', 'Typical Distance (300 bp)');
plot([0, typical_dist], [typical_usage, typical_usage], 'k--', 'LineWidth', 1.5, 'HandleVisibility', 'off');
plot(typical_dist, typical_usage, 'ko', 'MarkerFaceColor', 'k', 'MarkerSize', 8, 'HandleVisibility', 'off');

text(typical_dist + 50, typical_usage - 5, sprintf('%.1f%%', typical_usage), 'FontSize', 12, 'FontWeight', 'bold');

% Formatting
xlabel('Inter-PAS Distance (bp)', 'FontSize', 12);
ylabel('Proximal Site Usage (%)', 'FontSize', 12);
title('Proximal PAS Sequence Usage vs. Downstream Distance', 'FontSize', 14, 'FontWeight', 'bold');
legend('Location', 'southeast');
grid on;
xlim([0, 3000]);
ylim([0, 100]);

hold off;

% Save figure
saveas(gcf, fullfile(outputDir, 'Proximal_Usage_vs_InterPAS.png'));
fprintf('Figure saved: Proximal_Usage_vs_InterPAS.png\n');

%% --- ANALYSIS SUMMARY ---
fprintf('\n=== Analysis Summary ===\n');
fprintf('Proximal Usage at key distances:\n');
fprintf('  100 bp:  %.1f%%\n', interp1(inter_pas_distances, proximal_usage*100, 100));
fprintf('  300 bp:  %.1f%% (Typical)\n', interp1(inter_pas_distances, proximal_usage*100, 300));
fprintf('  500 bp:  %.1f%%\n', interp1(inter_pas_distances, proximal_usage*100, 500));
fprintf('  1000 bp: %.1f%%\n', interp1(inter_pas_distances, proximal_usage*100, 1000));

fprintf('\nBiological interpretation:\n');
fprintf('  - Usage increases with distance to the next site (kinetic competition)\n');
fprintf('  - Closely spaced PAS sites (< 300 bp) show lower proximal usage\n');
fprintf('  - Distant sites allow full termination efficacy at the proximal site\n');

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
    nodes_post_pas = 1:length(REH_sol);
    distances_bp = nodes_post_pas * P_sim.L_a;
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
    kPon_vals = P.kPon_min + P.kPon_slope * (0:N-1);
    
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
    
    Ef_ss = P.E_free; % Initialize Ef_ss to the free E-factor pool
    X_temp = fsolve(@(xx) ode_dynamics_multipleE(xx, P), X_guess, options);
    
    avg_E_bound = P.RE_val_bind_E(Ef_ss);
    P.kHon = P.kHon * avg_E_bound(PAS);
    
    X_final = fsolve(@(xx) ode_dynamics_multipleE(xx, P), X_temp, options);
    
    R_sol = X_final(1:N);
    REH_sol = X_final((N+1):(N+N_PAS));
    
    P_sim = P;
end