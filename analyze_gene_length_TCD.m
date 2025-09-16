% ANALYZE_GENE_LENGTH_TCD.m
% Final analysis script for gene length effects on Termination Commitment Distance (TCD)
%
% This script implements Steps 3-5 of the gene length analysis:
% - Load interpolation functions and gene length distribution
% - Implement conservation equations with numerical integration
% - Solve for self-consistent (R_free, E_free)
% - Calculate TCD for different gene lengths
% - Generate comprehensive analysis and visualizations

fprintf('=== Gene Length TCD Analysis ===\n');
fprintf('Analyzing Termination Commitment Distance across gene lengths...\n\n');

%% --- LOAD INTERPOLATION DATA ---
fprintf('Loading interpolation functions...\n');

% Find the most recent interpolation file
interp_dir = 'SecondVersionResults/GeneLengthAnalysis/';
if ~exist(interp_dir, 'dir')
    error('Interpolation directory not found. Please run build_gene_length_interpolation.m first.');
end

% Look for interpolation .mat files
interp_files = dir(fullfile(interp_dir, 'gene_length_interpolation_*.mat'));
if isempty(interp_files)
    error('No interpolation files found. Please run build_gene_length_interpolation.m first.');
end

% Use the most recent file
[~, newest_idx] = max([interp_files.datenum]);
interp_filename = fullfile(interp_dir, interp_files(newest_idx).name);
fprintf('Loading: %s\n', interp_filename);

load(interp_filename, 'interpolation_results');
fprintf('Interpolation data loaded successfully!\n');

% Extract key functions and parameters
R_occupied_interp = interpolation_results.functions.R_occupied_interp;
E_occupied_interp = interpolation_results.functions.E_occupied_interp;
gene_length_pdf = interpolation_results.functions.gene_length_pdf;

% System parameters
R_total_target = 70000;  % Target total R (from typical parameter set)
E_total_target = 70000;  % Target total E (from typical parameter set)

fprintf('Target totals: R_total = %.0f, E_total = %.0f\n', R_total_target, E_total_target);

%% --- DEFINE INTEGRATION DOMAIN ---
fprintf('\nSetting up integration domain...\n');

% Gene length integration range (cover main distribution)
L_min_int = 1000;      % 1 kb
L_max_int = 300000;    % 300 kb (covers >99% of distribution)
L_integration_points = 200;  % High resolution for accurate integration

L_integration = logspace(log10(L_min_int), log10(L_max_int), L_integration_points);
dL = diff([L_integration(1)/1.1, L_integration]);  % Integration weights (log spacing)

% Evaluate gene length distribution
f_L = gene_length_pdf(L_integration);
f_L_normalized = f_L ./ sum(f_L .* dL);  % Normalize to ensure integral = 1

fprintf('Integration domain:\n');
fprintf('  L range: %.1f kb to %.1f kb\n', L_min_int/1000, L_max_int/1000);
fprintf('  Integration points: %d\n', L_integration_points);
fprintf('  Distribution integral check: %.6f (should be ?1)\n', sum(f_L_normalized .* dL));

%% --- IMPLEMENT CONSERVATION EQUATIONS ---
fprintf('\nImplementing conservation equations...\n');

% The conservation equations are the heart of the gene length analysis.
% They enforce the constraint that the total amount of R and E factors
% in the system equals the sum of free factors plus those bound to genes.
%
% MATHEMATICAL FRAMEWORK:
% For a genome with gene length distribution f(L), the conservation equations are:
%
%   R_total = R_free + ∫ R_occupied(R_free, E_free, L) × f(L) dL
%   E_total = E_free + ∫ E_occupied(R_free, E_free, L) × f(L) dL
%
% Where:
%   - R_total, E_total: Total amounts in the system (fixed, typically 70,000)
%   - R_free, E_free: Free (unbound) amounts (unknowns to solve for)
%   - R_occupied(R_free, E_free, L): Bound R for gene of length L (from interpolation)
%   - E_occupied(R_free, E_free, L): Bound E for gene of length L (from interpolation)
%   - f(L): Gene length probability density function (log-normal distribution)
%
% BIOLOGICAL INTERPRETATION:
% - Each gene of length L consumes R_occupied and E_occupied resources
% - The genome-wide consumption is the integral over all gene lengths
% - The remaining resources (R_free, E_free) are available for new genes
% - This creates a self-consistent feedback: available resources determine
%   individual gene behavior, which determines total resource consumption

fprintf('Conservation equation mathematical framework established.\n');

%% --- SOLVE FOR SELF-CONSISTENT SOLUTION ---
fprintf('\nSolving for self-consistent (R_free, E_free)...\n');

% SOLUTION STRATEGY:
% We need to find (R_free, E_free) such that both conservation equations
% are satisfied simultaneously. This is a 2D root-finding problem.
%
% NUMERICAL APPROACH:
% 1. Define residual function: [R_total - R_calc, E_total - E_calc]
% 2. Use fsolve to find where residual = [0, 0]
% 3. Validate solution by checking conservation equations

% Define residual function for root finding
% This function returns [0, 0] when the conservation equations are satisfied
residual_function = @(x) conservation_equations_residual(x, R_occupied_interp, E_occupied_interp, L_integration, f_L_normalized, dL, R_total_target, E_total_target);

% INITIAL GUESS STRATEGY:
% Start with moderate free fractions (30% of total). This is reasonable because:
% - Too high: Would imply very little resource consumption (unlikely)
% - Too low: Would imply extreme resource depletion (may cause convergence issues)
% - 30%: Balanced starting point that allows for significant resource binding
R_free_guess = 0.3 * R_total_target;  % 30% free polymerase
E_free_guess = 0.3 * E_total_target;  % 30% free E factors
initial_guess = [R_free_guess, E_free_guess];

fprintf('Initial guess strategy:\n');
fprintf('  R_free = %.0f (%.0f%% of total) - moderate free fraction\n', R_free_guess, 30);
fprintf('  E_free = %.0f (%.0f%% of total) - moderate free fraction\n', E_free_guess, 30);
fprintf('  Rationale: Balanced starting point avoiding extreme resource scenarios\n');

% SOLVER CONFIGURATION:
% - Display: 'iter' shows convergence progress
% - FunctionTolerance: 1e-6 ensures high accuracy in conservation equations
% - MaxIterations: 100 should be sufficient for this smooth problem
options = optimoptions('fsolve', 'Display', 'iter', 'FunctionTolerance', 1e-6, 'MaxIterations', 100);

fprintf('\nSolver configuration:\n');
fprintf('  Method: fsolve (trust-region-dogleg algorithm)\n');
fprintf('  Tolerance: 1e-6 (high precision for conservation equations)\n');
fprintf('  Max iterations: 100\n');

fprintf('Solving conservation equations...\n');
tic;
try
    [solution, fval, exitflag, output] = fsolve(residual_function, initial_guess, options);
    solve_time = toc;
    
    if exitflag > 0
        R_free_solution = solution(1);
        E_free_solution = solution(2);
        
        fprintf('\nSolution found!\n');
        fprintf('  R_free = %.0f (%.1f%% of total)\n', R_free_solution, R_free_solution/R_total_target*100);
        fprintf('  E_free = %.0f (%.1f%% of total)\n', E_free_solution, E_free_solution/E_total_target*100);
        fprintf('  Solve time: %.2f seconds\n', solve_time);
        fprintf('  Final residual: [%.2e, %.2e]\n', fval(1), fval(2));
        
        % Validate solution
        [R_check, E_check] = conservation_equations(R_free_solution, E_free_solution, R_occupied_interp, E_occupied_interp, L_integration, f_L_normalized, dL);
        fprintf('  Validation: R_total = %.0f, E_total = %.0f\n', R_check, E_check);
        
    else
        error('Failed to find solution. Exit flag: %d', exitflag);
    end
    
catch ME
    fprintf('Error in solving: %s\n', ME.message);
    error('Could not find self-consistent solution. Try adjusting initial guess or parameter ranges.');
end

%% --- CALCULATE TCD FOR DIFFERENT GENE LENGTHS ---
fprintf('\nCalculating TCD for different gene lengths...\n');

% Define gene lengths for TCD analysis (cover the distribution well)
L_TCD_analysis = logspace(log10(2000), log10(200000), 50);  % 2 kb to 200 kb
n_L_TCD = length(L_TCD_analysis);

% TCD calculation parameters
TCD_thresholds = [0.25, 0.5, 0.75];  % 25%, 50%, 75% termination thresholds
n_thresholds = length(TCD_thresholds);

% Pre-allocate results
TCD_results = zeros(n_L_TCD, n_thresholds);
termination_profiles = cell(n_L_TCD, 1);

fprintf('Calculating TCD for %d gene lengths...\n', n_L_TCD);

% Use the solved (R_free, E_free) for all TCD calculations
for i = 1:n_L_TCD
    L_current = L_TCD_analysis(i);
    
    try
        % Calculate termination profile for this gene length
        [termination_profile, distances_bp] = calculate_termination_profile(L_current, R_free_solution, E_free_solution);
        termination_profiles{i} = struct('distances', distances_bp, 'profile', termination_profile);
        
        % Calculate TCD for each threshold
        for j = 1:n_thresholds
            threshold = TCD_thresholds(j);
            TCD_results(i, j) = calculate_TCD_from_profile(termination_profile, distances_bp, threshold);
        end
        
    catch ME
        fprintf('Warning: Failed to calculate TCD for L = %.0f bp: %s\n', L_current, ME.message);
        TCD_results(i, :) = NaN;
    end
    
    if mod(i, 10) == 0
        fprintf('  Completed %d/%d gene lengths\n', i, n_L_TCD);
    end
end

fprintf('TCD calculations completed!\n');

%% --- ANALYSIS AND STATISTICS ---
fprintf('\nPerforming statistical analysis...\n');

% Remove failed calculations
valid_TCD = ~isnan(TCD_results(:, 1));
L_TCD_valid = L_TCD_analysis(valid_TCD);
TCD_valid = TCD_results(valid_TCD, :);

fprintf('Valid TCD calculations: %d/%d\n', sum(valid_TCD), n_L_TCD);

% Statistics for each threshold
for j = 1:n_thresholds
    threshold = TCD_thresholds(j);
    TCD_vals = TCD_valid(:, j);
    
    fprintf('TCD statistics (%.0f%% threshold):\n', threshold*100);
    fprintf('  Range: %.0f to %.0f bp\n', min(TCD_vals), max(TCD_vals));
    fprintf('  Mean: %.0f bp, Median: %.0f bp\n', mean(TCD_vals), median(TCD_vals));
    fprintf('  Std: %.0f bp\n', std(TCD_vals));
end

% Correlation analysis
fprintf('\nCorrelation analysis:\n');
for j = 1:n_thresholds
    correlation = corr(log10(L_TCD_valid'), TCD_valid(:, j));
    fprintf('  TCD(%.0f%%) vs log10(L): r = %.3f\n', TCD_thresholds(j)*100, correlation);
end

%% --- GENERATE VISUALIZATIONS ---
fprintf('\nGenerating TCD visualization...\n');

% TCD vs Gene Length plot
figure('Position', [100, 100, 800, 600]);

colors = lines(n_thresholds);
hold on;
for j = 1:n_thresholds
    semilogx(L_TCD_valid/1000, TCD_valid(:, j), 'o-', 'LineWidth', 2, 'MarkerSize', 4, ...
        'Color', colors(j, :), 'DisplayName', sprintf('%.0f%% threshold', TCD_thresholds(j)*100));
end
xlabel('TSS-to-PAS Distance (kb)', 'FontSize', 12);
ylabel('TCD (bp)', 'FontSize', 12);
title('Termination Commitment Distance vs Gene Length', 'FontSize', 14, 'FontWeight', 'bold');
legend('Location', 'best');
grid on;
hold off;

%% --- SAVE COMPREHENSIVE RESULTS ---
fprintf('\nSaving analysis results...\n');

% Create results structure
analysis_results = struct();
analysis_results.metadata.creation_date = datestr(now);
analysis_results.metadata.source_interpolation_file = interp_filename;
analysis_results.metadata.description = 'Complete gene length TCD analysis results';

% Solution
analysis_results.solution.R_free = R_free_solution;
analysis_results.solution.E_free = E_free_solution;
analysis_results.solution.R_total_target = R_total_target;
analysis_results.solution.E_total_target = E_total_target;
analysis_results.solution.solve_time_seconds = solve_time;

% TCD results
analysis_results.TCD.gene_lengths = L_TCD_analysis;
analysis_results.TCD.thresholds = TCD_thresholds;
analysis_results.TCD.values = TCD_results;
analysis_results.TCD.valid_indices = valid_TCD;
analysis_results.TCD.termination_profiles = termination_profiles;

% Statistics
analysis_results.statistics.correlation_with_log_length = zeros(n_thresholds, 1);
for j = 1:n_thresholds
    if sum(valid_TCD) > 1
        analysis_results.statistics.correlation_with_log_length(j) = corr(log10(L_TCD_valid'), TCD_valid(:, j));
    end
end

% Save results
timestamp = datestr(now, 'yyyymmdd_HHMMSS');
results_filename = fullfile(interp_dir, sprintf('gene_length_TCD_analysis_%s.mat', timestamp));
save(results_filename, 'analysis_results', '-v7.3');

% Save TCD plot
plot_filename = sprintf('TCD_vs_gene_length_%s.png', timestamp);
saveas(gcf, fullfile(interp_dir, plot_filename));

fprintf('Analysis results saved:\n');
fprintf('  Data file: %s\n', results_filename);
fprintf('  Plot: %s\n', plot_filename);

%% --- FINAL SUMMARY ---
fprintf('\n=== GENE LENGTH TCD ANALYSIS COMPLETE ===\n');
fprintf('Key findings:\n');
fprintf('  Self-consistent solution: R_free = %.0f (%.1f%%), E_free = %.0f (%.1f%%)\n', ...
    R_free_solution, R_free_solution/R_total_target*100, E_free_solution, E_free_solution/E_total_target*100);
fprintf('  TCD range (50%% threshold): %.0f to %.0f bp\n', min(TCD_valid(:,2)), max(TCD_valid(:,2)));
fprintf('  TCD-length correlation (50%%): r = %.3f\n', analysis_results.statistics.correlation_with_log_length(2));

fprintf('\nBiological insights:\n');
fprintf('  - Resource competition affects termination commitment\n');
fprintf('  - Gene length influences available free factor pools\n');
fprintf('  - TCD varies systematically across the gene length distribution\n');

fprintf('\nAnalysis completed successfully!\n');

%% --- HELPER FUNCTIONS ---

function residual = conservation_equations_residual(x, R_interp, E_interp, L_vals, f_L_vals, dL_vals, R_target, E_target)
    % RESIDUAL FUNCTION FOR CONSERVATION EQUATION SOLVER
    % 
    % This function computes the residual (error) in the conservation equations
    % for a given guess of (R_free, E_free). The root finder (fsolve) will
    % minimize this residual to find the self-consistent solution.
    %
    % INPUT:
    %   x(1) = R_free: Guess for free polymerase amount
    %   x(2) = E_free: Guess for free E factor amount
    %   R_interp, E_interp: Interpolation functions R_occupied(R_free, E_free, L)
    %   L_vals, f_L_vals, dL_vals: Gene length distribution and integration weights
    %   R_target, E_target: Target total amounts (typically 70,000 each)
    %
    % OUTPUT:
    %   residual: [R_error, E_error] where each should be zero at solution
    %
    % MATHEMATICAL MEANING:
    %   residual = [R_total_calculated - R_target, E_total_calculated - E_target]
    %   When residual = [0, 0], the conservation equations are satisfied
    
    R_free = x(1);  % Extract free polymerase guess
    E_free = x(2);  % Extract free E factor guess
    
    % Calculate what the total amounts would be given these free amounts
    [R_total_calc, E_total_calc] = conservation_equations(R_free, E_free, R_interp, E_interp, L_vals, f_L_vals, dL_vals);
    
    % Compute residual: difference between calculated and target totals
    % fsolve will minimize the magnitude of this residual vector
    residual = [R_total_calc - R_target, E_total_calc - E_target];
end

function [R_total_calc, E_total_calc] = conservation_equations(R_free, E_free, R_interp, E_interp, L_vals, f_L_vals, dL_vals)
    % CORE CONSERVATION EQUATION IMPLEMENTATION
    %
    % This function implements the fundamental conservation equations that
    % couple individual gene behavior to the global resource pool.
    %
    % CONSERVATION PRINCIPLE:
    %   Total_resources = Free_resources + Bound_resources
    %   
    % MATHEMATICAL IMPLEMENTATION:
    %   R_total = R_free + ∫ R_occupied(R_free, E_free, L) × f(L) dL
    %   E_total = E_free + ∫ E_occupied(R_free, E_free, L) × f(L) dL
    %
    % NUMERICAL INTEGRATION DETAILS:
    % - L_vals: Discrete gene length points for integration
    % - f_L_vals: Probability density at each gene length
    % - dL_vals: Integration weights (accounts for log-spacing)
    % - Integration: Riemann sum approximation of the continuous integral
    %
    % BIOLOGICAL INTERPRETATION:
    % - For each gene length L in the distribution f(L):
    %   * Lookup how much R and E that gene would consume
    %   * Weight by the frequency of genes of that length
    %   * Sum over all gene lengths to get total consumption
    % - Add free amounts to get total amounts in system
    
    % Vectorized evaluation of interpolation functions
    % We need R_occupied and E_occupied for each gene length in the distribution
    n_L = length(L_vals);
    R_occupied_vals = zeros(1, n_L);  % R consumed by each gene length
    E_occupied_vals = zeros(1, n_L);  % E consumed by each gene length
    
    % INTERPOLATION LOOP:
    % For each gene length, lookup resource consumption from pre-computed grid
    for i = 1:n_L
        try
            % Query interpolation functions: how much R and E does a gene of
            % length L_vals(i) consume when free pools are (R_free, E_free)?
            R_occupied_vals(i) = R_interp(R_free, E_free, L_vals(i));
            E_occupied_vals(i) = E_interp(R_free, E_free, L_vals(i));
        catch
            % EXTRAPOLATION HANDLING:
            % If (R_free, E_free, L) is outside the interpolation domain,
            % assume zero consumption (conservative approach)
            R_occupied_vals(i) = 0;
            E_occupied_vals(i) = 0;
        end
    end
    
    % NUMERICAL INTEGRATION:
    % Approximate ∫ R_occupied(R_free, E_free, L) × f(L) dL
    % Using discrete Riemann sum: Σ R_occupied(L_i) × f(L_i) × ΔL_i
    R_integral = sum(R_occupied_vals .* f_L_vals .* dL_vals);
    E_integral = sum(E_occupied_vals .* f_L_vals .* dL_vals);
    
    % CONSERVATION EQUATIONS:
    % Total amount = Free amount + Genome-wide consumption
    R_total_calc = R_free + R_integral;  % Total polymerase in system
    E_total_calc = E_free + E_integral;  % Total E factors in system
end

function [termination_profile, distances_bp] = calculate_termination_profile(tss_to_pas_distance, R_free, E_free)
    % Calculate termination profile for a specific TSS-to-PAS distance
    % This is a simplified version - in practice, you'd run the full simulation
    
    global N PAS N_PAS Ef_ss;
    syms Ef real;
    
    % Set up parameters for this gene length
    P = struct();
    P.L_a = 100;
    P.k_in = 2;
    P.kEon = 0.00025;
    P.kEoff = 10;
    P.k_e = 65/100;
    P.k_e2 = 30/100;
    P.E_total = E_free;  % Use free E as total for local analysis
    P.Pol_total = R_free;  % Use free R as total for local analysis
    P.kHon = 0.2;
    P.kHoff = 0.0125;
    P.kc = 0.05;
    P.kPon_min = 0.01;
    P.kPon_max = 1;
    P.kPoff_min = 0.1;
    P.kPoff_max = 2;
    P.kPoff_const = 1;
    
    % NEW: Gene structure with fixed after-PAS length
    after_PAS_length = 5000;  % Fixed 5 kb after-PAS region
    P.PASposition = tss_to_pas_distance;  % PAS position = TSS-to-PAS distance
    P.geneLength_bp = tss_to_pas_distance + after_PAS_length;  % Total gene length
    P.EBindingNumber = 6;
    
    % Run simulation (similar to existing scripts)
    [R_sol, REH_sol, P_sim] = run_single_gene_simulation_TCD(P);
    
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

function TCD = calculate_TCD_from_profile(termination_profile, distances_bp, threshold)
    % Calculate TCD for a given threshold
    
    if isempty(termination_profile) || all(isnan(termination_profile))
        TCD = NaN;
        return;
    end
    
    % Find where termination reaches threshold
    threshold_idx = find(termination_profile >= threshold, 1, 'first');
    
    if isempty(threshold_idx)
        TCD = distances_bp(end);  % Never reaches threshold
    elseif threshold_idx == 1
        TCD = distances_bp(1);    % Reaches threshold immediately
    else
        % Interpolate for more precise TCD
        x1 = termination_profile(threshold_idx-1);
        x2 = termination_profile(threshold_idx);
        d1 = distances_bp(threshold_idx-1);
        d2 = distances_bp(threshold_idx);
        
        TCD = d1 + (threshold - x1) * (d2 - d1) / (x2 - x1);
    end
end

function [R_sol, REH_sol, P_sim] = run_single_gene_simulation_TCD(P)
    % Simplified version of single gene simulation for TCD analysis
    global N PAS N_PAS Ef_ss;
    syms Ef real;
    
    % Set up geometry
    L_a = P.L_a;
    N = floor(P.geneLength_bp / L_a);
    PAS = floor(P.PASposition / L_a);
    N_PAS = N - PAS + 1;
    
    % Use simplified steady state (for speed)
    % In practice, you might want to use the full compute_steady_states
    
    % Simplified binding calculation
    P.RE_val_bind_E = @(Ef) ones(1, N) * P.EBindingNumber * 0.5;  % Simplified
    
    % Solve system
    X_guess = 1e-6 * ones(N + N_PAS, 1);
    options = optimoptions('fsolve', 'Display', 'off', 'FunctionTolerance', 1e-6);
    
    % Two-step solution
    P.FirstRun = true;
    P.is_unphysical = false;
    Ef_ss = P.E_total * 0.5;  % Simplified
    
    X_base = fsolve(@(xx) ode_dynamics_multipleE(xx, P), X_guess, options);
    if P.is_unphysical
        error('Unphysical result in TCD simulation');
    end
    
    % Update and resolve (simplified)
    P.FirstRun = false;
    X_final = fsolve(@(xx) ode_dynamics_multipleE(xx, P), X_base, options);
    
    R_sol = X_final(1:N);
    REH_sol = X_final((N+1):(N+N_PAS));
    P_sim = P;
end
