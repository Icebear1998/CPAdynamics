% GENERATE_GENE_LENGTH_GRID.m
% Generates lookup data for R_occupied and E_occupied as functions of
% (R_free, E_free, L) for the gene length analysis
%
% This script implements Step 1 of the gene length analysis:
% - Grid Generation: Create 3D grid of (R_free, E_free, L)
% - Parallel Computation: Run single-gene model for each grid point
% - Data Saving: Store results for later interpolation

% Start parallel pool if not already running
if isempty(gcp('nocreate')); parpool; end

fprintf('=== Gene Length Grid Generation ===\n');
fprintf('Generating lookup data for R_occupied and E_occupied...\n\n');

%% --- PARAMETER RANGES ---

% Base totals from typical simulations
R_total_base = 70000;
E_total_base = 70000;

% R_free range: should be much larger than 1/70000 of total
% Using range from ~10% to ~90% of total (well above the 1/70000 threshold)
R_free_min = 0.1 * R_total_base;   % 7,000
R_free_max = 0.9 * R_total_base;   % 63,000
R_free_points = 15;                % Resolution for R_free

% E_free range: similar logic to R_free
E_free_min = 0.1 * E_total_base;   % 7,000  
E_free_max = 0.9 * E_total_base;   % 63,000
E_free_points = 15;                % Resolution for E_free

% Gene length range: based on human gene distribution (including introns)
% From the histogram: spans ~10^2 to ~10^6 bp, with most genes 10^3 to 10^5
L_min = 1000;      % 1 kb (lower end of main distribution)
L_max = 200000;    % 200 kb (covers most genes, excludes extreme outliers)
L_points = 25;     % Resolution for gene length

% Create parameter grids
R_free_values = linspace(R_free_min, R_free_max, R_free_points);
E_free_values = linspace(E_free_min, E_free_max, E_free_points);
L_values = logspace(log10(L_min), log10(L_max), L_points); % Log spacing for gene lengths

fprintf('Parameter Ranges:\n');
fprintf('  R_free: %.0f to %.0f (%d points)\n', R_free_min, R_free_max, R_free_points);
fprintf('  E_free: %.0f to %.0f (%d points)\n', E_free_min, E_free_max, E_free_points);
fprintf('  L: %.0f to %.0f bp (%d points, log-spaced)\n', L_min, L_max, L_points);
fprintf('  Total grid points: %d\n\n', R_free_points * E_free_points * L_points);

%% --- BASE PARAMETERS ---
% Use standard parameter set from existing analyses
P_base.L_a = 100;
P_base.k_in = 2;
P_base.kEon = 0.00025;
P_base.kEoff = 10;
P_base.k_e = 65/100;
P_base.k_e2 = 30/100;
P_base.kHon = 0.2;
P_base.kHoff = 0.0125;
P_base.kc = 0.05;
P_base.kPon_min = 0.01;
P_base.kPon_max = 1;
P_base.kPoff_min = 0.1;
P_base.kPoff_max = 2;
P_base.kPoff_const = 1;
P_base.EBindingNumber = 6;  % Use standard value

%% --- GRID GENERATION ---
fprintf('Generating parameter grid...\n');

% Create all combinations of parameters
[R_grid, E_grid, L_grid] = meshgrid(R_free_values, E_free_values, L_values);

% Flatten to vectors for parallel processing
R_free_vec = R_grid(:);
E_free_vec = E_grid(:);
L_vec = L_grid(:);
n_points = length(R_free_vec);

% Pre-allocate results
R_occupied_vec = zeros(n_points, 1);
E_occupied_vec = zeros(n_points, 1);
success_flag = zeros(n_points, 1);  % Track successful simulations

fprintf('Grid generated: %d total points\n', n_points);
fprintf('Starting parallel computation...\n\n');

%% --- PARALLEL COMPUTATION ---
tic;

parfor i = 1:n_points
    try
        % Extract parameters for this grid point
        R_free_i = R_free_vec(i);
        E_free_i = E_free_vec(i);
        L_i = L_vec(i);
        
        % Set up parameters for this simulation
        P_i = P_base;
        P_i.Pol_total = R_free_i;  % Use R_free as Pol_total for local analysis
        P_i.E_total = E_free_i;    % Use E_free as E_total for local analysis
        P_i.geneLength_bp = L_i;
        P_i.PASposition = max(1000, L_i - 5000);  % PAS near end, but not too close
        
        % Run single-gene simulation
        [R_sol, REH_sol, ~] = run_single_gene_simulation(P_i);
        
        % Calculate occupied amounts
        R_occupied_i = sum(R_sol) + sum(REH_sol);  % Total bound polymerase
        
        % Calculate E occupied (need to run binding calculation)
        E_occupied_i = calculate_E_occupied(P_i, R_sol, REH_sol);
        
        % Store results
        R_occupied_vec(i) = R_occupied_i;
        E_occupied_vec(i) = E_occupied_i;
        success_flag(i) = 1;
        
    catch ME
        % Handle failed simulations
        R_occupied_vec(i) = NaN;
        E_occupied_vec(i) = NaN;
        success_flag(i) = 0;
        
        % Optionally log the error (commented out to avoid parfor issues)
        % fprintf('Failed at point %d: %s\n', i, ME.message);
    end
    
    % Progress indication (every 1000 points)
    if mod(i, 1000) == 0
        fprintf('Completed %d/%d points\n', i, n_points);
    end
end

computation_time = toc;
successful_points = sum(success_flag);
success_rate = successful_points / n_points * 100;

fprintf('\nParallel computation completed!\n');
fprintf('  Total time: %.1f minutes\n', computation_time/60);
fprintf('  Successful simulations: %d/%d (%.1f%%)\n', successful_points, n_points, success_rate);

%% --- DATA ORGANIZATION ---
fprintf('\nOrganizing results...\n');

% Create results structure
results = struct();
results.metadata.creation_date = datestr(now);
results.metadata.computation_time_minutes = computation_time/60;
results.metadata.success_rate = success_rate;
results.metadata.description = 'Lookup data for gene length analysis: R_occupied and E_occupied as functions of (R_free, E_free, L)';

% Parameter information
results.parameters.R_free_range = [R_free_min, R_free_max];
results.parameters.E_free_range = [E_free_min, E_free_max];
results.parameters.L_range = [L_min, L_max];
results.parameters.R_free_points = R_free_points;
results.parameters.E_free_points = E_free_points;
results.parameters.L_points = L_points;
results.parameters.base_parameters = P_base;

% Grid data
results.grid.R_free_values = R_free_values;
results.grid.E_free_values = E_free_values;
results.grid.L_values = L_values;

% Results data
results.data.R_free_vec = R_free_vec;
results.data.E_free_vec = E_free_vec;
results.data.L_vec = L_vec;
results.data.R_occupied_vec = R_occupied_vec;
results.data.E_occupied_vec = E_occupied_vec;
results.data.success_flag = success_flag;

%% --- SAVE RESULTS ---
fprintf('Saving results...\n');

% Create output directory
output_dir = 'SecondVersionResults/GeneLengthAnalysis/';
if ~exist(output_dir, 'dir')
    mkdir(output_dir);
end

% Save MATLAB data file
timestamp = datestr(now, 'yyyymmdd_HHMMSS');
mat_filename = fullfile(output_dir, sprintf('gene_length_grid_data_%s.mat', timestamp));
save(mat_filename, 'results', '-v7.3');  % Use v7.3 for large files

% Save text file with summary and data
txt_filename = fullfile(output_dir, sprintf('gene_length_grid_data_%s.txt', timestamp));
fid = fopen(txt_filename, 'w');

% Header
fprintf(fid, '%% Gene Length Analysis - Grid Data\n');
fprintf(fid, '%% Generated on: %s\n', results.metadata.creation_date);
fprintf(fid, '%% Computation time: %.1f minutes\n', results.metadata.computation_time_minutes);
fprintf(fid, '%% Success rate: %.1f%%\n', results.metadata.success_rate);
fprintf(fid, '%% \n');

% Parameters
fprintf(fid, '%% Parameters:\n');
fprintf(fid, '%% R_free range: %.0f to %.0f (%d points)\n', R_free_min, R_free_max, R_free_points);
fprintf(fid, '%% E_free range: %.0f to %.0f (%d points)\n', E_free_min, E_free_max, E_free_points);
fprintf(fid, '%% L range: %.0f to %.0f bp (%d points)\n', L_min, L_max, L_points);
fprintf(fid, '%% Total grid points: %d\n', n_points);
fprintf(fid, '%% \n');

% Base parameters
fprintf(fid, '%% Base Parameters:\n');
fprintf(fid, '%% L_a = %g, k_in = %g, kEon = %g, kEoff = %g\n', P_base.L_a, P_base.k_in, P_base.kEon, P_base.kEoff);
fprintf(fid, '%% k_e = %g, k_e2 = %g, kHon = %g, kHoff = %g\n', P_base.k_e, P_base.k_e2, P_base.kHon, P_base.kHoff);
fprintf(fid, '%% kc = %g, EBindingNumber = %g\n', P_base.kc, P_base.EBindingNumber);
fprintf(fid, '%% \n');

% Data columns
fprintf(fid, '%% Data Format (columns):\n');
fprintf(fid, '%% 1: R_free, 2: E_free, 3: L, 4: R_occupied, 5: E_occupied, 6: success_flag\n');
fprintf(fid, '%% \n');

% Write data
for i = 1:n_points
    fprintf(fid, '%.6f\t%.6f\t%.6f\t%.6f\t%.6f\t%d\n', ...
        R_free_vec(i), E_free_vec(i), L_vec(i), ...
        R_occupied_vec(i), E_occupied_vec(i), success_flag(i));
end

fclose(fid);

fprintf('Results saved:\n');
fprintf('  MATLAB file: %s\n', mat_filename);
fprintf('  Text file: %s\n', txt_filename);

%% --- SUMMARY STATISTICS ---
fprintf('\n=== SUMMARY ===\n');
valid_indices = success_flag == 1;
if sum(valid_indices) > 0
    fprintf('R_occupied statistics (valid points):\n');
    fprintf('  Min: %.2f, Max: %.2f, Mean: %.2f\n', ...
        min(R_occupied_vec(valid_indices)), max(R_occupied_vec(valid_indices)), mean(R_occupied_vec(valid_indices)));
    
    fprintf('E_occupied statistics (valid points):\n');
    fprintf('  Min: %.2f, Max: %.2f, Mean: %.2f\n', ...
        min(E_occupied_vec(valid_indices)), max(E_occupied_vec(valid_indices)), mean(E_occupied_vec(valid_indices)));
    
    fprintf('Gene length coverage:\n');
    fprintf('  Min: %.0f bp, Max: %.0f bp\n', min(L_vec), max(L_vec));
end

fprintf('\nGrid generation complete!\n');
fprintf('Next steps:\n');
fprintf('1. Use this data to build interpolation functions\n');
fprintf('2. Implement conservation equations with gene length distribution\n');
fprintf('3. Solve for self-consistent (R_free, E_free)\n');
fprintf('4. Calculate TCD(L) relationships\n');

%% --- HELPER FUNCTIONS ---

function [R_sol, REH_sol, P_sim] = run_single_gene_simulation(P)
    % Run single gene simulation similar to existing scripts
    global N PAS N_PAS Ef_ss;
    syms Ef real;
    
    % Set up geometry
    L_a = P.L_a;
    N = floor(P.geneLength_bp / L_a);
    PAS = floor(P.PASposition / L_a);
    N_PAS = N - PAS + 1;
    
    % Compute steady states
    [r_E_BeforePas] = compute_steady_states(P, P.EBindingNumber + 1);
    
    % Set up binding values
    kPon_vals = linspace(P.kPon_min, P.kPon_max, PAS);
    RE_vals = sym(zeros(P.EBindingNumber + 1, N));
    
    for e = 1:(P.EBindingNumber + 1)
        for idx = 1:length(kPon_vals)
            RE_vals(e, idx) = subs(r_E_BeforePas(e), {'kPon', 'kPoff'}, {kPon_vals(idx), P.kPoff_const});
        end
        for idx = (PAS+1):N
            RE_vals(e, idx) = subs(r_E_BeforePas(e), {'kPon', 'kPoff'}, {P.kPon_max, P.kPoff_min});
        end
    end
    
    P.RE_val_bind_E = matlabFunction(simplify(sum(sym(1:P.EBindingNumber)' .* RE_vals(2:end, :), 1)), 'Vars', {Ef});
    
    % Solve system
    X_guess = 1e-6 * ones(N + N_PAS, 1);
    options = optimoptions('fsolve', 'Display', 'off', 'FunctionTolerance', 1e-8);
    
    % Two-step solution
    P.FirstRun = true;
    P.is_unphysical = false;
    Ef_ss = 0;
    
    X_base = fsolve(@(xx) ode_dynamics_multipleE(xx, P), X_guess, options);
    if P.is_unphysical
        error('Unphysical result in step 1');
    end
    
    % Update kHon and resolve
    avg_E_bound = P.RE_val_bind_E(Ef_ss);
    P.FirstRun = false;
    P.kHon = P.kHon * avg_E_bound(end);
    X_final = fsolve(@(xx) ode_dynamics_multipleE(xx, P), X_base, options);
    
    R_sol = X_final(1:N);
    REH_sol = X_final((N+1):(N+N_PAS));
    P_sim = P;
end

function E_occupied = calculate_E_occupied(P, R_sol, REH_sol)
    % Calculate total E factors bound to this gene
    global Ef_ss PAS;
    
    % Get binding function values
    avg_E_bound_profile = P.RE_val_bind_E(Ef_ss);
    
    % Calculate E bound to R states (before and after PAS)
    E_bound_R = sum(R_sol .* avg_E_bound_profile');
    
    % Calculate E bound to REH states (after PAS)
    E_bound_REH = sum(REH_sol .* avg_E_bound_profile(PAS:end));
    
    E_occupied = E_bound_R + E_bound_REH;
end
