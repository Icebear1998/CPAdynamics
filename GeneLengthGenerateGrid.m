% GeneLengthGenerateGrid.m — full finite-rate model, fixed shared free pools
% Generates lookup data for R_occupied and E_occupied as functions of
% (R_free, E_free, L) for the gene length analysis
%
% This script implements Step 1 of the gene length analysis
% - Grid Generation: Create 3D grid of (R_free, E_free, L)
% - Parallel Computation: Run single-gene model for each grid point
% - Data Saving: Store results for later interpolation

% Start parallel pool if not already running
if isempty(gcp('nocreate')); parpool; end

fprintf('=== Gene Length Grid Generation ===\n');
fprintf('Generating lookup data for R_occupied and E_occupied...\n\n');

%% --- PARAMETER RANGES ---

% Base totals from typical simulations
P_base = default_parameters();
R_total_base = P_base.Pol_total;
E_total_base = P_base.E_total;

% Cover all physically possible free pools, including zero. At fixed E,
% occupancy is exactly linear in R_free, so a coarse R grid is sufficient.
R_free_min = 0;
R_free_max = R_total_base;
R_free_points = 5;
E_free_min = 0;
E_free_max = E_total_base;
E_free_points = 5;

% TSS-to-PAS length range: based on human gene distribution (including introns)
% From the histogram: spans ~10^2 to ~10^6 bp, with most genes 10^3 to 10^5
% NOTE: L now represents TSS-to-PAS distance, not total gene length
L_min = 2500;      % 2.5 kb (lower end of main distribution)
L_max = 200000;    % 200 kb (covers most genes, excludes extreme outliers)
L_points = 10;     % Resolution for TSS-to-PAS length

% Fixed after-PAS length for all genes
after_PAS_length = 5000;  % 5 kb constant after-PAS region

%% --- BASE PARAMETERS ---
% k_in must be rescaled for the local model. The original value (e.g., 2)
% was for a global model where one gene represented the whole system.
% Here, we scale it to represent a single gene's promoter strength.
k_in_global = P_base.k_in;
num_active_genes = 10000; % Estimated number of active genes in the system
P_base.k_in = k_in_global / num_active_genes;
P_base.EBindingNumber = 5;  % Use standard value

% Create parameter grids
R_free_values = linspace(R_free_min, R_free_max, R_free_points);
E_free_values = linspace(E_free_min, E_free_max, E_free_points);
L_values = logspace(log10(L_min), log10(L_max), L_points); % Log spacing for gene lengths
L_values([1 end]) = [L_min L_max];

fprintf('Parameter Ranges:\n');
fprintf('  R_free: %.0f to %.0f (%d points)\n', R_free_min, R_free_max, R_free_points);
fprintf('  E_free: %.0f to %.0f (%d points)\n', E_free_min, E_free_max, E_free_points);
fprintf('  L (TSS-to-PAS): %.0f to %.0f bp (%d points, log-spaced)\n', L_min, L_max, L_points);
fprintf('  After-PAS length: %.0f bp (fixed for all genes)\n', after_PAS_length);
fprintf('  kPon: linear increase with slope = %.4f\n', P_base.kPon_slope);
fprintf('  Total grid points: %d\n\n', R_free_points * E_free_points * L_points);

%% --- GRID GENERATION ---
fprintf('Generating parameter grid...\n');

% Create all combinations of parameters
[R_grid, E_grid, L_grid] = ndgrid(R_free_values, E_free_values, L_values);

% Flatten to vectors for parallel processing
R_free_vec = R_grid(:);
E_free_vec = E_grid(:);
L_vec = L_grid(:);
n_points = length(R_free_vec);

% Pre-allocate results
R_occupied_vec = NaN(n_points, 1);
E_occupied_vec = NaN(n_points, 1);
rhs_residual_vec = NaN(n_points, 1);
error_messages = repmat({''}, n_points, 1);
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
        
        % NEW: L_i now represents TSS-to-PAS distance
        P_i.PASposition = L_i;     % PAS position = TSS-to-PAS distance
        P_i.geneLength_bp = L_i + after_PAS_length;  % Total gene length = TSS-to-PAS + after-PAS
        
        % Shared pools are prescribed here; do not solve per-gene totals.
        [~, ~, ~, details] = run_full_termination_simulation(P_i, P_i.EBindingNumber, ...
            'FreePools', [R_free_i, E_free_i]);
        R_occupied_i = details.Pol_bound;
        E_occupied_i = details.E_bound; % Exact sum of E over all R/RHE microstates
        rhs_residual_vec(i) = details.rhs_max_abs;

        % Store results
        R_occupied_vec(i) = R_occupied_i;
        E_occupied_vec(i) = E_occupied_i;
        success_flag(i) = 1;
    catch ME
        success_flag(i) = 0;
        error_messages{i} = ME.message;
    end
    
    % Progress indication (every 100 points)
    if mod(i, 100) == 0
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
results.metadata.model_variant = 'full_kinetics_rapid_EH_disassembly';
results.metadata.pool_mode = 'fixed_free_pools';
results.metadata.grid_layout = 'ndgrid';
results.metadata.computation_time_minutes = computation_time/60;
results.metadata.success_rate = success_rate;
results.metadata.description = 'Lookup data for gene length analysis: R_occupied and E_occupied as functions of (R_free, E_free, L) where L = TSS-to-PAS distance';

% Parameter information
results.parameters.R_free_range = [R_free_min, R_free_max];
results.parameters.E_free_range = [E_free_min, E_free_max];
results.parameters.L_range = [L_min, L_max];
results.parameters.R_free_points = R_free_points;
results.parameters.E_free_points = E_free_points;
results.parameters.L_points = L_points;
results.parameters.after_PAS_length = after_PAS_length;
results.parameters.base_parameters = P_base;

% Add system-wide parameters for consistency
results.parameters.R_total_base = R_total_base;
results.parameters.E_total_base = E_total_base;
results.parameters.num_active_genes = num_active_genes;

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
results.data.rhs_residual_vec = rhs_residual_vec;
results.data.error_messages = error_messages;

%% --- SAVE RESULTS ---
fprintf('Saving results...\n');

% Create output directory
output_dir = cpad_analysis_output_dir('GeneLengthAnalysis', 'SecondVersionResults');

% Save MATLAB data file
mat_filename = fullfile(output_dir, sprintf('full_gene_length_grid_data_%d.mat', P_base.EBindingNumber));
save(mat_filename, 'results', '-v7.3');  % Use v7.3 for large files

% Save text file with summary and data
txt_filename = fullfile(output_dir, sprintf('full_gene_length_grid_data_%d.txt', P_base.EBindingNumber));
fid = fopen(txt_filename, 'w');

% Header
fprintf(fid, '%% Gene Length Analysis - Full finite-rate grid data, fixed free pools\n');
fprintf(fid, '%% kEoff_engaged = %g s^-1\n', P_base.kEoff_engaged);
fprintf(fid, '%% Generated on: %s\n', results.metadata.creation_date);
fprintf(fid, '%% Computation time: %.1f minutes\n', results.metadata.computation_time_minutes);
fprintf(fid, '%% Success rate: %.1f%%\n', results.metadata.success_rate);
fprintf(fid, '%% \n');

% Parameters
fprintf(fid, '%% Parameters:\n');
fprintf(fid, '%% R_free range: %.0f to %.0f (%d points)\n', R_free_min, R_free_max, R_free_points);
fprintf(fid, '%% E_free range: %.0f to %.0f (%d points)\n', E_free_min, E_free_max, E_free_points);
fprintf(fid, '%% L (TSS-to-PAS) range: %.0f to %.0f bp (%d points)\n', L_min, L_max, L_points);
fprintf(fid, '%% After-PAS length: %.0f bp (fixed for all genes)\n', after_PAS_length);
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
fprintf(fid, '%% 1: R_free, 2: E_free, 3: L (TSS-to-PAS), 4: R_occupied, 5: E_occupied, 6: success_flag\n');
fprintf(fid, '%% Note: Total gene length = L + %.0f bp (fixed after-PAS region)\n', after_PAS_length);
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
    
    fprintf('TSS-to-PAS distance coverage:\n');
    fprintf('  Min: %.0f bp, Max: %.0f bp\n', min(L_vec), max(L_vec));
    fprintf('  Total gene lengths: %.0f to %.0f bp (including %.0f bp after-PAS)\n', ...
        min(L_vec) + after_PAS_length, max(L_vec) + after_PAS_length, after_PAS_length);
end

fprintf('\nGrid generation complete!\n');
fprintf('Next steps:\n');
fprintf('1. Use this data to build interpolation functions\n');
fprintf('2. Implement conservation equations with gene length distribution\n');
fprintf('3. Solve for self-consistent (R_free, E_free)\n');
fprintf('4. Calculate TCD(L) relationships\n');
