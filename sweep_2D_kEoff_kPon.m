% --- CONFIGURATION ---
% Fixed E binding number for analysis
EBindingNumber = 1;

% Define sweep ranges
kEoff_values = [1, 5, 10, 20];      % kEoff sweep range
kPon_slope_values = [0.005, 0.01, 0.02, 0.03];  % kPon_slope sweep range

% Define output directory
outputDir = 'SecondVersionResults/2DSweep/';
if ~exist(outputDir, 'dir')
    mkdir(outputDir);
end

% --- COMMON PARAMETERS ---
P.L_a = 100;
P.k_in = 2;
P.kEon = 0.00025;
P.kEoff = 10;  % Default (will be swept)
P.k_e = 65 / P.L_a;
P.k_e2 = 30 / P.L_a;
P.E_total = 100000;
P.L_total = 100000;
P.Pol_total = 70000;
P.kc = 0.1;
P.kHon = 0.1;  % Initial kHon value
P.kHoff = 0.01;

P.kPon_min = 0.01;
P.kPon_slope = 0.02;  % Default (will be swept)
P.kPoff = 1;

P.geneLength_bp = 25000;
P.PASposition = 20000;

fprintf('Starting 2D sweep: kEoff vs kPon_slope for E Binding Number = %d...\n', EBindingNumber);
fprintf('kEoff range: [%.1f, %.1f], %d points\n', min(kEoff_values), max(kEoff_values), length(kEoff_values));
fprintf('kPon_slope range: [%.3f, %.3f], %d points\n', min(kPon_slope_values), max(kPon_slope_values), length(kPon_slope_values));

% Pre-allocate result matrix
nKEoff = length(kEoff_values);
nKPon = length(kPon_slope_values);
avg_E_at_PAS = zeros(nKPon, nKEoff);  % Rows: kPon_slope, Cols: kEoff

%% Nested loop over parameter space
totalIterations = nKEoff * nKPon;
currentIter = 0;

for i = 1:nKEoff
    for j = 1:nKPon
        currentIter = currentIter + 1;
        
        % Update parameters
        P_local = P;
        P_local.kEoff = kEoff_values(i);
        P_local.kPon_slope = kPon_slope_values(j);
        
        fprintf('  [%d/%d] Simulating kEoff=%.1f, kPon_slope=%.3f...\n', ...
            currentIter, totalIterations, P_local.kEoff, P_local.kPon_slope);
        
        try
            % Run simulation using existing helper function
            [R_sol, REH_sol, P_local, r_E_BeforePas] = run_termination_simulation(P_local, EBindingNumber);
            
            % Calculate Ef_ss from solution
            Ef_ss_local = calculate_Efree_from_solution(R_sol, REH_sol, P_local);
            
            % Get average E binding profile
            avg_E_bound = P_local.RE_val_bind_E(Ef_ss_local);
            
            % Extract value at PAS
            PAS = floor(P_local.PASposition / P_local.L_a);
            avg_E_at_PAS(j, i) = avg_E_bound(PAS);
            
        catch ME
            warning('Simulation failed for kEoff=%.1f, kPon_slope=%.3f: %s', ...
                P_local.kEoff, P_local.kPon_slope, ME.message);
            avg_E_at_PAS(j, i) = NaN;
        end
    end
end

%% Create Contour Plot
fig1 = figure('Name', '2D Sweep: Average E at PAS', 'Position', [100, 100, 900, 700]);

% Create meshgrid for contour
[X, Y] = meshgrid(kEoff_values, kPon_slope_values);

% Filled contour background
nLevels = 20;
contourf(X, Y, avg_E_at_PAS, nLevels, 'LineStyle', 'none');
colorbar;
colormap(jet);
hold on;

% Overlay labeled isolines
[C, h] = contour(X, Y, avg_E_at_PAS, 8, 'k', 'LineWidth', 1.5);
clabel(C, h, 'FontSize', 11, 'Color', 'k', 'FontWeight', 'bold');
hold off;

% Labels and formatting
xlabel('k_{Eoff}', 'FontSize', 14, 'FontWeight', 'bold');
ylabel('k_{Pon} slope', 'FontSize', 14, 'FontWeight', 'bold');
title(sprintf('Average E Binding at PAS (Contour)\n(E Binding Number = %d)', EBindingNumber), ...
    'FontSize', 16, 'FontWeight', 'bold');
set(gca, 'FontSize', 12);
grid on;

% Save figure
saveas(fig1, fullfile(outputDir, sprintf('AvgE_at_PAS_Contour_N%d.png', EBindingNumber)));
saveas(fig1, fullfile(outputDir, sprintf('AvgE_at_PAS_Contour_N%d.fig', EBindingNumber)));

%% Create Heatmap
fig2 = figure('Name', '2D Sweep: Average E at PAS (Heatmap)', 'Position', [150, 150, 900, 700]);

imagesc(kEoff_values, kPon_slope_values, avg_E_at_PAS);
set(gca, 'YDir', 'normal');
colorbar;
colormap(jet);

xlabel('k_{Eoff}', 'FontSize', 14, 'FontWeight', 'bold');
ylabel('k_{Pon} slope', 'FontSize', 14, 'FontWeight', 'bold');
title(sprintf('Average E Binding at PAS (Heatmap)\n(E Binding Number = %d)', EBindingNumber), ...
    'FontSize', 16, 'FontWeight', 'bold');
set(gca, 'FontSize', 12);

% Save figure
saveas(fig2, fullfile(outputDir, sprintf('AvgE_at_PAS_Heatmap_N%d.png', EBindingNumber)));
saveas(fig2, fullfile(outputDir, sprintf('AvgE_at_PAS_Heatmap_N%d.fig', EBindingNumber)));

%% Save data to file
dataFile = fullfile(outputDir, sprintf('2D_sweep_data_N%d.mat', EBindingNumber));
save(dataFile, 'kEoff_values', 'kPon_slope_values', 'avg_E_at_PAS', 'EBindingNumber', 'P');

% Also save as text file for easy inspection
txtFile = fullfile(outputDir, sprintf('2D_sweep_data_N%d.txt', EBindingNumber));
fid = fopen(txtFile, 'w');
fprintf(fid, '2D Parameter Sweep Results\n');
fprintf(fid, 'E Binding Number: %d\n\n', EBindingNumber);
fprintf(fid, 'Average E Binding at PAS\n');
fprintf(fid, 'Rows: kPon_slope values, Columns: kEoff values\n\n');

% Header row
fprintf(fid, 'kPon_slope\\kEoff');
for i = 1:nKEoff
    fprintf(fid, '\t%.2f', kEoff_values(i));
end
fprintf(fid, '\n');

% Data rows
for j = 1:nKPon
    fprintf(fid, '%.4f', kPon_slope_values(j));
    for i = 1:nKEoff
        fprintf(fid, '\t%.4f', avg_E_at_PAS(j, i));
    end
    fprintf(fid, '\n');
end
fclose(fid);

fprintf('\n2D sweep completed! (EBindingNumber = %d)\n', EBindingNumber);
fprintf('Results saved to: %s\n', outputDir);
fprintf('  - Contour plot: AvgE_at_PAS_Contour_N%d.png\n', EBindingNumber);
fprintf('  - Heatmap:      AvgE_at_PAS_Heatmap_N%d.png\n', EBindingNumber);
fprintf('  - Data (.mat): 2D_sweep_data_N%d.mat\n', EBindingNumber);
fprintf('  - Data (.txt): 2D_sweep_data_N%d.txt\n', EBindingNumber);

%% Helper function to calculate Ef_ss from solution
function Ef_ss = calculate_Efree_from_solution(R_sol, REH_sol, P)
    % Calculate free E based on solution and conservation
    N = floor(P.geneLength_bp / P.L_a);
    PAS = floor(P.PASposition / P.L_a);
    
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
