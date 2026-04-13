% PASUAGEVSINTERPAS DISTANCE.m
% Calculate and plot Proximal PAS Usage vs. Inter-PAS Distance
% Uses flux-based CDF approach for accurate termination probability

clear; clc;
fprintf('=== Proximal PAS Usage vs Inter-PAS Distance ===\n\n');

% --- BASE PARAMETERS ---
P.L_a = 100;
P.k_in    = 2;
P.k_e     = 65/P.L_a;
P.k_e2    = 30/P.L_a;
P.E_total = 100000;
P.L_total = 100000;
P.Pol_total = 70000;

% We are not confident here
P.kEon = 0.0000025;
P.kEoff = 0.1;
P.kHon = 2; 
P.kHoff = 1; 
P.kc = 0.1; 

P.kPon_min = 0.01; % at TSS
P.kPon_slope = 0.005; % determined how fast Sep2P increasing from TSS
P.kPoff = 1;
P.geneLength_bp = 25000;
P.PASposition = 20000;

% --- CONFIGURATION ---
EBindingNumber = 2;
save_result = false;

% --- RUN SIMULATION ---
fprintf('Running simulation for EBindingNumber = %d...\n', EBindingNumber);

try
    [R_sol, REH_sol, P_sim] = run_termination_simulation(P, EBindingNumber);
    fprintf('Simulation completed successfully.\n');
catch ME
    error('Simulation failed: %s', ME.message);
end

% --- CALCULATE TERMINATION PROFILE (CDF) ---
fprintf('Calculating termination profile...\n');

% Use standardized function to calculate PAS usage profile
[exit_cdf, distances_bp] = calculate_pas_usage_profile(R_sol, REH_sol, P_sim);

% Prepend zero for interpolation at distance = 0
distances_for_interp = [0; distances_bp(:)];
cdf_for_interp = [0; exit_cdf(:)];

% --- CALCULATE PROXIMAL PAS USAGE ---
fprintf('Calculating proximal PAS usage...\n');

% Define inter-PAS distances
inter_pas_distances_bp = 0:50:2000;

% Calculate proximal PAS usage for each inter-PAS distance
proximal_usage_prob = zeros(size(inter_pas_distances_bp));

for i = 1:length(inter_pas_distances_bp)
    dist_bp = inter_pas_distances_bp(i);
    
    % CDF at this distance tells us fraction that terminated by this point (at proximal PAS)
    % So proximal usage = CDF at that distance
    proximal_usage_prob(i) = interp1(distances_for_interp, cdf_for_interp, dist_bp, 'linear', 'extrap');
end

% --- PLOT RESULTS ---
fprintf('Generating plot...\n');

figure('Position', [100, 100, 800, 600]);
hold on;

plot(inter_pas_distances_bp, proximal_usage_prob * 100, '-', 'LineWidth', 2.5, ...
     'MarkerSize', 6, 'Color', 'r', 'DisplayName', 'Proximal PAS Usage');

% Add vertical line at 300 bp (median tandem distance)
line([300 300], [0 100], 'Color', 'k', 'LineStyle', '--', 'LineWidth', 1.5, ...
     'DisplayName', 'Typical inter-PAS distance (300 bp)');

% Customize plot
xlabel('Inter-PAS Distance (bp)', 'FontSize', 12, 'FontWeight', 'bold');
ylabel('Proximal Site Usage (%)', 'FontSize', 12, 'FontWeight', 'bold');
%title(sprintf('Proximal PAS Usage vs Inter-PAS Distance (EBindingNumber=%d)', EBindingNumber), ...
      %'FontSize', 14, 'FontWeight', 'bold');
grid on;
legend('show', 'Location', 'best', 'FontSize', 11);
set(gca, 'FontSize', 11);
box on;
ylim([0 60]);

hold off;

% --- SAVE RESULTS ---
if save_result
    fprintf('Saving results...\n');
    
    % Prepare data structure
    data.EBindingNumber = EBindingNumber;
    data.inter_pas_distances_bp = inter_pas_distances_bp;
    data.proximal_usage_prob = proximal_usage_prob;
    data.exit_cdf = exit_cdf;
    data.distances_bp = distances_bp;
    data.R_sol = R_sol;
    data.REH_sol = REH_sol;
    
    % Save using utility function
    save_analysis_results('PASUsagevsInterPASDistance', data, P);
end

fprintf('\n=== Analysis Complete ===\n');
fprintf('Results at key distances:\n');
fprintf('  100 bp: %.1f%% proximal usage\n', interp1(inter_pas_distances_bp, proximal_usage_prob*100, 100));
fprintf('  300 bp: %.1f%% proximal usage\n', interp1(inter_pas_distances_bp, proximal_usage_prob*100, 300));
fprintf('  500 bp: %.1f%% proximal usage\n', interp1(inter_pas_distances_bp, proximal_usage_prob*100, 500));
fprintf('  1000 bp: %.1f%% proximal usage\n', interp1(inter_pas_distances_bp, proximal_usage_prob*100, 1000));

