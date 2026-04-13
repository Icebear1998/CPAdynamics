% PASUsageAnalysis.m

% --- CONFIGURATION ---
save_result = false;
% 1. Define the ranges for the two sweep parameters
kHoff_values = logspace(log10(0.001), log10(0.1), 6);
E_total_values = 30000:40000:200000;

% 2. Define the fixed, biologically relevant distance for the metric
fixed_distance_bp = 300;

% --- BASE PARAMETERS (P struct) ---
% (Parameter definitions are unchanged)
P.L_a = 100; P.k_in = 2; P.kEon = 0.0000025; P.kEoff = 0.1;
P.k_e = 65/100; P.k_e2 = 30/100; P.E_total = 100000;
P.L_total = 100000; P.Pol_total = 70000; P.kHon = 2;
P.kHoff = 1; P.kc = 0.1; P.kPon_min = 0.01; P.kPon_slope = 0.005;
P.kPoff = 1;
P.geneLength_bp = 25000; P.PASposition = 20000; P.EBindingNumber = 1;

% --- Pre-allocate the results matrix ---
results_matrix = zeros(length(E_total_values), length(kHoff_values));

% --- 2D PARAMETER SWEEP LOOP ---
fprintf('Starting 2D sweep over kHoff and E_total...\n');

for i = 1:length(E_total_values)
    % Create a temporary row vector for the results of the inner loop
    temp_row = zeros(1, length(kHoff_values));
    
    for j = 1:length(kHoff_values)
        fprintf('  Running E_total = %d, kHoff = %.4g (%d,%d)...\n', E_total_values(i), kHoff_values(j), i, j);
        
        P_run = P;
        P_run.E_total = E_total_values(i);
        P_run.kHoff = kHoff_values(j);

        [R_sol, REH_sol, P_sim] = run_termination_simulation(P_run, P_run.EBindingNumber);
        
        if any(isnan(R_sol))
            temp_row(j) = NaN;
            continue;
        end

        % Use standardized function to calculate PAS usage profile
        try
            [exit_cdf, distances_bp] = calculate_pas_usage_profile(R_sol, REH_sol, P_sim);
            
            % Interpolate to get usage at fixed distance
            % Prepend zero for interpolation at distance = 0
            distances_for_interp = [0; distances_bp(:)];
            cdf_for_interp = [0; exit_cdf(:)];
            
            usage_at_dist = interp1(distances_for_interp, cdf_for_interp, fixed_distance_bp, 'linear', 'extrap');
            temp_row(j) = usage_at_dist * 100;
        catch
            temp_row(j) = NaN;
        end
    end
    % Assign the entire computed row to the results matrix
    results_matrix(i, :) = temp_row;
    fprintf('Finished row %d/%d for E_total = %d\n', i, length(E_total_values), E_total_values(i));
end
disp('All simulations complete.');

% % --- PLOT THE CONTOUR MAP ---
figure('Position', [100, 100, 900, 700]);
[X, Y] = meshgrid(kHoff_values, E_total_values);
contourf(X, Y, results_matrix, 10, 'LineColor', 'k');
set(gca, 'XScale', 'log');
xlabel('PAS Strength (kHoff)', 'FontSize', 12, 'FontWeight', 'bold');
ylabel('Global Factor (E_{total})', 'FontSize', 12, 'FontWeight', 'bold');
title(sprintf('Phase Diagram of Cumulative Polymerase Exit (at %d bp) of Enum = %d', fixed_distance_bp, P.EBindingNumber), 'FontSize', 14, 'FontWeight', 'bold');
c = colorbar;
c.Label.String = 'Cumulative Exit (%)';
c.Label.FontSize = 12;
colormap('parula');
set(gca, 'FontSize', 10);
box on;

% --- Plotting ---
figure('Position', [100, 100, 800, 600]);
hold on;
colors = lines(size(results_matrix, 1));
for i = 1:size(results_matrix, 1)
    plot(kHoff_values, results_matrix(i, :), 'o-', 'LineWidth', 2, ...
         'Color', colors(i, :), 'DisplayName', sprintf('%s = %.2g', 'E_{total}', E_total_values(i)));
end
set(gca, 'XScale', 'log');
xlabel('kHoff', 'FontSize', 12, 'FontWeight', 'bold');
ylabel('Cumulative Exit (%)', 'FontSize', 12, 'FontWeight', 'bold');
title(sprintf('Phase Diagram of Cumulative Polymerase Exit (at %d bp) of Enum = %d', fixed_distance_bp, P.EBindingNumber), 'FontSize', 14, 'FontWeight', 'bold');
grid on; legend('show', 'Location', 'best'); set(gca, 'FontSize', 10); box on;
hold off;

if save_result
    % --- SAVE RESULTS ---
    % Prepare data structure for saving
    data.results_matrix = results_matrix;
    data.x_values = kHoff_values;
    data.y_values = E_total_values;
    data.x_label = 'kHoff';
    data.y_label = 'E_total';
    data.metric_distance = fixed_distance_bp;

    % Save results using the utility function
    save_analysis_results('PASUsageAnalysis', data, P);
end
