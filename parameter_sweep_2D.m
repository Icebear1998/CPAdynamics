% --- Model parameters (base values) ---
save_result = false;
P = struct();
P.k_in = 2;
P.kEon = 0.0000025;
P.kEoff = 0.1;
P.k_e = 65/100;
P.k_e2 = 30/100;
P.E_total = 100000;
P.L_total = 100000;
P.Pol_total = 70000;
P.kHon = 2;
P.kHoff = 1;
P.kc = 0.1;
P.kPon_min = 0.01;
P.kPon_slope = 0.005;
P.kPoff = 1;

L_a = 100;
geneLength_bp = 25000;
PASposition = 20000;
EBindingNumber = 1;

% Add required fields for PAS usage calculation
P.L_a = L_a;
P.geneLength_bp = geneLength_bp;
P.PASposition = PASposition;

% --- Global variables (set internally by run_termination_simulation) ---
global N PAS N_PAS Ef_ss;
Ef_ss = 0;

% --- Parameter Sweep Setup ---
param_pairs = {{'E_total', 'kHoff'}};

for pair_idx = 1:length(param_pairs)
    param1 = param_pairs{pair_idx}{1};
    param2 = param_pairs{pair_idx}{2};
    
    param1_values = 30000:10000:100000;
    param2_values = logspace(-2,1,10);

    cutoff_matrix = zeros(length(param2_values), length(param1_values));

    % --- 2D Parameter Sweep Loop ---
    fprintf('Starting 2D sweep for %s vs %s...\n', param1, param2);
    for i = 1:length(param2_values)
        for j = 1:length(param1_values)
            fprintf('  Running %s = %.2g, %s = %.2g (%d,%d)\n', param2, param2_values(i), param1, param1_values(j), i, j);

            P_run = P;
            P_run.(param1) = param1_values(j);
            P_run.(param2) = param2_values(i);

            try
                [R_sol, REH_sol, P_sim] = run_termination_simulation(P_run, EBindingNumber);
                [~, ~, cutoff_matrix(i,j)] = calculate_pas_cleavage_profile(R_sol, REH_sol, P_sim, 'PercentCleavage', 50);
            catch
                cutoff_matrix(i,j) = NaN;
            end
        end
        fprintf('Completed row %d/%d for %s = %.2g\n', i, length(param2_values), param2, param2_values(i));
    end
    fprintf('2D sweep complete.\n');

    % --- Plotting ---
    figure('Position', [100, 100, 800, 600]);
    hold on;
    colors = lines(size(cutoff_matrix, 1));
    for i = 1:size(cutoff_matrix, 1)
        plot(param1_values, cutoff_matrix(i, :), 'o-', 'LineWidth', 2, ...
             'Color', colors(i, :), 'DisplayName', sprintf('%s = %.2g', param2, param2_values(i)));
    end
    xlabel([strrep(param1, '_', '\_'), ' Value'], 'FontSize', 12);
    ylabel('Position at which 50% termination (bp)', 'FontSize', 12);
    title(['75% Termination Position vs ', strrep(param1, '_', '\_'), ' by ', strrep(param2, '_', '\_')], 'FontSize', 14);
    grid on; legend('show', 'Location', 'best'); set(gca, 'FontSize', 10); box on;
    hold off;
    
    % --- SAVE RESULTS ---
    % Prepare data structure for saving
    if save_result
        data.EBindingNumber = EBindingNumber;
        data.param1 = param1;
        data.param2 = param2;
        data.param1_values = param1_values;
        data.param2_values = param2_values;
        data.cutoff_matrix = cutoff_matrix;

        % Save results using the utility function
        save_analysis_results('parameter_sweep_2D', data, P);
    end

end
