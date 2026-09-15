% Support Figures: Average E Binding and Ser2P Profiles
% For each EBindingNumber, generates:
%   - Average E binding profile (solid line)
%   - Ser2P profile (dashed line, same color)
% on a single combined figure from TSS to PAS (inclusive).
% The full gene is still simulated and exported for consistent free pools.
saveData = strcmpi(getenv('CPAD_FORCE_SAVE'), 'true');

% --- BASE PARAMETERS ---
P = default_parameters();


% Analysis Scenarios
BindingNumbers = [1, 5, 10];
Colors = {'r', 'b', 'g'};

fprintf('Starting simulations for Support Figures...\n');

%% Combined plot: Average E (solid) + Ser2P (dashed)
fig = figure('Name', 'Average E and Ser2P Profiles', 'Position', [100, 100, 900, 600]);
hold on;

for idx = 1 : length(BindingNumbers)
    nb = BindingNumbers(idx);
    fprintf('Simulating EBindingNumber = %d...\n', nb);

    % --- Full finite-rate microstate distributions ---
    [~, ~, P_out, full_details] = run_full_termination_simulation(P, nb);

    % --- Average E binding and Ser2P profiles in a single pass ---
    avg_E_bound = full_details.avg_E_bound;
    avg_Ser2P = full_details.avg_Ser2P;

    % --- Plot coordinates ---
    x_coords = ((1 - P_out.PAS):(P_out.N - P_out.PAS)) * P.L_a / 1000;  % Position relative to PAS in kb
    
    % --- Save raw data ---
    if saveData
        dataOutDir = cpad_analysis_output_dir('Ser2P_Eaverage_Profile', 'Results');
        filename = fullfile(dataOutDir, sprintf('ProfileData_N%d.txt', nb));
        T = table(x_coords(:), avg_E_bound(:), avg_Ser2P(:), ...
            'VariableNames', {'PositionRelPAS_kb', 'Avg_E_bound', 'Avg_Ser2P'});
        writetable(T, filename, 'Delimiter', '\t');
    end

    % --- Plot ---
    col = Colors{idx};
    plot_nodes = 1:P_out.PAS;
    plot_coords = plot_nodes * P.L_a / 1000; % Distance from TSS, kb

    plot(plot_coords, avg_E_bound(plot_nodes), '-',  'Color', col, 'LineWidth', 2, ...
        'DisplayName', sprintf('Avg E (N=%d)', nb));
    plot(plot_coords, avg_Ser2P(plot_nodes),  '--', 'Color', col, 'LineWidth', 2, ...
        'DisplayName', sprintf('Ser2P (N=%d)', nb));
end

% Finalize plot
xlabel('Distance from TSS (kb)', 'FontSize', 12);
ylabel('Average Bound Factors', 'FontSize', 12);
title('Full finite-rate E binding (solid) and Ser2P (dashed)', 'FontSize', 14);
legend('Location', 'northwest', 'FontSize', 10);
pas_kb = P_out.PAS * P.L_a / 1000;
xlim([0, pas_kb]);
xticks(linspace(0, pas_kb, 5));
xline(pas_kb, 'k--', 'PAS', 'LineWidth', 1.5, ...
    'LabelHorizontalAlignment', 'left', 'HandleVisibility', 'off');
grid on;

if saveData
    outputDir = cpad_analysis_output_dir('SupportFigures', 'Results');
    saveas(fig, fullfile(outputDir, 'Average_E_and_Ser2P_Comparison.png'));
end

fprintf('Support figures generated.\n');
