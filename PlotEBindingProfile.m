% Support Figures: Average E Binding and Ser2P Profiles
% For each EBindingNumber, generates:
%   - Average E binding profile (solid line)
%   - Ser2P profile (dashed line, same color)
% on a single combined figure
saveData = false;

% Define output directory
outputDir = 'Results/SupportFigures/';
if saveData && ~exist(outputDir, 'dir')
    mkdir(outputDir);
end

% --- BASE PARAMETERS ---
P = default_parameters();


% Analysis Scenarios
BindingNumbers = [1, 3, 5];
Colors = {'r', 'g', 'b'};

fprintf('Starting simulations for Support Figures...\n');

%% Combined plot: Average E (solid) + Ser2P (dashed)
fig = figure('Name', 'Average E and Ser2P Profiles', 'Position', [100, 100, 900, 600]);
hold on;

for idx = 1 : length(BindingNumbers)
    nb = BindingNumbers(idx);
    fprintf('Simulating EBindingNumber = %d...\n', nb);

    P.EBindingNumber = nb;

    % --- Run simulation (uses file cache for symbolic step) ---
    [R_sol, REH_sol, P_out, r_E_BeforePas, r_P] = run_termination_simulation(P, nb);

    % --- Average E binding and Ser2P profiles in a single pass ---
    [avg_E_bound, avg_Ser2P] = P_out.RE_val_bind_E(P_out.Ef_ss);

    % --- Plot coordinates ---
    x_coords = ((1 - P_out.PAS):(P_out.N - P_out.PAS)) * P.L_a / 1000;  % Position relative to PAS in kb
    
    % --- Save raw data ---
    if saveData
        dataOutDir = fullfile('Results', 'Ser2P_Eaverage_Profile');
        if ~exist(dataOutDir, 'dir')
            mkdir(dataOutDir);
        end
        filename = fullfile(dataOutDir, sprintf('ProfileData_N%d.txt', nb));
        T = table(x_coords(:), avg_E_bound(:), avg_Ser2P(:), ...
            'VariableNames', {'PositionRelPAS_kb', 'Avg_E_bound', 'Avg_Ser2P'});
        writetable(T, filename, 'Delimiter', '\t');
    end

    % --- Plot ---
    col = Colors{idx};

    plot(x_coords, avg_E_bound, '-',  'Color', col, 'LineWidth', 2, ...
        'DisplayName', sprintf('Avg E (N=%d)', nb));
    plot(x_coords, avg_Ser2P,  '--', 'Color', col, 'LineWidth', 2, ...
        'DisplayName', sprintf('Ser2P (N=%d)', nb));
end

% Finalize plot
xlabel('Position relative to PAS (kb)', 'FontSize', 12);
ylabel('Average Bound Factors', 'FontSize', 12);
title('Average E-Factor Binding (solid) and Ser2P (dashed) Profiles', 'FontSize', 14);
legend('Location', 'northwest', 'FontSize', 10);
xline(0, 'k--', 'PAS', 'LineWidth', 1.5);
grid on;

if saveData
    saveas(fig, fullfile(outputDir, 'Average_E_and_Ser2P_Comparison.png'));
end

fprintf('Support figures generated.\n');