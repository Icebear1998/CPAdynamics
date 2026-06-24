% Support Figures: Average E Binding and Ser2P Profiles
% For each EBindingNumber, generates:
%   - Average E binding profile (solid line)
%   - Ser2P profile (dashed line, same color)
% on a single combined figure

clear;
close all;
global N PAS N_PAS Ef_ss;

% Define output directory
outputDir = 'SecondVersionResults/SupportFigures/';
if ~exist(outputDir, 'dir')
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

    % Retrieve geometry and Ef_ss set by run_termination_simulation

    % --- Average E binding and Ser2P profiles in a single pass ---
    [avg_E_bound, avg_Ser2P] = P_out.RE_val_bind_E(Ef_ss);

        % --- Plot coordinates ---
    x_coords = ((1 - PAS):(N - PAS)) * P.L_a / 1000;  % Position relative to PAS in kb
    
    % --- Save raw data ---
    dataOutDir = fullfile('SecondVersionResults', 'Ser2P_Eaverage_Profile');
    if ~exist(dataOutDir, 'dir')
        mkdir(dataOutDir);
    end
    
    filename = fullfile(dataOutDir, sprintf('ProfileData_N%d.txt', nb));
    T = table(x_coords(:), avg_E_bound(:), avg_Ser2P(:), ...
        'VariableNames', {'PositionRelPAS_kb', 'Avg_E_bound', 'Avg_Ser2P'});
    writetable(T, filename, 'Delimiter', '\t');

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

saveas(fig, fullfile(outputDir, 'Average_E_and_Ser2P_Comparison.png'));

fprintf('Support figures generated in %s\n', outputDir);