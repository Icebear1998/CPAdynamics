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

% --- Common Parameters ---
P.L_a     = 100;
P.k_in    = 2;
P.kEon    = 0.00025;
P.kEoff   = 1;
P.k_e     = 65 / P.L_a;
P.k_e2    = 30 / P.L_a;
P.E_total = 100000;
P.L_total = 100000;
P.Pol_total = 70000;
P.kc      = 0.1;
P.kHon    = 0.1;
P.kHoff   = 0.01;

P.kPon_min   = 0.01;
P.kPon_slope = 0.005;
P.kPoff      = 1;

P.geneLength_bp = 25000;
P.PASposition   = 20000;

% Analysis Scenarios
BindingNumbers = [9];
Colors = {'r', 'g', 'b', 'm', 'c'};

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
    kPon_vals = P_out.kPon_min + P_out.kPon_slope * (0 : N - 1);

    % --- Average E binding profile ---
    avg_E_bound = P_out.RE_val_bind_E(Ef_ss);

    % --- Ser2P profile ---
    % Build P_vals by substituting kPon and kPoff into r_P (matches CPA_multipleE_main.m)
    P_vals = sym(zeros(nb + 1, N));
    for e = 1:nb+1
        for i = 1:N
            P_vals(e, i) = subs(r_P(e), {'kPon', 'kPoff'}, {kPon_vals(i), P_out.kPoff});
        end
    end

    avg_Ser2P = zeros(1, N);
    for i = 1:N
        total_P_bound = 0;
        total_P = 0;
        for e = 1:nb+1
            P_e = double(R_sol(i) * double(subs(P_vals(e, i), 'Ef', Ef_ss)));
            total_P_bound = total_P_bound + (e - 1) * P_e;
            total_P = total_P + P_e;
        end
        if total_P > 0
            avg_Ser2P(i) = total_P_bound / total_P;
        else
            avg_Ser2P(i) = 0;
        end
    end

    % --- Plot ---
    x_coords = ((1 - PAS):(N - PAS)) * P.L_a / 1000;  % Position relative to PAS in kb
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