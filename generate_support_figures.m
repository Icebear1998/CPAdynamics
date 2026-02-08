% Multi-Threaded Simulation for Support Figures
% This script runs simulations with different EBindingNumbers and generates:
% 1. Average E-binding profile comparison
% 2. Pol II state distribution heatmap (for N=5 or similar)
% 3. Fraction of "Competent" Pol II

clear;
close all;
global N PAS N_PAS Ef_ss;

% Define output directory
outputDir = 'SecondVersionResults/SupportFigures/';
if ~exist(outputDir, 'dir')
    mkdir(outputDir);
end

% --- Common Parameters (Match CPA_multipleE_main.m) ---
P.L_a = 100;
P.k_in = 2;
P.kEon = 0.00025;
P.kEoff = 10;
P.k_e = 65 / P.L_a;
P.k_e2 = 30 / P.L_a;
P.E_total = 70000;
P.L_total = 100000;
P.Pol_total = 70000;
P.kc = 0.05;

% Parameters for kHon calculation (from main: kHon approx 0.2, then renormalized)
% Note: The main script recalculates kHon based on E_bound at PAS. We will replicate this.
Initial_kHon = 0.2;
P.kHoff = 0.0125;

kPon_min = 0.01;
kPon_slope = 0.005;
kPoff = 1;

geneLength_bp = 25000;
PASposition = 20000;
N = floor(geneLength_bp / P.L_a);
PAS = floor(PASposition / P.L_a);
N_PAS = N - PAS + 1;
kPon_vals = kPon_min + kPon_slope * (0 : N - 1);

% Analysis Scenarios
BindingNumbers = [ 1, 3, 5];
Colors = {'r', 'g', 'b'};
avg_E_profiles = {};

fprintf('Starting simulations for Support Figures...\n');

%% 1. Loop over Binding Numbers for Profile Comparison
figure('Name', 'Average E Binding Profile', 'Position', [100, 100, 800, 600]);
hold on;

for idx = 1 : length(BindingNumbers)
    nb = BindingNumbers(idx);
    fprintf('Simulating EBindingNumber = %d...\n', nb);

    P.EBindingNumber = nb;
    P.kHon = Initial_kHon;
    
    % Reset global Ef (if needed, though usually handled in dynamics)
    % Ef_ss = 0;

    % Step A: Pre-compute Steady States (Symbolic)
    [r_E_BeforePas, r_P] = compute_steady_states(P, nb + 1);

    % Prepare RE_vals (Symbolic array)
    RE_vals = sym(zeros(nb + 1, N));
    
    % Construct RE_vals matching main script logic
    for e = 1:nb+1
        % Substitute kPoff once (scalar)
        r_E_e_fixed = subs(r_E_BeforePas(e), 'kPoff', kPoff);
        
        % Substitute kPon vector element-wise
        for i = 1:N
            RE_vals(e, i) = subs(r_E_e_fixed, 'kPon', kPon_vals(i));
        end
    end
    
    % Create handle for total E bound (function of Ef)
    % Sums weighted average: sum( (1:nb) * Probability )
    P.RE_val_bind_E = matlabFunction(simplify(sum(sym(1:nb)' .* RE_vals(2:end, :), 1)), 'Vars', {'Ef'});
    
    % Step B: Solve ODE
    P.FirstRun = true;
    P.is_unphysical = false;
    X0 = 1e-6 * ones(N + N_PAS, 1);
    options = optimoptions('fsolve', 'Display', 'off', 'FunctionTolerance', 1e-8);
    X = fsolve(@(xx) ode_dynamics_multipleE(xx, P), X0, options);
    
    % Recalculate kHon constraint (match main script)
    avg_E_bound_temp = P.RE_val_bind_E(Ef_ss);
    P.kHon = P.kHon * avg_E_bound_temp(PAS);
    P.FirstRun = false;
    X = fsolve(@(xx) ode_dynamics_multipleE(xx, P), X, options);
    
    % Step C: Extract Results
    avg_E_bound = P.RE_val_bind_E(Ef_ss);
    
    % Store for plotting (PAS-centered coordinates)
    x_coords = ((1-PAS):(N-PAS)) * P.L_a / 1000;  % Position relative to PAS in kb
    plot(x_coords, avg_E_bound, 'Color', Colors{idx}, 'LineWidth', 2, ...
        'DisplayName', sprintf('Max Sites = %d', nb));
    
    % --- HEATMAP & Ser2P GENERATION (Only for N=5 default case) ---
    if nb == 5
        fprintf('Generating Heatmap and Ser2P profile for N=5...\n');
        
        % 1. Calculate Pol II State Distribution (Heatmap Data)
        HeatmapData = zeros(nb+1, N);
        
        % Calculate actual State Probabilities at Ef_ss
        RE_vals_numeric = zeros(nb+1, N);
        for e = 1:nb+1
            for i = 1:N
                val = double(subs(RE_vals(e, i), 'Ef', Ef_ss));
                RE_vals_numeric(e, i) = val;
            end
        end
        
        HeatmapData = RE_vals_numeric;
        
        % Plot Heatmap
        hFig = figure('Name', 'Pol II State Heatmap', 'Position', [150, 150, 1000, 400]);
        x_coords_heatmap = ((1-PAS):(N-PAS)) * P.L_a / 1000;  % PAS-centered
        imagesc(x_coords_heatmap, 0:nb, HeatmapData);
        colorbar;
        colormap(jet);
        xlabel('Position relative to PAS (kb)');
        ylabel('Number of Bound E Factors');
        title(sprintf('Pol II State Distribution (Max Sites = %d)', nb));
        xline(0, 'w--', 'PAS', 'LineWidth', 2);
        saveas(hFig, fullfile(outputDir, 'PolII_State_Heatmap_N5.png'));
        
        % 2. Calculate Ser2P Profile (avg_P)
        avg_P = zeros(1, N);
        
        for i = 1:N
            num_sum = 0; 
            den_sum = 0;
            for e = 1:nb+1
                % Get P_val for this state and pos
                r_P_e = r_P(e);
                % Subs kPon, kPoff, Ef
                val = double(subs(r_P_e, {'kPon', 'kPoff', 'Ef'}, {kPon_vals(i), kPoff, Ef_ss}));
                
                num_sum = num_sum + (e-1) * val;
                den_sum = den_sum + val;
            end
            if den_sum > 0
                avg_P(i) = num_sum / den_sum;
            end
        end
        
        % Plot Ser2P
        hFig2 = figure('Name', 'Ser2P Profile', 'Position', [150, 150, 800, 400]);
        x_coords_ser2p = ((1-PAS):(N-PAS)) * P.L_a / 1000;  % PAS-centered
        plot(x_coords_ser2p, avg_P, 'm-', 'LineWidth', 2.5);
        xlabel('Position relative to PAS (kb)');
        ylabel('Average Ser2P (Available Sites)');
        title('Predicted Ser2P Phosphorylation Profile (N=5)');
        xline(0, 'k--', 'PAS');
        grid on;
        saveas(hFig2, fullfile(outputDir, 'Ser2P_Profile_N5.png'));
    end
end

% Finalize Profile Plot
xlabel('Position relative to PAS (kb)');
ylabel('Average Number of Bound E Factors');
title('Average E-Factor Binding Profile');
xline(0, 'k--', 'PAS');
legend('Location', 'northwest');
grid on;

saveas(gcf, fullfile(outputDir, 'Average_E_Binding_Comparison.png'));

fprintf('Support figures generated in %s\n', outputDir);