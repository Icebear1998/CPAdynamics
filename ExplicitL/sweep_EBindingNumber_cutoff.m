% SWEEP_EBINDINGNUMBER_CUTOFF.m (Explicit L-factor version)
% Plot CAD (50% cleavage distance) vs. E Binding Number
fprintf('=== E Binding Number Sweep Analysis (Explicit L) ===\n\n');

% --- BASE PARAMETERS ---
P.L_a = 100;
P.k_in    = 2;
P.kEon    = 0.0001;
P.kEoff   = 1;
P.k_e     = 65/P.L_a;
P.k_e2    = 30/P.L_a;
P.E_total = 100000;
P.L_total = 100000;
P.Pol_total = 70000;

% Reaction 1: Tethered E-hexamer encounter (per E-factor)
P.kHon  = 7;
P.kHoff = 0.03;

% Reaction 2: L-factor diffusional recruitment (per E-factor per L-molecule)
P.kLon  = 2.2e-6;
P.kLoff = 0.001;

% Cleavage (only from fully assembled REHL state)
P.kc = 0.1;

P.kPon_min = 0.01;
P.kPon_slope = 0.05;
P.kPoff = 1;
P.geneLength_bp = 25000;
P.PASposition = 20000;


% --- SWEEP CONFIGURATION ---
EBindingNumber_values = [1, 2, 3, 4, 5];
cutoff_threshold = 0.5;  % 50% termination threshold

% Pre-allocate results
cutoff_positions = zeros(size(EBindingNumber_values));

% --- SWEEP LOOP ---
fprintf('Sweeping E Binding Number: %s\n', mat2str(EBindingNumber_values));

for i = 1:length(EBindingNumber_values)
    EBindingNumber = EBindingNumber_values(i);
    fprintf('Running EBindingNumber = %d... ', EBindingNumber);
    
    try
        % Run simulation (explicit L version returns REHL_sol)
        [R_sol, REH_sol, REHL_sol, P_sim] = run_termination_simulation(P, EBindingNumber);
        
        % Calculate termination profile (CDF) using REHL for cleavage
        [exit_cdf, distances_bp] = calculate_pas_usage_profile(R_sol, REH_sol, REHL_sol, P_sim);
        
        % Find position where cutoff threshold is reached
        if max(exit_cdf) >= cutoff_threshold
            cutoff_positions(i) = interp1(exit_cdf, distances_bp, cutoff_threshold, 'linear');
        else
            cutoff_positions(i) = NaN;  % Threshold not reached
            fprintf('(CDF max = %.1f%%) ', max(exit_cdf)*100);
        end
        
        fprintf('Cutoff position = %.0f bp\n', cutoff_positions(i));
        
    catch ME
        fprintf('ERROR: %s\n', ME.message);
        cutoff_positions(i) = NaN;
    end
end

% --- PLOT RESULTS ---
figure('Position', [100, 100, 800, 600]);
plot(EBindingNumber_values, cutoff_positions, 'o', 'LineWidth', 2.5, 'MarkerSize', 10, ...
     'Color', [0, 0.4470, 0.7410], 'MarkerFaceColor', [0, 0.4470, 0.7410]);

% Add data labels
for i = 1:length(EBindingNumber_values)
    if ~isnan(cutoff_positions(i))
        text(EBindingNumber_values(i), cutoff_positions(i) + 50, ...
             sprintf('%.0f', cutoff_positions(i)), ...
             'FontSize', 11, 'HorizontalAlignment', 'center');
    end
end

xlabel('Number of maximum CPA factor binding sites', 'FontSize', 12, 'FontWeight', 'bold');
ylabel('CAD (bp)', 'FontSize', 12, 'FontWeight', 'bold');
xlim([0.5 7.5]);
ax = gca;
ax.XTick = unique(round(ax.XTick));
xtickformat('%d');
set(gca, 'FontSize', 11);
grid on;
box on;

fprintf('\n=== Analysis Complete ===\n');
fprintf('Results:\n');
for i = 1:length(EBindingNumber_values)
    fprintf('  EBindingNumber = %d: %.0f bp\n', EBindingNumber_values(i), cutoff_positions(i));
end
