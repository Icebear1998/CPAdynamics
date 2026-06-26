% SWEEP_EBINDINGNUMBER_CUTOFF.m
% Plot position at which 50% cutoff vs. E Binding Number
% Simple script to analyze how E binding number affects termination distance
saveData = false;
fprintf('=== E Binding Number Sweep Analysis ===\n\n');

% --- BASE PARAMETERS ---
P = default_parameters();


% --- SWEEP CONFIGURATION ---
EBindingNumber_values = [1, 2, 3, 4, 5, 6];
cutoff_threshold = 0.5;  % 50% termination threshold

% Pre-allocate results
cutoff_positions = zeros(size(EBindingNumber_values));

% --- SWEEP LOOP ---
fprintf('Sweeping E Binding Number: %s\n', mat2str(EBindingNumber_values));

for i = 1:length(EBindingNumber_values)
    EBindingNumber = EBindingNumber_values(i);
    fprintf('Running EBindingNumber = %d... ', EBindingNumber);
    
    try
        % Run simulation
        [R_sol, REH_sol, P_sim] = run_termination_simulation(P, EBindingNumber);
        
        % Find distance where cutoff_threshold% of polymerases have cleaved
        [~, ~, CAD] = calculate_pas_cleavage_profile(R_sol, REH_sol, P_sim, 'PercentCleavage', cutoff_threshold * 100);
        cutoff_positions(i) = CAD;
        
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

xlabel('Number of maxium CPA factor binding', 'FontSize', 12, 'FontWeight', 'bold');
ylabel('CAD (bp)', 'FontSize', 12, 'FontWeight', 'bold');
%xlim([0.5 7.5]);
% Set x-axis ticks to only integer values
ax = gca;
ax.XTick = unique(round(ax.XTick));
% Format x-axis tick labels to display as integers
xtickformat('%d');
set(gca, 'FontSize', 11);
grid on;
box on;

fprintf('\n=== Analysis Complete ===\n');
fprintf('Results:\n');
for i = 1:length(EBindingNumber_values)
    fprintf('  EBindingNumber = %d: %.0f bp\n', EBindingNumber_values(i), cutoff_positions(i));
end

if saveData
    data.EBindingNumber_values = EBindingNumber_values;
    data.cutoff_positions = cutoff_positions;
    save_analysis_results('EBindingNumber_vs_CAD', data, P);
end
