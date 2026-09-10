% Full finite-rate R/RHE model: maximum E capacity versus CAD50.
% Includes engaged-E detachment followed by rapid EH disassembly.
saveData = strcmpi(getenv('CPAD_FORCE_SAVE'), 'true');
fprintf('=== Full R/RHE E Binding Number Sweep ===\n\n');

% --- BASE PARAMETERS ---
P = default_parameters();
P.kEoff_engaged = 0.05; % s^-1; trial value, 0 recovers the protected-E limit
fprintf('Engaged-E off-rate: %.3g s^-1; ordinary E off-rate: %.3g s^-1\n', ...
    P.kEoff_engaged, P.kEoff);

% --- SWEEP CONFIGURATION ---
EBindingNumber_values = [1, 2, 3, 4, 5, 6];
% Use 1:7 to include M=7; both maximum p and maximum e change with M.
cutoff_threshold = 0.5;  % 50% termination threshold

% Pre-allocate results
cutoff_positions = NaN(size(EBindingNumber_values));
max_exit_cdf = NaN(size(EBindingNumber_values));
Ef_ss_values = NaN(size(EBindingNumber_values));
Pol_free_values = NaN(size(EBindingNumber_values));
rhs_residuals = NaN(size(EBindingNumber_values));

% --- SWEEP LOOP ---
fprintf('Sweeping E Binding Number: %s\n', mat2str(EBindingNumber_values));

for i = 1:length(EBindingNumber_values)
    EBindingNumber = EBindingNumber_values(i);
    fprintf('Running EBindingNumber = %d... ', EBindingNumber);
    
    try
        % Solve all finite-rate microstates with self-consistent free pools.
        [R_sol, RHE_sol, P_sim, full_details] = run_full_termination_simulation(P, EBindingNumber);
        
        % Find distance where the requested fraction of polymerases has cleaved
        [~, ~, CAD, cleavage_diagnostics] = calculate_full_pas_cleavage_profile( ...
            R_sol, RHE_sol, P_sim, 'PercentCleavage', cutoff_threshold * 100);
        cutoff_positions(i) = CAD;
        max_exit_cdf(i) = cleavage_diagnostics.max_exit_cdf;
        Ef_ss_values(i) = P_sim.Ef_ss;
        Pol_free_values(i) = P_sim.Pol_free_ss;
        rhs_residuals(i) = full_details.rhs_max_abs;
        
        fprintf('CAD50 = %.0f bp; within-window cleavage = %.2f%%\n', ...
            cutoff_positions(i), 100*max_exit_cdf(i));
        
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

xlabel('Maximum E binding and phosphorylation capacity (M)', 'FontSize', 12, 'FontWeight', 'bold');
ylabel('CAD_{50} (bp)', 'FontSize', 12, 'FontWeight', 'bold');
title(sprintf('Full R/RHE model: engaged E off-rate = %.3g s^{-1}', P.kEoff_engaged));
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
    data.max_exit_cdf = max_exit_cdf;
    data.Ef_ss_values = Ef_ss_values;
    data.Pol_free_values = Pol_free_values;
    data.rhs_residuals = rhs_residuals;
    data.model_variant = 'full_kinetics_rapid_EH_disassembly';
    output_dir = cpad_analysis_output_dir('Full_EBindingNumber_vs_CAD', ...
        fullfile(fileparts(mfilename('fullpath')), 'SecondVersionResults'));
    stem = sprintf('Full_EBinding_vs_CAD_engagedOff%.3g_%s', ...
        P.kEoff_engaged, datestr(now, 'yyyymmdd_HHMMSSFFF'));
    save(fullfile(output_dir, [stem '.mat']), 'data', 'P');
    summary_table = table(EBindingNumber_values(:), cutoff_positions(:), max_exit_cdf(:), ...
        Ef_ss_values(:), Pol_free_values(:), rhs_residuals(:), ...
        'VariableNames', {'M', 'CAD50_bp', 'max_exit_cdf', 'E_free', 'Pol_free', 'rhs_max_abs'});
    writetable(summary_table, fullfile(output_dir, [stem '.csv']));
    saveas(gcf, fullfile(output_dir, [stem '.png']));
    fprintf('Full-model results saved to: %s\n', output_dir);
end
