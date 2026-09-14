% Compare engaged-E off-rates in the full finite-rate R/RHE model across E capacity.
saveData = strcmpi(getenv('CPAD_FORCE_SAVE'), 'true');
fprintf('=== E Binding Number: Full Finite-Rate Model ===\n\n');

% --- PARAMETERS AND SWEEP CONFIGURATION ---
P = default_parameters();
EBindingNumber_values = 1:6; % Use 1:7 or [1,5]; both maximum p and e change.
engaged_E_off_rates = [0, 0.05, 0.5]; % s^-1; 0 is the protected-E limit
cutoff_threshold = 0.5;
percent_cleavage = 100*cutoff_threshold;

% Each column corresponds to one full-model engaged-E off-rate.
num_capacities = numel(EBindingNumber_values);
num_series = numel(engaged_E_off_rates);
series_labels = cell(1, num_series);
for j = 1:num_series
    series_labels{j} = sprintf('Full: engaged E off = %.3g s^{-1}', engaged_E_off_rates(j));
    if engaged_E_off_rates(j) == 0
        series_labels{j} = [series_labels{j} ' (protected-E limit)'];
    end
end
cutoff_positions = NaN(num_capacities, num_series);
max_exit_cdf = NaN(num_capacities, num_series);
Ef_ss_values = NaN(num_capacities, num_series);
Pol_free_values = NaN(num_capacities, num_series);
rhs_residuals = NaN(num_capacities, num_series);
error_messages = repmat({''}, num_capacities, num_series);

fprintf('Ordinary E off-rate: %.3g s^-1\n', P.kEoff);
fprintf('Full-model engaged-E off-rates: %s s^-1\n', mat2str(engaged_E_off_rates));
fprintf('Sweeping E Binding Number: %s\n', mat2str(EBindingNumber_values));

% --- COMPUTE EVERY CURVE WITH THE CURRENT PARAMETERS ---
for i = 1:num_capacities
    EBindingNumber = EBindingNumber_values(i);
    for j = 1:num_series
        fprintf('M = %d, %s... ', EBindingNumber, series_labels{j});
        try
            P_case = P; % Each run starts with the unchanged base kHon.
            P_case.kEoff_engaged = engaged_E_off_rates(j);
            [R_sol, RHE_sol, P_sim, full_details] = ...
                run_full_termination_simulation(P_case, EBindingNumber);
            Pol_free = P_sim.Pol_free_ss;
            rhs_max_abs = full_details.rhs_max_abs;
            % Same flux observable for all curves: include the first cleavage
            % bin and return NaN if the requested threshold is not reached.
            [~, ~, CAD, cleavage_diagnostics] = calculate_full_pas_cleavage_profile( ...
                R_sol, RHE_sol, P_sim, 'PercentCleavage', percent_cleavage);
            cutoff_positions(i, j) = CAD;
            max_exit_cdf(i, j) = cleavage_diagnostics.max_exit_cdf;
            Ef_ss_values(i, j) = P_sim.Ef_ss;
            Pol_free_values(i, j) = Pol_free;
            rhs_residuals(i, j) = rhs_max_abs;
            fprintf('CAD%.0f = %.0f bp; within-window cleavage = %.2f%%\n', ...
                percent_cleavage, CAD, 100*max_exit_cdf(i, j));
        catch ME
            error_messages{i, j} = ME.message;
            fprintf('ERROR: %s\n', ME.message);
        end
    end
end

% --- PLOT COMPARISON ---
fig = figure('Color', 'w', 'Position', [100, 100, 1200, 780]);
ax = axes('Parent', fig);
hold(ax, 'on');
full_colors = [21, 111, 138; 207, 110, 23; 143, 67, 139]/255;
if numel(engaged_E_off_rates) > size(full_colors, 1)
    full_colors = [full_colors; lines(numel(engaged_E_off_rates)-size(full_colors, 1))];
end
for j = 1:num_series
    plot(ax, EBindingNumber_values, cutoff_positions(:, j), '-o', ...
        'Color', full_colors(j, :), 'MarkerFaceColor', full_colors(j, :), ...
        'LineWidth', 2.5, 'MarkerSize', 8, 'DisplayName', series_labels{j});
end
hold(ax, 'off');
xlabel(ax, 'Maximum E binding and phosphorylation capacity (M)', 'FontSize', 15);
ylabel(ax, sprintf('CAD_{%.0f} (bp)', percent_cleavage), 'FontSize', 15);
title(ax, 'Rapid disassembly after engaged E detachment', 'FontSize', 20, 'FontWeight', 'normal');
legend(ax, 'Location', 'northeast', 'Box', 'off', 'FontSize', 12, 'Interpreter', 'tex');
set(ax, 'FontSize', 13, 'XTick', EBindingNumber_values, 'Box', 'off', ...
    'TickDir', 'out', 'LineWidth', 1, 'GridAlpha', 0.12);
xlim(ax, [min(EBindingNumber_values)-0.25, max(EBindingNumber_values)+0.25]);
finite_cad = cutoff_positions(isfinite(cutoff_positions));
if isempty(finite_cad)
    ylim(ax, [0, 1]);
else
    ylim(ax, [0, max(1, 1.08*max(finite_cad))]);
end
grid(ax, 'on');

% One row per (engaged-E off-rate, M), all using the full finite-rate model.
model_names = repmat({'full'}, num_capacities, num_series);
capacity_grid = repmat(EBindingNumber_values(:), 1, num_series);
rate_grid = repmat(engaged_E_off_rates, num_capacities, 1);
summary_table = table(model_names(:), capacity_grid(:), rate_grid(:), ...
    cutoff_positions(:), max_exit_cdf(:), Ef_ss_values(:), Pol_free_values(:), ...
    rhs_residuals(:), error_messages(:), 'VariableNames', ...
    {'Model', 'M', 'kEoff_engaged', 'CAD_bp', 'max_exit_cdf', 'E_free', ...
     'Pol_free', 'rhs_max_abs', 'Error'});
fprintf('\n=== CAD%.0f Comparison (bp) ===\n', percent_cleavage);
disp(summary_table(:, {'Model', 'M', 'kEoff_engaged', 'CAD_bp'}));

if saveData
    data = struct();
    data.EBindingNumber_values = EBindingNumber_values;
    data.engaged_E_off_rates = engaged_E_off_rates;
    data.percent_cleavage = percent_cleavage;
    data.series_labels = series_labels;
    data.cutoff_positions = cutoff_positions;
    data.max_exit_cdf = max_exit_cdf;
    data.Ef_ss_values = Ef_ss_values;
    data.Pol_free_values = Pol_free_values;
    data.rhs_residuals = rhs_residuals;
    data.error_messages = error_messages;
    data.model_variant = 'full_kinetics_rapid_EH_disassembly';
    output_dir = cpad_analysis_output_dir('Full_EBindingNumber_vs_CAD', ...
        fullfile(fileparts(mfilename('fullpath')), 'SecondVersionResults'));
    stem = sprintf('Engaged_E_off_comparison_%s', datestr(now, 'yyyymmdd_HHMMSSFFF'));
    save(fullfile(output_dir, [stem '.mat']), 'data', 'P', 'summary_table');
    writetable(summary_table, fullfile(output_dir, [stem '.csv']));
    exportgraphics(fig, fullfile(output_dir, [stem '.png']), 'Resolution', 200);
    fprintf('Comparison results saved to: %s\n', output_dir);
end
