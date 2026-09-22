%% Sweep2DCad.m
% Full finite-rate CAD maps (x-y): kHd-kEd, kc-kHd and kc-kEd.
% kHd = kHoff/kHon; kEd = kEoff/kEon. On-rates stay at their
% reference values; each unswept parameter stays at its baseline.
% Ranges come from default_parameters(). Ratio envelopes can produce
% off-rates outside their individual source intervals at fixed on-rates.

%% User configuration
save_result = true;
EBindingNumbers = 5;      % Use [1 5] to compare binding capacities
nPoints = 10;             % Log-spaced points per axis (at least 2)
percent_cleavage = 50;

P = default_parameters();

%% Compute, plot and save each pair
for EBindingNumber = EBindingNumbers(:)'
    sweeps = run_cad_parameter_sweeps(P, EBindingNumber, nPoints, percent_cleavage);
    for sweep_index = 1:numel(sweeps)
        data = sweeps(sweep_index);
        plot_cad_parameter_sweep(data);
        if save_result || strcmpi(getenv('CPAD_FORCE_SAVE'), 'true')
            analysis_type = sprintf('sweep_2D_%s_%s_CAD', data.x_name, data.y_name);
            save_analysis_results(analysis_type, data, P, ...
                'ExtraInfo', sprintf('Full_engagedOff%.3g', P.kEoff_engaged));
        end
    end
end
fprintf('All CAD sweeps completed.\n');

%% Local sweep and plotting functions
function sweeps = run_cad_parameter_sweeps(P, EBindingNumber, nPoints, percent_cleavage)
% RUN_CAD_PARAMETER_SWEEPS Compute all three pairwise kHd/kEd/kc CAD maps.
% sweeps = run_cad_parameter_sweeps(P, M, nPoints, percent_cleavage)
% returns a struct array ordered kHd-kEd, kc-kHd, kc-kEd. Matrices have
% rows = y_values, columns = x_values. No plotting or file I/O occurs here.
% P.ranges supplies positive increasing bounds for kHd, kEd and kc.
% On-rates and unswept parameters remain at P's reference values. Engaged
% E detachment is independent of kEd; its full-model default is zero.
% Unreached CAD thresholds remain NaN with finite cleavage diagnostics;
% solver failures additionally carry an error message and NaN diagnostics.

validateattributes(EBindingNumber, {'numeric'}, ...
    {'scalar', 'integer', 'positive', 'finite'}, mfilename, 'EBindingNumber');
validateattributes(nPoints, {'numeric'}, ...
    {'scalar', 'integer', '>=', 2, 'finite'}, mfilename, 'nPoints');
validateattributes(percent_cleavage, {'numeric'}, ...
    {'scalar', 'real', 'finite', '>=', 0, '<=', 100}, mfilename, 'percent_cleavage');
validateattributes(P.kHon, {'numeric'}, {'scalar', 'real', 'finite', 'positive'});
validateattributes(P.kEon, {'numeric'}, {'scalar', 'real', 'finite', 'positive'});
if ~isfield(P, 'kEoff_engaged')
    P.kEoff_engaged = 0;
end

names = {'kHd', 'kEd', 'kc'};
values = cell(size(names));
for axis_index = 1:numel(names)
    name = names{axis_index};
    bounds = P.ranges.(name);
    validateattributes(bounds, {'numeric'}, ...
        {'vector', 'numel', 2, 'real', 'finite', 'positive'}, mfilename, ['ranges.' name]);
    if bounds(1) >= bounds(2)
        error('run_cad_parameter_sweeps:InvalidRange', ...
            'P.ranges.%s must have increasing bounds.', name);
    end
    values{axis_index} = logspace(log10(bounds(1)), log10(bounds(2)), nPoints);
end
base_values = [P.kHoff/P.kHon, P.kEoff/P.kEon, P.kc];

fprintf('Full-model CAD_%g sweeps for M = %d; fixed engaged E off = %.3g s^-1\n', ...
    percent_cleavage, EBindingNumber, P.kEoff_engaged);
fprintf('Running baseline ... ');
baseline = evaluate_cad(P, EBindingNumber, percent_cleavage);

common.model_variant = 'full_kinetics_rapid_EH_disassembly';
common.EBindingNumber = EBindingNumber;
common.percent_cleavage = percent_cleavage;
common.kHon_ref = P.kHon;
common.kEon_ref = P.kEon;
common.kEoff_engaged = P.kEoff_engaged;
common.source_ranges = P.ranges;
common.kHd_base = base_values(1);
common.kEd_base = base_values(2);
common.kc_base = P.kc;
common.CAD_base = baseline.CAD;
common.max_exit_cdf_base = baseline.max_exit_cdf;
common.rhs_residual_base = baseline.rhs_residual;
common.base_error = baseline.error;

pairs = [1 2; 3 1; 3 2];
results = cell(1, size(pairs, 1));
for pair_index = 1:size(pairs, 1)
    x_index = pairs(pair_index, 1);
    y_index = pairs(pair_index, 2);
    data = common;
    data.x_name = names{x_index};
    data.y_name = names{y_index};
    data.x_values = values{x_index};
    data.y_values = values{y_index};
    data.x_base = base_values(x_index);
    data.y_base = base_values(y_index);
    data.CAD_matrix = NaN(nPoints);
    data.max_exit_cdf_matrix = NaN(nPoints);
    data.rhs_residual_matrix = NaN(nPoints);
    data.error_matrix = repmat({''}, nPoints, nPoints);
    fprintf('\nStarting %s vs %s (%d runs)\n', data.x_name, data.y_name, nPoints^2);
    timer_id = tic;
    for row = 1:nPoints
        for col = 1:nPoints
            % Start each cell from P so previous sweeps cannot leak rates.
            P_run = set_sweep_parameter(P, data.x_name, data.x_values(col));
            P_run = set_sweep_parameter(P_run, data.y_name, data.y_values(row));
            fprintf('  [%d/%d] %s=%.3g, %s=%.3g ... ', ...
                (row-1)*nPoints+col, nPoints^2, data.x_name, data.x_values(col), ...
                data.y_name, data.y_values(row));
            result = evaluate_cad(P_run, EBindingNumber, percent_cleavage);
            data.CAD_matrix(row, col) = result.CAD;
            data.max_exit_cdf_matrix(row, col) = result.max_exit_cdf;
            data.rhs_residual_matrix(row, col) = result.rhs_residual;
            data.error_matrix{row, col} = result.error;
        end
    end
    data.elapsed_seconds = toc(timer_id);
    fprintf('Sweep completed in %.1f s\n', data.elapsed_seconds);
    results{pair_index} = data;
end
sweeps = [results{:}];
end

function P = set_sweep_parameter(P, name, value)
switch name
    case 'kHd'
        P.kHoff = value * P.kHon;
    case 'kEd'
        P.kEoff = value * P.kEon;
    case 'kc'
        P.kc = value;
end
end

function result = evaluate_cad(P, EBindingNumber, percent_cleavage)
result = struct('CAD', NaN, 'max_exit_cdf', NaN, 'rhs_residual', NaN, 'error', '');
try
    [R, RHE, P_sim, details] = run_full_termination_simulation(P, EBindingNumber);
    [~, ~, cad, diagnostics] = calculate_full_pas_cleavage_profile(R, RHE, P_sim, ...
        'PercentCleavage', percent_cleavage);
    result.CAD = cad;
    result.max_exit_cdf = diagnostics.max_exit_cdf;
    result.rhs_residual = details.rhs_max_abs;
    fprintf('CAD = %.0f bp; within-window cleavage = %.2f%%\n', ...
        cad, 100*result.max_exit_cdf);
catch ME
    result.error = sprintf('%s: %s', ME.identifier, ME.message);
    fprintf('FAILED: %s\n', result.error);
end
end

function fig = plot_cad_parameter_sweep(data)
% PLOT_CAD_PARAMETER_SWEEP Plot one result from run_cad_parameter_sweeps.
% Display is capped at 1600 bp like pyCPA/viz/plot_all.py; raw CAD is kept.
% Gray regions have no measured CAD. Markers distinguish unreached
% thresholds (open circles) from solver failures (crosses).

display_cap = 2000;
fig = figure('Name', sprintf('CAD: %s vs %s, M = %d', ...
    data.x_name, data.y_name, data.EBindingNumber), 'Position', [100 100 900 700]);
ax = axes('Parent', fig);
hold(ax, 'on');
[X, Y] = meshgrid(data.x_values, data.y_values);
cad = data.CAD_matrix;
finite_cad = cad(isfinite(cad));
if ~isempty(finite_cad)
    displayed_cad = cad;
    displayed_cad(displayed_cad > display_cap) = display_cap;
    contourf(ax, X, Y, displayed_cad, linspace(0, display_cap, 21), ...
        'LineStyle', 'none');
    % Use uncapped values for target contours; never infer CAD from NaNs.
    for level = [400 800]
        if min(finite_cad) < level && max(finite_cad) > level
            [C, h] = contour(ax, X, Y, cad, [level level], 'r--', 'LineWidth', 2.5);
            clabel(C, h, 'FontSize', 10, 'Color', 'r');
        end
    end
end
set(ax, 'XScale', 'log', 'YScale', 'log', 'FontSize', 12, 'Color', [0.85 0.85 0.85]);
caxis(ax, [0 display_cap]);
if exist('turbo', 'builtin') || exist('turbo', 'file')
    colormap(ax, turbo(256));
else
    colormap(ax, jet(256));
end
cb = colorbar(ax);
cb.Label.String = sprintf('CAD_{%g} (bp; display capped at %g)', ...
    data.percent_cleavage, display_cap);

censored = isnan(cad) & isfinite(data.max_exit_cdf_matrix);
failed = ~cellfun(@isempty, data.error_matrix);
legend_handles = gobjects(0);
if any(censored(:))
    legend_handles(end+1) = plot(ax, X(censored), Y(censored), 'ko', ...
        'MarkerSize', 7, 'LineWidth', 1.2, 'DisplayName', 'Threshold not reached');
end
if any(failed(:))
    legend_handles(end+1) = plot(ax, X(failed), Y(failed), 'kx', ...
        'MarkerSize', 8, 'LineWidth', 1.5, 'DisplayName', 'Simulation failed');
end
plot(ax, data.x_base, data.y_base, 'p', ...
    'MarkerSize', 18, 'MarkerFaceColor', 'w', 'MarkerEdgeColor', 'k', ...
    'LineWidth', 1.5, 'HandleVisibility', 'off');
if isfinite(data.CAD_base)
    base_label = sprintf('Base parameters CAD_{%g} = %.0f bp', ...
        data.percent_cleavage, data.CAD_base);
elseif isempty(data.base_error)
    base_label = sprintf('Base parameters CAD_{%g} = not reached', data.percent_cleavage);
else
    base_label = sprintf('Base parameters CAD_{%g} = simulation failed', data.percent_cleavage);
end
% A text-only legend entry reports baseline CAD without repeating the star.
base_legend = plot(ax, NaN, NaN, 'LineStyle', 'none', 'Marker', 'none', ...
    'DisplayName', base_label);
legend(ax, [base_legend legend_handles], 'Location', 'northwest', 'Interpreter', 'tex');
xlabel(ax, parameter_label(data.x_name), 'FontSize', 14);
ylabel(ax, parameter_label(data.y_name), 'FontSize', 14);
title(ax, sprintf('Full model CAD_{%g} (M = %d, engaged E off = %.3g s^{-1})', ...
    data.percent_cleavage, data.EBindingNumber, data.kEoff_engaged));
hold(ax, 'off');
end

function label = parameter_label(name)
switch name
    case 'kHd'
        label = 'k_{Hd} = k_{Hoff} / k_{Hon}';
    case 'kEd'
        label = 'k_{Ed} = k_{Eoff} / k_{Eon} (molecules)';
    case 'kc'
        label = 'k_c (s^{-1})';
end
end
