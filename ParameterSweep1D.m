% Parameter sweep setup
save_result = strcmpi(getenv('CPAD_FORCE_SAVE'), 'true');

% Model parameters (base values)
P_default = default_parameters();
L_a = P_default.L_a;

% --- PARAMETERS TO SWEEP ---
% Available: every parameter with an entry in P.ranges (default_parameters.m):
%            'k_e', 'k_e2', 'E_total', 'Pol_total', 'kc', 'kEon', 'kEoff',
%            'kEoff_engaged', 'kHoff', 'kHon', 'kPon_slope', 'kPoff',
%            plus 'k_in' (exploratory range; none in default_parameters).
% All parameters use the full finite-rate solver.
param_list = {'kc'};

EBindingNumber = 1;
n_points = 8;   % grid points across each P.ranges interval

% Iterate over each parameter to sweep
for param_idx = 1:length(param_list)
    % Reset to default parameter set
    P = P_default;
        
    param_to_sweep = param_list{param_idx};
    default_value = P.(param_to_sweep);
        
    % Sweep range comes from P_default.ranges (default_parameters.m);
    % the default value is added so it appears on the curve.
    if strcmp(param_to_sweep, 'k_in')
        % default_parameters supplies no range for k_in; keep exploratory range
        base_range = logspace(log10(0.2), log10(10), 5);
    elseif isfield(P_default.ranges, param_to_sweep)
        base_range = sweep_values_from_ranges(P_default, param_to_sweep, n_points);
    else
        error('Invalid parameter selected');
    end
    param_values = sort(unique([base_range, default_value]));

    % Initialize arrays
    cutoff_values = NaN(1, length(param_values));
    max_exit_cdf = NaN(size(cutoff_values));
    rhs_residuals = NaN(size(cutoff_values));
    error_messages = repmat({''}, size(cutoff_values));

    % --- PARAMETER SWEEP ---
    fprintf('Starting sweep for %s...\n', param_to_sweep);
    for k = 1:length(param_values)
        fprintf('  Running %s = %.4g (%d/%d)...\n', param_to_sweep, param_values(k), k, length(param_values));
        
        % Create local copy of parameters and update for this iteration
        P_run = P_default;
        P_run.(param_to_sweep) = param_values(k);

        try
            [R_sol, REH_sol, P_sim, full_details] = run_full_termination_simulation(P_run, EBindingNumber);
        catch ME
            fprintf('    Error in run_full_termination_simulation: %s\n', ME.message);
            error_messages{k} = ME.message;
            continue;
        end
        
        % --- CALCULATE TERMINATION PROFILE ---
        try
            [~, ~, cutoff_values(k), diagnostics] = calculate_full_pas_cleavage_profile( ...
                R_sol, REH_sol, P_sim, 'PercentCleavage', 50);
            max_exit_cdf(k) = diagnostics.max_exit_cdf;
            rhs_residuals(k) = full_details.rhs_max_abs;
            cutoff_values(k) = round(cutoff_values(k));
        catch ME
            cutoff_values(k) = NaN;
            error_messages{k} = ME.message;
        end
    end
    fprintf('Sweep for %s complete.\n\n', param_to_sweep);

    % --- PLOT RESULTS ---
    figure('Position', [100, 100, 800, 600]);
    hold on;
    
    % Use log scale for rate constant parameters
    log_scale_params = {'k_in', 'kEon', 'kEoff', 'kEoff_engaged', 'kHon', 'kHoff', 'kPon_slope'};
    if ismember(param_to_sweep, log_scale_params)
        set(gca, 'XScale', 'log');
    end
    
    % Plot sweep results
    plot(param_values, cutoff_values, 'o-', 'LineWidth', 2, 'MarkerSize', 8, ...
         'Color', [0, 0.4470, 0.7410], 'DisplayName', '50% Termination');
    
    % Highlight default parameter value
    default_cutoff = cutoff_values(param_values == default_value);
    plot(default_value, default_cutoff, 'ro', 'MarkerSize', 10, 'LineWidth', 2, ...
         'DisplayName', 'Default Value');
    
    % Customize plot appearance
    xlabel([strrep(param_to_sweep, '_', '\_'), ' Value'], 'FontSize', 12);
    ylabel('CAD (bp)', 'FontSize', 12);
    title(['Cleavage and termination distance (CAD) vs ', strrep(param_to_sweep, '_', '\_'), ...
           ' (EBindingNumber=', num2str(EBindingNumber), ')'], 'FontSize', 14, 'FontWeight', 'bold');
    grid on;
    legend('show', 'Location', 'best', 'FontSize', 10);
    set(gca, 'FontSize', 10);
    box on;
    
    % --- SAVE RESULTS ---
    if save_result
        data = struct();
        data.model_variant = 'full_kinetics_rapid_EH_disassembly';
        data.max_exit_cdf = max_exit_cdf;
        data.rhs_residuals = rhs_residuals;
        data.error_messages = error_messages;
        data.EBindingNumber = EBindingNumber;
        data.sweep_param = param_to_sweep;
        data.param_values = param_values;
        data.cutoff_values = cutoff_values;
        save_analysis_results('parameter_sweep_1D', data, P_default, ...
            'Folder', 'ParameterSweep1D', ...
            'ParamSet', sprintf('M%d_%s', EBindingNumber, param_to_sweep));
    end
end

