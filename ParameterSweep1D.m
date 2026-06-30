% Parameter sweep setup
save_result = strcmpi(getenv('CPAD_FORCE_SAVE'), 'true');

% Model parameters (base values)
P_default = default_parameters();
L_a = P_default.L_a;

% --- PARAMETERS TO SWEEP ---
% Available: 'k_e', 'k_e2', 'E_total', 'Pol_total', 'kc', 'kEon', 'kEoff', 
%            'k_in', 'kHoff', 'kHon', 'kPon_slope'
% Note: kEon and kEoff affect symbolic steady states and will be slower
param_list = {'kPon_slope'};

% % Ensure ode_dynamics_multipleE is available (assumed from context)
% if ~exist('ode_dynamics_multipleE', 'file')
%     error('ode_dynamics_multipleE function not found. Please ensure it is defined.');
% end
% if ~exist('compute_steady_states', 'file')
%     error('compute_steady_states function not found. Please ensure it is defined.');
% end

EBindingNumber = 5;

% Iterate over each parameter to sweep
for param_idx = 1:length(param_list)
    % Reset to default parameter set
    P = P_default;
        
    param_to_sweep = param_list{param_idx};
    default_value = P.(param_to_sweep);
    kHon_default = P.kHon;
        
    % Define sweep range (linear for k_e, k_e2, E_total, kc; log for others)
    switch param_to_sweep
        case 'k_e'
            % Range: 1–5 kb/min = 17–83 bp/s
            % Model unit: k_e = velocity / L_a
            base_range = [17, 25, 33, 50, 67, 83] / L_a;
            param_values = sort(unique([base_range, default_value]));

        case 'k_e2'
            % Range: 0.1–0.5 s^-1 = 10–50 bp/s for L_a = 100 bp
            base_range = [10, 20, 30, 40, 50] / L_a;
            param_values = sort(unique([base_range, default_value]));

        case 'E_total'
            % Range: 50,000–200,000 molecules
            base_range = 50000:25000:200000;
            param_values = sort(unique([base_range, default_value]));

        case 'Pol_total'
            % Range: 40,000–150,000 molecules
            base_range = [40000, 60000, 80000, 100000, 120000, 150000];
            param_values = sort(unique([base_range, default_value]));

        case 'kc'
            % Range: 0.007–0.09 s^-1 after productive-assembly/steric commitment correction
            base_range = logspace(log10(0.015), log10(0.4), 8);
            param_values = sort(unique([base_range, default_value]));

        case 'k_in'
            % Not explicitly re-estimated in the parameter appendix; keep exploratory range
            base_range = logspace(log10(0.2), log10(10), 5);
            param_values = sort(unique([base_range, default_value]));

        case 'kEon'
            % Range: 6.6e-7–1.7e-6 molecule^-1 s^-1
            base_range = logspace(log10(1.5e-7), log10(5e-6), 8);
            param_values = sort(unique([base_range, default_value]));

        case 'kEoff'
            % Range: 0.04–0.4 s^-1 after multivalent CTD-binding correction
            base_range = logspace(log10(0.2), log10(2), 8);
            param_values = sort(unique([base_range, default_value]));

        case 'kHon'
            % Range: 0.04–2 s^-1
            base_range = logspace(log10(0.04), log10(2), 8);
            param_values = sort(unique([base_range, default_value]));

        case 'kHoff'
            % Range: 0.04–4 s^-1
            base_range = logspace(log10(0.05), log10(5), 8);
            param_values = sort(unique([base_range, default_value]));

        case 'kPon_slope'
            % Range: 0.001–0.01 s^-1 per node
            base_range = logspace(log10(0.001), log10(0.05), 8);
            param_values = sort(unique([base_range, default_value]));

        otherwise
            error('Invalid parameter selected');
    end

    % Initialize arrays
    cutoff_values = zeros(1, length(param_values));

    % --- PARAMETER SWEEP ---
    fprintf('Starting sweep for %s...\n', param_to_sweep);
    for k = 1:length(param_values)
        fprintf('  Running %s = %.4g (%d/%d)...\n', param_to_sweep, param_values(k), k, length(param_values));
        
        % Create local copy of parameters and update for this iteration
        P_run = P_default;
        P_run.kHon = kHon_default;
        P_run.(param_to_sweep) = param_values(k);

        try
            [R_sol, REH_sol, P_sim] = run_termination_simulation(P_run, EBindingNumber);
        catch ME
            fprintf('    Error in run_termination_simulation: %s\n', ME.message);
            cutoff_values(k) = NaN;
            continue;
        end
        
        % --- CALCULATE TERMINATION PROFILE ---
        try
            [exit_cdf, distances_bp, cutoff_values(k)] = calculate_pas_cleavage_profile(R_sol, REH_sol, P_sim, 'PercentCleavage', 50);
            cutoff_values(k) = round(cutoff_values(k));
        catch
            cutoff_values(k) = -1;
        end
    end
    fprintf('Sweep for %s complete.\n\n', param_to_sweep);

    % --- PLOT RESULTS ---
    figure('Position', [100, 100, 800, 600]);
    hold on;
    
    % Use log scale for rate constant parameters
    log_scale_params = {'k_in', 'kEon', 'kEoff', 'kHon', 'kHoff', 'kPon_slope'};
    if ismember(param_to_sweep, log_scale_params)
        set(gca, 'XScale', 'log');
    end
    
    % Plot sweep results
    plot(param_values, cutoff_values, 'o-', 'LineWidth', 2, 'MarkerSize', 8, ...
         'Color', [0, 0.4470, 0.7410], 'DisplayName', '50% Termination');
    
    % Highlight default parameter value
    default_cutoff = interp1(param_values, cutoff_values, default_value, 'linear', 'extrap');
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
        data.EBindingNumber = EBindingNumber;
        data.sweep_param = param_to_sweep;
        data.param_values = param_values;
        data.cutoff_values = cutoff_values;
        save_analysis_results('parameter_sweep_1D', data, P_default);
    end
end

