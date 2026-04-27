% Parameter sweep setup
global N PAS N_PAS Ef_ss;
syms Ef real;
Ef_ss = 0;

save_result = false;

% Model parameters (base values)
P_default.L_a = 100;
P_default.k_in    = 2;
P_default.k_e     = 65/P_default.L_a;
P_default.k_e2    = 30/P_default.L_a;
P_default.E_total = 100000;
P_default.L_total = 100000;
P_default.Pol_total = 70000;

% We are not confident here
P_default.kEon = 0.0000025;
P_default.kEoff = 0.1;
P_default.kHon = 1; 
P_default.kHoff = 1; 
P_default.kc = 0.1; 

P_default.kPon_min = 0.01; % at TSS
P_default.kPon_slope = 0.005; % determined how fast Sep2P increasing from TSS
P_default.kPoff = 1;

P_default.geneLength_bp = 25000;
P_default.PASposition   = 20000;

N = floor(P_default.geneLength_bp / P_default.L_a);
PAS = floor(P_default.PASposition / P_default.L_a);
N_PAS = N - PAS + 1;

% --- PARAMETERS TO SWEEP ---
% Available: 'k_e', 'k_e2', 'E_total', 'Pol_total', 'kc', 'kEon', 'kEoff', 
%            'k_in', 'kHoff', 'kHon', 'kPon_slope'
% Note: kEon and kEoff affect symbolic steady states and will be slower
param_list = {'kHon'};%, 'kEon', 'kEoff'};

% % Ensure ode_dynamics_multipleE is available (assumed from context)
% if ~exist('ode_dynamics_multipleE', 'file')
%     error('ode_dynamics_multipleE function not found. Please ensure it is defined.');
% end
% if ~exist('compute_steady_states', 'file')
%     error('compute_steady_states function not found. Please ensure it is defined.');
% end

% Iterate over EBindingNumber
for EBindingNumber = 1:1
    fprintf('Running for EBindingNumber = %d\n', EBindingNumber);
    
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
                param_values = 10/L_a:10/L_a:600/L_a; % Linear range
            case 'k_e2'
                base_range = 15/L_a:10/L_a:65/L_a; % Linear range
                param_values = sort(unique([base_range, default_value]));
            case 'E_total'
                param_values = 30000:20000:150000; % Linear range
            case 'Pol_total'
                param_values = 50000:20000:240000;
            case 'kc'
                param_values = 0.1:0.1:0.9; % Linear range
            case 'k_in'
                base_range = logspace(log10(0.2), log10(10), 5); % Log range
                param_values = sort(unique([base_range, default_value]));
            case 'kEon'
                base_range = logspace(-4, -2, 5); % Log range
                param_values = sort(unique([base_range, default_value]));
            case 'kEoff'
                base_range = logspace(0, log10(100), 5); % Log range
                param_values = sort(unique([base_range, default_value]));
            case 'kHon'
                param_values = 1:1:10;
%                 base_range = logspace(-1, 1, 6); % Log range
%                 param_values = sort(unique([base_range, default_value]));
            case 'kHoff'
                base_range = logspace(-3, 0, 8); % Log range
                param_values = sort(unique([base_range, default_value]));
            case 'kPon_slope'
                param_values = 0.01:0.01:0.1; % Linear range
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
            
            % Update kHon_default for this iteration if sweeping kHon
            if strcmp(param_to_sweep, 'kHon')
                kHon_default = param_values(k);
            end

            try
                [R_sol, REH_sol, P_sim] = run_termination_simulation(P_run, EBindingNumber);
            catch ME
                fprintf('    Error in run_termination_simulation: %s\n', ME.message);
                cutoff_values(k) = NaN;
                continue;
            end
            
            % --- CALCULATE TERMINATION PROFILE ---
            try
                [exit_cdf, distances_bp] = calculate_pas_usage_profile(R_sol, REH_sol, P_sim);
                cutoff_values(k) = round(interp1(exit_cdf, distances_bp, 0.5, 'linear', 'extrap'));
            catch
                cutoff_values(k) = -1;
            end
        end
        fprintf('Sweep for %s complete.\n\n', param_to_sweep);

        % --- PLOT RESULTS ---
        figure('Position', [100, 100, 800, 600]);
        hold on;
        
        % Use log scale for rate constant parameters
        log_scale_params = {'k_in', 'kEon', 'kEoff', 'kHon', 'kHoff'};
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
end
