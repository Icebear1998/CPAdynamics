% Parameter sweep setup
global N PAS N_PAS Ef_ss;
syms Ef real;
Ef_ss = 0;

save_result = false;

% --- BASE PARAMETER VALUES (DEFINE ONCE) ---
L_a = 100;
geneLength_bp = 25000;
PASposition = 20000;

% Calculate discretization
N = floor(geneLength_bp / L_a);
PAS = floor(PASposition / L_a);
N_PAS = N - PAS + 1;

% Define all default parameter values
P_default = struct();
P_default.L_a = L_a;
P_default.geneLength_bp = geneLength_bp;
P_default.PASposition = PASposition;
P_default.k_in = 2;
P_default.kEon = 0.00025;
P_default.kEoff = 10;
P_default.k_e = 65/L_a;
P_default.k_e2 = 30/L_a;
P_default.E_total = 70000;
P_default.L_total = 100000;
P_default.Pol_total = 70000;
P_default.kHon = 0.2;
P_default.kHoff = 0.0125;
P_default.kc = 0.05;
P_default.kPon_min = 0.01;
P_default.kPon_slope = 0.02;
P_default.kPoff = 1;

% --- PARAMETERS TO SWEEP ---
% Available: 'k_e', 'k_e2', 'E_total', 'Pol_total', 'kc', 'kEon', 'kEoff', 
%            'k_in', 'kHoff', 'kHon', 'kPon_slope'
% Note: kEon and kEoff affect symbolic steady states and will be slower
param_list = {'E_total', 'Pol_total'};%, 'kEon', 'kEoff'};

% Ensure ode_dynamics_multipleE is available (assumed from context)
if ~exist('ode_dynamics_multipleE', 'file')
    error('ode_dynamics_multipleE function not found. Please ensure it is defined.');
end
if ~exist('compute_steady_states', 'file')
    error('compute_steady_states function not found. Please ensure it is defined.');
end

% Iterate over EBindingNumber
for EBindingNumber = 2:2
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
                base_range = logspace(-2, 0, 5); % Log range
                param_values = sort(unique([base_range, default_value]));
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

        % --- DETERMINE IF SYMBOLIC PRE-COMPUTATION IS POSSIBLE ---
        % Parameters that affect symbolic steady states must be recomputed each iteration
        params_affecting_symbolic = {'kEon', 'kEoff'};
        recompute_symbolic = ismember(param_to_sweep, params_affecting_symbolic);
        
        if ~recompute_symbolic
            % --- PRE-COMPUTE SYMBOLIC EXPRESSIONS ONCE ---
            fprintf('Pre-computing symbolic steady states once (parameter %s does not affect them)...\n', param_to_sweep);
            try
                RE_val_bind_E_func = compute_symbolic_function(P, EBindingNumber);
                fprintf('Symbolic pre-computation complete.\n');
            catch ME
                fprintf('Error in compute_steady_states: %s\n', ME.message);
                continue;
            end
        else
            fprintf('Parameter %s affects symbolic steady states - will recompute for each value.\n', param_to_sweep);
        end

        % --- PARAMETER SWEEP ---
        fprintf('Starting sweep for %s...\n', param_to_sweep);
        for k = 1:length(param_values)
            fprintf('  Running %s = %.4g (%d/%d)...\n', param_to_sweep, param_values(k), k, length(param_values));
            
            % Create local copy of parameters and update for this iteration
            P_run = P;
            P_run.kHon = kHon_default;
            P_run.(param_to_sweep) = param_values(k);
            
            % Update kHon_default for this iteration if sweeping kHon
            kHon_iter = kHon_default;
            if strcmp(param_to_sweep, 'kHon')
                kHon_iter = param_values(k);
            end
            
            % Compute symbolic expressions if needed, otherwise use pre-computed
            if recompute_symbolic
                fprintf('    Computing symbolic expressions for this parameter value...\n');
                try
                    RE_val_bind_E_func = compute_symbolic_function(P_run, EBindingNumber);
                catch ME
                    fprintf('    Error in compute_steady_states: %s\n', ME.message);
                    cutoff_values(k) = NaN;
                    continue;
                end
            end
            
            % Assign the function handle to P_run
            P_run.RE_val_bind_E = RE_val_bind_E_func;

            % --- SOLVE ODE SYSTEM (Two-step approach) ---
            options = optimoptions('fsolve', 'Display', 'off', 'FunctionTolerance', 1e-8, 'MaxIterations', 1000);
            X0 = 1e-6 * ones(N + N_PAS, 1);
            
            % Step 1: Solve with initial kHon
            P_run.FirstRun = true;
            P_run.is_unphysical = false;
            X = fsolve(@(xx) ode_dynamics_multipleE(xx, P_run), X0, options);
            
            if P_run.is_unphysical || any(isnan(X)) || any(isinf(X))
                cutoff_values(k) = -1;
                continue;
            end

            % Step 2: Adjust kHon based on E binding and re-solve
            avg_E_bound = P_run.RE_val_bind_E(Ef_ss);
            P_run.FirstRun = false;
            P_run.kHon = kHon_iter * avg_E_bound(PAS);
            X_final = fsolve(@(xx) ode_dynamics_multipleE(xx, P_run), X, options);

            R_sol = X_final(1:N);
            REH_sol = X_final(N+1 : N+N_PAS);
            
            % --- CALCULATE TERMINATION PROFILE ---
            try
                [exit_cdf, distances_bp] = calculate_pas_usage_profile(R_sol, REH_sol, P_run);
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
            save_analysis_results('parameter_sweep_1D', data, P);
        end
    end
end

%% --- HELPER FUNCTION ---
function RE_val_bind_E_func = compute_symbolic_function(P_in, EBindingNumber)
    % Compute symbolic function for E binding
    % This function computes the symbolic steady states and creates a function handle
    % for evaluating the average number of E factors bound to Pol II
    
    global N;
    syms Ef real;
    
    % Compute symbolic steady states
    [r_E_BeforePas] = compute_steady_states(P_in, EBindingNumber + 1);
    
    % Set up kPon values with linear increase
    kPon_vals = P_in.kPon_min + P_in.kPon_slope * (0:N-1);
    RE_vals = sym(zeros(EBindingNumber+1, N));
    
    % Substitute kPon and kPoff values into steady state expressions
    for e = 1:EBindingNumber+1
        for i = 1:N
            kPon_val = kPon_vals(i);
            RE_vals(e, i) = subs(r_E_BeforePas(e), {'kPon', 'kPoff'}, {kPon_val, P_in.kPoff});
        end
    end
    
    % Create function handle for average E binding
    % Average E = sum of (number of E molecules) * (fraction with that many E)
    RE_val_bind_E_func = matlabFunction(simplify(sum(sym(1:EBindingNumber)' .* RE_vals(2:end, :), 1)), 'Vars', {Ef});
end