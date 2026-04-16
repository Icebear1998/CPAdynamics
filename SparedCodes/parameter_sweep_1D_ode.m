% Parameter sweep setup using ODE solver to find steady state
global N PAS N_PAS Ef_ss;
syms Ef real;
Ef_ss = 0;

save_result = true;

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
param_list = {'E_total'};

% --- SCRIPT SETUP ---
if ~exist('compute_steady_states', 'file')
    error('compute_steady_states function not found. Please ensure it is defined.');
end

% Iterate over EBindingNumber
for EBindingNumber = 3:3
    fprintf('Running for EBindingNumber = %d\n', EBindingNumber);
    
    % Iterate over each parameter to sweep
    for param_idx = 1:length(param_list)
        P = P_default;
        
        param_to_sweep = param_list{param_idx};
        default_value = P.(param_to_sweep);
        kHon_default = P.kHon;
        
        % Define sweep range
        switch param_to_sweep
            case 'k_e'
                param_values = 45/L_a:10/L_a:85/L_a;
            case 'k_e2'
                base_range = 15/L_a:10/L_a:65/L_a;
                param_values = sort(unique([base_range, default_value]));
            case 'E_total'
                param_values = 50000:20000:260000;
            case 'Pol_total'
                param_values = 30000:20000:260000;
            case 'kc'
                param_values = 0.01:0.1:0.8;
            case 'k_in'
                base_range = logspace(log10(0.2), log10(10), 5);
                param_values = sort(unique([base_range, default_value]));
            case 'kEon'
                base_range = logspace(-4, -2, 5);
                param_values = sort(unique([base_range, default_value]));
            case 'kEoff'
                base_range = logspace(0, log10(100), 5);
                param_values = sort(unique([base_range, default_value]));
            case 'kHon'
                base_range = logspace(-2, 0, 5);
                param_values = sort(unique([base_range, default_value]));
            case 'kHoff'
                base_range = logspace(-3, 0, 8);
                param_values = sort(unique([base_range, default_value]));
            case 'kPon_slope'
                param_values = 0.01:0.01:0.1;
            otherwise
                error('Invalid parameter selected');
        end

        cutoff_values = zeros(1, length(param_values));

        % Determine if symbolic pre-computation is needed
        params_affecting_symbolic = {'kEon', 'kEoff'};
        recompute_symbolic = ismember(param_to_sweep, params_affecting_symbolic);
        
        if ~recompute_symbolic
            fprintf('Pre-computing symbolic steady states once...\n');
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
            
            P_run = P;
            P_run.(param_to_sweep) = param_values(k);
            
            kHon_iter = kHon_default;
            if strcmp(param_to_sweep, 'kHon')
                kHon_iter = param_values(k);
            end
            
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
            
            P_run.RE_val_bind_E = RE_val_bind_E_func;

            % --- SOLVE ODE SYSTEM TO STEADY STATE ---
            tspan = [0 1e6];
            X0 = 1e-6 * ones(N + N_PAS, 1);
            ode_options = odeset('RelTol', 1e-6, 'AbsTol', 1e-8, 'NonNegative', 1:(N + N_PAS));
            
            % Step 1: Solve with initial kHon to find Ef_ss
            P_run.FirstRun = true;
            [~, X_sol] = ode15s(@(t, xx) ode_dynamics_for_solver(t, xx, P_run), tspan, X0, ode_options);
            X = X_sol(end, :)';
            
            if Ef_ss < 0 || any(isnan(X)) || any(isinf(X))
                cutoff_values(k) = -1; % Mark as failed
                continue;
            end

            % Step 2: Adjust kHon based on E binding and re-solve
            avg_E_bound = P_run.RE_val_bind_E(Ef_ss);
            P_run.FirstRun = false;
            P_run.kHon = kHon_iter * avg_E_bound(PAS);
            
            [~, X_sol_final] = ode15s(@(t, xx) ode_dynamics_for_solver(t, xx, P_run), tspan, X, ode_options);
            X_final = X_sol_final(end, :)';

            R_sol = X_final(1:N);
            REH_sol = X_final(N+1 : N+N_PAS);
            
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
        
        log_scale_params = {'k_in', 'kEon', 'kEoff', 'kHon', 'kHoff'};
        if ismember(param_to_sweep, log_scale_params)
            set(gca, 'XScale', 'log');
        end
        
        plot(param_values, cutoff_values, 'o-', 'LineWidth', 2, 'MarkerSize', 8);
        
        default_cutoff = interp1(param_values, cutoff_values, default_value, 'linear', 'extrap');
        plot(default_value, default_cutoff, 'ro', 'MarkerSize', 10, 'LineWidth', 2, 'DisplayName', 'Default Value');
        
        xlabel([strrep(param_to_sweep, '_', '\_'), ' Value'], 'FontSize', 12);
        ylabel('CAD (bp)', 'FontSize', 12);
        title(['CAD vs ', strrep(param_to_sweep, '_', '\_'), ...
               ' (EBindingNumber=', num2str(EBindingNumber), ')'], 'FontSize', 14, 'FontWeight', 'bold');
        grid on;
        legend('show', 'Location', 'best', 'FontSize', 10);
        set(gca, 'FontSize', 10);
        box on;
        
        if save_result
            data.EBindingNumber = EBindingNumber;
            data.sweep_param = param_to_sweep;
            data.param_values = param_values;
            data.cutoff_values = cutoff_values;
            save_analysis_results('parameter_sweep_1D_ode', data, P);
        end
    end
end

%% --- HELPER FUNCTIONS ---
function RE_val_bind_E_func = compute_symbolic_function(P_in, EBindingNumber)
    global N;
    syms Ef real;
    
    [r_E_BeforePas] = compute_steady_states(P_in, EBindingNumber + 1);
    
    kPon_vals = P_in.kPon_min + P_in.kPon_slope * (0:N-1);
    RE_vals = sym(zeros(EBindingNumber+1, N));
    
    for e = 1:EBindingNumber+1
        for i = 1:N
            kPon_val = kPon_vals(i);
            RE_vals(e, i) = subs(r_E_BeforePas(e), {'kPon', 'kPoff'}, {kPon_val, P_in.kPoff});
        end
    end
    
    RE_val_bind_E_func = matlabFunction(simplify(sum(sym(1:EBindingNumber)' .* RE_vals(2:end, :), 1)), 'Vars', {Ef});
end

function dxdt = ode_dynamics_for_solver(t, X, P)
    global N PAS Ef_ss;
    
    % Silence the unused 't' variable warning
    %#ok<INUSD>
    
    k_in   = P.k_in;
    k_e    = P.k_e;
    k_e2   = P.k_e2;
    kHoff_t= P.kHoff;
    kc_t   = P.kc;
    RE_val_bind_E = P.RE_val_bind_E;
    kHon_t = P.kHon;
    
    R   = X(1:N);
    REH = X(N+1:end);
    
    if P.FirstRun
        if Ef_ss == 0
            Ef_ss = P.E_total; % Initial guess
        end
        
        REvalbindEAfterPas = RE_val_bind_E(Ef_ss);
        E_used = sum(R(1:N)'.* RE_val_bind_E(Ef_ss)) + sum(REH' .* REvalbindEAfterPas(PAS:N));
        E_f = P.E_total - E_used;
        Ef_ss = E_f;
    end
    
    Pol_f = P.Pol_total - sum(R) - sum(REH);
    
    dxdt = zeros(length(X),1);
    
    n = 1;
    dxdt(n) = Pol_f*k_in - k_e*R(n);
    
    for n = 2:(PAS-1)
        dxdt(n) = k_e*R(n-1) - k_e*R(n);
    end
    
    n = PAS;
    j = n - PAS + 1;
    dxdt(n) = k_e*R(n-1) - k_e*R(n) - kHon_t*R(n) + kHoff_t*REH(j);
    dxdt(N+j) = -k_e2*REH(j) + kHon_t*R(n) - kHoff_t*REH(j);
    
    for n = (PAS+1):N
        j = n - PAS + 1;
        dxdt(n) = k_e*R(n-1) - k_e*R(n) - kHon_t*R(n) + kHoff_t*REH(j);
        dxdt(N+j) = k_e2*REH(j-1) - k_e2*REH(j) + kHon_t*R(n) - kHoff_t*REH(j) - kc_t*REH(j);
    end
end
