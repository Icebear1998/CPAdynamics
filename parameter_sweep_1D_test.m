% Parameter sweep setup
global N PAS N_PAS Ef_ss;
syms Ef real;
 
% Model parameters (base values)
L_a = 100;
P.k_in = 2;
P.kEon = 0.00025;
P.kEoff = 10;
P.k_e = 65/L_a;
P.k_e2 = 30/L_a;
P.E_total = 70000;
P.L_total = 100000;
P.Pol_total = 70000;
P.kHon = 0.2;
P.kHoff = 0.0125;
P.kc = 0.05;
P.kPon_min = 0.01;
P.kPon_max = 1;
P.kPoff_min = 0.1;
P.kPoff_max = 2;
P.kPoff_const = 1;
 
P.kPon_min = 0.01; % at TSS
P.kPon_max = 1.5; % at PAS
P.kPoff_min = 0.1; % at PAS
P.kPoff_max = 20; % at TSS
P.kPoff_const = 1;
P.kPon_const = 1;
 
geneLength_bp = 25000;
PASposition = 20000;
N = floor(geneLength_bp / L_a);
PAS = floor(PASposition / L_a);
N_PAS = N - PAS + 1;
 
% List of parameters to sweep
% 'k_e2', 'k_in', 'kHoff', 'kPmax'
% important: 'E_total', 'k_e', 'kPmin', 'kc'
% 'E_total', 'k_e', 'kPmin', 'kc', 
% 'kEon',
param_list = {'E_total'};%'kEoff', 'kHon', 'E_total', 'k_e', 'kc'};
 
% Ensure ode_dynamics_multipleE is available (assumed from context)
if ~exist('ode_dynamics_multipleE', 'file')
    error('ode_dynamics_multipleE function not found. Please ensure it is defined.');
end
if ~exist('compute_steady_states', 'file')
    error('compute_steady_states function not found. Please ensure it is defined.');
end
 
% Iterate over EBindingNumber
for EBindingNumber = 5:5
    
    fprintf('Running for EBindingNumber = %d\n', EBindingNumber);
    
    % Iterate over each parameter to sweep
    for param_idx = 1:length(param_list)
        % Reset standard parameter set
        P.k_in = 2;
        P.kEon = 0.00025;
        P.kEoff = 10;
        P.k_e = 65/L_a;
        P.k_e2 = 30/L_a;
        P.E_total = 70000;
        P.L_total = 100000;
        P.Pol_total = 70000;
        P.kHon = 0.2;
        P.kHoff = 0.0125;
        P.kc = 0.05;
        P.kPon_min = 0.01;
        P.kPon_max = 1;
        P.kPoff_min = 0.1;
        P.kPoff_max = 2;
        P.kPoff_const = 1;
        
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
                param_values = 30000:10000:200000; % Linear range
            case 'kc'
                param_values = 0.1:0.1:0.9; % Linear range
            case 'k_in'
                base_range = logspace(log10(0.2), log10(10), 5); % Log range
                param_values = sort(unique([base_range, default_value]));
            case 'kEon'
                base_range = logspace(-4, -2, 6); % Log range
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
            case 'kPmin'
                base_range = logspace(log10(0.05), log10(0.4), 5); % Log range
                param_values = sort(unique([base_range, default_value]));
            case 'kPmax'
                base_range = logspace(1, 2, 5); % Log range
                param_values = sort(unique([base_range, default_value]));
            case 'kPon_max'
                param_values = logspace(-2, 1, 6); % Log range
            case 'Pol_total'
                param_values = 50000:20000:200000;
            otherwise
                error('Invalid parameter selected');
        end
 
        % Initialize arrays
        cutoff_values = zeros(1, length(param_values));
 
        % Parameter sweep
        for k = 1:length(param_values)
            Ef_ss = 0;
            P.kHon = kHon_default;
            % Update the parameter
            P.(param_to_sweep) = param_values(k);
            
            if strcmp(param_to_sweep, 'kHon')
                kHon_default = param_values(k);
            end
            
            try
                [r_E_BeforePas] = compute_steady_states(P, EBindingNumber + 1);
                disp('done compute steady states');
            catch ME
                fprintf('Error in compute_steady_states for %s = %.4g: %s\n', param_to_sweep, param_values(k), ME.message);
                cutoff_values(k) = NaN;
                continue;
            end
            
            kPon_vals = linspace(P.kPon_min, P.kPon_max, PAS); % Range of Kp for kPon increases linearly 
            RE_vals = sym(zeros(EBindingNumber+1, N));
 
            for e = 1:EBindingNumber+1
                for i = 1:PAS
                    kPon_val = kPon_vals(i);
                    RE_vals(e, i) = subs(r_E_BeforePas(e), {'kPon', 'kPoff'}, {kPon_val, P.kPoff_const});
                end
                for i = PAS+1:N
                    RE_vals(e, i) = subs(r_E_BeforePas(e), {'kPon', 'kPoff'}, {P.kPon_max, P.kPoff_const});
                end
            end
            disp('done compute EBindingNumber');
 
            P.RE_val_bind_E = matlabFunction(simplify(sum(sym(1:EBindingNumber)' .* RE_vals(2:end, :), 1)), 'Vars', {Ef});
            disp('done compute RE_val_bind_E');
 
            P.FirstRun = true;
            P.is_unphysical = false; % Reset flag
            X0 = 1e-6 * ones(N + N_PAS, 1); % Small positive initial guess
 
            options = optimoptions('fsolve', 'Display', 'off', 'FunctionTolerance', 1e-8, 'MaxIterations', 1000);
            X = fsolve(@(xx) ode_dynamics_multipleE_FIXED(xx, P), X0, options);
            if P.is_unphysical || any(isnan(X)) || any(isinf(X))
                cutoff_values(k) = -1;
                continue;
            end
 
            avg_E_bound = P.RE_val_bind_E(Ef_ss);
 
            % Recalculate kHon
            P.FirstRun = false;
            X1 = X;
            P.kHon = kHon_default * avg_E_bound(end);
            X_adj =  fsolve(@(xx) ode_dynamics_multipleE_FIXED(xx, P), X0, options);
 
            R_sol = X_adj(1:N);
            REH_sol = X_adj(N+1 : N+N_PAS);
            
            % --- CALCULATE TERMINATION PROFILE (FLUX-BASED) ---
            try
                [exit_cdf, distances_bp] = calculate_pas_usage_profile(R_sol, REH_sol, P);
                cutoff_values(k) = round(interp1(exit_cdf, distances_bp, 0.5, 'linear', 'extrap'));
            catch
                cutoff_values(k) = -1;
            end
        end
 
        % Plot
        figure('Position', [100, 100, 800, 600]); % Set figure size
        hold on;
        if ismember(param_to_sweep, {'k_in', 'kEon', 'kEoff', 'kHon', 'kHoff', 'kPmin', 'kPmax'})
            set(gca, 'XScale', 'log');
        end
        plot(param_values, cutoff_values, 'o-', 'LineWidth', 2, 'MarkerSize', 8, ...
             'Color', [0, 0.4470, 0.7410], 'DisplayName', '50% Termination');
        
        % Highlight default parameter value with interpolation
        default_cutoff = interp1(param_values, cutoff_values, default_value, 'linear', 'extrap');
        plot(default_value, default_cutoff, 'ro', 'MarkerSize', 10, 'LineWidth', 2, ...
             'DisplayName', 'Default Value');
        
        % Customize plot
        xlabel([strrep(param_to_sweep, '_', '\_'), ' Value'], 'FontSize', 12);
        ylabel('CAD (bp)', 'FontSize', 12);
        title(['Cleavage and termination distance (CAD) vs ', strrep(param_to_sweep, '_', '\_'), ' (EBindingNumber=', num2str(EBindingNumber), ')'], ...
              'FontSize', 14, 'FontWeight', 'bold');
        grid on;
        legend('show', 'Location', 'best', 'FontSize', 10);
        set(gca, 'FontSize', 10);
        box on; % Add border around plot
        
        % Add current date and time to title
        current_time = datestr(now, 'HH:MM AM dd-mmm-yyyy');
        annotation('textbox', [0.1, 0.95, 0.8, 0.05], ...
                   'String', ['Generated at ', current_time], ...
                   'EdgeColor', 'none', 'HorizontalAlignment', 'center', 'FontSize', 8);
 
%         % --- SAVE RESULTS ---
%         % Prepare data structure for saving
%         data.EBindingNumber = EBindingNumber;
%         data.sweep_param = param_to_sweep;
%         data.param_values = param_values;
%         data.cutoff_values = cutoff_values;
%         
%         % Save results using the utility function
%         save_analysis_results('parameter_sweep_1D', data, P);
    end
end