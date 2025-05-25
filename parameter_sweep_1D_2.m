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
P.kHon = 0.05;
P.kHoff = 0.0025;
P.kc = 0.8;
P.kPmin = 0.1;
P.kPmax = 40;

geneLength_bp = 25000;
PASposition = 20000;
N = floor(geneLength_bp / L_a);
PAS = floor(PASposition / L_a);
N_PAS = N - PAS + 1;

% List of parameters to sweep
param_list = {'k_e', 'k_e2', 'E_total', 'kc', 'k_in', 'kEon', 'kEoff', 'kHon', 'kHoff', 'kPmin', 'kPmax'};

% Iterate over EBindingNumber
for EBindingNumber = 2:5
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
        P.kHon = 0.05;
        P.kHoff = 0.0025;
        P.kc = 0.8;
        P.kPmin = 0.1;
        P.kPmax = 40;
        
        param_to_sweep = param_list{param_idx};
        default_value = P.(param_to_sweep);
        
        % Define sweep range (linear for k_e, k_e2, E_total, kc; log for others)
        switch param_to_sweep
            case 'k_e'
                param_values = 55/L_a:10/L_a:75/L_a; % Linear range
            case 'k_e2'
                param_values = 25/L_a:10/L_a:55/L_a; % Linear range
            case 'E_total'
                param_values = 60000:10000:80000; % Linear range
            case 'kc'
                param_values = 0.2:0.2:1.0; % Linear range
            case 'k_in'
                param_values = logspace(log10(0.2), log10(10), 8); % Log range
            case 'kEon'
                param_values = logspace(-5, -2, 8); % Log range
            case 'kEoff'
                param_values = logspace(0, log10(100), 8); % Log range
            case 'kHon'
                param_values = logspace(-2, 0, 8); % Log range
            case 'kHoff'
                param_values = logspace(-3, -1.5, 8); % Log range
            case 'kPmin'
                param_values = logspace(-2, 0, 8); % Log range
            case 'kPmax'
                param_values = logspace(1, 2, 5); % Log range
            otherwise
                error('Invalid parameter selected');
        end

        % Initialize arrays
        cutoff_values = zeros(1,length(param_values));
        cutoff2 = zeros(1, length(param_values));

        % Parameter sweep
        for k = 1:length(param_values)
            Ef_ss = 0;
            % Update the parameter
            P.(param_to_sweep) = param_values(k);
            
            [r_E_BeforePas] = compute_steady_states(P, EBindingNumber + 1);
            disp('done compute steady states');

            Kp_vals = linspace(P.kPmax, P.kPmin, PAS);
            RE_vals = sym(zeros(EBindingNumber + 1, PAS));

            for e = 1:EBindingNumber+1
                for idx = 1:length(Kp_vals)
                    kP_val = Kp_vals(idx);
                    RE_vals(e, idx) = subs(r_E_BeforePas(e), {'kP'}, {kP_val});
                end
            end

            P.RE_val_bind_E = matlabFunction(simplify(sum(sym(1:EBindingNumber)' .* RE_vals(2:end, :), 1)), 'Vars', {Ef});

            P.FirstRun = true;
            P.is_unphysical = false; % Reset flag
            X0 = zeros(N + N_PAS, 1);

            X = fsolve(@(xx) ode_dynamics_multipleE(xx, P), X0);
            if P.is_unphysical
                cutoff_values(k) = -1;
                continue;
            end

            avg_E_bound = P.RE_val_bind_E(Ef_ss);

            % Recalculate kHon
            P.FirstRun = false;
            P.kHon = P.kHon * avg_E_bound(end);
            X = fsolve(@(xx) ode_dynamics_multipleE(xx, P), X);

            R_sol = X(1:N);
            REH_sol = X(N+1 : N+N_PAS);
            
            % Calculate cutoff and REH at 800 bp
            ratio = (REH_sol(1:end) + R_sol(PAS:end)) / R_sol(PAS-1);
            idx_L50 = find(ratio < 0.5, 1, 'first');
            if isempty(idx_L50)
                cutoff_values(k) = - 1;
            else
                cutoff_values(k) = idx_L50;
            end
            
            idx_L90 = find(ratio < 0.1, 1, 'first');
            if isempty(idx_L90)
                cutoff2(k) = -1;
            else
                cutoff2(k) = idx_L90;
            end
        end

        % Plot
        figure;
        plot(param_values, cutoff_values, 'o-', 'LineWidth', 2, 'MarkerSize', 8);
        hold on;
        % Highlight default parameter value
        
        plot(default_value, interp1(param_values, cutoff_values, default_value, 'linear', 'extrap'), ...
             'ro', 'MarkerSize', 8, 'LineWidth', 2, 'DisplayName', 'Default');
        for i = 1:length(param_values)
            text(param_values(i), cutoff_values(i), sprintf('REH=%.2f', cutoff2(i)), ...
                'FontSize', 5, 'HorizontalAlignment', 'right', 'VerticalAlignment', 'bottom');
        end
        xlabel([strrep(param_to_sweep, '_', '\_'), ' Value']);
        ylabel('Fraction of Polymerases Processed by 800 bp');
        title(['CPA Cutoff vs ', strrep(param_to_sweep, '_', '\_'), ' (EBindingNumber=', num2str(EBindingNumber), ')']);
        grid on;
        
        % Save the plot
        filename = sprintf('CPA_Cutoff_EBinding%d_%s.png', EBindingNumber, param_to_sweep);
        saveas(gcf, filename);
        close(gcf); % Close the figure to save memory
    end
end