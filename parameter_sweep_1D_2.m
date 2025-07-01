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
% 'k_e2', 'k_in', 'kHoff', 'kPmax'
% important: 'E_total', 'k_e', 'kPmin', 'kc'
% 'E_total', 'k_e', 'kPmin', 'kc', 
param_list = {'kEon','kEoff', 'kHon', 'E_total', 'k_e', 'kPmin', 'kc'};

% Iterate over EBindingNumber
for EBindingNumber = 2:2
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
        kHon_default = P.kHon;
        
        % Define sweep range (linear for k_e, k_e2, E_total, kc; log for others)
        switch param_to_sweep
            case 'k_e'
                param_values = 45/L_a:10/L_a:85/L_a; % Linear range
            case 'k_e2'
                base_range = 15/L_a:10/L_a:65/L_a; % Linear range
                param_values = sort(unique([base_range, default_value]));
            case 'E_total'
                param_values = 50000:10000:100000; % Linear range
            case 'kc'
                param_values = 0.2:0.1:1.0; % Linear range
            case 'k_in'
                base_range = logspace(log10(0.2), log10(10), 5); % Log range
                param_values = sort(unique([base_range, default_value]));
            case 'kEon'
                base_range = logspace(-5, -2, 5); % Log range
                param_values = sort(unique([base_range, default_value]));
            case 'kEoff'
                base_range = logspace(0, log10(100), 5); % Log range
                param_values = sort(unique([base_range, default_value]));
            case 'kHon'
                base_range = logspace(-2, 0, 5); % Log range
                param_values = sort(unique([base_range, default_value]));
            case 'kHoff'
                base_range = logspace(-3, -1.5, 5); % Log range
                param_values = sort(unique([base_range, default_value]));
            case 'kPmin'
                base_range = logspace(log10(0.05), log10(0.4), 5); % Log range
                param_values = sort(unique([base_range, default_value]));
            case 'kPmax'
                base_range = logspace(1, 2, 5); % Log range
                param_values = sort(unique([base_range, default_value]));
            otherwise
                error('Invalid parameter selected');
        end

        % Initialize arrays
        cutoff_values = zeros(1,length(param_values));
        cutoff2 = zeros(1, length(param_values));

        % Parameter sweep
        for k = 1:length(param_values)
            Ef_ss = 0;
            P.kHon = kHon_default;
            % Update the parameter
            P.(param_to_sweep) = param_values(k);
            
            if strcmp(param_to_sweep,'kHon')
                kHon_default = param_values(k);
            end
           
                  
            [r_E_BeforePas] = compute_steady_states(P, EBindingNumber + 1);
            disp('done compute steady states');
            
            KpOver_vals = linspace(1/P.kPmax, 1/P.kPmin, PAS); % Range of Kp for kPon increases linearly
            Kp_vals = 1./KpOver_vals;
            %Kp_vals = linspace(P.kPmax, P.kPmin, PAS);
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
            P.kHon = kHon_default * avg_E_bound(end);
            X = fsolve(@(xx) ode_dynamics_multipleE(xx, P), X);

            R_sol = X(1:N);
            REH_sol = X(N+1 : N+N_PAS);
            
            % Calculate cutoff position in the gene using interpolation
            ratio = (REH_sol(1:end) + R_sol(PAS:end)) / R_sol(PAS-1);
            node_indices = 1:length(ratio);
            if all(ratio >= 0.5)
                cutoff_values(k) = -1;
            else
                cutoff_values(k) = interp1(ratio, node_indices, 0.5, 'linear') * L_a;% + (PAS - 1) * L_a;
            end

            if all(ratio >= 0.1)
                cutoff2(k) = -1;
            else
                cutoff2(k) = interp1(ratio, node_indices, 0.1, 'linear') * L_a;% + (PAS - 1) * L_a;
            end
        end

        % Plot
        figure;
        if ismember(param_to_sweep, {'k_in', 'kEon', 'kEoff', 'kHon', 'kHoff', 'kPmin', 'kPmax'})
            set(gca, 'XScale', 'log');
        end
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
        ylabel('Position at which 50% cleavage (*100 bp)');
        title(['CPA Cutoff vs ', strrep(param_to_sweep, '_', '\_'), ' (EBindingNumber=', num2str(EBindingNumber), ')']);
        grid on;
        
        % Save the plot
        filename = sprintf('CPA_Cutoff_EBinding%d_%s_2.png', EBindingNumber, param_to_sweep);
        saveas(gcf, filename);
        close(gcf); % Close the figure to save memory
    end
end