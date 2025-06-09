% 2D Parameter Sweep Setup
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
EBindingNumber = 2;
kHon_default = P.kHon;

% Define parameter pairs for 2D sweep
param_pairs = {
    {'E_total', 'kEon'}, ...
    {'E_total', 'kEoff'}, ...
    {'kc', 'kHoff'}, ...
    {'kPmin', 'kPmax'}, ...
    {'k_e', 'kHon'}%, ...
    %{'EBindingNumber', 'kPmin'}
};

% Loop through all parameter pairs
for pair_idx = 1:length(param_pairs)
    % Reset parameters
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
            
    param1 = param_pairs{pair_idx}{1};
    param2 = param_pairs{pair_idx}{2};
    default1 = P.(param1);
    default2 = P.(param2);

    % Define ranges for each parameter
    switch param1
        case 'k_e'
            param1_values = 35/L_a:10/L_a:95/L_a;
        case 'E_total'
            param1_values = 50000:10000:100000;
        case 'kc'
            param1_values = 0.2:0.1:1.0;
        case 'kPmin'
            param1_values = logspace(log10(0.05), log10(0.4), 6);
        case 'EBindingNumber'
            param1_values = 2:5; % Range for EBindingNumber
        otherwise
            error('Invalid parameter1 selected');
    end

    switch param2
        case 'k_e'
            param2_values = 35/L_a:10/L_a:95/L_a;
        case 'kEon'
            param2_values = logspace(-5, -2, 6);
        case 'kEoff'
            param2_values = logspace(0, log10(100), 8);
        case 'kHon'
            param2_values = logspace(-2, 0, 8);
        case 'kHoff'
            param2_values = logspace(-3, -1.5, 8);
        case 'kPmax'
            param2_values = logspace(1, 2, 8);
        case 'kPmin'
            param2_values = logspace(log10(0.05), log10(0.4), 6);
        otherwise
            error('Invalid parameter2 selected');
    end

    % Initialize matrix for cutoff values
    cutoff_matrix = zeros(length(param2_values), length(param1_values));

    % 2D Parameter Sweep
    for i = 1:length(param2_values)
        for j = 1:length(param1_values)
            Ef_ss = 0;
            P.kHon = kHon_default;
            
            P.(param1) = param1_values(j);
            P.(param2) = param2_values(i);
            
            % Update kHon_default if either param is kHon
            if strcmp(param1, 'kHon')
                kHon_default = param1_values(j);
            elseif strcmp(param2, 'kHon')
                kHon_default = param2_values(i);
            end
            
            [r_E_BeforePas] = compute_steady_states(P, EBindingNumber+1);
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
            P.is_unphysical = false;
            X0 = zeros(N + N_PAS, 1);

            X = fsolve(@(xx) ode_dynamics_multipleE(xx, P), X0);
            if P.is_unphysical
                cutoff_matrix(i,j) = 0;
                continue;
            end

            avg_E_bound = P.RE_val_bind_E(Ef_ss);

            P.FirstRun = false;
            P.kHon = kHon_default * avg_E_bound(end);
            X = fsolve(@(xx) ode_dynamics_multipleE(xx, P), X);

            R_sol = X(1:N);
            REH_sol = X(N+1 : N+N_PAS);
            
            % Calculate cutoff position in the gene using interpolation
            ratio = (REH_sol(1:end) + R_sol(PAS:end)) / R_sol(PAS-1);
            node_indices = 1:length(ratio);
            cutoff_matrix(i,j) = interp1(ratio, node_indices, 0.25, 'linear') * L_a;
%             if all(ratio >= 0.5)
%                 cutoff_matrix(i,j) = -1;
%             else
%                 cutoff_matrix(i,j) = interp1(ratio, node_indices, 0.5, 'linear') * L_a;
%             end  
        end
    end

    % Plot as contour plot
    figure;
    [XX, YY] = meshgrid(param1_values, param2_values);
    contourf(XX, YY, cutoff_matrix, 20, 'LineColor', 'none');
    if ismember(param1, {'kPmin', 'kPmax', 'kEon', 'kEoff', 'kHon', 'kHoff'})
        set(gca, 'XScale', 'log');
    end
    if ismember(param2, {'kPmin', 'kPmax', 'kEon', 'kEoff', 'kHon', 'kHoff'})
        set(gca, 'YScale', 'log');
    end
    colorbar;
    xlabel([strrep(param1, '_', '\_'), ' Value']);
    ylabel([strrep(param2, '_', '\_'), ' Value']);
    title(['2D Contour Plot of CPA Cutoff (EBindingNumber=', num2str(EBindingNumber), ')']);
    hold on;
    plot(default1, default2, 'ro', 'MarkerSize', 10, 'LineWidth', 2, 'DisplayName', 'Default');
    legend('show');
    grid on;

    % Save the plot
    filename = sprintf('CPA_Cutoff_2D_%s_%s.png', param1, param2);
    saveas(gcf, filename);
    close(gcf);
end