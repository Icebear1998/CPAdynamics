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

% Define parameters for 2D sweep
% Parameters to sweep: (kc,kHoff); (kPmin, kPmax); (ke, kHon); 
% (E_tot, kEon/kEoff); (Emax, kPmin)
param1 = 'kPmin';
param2 = 'EBindingNumber';
default1 = P.(param1);
default2 = 3;

% Define ranges for 2D sweep
base_range = logspace(log10(0.05), log10(0.4), 6); % Log range
param1_values = sort(unique([base_range, default_value]));
param2_values = 1:1:4; % Log range

% Initialize matrix for cutoff values
cutoff_matrix = zeros(length(param2_values), length(param1_values));

% 2D Parameter Sweep
for i = 1:length(param2_values)
    for j = 1:length(param1_values)
        Ef_ss = 0;
        P.kHon = kHon_default;
        
        P.(param1) = param1_values(j);
        %P.(param2) = param2_values(i);
        EBindingNumber = param2_values(i);
        
        if strcmp(param1,'kHon')
                kHon_default = param_values1(k);
        end
        if strcmp(param2,'kHon')
                kHon_default = param_values2(k);
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
        X = fsolve(@(xx) ode_dynamics_multipleE(xx, P), X0);

        R_sol = X(1:N);
        REH_sol = X(N+1 : N+N_PAS);
        
        % Calculate cutoff position in the gene using interpolation
        ratio = (REH_sol(1:end) + R_sol(PAS:end)) / R_sol(PAS-1);
        node_indices = 1:length(ratio);
        if all(ratio >= 0.5)
            cutoff_matrix(i,j) = -1;
        else
            cutoff_matrix(i,j) = interp1(ratio, node_indices, 0.5, 'linear') * L_a;% + (PAS - 1) * L_a;
        end  
    end
end

% Plot as contour plot
figure;
[XX, YY] = meshgrid(param1_values, param2_values);
contourf(XX, YY, cutoff_matrix, 20, 'LineColor', 'none');
set(gca, 'XScale', 'log');
colorbar;
xlabel([strrep(param1, '_', '\_'), ' Value']);
ylabel([strrep(param2, '_', '\_'), ' Value']);
title(['2D Contour Plot of CPA Cutoff (EBindingNumber=', num2str(EBindingNumber), ')']);
hold on;
plot(default1, default2, 'ro', 'MarkerSize', 10, 'LineWidth', 2, 'DisplayName', 'Default');
legend('show');
grid on;

% % Save the plot
% filename = 'CPA_Cutoff_2D_k_in_k_e.png';
% saveas(gcf, filename);
% close(gcf);