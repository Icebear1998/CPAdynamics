% Parameter setup
global N PAS N_PAS Ef_ss;
syms Ef real;
 
% Model parameters (fixed values)
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
 
P.kPon_min = 0.01; % at TSS
P.kPon_slope = 0.05; % slope of linear increase
P.kPoff = 1;
 
geneLength_bp = 25000;
PASposition = 20000;
N = floor(geneLength_bp / L_a);
PAS = floor(PASposition / L_a);
N_PAS = N - PAS + 1;
 
% Ensure ode_dynamics_multipleE is available
if ~exist('ode_dynamics_multipleE', 'file')
    error('ode_dynamics_multipleE function not found. Please ensure it is defined.');
end
if ~exist('compute_steady_states', 'file')
    error('compute_steady_states function not found. Please ensure it is defined.');
end
 
% Initialize storage for average E binding
avg_E_bind_data = zeros(1, 4); % 4 values for EBindingNumber 2 to 5
 
% Iterate over EBindingNumber
for EBindingNumber = 1:3
    fprintf('Running for EBindingNumber = %d\n', EBindingNumber);
    
    Ef_ss = 0;
    P.kHon = 0.05; % Fixed kHon value
    
    try
        [r_E_BeforePas] = compute_steady_states(P, EBindingNumber + 1);
        disp('done compute steady states');
    catch ME
        fprintf('Error in compute_steady_states for EBindingNumber = %d: %s\n', EBindingNumber, ME.message);
        continue;
    end
    
    kPon_vals = P.kPon_min + P.kPon_slope * (0:N-1); % Linear increase with slope
    RE_vals = sym(zeros(EBindingNumber + 1, N));
 
    for e = 1:EBindingNumber + 1
        for idx = 1:N
            kPon_val = kPon_vals(idx);
            RE_vals(e, idx) = subs(r_E_BeforePas(e), {'kPon', 'kPoff'}, {kPon_val, P.kPoff});
        end
    end
 
    P.RE_val_bind_E = matlabFunction(simplify(sum(sym(1:EBindingNumber)' .* RE_vals(2:end, :), 1)), 'Vars', {Ef});
 
    P.FirstRun = true;
    P.is_unphysical = false;
    X0 = 1e-6 * ones(N + N_PAS, 1);
 
    options = optimoptions('fsolve', 'Display', 'off', 'FunctionTolerance', 1e-8, 'MaxIterations', 1000);
    X = fsolve(@(xx) ode_dynamics_multipleE(xx, P), X0, options);
    if P.is_unphysical || any(isnan(X)) || any(isinf(X))
        continue;
    end
 
    avg_E_bound_kPon = P.RE_val_bind_E(Ef_ss);
 
    % Recalculate kHon
    P.FirstRun = false;
    X1 = X;
    P.kHon = P.kHon * avg_E_bound_kPon(end);
    X_adj = fsolve(@(xx) ode_dynamics_multipleE(xx, P), X1, options);
 
    R_sol = X_adj(1:N);
    REH_sol = X_adj(N+1 : N+N_PAS);
 
    % Compute average E binding at PAS
    avg_E_bind = avg_E_bound_kPon(end); % Use the last value (at PAS)
    avg_E_bind_data(EBindingNumber) = avg_E_bind; % Store for EBindingNumber 2 to 5
end
 
% Plot Average E binding vs EBindingNumber
figure('Position', [100, 100, 800, 600]);
hold on;
 
EBindingNumbers = 2:5;
plot(EBindingNumbers, avg_E_bind_data, 'o-', 'LineWidth', 2, 'MarkerSize', 8, ...
     'Color', [0, 0.4470, 0.7410], 'DisplayName', 'Average E Binding');
 
% Customize plot
xlabel('EBindingNumber', 'FontSize', 12);
ylabel('Average E Binding at PAS', 'FontSize', 12);
title('Average E Binding vs EBindingNumber (Generated at 03:55 PM EDT 24-Jun-2025)', ...
      'FontSize', 14, 'FontWeight', 'bold');
grid on;
legend('show', 'Location', 'best', 'FontSize', 10);
set(gca, 'FontSize', 10);
box on;
 
% Save the plot
% filename = 'Avg_E_Binding_vs_EBindingNumber.png';
% saveas(gcf, filename, 'png');
% print(gcf, filename, '-dpng', '-r300'); % High-resolution save
% close(gcf);

