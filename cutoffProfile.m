% Model setup for cutoff profile
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
P.kPon_min = 0.01;
P.kPon_max = 1;
P.kPoff_min = 0.1;
P.kPoff_max = 2;
P.kPoff_const = 1;
P.kPon_const = 1;

geneLength_bp = 25000;
PASposition = 20000;
N = floor(geneLength_bp / L_a);
PAS = floor(PASposition / L_a);
N_PAS = N - PAS + 1;

% Ensure required functions are available
if ~exist('ode_dynamics_multipleE', 'file')
    error('ode_dynamics_multipleE function not found. Please ensure it is defined.');
end
if ~exist('compute_steady_states', 'file')
    error('compute_steady_states function not found. Please ensure it is defined.');
end

% Fixed EBindingNumber
EBindingNumber = 5;

% Compute steady-state solution
Ef_ss = 0;
P.kHon = 0.05;

try
    [r_E_BeforePas] = compute_steady_states(P, EBindingNumber + 1);
    disp('Done computing steady states');
catch ME
    error('Error in compute_steady_states: %s', ME.message);
end

kPon_vals = linspace(P.kPon_min, P.kPon_max, PAS);
RE_kPon_vals = sym(zeros(EBindingNumber + 1, PAS));

for e = 1:EBindingNumber + 1
    for idx = 1:length(kPon_vals)
        kPon_val = kPon_vals(idx);
        RE_kPon_vals(e, idx) = subs(r_E_BeforePas(e), {'kPon', 'kPoff'}, {kPon_val, P.kPoff_const});
    end
end

P.RE_val_bind_E_kPon = matlabFunction(simplify(sum(sym(1:EBindingNumber)' .* RE_kPon_vals(2:end, :), 1)), 'Vars', {Ef});

P.FirstRun = true;
P.is_unphysical = false;
X0 = 1e-6 * ones(N + N_PAS, 1);

options = optimoptions('fsolve', 'Display', 'off', 'FunctionTolerance', 1e-8, 'MaxIterations', 1000);
X = fsolve(@(xx) ode_dynamics_multipleE(xx, P), X0, options);
if P.is_unphysical || any(isnan(X)) || any(isinf(X))
    error('Unphysical solution or convergence failure.');
end

avg_E_bound_kPon = P.RE_val_bind_E_kPon(Ef_ss);

% Recalculate kHon
P.FirstRun = false;
X1 = X;
P.kHon = P.kHon * avg_E_bound_kPon(end);
X_adj = fsolve(@(xx) ode_dynamics_multipleE(xx, P), X1, options);

R_sol = X_adj(1:N);
REH_sol = X_adj(N+1 : N+N_PAS);

% Compute ratio at each node after PAS
ratio = (REH_sol(1:end) + R_sol(PAS:end)) / R_sol(PAS-1);

% Define percentage thresholds (90% to 10% in steps of 10%)
percentages = 0:10:100;
thresholds = 1 - (percentages / 100); % Convert to ratio thresholds (0.9, 0.8, ..., 0.1)

% Initialize arrays for positions
positions = zeros(size(thresholds));

% Compute position for each threshold
for idx = 1:length(thresholds)
    t = thresholds(idx);
    node_indices = 1:length(ratio);
    pos = interp1(ratio, node_indices, t, 'linear', 'extrap') * L_a;
    positions(idx) = pos;
end


% Plot the profile
figure('Position', [100, 100, 800, 600]);
hold on;

plot(positions, percentages, 'o-', 'LineWidth', 2, 'MarkerSize', 8, ...
     'Color', [0, 0.4470, 0.7410], 'DisplayName', 'Cutoff Profile');

% Customize plot
xlabel('Position (bp)', 'FontSize', 12);
ylabel('Percentage Cutoff (%)', 'FontSize', 12);
title(['Cutoff Profile (EBindingNumber=', num2str(EBindingNumber), ')'], 'FontSize', 14, 'FontWeight', 'bold');
grid on;
legend('show', 'Location', 'best', 'FontSize', 10);
set(gca, 'FontSize', 10);
box on;

% % Save the plot
% saveas(gcf, 'Cutoff_Profile.png', 'png');
% close(gcf);