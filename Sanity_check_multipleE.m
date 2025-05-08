global N PAS N_PAS Ef_ss;
syms Ef real;

%% 1. Define Model Parameters
L_a = 100;
P.k_in    = 2;
P.kEon    = 0.00025;
P.kEoff   = 10;
P.k_e     = 65/L_a;
P.k_e2    = 30/L_a;
P.E_total = 70000;
P.L_total = 100000; % Not used in this model, but kept for consistency
P.Pol_total = 70000;
P.kHon    = 0.1;
P.kHoff   = 0.0001;
P.kc      = 0.6;
P.kPmin   = 1;
P.kPmax   = 8;

geneLength_bp = 25000;
PASposition   = 20000;
N      = floor(geneLength_bp / L_a);  % total nodes
PAS    = floor(PASposition / L_a);    % node index of PAS
N_PAS  = N - PAS + 1;                 % number of nodes at/after PAS

EBindingNumber = 1;

% Compute steady-state E-binding rates before PAS
[r_E_BeforePas, r_k_AfterPas] = compute_steady_states(P, EBindingNumber + 1);

Kp_vals = linspace(P.kPmax, P.kPmin, PAS); % Range of Kp
RE_vals = sym(zeros(EBindingNumber+1, PAS));

for e = 1:EBindingNumber+1
    for i = 1:length(Kp_vals)
        kP_val = Kp_vals(i);
        RE_vals(e, i) = subs(r_E_BeforePas(e), {'kP'}, {kP_val});
    end
end

P.RE_val_bind_E = simplify(sum(sym(1:EBindingNumber)' .* RE_vals(2:end, :), 1));
P.kHon_t_Ef = matlabFunction(r_k_AfterPas(1), 'Vars', {Ef});

%% 2. Solve the Steady State with fsolve
tspan = [0 5000]; % Increased time span for better steady-state approximation
X0 = zeros(N + N_PAS, 1);
[t, X_ode] = ode45(@(t, x) ode_backbone(t, x, P), tspan, X0);
X_init = X_ode(end,:);

% Run fsolve
options = optimoptions('fsolve', 'Display', 'iter', 'MaxIterations', 1000, 'MaxFunctionEvaluations', 100000);
[X, fval, exitflag] = fsolve(@(xx) ode_dynamics(0, xx, P), X_init);

% Check convergence
if exitflag > 0
    disp('fsolve converged successfully.');
else
    disp('fsolve did not converge. Check the initial guess or system equations.');
end

% Extract solutions
R_sol   = X(1:N);
REH_sol = X(N+1:N+N_PAS);

% Compute RE_vals at steady state
for e = 1:EBindingNumber+1
    for i = 1:PAS
        RE_vals(e, i) = R_sol(i) * double(subs(RE_vals(e, i), {'Ef'}, {Ef_ss}));
    end
end

%% 3. Plot Results Using a Function
plot_spatial_distribution(R_sol, REH_sol, RE_vals, N, PAS, EBindingNumber, 'Main Run: Spatial Distribution of R(l) and REH(l)');

%% 4. Perform Sanity Checks
fprintf('\n--- Main Run Checks ---\n');

% 4.1 Pol II Conservation Check
Pol_f = P.Pol_total - sum(R_sol) - sum(REH_sol);
Pol_used = sum(R_sol) + sum(REH_sol);
Pol_cons = Pol_used + Pol_f;
fprintf('Pol II conservation check:\n');
fprintf('  Total Pol II (expected): %d\n', P.Pol_total);
fprintf('  Sum of R + REH + Pol_f at final time: %g\n', Pol_cons);

% 4.2 E Conservation Check
E_used_sym = sum(R_sol(1:PAS)' .* P.RE_val_bind_E);
E_used = double(subs(E_used_sym, 'Ef', Ef_ss));
E_bound = E_used + sum(REH_sol); % REH has 1 E per complex
E_f = P.E_total - E_bound;
E_cons = E_f + E_bound;
fprintf('E conservation check at final time:\n');
fprintf('  E_total (expected): %g\n', P.E_total);
fprintf('  E_free + E_bound = %g\n', E_cons);

% 4.3 Non-negativity Check
fprintf('\n--- Non-Negativity Check at Final Time ---\n');
if any(R_sol < 0)
    neg_idx = find(R_sol < 0);
    neg_values = R_sol(neg_idx);
    fprintf('ERROR: Negative R at indices: %s, values: %s\n',...
         mat2str(neg_idx), mat2str(neg_values));
end
if any(REH_sol < 0)
    neg_idx = find(REH_sol < 0);
    neg_values = REH_sol(neg_idx);
    fprintf('ERROR: Negative REH at indices: %s, values: %s\n',...
         mat2str(neg_idx), mat2str(neg_values));
end
if ~any(R_sol < 0) && ~any(REH_sol < 0)
    fprintf('All state variables non-negative.\n');
end

% 4.4 E-Binding Before PAS Check
fprintf('\n--- E-Binding Before PAS Check ---\n');
% Check if the number of bound E molecules increases as kP increases
bound_E = zeros(1, PAS);
for i = 1:PAS
    bound_E(i) = sum((0:EBindingNumber) .* double(RE_vals(:, i)));
end
% Since kP decreases from kPmax to kPmin, bound_E should decrease
if all(diff(bound_E) <= 0)
    fprintf('E-binding decreases as kP decreases (kPmax to kPmin). [PASS]\n');
else
    fprintf('E-binding does not decrease monotonically with kP. [FAIL]\n');
end

% 4.5 Consistency of E_f
fprintf('\n--- E_f Consistency Check ---\n');
E_f_computed = P.E_total - E_bound;
if abs(E_f_computed - Ef_ss) < 1e-6
    fprintf('E_f from computation matches Ef_ss. [PASS]\n');
else
    fprintf('E_f mismatch: Computed E_f = %g, Ef_ss = %g. [FAIL]\n', E_f_computed, Ef_ss);
end

%% 5. Edge Case Tests
fprintf('\n--- Edge Case Tests ---\n');

% Edge Case 1: k_in = 0 ? No Pol II should initiate
P_temp = P; P_temp.k_in = 0;
[X_test] = fsolve(@(x) ode_dynamics(x, P_temp), X0, options);
if max(abs(X_test)) < 1e-6
    fprintf('Edge Case 1 (k_in = 0): No Pol II initiated as expected. [PASS]\n');
else
    fprintf('Edge Case 1 (k_in = 0): Pol II is present. [FAIL]\n');
    plot_spatial_distribution(X_test(1:N), X_test(N+1:end), RE_vals*0, N, PAS, EBindingNumber, 'Edge Case 1 (k_in = 0): Pol II Distribution');
end

% Edge Case 2: k_e = k_e2 = 0 ? Pol II should stay at initiation site
P_temp = P; P_temp.k_e = 0; P_temp.k_e2 = 0;
X_init_case2 = zeros(N + N_PAS, 1);
X_init_case2(1) = P.Pol_total; % All Pol II at l = 1 in R(1)
[X_test] = fsolve(@(x) ode_dynamics(x, P_temp), X_init_case2, options);
if all(X_test(2:PAS) < 1e-6) && all(X_test(PAS+1:end) < 1e-6)
    fprintf('Edge Case 2 (k_e = k_e2 = 0): Pol II does not elongate. [PASS]\n');
else
    fprintf('Edge Case 2 (k_e = k_e2 = 0): Pol II spreads along gene. [FAIL]\n');
    plot_spatial_distribution(X_test(1:N), X_test(N+1:end), RE_vals*0, N, PAS, EBindingNumber, 'Edge Case 2 (k_e = k_e2 = 0): Pol II Distribution');
end

% Edge Case 3: k_e2 = 0 ? Pol II should accumulate at PAS
P_temp = P; P_temp.k_e2 = 0;
X_init_case3 = zeros(N + N_PAS, 1);
X_init_case3(PAS) = P.Pol_total; % All Pol II at l = PAS in R(PAS)
[X_test] = fsolve(@(x) ode_dynamics(x, P_temp), X_init_case3, options);
R_test = X_test(1:N);
REH_test = X_test(N+1:end);
accumPAS = R_test(PAS) + REH_test(1);
otherSum = sum(R_test) + sum(REH_test) - accumPAS;
if otherSum < 1e-3
    fprintf('Edge Case 3 (k_e2 = 0): Pol II accumulates at PAS. [PASS]\n');
else
    fprintf('Edge Case 3 (k_e2 = 0): Pol II also present beyond PAS? [FAIL]\n');
    fprintf('Pol II at PAS: %g\n', accumPAS);
    fprintf('Pol II elsewhere: %g\n', otherSum);
    plot_spatial_distribution(R_test, REH_test, RE_vals*0, N, PAS, EBindingNumber, 'Edge Case 3 (k_e2 = 0): Pol II Distribution');
end

% Edge Case 4: k_c = 0 ? Pol II should accumulate at last node
P_temp = P; P_temp.kc = 0;
X_init_case4 = zeros(N + N_PAS, 1);
X_init_case4(N) = P.Pol_total; % All Pol II at l = N in R(N)
[X_test] = fsolve(@(x) ode_dynamics(x, P_temp), X_init_case4, options);
R_test = X_test(1:N);
REH_test = X_test(N+1:end);
lastPosPop = R_test(N) + REH_test(N_PAS);
othersPop = sum(R_test) + sum(REH_test) - lastPosPop;
if othersPop < 1e-3
    fprintf('Edge Case 4 (k_c = 0): Pol II accumulates at last node. [PASS]\n');
else
    fprintf('Edge Case 4 (k_c = 0): Pol II found in earlier nodes? [FAIL]\n');
    fprintf('Pol II at last node: %g\n', lastPosPop);
    fprintf('Pol II elsewhere: %g\n', othersPop);
    plot_spatial_distribution(R_test, REH_test, RE_vals*0, N, PAS, EBindingNumber, 'Edge Case 4 (k_c = 0): Pol II Distribution');
end

% Edge Case 5: E_total = 0 ? No E-binding should occur
P_temp = P; P_temp.E_total = 0;
[X_test] = fsolve(@(x) ode_dynamics(x, P_temp), X0, options);
R_test = X_test(1:N);
REH_test = X_test(N+1:end);
% Recompute RE_vals with E_total = 0
RE_vals_temp = sym(zeros(EBindingNumber+1, PAS));
for e = 1:EBindingNumber+1
    for i = 1:PAS
        RE_vals_temp(e, i) = R_test(i) * double(subs(r_E_BeforePas(e), {'kP', 'Ef'}, {Kp_vals(i), 0}));
    end
end
if max(abs(REH_test)) < 1e-6 && all(abs(double(RE_vals_temp(2:end, :))) < 1e-6, 'all')
    fprintf('Edge Case 5 (E_total = 0): No E-binding occurs. [PASS]\n');
else
    fprintf('Edge Case 5 (E_total = 0): E-binding occurs. [FAIL]\n');
    plot_spatial_distribution(R_test, REH_test, RE_vals_temp, N, PAS, EBindingNumber, 'Edge Case 5 (E_total = 0): Pol II Distribution');
end

% Edge Case 6: Large k_c, kHon ? Immediate exit at PAS
P_temp = P; P_temp.kc = 1; P_temp.kHon = 10; P_temp.k_in = 0.01;
[X_test] = fsolve(@(x) ode_dynamics(x, P_temp), X0, options);
R_test = X_test(1:N);
REH_test = X_test(N+1:end);
finalPol = sum(R_test) + sum(REH_test);
if finalPol < 1e-3
    fprintf('Edge Case 6 (k_c, kHon large): Immediate exit achieved. [PASS]\n');
else
    fprintf('Edge Case 6 (k_c, kHon large): Pol II remains in system? [FAIL]\n');
    fprintf('Total Pol II in system: %g\n', finalPol);
    plot_spatial_distribution(R_test, REH_test, RE_vals*0, N, PAS, EBindingNumber, 'Edge Case 6 (k_c, kHon large): Pol II Distribution');
end

fprintf('\n--- All Edge Cases Completed ---\n');

%% Plotting Function
function plot_spatial_distribution(R, REH, RE_vals, N, PAS, EBindingNumber, plot_title)
    l_values = (1-PAS):(N-PAS);
    REH_plot = [zeros(PAS-1,1); REH(:)];  % Assume REH=0 for nodes < PAS

    figure; hold on;
    % Plot RE with different numbers of bound E
    for e = 1:EBindingNumber+1
        RE_e = double(RE_vals(e, :));
        plot((1-PAS):0, RE_e, 'LineWidth', 2, 'DisplayName', sprintf('RE with %d E', e-1));
    end
    plot(l_values, R, 'b-', 'LineWidth', 2.5, 'DisplayName', 'R(l)');
    plot(l_values, REH_plot, 'r-', 'LineWidth', 2.5, 'DisplayName', 'REH(l)');
    xlabel('Distance from PAS (bp)', 'FontSize', 14);
    ylabel('Pol II (scaled)', 'FontSize', 14);
    legend('show', 'Location', 'best');
    title(plot_title);
    hold off;
end

%% ODE Dynamics Function (unchanged, included for completeness)
function dxdt = ode_dynamics(X, P)
global N PAS Ef_ss

k_in   = P.k_in;
k_e    = P.k_e;
k_e2   = P.k_e2;
kHoff_t= P.kHoff;
kc_t   = P.kc;
RE_val_bind_E = P.RE_val_bind_E;
kHon_t_Ef  = P.kHon_t_Ef;

R   = X(1:N);
REH = X(N+1:end);
E_f = double(P.E_total);

E_used_sym = sum(R(1:PAS)' .* RE_val_bind_E);
E_used = double(subs(E_used_sym, 'Ef', E_f));

E_f = P.E_total - E_used;
Ef_ss = E_f;
% If E_f < 0, throw error and stop solver
if E_f < 0
    error('Negative E_f at t = %g (E_f = %g). Stopping simulation.', t, E_f);
end

Pol_f = P.Pol_total - sum(R) - sum(REH);

kHon_t = kHon_t_Ef(E_f);

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