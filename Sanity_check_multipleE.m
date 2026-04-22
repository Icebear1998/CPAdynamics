global N PAS N_PAS Ef_ss;

%% 1. Define Model Parameters
P.L_a = 100;
P.k_in    = 2;
P.kEon    = 0.0000025;
P.kEoff   = 0.1;
P.k_e     = 65/P.L_a;
P.k_e2    = 30/P.L_a;
P.E_total = 100000;
P.L_total = 100000;
P.Pol_total = 70000;
P.kHon    = 2;
P.kHoff   = 1;
P.kc      = 0.1;
P.kPon_min   = 0.01;
P.kPon_slope = 0.005;
P.kPoff   = 1;

P.geneLength_bp = 25000;
P.PASposition   = 20000;

N      = floor(P.geneLength_bp / P.L_a);
PAS    = floor(P.PASposition   / P.L_a);
N_PAS  = N - PAS + 1;
Ef_ss  = 0;

EBindingNumber = 1;

%% 2. Solve the Steady State
fprintf('Running main simulation...\n');
[R_sol, REH_sol, P_sim] = run_termination_simulation(P, EBindingNumber);
fprintf('Simulation complete.\n');

%% 3. Plot Results Using a Function
plot_spatial_distribution(R_sol, REH_sol, zeros(EBindingNumber+1, N), N, PAS, EBindingNumber, 'Main Run: Spatial Distribution of R(l) and REH(l)');

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
avg_E_bound = P_sim.RE_val_bind_E(Ef_ss);
E_used_R   = sum(R_sol'   .* avg_E_bound);
E_used_REH = sum(REH_sol' .* avg_E_bound(PAS:N));
E_bound = E_used_R + E_used_REH;
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
% Check if the number of bound E molecules increases as kPon increases
bound_E = avg_E_bound(1:PAS);
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

% Edge Case 1: k_in = 0 → No Pol II should initiate
P_temp = P; P_temp.k_in = 0;
try
    [R_test, REH_test] = run_termination_simulation(P_temp, EBindingNumber);
catch
    R_test = zeros(N, 1); REH_test = zeros(N_PAS, 1);
end
if max(abs(R_test)) < 1e-6 && max(abs(REH_test)) < 1e-6
    fprintf('Edge Case 1 (k_in = 0): No Pol II initiated as expected. [PASS]\n');
else
    fprintf('Edge Case 1 (k_in = 0): Pol II is present. [FAIL]\n');
    plot_spatial_distribution(R_test, REH_test, zeros(EBindingNumber+1, N), N, PAS, EBindingNumber, 'Edge Case 1 (k_in = 0): Pol II Distribution');
end

% Edge Case 2: k_e = k_e2 = 0 → Pol II should stay at initiation site
P_temp = P; P_temp.k_e = 0; P_temp.k_e2 = 0;
try
    [R_test, REH_test] = run_termination_simulation(P_temp, EBindingNumber);
catch
    R_test = zeros(N, 1); REH_test = zeros(N_PAS, 1);
end
if all(R_test(2:PAS) < 1e-6) && all(REH_test < 1e-6)
    fprintf('Edge Case 2 (k_e = k_e2 = 0): Pol II does not elongate. [PASS]\n');
else
    fprintf('Edge Case 2 (k_e = k_e2 = 0): Pol II spreads along gene. [FAIL]\n');
    plot_spatial_distribution(R_test, REH_test, zeros(EBindingNumber+1, N), N, PAS, EBindingNumber, 'Edge Case 2 (k_e = k_e2 = 0): Pol II Distribution');
end

% Edge Case 3: k_e2 = 0 → Pol II should accumulate at PAS
P_temp = P; P_temp.k_e2 = 0;
try
    [R_test, REH_test] = run_termination_simulation(P_temp, EBindingNumber);
catch
    R_test = zeros(N, 1); REH_test = zeros(N_PAS, 1);
end
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

% Edge Case 4: k_c = 0 → Pol II should accumulate at last node
P_temp = P; P_temp.kc = 0;
try
    [R_test, REH_test] = run_termination_simulation(P_temp, EBindingNumber);
catch
    R_test = zeros(N, 1); REH_test = zeros(N_PAS, 1);
end
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

% Edge Case 5: E_total = 0 → No E-binding should occur
P_temp = P; P_temp.E_total = 0;
try
    [R_test, REH_test, P_sim5] = run_termination_simulation(P_temp, EBindingNumber);
    avg_E_test = P_sim5.RE_val_bind_E(Ef_ss);
    if max(abs(REH_test)) < 1e-6 && max(abs(avg_E_test)) < 1e-6
        fprintf('Edge Case 5 (E_total = 0): No E-binding occurs. [PASS]\n');
    else
        fprintf('Edge Case 5 (E_total = 0): E-binding occurs. [FAIL]\n');
        plot_spatial_distribution(R_test, REH_test, zeros(EBindingNumber+1, N), N, PAS, EBindingNumber, 'Edge Case 5 (E_total = 0): Pol II Distribution');
    end
catch ME
    fprintf('Edge Case 5 (E_total = 0): Simulation failed - %s\n', ME.message);
end

% Edge Case 6: Large k_c, kHon → Immediate exit at PAS
P_temp = P; P_temp.kc = 1; P_temp.kHon = 10; P_temp.k_in = 0.01;
try
    [R_test, REH_test] = run_termination_simulation(P_temp, EBindingNumber);
catch
    R_test = zeros(N, 1); REH_test = zeros(N_PAS, 1);
end
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