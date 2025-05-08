global N PAS N_PAS Ef_ss;

%% 1. Define Model Parameters
L_a = 100;
P.k_in    = 2;
P.kEon    = 0.00025;
P.kEoff   = 10;
P.k_e     = 65/L_a;
P.k_e2    = 30/L_a;
P.E_total = 70000;
P.L_total = 100000; % Not used, kept for consistency
P.Pol_total = 70000;
P.kHon    = 4;
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

%% 2. Compute Steady-State E-Binding Rates (Optimized)
% Instead of symbolic computation, compute numerically
[r_E_BeforePas, r_k_AfterPas] = compute_steady_states(P, EBindingNumber + 1);
P.r_E_BeforePas = r_E_BeforePas;
P.r_k_AfterPas = r_k_AfterPas;
disp('Done computing steady states');

% Precompute kP values
Kp_vals = linspace(P.kPmax, P.kPmin, PAS); % Range of kP
P.Kp_vals = Kp_vals;
% We’ll compute RE_vals numerically after solving for Ef_ss

% Convert kHon_t_Ef to a function
P.kHon_t_Ef = @(Ef) r_k_AfterPas(1) * Ef; % Simplified, assuming linear dependence

%% 3. Solve the Steady State with fsolve
% Initial guess using a short ode45 run
X0 = zeros(N + N_PAS, 1);
tspan = [0 1000];
[t, X_ode] = ode45(@(t, x) ode_dynamics(t, x, P), tspan, X0);
X_init = X_ode(end,:)';

% Run fsolve
options = optimoptions('fsolve', 'Display', 'iter', 'MaxIterations', 1000, 'MaxFunctionEvaluations', 100000);
[X, fval, exitflag] = fsolve(@(xx) ode_dynamics(0, xx, P), X_init, options);
disp('Done computing fsolve');

% Check convergence
if exitflag > 0
    disp('fsolve converged successfully.');
else
    disp('fsolve did not converge. Check the initial guess or system equations.');
end

% Extract solutions
R_sol   = X(1:N);
REH_sol = X(N+1:N+N_PAS);

%% 4. Compute RE_vals Numerically
% Now that we have Ef_ss, compute RE_vals numerically
RE_vals = zeros(EBindingNumber+1, PAS);
for e = 1:EBindingNumber+1
    % Evaluate r_E_BeforePas for each kP value
    r_e = r_E_BeforePas(e);
    for i = 1:PAS
        kP_val = Kp_vals(i);
        % Substitute kP and Ef numerically
        r_e_val = r_e(kP_val, Ef_ss);
        RE_vals(e, i) = R_sol(i) * r_e_val;
    end
end
disp('Done computing RE_vals');

% Compute the average number of E molecules bound at each position (Ser2P)
avg_E_bound = zeros(1, PAS);
num_E = (0:EBindingNumber)';
total_E_bound = sum(num_E .* RE_vals, 1);
total_R = sum(RE_vals, 1);
nonzero_idx = total_R > 0;
avg_E_bound(nonzero_idx) = total_E_bound(nonzero_idx) ./ total_R(nonzero_idx);
Ser2P = avg_E_bound;

%% 5. Plot Results
l_values = (1-PAS):(N-PAS);
REH_plot = [zeros(PAS-1,1); REH_sol];

figure; hold on;
% Plot RE with different numbers of bound E
for e = 1:EBindingNumber+1
    plot((1-PAS):0, RE_vals(e,:), 'LineWidth', 2, 'DisplayName', sprintf('RE with %d E', e-1));
end
plot(l_values, R_sol, 'b-', 'LineWidth', 2.5, 'DisplayName', 'R(l)');
plot(l_values, REH_plot, 'r-', 'LineWidth', 2.5, 'DisplayName', 'REH(l)');
plot((1-PAS):0, Ser2P, 'g-', 'LineWidth', 2.5, 'DisplayName', 'Ser2P (Avg E bound)');
xlabel('Distance from PAS (bp)', 'FontSize', 14);
ylabel('Pol II (scaled)', 'FontSize', 14);
legend('show', 'Location', 'best');
title('Spatial Distribution of R, REH, and Ser2P');
hold off;

% %% 6. Perform Sanity Checks
% fprintf('\n--- Main Run Checks ---\n');
% 
% % 6.1 Pol II Conservation Check
% Pol_f = P.Pol_total - sum(R_sol) - sum(REH_sol);
% Pol_used = sum(R_sol) + sum(REH_sol);
% Pol_cons = Pol_used + Pol_f;
% fprintf('Pol II conservation check:\n');
% fprintf('  Total Pol II (expected): %d\n', P.Pol_total);
% fprintf('  Sum of R + REH + Pol_f: %g\n', Pol_cons);
% 
% % 6.2 E Conservation Check
% E_used = sum(R_sol(1:PAS)' .* sum((1:EBindingNumber) .* RE_vals(2:end, :), 1));
% E_bound = E_used + sum(REH_sol); % REH has 1 E per complex
% E_f = P.E_total - E_bound;
% E_cons = E_f + E_bound;
% fprintf('E conservation check:\n');
% fprintf('  E_total (expected): %g\n', P.E_total);
% fprintf('  E_free + E_bound = %g\n', E_cons);
% 
% % 6.3 Non-negativity Check
% fprintf('\n--- Non-Negativity Check ---\n');
% if any(R_sol < 0) || any(REH_sol < 0)
%     fprintf('ERROR: Negative values found in R or REH.\n');
% else
%     fprintf('All state variables non-negative. [PASS]\n');
% end
% 
% % 6.4 E-Binding Before PAS Check
% fprintf('\n--- E-Binding Before PAS Check ---\n');
% if all(diff(Ser2P) <= 0)
%     fprintf('Average E-binding decreases as kP decreases (kPmax to kPmin). [PASS]\n');
% else
%     fprintf('Average E-binding does not decrease monotonically with kP. [FAIL]\n');
% end
% 
% % 6.5 E_f Consistency Check
% fprintf('\n--- E_f Consistency Check ---\n');
% if abs(E_f - Ef_ss) < 1e-6
%     fprintf('E_f matches Ef_ss. [PASS]\n');
% else
%     fprintf('E_f mismatch: Computed E_f = %g, Ef_ss = %g. [FAIL]\n', E_f, Ef_ss);
% end

% %% 7. Edge Case Tests
% fprintf('\n--- Edge Case Tests ---\n');
% 
% % Edge Case 1: k_in = 0
% P_temp = P; P_temp.k_in = 0;
% [X_test] = fsolve(@(x) ode_dynamics(0, x, P_temp), X0, options);
% if max(abs(X_test)) < 1e-6
%     fprintf('Edge Case 1 (k_in = 0): No Pol II initiated. [PASS]\n');
% else
%     fprintf('Edge Case 1 (k_in = 0): Pol II is present. [FAIL]\n');
% end
% 
% % Edge Case 2: k_e = k_e2 = 0
% P_temp = P; P_temp.k_e = 0; P_temp.k_e2 = 0;
% X_init_case2 = zeros(N + N_PAS, 1); X_init_case2(1) = P.Pol_total;
% [X_test] = fsolve(@(x) ode_dynamics(0, x, P_temp), X_init_case2, options);
% if all(X_test(2:end) < 1e-6)
%     fprintf('Edge Case 2 (k_e = k_e2 = 0): Pol II stays at initiation. [PASS]\n');
% else
%     fprintf('Edge Case 2 (k_e = k_e2 = 0): Pol II spreads. [FAIL]\n');
% end
% 
% % Edge Case 3: k_e2 = 0
% P_temp = P; P_temp.k_e2 = 0;
% X_init_case3 = zeros(N + N_PAS, 1); X_init_case3(PAS) = P.Pol_total;
% [X_test] = fsolve(@(x) ode_dynamics(0, x, P_temp), X_init_case3, options);
% R_test = X_test(1:N); REH_test = X_test(N+1:end);
% accumPAS = R_test(PAS) + REH_test(1);
% otherSum = sum(R_test) + sum(REH_test) - accumPAS;
% if otherSum < 1e-3
%     fprintf('Edge Case 3 (k_e2 = 0): Pol II accumulates at PAS. [PASS]\n');
% else
%     fprintf('Edge Case 3 (k_e2 = 0): Pol II beyond PAS. [FAIL]\n');
% end
% 
% % Edge Case 4: k_c = 0
% P_temp = P; P_temp.kc = 0;
% X_init_case4 = zeros(N + N_PAS, 1); X_init_case4(N) = P.Pol_total;
% [X_test] = fsolve(@(x) ode_dynamics(0, x, P_temp), X_init_case4, options);
% R_test = X_test(1:N); REH_test = X_test(N+1:end);
% lastPosPop = R_test(N) + REH_test(N_PAS);
% othersPop = sum(R_test) + sum(REH_test) - lastPosPop;
% if othersPop < 1e-3
%     fprintf('Edge Case 4 (k_c = 0): Pol II accumulates at last node. [PASS]\n');
% else
%     fprintf('Edge Case 4 (k_c = 0): Pol II in earlier nodes. [FAIL]\n');
% end
% 
% % Edge Case 5: E_total = 0
% P_temp = P; P_temp.E_total = 0;
% [X_test] = fsolve(@(x) ode_dynamics(0, x, P_temp), X0, options);
% if max(abs(X_test(N+1:end))) < 1e-6
%     fprintf('Edge Case 5 (E_total = 0): No E-binding occurs. [PASS]\n');
% else
%     fprintf('Edge Case 5 (E_total = 0): E-binding occurs. [FAIL]\n');
% end

%% Optimized Functions

function [r_E_BeforePas, r_k_AfterPas] = compute_steady_states(P, n)
    % Compute steady-state distributions numerically
    kEon = P.kEon;
    kEoff = P.kEoff;
    kHon = P.kHon;
    kHoff = P.kHoff;
    kc = P.kc;
    kPmax = P.kPmax;

    % Construct rate matrices numerically
    A_before = construct_rate_matrix(n, kEon, kEoff);
    A_after = construct_rate_matrix3D(n, kPmax, kEon, kEoff, kHon, kHoff, kc);

    % Compute steady-state distributions
    r_Before = null(A_before);
    r_After = null(A_after);

    % Normalize
    r_Before = r_Before / sum(r_Before);
    r_After = r_After / sum(r_After);

    % Compute cumulative sums for r_E_BeforePas
    RE_index = cumsum(0:n-1) + 1;
    r_E_BeforePas = zeros(1, n);
    for i = 1:n
        RE_index_i = RE_index(i:end);
        RE_index_i = RE_index_i(RE_index_i <= length(r_Before));
        r_E_BeforePas(i) = sum(r_Before(RE_index_i));
    end

    % Convert r_E_BeforePas to a function of kP and Ef
    r_E_BeforePas = @(kP, Ef) r_E_BeforePas; % Simplified for numeric use

    % Compute after-PAS parameters
    RE_index = cumsum(0:n-1) + 1;
    result_index = [];
    for i = 2:length(RE_index)
        result_index = [result_index, (RE_index(i:end) + (i-1))];
    end
    result_index = result_index(result_index <= length(r_After));
    RI_index = sum(1:n);
    Ri = sum(r_After(1:RI_index));
    Rii = sum(r_After(RI_index:end));

    kon_t = sum(r_After(result_index)) * kHon / Ri;
    r_k_AfterPas = kon_t; % Simplified, assuming kHon_t_Ef is linear in Ef
end

function A = construct_rate_matrix(n, kEon, kEoff)
    A = zeros(n*n, n*n);
    sz = [n n];

    % kPon indices (simplified, assuming kPon = 1, kPoff = kP)
    rows_kPon = repelem(1:n, fliplr(1:n));
    cellArray1 = arrayfun(@(x) x:n, 1:n, 'UniformOutput', false);
    cols_kPon = [cellArray1{:}];
    ind_kPon = sub2ind(sz, rows_kPon, cols_kPon);

    % kEon indices
    cellArray2 = arrayfun(@(x) 1:x, 1:n, 'UniformOutput', false);
    rows_kEon = [cellArray2{:}];
    cols_kEon = repelem(1:n, 1:n);
    ind_kEon = sub2ind(sz, rows_kEon, cols_kEon);

    % Vectorized assignment for kPon
    valid_pairs = ind_kPon(1:end-2) < ind_kPon(2:end-1);
    current_idx = ind_kPon(1:end-2);
    next_idx = ind_kPon(2:end-1);
    valid_idx = current_idx(valid_pairs);
    valid_next = next_idx(valid_pairs);
    A(sub2ind(size(A), valid_next, valid_idx)) = 1;
    A(sub2ind(size(A), valid_idx, valid_next)) = 1; % kP will be substituted later

    % Vectorized assignment for kEon
    a1 = 1;
    a2 = 0;
    for i = 1:length(ind_kEon)-1
        current_idx = ind_kEon(i);
        next_idx = ind_kEon(i+1);
        if next_idx == current_idx + 1
            A(next_idx, current_idx) = a2 * kEon * 1; % Ef will be substituted later
            A(current_idx, next_idx) = a1 * kEoff;
            a1 = a1 + 1;
            a2 = a2 - 1;
        else
            a2 = a1;
            a1 = 1;
        end
    end

    % Remove all-zero rows and columns
    row_mask = any(A ~= 0, 2);
    col_mask = any(A ~= 0, 1);
    A = A(row_mask, col_mask);

    % Set diagonal elements
    for i = 1:size(A,1)
        A(i,i) = -sum(A(:,i));
    end
end

function A = construct_rate_matrix3D(n, kP, kEon, kEoff, kHon, kHoff, kc)
    A = zeros(n*n + n*n, n*n + n*n);
    sz = [n n 2];

    % kPon indices
    rows_kPon = repelem(1:n, fliplr(1:n));
    cellArray1 = arrayfun(@(x) x:n, 1:n, 'UniformOutput', false);
    cols_kPon = [cellArray1{:}];
    pages_kPon = ones(1, length(rows_kPon));
    ind_kPon = sub2ind(sz, rows_kPon, cols_kPon, pages_kPon);

    % kEon indices
    cellArray2 = arrayfun(@(x) 1:x, 1:n, 'UniformOutput', false);
    rows_kEon = [cellArray2{:}];
    cols_kEon = repelem(1:n, 1:n);
    pages_kEon = ones(1, length(rows_kEon));
    ind_kEon = sub2ind(sz, rows_kEon, cols_kEon, pages_kEon);

    shift_index = n*n;

    % Vectorized assignment for kPon
    valid_pairs = ind_kPon(1:end-2) < ind_kPon(2:end-1);
    current_idx = ind_kPon(1:end-2);
    next_idx = ind_kPon(2:end-1);
    valid_idx = current_idx(valid_pairs);
    valid_next = next_idx(valid_pairs);
    A(sub2ind(size(A), valid_next, valid_idx)) = 1;
    A(sub2ind(size(A), valid_idx, valid_next)) = kP;

    valid_idx_shift = valid_idx(valid_idx >= n+1);
    valid_next_shift = valid_next(valid_idx >= n+1);
    A(sub2ind(size(A), valid_next_shift + shift_index, valid_idx_shift + shift_index)) = 1;
    A(sub2ind(size(A), valid_idx_shift + shift_index, valid_next_shift + shift_index)) = kP;

    % Vectorized assignment for kEon
    a1 = 1;
    a2 = 0;
    for i = 1:length(ind_kEon)-1
        current_idx = ind_kEon(i);
        next_idx = ind_kEon(i+1);
        if next_idx == current_idx + 1
            A(next_idx, current_idx) = a2 * kEon * 1; % Ef will be substituted later
            A(current_idx, next_idx) = a1 * kEoff;
            if ~ismember(current_idx, ind_kPon(1:n+1))
                A(next_idx + shift_index, current_idx + shift_index) = a2 * kEon * 1;
                A(current_idx + shift_index, next_idx + shift_index) = a1 * kEoff;
            end
            a1 = a1 + 1;
            a2 = a2 - 1;
        else
            a2 = a1;
            a1 = 1;
        end
    end

    % kHon and kHoff
    b = 1;
    for i = n+1:length(ind_kPon)
        current_idx = ind_kPon(i);
        next_idx = ind_kPon(min(i+1, length(ind_kPon)));
        A(current_idx + shift_index, current_idx) = b * kHon;
        A(current_idx, current_idx + shift_index) = kHoff;
        if i == length(ind_kPon)
            b = b + 1;
        end
        if current_idx > next_idx
            b = b + 1;
        end
    end

    % Remove all-zero rows and columns
    row_mask = any(A ~= 0, 2);
    col_mask = any(A ~= 0, 1);
    A = A(row_mask, col_mask);

    % Set diagonal elements
    for i = 1:(n*n)
        A(i,i) = -sum(A(:,i));
    end
end

function dxdt = ode_dynamics(t, X, P)
    global N PAS Ef_ss

    k_in   = P.k_in;
    k_e    = P.k_e;
    k_e2   = P.k_e2;
    kHoff_t= P.kHoff;
    kc_t   = P.kc;
    kHon_t_Ef = P.kHon_t_Ef;
    r_E_BeforePas = P.r_E_BeforePas;
    Kp_vals = P.Kp_vals;

    R   = X(1:N);
    REH = X(N+1:end);
    E_f = double(P.E_total);

    % Compute E_used numerically
    RE_vals = zeros(2, PAS); % For EBindingNumber = 1
    for e = 1:2
        r_e = r_E_BeforePas(e);
        for i = 1:PAS
            RE_vals(e, i) = R(i) * r_e(Kp_vals(i), E_f);
        end
    end
    E_used = sum(R(1:PAS)' .* sum((1:1) .* RE_vals(2:end, :), 1));

    E_f = P.E_total - E_used - sum(REH);
    Ef_ss = E_f;
    if E_f < 0
        error('Negative E_f at t = %g (E_f = %g).', t, E_f);
    end

    Pol_f = P.Pol_total - sum(R) - sum(REH);
    kHon_t = kHon_t_Ef(E_f);

    dxdt = zeros(length(X),1);

    n = 1;
    dxdt(n) = Pol_f * k_in - k_e * R(n);

    for n = 2:(PAS-1)
        dxdt(n) = k_e * R(n-1) - k_e * R(n);
    end

    n = PAS;
    j = n - PAS + 1;
    dxdt(n) = k_e * R(n-1) - k_e * R(n) - kHon_t * R(n) + kHoff_t * REH(j);
    dxdt(N+j) = -k_e2 * REH(j) + kHon_t * R(n) - kHoff_t * REH(j);

    for n = (PAS+1):N
        j = n - PAS + 1;
        dxdt(n) = k_e * R(n-1) - k_e * R(n) - kHon_t * R(n) + kHoff_t * REH(j);
        dxdt(N+j) = k_e2 * REH(j-1) - k_e2 * REH(j) + kHon_t * R(n) - kHoff_t * REH(j) - kc_t * REH(j);
    end
end