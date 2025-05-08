global N PAS N_PAS Ef_ss kHon_ss is_unphysical;

% ------------ MODEL PARAMETERS ------------
L_a = 100;
P.k_in    = 2;
P.kEon    = 0.00025;
P.kEoff   = 10;
P.k_e     = 65 / L_a;
P.k_e2    = 30 / L_a;
P.E_total = 70000;
P.L_total = 100000;
P.Pol_total = 70000;
P.kHon    = 0.1; 
P.kHoff   = 0.0025; 
P.kc      = 0.8;

% Gene setup
geneLength_bp = 25000;
PASposition   = 20000;
N      = floor(geneLength_bp / L_a);  % total nodes
PAS    = floor(PASposition / L_a);    % node index of PAS
N_PAS  = N - PAS + 1;                 % number of nodes at/after PAS

EBindingNumber = 4;

% Parameter ranges for analysis
kPmin_range = linspace(0.1, 3.9, 3); % From 0.1 to 20
kPmax_range = linspace(30, 50, 2);  % From 1 to 100
[kPmax_grid, kPmin_grid] = meshgrid(kPmax_range, kPmin_range);
Ef_ss_values = nan(size(kPmin_grid)); % Initialize with NaN
avg_E_bound_values = nan(size(kPmin_grid)); % Initialize with NaN

% Compute steady-state probabilities
[r_E_BeforePas] = compute_steady_states(P, EBindingNumber + 1);
% Loop over parameter combinations
for i = 1:length(kPmin_range)
    for j = 1:length(kPmax_range)
        %disp(['kPmin = ', num2str(kPmin_grid(i, j)), ', kPmax = ', num2str(kPmax_grid(i, j))]);
        if kPmin_grid(i, j) >= kPmax_grid(i, j)
            continue; % Skip if kPmin >= kPmax
        end
        P.kPmin = kPmin_grid(i, j);
        P.kPmax = kPmax_grid(i, j);

        % Precompute Kp values
        Kp_vals = linspace(P.kPmax, P.kPmin, PAS);

        % Precompute RE_vals numerically
        RE_vals = zeros(EBindingNumber + 1, PAS);
        for e = 1:EBindingNumber + 1
            RE_func = matlabFunction(r_E_BeforePas(e), 'Vars', {'kP', 'Ef'});
            for k = 1:PAS
                kP_val = Kp_vals(k);
                Ef_guess = P.E_total;
                RE_vals(e, k) = RE_func(kP_val, Ef_guess);
            end
        end

        % Precompute average E bound per Pol II
        P.Avg_E = sum((0:EBindingNumber)' .* RE_vals, 1);

        % Solve steady state
        X0 = zeros(N + N_PAS, 1);
        is_unphysical = false; % Reset flag
        X = fsolve(@(xx) ode_dynamics(xx, P), X0, optimoptions('fsolve', 'Display', 'off'));

        % Check if solution is unphysical
        if is_unphysical
            Ef_ss_values(i, j) = NaN;
            avg_E_bound_values(i, j) = NaN;
            continue;
        end

         % Extract solutions
        R_sol = X(1:N);

        % Compute avg_E_bound at PAS
        avg_E_bound = zeros(1, PAS);
        for k = 1:PAS
            total_R = R_sol(k);
            if total_R > 0
                for e = 1:EBindingNumber + 1
                    RE_func = matlabFunction(r_E_BeforePas(e), 'Vars', {'kP', 'Ef'});
                    RE_vals(e, k) = RE_func(Kp_vals(k), Ef_ss);
                end
                avg_E_bound(k) = sum((0:EBindingNumber) .* RE_vals(:, k)');
            else
                avg_E_bound(k) = 0;
            end
        end

        % Store Ef_ss and avg_E_bound at PAS
        Ef_ss_values(i, j) = Ef_ss;
        avg_E_bound_values(i, j) = avg_E_bound(end); % At PAS
    end
end

% Plot contour of Ef_ss
figure;
[C, h] = contourf(kPmin_grid, kPmax_grid, Ef_ss_values, 20);
colorbar;
hold on;
contour(kPmin_grid, kPmax_grid, Ef_ss_values, [0 0], 'r-', 'LineWidth', 2, 'DisplayName', 'Ef_ss = 0');
xlabel('kPmin');
ylabel('kPmax');
title('E_f (Ef_ss) as a Function of kPmin and kPmax (kPmin < kPmax)');
legend('Ef_ss = 0', 'Location', 'best');
hold off;

% Plot contour of avg_E_bound at PAS
figure;
[C, h] = contourf(kPmin_grid, kPmax_grid, avg_E_bound_values, 20);
colorbar;
hold on;
contour(kPmin_grid, kPmax_grid, Ef_ss_values, [0 0], 'r-', 'LineWidth', 2, 'DisplayName', 'Ef_ss = 0');
xlabel('kPmin');
ylabel('kPmax');
title('Average E Bound at PAS(kPmin < kPmax)');
legend('Ef_ss = 0', 'Location', 'best');
hold off;

% ------------ ODE DYNAMICS FUNCTION ------------
function dxdt = ode_dynamics(X, P)
    global N PAS Ef_ss kHon_ss is_unphysical;

    k_in   = P.k_in;
    k_e    = P.k_e;
    k_e2   = P.k_e2;
    kHoff_t= P.kHoff;
    kc_t   = P.kc;
    Avg_E  = P.Avg_E;
    kHon_t = P.kHon;

    R   = X(1:N);
    REH = X(N+1:end);

    % Compute E_f numerically
    E_used = sum(R(1:PAS)' .* Avg_E);
    E_f = P.E_total - E_used;
    Ef_ss = E_f;
    if E_f(end) < 0
        is_unphysical = true; % Set flag
        dxdt = 1e6 * ones(length(X), 1); % Large residual to signal fsolve to stop
        return;
    end

    Pol_f = P.Pol_total - sum(R) - sum(REH);

    dxdt = zeros(length(X), 1);

    % Node 1
    n = 1;
    dxdt(n) = Pol_f * k_in - k_e * R(n);

    % Nodes 2 to PAS-1
    for n = 2:(PAS-1)
        dxdt(n) = k_e * R(n-1) - k_e * R(n);
    end

    % Node PAS
    n = PAS;
    j = n - PAS + 1;
    dxdt(n) = k_e * R(n-1) - k_e * R(n) - kHon_t * R(n) + kHoff_t * REH(j);
    dxdt(N+j) = -k_e2 * REH(j) + kHon_t * R(n) - kHoff_t * REH(j);

    % Nodes PAS+1 to N
    for n = (PAS+1):N
        j = n - PAS + 1;
        dxdt(n) = k_e * R(n-1) - k_e * R(n) - kHon_t * R(n) + kHoff_t * REH(j);
        dxdt(N+j) = k_e2 * REH(j-1) - k_e2 * REH(j) + kHon_t * R(n) - kHoff_t * REH(j) - kc_t * REH(j);
    end
end