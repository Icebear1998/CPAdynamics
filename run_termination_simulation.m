function [R_sol, REH_sol, P, r_E_BeforePas, r_P] = run_termination_simulation(P, EBindingNumber)
    % RUN_TERMINATION_SIMULATION Run single simulation for given E binding number
    %
    % Inputs:
    %   P - Parameter structure containing all simulation parameters
    %   EBindingNumber - Number of E factors that can bind to polymerase
    %
    % Outputs:
    %   R_sol - Solution vector for active polymerase concentrations
    %   REH_sol - Solution vector for terminating polymerase concentrations
    %   P - Updated parameter structure with computed values
    %   r_E_BeforePas - Unused (kept for call-site compatibility); always []
    %   r_P - Unused (kept for call-site compatibility); always []

    % Setup geometry and store in P
    L_a = P.L_a;
    P.N = floor(P.geneLength_bp / L_a);
    P.PAS = floor(P.PASposition / L_a);
    P.N_PAS = P.N - P.PAS + 1;
    kHon_base = P.kHon;

    % Setup kPon values with linear increase
    kPon_vals = P.kPon_min + P.kPon_slope * (0:P.N-1);
    n_states = EBindingNumber + 1;

    % Symbolic outputs are no longer computed; numerical approach is used throughout.
    r_E_BeforePas = [];
    r_P = [];

    % OPTIMIZATION: Pre-compute the SVDs over a grid of Ef_val and interpolate
    num_grid_points = 100;
    Ef_grid = linspace(0, P.E_total, num_grid_points);
    avg_E_bound_grid = zeros(num_grid_points, length(kPon_vals));
    avg_Ser2P_grid = zeros(num_grid_points, length(kPon_vals));

    for i = 1:num_grid_points
        [avg_E_bound_grid(i, :), avg_Ser2P_grid(i, :)] = compute_avg_E_bound_numerical(Ef_grid(i), kPon_vals, P.kPoff, P.kEon, P.kEoff, n_states);
    end

    if P.E_total <= 0
        P.RE_val_bind_E = @(Ef_val) constant_E_bound(Ef_val, avg_E_bound_grid, avg_Ser2P_grid);
    else
        P.RE_val_bind_E = @(Ef_val) interpolate_E_bound(Ef_val, Ef_grid, avg_E_bound_grid, avg_Ser2P_grid);
    end

    options = optimoptions('fsolve', 'Display', 'off', 'FunctionTolerance', 1e-8);

    % Warm-start the Pol II solution with the base kHon.
    P.kHon = kHon_base;
    X_iter = fsolve(@(xx) ode_dynamics_multipleE(xx, P), ...
                    1e-6 * ones(P.N + P.N_PAS, 1), options);

    % Iterate until the free E used to set kHon matches the free E implied
    % by E conservation for the resulting Pol II distribution.
    Ef_current = max(0, P.E_total);
    max_iterations = 100;
    tolerance = 1e-7;
    damping = 0.5;

    for iteration = 1:max_iterations
        avg_E_bound = P.RE_val_bind_E(Ef_current);
        P.kHon = kHon_base * avg_E_bound(P.PAS);

        X_iter = fsolve(@(xx) ode_dynamics_multipleE(xx, P), X_iter, options);
        R_sol = max(0, X_iter(1:P.N));
        REH_sol = max(0, X_iter(P.N+1:P.N+P.N_PAS));

        E_bound = sum(R_sol(:)' .* avg_E_bound) + ...
                  sum(REH_sol(:)' .* avg_E_bound(P.PAS:P.N));
        Ef_implied = P.E_total - E_bound;

        if abs(Ef_implied - Ef_current) <= tolerance
            P.Ef_ss = Ef_current;
            return;
        end

        Ef_current = Ef_current + damping * (Ef_implied - Ef_current);
        Ef_current = min(P.E_total, max(0, Ef_current));
    end

    error('run_termination_simulation:EfreeNotConverged', ...
          'Self-consistent free E did not converge after %d iterations.', max_iterations);
end

function [avg_E, avg_S] = interpolate_E_bound(Ef_val, Ef_grid, E_grid, S_grid)
    % Fast 1D interpolation over pre-computed grid
    avg_E = interp1(Ef_grid, E_grid, Ef_val, 'linear', 'extrap');
    if nargout > 1
        avg_S = interp1(Ef_grid, S_grid, Ef_val, 'linear', 'extrap');
    end
end

function [avg_E, avg_S] = constant_E_bound(~, E_grid, S_grid)
    % Degenerate E_total <= 0 case: avoid interp1 on a repeated zero grid.
    avg_E = E_grid(1, :);
    if nargout > 1
        avg_S = S_grid(1, :);
    end
end
