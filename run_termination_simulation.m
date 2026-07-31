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

    % Warm-start with base kHon, then solve the free-E fixed point.
    P.kHon = kHon_base;
    [~, ~, X_base] = solve_ode_checked(P, 1e-6 * ones(P.N + P.N_PAS, 1), options, 'base kHon warm start');

    % Find Ef such that the free E implied by conservation is the same Ef
    % that sets avg E at the PAS and the effective kHon used by the ODE.
    [R_sol, REH_sol, P] = solve_selfconsistent_Efree(P, kHon_base, X_base, options);
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

function [R_sol, REH_sol, P] = solve_selfconsistent_Efree(P, kHon_base, X0, options)
    if P.E_total <= 0
        P.Ef_ss = 0;
        avg_E_bound = P.RE_val_bind_E(P.Ef_ss);
        P.kHon = kHon_base * avg_E_bound(P.PAS);
        [R_sol, REH_sol] = solve_ode_checked(P, X0, options, 'self-consistent E_total <= 0');
        return;
    end

    X_warm = X0;
    constraint_Ef = @(Ef_cand) selfconsistent_efree_constraint(Ef_cand);
    fzero_options = optimset('Display', 'off', 'TolX', 1e-8);

    try
        P.Ef_ss = fzero(constraint_Ef, [0, P.E_total], fzero_options);
    catch
        P.Ef_ss = fzero(constraint_Ef, P.E_total * 0.5, fzero_options);
    end

    % Final solve at the converged Ef. This makes the returned R/REH, kHon,
    % avg-E-at-PAS, and E conservation mutually consistent.
    [R_sol, REH_sol, P] = solve_ode_at_Ef(P.Ef_ss, X_warm);

    function val = selfconsistent_efree_constraint(Ef_cand)
        [R_tmp, REH_tmp, ~] = solve_ode_at_Ef(Ef_cand, X_warm);
        X_warm = [R_tmp; REH_tmp];
        re_all = P.RE_val_bind_E(Ef_cand);
        E_bound = sum(R_tmp(:)' .* re_all) + sum(REH_tmp(:)' .* re_all(P.PAS:P.N));
        val = Ef_cand - (P.E_total - E_bound);
    end

    function [R_tmp, REH_tmp, P_tmp] = solve_ode_at_Ef(Ef_cand, x0)
        P_tmp = P;
        avg_E_bound = P_tmp.RE_val_bind_E(Ef_cand);
        P_tmp.Ef_ss = Ef_cand;
        P_tmp.kHon = kHon_base * avg_E_bound(P_tmp.PAS);
        [R_tmp, REH_tmp] = solve_ode_checked(P_tmp, x0, options, 'self-consistent Ef candidate');
    end
end

function [R, REH, X] = solve_ode_checked(P, x0, options, context)
    [X, ~, exitflag] = fsolve(@(xx) ode_dynamics_multipleE(xx, P), x0, options);
    if exitflag <= 0
        error('run_termination_simulation:FsolveFailed', ...
              'fsolve failed during %s (exitflag %d).', context, exitflag);
    end

    tol = 1e-8;
    min_raw = min(X(:));
    if min_raw < -tol
        error('run_termination_simulation:NegativeConcentration', ...
              'Raw fsolve solution during %s has negative concentration %.3g.', ...
              context, min_raw);
    end

    X = max(0, X);
    R = X(1:P.N);
    REH = X(P.N+1:P.N+P.N_PAS);
end
