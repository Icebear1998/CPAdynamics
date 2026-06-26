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

    P.RE_val_bind_E = @(Ef_val) interpolate_E_bound(Ef_val, Ef_grid, avg_E_bound_grid, avg_Ser2P_grid);

    % Solve system - Step 1
    X_guess = 1e-6 * ones(P.N + P.N_PAS, 1);
    options = optimoptions('fsolve', 'Display', 'off', 'FunctionTolerance', 1e-8);

    P.FirstRun = true;
    P.is_unphysical = false;

    X_base = fsolve(@(xx) ode_dynamics_multipleE(xx, P), X_guess, options);

    % Solve for Ef_ss at the end of Step 1
    R_base = X_base(1:P.N);
    REH_base = X_base(P.N+1:P.N+P.N_PAS);
    P.Ef_ss = solve_Efree_steady_state(R_base, REH_base, P);

    % Solve system - Step 2 with adjusted kHon
    avg_E_bound = P.RE_val_bind_E(P.Ef_ss);
    P.FirstRun = false;
    P.kHon = kHon_base * avg_E_bound(end);

    X_final = fsolve(@(xx) ode_dynamics_multipleE(xx, P), X_base, options);

    R_sol = X_final(1:P.N);
    REH_sol = X_final(P.N+1:P.N+P.N_PAS);

    % Update P.Ef_ss to be consistent with final steady-state
    P.Ef_ss = solve_Efree_steady_state(R_sol, REH_sol, P);
end

function [avg_E, avg_S] = interpolate_E_bound(Ef_val, Ef_grid, E_grid, S_grid)
    % Fast 1D interpolation over pre-computed grid
    avg_E = interp1(Ef_grid, E_grid, Ef_val, 'linear', 'extrap');
    if nargout > 1
        avg_S = interp1(Ef_grid, S_grid, Ef_val, 'linear', 'extrap');
    end
end

function Ef_ss = solve_Efree_steady_state(R, REH, P)
    constraint_Ef = @(Ef_cand) efree_constraint(Ef_cand, R, REH, P);
    options = optimset('Display', 'off', 'TolX', 1e-8);
    try
        Ef_ss = fzero(constraint_Ef, [0, P.E_total], options);
    catch
        % Fallback if bounds are inconsistent
        Ef_ss = fzero(constraint_Ef, P.E_total * 0.5, options);
    end
end

function val = efree_constraint(Ef_cand, R, REH, P)
    re_all = P.RE_val_bind_E(Ef_cand);
    val = Ef_cand - (P.E_total - (sum(R(:)' .* re_all) + sum(REH(:)' .* re_all(P.PAS:P.N))));
end
