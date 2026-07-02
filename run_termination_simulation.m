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

    % Two-step solution strategy (matches the original approach):
    %
    % Step 1: Solve ODE with base kHon. This gives the Pol II distribution
    %         and the free E concentration Ef_ss_1.  The average E bound at
    %         the PAS at that free-E level determines the effective kHon.
    %
    % Step 2: Update kHon = kHon_base * avg_E_bound(PAS, Ef_ss_1), then
    %         re-solve the ODE from the Step 1 solution.  Update Ef_ss once
    %         more from the final distribution.
    %
    % This is NOT a full fixed-point iteration but faithfully reproduces the
    % original two-step logic without global variables or a nested solver.

    options = optimoptions('fsolve', 'Display', 'off', 'FunctionTolerance', 1e-8);

    % --- Step 1: solve with base kHon ---
    P.kHon = kHon_base;
    X_base = fsolve(@(xx) ode_dynamics_multipleE(xx, P), ...
                    1e-6 * ones(P.N + P.N_PAS, 1), options);
    R_base  = max(0, X_base(1:P.N));
    REH_base = max(0, X_base(P.N+1:P.N+P.N_PAS));

    % Compute Ef_ss from the Step 1 distribution
    P.Ef_ss = solve_Efree_steady_state(R_base, REH_base, P);

    % Update kHon using avg E bound at PAS evaluated at Ef_ss_1
    avg_E_bound_1 = P.RE_val_bind_E(P.Ef_ss);
    P.kHon = kHon_base * avg_E_bound_1(P.PAS);

    % --- Step 2: re-solve with updated kHon, warm-started ---
    X_iter = fsolve(@(xx) ode_dynamics_multipleE(xx, P), X_base, options);
    R_iter  = max(0, X_iter(1:P.N));
    REH_iter = max(0, X_iter(P.N+1:P.N+P.N_PAS));

    % Final Ef_ss from the Step 2 distribution
    P.Ef_ss = solve_Efree_steady_state(R_iter, REH_iter, P);

    R_sol = R_iter;
    REH_sol = REH_iter;
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
