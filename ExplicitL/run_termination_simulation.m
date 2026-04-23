function [R_sol, REH_sol, REHL_sol, P, r_E_BeforePas, r_P] = run_termination_simulation(P, EBindingNumber)
    % RUN_TERMINATION_SIMULATION Run single simulation with explicit L-factor
    %
    % Inputs:
    %   P - Parameter structure containing all simulation parameters
    %   EBindingNumber - Number of E factors that can bind to polymerase
    %
    % Outputs:
    %   R_sol    - Solution vector for active polymerase concentrations
    %   REH_sol  - Solution vector for E-hexamer bound concentrations
    %   REHL_sol - Solution vector for fully committed (E+L) concentrations
    %   P        - Updated parameter structure with computed values
    %   r_E_BeforePas - Symbolic steady state expressions for E binding
    %   r_P      - Symbolic steady state expressions for phosphorylation

    global N PAS N_PAS Ef_ss Lf_ss;
    syms Ef real;

    % Setup geometry
    L_a = P.L_a;
    N = floor(P.geneLength_bp / L_a);
    PAS = floor(P.PASposition / L_a);
    N_PAS = N - PAS + 1;
    kHon_base = P.kHon;

    % Compute symbolic steady states
    [r_E_BeforePas, r_P] = compute_steady_states(P, EBindingNumber + 1);

    % Setup kPon values with linear increase
    kPon_vals = P.kPon_min + P.kPon_slope * (0:N-1);
    RE_vals = sym(zeros(EBindingNumber + 1, N));

    for e = 1:(EBindingNumber + 1)
        for idx = 1:N
            RE_vals(e, idx) = subs(r_E_BeforePas(e), {'kPon', 'kPoff'}, {kPon_vals(idx), P.kPoff});
        end
    end

    % Create E binding function: avg number of E-factors bound at each node as function of Ef
    P.RE_val_bind_E = matlabFunction(simplify(sum(sym(1:EBindingNumber)' .* RE_vals(2:end, :), 1)), 'Vars', {Ef});

    % ============================================================
    % STEP 1: Solve for self-consistent E_free
    % ============================================================
    % Use N + 2*N_PAS state vector (R, REH, REHL)
    X_guess = 1e-6 * ones(N + 2*N_PAS, 1);
    options = optimoptions('fsolve', 'Display', 'off', 'FunctionTolerance', 1e-8);

    P.FirstRun = true;
    P.is_unphysical = false;
    Ef_ss = 0;
    Lf_ss = P.L_total;

    % For Step 1, use placeholder rates (kHon_base uniform, kLon uniform)
    % This gets us a reasonable E_free estimate
    avg_E_placeholder = ones(1, N);  % placeholder: uniform
    P.kHon_vals = kHon_base * avg_E_placeholder;
    P.kLon_base_vals = P.kLon * avg_E_placeholder;

    X_base = fsolve(@(xx) ode_dynamics_multipleE(xx, P), X_guess, options);

    if P.is_unphysical
        error('Unphysical result in Step 1');
    end

    % ============================================================
    % STEP 2: Solve with position-dependent rates
    % ============================================================
    % Compute avg_E_bound profile using the solved E_free
    avg_E_bound = P.RE_val_bind_E(Ef_ss);
    P.avg_E_bound = avg_E_bound;

    % Position-dependent kHon: kHon_base * avg_E_bound(l) at each node
    P.kHon_vals = kHon_base * avg_E_bound;

    % Position-dependent kLon_base: kLon * avg_E_bound(l)
    % (the actual rate in the ODE will be kLon_base_vals(n) * L_free)
    % P.kLon_base_vals = P.kLon * avg_E_bound;

    P.FirstRun = false;

    X_final = fsolve(@(xx) ode_dynamics_multipleE(xx, P), X_base, options);

    % Unpack solution
    R_sol    = X_final(1:N);
    REH_sol  = X_final(N+1 : N+N_PAS);
    REHL_sol = X_final(N+N_PAS+1 : N+2*N_PAS);
end
