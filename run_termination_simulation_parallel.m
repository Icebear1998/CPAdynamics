function [R_sol, REH_sol, P] = run_termination_simulation_parallel(P, EBindingNumber)
    % This function is now self-contained and safe for parallel execution.
    syms Ef real;
    
    % Calculate node information and store it in the P struct
    P.N = floor(P.geneLength_bp / P.L_a);
    P.PAS = floor(P.PASposition / P.L_a);
    P.N_PAS = P.N - P.PAS + 1;
    
    kHon_base = P.kHon;
    
    % --- Symbolic Pre-computation ---
    [r_E_BeforePas] = compute_steady_states(P, EBindingNumber + 1);
    kPon_vals = linspace(P.kPon_min, P.kPon_max, P.PAS);
    RE_vals = sym(zeros(EBindingNumber + 1, P.N));
    for e = 1:EBindingNumber + 1
        for idx = 1:length(kPon_vals)
            RE_vals(e, idx) = subs(r_E_BeforePas(e), {'kPon', 'kPoff'}, {kPon_vals(idx), P.kPoff_const});
        end
        for idx = P.PAS+1:P.N
            RE_vals(e, idx) = subs(r_E_BeforePas(e), {'kPon', 'kPoff'}, {P.kPon_max, P.kPoff_min});
        end
    end
    P.RE_val_bind_E = matlabFunction(simplify(sum(sym(1:EBindingNumber)' .* RE_vals(2:end, :), 1)), 'Vars', {Ef});

    % --- STABLE TWO-STEP SOLVER ---
    X_guess = 1e-6 * ones(P.N + P.N_PAS, 1);
    options = optimoptions('fsolve', 'Display', 'off', 'FunctionTolerance', 1e-8);
    
    % Create an anonymous function for fsolve that passes P along
    ode_handle = @(X, P_inst) ode_dynamics_multipleE_parallel(X, P_inst);

    % Step 1: Solve for the initial steady state
    P.FirstRun = true;
    try
        X_base = fsolve(@(X) ode_handle(X, P), X_guess, options);
    catch
        error('Solver failed in Step 1.');
    end

    % Step 2: Update kHon and solve for the final steady state
    avg_E_bound = P.RE_val_bind_E(P.Ef_ss_sol);
    P.FirstRun = false;
    P.kHon = kHon_base * avg_E_bound(end);
    
    X_final = fsolve(@(X) ode_handle(X, P), X_base, options);
    
    R_sol = X_final(1:P.N);
    REH_sol = X_final(P.N+1 : P.N+P.N_PAS);
end
