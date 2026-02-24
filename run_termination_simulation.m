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
    %   r_E_BeforePas - (Optional) Symbolic steady state expressions for E binding
    
    global N PAS N_PAS Ef_ss;
    syms Ef real;
    
    % Setup geometry
    L_a = P.L_a;
    N = floor(P.geneLength_bp / L_a);
    PAS = floor(P.PASposition / L_a);
    N_PAS = N - PAS + 1;
    kHon_base = P.kHon;
    
    % Compute steady states (with file-based cache keyed by EBindingNumber, kEon, kEoff)
    % kEon and kEoff are baked numerically into the symbolic result, so all three
    % parameters are needed to uniquely identify the cached expression.
    cacheDir = 'SymbolicCache';
    if ~exist(cacheDir, 'dir'), mkdir(cacheDir); end
    
    cacheFile = fullfile(cacheDir, sprintf('ss_N%d_kEon%s_kEoff%s.mat', ...
        EBindingNumber, ...
        strrep(num2str(P.kEon,  '%.6g'), '.', 'p'), ...
        strrep(num2str(P.kEoff, '%.6g'), '.', 'p')));
    
    if exist(cacheFile, 'file')
        fprintf('  [Cache] Loading symbolic steady states from %s\n', cacheFile);
        loaded = load(cacheFile, 'r_E_BeforePas', 'r_P');
        r_E_BeforePas = loaded.r_E_BeforePas;
        r_P = loaded.r_P;
    else
        fprintf('  [Cache] Computing symbolic steady states (EBindingNumber=%d, kEon=%.4g, kEoff=%.4g)...\n', ...
            EBindingNumber, P.kEon, P.kEoff);
        [r_E_BeforePas, r_P] = compute_steady_states(P, EBindingNumber + 1);
        save(cacheFile, 'r_E_BeforePas', 'r_P');
        fprintf('  [Cache] Saved to %s\n', cacheFile);
    end
    
    % Setup kPon values with linear increase
    kPon_vals = P.kPon_min + P.kPon_slope * (0:N-1);
    RE_vals = sym(zeros(EBindingNumber + 1, N));
    
    for e = 1:(EBindingNumber + 1)
        for idx = 1:N
            RE_vals(e, idx) = subs(r_E_BeforePas(e), {'kPon', 'kPoff'}, {kPon_vals(idx), P.kPoff});
        end
    end
    
    % Create E binding function
    P.RE_val_bind_E = matlabFunction(simplify(sum(sym(1:EBindingNumber)' .* RE_vals(2:end, :), 1)), 'Vars', {Ef});
    
    % Solve system - Step 1
    X_guess = 1e-6 * ones(N + N_PAS, 1);
    options = optimoptions('fsolve', 'Display', 'off', 'FunctionTolerance', 1e-8);
    
    P.FirstRun = true;
    P.is_unphysical = false;
    Ef_ss = 0;
    
    X_base = fsolve(@(xx) ode_dynamics_multipleE(xx, P), X_guess, options);
    
    if P.is_unphysical
        error('Unphysical result in Step 1');
    end
    
    % Solve system - Step 2 with adjusted kHon
    avg_E_bound = P.RE_val_bind_E(Ef_ss);
    P.FirstRun = false;
    P.kHon = kHon_base * avg_E_bound(end);
    
    
    X_final = fsolve(@(xx) ode_dynamics_multipleE(xx, P), X_base, options);
    
    R_sol = X_final(1:N);
    REH_sol = X_final(N+1:N+N_PAS);
end
