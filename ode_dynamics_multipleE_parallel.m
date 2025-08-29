function dxdt = ode_dynamics_multipleE_parallel(X, P)
    % This ODE function is now self-contained and safe for parallel execution.
    
    % Unpack parameters from the struct
    k_in   = P.k_in;
    k_e    = P.k_e;
    k_e2   = P.k_e2;
    kHoff  = P.kHoff;
    kc     = P.kc;
    kHon   = P.kHon;
    RE_val_bind_E = P.RE_val_bind_E;
    
    N = P.N;
    PAS = P.PAS;
    
    R   = X(1:N);
    REH = X(N+1:end);

    % Self-consistent Ef_ss calculation
    if P.FirstRun
        Ef_ss_guess = P.E_total; % Initial guess for free E
        E_used = sum(R(1:N)' .* RE_val_bind_E(Ef_ss_guess)) + sum(REH' .* RE_val_bind_E(Ef_ss_guess)(PAS:N));
        Ef_ss = P.E_total - E_used;
        
        % Store the solved Ef_ss back into the struct for use outside fsolve
        P.Ef_ss_sol = Ef_ss;
    else
        Ef_ss = P.Ef_ss_sol; % Use the previously solved value
    end

    Pol_f = P.Pol_total - sum(R) - sum(REH);
    Pol_f = max(0, Pol_f); % Ensure non-negative
    
    dxdt = zeros(length(X),1);
    
    dxdt(1) = Pol_f*k_in - k_e*R(1);
    for n = 2:(PAS-1)
        dxdt(n) = k_e*R(n-1) - k_e*R(n);
    end

    n = PAS; j = 1;
    dxdt(n) = k_e*R(n-1) - k_e*R(n) - kHon*R(n) + kHoff*REH(j);
    dxdt(N+j) = -k_e2*REH(j) + kHon*R(n) - kHoff*REH(j) - kc*REH(j);

    for n = (PAS+1):N
        j = n - PAS + 1;
        dxdt(n) = k_e*R(n-1) - k_e*R(n) - kHon*R(n) + kHoff*REH(j);
        dxdt(N+j) = k_e2*REH(j-1) - k_e2*REH(j) + kHon*R(n) - kHoff*REH(j) - kc*REH(j);
    end
end
