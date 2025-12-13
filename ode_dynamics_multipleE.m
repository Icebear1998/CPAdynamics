%% ------------ ODE DYNAMICS FUNCTION ------------
function dxdt = ode_dynamics_multipleE(X, P)
global N PAS Ef_ss

k_in   = P.k_in;
k_e    = P.k_e;
k_e2   = P.k_e2;
kHoff_t= P.kHoff;
kc_t   = P.kc;
RE_val_bind_E = P.RE_val_bind_E;
kHon_t = P.kHon;

R   = X(1:N);
REH = X(N+1:end);

if P.FirstRun
    %% FIXED: Self-consistent E_free calculation
    
    % Initialize Ef_ss if needed (warm start for subsequent calls)
    if Ef_ss == 0
        Ef_ss = P.E_total; % Initial guess
    end
    
    % Define constraint equation: Ef = E_total - E_used(Ef)
    % This enforces conservation: E_total = E_used + E_free
    constraint_Ef = @(Ef_candidate) calculate_Efree_constraint(Ef_candidate, R, REH, P);
    
    % Solve for self-consistent E_free using fsolve
    options_Ef = optimoptions('fsolve', 'Display', 'off', ...
        'FunctionTolerance', 1e-8, 'MaxIterations', 50, ...
        'OptimalityTolerance', 1e-8);
    
    try
        % Solve with warm start (using previous Ef_ss as initial guess)
        [Ef_solution, ~, exitflag] = fsolve(constraint_Ef, Ef_ss, options_Ef);
        
        % Check for valid solution
        if exitflag <= 0 || Ef_solution < 0 || isnan(Ef_solution) || isinf(Ef_solution)
            % Try again with different initial guess
            Ef_solution = fsolve(constraint_Ef, P.E_total * 0.5, options_Ef);
        end
        
        % Final check
        if Ef_solution < 0 || isnan(Ef_solution) || isinf(Ef_solution)
            P.is_unphysical = true;
            dxdt = 1e6 * ones(length(X), 1);
            return;
        end
        
        % Update global E_free
        Ef_ss = Ef_solution;
        
    catch ME
        % If self-consistent solver fails, flag as unphysical
        warning('E_free self-consistent solver failed: %s', ME.message);
        P.is_unphysical = true;
        dxdt = 1e6 * ones(length(X), 1);
        return;
    end
end

Pol_f = P.Pol_total - sum(R) - sum(REH);
% L_f = P.L_total - sum()
%kHon_t = REvalbindEAfterPas(PAS:N)*P.kHon;

dxdt = zeros(length(X),1);

n = 1;
dxdt(n) = Pol_f*k_in - k_e*R(n);

for n = 2:(PAS-1)
    dxdt(n) = k_e*R(n-1) - k_e*R(n);
end

n = PAS;
j = n - PAS + 1;
dxdt(n) = k_e*R(n-1) - k_e*R(n) - kHon_t*R(n) + kHoff_t*REH(j);
dxdt(N+j) = -k_e2*REH(j) + kHon_t*R(n) - kHoff_t*REH(j);

for n = (PAS+1):N
    j = n - PAS + 1;
    dxdt(n) = k_e*R(n-1) - k_e*R(n) - kHon_t*R(n) + kHoff_t*REH(j);
    dxdt(N+j) = k_e2*REH(j-1) - k_e2*REH(j) + kHon_t*R(n) - kHoff_t*REH(j) - kc_t*REH(j);
end

end

%% Helper function: Calculate E_free constraint
function residual = calculate_Efree_constraint(Ef_candidate, R, REH, P)
    % Constraint: Ef_candidate = E_total - E_used(Ef_candidate)
    % Residual: Ef_candidate - (E_total - E_used(Ef_candidate)) should be zero
    
    global N PAS
    
    % Evaluate binding function at candidate E_free
    REvalbindE_candidate = P.RE_val_bind_E(Ef_candidate);
    
    % Calculate E_used with this candidate
    E_used = sum(R(1:N)'.* REvalbindE_candidate) + ...
             sum(REH' .* REvalbindE_candidate(PAS:N));
    
    % Calculate what E_free should be given this E_used
    E_free_implied = P.E_total - E_used;
    
    % Residual: difference between candidate and implied value
    residual = Ef_candidate - E_free_implied;
end
