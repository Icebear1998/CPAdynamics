%% ------------ ODE DYNAMICS FUNCTION WITH DIAGNOSTICS ------------
function dxdt = ode_dynamics_multipleE_diagnostic(X, P)
global N PAS Ef_ss diagnostic_log

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
    % Set initial Ef_ss once (using P.E_total as a starting guess)
    if Ef_ss == 0
        Ef_ss = P.E_total; % Initial guess
    end

    % Convert symbolic expression to a numerical function
    REvalbindEAfterPas = RE_val_bind_E(Ef_ss);
    
    % Calculate E_used - THIS IS WHERE THE ISSUE OCCURS
    E_used = sum(R(1:N)'.* RE_val_bind_E(Ef_ss)) + sum(REH' .* REvalbindEAfterPas(PAS:N));
    
    % Calculate E_free
    E_f_calculated = P.E_total - E_used;
    
    % DIAGNOSTIC LOGGING
    if isfield(P, 'diagnostic_mode') && P.diagnostic_mode
        iter = length(diagnostic_log.iteration) + 1;
        diagnostic_log.iteration(iter) = iter;
        diagnostic_log.E_total(iter) = P.E_total;
        diagnostic_log.E_used(iter) = E_used;
        diagnostic_log.E_free(iter) = E_f_calculated;
        diagnostic_log.E_free_negative(iter) = (E_f_calculated < 0);
        diagnostic_log.max_R(iter) = max(R);
        diagnostic_log.max_REH(iter) = max(REH);
        
        % Calculate average E bound per polymerase
        if sum(R) > 0
            avg_E_per_R = sum(R'.* RE_val_bind_E(Ef_ss)) / sum(R);
        else
            avg_E_per_R = 0;
        end
        if sum(REH) > 0
            avg_E_per_REH = sum(REH' .* REvalbindEAfterPas(PAS:N)) / sum(REH);
        else
            avg_E_per_REH = 0;
        end
        diagnostic_log.avg_E_bound_per_R(iter) = avg_E_per_R;
        diagnostic_log.avg_E_bound_per_REH(iter) = avg_E_per_REH;
        
        % Print warning if violation detected
        if E_f_calculated < 0
            fprintf('  ⚠️  Iter %d: E_used (%.2f) > E_total (%.2f), violation = %.2f\n', ...
                iter, E_used, P.E_total, abs(E_f_calculated));
        end
    end
    
    % Original problematic code with abs()
    E_f = abs(E_f_calculated);
    Ef_ss = E_f;

    % If E_f < 0, throw error and stop solver
    if E_f_calculated < 0
        P.is_unphysical = true; % Set flag instead of throwing error
        dxdt = 1e6 * ones(length(X), 1); % Large residual to signal fsolve to stop
        return;
    end
end

Pol_f = P.Pol_total - sum(R) - sum(REH);

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

%% Helper function for tracking iterations
function stop = trackIteration(x, optimValues, state, P)
    stop = false;
    
    if strcmp(state, 'iter')
        fprintf('  fsolve iteration %d, residual norm = %.4e\n', ...
            optimValues.iteration, optimValues.resnorm);
    end
end

