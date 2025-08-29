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
    % Set initial Ef_ss once (using P.E_total as a starting guess)
    if Ef_ss == 0
        Ef_ss = P.E_total; % Initial guess
    end

    % Convert symbolic expression to a numerical function
    REvalbindEAfterPas = RE_val_bind_E(Ef_ss);
    E_used = sum(R(1:N)'.* RE_val_bind_E(Ef_ss)) + sum(REH' .* REvalbindEAfterPas(PAS:N)); %sum(REH, 1); 
    %E_f = abs(P.E_total - E_used);
    E_f = max(1, P.E_total - E_used);
    Ef_ss = E_f;
    %disp({sum(RE_val_bind_E(Ef_ss)), sum(R(1:PAS)), Ef_ss});

    % If E_f < 0, throw error and stop solver
    if Ef_ss < 0
        P.is_unphysical = true; % Set flag instead of throwing error
        dxdt = 1e6 * ones(length(X), 1); % Large residual to signal fsolve to stop
        return;
    end
end

Pol_f = P.Pol_total - sum(R) - sum(REH);
% L_f = P.L_total - sum()

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
