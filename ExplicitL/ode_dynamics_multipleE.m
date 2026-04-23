%% ------------ ODE DYNAMICS FUNCTION (Explicit L-factor version) ------------
% State vector: X = [R(1:N), REH(1:N_PAS), REHL(1:N_PAS)]
%   R    = elongating Pol II
%   REH  = Pol II with E-hexamer binding (Reaction 1 complete)
%   REHL = Pol II with E-hexamer + L-factor (fully committed, can cleave)
function dxdt = ode_dynamics_multipleE(X, P)
global N PAS Ef_ss Lf_ss

k_in   = P.k_in;
k_e    = P.k_e;
k_e2   = P.k_e2;
kHoff_t= P.kHoff;
kLoff_t= P.kLoff;
kc_t   = P.kc;
RE_val_bind_E = P.RE_val_bind_E;

N_PAS = N - PAS + 1;

% Unpack state vector
R    = X(1:N);
REH  = X(N+1 : N+N_PAS);
REHL = X(N+N_PAS+1 : N+2*N_PAS);

if P.FirstRun
    %% Self-consistent E_free calculation (same as original)
    if Ef_ss == 0
        Ef_ss = P.E_total;
    end

    constraint_Ef = @(Ef_candidate) calculate_Efree_constraint(Ef_candidate, R, REH, REHL, P);

    options_Ef = optimoptions('fsolve', 'Display', 'off', ...
        'FunctionTolerance', 1e-8, 'MaxIterations', 50, ...
        'OptimalityTolerance', 1e-8);

    try
        [Ef_solution, ~, exitflag] = fsolve(constraint_Ef, Ef_ss, options_Ef);

        if exitflag <= 0 || Ef_solution < 0 || isnan(Ef_solution) || isinf(Ef_solution)
            Ef_solution = fsolve(constraint_Ef, P.E_total * 0.5, options_Ef);
        end

        if Ef_solution < 0 || isnan(Ef_solution) || isinf(Ef_solution)
            P.is_unphysical = true;
            dxdt = 1e6 * ones(length(X), 1);
            return;
        end

        Ef_ss = Ef_solution;

    catch ME
        warning('E_free self-consistent solver failed: %s', ME.message);
        P.is_unphysical = true;
        dxdt = 1e6 * ones(length(X), 1);
        return;
    end
end

% Compute free Pol II and free L-factor pools
Pol_f = P.Pol_total - sum(R) - sum(REH) - sum(REHL);

% L_free: global conservation (each REHL sequesters one L-factor)
Lf_ss = P.L_total - sum(REHL);
if Lf_ss < 0, Lf_ss = 0; end

% Get position-dependent rates
% kHon_eff(l) = kHon_base * avg_E_bound(l)   [Reaction 1: tethered E-hexamer encounter]
% kLon_eff(l) = kLon * avg_E_bound(l) * Lf   [Reaction 2: L-factor diffusion, antenna effect]
kHon_vals = P.kHon_vals;   % Precomputed: kHon_base * avg_E_bound (vector, N nodes)
kLon_base_vals = P.kLon_base_vals; % Precomputed: kLon * avg_E_bound (vector, N nodes)

% ------------ Build ODE equations ------------
dxdt = zeros(length(X), 1);

% --- R equations (elongating Pol II) ---
% Node 1: entry
n = 1;
dxdt(n) = Pol_f*k_in - k_e*R(n);

% Nodes 2 to PAS-1: pure advection
for n = 2:(PAS-1)
    dxdt(n) = k_e*R(n-1) - k_e*R(n);
end

% --- At PAS node (j=1) ---
n = PAS;
j = 1;
kHon_eff = kHon_vals(n);
kLon_eff = kLon_base_vals(n) * Lf_ss;

dxdt(n)       = k_e*R(n-1) - k_e*R(n) - kHon_eff*R(n) + kHoff_t*REH(j);
dxdt(N+j)     = kHon_eff*R(n) - kHoff_t*REH(j) - k_e2*REH(j) - kLon_eff*REH(j) + kLoff_t*REHL(j);
dxdt(N+N_PAS+j) = kLon_eff*REH(j) - kLoff_t*REHL(j) - k_e2*REHL(j);

% --- After PAS (j>=2): advection + reactions + cleavage ---
for n = (PAS+1):N
    j = n - PAS + 1;
    kHon_eff = kHon_vals(n);
    kLon_eff = kLon_base_vals(n) * Lf_ss;

    % R: elongation in/out + hexamer binding/unbinding
    dxdt(n)       = k_e*R(n-1) - k_e*R(n) - kHon_eff*R(n) + kHoff_t*REH(j);

    % REH: advection + hexamer binding from R - hexamer unbinding - L binding + L unbinding
    dxdt(N+j)     = k_e2*REH(j-1) - k_e2*REH(j) + kHon_eff*R(n) - kHoff_t*REH(j) ...
                    - kLon_eff*REH(j) + kLoff_t*REHL(j);

    % REHL: advection + L binding - L unbinding - cleavage
    dxdt(N+N_PAS+j) = k_e2*REHL(j-1) - k_e2*REHL(j) + kLon_eff*REH(j) ...
                      - kLoff_t*REHL(j) - kc_t*REHL(j);
end

end

%% Helper function: Calculate E_free constraint
function residual = calculate_Efree_constraint(Ef_candidate, R, REH, REHL, P)
    % Constraint: Ef_candidate = E_total - E_used(Ef_candidate)
    % E factors are bound on R, REH, and REHL (all carry CTD-loaded E-factors)

    global N PAS

    % Evaluate binding function at candidate E_free
    REvalbindE_candidate = P.RE_val_bind_E(Ef_candidate);

    % E_used: E-factors bound on all polymerase species at all positions
    % R(1:N) contributes across all nodes; REH and REHL at PAS:N nodes
    E_used = sum(R(1:N)' .* REvalbindE_candidate) + ...
             sum(REH' .* REvalbindE_candidate(PAS:N)) + ...
             sum(REHL' .* REvalbindE_candidate(PAS:N));

    E_free_implied = P.E_total - E_used;
    residual = Ef_candidate - E_free_implied;
end
