%% ------------ ODE DYNAMICS FUNCTION ------------
function dxdt = ode_dynamics_multipleE(X, P)
% ODE_DYNAMICS_MULTIPLEE  Computes the RHS of the Pol II transcription/termination ODEs.
%
% Inputs:
%   X - State vector: [R (1:N); REH (1:N_PAS)]
%   P - Parameter structure containing geometry and rates:
%       - P.N, P.PAS, P.k_in, P.k_e, P.k_e2, P.kHon, P.kHoff, P.kc, P.Pol_total
%
% Output:
%   dxdt - Derivatives vector

k_in   = P.k_in;
k_e    = P.k_e;
k_e2   = P.k_e2;
kHoff_t= P.kHoff;
kc_t   = P.kc;
kHon_t = P.kHon;

N   = P.N;
PAS = P.PAS;

R   = X(1:N);
REH = X(N+1:end);

if isfield(P, 'Pol_free')
    Pol_f = P.Pol_free;
else
    Pol_f = P.Pol_total - sum(R) - sum(REH);
end

dxdt = zeros(length(X), 1);

% Node 1 (TSS)
n = 1;
dxdt(n) = Pol_f*k_in - k_e*R(n);

% Nodes 2 to PAS-1 (Pre-PAS elongation)
for n = 2:(PAS-1)
    dxdt(n) = k_e*R(n-1) - k_e*R(n);
end

% Node PAS (PAS recognition site boundary)
n = PAS;
j = n - PAS + 1;
dxdt(n) = k_e*R(n-1) - k_e*R(n) - kHon_t*R(n) + kHoff_t*REH(j);
dxdt(N+j) = -k_e2*REH(j) + kHon_t*R(n) - kHoff_t*REH(j);

% Nodes PAS+1 to N (Post-PAS elongation and termination)
for n = (PAS+1):N
    j = n - PAS + 1;
    dxdt(n) = k_e*R(n-1) - k_e*R(n) - kHon_t*R(n) + kHoff_t*REH(j);
    dxdt(N+j) = k_e2*REH(j-1) - k_e2*REH(j) + kHon_t*R(n) - kHoff_t*REH(j) - kc_t*REH(j);
end

end
