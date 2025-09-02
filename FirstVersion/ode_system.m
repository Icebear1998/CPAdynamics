function dxdt = ode_system(t,X, P)
%global k_in k_c kE_off kL_on kL_off k_e k_e2 E_total L_total N PAS Pol_total kE_on_f;
global N PAS Pol_total;

k_in = P.k_in;      % Define k_in
k_c = P.k_c;       % Define k_c
kE_off = P.kE_off;    % Define kE_off
kL_on = P.kL_on;     % Define kL_on
kL_off = P.kL_off;    % Define kL_off
k_e = P.k_e;      % Define k_e
k_e2 = P.k_e2;      % Define k_e2
E_total = P.E_total;   % Define E_total
L_total = P.L_total;   % Define L_total
kE_on_f = P.kE_on_f;

% Extract variables from X
R = X(1:N);
RE = X(N+1:2*N);
REL = X(2*N+1:end);

E_f = E_total - sum(RE) - sum(REL);
L_f = L_total - sum(REL);
Pol_f = Pol_total - sum(R) - sum(RE) - sum(REL);


dxdt = zeros(length(X), 1);  % Initialize the output vector for the ODEs

% ODE for l = 1 (initial condition)
n = 1;
    dxdt(n) = k_in * Pol_f - k_e * R(n) - kE_on_f(n) * R(n) * E_f + kE_off * RE(n);  % d/dt(R)
    dxdt(N+n) = - k_e * RE(n) + kE_on_f(n) * R(n) * E_f - kE_off * RE(n);    % d/dt(RE)

for n = 2:PAS-1
    dxdt(n) = k_e * R(n-1) - k_e * R(n) - kE_on_f(n) * R(n) * E_f + kE_off * RE(n);  % d/dt(R)
    dxdt(N+n) = k_e * RE(n-1) - k_e * RE(n) + kE_on_f(n) * R(n) * E_f - kE_off * RE(n);  % d/dt(RE)
end

n = PAS;
dxdt(n) = k_e * R(n-1) - k_e * R(n) - kE_on_f(n) * R(n) * E_f + kE_off * RE(n);  % d/dt(R)
dxdt(N+n) = k_e2 * RE(n-1) - k_e2 * RE(n) + kE_on_f(n) * R(n) * E_f - kE_off * RE(n) - kL_on * RE(n) * L_f + kL_off * REL(n-PAS+1);  % d/dt(RE)
dxdt(2*N+n-PAS+1) = - k_e2 * REL(n-PAS+1) - k_c * REL(n-PAS+1) + kL_on * RE(n) * L_f - kL_off * REL(n-PAS+1); %d/dt(REL)

for n = PAS+1:N
    dxdt(n) = k_e * R(n-1) - k_e * R(n) - kE_on_f(n) * R(n) * E_f + kE_off * RE(n);  % d/dt(R)
    dxdt(N+n) = k_e2 * RE(n-1) - k_e2 * RE(n) + kE_on_f(n) * R(n) * E_f - kE_off * RE(n) - kL_on * RE(n) * L_f + kL_off * REL(n-PAS+1);  % d/dt(RE)
    dxdt(2*N+n-PAS+1) = k_e2 * REL(n-PAS) - k_e2 * REL(n-PAS+1) - k_c * REL(n-PAS+1) + kL_on * RE(n)* L_f- kL_off * REL(n-PAS+1); %d/dt(REL)
end

end
