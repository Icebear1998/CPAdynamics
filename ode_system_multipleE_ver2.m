function dxdt = ode_system_multipleE_ver2(t, X, P)
%global k_in k_c kE_off kL_on kL_off k_e k_e2 E_total L_total N PAS Pol_total kE_on;
global N PAS N_PAS Pol_total;

k_in = P.k_in;      % Define k_in
k_c = P.k_c;       % Define k_c
kE_off = P.kE_off;    % Define kE_off
kL_on = P.kL_on;     % Define kL_on
%kL_off = P.kL_off;    % Define kL_off
kH_on = P.kH_on;    % Define RE + Hexamer binding rate 
kH_off = P.kH_off;  % Define RE + Hex unbinding rate
k_e = P.k_e;      % Define k_e
k_e2 = P.k_e2;      % Define k_e2
E_total = P.E_total;   % Define E_total
L_total = P.L_total;   % Define L_total
kE_on = P.kE_on;
RTotalAtPas = P.RTotalAtPas;

% Extract variables from X
% R = X(1:PAS); % length of N
% RE = X(PAS+1:2*PAS); % length of PAS


RE1 = X(1: N_PAS); n1 = 1; % length of N - PAS;
RE2 = X(N_PAS+1: 2*N_PAS); n2 = N_PAS+1; % length of N - PAS
RE3 = X(2*N_PAS+1: 3*N_PAS); n3 = 2*N_PAS+1; % length of N - PAS
RE1H = X(3*N_PAS+1: 4*N_PAS); n4 = 3*N_PAS+1; % length of N - PAS
RE2H = X(4*N_PAS+1: 5*N_PAS); n5 = 4*N_PAS+1; % length of N - PAS
RE3H = X(5*N_PAS+1: 6*N_PAS); n6 = 5*N_PAS+1; % length of N - PAS
REHL = X(6*N_PAS+1: 7*N_PAS); n7 = 6*N_PAS+1; % length of N - PAS

% E_f = E_total - sum(RE); % E free
% L_f = L_total - sum(REHL); % L free
% Pol_f = Pol_total - sum(R) - sum(RE); % Pol free

%E_f = E_total - sum(RE);
L_f = L_total - sum(REHL);
%Pol_f = Pol_total - sum(R) - sum(RE);


dxdt = zeros(length(X), 1);  % Initialize the output vector for the ODEs


REDis = P.EBindingDisAtPas*RTotalAtPas; % E binding equlabrium distribution at PAS
RE1(1) = REDis(2);
RE2(1) = REDis(3);
RE3(1) = REDis(4);

n = PAS + 1;
dxdt(n1 + (n - 1 - PAS)) = 0; % d/dt(RE1)
dxdt(n2 + (n - 1 - PAS)) = 0; % d/dt(RE2)
dxdt(n3 + (n - 1 - PAS)) = 0; % d/dt(RE3)

dxdt(n4 + (n - 1 - PAS)) = 0; % d/dt(RE1H)
dxdt(n5 + (n - 1 - PAS)) = 0; % d/dt(RE2H)
dxdt(n6 + (n - 1 - PAS)) = 0; % d/dt(RE3H)

dxdt(n7 + (n - 1 - PAS)) = 0;
    
for n = PAS+2:N
    dxdt(n1 + (n - 1 - PAS)) = k_e * RE1(n-PAS-1) - k_e * RE1(n-PAS) - kH_on * RE1(n-PAS) + kH_off * RE1H(n-PAS); % d/dt(RE1)
    dxdt(n2 + (n - 1 - PAS)) = k_e * RE2(n-PAS-1) - k_e * RE2(n-PAS) - 2 * kH_on * RE2(n-PAS) + kH_off * RE2H(n-PAS); % d/dt(RE2)
    dxdt(n3 + (n - 1 - PAS)) = k_e * RE3(n-PAS-1) - k_e * RE3(n-PAS) - 3 * kH_on * RE3(n-PAS) + kH_off * RE3H(n-PAS); % d/dt(RE3)
    
    dxdt(n4 + (n - 1 - PAS)) =  k_e2 * RE1H(n-PAS-1) - k_e2 * RE1H(n-PAS) + kH_on * RE1(n-PAS) - kH_off * RE1H(n-PAS) - kL_on * RE1H(n-PAS) * L_f; % d/dt(RE1H)
    dxdt(n5 + (n - 1 - PAS)) =  k_e2 * RE2H(n-PAS-1) - k_e2 * RE2H(n-PAS) + 2 * kH_on * RE2(n-PAS) - kH_off * RE2H(n-PAS) - kL_on * RE2H(n-PAS) * L_f; % d/dt(RE2H)
    dxdt(n6 + (n - 1 - PAS)) =  k_e2 * RE3H(n-PAS-1) - k_e2 * RE3H(n-PAS) + 3 * kH_on * RE3(n-PAS) - kH_off * RE1H(n-PAS) - kL_on * RE3H(n-PAS) * L_f; % d/dt(RE3H)
    
    dxdt(n7 + (n - 1 - PAS)) = k_e2 * REHL(n-PAS-1) - k_e2 * REHL(n-PAS) + kL_on * RE1H(n-PAS) * L_f + kL_on * RE2H(n-PAS) * L_f + kL_on * RE3H(n-PAS) * L_f - k_c * REHL(n-PAS); % d/dt(REHL)
end

end
