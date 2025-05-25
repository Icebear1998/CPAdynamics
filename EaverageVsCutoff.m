global N PAS N_PAS;

% Model parameters
L_a = 100;
P.k_in = 2;
P.k_e = 65/L_a;
P.k_e2 = 30/L_a;
P.E_total = 70000;
P.L_total = 100000;
P.Pol_total = 70000;
kHon = 0.1;
P.kHoff = 0.0025;
P.kc = 0.8;
P.kPmin = 0.05;
P.kPmax = 40;
geneLength_bp = 25000;
PASposition = 20000;
N = floor(geneLength_bp / L_a);
PAS = floor(PASposition / L_a);
N_PAS = N - PAS + 1;
EBindingNumber = 3;

% Parameters to vary
avg_E_bound_PAS_values = 0.5:0.5:10;
cutoff_values = zeros(size(avg_E_bound_PAS_values));
REH_at_800bp = zeros(size(avg_E_bound_PAS_values));

for i = 1:length(avg_E_bound_PAS_values)
    P.kHon = kHon * avg_E_bound_PAS_values(i); % Adjust kHon to achieve desired avg_E_bound
    
    X0 = zeros(N + N_PAS, 1);
    X = fsolve(@(xx) ode_dynamics(xx, P), X0);
    
    R_sol = X(1:N);
    REH_sol = X(N+1 : N+N_PAS);
    
    cutoff = sum(REH_sol(1:9)) / (sum(R_sol(PAS:PAS+8)) + sum(REH_sol(1:9)));
    cutoff_values(i) = cutoff;
    REH_at_800bp(i) = REH_sol(9); % REH at PAS + 8 (800 bp after PAS)
end

% Plot
figure;
plot(avg_E_bound_PAS_values, cutoff_values, 'o-', 'LineWidth', 2, 'MarkerSize', 8);
hold on;
for i = 1:length(avg_E_bound_PAS_values)
    text(avg_E_bound_PAS_values(i), cutoff_values(i), sprintf('Sum REH=%.2f', REH_at_800bp(i)), ...
        'FontSize', 8, 'HorizontalAlignment', 'right', 'VerticalAlignment', 'bottom');
end
xlabel('Average E Bound at PAS');
ylabel('Fraction of Polymerases Processed by 800 bp');
title('CPA Cutoff vs Average E Bound at PAS');
grid on;
hold off;

%% ------------ ODE DYNAMICS FUNCTION ------------
function dxdt = ode_dynamics(X, P)
global N PAS

k_in   = P.k_in;
k_e    = P.k_e;
k_e2   = P.k_e2;
kHoff_t= P.kHoff;
kc_t   = P.kc;
kHon_t = P.kHon;

R   = X(1:N);
REH = X(N+1:end);

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