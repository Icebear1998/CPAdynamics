global N PAS N_PAS Ef_ss kHon_ss;
syms Ef real;

% ------------ MODEL PARAMETERS ------------
L_a = 100;
P.k_in    = 2;
P.kEon    = 0.00025;
P.kEoff   = 10;
P.k_e     = 65/L_a;
P.k_e2    = 30/L_a;
P.E_total = 70000;
P.L_total = 100000;
P.Pol_total = 70000;
P.kHon_max = 4;
P.kHoff   = 0.0001;
P.kc      = 0.2;
P.kPmin   = 0.5; % Adjusted
P.kPmax   = 10;  % Adjusted
P.kLon    = 0.001;
P.kLoff   = 0.01;

geneLength_bp = 25000;
PASposition   = 20000;
N      = floor(geneLength_bp / L_a);  % total nodes
PAS    = floor(PASposition   / L_a);  % node index of PAS
N_PAS  = N - PAS + 1;                 % number of nodes at/after PAS

% Define kP profile across all nodes
Kp_vals = zeros(1, N);
Kp_vals(1:PAS) = linspace(P.kPmax, P.kPmin, PAS); % Linear decrease to PAS
Kp_vals(PAS:230) = linspace(P.kPmin, 0.1, 230-PAS+1); % Decrease to peak Ser2P
Kp_vals(231:N) = linspace(0.1, 0.5, N-230); % Increase to gene end

EBindingNumber = 3; 
[r_E_BeforePas, r_k_AfterPas, r_states] = compute_steady_states(P, EBindingNumber + 1); 
disp('done compute steady states');

% Compute r_states for all nodes
num_states = 10; % Number of states (4x4 grid, but some are zero due to constraints)
r_states_all = sym(zeros(num_states, N));
for i = 1:N
    P_temp = P;
    P_temp.kPmax = Kp_vals(i); % Update kP for this node
    [~, ~, r_states_node] = compute_steady_states(P_temp, EBindingNumber + 1);
    r_states_all(:, i) = r_states_node;
end
disp('done compute r_states_all');

% Compute average Ser2P at each node
Ser2P_avg = zeros(1, N);
for i = 1:N
    r = double(r_states_all(:, i));
    % Ser2P: 0 for state 1, 1 for states 2-3, 2 for states 4-6, 3 for states 7-10
    Ser2P_avg(i) = (0 * r(1) + 1 * (r(2) + r(3)) + 2 * (r(4) + r(5) + r(6)) + 3 * (r(7) + r(8) + r(9) + r(10)));
end

% Normalize Ser2P to match experimental scale (max 3 Ser2P -> 1)
Ser2P_avg = Ser2P_avg / 3;

% Plot Ser2P
figure;
plot(1:N, Ser2P_avg, 'b-', 'LineWidth', 2);
hold on;
xline(1, '--k', 'Sense TSS');
xline(PAS, '--k', 'PAS');
xlabel('Node (Position along gene)');
ylabel('Normalized Ser2P Level');
title('Average Ser2P Phosphorylation along Gene');
grid on;

% ------------ SOLVE THE STEADY STATE ------------
X0 = zeros(N + 2*N_PAS, 1);
X = fsolve(@(xx) ode_dynamics(xx, P), X0, optimoptions('fsolve', 'Display', 'iter'));
disp('done compute fsolve');

% Extract solutions
R_sol   = X(1:N);
REH_sol = X(N+1 : N+N_PAS);
REHL_sol = X(N+N_PAS+1 : N+2*N_PAS);

%% ------------ ODE DYNAMICS FUNCTION ------------
function dxdt = ode_dynamics(X, P)
global N PAS Ef_ss kHon_ss

k_in   = P.k_in;
k_e    = P.k_e;
k_e2   = P.k_e2;
kHoff  = P.kHoff;
kc     = P.kc;
kLon   = P.kLon;
kLoff  = P.kLoff;
RE_val_bind_E = P.RE_val_bind_E;
kHon_t_Ef  = P.kHon_t_Ef;

R   = X(1:N);
REH = X(N+1 : N+N_PAS);
REHL = X(N+N_PAS+1 : N+2*N_PAS);

E_f = double(P.E_total);
E_used_sym = sum(R(1:PAS)' .* RE_val_bind_E);
E_used = double(subs(E_used_sym, 'Ef', E_f));
E_f = P.E_total - E_used;
Ef_ss = E_f;

if E_f < 0
    error('Negative E_f (E_f = %g). Stopping simulation.', E_f);
end

Pol_f = P.Pol_total - sum(R) - sum(REH) - sum(REHL);
L_f = P.L_total - sum(REHL);

if Pol_f < 0 || L_f < 0
    error('Negative Pol_f (%g) or L_f (%g).', Pol_f, L_f);
end

kHon_t = kHon_t_Ef(E_f);
kHon_ss = kHon_t;
dxdt = zeros(length(X), 1);

% Before PAS
n = 1;
dxdt(n) = Pol_f*k_in - k_e*R(n);

for n = 2:(PAS-1)
    dxdt(n) = k_e*R(n-1) - k_e*R(n);
end

% At PAS
n = PAS;
j = n - PAS + 1;
dxdt(n) = k_e*R(n-1) - k_e*R(n) - kHon_t*R(n) + kHoff*REH(j);
dxdt(N+j) = -k_e2*REH(j) + kHon_t*R(n) - kHoff*REH(j) - kLon*L_f*REH(j) + kLoff*REHL(j);
dxdt(N+N_PAS+j) = kLon*L_f*REH(j) - kLoff*REHL(j) - kc*REHL(j);

% After PAS
for n = (PAS+1):N
    j = n - PAS + 1;
    dxdt(n) = k_e*R(n-1) - k_e*R(n) - kHon_t*R(n) + kHoff*REH(j) - kc*REHL(j);
    dxdt(N+j) = k_e2*REH(j-1) - k_e2*REH(j) + kHon_t*R(n) - kHoff*REH(j) - kLon*L_f*REH(j) + kLoff*REHL(j);
    dxdt(N+N_PAS+j) = kLon*L_f*REH(j) - kLoff*REHL(j) - kc*REHL(j);
end
end
