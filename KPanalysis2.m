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
P.kHon = 0.1; % based on typical k bind and estimated J factor for H.
P.kHoff = 0.0025; 
P.kc = 0.8; %not sure
P.kPmin   = 1; %not sure
P.kPmax   = 40; %not sure

geneLength_bp = 25000;
PASposition   = 20000;
N      = floor(geneLength_bp / L_a);  % total nodes
PAS    = floor(PASposition   / L_a);  % node index of PAS
N_PAS  = N - PAS + 1;                 % number of nodes at/after PAS

EBindingNumber = 3; 
[r_E_BeforePas] = compute_steady_states(P, EBindingNumber + 1); 
disp('done compute steady states');

KpOver_vals = linspace(1/P.kPmax, 1/P.kPmin, PAS); % Range of Kp
Kp_vals = 1./KpOver_vals;
%Kp_vals = linspace(P.kPmax, P.kPmin, PAS); % Range of Kp
RE_vals = sym(zeros(EBindingNumber, PAS));

for e = 1:EBindingNumber+1
    for i = 1:length(Kp_vals)
        kP_val = Kp_vals(i);
        RE_vals(e, i) = subs(r_E_BeforePas(e), {'kP'}, {kP_val});
    end
end
disp('done compute EBindingNumber');

P.RE_val_bind_E = matlabFunction(simplify(sum(sym(1:EBindingNumber)' .* RE_vals(2:end, :), 1)), 'Vars', {Ef});
disp('done compute RE_val_bind_E');
% ------------ SOLVE THE STEADY STATE ------------
%tspan = [0 1000]; % Increased time span for better steady-state approximation
X0 = zeros(N + N_PAS, 1);
% [t, X_ode] = ode45(@(t, x) ode_backbone(t, x, P), tspan, X0);
% X_init = X_ode(end,:);

X = fsolve(@(xx) ode_dynamics(xx, P), X0);
disp('done compute fsolve');
% Extract solutions
R_sol   = X(1:N);
REH_sol = X(N+1 : N+N_PAS);

avg_E_bound = zeros(1, PAS); % Row vector for positions 1 to PAS
for e = 1:EBindingNumber+1
    for i = 1:PAS
        RE_vals(e, i) = double(R_sol(i)*double(subs(RE_vals(e, i), {'Ef'}, {Ef_ss})));
        % Compute the average number of E molecules bound at each position
        avg_E_bound(i) = avg_E_bound(i) + (e-1)*(RE_vals(e, i)/R_sol(i));
    end
end

%disp(P.RE_val_bind_E(Ef_ss));

disp('done compute RE_vals');
% ------------ PLOT RESULTS ------------
% 1. Time evolution plot
l_values =  (1-PAS):(N-PAS);

figure; hold on;
for e = 1:EBindingNumber+1
    plot((1-PAS):0, RE_vals(e,:), 'LineWidth',2);
end

plot(l_values, R_sol, 'b-','LineWidth',2.5, 'DisplayName','R(l)');
plot(l_values, [zeros(PAS-1,1);REH_sol], 'r-','LineWidth',2.5, 'DisplayName','REH(l)');
xlabel('Time'); ylabel('Total Pol II');
legend({'Total R', 'Total REH'}, 'Location', 'best');
title('Time Evolution of R and REH');



%% ------------ ODE DYNAMICS FUNCTION ------------
function dxdt = ode_dynamics(X, P)
global N PAS Ef_ss kHon_ss

k_in   = P.k_in;
k_e    = P.k_e;
k_e2   = P.k_e2;
kHoff_t= P.kHoff;
kc_t   = P.kc;
RE_val_bind_E = P.RE_val_bind_E;
kHon_t = P.kHon;

R   = X(1:N);
REH = X(N+1:end);
E_f = double(P.E_total);

% Convert symbolic expression to a numerical function

E_used_Ef = RE_val_bind_E(E_f);
E_used = sum(R(1:PAS)'.* E_used_Ef);
%disp(size(E_used));

E_f = P.E_total - E_used;
Ef_ss = E_f;

% If E_f < 0, throw error and stop solver
if Ef_ss < 0
    error('Negative E_f at t = %g (E_f = %g). Stopping simulation.', t, E_f);
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
