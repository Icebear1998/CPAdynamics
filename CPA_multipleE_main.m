global N PAS N_PAS Ef_ss FirstRun;
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
P.kHon = 0.05; % based on typical k bind and estimated J factor for H.
P.kHoff = 0.0025; 
P.kc = 0.8; %not sure
% P.kPmin   = 0.1; %not sure
% P.kPmax   = 40; %not sure
kPon_min = 0.01; % at TSS
kPon_max = 0.5; % at PAS
kPoff_min = 0.001; % at PAS
kPoff_max = 0.4; % at TSS

geneLength_bp = 25000;
PASposition   = 20000;
N      = floor(geneLength_bp / L_a);  % total nodes
PAS    = floor(PASposition   / L_a);  % node index of PAS
N_PAS  = N - PAS + 1;                 % number of nodes at/after PAS
Ef_ss = 0;


EBindingNumber = 3; 
[r_E_BeforePas, r_P] = compute_steady_states(P, EBindingNumber + 1); 
disp('done compute steady states');

% KpOver_vals = linspace(1/P.kPmax, 1/P.kPmin, PAS); % Range of Kp for kPon increases linearly
% Kp_vals = 1./KpOver_vals;
kPon_vals = linspace(kPon_min, kPon_max, PAS); % Range of Kp for kPon increases linearly 
kPoff_vals = linspace(kPoff_max, kPoff_min, PAS); % Range of Kp for kPoff decreases linearly 
RE_vals = sym(zeros(EBindingNumber, PAS));
P_vals = sym(zeros(EBindingNumber, PAS));

for e = 1:EBindingNumber+1
    for i = 1:length(kPon_vals)
        kPon_val = kPon_vals(i);
        %kPoff_val = kPoff_vals(i);
        RE_vals(e, i) = subs(r_E_BeforePas(e), {'kPon', 'kPoff'}, {kPon_val, kPoff_max});
        P_vals(e, i) = subs(r_P(e), {'kPon', 'kPoff'}, {kPon_val, kPoff_max});
    end
end
disp('done compute EBindingNumber');

P.RE_val_bind_E = matlabFunction(simplify(sum(sym(1:EBindingNumber)' .* RE_vals(2:end, :), 1)), 'Vars', {Ef});
disp('done compute RE_val_bind_E');
% ------------ SOLVE THE STEADY STATE ------------
FirstRun = true;
%tspan = [0 1000]; % Increased time span for better steady-state approximation
X0 = zeros(N + N_PAS, 1);
% [t, X_ode] = ode45(@(t, x) ode_backbone(t, x, P), tspan, X0);
% X_init = X_ode(end,:);

X = fsolve(@(xx) ode_dynamics(xx, P), X0);
disp('done compute fsolve');

disp(Ef_ss);
avg_E_bound = P.RE_val_bind_E(Ef_ss);
disp(avg_E_bound(end));

disp('Recalculate kHon');
% Recalculate kHon (calculate kHon_tt)
FirstRun = false;
P.kHon = P.kHon * avg_E_bound(end);
X = fsolve(@(xx) ode_dynamics(xx, P), X);

%Extract solutions
R_sol   = X(1:N);
REH_sol = X(N+1 : N+N_PAS);

avg_E_bound = zeros(1, PAS); % Row vector for positions 1 to PAS
for e = 1:EBindingNumber+1
    for i = 1:PAS
        RE_vals(e, i) = double(R_sol(i)*double(subs(RE_vals(e, i), {'Ef'}, {Ef_ss})));
        P_vals(e, i) = double(R_sol(i)*double(subs(P_vals(e, i), {'Ef'}, {Ef_ss})));
        % Compute the average number of E molecules bound at each position
        avg_E_bound(i) = avg_E_bound(i) + (e-1)*(RE_vals(e, i)/R_sol(i));
    end
end
disp('done compute RE_vals');
% Compute the average number Ser2P at each position
avg_P_bound = zeros(1, PAS); % Row vector for positions 1 to PAS
for i = 1:PAS
    total_P_bound = 0;
    total_P = 0;
    for e = 1:EBindingNumber+1
        num_E = e - 1; % Number of E molecules bound (0 to EBindingNumber)       
        P_e = double(P_vals(e, i)); % Amount of Pol II with num_E E molecules
        total_P_bound = total_P_bound + num_E * P_e;
        total_P = total_P + P_e;
    end
    if total_P > 0
        avg_P(i) = total_P_bound / total_P;
    else
        avg_P(i) = 0; % Avoid division by zero
    end
end

Ser2P = avg_P;
hold on;
plot((1-PAS):0, Ser2P, 'g-','LineWidth',2.5, 'DisplayName','Ser2P');
plot((1-PAS):0, avg_E_bound, 'b-','LineWidth',2.5, 'DisplayName','AverageE');
legend({'Ser2P', 'AverageE'}, 'Location', 'best');
xlabel('position'); ylabel('AverageE');
hold off;
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
global N PAS Ef_ss FirstRun

k_in   = P.k_in;
k_e    = P.k_e;
k_e2   = P.k_e2;
kHoff_t= P.kHoff;
kc_t   = P.kc;
RE_val_bind_E = P.RE_val_bind_E;
kHon_t = P.kHon;

R   = X(1:N);
REH = X(N+1:end);

if FirstRun
    % Set initial Ef_ss once (using P.E_total as a starting guess)
    if Ef_ss == 0
        Ef_ss = P.E_total; % Initial guess
    end

    % Convert symbolic expression to a numerical function
    E_used = sum(R(1:PAS)'.* RE_val_bind_E(Ef_ss));
    E_f = P.E_total - E_used;
    Ef_ss = E_f;
    %disp({sum(RE_val_bind_E(Ef_ss)), sum(R(1:PAS)), Ef_ss});

    % If E_f < 0, throw error and stop solver
    if Ef_ss < 0
        error('Negative E_f = %g. Stopping simulation.', E_f);
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
