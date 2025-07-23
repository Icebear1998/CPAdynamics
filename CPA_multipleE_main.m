global N PAS N_PAS N_PAUSE Ef_ss FirstRun;
syms Ef real;
% ------------ MODEL PARAMETERS ------------
L_a = 100;
P.k_in    = 2;
P.kEon    = 0.00025;
P.kEoff   = 10;
P.k_e     = 65/L_a;
P.k_e2    = 30/L_a;
P.E_total = 10000;
P.L_total = 100000;
P.Pol_total = 70000;
P.kHon = 0.05; % based on typical k bind and estimated J factor for H.
P.kHoff = 0.0025; 
P.kc = 0.8; %not sure

kPon_min = 0.01; % at TSS
kPon_max = 1; % at PAS
kPoff_min = 0.1; % at PAS
kPoff_max = 2; % at TSS
kPoff_const = 1;
kPon_const = 1;

geneLength_bp = 25000;
PASposition   = 20000;
N      = floor(geneLength_bp / L_a);  % total nodes
PAS    = floor(PASposition   / L_a);  % node index of PAS
N_PAS  = N - PAS + 1;                 % number of nodes at/after PAS
PAUSE_LENGTH = 20; % length of pause region in number of nodes
N_PAUSE = PAS+PAUSE_LENGTH;
Ef_ss = 0;
P.kHon_afterPAS = P.kHon*ones(1,N_PAS);

EBindingNumber = 3; 
[r_E_BeforePas, r_P] = compute_steady_states(P, EBindingNumber + 1); 
disp('done compute steady states');


kPon_vals = linspace(kPon_min, kPon_max, PAS); % Range of Kp for kPon increases linearly 
%kPoff_vals = linspace(kPoff_min, kPoff_min, PAUSE_LENGTH); % Range of Kp for kPoff decreases linearly 

RE_vals = sym(zeros(EBindingNumber+1, N_PAUSE));
P_vals = sym(zeros(EBindingNumber+1, N_PAUSE));

for e = 1:EBindingNumber+1
    for i = 1:PAS
        kPon_val = kPon_vals(i);
        %kPoff_val = kPval_vals(i);
        RE_vals(e, i) = subs(r_E_BeforePas(e), {'kPon', 'kPoff'}, {kPon_val, kPoff_const});
        P_vals(e, i) = subs(r_P(e), {'kPon', 'kPoff'}, {kPon_val, kPoff_const});
    end
    for i = PAS+1:N
        %kPoff_val = kPoff_vals(i-PAS);
        RE_vals(e, i) = subs(r_E_BeforePas(e), {'kPon', 'kPoff'}, {kPon_max, kPoff_min});
        P_vals(e, i) = subs(r_P(e), {'kPon', 'kPoff'}, {kPon_max, kPoff_min});
    end
%     for i = N_PAUSE+1:N
%         RE_vals(e, i) = subs(r_E_BeforePas(e), {'kPon', 'kPoff'}, {kPon_max, kPoff_min});
%         P_vals(e, i) = subs(r_P(e), {'kPon', 'kPoff'}, {kPon_max, kPoff_min});
%     end
end
disp('done compute EBindingNumber');

P.RE_val_bind_E = matlabFunction(simplify(sum(sym(1:EBindingNumber)' .* RE_vals(2:end, :), 1)), 'Vars', {Ef});
disp('done compute RE_val_bind_E');

% ------------ SOLVE THE STEADY STATE ------------ %
FirstRun = true;
X0 = zeros(N + N_PAS, 1);
X = fsolve(@(xx) ode_dynamics(xx, P), X0);
disp('done compute fsolve');

disp(Ef_ss);
avg_E_bound = P.RE_val_bind_E(Ef_ss);
disp(avg_E_bound(PAS));

disp('Recalculate kHon');
% Recalculate kHon (calculate kHon_tt)
FirstRun = false;
P.kHon_afterPAS = P.kHon * avg_E_bound(PAS:end);
X = fsolve(@(xx) ode_dynamics(xx, P), X);

%Extract solutions
R_sol   = X(1:N);
REH_sol = X(N+1 : N+N_PAS);


for i = 1:N
    total_P_bound = 0;
    total_P = 0;
    for e = 1:EBindingNumber+1
        RE_vals(e, i) = double(R_sol(i)*double(subs(RE_vals(e, i), {'Ef'}, {Ef_ss})));
        P_e = double(R_sol(i)*double(subs(P_vals(e, i), {'Ef'}, {Ef_ss}))); % Amount of Pol II with num_E E molecules
        total_P_bound = total_P_bound + (e - 1) * P_e;
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
plot((1-PAS):N_PAS-1, Ser2P, 'g-','LineWidth',2.5, 'DisplayName','Ser2P');
plot((1-PAS):N_PAS-1, avg_E_bound, 'b-','LineWidth',2.5, 'DisplayName','AverageE');
legend({'Ser2P', 'AverageE'}, 'Location', 'best');
xlabel('position'); ylabel('AverageE');
hold off;
% ------------ PLOT RESULTS ------------
% 1. Time evolution plot
l_values =  (1-PAS):(N-PAS);

figure; hold on;
for e = 1:EBindingNumber+1
    plot((1-PAS):N_PAS-1, RE_vals(e,:), 'LineWidth',2);
end

plot(l_values, R_sol, 'b-','LineWidth',2.5, 'DisplayName','R(l)');
plot(l_values, [zeros(PAS-1,1);REH_sol], 'r-','LineWidth',2.5, 'DisplayName','REH(l)');
xlabel('Time'); ylabel('Total Pol II');
legend({'Total R', 'Total REH'}, 'Location', 'best');
title('Time Evolution of R and REH');



%% ------------ ODE DYNAMICS FUNCTION ------------
function dxdt = ode_dynamics(X, P)
global N PAS Ef_ss N_PAUSE FirstRun

k_in   = P.k_in;
k_e    = P.k_e;
k_e2   = P.k_e2;
kHoff_t= P.kHoff;
kc_t   = P.kc;
RE_val_bind_E = P.RE_val_bind_E;
kHon_t = P.kHon_afterPAS;

R   = X(1:N);
REH = X(N+1:end);

if FirstRun
    if Ef_ss == 0
        Ef_ss = P.E_total; % Initial guess for free E
    end

    % Calculate E_used based on current Ef_ss 
    %REvalbindE = RE_val_bind_E(Ef_ss);
    E_used = sum(R(1:N)' .* RE_val_bind_E(Ef_ss)) + sum(REH, 1);
    
    % Compute E_f as the difference between total and used E
    E_f = abs(P.E_total - E_used); % Ensure non-negative E_f
    
    % Update Ef_ss with the new E_f value
    Ef_ss = E_f;

    % Check for unphysical conditions
    if Ef_ss < 0
        error('Negative E_f = %g. Stopping simulation.', E_f);
        %Ef_ss = 0.01;
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
dxdt(n) = k_e*R(n-1) - k_e*R(n) - kHon_t(1)*R(n) + kHoff_t*REH(j);
dxdt(N+j) = -k_e2*REH(j) + kHon_t(1)*R(n) - kHoff_t*REH(j);

for n = (PAS+1):N
    j = n - PAS + 1;
    dxdt(n) = k_e*R(n-1) - k_e*R(n) - kHon_t(n-PAS+1)*R(n) + kHoff_t*REH(j);
    dxdt(N+j) = k_e2*REH(j-1) - k_e2*REH(j) + kHon_t(n-PAS+1)*R(n) - kHoff_t*REH(j) - kc_t*REH(j);
end

end
