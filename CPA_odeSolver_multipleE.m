clear all;
%global k_in k_c kE_off kL_on kL_off k_e k_e2 E_total L_total N PAS Pol_total L_a kE_on_f;
global N PAS N_PAS Pol_total;
timescale = 65;
% Define constants
L_a = 100;

P.k_in = 2/timescale;      % Define k_in
P.k_c = 0.8/timescale;       % Define k_c 
kP_on_min = 0.001/timescale;     % Define phosphorylation min value
kP_on_max = 0.1/timescale;       % Define phosphorylation max value  
P.kE_on = 0.0001/timescale;
P.kE_off = 10/timescale;    % Define kE_off
P.kL_on = 0.00025/timescale;     % Define kL_on
P.kL_off = 0.00001/timescale;    % Define kL_off
P.kH_on = 0.2/timescale;       % Define RE + Hexamer binding rate   
P.kH_off = 0.1/timescale;     % Define RE + Hex unbinding rate
P.k_e = 65/timescale/L_a;      % Define k_e
P.k_e2 = 30/timescale/L_a;      % Define k_e2
P.E_total = 70000;   % Define E_total
P.L_total = 100000;   % Define L_total

Pol_total = 70000;   % Pol II total number
N = floor(25000/L_a);        % Gene length
PAS = floor(20000/L_a);      % PAS site position
N_PAS = N - PAS +1;

EBindingNumber = 3; % Number of possible E binding

P.EBindingDisAtPas = compute_normalized_ratios(kP_on_max, P.kE_on, P.kE_off, P.E_total, Pol_total); %calculate the equlibrium distribution of E binding

% Time span for the simulation
tspan = [0 1000000];

% Solve the ODE system
X0 = zeros(2*N + N_PAS + 2*EBindingNumber*N_PAS,1); % Number of equations
[t, X] = ode45(@(t, x) ode_system_multipleE(t, x, P), tspan, X0);

% Extract the solutions for R(l), RE(l), REL(l), E_f, and L_f
R_sol = X(:, 1:N);
RE_sol = X(:, N+1:2*N);
RE1_sol = X(:,2*N+1: 2*N+N_PAS);
RE2_sol = X(:,2*N+N_PAS+1: 2*N+2*N_PAS); 
RE3_sol = X(:,2*N+2*N_PAS+1: 2*N+3*N_PAS);
RE1H_sol = X(:,2*N+3*N_PAS+1: 2*N+4*N_PAS); 
RE2H_sol = X(:,2*N+4*N_PAS+1: 2*N+5*N_PAS);
RE3H_sol = X(:,2*N+5*N_PAS+1: 2*N+6*N_PAS);
REHL_sol = X(:,2*N+6*N_PAS+1: 2*N+7*N_PAS);

% Define l values corresponding to the range of gene length
l_values = (1-PAS):(N-PAS);
% Plot R(l), RE(l), REL(l) and R(l) + RE(l) + REL(l) vs l at the final time step
figure;
hold on;
plot(100*l_values, 100*R_sol(end, :), 'b-','LineWidth',2.5, 'DisplayName', 'R');
plot(100*l_values, 100*RE_sol(end, :), 'r-','LineWidth',2.5, 'DisplayName', 'RE');
plot(100*l_values, 100*[zeros(PAS-1,1);RE1_sol(end, :)'],'LineWidth',2.5, 'DisplayName', 'RE1');
plot(100*l_values, 100*[zeros(PAS-1,1);RE2_sol(end, :)'],'LineWidth',2.5, 'DisplayName', 'RE2');
plot(100*l_values, 100*[zeros(PAS-1,1);RE3_sol(end, :)'],'LineWidth',2.5, 'DisplayName', 'RE3');
plot(100*l_values, 100*[zeros(PAS-1,1);(RE1H_sol(end, :)+RE2H_sol(end, :)+RE3H_sol(end, :))'],'LineWidth',2.5, 'DisplayName', 'REH');
plot(100*l_values, 100*[zeros(PAS-1,1); REHL_sol(end, :)'], 'g-','LineWidth',3, 'DisplayName', 'REHL');
%plot(l_values, R_sol(end, :) + RE_sol(end, :) + [zeros(PAS-1,1)'; REL_sol(end, :)]), 'k--', 'DisplayName', 'R(l) + RE(l) + REL(l)');

xlabel('Distance from PAS (Bp)', 'FontSize', 14);
ylabel('Average Number of Copy',  'FontSize', 14);
legend('show', 'Location', 'northwest');
title('Plot of R(l), RE(l), REL(l) and R(l)');
hold off;

% E_f = P.E_total - sum(RE_sol,2);% - sum(RE1_sol,2) - sum(RE2_sol,2) - sum(RE3_sol,2) - sum(RE1H_sol,2)-sum(RE2H_sol,2)-sum(RE3H_sol,2)-sum(REHL_sol,2);
% L_f = P.L_total - sum(REHL_sol,2);
% Pol_f = Pol_total - sum(R_sol,2) - sum(RE_sol,2);% - sum(RE1_sol,2) - sum(RE2_sol,2) - sum(RE3_sol,2) - sum(RE1H_sol,2)-sum(RE2H_sol,2)-sum(RE3H_sol,2)-sum(REHL_sol,2);
% 
% figure; plot(t, Pol_f, t, E_f, t, L_f);
% xlabel('time');
% ylabel('Free conc')
% legend({'Pol';'E';'L'})