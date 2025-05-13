clear all;
global N PAS Pol_total;
timescale = 65;
% Define constants
L_a = 100;

% Rate parameters scaled by the timescale
P.k_in    = 2/timescale;           % Initiation rate
P.k_c     = 0.2/timescale;         % Cleavage rate
kE_on_min = 0.00001/timescale;     % Minimum E binding rate
kE_on_max = 0.00025/timescale;     % Maximum E binding rate
P.kE_off  = 10/timescale;          % E unbinding rate
P.kL_on   = 0.00025/timescale;     % L binding rate
P.kL_off  = 0.001/timescale;       % L unbinding rate
P.k_e     = 65/timescale/L_a;      % Elongation rate before PAS
P.k_e2    = 30/timescale/L_a;      % Elongation rate after PAS
P.E_total = 70000;                 % Total E molecules
P.L_total = 100000;                % Total L molecules


N = floor(25000/L_a);        % Gene length
PAS = floor(20000/L_a);      % PAS site position
Pol_total = 70000;   % Pol II total number 
% lp = 1; % Persistence length of RNA (nm)

% Function for kE_on as a function of time
P.kE_on_f = @(n) (n<PAS)*((kE_on_max-kE_on_min)/PAS*n+kE_on_min)+(n>=PAS)*kE_on_max; % Example: linear increase with time
P.j_factor = @(n) ((4*lp/((10^4)*(n*100)))^1.5)*exp(-460*(lp^2)/(6.25*((n*100)^2)))*((1.25*10^5)/(lp^3)); %J-factor function for RNA looping

% Initial guesses for R(l), RE(l), REL(l), E_f, and L_f
% R0 = [(Pol/(2*PAS))*linspace(0, 1, PAS)'; zeros(N-PAS+1,1)];    % Initial guess for R(l)
% RE0 = [((Pol*(3/8))/PAS)*linspace(0, 1, PAS+1)'; zeros(N-PAS, 1)];   % Initial guess for RE(l)
% REL0 = [0; (Pol*(1/8)/(N-PAS+1))*ones(N-PAS+1, 1)];  % Initial guess for REL(l)

% Time span for the simulation
tspan = [0 1000000];

% Solve the ODE system
X0 = zeros(3*N-PAS+1,1);
[t, X] = ode45(@(t, x) ode_system(t, x, P), tspan, X0);

% Extract the solutions for R(l), RE(l), REL(l), E_f, and L_f
R_sol = X(:, 1:N);
RE_sol = X(:, N+1:2*N);
REL_sol = X(:, 2*N+1:end);

% Define l values corresponding to the range of gene length
l_values = (1-PAS):(N-PAS);

% Plot R(l), RE(l), REL(l) and R(l) + RE(l) + REL(l) vs l at the final time step
figure;
hold on;
plot(100*l_values, 100*R_sol(end, :), 'b-','LineWidth',2.5, 'DisplayName', 'R(l)');
plot(100*l_values, 100*RE_sol(end, :), 'r-','LineWidth',2.5, 'DisplayName', 'RE(l)');
plot(100*l_values, 100*[zeros(PAS-1,1); REL_sol(end, :)'], 'g-','LineWidth',3, 'DisplayName', 'REL(l)');
%plot(l_values, R_sol(end, :) + RE_sol(end, :) + [zeros(PAS-1,1)'; REL_sol(end, :)]), 'k--', 'DisplayName', 'R(l) + RE(l) + REL(l)');

xlabel('Distance from PAS (Bp)', 'FontSize', 14);
ylabel('Average Number of Copy',  'FontSize', 14);
legend('show', 'Location', 'northwest');
title('Plot of R(l), RE(l), REL(l) and R(l)');
hold off;

E_f = P.E_total - sum(RE_sol,2) - sum(REL_sol,2);
L_f = P.L_total - sum(REL_sol,2);
Pol_f = Pol_total - sum(R_sol,2) - sum(RE_sol,2) - sum(REL_sol,2);

figure; plot(t, Pol_f, t, E_f, t, L_f);
xlabel('time');
ylabel('Free conc')
legend({'Pol';'E';'L'})