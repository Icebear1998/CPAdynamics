clear all;
global N PAS Pol_total;
% Define constants
L_a = 100;

% Rate parameters scaled by the timescale
P.k_in    = 2;           % Initiation rate
P.k_c     = 0.6;         % Cleavage rate
kE_on_min = 0.00001;     % Minimum E binding rate
kE_on_max = 0.00025;     % Maximum E binding rate
P.kE_off  = 10;          % E unbinding rate
P.kL_on   = 0.00025;     % L binding rate
P.kL_off  = 0.001;       % L unbinding rate
P.k_e     = 65/L_a;      % Elongation rate before PAS
P.k_e2    = 30/L_a;      % Elongation rate after PAS
P.E_total = 70000;                 % Total E molecules
P.L_total = 100000;                % Total L molecules


N = floor(25000/L_a);        % Gene length
PAS = floor(20000/L_a);      % PAS site position
Pol_total = 70000;   % Pol II total number 
% lp = 1; % Persistence length of RNA (nm)

% Function for kE_on as a function of time
P.kE_on_f = @(n) (n<PAS)*((kE_on_max-kE_on_min)/PAS*n+kE_on_min)+(n>=PAS)*kE_on_max; % Example: linear increase with time
%P.j_factor = @(n) ((4*lp/((10^4)*(n*100)))^1.5)*exp(-460*(lp^2)/(6.25*((n*100)^2)))*((1.25*10^5)/(lp^3)); %J-factor function for RNA looping

% Initial guesses for R(l), RE(l), REL(l), E_f, and L_f
% R0 = [(Pol/(2*PAS))*linspace(0, 1, PAS)'; zeros(N-PAS+1,1)];    % Initial guess for R(l)
% RE0 = [((Pol*(3/8))/PAS)*linspace(0, 1, PAS+1)'; zeros(N-PAS, 1)];   % Initial guess for RE(l)
% REL0 = [0; (Pol*(1/8)/(N-PAS+1))*ones(N-PAS+1, 1)];  % Initial guess for REL(l)

% Time span for the simulation
tspan = [0 10000];

% Solve the ODE system
X0 = zeros(3*N-PAS+1,1);
[t, X_ode] = ode45(@(t, x) ode_system(t, x, P), tspan, X0);
X_init = X_ode(end,:);

options = optimoptions('fsolve', 'Display', 'iter', 'MaxIterations', 1000, 'MaxFunctionEvaluations', 100000);
[X, fval, exitflag] = fsolve(@(x) ode_system(0, x, P), X_init, options);

% Check convergence
if exitflag > 0
    disp('fsolve converged successfully.');
else
    disp('fsolve did not converge. Check the initial guess or system equations.');
end

% Extract the solutions
R_sol   = X(1:N);
RE_sol  = X(N+1:2*N);
REL_sol = X(2*N+1:end);

plot_spatial_distribution(R_sol, RE_sol, REL_sol, N, PAS, 'Main Run: Final Spatial Distribution: R(l), RE(l), and REL(l)');

% Compute free species for mass balance
E_f_time = P.E_total - sum(RE_sol) - sum(REL_sol);
L_f_time = P.L_total - sum(REL_sol);
Pol_f_time = Pol_total - sum(R_sol) - sum(RE_sol) - sum(REL_sol);


function plot_spatial_distribution(R, RE, REL, N, PAS, plot_title)
    l_values = (1-PAS):(N-PAS); 
    REL_plot = [zeros(PAS-1,1); REL(:)];  % Force REL into a column vector, assume REL=0 for nodes < PAS

    figure;
    hold on;
    plot(l_values, R, 'b-', 'LineWidth', 2.5, 'DisplayName', 'R(l)');
    plot(l_values, RE, 'r-', 'LineWidth', 2.5, 'DisplayName', 'RE(l)');
    plot(l_values, REL_plot, 'g-', 'LineWidth', 3, 'DisplayName', 'REL(l)');
    xlabel('Distance from PAS (bp)', 'FontSize', 14);
    ylabel('Pol II (scaled)', 'FontSize', 14);
    legend('show', 'Location', 'northwest');
    title(plot_title);
    hold off;
end