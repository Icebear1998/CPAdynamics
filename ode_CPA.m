% Define constants
global k_in k_c kE_on_initial kE_off kL_on kL_off k_e k_e2 E_total L_total N PAS Pol L_a kE_on_f;

k_in = 0.01;      % Define k_in
k_c = 0.001;       % Define k_c
kE_on_initial = 0.001;     % Define kE_on
kE_off = 0.05;    % Define kE_off
kL_on = 0;     % Define kL_on
kL_off = 0.001;    % Define kL_off
k_e = 0.1*10;      % Define k_e
k_e2 = 0.02;      % Define k_e2
E_total = 70000;   % Define E_total
L_total = 100000;   % Define L_total
L_a = 1;

N = 1000;        % Gene length
PAS = 800;      % PAS site position
Pol = 70000;   % Pol II total number 

% Function for kE_on as a function of time
kE_on_f = @(l) kE_on_initial + 0.001 * l; % Example: linear increase with time

% Initial guesses for R(l), RE(l), REL(l), E_f, and L_f
R0 = [(Pol/(2*PAS))*linspace(0, 1, PAS+1)'; zeros(N-PAS+1,1)];    % Initial guess for R(l)
RE0 = [((Pol*(3/8))/PAS)*linspace(0, 1, PAS+1)'; zeros(N-PAS, 1)];   % Initial guess for RE(l)
REL0 = [zeros(PAS, 1); (Pol*(1/8)/(N-PAS+1))*ones(N-PAS+1, 1)];  % Initial guess for REL(l)
E_f0 = 0.5*E_total;              % Initial guess for E_f
L_f0 = 0.5*L_total;              % Initial guess for L_f
X0 = [R0; RE0; REL0; E_f0; L_f0];  % Combine all initial guesses into a single vector

% Solve the system using fsolve
options = optimoptions('fsolve', 'Display', 'iter','Algorithm', 'trust-region');  % Optionally display iterations
[X_sol, fval, exitflag] = fsolve(@mySystem, X0, options);

% Extract the solutions for R(l), RE(l), REL(l), E_f, and L_f
R_sol = X_sol(1:length(R0));
RE_sol = X_sol(length(R0)+1:length(R0)+length(RE0));
REL_sol = [X_sol(length(R0)+length(RE0)+1:length(R0)+length(RE0)+length(REL0))];
E_f_sol = X_sol(end-1);
L_f_sol = X_sol(end);

disp('Solution for E_f:');
disp(E_f_sol);
disp('Solution for L_f:');
disp(L_f_sol);

% Define l values corresponding to the range of gene length
l_values = -PAS:(N-PAS);

% Plot R(l), RE(l), REL(l) and R(l) + RE(l) + REL(l) vs l
figure;
hold on;
plot(l_values, R_sol, 'b-', 'DisplayName', 'R(l)');
plot(l_values, RE_sol, 'r-', 'DisplayName', 'RE(l)');
plot(l_values, REL_sol, 'g-', 'DisplayName', 'REL(l)');
plot(l_values, R_sol + RE_sol + [REL_sol; zeros(length(R_sol) - length(REL_sol), 1)], 'k--', 'DisplayName', 'R(l) + RE(l) + REL(l)');

xlabel('l');
ylabel('Concentration');
legend('show', 'Location','northwest');
title('Plot of R(l), RE(l), REL(l) and R(l) + RE(l) + REL(l) vs l');
hold off;

% Define the system of equations as a function
function F = mySystem(X)
    global k_in k_c kE_on_initial kE_off kL_on kL_off k_e k_e2 E_total L_total N PAS Pol L_a kE_on_f;
    
    % Extract variables from X
    R = X(1:N+1);          
    RE = X(N+2:2*N+2);     
    REL = X(2*N+3:3*N+3);     
    E_f = X(end-1);          
    L_f = X(end);          

    F = zeros(length(X), 1);     % Initialize the output vector for all equations
    
    % Define the equation for l = 1 (special case for initial condition)
    F(1) = k_in - k_e * R(1) - kE_on_initial * R(1) * E_f + kE_off * RE(1);
    
    % Define the equations for l = 2 to N+1
    for l = 2:N+1
        if l <= PAS  
            kE_on = kE_on_f(l);
            k_e_current = k_e;
            kL_on = 0;
        else  % l >= 0, corresponding to indices 1002 to 2001
            k_e_current = k_e2;
            kL_on = 0.001;
            F(2*N+2 + l - PAS) = k_e2 * REL(l-1) - k_e2 * REL(l) - k_c * REL(l) - kL_off * REL(l) + kL_on * RE(l) * L_f;
        end
        
        % Index in MATLAB: l = 1 corresponds to l = -PAS
        F(l) = k_e_current * (R(l-1)/L_a) - k_e_current * (R(l)/L_a) - kE_on * R(l) * E_f + kE_off * RE(l);
        F(N+1 + l) = k_e_current * (RE(l-1)/L_a) - k_e_current * (RE(l)/L_a) + kE_on * R(l) * E_f - kE_off * RE(l) + kL_off * REL(min(l, length(REL))) - kL_on * RE(l) * L_f;
        
    end
    
    % Equations for E_f and L_f
    F(end-1) = E_f - (E_total - sum(RE) - sum(REL));
    F(end) = L_f - (L_total - sum(REL));
end
