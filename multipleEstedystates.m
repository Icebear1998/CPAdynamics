% Parameters
IF = 10;
K_on_D = 0.2;
K_off_D = 0.1;
k_on_E = 0.05;      % Association rate for binding E
k_off_E = 0.02;     % Dissociation rate for E-bound states

% Define the symbolic variables for each state
syms R0 R1 R2 R3 E R1E R2E R3E R2E2 R3E2 R3E3

% Define the system of equations based on the reactions
eq0 = R0 == IF;
eq1 = R0 + K_off_D * R1 - K_on_D * R0 == 0;
eq2 = - K_off_D * R1 + K_on_D * R0 + K_off_D * R2 - K_on_D * R1 - k_on_E * R1 * E - k_off_E * R1E == 0;
eq3 = - K_off_D * R2 + K_on_D * R1 + K_off_D * R3 - K_on_D * R2 - k_on_E * R2 * E - k_off_E * R2E == 0;
eq4 = - K_off_D * R3 + K_on_D * R2 - k_on_E * R3 * E - k_off_E * R3E == 0;
eq5 = k_on_E * R1 * E - k_off_E * R1E;
eq6 = k_on_E * R1E - KD * R2E == 0;
eq7 = 2 * k_on_E * R2 * E - k_off_E * R2E == 0;
eq8 = k_on_E * R2E - KD * R3E == 0;
eq9 = 3 * k_on_E * R3 * E - k_off_E * R3E == 0;
eq10 = 2 * k_on_E * R2E * E - k_off_E * R2E2 == 0;
eq11 = 2 * k_on_E * R3E * E - k_off_E * R3E2 == 0;
eq12 = k_on_E * R3E2 * E - k_off_E * R3E3 == 0;

% Solve the system of equations
sol = solve([eq1, eq2, eq3, eq4, eq5, eq6, eq7, eq8, eq9, eq10, eq11, eq12], ...
            [R0, R1, R2, R3, R0E, R1E, R2E, R3E, R2E2, R3E2, R3E3]);

% Display the solutions for each state
fprintf('Steady-state concentrations:\n');
fprintf('R0: %f\n', sol.R0);
fprintf('R1: %f\n', sol.R1);
fprintf('R2: %f\n', sol.R2);
fprintf('R3: %f\n', sol.R3);
fprintf('R0.E: %f\n', sol.R0E);
fprintf('R1.E: %f\n', sol.R1E);
fprintf('R2.E: %f\n', sol.R2E);
fprintf('R3.E: %f\n', sol.R3E);
fprintf('R2.E^2: %f\n', sol.R2E2);
fprintf('R3.E^2: %f\n', sol.R3E2);
fprintf('R3.E^3: %f\n', sol.R3E3);
