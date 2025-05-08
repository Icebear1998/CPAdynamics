% Parameters
KD = 1;          % Dissociation constant for unbound states
k_on_E = 1;      % Association rate for binding E
k_off_E = 1;     % Dissociation rate for E-bound states

% Define the symbolic variables for each state
syms R0 R1 R2 R3 R0E R1E R2E R3E R2E2 R3E2 R3E3

% Define the system of equations based on the reactions

% Unbound state reactions
eq1 = KD * R1 - R0 == 0;
eq2 = KD * R2 - R1 == 0;
eq3 = KD * R3 - R2 == 0;

% Binding reactions with E
eq4 = k_on_E * R0 * E - k_off_E * R0E == 0;
eq5 = k_on_E * R1 * E - k_off_E * R1E == 0;
eq6 = k_on_E * R1E - KD * R2E == 0;
eq7 = 2 * k_on_E * R2 * E - k_off_E * R2E == 0;
eq8 = k_on_E * R2E - KD * R3E == 0;
eq9 = 3 * k_on_E * R3 * E - k_off_E * R3E == 0;

% Further binding reactions for bound states
eq10 = 2 * k_on_E * R2E * E - k_off_E * R2E2 == 0;
eq11 = 2 * k_on_E * R3E * E - k_off_E * R3E2 == 0;
eq12 = k_on_E * R3E2 * E - k_off_E * R3E3 == 0;

% Solve the system of equations
sol = solve([eq1, eq2, eq3, eq4, eq5, eq6, eq7, eq8, eq9, eq10, eq11, eq12], ...
            [R0, R1, R2, R3, R0E, R1E, R2E, R3E, R2E2, R3E2, R3E3]);

% Display the solutions for each state
R0_sol = double(sol.R0);
R1_sol = double(sol.R1);
R2_sol = double(sol.R2);
R3_sol = double(sol.R3);
R0E_sol = double(sol.R0E);
R1E_sol = double(sol.R1E);
R2E_sol = double(sol.R2E);
R3E_sol = double(sol.R3E);
R2E2_sol = double(sol.R2E2);
R3E2_sol = double(sol.R3E2);
R3E3_sol = double(sol.R3E3);

fprintf('Steady-state concentrations:\n');
fprintf('R0: %f\n', R0_sol);
fprintf('R1: %f\n', R1_sol);
fprintf('R2: %f\n', R2_sol);
fprintf('R3: %f\n', R3_sol);
fprintf('R0.E: %f\n', R0E_sol);
fprintf('R1.E: %f\n', R1E_sol);
fprintf('R2.E: %f\n', R2E_sol);
fprintf('R3.E: %f\n', R3E_sol);
fprintf('R2.E^2: %f\n', R2E2_sol);
fprintf('R3.E^2: %f\n', R3E2_sol);
fprintf('R3.E^3: %f\n', R3E3_sol);
