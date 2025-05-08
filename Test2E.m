% Define symbolic variables
syms R0 R1 R2 R3 R1E R2E R2E2 R3E R3E2 R3E3 R4 R4E R4E2 R4E3 R4E4 Kd Kon Koff Ef a

% Add assumptions to avoid trivial solutions
assume(R0 ~= 0);
assume(R1 ~= 0);
assume(R2 ~= 0);
assume(Kd ~= 0);
assume(Kon ~= 0);
assume(Koff ~= 0);
assume(Ef ~= 0);

% Define the equations
eq0 = R0 == a;
eq1 = Kd*R1 - R0 == 0;
eq2 = - Kd*R1 + R0 + Kd*R2 - R1 + Koff*R1E - Kon*R1*Ef == 0;
eq3 = Kon*Ef*R1 - Koff*R1E + Kd*R2E - R1E == 0;
eq4 = R1 + Kd*R3 + Koff*R2E - (Kd + 1 + 2*Kon*Ef)*R2 == 0;
eq5 = R1E + 2*Kon*Ef*R2 + Kd*R3E + Koff*R2E2 - (Kd + Koff + 1 + Kon*Ef)*R2E == 0;
eq6 = Kon*Ef*R2E + Kd*R3E2 - (Koff + 1)*R2E2 == 0;

% Adjusted equations for R3 and its related species
eq7 = R2 + Koff*R3E + Kd*R4 - (Kd + 3*Kon*Ef)*R3 == 0; % Adjusted for R4
eq8 = 3*Kon*Ef*R3 + R2E + Koff*R3E2 + Kd*R4E - (Koff + Kd + 2*Kon*Ef)*R3E == 0; % Adjusted for R4E
eq9 = R2E2 + 2*Kon*Ef*R3E + Koff*R3E3 + Kd*R4E2 - (Kd + Koff + Kon*Ef)*R3E2 == 0; % Adjusted for R4E2
eq10 = Kon*Ef*R3E2 - Koff*R3E3 + Kd*R4E3 == 0; % Adjusted for R4E3

% Extended equations for R4 and its related species
eq11 = R3 + Koff*R4E - (Kd + 4*Kon*Ef)*R4 == 0; % Equation for R4
eq12 = 4*Kon*Ef*R4 + R3E + Koff*R4E2 - (Koff + Kd + 3*Kon*Ef)*R4E == 0; % Equation for R4E
eq13 = R3E2 + 3*Kon*Ef*R4E + Koff*R4E3 - (Kd + Koff + 2*Kon*Ef)*R4E2 == 0; % Equation for R4E2
eq14 = R3E3 + 2*Kon*Ef*R4E2 + Koff*R4E4 - (Kd + Koff + Kon*Ef)*R4E3 == 0; % Equation for R4E3
eq15 = Kon*Ef*R4E3 - Koff*R4E4 == 0; % Equation for R4E4


% Solve the equations
solutions = solve([eq0, eq1, eq2, eq3, eq4, eq5, eq6, eq7, eq8, eq9, eq10, eq11, eq12, eq13, eq14, eq15], ...
                  [R0, R1, R2, R3, R1E, R2E, R2E2, R3E, R3E2, R3E3, R4, R4E, R4E2, R4E3, R4E4]);

% Simplify the solutions
solutions = structfun(@simplify, solutions, 'UniformOutput', false);

% Define new symbolic variables based on the solutions
R0_t = solutions.R0 + solutions.R1 + solutions.R2 + solutions.R3;
RE_t = solutions.R1E + solutions.R2E + solutions.R3E;
RE2_t = solutions.R2E2 + solutions.R3E2;
RE3_t = solutions.R3E3;

% Simplify the new symbolic variables
R0_t_simplified = simplify(R0_t);
RE_t_simplified = simplify(RE_t);
RE2_t_simplified = simplify(RE2_t);
RE3_t_simplified = simplify(RE3_t);

% Substitute specific parameter values (e.g., Kon = 1, Koff = 2)
parameter_values = [a, Ef, Kd, Kon, Koff];  % List of parameters
specific_values = [1, 3,1, 1, 2];       % Corresponding values

R0_t_numeric = subs(R0_t_simplified, parameter_values, specific_values);
RE_t_numeric = subs(RE_t_simplified, parameter_values, specific_values);
RE2_t_numeric = subs(RE2_t_simplified, parameter_values, specific_values);
RE3_t_numeric = subs(RE3_t_simplified, parameter_values, specific_values);
R_t = R0_t_numeric + RE_t_numeric + RE2_t_numeric + RE3_t_numeric;

%Display the results
disp('Solution for R0:');
disp(solutions.R0);
disp('Solution for R1:');
disp(solutions.R1);
disp('Solution for R2:');
disp(solutions.R2);
disp('Solution for R3:');
disp(solutions.R3);
disp('Solution for R1E:');
disp(solutions.R1E);
disp('Solution for R2E:');
disp(solutions.R2E);
disp('Solution for R2E2:');
disp(solutions.R2E2);
disp('Solution for R3E:');
disp(solutions.R3E);
disp('Solution for R3E2:');
disp(solutions.R3E2);
disp('Solution for R3E3:');
disp(solutions.R3E3);

% disp('Solution for R:');
% disp(double(R0_t_numeric/R_t));
% disp('Solution for RE:');
% disp(double(RE_t_numeric/R_t));
% disp('Solution for RE2:');
% disp(double(RE2_t_numeric/R_t));
% disp('Solution for RE3:');
% disp(double(RE3_t_numeric/R_t));