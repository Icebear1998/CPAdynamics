% Define symbolic variables
syms R0 R1 R2 R3 R1E R2E R2E2 R3E R3E2 R3E3 Kp Kon Koff Ef a

%%
% Add assumptions to avoid trivial solutions
assume(R0 ~= 0);
assume(R1 ~= 0);
assume(R2 ~= 0);
assume(Kp ~= 0);
assume(Kon ~= 0);
assume(Koff ~= 0);
assume(Ef ~= 0);
%%
% Define the equations
eq1 = Kp*R1 - R0 == 0;
eq2 = R0 + Koff*R1E + Kp*R2 - (Kon*Ef + Kp + 1)*R1 == 0;
eq3 = Kon*Ef*R1 + Kp*R2E - (Koff + 1)*R1E == 0;
eq4 = R1 + Koff*R2E + Kp*R3 - (2*Kon*Ef + Kp + 1)*R2 == 0;
eq5 = 2*Kon*Ef*R2 + R1E + Koff*R2E2 + Kp*R3E - (Kon*Ef + Kp + Koff + 1)*R2E == 0;
eq6 = Kon*Ef*R2E + Kp*R3E2 - (Koff + 1)*R2E2 == 0;
eq7 = R2 + Koff*R3E - (3*Kon*Ef + Kp)*R3 == 0;
eq8 = 3*Kon*Ef*R3 + R2E + Koff*R3E2 - (2*Kon*Ef + Kp + Koff)*R3E == 0;
eq9 = 2*Kon*Ef*R3E + R2E2 + Koff*R3E3 - (Kon*Ef + Kp + Koff)*R3E2 == 0;
eq10 = Kon*Ef*R3E2 - Koff*R3E3 == 0;

% Solve the equations
solutions = solve([eq1, eq2, eq3, eq4, eq5, eq6, eq7, eq8, eq9, eq10], [R0, R1, R2, R3, R1E, R2E, R2E2, R3E, R3E2, R3E3]);
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
parameter_values = [Ef, Kon, Koff, Kp];  % List of parameters
specific_values = [1, 0.5, 2, 1];       % Corresponding values

R0_t_numeric = subs(R0_t_simplified, parameter_values, specific_values);
RE_t_numeric = subs(RE_t_simplified, parameter_values, specific_values);
RE2_t_numeric = subs(RE2_t_simplified, parameter_values, specific_values);
RE3_t_numeric = subs(RE3_t_simplified, parameter_values, specific_values);
R_t = R0_t_numeric + RE_t_numeric + RE2_t_numeric + RE3_t_numeric;

%%
% Display the results
% disp('Solution for R0:');
% disp(solutions.R0);
% disp('Solution for R1:');
% disp(solutions.R1);
% disp('Solution for R2:');
% disp(solutions.R2);
% disp('Solution for R3:');
% disp(solutions.R3);
% disp('Solution for R1E:');
% disp(solutions.R1E);
% disp('Solution for R2E:');
% disp(solutions.R2E);
% disp('Solution for R2E2:');
% disp(solutions.R2E2);
% disp('Solution for R3E:');
% disp(solutions.R3E);
% disp('Solution for R3E2:');
% disp(solutions.R3E2);
% disp('Solution for R3E3:');
% disp(solutions.R3E3);
%%

disp('Solution for R0/Rt:');
disp(double(R0_t_numeric/R_t));
disp('Solution for RE/Rt:');
disp(double(RE_t_numeric/R_t));
disp('Solution for RE2/Rt:');
disp(double(RE2_t_numeric/R_t));
disp('Solution for RE3/Rt:');
disp(double(RE3_t_numeric/R_t));