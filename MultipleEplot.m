% Define symbolic variables
syms R0 R1 R2 R3 R1E R2E R2E2 R3E R3E2 R3E3 Kp Kon Koff Ef

% Add assumptions to avoid trivial solutions
assume(R0 ~= 0);
assume(R1 ~= 0);
assume(R2 ~= 0);
assume(Kp ~= 0);
assume(Kon ~= 0);
assume(Koff ~= 0);
assume(Ef ~= 0);

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
solutions = structfun(@simplify, solutions, 'UniformOutput', false);

% Define normalized variables
R0_t = solutions.R0 + solutions.R1 + solutions.R2 + solutions.R3;
RE_t = solutions.R1E + solutions.R2E + solutions.R3E;
RE2_t = solutions.R2E2 + solutions.R3E2;
RE3_t = solutions.R3E3;
R_t = R0_t + RE_t + RE2_t + RE3_t;

R0_norm = simplify(R0_t / R_t);
RE_norm = simplify(RE_t / R_t);
RE2_norm = simplify(RE2_t / R_t);
RE3_norm = simplify(RE3_t / R_t);

% % Display the results
% disp('Solution for R0/Rt:');
% disp(R0_norm);
% disp('Solution for R1/Rt:');
% disp(RE_norm);
% disp('Solution for R2/Rt:');
% disp(RE2_norm);
% disp('Solution for R3/Rt:');
% disp(RE3_norm);

% Parameter values
Ef_val = 1;
Kon_val = 1;
Koff_val = 1.5;

R0_vals = (subs(R0_norm, [Ef, Kon, Koff], [Ef_val, Kon_val, Koff_val]));
RE_vals = (subs(RE_norm, [Ef, Kon, Koff], [Ef_val, Kon_val, Koff_val]));
RE2_vals = (subs(RE2_norm, [Ef, Kon, Koff], [Ef_val, Kon_val, Koff_val]));
RE3_vals = (subs(RE3_norm, [Ef, Kon, Koff], [Ef_val, Kon_val, Koff_val]));

% Display the results
disp('Solution for R0/Rt:');
disp(R0_vals);
disp('Solution for R1/Rt:');
disp(RE_vals);
disp('Solution for R2/Rt:');
disp(RE2_vals);
disp('Solution for R3/Rt:');
disp(RE3_vals);

% Generate data for plotting
Kp_vals = linspace(0.1, 10, 100); % Range of Kp
R0_vals = zeros(size(Kp_vals));
RE_vals = zeros(size(Kp_vals));
RE2_vals = zeros(size(Kp_vals));
RE3_vals = zeros(size(Kp_vals));

for i = 1:length(Kp_vals)
    Kp_val = Kp_vals(i);
    R0_vals(i) = double(subs(R0_norm, [Ef, Kon, Koff, Kp], [Ef_val, Kon_val, Koff_val, Kp_val]));
    RE_vals(i) = double(subs(RE_norm, [Ef, Kon, Koff, Kp], [Ef_val, Kon_val, Koff_val, Kp_val]));
    RE2_vals(i) = double(subs(RE2_norm, [Ef, Kon, Koff, Kp], [Ef_val, Kon_val, Koff_val, Kp_val]));
    RE3_vals(i) = double(subs(RE3_norm, [Ef, Kon, Koff, Kp], [Ef_val, Kon_val, Koff_val, Kp_val]));
end

% Plot the results
figure;
plot(Kp_vals, R0_vals, 'r', 'LineWidth', 1.5); hold on;
plot(Kp_vals, RE_vals, 'b', 'LineWidth', 1.5);
plot(Kp_vals, RE2_vals, 'g', 'LineWidth', 1.5);
plot(Kp_vals, RE3_vals, 'k', 'LineWidth', 1.5);
hold off;

% Customize the plot
set(gca, 'XDir', 'reverse');
xlabel('K_p');
ylabel('Normalized Concentrations');
title('Normalized Concentrations vs K_p');
legend({'R_0/R_t', 'R_E/R_t', 'R_{E2}/R_t', 'R_{E3}/R_t'}, 'Location', 'best');
grid on;
