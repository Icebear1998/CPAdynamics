% Define symbolic variables
syms R0 R1 R1E R2 R2E R2E2 Kp Kon Koff Ef a

% Add assumptions to avoid trivial solutions
assume(R0 ~= 0);
assume(R1 ~= 0);
assume(Kp ~= 0);
assume(Kon ~= 0);
assume(Koff ~= 0);
assume(Ef ~= 0);

% Define the equations
% Equations for two Es:
eq1 = Kp*R1 - R0 == 0;
eq2 = R0 + Koff*R1E + Kp*R2 - (Kon*Ef + Kp + 1)*R1 == 0;
eq3 = Kon*Ef*R1 + Kp*R2E - (Koff + 1)*R1E == 0;
eq4 = R1 + Koff*R2E - (2*Kon*Ef + Kp)*R2 == 0;
eq5 = 2*Kon*Ef*R2 + R1E + 2*Koff*R2E2 - (Kon*Ef + Kp + Koff)*R2E == 0;
eq6 = Kon*Ef*R2E - 2*Koff*R2E2 == 0;

% % Equations for 1 E:
% eq1 = Kp*R1 - R0 == 0;
% eq2 = R0 + Koff*R1E - (Kon*Ef + Kp)*R1 == 0;
% eq3 = Kon*Ef*R1 - Koff*R1E == 0;

% Solve the equations
solutions = solve([eq1, eq2, eq3, eq4, eq5, eq6], [R0, R1, R1E, R2, R2E, R2E2]);
% Simplify the solutions
solutions = structfun(@simplify, solutions, 'UniformOutput', false);

% Define new symbolic variables based on the solutions
R0_t = solutions.R0 + solutions.R1 + solutions.R2;
RE_t = solutions.R1E + solutions.R2E;
RE2_t = solutions.R2E2;
R_t = R0_t + RE_t + RE2_t;

% Simplify the new symbolic variables
R0_t_simplified = simplify(R0_t);
RE_t_simplified = simplify(RE_t);
RE2_t_simplified = simplify(RE2_t);

R0_norm = simplify(R0_t / R_t);
RE_norm = simplify(RE_t / R_t);
RE2_norm = simplify(RE2_t / R_t);

% % Substitute specific parameter values (e.g., Kon = 1, Koff = 2)
% parameter_values = [Ef, Kon, Koff, Kp];  % List of parameters
% specific_values = [1, 1, 2, 4];       % Corresponding values
% 
% R0_t_numeric = subs(R0_t_simplified, parameter_values, specific_values);
% RE_t_numeric = subs(RE_t_simplified, parameter_values, specific_values);
% RE2_t_numeric = subs(RE2_t_simplified, parameter_values, specific_values);
% R_t = R0_t_numeric + RE_t_numeric + RE2_t_numeric;

% Parameter values
Ef_val = 1;
Kon_val = 1;
Koff_val = 1;

% Generate data for plotting
%Kp_vals = logspace(-5, 2, 100); % Range of Kp
Kp_vals = linspace(0.01, 10, 100); % Range of Kp
R0_vals = zeros(size(Kp_vals));
RE_vals = zeros(size(Kp_vals));
RE2_vals = zeros(size(Kp_vals));



% % Display the results
% disp('Solution for R0:');
% disp(solutions.R0);
% disp('Solution for R1:');
% disp(solutions.R1);
% disp('Solution for R1E:');
% disp(solutions.R1E);
% disp('Solution for R2:');
% disp(solutions.R2);
% disp('Solution for R2E:');
% disp(solutions.R2E);
% disp('Solution for R2E2:');
% disp(solutions.R2E2);

% disp('Solution for R0/Rt:');
% disp(double(R0_t_numeric/R_t));
% disp('Solution for RE/Rt:');
% disp(double(RE_t_numeric/R_t));
% disp('Solution for RE2/Rt:');
% disp(double(RE2_t_numeric/R_t));

for i = 1:length(Kp_vals)
    Kp_val = Kp_vals(i);
    R0_vals(i) = double(subs(R0_norm, [Ef, Kon, Koff, Kp], [Ef_val, Kon_val, Koff_val, Kp_val]));
    RE_vals(i) = double(subs(RE_norm, [Ef, Kon, Koff, Kp], [Ef_val, Kon_val, Koff_val, Kp_val]));
    RE2_vals(i) = double(subs(RE2_norm, [Ef, Kon, Koff, Kp], [Ef_val, Kon_val, Koff_val, Kp_val]));
end

% Plot the results
figure;
plot(Kp_vals, R0_vals, 'r', 'LineWidth', 1.5); hold on;
plot(Kp_vals, RE_vals, 'b', 'LineWidth', 1.5);
plot(Kp_vals, RE2_vals, 'g', 'LineWidth', 1.5);
hold off;

% Customize the plot
set(gca, 'XDir', 'reverse');
title_str = sprintf('Normalized Concentrations vs K_p\n(Ef = %.2f, Kon = %.2f, Koff = %.2f)', Ef_val, Kon_val, Koff_val);
title(title_str);
xlabel('K_p');
ylabel('Normalized Concentrations');
legend({'R_0/R_t', 'R_E/R_t', 'R_{E2}/R_t'}, 'Location', 'best');
grid on;
