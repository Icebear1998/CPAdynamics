function result = compute_normalized_ratios(KpVal, KonVal, KoffVal, EfVal, PolVal)
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
    solutions = solve([eq1, eq2, eq3, eq4, eq5, eq6, eq7, eq8, eq9, eq10], ...
                      [R0, R1, R2, R3, R1E, R2E, R2E2, R3E, R3E2, R3E3]);
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

    % Parameter values
    Ef_val = EfVal/PolVal;
    Kon_val = KonVal;
    Koff_val = KoffVal;

%     % Generate data for plotting
%     % Kp_vals = linspace(kP_off_max, kP_off_min, 100); % Range of Kp
%     result = zeros(4, length(Kp_vals)); % 4x100 matrix to store results
% 
%     % Compute normalized values for each Kp and store in result matrix
%     for i = 1:length(Kp_vals)
%         Kp_val = Kp_vals(i);
%         result(1, i) = double(subs(R0_norm, [Ef, Kon, Koff, Kp], [Ef_val, Kon_val, Koff_val, Kp_val]));
%         result(2, i) = double(subs(RE_norm, [Ef, Kon, Koff, Kp], [Ef_val, Kon_val, Koff_val, Kp_val]));
%         result(3, i) = double(subs(RE2_norm, [Ef, Kon, Koff, Kp], [Ef_val, Kon_val, Koff_val, Kp_val]));
%         result(4, i) = double(subs(RE3_norm, [Ef, Kon, Koff, Kp], [Ef_val, Kon_val, Koff_val, Kp_val]));
%     end

    result = zeros(4, 1); % Store the result for Kp value
    Kp_val = KpVal;
    result(1, 1) = double(subs(R0_norm, [Ef, Kon, Koff, Kp], [Ef_val, Kon_val, Koff_val, Kp_val]));
    result(2, 1) = double(subs(RE_norm, [Ef, Kon, Koff, Kp], [Ef_val, Kon_val, Koff_val, Kp_val]));
    result(3, 1) = double(subs(RE2_norm, [Ef, Kon, Koff, Kp], [Ef_val, Kon_val, Koff_val, Kp_val]));
    result(4, 1) = double(subs(RE3_norm, [Ef, Kon, Koff, Kp], [Ef_val, Kon_val, Koff_val, Kp_val]));
end
