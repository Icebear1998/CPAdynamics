% Define symbolic variables
syms RE1 RE2 RE3 RE1H RE2H RE3H REHL k_e k_e2 kH_on kH_off kL_on L_f k_c;

% Add assumptions to avoid trivial solutions
assume(RE1 ~= 0);
assume(RE2 ~= 0);
assume(RE3 ~= 0);
assume(RE1H ~= 0);
assume(RE2H ~= 0);
assume(RE3H ~= 0);
assume(REHL ~= 0);
assume(k_e ~= 0);
assume(k_e2 ~= 0);
assume(kH_on ~= 0);
assume(kH_off ~= 0);
assume(kL_on ~= 0);
assume(L_f ~= 0);
assume(k_c ~= 0);

% Define the equations
eq1 = k_e * RE1 - k_e * RE2 - kH_on * RE2 + kH_off * RE1H; % d/dt(RE1)
eq2 = k_e * RE2 - k_e * RE3 - 2 * kH_on * RE3 + kH_off * RE2H; % d/dt(RE2)
eq3 = k_e * RE3 - 3 * kH_on * RE3 + kH_off * RE3H; % d/dt(RE3)

eq4 = k_e2 * RE1H - k_e2 * RE2H + kH_on * RE2 - kH_off * RE1H - kL_on * RE1H * L_f; % d/dt(RE1H)
eq5 = k_e2 * RE2H - k_e2 * RE3H + 2 * kH_on * RE3 - kH_off * RE2H - kL_on * RE2H * L_f; % d/dt(RE2H)
eq6 = k_e2 * RE3H + 3 * kH_on * RE3 - kH_off * RE3H - kL_on * RE3H * L_f; % d/dt(RE3H)

eq7 = k_e2 * REHL - k_c * REHL + kL_on * RE1H * L_f + kL_on * RE2H * L_f + kL_on * RE3H * L_f; % d/dt(REHL)

% Solve the equations
solutions = solve([eq1, eq2, eq3, eq4, eq5, eq6, eq7], [RE1, RE2, RE3, RE1H, RE2H, RE3H, REHL]);
solutions = structfun(@simplify, solutions, 'UniformOutput', false);

% Display the solutions
disp('Simplified solutions:');
disp(solutions);

% Generate normalized concentrations
RE_total = solutions.RE1 + solutions.RE2 + solutions.RE3 + ...
           solutions.RE1H + solutions.RE2H + solutions.RE3H + solutions.REHL;

RE1_norm = simplify(solutions.RE1 / RE_total);
RE2_norm = simplify(solutions.RE2 / RE_total);
RE3_norm = simplify(solutions.RE3 / RE_total);
RE1H_norm = simplify(solutions.RE1H / RE_total);
RE2H_norm = simplify(solutions.RE2H / RE_total);
RE3H_norm = simplify(solutions.RE3H / RE_total);
REHL_norm = simplify(solutions.REHL / RE_total);

% Display normalized concentrations
disp('Normalized concentrations:');
disp(struct('RE1_norm', RE1_norm, 'RE2_norm', RE2_norm, 'RE3_norm', RE3_norm, ...
    'RE1H_norm', RE1H_norm, 'RE2H_norm', RE2H_norm, 'RE3H_norm', RE3H_norm, 'REHL_norm', REHL_norm));

% Note: For specific parameter values, substitute them into the normalized solutions.
% Example:
% k_e_val = 1;
% kH_on_val = 0.1;
% ...
% RE1_val = double(subs(RE1_norm, [k_e, kH_on, ...], [k_e_val, kH_on_val, ...]));
