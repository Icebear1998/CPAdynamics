global N PAS N_PAS Ef_ss;
syms Ef real;
saveData = false;
% ------------ MODEL PARAMETERS ------------
L_a = 100;
P.k_in    = 2;
P.kEon    = 0.00025;
P.kEoff   = 10;
P.k_e     = 65/L_a;
P.k_e2    = 30/L_a;
P.E_total = 70000;
P.L_total = 100000;
P.Pol_total = 70000;
P.kHon = 0.2; % based on typical k bind and estimated J factor for H.
P.kHoff = 0.0125; 
P.kc = 0.05; %not sure

kPon_min = 0.01; % at TSS
kPon_max = 1; % at PAS
kPoff_min = 0.1; % at PAS
kPoff_max = 1; % at TSS
kPoff_const = 1;
kPon_const = 1;

geneLength_bp = 15000;
PASposition   = 10000;
N      = floor(geneLength_bp / L_a);  % total nodes
PAS    = floor(PASposition   / L_a);  % node index of PAS
N_PAS  = N - PAS + 1;                 % number of nodes at/after PAS
Ef_ss = 0;

EBindingNumber = 1; 
[r_E_BeforePas, r_P] = compute_steady_states(P, EBindingNumber + 1); 
disp('done compute steady states');


kPon_vals = linspace(kPon_min, kPon_max, 200); % Range of Kp for kPon increases linearly 
kPon_vals = kPon_vals(1:PAS);
%kPoff_vals = linspace(kPoff_min, kPoff_min, PAUSE_LENGTH); % Range of Kp for kPoff decreases linearly 

RE_vals = sym(zeros(EBindingNumber+1, N));
P_vals = sym(zeros(EBindingNumber+1, N));
avg_P = zeros(1,N);

for e = 1:EBindingNumber+1
    for i = 1:PAS
        kPon_val = kPon_vals(i);
        %kPoff_val = kPval_vals(i);
        RE_vals(e, i) = subs(r_E_BeforePas(e), {'kPon', 'kPoff'}, {kPon_val, kPoff_const});
        P_vals(e, i) = subs(r_P(e), {'kPon', 'kPoff'}, {kPon_val, kPoff_const});
    end
    for i = PAS+1:N
        kPon_val = kPon_vals(PAS);
        RE_vals(e, i) = subs(r_E_BeforePas(e), {'kPon', 'kPoff'}, {kPon_val, kPoff_const});
        P_vals(e, i) = subs(r_P(e), {'kPon', 'kPoff'}, {kPon_val, kPoff_const});
    end
end
disp('done compute EBindingNumber');

P.RE_val_bind_E = matlabFunction(simplify(sum(sym(1:EBindingNumber)' .* RE_vals(2:end, :), 1)), 'Vars', {Ef});
disp('done compute RE_val_bind_E');

% ------------ SOLVE THE STEADY STATE ------------ %
P.FirstRun = true;
P.is_unphysical = false; % Reset flag
X0 = 1e-6 * ones(N + N_PAS, 1); % Small positive initial guess
options = optimoptions('fsolve', 'Display', 'off', 'FunctionTolerance', 1e-8, 'MaxIterations', 1000);
X = fsolve(@(xx) ode_dynamics_multipleE(xx, P), X0, options);
if P.is_unphysical || any(isnan(X)) || any(isinf(X))
    error('Simulation failed: unphysical solution detected');
end
disp('done compute fsolve');

disp(Ef_ss);
avg_E_bound = P.RE_val_bind_E(Ef_ss);
disp(avg_E_bound(PAS));

disp('Recalculate kHon');
% Recalculate kHon (calculate kHon_tt)
P.FirstRun = false;
P.kHon = P.kHon * avg_E_bound(end);
X = fsolve(@(xx) ode_dynamics_multipleE(xx, P), X, options);

%Extract solutions
R_sol   = X(1:N);
REH_sol = X(N+1 : N+N_PAS);


for i = 1:N
    total_P_bound = 0;
    total_P = 0;
    for e = 1:EBindingNumber+1
        RE_vals(e, i) = double(R_sol(i)*double(subs(RE_vals(e, i), {'Ef'}, {Ef_ss})));
        P_e = double(R_sol(i)*double(subs(P_vals(e, i), {'Ef'}, {Ef_ss}))); % Amount of Pol II with num_E E molecules
        total_P_bound = total_P_bound + (e - 1) * P_e;
        total_P = total_P + P_e;
    end
    if total_P > 0
        avg_P(i) = total_P_bound / total_P;
    else
        avg_P(i) = 0; % Avoid division by zero
    end
end

Ser2P = avg_P;
hold on;
plot((1-PAS):N_PAS-1, Ser2P, 'g-','LineWidth',2.5, 'DisplayName','Ser2P');
plot((1-PAS):N_PAS-1, avg_E_bound, 'b-','LineWidth',2.5, 'DisplayName','AverageE');
legend({'Ser2P', 'AverageE'}, 'Location', 'best');
xlabel('position'); ylabel('AverageE');
hold off;
% ------------ PLOT RESULTS ------------
% 1. Time evolution plot
l_values =  (1-PAS):(N-PAS);

figure; hold on;
for e = 1:EBindingNumber+1
    plot((1-PAS):N_PAS-1, RE_vals(e,:), 'LineWidth',2);
end

plot(l_values, R_sol, 'b-','LineWidth',2.5, 'DisplayName','R(l)');
plot(l_values, [zeros(PAS-1,1);REH_sol], 'r-','LineWidth',2.5, 'DisplayName','REH(l)');
xlabel('Time'); ylabel('Total Pol II');
legend({'Total R', 'Total REH'}, 'Location', 'best');
title('Time Evolution of R and REH');
hold off;

% 1. Calculate the final, steady-state concentration of free Pol II
% This is the total Pol II minus all polymerases bound to the gene (R and REH).
bound_pol_II = sum(R_sol) + sum(REH_sol);
Pol_f_final = P.Pol_total - bound_pol_II;

% 2. Display the values in the command window for clarity
fprintf('\n--- Polymerase Distribution at Steady State ---\n');
fprintf('Total Pol II in system: %d\n', P.Pol_total);
fprintf('Total Bound Pol II (on gene): %.2f\n', bound_pol_II);
fprintf('Total Free Pol II (Pol_f):    %.2f\n', Pol_f_final);

if saveData
    % --- SAVE RESULTS ---
    % Prepare data structure for saving
    data.EBindingNumber = EBindingNumber;
    data.R_sol = R_sol;
    data.REH_sol = REH_sol;
    data.Ser2P = Ser2P;
    data.avg_E_bound = avg_E_bound;
    data.Ef_ss = Ef_ss;
    data.Pol_f_final = Pol_f_final;

    % Save results using the utility function
    save_analysis_results('CPA_multipleE_main', data, P);
end


