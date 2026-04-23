global N PAS N_PAS Ef_ss Lf_ss;
syms Ef real;
saveData = false;

% ------------ MODEL PARAMETERS ------------
P.L_a = 100;
P.k_in    = 2;
P.kEon    = 0.000002;
P.kEoff   = 1;
P.k_e     = 65/P.L_a;
P.k_e2    = 30/P.L_a;
P.E_total = 100000;
P.L_total = 150000;
P.Pol_total = 70000;

% Reaction 1: Tethered E-hexamer encounter (per E-factor)
P.kHon  = 2;        % s^{-1} per E-factor (from Appendix A, k1_single)
P.kHoff = 0.01;    % s^{-1} (hexamer dissociation)

% Reaction 2: L-factor diffusional recruitment (per E-factor per L-molecule)
P.kLon  = 4.4e-6;   % molecule^{-1} s^{-1} per E-factor (antenna effect)
P.kLoff = 0.001;        % s^{-1} (CstF-64 weak binding)

% Cleavage (only from fully assembled REHL state)
P.kc = 0.1;

P.kPon_min = 0.01;
P.kPon_slope = 0.007;
P.kPoff = 1;

P.geneLength_bp = 25000;
P.PASposition   = 20000;

N      = floor(P.geneLength_bp / P.L_a);
PAS    = floor(P.PASposition   / P.L_a);
N_PAS  = N - PAS + 1;
Ef_ss = 0;
Lf_ss = P.L_total;

EBindingNumber = 5;

% Run termination simulation (explicit L version)
[R_sol, REH_sol, REHL_sol, P, r_E_BeforePas, r_P] = run_termination_simulation(P, EBindingNumber);
disp('done compute simulation');

% Set up kPon values with linear increase for P_vals calculation
kPon_vals = P.kPon_min + P.kPon_slope * (0:N-1);

RE_vals = sym(zeros(EBindingNumber+1, N));
P_vals = sym(zeros(EBindingNumber+1, N));
avg_P = zeros(1,N);

for e = 1:EBindingNumber+1
    for i = 1:N
        kPon_val = kPon_vals(i);
        RE_vals(e, i) = subs(r_E_BeforePas(e), {'kPon', 'kPoff'}, {kPon_val, P.kPoff});
        P_vals(e, i) = subs(r_P(e), {'kPon', 'kPoff'}, {kPon_val, P.kPoff});
    end
end
disp('done compute P_vals for Ser2P');

% Get average E bound values
fprintf('E free: %d\n', Ef_ss);
fprintf('L free: %d\n', Lf_ss);
avg_E_bound = P.RE_val_bind_E(Ef_ss);

% Compute PAS usage profile using cleavage from REHL
[exit_cdf, distances_bp] = calculate_pas_usage_profile(R_sol, REH_sol, REHL_sol, P);
max_cdf = max(exit_cdf);
if max_cdf >= 0.5
    CAD = interp1(exit_cdf, distances_bp, 0.5, 'linear');
else
    CAD = Inf;
    fprintf('WARNING: CDF only reaches %.2f%% — 50%% cleavage not achieved within gene.\n', max_cdf*100);
    fprintf('  Max cleavage distance: %d bp (last node)\n', max(distances_bp));
end

% Compute Ser2P profile
for i = 1:N
    total_P_bound = 0;
    total_P = 0;
    for e = 1:EBindingNumber+1
        RE_vals(e, i) = double(R_sol(i)*double(subs(RE_vals(e, i), {'Ef'}, {Ef_ss})));
        P_e = double(R_sol(i)*double(subs(P_vals(e, i), {'Ef'}, {Ef_ss})));
        total_P_bound = total_P_bound + (e - 1) * P_e;
        total_P = total_P + P_e;
    end
    if total_P > 0
        avg_P(i) = total_P_bound / total_P;
    else
        avg_P(i) = 0;
    end
end

Ser2P = avg_P;

% ------------ PLOT 1: Ser2P and Average E bound ------------
figure;
hold on;
plot((1-PAS):N_PAS-1, Ser2P, 'g-','LineWidth',2.5, 'DisplayName','Ser2P');
plot((1-PAS):N_PAS-1, avg_E_bound, 'b-','LineWidth',2.5, 'DisplayName','AverageE');
legend({'Ser2P', 'AverageE'}, 'Location', 'best');
xlabel('Position relative to PAS (x100 bp)'); ylabel('Value');
title('Ser2P and Average E-factor Loading');
hold off;

% ------------ PLOT 2: Polymerase species along gene ------------
l_values = (1-PAS):(N-PAS);

figure; hold on;
plot(l_values, R_sol, 'b-','LineWidth',2.5, 'DisplayName','R(l) - Elongating');
plot(l_values, [zeros(PAS-1,1);REH_sol], 'r-','LineWidth',2.5, 'DisplayName','REH(l) - E+Hexamer');
plot(l_values, [zeros(PAS-1,1);REHL_sol], 'm-','LineWidth',2.5, 'DisplayName','REHL(l) - Committed');

CAD_nodes = CAD/100;
line([CAD_nodes CAD_nodes], ylim, 'Color', 'k', 'LineStyle', '--', 'LineWidth', 1.5, 'DisplayName', '50% CAD');
text(CAD_nodes, max(R_sol)*0.8, sprintf('CAD=%d bp', int32(CAD)));

xlabel('Position relative to PAS (x100 bp)'); ylabel('Concentration');
legend({'R (Elongating)', 'REH (E+Hexamer)', 'REHL (Committed)', '50% CAD'}, 'Location', 'best');
title('Explicit L-factor CPA Model');
hold off;

% ------------ SUMMARY STATISTICS ------------
bound_pol_II = sum(R_sol) + sum(REH_sol) + sum(REHL_sol);
Pol_f_final = P.Pol_total - bound_pol_II;
L_bound = sum(REHL_sol);
Lf_final = P.L_total - L_bound;

fprintf('\n--- Polymerase Distribution at Steady State ---\n');
fprintf('Total Pol II in system:       %d\n', P.Pol_total);
fprintf('Total Bound Pol II (on gene): %.2f\n', bound_pol_II);
fprintf('Total Free Pol II (Pol_f):    %.2f\n', Pol_f_final);
fprintf('\n--- L-factor Distribution ---\n');
fprintf('Total L-factors:     %d\n', P.L_total);
fprintf('Bound L (on REHL):   %.4f\n', L_bound);
fprintf('Free L (Lf):         %.2f\n', Lf_final);
fprintf('\n--- Assembly Metric ---\n');
fprintf('CAD (50%% cleavage distance): %.0f bp\n', CAD);

if saveData
    data.EBindingNumber = EBindingNumber;
    data.R_sol = R_sol;
    data.REH_sol = REH_sol;
    data.REHL_sol = REHL_sol;
    data.Ser2P = Ser2P;
    data.avg_E_bound = avg_E_bound;
    data.Ef_ss = Ef_ss;
    data.Lf_ss = Lf_ss;
    data.Pol_f_final = Pol_f_final;
    save_analysis_results('CPA_explicitL_main', data, P);
end
