% Main full finite-rate CPA simulation and microstate-derived profiles.
saveData = strcmpi(getenv('CPAD_FORCE_SAVE'), 'true');
P = default_parameters();
EBindingNumber = 5;

[R_sol, REH_sol, P, full_details] = run_full_termination_simulation(P, EBindingNumber);
avg_E_bound = full_details.avg_E_bound;
Ser2P = full_details.avg_Ser2P;
[~, ~, CAD, diagnostics] = calculate_full_pas_cleavage_profile( ...
    R_sol, REH_sol, P, 'PercentCleavage', 50);
positions_bp = ((1:P.N)-P.PAS)*P.L_a;

figure;
hold on;
plot(positions_bp, Ser2P, 'g-', 'LineWidth', 2.5, 'DisplayName', 'Ser2P');
plot(positions_bp, avg_E_bound, 'b-', 'LineWidth', 2.5, 'DisplayName', 'Average E');
legend('Location', 'best');
xlabel('Position relative to PAS (bp)');
ylabel('Average count per polymerase');
title('Full finite-rate E binding and phosphorylation');
hold off;

figure;
hold on;
plot(positions_bp, R_sol, 'b-', 'LineWidth', 2.5, 'DisplayName', 'R');
plot(positions_bp, [zeros(P.PAS-1, 1); REH_sol], 'r-', ...
    'LineWidth', 2.5, 'DisplayName', 'RHE');
if isfinite(CAD)
    xline(CAD, 'r--', 'LineWidth', 1.5, 'DisplayName', 'CAD_{50}');
end
xlabel('Position relative to PAS (bp)');
ylabel('Polymerase count');
legend('Location', 'best');
title('Full finite-rate CPA steady state');
hold off;

Pol_f_final = P.Pol_free_ss;
fprintf('CAD_50 = %.2f bp; within-window cleavage = %.2f%%\n', CAD, 100*diagnostics.max_exit_cdf);
fprintf('Free E = %.2f; bound E = %.2f\n', P.Ef_ss, full_details.E_bound);
fprintf('Free Pol II = %.2f; bound Pol II = %.2f\n', Pol_f_final, full_details.Pol_bound);

if saveData
    data = struct();
    data.model_variant = P.model_variant;
    data.EBindingNumber = EBindingNumber;
    data.R_sol = R_sol;
    data.REH_sol = REH_sol;
    data.Ser2P = Ser2P;
    data.avg_E_bound = avg_E_bound;
    data.Ef_ss = P.Ef_ss;
    data.Pol_f_final = Pol_f_final;
    data.max_exit_cdf = diagnostics.max_exit_cdf;
    data.rhs_residuals = full_details.rhs_max_abs;
    save_analysis_results('CPA_multipleE_main', data, P);
end
