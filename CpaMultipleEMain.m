% Main full finite-rate CPA simulation and microstate-derived profiles.
saveData = strcmpi(getenv('CPAD_FORCE_SAVE'), 'true');
P = default_parameters();
EBindingNumber = 5;

% Record the actual files resolved on the MATLAB path. A different checkout
% or shadowed parameter/helper file can change profiles without changing this script.
source_names = {'default_parameters', 'run_full_termination_simulation', ...
    'build_full_rate_matrices', 'build_full_internal_rate_matrices', 'build_full_state_map'};
source_paths = struct();
fprintf('=== Full finite-rate CPA: runtime sources ===\n');
for source_idx = 1:numel(source_names)
    source_name = source_names{source_idx};
    source_paths.(source_name) = which(source_name);
    fprintf('  %s: %s\n', source_name, source_paths.(source_name));
end

[R_sol, REH_sol, P, full_details] = run_full_termination_simulation(P, EBindingNumber);
fprintf('M=%d; kEon=%.8g, kEoff=%.8g, kEoff_engaged=%.8g\n', ...
    EBindingNumber, P.kEon, P.kEoff, P.kEoff_engaged);
fprintf('kHon=%.8g, kHoff=%.8g, kc=%.8g; k_e=%.8g, k_e2=%.8g\n', ...
    P.kHon, P.kHoff, P.kc, P.k_e, P.k_e2);
fprintf('kPon_min=%.8g, kPon_slope=%.8g, kPoff=%.8g; E_total=%.8g, Pol_total=%.8g\n', ...
    P.kPon_min, P.kPon_slope, P.kPoff, P.E_total, P.Pol_total);
avg_E_bound = full_details.avg_E_bound;
Ser2P = full_details.avg_Ser2P;
[exit_cdf, distances_bp, CAD, diagnostics] = calculate_full_pas_cleavage_profile( ...
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
sample_nodes = unique([max(1, P.PAS-1), P.PAS, min(P.N, P.PAS+round(500/P.L_a)), P.N]);
fprintf('  Position relative to PAS (bp)    Average E    Ser2P\n');
for sample_idx = sample_nodes
    fprintf('  %28g    %9.6f    %9.6f\n', ...
        positions_bp(sample_idx), avg_E_bound(sample_idx), Ser2P(sample_idx));
end

if saveData
    data = struct();
    data.model_variant = P.model_variant;
    data.source_paths = source_paths;
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
