% Physical checks for the full finite-rate model. Failures are assertions;
% solver errors are never converted into artificial zero-population passes.
P = default_parameters();
EBindingNumber = 1;
[R, RHE, P_sim, details] = run_full_termination_simulation(P, EBindingNumber);
check_full_solution(P_sim, details);
fprintf('Full-model baseline: resource conservation, nonnegativity and steady flux [PASS]\n');

positions_bp = ((1:P_sim.N)-P_sim.PAS)*P_sim.L_a;
figure;
plot(positions_bp, R, 'b-', positions_bp, [zeros(P_sim.PAS-1, 1); RHE], 'r-', 'LineWidth', 2);
xlabel('Position relative to PAS (bp)');
ylabel('Polymerase count');
legend('R', 'RHE', 'Location', 'best');
title('Full finite-rate steady state');
grid on;

% E is required for recognition; disabling E supply or either on-rate
% removes RHE while preserving transport and total-resource conservation.
for field = {'E_total', 'kEon', 'kHon'}
    P_case = P;
    P_case.(field{1}) = 0;
    [~, RHE_case, P_case, D] = run_full_termination_simulation(P_case, EBindingNumber);
    check_full_solution(P_case, D);
    assert(max(abs(RHE_case)) < 1e-8, 'RHE persists with %s = 0.', field{1});
    if ~strcmp(field{1}, 'kHon')
        assert(abs(D.E_bound) < 1e-8, 'Bound E persists without E supply or binding.');
    end
    fprintf('%s = 0: no PAS recognition [PASS]\n', field{1});
end

% With no cleavage, polymerases still leave the finite gene by elongation.
P_case = P;
P_case.kc = 0;
[R_case, RHE_case, P_case, D] = run_full_termination_simulation(P_case, EBindingNumber);
check_full_solution(P_case, D);
[cdf, ~] = calculate_full_pas_cleavage_profile(R_case, RHE_case, P_case, 'PercentCleavage', 0);
assert(all(cdf == 0), 'Cleavage must vanish when kc = 0.');
assert(abs(P_case.k_e*R_case(end)+P_case.k_e2*RHE_case(end)-P_case.k_in*P_case.Pol_free_ss) < 1e-6, ...
    'Initiation and terminal elongation flux must balance when kc = 0.');
fprintf('kc = 0: no cleavage, terminal elongation balances initiation [PASS]\n');

% The algebraic steady solver explicitly requires positive transport rates.
for field = {'k_in', 'k_e', 'k_e2'}
    P_case = P;
    P_case.(field{1}) = 0;
    rejected = false;
    try
        run_full_termination_simulation(P_case, EBindingNumber);
    catch ME
        if ~strcmp(ME.identifier, 'run_full_termination_simulation:InvalidTransport')
            rethrow(ME);
        end
        rejected = true;
    end
    assert(rejected, 'Zero %s must be rejected by the steady solver.', field{1});
    fprintf('%s = 0: documented unsupported steady-solver input [PASS]\n', field{1});
end
fprintf('All full-model sanity checks passed.\n');

function check_full_solution(P, D)
    E_from_microstates = sum(D.R_micro*D.R_states(:, 2)) ...
        + sum(D.RHE_micro*D.RHE_states(:, 2));
    Pol_from_microstates = sum(D.R_micro(:))+sum(D.RHE_micro(:));
    assert(abs(P.Ef_ss+E_from_microstates-P.E_total) < 1e-6, 'E conservation failed.');
    assert(abs(P.Pol_free_ss+Pol_from_microstates-P.Pol_total) < 1e-6, 'Pol II conservation failed.');
    assert(min(D.state_ss) >= -1e-9, 'Negative microstate population.');
    assert(D.rhs_max_abs < 1e-6 && abs(D.flux_residual) < 1e-6, 'Steady-state balance failed.');
end
