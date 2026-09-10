function [R_sol, RHE_sol, P, details] = run_full_termination_simulation(P, EBindingNumber)
% RUN_FULL_TERMINATION_SIMULATION Full finite-rate steady state (no fast closure).
%
% [R,RHE,P_out,D] = run_full_termination_simulation(P,M)
% R and RHE are aggregate column profiles compatible with cleavage analysis.
% P_out includes N, PAS, N_PAS, Ef_ss, Pol_free_ss and kEoff_engaged.
% P.kHon remains the BASE per-E recognition rate; it is never rescaled.
% D contains microstates, state maps, the full steady vector and diagnostics.
%
% At fixed E_free, ALL microstate/transport equations are linear. Solve them
% per unit initiation flux, scale by Pol conservation, and use fzero for E
% conservation. This algebraic steady solve does not assume fast binding.
[model, P] = build_full_rate_matrices(P, EBindingNumber);
if min([P.k_in, P.k_e, P.k_e2]) <= 0
    error('run_full_termination_simulation:InvalidTransport', ...
        'The steady solver requires positive initiation and elongation rates.');
end
if P.E_total == 0
    Ef = 0;
else
    % fzero's stopping tolerance is on E_free, not the E-balance residual.
    % With E_free ~1e5, TolX=1e-9 can stop before the absolute 1e-6
    % conservation check passes. Use fzero's double-precision default.
    options = optimset('TolX', eps, 'FunValCheck', 'on');
    [Ef, ~, exitflag] = fzero(@(value) E_residual(value, P, model), [0, P.E_total], options);
    if exitflag <= 0
        error('run_full_termination_simulation:EfreeNotConverged', ...
            'Free-E conservation solve failed (exitflag %d).', exitflag);
    end
end
[x, Pol_f] = conditional_state(Ef, P, model);
X = [x; Ef; Pol_f];
R_micro = reshape(x(1:model.r_size), model.sr, model.N)';
RHE_micro = reshape(x(model.r_size+1:end), model.sh, model.N_PAS)';
R_sol = sum(R_micro, 2);
RHE_sol = sum(RHE_micro, 2);
P.Ef_ss = Ef;
P.Pol_free_ss = Pol_f;
P.EBindingNumber = EBindingNumber;
P.model_variant = 'full_kinetics_rapid_EH_disassembly';

details.R_micro = R_micro;
details.RHE_micro = RHE_micro;
details.R_states = model.states.R;
details.RHE_states = model.states.RHE;
details.state_ss = X;
details.E_conservation_residual = Ef+model.e_counts'*x-P.E_total;
details.Pol_conservation_residual = Pol_f+sum(x)-P.Pol_total;
details.rhs_max_abs = norm(ode_dynamics_full_multipleE(0, X, model), inf);
details.flux_residual = P.k_e*R_sol(end)+P.k_e2*RHE_sol(end) ...
                      + P.kc*sum(RHE_sol)-P.k_in*Pol_f;
if any(~isfinite(X)) || min(X) < -1e-9 || ...
        max(abs([details.E_conservation_residual, details.Pol_conservation_residual, ...
                 details.rhs_max_abs, details.flux_residual])) > 1e-6
    error('run_full_termination_simulation:InvalidSolution', ...
        ['Full steady solution failed validation for M=%d, kEoff_engaged=%.6g. ' ...
         'E_free=%.16g; nonfinite populations=%d; min population=%.3e ' ...
         '(must be >= -1e-9); E balance=%.3e; Pol balance=%.3e; ' ...
         'max |dX/dt|=%.3e; flux balance=%.3e ' ...
         '(each absolute residual must be <= 1e-6).'], ...
        EBindingNumber, P.kEoff_engaged, Ef, nnz(~isfinite(X)), min(X), ...
        details.E_conservation_residual, details.Pol_conservation_residual, ...
        details.rhs_max_abs, details.flux_residual);
end
details.model = model; % reusable sparse operator for optional time integration
end

function residual = E_residual(Ef, P, model)
[x, ~] = conditional_state(Ef, P, model);
residual = Ef+model.e_counts'*x-P.E_total;
end

function [x, Pol_f] = conditional_state(Ef, P, model)
occupation = -(model.A0+Ef*model.AE) \ model.unit_source;
if any(~isfinite(occupation)) || min(occupation) < -1e-9
    error('run_full_termination_simulation:InvalidConditionalState', ...
        'Conditional linear solve produced invalid microstate populations.');
end
flux = P.Pol_total/(1/P.k_in+sum(occupation));
x = flux*occupation;
Pol_f = flux/P.k_in;
end
