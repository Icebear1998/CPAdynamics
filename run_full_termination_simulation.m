function [R_sol, RHE_sol, P, details] = run_full_termination_simulation(P, EBindingNumber, varargin)
% RUN_FULL_TERMINATION_SIMULATION Full finite-rate steady state (no fast closure).
%
% [R,RHE,P_out,D] = run_full_termination_simulation(P,M)
% R and RHE are aggregate column profiles compatible with cleavage analysis.
% P_out includes N, PAS, N_PAS, Ef_ss, Pol_free_ss and kEoff_engaged.
% P.kHon remains the BASE per-E recognition rate; it is never rescaled.
% D contains microstates, state maps, the full steady vector and diagnostics.
% D.avg_E_bound and D.avg_Ser2P average over all Pol II at each node,
% including RHE after PAS. Empty nodes have NaN averages and zero occupancy.
%
% [R,RHE,P_out,D] = run_full_termination_simulation(P,M,'FreePools',[R_free E_free])
% holds BOTH free pools fixed for one gene in a genome-wide calculation.
% Initiation is k_in*R_free; P.Pol_total and P.E_total do not constrain
% this local solution. D.Pol_bound and D.E_bound are its resource demand.
% Conservation against P's totals is only tested in the default closed-pool
% mode; fixed-pool results instead validate the steady RHS and flux balance.
%
% At fixed E_free, ALL microstate/transport equations are linear. Solve them
% per unit initiation flux, scale by Pol conservation, and use fzero for E
% conservation. This algebraic steady solve does not assume fast binding.
parser = inputParser;
addParameter(parser, 'FreePools', [], @(v) isempty(v) || ...
    (isnumeric(v) && isreal(v) && isvector(v) && numel(v) == 2 && ...
     all(isfinite(v)) && all(v >= 0)));
parse(parser, varargin{:});
free_pools = parser.Results.FreePools;
[model, P] = build_full_rate_matrices(P, EBindingNumber);
if min([P.k_in, P.k_e, P.k_e2]) <= 0
    error('run_full_termination_simulation:InvalidTransport', ...
        'The steady solver requires positive initiation and elongation rates.');
end
if ~isempty(free_pools)
    Ef = free_pools(2);
elseif P.E_total == 0
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
if isempty(free_pools)
    [x, Pol_f] = conditional_state(Ef, P, model);
    details.pool_mode = 'self_consistent';
else
    [x, Pol_f] = conditional_state(Ef, P, model, free_pools(1));
    details.pool_mode = 'fixed_free_pools';
end
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
details.E_bound_per_node = R_micro*model.states.R(:, 2);
details.Ser2P_per_node = R_micro*model.states.R(:, 1);
details.E_bound_per_node(P.PAS:end) = details.E_bound_per_node(P.PAS:end) ...
    + RHE_micro*model.states.RHE(:, 2);
details.Ser2P_per_node(P.PAS:end) = details.Ser2P_per_node(P.PAS:end) ...
    + RHE_micro*model.states.RHE(:, 1);
pol_per_node = R_sol;
pol_per_node(P.PAS:end) = pol_per_node(P.PAS:end)+RHE_sol;
details.avg_E_bound = NaN(P.N, 1);
details.avg_Ser2P = NaN(P.N, 1);
occupied = pol_per_node > 0;
details.avg_E_bound(occupied) = details.E_bound_per_node(occupied)./pol_per_node(occupied);
details.avg_Ser2P(occupied) = details.Ser2P_per_node(occupied)./pol_per_node(occupied);
details.E_bound = model.e_counts'*x;
details.Pol_bound = sum(x);
details.E_conservation_residual = NaN;
details.Pol_conservation_residual = NaN;
conservation_residuals = [];
if isempty(free_pools)
    details.E_conservation_residual = Ef+details.E_bound-P.E_total;
    details.Pol_conservation_residual = Pol_f+details.Pol_bound-P.Pol_total;
    conservation_residuals = [details.E_conservation_residual, details.Pol_conservation_residual];
end
details.rhs_max_abs = norm(ode_dynamics_full_multipleE(0, X, model), inf);
details.flux_residual = P.k_e*R_sol(end)+P.k_e2*RHE_sol(end) ...
                      + P.kc*sum(RHE_sol)-P.k_in*Pol_f;
if any(~isfinite(X)) || min(X) < -1e-9 || ...
        max(abs([conservation_residuals, details.rhs_max_abs, details.flux_residual])) > 1e-6
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

function [x, Pol_f] = conditional_state(Ef, P, model, Pol_f)
occupation = -(model.A0+Ef*model.AE) \ model.unit_source;
if any(~isfinite(occupation)) || min(occupation) < -1e-9
    error('run_full_termination_simulation:InvalidConditionalState', ...
        'Conditional linear solve produced invalid microstate populations.');
end
if nargin < 4
    flux = P.Pol_total/(1/P.k_in+sum(occupation));
    Pol_f = flux/P.k_in;
else
    flux = P.k_in*Pol_f;
end
x = flux*occupation;
end
