function [exit_cdf, distances_bp, usage_at_distance] = calculate_pas_usage_profile(R_sol, REH_sol, REHL_sol, P_sim, varargin)
% CALCULATE_PAS_USAGE_PROFILE - Calculate polymerase exit profile downstream of PAS
%   (Explicit L-factor version: cleavage occurs from REHL, not REH)
%
% USAGE:
%   [exit_cdf, distances_bp] = calculate_pas_usage_profile(R_sol, REH_sol, REHL_sol, P_sim)
%   [~, ~, usage_at_distance] = calculate_pas_usage_profile(R_sol, REH_sol, REHL_sol, P_sim, 'Distance', target_distance)
%
% INPUTS:
%   R_sol      - Steady-state polymerase concentrations (all N nodes)
%   REH_sol    - Steady-state REH complex concentrations (N_PAS nodes)
%   REHL_sol   - Steady-state REHL complex concentrations (N_PAS nodes)
%   P_sim      - Parameter structure with fields: kc, k_e, k_e2, L_a
%   'Distance' - (Optional) Target distance for interpolation (bp)
%
% OUTPUTS:
%   exit_cdf           - Cumulative exit probability at each node
%   distances_bp       - Distances downstream of PAS (bp)
%   usage_at_distance  - Usage probability at target distance (percentage)

% Parse inputs
p = inputParser;
addParameter(p, 'Distance', [], @isnumeric);
parse(p, varargin{:});
target_distance = p.Results.Distance;

% Handle NaN inputs
if any(isnan(R_sol)) || any(isnan(REH_sol)) || any(isnan(REHL_sol))
    exit_cdf = NaN(size(REHL_sol));
    distances_bp = NaN(size(REHL_sol));
    usage_at_distance = NaN;
    return;
end

% Calculate fluxes
% Cleavage only from REHL (fully committed state)
flux_cleavage_per_node = P_sim.kc * REHL_sol;

% Run-off flux: R, REH, and REHL at terminal node
flux_R_exit   = P_sim.k_e  * R_sol(end);
flux_REH_exit = P_sim.k_e2 * REH_sol(end);
flux_REHL_exit = P_sim.k_e2 * REHL_sol(end);
total_outflux = sum(flux_cleavage_per_node) + flux_R_exit + flux_REH_exit + flux_REHL_exit;

% Calculate CDF
if total_outflux > 1e-9
    cumulative_exit_flux = cumsum(flux_cleavage_per_node);
    exit_cdf = cumulative_exit_flux / total_outflux;
else
    exit_cdf = zeros(size(REHL_sol));
end

% Calculate distances
nodes_post_pas = 1:length(REHL_sol);
distances_bp = nodes_post_pas * P_sim.L_a;

% Interpolate if distance specified
if ~isempty(target_distance) && nargout > 2
    bp_for_interp = [0; distances_bp(:)];
    cdf_for_interp = [0; exit_cdf(:)];
    usage_at_distance = interp1(bp_for_interp, cdf_for_interp, target_distance, 'linear', 'extrap') * 100;
else
    usage_at_distance = [];
end

% Ensure column vectors
exit_cdf = exit_cdf(:);
distances_bp = distances_bp(:);

end
