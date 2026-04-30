function [exit_cdf, distances_bp, CAD] = calculate_pas_cleavage_profile(R_sol, REH_sol, P_sim, varargin)
% CALCULATE_PAS_CLEAVAGE_PROFILE - Calculate polymerase cleavage profile downstream of PAS
%
% USAGE:
%   [exit_cdf, distances_bp] = calculate_pas_cleavage_profile(R_sol, REH_sol, P_sim)
%   [~, ~, CAD] = calculate_pas_cleavage_profile(R_sol, REH_sol, P_sim, 'PercentCleavage', 50)
%
% INPUTS:
%   R_sol            - Steady-state polymerase concentrations
%   REH_sol          - Steady-state REH complex concentrations
%   P_sim            - Parameter structure with fields: kc, k_e, k_e2, L_a
%   'PercentCleavage'- (Optional) CDF threshold (0-100). Returns the downstream
%                      distance (bp) at which this % of polymerases have cleaved.
%
% OUTPUTS:
%   exit_cdf     - Cumulative exit probability at each node
%   distances_bp - Distances downstream of PAS (bp)
%   CAD          - Cleavage/Arrest Distance at the requested PercentCleavage threshold (bp)

% Parse inputs
p = inputParser;
addParameter(p, 'PercentCleavage', [], @isnumeric);
parse(p, varargin{:});
percent_cleavage = p.Results.PercentCleavage;

% Handle NaN inputs
if any(isnan(R_sol)) || any(isnan(REH_sol))
    exit_cdf = NaN(size(REH_sol));
    distances_bp = NaN(size(REH_sol));
    CAD = NaN;
    return;
end

% Calculate fluxes
flux_cleavage_per_node = P_sim.kc * REH_sol;
flux_R_exit = P_sim.k_e * R_sol(end);
flux_REH_exit = P_sim.k_e2 * REH_sol(end);
total_outflux = sum(flux_cleavage_per_node) + flux_R_exit + flux_REH_exit;

% Calculate CDF
if total_outflux > 1e-9
    cumulative_exit_flux = cumsum(flux_cleavage_per_node);
    exit_cdf = cumulative_exit_flux / total_outflux;
else
    exit_cdf = zeros(size(REH_sol));
end

% Calculate distances
nodes_post_pas = 1:length(REH_sol);
distances_bp = nodes_post_pas * P_sim.L_a;

% Interpolate CDF to find distance at requested cleavage percentage
if ~isempty(percent_cleavage)
    cdf_threshold = percent_cleavage / 100;
    CAD = interp1([0; exit_cdf(:)], [0; distances_bp(:)], cdf_threshold, 'linear', 'extrap');
else
    CAD = [];
end

% Ensure column vectors
exit_cdf = exit_cdf(:);
distances_bp = distances_bp(:);

end
