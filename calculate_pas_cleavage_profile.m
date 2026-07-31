function [exit_cdf, distances_bp, CAD, diagnostics] = calculate_pas_cleavage_profile(R_sol, REH_sol, P_sim, varargin)
% CALCULATE_PAS_CLEAVAGE_PROFILE - Calculate polymerase cleavage profile downstream of PAS
%
% USAGE:
%   [exit_cdf, distances_bp] = calculate_pas_cleavage_profile(R_sol, REH_sol, P_sim)
%   [~, ~, CAD] = calculate_pas_cleavage_profile(R_sol, REH_sol, P_sim, 'PercentCleavage', 50)
%   [~, ~, CAD] = calculate_pas_cleavage_profile(..., 'AllowExtrapolation', true)
%
% INPUTS:
%   R_sol            - Steady-state polymerase concentrations
%   REH_sol          - Steady-state REH complex concentrations
%   P_sim            - Parameter structure with fields: kc, k_e, k_e2, L_a
%   'PercentCleavage'  - (Optional) CDF threshold (0-100). Returns the downstream
%                        distance (bp) at which this % of polymerases have cleaved.
%   'AllowExtrapolation' - (Optional, default false). If false, CAD is NaN when
%                        the threshold is not reached inside the simulated window.
%
% OUTPUTS:
%   exit_cdf     - Cumulative exit probability at each node
%   distances_bp - Distances downstream of PAS (bp)
%   CAD          - Cleavage/Arrest Distance at the requested PercentCleavage threshold (bp)
%   diagnostics  - Struct with fluxes, threshold status, and mass-balance residual

% Parse inputs
p = inputParser;
addParameter(p, 'PercentCleavage', [], @isnumeric);
addParameter(p, 'AllowExtrapolation', false, @(x) islogical(x) || isnumeric(x));
addParameter(p, 'MassBalanceTolerance', 1e-6, @isnumeric);
parse(p, varargin{:});
percent_cleavage = p.Results.PercentCleavage;
allow_extrapolation = logical(p.Results.AllowExtrapolation);
mass_balance_tolerance = p.Results.MassBalanceTolerance;

% Handle NaN inputs
if any(isnan(R_sol)) || any(isnan(REH_sol))
    exit_cdf = NaN(size(REH_sol));
    distances_bp = NaN(size(REH_sol));
    CAD = NaN;
    diagnostics = struct('threshold_reached', false, 'mass_balance_residual', NaN, ...
                         'mass_balance_relative_residual', NaN);
    return;
end

% Calculate fluxes
R_sol = R_sol(:);
REH_sol = REH_sol(:);

% The ODE removes cleavage strictly downstream of the PAS node. Keep the flux
% bookkeeping consistent with ode_dynamics_multipleE.m by excluding REH(1).
flux_cleavage_per_node = P_sim.kc * REH_sol;
if ~isempty(flux_cleavage_per_node)
    flux_cleavage_per_node(1) = 0;
end
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

% Global steady-state mass balance: initiation influx should equal all exit
% fluxes when R/REH are an actual steady-state solution of the ODE.
if isfield(P_sim, 'Pol_free')
    Pol_f = P_sim.Pol_free;
elseif isfield(P_sim, 'Pol_total')
    Pol_f = P_sim.Pol_total - sum(R_sol) - sum(REH_sol);
else
    Pol_f = NaN;
end

if isfield(P_sim, 'k_in') && ~isnan(Pol_f)
    influx = Pol_f * P_sim.k_in;
    mass_balance_residual = influx - total_outflux;
    mass_balance_relative_residual = abs(mass_balance_residual) / max(abs(influx), 1);
else
    influx = NaN;
    mass_balance_residual = NaN;
    mass_balance_relative_residual = NaN;
end

diagnostics = struct();
diagnostics.flux_cleavage_per_node = flux_cleavage_per_node;
diagnostics.flux_R_exit = flux_R_exit;
diagnostics.flux_REH_exit = flux_REH_exit;
diagnostics.total_outflux = total_outflux;
diagnostics.influx = influx;
diagnostics.mass_balance_residual = mass_balance_residual;
diagnostics.mass_balance_relative_residual = mass_balance_relative_residual;
diagnostics.max_exit_cdf = max(exit_cdf);
diagnostics.threshold_reached = true;

if ~isnan(mass_balance_relative_residual) && mass_balance_relative_residual > mass_balance_tolerance
    warning('calculate_pas_cleavage_profile:MassBalanceResidual', ...
            ['Steady-state mass-balance residual is %.3g (relative %.3g). ' ...
             'Check R/REH steady state and flux bookkeeping.'], ...
            mass_balance_residual, mass_balance_relative_residual);
end

% Interpolate CDF to find distance at requested cleavage percentage
if ~isempty(percent_cleavage)
    cdf_threshold = percent_cleavage / 100;
    diagnostics.threshold_reached = max(exit_cdf) >= cdf_threshold;
    x_cdf = [0; exit_cdf(:)];
    y_dist = [0; distances_bp(:)];
    [x_cdf, unique_idx] = unique(x_cdf, 'stable');
    y_dist = y_dist(unique_idx);

    if diagnostics.threshold_reached
        CAD = interp1(x_cdf, y_dist, cdf_threshold, 'linear');
    elseif allow_extrapolation
        warning('calculate_pas_cleavage_profile:ExtrapolatedCAD', ...
                ['Requested %.1f%% cleavage is not reached in the simulated %.0f bp window ' ...
                 '(max CDF %.3f); extrapolating because AllowExtrapolation=true.'], ...
                percent_cleavage, distances_bp(end), max(exit_cdf));
        if numel(x_cdf) >= 2
            CAD = interp1(x_cdf, y_dist, cdf_threshold, 'linear', 'extrap');
        else
            CAD = NaN;
        end
    else
        warning('calculate_pas_cleavage_profile:UnreachedThreshold', ...
                ['Requested %.1f%% cleavage is not reached in the simulated %.0f bp window ' ...
                 '(max CDF %.3f); returning CAD=NaN.'], ...
                percent_cleavage, distances_bp(end), max(exit_cdf));
        CAD = NaN;
    end
else
    CAD = [];
end

% Ensure column vectors
exit_cdf = exit_cdf(:);
distances_bp = distances_bp(:);

end
