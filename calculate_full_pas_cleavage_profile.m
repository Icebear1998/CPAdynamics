function [exit_cdf, distances_bp, CAD, diagnostics] = calculate_full_pas_cleavage_profile(R_sol, RHE_sol, P, varargin)
% CALCULATE_FULL_PAS_CLEAVAGE_PROFILE Full-model flux CDF and cleavage distance.
% Cleavage includes RHE(1), representing (0,L_a] bp; its CDF is at L_a.
% Interpolation includes (0,0). Unreached thresholds return NaN, not an
% extrapolated distance. Returned vectors exclude the prepended origin for
% compatibility with the original MATLAB profile interface.
parser = inputParser;
addParameter(parser, 'PercentCleavage', 50, ...
    @(v) isnumeric(v) && isscalar(v) && isreal(v) && isfinite(v) && v >= 0 && v <= 100);
parse(parser, varargin{:});
validateattributes(R_sol, {'numeric'}, {'vector', 'nonempty', 'real', 'finite', 'nonnegative'});
validateattributes(RHE_sol, {'numeric'}, {'vector', 'nonempty', 'real', 'finite', 'nonnegative'});
R_sol = R_sol(:);
RHE_sol = RHE_sol(:);
cleavage_flux = P.kc*RHE_sol;
total_outflux = sum(cleavage_flux)+P.k_e*R_sol(end)+P.k_e2*RHE_sol(end);
if total_outflux > 1e-9
    exit_cdf = cumsum(cleavage_flux)/total_outflux;
else
    exit_cdf = zeros(size(RHE_sol));
end
distances_bp = (1:numel(RHE_sol))'*P.L_a;
threshold = parser.Results.PercentCleavage/100;
CAD = NaN;
if threshold == 0
    CAD = 0;
else
    j = find(exit_cdf >= threshold, 1, 'first');
    if isempty(j)
        warning('calculate_full_pas_cleavage_profile:ThresholdNotReached', ...
            'Requested %.1f%% cleavage is not reached; the window reaches %.2f%%.', ...
            parser.Results.PercentCleavage, 100*exit_cdf(end));
    else
        cdf_with_origin = [0; exit_cdf];
        distances_with_origin = [0; distances_bp];
        fraction = (threshold-cdf_with_origin(j))/(cdf_with_origin(j+1)-cdf_with_origin(j));
        CAD = distances_with_origin(j)+fraction*(distances_with_origin(j+1)-distances_with_origin(j));
    end
end
diagnostics.max_exit_cdf = exit_cdf(end);
diagnostics.total_outflux = total_outflux;
diagnostics.threshold_reached = ~isnan(CAD);
end
