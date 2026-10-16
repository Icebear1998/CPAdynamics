function values = sweep_values_from_ranges(P, param, n_points)
% SWEEP_VALUES_FROM_RANGES  Sweep grid spanning P.ranges.(param).
%
%   values = sweep_values_from_ranges(P, param, n_points)
%
% P is the struct from default_parameters(). Elongation rates and pool
% sizes are spaced linearly; every rate constant is spaced logarithmically.
% Errors if default_parameters supplies no range for param.

    linear_params = {'k_e', 'k_e2', 'E_total', 'Pol_total'};

    if ~isfield(P, 'ranges') || ~isfield(P.ranges, param)
        error('sweep_values_from_ranges:NoRange', ...
              'default_parameters supplies no range for "%s".', param);
    end
    bounds = P.ranges.(param);

    if ismember(param, linear_params)
        values = linspace(bounds(1), bounds(2), n_points);
    else
        values = logspace(log10(bounds(1)), log10(bounds(2)), n_points);
    end
end
