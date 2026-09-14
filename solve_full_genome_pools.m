function solution = solve_full_genome_pools(interpolation_results, L_values, weights, R_total, E_total, num_genes)
% SOLVE_FULL_GENOME_POOLS Shared free pools using full-model gene occupancies.
% weights are normalized integration weights for L_values. Grid data must
% come from the fixed-free-pool, finite-rate model and cover [0,R_total],
% [0,E_total] and every requested length. No extrapolation is permitted.
% At fixed E_free, occupancy is exactly proportional to R_free. Eliminate
% R_free analytically and solve only E conservation on [0,E_total].

validateattributes(R_total, {'numeric'}, {'scalar', 'real', 'finite', 'positive'});
validateattributes(E_total, {'numeric'}, {'scalar', 'real', 'finite', 'nonnegative'});
validateattributes(num_genes, {'numeric'}, {'scalar', 'real', 'finite', 'integer', 'positive'});
validateattributes(L_values, {'numeric'}, {'vector', 'nonempty', 'real', 'finite', 'positive'});
validateattributes(weights, {'numeric'}, {'vector', 'numel', numel(L_values), 'real', 'finite', 'nonnegative'});
L_values = L_values(:);
weights = weights(:);
if abs(sum(weights)-1) > 1e-8
    error('solve_full_genome_pools:InvalidWeights', 'Integration weights must sum to one.');
end
metadata = interpolation_results.metadata;
if ~isfield(metadata, 'model_variant') || ...
        ~strcmp(metadata.model_variant, 'full_kinetics_rapid_EH_disassembly') || ...
        ~isfield(metadata, 'pool_mode') || ~strcmp(metadata.pool_mode, 'fixed_free_pools')
    error('solve_full_genome_pools:IncompatibleModel', ...
        'Rebuild the gene-length grid and interpolation using the full finite-rate model.');
end
grid = interpolation_results.original_grid;
if grid.R_free_range(1) > 0 || grid.R_free_range(2) < R_total || ...
        grid.E_free_range(1) > 0 || grid.E_free_range(2) < E_total || ...
        any(L_values < grid.L_range(1) | L_values > grid.L_range(2))
    error('solve_full_genome_pools:OutsideGrid', ...
        'Regenerate a grid covering zero through the target totals and all integration lengths.');
end
R_interp = interpolation_results.functions.R_occupied_interp;
E_interp = interpolation_results.functions.E_occupied_interp;

if E_total == 0
    E_free = 0;
else
    options = optimset('TolX', eps, 'FunValCheck', 'on');
    [E_free, ~, exitflag] = fzero(@evaluate_pools, [0, E_total], options);
    if exitflag <= 0
        error('solve_full_genome_pools:NotConverged', 'Shared free-E conservation solve did not converge.');
    end
end
[E_residual, R_free, mean_R, mean_E] = evaluate_pools(E_free);
residuals = [R_free+num_genes*mean_R-R_total, E_residual];
if any(~isfinite(residuals)) || max(abs(residuals)) > 1e-6 || ...
        R_free < 0 || R_free > R_total || E_free < 0 || E_free > E_total
    error('solve_full_genome_pools:InvalidSolution', ...
        'Shared pool solution failed conservation: R residual %.3g, E residual %.3g.', residuals(1), residuals(2));
end
solution.R_free = R_free;
solution.E_free = E_free;
solution.R_total_target = R_total;
solution.E_total_target = E_total;
solution.num_active_genes = num_genes;
solution.mean_R_occupied = mean_R;
solution.mean_E_occupied = mean_E;
solution.conservation_residuals = residuals;

    function [residual, trial_R, mean_R_at_trial, mean_E_at_trial] = evaluate_pools(trial_E)
        query_shape = ones(size(L_values));
        R_at_reference = R_interp(R_total*query_shape, trial_E*query_shape, L_values);
        E_at_reference = E_interp(R_total*query_shape, trial_E*query_shape, L_values);
        if any(~isfinite(R_at_reference)) || any(~isfinite(E_at_reference)) || ...
                any(R_at_reference < 0) || any(E_at_reference < 0)
            error('solve_full_genome_pools:InvalidInterpolation', ...
                'Occupancy interpolation returned invalid values; rebuild the grid.');
        end
        mean_R_per_free = sum(weights.*R_at_reference(:))/R_total;
        trial_R = R_total/(1+num_genes*mean_R_per_free);
        mean_R_at_trial = trial_R*mean_R_per_free;
        mean_E_at_trial = (trial_R/R_total)*sum(weights.*E_at_reference(:));
        residual = trial_E+num_genes*mean_E_at_trial-E_total;
    end
end
