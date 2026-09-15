function [R_occupied, E_occupied, diagnostics] = validate_full_gene_length_grid(data, EBindingNumber)
% VALIDATE_FULL_GENE_LENGTH_GRID Accept solver roundoff, reject invalid rows.
% Sparse steady solves can leave tiny signed populations at exact-zero
% boundaries. The solver accepts roundoff down to -1e-9; an exact >= 0
% occupancy check must not subsequently discard those successful rows.
% Only roundoff is corrected. Failed/nonfinite/materially negative rows
% still prevent construction of a complete Cartesian interpolant.

validateattributes(EBindingNumber, {'numeric'}, ...
    {'scalar', 'real', 'finite', 'integer', 'positive'});
fields = {'R_free_vec', 'E_free_vec', 'L_vec', ...
          'R_occupied_vec', 'E_occupied_vec', 'success_flag'};
n = numel(data.success_flag);
for k = 1:numel(fields)
    validateattributes(data.(fields{k}), {'numeric'}, ...
        {'vector', 'real', 'nonempty', 'numel', n}, mfilename, fields{k});
    data.(fields{k}) = data.(fields{k})(:);
end
R_occupied = data.R_occupied_vec;
E_occupied = data.E_occupied_vec;

% Absolute solver tolerance plus a small floating-point allowance for sums
% of large populations. E sums are bounded in scale by M times Pol occupancy.
R_scale = max(ones(n, 1), abs(R_occupied));
E_scale = max([ones(n, 1), EBindingNumber*abs(R_occupied), abs(E_occupied)], [], 2);
R_tolerance = 1e-9 + 64*eps(R_scale);
E_tolerance = 1e-9 + 64*eps(E_scale);
failed = data.success_flag ~= 1;
nonfinite = ~isfinite(R_occupied) | ~isfinite(E_occupied);
negative = R_occupied < -R_tolerance | E_occupied < -E_tolerance;
bad_coordinates = ~isfinite(data.R_free_vec) | data.R_free_vec < 0 ...
    | ~isfinite(data.E_free_vec) | data.E_free_vec < 0 ...
    | ~isfinite(data.L_vec) | data.L_vec <= 0;
zero_pol = data.R_free_vec == 0;
zero_e = zero_pol | data.E_free_vec == 0;
bad_boundary = (zero_pol & abs(R_occupied) > R_tolerance) ...
    | (zero_e & abs(E_occupied) > E_tolerance);
invalid = failed | nonfinite | negative | bad_coordinates | bad_boundary;
if any(invalid)
    first = find(invalid, 1);
    solver_message = '';
    if isfield(data, 'error_messages') && numel(data.error_messages) >= first
        solver_message = data.error_messages{first};
    end
    error('GeneLengthBuildInterpolation:IncompleteGrid', ...
        ['Full-model interpolation requires a complete grid. Invalid rows: %d/%d ' ...
         '(failed=%d, nonfinite=%d, materially negative=%d, bad coordinates=%d, ' ...
         'inconsistent zero-pool boundary=%d; categories may overlap).\n' ...
         'First invalid row %d: R_free=%.16g, E_free=%.16g, L=%.16g, ' ...
         'R_occupied=%.16g, E_occupied=%.16g; tolerances R=%.3g, E=%.3g.\n' ...
         'Solver message: %s\nInspect these rows in results.data and regenerate if needed.'], ...
        nnz(invalid), n, nnz(failed), nnz(nonfinite), nnz(negative), ...
        nnz(bad_coordinates), nnz(bad_boundary), first, ...
        data.R_free_vec(first), data.E_free_vec(first), data.L_vec(first), ...
        R_occupied(first), E_occupied(first), R_tolerance(first), E_tolerance(first), solver_message);
end

correct_R = R_occupied < 0 | (zero_pol & R_occupied ~= 0);
correct_E = E_occupied < 0 | (zero_e & E_occupied ~= 0);
diagnostics.corrected_rows = nnz(correct_R | correct_E);
diagnostics.corrected_R_values = nnz(correct_R);
diagnostics.corrected_E_values = nnz(correct_E);
diagnostics.minimum_raw_R_occupied = min(R_occupied);
diagnostics.minimum_raw_E_occupied = min(E_occupied);
diagnostics.max_R_tolerance = max(R_tolerance);
diagnostics.max_E_tolerance = max(E_tolerance);
diagnostics.tolerance_rule = '1e-9 + 64*eps(population scale); E scale includes M*Pol';
R_occupied(correct_R) = 0;
E_occupied(correct_E) = 0;
end
