function [r_E_num, r_P_num] = compute_steady_states_numerical(kPon_val, kPoff_val, kEon_val, kEoff_val, Ef_val, n)
    % COMPUTE_STEADY_STATES_NUMERICAL Compute steady-state distributions numerically
    %
    % Numerical version of the symbolic compute_steady_states + substitution.
    % Returns the marginalized probability distributions for E-binding and
    % P-phosphorylation states at a specific set of parameter values.
    %
    % Inputs:
    %   kPon_val  - kPon rate constant (scalar)
    %   kPoff_val - kPoff rate constant (scalar)
    %   kEon_val  - kEon rate constant (scalar)
    %   kEoff_val - kEoff rate constant (scalar)
    %   Ef_val    - Free E concentration (scalar)
    %   n         - Number of states = EBindingNumber + 1
    %
    % Outputs:
    %   r_E_num - Marginalized E-binding distribution (1 x n), r_E_num(i) = P(E = i-1)
    %   r_P_num - Marginalized P-phosphorylation distribution (1 x n), r_P_num(i) = P(P-block = i-1)

    % Build the numerical rate matrix
    A = build_rate_matrix_numerical(n, kPon_val, kPoff_val, kEon_val, kEoff_val, Ef_val);

    % Compute null space via SVD
    [~, S, V] = svd(A);
    s = diag(S);
    [~, min_idx] = min(abs(s));
    ss_dist = V(:, min_idx);

    % Ensure non-negative
    if sum(ss_dist) < 0
        ss_dist = -ss_dist;
    end

    % Normalize
    ss_dist = ss_dist / sum(ss_dist);

    total_states = length(ss_dist);

    % Compute r_E: marginalized E distribution using same indexing as
    % compute_cumulative_sum in compute_steady_states.m
    RE_index = cumsum(0:n-1) + 1;
    r_E_num = zeros(1, n);

    for i = 1:n
        RE_index_filtered = RE_index(RE_index <= total_states);
        RE_index_copy = RE_index_filtered(i:end);
        if ~isempty(RE_index_copy)
            r_E_num(i) = sum(ss_dist(RE_index_copy));
        end
        RE_index = RE_index + 1;
    end

    % Normalize r_E
    r_E_num = r_E_num / sum(r_E_num);

    % Compute r_P: marginalized P distribution using same indexing as
    % compute_block_sums in compute_steady_states.m
    block_lengths = 1:n;
    start_indices = cumsum([1 block_lengths(1:end-1)]);
    r_P_num = zeros(1, n);

    for i = 1:n
        block_start = start_indices(i);
        block_end = block_start + block_lengths(i) - 1;
        if block_start > total_states
            r_P_num(i) = 0;
            continue;
        end
        block_end = min(block_end, total_states);
        r_P_num(i) = sum(ss_dist(block_start:block_end));
    end

    % Normalize r_P
    r_P_num = r_P_num / sum(r_P_num);
end


function A = build_rate_matrix_numerical(n, kPon, kPoff, kEon, kEoff, Ef)
    % BUILD_RATE_MATRIX_NUMERICAL Numerical rate matrix construction
    % Replicates the logic of construct_rate_matrix.m with numerical values.

    sz = [n, n];
    A_full = zeros(n*n, n*n);

    % --- kPon/kPoff transitions (phosphorylation) ---
    rows_kPon = repelem(1:n, fliplr(1:n));
    cellArray1 = arrayfun(@(x) x:n, 1:n, 'UniformOutput', false);
    cols_kPon = [cellArray1{:}];
    ind_kPon = sub2ind(sz, rows_kPon, cols_kPon);

    for i = 1:length(ind_kPon)-2
        current_idx = ind_kPon(i);
        next_idx = ind_kPon(i+1);
        if next_idx > current_idx
            A_full(next_idx, current_idx) = kPon;
            A_full(current_idx, next_idx) = kPoff;
        end
    end

    % --- kEon/kEoff transitions (E binding) ---
    cellArray2 = arrayfun(@(x) 1:x, 1:n, 'UniformOutput', false);
    rwos_kEon = [cellArray2{:}];
    cols_kEon = repelem(1:n, 1:n);
    ind_kEon = sub2ind(sz, rwos_kEon, cols_kEon);

    a1 = 1;
    a2 = 0;
    for i = 1:length(ind_kEon)-1
        current_idx = ind_kEon(i);
        next_idx = ind_kEon(i+1);
        if next_idx == current_idx + 1
            A_full(next_idx, current_idx) = a2 * kEon * Ef;
            A_full(current_idx, next_idx) = a1 * kEoff;
            a1 = a1 + 1;
            a2 = a2 - 1;
        else
            a2 = a1;
            a1 = 1;
        end
    end

    % Remove all-zero rows and columns
    row_mask = any(A_full ~= 0, 2);
    col_mask = any(A_full ~= 0, 1);
    A = A_full(row_mask, col_mask);

    % Set diagonal: A(i,i) = -sum(column i)
    for i = 1:size(A, 1)
        A(i, i) = -sum(A(:, i)) + A(i, i);
    end
end
