function avg_E_bound = compute_avg_E_bound_numerical(Ef_val, kPon_vals, kPoff_val, kEon_val, kEoff_val, n)
    % COMPUTE_AVG_E_BOUND_NUMERICAL Compute average E bound at each position numerically
    %
    % Instead of using symbolic null-space computation (which becomes numerically
    % unstable for large n due to extremely complex rational expressions), this
    % function builds the rate matrix numerically for each position and computes
    % the null space directly using SVD, which is numerically robust.
    %
    % Inputs:
    %   Ef_val     - Free E concentration (scalar)
    %   kPon_vals  - kPon values for each position (1 x N vector)
    %   kPoff_val  - kPoff rate constant (scalar)
    %   kEon_val   - kEon rate constant (scalar)
    %   kEoff_val  - kEoff rate constant (scalar)
    %   n          - Number of states = EBindingNumber + 1
    %
    % Outputs:
    %   avg_E_bound - Average number of E factors bound at each position (1 x N)

    N_positions = length(kPon_vals);
    avg_E_bound = zeros(1, N_positions);

    % Pre-compute the state index mapping
    % States are indexed in blocks by P-level:
    %   P=0: E=0                        -> state 1
    %   P=1: E=0, E=1                   -> states 2, 3
    %   P=2: E=0, E=1, E=2              -> states 4, 5, 6
    %   ...
    %   P=p: E=0, ..., E=p              -> states (p*(p+1)/2 + 1) to ((p+1)*(p+2)/2)
    total_states = n * (n + 1) / 2;

    % Build index arrays for E-count of each state
    E_count = zeros(total_states, 1);
    idx = 0;
    for p = 0:(n-1)
        for e = 0:p
            idx = idx + 1;
            E_count(idx) = e;
        end
    end

    for pos = 1:N_positions
        kPon = kPon_vals(pos);

        % Build numerical rate matrix
        A = build_numerical_rate_matrix(n, kPon, kPoff_val, kEon_val, kEoff_val, Ef_val);

        % Compute numerical null space via SVD
        % The null space of the rate matrix gives the steady-state distribution
        [~, S, V] = svd(A);
        s = diag(S);
        % Find the smallest singular value (should be ~0 for the null space)
        [~, min_idx] = min(abs(s));
        ss_dist = V(:, min_idx);

        % Ensure non-negative (flip sign if needed)
        if sum(ss_dist) < 0
            ss_dist = -ss_dist;
        end

        % Normalize to probability distribution
        ss_dist = ss_dist / sum(ss_dist);

        % Compute average E bound = sum(E_count .* probability)
        avg_E_bound(pos) = sum(E_count .* ss_dist);
    end
end


function A = build_numerical_rate_matrix(n, kPon, kPoff, kEon, kEoff, Ef)
    % BUILD_NUMERICAL_RATE_MATRIX Build the rate matrix numerically
    %
    % This replicates the logic of construct_rate_matrix.m but works with
    % numerical values directly, avoiding symbolic computation.

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

    % Set diagonal: A(i,i) = -sum of column i (excluding diagonal)
    for i = 1:size(A, 1)
        A(i, i) = -sum(A(:, i)) + A(i, i); % subtract existing diagonal if any
    end
end
