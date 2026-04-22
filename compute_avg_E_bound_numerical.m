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
        A = build_rate_matrix_numerical(n, kPon, kPoff_val, kEon_val, kEoff_val, Ef_val);

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
