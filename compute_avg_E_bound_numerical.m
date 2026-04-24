function [avg_E_bound, avg_Ser2P] = compute_avg_E_bound_numerical(Ef_val, kPon_vals, kPoff_val, kEon_val, kEoff_val, n)
    % COMPUTE_AVG_E_BOUND_NUMERICAL Compute average E bound and Ser2P at each position
    %
    % Builds the rate matrix numerically for each position and computes the
    % null space via SVD. Both avg_E_bound and avg_Ser2P are derived from the
    % same SVD so requesting both outputs adds no extra cost.
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
    %   avg_Ser2P   - Average Ser2P (P-level) at each position (1 x N)

    N_positions = length(kPon_vals);
    avg_E_bound = zeros(1, N_positions);
    avg_Ser2P   = zeros(1, N_positions);

    % Pre-compute per-state E and P counts (same indexing: block p contains E=0..p)
    total_states = n * (n + 1) / 2;
    E_count = zeros(total_states, 1);
    P_count = zeros(total_states, 1);
    idx = 0;
    for p = 0:(n-1)
        for e = 0:p
            idx = idx + 1;
            E_count(idx) = e;
            P_count(idx) = p;
        end
    end

    for pos = 1:N_positions
        kPon = kPon_vals(pos);

        A = build_rate_matrix_numerical(n, kPon, kPoff_val, kEon_val, kEoff_val, Ef_val);

        [~, S, V] = svd(A);
        s = diag(S);
        [~, min_idx] = min(abs(s));
        ss_dist = V(:, min_idx);

        if sum(ss_dist) < 0
            ss_dist = -ss_dist;
        end
        ss_dist = ss_dist / sum(ss_dist);

        avg_E_bound(pos) = sum(E_count .* ss_dist);
        avg_Ser2P(pos)   = sum(P_count .* ss_dist);
    end
end
