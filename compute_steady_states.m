function [r_E_BeforePas, r_P] = compute_steady_states(P, n)
    % Define symbolic variables
    syms kP Ef kEon kEoff kHon kHoff kc
   
    
    % Construct transition matrices symbolically
    A_before = construct_rate_matrix(n);
    %A_after = construct_rate_matrix3D(n);
    
    A_before = subs(A_before, [kEon, kEoff], [P.kEon, P.kEoff]);
    %A_after = subs(A_after, [kP, kEon, kEoff, kHon, kHoff, kc], [P.kPmax, P.kEon, P.kEoff, P.kHon, P.kHoff, P.kc]);

    % Compute steady-state distributions symbolically
    r_Before = null(A_before);  
    %r_After = null(A_after);
    
    % Compute cumulative sums of steady-state probabilities symbolically
    r_E_BeforePas = compute_cumulative_sum(r_Before, n);
    r_P = compute_block_sums(r_Before, n);

    % Normalize the steady-state probabilities
    r_E_BeforePas = simplify(r_E_BeforePas / sum(r_E_BeforePas));
    r_P = simplify(r_P / sum(r_P));

end

function r_E = compute_cumulative_sum(r, n)
    RE_index = cumsum(0:n-1) + 1;
    r_E = sym(zeros(1, n)); % Initialize symbolic array
    
    for i = 1:n
        RE_index = RE_index(RE_index <= length(r));
        RE_index_copy = RE_index(i:end);
        r_E(i) = sum(r(RE_index_copy));
        RE_index = RE_index + 1;
    end
end

function r_P = compute_block_sums(r, n)
    block_lengths = 1:n; % Block sizes: 1, 2, 3, ..., n
    start_indices = cumsum([1 block_lengths(1:end-1)]); % Start indices: 1, 2, 4, 7, ...

    r_P = sym(zeros(1, n));

    % Compute the sum for each block
    for i = 1:n
        % Determine the range of indices for the current block
        block_start = start_indices(i);
        block_end = block_start + block_lengths(i) - 1;
        
        % Ensure we don’t exceed the length of r
        if block_start > length(r)
            r_P(i) = 0; % If there are no elements to sum, set to 0
            continue;
        end
        block_end = min(block_end, length(r)); % Cap at the length of r
        
        % Sum the elements in the current block
        block_indices = block_start:block_end;
        r_P(i) = sum(r(block_indices));
    end
end

