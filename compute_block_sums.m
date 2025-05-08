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