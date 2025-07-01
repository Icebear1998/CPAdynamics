function A = construct_rate_matrix(n)
    syms kP kPon kPoff kEon kEoff Ef real
    A = sym(zeros(n*n, n*n));
    sz = [n n]; 
    
    
    rows_kPon = repelem(1:n, fliplr(1:n)); % Repeat each number decreasingly
    cellArray1 = arrayfun(@(x) x:n, 1:n, 'UniformOutput', false); % Generate sequences
    cols_kPon = [cellArray1{:}]; % Concatenate horizontally
    ind_kPon = sub2ind(sz, rows_kPon, cols_kPon);
    %disp(ind_kPon);
    
    cellArray2 = arrayfun(@(x) 1:x, 1:n, 'UniformOutput', false); % Generate sequences
    rwos_kEon = [cellArray2{:}]; % Concatenate all sequences horizontally
    cols_kEon = repelem(1:n, 1:n); % Repeat each number i, i times
    ind_kEon = sub2ind(sz, rwos_kEon, cols_kEon);
    %disp(ind_kEon);
    
    % Loop through the array r
    for i = 1:length(ind_kPon)-2
        % Get the current and next index from r
        current_idx = ind_kPon(i);
        next_idx = ind_kPon(i+1);

         % Only assign values if the next index is valid (skip invalid pairs)
        if next_idx > current_idx
%             A(next_idx, current_idx) = 1; % 1 or kPon
%             A(current_idx, next_idx) = kP; % kP or kPoff
            A(next_idx, current_idx) = kPon; % 1 or kPon
            A(current_idx, next_idx) = kPoff; % kP or kPoff
        end
    end
    
    a1 = 1;
    a2 = 0;
    % Loop through the array r
    for i = 1:length(ind_kEon)-1
        % Get the current and next index from r
        current_idx = ind_kEon(i);
        next_idx = ind_kEon(i+1);

         % Only assign values if the next index is valid (skip invalid pairs)
        if next_idx == current_idx + 1
            A(next_idx, current_idx) = a2*kEon*Ef; % Assign to matrix B
            A(current_idx, next_idx) = a1*kEoff; % If you want symmetry
            a1 = a1 + 1;
            a2 = a2 - 1;
        else
            a2 = a1;
            a1 = 1;
        end
    end
    
    % Find rows and columns that are not all zeros (symbolically)
    row_mask = any(A ~= 0, 2);  % Check row-wise nonzero elements
    col_mask = any(A ~= 0, 1);  % Check column-wise nonzero elements

    % Convert logical masks using isAlways for symbolic variables
    row_mask = isAlways(row_mask);
    col_mask = isAlways(col_mask);

    % Apply the mask to remove all-zero rows and columns
    A = A(row_mask, col_mask);
    
    for i=1:size(A,1)
        A(i,i) = -sum(A(:,i));
    end
    
end
