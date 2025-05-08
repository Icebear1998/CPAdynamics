function A = construct_rate_matrix3D(n)
    syms Ef kP kPon kPoff kEon kEoff kHon kHoff kc real
    A = sym(zeros(n*n + n*n, n*n + n*n));
    sz = [n n 2]; 
    
%     kPon = 1;
%     kPoff = -1;
%     kEon = 2;
%     kEoff = -2;
   
    rows_kPon = repelem(1:n, fliplr(1:n)); % Repeat each number decreasingly
    cellArray1 = arrayfun(@(x) x:n, 1:n, 'UniformOutput', false); % Generate sequences
    cols_kPon = [cellArray1{:}]; % Concatenate horizontally
    pages_kPon = ones(1,length(rows_kPon));
    ind_kPon = sub2ind(sz, rows_kPon, cols_kPon, pages_kPon);
    
    cellArray2 = arrayfun(@(x) 1:x, 1:n, 'UniformOutput', false); % Generate sequences
    rwos_kEon = [cellArray2{:}]; % Concatenate all sequences horizontally
    cols_kEon = repelem(1:n, 1:n); % Repeat each number i, i times
    pages_kEon = ones(1,length(rwos_kEon));
    ind_kEon = sub2ind(sz, rwos_kEon, cols_kEon, pages_kEon);
    %disp(ind_kEon);
    
    % index to remove
    
    
    shift_index = n*n;
    % Loop through the array r
    for i = 1:length(ind_kPon)-2
        % Get the current and next index from r
        current_idx = ind_kPon(i);
        next_idx = ind_kPon(i+1);

         % Only assign values if the next index is valid (skip invalid pairs)
        if next_idx > current_idx
            A(next_idx, current_idx) = 1; % Assign to matrix B
            A(current_idx, next_idx) = kP; % If you want symmetry
            
            if (i >= n+1)
                A(next_idx+shift_index, current_idx+shift_index) = 1; % Assign to matrix B
                A(current_idx+shift_index, next_idx+shift_index) = kP; % If you want symmetry
            end
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
           
            if (~ismember(current_idx, ind_kPon(1:n+1)))
                A(next_idx+shift_index, current_idx+shift_index) = a2*kEon*Ef; % Assign to matrix B
                A(current_idx+shift_index, next_idx+shift_index) = a1*kEoff; % If you want symmetry
            end
            
            a1 = a1 + 1;
            a2 = a2 - 1;
        else
            a2 = a1;
            a1 = 1;
        end
    end
    
    b = 1;
    for i = n+1:length(ind_kPon)
        current_idx = ind_kPon(i);
        next_idx = ind_kPon(min(i+1, length(ind_kPon)));
        
        if (i == length(ind_kPon))
            b = b+1;
        end
        
        A(current_idx+shift_index, current_idx) = b*kHon; % Assign to matrix B
        A(current_idx, current_idx+shift_index) = kHoff; % If you want symmetry
        
        if (current_idx > next_idx)
            b = b + 1;
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
    
    for i=1:(n*n)
        A(i,i) = -sum(A(:,i));
              
%         if (i > sum(1:n))
%            A(i,i) =  A(i,i) - kc;
%         end
    end
end
