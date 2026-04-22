function A = build_rate_matrix_numerical(n, kPon, kPoff, kEon, kEoff, Ef)
% BUILD_RATE_MATRIX_NUMERICAL  Numerical rate matrix for E-binding / P-phosphorylation states
%
% Replicates the logic of construct_rate_matrix.m using numerical values
% directly, avoiding symbolic computation.
%
% Inputs:
%   n     - Number of states = EBindingNumber + 1
%   kPon  - Ser2P phosphorylation on-rate (scalar)
%   kPoff - Ser2P phosphorylation off-rate (scalar)
%   kEon  - E factor on-rate (scalar)
%   kEoff - E factor off-rate (scalar)
%   Ef    - Free E concentration (scalar)
%
% Output:
%   A - Rate matrix with all-zero rows/columns removed and diagonal set

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

% Set diagonal: A(i,i) = -sum of off-diagonal column entries
for i = 1:size(A, 1)
    A(i, i) = -sum(A(:, i)) + A(i, i);
end

end
