function [A_const, A_Pon, A_Ef] = build_full_internal_rate_matrices(P, S, recognized)
% BUILD_FULL_INTERNAL_RATE_MATRICES Local column-oriented reaction generators.
% A = A_const + kPon*A_Pon + E_free*A_Ef. Columns sum to zero.
% Phosphorylation has no site multiplicity, matching the original MATLAB model.
% An RHE state's e-1 unengaged E factors use kEoff. Its engaged-E detachment
% changes RHE to R and is added separately by build_full_rate_matrices.
if recognized
    states = S.RHE;
    lookup = S.RHE_index;
else
    states = S.R;
    lookup = S.R_index;
end
count = size(states, 1);
A_const = zeros(count);
A_Pon = zeros(count);
A_Ef = zeros(count);
for j = 1:count
    p = states(j, 1);
    e = states(j, 2);
    if p < S.M
        i = lookup(p+2, e+1);
        A_Pon(i, j) = A_Pon(i, j)+1;
        A_Pon(j, j) = A_Pon(j, j)-1;
    end
    if p > e
        i = lookup(p, e+1);
        A_const(i, j) = A_const(i, j)+P.kPoff;
        A_const(j, j) = A_const(j, j)-P.kPoff;
        i = lookup(p+1, e+2);
        rate = (p-e)*P.kEon;
        A_Ef(i, j) = A_Ef(i, j)+rate;
        A_Ef(j, j) = A_Ef(j, j)-rate;
    end
    exchangeable = e-double(recognized);
    if exchangeable > 0
        i = lookup(p+1, e);
        rate = exchangeable*P.kEoff;
        A_const(i, j) = A_const(i, j)+rate;
        A_const(j, j) = A_const(j, j)-rate;
    end
end
A_const = sparse(A_const);
A_Pon = sparse(A_Pon);
A_Ef = sparse(A_Ef);
end
