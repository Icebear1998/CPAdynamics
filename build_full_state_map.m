function S = build_full_state_map(M)
% BUILD_FULL_STATE_MAP Fixed (p,e) maps for the finite-rate R/RHE model.
% R: 0 <= e <= p <= M. RHE: 1 <= e <= p <= M (one E engaged with H).
% The state order matches Python: increasing p, then increasing e.
% Lookup arrays use (p+1,e+1), since biochemical counts start at zero.
validateattributes(M, {'numeric'}, {'scalar', 'real', 'finite', 'integer', 'positive'});
S.M = M;
S.sr = (M+1)*(M+2)/2;
S.sh = M*(M+1)/2;
S.R = zeros(S.sr, 2);
S.RHE = zeros(S.sh, 2);
S.R_index = zeros(M+1);
S.RHE_index = zeros(M+1);
r = 0;
h = 0;
for p = 0:M
    for e = 0:p
        r = r+1;
        S.R(r, :) = [p, e];
        S.R_index(p+1, e+1) = r;
        if e >= 1
            h = h+1;
            S.RHE(h, :) = [p, e];
            S.RHE_index(p+1, e+1) = h;
        end
    end
end
end
