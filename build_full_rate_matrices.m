function [model, P] = build_full_rate_matrices(P, EBindingNumber)
% BUILD_FULL_RATE_MATRICES Sparse finite-rate R/RHE transport and reactions.
% x contains all R nodes, then all RHE nodes, with each node's (p,e) states
% contiguous. A(E_free) = model.A0 + E_free*model.AE acts on x.
% No local equilibrium, symbolic calculations or global variables are used.
% kEoff_engaged defaults to zero if absent (the earlier full-model limit).
if ~isfield(P, 'kEoff_engaged')
    P.kEoff_engaged = 0;
end
rateFields = {'k_in', 'k_e', 'k_e2', 'kEon', 'kEoff', 'kHon', 'kHoff', ...
              'kc', 'kPon_min', 'kPon_slope', 'kPoff', 'kEoff_engaged', ...
              'E_total', 'Pol_total'};
for k = 1:numel(rateFields)
    name = rateFields{k};
    validateattributes(P.(name), {'numeric'}, {'scalar', 'real', 'finite', 'nonnegative'}, mfilename, name);
end
validateattributes(P.L_a, {'numeric'}, {'scalar', 'real', 'finite', 'positive'});
validateattributes(P.geneLength_bp, {'numeric'}, {'scalar', 'real', 'finite', 'positive'});
validateattributes(P.PASposition, {'numeric'}, {'scalar', 'real', 'finite', 'positive'});
P.N = floor(P.geneLength_bp/P.L_a);
P.PAS = floor(P.PASposition/P.L_a); % MATLAB 1-based PAS node
P.N_PAS = P.N-P.PAS+1;
if P.N < 1 || P.PAS < 1 || P.PAS > P.N
    error('build_full_rate_matrices:InvalidGeometry', 'PAS must be a node inside the simulated gene.');
end
S = build_full_state_map(EBindingNumber);
model.states = S;
model.N = P.N;
model.PAS = P.PAS;
model.N_PAS = P.N_PAS;
model.sr = S.sr;
model.sh = S.sh;
model.r_size = P.N*S.sr;
model.size = model.r_size+P.N_PAS*S.sh;
model.k_in = P.k_in;
model.e_counts = [repmat(S.R(:, 2), P.N, 1); repmat(S.RHE(:, 2), P.N_PAS, 1)];

[R_const, R_Pon, R_Ef] = build_full_internal_rate_matrices(P, S, false);
[H_const, H_Pon, H_Ef] = build_full_internal_rate_matrices(P, S, true);
kPon = P.kPon_min+P.kPon_slope*(0:P.N-1)';
T_R = transport_matrix(P.N, P.k_e);
T_H = transport_matrix(P.N_PAS, P.k_e2);
A_R = kron(speye(P.N), R_const) + kron(spdiags(kPon, 0, P.N, P.N), R_Pon) ...
      + kron(T_R, speye(S.sr));
A_H = kron(speye(P.N_PAS), H_const) ...
      + kron(spdiags(kPon(P.PAS:end), 0, P.N_PAS, P.N_PAS), H_Pon) ...
      + kron(T_H, speye(S.sh)) - P.kc*speye(P.N_PAS*S.sh);
model.A0 = blkdiag(A_R, A_H);
model.AE = blkdiag(kron(speye(P.N), R_Ef), kron(speye(P.N_PAS), H_Ef));

% Each RHE state has three cross-population reactions (two entries each).
count = 6*P.N_PAS*S.sh;
rows = zeros(count, 1);
cols = zeros(count, 1);
rates = zeros(count, 1);
offset = 0;
for node = P.PAS:P.N
    for hLocal = 1:S.sh
        pe = S.RHE(hLocal, :);
        r = (node-1)*S.sr+S.R_index(pe(1)+1, pe(2)+1);
        h = model.r_size+(node-P.PAS)*S.sh+hLocal;
        rAfterLoss = (node-1)*S.sr+S.R_index(pe(1)+1, pe(2));
        take = offset+(1:6);
        rows(take) = [h; r; r; h; rAfterLoss; h];
        cols(take) = [r; r; h; h; h; h];
        % Exactly one E is H-engaged: its off-rate has no e multiplicity.
        rates(take) = [pe(2)*P.kHon; -pe(2)*P.kHon; P.kHoff; -P.kHoff; ...
                       P.kEoff_engaged; -P.kEoff_engaged];
        offset = offset+6;
    end
end
model.A0 = model.A0+sparse(rows, cols, rates, model.size, model.size);
model.unit_source = zeros(model.size, 1);
model.unit_source(1) = 1; % newly initiated R(p=0,e=0)
model.initiation = P.k_in*model.unit_source;
weights = sparse(model.e_counts');
model.E_release0 = -weights*model.A0;
model.E_release1 = -weights*model.AE;
model.Pol_release0 = -sum(model.A0, 1);
model.Pol_release1 = -sum(model.AE, 1);
end

function T = transport_matrix(nodes, rate)
T = -rate*speye(nodes);
if nodes > 1
    T = T+sparse(2:nodes, 1:nodes-1, rate*ones(1, nodes-1), nodes, nodes);
end
end

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

function [A_const, A_Pon, A_Ef] = build_full_internal_rate_matrices(P, S, recognized)
% BUILD_FULL_INTERNAL_RATE_MATRICES Local column-oriented reaction generators.
% A = A_const + kPon*A_Pon + E_free*A_Ef. Columns sum to zero.
% Phosphorylation has no site multiplicity.
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
