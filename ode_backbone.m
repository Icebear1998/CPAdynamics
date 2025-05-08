function dxdt = ode_backbone(t, X, P)
global N PAS Ef_ss;

% Assign constant parameters every call
k_in = P.k_in; 
k_e  = P.k_e; 
k_e2 = P.k_e2;
RE_val_bind_E = P.RE_val_bind_E;
kHon_t_Ef  = P.kHon_t_Ef;
kHoff_t = P.kHoff;
kc_t    = P.kc;

% Declare E_f as persistent so it is assigned only once
persistent E_f_initialized E_f

if isempty(E_f_initialized)  % Runs only once
    E_f = P.E_total;         % Initialize Ef to E_total on the first call
    E_f_initialized = true;
end

% Extract state variables
R   = X(1:N);   
REH = X(N+1:end);

% Compute the used E, update E_f every call
E_used_sym = sum(R(1:PAS)' .* RE_val_bind_E);
E_used = double(subs(E_used_sym, 'Ef', E_f));  % Use the current Ef

E_f = (P.E_total - E_used);  % Update E_f dynamically in every call
%E_f = max(P.E_total - E_used, 0);  % Ensure free E does not become negative
Ef_ss = E_f;

Pol_f = (P.Pol_total - sum(R) - sum(REH));

% If E_f < 0, throw error and stop solver
if E_f < 0
    error('Negative E_f at t = %g (E_f = %g). Stopping simulation.', t, E_f);
end
%Update kHon_t based on the new E_f
kHon_t = kHon_t_Ef(E_f);

dxdt = zeros(length(X),1);

n = 1;
dxdt(n) = Pol_f*k_in - k_e*R(n);

for n = 2:(PAS-1)
    dxdt(n) = k_e*R(n-1) - k_e*R(n);
end

n = PAS;
j = n - PAS + 1;
dxdt(n) = k_e*R(n-1) - k_e*R(n) - kHon_t*R(n) + kHoff_t*REH(j);
dxdt(N+j) = -k_e2*REH(j) + kHon_t*R(n) - kHoff_t*REH(j);

for n = (PAS+1):N
    j = n - PAS + 1;
    dxdt(n) = k_e*R(n-1) - k_e*R(n) - kHon_t*R(n) + kHoff_t*REH(j);
    dxdt(N+j) = k_e2*REH(j-1) - k_e2*REH(j) + kHon_t*R(n) - kHoff_t*REH(j) - kc_t*REH(j);
end

end
