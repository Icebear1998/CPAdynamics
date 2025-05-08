global N PAS N_PAS Ef_ss kHon_ss is_unphysical;

% ------------ MODEL PARAMETERS ------------
L_a = 100;
P.k_in    = 2;
P.kEon    = 0.00025;
P.kEoff   = 10;
P.k_e     = 65/L_a;
P.k_e2    = 30/L_a;
P.E_total = 70000;
P.L_total = 100000;
P.Pol_total = 70000;
P.kHon    = 0.1; % based on typical k bind and estimated J factor for H.
P.kHoff   = 0.0025; 
P.kc      = 0.8; % not sure
P.kPmax   = 40; % Fixed value

% Gene setup
geneLength_bp = 25000;
PASposition   = 20000;
N      = floor(geneLength_bp / L_a);  % total nodes
PAS    = floor(PASposition / L_a);    % node index of PAS
N_PAS  = N - PAS + 1;                 % number of nodes at/after PAS

EBindingNumber = 3;
[r_E_BeforePas] = compute_steady_states(P, EBindingNumber + 1); 
disp('done compute steady states');

% Linear range for kPmin
kPmin_min = 0.1;  % Minimum kPmin value
kPmin_max = 8;   % Maximum kPmin value
num_kPmin = 3;   % Number of kPmin points
kPmin_range = linspace(kPmin_min, kPmin_max, num_kPmin);

Kp_vals = linspace(P.kPmax, P.kPmin, PAS); % Range of Kp
RE_vals = sym(zeros(EBindingNumber + 1, PAS));

for e = 1:EBindingNumber+1
    for i = 1:length(Kp_vals)
        kP_val = Kp_vals(i);
        RE_vals(e, i) = subs(r_E_BeforePas(e), {'kP'}, {kP_val});
    end
end
disp('done compute EBindingNumber');

% Efficient conversion to function
P.RE_val_bind_E = matlabFunction(simplify(sum((1:EBindingNumber)' .* RE_vals(2:end, :), 1)), 'Vars', {Ef});
disp('done compute RE_val_bind_E');

% Arrays to store results
avg_E_bound_values = zeros(length(kPmin_range), 1);
Ef_ss_values = zeros(length(kPmin_range), 1);

% Loop over kPmin values
for idx = 1:length(kPmin_range)
    P.kPmin = kPmin_range(idx);
    disp(['kPmin = ', num2str(P.kPmin)]);
    
    if P.kPmin >= P.kPmax
        avg_E_bound_values(idx) = NaN;
        Ef_ss_values(idx) = NaN;
        continue; % Skip if kPmin >= kPmax
    end
    
    % Update Kp_vals for current kPmin
    Kp_vals = linspace(P.kPmax, P.kPmin, PAS);
    
    % Recompute RE_vals for current kPmin
    for e = 1:EBindingNumber+1
        for i = 1:length(Kp_vals)
            kP_val = Kp_vals(i);
            RE_vals(e, i) = subs(r_E_BeforePas(e), {'kP'}, {kP_val});
        end
    end
    
    % ------------ SOLVE THE STEADY STATE ------------
    X0 = zeros(N + N_PAS, 1);
    is_unphysical = false; % Reset flag
    X = fsolve(@(xx) ode_dynamics(xx, P), X0, optimoptions('fsolve', 'Display', 'off'));
    disp('done compute fsolve');
    
    % Extract solutions
    R_sol   = X(1:N);
    REH_sol = X(N+1 : N+N_PAS);
    
    % Compute avg_E_bound at PAS
    avg_E_bound = zeros(1, PAS);
    for e = 1:EBindingNumber+1
        for i = 1:PAS
            RE_vals(e, i) = double(R_sol(i) * double(subs(RE_vals(e, i), {'Ef'}, {Ef_ss})));
            % Compute the average number of E molecules bound at each position
            if R_sol(i) > 0
                avg_E_bound(i) = avg_E_bound(i) + (e-1) * (RE_vals(e, i) / R_sol(i));
            end
        end
    end
    disp('done compute RE_vals');
    
    % Store results
    avg_E_bound_values(idx) = avg_E_bound(end); % At PAS
    Ef_ss_values(idx) = Ef_ss;
end

% ------------ PLOT RESULTS ------------
figure; hold on;
plot(kPmin_range, avg_E_bound_values, 'o-', 'LineWidth', 2, 'MarkerSize', 8, 'DisplayName', 'Avg E Bound');
grid on;
xlabel('kPmin');
ylabel('Average E Bound at PAS');
title('Average E Bound at PAS vs. kPmin (kPmax = 20)');
legend('Location', 'best');
xlim([min(kPmin_range)*0.5, max(kPmin_range)*2]);
ylim([0, max(avg_E_bound_values, [], 'omitnan')*1.1]);
hold off;

% ------------ ODE DYNAMICS FUNCTION ------------
function dxdt = ode_dynamics(X, P)
global N PAS Ef_ss kHon_ss is_unphysical;

k_in   = P.k_in;
k_e    = P.k_e;
k_e2   = P.k_e2;
kHoff_t= P.kHoff;
kc_t   = P.kc;
RE_val_bind_E = P.RE_val_bind_E;
kHon_t = P.kHon;

R   = X(1:N);
REH = X(N+1:end);

% Compute E_f numerically
E_used = sum(R(1:PAS)' .* RE_val_bind_E(Ef_ss));
E_f = P.E_total - E_used;
Ef_ss = E_f;
if E_f < 0
    is_unphysical = true; % Set flag instead of throwing error
    dxdt = 1e6 * ones(length(X), 1); % Large residual to signal fsolve to stop
    return;
end

Pol_f = P.Pol_total - sum(R) - sum(REH);
dxdt = zeros(length(X), 1);

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