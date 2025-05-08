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

% Range of EBindingNumber to test
EBindingNumbers = 1:5;

% Logarithmic range for kPmin
kPmin_min = 0.001;  % Minimum kPmin value
kPmin_max = 0.1;    % Maximum kPmin value
num_kPmin = 3;      % Number of kPmin points
kPmin_range = logspace(log10(kPmin_min), log10(kPmin_max), num_kPmin);

% Arrays to store results
max_avg_E_bound = zeros(length(EBindingNumbers), 1);
optimal_kPmin = zeros(length(EBindingNumbers), 1);
optimal_Ef_ss = zeros(length(EBindingNumbers), 1);

% Loop over EBindingNumber
for idx_E = 1:length(EBindingNumbers)
    EBindingNumber = EBindingNumbers(idx_E);
    disp(['EBindingNumber = ', num2str(EBindingNumber)]);
    
    % Compute steady-state probabilities
    [r_E_BeforePas] = compute_steady_states(P, EBindingNumber + 1); 
    disp('done compute steady states');
    
    % Initialize maximum avg_E_bound for this EBindingNumber
    max_E_bound = -Inf;
    best_kPmin = NaN;
    best_Ef_ss = NaN;
    
    % Loop over kPmin to find maximum avg_E_bound
    for idx = 1:length(kPmin_range)
        P.kPmin = kPmin_range(idx);
        disp(['kPmin = ', num2str(P.kPmin)]);
        
        if P.kPmin >= P.kPmax
            continue; % Skip if kPmin >= kPmax
        end
        
        % Update Kp_vals for current kPmin
%         KpOver_vals = linspace(1/P.kPmax, 1/P.kPmin, PAS); % Range of Kp
%         Kp_vals = 1./KpOver_vals;
        Kp_vals = linspace(P.kPmax, P.kPmin, PAS); % Range of Kp
        % Recompute RE_vals for current kPmin
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
        
        % ------------ SOLVE THE STEADY STATE ------------
        X0 = zeros(N + N_PAS, 1);
        is_unphysical = false; % Reset flag
        X = fsolve(@(xx) ode_dynamics(xx, P), X0, optimoptions('fsolve', 'Display', 'off'));
        disp('done compute fsolve');
        
        % Check if solution is unphysical
        if is_unphysical
            continue;
        end
        
        % Extract solutions
        R_sol   = X(1:N);
        REH_sol = X(N+1 : N+N_PAS);
        
        % Compute avg_E_bound at PAS
        avg_E_bound = P.RE_val_bind_E(Ef_ss);
        disp('done compute RE_vals');
        
        % Update maximum avg_E_bound if current is higher
        if avg_E_bound(end) > max_E_bound
            max_E_bound = avg_E_bound(end);
            best_kPmin = P.kPmin;
            best_Ef_ss = Ef_ss;
        end
    end
    
    % Store the maximum avg_E_bound and optimal parameters for this EBindingNumber
    max_avg_E_bound(idx_E) = max_E_bound;
    optimal_kPmin(idx_E) = best_kPmin;
    optimal_Ef_ss(idx_E) = best_Ef_ss;
end

% ------------ PLOT RESULTS ------------
figure; hold on;
plot(EBindingNumbers, max_avg_E_bound, 'o-', 'LineWidth', 2, 'MarkerSize', 8, 'DisplayName', 'Max Avg E Bound');
grid on;
xlabel('Number of E Bindings (EBindingNumber)');
ylabel('Maximum Average E Bound at PAS');
title('Maximum Average E Bound vs. EBindingNumber');

% Add labels for optimal kPmin and Ef_ss
for idx_E = 1:length(EBindingNumbers)
    if ~isnan(optimal_kPmin(idx_E)) && ~isnan(max_avg_E_bound(idx_E))
        text(EBindingNumbers(idx_E), max_avg_E_bound(idx_E), ...
            sprintf('kPmin=%.3e, Ef_ss=%.1f', optimal_kPmin(idx_E), optimal_Ef_ss(idx_E)), ...
            'FontSize', 8, 'HorizontalAlignment', 'right', 'VerticalAlignment', 'bottom');
    end
end

legend('Location', 'best');
xlim([min(EBindingNumbers)-1, max(EBindingNumbers)+1]);
ylim([0, max(max_avg_E_bound, [], 'omitnan')*1.1]);
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
Ef_ss = E_f(end);
if Ef_ss < 0
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