% --- Model parameters (base values) ---
P = struct();
P.k_in = 2;
P.kEon = 0.00025;
P.kEoff = 10;
P.k_e = 65/100;
P.k_e2 = 30/100;
P.E_total = 70000;
P.L_total = 100000;
P.Pol_total = 70000;
P.kHon = 0.05;
P.kHoff = 0.0025;
P.kc = 0.8;
P.kPon_min = 0.01;
P.kPon_max = 1;
P.kPoff_min = 0.1;
P.kPoff_max = 2;
P.kPoff_const = 1;

L_a = 100;
geneLength_bp = 25000;
PASposition = 20000;
EBindingNumber = 3;

% --- Global variables for ODE function ---
global N PAS N_PAS Ef_ss;
N = floor(geneLength_bp / L_a);
PAS = floor(PASposition / L_a);
N_PAS = N - PAS + 1;

% --- OPTIMIZATION: Perform Symbolic Pre-computation ONCE ---
syms Ef real;
fprintf('Performing one-time symbolic pre-computation...\n');
[r_E_BeforePas] = compute_steady_states(P, EBindingNumber+1);
kPon_vals = linspace(P.kPon_min, P.kPon_max, PAS);
RE_vals = sym(zeros(EBindingNumber + 1, N));
for e = 1:EBindingNumber + 1
    for idx = 1:PAS
        kPon_val = kPon_vals(idx);
        RE_vals(e, idx) = subs(r_E_BeforePas(e), {'kPon', 'kPoff'}, {kPon_val, P.kPoff_const});
    end
    for idx = PAS+1:N
        RE_vals(e, idx) = subs(r_E_BeforePas(e), {'kPon', 'kPoff'}, {P.kPon_max, P.kPoff_min});
    end
end
P.RE_val_bind_E = matlabFunction(simplify(sum(sym(1:EBindingNumber)' .* RE_vals(2:end, :), 1)), 'Vars', {Ef});
fprintf('Pre-computation complete.\n\n');

% --- Parameter Sweep Setup ---
param_pairs = {{'E_total', 'kHon'}};

for pair_idx = 1:length(param_pairs)
    param1 = param_pairs{pair_idx}{1};
    param2 = param_pairs{pair_idx}{2};
    
    param1_values = 40000:-10000:20000;
    param2_values = 0.05:-0.01:0.02;

    cutoff_matrix = zeros(length(param2_values), length(param1_values));

    % --- SOLVER OPTIONS for ode15s ---
    % We set tolerances and enforce that all concentrations are non-negative.
    ode_options = odeset('RelTol', 1e-6, 'AbsTol', 1e-8, 'NonNegative', 1:(N + N_PAS));
    t_span = [0 1e6]; % Simulate for a long time to ensure steady state is reached

    % --- 2D Parameter Sweep Loop ---
    for i = 1:length(param2_values)
        X_guess = 1e-2 * ones(N + N_PAS, 1);
        for j = 1:length(param1_values)
            fprintf('Running sweep: %s = %.2g, %s = %.2g\n', param2, param2_values(i), param1, param1_values(j));
            
            P_run = P;
            P_run.(param1) = param1_values(j);
            P_run.(param2) = param2_values(i);
            kHon_original = P_run.kHon;
            
            % --- STABLE TWO-STEP SOLVER (ode15s) ---
            X_final = NaN(size(X_guess)); % Default to NaN
            try
                % Step 1: Find the steady state for the BASE kHon value.
                P_run.FirstRun = true;
                Ef_ss = 0; % Reset Ef_ss for the self-consistency calculation
                [~, X_t_base] = ode15s(@(t,X) ode_dynamics_multipleE(t, X, P_run), t_span, X_guess, ode_options);
                X_base = X_t_base(end, :)'; % The steady state is the last time point

                % Step 2: Update kHon and solve for the final steady state.
                avg_E_bound = P_run.RE_val_bind_E(Ef_ss);
                P_run.FirstRun = false; 
                P_run.kHon = kHon_original * avg_E_bound(end);
                
                % Use the result of the first run as the initial condition for the second
                [~, X_t_final] = ode15s(@(t,X) ode_dynamics_multipleE(t, X, P_run), t_span, X_base, ode_options);
                X_final = X_t_final(end, :)';

            catch ME
                warning('ode15s failed for point (%d,%d): %s', i, j, ME.message);
            end

            X_guess = X_final;
            
            % --- Process and Store Result ---
            if any(isnan(X_final)); cutoff_matrix(i,j) = NaN; continue; end
            
            R_sol = X_final(1:N);
            REH_sol = X_final(N+1 : N+N_PAS);
            
            ratio = (REH_sol(1:end) + R_sol(PAS:end)) / (R_sol(PAS-1) + 1e-9);
            node_indices = 1:(length(ratio));
            first_idx_below = find(ratio <= 0.85, 1, 'first');
            
            if isempty(first_idx_below); cutoff_node = length(node_indices);
            elseif first_idx_below == 1; cutoff_node = interp1([1, ratio(1)], [0, 1], 0.85);
            else; cutoff_node = interp1([ratio(first_idx_below-1), ratio(first_idx_below)], [node_indices(first_idx_below-1), node_indices(first_idx_below)], 0.85);
            end
            if isnan(cutoff_node); cutoff_node = length(node_indices); end
            cutoff_matrix(i,j) = cutoff_node * L_a;
        end
    end

    % --- Plotting ---
    figure('Position', [100, 100, 800, 600]);
    hold on;
    colors = lines(size(cutoff_matrix, 1));
    for i = 1:size(cutoff_matrix, 1)
        plot(param1_values, cutoff_matrix(i, :), 'o-', 'LineWidth', 2, ...
             'Color', colors(i, :), 'DisplayName', sprintf('%s = %.2g', param2, param2_values(i)));
    end
    xlabel([strrep(param1, '_', '\_'), ' Value'], 'FontSize', 12);
    ylabel('CPA Cutoff Position (bp)', 'FontSize', 12);
    title(['CPA Cutoff vs ', strrep(param1, '_', '\_'), ' by ', strrep(param2, '_', '\_')], 'FontSize', 14);
    grid on; legend('show', 'Location', 'best'); set(gca, 'FontSize', 10); box on;
    hold off;
end

%% --- Helper Functions ---

function dxdt = ode_dynamics_multipleE(t, X, P)
    % MODIFIED to accept the time argument 't', as required by ODE solvers.
    global N PAS N_PAS Ef_ss;
    
    k_in = P.k_in; k_e = P.k_e; k_e2 = P.k_e2;
    kHoff = P.kHoff; kc = P.kc; kHon = P.kHon;
    RE_val_bind_E = P.RE_val_bind_E;

    R   = X(1:N);
    REH = X(N+1:end);

    if P.FirstRun
        if Ef_ss == 0; Ef_ss = P.E_total; end
        E_used = sum(R(1:N)' .* RE_val_bind_E(Ef_ss))+ sum(REH, 1);% + sum(REH' .* RE_val_bind_E(Ef_ss)(PAS:N));
        Ef_ss = P.E_total - E_used;
    end

    Pol_f = P.Pol_total - sum(R) - sum(REH);
    dxdt = zeros(length(X),1);
    
    dxdt(1) = Pol_f*k_in - k_e*R(1);
    for n = 2:(PAS-1)
        dxdt(n) = k_e*R(n-1) - k_e*R(n);
    end

    n = PAS; j = 1;
    dxdt(n) = k_e*R(n-1) - k_e*R(n) - kHon*R(n) + kHoff*REH(j);
    dxdt(N+j) = -k_e2*REH(j) + kHon*R(n) - kHoff*REH(j) - kc*REH(j);

    for n = (PAS+1):N
        j = n - PAS + 1;
        dxdt(n) = k_e*R(n-1) - k_e*R(n) - kHon*R(n) + kHoff*REH(j);
        dxdt(N+j) = k_e2*REH(j-1) - k_e2*REH(j) + kHon*R(n) - kHoff*REH(j) - kc*REH(j);
    end
end