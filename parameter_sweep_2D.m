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
P.kHon = 0.2;
P.kHoff = 0.0125;
P.kc = 0.05;
P.kPon_min = 0.01;
P.kPon_max = 1;
P.kPoff_min = 0.1;
P.kPoff_max = 2;
P.kPoff_const = 1;

L_a = 100;
geneLength_bp = 25000;
PASposition = 20000;
EBindingNumber = 1;

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
param_pairs = {{'E_total', 'kHoff'}};

for pair_idx = 1:length(param_pairs)
    param1 = param_pairs{pair_idx}{1};
    param2 = param_pairs{pair_idx}{2};
    
    param1_values = 40000:-10000:20000;
    param2_values = 0.06:-0.01:0.02;

    cutoff_matrix = zeros(length(param2_values), length(param1_values));
    options = optimoptions('fsolve', 'Algorithm', 'trust-region', 'Display', 'none', 'FunctionTolerance', 1e-10);

    % --- 2D Parameter Sweep Loop ---
    for i = 1:length(param2_values)
        
        % --- NEW: Initialize history buffer for each row of the sweep ---
        X_history = cell(1, length(param1_values));
        
        for j = 1:length(param1_values)
            fprintf('Running sweep: %s = %.2g, %s = %.2g\n', param2, param2_values(i), param1, param1_values(j));
            
            % --- NEW: EXTRAPOLATION-BASED GUESS LOGIC ---
            if j == 1
                % For the first point, use a default guess
                X_guess = 1e-2 * ones(N + N_PAS, 1);
            elseif j == 2
                % For the second point, use the previous solution (simple continuation)
                X_guess = X_history{j-1};
            else
                % For all subsequent points, use linear extrapolation
                X1 = X_history{j-2}; % Penultimate solution
                X2 = X_history{j-1}; % Last solution
                p1 = param1_values(j-2);
                p2 = param1_values(j-1);
                p_next = param1_values(j);
                
                % Linear extrapolation formula: X_guess = X2 + (X2-X1)*(p_next-p2)/(p2-p1)
                X_guess = X2 + (X2 - X1) * ((p_next - p2) / (p2 - p1));
                
                % Ensure the guess remains non-negative
                X_guess(X_guess < 0) = 0;
            end

            P_run = P;
            P_run.(param1) = param1_values(j);
            P_run.(param2) = param2_values(i);
            kHon_original = P_run.kHon;
            
            % --- DIRECT TWO-STEP SOLVER ---
            P_run.is_unphysical = false;
            
            P_run.FirstRun = true; Ef_ss = 0;
            try
                X_base = fsolve(@(xx) ode_dynamics_multipleE(xx, P_run), X_guess, options);
            catch
                P_run.is_unphysical = true;
            end
            
            if P_run.is_unphysical || any(isnan(X_base))
                cutoff_matrix(i,j) = NaN; continue;
            end

            avg_E_bound = P_run.RE_val_bind_E(Ef_ss);
            P_run.FirstRun = false;
            P_run.kHon = kHon_original * avg_E_bound(end);
            X_final = fsolve(@(xx) ode_dynamics_multipleE(xx, P_run), X_base, options);
            
            % --- Store result in history buffer ---
            X_history{j} = X_final;
            
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
    
    % --- SAVE RESULTS ---
    % Prepare data structure for saving
    data.EBindingNumber = EBindingNumber;
    data.param1 = param1;
    data.param2 = param2;
    data.param1_values = param1_values;
    data.param2_values = param2_values;
    data.cutoff_matrix = cutoff_matrix;
    
    % Save results using the utility function
    save_analysis_results('parameter_sweep_2D', data, P);
end

% --- Helper Functions (ode_dynamics_multipleE, compute_steady_states, etc. would be here) ---