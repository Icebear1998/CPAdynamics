% --- CONFIGURATION: CHOOSE THE PARAMETER TO SWEEP ---
sweep_param_name = 'kHoff';

switch sweep_param_name
    case 'E_total'
        sweep_param_values = [10000, 30000, 70000, 100000, 120000, 150000];
    case 'kHoff'
        sweep_param_values = [0.001, 0.01, 0.1, 1, 5, 10];
    otherwise
        error('Selected sweep parameter is not defined.');
end

% --- BASE PARAMETERS ---
P.k_in = 2; P.k_e = 65/100; P.k_e2 = 30/100; P.E_total = 70000;
P.Pol_total = 70000; P.kEon = 0.00025; P.kEoff = 10; P.kHon = 0.2;
P.kHoff = 0.0125; P.kc = 0.05; P.kPon_min = 0.01; P.kPon_max = 1;
P.kPoff_const = 1; P.kPoff_max = 2; P.kPoff_min = 0.1; P.L_a = 100;
P.geneLength_bp = 25000; P.PASposition = 20000; EBindingNumber = 1;

% --- SIMULATION SETUP ---
inter_pas_distances_bp = 0:100:2500;
proximal_usage_results_cdf = zeros(length(inter_pas_distances_bp), length(sweep_param_values));

% --- Start a parallel pool of workers ---
if isempty(gcp('nocreate')); parpool; end
p = gcp; % Get the current parallel pool object
fprintf('Parallel pool started with %d workers...\n', p.NumWorkers);


% --- PARALLEL PARAMETER SWEEP LOOP ---
fprintf('Starting parallel sweep over parameter: %s\n', sweep_param_name);
parfor p_idx = 1:length(sweep_param_values)
    P_run = P; 
    P_run.(sweep_param_name) = sweep_param_values(p_idx);
    
    [R_sol, REH_sol, P_sim] = run_termination_simulation_parallel(P_run, EBindingNumber);

    % --- Flux-based CDF calculation ---
    flux_cleavage_per_node = P_sim.kc * REH_sol;
    total_cleavage_flux    = sum(flux_cleavage_per_node);
    
    if total_cleavage_flux > 1e-9
        cumulative_cleavage_flux = cumsum(flux_cleavage_per_node);
        termination_cdf = cumulative_cleavage_flux / total_cleavage_flux;
    else
        termination_cdf = zeros(size(REH_sol));
    end

    nodes_post_pas = 1:length(REH_sol);
    bp_post_pas = nodes_post_pas * P_sim.L_a;
    bp_for_interp = [0; bp_post_pas(:)];
    cdf_for_interp = [0; termination_cdf(:)];
    
    proximal_usage_results_cdf(:, p_idx) = interp1(bp_for_interp, cdf_for_interp, inter_pas_distances_bp, 'linear', 'extrap');
end
disp('All simulations complete.');

% --- PLOT THE RESULTS ---
figure('Position', [100, 100, 800, 600]);
hold on;
colors = lines(length(sweep_param_values));

for p_idx = 1:length(sweep_param_values)
    plot(inter_pas_distances_bp, proximal_usage_results_cdf(:, p_idx) * 100, ...
         'LineWidth', 2.5, 'Color', colors(p_idx,:), ...
         'DisplayName', sprintf('%s = %.3g', strrep(sweep_param_name, '_', '\_'), sweep_param_values(p_idx)));
end

plot_title = sprintf('CDF of Proximal PAS Usage vs. %s (EBindingNumber=%d)', strrep(sweep_param_name, '_', ' '), EBindingNumber);
xlabel('Inter-PAS Distance (bp)', 'FontSize', 12);
ylabel('Cumulative Proximal Site Usage (%)', 'FontSize', 12);
title(plot_title, 'FontSize', 14, 'FontWeight', 'bold');
grid on; legend('show', 'Location', 'best'); set(gca, 'FontSize', 10); box on;
ylim([0 100]);

function [R_sol, REH_sol, P] = run_termination_simulation_parallel(P, EBindingNumber)
    % This function is now self-contained and safe for parallel execution.
    syms Ef real;
    
    % Calculate node information and store it in the P struct
    P.N = floor(P.geneLength_bp / P.L_a);
    P.PAS = floor(P.PASposition / P.L_a);
    P.N_PAS = P.N - P.PAS + 1;
    
    kHon_base = P.kHon;
    
    % --- Symbolic Pre-computation ---
    [r_E_BeforePas] = compute_steady_states(P, EBindingNumber + 1);
    kPon_vals = linspace(P.kPon_min, P.kPon_max, P.PAS);
    RE_vals = sym(zeros(EBindingNumber + 1, P.N));
    for e = 1:EBindingNumber + 1
        for idx = 1:length(kPon_vals)
            RE_vals(e, idx) = subs(r_E_BeforePas(e), {'kPon', 'kPoff'}, {kPon_vals(idx), P.kPoff_const});
        end
        for idx = P.PAS+1:P.N
            RE_vals(e, idx) = subs(r_E_BeforePas(e), {'kPon', 'kPoff'}, {P.kPon_max, P.kPoff_min});
        end
    end
    P.RE_val_bind_E = matlabFunction(simplify(sum(sym(1:EBindingNumber)' .* RE_vals(2:end, :), 1)), 'Vars', {Ef});

    % --- STABLE TWO-STEP SOLVER ---
    X_guess = 1e-6 * ones(P.N + P.N_PAS, 1);
    options = optimoptions('fsolve', 'Display', 'off', 'FunctionTolerance', 1e-8);
    
    % Create an anonymous function for fsolve that passes P along
    ode_handle = @(X, P_inst) ode_dynamics_multipleE_parallel(X, P_inst);

    % Step 1: Solve for the initial steady state
    P.FirstRun = true;
    P.Ef_ss_sol = 0;
    X_base = fsolve(@(X) ode_handle(X, P), X_guess, options);
%     try
%         X_base = fsolve(@(X) ode_handle(X, P), X_guess, options);
%     catch
%         error('Solver failed in Step 1.');
%     end

    % Step 2: Update kHon and solve for the final steady state
    avg_E_bound = P.RE_val_bind_E(P.Ef_ss_sol);
    P.FirstRun = false;
    P.kHon = kHon_base * avg_E_bound(end);
    
    X_final = fsolve(@(X) ode_handle(X, P), X_base, options);
    
    R_sol = X_final(1:P.N);
    REH_sol = X_final(P.N+1 : P.N+P.N_PAS);
end

function dxdt = ode_dynamics_multipleE_parallel(X, P)
    % This ODE function is now self-contained and safe for parallel execution.
    
    % Unpack parameters from the struct
    k_in   = P.k_in;
    k_e    = P.k_e;
    k_e2   = P.k_e2;
    kHoff  = P.kHoff;
    kc     = P.kc;
    kHon   = P.kHon;
    RE_val_bind_E = P.RE_val_bind_E;
    
    N = P.N;
    PAS = P.PAS;
    
    R   = X(1:N);
    REH = X(N+1:end);

    % Self-consistent Ef_ss calculation
    if P.FirstRun
        if P.Ef_ss_sol == 0
            P.Ef_ss_sol = P.E_total; % Initial guess for free E
        end
        
        REvalbindEAfterPas = RE_val_bind_E(P.Ef_ss_sol);
        E_used = sum(R(1:N)' .* REvalbindEAfterPas) + sum(REH' .* REvalbindEAfterPas(PAS:N));
        Ef_ss = P.E_total - E_used;
        
        % Store the solved Ef_ss back into the struct for use outside fsolve
        P.Ef_ss_sol = Ef_ss;
    else
        Ef_ss = P.Ef_ss_sol; % Use the previously solved value
    end

    Pol_f = P.Pol_total - sum(R) - sum(REH);
    Pol_f = max(0, Pol_f); % Ensure non-negative
    
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
