% PASUsageAnalysis_parallel.m

% MODIFIED: Start a parallel pool of workers. MATLAB will use the default number of cores.
if isempty(gcp('nocreate')); parpool; end

% --- CONFIGURATION ---
save_result = false;
% 1. Define the ranges for the two sweep parameters
kHoff_values = logspace(log10(0.001), log10(0.1), 4);
E_total_values = 60000:40000:220000;

% 2. Define the fixed, biologically relevant distance for the metric
fixed_distance_bp = 300;

% --- BASE PARAMETERS (P struct) ---
% (Parameter definitions are unchanged)
P.L_a = 100; P.k_in = 2; P.kEon = 0.00025; P.kEoff = 10;
P.k_e = 65/100; P.k_e2 = 30/100; P.E_total = 70000;
P.L_total = 100000; P.Pol_total = 70000; P.kHon = 0.25;
P.kHoff = 0.0125; P.kc = 0.05; P.kPon_min = 0.01; P.kPon_max = 0.1;
P.kPoff_min = 0.1; P.kPoff_max = 1; P.kPoff_const = 1;
P.geneLength_bp = 25000; P.PASposition = 20000; P.EBindingNumber = 5;

% --- Pre-allocate the results matrix ---
% MODIFIED: It's good practice to pre-allocate before the parallel loop.
results_matrix = zeros(length(E_total_values), length(kHoff_values));

% --- 2D PARAMETER SWEEP LOOP ---
fprintf('Starting 2D sweep over kHoff and E_total in parallel...\n');

% MODIFIED: Changed the outer 'for' loop to 'parfor'.
% Each worker will process one value of E_total_values (one 'i').
parfor i = 1:length(E_total_values)
    % MODIFIED: Create a temporary row vector for the results of the inner loop.
    % This is a cleaner pattern for parfor.
    temp_row = zeros(1, length(kHoff_values));
    
    % The inner loop remains a standard 'for' loop, executed by the worker.
    for j = 1:length(kHoff_values)
        P_run = P;
        P_run.E_total = E_total_values(i);
        P_run.kHoff = kHoff_values(j);

        % We can't print from inside a parfor loop easily, so it's best to remove this line
        % or use a more advanced method if progress monitoring is essential.
        % fprintf('Running: E_total = %d, kHoff = %.4g\n', E_total_values(i), kHoff_values(j));

        [R_sol, REH_sol, P_sim] = run_termination_simulation(P_run, P_run.EBindingNumber);
        
        if any(isnan(R_sol))
            temp_row(j) = NaN;
            continue;
        end

        flux_cleavage_per_node = P_sim.kc * REH_sol;
        flux_R_exit = P_sim.k_e * R_sol(end);
        flux_REH_exit = P_sim.k_e2 * REH_sol(end);
        total_outflux = sum(flux_cleavage_per_node) + flux_R_exit + flux_REH_exit;

        if total_outflux > 1e-9
            cumulative_exit_flux = cumsum(flux_cleavage_per_node);
            exit_cdf = cumulative_exit_flux / total_outflux;
        else
            exit_cdf = zeros(size(REH_sol));
        end

        nodes_post_pas = 1:length(REH_sol);
        bp_post_pas = nodes_post_pas * P_sim.L_a;
        bp_for_interp = [0; bp_post_pas(:)];
        cdf_for_interp = [0; exit_cdf(:)];
        
        usage_at_dist = interp1(bp_for_interp, cdf_for_interp, fixed_distance_bp, 'linear', 'extrap');
        
        temp_row(j) = usage_at_dist * 100;
    end
    % MODIFIED: Assign the entire computed row to the results matrix.
    results_matrix(i, :) = temp_row;
    fprintf('Finished row for E_total = %d\n', E_total_values(i));
end
disp('All parallel simulations complete.');

% % --- PLOT THE CONTOUR MAP ---
figure('Position', [100, 100, 900, 700]);
[X, Y] = meshgrid(kHoff_values, E_total_values);
contourf(X, Y, results_matrix, 10, 'LineColor', 'k');
set(gca, 'XScale', 'log');
xlabel('PAS Strength (kHoff)', 'FontSize', 12, 'FontWeight', 'bold');
ylabel('Global Factor (E_{total})', 'FontSize', 12, 'FontWeight', 'bold');
title(sprintf('Phase Diagram of Cumulative Polymerase Exit (at %d bp) of Enum = %d', fixed_distance_bp, P.EBindingNumber), 'FontSize', 14, 'FontWeight', 'bold');
c = colorbar;
c.Label.String = 'Cumulative Exit (%)';
c.Label.FontSize = 12;
colormap('parula');
set(gca, 'FontSize', 10);
box on;

% --- Plotting ---
figure('Position', [100, 100, 800, 600]);
hold on;
colors = lines(size(results_matrix, 1));
for i = 1:size(results_matrix, 1)
    plot(kHoff_values, results_matrix(i, :), 'o-', 'LineWidth', 2, ...
         'Color', colors(i, :), 'DisplayName', sprintf('%s = %.2g', 'E_{total}', E_total_values(i)));
end
set(gca, 'XScale', 'log');
xlabel('kHoff', 'FontSize', 12, 'FontWeight', 'bold');
ylabel('Cumulative Exit (%)', 'FontSize', 12, 'FontWeight', 'bold');
title(sprintf('Phase Diagram of Cumulative Polymerase Exit (at %d bp) of Enum = %d', fixed_distance_bp, P.EBindingNumber), 'FontSize', 14, 'FontWeight', 'bold');
grid on; legend('show', 'Location', 'best'); set(gca, 'FontSize', 10); box on;
hold off;

if save_result
    % --- SAVE RESULTS ---
    % Prepare data structure for saving
    data.results_matrix = results_matrix;
    data.x_values = kHoff_values;
    data.y_values = E_total_values;
    data.x_label = 'kHoff';
    data.y_label = 'E_total';
    data.metric_distance = fixed_distance_bp;

    % Save results using the utility function
    save_analysis_results('PASUsageAnalysis', data, P);
end

%% --- Helper function to run the simulation ---
function [R_sol, REH_sol, P] = run_termination_simulation(P, EBindingNumber)
    % MODIFIED: Now returns the P struct as well, as it contains k_e, k_e2 etc.
    global N PAS N_PAS Ef_ss;
    syms Ef real;
    
    L_a = P.L_a; N = floor(P.geneLength_bp / L_a); PAS = floor(P.PASposition / L_a); N_PAS = N - PAS + 1;
    kHon_base = P.kHon;
    
    [r_E_BeforePas] = compute_steady_states(P, EBindingNumber + 1);
    kPon_vals = linspace(P.kPon_min, P.kPon_max, PAS);
    RE_vals = sym(zeros(EBindingNumber + 1, N));
    for e = 1:EBindingNumber + 1
        for idx = 1:length(kPon_vals); RE_vals(e, idx) = subs(r_E_BeforePas(e), {'kPon', 'kPoff'}, {kPon_vals(idx), P.kPoff_const}); end
        for idx = PAS+1:N; RE_vals(e, idx) = subs(r_E_BeforePas(e), {'kPon', 'kPoff'}, {P.kPon_max, P.kPoff_min}); end
    end
    P.RE_val_bind_E = matlabFunction(simplify(sum(sym(1:EBindingNumber)' .* RE_vals(2:end, :), 1)), 'Vars', {Ef});

    X_guess = 1e-6 * ones(N + N_PAS, 1);
    options = optimoptions('fsolve', 'Display', 'off', 'FunctionTolerance', 1e-8);

    P.FirstRun = true; P.is_unphysical = false; Ef_ss = 0;
    try; X_base = fsolve(@(xx) ode_dynamics_multipleE(xx, P), X_guess, options);
    catch; error('Solver failed in Step 1.'); end
    if P.is_unphysical; error('Solver returned unphysical result in Step 1.'); end

    avg_E_bound = P.RE_val_bind_E(Ef_ss);
    P.FirstRun = false; P.kHon = kHon_base * avg_E_bound(end);
    X_final = fsolve(@(xx) ode_dynamics_multipleE(xx, P), X_base, options);
    
    R_sol = X_final(1:N);
    REH_sol = X_final(N+1 : N+N_PAS);
end


%% ------------ ODE DYNAMICS FUNCTION ------------
function dxdt = ode_dynamics_multipleE(X, P)
global N PAS Ef_ss

k_in   = P.k_in;
k_e    = P.k_e;
k_e2   = P.k_e2;
kHoff_t= P.kHoff;
kc_t   = P.kc;
RE_val_bind_E = P.RE_val_bind_E;
kHon_t = P.kHon;

R   = X(1:N);
REH = X(N+1:end);

if P.FirstRun
    % Set initial Ef_ss once (using P.E_total as a starting guess)
    if Ef_ss == 0
        Ef_ss = P.E_total; % Initial guess
    end

    % Convert symbolic expression to a numerical function
    REvalbindEAfterPas = RE_val_bind_E(Ef_ss);
    E_used = sum(R(1:N)'.* RE_val_bind_E(Ef_ss)) + sum(REH' .* REvalbindEAfterPas(PAS:N)); %sum(REH, 1); 
    E_f = abs(P.E_total - E_used);
    %E_f = max(1, P.E_total - E_used);
    Ef_ss = E_f;
    %disp({sum(RE_val_bind_E(Ef_ss)), sum(R(1:PAS)), Ef_ss});

    % If E_f < 0, throw error and stop solver
    if Ef_ss < 0
        P.is_unphysical = true; % Set flag instead of throwing error
        dxdt = 1e6 * ones(length(X), 1); % Large residual to signal fsolve to stop
        return;
    end
end

Pol_f = P.Pol_total - sum(R) - sum(REH);
% L_f = P.L_total - sum()

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