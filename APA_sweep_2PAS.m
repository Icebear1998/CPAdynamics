% SCRIPT for 2D parameter sweep of APA site choice
% Varies inter-PAS distance and relative PAS strength

% clear all;
% close all;
% clc;

% --- BASE PARAMETERS ---
% These are the fixed parameters for all simulations
P_base.k_in = 2;
P_base.k_e = 65/100;
P_base.k_e2 = 30/100;
P_base.E_total = 70000;
P_base.Pol_total = 70000;
P_base.kHoff_prox = 0.0025;
P_base.kc_prox = 0.8;
P_base.kHoff_dist = 0.0025;
P_base.kc_dist = 0.8;
P_base.L_a = 100;
P_base.geneLength_bp = 25000;
P_base.PAS_prox_pos = 20000;

% For this simulation, we use a simplified constant for average E-binding
% A full model would calculate this dynamically
P_base.avg_E_bound_value = 1.0;

% --- SWEEP PARAMETERS ---
% Define the ranges for the two axes of the parameter sweep

% Axis 1: Distance between proximal and distal PAS (bp)
inter_pas_distances = 100:100:1500;

% Axis 2: Relative strength of the proximal PAS
% This is the ratio: kHon_prox / kHon_distal
% A value of 0.2 means the proximal site is 20% as strong as the distal one.
strength_ratios = 0.1:0.1:1.0;

% Set a fixed, strong value for the distal PAS. We will vary the proximal one.
P_base.kHon_dist = 0.1;

% --- INITIALIZE RESULTS MATRIX ---
% This matrix will store the proximal usage percentage for each parameter pair
proximal_usage_matrix = zeros(length(strength_ratios), length(inter_pas_distances));

% --- 2D PARAMETER SWEEP LOOP ---
disp('Starting 2D parameter sweep for APA...');
total_sims = length(strength_ratios) * length(inter_pas_distances);
sim_count = 0;

% Outer loop: Iterate over relative strengths
for i = 1:length(strength_ratios)
    % Inner loop: Iterate over distances
    for j = 1:length(inter_pas_distances)
        sim_count = sim_count + 1;
        fprintf('Running simulation %d of %d...\n', sim_count, total_sims);
        
        % 1. Create a fresh parameter set for this specific simulation
        P = P_base;
        
        % 2. Set the parameters for this specific run
        P.kHon_prox = P.kHon_dist * strength_ratios(i); % Set proximal strength
        PAS_dist_pos = P.PAS_prox_pos + inter_pas_distances(j); % Set distal position
        
        % 3. Calculate node indices and vector lengths for this run
        P.N = floor(P.geneLength_bp / P.L_a);
        P.PAS_prox = floor(P.PAS_prox_pos / P.L_a);
        P.PAS_dist = floor(PAS_dist_pos / P.L_a);
        P.N_PAS_prox = P.N - P.PAS_prox + 1;
        P.N_PAS_dist = P.N - P.PAS_dist + 1;
        
        % Ensure distal PAS is within the gene
        if P.PAS_dist >= P.N
            proximal_usage_matrix(i, j) = NaN; % Mark as invalid if distal PAS is outside
            continue;
        end
        
        % 4. Set up the initial guess vector with the correct size
        total_vars = P.N + P.N_PAS_prox + P.N_PAS_dist;
        X0 = 1e-6 * ones(total_vars, 1);
        
        % 5. Run the solver
        options = optimoptions('fsolve', 'Display', 'off', 'FunctionTolerance', 1e-9);
        [X_sol, ~, exitflag] = fsolve(@(xx) ode_dynamics_2PAS_sweep(xx, P), X0, options);
        
        % 6. Process the result
        if exitflag > 0 % Check if solver converged
            REH_prox_sol = X_sol(P.N+1 : P.N+P.N_PAS_prox);
            REH_dist_sol = X_sol(P.N+P.N_PAS_prox+1 : end);
            
            total_prox_term = sum(REH_prox_sol);
            total_dist_term = sum(REH_dist_sol);
            total_term = total_prox_term + total_dist_term;
            
            if total_term > 0
                proximal_usage_matrix(i, j) = (total_prox_term / total_term) * 100;
            else
                proximal_usage_matrix(i, j) = 0;
            end
        else
            proximal_usage_matrix(i, j) = NaN; % Mark as failed convergence
        end
    end
end
disp('Parameter sweep complete.');

% --- PLOT THE RESULTS ---
disp('Generating plot...');
figure('Position', [100, 100, 900, 700]);
hold on;

% Use imagesc for a clean heatmap
imagesc(inter_pas_distances, strength_ratios, proximal_usage_matrix);
set(gca, 'YDir', 'normal'); % Place 0.1 ratio at the bottom

% Add contour lines for clarity
contour(inter_pas_distances, strength_ratios, proximal_usage_matrix, 'ShowText', 'on', 'LineColor', 'k');

% Customize the plot
xlabel('Inter-PAS Distance (bp)', 'FontSize', 12, 'FontWeight', 'bold');
ylabel('Relative Proximal PAS Strength (kHon_{prox} / kHon_{dist})', 'FontSize', 12, 'FontWeight', 'bold');
title('APA Site Choice: Distance vs. Relative PAS Strength', 'FontSize', 14, 'FontWeight', 'bold');
c = colorbar;
c.Label.String = 'Proximal Site Usage (%)';
c.Label.FontSize = 12;
c.Label.FontWeight = 'bold';
colormap('parula'); % A nice, perceptually uniform colormap
grid on;
set(gca, 'FontSize', 10);
box on;

%% --- ODE DYNAMICS FUNCTION (MODIFIED FOR SWEEP) ---
% This version takes all node info from the 'P' struct
function dxdt = ode_dynamics_2PAS_sweep(X, P)
    % Unpack node indices from P struct
    N = P.N;
    PAS_prox = P.PAS_prox;
    PAS_dist = P.PAS_dist;
    N_PAS_prox = P.N_PAS_prox;
    
    % Unpack parameters from P
    k_in = P.k_in;
    k_e = P.k_e; k_e2 = P.k_e2;
    kHon_prox = P.kHon_prox; kHoff_prox = P.kHoff_prox; kc_prox = P.kc_prox;
    kHon_dist = P.kHon_dist; kHoff_dist = P.kHoff_dist; kc_dist = P.kc_dist;

    % Unpack state vector X
    R = X(1:N);
    REH_prox = X(N+1 : N+N_PAS_prox);
    REH_dist = X(N+N_PAS_prox+1 : end);

    % Simplified Ef_ss calculation
    E_used = (sum(R) + sum(REH_prox) + sum(REH_dist)) * P.avg_E_bound_value;
    Ef_ss = max(0, P.E_total - E_used);

    Pol_f = P.Pol_total - sum(R) - sum(REH_prox) - sum(REH_dist);
    dxdt = zeros(length(X), 1);

    % Region 1: Before any PAS
    dxdt(1) = Pol_f*k_in - k_e*R(1);
    for n = 2:(PAS_prox-1)
        dxdt(n) = k_e*R(n-1) - k_e*R(n);
    end

    % Region 2: Between Proximal and Distal PAS
    for n = PAS_prox:(PAS_dist-1)
        j = n - PAS_prox + 1;
        dxdt_R = k_e*R(n-1) - k_e*R(n) - kHon_prox*R(n) + kHoff_prox*REH_prox(j);
        if j == 1
            dxdt_REH_prox = -k_e2*REH_prox(j) + kHon_prox*R(n) - kHoff_prox*REH_prox(j) - kc_prox*REH_prox(j);
        else
            dxdt_REH_prox = k_e2*REH_prox(j-1) - k_e2*REH_prox(j) + kHon_prox*R(n) - kHoff_prox*REH_prox(j) - kc_prox*REH_prox(j);
        end
        dxdt(n) = dxdt_R;
        dxdt(N+j) = dxdt_REH_prox;
    end

    % Region 3: At and After Distal PAS
    for n = PAS_dist:N
        j = n - PAS_prox + 1;
        k = n - PAS_dist + 1;
        dxdt_R = k_e*R(n-1) - k_e*R(n) - kHon_prox*R(n) + kHoff_prox*REH_prox(j) - kHon_dist*R(n) + kHoff_dist*REH_dist(k);
        dxdt_REH_prox = k_e2*REH_prox(j-1) - k_e2*REH_prox(j) + kHon_prox*R(n) - kHoff_prox*REH_prox(j) - kc_prox*REH_prox(j);
        if k == 1
            dxdt_REH_dist = -k_e2*REH_dist(k) + kHon_dist*R(n) - kHoff_dist*REH_dist(k) - kc_dist*REH_dist(k);
        else
            dxdt_REH_dist = k_e2*REH_dist(k-1) - k_e2*REH_dist(k) + kHon_dist*R(n) - kHoff_dist*REH_dist(k) - kc_dist*REH_dist(k);
        end
        dxdt(n) = dxdt_R;
        dxdt(N+j) = dxdt_REH_prox;
        dxdt(N+N_PAS_prox+k) = dxdt_REH_dist;
    end
end