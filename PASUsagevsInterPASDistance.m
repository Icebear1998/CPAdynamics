% SCRIPT to calculate and plot Proximal PAS Usage vs. Inter-PAS Distance

clear all;
global N PAS N_PAS Ef_ss;
syms Ef real;

% --- Model parameters ---
L_a = 100; % bp per node
P.k_in = 2;
P.kEon = 0.00025;
P.kEoff = 10;
P.k_e = 65/L_a;
P.k_e2 = 30/L_a;
P.E_total = 70000;
P.Pol_total = 70000;
P.kHon = 0.05;
P.kHoff = 0.0025;
P.kc = 0.8;
P.kPon_min = 0.01;
P.kPon_max = 1;
P.kPoff_const = 1;
P.kPoff_max = 2;
P.kPoff_min = 0.1;
EBindingNumber = 2; % Set the desired EBindingNumber

geneLength_bp = 25000;
PASposition = 20000;
N = floor(geneLength_bp / L_a);
PAS = floor(PASposition / L_a);
N_PAS = N - PAS + 1;

% --- Run Simulation to Get Termination Profile ---
% (This section uses the robust iterative solver to get a stable result)
disp('Starting simulation...');
kHon_base = P.kHon; % Store the base kHon for this run
    
% --- Symbolic Pre-computation ---
[r_E_BeforePas] = compute_steady_states(P, EBindingNumber + 1);
kPon_vals = linspace(P.kPon_min, P.kPon_max, PAS);
RE_vals = sym(zeros(EBindingNumber + 1, N));
for e = 1:EBindingNumber + 1
    for idx = 1:length(kPon_vals)
        kPon_val = kPon_vals(idx);
        RE_vals(e, idx) = subs(r_E_BeforePas(e), {'kPon', 'kPoff'}, {kPon_val, P.kPoff_const});
    end
    for idx = PAS+1:N
        RE_vals(e, idx) = subs(r_E_BeforePas(e), {'kPon', 'kPoff'}, {P.kPon_max, P.kPoff_min});
    end
end
P.RE_val_bind_E = matlabFunction(simplify(sum(sym(1:EBindingNumber)' .* RE_vals(2:end, :), 1)), 'Vars', {Ef});

% --- STABLE TWO-STEP SOLVER ---
X_guess = 1e-6 * ones(N + N_PAS, 1);
options = optimoptions('fsolve', 'Display', 'off', 'FunctionTolerance', 1e-8);

% Step 1: Solve for the initial steady state and self-consistent Ef_ss
P.FirstRun = true; 
P.is_unphysical = false; 
Ef_ss = 0;
try
    X_base = fsolve(@(xx) ode_dynamics_multipleE(xx, P), X_guess, options);
catch
    error('Solver failed in Step 1.');
end
if P.is_unphysical; error('Solver returned unphysical result in Step 1.'); end

% Step 2: Update kHon and solve for the final steady state
avg_E_bound = P.RE_val_bind_E(Ef_ss);
P.FirstRun = false; 
P.kHon = kHon_base * avg_E_bound(end);

% Use the result of Step 1 as the initial guess for Step 2
X_final = fsolve(@(xx) ode_dynamics_multipleE(xx, P), X_base, options);

% --- Calculate and return the final ratio ---
R_sol = X_final(1:N);
REH_sol = X_final(N+1 : N+N_PAS);
ratio_vector = (REH_sol(1:end) + R_sol(PAS:end)) / (R_sol(PAS-1) + 1e-9);

% --- GENERATE THE REQUESTED PLOT ---
disp('Generating APA plot...');

% 1. Define the relevant distances between a proximal and a distal PAS
inter_pas_distances_bp = 0:100:2500; % e.g., from 0 to 2.5 kb

% 2. Calculate the proximal site usage for each distance
proximal_usage_prob = zeros(size(inter_pas_distances_bp));
for i = 1:length(inter_pas_distances_bp)
    dist_bp = inter_pas_distances_bp(i);
    % Convert distance in bp to a node index
    node_idx = round(dist_bp / L_a);

    % Ensure node_idx is within the bounds of our 'ratio' vector
    if node_idx < 1
        node_idx = 1;
    end
    if node_idx > length(ratio_vector)
        node_idx = length(ratio_vector);
    end

    % The probability of using the DISTAL site is the read-through 'ratio'
    distal_prob = ratio_vector(node_idx);
    
    % The choice is binary, so the probability of using the PROXIMAL site is 1 minus the distal probability
    proximal_usage_prob(i) = 1 - distal_prob;
end

% 3. Create the plot
figure('Position', [100, 100, 800, 600]);
hold on;

plot(inter_pas_distances_bp, proximal_usage_prob * 100, 'o-', 'LineWidth', 2, 'MarkerSize', 6);

% Add a vertical line at 300 nt for biological context
line([300 300], ylim, 'Color', 'r', 'LineStyle', '--', 'LineWidth', 1.5, 'DisplayName', 'Median Tandem Distance');

% Customize the plot
xlabel('Inter-PAS Distance (bp)', 'FontSize', 12);
ylabel('Proximal Site Usage (%)', 'FontSize', 12);
title(['Predicted Proximal PAS Usage (EBindingNumber=', num2str(EBindingNumber), ')'], 'FontSize', 14, 'FontWeight', 'bold');
grid on;
legend('show', 'Location', 'best');
set(gca, 'FontSize', 10);
box on;
ylim([0 100]); % Ensure y-axis is from 0% to 100%

disp('Done.');