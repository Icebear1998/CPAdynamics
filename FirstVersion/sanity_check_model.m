global N PAS Pol_total;

%% 1. Define Model Parameters
L_a = 100;  % segment length in bp

% Rate parameters scaled by the timescale
P.k_in    = 2;           % Initiation rate
P.k_c     = 0.6;         % Cleavage rate
kE_on_min = 0.00001;     % Minimum E binding rate
kE_on_max = 0.00025;     % Maximum E binding rate
P.kE_off  = 10;          % E unbinding rate
P.kL_on   = 0.00025;     % L binding rate
P.kL_off  = 0.001;       % L unbinding rate
P.k_e     = 65/L_a;      % Elongation rate before PAS
P.k_e2    = 30/L_a;      % Elongation rate after PAS
P.E_total = 70000;       % Total E molecules
P.L_total = 100000;      % Total L molecules

% Gene and PAS specifications
geneLength_bp = 25000;     % Total gene length in bp
PASposition   = 20000;     % PAS position in bp
N = floor(geneLength_bp / L_a);    % Number of nodes along the gene
PAS = floor(PASposition / L_a);    % Node index corresponding to PAS
Pol_total = 70000;         % Total Pol II number

% Define a position-dependent function for E binding rate (kE_on)
P.kE_on_f = @(n) (n<PAS)*((kE_on_max-kE_on_min)/PAS*n+kE_on_min)+(n>=PAS)*kE_on_max;

%% 2. Set Up the Initial Guess Using ode45
tspan = [0 10000]; % Increased time span for better steady-state approximation
X0 = zeros(3*N-PAS+1,1);
[t, X_ode] = ode45(@(t, x) ode_system(t, x, P), tspan, X0);
X_init = X_ode(end,:);

%% 3. Run fsolve with the Initial Guess
options = optimoptions('fsolve', 'Display', 'iter', 'MaxIterations', 1000, 'MaxFunctionEvaluations', 100000);
[X, fval, exitflag] = fsolve(@(x) ode_system(0, x, P), X_init, options);

% Check convergence
if exitflag > 0
    disp('fsolve converged successfully.');
else
    disp('fsolve did not converge. Check the initial guess or system equations.');
end

% Extract the solutions
R_sol   = X(1:N);
RE_sol  = X(N+1:2*N);
REL_sol = X(2*N+1:end);

%% 4. Plot Main Run Using the Function
plot_spatial_distribution(R_sol, RE_sol, REL_sol, N, PAS, 'Main Run: Final Spatial Distribution: R(l), RE(l), and REL(l)');

% Compute free species for mass balance
E_f_time = P.E_total - sum(RE_sol) - sum(REL_sol);
L_f_time = P.L_total - sum(REL_sol);
Pol_f_time = Pol_total - sum(R_sol) - sum(RE_sol) - sum(REL_sol);

% %% 5. Perform Sanity Checks and Print Key Outputs
% fprintf('\n--- Main Run Checks ---\n');
% R_final   = R_sol(:);
% RE_final  = RE_sol(:);
% REL_final = REL_sol(:);
% 
% % 5.1 Pol II Conservation Check
% Pol_used = sum(R_final) + sum(RE_final) + sum(REL_final);
% Pol_cons = Pol_used + Pol_f_time;
% fprintf('Pol II conservation check:\n');
% fprintf('  Total Pol II (expected): %d\n', Pol_total);
% fprintf('  Sum of R + RE + REL at final time: %g\n', Pol_cons);
% 
% % 5.2 E and L Conservation
% E_bound_final = sum(RE_final) + sum(REL_final);
% E_cons_final  = E_f_time + E_bound_final;
% fprintf('E conservation check at final time:\n');
% fprintf('  E_total (expected): %g\n', P.E_total);
% fprintf('  E_free + E_bound = %g\n', E_cons_final);
% 
% L_bound_final = sum(REL_final);
% L_cons_final  = L_f_time + L_bound_final;
% fprintf('L conservation check at final time:\n');
% fprintf('  L_total (expected): %g\n', P.L_total);
% fprintf('  L_free + L_bound = %g\n', L_cons_final);
% 
% % 5.3 Non-negativity Check at Final Time
% fprintf('\n--- Non-Negativity Check at Final Time ---\n');
% if any(R_final < 0)
%     neg_idx = find(R_final < 0);
%     neg_values = R_final(neg_idx);
%     fprintf('ERROR: Negative R at indices: %s, values: %s\n',...
%          mat2str(neg_idx), mat2str(neg_values));
% end
% if any(RE_final < 0)
%     neg_idx = find(RE_final < 0);
%     neg_values = RE_final(neg_idx);
%     fprintf('ERROR: Negative RE at indices: %s, values: %s\n',...
%          mat2str(neg_idx), mat2str(neg_values));
% end
% if any(REL_final < 0)
%     neg_idx = find(REL_final < 0);
%     neg_values = REL_final(neg_idx);
%     fprintf('ERROR: Negative REL at indices: %s, values: %s\n',...
%          mat2str(neg_idx), mat2str(neg_values));
% end
% if ~any(R_final < 0) && ~any(RE_final < 0) && ~any(REL_final < 0)
%     fprintf('All state variables non-negative.\n');
% end
% 
% %% 6. Edge Case Tests
% fprintf('\n--- Edge Case Tests ---\n');
% 
% % Edge Case 1: kE_on_f = 0 ? RE should be zero
% P_temp = P; P_temp.kE_on_f = @(n) 0;
% [X_test] = fsolve(@(x) ode_system(0, x, P_temp), X_init);
% RE_test = X_test(N+1:2*N);
% if max(max(RE_test)) < 1e-6
%     fprintf('Edge Case 1 (kE_on_f = 0): RE is zero. [PASS]\n');
% else
%     fprintf('Edge Case 1 (kE_on_f = 0): RE is nonzero. [FAIL]\n');
% end
% 
% % Edge Case 2: k_c = 0 ? REL should accumulate and never be cleaved
% P_temp = P; P_temp.k_c = 0;
% [X_test] = fsolve(@(x) ode_system(0, x, P_temp), X_init);
% REL_test = X_test(2*N+1:end);
% if min(diff(REL_test(end, :))) >= -1e-10 
%     fprintf('Edge Case 2 (k_c = 0): REL accumulates as expected. [PASS]\n');
% else
%     fprintf('Edge Case 2 (k_c = 0): REL decreases unexpectedly. [FAIL]\n');
% end
% 
% % Edge Case 3: k_in = 0 ? No Pol II should initiate
% P_temp = P; P_temp.k_in = 0;
% [X_test] = fsolve(@(x) ode_system(0, x, P_temp), X_init);
% if max(max(X_test)) < 1e-6
%     fprintf('Edge Case 3 (k_in = 0): No Pol II initiated as expected. [PASS]\n');
% else
%     fprintf('Edge Case 3 (k_in = 0): Pol II is present. [FAIL]\n');
% end
% 
% % Edge Case 4: k_e = k_e2 = 0 ? Pol II should stay at initiation site
% P_temp = P; P_temp.k_e = 0; P_temp.k_e2 = 0;
% X_init_case4 = zeros(3*N - PAS + 1, 1);
% X_init_case4(1) = Pol_total; % All Pol II at l = 1 in R(1)
% 
% [t, X_test] = ode45(@(t, x) ode_system(t, x, P_temp), tspan, X_init_case4);
% R_test = X_test(end,1:N);
% RE_test = X_test(end,N+1:2*N);
% REL_test = X_test(end,2*N+1:end);
% if all(R_test(2:end) < 1e-6) && all(RE_test(2:end) < 1e-6)
%     fprintf('Edge Case 4 (k_e = 0): Pol II does not elongate. [PASS]\n');
% else
%     fprintf('Edge Case 4 (k_e = 0): Pol II spreads along gene. [FAIL]\n');
%     % Plot the distribution
%     plot_spatial_distribution(R_test, RE_test, REL_test, N, PAS, 'Edge Case 4 (k_e = k_e2 = 0): Pol II Distribution');
% end
% 
% % Edge Case 5: E_total = 0 ? No RE should form
% P_temp = P; P_temp.E_total = 0;
% [X_test] = fsolve(@(x) ode_system(0, x, P_temp), X_init);
% if max(max(X_test(N+1:2*N))) < 1e-6
%     fprintf('Edge Case 5 (E_total = 0): RE remains zero. [PASS]\n');
% else
%     fprintf('Edge Case 5 (E_total = 0): RE is nonzero. [FAIL]\n');
% end
% 
% % Edge Case 6: L_total = 0 ? No REL should form
% P_temp = P; P_temp.L_total = 0;
% [X_test] = fsolve(@(x) ode_system(0, x, P_temp), X_init);
% if max(max(X_test(2*N+1:end))) < 1e-6
%     fprintf('Edge Case 6 (L_total = 0): REL remains zero. [PASS]\n');
% else
%     fprintf('Edge Case 6 (L_total = 0): REL is nonzero. [FAIL]\n');
% end
% 
% % Case 7: k_e2 = 0 => Pol II accumulates at PAS, no exit
% P_temp = P; 
% P_temp.k_e2 = 0;
% X_init_case7 = zeros(3*N - PAS + 1, 1);
% X_init_case7(PAS) = Pol_total; % All Pol II at l = PAS in R(PAS)
% [X_test] = fsolve(@(x) ode_system(0, x, P_temp), X_init_case7);
% R_endA   = X_test(1:N);
% RE_endA  = X_test(N+1:2*N);
% REL_endA = X_test(2*N+1:end);
% accumPAS = R_endA(PAS) + RE_endA(PAS) + REL_endA(PAS - PAS + 1);
% otherSum = (sum(R_endA) + sum(RE_endA) + sum(REL_endA)) - accumPAS;
% if otherSum < 1e-3
%     fprintf('Case 7 (k_e2=0): Pol II accumulates at PAS. [PASS]\n');
% else
%     fprintf('Case 7 (k_e2=0): Pol II also present beyond PAS? [FAIL]\n');
%     fprintf('Pol II at PAS: %g\n', accumPAS);
%     fprintf('Pol II elsewhere: %g\n', otherSum);
%     % Plot the distribution
%     plot_spatial_distribution(R_endA, RE_endA, REL_endA, N, PAS, 'Case 7 (k_e2 = 0): Pol II Distribution');
% end
% 
% % Case 8: k_c = 0 => no cleavage, Pol II accumulates at last position
% P_temp = P;
% P_temp.k_c = 0;
% X_init_case8 = zeros(3*N - PAS + 1, 1);
% X_init_case8(N) = Pol_total; % All Pol II at l = N in R(N)
% [t, X_test] = ode45(@(t, x) ode_system(t, x, P_temp), tspan, X_init_case8);
% R_endB   = X_test(end,1:N);
% RE_endB  = X_test(end,N+1:2*N);
% REL_endB = X_test(end,2*N+1:end);
% lastPosPop = R_endB(N) + RE_endB(N) + REL_endB(N - PAS + 1);
% othersPop  = sum(R_endB) + sum(RE_endB) + sum(REL_endB) - lastPosPop;
% if othersPop < 1e-3
%     fprintf('Case 8 (k_c=0): Pol II accumulates at last node, no cleavage. [PASS]\n');
% else
%     fprintf('Case 8 (k_c=0): Pol II found in earlier nodes? [FAIL]\n');
%     fprintf('Pol II at last node: %g\n', lastPosPop);
%     fprintf('Pol II elsewhere: %g\n', othersPop);
%     % Plot the distribution
%     plot_spatial_distribution(R_endB, RE_endB, REL_endB, N, PAS, 'Case 8 (k_c = 0): Pol II Distribution');
% end
% 
% % Case 9: k_c and kL_on very large => immediate exit
% P_temp = P;
% P_temp.k_c   = 1;
% P_temp.kL_on = 0.1/65; % Assuming timescale = 65 as in previous scripts
% P_temp.k_in  = 0.01; % Reduce initiation rate to achieve near-zero Pol II
% [X_test] = fsolve(@(x) ode_system(0, x, P_temp), X_init);
% R_endC   = X_test(1:N);
% RE_endC  = X_test(N+1:2*N);
% REL_endC = X_test(2*N+1:end);
% finalPolC = sum(R_endC) + sum(RE_endC) + sum(REL_endC);
% if finalPolC < 1e-3
%     fprintf('Case 9 (k_c,kL_on large): Immediate exit achieved. [PASS]\n');
% else
%     fprintf('Case 9 (k_c,kL_on large): Pol II remains in system? [FAIL]\n');
%     fprintf('Total Pol II in system: %g\n', finalPolC);
%     % Plot the distribution
%     plot_spatial_distribution(R_endC, RE_endC, REL_endC, N, PAS, 'Case 9 (k_c, kL_on large): Pol II Distribution');
% end
% 
% fprintf('\n--- All Edge Cases Completed ---\n');
% 
function plot_spatial_distribution(R, RE, REL, N, PAS, plot_title)
    % Plot the spatial distribution of R(l), RE(l), and REL(l)
    % Inputs:
    %   R: Vector of R(l) values (length N)
    %   RE: Vector of RE(l) values (length N)
    %   REL: Vector of REL(l) values (length N - PAS + 1)
    %   N: Total number of nodes along the gene
    %   PAS: Node index of the PAS
    %   plot_title: Title for the plot

    l_values = (1-PAS):(N-PAS); 
    REL_plot = [zeros(PAS-1,1); REL(:)];  % Force REL into a column vector, assume REL=0 for nodes < PAS

    figure;
    hold on;
    plot(l_values, R, 'b-', 'LineWidth', 2.5, 'DisplayName', 'R(l)');
    plot(l_values, RE, 'r-', 'LineWidth', 2.5, 'DisplayName', 'RE(l)');
    plot(l_values, REL_plot, 'g-', 'LineWidth', 3, 'DisplayName', 'REL(l)');
    xlabel('Distance from PAS (bp)', 'FontSize', 14);
    ylabel('Pol II (scaled)', 'FontSize', 14);
    legend('show', 'Location', 'northwest');
    title(plot_title);
    hold off;
end