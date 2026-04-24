global N PAS N_PAS Ef_ss;
syms Ef real;
saveData = false;
% ------------ MODEL PARAMETERS ------------
P.L_a = 100;
P.k_in    = 2;
P.k_e     = 65/P.L_a;
P.k_e2    = 30/P.L_a;
P.E_total = 100000;
P.L_total = 100000;
P.Pol_total = 70000;

% We are not confident here
P.kEon = 0.0000025;
P.kEoff = 0.1;
P.kHon = 2; 
P.kHoff = 1; % for early poly A site with kd ~ 2000
P.kc = 0.1; 

P.kPon_min = 0.01; % at TSS
P.kPon_slope = 0.005; % determined how fast Sep2P increasing from TSS
P.kPoff = 1;

P.geneLength_bp = 25000;
P.PASposition   = 20000;

N      = floor(P.geneLength_bp / P.L_a);  % total nodes
PAS    = floor(P.PASposition   / P.L_a);  % node index of PAS
N_PAS  = N - PAS + 1;                 % number of nodes at/after PAS
Ef_ss = 0;

EBindingNumber = 1; 

% Run termination simulation using the new function
[R_sol, REH_sol, P, r_E_BeforePas, r_P] = run_termination_simulation(P, EBindingNumber);
disp('done compute simulation');

% Get average E bound and Ser2P profiles in a single pass
fprintf('E free %d\n', Ef_ss);
kPon_vals = P.kPon_min + P.kPon_slope * (0:N-1);
[avg_E_bound, Ser2P] = P.RE_val_bind_E(Ef_ss);

[exit_cdf, distances_bp] = calculate_pas_usage_profile(R_sol, REH_sol, P);
CAD = interp1(exit_cdf, distances_bp, 0.5, 'linear', 'extrap');

hold on;
plot((1-PAS):N_PAS-1, Ser2P, 'g-','LineWidth',2.5, 'DisplayName','Ser2P');
plot((1-PAS):N_PAS-1, avg_E_bound, 'b-','LineWidth',2.5, 'DisplayName','AverageE');
legend({'Ser2P', 'AverageE'}, 'Location', 'best');
xlabel('position'); ylabel('AverageE');
hold off;
% ------------ PLOT RESULTS ------------
% 1. Time evolution plot
l_values =  (1-PAS):(N-PAS);

figure; hold on;
% for e = 1:EBindingNumber+1
%     plot((1-PAS):N_PAS-1, RE_vals(e,:), 'LineWidth',2);
% end

plot(l_values, R_sol, 'b-','LineWidth',2.5, 'DisplayName','R(l)');
plot(l_values, [zeros(PAS-1,1);REH_sol], 'r-','LineWidth',2.5, 'DisplayName','REH(l)');
CAD = CAD/100;
line([CAD CAD], ylim, 'Color', 'r', 'LineStyle', '--', 'LineWidth', 1.5, 'DisplayName', 'Median Tandem Distance');
text(CAD, 8, num2str(int32(CAD)));
xlabel('Position relative to PAS (x100 basepair)'); ylabel('Concentration');
legend({'Total R', 'Total REH', '50% TCD'}, 'Location', 'best');
title('CPA model dynamic simulation');
hold off;

% 1. Calculate the final, steady-state concentration of free Pol II
% This is the total Pol II minus all polymerases bound to the gene (R and REH).
bound_pol_II = sum(R_sol) + sum(REH_sol);
Pol_f_final = P.Pol_total - bound_pol_II;

% 2. Display the values in the command window for clarity
fprintf('\n--- Polymerase Distribution at Steady State ---\n');
fprintf('Total Pol II in system: %d\n', P.Pol_total);
fprintf('Total Bound Pol II (on gene): %.2f\n', bound_pol_II);
fprintf('Total Free Pol II (Pol_f):    %.2f\n', Pol_f_final);

if saveData
    % --- SAVE RESULTS ---
    % Prepare data structure for saving
    data.EBindingNumber = EBindingNumber;
    data.R_sol = R_sol;
    data.REH_sol = REH_sol;
    data.Ser2P = Ser2P;
    data.avg_E_bound = avg_E_bound;
    data.Ef_ss = Ef_ss;
    data.Pol_f_final = Pol_f_final;

    % Save results using the utility function
    save_analysis_results('CPA_multipleE_main', data, P);
end


