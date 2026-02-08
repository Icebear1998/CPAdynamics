%% DIAGNOSTIC SCRIPT: E_free Calculation Issue
% This script tracks and visualizes the E_free calculation problem during fsolve iterations
% Issue: E_used can exceed E_total, especially for EBindingNumber > 1

% clear all;
% close all;

%% Setup logging mechanism
global N PAS N_PAS Ef_ss diagnostic_log;

% Initialize diagnostic log
diagnostic_log = struct();
diagnostic_log.iteration = [];
diagnostic_log.E_total = [];
diagnostic_log.E_used = [];
diagnostic_log.E_free = [];
diagnostic_log.E_free_negative = [];
diagnostic_log.max_R = [];
diagnostic_log.max_REH = [];
diagnostic_log.avg_E_bound_per_R = [];
diagnostic_log.avg_E_bound_per_REH = [];

%% Model parameters (base values)
L_a = 100;
P.k_in = 2;
P.kEon = 0.00025;
P.kEoff = 10;
P.k_e = 65/L_a;
P.k_e2 = 30/L_a;
P.E_total = 70000;
P.L_total = 100000;
P.Pol_total = 70000;
P.kHon = 0.2;
P.kHoff = 0.0125;
P.kc = 0.05;
P.kPon_min = 0.01;
P.kPon_max = 1.5;
P.kPoff_min = 0.1;
P.kPoff_max = 20;
P.kPoff_const = 1;
P.kPon_const = 1;

geneLength_bp = 25000;
PASposition = 20000;
N = floor(geneLength_bp / L_a);
PAS = floor(PASposition / L_a);
N_PAS = N - PAS + 1;

%% Test multiple EBindingNumber values
EBindingNumbers_to_test = [1, 2, 3];
results_summary = struct();

for idx = 1:length(EBindingNumbers_to_test)
    EBindingNumber = EBindingNumbers_to_test(idx);
    
    fprintf('\n========================================\n');
    fprintf('Testing EBindingNumber = %d\n', EBindingNumber);
    fprintf('========================================\n');
    
    % Reset diagnostic log
    diagnostic_log = struct();
    diagnostic_log.iteration = [];
    diagnostic_log.E_total = [];
    diagnostic_log.E_used = [];
    diagnostic_log.E_free = [];
    diagnostic_log.E_free_negative = [];
    diagnostic_log.max_R = [];
    diagnostic_log.max_REH = [];
    diagnostic_log.avg_E_bound_per_R = [];
    diagnostic_log.avg_E_bound_per_REH = [];
    
    Ef_ss = 0;
    
    try
        % Compute steady states for E binding
        syms Ef real;
        [r_E_BeforePas] = compute_steady_states(P, EBindingNumber + 1);
        
        % Setup position-dependent kPon
        kPon_vals = linspace(P.kPon_min, P.kPon_max, PAS);
        RE_vals = sym(zeros(EBindingNumber+1, N));
        
        for e = 1:EBindingNumber+1
            for i = 1:PAS
                kPon_val = kPon_vals(i);
                RE_vals(e, i) = subs(r_E_BeforePas(e), {'kPon', 'kPoff'}, {kPon_val, P.kPoff_const});
            end
            for i = PAS+1:N
                RE_vals(e, i) = subs(r_E_BeforePas(e), {'kPon', 'kPoff'}, {P.kPon_max, P.kPoff_const});
            end
        end
        
        % Create average E binding function
        P.RE_val_bind_E = matlabFunction(simplify(sum(sym(1:EBindingNumber)' .* RE_vals(2:end, :), 1)), 'Vars', {Ef});
        
        % Create probability distribution function for E binding states
        P.RE_vals_func = cell(EBindingNumber+1, 1);
        for e = 1:EBindingNumber+1
            P.RE_vals_func{e} = matlabFunction(RE_vals(e, :), 'Vars', {Ef});
        end
        
        % Initial condition
        P.FirstRun = true;
        P.is_unphysical = false;
        P.diagnostic_mode = true;  % Enable diagnostic tracking
        X0 = 1e-6 * ones(N + N_PAS, 1);
        
        % Solve with tracking
        options = optimoptions('fsolve', 'Display', 'iter', ...
            'FunctionTolerance', 1e-8, 'MaxIterations', 1000, ...
            'OutputFcn', @(x, optimValues, state) trackIteration(x, optimValues, state, P));
        
        fprintf('\nSolving ODE system...\n');
        X = fsolve(@(xx) ode_dynamics_multipleE_diagnostic(xx, P), X0, options);
        
        % Store results
        results_summary(idx).EBindingNumber = EBindingNumber;
        results_summary(idx).converged = ~P.is_unphysical && ~any(isnan(X));
        results_summary(idx).final_Ef_ss = Ef_ss;
        results_summary(idx).diagnostic_log = diagnostic_log;
        
        % Print summary statistics
        fprintf('\n--- Summary Statistics ---\n');
        fprintf('Converged: %s\n', mat2str(results_summary(idx).converged));
        fprintf('Final E_free: %.2f\n', Ef_ss);
        fprintf('E_total: %.2f\n', P.E_total);
        fprintf('Number of iterations: %d\n', length(diagnostic_log.iteration));
        fprintf('Times E_free went negative: %d\n', sum(diagnostic_log.E_free_negative));
        fprintf('Max E_used: %.2f (%.1f%% of E_total)\n', ...
            max(diagnostic_log.E_used), max(diagnostic_log.E_used)/P.E_total*100);
        fprintf('Min E_free: %.2f\n', min(diagnostic_log.E_free));
        
    catch ME
        fprintf('\nERROR: %s\n', ME.message);
        results_summary(idx).EBindingNumber = EBindingNumber;
        results_summary(idx).converged = false;
        results_summary(idx).error_message = ME.message;
        results_summary(idx).diagnostic_log = diagnostic_log;
    end
end

%% Visualization
fprintf('\n========================================\n');
fprintf('Generating diagnostic plots...\n');
fprintf('========================================\n');

for idx = 1:length(results_summary)
    if ~results_summary(idx).converged
        continue;
    end
    
    log = results_summary(idx).diagnostic_log;
    EBindingNumber = results_summary(idx).EBindingNumber;
    
    figure('Position', [100, 100, 1400, 900]);
    sgtitle(sprintf('E_{free} Calculation Diagnostics (EBindingNumber = %d)', EBindingNumber), ...
        'FontSize', 14, 'FontWeight', 'bold');
    
    % Plot 1: E_total vs E_used
    subplot(2, 3, 1);
    plot(log.iteration, log.E_total, 'k-', 'LineWidth', 2, 'DisplayName', 'E_{total}');
    hold on;
    plot(log.iteration, log.E_used, 'b-', 'LineWidth', 1.5, 'DisplayName', 'E_{used}');
    plot(log.iteration, log.E_free, 'r-', 'LineWidth', 1.5, 'DisplayName', 'E_{free}');
    % Highlight when E_used > E_total
    violation_idx = log.E_used > log.E_total;
    if any(violation_idx)
        scatter(log.iteration(violation_idx), log.E_used(violation_idx), 100, 'rx', ...
            'LineWidth', 3, 'DisplayName', 'E_{used} > E_{total}');
    end
    xlabel('Iteration');
    ylabel('Concentration');
    title('E Conservation Tracking');
    legend('Location', 'best');
    grid on;
    
    % Plot 2: E_free over iterations
    subplot(2, 3, 2);
    plot(log.iteration, log.E_free, 'r-', 'LineWidth', 2);
    hold on;
    yline(0, 'k--', 'LineWidth', 1.5, 'DisplayName', 'Zero line');
    xlabel('Iteration');
    ylabel('E_{free}');
    title('Free E Factor Evolution');
    grid on;
    ylim([min(log.E_free)*1.2, max(log.E_free)*1.2]);
    
    % Plot 3: E_used / E_total ratio
    subplot(2, 3, 3);
    ratio = log.E_used ./ log.E_total;
    plot(log.iteration, ratio, 'b-', 'LineWidth', 2);
    hold on;
    yline(1, 'r--', 'LineWidth', 2, 'DisplayName', 'Ratio = 1');
    xlabel('Iteration');
    ylabel('E_{used} / E_{total}');
    title('E Usage Ratio');
    legend('Location', 'best');
    grid on;
    
    % Plot 4: Max polymerase concentrations
    subplot(2, 3, 4);
    yyaxis left;
    plot(log.iteration, log.max_R, 'b-', 'LineWidth', 1.5);
    ylabel('max(R)');
    yyaxis right;
    plot(log.iteration, log.max_REH, 'r-', 'LineWidth', 1.5);
    ylabel('max(REH)');
    xlabel('Iteration');
    title('Maximum Polymerase Concentrations');
    grid on;
    
    % Plot 5: Average E bound per polymerase
    subplot(2, 3, 5);
    plot(log.iteration, log.avg_E_bound_per_R, 'b-', 'LineWidth', 2, 'DisplayName', 'per R');
    hold on;
    plot(log.iteration, log.avg_E_bound_per_REH, 'r-', 'LineWidth', 2, 'DisplayName', 'per REH');
    xlabel('Iteration');
    ylabel('Avg E bound');
    title('Average E Factors Bound per Polymerase');
    legend('Location', 'best');
    grid on;
    
    % Plot 6: Phase diagram (E_free vs E_used)
    subplot(2, 3, 6);
    scatter(log.E_free, log.E_used, 50, log.iteration, 'filled');
    hold on;
    plot([0 max(log.E_free)], [P.E_total P.E_total-max(log.E_free)], 'r--', ...
        'LineWidth', 2, 'DisplayName', 'E_{free} + E_{used} = E_{total}');
    xlabel('E_{free}');
    ylabel('E_{used}');
    title('Phase Diagram');
    colorbar;
    ylabel(colorbar, 'Iteration');
    legend('Location', 'best');
    grid on;
    axis equal;
    xlim([min(log.E_free)*1.1, max(log.E_free)*1.1]);
    ylim([min(log.E_used)*0.9, max(log.E_used)*1.1]);
end

%% Analysis and Recommendations
fprintf('\n========================================\n');
fprintf('DIAGNOSTIC ANALYSIS\n');
fprintf('========================================\n');

for idx = 1:length(results_summary)
    fprintf('\n--- EBindingNumber = %d ---\n', results_summary(idx).EBindingNumber);
    
    if ~results_summary(idx).converged
        fprintf('FAILED TO CONVERGE\n');
        if isfield(results_summary(idx), 'error_message')
            fprintf('Error: %s\n', results_summary(idx).error_message);
        end
        continue;
    end
    
    log = results_summary(idx).diagnostic_log;
    
    % Check for issues
    n_violations = sum(log.E_used > log.E_total);
    max_violation = max(log.E_used - log.E_total);
    
    fprintf('Conservation violations: %d times\n', n_violations);
    fprintf('Max violation: %.2f (%.1f%% of E_total)\n', ...
        max_violation, max_violation/P.E_total*100);
    
    if n_violations > 0
        fprintf('⚠�?  ISSUE CONFIRMED: E_used exceeded E_total during iterations\n');
    else
        fprintf('✓ No conservation violations detected\n');
    end
end

fprintf('\n========================================\n');
fprintf('KEY INSIGHTS\n');
fprintf('========================================\n');
fprintf('1. The issue is a circular dependency:\n');
fprintf('   - E_free is calculated as: E_free = E_total - E_used\n');
fprintf('   - But E_used depends on E_free: E_used = f(R, REH, E_free)\n');
fprintf('   - During fsolve iterations, R and REH are not self-consistent with E_free\n\n');

fprintf('2. Why EBindingNumber > 1 is worse:\n');
fprintf('   - Higher E binding means more E factors bound per polymerase\n');
fprintf('   - RE_val_bind_E(Ef) becomes more sensitive to Ef\n');
fprintf('   - Small errors in Ef estimate lead to large errors in E_used\n\n');

fprintf('3. The abs() "fix" is incorrect because:\n');
fprintf('   - It hides the conservation violation\n');
fprintf('   - It allows unphysical states (negative E_free becomes positive)\n');
fprintf('   - The solver converges to wrong solution\n\n');

fprintf('========================================\n');
fprintf('INFORMATION NEEDED TO FIX\n');
fprintf('========================================\n');
fprintf('To properly fix this issue, I need:\n\n');
fprintf('1. ✓ RE_val_bind_E function structure - HAVE THIS\n');
fprintf('2. ✓ How E factors bind to polymerases - HAVE THIS\n');
fprintf('3. ⚠�?  NEED: Should we solve for E_free self-consistently WITHIN each ODE evaluation?\n');
fprintf('4. ⚠�?  NEED: Or should we iterate E_free OUTSIDE the fsolve loop?\n');
fprintf('5. ⚠�?  NEED: Typical convergence criteria for E_free iteration\n');
fprintf('6. ✓ Conservation laws - HAVE THIS\n\n');

fprintf('PROPOSED SOLUTIONS:\n');
fprintf('A) Self-consistent solver (RECOMMENDED):\n');
fprintf('   - Within each ode_dynamics call, solve for E_free that satisfies:\n');
fprintf('     E_free = E_total - sum(R .* RE_val_bind_E(E_free)) - sum(REH .* RE_val_bind_E(E_free)[PAS:N])\n');
fprintf('   - Use iterative method (fixed point or Newton) to find E_free\n\n');

fprintf('B) Outer iteration loop:\n');
fprintf('   - Iterate: solve ODE -> update E_free -> solve ODE -> ...\n');
fprintf('   - Until E_free converges\n\n');

fprintf('C) Coupled system:\n');
fprintf('   - Add E_free as an additional variable in the ODE system\n');
fprintf('   - Add conservation equation as additional constraint\n\n');

fprintf('Please review the diagnostic plots and confirm which approach to use.\n');
fprintf('========================================\n');

%% Helper function for tracking iterations
function stop = trackIteration(x, optimValues, state, P)
    stop = false;
    
    if strcmp(state, 'iter')
        % For fsolve, the available fields are: iteration, funccount, fval
        fprintf('  fsolve iteration %d, function count = %d, |F(x)| = %.4e\n', ...
            optimValues.iteration, optimValues.funccount, norm(optimValues.fval));
    end
end

