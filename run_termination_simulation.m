function [R_sol, REH_sol, P, r_E_BeforePas, r_P] = run_termination_simulation(P, EBindingNumber)
    % RUN_TERMINATION_SIMULATION Run single simulation for given E binding number
    %
    % Inputs:
    %   P - Parameter structure containing all simulation parameters
    %   EBindingNumber - Number of E factors that can bind to polymerase
    %
    % Outputs:
    %   R_sol - Solution vector for active polymerase concentrations
    %   REH_sol - Solution vector for terminating polymerase concentrations
    %   P - Updated parameter structure with computed values
    %   r_E_BeforePas - Unused (kept for call-site compatibility); always []
    %   r_P - Unused (kept for call-site compatibility); always []
    
    global N PAS N_PAS Ef_ss;
    
    % Setup geometry
    L_a = P.L_a;
    N = floor(P.geneLength_bp / L_a);
    PAS = floor(P.PASposition / L_a);
    N_PAS = N - PAS + 1;
    kHon_base = P.kHon;
    
    % Symbolic outputs are no longer computed; numerical approach is used throughout.
    r_E_BeforePas = [];
    r_P = [];
    
    % Setup kPon values with linear increase
    kPon_vals = P.kPon_min + P.kPon_slope * (0:N-1);
    
    % Create E binding function using numerical computation
    % For high EBindingNumber (>=5), symbolic expressions become too complex
    % and cause numerical overflow/NaN. Use numerical null-space instead.
    n_states = EBindingNumber + 1;
    
    % OPTIMIZATION: Pre-compute the SVDs over a grid of Ef_val and interpolate
    % Doing SVDs inside the ODE fsolve loop is millions of times too slow
    % because of the nested solvers.
    num_grid_points = 100;
    Ef_grid = linspace(0, P.E_total, num_grid_points);
    avg_E_bound_grid = zeros(num_grid_points, length(kPon_vals));
    avg_Ser2P_grid = zeros(num_grid_points, length(kPon_vals));
    
    for i = 1:num_grid_points
        [avg_E_bound_grid(i, :), avg_Ser2P_grid(i, :)] = compute_avg_E_bound_numerical(Ef_grid(i), kPon_vals, P.kPoff, P.kEon, P.kEoff, n_states);
    end
    
    P.RE_val_bind_E = @(Ef_val) interpolate_E_bound(Ef_val, Ef_grid, avg_E_bound_grid, avg_Ser2P_grid);
    
    % Solve system - Step 1
    X_guess = 1e-6 * ones(N + N_PAS, 1);
    options = optimoptions('fsolve', 'Display', 'off', 'FunctionTolerance', 1e-8);
    
    P.FirstRun = true;
    P.is_unphysical = false;
    Ef_ss = 0;
    
    X_base = fsolve(@(xx) ode_dynamics_multipleE(xx, P), X_guess, options);
    
    if P.is_unphysical
        error('Unphysical result in Step 1');
    end
    
    % Solve system - Step 2 with adjusted kHon
    avg_E_bound = P.RE_val_bind_E(Ef_ss);
    P.FirstRun = false;
    P.kHon = kHon_base * avg_E_bound(end);
    
    
    X_final = fsolve(@(xx) ode_dynamics_multipleE(xx, P), X_base, options);
    
    R_sol = X_final(1:N);
    REH_sol = X_final(N+1:N+N_PAS);
end

function [avg_E, avg_S] = interpolate_E_bound(Ef_val, Ef_grid, E_grid, S_grid)
    % Fast 1D interpolation over pre-computed grid
    avg_E = interp1(Ef_grid, E_grid, Ef_val, 'linear', 'extrap');
    if nargout > 1
        avg_S = interp1(Ef_grid, S_grid, Ef_val, 'linear', 'extrap');
    end
end
