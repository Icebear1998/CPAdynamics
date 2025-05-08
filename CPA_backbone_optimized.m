global N PAS N_PAS Ef_ss kHon_ss;

% ------------ MODEL PARAMETERS ------------
L_a = 100;
P.k_in    = 2;
P.kEon    = 0.00025;
P.kEoff   = 10;
P.k_e     = 65 / L_a;
P.k_e2    = 30 / L_a;
P.E_total = 70000;
P.L_total = 100000;
P.Pol_total = 70000;
P.kHon    = 0.1; 
P.kHoff   = 0.0025; 
P.kc      = 0.8;
P.kPmin   = 2; 
P.kPmax   = 15; 

geneLength_bp = 25000;
PASposition   = 20000;
N      = floor(geneLength_bp / L_a);  % total nodes
PAS    = floor(PASposition / L_a);    % node index of PAS
N_PAS  = N - PAS + 1;                 % number of nodes at/after PAS

EBindingNumber = 4;

% Compute steady-state probabilities (assumed to return symbolic expressions)
[r_E_BeforePas] = compute_steady_states(P, EBindingNumber + 1);
disp('Done computing steady states');

% Precompute Kp values for positions 1 to PAS
Kp_vals = linspace(P.kPmax, P.kPmin, PAS); % Range of Kp

% Precompute RE_vals and P_vals numerically for all positions
RE_vals = zeros(EBindingNumber + 1, PAS);
%P_vals = zeros(EBindingNumber + 1, PAS);

% Convert symbolic expressions to numerical functions
for e = 1:EBindingNumber + 1
    RE_func = matlabFunction(r_E_BeforePas(e), 'Vars', {'kP', 'Ef'});
    %P_func = matlabFunction(r_P(e), 'Vars', {'kP', 'Ef'});
    for i = 1:PAS
        kP_val = Kp_vals(i);
        % Initial Ef guess (will be updated later)
        Ef_guess = P.E_total;
        RE_vals(e, i) = RE_func(kP_val, Ef_guess);
        %P_vals(e, i) = P_func(kP_val, Ef_guess);
    end
end
disp('Done computing EBindingNumber');

% Precompute average E bound per Pol II at each position
P.Avg_E = sum((0:EBindingNumber)' .* RE_vals, 1); % Avg E at each position
disp('Done computing Avg_E');

% Convert kHon_t_Ef to a numerical function
P.kHon_t_func = matlabFunction(r_k_AfterPas(1), 'Vars', {'Ef'});
disp('Done converting kHon_t_Ef to numerical function');

% ------------ SOLVE THE STEADY STATE ------------
X0 = zeros(N + N_PAS, 1);
X = fsolve(@(xx) ode_dynamics(xx, P), X0);
disp('Done solving steady state with fsolve');

% Extract solutions
R_sol   = X(1:N);
REH_sol = X(N+1 : N+N_PAS);

Ef_ss = Ef_ss(end);
%kHon_ss = kHon_ss;
disp(['Ef_ss: ', num2str(Ef_ss)]);
%disp(['kHon_ss: ', num2str(kHon_ss)]);

% Update RE_vals and P_vals with steady-state Ef_ss and R_sol
avg_E_bound = zeros(1, PAS);
%avg_P = zeros(1, PAS);

for i = 1:PAS
    total_R = R_sol(i);
    if total_R > 0
        % Update RE_vals and P_vals at position i
        for e = 1:EBindingNumber + 1
            RE_func = matlabFunction(r_E_BeforePas(e), 'Vars', {'kP', 'Ef'});
            RE_vals(e, i) = total_R * RE_func(Kp_vals(i), Ef_ss);
            %P_func = matlabFunction(r_P(e), 'Vars', {'kP', 'Ef'});
            %P_vals(e, i) = total_R * P_func(Kp_vals(i), Ef_ss);
        end
        % Compute average E bound
        avg_E_bound(i) = sum((0:EBindingNumber) .* RE_vals(:, i)') / total_R;
        % Compute average Ser2P
%         total_P_bound = sum((0:EBindingNumber) .* P_vals(:, i)');
%         total_P = sum(P_vals(:, i));
%         if total_P > 0
%             avg_P(i) = total_P_bound / total_P;
%         else
%             avg_P(i) = 0;
%         end
     else
         avg_E_bound(i) = 0;
         %avg_P(i) = 0;
     end
end
disp(['kHon_ss: ', num2str(P.kHon*avg_E_bound(end))])
disp('Done computing RE_vals and Ser2P');

% Plot Ser2P and Average E
figure;
hold on;
plot((1-PAS):0, avg_E_bound, 'b-', 'LineWidth', 2.5, 'DisplayName', 'AverageE');
%plot((1-PAS):0, avg_P, 'g-', 'LineWidth', 2.5, 'DisplayName', 'Ser2P');
legend({'AverageE', 'Ser2P'}, 'Location', 'best');
xlabel('Position');
ylabel('Average');
title('Ser2P and Average E Along Gene');
hold off;

% Re-solve the model system with new calculated kHon value
P.kHon = P.kHon * avg_E_bound(end);
X = fsolve(@(xx) ode_dynamics(xx, P), X);
disp('Done solving steady state with fsolve with calculated kHon value');

% Extract solutions
R_sol   = X(1:N);
REH_sol = X(N+1 : N+N_PAS);


% ------------ PLOT RESULTS ------------
% Time evolution plot
l_values = (1-PAS):(N-PAS);

figure;
hold on;
for e = 1:EBindingNumber+1
    plot((1-PAS):0, RE_vals(e,:), 'LineWidth',2);
end
plot(l_values, R_sol, 'b-', 'LineWidth', 2.5, 'DisplayName', 'R(l)');
plot(l_values, [zeros(PAS-1,1); REH_sol], 'g-', 'LineWidth', 2.5, 'DisplayName', 'REH(l)');
xlabel('Position');
ylabel('Total Pol II');
legend({'Total R', 'Total REH'}, 'Location', 'best');
title('Spatial Distribution of R and REH');
hold off;

% ------------ ODE DYNAMICS FUNCTION ------------
function dxdt = ode_dynamics(X, P)
    global N PAS Ef_ss kHon_ss;

    k_in   = P.k_in;
    k_e    = P.k_e;
    k_e2   = P.k_e2;
    kHoff_t= P.kHoff;
    kc_t   = P.kc;
    Avg_E  = P.Avg_E;
    %kHon_t_func = P.kHon_t_func;
    kHon_t = P.kHon; % Use target steady-state value for faster dynamics

    R   = X(1:N);
    REH = X(N+1:end);

    % Compute E_f numerically
    E_used = sum(R(1:PAS)'.* Avg_E);
    E_f = P.E_total - E_used;
    Ef_ss = E_f;
    if E_f < 0
        error('Negative E_f at t = %g (E_f = %g). Stopping simulation.', 0, E_f);
    end

    Pol_f = P.Pol_total - sum(R) - sum(REH);
%     kHon_t = kHon_t_func(E_f(end)); % Or use kHon_t = 0.0693 directly
%     kHon_ss = kHon_t;

    dxdt = zeros(length(X), 1);

    % Node 1
    n = 1;
    dxdt(n) = Pol_f * k_in - k_e * R(n);

    % Nodes 2 to PAS-1
    for n = 2:(PAS-1)
        dxdt(n) = k_e * R(n-1) - k_e * R(n);
    end

    % Node PAS
    n = PAS;
    j = n - PAS + 1;
    dxdt(n) = k_e * R(n-1) - k_e * R(n) - kHon_t * R(n) + kHoff_t * REH(j);
    dxdt(N+j) = -k_e2 * REH(j) + kHon_t * R(n) - kHoff_t * REH(j);

    % Nodes PAS+1 to N
    for n = (PAS+1):N
        j = n - PAS + 1;
        dxdt(n) = k_e * R(n-1) - k_e * R(n) - kHon_t * R(n) + kHoff_t * REH(j);
        dxdt(N+j) = k_e2 * REH(j-1) - k_e2 * REH(j) + kHon_t * R(n) - kHoff_t * REH(j) - kc_t * REH(j);
    end
end