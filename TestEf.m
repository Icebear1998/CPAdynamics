global N PAS N_PAS Ef_ss kHon_ss is_unphysical;
syms Ef real;
% ------------ MODEL PARAMETERS ------------
L_a = 100;
P.k_in    = 2;
P.kEon    = 0.00025;
P.kEoff   = 10;
P.k_e     = 65/L_a;
P.k_e2    = 30/L_a;
P.E_total = 200000;
P.L_total = 100000;
P.Pol_total = 70000;
P.kHon = 0.1; % based on typical k bind and estimated J factor for H.
P.kHoff = 0.0025; 
P.kc = 0.8; %not sure
P.kPmin   = 0.1; %not sure
P.kPmax   = 40; %not sure

geneLength_bp = 25000;
PASposition   = 20000;
N      = floor(geneLength_bp / L_a);  % total nodes
PAS    = floor(PASposition / L_a);  % node index of PAS
N_PAS  = N - PAS + 1;                 % number of nodes at/after PAS
Ef_ss = 0;

EBindingNumber = 3; 
[r_E_BeforePas] = compute_steady_states(P, EBindingNumber + 1); 
disp('done compute steady states');

KpOver_vals = linspace(1/P.kPmax, 1/P.kPmin, PAS); % Range of Kp for kPon increases linearly
Kp_vals = 1./KpOver_vals;
RE_vals = sym(zeros(EBindingNumber + 1, PAS));

for e = 1:EBindingNumber+1
    for i = 1:length(Kp_vals)
        kP_val = Kp_vals(i);
        RE_vals(e, i) = subs(r_E_BeforePas(e), {'kP'}, {kP_val});
    end
end
disp('done compute EBindingNumber');

P.RE_val_bind_E = matlabFunction(simplify(sum(sym(1:EBindingNumber)' .* RE_vals(2:end, :), 1)), 'Vars', {Ef});
disp('done compute RE_val_bind_E');
% ------------ SOLVE THE STEADY STATE ------------
X0 = zeros(N + N_PAS, 1);
% Clear persistent variables before running fsolve
clear ode_dynamics; % Ensures persistent variables are reset
X = fsolve(@(xx) ode_dynamics(xx, P), X0, optimoptions('fsolve', 'Display', 'off'));
disp('done compute fsolve');
%Extract solutions
R_sol   = X(1:N);
REH_sol = X(N+1 : N+N_PAS);

% Compute avg_E_bound at PAS
avg_E_bound = P.RE_val_bind_E(Ef_ss);
disp('done compute RE_vals');

% Retrieve the logged values from ode_dynamics
[total_RE_vals, total_R, E_used_vals] = ode_dynamics('get_logs');
iterations = 1:length(total_RE_vals);

% Normalize quantities for plotting
max_total_RE = max(total_RE_vals, [], 'omitnan');
max_total_R = max(total_R, [], 'omitnan');
max_E_used = max(E_used_vals, [], 'omitnan');

if max_total_RE == 0, max_total_RE = 1; end
if max_total_R == 0, max_total_R = 1; end
if max_E_used == 0, max_E_used = 1; end

normalized_total_RE = total_RE_vals / max_total_RE;
normalized_total_R = total_R / max_total_R;
normalized_E_used = E_used_vals / max_E_used;

% Find indices of maximum values for labeling
[~, max_idx_total_RE] = max(normalized_total_RE);
[~, max_idx_total_R] = max(normalized_total_R);
[~, max_idx_E_used] = max(normalized_E_used);

% ------------ PLOT RESULTS ------------
% % Plot 1: Original Time Evolution Plot
% l_values = (1-PAS):(N-PAS);
% 
% figure; hold on;
% for e = 1:EBindingNumber+1
%     plot((1-PAS):0, RE_vals(e,:), 'LineWidth',2);
% end
% 
% plot(l_values, R_sol, 'b-','LineWidth',2.5, 'DisplayName','R(l)');
% plot(l_values, [zeros(PAS-1,1);REH_sol], 'r-','LineWidth',2.5, 'DisplayName','REH(l)');
% xlabel('Time'); ylabel('Total Pol II');
% legend({'Total R', 'Total REH'}, 'Location', 'best');
% title('Time Evolution of R and REH');
% hold off;

% Plot 2: Normalized quantities over fsolve iterations
figure; hold on;
plot(iterations, normalized_total_R, 'o-', 'LineWidth', 2, 'MarkerSize', 8, 'DisplayName', sprintf('Total R (max=%.1f)', max_total_R));
plot(iterations, normalized_total_RE, 's-', 'LineWidth', 2, 'MarkerSize', 8, 'DisplayName', sprintf('Total RE_val_bind_E (max=%.1f)', max_total_RE));
plot(iterations, normalized_E_used, 'd-', 'LineWidth', 2, 'MarkerSize', 8, 'DisplayName', sprintf('E_used (max=%.1f)', max_E_used));

% Label maximum values
text(iterations(max_idx_total_R), normalized_total_R(max_idx_total_R), sprintf('%.1f', max_total_R), 'FontSize', 8, 'HorizontalAlignment', 'center', 'VerticalAlignment', 'bottom');
text(iterations(max_idx_total_RE), normalized_total_RE(max_idx_total_RE), sprintf('%.1f', max_total_RE), 'FontSize', 8, 'HorizontalAlignment', 'center', 'VerticalAlignment', 'bottom');
text(iterations(max_idx_E_used), normalized_E_used(max_idx_E_used), sprintf('%.1f', max_E_used), 'FontSize', 8, 'HorizontalAlignment', 'center', 'VerticalAlignment', 'bottom');

grid on;
xlabel('fsolve Iteration');
ylabel('Normalized Values');
title('Normalized Quantities Over fsolve Iterations');
legend('Location', 'best');
xlim([min(iterations)-1, max(iterations)+1]);
ylim([-0.1, 1.1]);
hold off;

%% ------------ ODE DYNAMICS FUNCTION ------------
function varargout = ode_dynamics(X, P)
global N PAS Ef_ss kHon_ss is_unphysical;

persistent total_RE_vals total_R E_used_vals iteration_counter;

% Initialize persistent variables
if isempty(iteration_counter)
    iteration_counter = 0;
    total_RE_vals = [];
    total_R = [];
    E_used_vals = [];
end

% Handle the 'get_logs' command to return logged data
if ischar(X) && strcmp(X, 'get_logs')
    varargout = {total_RE_vals, total_R, E_used_vals};
    return;
end

k_in   = P.k_in;
k_e    = P.k_e;
k_e2   = P.k_e2;
kHoff_t= P.kHoff;
kc_t   = P.kc;
RE_val_bind_E = P.RE_val_bind_E;
kHon_t = P.kHon;

R   = X(1:N);
REH = X(N+1:end);

% Set initial Ef_ss once (using P.E_total as a starting guess)
if Ef_ss == 0
    Ef_ss = P.E_total; % Initial guess
end

% Convert symbolic expression to a numerical function
E_used = sum(R(1:PAS)' .* RE_val_bind_E(Ef_ss));
E_f = P.E_total - E_used;

% Log the values before updating Ef_ss
iteration_counter = iteration_counter + 1;
total_RE_vals(iteration_counter) = sum(RE_val_bind_E(Ef_ss));
total_R(iteration_counter) = sum(R(1:PAS));
E_used_vals(iteration_counter) = E_used;
disp({total_RE_vals(iteration_counter), total_R(iteration_counter), E_used_vals(iteration_counter)});

Ef_ss = E_f;

% Check if E_f < 0 without reassigning Ef_ss in the output
if E_f < 0
    is_unphysical = true; % Set flag instead of throwing error
    dxdt = 1e6 * ones(length(X), 1); % Large residual to signal fsolve to stop
    varargout = {dxdt};
    return;
end

Pol_f = P.Pol_total - sum(R) - sum(REH);
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

varargout = {dxdt};
end