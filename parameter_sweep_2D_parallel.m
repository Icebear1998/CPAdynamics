% 2D Parameter Sweep Setup
global N PAS N_PAS Ef_ss;
syms Ef real;

% Model parameters (base values) - Define these once at the top
BASE_PARAMS = struct();
BASE_PARAMS.k_in = 2;
BASE_PARAMS.kEon = 0.00025;
BASE_PARAMS.kEoff = 10;
BASE_PARAMS.k_e = 65/100;  % Using L_a = 100
BASE_PARAMS.k_e2 = 30/100;
BASE_PARAMS.E_total = 70000;
BASE_PARAMS.L_total = 100000;
BASE_PARAMS.Pol_total = 70000;
BASE_PARAMS.kHon = 0.05;
BASE_PARAMS.kHoff = 0.0025;
BASE_PARAMS.kc = 0.8;
BASE_PARAMS.kPon_min = 0.01; % at TSS
BASE_PARAMS.kPon_max = 1; % at PAS
BASE_PARAMS.kPoff_min = 0.1; % at PAS
BASE_PARAMS.kPoff_max = 2; % at TSS
BASE_PARAMS.kPoff_const = 1;
BASE_PARAMS.kPon_const = 1;

L_a = 100;
geneLength_bp = 25000;
PASposition = 20000;
N = floor(geneLength_bp / L_a);
PAS = floor(PASposition / L_a);
N_PAS = N - PAS + 1;
EBindingNumber = 3;

% Define parameter pairs for 2D sweep
param_pairs = {
%     {'E_total', 'kEoff'}, ...
%     {'E_total', 'kEon'}, ...
%    {'kc', 'kHoff'}
%     {'kPmin', 'kPmax'}, ...
%     {'k_e', 'E_total'}
%     {'k_e', 'kHon'}
%     {'kPon_min', 'kPon_max'}
    {'E_total', 'kHon'}
};

% Loop through all parameter pairs
for pair_idx = 1:length(param_pairs)
    param1 = param_pairs{pair_idx}{1};
    param2 = param_pairs{pair_idx}{2};
    
    % Define ranges for each parameter
    switch param1
        case 'k_e'
            param1_values = 35/L_a:10/L_a:95/L_a;
        case 'E_total'
            param1_values = 40000:-10000:10000;
        case 'kc'
            param1_values = 0.2:0.1:1.0;
        case 'kPmin'
            param1_values = logspace(log10(0.05), log10(0.4), 6);
        case 'EBindingNumber'
            param1_values = 2:5; % Range for EBindingNumber
        case 'kPon_min'
            param1_values = logspace(log10(0.001), log10(1), 6);
        otherwise
            error('Invalid parameter1 selected');
    end

    switch param2
        case 'k_e'
            param2_values = 35/L_a:10/L_a:95/L_a;
        case 'kEon'
            param2_values = logspace(-5, -2, 6);
        case 'kEoff'
            param2_values = logspace(0, log10(100), 6);
        case 'kHon'
            %param2_values = logspace(log10(0.01), log10(0.05), 5);
            param2_values = 0.04:-0.01:0.01;
        case 'kHoff'
            param2_values = logspace(-3, -1.5, 6);
        case 'kPmax'
            param2_values = logspace(1, 2, 6);
        case 'kPmin'
            param2_values = logspace(log10(0.05), log10(0.4), 6);
        case 'E_total'
            param2_values = 50000:10000:100000;
        case 'kPon_max'
            param2_values = logspace(log10(1), log10(20), 6);
        otherwise
            error('Invalid parameter2 selected');
    end

    % Get default values for this parameter pair
    default1 = BASE_PARAMS.(param1);
    default2 = BASE_PARAMS.(param2);

    % Initialize matrix for cutoff values
    cutoff_matrix = zeros(length(param2_values), length(param1_values));

    % Define options before the loop
    options = optimoptions('fsolve', ...
    'Algorithm', 'trust-region', ... % 'trust-region' is often more robust
    'Display', 'none', ...           % Suppress solver output in the command window
    'FunctionTolerance', 1e-10, 'TolX', 1e-10);      % A tighter tolerance for accuracy
    % 2D Parameter Sweep
    for i = 1:length(param2_values)
        % Initialize the guess for the first point in the inner loop
        X_guess = 1e-2 *ones(N + N_PAS, 1);  
        for j = 1:length(param1_values)
            % CRITICAL: Reset ALL parameters and global variables for each iteration
            P = BASE_PARAMS;  % Reset to base parameters
            Ef_ss = 0;        % Reset global variable
            
            % Set the two parameters being varied
            P.(param1) = param1_values(j);
            P.(param2) = param2_values(i);
            
            % Store original kHon value for this parameter combination
            kHon_original = P.kHon;
            
            [r_E_BeforePas] = compute_steady_states(P, EBindingNumber+1);
            disp('done compute steady states');

            kPon_vals = linspace(P.kPon_min, P.kPon_max, PAS); % Range of kPon increases linearly
            %kPoff_vals = linspace(P.kPoff_max, P.kPoff_min, PAS); % Range of kPoff decreases linear
            
            RE_vals = sym(zeros(EBindingNumber + 1, N));

            for e = 1:EBindingNumber + 1
                for idx = 1:PAS
                    kPon_val = kPon_vals(idx);
                    RE_vals(e, idx) = subs(r_E_BeforePas(e), {'kPon', 'kPoff'}, {kPon_val, P.kPoff_const});
                end
                for idx = PAS+1:N
                    RE_vals(e, idx) = subs(r_E_BeforePas(e), {'kPon', 'kPoff'}, {P.kPon_max, P.kPoff_min});
                end
            end  

            P.RE_val_bind_E = matlabFunction(simplify(sum(sym(1:EBindingNumber)' .* RE_vals(2:end, :), 1)), 'Vars', {Ef});

            P.FirstRun = true;
            P.is_unphysical = false;
            %X0 = zeros(N + N_PAS, 1);

            X = fsolve(@(xx) ode_dynamics_multipleE(xx, P), X_guess, options);
            if P.is_unphysical
                cutoff_matrix(i,j) = 0;
                continue;
            end

            avg_E_bound = P.RE_val_bind_E(Ef_ss);

            P.FirstRun = false;
            % Update kHon based on the binding, but use the original value as base
            disp(kHon_original);
            P.kHon = kHon_original * avg_E_bound(end);
            X = fsolve(@(xx) ode_dynamics_multipleE(xx, P), X, options);
            
            % Store the final solution to use as the next initial guess
            X_guess = X;

            R_sol = X(1:N);
            REH_sol = X(N+1 : N+N_PAS);
            
            % Calculate cutoff position in the gene using interpolation
            ratio = (REH_sol(1:end) + R_sol(PAS:end)) / R_sol(PAS-1);
            node_indices = 1:(length(ratio));
            smoothed_ratio = smooth(ratio, 5, 'moving'); % 5-point moving average
            cutoff_matrix(i,j) = interp1(smoothed_ratio, node_indices, 0.85, 'linear', 'extrap') * L_a;
        end
    end

%     % Plot as contour plot
%     figure;
%     hold on;
%     [XX, YY] = meshgrid(param1_values, param2_values);
%     contourf(XX, YY, cutoff_matrix, 20, 'LineColor', 'none');
%     if ismember(param1, {'kPmin', 'kPmax', 'kEon', 'kEoff', 'kHon', 'kHoff', 'kPon_min', 'kPon_max'})
%         set(gca, 'XScale', 'log');
%     end
%     if ismember(param2, {'kPmin', 'kPmax', 'kEon', 'kEoff', 'kHon', 'kHoff', 'kPon_min', 'kPon_max'})
%         set(gca, 'YScale', 'log');
%     end
%     colorbar;
%     xlabel([strrep(param1, '_', '\_'), ' Value']);
%     ylabel([strrep(param2, '_', '\_'), ' Value']);
%     title(['2D Contour Plot of CPA Cutoff (EBindingNumber=', num2str(EBindingNumber), ')']);
%     hold on;
%     %plot(default1, default2, 'ro', 'MarkerSize', 10, 'LineWidth', 2, 'DisplayName', 'Default');
%     legend('show');
%     grid on;

%     % Save the plot
%     filename = sprintf('CPA_Cutoff_2D_%s_%s_2.png', param1, param2);
%     saveas(gcf, filename);
%     close(gcf);
%     
%     % Display the cutoff value at default parameters for verification
%     [~, idx1] = min(abs(param1_values - default1));
%     [~, idx2] = min(abs(param2_values - default2));
%     fprintf('Default cutoff for %s vs %s: %.2f\n', param1, param2, cutoff_matrix(idx2, idx1));

    % Plot as line plot
    figure('Position', [100, 100, 800, 600]); % Set figure size for clarity
    hold on;
    
    % Define a set of colors from a colormap
    colors = lines(size(cutoff_matrix, 1)); % Use 'lines' colormap for distinct colors
    
    % Plot each row of cutoff_matrix as a line (one line per param2_values)
    for i = 1:size(cutoff_matrix, 1)
        % Replace NaN with a sentinel value (e.g., -1) for plotting
        plot_data = cutoff_matrix(i, :);
        plot_data(isnan(plot_data)) = -1; % Use -1 as a placeholder (adjust if needed)
        
        plot(param1_values, plot_data, 'o-', 'LineWidth', 2, 'MarkerSize', 6, ...
             'Color', colors(i, :), ... % Use colormap-based color
             'DisplayName', sprintf('%s = %.2g', param2, param2_values(i)));
    end
    
    % Apply log scale if applicable
    if ismember(param1, {'kPmin', 'kPmax', 'kEon', 'kEoff', 'kHon', 'kHoff', 'kPon_min', 'kPon_max'})
        set(gca, 'XScale', 'log');
    end
    % Note: Y-scale log may not be ideal for cutoff values; comment out if unnecessary
    % if ismember(param2, {'kPmin', 'kPmax', 'kEon', 'kEoff', 'kHon', 'kHoff', 'kPon_min', 'kPon_max'})
    %     set(gca, 'YScale', 'log');
    % end
    
    % Customize plot
    xlabel([strrep(param1, '_', '\_'), ' Value'], 'FontSize', 12);
    ylabel('CPA Cutoff Position (bp)', 'FontSize', 12);
    title(['CPA Cutoff vs ', strrep(param1, '_', '\_'), ' by ', strrep(param2, '_', '\_'), ...
           ' (EBindingNumber=', num2str(EBindingNumber), ')'], 'FontSize', 14, 'FontWeight', 'bold');
    grid on;
    legend('show', 'Location', 'best', 'FontSize', 10);
    set(gca, 'FontSize', 10);
    box on;
    hold off;
    
%     % Add current date and time (10:42 AM EDT 29-Jul-2025)
%     annotation('textbox', [0.1, 0.95, 0.8, 0.05], ...
%                'String', 'Generated at 10:42 AM EDT 29-Jul-2025', ...
%                'EdgeColor', 'none', 'HorizontalAlignment', 'center', 'FontSize', 8);
%     
%     % Save the plot
%     filename = sprintf('CPA_Cutoff_LinePlot_EBinding%d_%s_%s.png', EBindingNumber, param1, param2);
%     saveas(gcf, filename, 'png');
%     print(gcf, filename, '-dpng', '-r300'); % High-resolution save
%     close(gcf);
end