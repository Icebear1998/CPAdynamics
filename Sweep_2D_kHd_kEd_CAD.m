%% sweep_2D_kHd_kEd_CAD.m
%  2D contour map of CAD_50 as a function of:
%     kHd = kHoff / kHon   (PAS recognition dissociation constant)
%     kEd = kEoff / kEon   (E-factor binding dissociation constant)
%
%  Purpose: Parameter sensitivity analysis for n = 1 (EBindingNumber = 1)
%           to show that CAD remains above the 400–800 bp experimental
%           window across literature-plausible parameter ranges.
%
%  How kHd and kEd are constructed:
%     Specify literature ranges for kHon, kHoff, kEon, kEoff.
%     kHd range = [min_kHoff/max_kHon,  max_kHoff/min_kHon]
%     kEd range = [min_kEoff/max_kEon,  max_kEoff/min_kEon]
%     For the actual simulation, on-rates are fixed at reference values
%     and off-rates are computed as kHoff = kHd * kHon_ref, etc.

%% ========== USER-CONFIGURABLE SECTION ==========
save_result = true;
EBindingNumber = 1;       % n = 1 case
nPoints = 4;              % Number of grid points per axis

% --- Literature ranges for individual rates (UPDATE THESE) ---
kHon_range  = [0.5,  5];       % PAS recognition on-rate   [1/s]
kHoff_range = [0.5,  5];        % PAS recognition off-rate  [1/s]
kEon_range  = [1e-6, 5e-6];     % E-factor on-rate          [1/(s·molecule)]
kEoff_range = [0.1, 5];        % E-factor off-rate         [1/s]

% --- Compute kHd and kEd ranges from individual rate ranges ---
kHd_min = kHoff_range(1) / kHon_range(2);   % most favorable (tight binding)
kHd_max = kHoff_range(2) / kHon_range(1);   % most unfavorable (weak binding)
kEd_min = kEoff_range(1) / kEon_range(2);
kEd_max = kEoff_range(2) / kEon_range(1);

kHd_values = logspace(log10(kHd_min), log10(kHd_max), nPoints);
kEd_values = logspace(log10(kEd_min), log10(kEd_max), nPoints);

fprintf('Derived sweep ranges:\n');
fprintf('  kHd = kHoff/kHon : [%.2g, %.2g] \n', kHd_min, kHd_max);
fprintf('  kEd = kEoff/kEon : [%.2g, %.2g] \n', kEd_min, kEd_max);

% --- CAD threshold ---
percent_cleavage = 50;    % CAD_50

%% ========== BASE PARAMETERS ==========
P = default_parameters();



%% ========== GLOBAL VARIABLES ==========
global N PAS N_PAS Ef_ss;
Ef_ss = 0;

%% ========== 2D SWEEP ==========
% --- Reference on-rates (held fixed; off-rates vary as kHoff = kHd * kHon_ref) ---
kHon_ref = P.kHon;   % reference PAS recognition on-rate  [1/s]
kEon_ref = P.kEon;   % reference E-factor on-rate          [1/(s·molecule)]

nH = length(kHd_values);
nE = length(kEd_values);
CAD_matrix = NaN(nE, nH);   % Rows: kEd, Cols: kHd

totalIter = nH * nE;
currentIter = 0;

fprintf('Starting 2D sweep: kHd vs kEd for n = %d\n', EBindingNumber);
fprintf('  kHd range : [%.2g, %.2g], %d points  (kHon_ref = %.2g)\n', ...
    min(kHd_values), max(kHd_values), nH, kHon_ref);
fprintf('  kEd range : [%.2g, %.2g], %d points  (kEon_ref = %.2g)\n', ...
    min(kEd_values), max(kEd_values), nE, kEon_ref);
fprintf('  Total runs: %d\n\n', totalIter);

tic;
for i = 1:nE
    for j = 1:nH
        currentIter = currentIter + 1;

        % Compute actual off-rates from ratios and reference on-rates
        P_run = P;
        P_run.kEon  = kEon_ref;
        P_run.kEoff = kEd_values(i) * kEon_ref;    % kEoff = kEd * kEon_ref
        P_run.kHon  = kHon_ref;
        P_run.kHoff = kHd_values(j) * kHon_ref;    % kHoff = kHd * kHon_ref

        fprintf('  [%d/%d] kEd=%.2g (kEoff=%.2g), kHd=%.2g (kHoff=%.2g) ... ', ...
            currentIter, totalIter, ...
            kEd_values(i), P_run.kEoff, ...
            kHd_values(j), P_run.kHoff);

        try
            [R_sol, REH_sol, P_sim] = run_termination_simulation(P_run, EBindingNumber);
            [~, ~, cad] = calculate_pas_cleavage_profile(R_sol, REH_sol, P_sim, ...
                'PercentCleavage', percent_cleavage);
            CAD_matrix(i, j) = cad;
            fprintf('CAD_%d = %.0f bp\n', percent_cleavage, cad);
        catch ME
            fprintf('FAILED: %s\n', ME.message);
            CAD_matrix(i, j) = NaN;
        end
    end
end
elapsed = toc;
fprintf('\nSweep completed in %.1f s\n\n', elapsed);

%% ========== BASE-PARAMETER CAD ==========
% Dissociation constants at the reference (base) parameter values
kHd_base = P.kHoff / P.kHon;   % = kHoff_ref / kHon_ref
kEd_base = P.kEoff / P.kEon;   % = kEoff_ref / kEon_ref

fprintf('Running simulation at base parameters (kHd=%.3g, kEd=%.3g) ...\n', kHd_base, kEd_base);
try
    [R_base, REH_base, P_base_sim] = run_termination_simulation(P, EBindingNumber);
    [~, ~, CAD_base] = calculate_pas_cleavage_profile(R_base, REH_base, P_base_sim, ...
        'PercentCleavage', percent_cleavage);
    fprintf('  CAD_%d (base) = %.0f bp\n\n', percent_cleavage, CAD_base);
catch ME
    fprintf('  FAILED: %s\n\n', ME.message);
    CAD_base = NaN;
end

%% ========== PLOTTING ==========
% --- Version-safe colormap (turbo requires R2020b+) ---
if exist('turbo', 'builtin') || exist('turbo', 'file')
    cmap = turbo(256);
else
    cmap = jet(256);
end

% --- Figure 1: Filled contour ---
fig1 = figure('Name', 'CAD_{50} Contour: kHd vs kEd', ...
    'Position', [100 100 900 700]);

[X, Y] = meshgrid(kHd_values, kEd_values);

nLevels = 20;
contourf(X, Y, CAD_matrix, nLevels, 'LineStyle', 'none');
hold on;

% Draw experimental target band (400–800 bp) as contour lines
[C400, h400] = contour(X, Y, CAD_matrix, [400 400], 'r--', 'LineWidth', 2.5);
[C800, h800] = contour(X, Y, CAD_matrix, [800 800], 'r--', 'LineWidth', 2.5);

% Mark base-parameter point
plot(kHd_base, kEd_base, 'pw', 'MarkerSize', 18, 'MarkerFaceColor', 'w', 'LineWidth', 1.5);
if ~isnan(CAD_base)
    text(kHd_base, kEd_base, sprintf('  CAD_{%d} = %.0f bp', percent_cleavage, CAD_base), ...
        'Color', 'w', 'FontSize', 11, 'FontWeight', 'bold', ...
        'VerticalAlignment', 'middle');
end
hold off;

set(gca, 'XScale', 'log', 'YScale', 'log', 'FontSize', 12);
xlabel('k_{Hd} = k_{Hoff} / k_{Hon}', 'FontSize', 14, 'FontWeight', 'bold');
ylabel('k_{Ed} = k_{Eoff} / k_{Eon}', 'FontSize', 14, 'FontWeight', 'bold');
title(sprintf('CAD_{%d} Contour Map  (n = %d)', ...
    percent_cleavage, EBindingNumber), ...
    'FontSize', 15, 'FontWeight', 'bold');
cb = colorbar;
cb.Label.String = sprintf('CAD_{%d}  (bp)', percent_cleavage);
cb.Label.FontSize = 13;
colormap(cmap);

%% ========== SAVE ==========
if save_result
    data.EBindingNumber = EBindingNumber;
    data.kHd_values = kHd_values;
    data.kEd_values = kEd_values;
    data.kHon_ref = kHon_ref;
    data.kEon_ref = kEon_ref;
    data.kHd_base = kHd_base;
    data.kEd_base = kEd_base;
    data.CAD_base = CAD_base;
    data.kHon_range  = kHon_range;
    data.kHoff_range = kHoff_range;
    data.kEon_range  = kEon_range;
    data.kEoff_range = kEoff_range;
    data.CAD_matrix = CAD_matrix;
    data.percent_cleavage = percent_cleavage;
    
    save_analysis_results('sweep_2D_kHd_kEd_CAD', data, P);
    fprintf('Results saved via save_analysis_results.\n');
end

fprintf('Done.\n');
