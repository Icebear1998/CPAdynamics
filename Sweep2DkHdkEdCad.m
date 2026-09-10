%% Sweep2DkHdkEdCad.m
%  Full finite-rate R/RHE model: 2D contour map of CAD_50 as a function of:
%     kHd = kHoff / kHon   (PAS recognition dissociation constant)
%     kEd = kEoff / kEon   (E-factor binding dissociation constant)
%
%  Purpose: Compare CAD with the 400-800 bp experimental window across
%           the revised appendix's sensitivity ranges for the chosen M.
%
%  How kHd and kEd are constructed:
%     Read defaults and ranges from default_parameters.m (P.ranges).
%     kHd range = [min_kHoff/max_kHon,  max_kHoff/min_kHon]
%     kEd range = [min_kEoff/max_kEon,  max_kEoff/min_kEon]
%     For the actual simulation, on-rates are fixed at reference values
%     and off-rates are computed as kHoff = kHd * kHon_ref, etc.
%     This sweeps the full RATIO envelope: reconstructed off-rates can lie
%     outside their individual intervals because the on-rates stay fixed.
%     Full-model kinetics also depend on absolute rates. kEoff_engaged stays
%     fixed across this map, independently of the ordinary E off-rate.

%% ========== USER-CONFIGURABLE SECTION ==========
save_result = true;
EBindingNumber = 5;       % maximum E binding / phosphorylation capacity
nPoints = 4;              % Number of grid points per axis

% --- Base parameters and centralized appendix ranges ---
P = default_parameters();
P.kEoff_engaged = 0.05; % s^-1; trial engaged-E detachment with rapid EH disassembly
kHon_range  = P.ranges.kHon;
kHoff_range = P.ranges.kHoff;
kEon_range  = P.ranges.kEon;
kEoff_range = P.ranges.kEoff;

% --- Compute kHd and kEd ranges from individual rate ranges ---
kHd_min = P.ranges.kHd(1);
kHd_max = P.ranges.kHd(2);
kEd_min = P.ranges.kEd(1);
kEd_max = P.ranges.kEd(2);

kHd_values = logspace(log10(kHd_min), log10(kHd_max), nPoints);
kEd_values = logspace(log10(kEd_min), log10(kEd_max), nPoints);

fprintf('Derived sweep ranges:\n');
fprintf('  kHd = kHoff/kHon : [%.2g, %.2g] \n', kHd_min, kHd_max);
fprintf('  kEd = kEoff/kEon : [%.2g, %.2g] \n', kEd_min, kEd_max);

% --- CAD threshold ---
percent_cleavage = 50;    % CAD_50

%% ========== 2D SWEEP ==========
% --- Reference on-rates (held fixed; off-rates vary as kHoff = kHd * kHon_ref) ---
kHon_ref = P.kHon;   % reference PAS recognition on-rate  [1/s]
kEon_ref = P.kEon;   % reference E-factor on-rate          [1/(s·molecule)]

nH = length(kHd_values);
nE = length(kEd_values);
CAD_matrix = NaN(nE, nH);   % Rows: kEd, Cols: kHd
max_exit_cdf_matrix = NaN(nE, nH);
rhs_residual_matrix = NaN(nE, nH);

totalIter = nH * nE;
currentIter = 0;

fprintf('Starting full-model 2D sweep: kHd vs kEd for M = %d\n', EBindingNumber);
fprintf('  Engaged-E off-rate held fixed: %.3g s^-1\n', P.kEoff_engaged);
fprintf('  kHd range : [%.2g, %.2g], %d points  (kHon_ref = %.2g)\n', ...
    min(kHd_values), max(kHd_values), nH, kHon_ref);
fprintf('  kEd range : [%.2g, %.2g], %d points  (kEon_ref = %.2g)\n', ...
    min(kEd_values), max(kEd_values), nE, kEon_ref);
fprintf('  Total runs: %d\n\n', totalIter);
fprintf('  Realized kHoff at fixed kHon: [%.3g, %.3g] s^-1\n', ...
    kHd_min*kHon_ref, kHd_max*kHon_ref);
fprintf('  Realized kEoff at fixed kEon: [%.3g, %.3g] s^-1\n\n', ...
    kEd_min*kEon_ref, kEd_max*kEon_ref);

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
            [R_sol, RHE_sol, P_sim, full_details] = ...
                run_full_termination_simulation(P_run, EBindingNumber);
            [~, ~, cad, cleavage_diagnostics] = calculate_full_pas_cleavage_profile(R_sol, RHE_sol, P_sim, ...
                'PercentCleavage', percent_cleavage);
            CAD_matrix(i, j) = cad;
            max_exit_cdf_matrix(i, j) = cleavage_diagnostics.max_exit_cdf;
            rhs_residual_matrix(i, j) = full_details.rhs_max_abs;
            fprintf('CAD_%d = %.0f bp; within-window cleavage = %.2f%%\n', ...
                percent_cleavage, cad, 100*cleavage_diagnostics.max_exit_cdf);
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
max_exit_cdf_base = NaN;
try
    [R_base, RHE_base, P_base_sim] = run_full_termination_simulation(P, EBindingNumber);
    [~, ~, CAD_base, base_diagnostics] = calculate_full_pas_cleavage_profile(R_base, RHE_base, P_base_sim, ...
        'PercentCleavage', percent_cleavage);
    max_exit_cdf_base = base_diagnostics.max_exit_cdf;
    fprintf('  CAD_%d (base) = %.0f bp; within-window cleavage = %.2f%%\n\n', ...
        percent_cleavage, CAD_base, 100*max_exit_cdf_base);
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
fig1 = figure('Name', 'Full model CAD_{50} Contour: kHd vs kEd', ...
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
title(sprintf('Full model CAD_{%d} (M = %d, engaged E off = %.3g s^{-1})', ...
    percent_cleavage, EBindingNumber, P.kEoff_engaged), ...
    'FontSize', 15, 'FontWeight', 'bold');
cb = colorbar;
cb.Label.String = sprintf('CAD_{%d}  (bp)', percent_cleavage);
cb.Label.FontSize = 13;
colormap(cmap);

%% ========== SAVE ==========
if save_result
    data = struct();
    data.model_variant = 'full_kinetics_rapid_EH_disassembly';
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
    data.max_exit_cdf_matrix = max_exit_cdf_matrix;
    data.max_exit_cdf_base = max_exit_cdf_base;
    data.rhs_residual_matrix = rhs_residual_matrix;
    data.percent_cleavage = percent_cleavage;
    
    save_analysis_results('sweep_2D_kHd_kEd_CAD', data, P, ...
        'ExtraInfo', sprintf('Full_engagedOff%.3g', P.kEoff_engaged));
    fprintf('Results saved via save_analysis_results.\n');
end

fprintf('Done.\n');
