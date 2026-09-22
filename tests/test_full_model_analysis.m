function tests = test_full_model_analysis
% Full-model migration: shared pools, microstate observables and genome solve.
tests = functiontests(localfunctions);
end

function setupOnce(testCase)
testCase.TestData.originalPath = path;
addpath(fileparts(fileparts(mfilename('fullpath'))));
P = default_parameters();
P.geneLength_bp = 800;
P.PASposition = 400;
P.k_in = 0.01;
P.kEon = 0.01;
P.Pol_total = 100;
P.E_total = 200;
testCase.TestData.P = P;
end

function teardownOnce(testCase)
path(testCase.TestData.originalPath);
end

function testFixedPoolsRecoverClosedSolution(testCase)
P = testCase.TestData.P;
for M = [1 3]
    [R, RHE, P_closed, D_closed] = run_full_termination_simulation(P, M);
    pools = [P_closed.Pol_free_ss P_closed.Ef_ss];
    [R_fixed, RHE_fixed, P_fixed, D_fixed] = run_full_termination_simulation(P, M, 'FreePools', pools);
    verifyEqual(testCase, R_fixed, R, 'AbsTol', 1e-9);
    verifyEqual(testCase, RHE_fixed, RHE, 'AbsTol', 1e-9);
    verifyEqual(testCase, D_fixed.E_bound, D_closed.E_bound, 'AbsTol', 1e-9);
    verifyEqual(testCase, P_fixed.kHon, P.kHon);
    verifyEqual(testCase, [P_fixed.Pol_free_ss P_fixed.Ef_ss], pools);
    verifyEqual(testCase, D_fixed.pool_mode, 'fixed_free_pools');
    verifyTrue(testCase, isnan(D_fixed.E_conservation_residual));
    verifyLessThan(testCase, D_fixed.rhs_max_abs, 1e-6);
end
end

function testFixedPoolsIgnoreLocalTotalsAndScaleWithFreePol(testCase)
P = testCase.TestData.P;
[R, RHE, ~, D] = run_full_termination_simulation(P, 3, 'FreePools', [10 20]);
P.Pol_total = 1;
P.E_total = 1; % Prescribed shared pools may exceed a single gene's totals.
[R_same, RHE_same] = run_full_termination_simulation(P, 3, 'FreePools', [10 20]);
verifyEqual(testCase, R_same, R, 'AbsTol', 1e-10);
verifyEqual(testCase, RHE_same, RHE, 'AbsTol', 1e-10);
[R_twice, RHE_twice, ~, D_twice] = run_full_termination_simulation(P, 3, 'FreePools', [20 20]);
verifyEqual(testCase, R_twice, 2*R, 'AbsTol', 1e-10);
verifyEqual(testCase, RHE_twice, 2*RHE, 'AbsTol', 1e-10);
verifyEqual(testCase, D_twice.E_bound, 2*D.E_bound, 'AbsTol', 1e-10);
verifyEqual(testCase, D_twice.avg_E_bound, D.avg_E_bound, 'AbsTol', 1e-10);
end

function testMicrostateAveragesConserveBoundResources(testCase)
P = testCase.TestData.P;
[R, RHE, P_sim, D] = run_full_termination_simulation(P, 3);
populations = R+[zeros(P_sim.PAS-1, 1); RHE];
verifySize(testCase, D.avg_E_bound, [P_sim.N 1]);
verifyEqual(testCase, sum(D.avg_E_bound.*populations), D.E_bound, 'AbsTol', 1e-9);
verifyEqual(testCase, sum(D.avg_Ser2P.*populations), ...
    sum(D.R_micro*D.R_states(:, 1))+sum(D.RHE_micro*D.RHE_states(:, 1)), 'AbsTol', 1e-9);
verifyGreaterThanOrEqual(testCase, D.avg_E_bound, zeros(P_sim.N, 1));
verifyLessThanOrEqual(testCase, D.avg_E_bound, D.avg_Ser2P+1e-12);
verifyLessThanOrEqual(testCase, D.avg_Ser2P, 3*ones(P_sim.N, 1)+1e-12);

[~, RHE_zero, ~, D_zero] = run_full_termination_simulation(P, 3, 'FreePools', [10 0]);
verifyEqual(testCase, RHE_zero, zeros(size(RHE_zero)), 'AbsTol', 1e-12);
verifyEqual(testCase, D_zero.E_bound, 0, 'AbsTol', 1e-12);
[R_empty, RHE_empty, ~, D_empty] = run_full_termination_simulation(P, 3, 'FreePools', [0 20]);
verifyEqual(testCase, [R_empty; RHE_empty], zeros(size([R_empty; RHE_empty])));
verifyTrue(testCase, all(isnan(D_empty.avg_E_bound)));
verifyTrue(testCase, all(isnan(D_empty.avg_Ser2P)));
verifyEqual(testCase, D_empty.E_bound, 0);
end

function testGenomeSolveMatchesAnalyticConservation(testCase)
I = analytic_interpolation();
A = run_genome_analysis(I);
S = A.solution;
R_expected = 100/(1+10*sum(A.integration.weights.*(1+A.integration.lengths/1000)));
E_expected = 200/(1+10*R_expected*0.01);
verifyEqual(testCase, S.R_free, R_expected, 'AbsTol', 1e-10);
verifyEqual(testCase, S.E_free, E_expected, 'AbsTol', 1e-9);
verifyEqual(testCase, S.conservation_residuals, [0 0], 'AbsTol', 1e-8);
verifyTrue(testCase, all(cellfun(@isempty, A.TCD.error_messages)));
% Check one of the actual CAD profiles against the fixed-pool full model.
P = I.original_grid.base_parameters;
P.PASposition = A.TCD.gene_lengths(1);
P.geneLength_bp = P.PASposition+I.original_grid.after_PAS_length;
[R, H, P] = run_full_termination_simulation(P, P.EBindingNumber, ...
    'FreePools', [S.R_free S.E_free]);
[cdf, ~] = calculate_full_pas_cleavage_profile(R, H, P, 'PercentCleavage', 0);
verifyEqual(testCase, A.TCD.termination_profiles{1}.profile, cdf, 'AbsTol', 1e-10);

I.original_grid.E_total_base = 0;
A_zero = run_genome_analysis(I);
verifyEqual(testCase, A_zero.solution.E_free, 0);
I.original_grid.E_total_base = 200;
I.functions.E_occupied_interp = @(r, e, L) zeros(size(r));
A_unbound = run_genome_analysis(I);
verifyEqual(testCase, A_unbound.solution.E_free, 200);
end

function testGenomeRejectsIncompatibleOrOutOfRangeData(testCase)
I = analytic_interpolation();
I.metadata.model_variant = 'incompatible';
verifyError(testCase, @() run_genome_analysis(I), ...
    'GeneLengthAnalyze:IncompatibleInterpolation');
I = analytic_interpolation();
I.original_grid.R_total_base = 101;
verifyError(testCase, @() run_genome_analysis(I), ...
    'solve_full_genome_pools:OutsideGrid');
end

function I = analytic_interpolation()
I.metadata.model_variant = 'full_kinetics_rapid_EH_disassembly';
I.metadata.pool_mode = 'fixed_free_pools';
I.original_grid.R_free_range = [0 100];
I.original_grid.E_free_range = [0 200];
I.original_grid.L_range = [100 1000];
I.original_grid.R_total_base = 100;
I.original_grid.E_total_base = 200;
I.original_grid.num_active_genes = 10;
I.original_grid.after_PAS_length = 500;
I.original_grid.base_parameters = default_parameters();
I.original_grid.base_parameters.EBindingNumber = 1;
I.functions.R_occupied_interp = @(r, e, L) r.*(1+L/1000);
I.functions.E_occupied_interp = @(r, e, L) r.*e*0.01;
I.functions.gene_length_pdf = @(L) ones(size(L))/900;
end

function A = run_genome_analysis(interpolation_results)
% Test the consolidated solver through the normal analysis entry point.
output_root = tempname;
mkdir(output_root);
old_root = getenv('CPAD_RESULTS_ROOT');
old_visibility = get(groot, 'defaultFigureVisible');
old_figures = findall(groot, 'Type', 'figure');
cleanup = onCleanup(@() restore_environment(output_root, old_root, old_visibility, old_figures)); %#ok<NASGU>
setenv('CPAD_RESULTS_ROOT', output_root);
set(groot, 'defaultFigureVisible', 'off');
output_dir = cpad_analysis_output_dir('GeneLengthAnalysis');
save(fullfile(output_dir, 'full_gene_length_interpolation_test.mat'), 'interpolation_results');
A = run_analysis();
end

function A = run_analysis()
evalc('GeneLengthAnalyze');
A = analysis_results;
end

function restore_environment(output_root, old_root, old_visibility, old_figures)
setenv('CPAD_RESULTS_ROOT', old_root);
set(groot, 'defaultFigureVisible', old_visibility);
close(setdiff(findall(groot, 'Type', 'figure'), old_figures));
rmdir(output_root, 's');
end
