function tests = test_gene_length_grid_validation
% Successful solver rows must survive harmless signed-zero roundoff.
tests = functiontests(localfunctions);
end

function setupOnce(testCase)
testCase.TestData.originalPath = path;
addpath(fileparts(fileparts(mfilename('fullpath'))));
end

function teardownOnce(testCase)
path(testCase.TestData.originalPath);
end

function testSuccessful250PointGridWith32NegativeRoundoffRows(testCase)
D = example_grid();
indices = find(D.R_free_vec > 0 & D.E_free_vec == 0);
indices = indices(1:32);
D.E_occupied_vec(indices) = -1e-12;
verifyEqual(testCase, nnz(D.E_occupied_vec >= 0), 218);
[R, E, diagnostics] = build_grid(D);
verifySize(testCase, R, [250 1]);
verifySize(testCase, E, [250 1]);
verifyEqual(testCase, diagnostics.corrected_rows, 32);
verifyEqual(testCase, E(indices), zeros(32, 1));
verifyEqual(testCase, R, D.R_occupied_vec);
unchanged = true(250, 1);
unchanged(indices) = false;
verifyEqual(testCase, E(unchanged), D.E_occupied_vec(unchanged));
F = griddedInterpolant({linspace(0, 70000, 5), linspace(0, 100000, 5), ...
    linspace(2500, 200000, 10)}, reshape(E, [5 5 10]), 'linear', 'none');
verifyEqual(testCase, F(D.R_free_vec, D.E_free_vec, D.L_vec), E, 'AbsTol', 1e-10);
end

function testRoundoffBoundariesAndMixedVectorOrientations(testCase)
D = example_grid();
D.R_occupied_vec(1) = -1e-12;
D.E_occupied_vec(1) = 1e-12;
D.E_occupied_vec(2) = -1e-9; % Within absolute tolerance.
D.E_occupied_vec(7) = 1e-12; % Positive occupancy away from a zero-pool boundary stays positive.
D.success_flag = D.success_flag';
D.R_free_vec = D.R_free_vec';
[R, E] = build_grid(D);
verifyEqual(testCase, R(1), 0);
verifyEqual(testCase, E(1:2), [0; 0]);
verifyEqual(testCase, E(7), 1e-12);
verifySize(testCase, R, [250 1]);
end

function testInvalidRowsStillRejected(testCase)
for failure = 1:5
    D = example_grid();
    switch failure
        case 1, D.success_flag(2) = 0;
        case 2, D.E_occupied_vec(2) = NaN;
        case 3, D.E_occupied_vec(2) = -1e-3;
        case 4, D.R_occupied_vec(2) = Inf;
        case 5, D.E_occupied_vec(2) = 1; % E_free=0 cannot support bound E.
    end
    verifyError(testCase, @() build_grid(D), ...
        'GeneLengthBuildInterpolation:IncompleteGrid');
end
end

function D = example_grid()
[r, e, L] = ndgrid(linspace(0, 70000, 5), linspace(0, 100000, 5), linspace(2500, 200000, 10));
D.R_free_vec = r(:);
D.E_free_vec = e(:);
D.L_vec = L(:);
D.R_occupied_vec = r(:).*L(:)/1e5;
D.E_occupied_vec = D.R_occupied_vec.*e(:)/1e5;
D.success_flag = ones(250, 1);
D.error_messages = repmat({''}, 250, 1);
end

function [R, E, diagnostics] = build_grid(D)
% Exercise validation through the public builder, including saved interpolants.
output_root = tempname;
mkdir(output_root);
old_root = getenv('CPAD_RESULTS_ROOT');
old_visibility = get(groot, 'defaultFigureVisible');
old_figures = findall(groot, 'Type', 'figure');
cleanup = onCleanup(@() restore_environment(output_root, old_root, old_visibility, old_figures)); %#ok<NASGU>
setenv('CPAD_RESULTS_ROOT', output_root);
set(groot, 'defaultFigureVisible', 'off');
grid_dir = cpad_analysis_output_dir('GeneLengthAnalysis');
results.data = D;
results.metadata.model_variant = 'full_kinetics_rapid_EH_disassembly';
results.metadata.pool_mode = 'fixed_free_pools';
results.metadata.grid_layout = 'ndgrid';
results.metadata.success_rate = 100;
results.parameters.base_parameters = default_parameters();
results.parameters.base_parameters.EBindingNumber = 5;
results.grid.R_free_values = linspace(0, 70000, 5);
results.grid.E_free_values = linspace(0, 100000, 5);
results.grid.L_values = linspace(2500, 200000, 10);
save(fullfile(grid_dir, 'full_gene_length_grid_data_test.mat'), 'results');
I = run_builder();
R = I.functions.R_occupied_interp(D.R_free_vec(:), D.E_free_vec(:), D.L_vec(:));
E = I.functions.E_occupied_interp(D.R_free_vec(:), D.E_free_vec(:), D.L_vec(:));
diagnostics = I.validation.grid_roundoff;
end

function I = run_builder()
evalc('GeneLengthBuildInterpolation');
I = interpolation_results;
end

function restore_environment(output_root, old_root, old_visibility, old_figures)
setenv('CPAD_RESULTS_ROOT', old_root);
set(groot, 'defaultFigureVisible', old_visibility);
close(setdiff(findall(groot, 'Type', 'figure'), old_figures));
rmdir(output_root, 's');
end
