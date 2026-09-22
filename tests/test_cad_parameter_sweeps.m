function tests = test_cad_parameter_sweeps
% Exercise the public script, including its local helpers and saved output.
% Run with: runtests('tests/test_cad_parameter_sweeps.m')
tests = functiontests(localfunctions);
end

function setupOnce(testCase)
testCase.TestData.originalPath = path;
project_root = fileparts(fileparts(mfilename('fullpath')));
addpath(project_root);
testCase.TestData.warningState = warning('off', ...
    'calculate_full_pas_cleavage_profile:ThresholdNotReached');
testCase.TestData.originalVisibility = get(groot, 'DefaultFigureVisible');
testCase.TestData.originalFigures = findall(groot, 'Type', 'figure');
set(groot, 'DefaultFigureVisible', 'off');
testCase.TestData.originalResultsRoot = getenv('CPAD_RESULTS_ROOT');
testCase.TestData.originalForceSave = getenv('CPAD_FORCE_SAVE');
testCase.TestData.outputRoot = tempname;
mkdir(testCase.TestData.outputRoot);
setenv('CPAD_RESULTS_ROOT', testCase.TestData.outputRoot);
setenv('CPAD_FORCE_SAVE', 'true');

% Run the actual entry point so tests do not depend on private functions.
run(fullfile(project_root, 'Sweep2DCad.m'));
testCase.TestData.P = P;
testCase.TestData.sweeps = sweeps;
testCase.TestData.nPoints = nPoints;
testCase.TestData.percentCleavage = percent_cleavage;
testCase.TestData.bindingNumber = EBindingNumber;
end

function teardownOnce(testCase)
figures = findall(groot, 'Type', 'figure');
close(setdiff(figures, testCase.TestData.originalFigures));
set(groot, 'DefaultFigureVisible', testCase.TestData.originalVisibility);
setenv('CPAD_RESULTS_ROOT', testCase.TestData.originalResultsRoot);
setenv('CPAD_FORCE_SAVE', testCase.TestData.originalForceSave);
rmdir(testCase.TestData.outputRoot, 's');
path(testCase.TestData.originalPath);
warning(testCase.TestData.warningState);
end

function testPairsMatchIndependentFullModelSolves(testCase)
P = testCase.TestData.P;
sweeps = testCase.TestData.sweeps;
verifyEqual(testCase, {sweeps.x_name}, {'kHd', 'kc', 'kc'});
verifyEqual(testCase, {sweeps.y_name}, {'kEd', 'kHd', 'kEd'});
for k = 1:3
    data = sweeps(k);
    verifySize(testCase, data.CAD_matrix, repmat(testCase.TestData.nPoints, 1, 2));
    verifyEqual(testCase, data.x_values([1 end]), P.ranges.(data.x_name), 'RelTol', 1e-12);
    verifyEqual(testCase, data.y_values([1 end]), P.ranges.(data.y_name), 'RelTol', 1e-12);
    % Opposite corners check orientation and that unswept rates stay fixed.
    for cell_index = [1 numel(data.x_values); numel(data.y_values) 1]'
        row = cell_index(1);
        col = cell_index(2);
        expected = P;
        switch k
            case 1
                expected.kHoff = data.x_values(col)*P.kHon;
                expected.kEoff = data.y_values(row)*P.kEon;
            case 2
                expected.kc = data.x_values(col);
                expected.kHoff = data.y_values(row)*P.kHon;
            case 3
                expected.kc = data.x_values(col);
                expected.kEoff = data.y_values(row)*P.kEon;
        end
        [R, RHE, P_sim, details] = run_full_termination_simulation( ...
            expected, testCase.TestData.bindingNumber);
        [~, ~, cad, diagnostics] = calculate_full_pas_cleavage_profile( ...
            R, RHE, P_sim, 'PercentCleavage', testCase.TestData.percentCleavage);
        verifyEmpty(testCase, data.error_matrix{row, col});
        if isnan(cad)
            verifyTrue(testCase, isnan(data.CAD_matrix(row, col)));
        else
            verifyEqual(testCase, data.CAD_matrix(row, col), cad, 'AbsTol', 1e-8);
        end
        verifyEqual(testCase, data.max_exit_cdf_matrix(row, col), ...
            diagnostics.max_exit_cdf, 'AbsTol', 1e-12);
        verifyEqual(testCase, data.rhs_residual_matrix(row, col), ...
            details.rhs_max_abs, 'AbsTol', 1e-10);
    end
end
verifyEqual(testCase, [sweeps.x_base], [P.kHoff/P.kHon, P.kc, P.kc]);
verifyEqual(testCase, [sweeps.y_base], [P.kEoff/P.kEon, P.kHoff/P.kHon, P.kEoff/P.kEon]);
end

function testPlotAndSaveEveryPair(testCase)
for data = testCase.TestData.sweeps
    figure_name = sprintf('CAD: %s vs %s, M = %d', ...
        data.x_name, data.y_name, data.EBindingNumber);
    fig = findall(groot, 'Type', 'figure', 'Name', figure_name);
    fig = setdiff(fig, testCase.TestData.originalFigures);
    assertNumElements(testCase, fig, 1);
    ax = findobj(fig, 'Type', 'axes');
    verifyEqual(testCase, ax.XScale, 'log');
    verifyEqual(testCase, ax.YScale, 'log');
    verifyEqual(testCase, ax.CLim, [0 1600]);
    plot_legend = findobj(fig, 'Type', 'legend');
    assertNumElements(testCase, plot_legend, 1);
    verifyEqual(testCase, plot_legend.Location, 'northwest');
    verifyTrue(testCase, startsWith(plot_legend.String{1}, ...
        sprintf('Base parameters CAD_{%g} = ', data.percent_cleavage)));
    contours = findobj(ax, 'Type', 'contour', 'Fill', 'on');
    if any(isfinite(data.CAD_matrix(:)))
        assertNumElements(testCase, contours, 1);
        expected_display = data.CAD_matrix;
        expected_display(expected_display > 1600) = 1600;
        verifyEqual(testCase, contours.ZData, expected_display);
    end

    analysis_type = sprintf('sweep_2D_%s_%s_CAD', data.x_name, data.y_name);
    output_dir = fullfile(testCase.TestData.outputRoot, analysis_type);
    file_prefix = sprintf('sweep2D_%s_%s_CAD_EBinding%d_', ...
        data.x_name, data.y_name, data.EBindingNumber);
    files = dir(fullfile(output_dir, [file_prefix '*_data.txt']));
    assertNumElements(testCase, files, 1);
    content = fileread(fullfile(files(1).folder, files(1).name));
    verifyTrue(testCase, contains(content, sprintf('rows = %s, columns = %s', ...
        data.y_name, data.x_name)));
    verifyTrue(testCase, contains(content, 'full_kinetics_rapid_EH_disassembly'));
    verifyTrue(testCase, contains(content, 'max_exit_cdf_matrix'));
    verifyTrue(testCase, contains(content, 'rhs_residual_matrix'));
    verifyNumElements(testCase, dir(fullfile(output_dir, [file_prefix '*.png'])), 1);
end
end
