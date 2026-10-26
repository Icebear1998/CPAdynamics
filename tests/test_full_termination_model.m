function tests = test_full_termination_model
% Full finite-rate model: Python parity, reaction bookkeeping and integration.
tests = functiontests(localfunctions);
end

function setupOnce(testCase)
testCase.TestData.originalPath = path;
testCase.TestData.root = fileparts(fileparts(mfilename('fullpath')));
addpath(testCase.TestData.root);
addpath(fullfile(testCase.TestData.root, 'tests'));
end

function teardownOnce(testCase)
path(testCase.TestData.originalPath);
end

function testStateMaps(testCase)
for M = 1:7
    P = default_parameters();
    P.geneLength_bp = 100;
    P.PASposition = 100;
    model = build_full_rate_matrices(P, M);
    S = model.states;
    verifySize(testCase, S.R, [(M+1)*(M+2)/2, 2]);
    verifySize(testCase, S.RHE, [M*(M+1)/2, 2]);
    verifyTrue(testCase, all(S.RHE(:, 2) >= 1));
end
end

function testPythonOperatorParity(testCase)
reference = load(fullfile(testCase.TestData.root, 'tests', 'fixtures', 'full_model_python_reference.mat'));
for M = [1, 5]
    [model, ~] = build_full_rate_matrices(reference.parameters, M);
    prefix = sprintf('M%d_', M);
    verifyEqual(testCase, model.states.R, double(reference.([prefix 'R_states'])));
    verifyEqual(testCase, model.states.RHE, double(reference.([prefix 'RHE_states'])));
    verifyEqual(testCase, full(model.A0), full(reference.([prefix 'A0'])), 'AbsTol', 1e-12);
    verifyEqual(testCase, full(model.AE), full(reference.([prefix 'AE'])), 'AbsTol', 1e-15);
    X = reference.([prefix 'state']);
    verifyEqual(testCase, ode_dynamics_full_multipleE(0, X, model), ...
        reference.([prefix 'rhs']), 'AbsTol', 1e-9);
    verifyEqual(testCase, full(full_model_jacobian(0, X, model)), ...
        full(reference.([prefix 'jacobian'])), 'AbsTol', 1e-12);
end
end

function testReferenceParameterSweepMatchesPython(testCase)
P = python_reference_parameters();
P.kEoff_engaged = 0.05;
expected = [1698.670, 1127.246, 975.185, 926.894, 912.346, 908.631];
for M = 1:6
    [R, RHE, P_out, D] = run_full_termination_simulation(P, M);
    [~, ~, cad] = calculate_full_pas_cleavage_profile(R, RHE, P_out);
    verifyEqual(testCase, cad, expected(M), 'AbsTol', 0.002);
    verifyLessThan(testCase, abs(D.E_conservation_residual), 1e-6);
    verifyLessThan(testCase, abs(D.Pol_conservation_residual), 1e-6);
    verifyLessThan(testCase, D.rhs_max_abs, 1e-7);
    verifyGreaterThanOrEqual(testCase, min(D.state_ss), -1e-10);
    verifyEqual(testCase, P_out.kHon, P.kHon);
end
end

function testDisabledAndNormalRateMatchPython(testCase)
P = python_reference_parameters();
for row = [1, 5; 1535.653580763, 845.209693107; 3169.511, 1581.043]
    M = row(1);
    for rateIndex = 1:2
        P.kEoff_engaged = (rateIndex-1)*0.5;
        [R, RHE, P_out] = run_full_termination_simulation(P, M);
        [~, ~, cad] = calculate_full_pas_cleavage_profile(R, RHE, P_out);
        verifyEqual(testCase, cad, row(rateIndex+1), 'AbsTol', 0.002);
    end
end
end

function testZeroEngagedOffSweepConservesPools(testCase)
% Regression: M=1 and M=2 were rejected with the looser fzero tolerance.
P = python_reference_parameters();
P.kEoff_engaged = 0;
expected = [1535.653580763, 1032.912272151, 898.830651653, ...
            857.555743681, 845.209693107, 842.071228735];
for M = 1:6
    [R, RHE, P_out, D] = run_full_termination_simulation(P, M);
    [~, ~, cad] = calculate_full_pas_cleavage_profile(R, RHE, P_out);
    verifyEqual(testCase, cad, expected(M), 'AbsTol', 0.002);
    verifyLessThan(testCase, abs(D.E_conservation_residual), 1e-6);
    verifyLessThan(testCase, abs(D.Pol_conservation_residual), 1e-6);
    verifyLessThan(testCase, D.rhs_max_abs, 1e-6);
    verifyLessThan(testCase, abs(D.flux_residual), 1e-6);
end
end

function testIndependentHAndEDissociation(testCase)
P = isolated_parameters();
P.kHoff = 0.2;
P.kEoff_engaged = 0.05;
[model, ~] = build_full_rate_matrices(P, 1);
X = zeros(model.size+2, 1);
X(model.r_size+1) = 1; % RHE(p=1,e=1)
dX = ode_dynamics_full_multipleE(0, X, model);
verifyEqual(testCase, dX(model.states.R_index(2, 2)), 0.2, 'AbsTol', 1e-14);
verifyEqual(testCase, dX(model.states.R_index(2, 1)), 0.05, 'AbsTol', 1e-14);
verifyEqual(testCase, dX(model.r_size+1), -0.25, 'AbsTol', 1e-14);
verifyEqual(testCase, dX(end-1), 0.05, 'AbsTol', 1e-14);
verifyEqual(testCase, dX(end), 0, 'AbsTol', 1e-14);
end

function testEngagedLossHasNoMultiplicityAndPreservesP(testCase)
P = isolated_parameters();
P.kEoff_engaged = 0.07;
[model, ~] = build_full_rate_matrices(P, 5);
for h = 1:model.sh
    pe = model.states.RHE(h, :);
    X = zeros(model.size+2, 1);
    X(model.r_size+h) = 3;
    expected = zeros(size(X));
    expected(model.r_size+h) = -0.21;
    expected(model.states.R_index(pe(1)+1, pe(2))) = 0.21;
    expected(end-1) = 0.21;
    verifyEqual(testCase, ode_dynamics_full_multipleE(0, X, model), expected, 'AbsTol', 1e-13);
end
end

function testFreeSiteLossScalesWithFreeSites(testCase)
P = isolated_parameters();
P.kPoff = 0.3;
[model, ~] = build_full_rate_matrices(P, 5);
for pe = [4, 1; 2, 2; 5, 0]'
    X = zeros(model.size+2, 1);
    source = model.states.R_index(pe(1)+1, pe(2)+1);
    X(source) = 2;
    expected = zeros(size(X));
    rate = 2*(pe(1)-pe(2))*P.kPoff;
    if rate > 0
        expected(source) = -rate;
        expected(model.states.R_index(pe(1), pe(2)+1)) = rate;
    end
    verifyEqual(testCase, ode_dynamics_full_multipleE(0, X, model), expected, 'AbsTol', 1e-13);
end
end

function testFirstBinCleavageAndCoordinates(testCase)
P = isolated_parameters();
P.kc = 2;
[model, ~] = build_full_rate_matrices(P, 1);
X = zeros(model.size+2, 1);
X(model.r_size+1) = 3;
dX = ode_dynamics_full_multipleE(0, X, model);
verifyEqual(testCase, dX(model.r_size+1), -6, 'AbsTol', 1e-14);
verifyEqual(testCase, dX(end-1:end), [6; 6], 'AbsTol', 1e-14);
[cdf, distances, cad] = calculate_full_pas_cleavage_profile([0; 0], [3; 0], P);
verifyEqual(testCase, cdf, [1; 1], 'AbsTol', 1e-14);
verifyEqual(testCase, distances, [100; 200]);
verifyEqual(testCase, cad, 50, 'AbsTol', 1e-14);
end

function testCommitmentHoldsResourcesUntilRelease(testCase)
P = isolated_parameters();
P.kc = 2;
P.k_rel = 0.5;
[model, ~] = build_full_rate_matrices(P, 5);
for h = 1:model.sh
    e = model.states.RHE(h, 2);
    X = zeros(model.size+2, 1);
    X(model.r_size+h) = 3;
    expected = zeros(size(X));
    expected(model.r_size+h) = -6;
    expected(model.c_offset+h) = 6; % same (p,e); nothing released yet
    verifyEqual(testCase, ode_dynamics_full_multipleE(0, X, model), expected, 'AbsTol', 1e-13);
    X(:) = 0;
    X(model.c_offset+h) = 3;
    expected(:) = 0;
    expected(model.c_offset+h) = -1.5;
    expected(end-1:end) = [1.5*e; 1.5];
    verifyEqual(testCase, ode_dynamics_full_multipleE(0, X, model), expected, 'AbsTol', 1e-13);
end
end

function testCommittedExchangeFreezesP(testCase)
P = isolated_parameters();
P.k_rel = 1;
P.kEoff = 0.3;
P.kEon = 0.01;
P.kPoff = 5;
P.kPon_min = 5;
[model, ~] = build_full_rate_matrices(P, 3);
idx = @(p, e) model.c_offset+model.states.RHE_index(p+1, e+1);
X = zeros(model.size+2, 1);
X(idx(3, 2)) = 1;
X(end-1) = 10; % E_free
dX = ode_dynamics_full_multipleE(0, X, model);
verifyEqual(testCase, dX(idx(3, 1)), 0.3, 'AbsTol', 1e-13); % one unengaged E leaves
verifyEqual(testCase, dX(idx(3, 3)), 0.1, 'AbsTol', 1e-13); % one free site binds E
verifyEqual(testCase, dX(idx(3, 2)), -1.4, 'AbsTol', 1e-13);
verifyEqual(testCase, nnz(abs(dX(1:model.c_offset)) > 0), 0); % no p change, no R/RHE
end

function testFastReleaseRecoversCommitmentRelease(testCase)
P = default_parameters();
for M = [1, 5]
    P.k_rel = Inf;
    [R0, RHE0, P0] = run_full_termination_simulation(P, M);
    [~, ~, cad0] = calculate_full_pas_cleavage_profile(R0, RHE0, P0);
    P.k_rel = 1e6;
    [R1, RHE1, P1, D1] = run_full_termination_simulation(P, M);
    [~, ~, cad1] = calculate_full_pas_cleavage_profile(R1, RHE1, P1);
    verifyEqual(testCase, cad1, cad0, 'RelTol', 1e-4);
    verifyEqual(testCase, P1.Ef_ss, P0.Ef_ss, 'RelTol', 1e-4);
    verifyLessThan(testCase, abs(D1.flux_residual), 1e-6);
end
end

function testCommittedPoolCountsAsDemand(testCase)
P = default_parameters();
P.k_rel = 0.01;
[~, RHE, P_out, D] = run_full_termination_simulation(P, 5);
verifyEqual(testCase, D.C_total, P.kc*sum(RHE)/P.k_rel, 'RelTol', 1e-8);
verifyEqual(testCase, D.Pol_bound, sum(D.R_micro(:))+sum(RHE)+D.C_total, 'RelTol', 1e-12);
verifyEqual(testCase, P_out.Pol_free_ss+D.Pol_bound, P.Pol_total, 'AbsTol', 1e-6);
verifyEqual(testCase, P_out.Ef_ss+D.E_bound, P.E_total, 'AbsTol', 1e-6);
verifyGreaterThan(testCase, D.C_E_bound, D.C_total); % engaged E plus exchanged E
end

function testZeroEHasNoRecognition(testCase)
P = default_parameters();
P.E_total = 0;
P.kEoff_engaged = 0.05;
[~, RHE, P_out] = run_full_termination_simulation(P, 5);
verifyEqual(testCase, RHE, zeros(size(RHE)), 'AbsTol', 1e-10);
verifyEqual(testCase, P_out.Ef_ss, 0);
end

function testUnreachedCleavageIsNotExtrapolated(testCase)
P = isolated_parameters();
P.kc = 1;
P.k_e = 1;
verifyWarning(testCase, @() calculate_full_pas_cleavage_profile([0; 9], [1; 0], P), ...
    'calculate_full_pas_cleavage_profile:ThresholdNotReached');
warningState = warning('off', 'calculate_full_pas_cleavage_profile:ThresholdNotReached');
cleanup = onCleanup(@() warning(warningState)); %#ok<NASGU>
[cdf, ~, cad, diagnostics] = calculate_full_pas_cleavage_profile([0; 9], [1; 0], P);
verifyEqual(testCase, cdf, [0.1; 0.1], 'AbsTol', 1e-14);
verifyTrue(testCase, isnan(cad));
verifyFalse(testCase, diagnostics.threshold_reached);
end

function testJacobianAndConservation(testCase)
P = default_parameters();
P.geneLength_bp = 500;
P.PASposition = 300;
P.kEoff_engaged = 0.05;
[model, ~] = build_full_rate_matrices(P, 3);
x = 0.5+0.1*sin((1:model.size)');
X = [x; P.E_total-model.e_counts'*x; P.Pol_total-sum(x)];
v = cos((1:numel(X))');
h = 0.001;
dX = ode_dynamics_full_multipleE(0, X, model);
finiteDifference = (ode_dynamics_full_multipleE(0, X+h*v, model) - ...
                    ode_dynamics_full_multipleE(0, X-h*v, model))/(2*h);
verifyEqual(testCase, full_model_jacobian(0, X, model)*v, finiteDifference, 'AbsTol', 1e-7);
verifyEqual(testCase, model.e_counts'*dX(1:model.size)+dX(end-1), 0, 'AbsTol', 1e-8);
verifyEqual(testCase, sum(dX(1:model.size))+dX(end), 0, 'AbsTol', 1e-8);
end

function testTimeIntegrationMatchesSteadyState(testCase)
P = default_parameters();
P.geneLength_bp = 1000;
P.PASposition = 600;
P.Pol_total = 100;
P.E_total = 200;
P.kEon = 0.001;
P.kEoff_engaged = 0.05;
for M = [1, 5]
    [~, ~, ~, D] = run_full_termination_simulation(P, M);
    model = D.model;
    X0 = [zeros(model.size, 1); P.E_total; P.Pol_total];
    options = odeset('RelTol', 1e-8, 'AbsTol', 1e-10, ...
        'Jacobian', @(t, X) full_model_jacobian(t, X, model));
    % Long horizon: the committed pool relaxes on the 1/k_rel ~ 100 s scale.
    [~, trajectory] = ode15s(@(t, X) ode_dynamics_full_multipleE(t, X, model), [0, 3000], X0, options);
    verifyEqual(testCase, trajectory(end, :)', D.state_ss, 'AbsTol', 1e-5);
end
end

function P = python_reference_parameters()
% Python CAD fixtures predate the revised appendix defaults. Keep their
% original inputs explicit so future default changes do not invalidate them.
P = default_parameters();
P.k_e = 65/P.L_a;
P.kEon = 2.5e-6;
P.kHon = 4;
P.kHoff = 2;
P.kc = 0.13;
P.kPon_min = 0.01;
P.kPon_slope = 0.005;
P.kPoff = 1;
P.Pol_total = 70000; % default_parameters later moved to 100000
P.k_rel = Inf; % fixtures release resources at commitment
end

function P = isolated_parameters()
P = default_parameters();
P.geneLength_bp = 100;
P.PASposition = 100;
zeroFields = {'k_in', 'k_e', 'k_e2', 'kEon', 'kEoff', 'kHon', 'kHoff', ...
              'kc', 'kPon_min', 'kPon_slope', 'kPoff', 'kEoff_engaged'};
for i = 1:numel(zeroFields)
    P.(zeroFields{i}) = 0;
end
P.k_rel = Inf;
end
