function tests = test_pas_cleavage_at_pas
% Legacy equilibrium-model regression tests for cleavage at the PAS node.
tests = functiontests(localfunctions);
end

function setupOnce(testCase)
testCase.TestData.originalPath = path;
archive_root = fileparts(fileparts(mfilename('fullpath')));
addpath(archive_root);
end

function teardownOnce(testCase)
path(testCase.TestData.originalPath);
end

function testPasStateHasCleavageLoss(testCase)
P.N = 3;
P.PAS = 2;
P.k_in = 0;
P.k_e = 0;
P.k_e2 = 0;
P.kHon = 0;
P.kHoff = 0;
P.kc = 3;
P.Pol_total = 0;
P.Pol_free = 0;

X = zeros(P.N + (P.N - P.PAS + 1), 1);
X(P.N + 1) = 2;  % REH(1), the PAS node
dxdt = ode_dynamics_multipleE(X, P);

verifyEqual(testCase, dxdt(P.N + 1), -6, 'AbsTol', 1e-12);
end

function testPasCleavageIsAtFirstBinRightEdge(testCase)
R_sol = zeros(2, 1);
REH_sol = [3; 0];
P_sim.kc = 2;
P_sim.k_e = 0;
P_sim.k_e2 = 0;
P_sim.L_a = 100;

[exit_cdf, distances_bp, CAD] = calculate_pas_cleavage_profile( ...
    R_sol, REH_sol, P_sim, 'PercentCleavage', 50);

verifyEqual(testCase, distances_bp(1), 100);
verifyEqual(testCase, exit_cdf(1), 1, 'AbsTol', 1e-12);
verifyEqual(testCase, CAD, 100);
end
