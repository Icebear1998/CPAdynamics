%% SANITYCHECKFULLMODEL  Run and report the full finite-rate model checks.
%
% Companion to SanityCheckMultipleE.m (which covers the equilibrium model).
% Runs tests/test_full_termination_model.m and prints a readable PASS/FAIL
% report, then reports the headline numbers for the current default_parameters.
%
% Usage (from the repository root):
%   SanityCheckFullModel
%
% Exit behaviour: the script errors at the end if any test failed, so it can
% also be used from a batch run:
%   matlab -batch "SanityCheckFullModel"

repoRoot = fileparts(mfilename('fullpath'));
testFile = fullfile(repoRoot, 'tests', 'test_full_termination_model.m');
if exist(testFile, 'file') ~= 2
    error('SanityCheckFullModel:MissingTests', 'Cannot find %s', testFile);
end

fprintf('=== Full finite-rate R/RHE model: sanity checks ===\n');
fprintf('MATLAB %s\n\n', version);

fixture = fullfile(repoRoot, 'tests', 'fixtures', 'full_model_python_reference.mat');
if exist(fixture, 'file') ~= 2
    fprintf(['NOTE: the Python cross-check fixture is missing, so those tests\n' ...
             '      will be skipped. Regenerate it with:\n' ...
             '        cd python && python export_matlab_full_reference.py\n\n']);
end

startTime = tic;
results = runtests(testFile);
elapsed = toc(startTime);

fprintf('\n--- Individual checks ---\n');
for k = 1:numel(results)
    shortName = regexprep(results(k).Name, '^.*/', '');
    if results(k).Passed
        verdict = 'PASS';
    elseif results(k).Incomplete
        verdict = 'SKIP';
    else
        verdict = 'FAIL';
    end
    fprintf('  [%s] %-55s %6.2f s\n', verdict, shortName, results(k).Duration);
end

nPassed = sum([results.Passed]);
nFailed = sum([results.Failed]);
nSkipped = sum([results.Incomplete]);
fprintf('\n--- Summary ---\n');
fprintf('  %d passed, %d failed, %d skipped (%.1f s total)\n', ...
    nPassed, nFailed, nSkipped, elapsed);

%% Headline numbers for the parameters this repository currently ships.
fprintf('\n--- Current default_parameters() with the full model ---\n');
P = default_parameters();
P.kEoff_engaged = 0.05;
fprintf('  kHoff = %.4g, kEoff = %.4g, kEoff_engaged = %.4g\n', ...
    P.kHoff, P.kEoff, P.kEoff_engaged);
fprintf('   %-3s %10s %10s %12s %12s %12s\n', ...
    'M', 'CAD50 bp', 'max CDF', 'E_free', 'Pol_free', 'max|dX/dt|');
for M = 1:5
    [R_sol, RHE_sol, P_out, details] = run_full_termination_simulation(P, M);
    warnState = warning('off', 'calculate_full_pas_cleavage_profile:ThresholdNotReached');
    [~, ~, CAD, diagnostics] = calculate_full_pas_cleavage_profile(R_sol, RHE_sol, P_out);
    warning(warnState);
    fprintf('   %-3d %10.2f %10.5f %12.4g %12.4g %12.2e\n', ...
        M, CAD, diagnostics.max_exit_cdf, P_out.Ef_ss, P_out.Pol_free_ss, ...
        details.rhs_max_abs);
end

if nFailed > 0
    error('SanityCheckFullModel:Failed', '%d full-model sanity check(s) failed.', nFailed);
end
fprintf('\nAll full-model sanity checks passed.\n');
