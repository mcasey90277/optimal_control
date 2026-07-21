function test_pmp_refine_indicator()
% TEST_PMP_REFINE_INDICATOR  Indicator sanity vs the known-good 1.12x file.
%
% INPUTS:  none
% OUTPUTS: none (prints ALL PASS or asserts)
here = fileparts(mfilename('fullpath'));
addpath(here);
addpath(fullfile(here, '..', '..', 'ms_band'));   run_setup_paths(here);
seed = fullfile(here, '..', 'results', 'minfuel', 'legacy_ms_f1120.mat');
assert(isfile(seed), 'missing seed file %s', seed);

opts = struct('M', 40, 'epsEval', 1e-4, 'mode', 'd', 'nbr', 3);
[score, tauSwitch, diag] = pmp_refine_indicator(seed, opts);

nN = numel(diag.Snode);  N = nN - 1;
assert(isequal(size(score), [1 N]), 'score must be 1xN');
assert(all(score >= 0), 'score must be nonnegative');
assert(any(score > 0), 'some intervals must score > 0 (switches present)');
assert(~isempty(tauSwitch) && issorted(tauSwitch), 'tauSwitch nonempty & sorted');
assert(numel(tauSwitch) >= 8 && numel(tauSwitch) <= 14, ...
       'expect ~10-12 direct switches at 1.12x, got %d', numel(tauSwitch));
assert(diag.betaSpread < 5, 'beta spread should be small (mode-d), got %.2f%%', diag.betaSpread);
assert(isequal(size(diag.Hres), [1 nN]), 'Hres must be 1xnN');
assert(diag.nViol >= 0, 'nViol counts');
fprintf('ALL PASS (nSwitch=%d, betaSpread=%.2f%%, HresMax=%.2e, nViol=%d)\n', ...
        numel(tauSwitch), diag.betaSpread, diag.HresMax, diag.nViol);
end

function run_setup_paths(here)
% add ms_band + its dependency chain (sundman_minfuel, lowThrust, pumpkyn)
old = cd(fullfile(here, '..', '..', 'ms_band'));  c = onCleanup(@() cd(old));
setup_paths();
end
