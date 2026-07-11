function test_refine_loop_smoke()
% TEST_REFINE_LOOP_SMOKE  2-round loop on the fast 1.12x file runs & records.
%
% NOTE: runs up to 2 eps=0 re-solves of the 10-switch 1.12x problem (~min each).
% INPUTS:  none
% OUTPUTS: none (prints ALL PASS or asserts)
here = fileparts(mfilename('fullpath'));  addpath(here);
addpath(fullfile(here, '..', '..', 'ms_band'));   % ms_band/setup_paths does NOT add ms_band itself
old = cd(fullfile(here, '..', '..', 'ms_band'));  c = onCleanup(@() cd(old));
setup_paths();  cd(old);

seed = fullfile(here, '..', 'results', 'minfuel', 'legacy_ms_f1120.mat');
opts = struct('maxRounds', 2, 'tag', 'smoke_1p12', 'K', 6, 'maxAdd', 30);
history = refine_loop(seed, opts);

assert(numel(history) >= 2, 'expect >= 2 measured rounds, got %d', numel(history));
assert(isnan(history(1).maxSwitchMove), 'round 1 has no previous move');
assert(all([history.nNodes] == sort([history.nNodes])), 'nodes non-decreasing');
assert(history(end).maxDefect < 1e-6 || ~history(end).converged, 'tight or flagged');
assert(isfile(fullfile(here, 'refine_smoke_1p12.png')), 'figure written');
fprintf('ALL PASS (rounds=%d, final switches=%d, maxMove=%.2e)\n', ...
        numel(history), history(end).switches, history(end).maxSwitchMove);
end
