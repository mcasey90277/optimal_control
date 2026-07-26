function test_refine_loop_smoke()
% TEST_REFINE_LOOP_SMOKE  2-round loop on the fast 1.12x file runs & records.
%
% NOTE: runs up to 2 eps=0 re-solves of the 10-switch 1.12x problem (~min each).
% INPUTS:  none
% OUTPUTS: none (prints ALL PASS or asserts)
here = fileparts(mfilename('fullpath'));
% One campaign setup: direct/setup_paths adds lib, certify, viz, tests,
% cr3bp_common AND the indirect ms_band verifier. Before the 2026-07-26 flatten
% this cd'd into ms_band to borrow its setup_paths -- via '../../ms_band', a
% path that stopped resolving when ms_band moved to indirect/ on 2026-07-21, so
% this test had been failing at startup for five days.
addpath(fullfile(here, '..'));  setup_paths();

seed = fullfile(here, '..', 'results', 'minfuel', 'legacy_ms_f1120.mat');
opts = struct('maxRounds', 2, 'tag', 'smoke_1p12', 'K', 6, 'maxAdd', 30);
history = refine_loop(seed, opts);

assert(numel(history) >= 2, 'expect >= 2 measured rounds, got %d', numel(history));
assert(isnan(history(1).maxSwitchMove), 'round 1 has no previous move');
assert(all([history.nNodes] == sort([history.nNodes])), 'nodes non-decreasing');
assert(history(end).maxDefect < 1e-6 || ~history(end).converged, 'tight or flagged');
rcfgT = minfuel_config();   % refine_loop writes to the results tree, not next to this test
assert(isfile(fullfile(rcfgT.dirs.root, 'refine_smoke_1p12.png')), ...
    'figure not written to %s', rcfgT.dirs.root);
fprintf('ALL PASS (rounds=%d, final switches=%d, maxMove=%.2e)\n', ...
        numel(history), history(end).switches, history(end).maxSwitchMove);
end
