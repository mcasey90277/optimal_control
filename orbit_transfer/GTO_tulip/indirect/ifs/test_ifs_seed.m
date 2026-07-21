function test_ifs_seed()
% TEST_IFS_SEED  Full-problem and window seeds build a bounded-residual Z.
%
% INPUTS:  none
% OUTPUTS: none (prints ALL PASS or asserts)
here = fileparts(mfilename('fullpath'));  addpath(here);  setup_paths();
seed = fullfile(here,'..','sundman_minfuel','results','minfuel','legacy_ms_f1120.mat');
assert(isfile(seed), 'missing %s', seed);

% window mode (rung 1): exactly one switch inside, fixedState terminal
[Zw, probW, metaW] = ifs_seed(seed, struct('mode','window','winSwitch',5,'winPad',60));
assert(probW.k == 1, 'window must isolate ONE switch, got k=%d', probW.k);
assert(strcmp(probW.termMode,'fixedState'), 'window is fixedState');
Rw = ifs_residual(Zw, probW);
assert(numel(Rw) == 8+17, 'window R square');
assert(metaW.seedResNorm < 1e-1, 'window seed residual should be small, got %.2e', metaW.seedResNorm);

% full mode (rung 2/3): rendezvous terminal, many switches
[Zf, probF, metaF] = ifs_seed(seed, struct('mode','full'));
assert(probF.k >= 8, 'full 1.12x should have ~10 switches, got %d', probF.k);
assert(strcmp(probF.termMode,'rendezvous'), 'full is rendezvous');
Rf = ifs_residual(Zf, probF);
assert(numel(Rf) == 8 + 17*probF.k, 'full R square');
fprintf('ALL PASS (window k=1 seedRes=%.2e; full k=%d seedRes=%.2e)\n', ...
        metaW.seedResNorm, probF.k, norm(Rf));
end
