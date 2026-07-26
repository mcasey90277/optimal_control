% TEST_INSERTION_STATES  Pin the six endpoint states against silent drift.
%
% insertion_states is the single source of truth for where every low-thrust
% transfer in this repo STARTS and ENDS. If it moves, every certified result
% moves with it -- silently, because a transfer to slightly-wrong endpoints
% still converges and still looks like a valid trajectory. Nothing downstream
% would flag it.
%
% The numbers below were measured on 2026-07-26 from the pre-move file (then in
% GTO_tulip/direct/sundman_minfuel) and are reproduced here to full double
% precision. They are the gate on the move into cr3bp_common: the move deleted
% two stale addpath calls, and this test is the evidence that deleting them
% changed no endpoint.
%
% Two independent checks, because neither alone is enough:
%   1. PINNED LITERALS (all six target/criterion pairs) -- catches any change,
%      including in the alternate criteria no seed exists for. A regression pin,
%      not a derivation: it says the endpoints are what they have always been,
%      not that they are physically right.
%   2. AGAINST THE CERTIFIED SEEDS -- the default criteria must reproduce what
%      the banked energy backbones actually hold. This is the check with real
%      authority: it ties the helper to the .mat files the campaigns solved
%      from, so the pin above cannot drift into self-consistent nonsense.
%
% A deliberate endpoint change SHOULD fail this test; update the pinned values
% in the same commit, and re-run the affected campaign.
%
% MERGED 2026-07-26 from GTO_tulip/direct/sundman_minfuel/test_insertion_states.m,
% which followed the function here. That file had been BROKEN since the
% 2026-07-21 reorg -- it loaded '../elfo/results/energy_elfo_freetf.mat', a path
% that stopped existing when ELFO moved to GTO_ELFO/, so it failed at load time
% and had not run for five days. Its seed comparison (check 2) is preserved
% here with correct paths.
%
% INPUTS:  none
% OUTPUTS: none (throws on any change)
%
% REFERENCES:
%   [1] cr3bp_common/insertion_states.m (the function under test).
%   [2] orbit_transfer/CODE_STRUCTURE.md (Tier-1: the move and its rationale).
here = fileparts(mfilename('fullpath'));
addpath(fullfile(here, '..'));
setup_cr3bp_common();          % gto_{tulip,elfo}_endpoints need pumpkyn

% target, criterion, expected label, expected rvf (pre-move, full precision)
CASES = {
  'tulip','campaign','tulipCampaign', [ 1.00658107295709,   0.0425745746906059,  -0.0557780910480905, ...
                                       -0.16004281347248,   0.0665702939657711,  -0.260455693516549  ]
  'tulip','maxydot', 'tulipMaxYdot',  [ 0.996283777921214,  7.32921748798876e-05, 0.013002103551144,  ...
                                       -0.0015994690531156, 1.15453423647664,    -0.00269486593903987]
  'tulip','apoapsis','tulipApoapsis', [ 1.04379468739791,   0,                   -0.0798645605909944, ...
                                        0,                  0.1334177677339,      0                  ]
  'elfo', 'nearest', 'elfoNearest',   [ 0.961765053808788, -0.0124457238162329,  -0.0319836796962375, ...
                                        0.0659987089188916,-0.230731759329445,   -0.31227600841443   ]
  'elfo', 'apolune', 'elfoApolune',   [ 0.972899317200594, -0.0250289679097755,  -0.0434389570168043, ...
                                        0.204972499987968, -0.107458244318226,   -0.00860189457954956]
  'elfo', 'perilune','elfoPerilune',  [ 0.989695156350785,  0.0048622726871064,   0.00781632222450811, ...
                                       -1.39977615440419,   0.466722684996873,    0.0396184262575154 ]
};

% The GTO departure is shared by every pipeline -- one value, checked once.
RV0 = [0.00349629072294633, -0.0072962582600817, 0, ...
       4.19147893803368,     8.98865558978329,   0];

fprintf('\n=== insertion_states endpoint gate ===\n');

%% ------------------------------------------------------------------------
%% Check 1: pinned literals, all six target/criterion pairs
%% ------------------------------------------------------------------------

% Tolerance is not zero: gto_{tulip,elfo}_endpoints run pumpkyn propagations,
% so the last bits can move with a library or platform change. It is tight
% enough that any real endpoint change fails loudly.
TOL = 1e-12;

for k = 1:size(CASES,1)
    tgt = CASES{k,1};  crit = CASES{k,2};  wantLabel = CASES{k,3};  wantRvf = CASES{k,4};

    % No warnings expected: the stale addpath that warned on every ELFO call
    % pre-move is gone, and this asserts it stays gone.
    lastwarn('');
    [rv0, rvf, meta] = insertion_states(tgt, crit);
    [wmsg, wid] = lastwarn;

    assert(isequal(size(rvf), [1 6]), 'rvf must be 1x6 for %s/%s', tgt, crit);
    assert(strcmp(meta.label, wantLabel), ...
        'label changed for %s/%s: expected %s, got %s', tgt, crit, wantLabel, meta.label);

    dRv0 = max(abs(rv0(:).' - RV0));
    assert(dRv0 <= TOL, ...
        'GTO departure rv0 CHANGED (max |d| = %.3g) on %s/%s -- every pipeline starts here.', ...
        dRv0, tgt, crit);

    d = max(abs(rvf - wantRvf));
    assert(d <= TOL, ...
        ['insertion state CHANGED for %s/%s (max |d| = %.3g > %.3g).\n' ...
         '  got : %s\n  want: %s\n' ...
         'If this change is deliberate, update the pinned value here in the same ' ...
         'commit and re-run the affected campaign -- certified results depend on it.'], ...
        tgt, crit, d, TOL, mat2str(rvf, 15), mat2str(wantRvf, 15));

    assert(isempty(wmsg), ...
        ['insertion_states(%s,%s) raised a warning: [%s] %s\n' ...
         'A stale addpath warned on every ELFO call before the cr3bp_common move; ' ...
         'this asserts no path warning has crept back in.'], tgt, crit, wid, wmsg);

    fprintf('  %-6s %-9s %-15s max|d| = %.2e   (no warnings)\n', tgt, crit, meta.label, d);
end

% Defaults must map to the criteria the certified seeds were built with.
[~, rvfTd, mT] = insertion_states('tulip');
[~, rvfEd, mE] = insertion_states('elfo');
assert(strcmp(mT.criterion, 'campaign'), 'tulip default criterion changed: %s', mT.criterion);
assert(strcmp(mE.criterion, 'nearest'),  'elfo default criterion changed: %s',  mE.criterion);
fprintf('  defaults: tulip->campaign, elfo->nearest\n');

%% ------------------------------------------------------------------------
%% Check 2: against the certified seeds
%% ------------------------------------------------------------------------
% The authority behind the pinned literals. If the helper and the banked
% backbones ever disagree, every drift guard downstream fires -- so catch it
% here, where the message can say why.
ot     = fullfile(here, '..', '..');                     % orbit_transfer/
% Derive the tulip seed location from minfuel_config rather than hardcoding it.
% This line named .../direct/sundman_minfuel/results/... and broke silently when
% the 2026-07-26 flatten moved that tree to direct/results -- exactly the class of
% breakage this file exists to catch, in this file. cfg.dirs follows the move.
cfgS   = minfuel_config();
seedT  = fullfile(cfgS.dirs.energy, cfgS.fname('energy', 1.12));
seedE  = fullfile(ot, 'GTO_ELFO',  'direct', 'elfo', 'results', 'energy_elfo_freetf.mat');

for s = {seedT, 'tulip'; seedE, 'elfo'}.'
    assert(isfile(s{1}), ...
        ['seed not found: %s\nThis test compares endpoints against the certified ' ...
         'backbones; a moved seed tree must be fixed here, not skipped.'], s{1});
end

E  = load(seedT);
Ee = load(seedE);

dRv0 = norm(RV0 - E.rv0(:).');
assert(dRv0 < 1e-12, 'rv0 disagrees with the tulip backbone seed (|d| = %.3g)', dRv0);

dT = norm(rvfTd - E.rvf(:).');
assert(dT < 1e-12, ...
    ['tulip ''campaign'' endpoint disagrees with the backbone seed %s (|d| = %.3g).\n' ...
     'The helper and the solved-from data have diverged -- every downstream drift ' ...
     'guard will fire. Fix the disagreement, do not relax this tolerance.'], seedT, dT);

dE = norm(rvfEd - Ee.rvf(:).');
assert(dE < 1e-12, ...
    ['elfo ''nearest'' endpoint disagrees with the ELFO seed %s (|d| = %.3g).'], seedE, dE);

fprintf('  vs certified seeds: tulip |d| = %.2e, elfo |d| = %.2e\n', dT, dE);

% Unknown inputs must fail loudly rather than silently returning a default.
for bad = {{'banana','campaign'}, {'tulip','banana'}, {'elfo','banana'}}
    try
        insertion_states(bad{1}{1}, bad{1}{2});
        error('test_insertion_states:noThrow', ...
            'insertion_states(%s,%s) should have thrown', bad{1}{1}, bad{1}{2});
    catch ME
        assert(startsWith(ME.identifier, 'insertion_states:'), ...
            'wrong error id for %s/%s: %s', bad{1}{1}, bad{1}{2}, ME.identifier);
    end
end
fprintf('  unknown target/criterion reject loudly\n');

fprintf('\ntest_insertion_states: ALL PASS (%d endpoints pinned)\n', size(CASES,1));
