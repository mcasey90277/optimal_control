function ok = test_run_phase_torus_p0()
% TEST_RUN_PHASE_TORUS_P0  The P0 repairs to the torus driver (review of
% 2026-09-17, FINDINGS 78), one check per behaviour, through the driver's
% local-function seam `run_phase_torus('localfunctions')`:
%
%   1. commonOpts carries the campaign's departure phase and its first
%      anchor (file + phase) to the sheet stage; physicsOpts asks for the
%      closures only (no anchor to polish) for the filler and discovery.
%   2. spawnArc clears stale verdict files, and its shell wrapper writes a
%      .fail when MATLAB exits without writing a verdict.
%   3. campaignEnded lets a 'budget' campaign resume under a larger bound.
%   4. roundExtras offers this round's own holes file to the packager.
%   5. adoptLiveJobs adopts a recorded job that is alive and retires a dead
%      one, so a resume never spawns a second writer.
%   6. promotionVerdict promotes only a NEWLY registered root that beats
%      the spine and lies on no known family (both paths, one rule).
%
% INPUTS:  none
% OUTPUTS: ok [logical]  every check passed (prints PASS/FAIL per check)

ok = true;
here = fileparts(fileparts(mfilename('fullpath')));            % DRO_tulip/indirect
addpath(here, fullfile(fileparts(fileparts(here)), 'costate_common'));
H = run_phase_torus('localfunctions');
tmp = fullfile(tempdir, sprintf('torus_p0_%s', char(java.util.UUID.randomUUID())));
mkdir(tmp);  cleanup = onCleanup(@() rmdir(tmp, 's'));
lg = @(varargin) [];                                            % a silent log

% ---- 1. the operating point reaches the sheet, the closures the filler ----
in = struct('sD', [0.2 0.5], 'sA', [0.1 0.6], 'tag', 't', ...
            'engine', struct('thrustN', 0.07, 'ispS', 900, 'm0kg', 150), ...
            'orbits', struct('tauDRO', 1, 'NpTulip', 7, 'pmTulip', -1), ...
            'librarySeeds', false, 'rib', struct('wallSec', 900));
anchors = {'a', '/x/anchor_a.mat', 0.61, 'A'; 'b', '/x/anchor_b.mat', 0.12, 'B'};
c = H.commonOpts(in, anchors, '/x/arcs', 'arrival_arc_t_*.mat');
ok = chk(ok, isfield(c, 'anchorMat') && strcmp(c.anchorMat, '/x/anchor_a.mat'), 'commonOpts: the FIRST anchor''s file goes to the sheet');
ok = chk(ok, isfield(c, 'anchorSA') && c.anchorSA == 0.61, 'commonOpts: ... and its arrival phase');
ok = chk(ok, isequal(c.sD, [0.2 0.5]), 'commonOpts: the departure list travels whole (sD(1) is the spine)');
so = H.physicsOpts(in.engine, in.orbits, 0.2);
ok = chk(ok, isfield(so, 'physicsOnly') && so.physicsOnly && ~isfield(so, 'anchorMat') && so.sD == 0.2, ...
         'physicsOpts: closures only, no anchor, at the spine''s departure phase');

% ---- 2. spawnArc: stale verdicts cleared; a verdict is guaranteed --------
jobDir = fullfile(tmp, 'jobs');  mkdir(jobDir);
outMat = fullfile(tmp, 'arrival_arc_t_a_dn_long.mat');
fclose(fopen([outMat '.fail'], 'w'));  fclose(fopen([outMat '.done'], 'w'));
arc = struct('nStep', 10, 'deadlineSec', 60, 'span', 1.15);
[pid, jobFile] = H.spawnArc('a', '/x/anchor_a.mat', 0.61, -1, outMat, in.engine, in.orbits, 0.2, in.sA, ...
                            arc, jobDir, '/usr/bin/false', tmp, here, lg);
ok = chk(ok, isfinite(pid) && isfile(jobFile), 'spawnArc: launched, job file written');
ok = chk(ok, ~isfile([outMat '.done']), 'spawnArc: a stale .done is removed before launch');
t0 = tic;  while ~isfile([outMat '.fail']) && toc(t0) < 15, pause(0.2); end
ok = chk(ok, isfile([outMat '.fail']) && contains(fileread([outMat '.fail']), 'exited'), ...
         'spawnArc: MATLAB exiting without a verdict yields a .fail that says so');
ok = chk(ok, isfile([outMat '.pid']), 'spawnArc: the wrapper records the MATLAB pid beside the arc');

% ---- 3. a budget campaign resumes under a larger bound -------------------
st = struct('status', 'budget', 'roundsDone', 2, 'reason', 'r');
[e1, ~] = H.campaignEnded(st, 6);   [e2, ~] = H.campaignEnded(st, 2);
[e3, ~] = H.campaignEnded(setfield(st, 'status', 'done'), 6); %#ok<SFLD>
[e4, ~] = H.campaignEnded(setfield(st, 'status', 'running'), 6); %#ok<SFLD>
ok = chk(ok, ~e1 && e2 && e3 && ~e4, 'campaignEnded: budget resumes when maxRounds > roundsDone; done stays done');

% ---- 4. this round's holes file is offered to the packager ---------------
mkdir(fullfile(tmp, 'round_01'));  mkdir(fullfile(tmp, 'round_02'));
fclose(fopen(fullfile(tmp, 'round_01', 'fine_rib_col01.mat'), 'w'));
fclose(fopen(fullfile(tmp, 'round_02', 'fine_rib_direct_holes.mat'), 'w'));
x = H.roundExtras(tmp, 2);
ok = chk(ok, any(endsWith(x, fullfile('round_01', 'fine_rib_col01.mat'))) && ...
             any(endsWith(x, fullfile('round_02', 'fine_rib_direct_holes.mat'))), ...
         'roundExtras: earlier ribs AND this round''s existing holes file');

% ---- 5. resume adopts live jobs and retires dead ones --------------------
need = {'a', '/x/a.mat', 0.61, -1, fullfile(tmp, 'arc_a.mat'); 'b', '/x/b.mat', 0.12, +1, fullfile(tmp, 'arc_b.mat')};
st = struct('jobs', struct('pid', {feature('getpid'), 999999}, 'file', {need{1,5}, need{2,5}}));
[adopted, st] = H.adoptLiveJobs(need, st, lg);
ok = chk(ok, isequal(adopted(:).', [true false]), 'adoptLiveJobs: a live recorded job is adopted, a dead one is not');
ok = chk(ok, numel(st.jobs) == 1 && st.jobs(1).pid == feature('getpid'), 'adoptLiveJobs: the dead job leaves the record');

% ---- 6. one promotion rule for both paths --------------------------------
ok = chk(ok, ~H.promotionVerdict(false, 17, 20, 0.05, 0), 'promotionVerdict: a root already registered is never re-anchored');
ok = chk(ok, ~H.promotionVerdict(true, 17, 20, 0.05, 2), 'promotionVerdict: a root on a known family is not anchored');
ok = chk(ok, ~H.promotionVerdict(true, 19.98, 20, 0.05, 0), 'promotionVerdict: a gain below acceptDays is not enough');
ok = chk(ok,  H.promotionVerdict(true, 17, 20, 0.05, 0), 'promotionVerdict: new, faster, unattached -> anchor');
ok = chk(ok,  H.promotionVerdict(true, 17, NaN, 0.05, -1), 'promotionVerdict: an empty column takes any new unattached root');

% ---- 7. 'packaged' is not enough: every stage that ran must have passed ---
good = struct('state', 'packaged', 'blockers', {{}}, 'stages', struct('package', true, 'audit', true, 'sweep', NaN));
badA = struct('state', 'packaged', 'blockers', {{'audit: 1 of 2 entries bad'}}, 'stages', struct('package', true, 'audit', false, 'sweep', true));
ok = chk(ok, ~throws(@() H.assertPackaged(good, 'round_01', 'test')), 'assertPackaged: all stages ok (one not run) passes');
ok = chk(ok,  throws(@() H.assertPackaged(badA, 'round_01', 'test')), 'assertPackaged: a catalog whose audit failed is refused');
ok = chk(ok,  throws(@() H.assertPackaged(setfield(good, 'state', 'pending'), 'round_01', 'test')), 'assertPackaged: a call that did not package is refused'); %#ok<SFLD>

% ---- 8. the registry knows a root by its COSTATES, not its flight time ------
% Two roots at one arrival phase whose flight times are 1e-4 d apart are one
% root only if their costates agree; the old rule (t_f within 1e-3 d) merged them.
reg = fullfile(tmp, 'registry.mat');
zA = [1; 2; 3; 4; 5; 6; 7; 4.0];  zB = [7; 6; 5; 4; 3; 2; 1; 4.0];
root = @(z, tf) struct('sA', 0.25, 'tfDays', tf, 'z', z);
ok = chk(ok,  H.registerRoot(reg, root(zA, 17.7050), 0, 'first'), 'a root is registered');
ok = chk(ok, ~H.registerRoot(reg, root(zA*(1 + 1e-9), 17.7050 + 2e-7), 0, 'a re-polish'), 'the same costates again: a duplicate');
ok = chk(ok,  H.registerRoot(reg, root(zB, 17.7051), 0, 'another root'), 'other costates at nearly the same flight time: a DISTINCT root, registered');
L = load(reg);  ok = chk(ok, numel(L.direct) == 2, 'the registry holds two roots');
ok = chk(ok,  H.sameRoot(zA, zA*(1 + 1e-9)) && H.sameRoot(zA*(1 + 1e-9), zA), 'sameRoot is symmetric');
ok = chk(ok, ~H.sameRoot(zA, [zA(1:7); 4.0 + 1e-3]), 'the same costates with another flight time are not the same root');
ok = chk(ok, ~H.sameRoot(zA, [NaN; zA(2:8)]) && ~H.sameRoot(zA, [zeros(7, 1); 4.0]), 'a non-finite or all-zero costate vector is never "the same root"');
ok = chk(ok, throws(@() H.registerRoot(reg, root([NaN; zA(2:8)], 17.7), 0, 'bad')), 'a root with non-finite costates is refused, not registered');
ok = chk(ok,  H.registerRoot(reg, root(zA, 17.7050), 0.25, 'another spine'), 'the same costates at ANOTHER departure phase are a distinct record');

% ---- 9. the wait for a round scales with the round -------------------------------
% (48 h was sized for 24 x 24; a 48 x 48 round of ribs + audit + sweep is about
% 55 h of work, and the driver used to give up on it)
h = @(nD, nA, nW, user) H.roundDeadlineSec(nD, nA, nW, user)/3600;
ok = chk(ok, h(24, 24, 4, []) == 48, 'a 24 x 24 round keeps the 48 h wait');
ok = chk(ok, h(48, 48, 4, []) > 100 && h(24, 96, 4, []) > 100, sprintf('a 48 x 48 round waits %.0f h, 24 x 96 waits %.0f h', h(48, 48, 4, []), h(24, 96, 4, [])));
ok = chk(ok, h(48, 48, 8, []) < h(48, 48, 4, []), 'more rib workers, a shorter wait');
ok = chk(ok, h(48, 48, 4, 200) == 200, 'an explicit .roundDeadlineHours is honoured');

if ok, fprintf('test_run_phase_torus_p0: ALL PASS\n'); else, fprintf('test_run_phase_torus_p0: FAIL\n'); end
end

function ok = chk(ok, c, msg)
% CHK  Print one PASS/FAIL line and fold it into ok.  INPUTS: ok; c; msg.
% OUTPUTS: ok.
if c, fprintf('  PASS  %s\n', msg); else, fprintf('  FAIL  %s\n', msg); end
ok = ok && c;
end

function tf = throws(f)
% THROWS  Does calling f throw?  INPUTS: f (handle).  OUTPUTS: tf.
tf = false;
try, f(); catch, tf = true; end
end
