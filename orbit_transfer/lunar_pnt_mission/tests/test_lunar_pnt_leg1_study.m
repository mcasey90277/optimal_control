function [ok, skipped] = test_lunar_pnt_leg1_study(resultsFile)
%% Purpose:
%
%   Verifies the SAVED results of lunar_pnt_leg1_study (it does not re-run
%   the solve): (1) structure: S, P1 [1 x 6], T0, knobs, run_meta (git
%   head + dirty flag, package identity with SHA-256, droMd5, droPhase0,
%   regime), failed; PRODUCTION results: both test-seam lists empty (and
%   no final-T0 seam) and run_meta.verdict exactly 'PASS'; (2) six accepted transfers, or a non-empty failed list
%   whose every entry carries a reason; (3) every accepted transfer's gates
%   pass (G.pass), its hand-off target is (k-1)/6 with |errPhase| below the
%   gate tolerance, its departure is the one GTO's state at its tDep
%   (recomputed here from S.orient0, S.gtoPhase0 and leg1_inputs), its
%   stored history starts at rv0 at tDep and ends at the flown arrival at
%   tArr, its thrust directions are unit vectors equal to -lambda_v/
%   |lambda_v| of the stored flight and re-fly it (stored alphaMissKm <
%   the study's alphaTolKm); (4) T0 >= max(tArr), and T0 is a slot
%   instant: mod(delta + T0/tauD, 1/6) = 0; (5) HEADLINE PINS measured on
%   the first run: T0, the six flight times, total stage propellant, total
%   stack wet mass, largest jet power, minimum located Earth and Moon
%   altitudes, achieved arrival maxima, number refuted (see the numbers
%   below and the reproducibility note); (6) FORCED FINAL-GATE FAILURE
%   (deliberate, through the seamFinalT0 injection knob; this section runs
%   the study, about 16 s: doSecond false, no figures, into a temporary
%   folder): with seams set and no outDir the study refuses to write into
%   results/ (lunar_pnt_leg1_study:seam); with T0 offset by 0.01 ND for
%   the final gates every satellite fails 'handoff at final T0', none is
%   accepted, failed lists six reasons, the verdict is not PASS, the
%   summary does not claim six closed, no sidecar is exported, and the
%   checks of (1) refuse the file (verdict not PASS). Skipped when the
%   function is given a resultsFile.
%   Returns skipped = true (and ok = true) when the results file does not
%   exist.
%
%% Inputs:
%
%  resultsFile              char (optional)         results to check
%                                                   [results/lunar_pnt_leg1_
%                                                   study.mat] (seam for the
%                                                   deliberate-failure runs on
%                                                   tampered copies)
%
%% Outputs:
%
%  ok                       logical                 true if every assertion
%                                                   held (else it throws)
%  skipped                  logical                 true if no results file
%
%% Revision History:
%  M. Casey                                                   (c) 11/11/2026
%  M. Casey  fix round 1: seams empty, verdict PASS, alpha check,
%            miss pins 1e-4 km / 1e-6 m/s                        (c) 11/11/2026
%  Copyright Coorbital Inc.
%% ------------------------ Begin Code Sequence ---------------------------

here = fileparts(mfilename('fullpath'));
root = fileparts(here);
run(fullfile(root, 'setup_lunar_pnt_paths.m'));
if nargin < 1 || isempty(resultsFile)
    resultsFile = fullfile(root, 'results', 'lunar_pnt_leg1_study.mat');
end
f = resultsFile;
skipped = false;
if ~isfile(f)
    fprintf('  SKIPPED: no results file (%s); run lunar_pnt_leg1_study first\n', f);
    ok = true;  skipped = true;
    return
end
R = load(f);

% ---- 1. structure ----
need = {'S', 'P1', 'T0', 'knobs', 'run_meta', 'failed'};
miss = need(~isfield(R, need));
assert(isempty(miss), 'results lack %s', strjoin(miss, ', '));
assert(numel(R.P1) == 6, 'P1 has %d entries', numel(R.P1));
rm = R.run_meta;
assert(isfield(rm, 'optimal_control') && isfield(rm, 'optimal_controlDirty') && ...
       numel(rm.identity.sha256) == 64 && numel(rm.droMd5) == 32 && ...
       isequal(rm.droPhase0, (0:5)/6) && abs(rm.tauD - 1) < 1e-9, 'run_meta incomplete');
kn = R.knobs;
assert(isempty(kn.seamRefuse) && isempty(kn.seamRefute) && ...
       (~isfield(kn, 'seamFinalT0') || isequal(kn.seamFinalT0, 0)), ...
       'production results were written with test seams set');
assert(strcmp(rm.verdict, 'PASS'), 'saved verdict is "%s", not PASS', rm.verdict);
pf = {'z8', 'tf', 'tDep', 'tArr', 'rv0', 'rvf', 'G', 't', 'X', 'alpha', 'route', 'attempts', ...
      'accepted', 'k', 'slot', 'iD', 'iA', 'iS', 'n', 'droPhaseArr'};
miss = pf(~isfield(R.P1, pf));
assert(isempty(miss), 'P1 lacks %s', strjoin(miss, ', '));

% ---- 2. six accepted, or failures with reasons ----
acc = [R.P1.accepted];
if ~all(acc)
    assert(~isempty(R.failed) && all(arrayfun(@(x) ~isempty(x.reason), R.failed)), ...
           'satellites not accepted but failed lists no reasons');
    fprintf('  NOT ALL ACCEPTED: failed %s\n', strjoin(arrayfun(@(x) sprintf('sat %d (%s)', x.k, x.reason), ...
                                                     R.failed, 'UniformOutput', false), '; '));
else
    assert(isempty(R.failed), 'all accepted but failed is not empty');
end

% ---- 3. per accepted transfer ----
L = leg1_inputs(struct('catalogMat', rm.identity.file));
tol = R.P1(find(acc, 1)).G.opts.handoffTolPhase;
for k = find(acc)
    q = R.P1(k);
    assert(q.k == k && isequal(q.G.pass, true), 'sat %d: k %d, G.pass %d', k, q.k, q.G.pass);
    assert(abs(q.G.handoff.target - (k - 1)/6) < 1e-15 && abs(q.G.handoff.errPhase) < tol, ...
           'sat %d: hand-off target %.6f err %.3g', k, q.G.handoff.target, q.G.handoff.errPhase);
    o = R.S.orient0 + L.sgn*(180/pi)*q.tDep;
    s = mod(R.S.gtoPhase0 + q.tDep/L.tauG, 1);
    dG = max(abs(L.gtoState(o, s) - q.rv0(:)));
    assert(dG < 1e-9, 'sat %d: departure is %.3g ND off the one GTO', k, dG);
    assert(abs(q.t(1) - q.tDep) < 1e-15 && abs(q.t(end) - q.tArr) < 1e-9 && ...
           max(abs(q.X(1:6, 1) - q.rv0(:))) < 1e-14, 'sat %d: stored history endpoints', k);
    un = abs(sqrt(sum(q.alpha.^2, 1)) - 1);
    assert(max(un) < 1e-12, 'sat %d: thrust directions not unit (%.3g)', k, max(un));
    assert(q.alphaMissKm < kn.alphaTolKm, 'sat %d: stored alpha re-flies to %.3g km (tol %g)', k, ...
           q.alphaMissKm, kn.alphaTolKm);
    fprintf('  sat %d: cell (%d,%d) n %d, orient %+.3f deg, hand-off err %+.2e, Earth %.1f km, Moon %.1f km, %s\n', ...
            k, q.iD, q.iA, q.n, q.orientDeg, q.G.handoff.errPhase, q.G.earth.hMinKm, q.G.moon.hMinKm, ...
            q.G.second.status);
end

% ---- 4. T0 ----
tArr = [R.P1(acc).tArr];
assert(R.T0 >= max(tArr) - 1e-12, 'T0 %.12g before the last arrival %.12g', R.T0, max(tArr));
ph = mod(R.S.delta + R.T0/L.tauD, 1/6);
assert(min(ph, 1/6 - ph) < 1e-9, 'T0 is not a slot instant (%.3g)', ph);
dayNd = rm.regime.tStar/86400;
fprintf('  T0 %.6f d (last arrival %.6f d)\n', R.T0*dayNd, max(tArr)*dayNd);

% ---- 5. headline pins ----
% Measured 2026-11-11 (package SHA-256 e57d9236..., one 0 deg sheet), two
% runs in one session: every number below was BIT-IDENTICAL between them
% (z8, tf, T0, masses, clearances, misses: difference 0). The tolerances
% only allow solver / integrator noise across sessions or MATLAB updates:
% the polish root stops at |g| < 1e-10 ND, so tf to 1e-6 d; masses follow
% tf; clearances are located minima (fzero, 1e-15) to 0.01 km; the
% arrival misses to 1e-4 km and 1e-6 m/s.
pin = struct('T0Days', 2.647279358370, ...
             'tfDays', [2.434151938 2.343373057 2.386896438 2.571656464 2.589840788 2.618790470], ...
             'propKg', 760.417943043, 'wetKg', 1774.480634499, 'maxJetKW', 86.837324822, ...
             'minEarthKm', 349.863000, 'minMoonKm', 841.342778, ...
             'maxMissKm', 0.0473016831, 'maxMissMps', 0.000663426206, 'nRefuted', 0);
assert(all(acc), 'the pins are for a run with six accepted transfers');
G = [R.P1.G];  ms = [G.mass];  ea = [G.earth];  mo = [G.moon];  ar = [G.arrival];
got = struct('T0Days', R.T0*dayNd, 'tfDays', [R.P1.tf]*dayNd, 'propKg', sum([ms.stagePropKg]), ...
             'wetKg', sum([ms.wetKg]), 'maxJetKW', max([ms.jetKW]), 'minEarthKm', min([ea.hMinKm]), ...
             'minMoonKm', min([mo.hMinKm]), 'maxMissKm', max([ar.missKm]), 'maxMissMps', max([ar.missMps]), ...
             'nRefuted', nnz([G.refuted]));
tolP = struct('T0Days', 1e-9, 'tfDays', 1e-6, 'propKg', 1e-3, 'wetKg', 1e-3, 'maxJetKW', 1e-4, ...
              'minEarthKm', 0.01, 'minMoonKm', 0.01, 'maxMissKm', 1e-4, 'maxMissMps', 1e-6, 'nRefuted', 0);
for f = fieldnames(pin).'
    e = max(abs(got.(f{1}) - pin.(f{1})));
    assert(e <= tolP.(f{1}), 'pin %s: measured %s, pinned %s (tol %g)', f{1}, mat2str(got.(f{1}), 10), ...
           mat2str(pin.(f{1}), 10), tolP.(f{1}));
    fprintf('  pin %-11s %s  (|diff| %.2g, tol %g)\n', f{1}, mat2str(got.(f{1}), 10), e, tolP.(f{1}));
end

% ---- 6. forced final-gate failure (seamFinalT0) ----
if nargin < 1 || strcmp(resultsFile, fullfile(root, 'results', 'lunar_pnt_leg1_study.mat'))
    study = fullfile(root, 'lunar_pnt_leg1_study.m');
    id = '';
    try
        exec_study(study, struct('headless', true, 'seamFinalT0', 0.01));
    catch err
        id = err.identifier;
    end
    assert(strcmp(id, 'lunar_pnt_leg1_study:seam'), 'seams with default outDir: expected a refusal, got "%s"', id);
    fprintf('  seams set, outDir default: refused (%s)\n', id);
    scratch = '/private/tmp/claude-501/-Users-msc-Desktop-optimal-control/70650963-d331-4c2e-9b00-6156c2bd5033/scratchpad/leg1';
    if ~isfolder(scratch), mkdir(scratch); end
    tmp = tempname(scratch);
    mkdir(tmp);
    cleanTmp = onCleanup(@() rmdir(tmp, 's'));
    t0 = tic;
    txt = exec_study(study, struct('headless', true, 'seamFinalT0', 0.01, 'doSecond', false, ...
                                   'makeFigures', false, 'outDir', tmp));
    Ff = load(fullfile(tmp, 'lunar_pnt_leg1_study.mat'), 'P1', 'failed', 'run_meta');
    reasons = {Ff.failed.reason};
    assert(~any([Ff.P1.accepted]) && numel(Ff.failed) == 6 && ...
           all(contains(reasons, 'handoff at final T0')), 'forced: acceptance did not follow the final gates');
    assert(~strcmp(Ff.run_meta.verdict, 'PASS') && contains(Ff.run_meta.verdict, 'TEST SEAMS'), ...
           'forced: verdict "%s"', Ff.run_meta.verdict);
    assert(~contains(txt, '6 of 6 closed') && contains(txt, '0 of 6 closed'), 'forced: summary claims closure');
    assert(~isfile(fullfile(tmp, 'leg1_transfers_gto_dro_5N.mat')) && contains(txt, 'CAPTURE: skipped (test seams'), ...
           'forced: a sidecar was exported from a seam run');
    id = '';
    try
        evalc('test_lunar_pnt_leg1_study(fullfile(tmp, ''lunar_pnt_leg1_study.mat''))');
    catch err
        id = err.message;
    end
    assert(contains(id, 'test seams') || contains(id, 'not PASS'), 'forced results were not refused by the checks');
    fprintf('  forced (T0 + 0.01 ND at the final gates, %.0f s): 0 of 6 accepted, reasons "%s", verdict "%s", no sidecar; checks refuse it: %s\n', ...
            toc(t0), reasons{1}, Ff.run_meta.verdict, id);
end
fprintf('  PASS\n');
ok = true;
end

% ------------------------------------------------------------------------
function txt = exec_study(file, leg1_opts)
% EXEC_STUDY  Run the study script in this function's own workspace (its
%   clearvars cannot reach the caller) with leg1_opts set (the script
%   reads it from this workspace); output captured.
%   INPUTS: file char; leg1_opts struct.  OUTPUTS: txt char.
assert(isfile(file) && isstruct(leg1_opts), 'exec_study: bad arguments');
txt = evalc('run(file)');
end
