%% LUNAR_PNT_LEG1_STUDY  Six 5 N stages from one inertial GTO onto the six DRO slots leg 2 starts from
%   Edit the knobs, press Run. One launch puts six stacks on one GTO fixed
%   in inertial space (350 x 35,786 km); in the Earth-Moon rotating frame
%   its orientation turns at the frame rate. Each stack leaves the GTO at a
%   catalog departure phase after whole GTO revolutions, flies a 5 N per
%   150 kg minimum-time transfer (engine always on) and arrives EXACTLY on
%   its DRO slot; the six slots are one sixth of a DRO period apart, and
%   the hand-off epoch T0 is the first instant after the last arrival at
%   which the satellites sit at the DRO phases (k-1)/6 that leg 2
%   (lunar_pnt_mission_study) starts from. Every step calls a tested
%   library function (lib/); nothing is reimplemented here.
%
%     0  knobs                              (one source; leg1_opts overrides)
%     1  inputs and identity                I1 leg1_inputs checks; package
%                                           snapshot in <outDir>/inputs
%     2  schedule                           C1 a complete assignment exists
%     3  polish each satellite: polish ->   per attempt: polish, flight, arrival,
%        gates -> exclude -> re-schedule    one GTO, Earth, Moon, mass (blocking);
%        that satellite only (same triple)  second order (prefer non-refuted)
%     4  gates table                        GA GF GO GH GE GM GW G2 GT GD, VERDICT
%     5  mass and power table + caveat
%     6  save                               results/lunar_pnt_leg1_study.mat
%     7  capture                            results/leg1_transfers_gto_dro_5N.mat
%     8  figures                            leg1_schedule, leg1_3d, leg1_mass,
%                                           mission_timeline (combined_timeline)
%
%   HOW A SATELLITE IS CLOSED. The scheduler (schedule_gto_dro) chooses a
%   global triple (orient0, gtoPhase0, delta) and, per slot, the earliest
%   admissible catalog cell. Slot by slot the chosen cell is polished
%   (polish_gto_dro: orientation walk, then the arrival root with the
%   departure pinned) and gated (leg1_gates). A polish failure or a failed
%   BLOCKING gate (flight admissibility, arrival, one GTO, located Earth
%   and Moon clearance, mass feasibility) adds the candidate [iD iA iS n] to
%   the exclusion list and that satellite is re-scheduled with the triple
%   LOCKED (opts.lockTriple), so the other five keep their transfers. A
%   REFUTED transfer (second order) is not blocking but the next candidate
%   is tried; when none within maxAttempts is non-refuted the earliest
%   refuted one is kept and flagged. After maxAttempts (or when the slot
%   has no candidate left) knobs.offgridSolver is called once when it is a
%   function handle E = f(L, orientDeg, sD, sA) (E.ok, E.z8, E.tf): its
%   cell is polished as a one-cell sheet at that orientation (route
%   'offgrid'); with no solver configured (the default) the satellite is
%   marked FAILED with its reason. Nothing is loosened and no unpolished
%   catalog entry is ever substituted. A REFUTED-but-otherwise-accepted
%   candidate kept by the fallback rule does NOT trigger the off-grid
%   solver: the solver is called only when no candidate passed the
%   blocking gates.
%
%   T0 AND LABELS. Polished arrival epochs differ from the table's, so the
%   hand-off epoch is recomputed from the POLISHED arrivals with the
%   schedule's slot grid: T0 = the first t >= the last accepted arrival
%   with mod(delta + t/tauD, 1/6) = 0; slot j's satellite is labelled k =
%   6*mod(delta + (j-1)/6 + T0/tauD, 1) + 1 (mod 6), the leg-2 spacecraft
%   number. The T0-dependent hand-off gate is evaluated only on this final
%   T0 (every gate re-run without the second order; the second-order block
%   does not depend on T0 and is carried over). A satellite is ACCEPTED
%   only if those final gates pass (G.pass, incl. the hand-off at the final
%   T0 and k) AND its stored thrust direction re-flies the flight below
%   alphaTolKm; otherwise it joins failed with the failing gates named. A
%   failed slot's stored G.handoff is blanked ('not evaluated at final T0').
%
%   STORED HISTORIES (leg-1 ND clock, t = 0 at the earliest permitted
%   departure): per satellite the flown history t [1 x M], X [7 x M]
%   (position, velocity, mass fraction) from departure to arrival, the
%   thrust direction alpha [3 x M] = -lambda_v/|lambda_v| (the tfMinProp
%   convention: the thrust acceleration is (Tnd/m)*alpha; verified in the
%   study by re-integrating the CR3BP with it, gate GD), plus orient0,
%   gtoPhase0, tDep (the GTO ride before departure: L.gtoState(orient0 +
%   sgn*(180/pi)*t, mod(gtoPhase0 + t/tauG, 1))) and droPhaseArr (DRO
%   phase at tArr; afterwards L.slotState(mod(droPhaseArr + (t -
%   tArr)/tauD, 1))).
%
%   leg1_opts: a struct in the workspace before running overrides any knob
%   (unknown names are refused). The script CLEARS its workspace
%   (clearvars -except leg1_opts). TEST SEAMS seamRefuse / seamRefute
%   (candidate keys 'iD_iA_iS_n') force a refusal or a refuted verdict on
%   a really polished and gated candidate, to exercise the exclusion,
%   re-schedule, refuted-fallback and failed branches; seamFinalT0 (ND)
%   offsets the T0 the final gates use (forces the final hand-off to fail).
%   Empty / 0 by default. When any is set the study refuses to write into
%   the default results/ folder (outDir must point elsewhere), skips the
%   sidecar export, and names the seams in the VERDICT.
%
%   PACKAGE SNAPSHOT. The package at the default path will be replaced by
%   the fleet's; section 1 keeps the one this run used: a file in
%   <outDir>/inputs whose own SHA-256 equals the package's IS the snapshot
%   (any name), else the package is copied there as
%   costate_catalog_gto_dro_<first 8 hex>.mat (never overwriting;
%   leg1_snapshot_package). run_meta.packageSnapshot and
%   .packageSnapshotStatus record it; leg1_package_for finds it by hash.
%
%   SCOPE. The second-order verdict and "minimum-time" refer to the
%   fixed-endpoint problem to the slot's arrival point (a local test at the
%   root); the moving-target rendezvous condition is not checked.
%
%  M. Casey                                                   (c) 11/11/2026
%  M. Casey  final fix wave: package snapshot, scope line     (c) 11/11/2026
%  Copyright Coorbital Inc.

clearvars -except leg1_opts
if ~exist('leg1_opts', 'var'), leg1_opts = struct(); end
here = fileparts(mfilename('fullpath'));

%% ========================================================================
%  0. KNOBS -- every adjustable parameter, with units; one source.
%% ========================================================================
catalogMat      = '';                % GTO -> DRO package ('' = leg1_inputs default)
leg2Mat         = '';                % leg-2 results ('' = results/lunar_pnt_mission_study.mat)
maxAttempts     = 6;                 % catalog candidates tried per satellite
offgridSolver   = [];                % E = f(L, orientDeg, sD, sA) fallback; [] = none
dryFrac         = 0.15;              % stage dry mass per kg of stage propellant
doSecond        = true;              % second-order verdicts (conjugate + hypothesis gates)
schedOpts       = struct();          % extra schedule_gto_dro options (grids, tolPhase, ...)
polishOpts      = struct();          % polish_gto_dro options
alphaTolKm      = 1;                 % GD: re-integration with alpha vs the flight [km]
makeFigures     = true;              % section 8
exportTransfers = true;              % section 7
headless        = false;             % figures invisible and closed after export
outDir          = '';                % results folder ('' = results/ beside this script)
seamRefuse      = {};                % TEST SEAM: candidate keys 'iD_iA_iS_n' treated as refused
seamRefute      = {};                % TEST SEAM: candidate keys treated as refuted (second order)
seamFinalT0     = 0;                 % TEST SEAM: offset [ND] added to the T0 the final gates use

knobs = struct('catalogMat', catalogMat, 'leg2Mat', leg2Mat, 'maxAttempts', maxAttempts, ...
               'offgridSolver', offgridSolver, 'dryFrac', dryFrac, 'doSecond', doSecond, ...
               'schedOpts', schedOpts, 'polishOpts', polishOpts, 'alphaTolKm', alphaTolKm, ...
               'makeFigures', makeFigures, 'exportTransfers', exportTransfers, ...
               'headless', headless, 'outDir', outDir, 'seamRefuse', {seamRefuse}, ...
               'seamRefute', {seamRefute}, 'seamFinalT0', seamFinalT0);
for f = fieldnames(leg1_opts).'
    if ~isfield(knobs, f{1})
        error('lunar_pnt_leg1_study:knob', 'unknown leg1_opts field ''%s''', f{1});
    end
end
for f = fieldnames(leg1_opts).'
    knobs.(f{1}) = leg1_opts.(f{1});
end
if isempty(knobs.outDir), knobs.outDir = fullfile(here, 'results'); end
if isempty(knobs.leg2Mat), knobs.leg2Mat = fullfile(here, 'results', 'lunar_pnt_mission_study.mat'); end
maxAttempts = knobs.maxAttempts;      offgridSolver = knobs.offgridSolver;
dryFrac     = knobs.dryFrac;          doSecond      = knobs.doSecond;
schedOpts   = knobs.schedOpts;        polishOpts    = knobs.polishOpts;
alphaTolKm  = knobs.alphaTolKm;       makeFigures   = knobs.makeFigures;
exportTransfers = knobs.exportTransfers;
headless    = knobs.headless;         outDir        = knobs.outDir;
seam = struct('refuse', {knobs.seamRefuse}, 'refute', {knobs.seamRefute}, 'finalT0', knobs.seamFinalT0);
seamUsed = ~isempty(seam.refuse) || ~isempty(seam.refute) || seam.finalT0 ~= 0;
if seamUsed && strcmp(canon(outDir), canon(fullfile(here, 'results')))
    error('lunar_pnt_leg1_study:seam', ['test seams are set (seamRefuse / seamRefute / seamFinalT0): refusing to ' ...
          'write into the production results folder %s; set leg1_opts.outDir elsewhere'], outDir);
end
if ~isfolder(outDir), mkdir(outDir); end
figVis = 'on';
if headless, figVis = 'off'; end

gates = cell(0, 4);                                    % {id, status, value, detail}
secName = {'1 inputs', '2 schedule', '3 polish', '4 gates', '5 mass', '6 save', '7 capture', '8 figures'};
secWall = nan(1, numel(secName));
tStudy  = tic;
gitBin = '/Library/Developer/CommandLineTools/usr/bin/git';
gitOC  = git_state(gitBin, here);
fprintf('lunar_pnt_leg1_study: optimal_control %s; MATLAB %s\n', gitOC.text, version);

%% ========================================================================
%  1. INPUTS AND IDENTITY
%% ========================================================================
tSec = tic;
fprintf('1. INPUTS\n');
run(fullfile(here, 'setup_lunar_pnt_paths.m'));
lOpts = struct();
if ~isempty(knobs.catalogMat), lOpts.catalogMat = knobs.catalogMat; end
L = leg1_inputs(lOpts);
dayNd = L.tStar/86400;
propKgPerNd = L.m0kg*L.phys.Tnd/L.phys.cnd;
fprintf('   package %s\n   SHA-256 %s, created %s, %d entries, sheet orientations [%s] deg, rung %d\n', ...
        L.identity.file, L.identity.sha256, L.identity.created, L.identity.n_entries, ...
        num2str(L.sheetOrient), L.identity.rung);
fprintf('   regime %g N, Isp %g s, m0 %g kg; GTO period %.6f ND (%.4f h), DRO period %.6f ND (%.4f d), sgn %+d\n', ...
        L.thrustN, L.ispS, L.m0kg, L.tauG, L.tauG*L.tStar/3600, L.tauD, L.tauD*dayNd, L.sgn);
fprintf('   clean 5 N cells per sheet: [%s]; propellant %.4f kg per ND of burn at 150 kg\n', ...
        num2str(squeeze(sum(sum(L.clean, 1), 2)).'), propKgPerNd);
% keep the package this run used: any file in <outDir>/inputs with its SHA-256
% is the snapshot; else copy it there as costate_catalog_gto_dro_<sha8>.mat
% (never overwriting)
[pkgSnapshot, pkgSnapshotStatus, snapMsg] = leg1_snapshot_package(L.identity, fullfile(outDir, 'inputs'));
fprintf('   %s\n', snapMsg);
nCk = numel(L.checks);
gates = addg(gates, 'I1', pass3(all([L.checks.pass])), nCk, sprintf('%d of %d leg1_inputs checks', ...
             nnz([L.checks.pass]), nCk));
secWall(1) = toc(tSec);

%% ========================================================================
%  2. SCHEDULE (catalog table, whole GTO revolutions; minimise T0)
%% ========================================================================
tSec = tic;
fprintf('\n2. SCHEDULE\n');
sOpts = schedOpts;
sOpts.propKgPerNd = propKgPerNd;
S0 = [];
try
    S0 = schedule_gto_dro(L.tf, L.sD, L.sA, L.sheetOrient, L.tauG, L.tauD, L.sgn, sOpts);
catch err
    if ~strcmp(err.identifier, 'schedule_gto_dro:noAssignment'), rethrow(err); end
    fprintf('   NO COMPLETE ASSIGNMENT: %s\n', err.message);
end
if isempty(S0)
    error('lunar_pnt_leg1_study:schedule', 'STOPPED AT SCHEDULE (C1 FAIL: no complete assignment)');
end
fprintf('   triple: orient0 %g deg, gtoPhase0 %.6f, delta %.6f; %d feasible triples; max orientation offset %.3f deg\n', ...
        S0.orient0, S0.gtoPhase0, S0.delta, S0.nFeasible, S0.maxOrientOffsetDeg);
fprintf('   T0 %.4f d, last arrival %.4f d, propellant %.2f kg (table, unpolished)\n', S0.T0*dayNd, ...
        S0.lastArr*dayNd, S0.totalPropKg);
print_schedule(S0, dayNd);
gates = addg(gates, 'C1', 'PASS', S0.T0*dayNd, 'complete assignment; value = table T0 [d]');
secWall(2) = toc(tSec);

%% ========================================================================
%  3. POLISH each satellite; on failure exclude and re-schedule that one
%% ========================================================================
tSec = tic;
if isa(offgridSolver, 'function_handle'), ogTxt = func2str_(offgridSolver); else, ogTxt = 'not configured'; end
fprintf('\n3. POLISH (at most %d candidates per satellite; off-grid solver %s)\n', maxAttempts, ogTxt);
sLock = sOpts;
sLock.lockTriple = struct('orient0', S0.orient0, 'gtoPhase0', S0.gtoPhase0, 'delta', S0.delta);
gOpts = struct('dryFrac', dryFrac, 'doSecond', doSecond);
S = S0;
exclude = zeros(6*(maxAttempts + 1), 4);  nEx = 0;    % rows [iD iA iS n], used 1:nEx
rec = containers.Map('KeyType', 'char', 'ValueType', 'any');
chosen = cell(1, 6);  slotOutcome = cell(1, 6);  attempts = cell(1, 6);
nCandSlot = zeros(1, 6);
for jj = 1:6
    kNow = find(S.slot == jj);
    nCandSlot(jj) = S.nCand(kNow);
    fprintf('   slot %d (now satellite %d, %d candidates):\n', jj, kNow, S.nCand(kNow));
    att = cell(1, maxAttempts + 1);  nAtt = 0;
    refKeys = cell(1, maxAttempts);  nRef = 0;
    while true
        kNow = find(S.slot == jj);
        c   = cand_of(S, kNow);
        key = ckey(c);
        if ~isKey(rec, key)
            rec(key) = apply_seam(evaluate(L, c, S, kNow, polishOpts, gOpts, 'catalog'), key, seam);
        end
        r = rec(key);
        nAtt = nAtt + 1;
        att{nAtt} = attempt_row(r, key);
        print_attempt(nAtt, r, dayNd);
        if r.blockOk && ~r.refuted
            outcome = 'accepted';
            break
        end
        if r.blockOk && r.refuted
            nRef = nRef + 1;  refKeys{nRef} = key;
        end
        nEx = nEx + 1;
        exclude(nEx, :) = [c.iD c.iA c.iS c.n];
        if nAtt >= maxAttempts
            outcome = 'exhausted';
            break
        end
        Snew = reschedule(L, sLock, exclude(1:nEx, :));
        if isempty(Snew)
            outcome = 'noCandidates';
            break
        end
        assert_others(S, Snew, jj);
        S = Snew;
    end
    if ~strcmp(outcome, 'accepted')
        % the last excluded row made the slot empty (noCandidates) or was never
        % re-scheduled (exhausted): S still points at it, so drop it again
        nEx = nEx - 1;
        if nRef > 0
            rr = rec(refKeys{1}).cand;
            keep = ~ismember(exclude(1:nEx, :), [rr.iD rr.iA rr.iS rr.n], 'rows');
            ex = exclude(1:nEx, :);  ex = ex(keep, :);
            nEx = size(ex, 1);  exclude(1:nEx, :) = ex;
            Snew = reschedule(L, sLock, exclude(1:nEx, :));
            assert(~isempty(Snew), 'lunar_pnt_leg1_study:refuted', 'restoring the refuted candidate failed');
            assert_others(S, Snew, jj);
            S = Snew;
            kNow = find(S.slot == jj);
            assert(strcmp(ckey(cand_of(S, kNow)), refKeys{1}), 'lunar_pnt_leg1_study:refuted', ...
                   'the re-schedule did not return the earliest refuted candidate');
            outcome = 'refuted';
            key = refKeys{1};
            fprintf('      no non-refuted candidate in %d attempt(s): keeping the earliest REFUTED one %s\n', ...
                    nAtt, key);
        elseif isa(offgridSolver, 'function_handle')
            r1 = rec(att{1}.key);
            c0 = r1.cand;
            fprintf('      %s after %d attempt(s): off-grid solver at %.3f deg, sD %.4f, sA %.4f\n', ...
                    outcome, nAtt, c0.orientDeg, L.sD(c0.iD), L.sA(c0.iA));
            [r, key] = offgrid_attempt(offgridSolver, L, c0, S, find(S.slot == jj), polishOpts, gOpts, r1.P);
            rec(key) = r;
            nAtt = nAtt + 1;
            att{nAtt} = attempt_row(r, key);
            print_attempt(nAtt, r, dayNd);
            if r.blockOk
                outcome = tern(r.refuted, 'offgridRefuted', 'offgrid');
            else
                outcome = 'failed';
            end
        else
            outcome = 'failed';
        end
        if strcmp(outcome, 'failed')
            key = att{nAtt}.key;
            fprintf('      SLOT %d FAILED: candidates exhausted (%d tried of %d; %s)\n', jj, nAtt, nCandSlot(jj), ...
                    tern(isa(offgridSolver, 'function_handle'), 'the off-grid solver did not close it', ...
                         'no off-grid solver configured'));
        end
    end
    chosen{jj} = key;
    slotOutcome{jj} = outcome;
    attempts{jj} = [att{1:nAtt}];
end
exclude = exclude(1:nEx, :);

% ---- T0 from the POLISHED arrivals, labels, final gates ----
okSlot = ~strcmp(slotOutcome, 'failed');
tArrPol = nan(1, 6);
for jj = find(okSlot), tArrPol(jj) = rec(chosen{jj}).P.tArr; end
xT = 6*(S.delta + max(tArrPol(okSlot))/L.tauD);
T0 = (ceil(xT - 1e-12)/6 - S.delta)*L.tauD;
kOf = mod(round(6*mod(S.delta + (0:5)/6 + T0/L.tauD, 1)), 6) + 1;
assert(isequal(sort(kOf), 1:6), 'lunar_pnt_leg1_study:labels', 'labels %s are not a permutation', mat2str(kOf));
fprintf('\n   T0 from the polished arrivals: %.6f d (table schedule %.6f d; last polished arrival %.6f d)\n', ...
        T0*dayNd, S.T0*dayNd, max(tArrPol(okSlot))*dayNd);
gates = addg(gates, 'T1', tern(abs(T0 - S.T0) < 1e-12, 'PASS', 'NOTE'), (T0 - S.T0)*dayNd, ...
             'hand-off epoch from the polished arrivals equals the table schedule''s; value = difference [d]');
P1c = cell(1, 6);
failedC = cell(1, 6);
for jj = 1:6
    r = rec(chosen{jj});
    k = kOf(jj);
    q = r.P;
    q.accepted = okSlot(jj);                % final: also Gf.pass and the alpha check below
    q.refuted = r.refuted;
    q.outcome = slotOutcome{jj};
    q.route = r.route;
    q.k = k;  q.slot = jj;
    q.iD = r.cand.iD;  q.iA = r.cand.iA;  q.iS = r.cand.iS;  q.n = r.cand.n;
    q.cand = r.cand;
    q.attempts = attempts{jj};
    q.orient0 = S.orient0;  q.gtoPhase0 = S.gtoPhase0;  q.delta = S.delta;
    q.G = r.G;
    q.t = [];  q.X = [];  q.alpha = [];  q.droPhaseArr = NaN;  q.alphaMissKm = NaN;  q.failReason = '';
    if okSlot(jj)
        ctx = struct('orient0', S.orient0, 'gtoPhase0', S.gtoPhase0, 'T0', T0 + seam.finalT0, 'kSat', k);
        Gf = leg1_gates(r.L, r.P, ctx, setfield_(gOpts, 'doSecond', false));
        Gf.second = r.G.second;
        Gf.refuted = r.G.refuted;
        Gf.secondUnknown = r.G.secondUnknown;
        q.G = Gf;
        F = r.P.flight;
        q.t = r.P.tDep + F.t(:).';
        q.X = F.Y(:, 1:7).';
        lv = F.Y(:, 11:13).';
        q.alpha = -lv./sqrt(sum(lv.^2, 1));
        q.droPhaseArr = Gf.handoff.phaseArr;
        q.alphaMissKm = fly_alpha(q.t, q.X, q.alpha, L);
        % ACCEPTANCE ON THE FINAL GATES: recomputed T0 and final k (hand-off),
        % and the stored thrust direction re-flying the flight
        names = {'flight', Gf.flight.ok; 'arrival', Gf.arrival.pass; 'oneGto', Gf.oneGto.pass; ...
                 'handoff at final T0', Gf.handoff.pass; 'earth', Gf.earth.pass; 'moon', Gf.moon.pass; ...
                 'mass', Gf.mass.feasible; 'thrust direction', q.alphaMissKm < alphaTolKm};
        bad = names(~[names{:, 2}], 1);
        q.accepted = isequal(Gf.pass, true) && isempty(bad);
        if ~q.accepted
            if isempty(bad), bad = {'G.pass false'}; end
            failedC{jj} = struct('k', k, 'slot', jj, 'reason', sprintf('final gates failed: %s', ...
                                 strjoin(bad.', ', ')), 'attempts', attempts{jj}, 'nCand', nCandSlot(jj));
            q.failReason = failedC{jj}.reason;
            fprintf('   *** satellite %d (slot %d) NOT ACCEPTED: %s ***\n', k, jj, q.failReason);
        end
    else
        % no final gates for a failed slot: its loop G carried a provisional
        % T0 and k, so the hand-off block is blanked (no stale k is saved)
        if isstruct(q.G) && isfield(q.G, 'handoff')
            q.G.handoff = struct('phaseArr', NaN, 'distKm', NaN, 'dVelMps', NaN, 'tfFlown', NaN, ...
                                 'tArrFlown', NaN, 'dtArr', NaN, 'phaseAtT0', NaN, 'target', NaN, ...
                                 'errPhase', NaN, 'reason', 'not evaluated at final T0', 'pass', false);
            q.G.pass = false;
        end
        failedC{jj} = struct('k', k, 'slot', jj, 'reason', ...
                             sprintf('candidates exhausted (%s): %s', ...
                                     tern(isa(offgridSolver, 'function_handle'), 'the off-grid solver did not close it', ...
                                          'no off-grid solver configured'), ...
                                     strjoin(arrayfun(@(a) sprintf('%s: %s', a.key, a.reason), attempts{jj}, ...
                                                      'UniformOutput', false), '; ')), ...
                             'attempts', attempts{jj}, 'nCand', nCandSlot(jj));
        q.failReason = failedC{jj}.reason;
    end
    P1c{k} = q;
end
P1 = [P1c{:}];
acc = [P1.accepted];
failed = [failedC{:}];
if isempty(failed), failed = struct('k', {}, 'slot', {}, 'reason', {}, 'attempts', {}, 'nCand', {}); end
secWall(3) = toc(tSec);

%% ========================================================================
%  4. GATES TABLE (each computed here, PASS / FAIL)
%% ========================================================================
tSec = tic;
fprintf('\n4. GATES per satellite (leg-2 numbering k; T0 = %.6f d)\n', T0*dayNd);
fprintf('   %-2s %-4s %-9s %-2s %-8s %8s %8s %8s %8s %9s %9s %9s %9s %10s %9s %-13s %s\n', 'k', 'slot', 'cell', 'n', ...
        'route', 'tDep d', 'tArr d', 'tf d', 'orient', 'miss km', 'miss m/s', 'dGTO ND', 'handoff', ...
        'Earth km', 'Moon km', 'second', 'pass');
for k = 1:6
    q = P1(k);
    if ~q.accepted
        fprintf('   %-2d %-4d %-9s %-2d %-8s FAILED: %s\n', k, q.slot, sprintf('(%d,%d,%d)', q.iD, q.iA, q.iS), ...
                q.n, q.route, q.failReason);
        continue
    end
    G = q.G;
    fprintf('   %-2d %-4d %-9s %-2d %-8s %8.4f %8.4f %8.4f %+8.3f %9.3g %9.3g %9.2g %+9.2e %10.1f %9.1f %-13s %d\n', ...
            k, q.slot, sprintf('(%d,%d,%d)', q.iD, q.iA, q.iS), q.n, q.route, q.tDep*dayNd, q.tArr*dayNd, ...
            q.tf*dayNd, q.orientDeg, G.arrival.missKm, G.arrival.missMps, G.oneGto.dStateNd, G.handoff.errPhase, ...
            G.earth.hMinKm, G.moon.hMinKm, G.second.status, G.pass);
end
fprintf('   located minima: Earth at t [%s] ND after departure, Moon at t [%s] ND\n', ...
        num2str(arrayfun(@(q) gnum(q.G, 'earth', 'tAtMin'), P1(acc)), 5), ...
        num2str(arrayfun(@(q) gnum(q.G, 'moon', 'tAtMin'), P1(acc)), 5));
Ga = [P1(acc).G];
gates = addg(gates, 'P6', pass3(all(acc)), nnz(acc), sprintf('%d of 6 satellites closed (failed: %s)', nnz(acc), ...
             tern(all(acc), 'none', strtrim(sprintf('%d ', find(~acc))))));
if any(acc)
    ar = [Ga.arrival];  og = [Ga.oneGto];  ho = [Ga.handoff];  ea = [Ga.earth];  mo = [Ga.moon];
    fl = [Ga.flight];   ms = [Ga.mass];    se = [Ga.second];
    gates = addg(gates, 'GA', pass3(all([ar.pass])), max([ar.missKm]), sprintf(['re-flown arrival miss < 100 km and ' ...
                 '< 10 m/s; value = max km (max %.3g m/s)'], max([ar.missMps])));
    gates = addg(gates, 'GF', pass3(all([fl.ok])), nnz([fl.ok]), 'flight admissible (validate_flight); value = count');
    gates = addg(gates, 'GO', pass3(all([og.pass])), max([og.dStateNd]), ...
                 'departure = the one GTO''s state at tDep (recomputed from orient0, gtoPhase0); value = max ND');
    gates = addg(gates, 'GH', pass3(all([ho.pass])), max(abs([ho.errPhase])), sprintf(['DRO phase of the flown arrival ' ...
                 'at T0 = (k-1)/6 (tol %.0e); value = max |error| [period]'], Ga(1).opts.handoffTolPhase));
    gates = addg(gates, 'GE', pass3(all([ea.pass])), min([ea.hMinKm]), sprintf('located minimum Earth altitude >= %g km; value = min km', ...
                 ea(1).floorKm));
    gates = addg(gates, 'GM', pass3(all([mo.pass])), min([mo.hMinKm]), sprintf('located minimum Moon altitude >= %g km; value = min km', ...
                 mo(1).floorKm));
    gates = addg(gates, 'GW', pass3(all([ms.feasible])), max([ms.stagePropKg]), 'stage mass closes; value = max stage propellant [kg]');
    st = {se.status};
    if all(strcmp(st, 'pass'))
        s2 = 'PASS';
    elseif any(strcmp(st, 'refuted'))
        s2 = 'REFUTED';
    else
        s2 = 'UNRESOLVED';
    end
    gates = addg(gates, 'G2', s2, nnz(strcmp(st, 'pass')), sprintf(['second order (reported, not blocking): ' ...
                 'statuses [%s]; value = count pass'], strjoin(st, ' ')));
end
tArrAcc = [P1(acc).tArr];
gtOk = isempty(tArrAcc) || (T0 >= max(tArrAcc) - 1e-12);
phT0 = mod(S.delta + T0/L.tauD, 1/6);
gtOk = gtOk && min(phT0, 1/6 - phT0) < 1e-9;
gates = addg(gates, 'GT', pass3(gtOk), (T0 - max([tArrAcc -Inf]))*dayNd, ...
             'T0 >= every arrival and a slot instant; value = T0 - last arrival [d]');
% GD: the stored thrust directions reproduce the flight (CR3BP + (Tnd/m)*alpha,
% mdot = -Tnd/cnd, ode113 1e-12), and the opposite sign does not (deliberate failure)
gdKm = nan(1, 6);
for k = find(acc)
    gdKm(k) = P1(k).alphaMissKm;           % computed with the final gates (section 3)
end
if any(acc)
    k1 = find(acc, 1);
    flipKm = fly_alpha(P1(k1).t, P1(k1).X, -P1(k1).alpha, L);
    fprintf('   GD: re-integration with alpha misses the flown arrival by [%s] km; with -alpha (sat %d) %.0f km\n', ...
            num2str(gdKm(acc), 3), k1, flipKm);
    gates = addg(gates, 'GD', pass3(max(gdKm(acc)) < alphaTolKm && flipKm > 1000), max(gdKm(acc)), ...
                 sprintf(['thrust direction alpha = -lambda_v/|lambda_v| re-flies the flight (< %g km); the opposite ' ...
                          'sign misses by %.0f km; value = max km'], alphaTolKm, flipKm));
end
fprintf('\n   %-4s %-10s %12s  %s\n', 'gate', 'status', 'value', 'detail');
for kg = 1:size(gates, 1)
    fprintf('   %-4s %-10s %12.5g  %s\n', gates{kg, 1}, gates{kg, 2}, gates{kg, 3}, gates{kg, 4});
end
blocking = ~ismember(gates(:, 1), {'G2', 'T1'});
failedIds = gates(blocking & ~strcmp(gates(:, 2), 'PASS'), 1);
if isempty(failedIds)
    verdict = 'PASS';
else
    verdict = ['not passed: ' strjoin(failedIds.', ', ')];
end
info2 = gates(strcmp(gates(:, 1), 'G2'), 2);
if ~isempty(info2) && ~strcmp(info2{1}, 'PASS')
    verdict = sprintf('%s (second order %s, reported)', verdict, info2{1});
end
if seamUsed
    verdict = sprintf('%s (TEST SEAMS: refused %s; refuted %s; final T0 offset %g ND)', verdict, ...
                      strjoin(seam.refuse, ' '), strjoin(seam.refute, ' '), seam.finalT0);
end
fprintf('   VERDICT: %s\n', verdict);
secWall(4) = toc(tSec);

%% ========================================================================
%  5. MASS AND POWER
%% ========================================================================
tSec = tic;
fprintf('\n5. MASS AND POWER (payload %g kg, dryFrac %g)\n', L.m0kg, dryFrac);
fprintf('   %-2s %8s %9s %9s %9s %9s %8s %8s %8s\n', 'k', 'tf d', 'mf/m0', 'prop kg', 'dry kg', 'wet kg', ...
        'thrust N', 'jet kW', 'dV km/s');
for k = find(acc)
    m = P1(k).G.mass;
    fprintf('   %-2d %8.4f %9.6f %9.3f %9.3f %9.3f %8.3f %8.2f %8.4f\n', k, P1(k).tf*dayNd, m.mfFrac, ...
            m.stagePropKg, m.stageDryKg, m.wetKg, m.thrustN, m.jetKW, m.dvKms);
end
if any(acc)
    mm = [Ga.mass];
    massTot = struct('stagePropKg', sum([mm.stagePropKg]), 'stageDryKg', sum([mm.stageDryKg]), ...
                     'wetKg', sum([mm.wetKg]), 'maxThrustN', max([mm.thrustN]), 'maxJetKW', max([mm.jetKW]), ...
                     'sumTfDays', sum([P1(acc).tf])*dayNd);
    fprintf('   total: stage propellant %.3f kg, stage dry %.3f kg, stack wet %.3f kg (launch mass of the %d stacks); ', ...
            massTot.stagePropKg, massTot.stageDryKg, massTot.wetKg, nnz(acc));
    fprintf('largest thrust %.3f N, largest jet power %.2f kW\n', massTot.maxThrustN, massTot.maxJetKW);
else
    massTot = struct();
end
caveat = sprintf(['   CAVEAT: these are minimum-time, engine-always-on transfers in the catalog''s numerical regime ' ...
                  '(%g N per %g kg, Isp %g s), scaled to a stack that arrives at %g kg with dryFrac %g (an ' ...
                  'assumption); they are not a vehicle design.'], L.thrustN, L.m0kg, L.ispS, L.m0kg, dryFrac);
fprintf('%s\n', caveat);
secWall(5) = toc(tSec);

%% ========================================================================
%  6. SAVE
%% ========================================================================
tSec = tic;
[~, hostName] = system('hostname');
s64 = (0:63)/64;
knobsSave = knobs;
if isa(knobsSave.offgridSolver, 'function_handle'), knobsSave.offgridSolver = func2str_(knobsSave.offgridSolver); end
run_meta = struct('date', char(datetime('now', 'Format', 'yyyy-MM-dd HH:mm')), 'matlab', version, ...
                  'host', strtrim(hostName), 'optimal_control', gitOC.head, ...
                  'optimal_controlDirty', gitOC.dirty, 'git', gitOC, 'knobs', knobsSave, ...
                  'identity', L.identity, 'packageSnapshot', pkgSnapshot, ...
                  'packageSnapshotStatus', pkgSnapshotStatus, 'droMd5', md5_of(L.slotState(s64)), 'droPhase0', (0:5)/6, ...
                  'tauD', L.tauD, 'tauG', L.tauG, 'sgn', L.sgn, ...
                  'regime', struct('thrustN', L.thrustN, 'ispS', L.ispS, 'm0kg', L.m0kg, 'muStar', L.muStar, ...
                                   'lStar', L.lStar, 'tStar', L.tStar, 'Tnd', L.phys.Tnd, 'cnd', L.phys.cnd), ...
                  'verdict', verdict, 'caveat', strtrim(caveat), 'sectionNames', {secName}, ...
                  'sectionWallS', secWall, 'studyWallS', toc(tStudy));
saveTarget = fullfile(outDir, 'lunar_pnt_leg1_study.mat');
T0sched = S.T0;
knobs = knobsSave;
save(saveTarget, 'S', 'S0', 'P1', 'T0', 'T0sched', 'knobs', 'run_meta', 'failed', 'gates', 'exclude', ...
     'massTot', 'slotOutcome', '-v7.3');
fprintf('\n6. SAVED %s\n', saveTarget);
secWall(6) = toc(tSec);

%% ========================================================================
%  7. CAPTURE (sidecar of the accepted transfers)
%% ========================================================================
tSec = tic;
sidecar = fullfile(outDir, 'leg1_transfers_gto_dro_5N.mat');
exportRec = struct('status', 'skipped', 'file', '', 'reason', '');
if ~exportTransfers
    exportRec.reason = 'exportTransfers = false';
    fprintf('\n7. CAPTURE: skipped (exportTransfers = false)\n');
elseif seamUsed
    exportRec.reason = 'test seams set: forced verdicts are not exported';
    fprintf('\n7. CAPTURE: skipped (test seams set; forced verdicts are never exported)\n');
elseif ~any(acc)
    exportRec.reason = 'no accepted transfer';
    fprintf('\n7. CAPTURE: nothing to export\n');
else
    fprintf('\n7. CAPTURE\n');
    try
        Ex = export_leg1_transfers(saveTarget, struct('outFile', sidecar, 'verbose', true));
        exportRec = struct('status', 'written', 'file', sidecar, 'reason', sprintf('%d entries', numel(Ex.entries)));
    catch xerr
        fprintf('   *** EXPORT FAILED: %s: %s ***\n', xerr.identifier, xerr.message);
        exportRec = struct('status', 'failed', 'file', '', 'reason', sprintf('%s: %s', xerr.identifier, xerr.message));
    end
end
secWall(7) = toc(tSec);

%% ========================================================================
%  8. FIGURES (slide-size fonts, PNG)
%% ========================================================================
tSec = tic;
timelineRec = struct('status', 'skipped', 'reason', '');
if makeFigures && any(acc)
    fprintf('\n8. FIGURES -> %s\n', outDir);
    figs = gobjects(0);
    try
        FS = 14;
        cGto = [0.65 0.65 0.65];  cXfer = [0.85 0.33 0.10];  cDro = [0.30 0.50 0.85];
        cSc = lines(6);
        lab = arrayfun(@(q) sprintf('sat %d', q.k), P1, 'UniformOutput', false);

        % leg1_schedule.png
        figs(end+1) = figure('Theme', 'light', 'Color', 'w', 'Position', [60 60 1280 600], 'Visible', figVis);
        ax = axes(figs(end));  hold(ax, 'on');
        for k = find(acc)
            q = P1(k);
            hv = tern(k == find(acc, 1), 'on', 'off');
            plot(ax, [0 q.tDep]*dayNd, [k k], '-', 'Color', cGto, 'LineWidth', 16, 'DisplayName', 'GTO wait', 'HandleVisibility', hv);
            plot(ax, [q.tDep q.tArr]*dayNd, [k k], '-', 'Color', cXfer, 'LineWidth', 16, 'DisplayName', 'leg 1 (5 N min-time)', 'HandleVisibility', hv);
            plot(ax, [q.tArr T0]*dayNd, [k k], '-', 'Color', cDro, 'LineWidth', 16, 'DisplayName', 'DRO coast to T0', 'HandleVisibility', hv);
            text(ax, q.tArr*dayNd - 0.05, k - 0.36, sprintf('cell (%d,%d) n %d, slot %d', q.iD, q.iA, q.n, q.slot), ...
                 'FontSize', FS - 3, 'HorizontalAlignment', 'right');
        end
        xline(ax, T0*dayNd, '--k', 'LineWidth', 1.4, 'DisplayName', sprintf('hand-off T0 = %.3f d', T0*dayNd));
        set(ax, 'YTick', 1:6, 'YTickLabel', lab, 'YDir', 'reverse', 'FontSize', FS, 'Color', 'w');
        xlim(ax, [0 T0*dayNd*1.04]);  ylim(ax, [0.4 6.6]);  grid(ax, 'on');
        xlabel(ax, 'leg-1 time [d]');
        legend(ax, 'Location', 'southoutside', 'Orientation', 'horizontal', 'FontSize', FS);
        title(ax, sprintf('Leg 1: six 5 N stages from one GTO to the DRO slots, hand-off at %.3f d', T0*dayNd), ...
              'FontSize', FS + 1);
        exportgraphics(figs(end), fullfile(outDir, 'leg1_schedule.png'), 'Resolution', 150);

        % leg1_3d.png
        figs(end+1) = figure('Theme', 'light', 'Color', 'w', 'Position', [60 60 1280 760], 'Visible', figVis);
        ax = axes(figs(end));  hold(ax, 'on');
        sPl = linspace(0, 1, 721);
        rD = L.slotState(sPl);
        plot3(ax, rD(1, :), rD(2, :), rD(3, :), '-', 'Color', cDro, 'LineWidth', 1.8, 'DisplayName', ...
              sprintf('DRO (%.2f d)', L.tauD*dayNd));
        oU = unique(round([P1(acc).orientDeg], 9));
        for ko = 1:numel(oU)
            rG = L.gtoState(oU(ko), sPl);
            plot3(ax, rG(1, :), rG(2, :), rG(3, :), '-', 'Color', cGto, 'LineWidth', 1.0, ...
                  'DisplayName', 'GTO at the departure orientations', 'HandleVisibility', tern(ko == 1, 'on', 'off'));
        end
        for k = find(acc)
            X = P1(k).X;
            plot3(ax, X(1, :), X(2, :), X(3, :), '-', 'Color', cSc(k, :), 'LineWidth', 1.1, ...
                  'DisplayName', sprintf('sat %d (%+.2f deg)', k, P1(k).orientDeg));
            plot3(ax, X(1, end), X(2, end), X(3, end), 'o', 'Color', cSc(k, :), 'MarkerFaceColor', cSc(k, :), ...
                  'MarkerSize', 7, 'HandleVisibility', 'off');
        end
        [sx, sy, sz] = sphere(48);
        rE = 6378.137/L.lStar;  rMn = 1737.4/L.lStar;
        surf(ax, -L.muStar + rE*sx, rE*sy, rE*sz, 'FaceColor', [0.20 0.45 0.85], 'EdgeColor', 'none', ...
             'DisplayName', 'Earth (to scale)');
        surf(ax, 1 - L.muStar + rMn*sx, rMn*sy, rMn*sz, 'FaceColor', [0.55 0.55 0.55], 'EdgeColor', 'none', ...
             'DisplayName', 'Moon (to scale)');
        axis(ax, 'equal');  grid(ax, 'on');  view(ax, [-15 55]);
        set(ax, 'FontSize', FS, 'Color', 'w');
        xlabel(ax, 'x [ND]');  ylabel(ax, 'y [ND]');  zlabel(ax, 'z [ND]');
        legend(ax, 'Location', 'eastoutside', 'FontSize', FS - 1);
        title(ax, 'Leg 1 in the Earth-Moon rotating frame: GTO -> DRO slots, 5 N min-time', 'FontSize', FS + 1);
        exportgraphics(figs(end), fullfile(outDir, 'leg1_3d.png'), 'Resolution', 150);

        % leg1_mass.png
        figs(end+1) = figure('Theme', 'light', 'Color', 'w', 'Position', [60 60 1280 640], 'Visible', figVis);
        ax = axes(figs(end));  hold(ax, 'on');
        mm = nan(6, 3);
        for k = find(acc)
            m = P1(k).G.mass;
            mm(k, :) = [L.m0kg, m.stageDryKg, m.stagePropKg];
        end
        hb = bar(ax, 1:6, mm, 0.6, 'stacked');
        hb(1).FaceColor = [0.30 0.60 0.35];  hb(1).DisplayName = sprintf('payload (%g kg)', L.m0kg);
        hb(2).FaceColor = [0.50 0.50 0.50];  hb(2).DisplayName = sprintf('stage dry (%g x propellant)', dryFrac);
        hb(3).FaceColor = cXfer;             hb(3).DisplayName = 'stage propellant';
        for k = find(acc)
            m = P1(k).G.mass;
            text(ax, k, m.wetKg + 8, sprintf('%.1f kg\n%.2f N\n%.1f kW', m.wetKg, m.thrustN, m.jetKW), ...
                 'HorizontalAlignment', 'center', 'VerticalAlignment', 'bottom', 'FontSize', FS - 3);
        end
        set(ax, 'XTick', 1:6, 'XTickLabel', lab, 'FontSize', FS, 'Color', 'w');
        ylim(ax, [0 1.25*max(sum(mm, 2))]);
        ylabel(ax, 'stack mass at departure [kg]');  grid(ax, 'on');
        legend(ax, 'Location', 'northoutside', 'Orientation', 'horizontal', 'FontSize', FS);
        title(ax, sprintf(['Stacks sized to arrive at %g kg: total wet %.0f kg, propellant %.0f kg ' ...
                           '(min-time, engine always on; not a vehicle design)'], L.m0kg, massTot.wetKg, ...
                          massTot.stagePropKg), 'FontSize', FS);
        exportgraphics(figs(end), fullfile(outDir, 'leg1_mass.png'), 'Resolution', 150);
        fprintf('   leg1_schedule.png, leg1_3d.png, leg1_mass.png written\n');
    catch err
        if headless, close(figs(isgraphics(figs))); end
        rethrow(err);
    end
    if headless, close(figs(isgraphics(figs))); end

    % mission_timeline.png (combined_timeline checks compatibility first)
    if all(acc) && isfile(knobs.leg2Mat)
        try
            [Tline, tlInfo] = combined_timeline(saveTarget, knobs.leg2Mat, struct('catalogMat', L.identity.file, ...
                                                'outPng', fullfile(outDir, 'mission_timeline.png')));
            fprintf('   mission_timeline.png written (%s)\n', strjoin(tlInfo.checks, '; '));
            fprintf('   combined timeline [d]: hand-off %.3f, leg-2 last arrival %.3f, service end %.3f\n', ...
                    tlInfo.T0Days, max(Tline.leg2(:, 2)), tlInfo.endDays);
            timelineRec = struct('status', 'written', 'reason', '');
        catch terr
            fprintf('   *** COMBINED TIMELINE REFUSED: %s: %s ***\n', terr.identifier, terr.message);
            timelineRec = struct('status', 'refused', 'reason', sprintf('%s: %s', terr.identifier, terr.message));
        end
    else
        timelineRec.reason = 'not all six accepted, or no leg-2 results';
        fprintf('   mission_timeline.png: not drawn (%s)\n', timelineRec.reason);
    end
end
secWall(8) = toc(tSec);
run_meta.export = exportRec;
run_meta.timeline = timelineRec;
run_meta.sectionWallS = secWall;
run_meta.studyWallS = toc(tStudy);
save(saveTarget, 'run_meta', '-append');

fprintf('\n== SUMMARY  lunar PNT leg 1: six 5 N stages, one GTO -> the DRO slots\n');
fprintf('   T0 %.4f d; %d of 6 closed; total stage propellant %.2f kg; VERDICT %s\n', T0*dayNd, nnz(acc), ...
        gfield(massTot, 'stagePropKg'), verdict);
fprintf(['   SCOPE: the second-order verdict and "minimum-time" refer to the fixed-endpoint problem to the ' ...
         'slot''s arrival point (a local test at the root); the moving-target rendezvous condition is not checked.\n']);
fprintf('   wall per section [s]:');
for ks = 1:numel(secName), fprintf('  %s %.1f', secName{ks}, secWall(ks)); end
fprintf('  (total %.0f)\n', toc(tStudy));


%% ------------------------------------------------------------------------
function c = cand_of(S, k)
% CAND_OF  The chosen candidate of satellite k of a schedule, as a polish
%   candidate.  INPUTS: S schedule; k index.  OUTPUTS: c struct.
c = struct('iD', S.iD(k), 'iA', S.iA(k), 'iS', S.iS(k), 'n', S.n(k), 'tDep', S.tDep(k), ...
           'tArr', S.tArr(k), 'tf', S.tf(k), 'orientDeg', S.orientDeg(k), ...
           'slotPhaseAtDep', S.slotPhaseAtDep(k));
end

function key = ckey(c)
% CKEY  Candidate identity 'iD_iA_iS_n'.  INPUTS: c struct.  OUTPUTS: key char.
key = sprintf('%d_%d_%d_%d', c.iD, c.iA, c.iS, c.n);
end

function r = evaluate(L, c, S, k, polishOpts, gOpts, route)
% EVALUATE  Polish one candidate and gate it with the current schedule's
%   context. Blocking: polish accepted, flight admissible, arrival, one
%   GTO, Earth and Moon (located), mass feasible. The hand-off gate depends
%   on T0 and is decided after the final schedule (not blocking here).
%   INPUTS: L leg-1 inputs; c candidate; S schedule; k its satellite index;
%   polishOpts, gOpts structs; route char.
%   OUTPUTS: r struct (.P .G .cand .route .blockOk .refuted .reason .L).
t0 = tic;
P = polish_gto_dro(L, c, polishOpts);
r = struct('P', P, 'G', [], 'cand', c, 'route', route, 'blockOk', false, 'refuted', false, ...
           'reason', '', 'second', 'n/a', 'L', L, 'wallS', NaN);
if ~P.ok
    r.reason = ['polish: ' P.reason];
else
    ctx = struct('orient0', S.orient0, 'gtoPhase0', S.gtoPhase0, 'T0', S.T0, 'kSat', k);
    G = leg1_gates(L, P, ctx, gOpts);
    r.G = G;
    r.second = G.second.status;
    names = {'flight', G.flight.ok; 'arrival', G.arrival.pass; 'oneGto', G.oneGto.pass; ...
             'earth', G.earth.pass; 'moon', G.moon.pass; 'mass', G.mass.feasible};
    bad = names(~[names{:, 2}], 1);
    r.blockOk = isempty(bad);
    r.refuted = G.refuted;
    if r.blockOk
        r.reason = sprintf('gates pass; second order %s', G.second.status);
    else
        r.reason = sprintf('blocking gate(s) failed: %s', strjoin(bad.', ', '));
        if any(strcmp(bad, 'earth'))
            r.reason = sprintf('%s (Earth %.1f km)', r.reason, G.earth.hMinKm);
        end
        if any(strcmp(bad, 'moon'))
            r.reason = sprintf('%s (Moon %.1f km)', r.reason, G.moon.hMinKm);
        end
    end
end
r.wallS = toc(t0);
end

function r = apply_seam(r, key, seam)
% APPLY_SEAM  TEST SEAM: a candidate listed in seam.refuse is treated as
%   refused (blocking), one in seam.refute as refuted (its gates' second
%   order marked refuted) when its real gates pass. Never used by default.
%   INPUTS: r evaluate output; key char; seam struct (.refuse .refute
%   cellstr).  OUTPUTS: r.
if any(strcmp(key, seam.refuse))
    r.blockOk = false;
    r.reason = ['TEST SEAM: refused (' r.reason ')'];
elseif any(strcmp(key, seam.refute)) && r.blockOk
    r.refuted = true;
    r.G.refuted = true;
    r.G.second.status = 'refuted';
    r.G.second.reason = ['TEST SEAM: marked refuted; real: ' r.G.second.reason];
    r.second = 'refuted';
    r.reason = 'gates pass; second order refuted (TEST SEAM)';
end
end

function [r, key] = offgrid_attempt(solver, L, c0, S, k, polishOpts, gOpts, Ptemplate)
% OFFGRID_ATTEMPT  The fallback: solve the cell at its true orientation
%   with solver(L, orientDeg, sD, sA), then polish it as a one-cell sheet
%   at that orientation (no walk) and gate it.  INPUTS: solver handle; L;
%   c0 first candidate of the slot; S; k; polishOpts, gOpts; Ptemplate a
%   polish struct (its field set, for a refusal record).
%   OUTPUTS: r (as evaluate, route 'offgrid'); key char.
key = sprintf('offgrid_%s', ckey(c0));
E = solver(L, c0.orientDeg, L.sD(c0.iD), L.sA(c0.iA));
if ~(isstruct(E) && isfield(E, 'ok') && isequal(E.ok, true) && numel(E.z8) == 8)
    why = 'off-grid solver did not return a solution';
    if isstruct(E) && isfield(E, 'reason'), why = sprintf('%s (%s)', why, E.reason); end
    P = Ptemplate;
    P.ok = false;  P.reason = why;  P.z8 = zeros(0, 1);  P.tf = NaN;  P.tArr = NaN;
    r = struct('P', P, 'G', [], 'cand', c0, 'route', 'offgrid', 'blockOk', false, 'refuted', false, ...
               'reason', why, 'second', 'n/a', 'L', L, 'wallS', NaN);
    return
end
L2 = L;
nS = numel(L.sheetOrient) + 1;
L2.sheetOrient(nS) = c0.orientDeg;
L2.tf(:, :, nS) = NaN;     L2.tf(c0.iD, c0.iA, nS) = E.z8(8);
L2.clean(:, :, nS) = false; L2.clean(c0.iD, c0.iA, nS) = true;
L2.z8(:, :, :, nS) = NaN;  L2.z8(:, c0.iD, c0.iA, nS) = E.z8(:);
c = c0;  c.iS = nS;
r = evaluate(L2, c, S, k, polishOpts, gOpts, 'offgrid');
r.cand.iS = c0.iS;                       % the seed sheet for the record
end

function a = attempt_row(r, key)
% ATTEMPT_ROW  One attempt as a flat struct.  INPUTS: r evaluate output;
%   key char.  OUTPUTS: a struct.
a = struct('key', key, 'iD', r.cand.iD, 'iA', r.cand.iA, 'iS', r.cand.iS, 'n', r.cand.n, ...
           'tDep', r.cand.tDep, 'orientDeg', r.cand.orientDeg, 'route', r.route, ...
           'polishOk', r.P.ok, 'blockOk', r.blockOk, 'refuted', r.refuted, 'second', r.second, ...
           'reason', r.reason, 'wallS', r.wallS);
end

function print_attempt(n, r, dayNd)
% PRINT_ATTEMPT  One line per attempt.  INPUTS: n attempt number; r;
%   dayNd [d per ND].  OUTPUTS: none (prints).
c = r.cand;
if r.P.ok
    G = r.G;
    fprintf(['      %d. cell (%d,%d,%d) n %d %s, tDep %.4f d, orient %+.3f deg: tf %.5f d (table %.5f), ' ...
             'miss %.3g km %.3g m/s, Earth %.1f km, Moon %.1f km, second %s -> %s (%.1f s)\n'], n, c.iD, c.iA, ...
            c.iS, c.n, r.route, c.tDep*dayNd, c.orientDeg, r.P.tf*dayNd, r.P.tfCell*dayNd, G.arrival.missKm, ...
            G.arrival.missMps, G.earth.hMinKm, G.moon.hMinKm, G.second.status, ...
            tern(r.blockOk, tern(r.refuted, 'REFUTED (non-blocking)', 'ACCEPTED'), 'REFUSED'), r.wallS);
else
    fprintf('      %d. cell (%d,%d,%d) n %d %s, tDep %.4f d, orient %+.3f deg: REFUSED: %s (%.1f s)\n', n, c.iD, ...
            c.iA, c.iS, c.n, r.route, c.tDep*dayNd, c.orientDeg, r.reason, r.wallS);
end
end

function S = reschedule(L, sLock, exclude)
% RESCHEDULE  The locked-triple schedule under an exclusion list, or []
%   when the slot has no candidate left.  INPUTS: L; sLock opts; exclude
%   [m x 4].  OUTPUTS: S schedule or [].
sLock.exclude = exclude;
try
    S = schedule_gto_dro(L.tf, L.sD, L.sA, L.sheetOrient, L.tauG, L.tauD, L.sgn, sLock);
catch err
    if ~strcmp(err.identifier, 'schedule_gto_dro:noAssignment'), rethrow(err); end
    S = [];
end
end

function assert_others(Sold, Snew, jSlot)
% ASSERT_OTHERS  A locked re-schedule must leave every other slot's chosen
%   candidate unchanged.  INPUTS: Sold, Snew schedules; jSlot the slot
%   re-scheduled.  OUTPUTS: none (throws lunar_pnt_leg1_study:lock).
for jj = setdiff(1:6, jSlot)
    a = cand_of(Sold, find(Sold.slot == jj));
    b = cand_of(Snew, find(Snew.slot == jj));
    if ~(strcmp(ckey(a), ckey(b)) && a.tDep == b.tDep)
        error('lunar_pnt_leg1_study:lock', 're-scheduling slot %d changed slot %d (%s -> %s)', jSlot, jj, ...
              ckey(a), ckey(b));
    end
end
end

function print_schedule(S, dayNd)
% PRINT_SCHEDULE  The schedule table.  INPUTS: S; dayNd.  OUTPUTS: none.
fprintf('   %-2s %-4s %-3s %-3s %-3s %-2s %9s %9s %9s %9s %8s %6s\n', 'k', 'slot', 'iD', 'iA', 'iS', 'n', ...
        'tDep d', 'tArr d', 'tf d', 'orient', 'resid', 'nCand');
for k = 1:6
    fprintf('   %-2d %-4d %-3d %-3d %-3d %-2d %9.4f %9.4f %9.4f %+9.3f %8.4f %6d\n', k, S.slot(k), S.iD(k), ...
            S.iA(k), S.iS(k), S.n(k), S.tDep(k)*dayNd, S.tArr(k)*dayNd, S.tf(k)*dayNd, S.orientDeg(k), ...
            S.resid(k), S.nCand(k));
end
end

function km = fly_alpha(t, X, alpha, L)
% FLY_ALPHA  Re-integrate the CR3BP with thrust acceleration (Tnd/m)*alpha
%   (alpha interpolated by pchip, renormalised) and mass rate -Tnd/cnd from
%   the first stored state; position miss at the last stored epoch [km].
%   INPUTS: t [1 x M] ND; X [7 x M]; alpha [3 x M]; L.  OUTPUTS: km.
mu = L.muStar;  Tnd = L.phys.Tnd;  cnd = L.phys.cnd;
pp = pchip(t, alpha);
rhs = @(tt, x) thrust_rhs(tt, x, pp, mu, Tnd, cnd);
op = odeset('RelTol', 1e-12, 'AbsTol', 1e-13);
[~, Y] = ode113(rhs, [t(1) t(end)], X(:, 1), op);
km = L.lStar*sqrt(sum((Y(end, 1:3).' - X(1:3, end)).^2));
end

function dx = thrust_rhs(tt, x, pp, mu, Tnd, cnd)
% THRUST_RHS  CR3BP (barycentric rotating, Earth at -mu, Moon at 1 - mu)
%   with thrust.  INPUTS: tt ND; x [7 x 1]; pp pchip of alpha; mu; Tnd;
%   cnd.  OUTPUTS: dx [7 x 1].
a = ppval(pp, tt);
a = a/sqrt(sum(a.^2));
r1 = sqrt((x(1) + mu)^2 + x(2)^2 + x(3)^2);
r2 = sqrt((x(1) - 1 + mu)^2 + x(2)^2 + x(3)^2);
ax = 2*x(5) + x(1) - (1 - mu)*(x(1) + mu)/r1^3 - mu*(x(1) - 1 + mu)/r2^3;
ay = -2*x(4) + x(2) - (1 - mu)*x(2)/r1^3 - mu*x(2)/r2^3;
az = -(1 - mu)*x(3)/r1^3 - mu*x(3)/r2^3;
dx = [x(4:6); [ax; ay; az] + (Tnd/x(7))*a; -Tnd/cnd];
end

function c = canon(p)
% CANON  Canonical absolute path (need not exist).  INPUTS: p char.
%   OUTPUTS: c char.
p = char(p);
if ~java.io.File(p).isAbsolute(), p = fullfile(pwd, p); end
c = char(java.io.File(p).getCanonicalPath());
end

function s = pass3(c)
% PASS3  PASS / FAIL from a logical.  INPUTS: c logical.  OUTPUTS: s char.
if c, s = 'PASS'; else, s = 'FAIL'; end
end

function s = tern(c, a, b)
% TERN  Inline conditional.  INPUTS: c logical; a, b.  OUTPUTS: s.
if c, s = a; else, s = b; end
end

function G = addg(G, id, status, value, detail)
% ADDG  Append one gate row {id, status, value, detail}.
%   INPUTS: G {n x 4}; id, status, detail char; value scalar.  OUTPUTS: G.
G(end+1, :) = {id, status, value, detail};
end

function s = setfield_(s, name, value)
% SETFIELD_  Copy of struct s with field name set.  INPUTS: s; name; value.  OUTPUTS: s.
s.(name) = value;
end

function v = gfield(s, name)
% GFIELD  s.(name), or NaN when absent.  INPUTS: s struct; name.  OUTPUTS: v.
if isstruct(s) && isfield(s, name) && ~isempty(s.(name)), v = s.(name); else, v = NaN; end
end

function v = gnum(G, a, b)
% GNUM  G.(a).(b) or NaN.  INPUTS: G struct; a, b char.  OUTPUTS: v.
if isstruct(G) && isfield(G, a) && isfield(G.(a), b), v = G.(a).(b); else, v = NaN; end
end

function s = func2str_(h)
% FUNC2STR_  Text of a function handle.  INPUTS: h handle.  OUTPUTS: s char.
s = char(func2str(h));
end

function h = md5_of(x)
% MD5_OF  Hex MD5 of a numeric array's bytes, size included (the recipe of
%   catalog_identity, so leg 2's fingerprint droMd5 compares directly).
%   INPUTS: x numeric array.  OUTPUTS: h char [1 x 32].
md = java.security.MessageDigest.getInstance('MD5');
md.update(typecast([double(size(x)), double(x(:)).'], 'uint8'));
h = lower(reshape(dec2hex(typecast(md.digest(), 'uint8'), 2).', 1, []));
end

function g = git_state(gitBin, folder)
% GIT_STATE  Read-only git state of the repository holding a folder.
%   INPUTS:  gitBin char (git executable); folder char.
%   OUTPUTS: g struct .head .root .dirty .nDirty .tag .text.
[s1, h]   = system(sprintf('"%s" -C "%s" rev-parse --short HEAD', gitBin, folder));
[~,  top] = system(sprintf('"%s" -C "%s" rev-parse --show-toplevel', gitBin, folder));
[s3, st]  = system(sprintf('"%s" -C "%s" status --porcelain --untracked-files=no', gitBin, folder));
g = struct('head', 'unknown', 'root', strtrim(top), 'dirty', false, 'nDirty', NaN, ...
           'tag', 'unknown', 'text', '');
if s1 == 0, g.head = strtrim(h); end
if s3 == 0
    ln = strsplit(strtrim(st), newline);
    g.nDirty = nnz(~cellfun(@isempty, strtrim(ln)));
    g.dirty  = g.nDirty > 0;
end
g.tag = [g.head tern(g.dirty, '-dirty', '')];
if s3 ~= 0
    g.text = sprintf('%s (status unknown)', g.head);
elseif g.dirty
    g.text = sprintf('%s (DIRTY: %d tracked path(s) modified)', g.head, g.nDirty);
else
    g.text = sprintf('%s (clean)', g.head);
end
end
