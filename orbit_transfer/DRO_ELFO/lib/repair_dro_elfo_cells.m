function R = repair_dro_elfo_cells(sheetMat, mode, opts)
%% Purpose:
%
%   Re-solve cells of a DRO -> ELFO sheet from external seeds, through
%   thrust_ladder_library's seedTop / redo options. A copy of
%   GTO_DRO/lib/repair_gto_dro_cells.m with the parameter block changed and
%   two rules changed (below: ACCEPTANCE and SHEET KEY):
%
%   'unsolved'  a cell the engine retired without any OK rung (ATT >= maxAtt)
%               is re-solved ONCE, seeded from the nearest cell in the
%               departure-phase direction that holds a top-rung entry (phase
%               continuation).
%   'basin'     every OK top-rung entry gets the conjugate re-solve
%               (seed_from_z8 -> ms_tfmin conjTest; verdicts cached in the
%               sidecar keyed to the entry's z8); a REFUTED top rung whose
%               second rung is OK is re-solved from that second-rung entry
%               (reverse warm start).
%
%   ACCEPTANCE (astra disposition C5, GTO_DRO/reviews/astra_stage2_2026-11-
%   11_dispositions.md). The re-solve is accepted when its new top rung is
%   OK, passes the campaign gate (flown miss < meta.gateKm, acceptance
%   |dz| < accTol) and is faster than the stored top rung (any time, when
%   the cell had none). It is accepted EVEN WHEN the re-walked ladder lost
%   lower rungs: the PREVIOUS lower rungs are kept -- every rung below the
%   top that was OK before the repair keeps its previous entry (an extremal
%   of its own branch); a rung that was empty before takes the re-solve's
%   entry when that one is OK and passes the gate; any other rung is empty.
%   CLEARANCE is part of the acceptance: every re-solve entry that could be
%   stored (the new top rung, and a re-solve entry for an empty lower rung)
%   is flown and its LOCATED minimum altitude above the Moon and the Earth
%   (costate_common/located_min_altitude on the true min-time dynamics, the
%   routine body_clearance_sweep uses) must clear opts.floorMoonKm and
%   opts.floorEarthKm; a top rung that dips below a floor is 'restored
%   (clearance)', so a re-solve below the floor never replaces a stored
%   entry, and a lower-rung re-solve entry below a floor is not taken.
%   Otherwise the cell is RESTORED from its backup. The sidecar row records
%   which lower rungs kept their previous entry (.keptOld) and which the
%   re-solve lost (.lostNew).
%
%   SHEET KEY. The engine rewrites the whole sheet file at every save, so
%   the top-level sheet_key variable the front door appended is gone after
%   an engine call. The tool reads the key before anything else (from the
%   sheet; from a pending cell backup or the snapshot when an interrupted
%   repair already lost it), checks it against the sheet's meta (target =
%   arrFamily, depSize = depParams.tau, raanDeg = arrParams.raanDeg for an
%   ELFO arrival) and writes it back with every sheet save; a sheet with no
%   recoverable key, or a key that disagrees with its meta, is refused
%   (repair_dro_elfo_cells:sheetKey) before the sheet is touched.
%
%   DATA SAFETY (copied). The engine clears a redo cell on disk before it
%   solves, so before every engine call the cell's slices are written to
%   <sheet>_repair_backup_iDD_aAA.mat (atomically), deleted only after the
%   attempt's sidecar row is saved; at the start of every call pending
%   backups of this sheet are restored into it (outcome 'restored
%   (interrupted)'). Before the first engine call of a call the whole sheet
%   is copied to <sheet>.bak_repair beside it (the recovery source for an
%   unreadable sheet: no sheet is overwritten without a backup beside it).
%   The tool's own saves of the sheet and the sidecar are temp file + move.
%
%   DIRECT SIDECAR (copied, per rung). With engine option keepDirect the
%   engine's redo rewrites the cell's slices of <sheet>_direct.mat with the
%   re-solve's direct solutions; every rung whose stored entry is NOT the
%   re-solve's (a restored cell, or the lower rungs that kept their previous
%   entry) is marked not recorded (NaN, DDONE false).
%
%   PHYSICS GUARD (copied). opts.engineOpts must agree with the sheet's
%   meta on ispS, m0kg, N, floorKm, floorEarthKm, gateKm, depFamily,
%   arrFamily, tauDRO and every field of depParams / arrParams the sheet
%   records (the ELFO elements and raanDeg included); accTol (1e-6) and
%   thrLock (true) are campaign constants. A disagreement is refused
%   (repair_dro_elfo_cells:engineOpts) before the sheet is touched.
%
%   FLEET LOCK. The batched driver's lock is <sheet folder>/
%   .dro_elfo_driver.lock (run_dro_elfo_batched.sh: OUTDIR/.dro_elfo_
%   driver.lock). If it exists the call refuses
%   (repair_dro_elfo_cells:fleetRunning) unless opts.ignoreLock; a call that
%   will write the sheet TAKES it (mkdir) and releases it through onCleanup.
%   A conjOnly call may run under the lock unless a pending backup needs
%   recovery. A foreign lock is never removed.
%
%   Every engine call passes directGateKm = inf (no direct gate) and keeps
%   the sheet's gateKm; a cell is not touched twice in the same mode unless
%   opts.force.
%
%% Inputs:
%
%  sheetMat                 char                    Sheet .mat written by
%                                                   run_dro_elfo_torus
%                                                   (thrust_ladder_library +
%                                                   sheet_key)
%
%  mode                     char                    'unsolved' | 'basin'
%
%  opts                     struct                  .engineOpts (required:
%                                                   the sheet's engine
%                                                   options, e.g. the
%                                                   planOnly optsA of
%                                                   run_dro_elfo_torus;
%                                                   rungs, nD, nA, sD0, sA0,
%                                                   cells, redo, seedTop,
%                                                   tf0 are overwritten from
%                                                   the sheet)
%                                                   .conjOnly [false] (basin
%                                                   only: verdicts)
%                                                   .force [false]
%                                                   .seedOverride [] (a
%                                                   seedTop struct used for
%                                                   every repair; tests)
%                                                   .maxCells [inf]
%                                                   .logFile ['']
%                                                   .directGateKm [inf]
%                                                   .ignoreLock [false]
%                                                   .floorMoonKm
%                                                   [meta.floorKm, else 500]
%                                                   .floorEarthKm [200]
%                                                   ([] = no floor)
%                                                   (.engineFcn TEST SEAM
%                                                   [@thrust_ladder_library]
%                                                   (sheetMat, engineOpts))
%                                                   (.conjFcn TEST SEAM []
%                                                   (z8, rv0, rvf, ndp,
%                                                   meta) -> 1 | 0 | NaN;
%                                                   [] = the ms_tfmin
%                                                   conjugate re-solve)
%                                                   (.clearFcn TEST SEAM []
%                                                   (z8, rv0, rvf, phys) ->
%                                                   [1 x 2] located [Moon
%                                                   Earth] altitude km, NaN
%                                                   = not located; [] = fly
%                                                   + located_min_altitude)
%
%% Outputs:
%
%  R                        struct                  .mode .nCandidates
%                                                   .nAttempted .nRepaired
%                                                   .nRestored .nFailed
%                                                   (engine threw)
%                                                   .nRecovered .nKeptOld
%                                                   (lower rungs that kept
%                                                   their previous entry in
%                                                   accepted repairs)
%                                                   .nConjRefuted .conjPass
%                                                   [nD x nA] (NaN not run /
%                                                   0 / 1) .rows [1 x nRows]
%                                                   .sideMat
%
%% Revision History:
%  M. Casey                                                   (c) 11/19/2026
%  M. Casey  fix round 1: clearance gate before acceptance    (c) 11/19/2026
%  Copyright Coorbital Inc.
%% ------------------------ Begin Code Sequence ---------------------------

if nargin == 0
    % self-demo: verdicts only (conjOnly) on a 1 x 1 cut of the compare sheet
    here = fileparts(fileparts(mfilename('fullpath')));
    addpath(fullfile(here, 'tests'), fullfile(here, 'catalog'));
    sheetMat = fullfile(tempname, 'dro_elfo_tau1_raan000.mat');
    fixture_dro_elfo_sheet(sheetMat, struct('iD', 2, 'iA', 2, 'nR', 2));
    plan = run_dro_elfo_torus(1, struct('planOnly', true));
    mode = 'basin';
    opts = struct('engineOpts', plan.optsA, 'conjOnly', true);
end
assert(any(strcmp(mode, {'unsolved', 'basin'})), 'repair_dro_elfo_cells:mode', ...
       'mode must be unsolved or basin');
assert(isfield(opts, 'engineOpts'), 'repair_dro_elfo_cells:engineOpts', ...
       'opts.engineOpts is required');
  accTolCampaign = 1e-6;                     % campaign constants: the sheet's meta does
 thrLockCampaign = true;                     % not record them
        conjOnly = fdef(opts, 'conjOnly', false);
           force = fdef(opts, 'force', false);
         seedOvr = fdef(opts, 'seedOverride', []);
        maxCells = fdef(opts, 'maxCells', inf);
         logFile = fdef(opts, 'logFile', '');
         dirGate = fdef(opts, 'directGateKm', inf);
      ignoreLock = fdef(opts, 'ignoreLock', false);
       engineFcn = fdef(opts, 'engineFcn', @thrust_ladder_library);
         conjFcn = fdef(opts, 'conjFcn', []);
        clearFcn = fdef(opts, 'clearFcn', @entry_clearance);
assert(~(conjOnly && strcmp(mode, 'unsolved')), 'repair_dro_elfo_cells:mode', ...
       'conjOnly applies to mode basin only');
lg       = @(varargin) logmsg(logFile, sprintf(varargin{:}));
sideMat  = [sheetMat(1:end-4) '_repairs.mat'];
bkPrefix = [sheetMat(1:end-4) '_repair_backup_'];
snapMat  = [sheetMat '.bak_repair'];

pending = dir([bkPrefix 'i*_a*.mat']);
pending = pending(~cellfun('isempty', regexp({pending.name}, '_i\d\d_a\d\d\.mat$', 'once')));

Q = load_sheet_or_throw(sheetMat, '', snapMat, lg);
meta = Q.meta;  rungs = Q.rungs;  [nD, nA, nR] = size(Q.OK);

% the sheet key: sheet, else a pending backup, else the snapshot; checked against meta
key = recover_key(Q, pending, snapMat);
check_key(key, meta, sheetMat);

% engine options: physics must agree with the sheet; the sheet owns the grid
eo = check_engine_opts(opts.engineOpts, meta, accTolCampaign, thrLockCampaign);
for f = {'cells', 'redo', 'seedTop', 'tf0', 'maxCells', 'batchSec'}
    if isfield(eo, f{1}), eo = rmfield(eo, f{1}); end
end
eo.rungs = rungs;  eo.nD = nD;  eo.nA = nA;  eo.sD0 = meta.sD0;  eo.sA0 = meta.sA0;
eo.resume = true;  eo.logFile = logFile;  eo.directGateKm = dirGate;
accTol = accTolCampaign;
floorMoonKm = 500;
if isfield(meta, 'floorKm') && ~isempty(meta.floorKm), floorMoonKm = meta.floorKm; end
floorMoonKm  = nofloor(fdef(opts, 'floorMoonKm', floorMoonKm));
floorEarthKm = nofloor(fdef(opts, 'floorEarthKm', 200));
gateKm = meta.gateKm;                          % the campaign gate (engine call keeps it)
sDgrid = mod(meta.sD0 + (0:nD-1)/nD, 1);
assert(max(abs(sDgrid - Q.sD(:).')) < 1e-12, 'repair_dro_elfo_cells:grid', ...
       'the sheet departure grid is not the engine grid (sD0 %g, nD %d)', meta.sD0, nD);

% fleet lock
lockPath = fullfile(fileparts(sheetMat), '.dro_elfo_driver.lock');
writing = ~conjOnly || ~isempty(pending);
if ~ignoreLock && writing
    if isfolder(lockPath) || isfile(lockPath)
        error('repair_dro_elfo_cells:fleetRunning', ...
              'fleet driver lock %s exists (pending backups: %d); pass opts.ignoreLock to override', ...
              lockPath, numel(pending));
    end
    [lockMade, ~, lockId] = mkdir(lockPath);
    if ~lockMade || strcmp(lockId, 'MATLAB:MKDIR:DirectoryExists')
        error('repair_dro_elfo_cells:fleetRunning', 'could not take the fleet driver lock %s', lockPath);
    end
    lockRelease = onCleanup(@() rmdir(lockPath));
end

rows = struct('mode', {}, 'iD', {}, 'iA', {}, 'seedFrom', {}, 'tfBefore', {}, 'tfAfter', {}, ...
              'okAfter', {}, 'outcome', {}, 'conjBefore', {}, 'conjAfter', {}, 'date', {}, ...
              'wallS', {}, 'directGateKm', {}, 'keptOld', {}, 'lostNew', {}, 'clearKm', {});
cache = struct('z8', nan(8, nD, nA), 'pass', nan(nD, nA));
if isfile(sideMat)
    S = load(sideMat);
    rows = S.rows;
    if ~isempty(rows) && ~isfield(rows, 'clearKm'), [rows.clearKm] = deal(nan(1, 2)); end
    if isfield(S, 'cache') && isequal(size(S.cache.pass), [nD nA]), cache = S.cache; end
end

% recover interrupted repairs: restore the backed-up cell, record, delete
nRecovered = 0;
for kp = 1:numel(pending)
    bfile = fullfile(pending(kp).folder, pending(kp).name);
    B = load(bfile);
    Qr = load_sheet_or_throw(sheetMat, bfile, snapMat, lg);
    Qr = restore_cell(Qr, B.bk, B.iD, B.iA);
    Qr.sheet_key = key;
    save_atomic(sheetMat, Qr);
    clear_direct_rungs(sheetMat, B.iD, B.iA, true(1, nR), lg);
    row = B.row;  row.outcome = 'restored (interrupted)';  row.tfBefore = B.bk.TF(1);
    if ~isfield(row, 'clearKm'), row.clearKm = nan(1, 2); end      % backup of an older version
    row.tfAfter = NaN;  row.okAfter = false;
    rows(end+1) = row;
    save_atomic(sideMat, struct('rows', rows, 'cache', cache));
    delete(bfile);
    lg('  (%2d,%2d) recovered an interrupted %s repair from %s', B.iD, B.iA, B.mode, pending(kp).name);
    nRecovered = nRecovered + 1;
end
if nRecovered > 0, Q = load(sheetMat); end
if ~isfield(Q, 'sheet_key') && writing          % a keyless sheet gets its key back now
    Q.sheet_key = key;  save_atomic(sheetMat, Q);
end

already = false(nD, nA);
for kr = 1:numel(rows)
    if strcmp(rows(kr).mode, mode) && ~strcmp(rows(kr).outcome, 'restored (interrupted)')
        already(rows(kr).iD, rows(kr).iA) = true;
    end
end

ndp = nd_propulsion(rungs(1), meta.ispS, meta.m0kg, meta.lStar, meta.tStar);
[tD, rvD] = get_family_orbit(meta.depFamily, meta.depParams);
[tA, rvA] = get_family_orbit(meta.arrFamily, meta.arrParams);
depAt = phase_state(tD, rvD);  arrAt = phase_state(tA, rvA);
if isempty(conjFcn)
    conjFcn = @(z8, rv0, rvf, ndpX, metaX) conj_verdict(z8, rv0, rvf, ndpX, metaX, lg);
end

conjPass = nan(nD, nA);
cand = zeros(0, 2);  seedOf = {};
switch mode
    case 'unsolved'
        maxAtt = fdef(eo, 'maxAtt', 2);
        for iD = 1:nD
            for iA = 1:nA
                if any(Q.OK(iD, iA, :)) || Q.ATT(iD, iA) < maxAtt, continue, end
                if already(iD, iA) && ~force, continue, end
                nb = [];
                for step = 1:floor(nD/2)
                    for sgn = [-1 1]
                        jD = mod(iD - 1 + sgn*step, nD) + 1;
                        if Q.OK(jD, iA, 1), nb = jD; break, end
                    end
                    if ~isempty(nb), break, end
                end
                if isempty(nb)
                    lg('  (%2d,%2d) unsolved: no solved neighbour in the departure direction -- skipped', iD, iA);
                    continue
                end
                cand(end+1, :) = [iD iA];
                seedOf{end+1} = struct('iD', nb, 'iA', iA, 'kr', 1);
            end
        end
    case 'basin'
        cacheDirty = false;
        for iD = 1:nD
            for iA = 1:nA
                if ~Q.OK(iD, iA, 1), continue, end
                z8 = Q.Z8(:, iD, iA, 1);
                if ~force && ~isnan(cache.pass(iD, iA)) && isequal(cache.z8(:, iD, iA), z8)
                    conjPass(iD, iA) = cache.pass(iD, iA);
                else
                    conjPass(iD, iA) = conjFcn(z8, depAt(Q.sD(iD)), arrAt(Q.sA(iA)), ndp, meta);
                    cache.z8(:, iD, iA) = z8;  cache.pass(iD, iA) = conjPass(iD, iA);
                    cacheDirty = true;
                end
                if conjPass(iD, iA) == 0 && nR >= 2 && Q.OK(iD, iA, 2) && (~already(iD, iA) || force)
                    cand(end+1, :) = [iD iA];
                    seedOf{end+1} = struct('iD', iD, 'iA', iA, 'kr', 2);
                end
            end
        end
        if cacheDirty, save_atomic(sideMat, struct('rows', rows, 'cache', cache)); end
end
R = struct('mode', mode, 'nCandidates', size(cand, 1), 'nAttempted', 0, 'nRepaired', 0, ...
           'nRestored', 0, 'nFailed', 0, 'nRecovered', nRecovered, 'nKeptOld', 0, ...
           'nConjRefuted', nnz(conjPass == 0), 'conjPass', conjPass, 'rows', rows, 'sideMat', sideMat);
lg('repair %s: %d candidate cell(s), %d refuted top rungs', mode, R.nCandidates, R.nConjRefuted);
if conjOnly || isempty(cand), return, end

for kc = 1:min(size(cand, 1), maxCells)
    iD = cand(kc, 1);  iA = cand(kc, 2);  sf = seedOf{kc};
    t0 = tic;
    R.nAttempted = R.nAttempted + 1;
    if R.nAttempted == 1                              % snapshot before the first engine call
        tmpSnap = sprintf('%s_tmp%06d', snapMat, randi(999999));
        copyfile(sheetMat, tmpSnap);  movefile(tmpSnap, snapMat, 'f');
    end
    bk = backup_cell(Q, iD, iA);
    row = struct('mode', mode, 'iD', iD, 'iA', iA, 'seedFrom', sf, 'tfBefore', bk.TF(1), ...
                 'tfAfter', NaN, 'okAfter', false, 'outcome', '', 'conjBefore', conjPass(iD, iA), ...
                 'conjAfter', NaN, 'date', char(datetime('now', 'Format', 'yyyy-MM-dd HH:mm')), ...
                 'wallS', NaN, 'directGateKm', dirGate, 'keptOld', zeros(1, 0), 'lostNew', zeros(1, 0), ...
                 'clearKm', nan(1, 2));
    bfile = sprintf('%si%02d_a%02d.mat', bkPrefix, iD, iA);
    save_atomic(bfile, struct('bk', bk, 'mode', mode, 'iD', iD, 'iA', iA, 'seedFrom', sf, ...
                              'date', row.date, 'row', row, 'sheet_key', key));
    try
        if isempty(seedOvr)
            seedTop = seed_top_from_entry(sheetMat, sf.iD, sf.iA, sf.kr, rungs(1));
        else
            seedTop = seedOvr;
        end
        e = eo;  e.cells = [iD iA];  e.redo = true;  e.seedTop = seedTop;
        engineFcn(sheetMat, e);
        Qn = load(sheetMat);
        Qn.sheet_key = key;                           % the engine dropped it
        okR    = squeeze(Qn.OK(iD, iA, :)).';
        gateR  = okR & squeeze(Qn.FLYKM(iD, iA, :)).' < gateKm & squeeze(Qn.ACCDZ(iD, iA, :)).' < accTol;
        tfNew  = Qn.TF(iD, iA, 1);
        row.tfAfter = tfNew;  row.okAfter = okR(1);
        faster = okR(1) && (isnan(row.tfBefore) || tfNew < row.tfBefore);
        % clearance of every re-solve entry that could be stored (top rung,
        % and re-solve entries for lower rungs that were empty)
        clrR = true(1, nR);
        for kq = find(gateR & ([true, ~bk.OK(2:end).']))
            ndq  = nd_propulsion(rungs(kq), meta.ispS, meta.m0kg, meta.lStar, meta.tStar);
            phq  = struct('Tnd', ndq.Tnd, 'cnd', ndq.cnd, 'muStar', meta.muStar, 'lStar', meta.lStar, ...
                          'tStar', meta.tStar, 'm0kg', meta.m0kg);
            hq   = clearFcn(Qn.Z8(:, iD, iA, kq), depAt(Q.sD(iD)), arrAt(Q.sA(iA)), phq);
            clrR(kq) = hq(1) >= floorMoonKm && hq(2) >= floorEarthKm;       % NaN fails
            if kq == 1, row.clearKm = hq(:).'; end
        end
        if faster && gateR(1) && ~clrR(1)
            row.outcome = 'restored (clearance)';
            Qn = restore_cell(Qn, bk, iD, iA);
            save_atomic(sheetMat, Qn);  Q = Qn;
            clear_direct_rungs(sheetMat, iD, iA, true(1, nR), lg);
            R.nRestored = R.nRestored + 1;
        elseif faster && gateR(1)
            % C5: the new top rung; the previous lower rungs are kept
            [Qn, takeNew, keptOld] = merge_lower_rungs(Qn, bk, iD, iA, gateR & clrR);
            row.keptOld = find(keptOld);
            row.lostNew = find(bk.OK(:).' & ~okR);
            row.outcome = 'repaired';
            save_atomic(sheetMat, Qn);
            clear_direct_rungs(sheetMat, iD, iA, ~takeNew, lg);
            if strcmp(mode, 'basin')
                row.conjAfter = conjFcn(Qn.Z8(:, iD, iA, 1), depAt(Q.sD(iD)), arrAt(Q.sA(iA)), ndp, meta);
                cache.z8(:, iD, iA) = Qn.Z8(:, iD, iA, 1);  cache.pass(iD, iA) = row.conjAfter;
            end
            Q = Qn;  R.nRepaired = R.nRepaired + 1;  R.nKeptOld = R.nKeptOld + numel(row.keptOld);
        else
            if ~okR(1),         row.outcome = 'restored (repair failed)';
            elseif ~faster,     row.outcome = 'restored (not faster)';
            else,               row.outcome = 'restored (campaign gate)';
            end
            Qn = restore_cell(Qn, bk, iD, iA);
            save_atomic(sheetMat, Qn);  Q = Qn;
            clear_direct_rungs(sheetMat, iD, iA, true(1, nR), lg);
            R.nRestored = R.nRestored + 1;
        end
    catch ME
        Qn = load_sheet_or_throw(sheetMat, bfile, snapMat, lg);
        Qn = restore_cell(Qn, bk, iD, iA);
        Qn.sheet_key = key;
        save_atomic(sheetMat, Qn);  Q = Qn;
        clear_direct_rungs(sheetMat, iD, iA, true(1, nR), lg);
        row.outcome = ['restored (error: ' ME.message ')'];
        R.nFailed = R.nFailed + 1;  R.nRestored = R.nRestored + 1;
    end
    row.wallS = toc(t0);
    rows(end+1) = row;
    save_atomic(sideMat, struct('rows', rows, 'cache', cache));
    delete(bfile);                                   % only after the row is saved
    lg('  (%2d,%2d) %s: %s  tf %.5f -> %.5f  kept previous rungs [%s]  (%.0f s)', iD, iA, mode, ...
       row.outcome, row.tfBefore, row.tfAfter, num2str(row.keptOld), row.wallS);
end
R.rows = rows;
if nargin == 0, disp(R), end
end

% ------------------------------------------------------------------------
function [Qn, takeNew, keptOld] = merge_lower_rungs(Qn, bk, iD, iA, gateR)
% MERGE_LOWER_RUNGS  Astra C5: the re-solve's top rung stands; every lower
% rung that was OK before keeps its previous entry; an empty lower rung
% takes the re-solve's entry when it is OK and passes the campaign gate;
% any other rung is emptied (OK = false, slices NaN).
% INPUTS: Qn [1 x 1] sheet struct after the engine call; bk [1 x 1] cell
%         backup (fields [nR x 1], Z8 [8 x nR]); iD, iA [1 x 1] int; gateR
%         [1 x nR] logical (re-solve OK, gated and clear per rung).
% OUTPUTS: Qn [1 x 1] merged sheet struct; takeNew [1 x nR] logical (the
%          stored entry is the re-solve's); keptOld [1 x nR] logical.
nR = numel(bk.TF);
takeNew = false(1, nR);  keptOld = false(1, nR);
takeNew(1) = true;
for kr = 2:nR
    if bk.OK(kr)
        keptOld(kr) = true;
        Qn.TF(iD, iA, kr)    = bk.TF(kr);     Qn.FLYKM(iD, iA, kr) = bk.FLYKM(kr);
        Qn.ACCDZ(iD, iA, kr) = bk.ACCDZ(kr);  Qn.RES(iD, iA, kr)   = bk.RES(kr);
        Qn.WALL(iD, iA, kr)  = bk.WALL(kr);   Qn.OK(iD, iA, kr)    = true;
        Qn.Z8(:, iD, iA, kr) = bk.Z8(:, kr);
    elseif gateR(kr)
        takeNew(kr) = true;
    else
        Qn.TF(iD, iA, kr) = NaN;     Qn.FLYKM(iD, iA, kr) = NaN;  Qn.ACCDZ(iD, iA, kr) = NaN;
        Qn.RES(iD, iA, kr) = NaN;    Qn.WALL(iD, iA, kr) = NaN;   Qn.OK(iD, iA, kr) = false;
        Qn.Z8(:, iD, iA, kr) = NaN;
    end
end
end

function key = recover_key(Q, pending, snapMat)
% RECOVER_KEY  The sheet_key from the sheet, else from a pending cell
% backup, else from the whole-sheet snapshot; refused when none has it.
% INPUTS: Q [1 x 1] sheet struct; pending [p x 1] dir() struct of cell
% backups; snapMat [1 x n] char.  OUTPUTS: key [1 x 1] struct (target,
% depSize, raanDeg).
if isfield(Q, 'sheet_key'), key = Q.sheet_key; return, end
for kp = 1:numel(pending)
    B = load(fullfile(pending(kp).folder, pending(kp).name));
    if isfield(B, 'sheet_key'), key = B.sheet_key; return, end
end
if isfile(snapMat)
    vars = who('-file', snapMat);
    if ismember('sheet_key', vars)
        P = load(snapMat, '-mat', 'sheet_key');  key = P.sheet_key;  return
    end
end
error('repair_dro_elfo_cells:sheetKey', ...
      'the sheet carries no sheet_key and none is recoverable (cell backups, %s) -- refusing', snapMat);
end

function check_key(key, meta, sheetMat)
% CHECK_KEY  Refuse a sheet_key that disagrees with the sheet's meta.
% INPUTS: key [1 x 1] struct; meta [1 x 1] struct; sheetMat [1 x n] char.
% OUTPUTS: none (throws repair_dro_elfo_cells:sheetKey).
ok = all(isfield(key, {'target', 'depSize', 'raanDeg'})) && strcmp(key.target, meta.arrFamily) && ...
     isfield(meta.depParams, 'tau') && abs(key.depSize - meta.depParams.tau) < 1e-12;
if ok && strcmp(meta.arrFamily, 'elfo')
    ok = isfield(meta.arrParams, 'raanDeg') && ...
         abs(mod(key.raanDeg - meta.arrParams.raanDeg + 180, 360) - 180) < 1e-9;
end
if ~ok
    error('repair_dro_elfo_cells:sheetKey', 'sheet %s: sheet_key disagrees with the sheet''s meta -- refusing', sheetMat);
end
end

function eo = check_engine_opts(eo, meta, accTolC, thrLockC)
% CHECK_ENGINE_OPTS  Refuse engine options that disagree with the sheet's
% meta; fill omitted ones from meta; accTol / thrLock must equal the
% campaign constants.  INPUTS: eo [1 x 1] struct; meta [1 x 1] struct;
% accTolC, thrLockC [1 x 1].  OUTPUTS: eo [1 x 1] struct.
cc = {'accTol', accTolC; 'thrLock', thrLockC};
for kk = 1:2
    k = cc{kk, 1};
    if ~isfield(eo, k)
        eo.(k) = cc{kk, 2};
    elseif ~isequal(eo.(k), cc{kk, 2})
        error('repair_dro_elfo_cells:engineOpts', 'engineOpts.%s = %s but the campaign constant is %s', ...
              k, mat2str_any(eo.(k)), mat2str_any(cc{kk, 2}));
    end
end
for f = {'ispS', 'm0kg', 'N', 'floorKm', 'floorEarthKm', 'gateKm', 'depFamily', 'arrFamily', 'tauDRO'}
    k = f{1};
    if ~isfield(meta, k), continue, end
    if ~isfield(eo, k)
        eo.(k) = meta.(k);
    elseif ~isequal(eo.(k), meta.(k))
        error('repair_dro_elfo_cells:engineOpts', 'engineOpts.%s = %s but the sheet has %s', ...
              k, mat2str_any(eo.(k)), mat2str_any(meta.(k)));
    end
end
for f = {'depParams', 'arrParams'}
    k = f{1};
    if ~isfield(meta, k), continue, end
    if isfield(eo, k)
        for g = fieldnames(eo.(k)).'
            if isfield(meta.(k), g{1}) && ~isequal(eo.(k).(g{1}), meta.(k).(g{1}))
                error('repair_dro_elfo_cells:engineOpts', 'engineOpts.%s.%s = %s but the sheet has %s', ...
                      k, g{1}, mat2str_any(eo.(k).(g{1})), mat2str_any(meta.(k).(g{1})));
            end
        end
    end
    eo.(k) = meta.(k);                         % the sheet's own parameters
end
end

function s = mat2str_any(v)
% MAT2STR_ANY  Printable form of an option.  INPUTS: v [any size] numeric,
% logical or char.  OUTPUTS: s [1 x n] char.
if ischar(v), s = v; elseif isnumeric(v) || islogical(v), s = mat2str(v); else, s = class(v); end
end

function Q = load_sheet_or_throw(sheetMat, bfile, snapMat, lg)
% LOAD_SHEET_OR_THROW  Load the sheet; on failure keep the backup and throw
% repair_dro_elfo_cells:sheetUnreadable naming the snapshot.
% INPUTS: sheetMat, bfile ('' at call start), snapMat [1 x n] char; lg
% [1 x 1] function handle.  OUTPUTS: Q [1 x 1] sheet struct.
try
    Q = load(sheetMat);
    assert(isfield(Q, 'OK') && isfield(Q, 'Z8'), 'sheet lacks OK/Z8');
catch ME
    lg('sheet %s unreadable (%s); cell backup %s kept; whole-sheet recovery source %s', ...
       sheetMat, ME.message, bfile, snapMat);
    error('repair_dro_elfo_cells:sheetUnreadable', ...
          'sheet %s unreadable (%s); cell backup "%s" kept; recover from the snapshot %s', ...
          sheetMat, ME.message, bfile, snapMat);
end
end

function save_atomic(path, S)
% SAVE_ATOMIC  Save the fields of struct S to path via temp file + movefile.
% INPUTS: path [1 x n] char; S [1 x 1] struct.  OUTPUTS: none (file written).
tmp = sprintf('%s_tmp%06d.mat', path(1:end-4), randi(999999));
save(tmp, '-struct', 'S');
movefile(tmp, path, 'f');
end

function cp = conj_verdict(z8, rv0, rvf, ndp, meta, lg)
% CONJ_VERDICT  The S4 junction test on one entry: NaN when the re-solve
% does not reproduce the entry (copied from repair_gto_dro_cells).
% INPUTS: z8 [8 x 1]; rv0, rvf [6 x 1]; ndp [1 x 1] struct (.Tnd .cnd);
% meta [1 x 1] struct; lg [1 x 1] handle.  OUTPUTS: cp [1 x 1] 1 / 0 / NaN.
cp = NaN;
rv0 = rv0(1:6);  rv0 = rv0(:).';  rvf = rvf(1:6);  rvf = rvf(:).';
try
    seed = seed_from_z8(z8, rv0, 24, ndp.Tnd, ndp.cnd, meta.muStar);
    [z, info] = ms_tfmin(rv0, rvf, seed, ndp.Tnd, ndp.cnd, meta.muStar, ...
                         struct('conjTest', true, 'wallSec', 120));
    if info.converged && max(abs(z - z8)) < 1e-6
        cp = double(info.conj.pass);
    else
        lg('  conjugate re-solve did not reproduce the entry (|dz| %.1e) -- verdict NaN', max(abs(z - z8)));
    end
catch ME
    lg('  conjugate re-solve error: %s -- verdict NaN', ME.message);
end
end

function seedTop = seed_top_from_entry(sheetMat, iD, iA, kr, thrustTargetN)
% SEED_TOP_FROM_ENTRY  Fly a sheet entry and resample it on the engine's
% N + 1 nodes as an opts.seedTop (thr = 1, t_f scaled by (T/Ttarget)^0.6);
% a copy of GTO_DRO/lib/seed_top_from_entry.m, local so no other
% campaign's folder goes on the path.
% INPUTS: sheetMat [1 x n] char; iD, iA, kr [1 x 1] int; thrustTargetN
% [1 x 1] N.  OUTPUTS: seedTop [1 x 1] struct .X [7 x N+1], .U [4 x N+1],
% .tf [1 x 1], .from [1 x 1] struct.
Q = load(sheetMat, 'Z8', 'OK', 'rungs', 'sD', 'meta');
if ~Q.OK(iD, iA, kr)
    error('repair_dro_elfo_cells:seedNotOK', 'seed entry (%d,%d,rung %g N) of %s is not OK', ...
          iD, iA, Q.rungs(kr), sheetMat);
end
meta = Q.meta;
z8   = Q.Z8(:, iD, iA, kr);
TN   = Q.rungs(kr);
ndp  = nd_propulsion(TN, meta.ispS, meta.m0kg, meta.lStar, meta.tStar);
[tD, rvD] = get_family_orbit(meta.depFamily, meta.depParams);
dep  = phase_state(tD, rvD);
rv0  = dep(Q.sD(iD));
[tt, yy] = pumpkyn.cr3bp.tfMinProp(z8(8), [rv0(:); 1; z8(1:7)], ndp.Tnd, ndp.cnd, meta.muStar);
N    = meta.N;
tg   = linspace(0, z8(8), N+1);
[tu, iu] = unique(tt(:), 'stable');
if numel(tu) < numel(tt), tt = tu;  yy = yy(iu, :); end
Y    = interp1(tt, yy, tg, 'spline');
if any(~isfinite(Y(:)))
    error('repair_dro_elfo_cells:seedNonFinite', 'resampled flight of entry (%d,%d,%g N) is not finite', iD, iA, TN);
end
lamV  = Y(:, 11:13).';
alpha = -lamV ./ max(vecnorm(lamV, 2, 1), eps);
seedTop = struct('X', Y(:, 1:7).', 'U', [alpha; ones(1, N+1)], ...
                 'tf', z8(8)*(TN/thrustTargetN)^0.6, ...
                 'from', struct('iD', iD, 'iA', iA, 'kr', kr, 'thrustN', TN, 'tf_nd', z8(8)));
end

function bk = backup_cell(Q, iD, iA)
% BACKUP_CELL  Copy one cell's slices.  INPUTS: Q [1 x 1] sheet struct;
% iD, iA [1 x 1] int.  OUTPUTS: bk [1 x 1] struct (TF..OK [nR x 1], Z8
% [8 x nR], ATT [1 x 1]).
bk = struct('TF', squeeze(Q.TF(iD, iA, :)), 'FLYKM', squeeze(Q.FLYKM(iD, iA, :)), ...
            'ACCDZ', squeeze(Q.ACCDZ(iD, iA, :)), 'RES', squeeze(Q.RES(iD, iA, :)), ...
            'WALL', squeeze(Q.WALL(iD, iA, :)), 'OK', squeeze(Q.OK(iD, iA, :)), ...
            'Z8', reshape(Q.Z8(:, iD, iA, :), 8, []), 'ATT', Q.ATT(iD, iA));
end

function Q = restore_cell(Q, bk, iD, iA)
% RESTORE_CELL  Put a backed-up cell back.  INPUTS: Q [1 x 1] sheet struct;
% bk [1 x 1] cell backup; iD, iA [1 x 1] int.  OUTPUTS: Q [1 x 1] struct.
nR = numel(bk.TF);
Q.TF(iD, iA, :)    = reshape(bk.TF, 1, 1, nR);
Q.FLYKM(iD, iA, :) = reshape(bk.FLYKM, 1, 1, nR);
Q.ACCDZ(iD, iA, :) = reshape(bk.ACCDZ, 1, 1, nR);
Q.RES(iD, iA, :)   = reshape(bk.RES, 1, 1, nR);
Q.WALL(iD, iA, :)  = reshape(bk.WALL, 1, 1, nR);
Q.OK(iD, iA, :)    = reshape(bk.OK, 1, 1, nR);
Q.Z8(:, iD, iA, :) = reshape(bk.Z8, 8, 1, 1, nR);
Q.ATT(iD, iA)      = bk.ATT;
end

function clear_direct_rungs(sheetMat, iD, iA, mask, lg)
% CLEAR_DIRECT_RUNGS  Mark the masked rungs of one cell of <sheet>_direct.mat
% not recorded (NaN, DDONE false); no sidecar or an empty mask = nothing to
% do. A sidecar that cannot be read or written is logged, never thrown.
% INPUTS: sheetMat [1 x n] char; iD, iA [1 x 1] int; mask [1 x nR]
% logical; lg [1 x 1] function handle.  OUTPUTS: none (file written).
dirMat = [sheetMat(1:end-4) '_direct.mat'];
if ~isfile(dirMat) || ~any(mask), return, end
try
    D = load(dirMat);
    if size(D.DDONE, 1) < iD || size(D.DDONE, 2) < iA || size(D.DDONE, 3) ~= numel(mask)
        lg('  (%2d,%2d) direct sidecar %s does not fit the sheet, left as is', iD, iA, dirMat);
        return
    end
    kk = find(mask);
    D.DX(:, :, iD, iA, kk) = NaN;   D.DU(:, :, iD, iA, kk) = NaN;
    D.DTF(iD, iA, kk) = NaN;        D.DEARTHKM(iD, iA, kk) = NaN;
    D.DMOONKM(iD, iA, kk) = NaN;    D.DDONE(iD, iA, kk) = false;
    save_atomic(dirMat, D);
    lg('  (%2d,%2d) direct-sidecar slices of rung(s) [%s] marked not recorded', iD, iA, num2str(kk));
catch ME
    lg('  (%2d,%2d) the direct sidecar %s could not be cleared: %s', iD, iA, dirMat, ME.message);
end
end

function h = entry_clearance(z8, rv0, rvf, phys)
% ENTRY_CLEARANCE  Located minimum altitude of one entry's flight above the
% Moon and the Earth (fly_transfer + located_min_altitude on the true
% min-time dynamics, as body_clearance_sweep measures it); NaN when the
% flight throws, is non-finite, or a minimum is not located.
% INPUTS: z8 [8 x 1]; rv0, rvf [6 x 1]; phys [1 x 1] struct (.Tnd .cnd
% .muStar .lStar .tStar .m0kg).  OUTPUTS: h [1 x 2] [Moon Earth] km.
h = [NaN NaN];
try
    F = fly_transfer(z8, rv0(:), rvf(:), phys);
    if ~all(isfinite(F.Y(:, 1:6)), 'all'), return, end
    lo = struct('lStar', phys.lStar, 'phys', phys);
    Cm = located_min_altitude(F.t(:), F.Y, [1 - phys.muStar 0 0], 1737.4, -inf, lo);
    Ce = located_min_altitude(F.t(:), F.Y, [-phys.muStar 0 0], 6378.137, -inf, lo);
    if Cm.nFailed == 0, h(1) = Cm.hMinKm; end
    if Ce.nFailed == 0, h(2) = Ce.hMinKm; end
catch
    h = [NaN NaN];
end
end

function f = nofloor(f)
% NOFLOOR  [] means no floor.  INPUTS: f [] or [1 x 1] km.  OUTPUTS: f
% [1 x 1] km (-inf for none).
if isempty(f), f = -inf; end
end

function v = fdef(s, f, v0)
% FDEF  s.(f) if present else v0.  INPUTS: s [1 x 1] struct; f [1 x n]
% char; v0 [any].  OUTPUTS: v [any].
if isfield(s, f), v = s.(f); else, v = v0; end
end

function logmsg(f, s)
% LOGMSG  Append-mode file log, stdout when f is empty.  INPUTS: f [1 x n]
% char log path ('' = stdout); s [1 x m] char.  OUTPUTS: none.
if isempty(f), fprintf('%s\n', s);
else, fid = fopen(f, 'a'); fprintf(fid, '%s\n', s); fclose(fid); end
end
