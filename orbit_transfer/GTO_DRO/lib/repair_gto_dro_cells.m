function R = repair_gto_dro_cells(sheetMat, mode, opts)
%% Purpose:
%
%   Re-solve cells of a GTO -> DRO sheet from external seeds, through
%   thrust_ladder_library's seedTop / redo options. Two repairs the fleet
%   cannot do for itself (stage 1 findings, GTO_DRO/README.md):
%
%   'unsolved'  a cell the engine retired without any OK rung (the perigee
%               column at 15 N burns its attempts at the iteration cap) is
%               re-solved ONCE, seeded from the nearest solved cell in the
%               departure-phase direction (phase continuation).
%   'basin'     every OK top-rung entry gets the conjugate re-solve
%               (seed_from_z8 -> ms_tfmin conjTest, the stage 1 S4 junction
%               test; verdicts are cached in the sidecar keyed to the
%               entry's z8, so an unchanged entry is not re-tested); a
%               refuted entry whose second rung is OK is re-solved from that
%               second-rung entry (reverse warm start: stage 1 measured
%               15 N <- 12 N landing on the faster basin). The repair is
%               kept only when the new top rung is OK and faster, no rung
%               that was OK is lost, and every OK rung passes the campaign
%               gate (flown miss < meta.gateKm, acceptance < accTol).
%               Otherwise the cell is RESTORED from its backup.
%
%   DATA SAFETY. The engine clears a redo cell on disk before it solves, so
%   before every engine call the cell's original slices are written to
%   <sheet>_repair_backup_iDD_aAA.mat (atomically). That file is deleted
%   only after the final sidecar row of the attempt is saved. At the START
%   of every call, pending backups of this sheet are restored into it
%   (outcome 'restored (interrupted)', R.nRecovered), which covers a kill,
%   a MUMPS segfault or Ctrl-C mid-repair. The tool's own saves of the
%   sheet and the sidecar are temp-file + movefile. If the sheet cannot be
%   loaded for the restore the backup is left in place and the call throws
%   repair_gto_dro_cells:sheetUnreadable.
%
%   PHYSICS GUARD. opts.engineOpts must agree with the sheet's meta on
%   ispS, m0kg, N, floorKm, floorEarthKm, gateKm, depFamily, arrFamily,
%   tauDRO and every field of depParams / arrParams the sheet records;
%   a disagreement is refused (repair_gto_dro_cells:engineOpts) before the
%   sheet is touched, an omitted field is filled from meta. accTol (1e-6)
%   and thrLock (true) are CAMPAIGN CONSTANTS of this tool, defined in the
%   knob block, because the sheet's meta does not record them: engineOpts
%   carrying a different value is refused, an omitted one is set, and the
%   tool's own acceptance gate uses the constant. Solver budgets
%   (maxCpuSec, maxIter, ...) are free.
%
%   FLEET LOCK. If <results>/.gtodrocat_driver.lock (the batched driver's
%   lock; results/ = parent of the sheet's folder) exists, the call refuses
%   (repair_gto_dro_cells:fleetRunning) unless opts.ignoreLock. A conjOnly
%   call may run under the lock (it writes only the sidecar) unless a
%   pending backup needs recovery. Any call that will write the sheet
%   TAKES the lock itself (mkdir, the same atomic acquisition the batched
%   driver uses, so the driver refuses to start during a repair) and
%   releases it through onCleanup, also on error and Ctrl-C. With
%   opts.ignoreLock the tool neither checks nor takes it. A foreign lock
%   is never removed.
%
%   SNAPSHOT. Before the first engine call of a call, the whole sheet is
%   copied to <sheet>.bak_repair (previous snapshot overwritten, copy to a
%   temp name then movefile) and stays on disk: the engine's own sheet
%   saves are not atomic, and the per-cell backups cannot rebuild other
%   cells, so this snapshot is the recovery source for an unreadable
%   sheet (repair_gto_dro_cells:sheetUnreadable names it).
%
%   Every engine call passes directGateKm = 1000 (the direct control of
%   the faster basin flies ~194 km off at N = 400; refinement closes it)
%   and leaves gateKm at the sheet's own value; this tool then enforces the
%   campaign gate itself. A refused or throwing engine call ends in a
%   restored cell and a sidecar row. A cell is not touched twice in the
%   same mode unless opts.force.
%
%% Inputs:
%
%  sheetMat                 char                    sheet .mat written by
%                                                   thrust_ladder_library
%
%  mode                     char                    'unsolved' | 'basin'
%
%  opts                     struct                  .engineOpts (required:
%                                                   the sheet's engine
%                                                   options, e.g. planOnly
%                                                   optsA of
%                                                   run_gto_dro_catalog;
%                                                   rungs, nD, nA, sD0, sA0,
%                                                   cells, redo, seedTop,
%                                                   tf0 are overwritten from
%                                                   the sheet),
%                                                   .conjOnly [false] (basin
%                                                   only: verdicts),
%                                                   .force [false],
%                                                   .seedOverride [] (a
%                                                   seedTop struct used for
%                                                   every repair; tests),
%                                                   .maxCells [inf],
%                                                   .logFile [''],
%                                                   .directGateKm [1000],
%                                                   .ignoreLock [false]
%
%% Outputs:
%
%  R                        struct                  .mode, .nCandidates,
%                                                   .nAttempted (cells
%                                                   tried under maxCells),
%                                                   .nRepaired, .nRestored,
%                                                   .nFailed (engine threw),
%                                                   .nRecovered (interrupted
%                                                   repairs restored at the
%                                                   start), .nConjRefuted,
%                                                   .conjPass [nD x nA]
%                                                   (NaN not run / 0 / 1),
%                                                   .rows [1 x nRows],
%                                                   .sideMat
%
%% Revision History:
%  M. Casey                                                   (c) 11/10/2026
%  Copyright Coorbital Inc.
%% ------------------------ Begin Code Sequence ---------------------------

if nargin == 0
    here = fileparts(fileparts(mfilename('fullpath')));
    addpath(fullfile(here, 'tests'), fullfile(here, 'catalog'));
    sheetMat = fixture_demo_sheet();  mode = 'basin';
    plan = run_gto_dro_catalog(2, [], [], [15 12], struct('planOnly', true));
    opts = struct('engineOpts', plan.optsA, 'conjOnly', true);
end
assert(any(strcmp(mode, {'unsolved', 'basin'})), 'repair_gto_dro_cells:mode', ...
       'mode must be unsolved or basin');
assert(isfield(opts, 'engineOpts'), 'repair_gto_dro_cells:engineOpts', ...
       'opts.engineOpts is required');
accTolCampaign  = 1e-6;                       % campaign constants: the sheet's meta does
thrLockCampaign = true;                       % not record them
conjOnly   = fdef(opts, 'conjOnly', false);
force      = fdef(opts, 'force', false);
seedOvr    = fdef(opts, 'seedOverride', []);
maxCells   = fdef(opts, 'maxCells', inf);
logFile    = fdef(opts, 'logFile', '');
dirGate    = fdef(opts, 'directGateKm', 1000);
ignoreLock = fdef(opts, 'ignoreLock', false);
assert(~(conjOnly && strcmp(mode, 'unsolved')), 'repair_gto_dro_cells:mode', ...
       'conjOnly applies to mode basin only');
lg       = @(varargin) logmsg(logFile, sprintf(varargin{:}));
sideMat  = [sheetMat(1:end-4) '_repairs.mat'];
bkPrefix = [sheetMat(1:end-4) '_repair_backup_'];

snapMat = [sheetMat '.bak_repair'];
Q = load_sheet_or_throw(sheetMat, '', snapMat, lg);
meta = Q.meta;  rungs = Q.rungs;  [nD, nA, nR] = size(Q.OK);

% engine options: physics must agree with the sheet; the sheet owns the grid
eo = check_engine_opts(opts.engineOpts, meta, accTolCampaign, thrLockCampaign);
for f = {'cells', 'redo', 'seedTop', 'tf0', 'maxCells', 'batchSec'}
    if isfield(eo, f{1}), eo = rmfield(eo, f{1}); end
end
eo.rungs = rungs;  eo.nD = nD;  eo.nA = nA;  eo.sD0 = meta.sD0;  eo.sA0 = meta.sA0;
eo.resume = true;  eo.logFile = logFile;  eo.directGateKm = dirGate;
accTol = accTolCampaign;
gateKm = meta.gateKm;                          % the campaign gate (engine call keeps it)
sDgrid = mod(meta.sD0 + (0:nD-1)/nD, 1);
assert(max(abs(sDgrid - Q.sD(:).')) < 1e-12, 'repair_gto_dro_cells:grid', ...
       'the sheet departure grid is not the engine grid (sD0 %g, nD %d)', meta.sD0, nD);

% fleet lock and pending backups
pending = dir([bkPrefix 'i*_a*.mat']);
pending = pending(~cellfun('isempty', regexp({pending.name}, '_i\d\d_a\d\d\.mat$', 'once')));
lockPath = fullfile(fileparts(fileparts(sheetMat)), '.gtodrocat_driver.lock');
writing = ~conjOnly || ~isempty(pending);
if ~ignoreLock && writing
    if isfolder(lockPath) || isfile(lockPath)
        error('repair_gto_dro_cells:fleetRunning', ...
              'fleet driver lock %s exists (pending backups: %d); pass opts.ignoreLock to override', ...
              lockPath, numel(pending));
    end
    [lockMade, ~, lockId] = mkdir(lockPath);
    if ~lockMade || strcmp(lockId, 'MATLAB:MKDIR:DirectoryExists')
        error('repair_gto_dro_cells:fleetRunning', 'could not take the fleet driver lock %s', lockPath);
    end
    lockRelease = onCleanup(@() rmdir(lockPath));
end

rows = struct('mode', {}, 'iD', {}, 'iA', {}, 'seedFrom', {}, 'tfBefore', {}, 'tfAfter', {}, ...
              'okAfter', {}, 'outcome', {}, 'conjBefore', {}, 'conjAfter', {}, 'date', {}, ...
              'wallS', {}, 'directGateKm', {});
cache = struct('z8', nan(8, nD, nA), 'pass', nan(nD, nA));
if isfile(sideMat)
    S = load(sideMat);
    rows = S.rows;
    if isfield(S, 'cache') && isequal(size(S.cache.pass), [nD nA]), cache = S.cache; end
end

% recover interrupted repairs: restore the backed-up cell, record, delete
nRecovered = 0;
for kp = 1:numel(pending)
    bfile = fullfile(pending(kp).folder, pending(kp).name);
    B = load(bfile);
    Qr = load_sheet_or_throw(sheetMat, bfile, snapMat, lg);
    Qr = restore_cell(Qr, B.bk, B.iD, B.iA);
    save_atomic(sheetMat, Qr);
    row = B.row;  row.outcome = 'restored (interrupted)';  row.tfBefore = B.bk.TF(1);
    row.tfAfter = NaN;  row.okAfter = false;
    rows(end+1) = row;
    save_atomic(sideMat, struct('rows', rows, 'cache', cache));
    delete(bfile);
    lg('  (%2d,%2d) recovered an interrupted %s repair from %s', B.iD, B.iA, B.mode, pending(kp).name);
    nRecovered = nRecovered + 1;
end
if nRecovered > 0, Q = load(sheetMat); end

already = false(nD, nA);
for kr = 1:numel(rows)
    if strcmp(rows(kr).mode, mode) && ~strcmp(rows(kr).outcome, 'restored (interrupted)'), already(rows(kr).iD, rows(kr).iA) = true; end
end

ndp = nd_propulsion(rungs(1), meta.ispS, meta.m0kg, meta.lStar, meta.tStar);
[tD, rvD] = get_family_orbit(meta.depFamily, meta.depParams);
[tA, rvA] = get_family_orbit(meta.arrFamily, meta.arrParams);
depAt = phase_state(tD, rvD);  arrAt = phase_state(tA, rvA);

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
                    conjPass(iD, iA) = conj_verdict(z8, depAt(Q.sD(iD)), arrAt(Q.sA(iA)), ndp, meta, lg, iD, iA);
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
           'nRestored', 0, 'nFailed', 0, 'nRecovered', nRecovered, 'nConjRefuted', nnz(conjPass == 0), ...
           'conjPass', conjPass, 'rows', rows, 'sideMat', sideMat);
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
                 'wallS', NaN, 'directGateKm', dirGate);
    bfile = sprintf('%si%02d_a%02d.mat', bkPrefix, iD, iA);
    save_atomic(bfile, struct('bk', bk, 'mode', mode, 'iD', iD, 'iA', iA, 'seedFrom', sf, ...
                              'date', row.date, 'row', row));
    try
        if isempty(seedOvr)
            seedTop = seed_top_from_entry(sheetMat, sf.iD, sf.iA, sf.kr, rungs(1));
        else
            seedTop = seedOvr;
        end
        e = eo;  e.cells = [iD iA];  e.redo = true;  e.seedTop = seedTop;
        thrust_ladder_library(sheetMat, e);
        Qn = load(sheetMat);
        tfNew = Qn.TF(iD, iA, 1);  okNew = Qn.OK(iD, iA, 1);
        okR = squeeze(Qn.OK(iD, iA, :));
        gateOK = all(squeeze(Qn.FLYKM(iD, iA, okR)) < gateKm) && ...
                 all(squeeze(Qn.ACCDZ(iD, iA, okR)) < accTol);
        row.tfAfter = tfNew;  row.okAfter = okNew;
        faster = okNew && (isnan(row.tfBefore) || tfNew < row.tfBefore);
        lostRungs = strcmp(mode, 'basin') && any(bk.OK(:) & ~okR(:));
        if faster && ~lostRungs && gateOK
            row.outcome = 'repaired';
            if strcmp(mode, 'basin')
                row.conjAfter = conj_verdict(Qn.Z8(:, iD, iA, 1), depAt(Q.sD(iD)), ...
                                             arrAt(Q.sA(iA)), ndp, meta, lg, iD, iA);
                cache.z8(:, iD, iA) = Qn.Z8(:, iD, iA, 1);  cache.pass(iD, iA) = row.conjAfter;
            end
            Q = Qn;  R.nRepaired = R.nRepaired + 1;
        else
            if ~okNew,          row.outcome = 'restored (repair failed)';
            elseif ~faster,     row.outcome = 'restored (not faster)';
            elseif lostRungs,   row.outcome = 'restored (lost rungs)';
            else,               row.outcome = 'restored (campaign gate)';
            end
            Qn = restore_cell(Qn, bk, iD, iA);
            save_atomic(sheetMat, Qn);  Q = Qn;
            R.nRestored = R.nRestored + 1;
        end
    catch ME
        Qn = load_sheet_or_throw(sheetMat, bfile, snapMat, lg);
        Qn = restore_cell(Qn, bk, iD, iA);
        save_atomic(sheetMat, Qn);  Q = Qn;
        row.outcome = ['restored (error: ' ME.message ')'];
        R.nFailed = R.nFailed + 1;  R.nRestored = R.nRestored + 1;
    end
    row.wallS = toc(t0);
    rows(end+1) = row;
    save_atomic(sideMat, struct('rows', rows, 'cache', cache));
    delete(bfile);                                   % only after the row is saved
    lg('  (%2d,%2d) %s: %s  tf %.5f -> %.5f  (%.0f s)', iD, iA, mode, row.outcome, ...
       row.tfBefore, row.tfAfter, row.wallS);
end
R.rows = rows;
if nargin == 0, disp(R), end
end

% ------------------------------------------------------------------------
function eo = check_engine_opts(eo, meta, accTolC, thrLockC)
% CHECK_ENGINE_OPTS  Refuse engine options that disagree with the sheet's
% meta; fill omitted ones from meta; accTol / thrLock must equal the
% campaign constants.  INPUTS: eo; meta; accTolC, thrLockC scalars.  OUTPUTS: eo.
cc = {'accTol', accTolC; 'thrLock', thrLockC};
for kk = 1:2
    k = cc{kk, 1};
    if ~isfield(eo, k)
        eo.(k) = cc{kk, 2};
    elseif ~isequal(eo.(k), cc{kk, 2})
        error('repair_gto_dro_cells:engineOpts', 'engineOpts.%s = %s but the campaign constant is %s', ...
              k, mat2str_any(eo.(k)), mat2str_any(cc{kk, 2}));
    end
end
for f = {'ispS', 'm0kg', 'N', 'floorKm', 'floorEarthKm', 'gateKm', 'depFamily', 'arrFamily', 'tauDRO'}
    k = f{1};
    if ~isfield(meta, k), continue, end
    if ~isfield(eo, k)
        eo.(k) = meta.(k);
    elseif ~isequal(eo.(k), meta.(k))
        error('repair_gto_dro_cells:engineOpts', 'engineOpts.%s = %s but the sheet has %s', ...
              k, mat2str_any(eo.(k)), mat2str_any(meta.(k)));
    end
end
for f = {'depParams', 'arrParams'}
    k = f{1};
    if ~isfield(meta, k), continue, end
    if isfield(eo, k)
        for g = fieldnames(eo.(k)).'
            if isfield(meta.(k), g{1}) && ~isequal(eo.(k).(g{1}), meta.(k).(g{1}))
                error('repair_gto_dro_cells:engineOpts', 'engineOpts.%s.%s = %s but the sheet has %s', ...
                      k, g{1}, mat2str_any(eo.(k).(g{1})), mat2str_any(meta.(k).(g{1})));
            end
        end
    end
    eo.(k) = meta.(k);                         % the sheet's own (possibly completed) parameters
end
end

function s = mat2str_any(v)
% MAT2STR_ANY  Printable form of a numeric / char option.  INPUTS: v.  OUTPUTS: s [char].
if ischar(v), s = v; elseif isnumeric(v) || islogical(v), s = mat2str(v); else, s = class(v); end
end

function Q = load_sheet_or_throw(sheetMat, bfile, snapMat, lg)
% LOAD_SHEET_OR_THROW  Load the sheet; on failure keep the backup (the
% recovery record) and throw repair_gto_dro_cells:sheetUnreadable, naming
% the whole-sheet snapshot snapMat as the recovery source.
% INPUTS: sheetMat; bfile ('' at call start); snapMat; lg.  OUTPUTS: Q.
try
    Q = load(sheetMat);
    assert(isfield(Q, 'OK') && isfield(Q, 'Z8'), 'sheet lacks OK/Z8');
catch ME
    lg('sheet %s unreadable (%s); cell backup %s kept; whole-sheet recovery source %s', ...
       sheetMat, ME.message, bfile, snapMat);
    error('repair_gto_dro_cells:sheetUnreadable', ...
          'sheet %s unreadable (%s); cell backup "%s" kept; recover from the snapshot %s', ...
          sheetMat, ME.message, bfile, snapMat);
end
end

function save_atomic(path, S)
% SAVE_ATOMIC  Save the fields of struct S to path via temp file + movefile.
% INPUTS: path [char]; S struct.  OUTPUTS: none.
tmp = sprintf('%s_tmp%06d.mat', path(1:end-4), randi(999999));
save(tmp, '-struct', 'S');
movefile(tmp, path, 'f');
end

function cp = conj_verdict(z8, rv0, rvf, ndp, meta, lg, iD, iA)
% CONJ_VERDICT  The stage 1 S4 junction test on one entry: NaN when the
% re-solve does not reproduce the entry.
% INPUTS: z8 [8x1]; rv0, rvf [6x1]; ndp; meta; lg; iD; iA.
% OUTPUTS: cp - 1 / 0 / NaN.
cp = NaN;
rv0 = rv0(1:6);  rv0 = rv0(:).';  rvf = rvf(1:6);  rvf = rvf(:).';
try
    seed = seed_from_z8(z8, rv0, 24, ndp.Tnd, ndp.cnd, meta.muStar);
    [z, info] = ms_tfmin(rv0, rvf, seed, ndp.Tnd, ndp.cnd, meta.muStar, ...
                         struct('conjTest', true, 'wallSec', 120));
    if info.converged && max(abs(z - z8)) < 1e-6
        cp = double(info.conj.pass);
    else
        lg('  (%2d,%2d) conjugate re-solve did not reproduce the entry (|dz| %.1e) -- verdict NaN', ...
           iD, iA, max(abs(z - z8)));
    end
catch ME
    lg('  (%2d,%2d) conjugate re-solve error: %s -- verdict NaN', iD, iA, ME.message);
end
end

function bk = backup_cell(Q, iD, iA)
% BACKUP_CELL  Copy one cell's slices.  INPUTS: Q; iD; iA.  OUTPUTS: bk struct.
bk = struct('TF', squeeze(Q.TF(iD, iA, :)), 'FLYKM', squeeze(Q.FLYKM(iD, iA, :)), ...
            'ACCDZ', squeeze(Q.ACCDZ(iD, iA, :)), 'RES', squeeze(Q.RES(iD, iA, :)), ...
            'WALL', squeeze(Q.WALL(iD, iA, :)), 'OK', squeeze(Q.OK(iD, iA, :)), ...
            'Z8', squeeze(Q.Z8(:, iD, iA, :)), 'ATT', Q.ATT(iD, iA));
end

function Q = restore_cell(Q, bk, iD, iA)
% RESTORE_CELL  Put a backed-up cell back.  INPUTS: Q; bk; iD; iA.  OUTPUTS: Q.
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

function v = fdef(s, f, v0)
% FDEF  s.(f) if present else v0.  INPUTS: s; f; v0.  OUTPUTS: v.
if isfield(s, f), v = s.(f); else, v = v0; end
end

function logmsg(f, s)
% LOGMSG  Append-mode file log, stdout when f is empty.  INPUTS: f; s.
if isempty(f), fprintf('%s\n', s);
else, fid = fopen(f, 'a'); fprintf(fid, '%s\n', s); fclose(fid); end
end
