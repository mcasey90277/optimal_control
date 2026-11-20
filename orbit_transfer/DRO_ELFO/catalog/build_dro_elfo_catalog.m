function [chainBlockers, report] = build_dro_elfo_catalog(runIn)
%% Purpose:
%
%   THE CHAIN for the DRO -> ELFO costate catalog AFTER the fleet
%   (run_dro_elfo_torus / run_dro_elfo_batched.sh write the sheets; this
%   function turns finished sheets into a packaged, verified catalog). A
%   copy of GTO_DRO/catalog/build_gto_dro_catalog.m with the parameter block
%   changed; kept from it: one switch per stage, every stage fenced, the
%   blocker list and the waiver rules, the fleet-lock refusal, the
%   staging-copy packaging and the preflight. There is NO deliverable stage
%   (spec section B): the end product is the verified catalog
%   DRO_ELFO/results/costate_catalog_dro_elfo.mat and a verdict file.
%
%     1. census           sheets found (by their exact plan names), cells
%                         solved per rung, retired cells, sheet_key checked
%                         against the plan row
%     2. clearance        body_clearance_sweep per sheet: LOCATED minimum
%                         altitude above the Moon and the Earth, demotes
%                         below the floors
%     3. repairUnsolved   repair_dro_elfo_cells 'unsolved' (phase
%                         continuation)
%     4. repairBasin      repair_dro_elfo_cells 'basin' (conjugate verdict on
%                         every top rung, reverse warm start on refuted ones;
%                         astra C5 acceptance)
%     2b. clearance re-sweep on every sheet a repair changed (keyed to each
%                         entry's z8: only re-solved entries fly)
%     5. package          build_costate_catalog_family on a CLEAN STAGING
%                         COPY of the official sheets, then the additive
%                         sheet_key per packaged sheet (target, depSize,
%                         raanDeg: spec section E), validated by
%                         catalog_schema
%     6. conj             conj_catalog_pass census -> writeback
%     7. gates            gates_catalog_pass census -> writeback
%     8. audit            3 entries replayed THROUGH dro_elfo_pick: the
%                         picker must return the sheet's own entry, which is
%                         then flown (flown miss < meta.gateKm), handed to
%                         pumpkyn tfMin (|dz| < 1e-6) and measured (located
%                         Moon and Earth clearance at the floors)
%     9. pictures         coverage per sheet and rung; t_f and m_f ranges
%    10. preflight        what is ON DISK, whatever ran in this call
%
%   SHEETS BY NAME. The sheet folder also holds <tag>_direct.mat,
%   <tag>_progress.txt, <tag>_clearance.mat, <tag>_repairs.mat, repair
%   backups, experiment files such as <tag>_conj.mat, and possibly a tulip
%   control sheet; a glob such as dro_elfo_*.mat matches several of them.
%   The official sheets are therefore listed by their EXACT names from the
%   front door's plan (run_dro_elfo_torus planOnly with runIn.planOpts),
%   ELFO targets only; the tulip control never enters.
%
%   RESUMABLE. Every stage resumes from its own record: the clearance sweep
%   and both catalog sweeps from their sidecars, the repairs from their
%   _repairs.mat rows (a cell is not re-tried in the same mode). The package
%   stage builds the new catalog to a temporary file and keeps the EXISTING
%   catalog (with its verdicts and sweep sidecars) when the content is
%   unchanged (same sheet keys, has_solution, entry_index and z8); only a
%   changed package replaces it, deleting the stale sweep sidecars and
%   keeping any old .bak_conj / .bak_gates under a dated name. So a killed
%   run is resumed by running the same call again.
%
%   WAIVERS (deliberate, off by default). A blocker is dropped from the
%   VERDICT test only, never from the returned list, and the number waived
%   is logged:
%     allowConjBlockers         'conj: ...' (refuted entries, shipped with
%                               their conj_pass = 0)
%     allowClearanceDemotions   'clearance-demoted: ...' (entries demoted to
%                               OK = false by the clearance sweep)
%   Nothing else is waivable: not 'clearance-incomplete:', 'clearance
%   threw', 'conj incomplete', 'stale:', 'unverified:', 'gates:', nor any
%   other preflight line. One exception the other way: a thrown pictures
%   stage ('pictures threw: ...') is reported in the returned list but does
%   NOT block the verdict (the figures are a view of the sheets, not a
%   check on them).
%
%   PREFLIGHT (on disk; each failure a blocker): 'catalog-missing:';
%   'stale:' the packaged catalog differs by CONTENT from the current sheets
%   (has_solution and the z8 of the OK entries, twin found by sheet_key) or
%   a sheet's clearance sidecar does not cover the current z8 of every OK
%   entry; 'unverified:' no conj / gate record, an accepted entry not run
%   (conj_pass = -1) or ungated, an OK entry below a floor, the audit
%   missing, failed or not produced from the current sheets (fingerprint:
%   sheet, cell, rung and z8 of each replayed entry); and the re-derived
%   'gates:', 'census:', 'clearance-demoted:' (from the sidecar record) and
%   'conj:' (refuted entries in the package), so switching a stage off is
%   not an implicit waiver.
%
%   FLEET LOCK. A call with a sheet-writing stage on (clearance,
%   repairUnsolved, repairBasin, package) refuses with
%   build_dro_elfo_catalog:fleetRunning when the batched driver's lock
%   <catDir>/.dro_elfo_driver.lock exists, else takes it (mkdir, released
%   by onCleanup) and passes ignoreLock to the repair tool. runIn.ignoreLock
%   skips the check and the take; a foreign lock is never removed.
%
%% Inputs:
%
%  runIn                    struct (optional)       Unknown fields refused
%                                                   (:runIn)
%   .stages                 struct | cellstr        Stage switches: a struct
%                                                   of logicals overriding
%                                                   the defaults (all true),
%                                                   or a cellstr naming the
%                                                   ONLY stages to run.
%                                                   Names: census clearance
%                                                   repairUnsolved
%                                                   repairBasin package conj
%                                                   gates audit pictures
%                                                   preflight
%   .catDir                 char                    Sheet folder
%                                                   [DRO_ELFO/results/
%                                                   catalog]
%   .outDir                 char                    Catalog, audit, verdict,
%                                                   pictures, logs
%                                                   [DRO_ELFO/results]
%   .planOpts               struct                  Options for the front
%                                                   door's plan (e.g.
%                                                   .sheets) [struct()]
%   .expectSheets           double                  [12]
%   .floorMoonKm            double                  [500]
%   .floorEarthKm           double                  [200]
%   .maxRepair              double                  Cells per sheet per
%                                                   repair mode [inf]
%   .allowConjBlockers      logical                 [false]
%   .allowClearanceDemotions logical                [false]
%   .ignoreLock             logical                 [false]
%
%% Outputs:
%
%  chainBlockers            cell                    Blocker strings (empty =
%                                                   clean); printed and
%                                                   logged
%
%  report                   struct                  .verified (no blocker
%                                                   left after the waivers)
%                                                   .verdictBlockers .nWaived
%                                                   .catalog .verdictFile
%                                                   .logFile .sheetMats
%
%% Revision History:
%  M. Casey                                                   (c) 11/19/2026
%  M. Casey  fix round 1: helper sizes, pictures-verdict note (c) 11/19/2026
%  Copyright Coorbital Inc.
%% ------------------------ Begin Code Sequence ---------------------------

if nargin < 1 || isempty(runIn), runIn = struct(); end
known = {'stages', 'catDir', 'outDir', 'planOpts', 'expectSheets', 'floorMoonKm', 'floorEarthKm', ...
         'maxRepair', 'allowConjBlockers', 'allowClearanceDemotions', 'ignoreLock'};
bad = setdiff(fieldnames(runIn), known);
if ~isempty(bad)
    error('build_dro_elfo_catalog:runIn', 'unknown field(s): %s', strjoin(bad, ', '));
end
here = fileparts(mfilename('fullpath'));
addpath(fileparts(here));
setup_dro_elfo_paths();
addpath(here);

%% ========================================================================
%  0. PARAMETERS AND STAGE SWITCHES
%% ========================================================================
names = {'census', 'clearance', 'repairUnsolved', 'repairBasin', 'package', 'conj', 'gates', ...
         'audit', 'pictures', 'preflight'};
run_ = cell2struct(num2cell(true(numel(names), 1)), names, 1);
if isfield(runIn, 'stages')
    st = runIn.stages;
    if iscell(st) || isstring(st)
        st = cellstr(st);
        badS = setdiff(st, names);
        if ~isempty(badS), error('build_dro_elfo_catalog:stages', 'unknown stage(s): %s', strjoin(badS, ', ')); end
        for f = names, run_.(f{1}) = ismember(f{1}, st); end
    else
        badS = setdiff(fieldnames(st), names);
        if ~isempty(badS), error('build_dro_elfo_catalog:stages', 'unknown stage(s): %s', strjoin(badS, ', ')); end
        for f = fieldnames(st).', run_.(f{1}) = logical(st.(f{1})); end
    end
end
               catDir = fdef(runIn, 'catDir', fullfile(fileparts(here), 'results', 'catalog'));
               outDir = fdef(runIn, 'outDir', fullfile(fileparts(here), 'results'));
             planOpts = fdef(runIn, 'planOpts', struct());
         expectSheets = fdef(runIn, 'expectSheets', 12);
          floorMoonKm = fdef(runIn, 'floorMoonKm', 500);
         floorEarthKm = fdef(runIn, 'floorEarthKm', 200);
            maxRepair = fdef(runIn, 'maxRepair', inf);
    allowConjBlockers = fdef(runIn, 'allowConjBlockers', false);
allowClearanceDemotions = fdef(runIn, 'allowClearanceDemotions', false);
           ignoreLock = fdef(runIn, 'ignoreLock', false);
if ~isfolder(outDir), mkdir(outDir); end
logFile = fullfile(outDir, ['chain_' char(datetime('now', 'Format', 'yyyyMMdd_HHmmss')) '.log']);
lg = @(varargin) logmsg(logFile, sprintf(varargin{:}));
files = struct('catalog', fullfile(outDir, 'costate_catalog_dro_elfo.mat'), ...
               'audit', fullfile(outDir, 'audit_dro_elfo.mat'), ...
               'verdict', fullfile(outDir, 'chain_verdict_dro_elfo.mat'), ...
               'pictures', fullfile(outDir, 'pictures'), ...
               'stage', fullfile(outDir, 'package_stage'));
floors = struct('moon', floorMoonKm, 'earth', floorEarthKm);
blockers = {};

% the sheets, by their exact plan names (never a glob), ELFO targets only
po = planOpts;  po.planOnly = true;
planAll = run_dro_elfo_torus([], po);
plan = planAll(strcmp({planAll.target}, 'elfo'));
if numel(plan) < numel(planAll)
    lg('0. plan: %d non-ELFO sheet(s) in the plan ignored (the chain is ELFO only)', numel(planAll) - numel(plan));
end
sheetMats = {};  planRow = [];
for kp = 1:numel(plan)
    p = fullfile(catDir, [plan(kp).tag '.mat']);
    if isfile(p), sheetMats{end+1} = p;  planRow(end+1) = kp; end
end
nS = numel(sheetMats);
lg('0. %s: %d of %d planned ELFO sheets on disk in %s', char(datetime('now', 'Format', 'yyyy-MM-dd HH:mm:ss')), ...
   nS, numel(plan), catDir);

% the fleet lock, for any call that writes sheets
lockPath = fullfile(catDir, '.dro_elfo_driver.lock');
if ~ignoreLock && (run_.clearance || run_.repairUnsolved || run_.repairBasin || run_.package)
    if isfolder(lockPath) || isfile(lockPath)
        error('build_dro_elfo_catalog:fleetRunning', ...
              'fleet driver lock %s exists; a fleet or a repair is running (runIn.ignoreLock overrides)', lockPath);
    end
    if ~isfolder(catDir), mkdir(catDir); end
    [lockMade, ~, lockId] = mkdir(lockPath);
    if ~lockMade || strcmp(lockId, 'MATLAB:MKDIR:DirectoryExists')
        error('build_dro_elfo_catalog:fleetRunning', 'could not take the fleet driver lock %s', lockPath);
    end
    lockRelease = onCleanup(@() rmdir(lockPath));
end

%% ========================================================================
%  1. CENSUS
%% ========================================================================
if run_.census
    try
        lg('1. census: %d of %d sheets on disk', nS, expectSheets);
        if nS < expectSheets, blockers{end+1} = sprintf('census: %d of %d sheets', nS, expectSheets); end
        for ks = 1:nS
            vars = who('-file', sheetMats{ks});
            Q = load(sheetMats{ks}, 'OK', 'ATT', 'rungs');
            [nD, nA, ~] = size(Q.OK);
            solved = squeeze(sum(sum(Q.OK, 1), 2)).';
            retired = nnz(~any(Q.OK, 3) & Q.ATT >= 2);
            [~, nm] = fileparts(sheetMats{ks});
            lg('   %-24s cells %3d  per rung [%s]  retired %2d', nm, nD*nA, num2str(solved), retired);
            keyOK = ismember('sheet_key', vars);
            if keyOK
                K = load(sheetMats{ks}, 'sheet_key');
                keyOK = isequaln(K.sheet_key, plan(planRow(ks)).sheet_key);
            end
            if ~keyOK
                blockers{end+1} = sprintf('census: %s carries no sheet_key or one that is not its plan row''s', nm);
            end
        end
    catch ME
        lg('1. census THREW: %s', ME.message);  blockers{end+1} = ['census threw: ' ME.message];
    end
end

%% ========================================================================
%  2. CLEARANCE SWEEP (sheets; demotes)
%% ========================================================================
clrBlk = repmat({''}, 1, nS);
if run_.clearance
    for ks = 1:nS
        clrBlk{ks} = sweep_sheet(sheetMats{ks}, floors, logFile, '2. clearance');
    end
end

%% ========================================================================
%  3./4. REPAIRS (sheets; sidecars)
%% ========================================================================
changed = false(1, nS);
for mode = {'unsolved', 'basin'}
    sw = pick(strcmp(mode{1}, 'unsolved'), run_.repairUnsolved, run_.repairBasin);
    if ~sw, continue, end
    for ks = 1:nS
        try
            [~, nm] = fileparts(sheetMats{ks});
            R = repair_dro_elfo_cells(sheetMats{ks}, mode{1}, ...
                    struct('engineOpts', plan(planRow(ks)).optsA, 'maxCells', maxRepair, ...
                           'logFile', logFile, 'ignoreLock', true));
            lg(['3/4. repair %-8s %-24s candidates %d attempted %d repaired %d restored %d ', ...
                'failed %d refuted %d recovered %d kept-previous-rungs %d'], mode{1}, nm, R.nCandidates, ...
               R.nAttempted, R.nRepaired, R.nRestored, R.nFailed, R.nConjRefuted, R.nRecovered, R.nKeptOld);
            if R.nRepaired > 0 || R.nRecovered > 0, changed(ks) = true; end
            if R.nFailed > 0, blockers{end+1} = sprintf('repair %s: %s %d failed', mode{1}, nm, R.nFailed); end
        catch ME
            lg('3/4. repair %s THREW on %s: %s', mode{1}, sheetMats{ks}, ME.message);
            blockers{end+1} = ['repair threw: ' ME.message];
        end
    end
end

%% ========================================================================
%  2b. CLEARANCE RE-SWEEP on every sheet a repair changed
%% ========================================================================
if run_.clearance
    for ks = find(changed)
        clrBlk{ks} = sweep_sheet(sheetMats{ks}, floors, logFile, '2b. clearance re-sweep');
    end
    blockers = [blockers, clrBlk(~cellfun(@isempty, clrBlk))];
end

%% ========================================================================
%  5. PACKAGE (clean staging copy of the official sheets; sheet_key per
%     packaged sheet; kept when the content is unchanged)
%% ========================================================================
if run_.package
    try
        spec = struct( ...
          'glob', 'dro_elfo_tau*_raan*.mat', ...
          'name', 'costate_catalog_dro_elfo', ...
          'description', ['Minimum-time low-thrust DRO -> ELFO costate CATALOG: one (DRO size x ', ...
            'ELFO plane orientation) phasing sheet per pair, thrust rungs walked down from 15 N. ', ...
            'COMPACT format: canonical ND quantities only; derived values via .derive. Each entry ', ...
            'z8 = [lambda_r(3); lambda_v(3); lambda_m; tf] is a converged PMP solution accepted by ', ...
            'pumpkyn.cr3bp.tfMin, solved with a 500 km lunar floor and swept for the LOCATED ', ...
            'minimum altitude above the Moon (500 km) and the Earth (200 km). SHEET IDENTITY: ', ...
            'sheets(k).sheet_key = (target, depSize = DRO period ND, raanDeg = ELFO plane ', ...
            'orientation, deg); tau_dep is the DRO period, tau_arr the ELFO Kepler period, the ', ...
            'orientation is arr_params.raanDeg; Np, pm are NaN and period_tulip_nd is the arrival ', ...
            'period, as the packager fills them for a non-tulip arrival. Pick with dro_elfo_pick. ', ...
            'conj_pass: 1 pass, 0 refuted (an extremal, not a minimum), -1 not run.'], ...
          'provenance', ['DRO -> ELFO fleet (run_dro_elfo_torus, thrust_ladder_library) and chain ', ...
            '(build_dro_elfo_catalog): located clearance sweep, phase-continuation and ', ...
            'reverse-warm-start repairs (astra C5 acceptance) recorded per sheet in ', ...
            '_clearance.mat / _repairs.mat sidecars. Casey/Koblick 2026-11.'], ...
          'depReconstruction', ['dro: [tau,rv] = get_family_orbit(''dro'', dep_params); elfo: ', ...
            '[tau,rv] = get_family_orbit(''elfo'', arr_params) (oe, raanDeg); state at phase f: ', ...
            'phase_state(tau, rv)(f)']);
        if isfolder(files.stage), rmdir(files.stage, 's'); end
        mkdir(files.stage);
        for ks = 1:nS, copyfile(sheetMats{ks}, files.stage); end
        newMat = [files.catalog(1:end-4) '_new.mat'];
        if isfile(newMat), delete(newMat); end
        catNew = build_costate_catalog_family(files.stage, newMat, spec);
        catNew = attach_sheet_keys(catNew, sheetMats);
        problems = catalog_schema('validate', catNew);
        assert(isempty(problems), 'build_dro_elfo_catalog:schema', 'the keyed catalog fails validation: %s', ...
               strjoin(problems, '; '));
        unchanged = false;
        if isfile(files.catalog)
            L = load(files.catalog);  fn = fieldnames(L);
            unchanged = same_package(catNew, L.(fn{1}));
        end
        if unchanged
            delete(newMat);
            lg('5. package: content unchanged (%d sheets, %d entries) -- the existing catalog and its verdicts kept', ...
               numel(catNew.sheets), catNew.n_entries);
        else
            for sc = {'_conjprog.mat', '_gatesprog.mat'}
                pp = [files.catalog(1:end-4) sc{1}];
                if isfile(pp), delete(pp); end
            end
            for bk = {'.bak_conj', '.bak_gates'}
                pp = [files.catalog bk{1}];
                if isfile(pp)
                    dated = [pp '_' char(datetime('now', 'Format', 'yyyyMMdd_HHmmss'))];
                    while isfile(dated), dated = [dated '_b']; end
                    movefile(pp, dated);
                    lg('5. package: previous %s kept as %s', bk{1}, dated);
                end
            end
            Sout = struct('costate_catalog_dro_elfo', catNew);
            save(newMat, '-struct', 'Sout');
            movefile(newMat, files.catalog, 'f');
            lg('5. package: %d sheets, %d entries -> %s', numel(catNew.sheets), catNew.n_entries, files.catalog);
        end
        rmdir(files.stage, 's');         % kept only on failure, for inspection
    catch ME
        lg('5. package THREW: %s', ME.message);  blockers{end+1} = ['package threw: ' ME.message];
    end
end

%% ========================================================================
%  6./7. CONJUGATE AND HYPOTHESIS SWEEPS (package; sidecars; writeback)
%% ========================================================================
if run_.conj && isfile(files.catalog)
    try
        Sc = conj_catalog_pass(files.catalog, struct('logFile', logFile));
        if Sc.done
            Sc = conj_catalog_pass(files.catalog, struct('logFile', logFile, 'writeback', true));
            lg('6. conj: %d pass, %d FAIL, %d not run', Sc.nPass, Sc.nFail, Sc.nNotrun);
            if Sc.nFail > 0
                blockers{end+1} = sprintf(['conj: %d refuted entries (consult conj_pass; ', ...
                    'repairBasin re-solves refuted TOP rungs)'], Sc.nFail);
            end
        else
            blockers{end+1} = sprintf('conj incomplete: %d to do', Sc.nTodo);
        end
    catch ME
        lg('6. conj THREW: %s', ME.message);  blockers{end+1} = ['conj threw: ' ME.message];
    end
end
if run_.gates && isfile(files.catalog)
    try
        Sg = gates_catalog_pass(files.catalog, struct('logFile', logFile));
        if Sg.done
            Sg = gates_catalog_pass(files.catalog, struct('logFile', logFile, 'writeback', true));
            lg('7. gates: %d done, H2 fail %d, H3 fail %d, abnormal %d', ...
               Sg.nDone, Sg.nH2fail, Sg.nH3fail, Sg.nAbnormal);
            if Sg.nH2fail + Sg.nH3fail + Sg.nAbnormal > 0
                blockers{end+1} = sprintf('gates: %d entries fail a hypothesis gate', ...
                                          Sg.nH2fail + Sg.nH3fail + Sg.nAbnormal);
            end
        else
            blockers{end+1} = sprintf('gates incomplete: %d to do', Sg.nTodo);
        end
    catch ME
        lg('7. gates THREW: %s', ME.message);  blockers{end+1} = ['gates threw: ' ME.message];
    end
end

%% ========================================================================
%  8. AUDIT (3 entries replayed through dro_elfo_pick)
%% ========================================================================
if run_.audit
    try
        results = audit_entries(files.catalog, sheetMats, floors, lg);
        save(files.audit, 'results');
        lg('8. audit: picked %d/%d, flown %d/%d, accepted %d/%d, clearance %d/%d', ...
           nnz([results.passPick]), numel(results), nnz([results.passFly]), numel(results), ...
           nnz([results.passAccept]), numel(results), nnz([results.passClear]), numel(results));
        if numel(results) < 3
            blockers{end+1} = sprintf('audit: only %d entr%s to replay (3 wanted)', numel(results), ...
                                      pick(isscalar(results), 'y', 'ies'));
        end
        if ~all([results.passPick] & [results.passFly] & [results.passAccept] & [results.passClear])
            blockers{end+1} = 'audit: a replayed entry failed a gate';
        end
    catch ME
        lg('8. audit THREW: %s', ME.message);  blockers{end+1} = ['audit threw: ' ME.message];
    end
end

%% ========================================================================
%  9. PICTURES (coverage per sheet and rung; t_f and m_f ranges)
%% ========================================================================
if run_.pictures
    visBefore = get(0, 'DefaultFigureVisible');
    figsBefore = findall(0, 'Type', 'figure');
    try
        if ~isfolder(files.pictures), mkdir(files.pictures); end
        set(0, 'DefaultFigureVisible', 'off');
        catalog_pictures(sheetMats, files.pictures);
        lg('9. pictures -> %s', files.pictures);
    catch ME
        lg('9. pictures THREW: %s', ME.message);  blockers{end+1} = ['pictures threw: ' ME.message];
    end
    set(0, 'DefaultFigureVisible', visBefore);
    figsNow = findall(0, 'Type', 'figure');
    close(figsNow(~ismember(figsNow, figsBefore)));
end

%% ========================================================================
%  10. PREFLIGHT (what is on disk; re-derives what switched-off stages
%      would report)
%% ========================================================================
if run_.preflight
    try
        [pre, nRefuted] = preflight(files, sheetMats, floors, expectSheets);
    catch ME
        pre = {['unverified: preflight threw: ' ME.message]};  nRefuted = 0;
    end
    keep = true(size(pre));
    for kb = 1:numel(pre)
        pfx = pre{kb}(1:find(pre{kb} == ':', 1));
        if any(strcmp(pfx, {'gates:', 'census:', 'clearance-demoted:'})) && any(startsWith(blockers, pfx))
            keep(kb) = false;
        end
    end
    pre = pre(keep);
    for kb = 1:numel(pre), lg('10. preflight: %s', pre{kb}); end
    blockers = [blockers, pre];
    if nRefuted > 0 && ~any(startsWith(blockers, 'conj:'))
        blockers{end+1} = sprintf('conj: %d refuted entries on disk (conj_pass = 0 in the package)', nRefuted);
    end
end

%% ========================================================================
%  VERDICT (waivers apply to the verdict only)
%% ========================================================================
verdictBlockers = blockers(~startsWith(blockers, 'pictures threw'));
nWaived = 0;
for wv = {{allowConjBlockers, 'conj:', 'allowConjBlockers'}, ...
          {allowClearanceDemotions, 'clearance-demoted:', 'allowClearanceDemotions'}}
    if wv{1}{1}
        isW = startsWith(verdictBlockers, wv{1}{2});
        if any(isW), lg('%s: %d blocker(s) waived for the verdict', wv{1}{3}, nnz(isW)); end
        nWaived = nWaived + nnz(isW);
        verdictBlockers = verdictBlockers(~isW);
    end
end
verified = run_.preflight && isempty(verdictBlockers);
if isempty(blockers)
    lg('CHAIN CLEAN: every stage that ran, passed.');
else
    lg('CHAIN COMPLETE WITH %d BLOCKER(S):', numel(blockers));
    for kb = 1:numel(blockers), lg('  - %s', blockers{kb}); end
end
if verified
    lg('CATALOG VERIFIED (%d blocker(s) waived): %s', nWaived, files.catalog);
else
    lg('CATALOG NOT VERIFIED: %d blocker(s) stand after the waivers%s', numel(verdictBlockers), ...
       pick(run_.preflight, '', ' (preflight off)'));
end
verdict = struct('date', char(datetime('now', 'Format', 'yyyy-MM-dd HH:mm:ss')), 'verified', verified, ...
                 'blockers', {blockers}, 'verdictBlockers', {verdictBlockers}, 'nWaived', nWaived, ...
                 'sheetMats', {sheetMats}, 'logFile', logFile);
save(files.verdict, '-struct', 'verdict');
chainBlockers = blockers;
report = struct('verified', verified, 'verdictBlockers', {verdictBlockers}, 'nWaived', nWaived, ...
                'catalog', files.catalog, 'verdictFile', files.verdict, 'logFile', logFile, ...
                'sheetMats', {sheetMats});
end

% ------------------------------------------------------------------------
function cat_ = attach_sheet_keys(cat_, sheetMats)
% ATTACH_SHEET_KEYS  The additive per-sheet sheet_key (spec section E):
% each packaged sheet is matched to its source sheet file by (DRO period,
% ELFO raanDeg) and takes that file's sheet_key, which must agree with the
% packaged recipes (target = arr_family, depSize = dep_params.tau, raanDeg =
% arr_params.raanDeg). Exactly one source per packaged sheet.
% INPUTS: cat_ [1 x 1] catalog struct; sheetMats {1 x nS} char paths.
% OUTPUTS: cat_ [1 x 1] catalog struct with sheets(k).sheet_key [1 x 1].
keys = cell(1, numel(sheetMats));
for ks = 1:numel(sheetMats)
    K = load(sheetMats{ks}, 'sheet_key');
    assert(isfield(K, 'sheet_key'), 'build_dro_elfo_catalog:sheetKey', '%s carries no sheet_key', sheetMats{ks});
    keys{ks} = K.sheet_key;
end
for kp = 1:numel(cat_.sheets)
    sh = cat_.sheets(kp);
    hit = false(1, numel(keys));
    for ks = 1:numel(keys)
        k = keys{ks};
        hit(ks) = strcmp(k.target, sh.arr_family) && abs(k.depSize - sh.dep_params.tau) < 1e-12 && ...
                  isfield(sh.arr_params, 'raanDeg') && ...
                  abs(mod(k.raanDeg - sh.arr_params.raanDeg + 180, 360) - 180) < 1e-9;
    end
    assert(nnz(hit) == 1, 'build_dro_elfo_catalog:sheetKey', ...
           'packaged sheet %d (DRO %g, raanDeg %g) matches %d source sheet keys, not 1', ...
           kp, sh.dep_params.tau, sh.arr_params.raanDeg, nnz(hit));
    cat_.sheets(kp).sheet_key = keys{hit};
end
end

function same = same_package(a, b)
% SAME_PACKAGE  Two catalogs hold the same entries: same sheet count, and
% per sheet (in order) the same sheet_key, has_solution, entry_index, z8
% and rungs.  INPUTS: a, b [1 x 1] catalog structs.  OUTPUTS: same [1 x 1]
% logical.
same = numel(a.sheets) == numel(b.sheets) && isfield(b.sheets, 'sheet_key') && ...
       isequal(a.rungs_N, b.rungs_N);
for ks = 1:numel(a.sheets) * same
    sa = a.sheets(ks);  sb = b.sheets(ks);
    same = same && isequaln(sa.sheet_key, sb.sheet_key) && isequal(sa.has_solution, sb.has_solution) ...
                && isequal(sa.entry_index, sb.entry_index) && isequaln(sa.z8, sb.z8);
end
end

function [blk, nRefuted] = preflight(files, sheetMats, floors, expectSheets)
% PREFLIGHT  What is ON DISK, independent of the stages that ran this call
% (header PREFLIGHT).  INPUTS: files [1 x 1] struct (.catalog .audit char);
% sheetMats {1 x nS} char; floors [1 x 1] struct (.moon .earth [1 x 1] km);
% expectSheets [1 x 1].
% OUTPUTS: blk {1 x b} char blocker strings; nRefuted [1 x 1] (accepted
% entries with conj_pass = 0 in the package).
blk = {};  nRefuted = 0;
nS = numel(sheetMats);
nDemoted = 0;  nStaleClr = 0;  nBelowOK = 0;  staleNames = {};
for ks = 1:nS
    Q = load(sheetMats{ks}, 'OK', 'Z8');
    [~, nm] = fileparts(sheetMats{ks});
    C = clearance_record(sheetMats{ks}, Q, floors);
    nDemoted = nDemoted + C.nDemoted;
    if nnz(Q.OK) > 0 && ~C.covers
        nStaleClr = nStaleClr + 1;  staleNames{end+1} = nm;
    end
    nBelowOK = nBelowOK + C.nBelowOK;
end
if nStaleClr > 0
    blk{end+1} = sprintf('stale: the clearance sidecar of %d sheet(s) does not cover the current z8 of every OK entry (%s)', ...
                         nStaleClr, strjoin(staleNames, ', '));
end
if nBelowOK > 0
    blk{end+1} = sprintf('unverified: %d OK entries are below a clearance floor (re-run the clearance stage)', nBelowOK);
end
if nDemoted > 0
    blk{end+1} = sprintf('clearance-demoted: %d entries demoted by the clearance sweep (OK = false, absent from the catalog)', nDemoted);
end
if ~isfile(files.catalog)
    blk{end+1} = sprintf('catalog-missing: %s does not exist', files.catalog);
    if nS < expectSheets
        blk{end+1} = sprintf('census: %d official sheets on disk, %d expected', nS, expectSheets);
    end
    return
end
L = load(files.catalog);  fn = fieldnames(L);  cat_ = L.(fn{1});
nonEmpty = 0;  nStaleSheet = 0;  nMissingRec = 0;  nNoConj = 0;  nNoGate = 0;
nH2 = 0;  nH3 = 0;  nAbn = 0;
lamVTol = 1e-6;
if isfield(cat_, 'hyp_gates') && isfield(cat_.hyp_gates, 'lamVTol'), lamVTol = cat_.hyp_gates.lamVTol; end
hasKeys = isfield(cat_.sheets, 'sheet_key');
for ks = 1:nS
    Q = load(sheetMats{ks}, 'OK', 'Z8', 'sheet_key');
    if nnz(Q.OK) == 0, continue, end
    nonEmpty = nonEmpty + 1;
    kp = [];
    if hasKeys && isfield(Q, 'sheet_key')
        for kk = 1:numel(cat_.sheets)
            if isequaln(cat_.sheets(kk).sheet_key, Q.sheet_key), kp = kk; break, end
        end
    end
    if isempty(kp), nStaleSheet = nStaleSheet + 1; continue, end
    sh = cat_.sheets(kp);
    same = isequal(sh.has_solution, Q.OK);
    if same
        Z = reshape(Q.Z8, 8, []);
        same = isequaln(sh.z8(:, sh.entry_index(Q.OK)), Z(:, Q.OK(:)));
    end
    if ~same, nStaleSheet = nStaleSheet + 1; end
    if ~(isfield(sh, 'conj_pass') && isfield(sh, 'gate_min_lamv') && isfield(sh, 'gate_dimS')) ...
            || isempty(sh.conj_pass) || isempty(sh.gate_min_lamv)
        nMissingRec = nMissingRec + 1;  continue
    end
    nNoConj  = nNoConj + nnz(sh.has_solution & sh.conj_pass == -1);
    nNoGate  = nNoGate + nnz(sh.has_solution & isnan(sh.gate_min_lamv));
    nRefuted = nRefuted + nnz(sh.has_solution & sh.conj_pass == 0);
    dn  = sh.has_solution & ~isnan(sh.gate_min_lamv);
    nH2 = nH2 + nnz(dn & sh.gate_min_lamv <= lamVTol);
    nH3 = nH3 + nnz(dn & sh.gate_min_qmt <= 0);
    nAbn = nAbn + nnz(dn & sh.gate_dimS ~= 1);
end
if nH2 + nH3 + nAbn > 0
    blk{end+1} = sprintf('gates: %d entries fail a hypothesis gate (on disk: H2 %d, H3 %d, abnormal %d)', ...
                         nH2 + nH3 + nAbn, nH2, nH3, nAbn);
end
if min(numel(cat_.sheets), nS) < expectSheets
    blk{end+1} = sprintf('census: %d packaged / %d official sheets on disk, %d expected', ...
                         numel(cat_.sheets), nS, expectSheets);
end
if ~hasKeys
    blk{end+1} = 'stale: the packaged catalog carries no sheet_key (re-run with package on)';
end
if numel(cat_.sheets) ~= nonEmpty
    blk{end+1} = sprintf('stale: the package holds %d sheets, the sheets on disk have accepted entries on %d', ...
                         numel(cat_.sheets), nonEmpty);
end
if nStaleSheet > 0
    blk{end+1} = sprintf('stale: %d sheet(s) differ from the package (has_solution or z8); re-run with package on', nStaleSheet);
end
if nMissingRec > 0
    blk{end+1} = sprintf('unverified: %d packaged sheet(s) carry no conj_pass / gate record', nMissingRec);
end
if nNoConj > 0 || nNoGate > 0
    blk{end+1} = sprintf('unverified: %d accepted entries without a conj verdict, %d without hypothesis gates', nNoConj, nNoGate);
end
% the audit: present, passed, and replayed THESE sheets
if ~isfile(files.audit)
    blk{end+1} = sprintf('unverified: audit file %s does not exist', files.audit);
else
    Au = load(files.audit);
    okAudit = isfield(Au, 'results') && ~isempty(Au.results) && isfield(Au.results, 'z8') && ...
              all([Au.results.passPick] & [Au.results.passFly] & [Au.results.passAccept] & [Au.results.passClear]);
    if okAudit
        nms = cell(size(sheetMats));
        for ks = 1:nS, [~, nms{ks}] = fileparts(sheetMats{ks}); end
        for kr = 1:numel(Au.results)
            r = Au.results(kr);
            ks = find(strcmp(nms, r.sheet), 1);
            if isempty(ks), okAudit = false; break, end
            Q = load(sheetMats{ks}, 'OK', 'Z8');
            if ~(Q.OK(r.iD, r.iA, r.kr) && isequal(Q.Z8(:, r.iD, r.iA, r.kr), r.z8(:)))
                okAudit = false;  break
            end
        end
    end
    if ~okAudit
        blk{end+1} = 'unverified: the audit did not pass or was not produced from the current sheets (fingerprint mismatch)';
    end
end
end

function C = clearance_record(sheetMat, Q, floors)
% CLEARANCE_RECORD  What the clearance sidecar says about the CURRENT
% sheet: covers (every OK entry has a record for its current z8), nDemoted
% (records for the current z8 below a floor whose entry is OK = false),
% nBelowOK (OK entries whose current record is below a floor).
% INPUTS: sheetMat [1 x n] char; Q [1 x 1] struct (.OK [nD x nA x nR],
% .Z8 [8 x nD x nA x nR]); floors [1 x 1] struct (.moon .earth km).
% OUTPUTS: C [1 x 1] struct .covers logical, .nDemoted, .nBelowOK [1 x 1].
C = struct('covers', false, 'nDemoted', 0, 'nBelowOK', 0);
sideMat = [sheetMat(1:end-4) '_clearance.mat'];
if ~isfile(sideMat), C.covers = nnz(Q.OK) == 0; return, end
E = load(sideMat);
if ~all(isfield(E, {'MOONKM', 'EARTHKM', 'done', 'Z8SWEPT'})) || ~isequal(size(E.Z8SWEPT), size(Q.Z8))
    return
end
A = reshape(E.Z8SWEPT, 8, []);  B = reshape(Q.Z8, 8, []);
cur = E.done(:).' & all(A == B | (isnan(A) & isnan(B)), 1);        % record for the current z8
cur = reshape(cur, size(Q.OK));
fm = floors.moon;   if isempty(fm), fm = -inf; end
fe = floors.earth;  if isempty(fe), fe = -inf; end
below = cur & ~(E.MOONKM >= fm & E.EARTHKM >= fe);
C.covers   = all(cur(Q.OK));
C.nDemoted = nnz(below & ~Q.OK);
C.nBelowOK = nnz(below & Q.OK);
end

function results = audit_entries(catMat, sheetMats, floors, lg)
% AUDIT_ENTRIES  Up to 3 OK entries (distinct rungs and sheets preferred,
% deterministic (sheet, iD, iA, kr) walk) replayed THROUGH dro_elfo_pick:
% the picker, asked for the entry's sheet key, cell phases and rung, must
% return exactly that entry (passPick); the returned z8 is flown from the
% sheet-local endpoints (passFly: flown miss < meta.gateKm), handed to
% pumpkyn tfMin (passAccept: |dz| < 1e-6) and measured with the shared
% located minimum (passClear: Moon and Earth at the floors).
% INPUTS: catMat [1 x n] char; sheetMats {1 x nS} char; floors [1 x 1]
% struct (.moon .earth km); lg [1 x 1] function handle.
% OUTPUTS: results [n x 1] struct, n <= 3 (.sheet .iD .iA .kr .rungN .z8 .tf_nd
% .flownKm .acceptDz .moonKm .earthKm .passPick .passFly .passAccept
% .passClear).
assert(isfile(catMat), 'build_dro_elfo_catalog:audit', 'no catalog %s to audit', catMat);
L = load(catMat);  fn = fieldnames(L);  cat_ = L.(fn{1});
cand = zeros(0, 5);                           % [ks iD iA kr rungN]
for ks = 1:numel(sheetMats)
    Q = load(sheetMats{ks}, 'OK', 'rungs');
    [nD, nA, nR] = size(Q.OK);
    for iD = 1:nD
        for iA = 1:nA
            for kr = 1:nR
                if Q.OK(iD, iA, kr), cand(end+1, :) = [ks iD iA kr Q.rungs(kr)]; end
            end
        end
    end
end
nPick = min(3, size(cand, 1));
pk = zeros(nPick, 1);
for kk = 1:nPick                              % greedy: new rung + new sheet first
    score = -inf(size(cand, 1), 1);
    for kc = 1:size(cand, 1)
        if any(pk(1:kk-1) == kc), continue, end
        score(kc) = double(~any(cand(pk(1:kk-1), 5) == cand(kc, 5))) + ...
                    double(~any(cand(pk(1:kk-1), 1) == cand(kc, 1)));
    end
    [~, pk(kk)] = max(score);
end
results = repmat(struct('sheet', '', 'iD', 0, 'iA', 0, 'kr', 0, 'rungN', 0, 'z8', zeros(8, 1), ...
    'tf_nd', NaN, 'flownKm', NaN, 'acceptDz', NaN, 'moonKm', NaN, 'earthKm', NaN, ...
    'passPick', false, 'passFly', false, 'passAccept', false, 'passClear', false), nPick, 1);
fm = floors.moon;   if isempty(fm), fm = -inf; end
fe = floors.earth;  if isempty(fe), fe = -inf; end
for kk = 1:nPick
    ks = cand(pk(kk), 1);  iD = cand(pk(kk), 2);  iA = cand(pk(kk), 3);  kr = cand(pk(kk), 4);
    Q = load(sheetMats{ks}, 'Z8', 'sD', 'sA', 'rungs', 'meta', 'sheet_key');
    meta = Q.meta;  [~, nm] = fileparts(sheetMats{ks});
    r = results(kk);
    r.sheet = nm;  r.iD = iD;  r.iA = iA;  r.kr = kr;  r.rungN = Q.rungs(kr);
    r.z8 = Q.Z8(:, iD, iA, kr);
    e = dro_elfo_pick(cat_, Q.sheet_key.depSize, Q.sheet_key.raanDeg, Q.sD(iD), Q.sA(iA), Q.rungs(kr));
    r.passPick = e.iD == iD && e.iA == iA && e.kr == kr && isequal(e.z8, r.z8);
    r.tf_nd = e.tf_nd;
    [tD, rvD] = get_family_orbit(meta.depFamily, meta.depParams);
    [tA, rvA] = get_family_orbit(meta.arrFamily, meta.arrParams);
    depAt = phase_state(tD, rvD);  arrAt = phase_state(tA, rvA);
    rv0 = depAt(e.sD);  rvf = arrAt(e.sA);
    ndp = nd_propulsion(e.thrustN, meta.ispS, meta.m0kg, meta.lStar, meta.tStar);
    phys = struct('Tnd', ndp.Tnd, 'cnd', ndp.cnd, 'muStar', meta.muStar, 'lStar', meta.lStar, ...
                  'tStar', meta.tStar, 'm0kg', meta.m0kg);
    F = fly_transfer(e.z8, rv0(:), rvf(:), phys);
    r.flownKm = F.flyKm;  r.passFly = F.flyKm < meta.gateKm;
    rv0r = rv0(1:6);  rvfr = rvf(1:6);
    zAcc = quiet_call(@() pumpkyn.cr3bp.tfMin(rv0r(:).', rvfr(:).', e.z8, ndp.Tnd, ndp.cnd, meta.muStar));
    r.acceptDz = max(abs(zAcc(:) - e.z8(:)));  r.passAccept = r.acceptDz < 1e-6;
    lo = struct('lStar', meta.lStar, 'phys', phys);
    Cm = located_min_altitude(F.t, F.Y, [1 - meta.muStar 0 0], 1737.4, fm, lo);
    Ce = located_min_altitude(F.t, F.Y, [-meta.muStar 0 0], 6378.137, fe, lo);
    r.moonKm = Cm.hMinKm;  r.earthKm = Ce.hMinKm;  r.passClear = Cm.pass && Ce.pass;
    results(kk) = r;
    lg('   audit %-24s (%d,%d) %5.2f N  pick %d  flown %.3f km  |dz| %.1e  Moon %.0f km  Earth %.0f km', ...
       nm, iD, iA, r.rungN, r.passPick, r.flownKm, r.acceptDz, r.moonKm, r.earthKm);
end
end

function out = quiet_call(f)
% QUIET_CALL  f() with its console output captured and dropped (pumpkyn
% tfMin prints its solver iterations).  INPUTS: f [1 x 1] function handle,
% no arguments.  OUTPUTS: out [any] - f's first output (here z [8 x 1]).
assert(isa(f, 'function_handle'), 'build_dro_elfo_catalog:quietCall', 'f must be a function handle');
out = [];
evalc('out = f();');
end

function blk = sweep_sheet(sheetMat, floors, logFile, label)
% SWEEP_SHEET  The clearance sweep (census + demotion) on one sheet,
% fenced.  INPUTS: sheetMat [1 x n] char; floors [1 x 1] struct (.moon
% .earth km); logFile, label [1 x m] char.  OUTPUTS: blk [1 x k] char
% blocker string, '' when clean.
blk = '';
[~, nm] = fileparts(sheetMat);
try
    S = body_clearance_sweep(sheetMat, struct('floorMoonKm', floors.moon, 'floorEarthKm', floors.earth, ...
                                              'logFile', logFile, 'demote', true));
    if S.done
        logmsg(logFile, sprintf('%s %-24s min Moon %.0f km, min Earth %.0f km, %d flown now, %d below (record), %d not located', ...
               label, nm, S.minMoonKm, S.minEarthKm, S.nFlown, S.nBelow, S.nNonFinite));
        if S.nBelow > 0
            blk = sprintf('clearance-demoted: %s demoted %d entries (Moon %d, Earth %d; OK = false, absent from the catalog)', ...
                          nm, S.nBelow, S.nBelowMoon, S.nBelowEarth);
        end
    else
        blk = sprintf('clearance-incomplete: %s incomplete (%d to do)', nm, S.nTodo);
        logmsg(logFile, sprintf('%s %s', label, blk));
    end
catch ME
    logmsg(logFile, sprintf('%s THREW on %s: %s', label, nm, ME.message));
    blk = ['clearance threw: ' ME.message];
end
end

function catalog_pictures(sheetMats, picDir)
% CATALOG_PICTURES  One coverage figure per sheet (OK per rung on the
% nD x nA torus) and one summary figure (t_f and m_f ranges per rung over
% every sheet).  INPUTS: sheetMats {1 x nS} char; picDir [1 x n] char.
% OUTPUTS: none (PNG files written).
nS = numel(sheetMats);
rungAll = [];  tfAll = [];  mfAll = [];
for ks = 1:nS
    Q = load(sheetMats{ks}, 'OK', 'TF', 'rungs', 'meta', 'sD', 'sA');
    [nD, nA, nR] = size(Q.OK);
    [~, nm] = fileparts(sheetMats{ks});
    fig = figure('Color', 'w', 'Position', [100 100 max(600, 160*nR) 320]);
    tl = tiledlayout(fig, 1, nR, 'TileSpacing', 'compact', 'Padding', 'compact');
    for kr = 1:nR
        ax = nexttile(tl);
        imagesc(ax, Q.sA, Q.sD, double(Q.OK(:, :, kr)));
        set(ax, 'YDir', 'normal');  clim(ax, [0 1]);  colormap(ax, [0.92 0.92 0.92; 0.1 0.5 0.2]);
        title(ax, sprintf('%g N: %d/%d', Q.rungs(kr), nnz(Q.OK(:, :, kr)), nD*nA), 'FontSize', 11);
        xlabel(ax, 'arrival phase');  if kr == 1, ylabel(ax, 'departure phase'); end
    end
    title(tl, strrep(nm, '_', ' '), 'FontSize', 13);
    exportgraphics(fig, fullfile(picDir, [nm '_coverage.png']), 'Resolution', 120);
    close(fig);
    for kr = 1:nR
        m = Q.OK(:, :, kr);
        tf = Q.TF(:, :, kr);  tf = tf(m);
        ndp = nd_propulsion(Q.rungs(kr), Q.meta.ispS, Q.meta.m0kg, Q.meta.lStar, Q.meta.tStar);
        mf = 1 - ndp.Tnd*tf/ndp.cnd;
        rungAll = [rungAll; Q.rungs(kr)*ones(numel(tf), 1)];
        tfAll   = [tfAll; tf(:)*Q.meta.tStar/86400];
        mfAll   = [mfAll; mf(:)*Q.meta.m0kg];
    end
end
if isempty(tfAll), return, end
fig = figure('Color', 'w', 'Position', [100 100 900 380]);
tl = tiledlayout(fig, 1, 2, 'TileSpacing', 'compact', 'Padding', 'compact');
ax1 = nexttile(tl);
scatter(ax1, rungAll, tfAll, 14, 'filled', 'MarkerFaceAlpha', 0.35, 'DisplayName', 'accepted entries');
set(ax1, 'XScale', 'log');  grid(ax1, 'on');
xlabel(ax1, 'thrust (N)');  ylabel(ax1, 't_f (days)');  title(ax1, 'transfer time per rung');
legend(ax1, 'Location', 'northeast');
ax2 = nexttile(tl);
scatter(ax2, rungAll, mfAll, 14, 'filled', 'MarkerFaceAlpha', 0.35, 'DisplayName', 'accepted entries');
set(ax2, 'XScale', 'log');  grid(ax2, 'on');
xlabel(ax2, 'thrust (N)');  ylabel(ax2, 'm_f (kg)');  title(ax2, 'final mass per rung');
legend(ax2, 'Location', 'southeast');
exportgraphics(fig, fullfile(picDir, 'catalog_ranges.png'), 'Resolution', 120);
close(fig);
end

function v = fdef(s, f, v0)
% FDEF  s.(f) if present else v0.  INPUTS: s [1 x 1] struct; f [1 x n]
% char; v0 [any].  OUTPUTS: v [any].
if isfield(s, f), v = s.(f); else, v = v0; end
end

function v = pick(c, a, b)
% PICK  Inline conditional.  INPUTS: c [1 x 1] logical; a, b [any].
% OUTPUTS: v [any] (a when c, else b).
if c, v = a; else, v = b; end
end

function logmsg(f, s)
% LOGMSG  Append-mode file log AND stdout.  INPUTS: f [1 x n] char log path
% ('' = stdout only); s [1 x m] char.  OUTPUTS: none.
fprintf('%s\n', s);
if ~isempty(f), fid = fopen(f, 'a'); fprintf(fid, '%s\n', s); fclose(fid); end
end
