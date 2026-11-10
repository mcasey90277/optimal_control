function chainBlockers = build_gto_dro_catalog(runIn)
%% Purpose:
%
%   THE CHAIN for the GTO -> DRO costate catalog AFTER the fleet
%   (run_gto_dro_catalog / run_gto_dro_batched.sh write the sheets; this
%   script turns sheets into a shippable catalog), build_70mN_library
%   style: parameter block first, one switch per stage, every stage
%   fenced, a blocker list the deliverable stage refuses on.
%
%     1. census           sheets found, cells solved per rung, retired cells
%     2. earth            earth_clearance_sweep per sheet, writeback (demotes)
%     3. repairUnsolved   repair_gto_dro_cells 'unsolved' (phase continuation)
%     4. repairBasin      repair_gto_dro_cells 'basin' (S4 on every top rung,
%                         reverse warm start on refuted ones)
%     2b. earth re-sweep  on every sheet a repair changed (the sweep is keyed
%                         to each entry's z8, so only re-solved entries fly)
%     5. package          build_costate_catalog_family on a CLEAN STAGING
%                         COPY of the official sheets (sheet key = orientDeg)
%     6. conj             conj_catalog_pass census -> writeback
%     7. gates            gates_catalog_pass census -> writeback
%     8. audit            audit_gto_dro_entries (flown, accepted, Earth floor)
%     9. pictures         coverage per sheet, t_f and m_f ranges per rung
%    10. deliverable      zip of catalog + pickers + README (refused on blockers)
%
%   Stages 2 to 4 modify the SHEETS (each keeps its own backup/sidecar);
%   stage 5 rebuilds the package from the sheets, and stages 6 to 7 run on
%   the package, so re-running the chain after a repair re-packages and
%   re-sweeps (the conj and gates sidecars and the conj writeback backup are
%   deleted when the package is rebuilt, since their entry indices would be
%   stale). With the package stage OFF the existing package and its sweep
%   sidecars are reused.
%
%   PACKAGE STAGING. The packager globs its folder, and the catalog folder
%   also holds <tag>_earth.mat, <tag>_repairs.mat, the repair tool's
%   <tag>_repair_backup_iDD_aDD.mat and the pilot's <tag>_deep.mat, all of
%   which match the sheet glob. The stage therefore copies only the
%   official sheets (explicit names, after complete_dep_params patched them
%   in place) into <outDir>/package_stage/ (emptied first) and packages
%   THAT folder. The packager records no source path in the catalog.
%
%   WAIVERS (deliberate decisions, off by default). A blocker is dropped from
%   the deliverable's refusal test only, never from the returned list, and
%   the number waived is logged:
%     allowConjBlockers     'conj: N refuted ...' (an extremal that is not a
%                           local minimum, shipped with its conj_pass verdict)
%     allowEarthDemotions   'earth-demoted: ...' (entries demoted to OK =
%                           false by the Earth sweep: absent from the
%                           catalog, an honest and visible state)
%   'conj incomplete', 'conj threw', 'earth-incomplete:' and 'earth threw'
%   are never waived, and neither is anything from the preflight.
%
%   DELIVERABLE PREFLIGHT. Whatever stages ran in THIS call, the deliverable
%   stage checks what is ON DISK, and each failure is a blocker that cannot
%   be waived: (a) 'catalog-missing:' the catalog file exists; (b) 'stale:'
%   the packaged catalog agrees with the current sheets by CONTENT: per
%   official sheet with accepted entries the packaged has_solution equals
%   the sheet's OK and the packaged z8 columns (picked through entry_index)
%   equal the sheet's Z8 at the OK entries; (c) 'unverified:' every packaged
%   sheet carries conj_pass, gate_min_lamv and gate_dimS and no accepted
%   entry is 'not run' (conj_pass = -1) or ungated, and every sheet's
%   _earth.mat sidecar covers the current z8 of every OK entry at or above
%   the floor; (d) the audit file exists, passed, and its fingerprint (the
%   tag, cell, rung and z8 of every replayed entry) still equals the sheets.
%   Refuted entries found on disk add the waivable 'conj:' blocker when this
%   call's conj stage did not already report one.
%
%   FLEET LOCK. A call with a sheet-writing stage on (earth, repairUnsolved,
%   repairBasin, package) takes the fleet driver lock
%   <catDir>/../.gtodrocat_driver.lock (mkdir, released by onCleanup) and
%   refuses with build_gto_dro_catalog:fleetRunning if it exists; it passes
%   ignoreLock = true to repair_gto_dro_cells because it already holds it.
%   runIn.ignoreLock = true skips the check and the take (a foreign lock is
%   never removed either way).
%
%% Inputs:
%
%  runIn                    struct (optional)       Stage switches
%                                                   .census .earth
%                                                   .repairUnsolved
%                                                   .repairBasin .package
%                                                   .conj .gates .audit
%                                                   .pictures .deliverable
%                                                   (all true); .catDir
%                                                   [../results/catalog],
%                                                   .outDir [../results],
%                                                   .expectSheets [16],
%                                                   .floorEarthKm [200],
%                                                   .maxRepair [inf],
%                                                   .allowConjBlockers
%                                                   [false],
%                                                   .allowEarthDemotions
%                                                   [false], .ignoreLock
%                                                   [false]
%
%% Outputs:
%
%  chainBlockers            cell                    Blocker strings (empty
%                                                   = clean); printed too
%
%% Revision History:
%  M. Casey                                                   (c) 11/10/2026
%  fix round 1: deliverable preflight, fleet lock, earth waiver, dated
%  sweep-backup rename, staging cleanup
%  Copyright Coorbital Inc.
%% ------------------------ Begin Code Sequence ---------------------------

if nargin < 1, runIn = struct(); end
here = fileparts(mfilename('fullpath'));
run(fullfile(fileparts(here), 'setup_gto_dro_paths.m'));
addpath(here);

%% ========================================================================
%  0. PARAMETERS AND STAGE SWITCHES
%% ========================================================================
run_ = struct('census', true, 'earth', true, 'repairUnsolved', true, 'repairBasin', true, ...
              'package', true, 'conj', true, 'gates', true, 'audit', true, 'pictures', true, ...
              'deliverable', true);
for f = fieldnames(run_).'
    if isfield(runIn, f{1}), run_.(f{1}) = runIn.(f{1}); end
end
          catDir = fdef(runIn, 'catDir', fullfile(fileparts(here), 'results', 'catalog'));
          outDir = fdef(runIn, 'outDir', fullfile(fileparts(here), 'results'));
    expectSheets = fdef(runIn, 'expectSheets', 16);
    floorEarthKm = fdef(runIn, 'floorEarthKm', 200);
       maxRepair = fdef(runIn, 'maxRepair', inf);
allowConjBlockers = fdef(runIn, 'allowConjBlockers', false);
allowEarthDemotions = fdef(runIn, 'allowEarthDemotions', false);
      ignoreLock = fdef(runIn, 'ignoreLock', false);
         orients = [0 90 180 270];
            taus = [0.5 1 2 3];
if ~isfolder(outDir), mkdir(outDir); end
logFile = fullfile(outDir, ['chain_' char(datetime('now', 'Format', 'yyyyMMdd_HHmm')) '.log']);
lg = @(varargin) logmsg(logFile, sprintf(varargin{:}));
files = struct('catalog', fullfile(outDir, 'costate_catalog_gto_dro.mat'), ...
               'audit', fullfile(outDir, 'audit_gto_dro.mat'), ...
               'pictures', fullfile(outDir, 'pictures'), ...
               'deliverable', fullfile(outDir, 'deliverable'), ...
               'stage', fullfile(outDir, 'package_stage'));
blockers = {};

% the sheets, by their official names (never a bare glob: _deep, _earth,
% _repairs, _repair_backup_* and .bak_* must not enter)
sheetMats = {};
for od = orients
    for ta = taus
        p = fullfile(catDir, sprintf('gto_dro_o%03d_tau%g.mat', od, ta));
        if isfile(p), sheetMats{end+1} = p; end
    end
end
nS = numel(sheetMats);

% the fleet lock, for any call that writes sheets (same directory, same
% atomic mkdir as repair_gto_dro_cells and run_gto_dro_batched.sh)
lockPath = fullfile(fileparts(catDir), '.gtodrocat_driver.lock');
if ~ignoreLock && (run_.earth || run_.repairUnsolved || run_.repairBasin || run_.package)
    if isfolder(lockPath) || isfile(lockPath)
        error('build_gto_dro_catalog:fleetRunning', ...
              'fleet driver lock %s exists; a fleet or a repair is running (runIn.ignoreLock overrides)', lockPath);
    end
    if ~isfolder(fileparts(lockPath)), mkdir(fileparts(lockPath)); end
    [lockMade, ~, lockId] = mkdir(lockPath);
    if ~lockMade || strcmp(lockId, 'MATLAB:MKDIR:DirectoryExists')
        error('build_gto_dro_catalog:fleetRunning', 'could not take the fleet driver lock %s', lockPath);
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
            Q = load(sheetMats{ks}, 'OK', 'ATT', 'rungs');
            [nD, nA, ~] = size(Q.OK);
            solved = squeeze(sum(sum(Q.OK, 1), 2)).';
            retired = nnz(~any(Q.OK, 3) & Q.ATT(:, :, 1) >= 2);
            [~, nm] = fileparts(sheetMats{ks});
            lg('   %-22s cells %3d  per rung [%s]  retired %2d', nm, nD*nA, num2str(solved), retired);
        end
    catch ME
        lg('1. census THREW: %s', ME.message);  blockers{end+1} = ['census threw: ' ME.message];
    end
end

%% ========================================================================
%  2. EARTH CLEARANCE SWEEP (sheets; writeback demotes)
%% ========================================================================
earthBlk = repmat({''}, 1, nS);
if run_.earth
    for ks = 1:nS
        earthBlk{ks} = sweep_sheet(sheetMats{ks}, floorEarthKm, logFile, '2. earth');
    end
end

%% ========================================================================
%  3./4. REPAIRS (sheets; sidecars)
%% ========================================================================
plan = [];
changed = false(1, nS);
for mode = {'unsolved', 'basin'}
    sw = pick(strcmp(mode{1}, 'unsolved'), run_.repairUnsolved, run_.repairBasin);
    if ~sw, continue, end
    for ks = 1:nS
        try
            [~, nm] = fileparts(sheetMats{ks});
            if isempty(plan)
                plan = run_gto_dro_catalog([], [], [], [], struct('planOnly', true));
            end
            kp = find(strcmp({plan.tag}, nm), 1);
            assert(~isempty(kp), 'build_gto_dro_catalog:plan', 'no plan row for %s', nm);
            R = repair_gto_dro_cells(sheetMats{ks}, mode{1}, ...
                    struct('engineOpts', plan(kp).optsA, 'maxCells', maxRepair, 'logFile', logFile, ...
                           'ignoreLock', true));
            lg(['3/4. repair %-8s %-22s candidates %d attempted %d repaired %d restored %d ', ...
                'failed %d refuted %d recovered %d'], mode{1}, nm, R.nCandidates, R.nAttempted, ...
               R.nRepaired, R.nRestored, R.nFailed, R.nConjRefuted, R.nRecovered);
            if R.nRepaired > 0 || R.nRecovered > 0, changed(ks) = true; end
            if R.nFailed > 0, blockers{end+1} = sprintf('repair %s: %s %d failed', mode{1}, nm, R.nFailed); end
        catch ME
            lg('3/4. repair %s THREW on %s: %s', mode{1}, sheetMats{ks}, ME.message);
            blockers{end+1} = ['repair threw: ' ME.message];
        end
    end
end

%% ========================================================================
%  2b. EARTH RE-SWEEP on every sheet a repair changed
%% ========================================================================
if run_.earth
    for ks = find(changed)
        earthBlk{ks} = sweep_sheet(sheetMats{ks}, floorEarthKm, logFile, '2b. earth re-sweep');
    end
    blockers = [blockers, earthBlk(~cellfun(@isempty, earthBlk))];
end

%% ========================================================================
%  5. PACKAGE (sheet key = orientDeg; dep_params completed on disk; clean
%     staging copy of the official sheets only)
%% ========================================================================
if run_.package
    try
        complete_dep_params(sheetMats);
        spec = struct( ...
          'glob', 'gto_dro_o*_tau*.mat', ...
          'name', 'costate_catalog_gto_dro', ...
          'description', ['Minimum-time low-thrust GTO->DRO costate CATALOG: one (GTO orientation ', ...
            'x DRO period) phasing sheet per pair, 12 x 6 phasing, thrust rungs walked from 15 N. ', ...
            'COMPACT format: canonical ND quantities only; derived values via .derive. Each entry ', ...
            'z8 = [lambda_r(3); lambda_v(3); lambda_m; tf] is a converged PMP solution accepted by ', ...
            'pumpkyn.cr3bp.tfMin, solved with a 200 km Earth-altitude floor and swept for the ', ...
            'flown Earth clearance. Departure is the GTO PSEUDO-family (algebraic Kepler locus, ', ...
            'see get_family_orbit), 350 x 35786 km, orientDeg = argument of perigee from the ', ...
            'Earth->Moon line, prograde; sheets(k).tauDRO/.tau_dep hold orientDeg IN DEGREES (the ', ...
            'GTO sheet key rule), the DRO period is arr_params.tau (ND). conj_pass marks entries ', ...
            'the conjugate test refutes: consult it, a refuted entry is an extremal, not a minimum.'], ...
          'provenance', ['GTO -> DRO catalog run (run_gto_dro_catalog), 2026-11; pilot sheet ', ...
            '(0 deg, tau 1) then the fleet; two-point top-rung multistart, first-success-wins; ', ...
            'Earth sweep, phase-continuation and reverse-warm-start repairs recorded per sheet ', ...
            'in _earth.mat / _repairs.mat sidecars. Casey/Koblick 2026-11.'], ...
          'depReconstruction', ['gto: [tau,rv] = get_family_orbit(''gto'', dep_params) -> algebraic ', ...
            'Kepler locus (NO propagation); state at phase f: phase_state(tau, rv)(f); dro: ', ...
            '[tau,rv] = get_family_orbit(''dro'', arr_params) -> getDRO seed, cont_np, prop; ', ...
            'state at phase f: phase_state(tau, rv)(f)']);
        % the sweep sidecars index entries of the OLD package: delete them
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
        % clean staging folder: ONLY the official sheets
        if isfolder(files.stage), rmdir(files.stage, 's'); end
        mkdir(files.stage);
        for ks = 1:nS, copyfile(sheetMats{ks}, files.stage); end
        cat_ = build_costate_catalog_family(files.stage, files.catalog, spec);
        lg('5. package: %d sheets, %d entries -> %s', ...
           numel(cat_.sheets), count_entries(cat_), files.catalog);
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
%  8. AUDIT (replay with the Earth floor as a gate)
%% ========================================================================
if run_.audit
    try
        results = audit_gto_dro_entries(catDir, sheetMats, floorEarthKm);
        save(files.audit, 'results');
        lg('8. audit (%s rule): fly %d/%d, accept %d/%d, earth %d/%d', results(1).rule, ...
           nnz([results.passFly]), numel(results), nnz([results.passAccept]), numel(results), ...
           nnz([results.passEarth]), numel(results));
        if ~all([results.passFly] & [results.passAccept] & [results.passEarth])
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
%  10. DELIVERABLE (preflight on what is on disk; refused on any blocker but
%      a pictures failure and the waived conj: / earth-demoted: blockers)
%% ========================================================================
if run_.deliverable
    try
        [pre, nRefuted] = deliverable_preflight(files, sheetMats, floorEarthKm, expectSheets);
    catch ME
        pre = {['unverified: preflight threw: ' ME.message]};  nRefuted = 0;
    end
    % gates:, census: and earth-demoted: are also reported by their own
    % stages; list the preflight's re-derivation only when no blocker with
    % that prefix stands already
    keep = true(size(pre));
    for kb = 1:numel(pre)
        pfx = pre{kb}(1:find(pre{kb} == ':', 1));
        if any(strcmp(pfx, {'gates:', 'census:', 'earth-demoted:'})) && any(startsWith(blockers, pfx))
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
shipBlockers = blockers(~startsWith(blockers, 'pictures threw'));
for wv = {{allowConjBlockers, 'conj:', 'allowConjBlockers'}, ...
          {allowEarthDemotions, 'earth-demoted:', 'allowEarthDemotions'}}
    if wv{1}{1}
        isW = startsWith(shipBlockers, wv{1}{2});
        if any(isW)
            lg('10. %s: %d blocker(s) waived for the deliverable', wv{1}{3}, nnz(isW));
        end
        shipBlockers = shipBlockers(~isW);
    end
end
if run_.deliverable
    if ~isempty(shipBlockers)
        lg('10. deliverable REFUSED -- %d blocker(s):', numel(shipBlockers));
        for kb = 1:numel(shipBlockers), lg('     - %s', shipBlockers{kb}); end
    else
        try
            if ~isfolder(files.deliverable), mkdir(files.deliverable); end
            stamp = char(datetime('now', 'Format', 'yyyyMMdd'));
            stage = fullfile(files.deliverable, ['costate_catalog_gto_dro_' stamp]);
            if isfolder(stage), rmdir(stage, 's'); end
            mkdir(stage);
            copyfile(files.catalog, stage);
            copyfile(fullfile(here, 'costate_catalog_pick.m'), stage);
            copyfile(fullfile(here, 'costate_lib_describe.m'), stage);
            copyfile(fullfile(fileparts(here), 'README.md'), fullfile(stage, 'README.md'));
            zipPath = [stage '.zip'];
            zip(zipPath, stage);
            lg('10. deliverable: %s', zipPath);
        catch ME
            lg('10. deliverable THREW: %s', ME.message);  blockers{end+1} = ['deliverable threw: ' ME.message];
        end
    end
end

%% ========================================================================
%  BLOCKER SUMMARY
%% ========================================================================
if isempty(blockers)
    lg('CHAIN CLEAN: every stage that ran, passed.');
else
    lg('CHAIN COMPLETE WITH %d BLOCKER(S):', numel(blockers));
    for kb = 1:numel(blockers), lg('  - %s', blockers{kb}); end
    lg('Nothing ships until these clear; the measurements above are still valid.');
end
chainBlockers = blockers;
end

% ------------------------------------------------------------------------
function [blk, nRefuted] = deliverable_preflight(files, sheetMats, floorEarthKm, expectSheets)
% DELIVERABLE_PREFLIGHT  What is ON DISK, independent of the stages that ran
% this call: catalog present, consistent in CONTENT with the current sheets
% (packaged has_solution == sheet OK; packaged z8(:, entry_index) == sheet
% Z8 at the OK entries), verdict and gate records complete, Earth sweep
% current, audit present and fingerprint-current; and the verdicts that
% the switchable stages also report are RE-DERIVED from disk so that
% switching a stage off is not an implicit waiver: gate failures (gates:,
% counted as gates_catalog_pass counts them, tolerance from
% hyp_gates.lamVTol), a sheet shortfall (census:), Earth demotions
% (earth-demoted:, from each sheet's meta.earthSweep.nDemoted).  INPUTS:
% files struct (.catalog, .audit); sheetMats cell; floorEarthKm double;
% expectSheets double.  OUTPUTS: blk cell of blocker strings (prefixes
% catalog-missing: stale: unverified: gates: census: earth-demoted:);
% nRefuted double (accepted entries with conj_pass = 0 in the package).
blk = {};  nRefuted = 0;
if ~isfile(files.catalog)
    blk{end+1} = sprintf('catalog-missing: %s does not exist', files.catalog);
    nDem0 = 0;                           % demotions are visible without a catalog too
    for ks = 1:numel(sheetMats)
        Q = load(sheetMats{ks}, 'meta');
        if isfield(Q.meta, 'earthSweep'), nDem0 = nDem0 + Q.meta.earthSweep.nDemoted; end
    end
    if nDem0 > 0
        blk{end+1} = sprintf('earth-demoted: %d entries demoted by the Earth sweep (meta.earthSweep; OK = false)', nDem0);
    end
    return
end
L = load(files.catalog);  fn = fieldnames(L);  cat_ = L.(fn{1});
nonEmpty = 0;  nStaleSheet = 0;  nMissingRec = 0;  nNoConj = 0;  nNoGate = 0;  nEarth = 0;
nDemoted = 0;  nH2 = 0;  nH3 = 0;  nAbn = 0;
lamVTol = 1e-6;                          % gates_catalog_pass default, used if the catalog does not record it
if isfield(cat_, 'hyp_gates') && isfield(cat_.hyp_gates, 'lamVTol'), lamVTol = cat_.hyp_gates.lamVTol; end
for ks = 1:numel(sheetMats)
    Q = load(sheetMats{ks}, 'OK', 'Z8', 'meta');
    if isfield(Q.meta, 'earthSweep'), nDemoted = nDemoted + Q.meta.earthSweep.nDemoted; end
    if nnz(Q.OK) == 0, continue, end
    nonEmpty = nonEmpty + 1;
    % the Earth sweep sidecar must cover the current z8 of every OK entry
    sideMat = [sheetMats{ks}(1:end-4) '_earth.mat'];
    okEarth = isfile(sideMat);
    if okEarth
        Se = load(sideMat);
        m = Q.OK(:);
        A = reshape(Se.Z8SWEPT, 8, []);  B = reshape(Q.Z8, 8, []);
        okEarth = all(Se.done(m)) && isequaln(A(:, m), B(:, m)) && all(isfinite(Se.EARTHKM(m))) ...
                  && all(Se.EARTHKM(m) >= floorEarthKm);
    end
    if ~okEarth, nEarth = nEarth + 1; end
    % the packaged twin of this sheet, by key
    kp = find([cat_.sheets.tau_dep] == Q.meta.tauDRO & [cat_.sheets.tau_arr] == Q.meta.arrParams.tau, 1);
    if isempty(kp), nStaleSheet = nStaleSheet + 1; continue, end
    sh = cat_.sheets(kp);
    same = isequal(sh.has_solution, Q.OK);
    if same
        Z = reshape(Q.Z8, 8, []);
        same = isequal(sh.z8(:, sh.entry_index(Q.OK)), Z(:, find(Q.OK(:))));
    end
    if ~same, nStaleSheet = nStaleSheet + 1; end
    if ~(isfield(sh, 'conj_pass') && isfield(sh, 'gate_min_lamv') && isfield(sh, 'gate_dimS'))
        nMissingRec = nMissingRec + 1;  continue
    end
    nNoConj = nNoConj + nnz(sh.has_solution & sh.conj_pass == -1);
    nNoGate = nNoGate + nnz(sh.has_solution & isnan(sh.gate_min_lamv));
    nRefuted = nRefuted + nnz(sh.has_solution & sh.conj_pass == 0);
    % the hypothesis gates, counted as gates_catalog_pass's census counts them
    dn = sh.has_solution & ~isnan(sh.gate_min_lamv);
    nH2 = nH2 + nnz(dn & sh.gate_min_lamv <= lamVTol);
    nH3 = nH3 + nnz(dn & sh.gate_min_qmt <= 0);
    nAbn = nAbn + nnz(dn & sh.gate_dimS ~= 1);
end
if nH2 + nH3 + nAbn > 0
    blk{end+1} = sprintf('gates: %d entries fail a hypothesis gate (on disk: H2 %d, H3 %d, abnormal %d)', nH2 + nH3 + nAbn, nH2, nH3, nAbn);
end
if min(numel(cat_.sheets), numel(sheetMats)) < expectSheets
    blk{end+1} = sprintf('census: %d packaged / %d official sheets on disk, %d expected', numel(cat_.sheets), numel(sheetMats), expectSheets);
end
if nDemoted > 0
    blk{end+1} = sprintf('earth-demoted: %d entries demoted by the Earth sweep (meta.earthSweep; OK = false, absent from the catalog)', nDemoted);
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
if nEarth > 0
    blk{end+1} = sprintf('unverified: %d sheet(s) without a current Earth sweep covering every accepted entry', nEarth);
end
% the audit: present, passed, and replayed THESE sheets
if ~isfile(files.audit)
    blk{end+1} = sprintf('unverified: audit file %s does not exist', files.audit);
else
    Au = load(files.audit);
    okAudit = isfield(Au, 'results') && ~isempty(Au.results) && isfield(Au.results, 'z8') ...
              && all([Au.results.passFly] & [Au.results.passAccept] & [Au.results.passEarth]);
    if okAudit
        names = cell(size(sheetMats));
        for ks = 1:numel(sheetMats), [~, names{ks}] = fileparts(sheetMats{ks}); end
        for kr = 1:numel(Au.results)
            r = Au.results(kr);
            ks = find(strcmp(names, r.sheet), 1);
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

% ------------------------------------------------------------------------
function blk = sweep_sheet(sheetMat, floorKm, logFile, label)
% SWEEP_SHEET  Earth sweep (census, then writeback if complete) on one
% sheet, fenced.  INPUTS: sheetMat char; floorKm double; logFile char; label char.
% OUTPUTS: blk char - blocker string, '' when clean.
blk = '';
[~, nm] = fileparts(sheetMat);
try
    S = earth_clearance_sweep(sheetMat, struct('floorKm', floorKm, 'logFile', logFile));
    nFlownNow = S.nFlown;
    if S.done
        S = earth_clearance_sweep(sheetMat, struct('floorKm', floorKm, 'writeback', true, 'logFile', logFile));
        logmsg(logFile, sprintf('%s %-22s min %.0f km, %d flown now, %d demoted (record), %d non-finite', ...
               label, nm, S.minAltKm, nFlownNow, S.nBelow, S.nNonFinite));
        if S.nBelow > 0, blk = sprintf('earth-demoted: %s demoted %d entries (OK = false, absent from the catalog)', nm, S.nBelow); end
    else
        blk = sprintf('earth-incomplete: %s incomplete (%d to do)', nm, S.nTodo);
        logmsg(logFile, sprintf('%s %s', label, blk));
    end
catch ME
    logmsg(logFile, sprintf('%s THREW on %s: %s', label, nm, ME.message));
    blk = ['earth threw: ' ME.message];
end
end

function complete_dep_params(sheetMats)
% COMPLETE_DEP_PARAMS  Patch each sheet's meta.depParams to the full
% {sma_km, ecc, orientDeg} triple, idempotently and atomically (copy of
% GTO_tulip/catalog/build_gto_catalog.m's subfunction, by explicit file
% list instead of a glob).  INPUTS: sheetMats cell.  OUTPUTS: none.
for kf = 1:numel(sheetMats)
    fpath = sheetMats{kf};
    S = load(fpath);
    dp = S.meta.depParams;
    if isfield(dp, 'sma_km') && isfield(dp, 'ecc'), continue, end
    assert(strcmpi(S.meta.depFamily, 'gto'), 'complete_dep_params:family', ...
        'sheet %s: depFamily is ''%s'', expected ''gto''', fpath, S.meta.depFamily);
    [~, ~, info] = get_family_orbit('gto', struct('orientDeg', dp.orientDeg));
    S.meta.depParams = struct('sma_km', info.sma_km, 'ecc', info.ecc, 'orientDeg', info.orientDeg);
    assert(abs(info.orientDeg - S.meta.tauDRO) < 1e-12, 'complete_dep_params:key', ...
        'sheet %s: orientDeg (%.6g) ~= tauDRO key (%.6g)', fpath, info.orientDeg, S.meta.tauDRO);
    tmpPath = [tempname(fileparts(fpath)) '.mat'];
    save(tmpPath, '-struct', 'S');
    Schk = load(tmpPath, 'meta');
    okTmp = isfield(Schk.meta.depParams, 'sma_km') && isfield(Schk.meta.depParams, 'ecc');
    if ~okTmp
        if isfile(tmpPath), delete(tmpPath); end
        error('complete_dep_params:tmpVerify', 'sheet %s: temp file failed verification', fpath);
    end
    movefile(tmpPath, fpath);
    fprintf('complete_dep_params: patched %s -> sma_km=%.6g ecc=%.6g orientDeg=%g\n', ...
        fpath, info.sma_km, info.ecc, info.orientDeg);
end
end

function n = count_entries(cat_)
% COUNT_ENTRIES  Accepted entries across sheets.  INPUTS: cat_.  OUTPUTS: n.
n = 0;
for ks = 1:numel(cat_.sheets), n = n + nnz(cat_.sheets(ks).has_solution); end
end

function catalog_pictures(sheetMats, picDir)
% CATALOG_PICTURES  One coverage figure per sheet (OK per rung on the
% nD x nA torus) and one summary figure (t_f and m_f ranges per rung over
% every sheet).  INPUTS: sheetMats cell; picDir char.  OUTPUTS: none (PNGs).
nS = numel(sheetMats);
rungAll = [];  tfAll = [];  mfAll = [];
for ks = 1:nS
    Q = load(sheetMats{ks}, 'OK', 'TF', 'rungs', 'meta', 'sD', 'sA');
    [nD, nA, nR] = size(Q.OK);
    [~, nm] = fileparts(sheetMats{ks});
    fig = figure('Color', 'w', 'Position', [100 100 1100 320]);
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
% FDEF  s.(f) if present else v0.  INPUTS: s; f; v0.  OUTPUTS: v.
if isfield(s, f), v = s.(f); else, v = v0; end
end

function v = pick(c, a, b)
% PICK  Inline conditional.  INPUTS: c; a; b.  OUTPUTS: v.
if c, v = a; else, v = b; end
end

function logmsg(f, s)
% LOGMSG  Append-mode file log AND stdout.  INPUTS: f; s.  OUTPUTS: none.
fprintf('%s\n', s);
if ~isempty(f), fid = fopen(f, 'a'); fprintf(fid, '%s\n', s); fclose(fid); end
end
