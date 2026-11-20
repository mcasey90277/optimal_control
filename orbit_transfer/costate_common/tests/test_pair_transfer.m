function [ok, wall] = test_pair_transfer(cfg)
%% Purpose:
%
%   Test of pair_transfer, the shared one-transfer solver and gate for any
%   orbit pair. Sections (run all by default, or name a subset):
%
%     'plan'    planOnly rung lists: 5 N -> [15 12 10 7 5]; 0.07 N -> the
%               eleven standard rungs then 0.07; 15 N -> [15]; opts.rungs
%               [15 5] used as given; a rung list not ending at thrustN and
%               an unknown option refuse by name.
%     'fail'    deliberate failures that must not throw: maxIter = 1 on the
%               ELFO cell gives ok = false, reason 'top rung failed';
%               msHardCapSec = -1 makes the ENGINE throw (its own option
%               check), caught as reason 'engine error: ...' with the
%               scratch folder removed. A bad family name errors BY DESIGN
%               before the engine (get_family_orbit:family). The mass guard
%               ('mass' mode) gives 'mass exhausted' at an impossible t_f.
%     'elfo'    DRO tau 1 -> ELFO raanDeg 0, sD 0, sA 0, 15 N: ok, t_f within
%               1e-6 of the stage 0 value 0.0900542599967 ND, flown miss
%               < 1e-3 km, gates.pass, located Moon minimum >= 500 km;
%               INDEPENDENTLY of located_min_altitude, the Moon distance of
%               every flight sample recomputed here: the located minimum is
%               <= that sample minimum + 1e-6 km and within 5 km of it, and
%               moonMinT lies in [0, t_f]; junctions [14 x 25], flight end
%               state within 1e-6 ND of rvf, flight.Y [M x 14], unit alpha.
%     'regate'  every blocking gate fails on the stored ELFO result through
%               pair_transfer('regate', T, th), without solving: floorKm
%               1e5 (moonMinKm < floorKm), gateKm 1e-12, accTol 0,
%               floorEarthKm 1e7, and tampered ok / admissible / moon
%               located-minimum failure; each regate is recorded in
%               entry.provenance.regate with the thresholds the verdicts
%               carry; an unknown threshold refuses.
%     'cap'     the second-order fence on the stored ELFO result
%               (pair_transfer('second', T, o), no ladder solve): cap 0.5 s
%               with the 5 s secondDelaySec seam gives secondStatus
%               'capped', conjPass NaN, conjReason 'capped after ...';
%               cap 300 s (fenced) and cap 0 (in-process) both reproduce
%               the stored conjPass; the pool state after each call equals
%               the state before it (no pool leaked, a pool already open is
%               kept).
%     'equiv'   GTO orient 0 -> DRO tau 1, sD 0.10, sA 0, 15 N, floorEarthKm
%               200 (gto_dro_entry's recipe and demo cell; directGateKm 100
%               as the wrapper's engine default): z8 equal to gto_dro_entry's
%               demo entry.z8 to 1e-6. GTO_DRO/lib is put on the path for
%               that call only, and the path after the call must equal the
%               path before it ([SKIP] when the wrapper result came from
%               the cache, since nothing was called).
%     'tulip'   DRO tau 1 -> tulip Np 7 pm -1, 15 N, at a cell of the
%               shipped catalog's OWN grid that has a 15 N entry (DRO_tulip/
%               direct/results/costate_catalog_dro_tulip.mat, tauDRO == 1
%               sheet, read only; the solved cell nearest (0, 0.075378)):
%               gates.pass, t_f printed beside the catalog's, no equality
%               asserted. Without the file: sD 0, sA 0.075378 and a skip
%               note.
%
%   CACHE: solved results are cached in tempdir/pair_transfer_test_cache,
%   keyed by a SHA-256 of (dep, arr, sD, sA, thrustN, opts, cacheVersion)
%   AND of the source text of pair_transfer, located_min_altitude,
%   thrust_ladder_library, fly_transfer and ms_tfmin (plus gto_dro_entry
%   for the equivalence file), so editing any of them re-solves. Use
%   cfg.useCache = false after a change to anything else the solve calls
%   (casadi_mintime_dro, harvest_ms_seed, pumpkyn, CasADi, ...), before a
%   commit, and whenever a result must be measured rather than recalled.
%   The uncached run takes several minutes (ELFO about 210 s, equivalence
%   about 200 s, tulip about 90 s): run it in the background.
%
%   PATH AND POOL: the path is snapshotted and restored on exit; the test
%   adds costate_common, oclib and CasADi, and runs pumpkynPie's startup
%   only when pumpkyn is missing. The test opens no pool itself:
%   pair_transfer opens and deletes its own fence pool when none is open.
%
%% Inputs:
%
%  cfg                      cellstr or struct       sections (cellstr), or
%                           (optional)              .sections [all]
%                                                   .useCache [true]
%
%% Outputs:
%
%  ok                       logical                 true iff every check
%                                                   passed
%
%  wall                     struct                  wall time per section, s
%
%% Revision History:
%  M. Casey                                                   (c) 11/19/2026
%  M. Casey  review round 1: source-hashed cache, useCache,    (c) 11/19/2026
%            cap / engine-error / mass sections, independent
%            Moon check, catalog-grid tulip cell
%  Copyright Coorbital Inc.
%% ------------------------ Begin Code Sequence ---------------------------
allSec = {'plan', 'fail', 'elfo', 'regate', 'cap', 'equiv', 'tulip'};
if nargin < 1 || isempty(cfg), cfg = struct(); end
if ischar(cfg), cfg = {cfg}; end
if iscell(cfg), cfg = struct('sections', {cfg}); end
sections = allSec;
if isfield(cfg, 'sections') && ~isempty(cfg.sections), sections = cfg.sections; end
useCache = ~isfield(cfg, 'useCache') || cfg.useCache;
cacheVersion = 'v2';

here    = fileparts(mfilename('fullpath'));          % costate_common/tests
cc      = fileparts(here);                           % costate_common
ot      = fileparts(cc);                             % orbit_transfer
home    = getenv('HOME');
oldPath = path;
restorePath = onCleanup(@() path(oldPath));
addpath(cc, fullfile(fileparts(ot), 'oclib'), fullfile(home, 'casadi-3.7.0'));
if isempty(which('pumpkyn.cr3bp.tfMin'))
    od = cd(fullfile(home, 'Desktop', 'proj7', 'external', 'pumpkynPie'));  startup();  cd(od);
end
cacheDir = fullfile(tempdir, 'pair_transfer_test_cache');
if ~isfolder(cacheDir), mkdir(cacheDir); end
srcHash = source_hash({'pair_transfer', 'located_min_altitude', 'thrust_ladder_library', ...
                       'fly_transfer', 'ms_tfmin'});
C = struct('dir', cacheDir, 'ver', [cacheVersion '_' srcHash], 'use', useCache);
fprintf('  (cache %s, key suffix %s)\n', onoff(useCache), C.ver);

ok = true;
wall = struct();
mu    = 0.012150585609624;
lStar = 389703.264829278;
dro1  = struct('family', 'dro',  'params', struct('tau', 1));
elfo0 = struct('family', 'elfo', 'params', struct('raanDeg', 0));

%% plan
if ismember('plan', sections)
    t0 = tic;
    po = struct('planOnly', true);
    T5 = pair_transfer(dro1, elfo0, 0, 0, 5, po);
    ok = check(sprintf('5 N rungs %s', mat2str(T5.rungs)), isequal(T5.rungs, [15 12 10 7 5])) && ok;
    T7 = pair_transfer(dro1, elfo0, 0, 0, 0.07, po);
    ok = check(sprintf('0.07 N rungs %s', mat2str(T7.rungs)), ...
               isequal(T7.rungs, [15 12 10 7 5 3 2 1.5 1 0.75 0.5 0.07])) && ok;
    T15 = pair_transfer(dro1, elfo0, 0, 0, 15, po);
    ok = check(sprintf('15 N rungs %s', mat2str(T15.rungs)), isequal(T15.rungs, 15)) && ok;
    po.rungs = [15 5];
    Tg = pair_transfer(dro1, elfo0, 0, 0, 5, po);
    ok = check(sprintf('opts.rungs [15 5] used as given: %s, engine rungs %s', mat2str(Tg.rungs), ...
               mat2str(Tg.engineOpts.rungs)), isequal(Tg.rungs, [15 5]) && ...
               isequal(Tg.engineOpts.rungs, [15 5])) && ok;
    ok = check('planOnly: ok false, reason ''planOnly'', recipe in engineOpts', ~T5.ok && ...
               strcmp(T5.reason, 'planOnly') && T5.engineOpts.thrLock == true && ...
               T5.engineOpts.tf0 == 0.3 && T5.engineOpts.maxIter == 6000 && ...
               T5.engineOpts.msHardCapSec == 300 && strcmp(T5.engineOpts.arrFamily, 'elfo')) && ok;
    ok = check_error('rungs [15 7] at 5 N refuses', 'pair_transfer:rungs', ...
                     @() pair_transfer(dro1, elfo0, 0, 0, 5, struct('planOnly', true, 'rungs', [15 7]))) && ok;
    ok = check_error('unknown option refuses', 'pair_transfer:opts', ...
                     @() pair_transfer(dro1, elfo0, 0, 0, 5, struct('planOnly', true, 'flooorKm', 1))) && ok;
    wall.plan = toc(t0);
    fprintf('  (plan: %.1f s)\n', wall.plan);
end

%% fail
if ismember('fail', sections)
    t0 = tic;
    try
        Tf = pair_transfer(dro1, elfo0, 0, 0, 15, struct('maxIter', 1, 'quiet', true));
        ok = check(sprintf('maxIter 1: ok %d, reason ''%s'', gates.pass %d', Tf.ok, Tf.reason, Tf.gates.pass), ...
                   ~Tf.ok && strcmp(Tf.reason, 'top rung failed') && ~Tf.gates.pass && ...
                   ~Tf.okPerRung(1)) && ok;
        ok = check('failed solve keeps the engine log for diagnosis', ...
                   ~isempty(Tf.scratch.log) && isfile(Tf.scratch.log)) && ok;
        if isfolder(Tf.scratch.outDir), rmdir(Tf.scratch.outDir, 's'); end
    catch err
        ok = check(['maxIter 1 threw: ' err.message], false) && ok;
    end
    try
        Tx = pair_transfer(dro1, elfo0, 0, 0, 15, struct('msHardCapSec', -1, 'quiet', true));
        ok = check(sprintf('engine throw caught: ok %d, reason ''%s''', Tx.ok, Tx.reason), ...
                   ~Tx.ok && startsWith(Tx.reason, 'engine error: ') && ...
                   contains(Tx.reason, 'msHardCapSec') && ~Tx.gates.pass) && ok;
        ok = check('engine throw: scratch sheet and the folder it made are gone', ...
                   ~isfile(Tx.scratch.sheet) && ~isfolder(Tx.scratch.outDir)) && ok;
        if isfolder(Tx.scratch.outDir), rmdir(Tx.scratch.outDir, 's'); end
    catch err
        ok = check(['engine throw escaped: ' err.message], false) && ok;
    end
    bad = struct('family', 'nofamily', 'params', struct());
    ok = check_error('bad family errors by design before the engine', 'get_family_orbit:family', ...
                     @() pair_transfer(dro1, bad, 0, 0, 15, struct('quiet', true))) && ok;
    M0 = pair_transfer('mass', 1e3, 15, struct());
    M1 = pair_transfer('mass', 0.0900542599967, 15, struct());
    ok = check(sprintf('mass guard: t_f 1e3 ND -> ok %d ''%s'' (mfFrac %.3g); t_f 0.0901 -> m_f %.4f kg', ...
               M0.ok, M0.reason, M0.mfFrac, M1.mf_kg), ~M0.ok && strcmp(M0.reason, 'mass exhausted') && ...
               isnan(M0.mf_kg) && M1.ok && abs(M1.mf_kg - 119.1499) < 1e-3) && ok;
    wall.fail = toc(t0);
    fprintf('  (fail: %.1f s)\n', wall.fail);
end

%% elfo
Te = [];
if any(ismember({'elfo', 'regate', 'cap'}, sections))
    t0 = tic;
    Te = cached_solve(C, dro1, elfo0, 0, 0, 15, struct());
    g = Te.gates;
    fprintf(['  ELFO 15 N: ok %d, t_f %.13f ND (%.4f d), m_f %.4f kg, dV %.4f km/s; flown %.3e km %.3e m/s, ' ...
             'accDz %.1e, |R| %.1e, admissible %d\n'], Te.ok, Te.tf_nd, Te.tf_d, Te.mf_kg, Te.dV_kms, ...
            g.flownKm, g.flownMps, g.accDz, g.msNormR, g.admissible);
    fprintf(['  ELFO 15 N: Moon located %.4f km at t %.6f (sample %.4f km), Earth located %.1f km; ' ...
             'conjPass %g (%s); minLamV %.4g minQmt %.4g dimS %g; second %s\n'], g.moonMinKm, ...
            g.moonMinT, g.moonSampleMinKm, g.earthMinKm, g.conjPass, g.conjReason, g.minLamV, ...
            g.minQmt, g.dimS, g.secondStatus);
    if ismember('elfo', sections)
        ok = check('ok', Te.ok) && ok;
        ok = check(sprintf('t_f %.13f within 1e-6 of 0.0900542599967 (diff %.2e)', Te.tf_nd, ...
                   Te.tf_nd - 0.0900542599967), abs(Te.tf_nd - 0.0900542599967) < 1e-6) && ok;
        ok = check(sprintf('flown %.3e km < 1e-3', g.flownKm), g.flownKm < 1e-3) && ok;
        ok = check(sprintf('gates.pass (failed: %s)', strjoin(g.failed, ', ')), g.pass) && ok;
        ok = check(sprintf('Moon located minimum %.4f km >= 500', g.moonMinKm), g.moonMinKm >= 500) && ok;
        hMoon = lStar*vecnorm(Te.flight.Y(:, 1:3) - [1 - mu 0 0], 2, 2) - 1737.4;   % independent
        sMin  = min(hMoon);
        ok = check(sprintf('independent sample minimum %.4f km: located %.4f km <= it + 1e-6 and within 5 km', ...
                   sMin, g.moonMinKm), g.moonMinKm <= sMin + 1e-6 && g.moonMinKm >= sMin - 5) && ok;
        ok = check(sprintf('moonMinT %.6f in [0, t_f %.6f]', g.moonMinT, Te.tf_nd), ...
                   g.moonMinT >= 0 && g.moonMinT <= Te.tf_nd) && ok;
        ok = check(sprintf('junctions.yj %s, tj %s', mat2str(size(Te.junctions.yj)), ...
                   mat2str(size(Te.junctions.tj))), isequal(size(Te.junctions.yj), [14 25]) && ...
                   isequal(size(Te.junctions.tj), [1 25])) && ok;
        dEnd = max(abs(Te.flight.Y(end, 1:6) - Te.rvf));
        ok = check(sprintf('flight end state within 1e-6 ND of rvf (%.2e)', dEnd), dEnd < 1e-6) && ok;
        ok = check(sprintf('flight.Y %s, alpha %s unit (%.1e)', mat2str(size(Te.flight.Y)), ...
                   mat2str(size(Te.flight.alpha)), max(abs(vecnorm(Te.flight.alpha, 2, 2) - 1))), ...
                   size(Te.flight.Y, 2) == 14 && isequal(size(Te.flight.alpha), [numel(Te.flight.t) 3]) && ...
                   max(abs(vecnorm(Te.flight.alpha, 2, 2) - 1)) < 1e-12) && ok;
        ok = check('entry in library format (sD sA z8 tf_nd junctions verdicts provenance)', ...
                   all(isfield(Te.entry, {'sD', 'sA', 'z8', 'tf_nd', 'junctions', 'verdicts', 'provenance'})) && ...
                   isequal(Te.entry.z8, Te.z8) && isa(Te.entry.provenance.created, 'datetime') && ...
                   isfield(Te.entry.provenance, 'engineOpts')) && ok;
    end
    wall.elfo = toc(t0);
    fprintf('  (elfo: %.1f s)\n', wall.elfo);
end

%% regate
if ismember('regate', sections)
    t0 = tic;
    T0 = pair_transfer('regate', Te, struct());
    ok = check('regate with no change reproduces pass', T0.gates.pass == Te.gates.pass && T0.gates.pass) && ok;
    R = pair_transfer('regate', Te, struct('floorKm', 1e5));
    ok = check(sprintf('floorKm 1e5: pass %d, moonMinKm %.1f < 1e5, failed {%s}', R.gates.pass, ...
               R.gates.moonMinKm, strjoin(R.gates.failed, ', ')), ~R.gates.pass && ...
               R.gates.moonMinKm < 1e5 && isequal(R.gates.failed, {'moon'})) && ok;
    rec = R.entry.provenance.regate{end};
    ok = check(sprintf('regate recorded in provenance (%d record, floorKm %g = verdict thresholds)', ...
               numel(R.entry.provenance.regate), rec.thresholds.floorKm), ...
               isscalar(R.entry.provenance.regate) && rec.thresholds.floorKm == 1e5 && ...
               isequal(rec.thresholds, R.entry.verdicts.thresholds) && isa(rec.at, 'datetime')) && ok;
    R2 = pair_transfer('regate', R, struct('floorKm', 500));
    ok = check(sprintf('second regate appends (%d records), pass restored %d', ...
               numel(R2.entry.provenance.regate), R2.gates.pass), ...
               numel(R2.entry.provenance.regate) == 2 && R2.gates.pass && ...
               isequal(R2.entry.provenance.regate{end}.thresholds, R2.entry.verdicts.thresholds)) && ok;
    R = pair_transfer('regate', Te, struct('gateKm', 1e-12));
    ok = check(sprintf('gateKm 1e-12: failed {%s}', strjoin(R.gates.failed, ', ')), ...
               ~R.gates.pass && isequal(R.gates.failed, {'flown'})) && ok;
    R = pair_transfer('regate', Te, struct('accTol', 0));
    ok = check(sprintf('accTol 0: failed {%s}', strjoin(R.gates.failed, ', ')), ...
               ~R.gates.pass && isequal(R.gates.failed, {'accDz'})) && ok;
    R = pair_transfer('regate', Te, struct('floorEarthKm', 1e7));
    ok = check(sprintf('floorEarthKm 1e7: failed {%s}', strjoin(R.gates.failed, ', ')), ...
               ~R.gates.pass && isequal(R.gates.failed, {'earth'})) && ok;
    Tt = Te;  Tt.ok = false;
    R = pair_transfer('regate', Tt, struct());
    ok = check(sprintf('tampered ok = false: failed {%s}', strjoin(R.gates.failed, ', ')), ...
               ~R.gates.pass && ismember('ok', R.gates.failed)) && ok;
    Tt = Te;  Tt.gates.admissible = false;
    R = pair_transfer('regate', Tt, struct());
    ok = check(sprintf('tampered admissible = false: failed {%s}', strjoin(R.gates.failed, ', ')), ...
               ~R.gates.pass && isequal(R.gates.failed, {'admissible'})) && ok;
    Tt = Te;  Tt.gates.moonNFailed = 1;
    R = pair_transfer('regate', Tt, struct());
    ok = check(sprintf('tampered Moon minimum not located: failed {%s}', strjoin(R.gates.failed, ', ')), ...
               ~R.gates.pass && isequal(R.gates.failed, {'moon'})) && ok;
    ok = check('regate updates entry.verdicts and thresholds', ...
               ~R.entry.verdicts.pass && R.thresholds.floorKm == Te.thresholds.floorKm) && ok;
    ok = check_error('regate unknown threshold refuses', 'pair_transfer:regate', ...
                     @() pair_transfer('regate', Te, struct('maxIter', 1))) && ok;
    wall.regate = toc(t0);
    fprintf('  (regate: %.1f s)\n', wall.regate);
end

%% cap
if ismember('cap', sections)
    t0 = tic;
    p0 = gcp('nocreate');
    Sc = pair_transfer('second', Te, struct('msHardCapSec', 0.5, 'secondDelaySec', 5));
    ok = check(sprintf('cap 0.5 s + 5 s delay: status %s, conjPass %g, ''%s''', Sc.gates.secondStatus, ...
               Sc.gates.conjPass, Sc.gates.conjReason), strcmp(Sc.gates.secondStatus, 'capped') && ...
               isnan(Sc.gates.conjPass) && startsWith(Sc.gates.conjReason, 'capped after') && ...
               strcmp(Sc.entry.verdicts.secondStatus, 'capped')) && ok;
    ok = check(sprintf('capped run: gates.pass unchanged (%d, second order is not blocking)', Sc.gates.pass), ...
               Sc.gates.pass == Te.gates.pass) && ok;
    ok = check('pool state after the capped run equals the state before', same_pool(p0)) && ok;
    Sf = pair_transfer('second', Te, struct('msHardCapSec', 300));
    ok = check(sprintf('fenced (cap 300 s): status %s, conjPass %g = stored %g', Sf.gates.secondStatus, ...
               Sf.gates.conjPass, Te.gates.conjPass), strcmp(Sf.gates.secondStatus, Te.gates.secondStatus) && ...
               isequal(Sf.gates.conjPass, Te.gates.conjPass) && Sf.gates.dimS == Te.gates.dimS) && ok;
    ok = check('pool state after the fenced run equals the state before', same_pool(p0)) && ok;
    Si = pair_transfer('second', Te, struct('msHardCapSec', 0));
    ok = check(sprintf('in-process (cap 0): status %s, conjPass %g', Si.gates.secondStatus, Si.gates.conjPass), ...
               strcmp(Si.gates.secondStatus, Te.gates.secondStatus) && ...
               isequal(Si.gates.conjPass, Te.gates.conjPass)) && ok;
    wall.cap = toc(t0);
    fprintf('  (cap: %.1f s)\n', wall.cap);
end

%% equiv
if ismember('equiv', sections)
    t0 = tic;
    gto0 = struct('family', 'gto', 'params', struct('orientDeg', 0));
    Tq = cached_solve(C, gto0, dro1, 0.10, 0, 15, ...
                      struct('floorEarthKm', 200, 'directGateKm', 100, 'secondOrder', false));
    eHash = source_hash({'gto_dro_entry', 'thrust_ladder_library', 'ms_tfmin'}, ...
                        fullfile(ot, 'GTO_DRO', 'lib'));
    eFile = fullfile(C.dir, sprintf('gto_dro_entry_demo_%s_%s.mat', cacheVersion, eHash));
    pathAtEquiv = path;
    fromCache = C.use && isfile(eFile);
    if fromCache
        S = load(eFile, 'entry');  entry = S.entry;
        fprintf('  (cached: %s)\n', eFile);
    else
        pBefore = path;
        restoreEq = onCleanup(@() path(pBefore));        % gto_dro_entry runs its setup
        addpath(fullfile(ot, 'GTO_DRO', 'lib'));
        tmpEq = tempname;  mkdir(tmpEq);
        entry = gto_dro_entry(0, 1, 0.10, 0, 15, struct('outMat', fullfile(tmpEq, 'demo.mat')));
        rmdir(tmpEq, 's');
        clear restoreEq                                  % path back now
        save(eFile, 'entry');
    end
    dz = max(abs(Tq.z8 - entry.z8));
    fprintf('  equiv: pair_transfer t_f %.13f (ok %d), gto_dro_entry t_f %.13f (ok %d), max|dz8| %.2e\n', ...
            Tq.tf_nd, Tq.ok, entry.tf_nd, entry.ok, dz);
    ok = check(sprintf('GTO demo cell z8 equals gto_dro_entry to 1e-6 (%.2e)', dz), ...
               Tq.ok && entry.ok && dz < 1e-6) && ok;
    if fromCache
        fprintf('  [SKIP] path after the gto_dro_entry call equals the path before it (result from cache: not called)\n');
    else
        % the shared session may already carry campaign folders: what must
        % hold is that the equivalence call leaves the path as it found it
        ok = check('path after the gto_dro_entry call equals the path before it', ...
                   isequal(path, pathAtEquiv)) && ok;
    end
    wall.equiv = toc(t0);
    fprintf('  (equiv: %.1f s)\n', wall.equiv);
end

%% tulip
if ismember('tulip', sections)
    t0 = tic;
    tul = struct('family', 'tulip', 'params', struct('Np', 7, 'pm', -1));
    catFile = fullfile(ot, 'DRO_tulip', 'direct', 'results', 'costate_catalog_dro_tulip.mat');
    sDt = 0;  sAt = 0.075378;  tfCat = NaN;
    if isfile(catFile)
        S  = load(catFile, 'costate_catalog_dro_tulip');
        sh = S.costate_catalog_dro_tulip.sheets;
        sh = sh([sh.tauDRO] == 1 & [sh.Np] == 7 & [sh.pm] == -1);
        [iDs, iAs] = find(sh.has_solution(:, :, 1));
        dW = @(a, b) abs(mod(a - b + 0.5, 1) - 0.5);
        [~, kBest] = min(dW(sh.sD_frac(iDs), 0).^2 + dW(sh.sA_frac(iAs), 0.075378).^2);
        iD = iDs(kBest);  iA = iAs(kBest);
        sDt = sh.sD_frac(iD);  sAt = sh.sA_frac(iA);  tfCat = sh.tf_nd(iD, iA, 1);
        fprintf('  tulip: catalog cell (%d, %d) of %d solved 15 N cells: sD %.6f, sA %.6f\n', ...
                iD, iA, numel(iDs), sDt, sAt);
    else
        fprintf('  tulip: shipped catalog not found (%s); solving at (0, 0.075378), comparison skipped\n', catFile);
    end
    Tu = cached_solve(C, dro1, tul, sDt, sAt, 15, struct());
    g = Tu.gates;
    fprintf(['  tulip 15 N: ok %d, t_f %.10f ND (%.4f d), flown %.2e km, Moon located %.1f km, ' ...
             'conjPass %g, second %s\n'], Tu.ok, Tu.tf_nd, Tu.tf_d, g.flownKm, g.moonMinKm, ...
            g.conjPass, g.secondStatus);
    ok = check(sprintf('tulip gates.pass (failed: %s)', strjoin(g.failed, ', ')), g.pass) && ok;
    if isfile(catFile)
        fprintf('  tulip: pair_transfer t_f %.10f ND beside the shipped catalog %.10f ND (diff %+.3e, %+.3f%%)\n', ...
                Tu.tf_nd, tfCat, Tu.tf_nd - tfCat, 100*(Tu.tf_nd - tfCat)/tfCat);
    end
    wall.tulip = toc(t0);
    fprintf('  (tulip: %.1f s)\n', wall.tulip);
end

if ok, fprintf('test_pair_transfer: ALL PASS\n');
else,  fprintf('test_pair_transfer: FAILURES\n'); end
end

% ------------------------------------------------------------------------
function T = cached_solve(C, dep, arr, sD, sA, thrustN, opts)
% CACHED_SOLVE  pair_transfer through a tempdir cache keyed by a SHA-256 of
%   every input and of the source hash in C.ver; C.use false re-solves
%   (and refreshes the cache).  INPUTS: C (.dir .ver .use); dep, arr;
%   sD, sA; thrustN; opts.  OUTPUTS: T.
key = jsonencode(struct('dep', dep, 'arr', arr, 'sD', sD, 'sA', sA, 'thrustN', thrustN, ...
                        'opts', opts, 'ver', C.ver));
f   = fullfile(C.dir, ['pt_' sha_hex(key) '.mat']);
if C.use && isfile(f)
    S = load(f, 'T');  T = S.T;
    fprintf('  (cached: %s)\n', f);
else
    opts.quiet = true;
    T = pair_transfer(dep, arr, sD, sA, thrustN, opts);
    save(f, 'T');
end
end

% ------------------------------------------------------------------------
function h = source_hash(names, extraDir)
% SOURCE_HASH  SHA-256 (first 16 hex) over the source text of the named
%   functions, resolved with which (extraDir searched first when given).
%   INPUTS: names cellstr; extraDir char (optional).  OUTPUTS: h char.
parts = cell(1, numel(names));
for k = 1:numel(names)
    f = '';
    if nargin > 1 && isfile(fullfile(extraDir, [names{k} '.m']))
        f = fullfile(extraDir, [names{k} '.m']);
    end
    if isempty(f), f = which(names{k}); end
    if isempty(f)
        error('test_pair_transfer:source', 'cannot resolve %s for the cache key', names{k});
    end
    parts{k} = [names{k}, ':', sha_hex(fileread(f))];
end
h = sha_hex(strjoin(parts, ';'));
h = h(1:16);
end

% ------------------------------------------------------------------------
function s = sha_hex(str)
% SHA_HEX  SHA-256 of a char array (UTF-8) as lower-case hex.
%   INPUTS: str char.  OUTPUTS: s char [1 x 64].
md = java.security.MessageDigest.getInstance('SHA-256');
b  = typecast(md.digest(unicode2native(str, 'UTF-8')), 'uint8');
s  = lower(reshape(dec2hex(b, 2).', 1, []));
end

% ------------------------------------------------------------------------
function tf = same_pool(p0)
% SAME_POOL  The pool state equals p0's: none then and none now, or the
%   same pool still open.  INPUTS: p0 gcp('nocreate') taken earlier.
%   OUTPUTS: tf logical.
p1 = gcp('nocreate');
if isempty(p0)
    tf = isempty(p1);
else
    tf = ~isempty(p1) && isvalid(p0) && p1 == p0;
end
end

% ------------------------------------------------------------------------
function s = onoff(c)
% ONOFF  'on' / 'off'.  INPUTS: c logical.  OUTPUTS: s char.
if c, s = 'on'; else, s = 'off'; end
end

% ------------------------------------------------------------------------
function ok = check_error(name, id, fcn)
% CHECK_ERROR  fcn must throw error id.  INPUTS: name; id; fcn handle.
%   OUTPUTS: ok.
try
    fcn();
    ok = check([name ' (no error thrown)'], false);
catch err
    ok = check(sprintf('%s (%s)', name, err.identifier), strcmp(err.identifier, id));
end
end

% ------------------------------------------------------------------------
function ok = check(name, cond)
% CHECK  One gate line.  INPUTS: name; cond.  OUTPUTS: ok.
ok = logical(cond);
if ok, tag = 'PASS'; else, tag = 'FAIL'; end
fprintf('  [%s] %s\n', tag, name);
end
