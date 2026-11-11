%% LUNAR_PNT_MISSION_STUDY  Six spacecraft from a DRO onto the six tulip slots, then the south-pole GDOP time series
%   Edit the knobs, press Run. Six spacecraft wait in the 4.43 d DRO at
%   phases (k-1)/6; each is scheduled onto one of the six slots of the
%   7-petal tulip constellation (the paper's phasing) through a certified
%   70 mN min-time catalog cell and polished onto its EXACT slot state at a
%   pinned arrival epoch. The product is ONE time series: GDOP at the lunar
%   south pole (latitude -90, 0 deg mask) from clock 0, through the
%   assembly, to serviceDays after the last arrival. Every step calls a
%   tested library function (lib/); nothing is reimplemented here.
%   The LunaNet R1-R3 requirement scoring (validated against the paper's
%   Table 1 floors) exists behind scoreRequirements = true and is NOT run by
%   default.
%
%     0  knobs                              (one source; study_opts overrides)
%     1  inputs and the orbit handles       I1 inputs, I2 interpolation vs propagation
%                                           (position and velocity)
%     2  schedule                           C1 a complete assignment exists
%     3  polish: schedule -> polish ->      P1..P6 accepted, PF flown, PR residual,
%        exclude -> re-schedule loop        PP pointwise, PH hypotheses, PC conjugate,
%                                           PA pumpkyn tfMin acceptance
%     4  histories                          H1 departure jumps, H2 arrival jumps
%     5  validation                         G0 pole series vs tulipConstellation
%                                           (G0R Table 1 when scoreRequirements);
%                                           anything but PASS, incl. NOT RUN, stops
%     6  the south-pole GDOP series         summary statistics (no verdicts);
%                                           R1, R2, R3 only when scoreRequirements
%     7  figures
%     8  movie                              mission_movie on the SAVED results,
%                                           rendered after the save in 9
%     9  gate table, VERDICT, save
%
%   HOW THE SCHEDULE IS FOUND. A polish can refuse its candidate (a root
%   that needs a departure before clock 0, a root not found, a failed or
%   erroring certification gate). A refusal changes the bottleneck, so the
%   assignment is re-optimised rather than walked per spacecraft: every
%   refused candidate [k j iD iA nRev] joins schedule_rendezvous's exclusion
%   list, CUMULATIVELY (it stays excluded in every later round), the
%   schedule is recomputed, and only spacecraft whose chosen candidate
%   changed are polished (polishes are cached by candidate key). The loop
%   ends when all six chosen candidates are accepted (P gates PASS), after
%   maxRounds, or when no complete assignment remains (C1 FAIL).
%
%   study_opts: a struct in the workspace before running overrides any knob
%   (e.g. study_opts.serviceDays = 3; study_opts.headless = true). Unknown
%   names are refused. The script CLEARS its workspace (clearvars -except
%   study_opts): run it from a function, or accept that the base workspace
%   is cleared. Used by tests/test_lunar_pnt_mission_study.
%
%   REUSE (knob reuse): polishes from an earlier results file seed the
%   cache only if (1) the file's fingerprint (catalog name, created date,
%   MD5 of sh.z8 and sh.tf_nd, regime constants, dtau6, the spacecraft DRO
%   phases droPhase0, MD5 of the DRO and tulip handles at 64 phases,
%   polish thresholds, polishGates) equals this run's -- otherwise the
%   whole reuse is refused (lunar_pnt_mission_study:reuse; a file written
%   before a fingerprint field existed is refused too) -- and (2) each
%   entry, when taken for a candidate, has the candidate's tArr (1e-12),
%   tDep + tf = tArr (1e-12), the current slot state as its rvf (1e-12),
%   the spacecraft's current DRO state at its tDep as its rv0 (1e-12) and
%   the current thresholds; an entry failing any of these is ignored and
%   re-polished (printed). Reused polishes are NOT re-certified: the P rows
%   are tagged REUSED, run_meta.polishSource names their source, and the
%   VERDICT says so; reused REFUSALS (which decide the exclusions) are
%   named in the VERDICT too, even when all six final polishes are fresh.
%
%   PROVENANCE: run_meta records the HEAD of this repository (found from
%   the script's folder) and of proj7, each with a dirty flag (tracked
%   files modified); a fresh polish's "certified at" head carries -dirty
%   when the working tree was dirty. run_meta.verdictSuffix keeps the
%   VERDICT's disclosure text.
%
%  M. Casey                                                   (c) 11/10/2026
%  M. Casey  rescope: south-pole GDOP series by default        (c) 11/10/2026
%  M. Casey  guarded reuse (fingerprint, per-entry checks)      (c) 11/10/2026
%  M. Casey  final-review fixes: git dirty flags, I2 velocity, orbit
%            fingerprint + rv0 entry check, pole_series_stats, dtSec
%            sets the pole cadence, H.r saved as double       (c) 11/10/2026
%  Copyright Coorbital Inc.

clearvars -except study_opts
if ~exist('study_opts', 'var'), study_opts = struct(); end
here = fileparts(mfilename('fullpath'));

%% ========================================================================
%  0. KNOBS -- every adjustable parameter, with units; one source.
%% ========================================================================
tol         = 0.005;                 % arrival phase tolerance [fraction of a tulip period]
nMaxRevs    = 10;                    % longest DRO wait [whole revolutions]
dtSec       = 300;                   % history grid step AND pole-series cadence [s]
serviceDays = 30;                    % service window after the last arrival [d]
validate    = true;                  % run the validation gate(s); false = STOPPED AT VALIDATION
maxRounds   = 12;                    % schedule -> polish -> exclude rounds before giving up
polishGates = true;                  % full certification on every polish (false: fly gate only)
i2Km        = 1;                     % I2 limit: orbit handle vs direct propagation, position [km]
i2Mps       = 0.1;                   % I2 limit: the same, velocity [m/s] (measured max 0.0016)
i2Phases    = 2400;                  % I2 uniform phases (plus 200 around the tulip perilune)
scoreRequirements = false;           % also score R1-R3 on the paper grid (validated vs Table 1)
scoreGrid   = 'paper';               % grid R1-R3 are scored on (scoreRequirements only)
chunkSel    = [];                    % score-grid chunks ([] = all five; scoreRequirements only)
reuse       = '';                    % an earlier lunar_pnt_mission_study.mat whose polishes seed the cache
makeFigures = true;                  % section 7
makeMovie   = false;                 % section 8: render results/lunar_pnt_mission.mp4 (~5 min)
headless    = false;                 % figures invisible and closed after export
outDir      = '';                    % results folder ('' = results/ beside this script)
g0SlotOverride = [];                 % TEST SEAM: [1 x 6] slot offsets [ND] for G0's own series (deliberate failure); [] = none

% study_opts overrides; unknown names are refused BEFORE anything is assigned
knobs = struct('tol', tol, 'nMaxRevs', nMaxRevs, 'dtSec', dtSec, 'serviceDays', serviceDays, ...
               'validate', validate, 'maxRounds', maxRounds, 'polishGates', polishGates, ...
               'i2Km', i2Km, 'i2Mps', i2Mps, 'i2Phases', i2Phases, 'scoreRequirements', scoreRequirements, ...
               'scoreGrid', scoreGrid, 'chunkSel', chunkSel, 'reuse', reuse, ...
               'makeFigures', makeFigures, 'makeMovie', makeMovie, 'headless', headless, ...
               'outDir', outDir, 'g0SlotOverride', g0SlotOverride);
knobs0 = knobs;
for f = fieldnames(study_opts).'
    if ~isfield(knobs, f{1})
        error('lunar_pnt_mission_study:knob', 'unknown study_opts field ''%s''', f{1});
    end
end
for f = fieldnames(study_opts).'
    knobs.(f{1}) = study_opts.(f{1});
end
if isempty(knobs.outDir), knobs.outDir = fullfile(here, 'results'); end
tol         = knobs.tol;          nMaxRevs    = knobs.nMaxRevs;
dtSec       = knobs.dtSec;        serviceDays = knobs.serviceDays;
validate    = knobs.validate;     maxRounds   = knobs.maxRounds;
polishGates = knobs.polishGates;  i2Km        = knobs.i2Km;
i2Mps       = knobs.i2Mps;
i2Phases    = knobs.i2Phases;     scoreRequirements = knobs.scoreRequirements;
scoreGrid   = knobs.scoreGrid;    chunkSel    = knobs.chunkSel;
reuse       = knobs.reuse;        makeFigures = knobs.makeFigures;
makeMovie   = knobs.makeMovie;    headless    = knobs.headless;
outDir      = knobs.outDir;       g0SlotOverride = knobs.g0SlotOverride;
if ~isfolder(outDir), mkdir(outDir); end
figVis = 'on';
if headless, figVis = 'off'; end
% knobs that change the result (display knobs excluded) and differ from the defaults
knNames = setdiff(fieldnames(knobs0), {'headless', 'outDir', 'makeFigures', 'makeMovie', 'reuse'}, 'stable');
isRed   = cellfun(@(n) ~isequal(knobs.(n), knobs0.(n)), knNames);
reduced = cellfun(@(n) sprintf('%s %s', n, knob_text(knobs.(n))), knNames(isRed), 'UniformOutput', false);

gates = cell(0, 4);                                    % {id, status, value, detail}
S = [];  S0 = [];  P = [];  Pattempted = {};  H = [];  V = [];  Vpole = [];
Gservice = [];  pole = [];
polishCache = struct('key', {{}}, 'P', {{}}, 'origin', {{}});
polishSource = struct('sc', num2cell(1:6), 'origin', 'none', 'file', '', 'date', '', 'head', '');
rounds = struct('round', {}, 'assign', {}, 'tLastDays', {}, 'polished', {}, ...
                'refused', {}, 'refusedRows', {}, 'wallS', {});
I2 = struct();  stopped = '';
secName = {'1 inputs', '2 schedule', '3 polish', '4 histories', '5 validation', ...
           '6 GDOP', '7 figures', '8 movie'};
secWall = nan(1, numel(secName));
tStudy  = tic;
reusedRefusals = cell(0, 2);                           % {key, source file} of reused REFUSED polishes
nRefLog = 0;                                           % rows used in refLog (section 3)

% provenance (read-only git through the Command Line Tools binary): the
% repository holding this script, located from the script's own folder,
% and proj7, located as setup_lunar_pnt_paths does; each with its HEAD and
% a dirty flag (tracked files modified, untracked files ignored)
gitBin = '/Library/Developer/CommandLineTools/usr/bin/git';
gitOC  = git_state(gitBin, here);
gitP7  = git_state(gitBin, fullfile(getenv('HOME'), 'Desktop', 'proj7'));
fprintf('lunar_pnt_mission_study: optimal_control %s; proj7 %s; MATLAB %s\n', ...
        gitOC.text, gitP7.text, version);

%% ========================================================================
%  1. INPUTS AND THE ORBIT HANDLES
%% ========================================================================
tSec = tic;
fprintf('1. INPUTS\n');
run(fullfile(here, 'setup_lunar_pnt_paths.m'));
M = mission_inputs();
% the identity a reused results file must match (thresholds read from the
% polish itself: a preflight refusal returns them without any solve)
thrNow = polish_rendezvous(M, 1, 1, struct('iD', 1, 'iA', 1, 'nRev', 0, 'tDep', 0, 'tArr', 0), ...
                           struct('gates', polishGates)).thresholds;
fingerprint = catalog_fingerprint(M, thrNow, polishGates);
dayNd = M.tStar/86400;
fprintf('   70 mN catalog %dx%d cells; DRO period %.6f d, tulip period %.6f d; propellant %.4f kg/ND\n', ...
        numel(M.sD), numel(M.sA), M.tauD*dayNd, M.tauT*dayNd, M.propKgPerNd);
nCk = numel(M.checks);
gates = addg(gates, 'I1', pass3(all([M.checks.pass])), nCk, ...
             sprintf('%d of %d mission_inputs checks', nnz([M.checks.pass]), nCk));

% I2: the periodic phase_state handles vs an independent tight propagation
% (ode113 on pumpkyn's CR3BP field, RelTol = AbsTol = 1e-13) from the
% handle's own phase-0 state, at the same orbit times
odeTight = odeset('RelTol', 1e-13, 'AbsTol', 1e-13);
sI2 = unique([linspace(0, 1, i2Phases + 1), linspace(0.45, 0.55, 201)]);
sI2 = sI2(sI2 < 1);
hand = {'tulip', M.arrState, M.tauT; 'DRO', M.depState, M.tauD};
for kh = 1:2
    [~, Yp] = ode113(@pumpkyn.cr3bp.eom, sI2*hand{kh, 3}, hand{kh, 2}(0).', odeTight, M.muStar, 1);
    Xh  = hand{kh, 2}(sI2).';
    eKm = M.lStar*sqrt(sum((Xh(:, 1:3) - Yp(:, 1:3)).^2, 2));
    eMs = 1000*M.lStar/M.tStar*sqrt(sum((Xh(:, 4:6) - Yp(:, 4:6)).^2, 2));
    [eMax, iMax] = max(eKm);
    [vMax, iVmx] = max(eMs);
    rMoon = M.lStar*sqrt(sum((Xh(iMax, 1:3) - [1 - M.muStar 0 0]).^2, 2));
    I2.(lower(hand{kh, 1})) = struct('maxKm', eMax, 'atPhase', sI2(iMax), 'rMoonKm', rMoon, ...
                                     'medianKm', median(eKm), 'p99Km', prctile(eKm, 99), 'nPhases', numel(sI2), ...
                                     'maxMps', vMax, 'atPhaseV', sI2(iVmx), 'medianMps', median(eMs));
    fprintf('   I2 %-5s handle vs propagation (%d phases): max %.4f km at phase %.4f (%.0f km from the Moon), median %.4f km, p99 %.4f km; velocity max %.2e m/s at phase %.4f\n', ...
            hand{kh, 1}, numel(sI2), eMax, sI2(iMax), rMoon, median(eKm), prctile(eKm, 99), vMax, sI2(iVmx));
end
% for information: the slot states vs pumpkyn's tulipConstellation (no
% continuation; the GDOP reference's own construction)
[tauC, rC] = pumpkyn.cr3bp.tulipConstellation(M.sh.Np, M.tauT, 1, M.sh.pm, M.dtau6);
dC = zeros(numel(tauC), 6);
for jj = 1:6
    Xs = M.arrState(mod(M.slotPhase0(jj) + tauC(:).'/M.tauT, 1)).';
    dC(:, jj) = M.lStar*sqrt(sum((Xs(:, 1:3) - rC(:, :, jj)).^2, 2));
end
[dMax, iC] = max(dC(:));
[iCt, iCj] = ind2sub(size(dC), iC);
phC = mod(M.slotPhase0(iCj) + tauC(iCt)/M.tauT, 1);
I2.vsTulipConstellation = struct('medianKm', median(dC(:)), 'p99Km', prctile(dC(:), 99), ...
                                 'maxKm', dMax, 'atPhase', phC);
fprintf('   (info) slot states vs tulipConstellation at its own samples: median %.3f km, p99 %.3f km, max %.3f km at orbit phase %.4f\n', ...
        median(dC(:)), prctile(dC(:), 99), dMax, phC);
if I2.tulip.maxKm < 0.1*dMax
    fprintf('          the tight propagation agrees with the phase_state handles, not with tulipConstellation (no continuation)\n');
end
i2Max  = max(I2.tulip.maxKm, I2.dro.maxKm);
i2VMax = max(I2.tulip.maxMps, I2.dro.maxMps);
i2Ok   = i2Max < i2Km && i2VMax < i2Mps;
gates = addg(gates, 'I2', pass3(i2Ok), i2Max, ...
             sprintf(['tulip max %.4f km @ phase %.4f, DRO max %.4f km @ phase %.4f (limit %g km); ' ...
                      'velocity max %.2e m/s (limit %g m/s); value = max km'], ...
                     I2.tulip.maxKm, I2.tulip.atPhase, I2.dro.maxKm, I2.dro.atPhase, i2Km, i2VMax, i2Mps));
if ~i2Ok
    stopped = sprintf('STOPPED AT INPUTS (I2: orbit handle error %.3g km, %.3g m/s)', i2Max, i2VMax);
end
secWall(1) = toc(tSec);

%% ========================================================================
%  2. SCHEDULE (catalog cells + whole-revolution waits; minimise the last arrival)
%% ========================================================================
tSec = tic;
if isempty(stopped)
    fprintf('\n2. SCHEDULE: tol %.4g of a tulip period, waits 0..%d DRO revolutions\n', tol, nMaxRevs);
    sOpts = struct('tol', tol, 'nMaxRevs', nMaxRevs, 'propKgPerNd', M.propKgPerNd);
    c1Detail = '';
    try
        S0 = schedule_rendezvous(M.tf, M.sD, M.sA, M.tauD, M.tauT, M.droPhase0, M.slotPhase0, sOpts);
    catch err
        if ~strcmp(err.identifier, 'schedule_rendezvous:noAssignment'), rethrow(err); end
        c1Detail = err.message;
    end
    if isempty(S0)
        fprintf('   NO COMPLETE ASSIGNMENT: %s\n', c1Detail);
        gates = addg(gates, 'C1', 'FAIL', NaN, 'no complete assignment at these knobs');
        stopped = 'STOPPED AT SCHEDULE (C1: no complete assignment)';
    else
        fprintf('   best arrival [d], rows spacecraft, columns slot (admissible candidates in brackets):\n');
        for kk = 1:6
            fprintf('   sc %d', kk);
            fprintf('  %8.3f (%3d)', [S0.best.tArr(kk, :)*dayNd; S0.nAdmissible(kk, :)]);
            fprintf('\n');
        end
        fprintf('   assignment (slot of spacecraft 1..6): [%s]; last arrival %.3f d; propellant %.3f kg (unpolished)\n', ...
                sprintf('%d ', S0.assign), S0.tLast*dayNd, S0.totalPropKg);
        fprintf('   sensitivity of the last arrival [d] (unpolished):\n   %10s', 'tol \ revs');
        sensRevs = [5 10 20];  sensTol = [0.0025 0.005 0.01];
        fprintf(' %8d', sensRevs);  fprintf('\n');
        sens = nan(numel(sensTol), numel(sensRevs));
        for a = 1:numel(sensTol)
            fprintf('   %10.4g', sensTol(a));
            for b = 1:numel(sensRevs)
                try
                    Ss = schedule_rendezvous(M.tf, M.sD, M.sA, M.tauD, M.tauT, M.droPhase0, M.slotPhase0, ...
                             struct('tol', sensTol(a), 'nMaxRevs', sensRevs(b), 'propKgPerNd', M.propKgPerNd));
                    sens(a, b) = Ss.tLast*dayNd;
                    fprintf(' %8.3f', sens(a, b));
                catch err
                    if ~strcmp(err.identifier, 'schedule_rendezvous:noAssignment'), rethrow(err); end
                    fprintf(' %8s', 'none');
                end
            end
            fprintf('\n');
        end
        gates = addg(gates, 'C1', 'PASS', S0.tLast*dayNd, 'complete assignment; value = unpolished last arrival [d]');
    end
end
secWall(2) = toc(tSec);

%% ========================================================================
%  3. POLISH -- schedule -> polish -> exclude refused (cumulative) -> re-schedule
%% ========================================================================
tSec = tic;
if isempty(stopped)
    fprintf('\n3. POLISH onto the exact slots (gates %s), at most %d rounds\n', pass3(polishGates), maxRounds);
    cache     = containers.Map('KeyType', 'char', 'ValueType', 'any');
    cacheOrig = containers.Map('KeyType', 'char', 'ValueType', 'any');   % key -> origin struct
    unchecked = containers.Map('KeyType', 'char', 'ValueType', 'logical'); % reused, not yet entry-checked
    freshOrig = struct('file', '', 'date', '', 'head', '');               % filled at save time
    refLog    = cell(6*maxRounds, 2);                                     % reused refusals {key, source}
    if ~isempty(reuse)
        L = load_reuse(reuse, fingerprint);                 % refuses on any identity mismatch
        for kc = 1:numel(L.polishCache.key)
            cache(L.polishCache.key{kc})     = L.polishCache.P{kc};
            cacheOrig(L.polishCache.key{kc}) = L.polishCache.origin{kc};
            unchecked(L.polishCache.key{kc}) = true;
        end
        fprintf('   reuse: fingerprint matches; %d cached polish(es) from %s (each checked when taken)\n', ...
                cache.Count, reuse);
    end
    exclude = nan(6*maxRounds, 5);  nEx = 0;              % rows [k j iD iA nRev], cumulative
    S = S0;  done = false;  Pc = cell(1, 6);
    roundsC = cell(1, maxRounds);
    for rr = 1:maxRounds
        tRound = tic;
        if rr > 1
            try
                S = schedule_rendezvous(M.tf, M.sD, M.sA, M.tauD, M.tauT, M.droPhase0, M.slotPhase0, ...
                                        setfield_(sOpts, 'exclude', exclude(1:nEx, :)));
            catch err
                if ~strcmp(err.identifier, 'schedule_rendezvous:noAssignment'), rethrow(err); end
                fprintf('   round %d: NO COMPLETE ASSIGNMENT after %d exclusions\n', rr, nEx);
                gates = setg(gates, 'C1', 'FAIL', NaN, sprintf('no complete assignment after %d exclusions (round %d)', ...
                                                                nEx, rr));
                stopped = 'STOPPED AT POLISH (C1: no complete assignment remains after exclusions)';
                break
            end
        end
        fprintf('   round %d: assignment [%s], last arrival %.3f d\n', rr, sprintf('%d ', S.assign), S.tLast*dayNd);
        isNew = false(1, 6);  refused = cell(1, 6);  refRows = nan(6, 5);
        for kk = 1:6
            jj  = S.assign(kk);
            key = sprintf('%d_%d_%d_%d_%d', kk, jj, S.iD(kk), S.iA(kk), S.nRev(kk));
            if isKey(unchecked, key)
                rvfNow = M.arrState(mod(M.slotPhase0(jj) + S.tArr(kk)/M.tauT, 1));
                rv0Now = M.depState(mod(M.droPhase0(kk) + cache(key).tDep/M.tauD, 1));
                why = reuse_entry_problem(cache(key), S.tArr(kk), rvfNow, rv0Now, thrNow);
                remove(unchecked, key);
                if isempty(why)
                    o = cacheOrig(key);
                    fprintf('      reused sc %d -> slot %d, cell (%d,%d), wait %d: %s (certified %s at %s)\n', ...
                            kk, jj, S.iD(kk), S.iA(kk), S.nRev(kk), cache(key).reason, o.date, o.head);
                    if ~cache(key).ok
                        nRefLog = nRefLog + 1;
                        refLog(nRefLog, :) = {key, sprintf('%s, certified %s at %s', o.file, o.date, o.head)};
                    end
                else
                    fprintf('      IGNORED cached sc %d -> slot %d, cell (%d,%d), wait %d: %s; re-polishing\n', ...
                            kk, jj, S.iD(kk), S.iA(kk), S.nRev(kk), why);
                    remove(cache, key);  remove(cacheOrig, key);
                end
            end
            if ~isKey(cache, key)
                cand = struct('iD', S.iD(kk), 'iA', S.iA(kk), 'nRev', S.nRev(kk), ...
                              'tDep', S.tDep(kk), 'tArr', S.tArr(kk));
                Pk = polish_rendezvous(M, kk, jj, cand, struct('gates', polishGates));
                cache(key) = Pk;
                cacheOrig(key) = freshOrig;
                isNew(kk) = true;
                fprintf('      polish sc %d -> slot %d, cell (%d,%d), wait %d: %s (%.1f s)\n', ...
                        kk, jj, S.iD(kk), S.iA(kk), S.nRev(kk), Pk.reason, Pk.wallS);
                if startsWith(Pk.reason, 'gate error')
                    fprintf('      *** GATE ERROR on sc %d: %s ***\n', kk, Pk.reason);
                end
            end
            Pc{kk} = cache(key);
            o = cacheOrig(key);
            if isempty(o.file)
                polishSource(kk) = struct('sc', kk, 'origin', 'fresh', 'file', '', 'date', '', 'head', '');
            else
                polishSource(kk) = struct('sc', kk, 'origin', 'reused', 'file', o.file, 'date', o.date, 'head', o.head);
            end
            if ~Pc{kk}.ok
                refused{kk} = sprintf('sc %d cell (%d,%d) n %d: %s', kk, S.iD(kk), S.iA(kk), S.nRev(kk), Pc{kk}.reason);
                refRows(kk, :) = [kk jj S.iD(kk) S.iA(kk) S.nRev(kk)];
            end
        end
        refused = refused(~cellfun(@isempty, refused));
        refRows = refRows(~isnan(refRows(:, 1)), :);
        roundsC{rr} = struct('round', rr, 'assign', S.assign, 'tLastDays', S.tLast*dayNd, ...
                             'polished', find(isNew), 'refused', {refused}, 'refusedRows', refRows, ...
                             'wallS', toc(tRound));
        if isempty(refRows)
            fprintf('   round %d: all six accepted (%d new polish(es), %.0f s)\n', rr, nnz(isNew), toc(tRound));
            done = true;
            break
        end
        fprintf('   round %d: %d refused -> excluded (cumulative %d): %s\n', rr, size(refRows, 1), ...
                nEx + size(refRows, 1), strjoin(refused, '; '));
        exclude(nEx + (1:size(refRows, 1)), :) = refRows;
        nEx = nEx + size(refRows, 1);
    end
    rounds = [roundsC{:}];
    reusedRefusals = refLog(1:nRefLog, :);
    ck = keys(cache);
    polishCache = struct('key', {ck}, 'P', {values(cache, ck)}, 'origin', {values(cacheOrig, ck)});
    Pattempted = Pc;
    if ~done && isempty(stopped)
        stopped = sprintf('STOPPED AT POLISH (not all six accepted after %d round(s))', maxRounds);
    end
    % the table and gates of the final (or last attempted) chosen candidates
    if all(~cellfun(@isempty, Pc))
        if done, label = 'realised'; else, label = 'LAST ATTEMPTED ROUND (not realised)'; end
        rank0 = nan(1, 6);
        for kk = 1:6
            R0 = S0.ranked{kk, S.assign(kk)};
            hit = find(R0(:, 2) == S.iD(kk) & R0(:, 3) == S.iA(kk) & R0(:, 4) == S.nRev(kk), 1);
            if ~isempty(hit), rank0(kk) = hit; end
        end
        fprintf('\n   polish table, %s:\n', label);
        fprintf('   %-3s %-4s %-8s %-4s %8s %8s %10s %10s %5s %9s %8s %8s %6s %5s %s\n', 'sc', 'slot', 'cell', 'wait', ...
                't_dep d', 't_arr d', 'tf cell d', 'tf pol d', 'evals', 'flown km', 'prop kg', 'dV km/s', 'famChg', 'rank', 'status');
        for kk = 1:6
            Pk = Pc{kk};
            fprintf('   %-3d %-4d %-8s %-4d %8.3f %8.3f %10.5f %10.5f %5d %9.2e %8.3f %8.4f %6d %5d %s\n', kk, S.assign(kk), ...
                    sprintf('(%d,%d)', S.iD(kk), S.iA(kk)), S.nRev(kk), Pk.tDep*dayNd, Pk.tArr*dayNd, ...
                    Pk.tfCell*dayNd, Pk.tf*dayNd, Pk.nEval, gfield(Pk.gates, 'flyKm'), Pk.propKg, Pk.dvKms, ...
                    Pk.familyChanged, rank0(kk), [tern(Pk.ok, 'accepted', 'REFUSED') ...
                    tern(strcmp(polishSource(kk).origin, 'reused'), ' (REUSED)', '')]);
        end
        fprintf('   rank = position of the chosen candidate in the round-1 ranked list of its (spacecraft, slot) pair\n');
        Pall = [Pc{:}];
        if done
            tRealised = max([Pall.tArr]);
            assert(abs(tRealised - S.tLast) <= 1e-12, 'lunar_pnt_mission_study:epoch', ...
                   'polished last arrival %.15g ND differs from the scheduler''s %.15g ND', tRealised, S.tLast);
            fprintf('   last arrival: %.3f d unpolished first round -> %.3f d realised (max of the polished arrival epochs); total propellant %.3f kg\n', ...
                    S0.tLast*dayNd, tRealised*dayNd, sum([Pall.propKg]));
            P = Pall;
        else
            fprintf('   last attempted round: scheduled last arrival %.3f d; NOT all six accepted\n', S.tLast*dayNd);
        end
        Gs = [Pall.gates];
        isReu = strcmp({polishSource.origin}, 'reused');
        tagAll = tern(any(isReu), 'REUSED (not re-certified here): ', '');
        for kk = 1:6
            gates = addg(gates, sprintf('P%d', kk), pass3(Pc{kk}.ok), double(Pc{kk}.ok), ...
                         sprintf('%ssc %d slot %d: %s', tern(isReu(kk), 'REUSED (not re-certified here): ', ''), ...
                                 kk, S.assign(kk), Pc{kk}.reason));
        end
        gates = addg(gates, 'PF', pass3(all(arrayfun(@(g) isequal(g.passFly, true), Gs))), max([Gs.flyKm]), ...
                     sprintf('%sflown miss < 100 km and < 10 m/s on all six; value = max km (max %.2e m/s)', tagAll, max([Gs.flyVms])));
        gates = addg(gates, 'PR', pass3(all(arrayfun(@(g) isequal(g.passResid, true), Gs))), max([Gs.resid]), ...
                     [tagAll 'ms_tfmin reproduces the root, |R| < 1e-8; value = max |R|']);
        gates = addg(gates, 'PP', pass3(all(arrayfun(@(g) isequal(g.passPW, true), Gs))), ...
                     max(arrayfun(@(g) gfield(g.pw, 'Hmax'), Gs)), [tagAll 'pmp_pointwise_checks on all six; value = max Hmax']);
        gates = addg(gates, 'PH', pass3(all(arrayfun(@(g) isequal(g.passHyp, true), Gs))), ...
                     min(arrayfun(@(g) gfield(g.hyp, 'minLamV'), Gs)), [tagAll 'minLamV > 0, minQmt > 0, dimS = 1; value = min minLamV']);
        gates = addg(gates, 'PC', conj_status(Gs), sum([Gs.conjPass] == 1), [tagAll 'conjugate test pass and covered; value = count passing']);
        gates = addg(gates, 'PA', pass3(all(arrayfun(@(g) isequal(g.passAcc, true), Gs))), max([Gs.accDz]), ...
                     [tagAll 'pumpkyn tfMin reproduces the polished z8 to 1e-6; value = max |dz|']);
    end
end
secWall(3) = toc(tSec);

%% ========================================================================
%  4. HISTORIES (one clock, one grid; hand-over jumps measured)
%% ========================================================================
tSec = tic;
if isempty(stopped)
    fprintf('\n4. HISTORIES at %g s, %g service days\n', dtSec, serviceDays);
    H = assemble_histories(M, Pc, S.assign, struct('dtSec', dtSec, 'serviceDays', serviceDays));
    fprintf('   %d epochs; arrivals %.3f .. %.3f d; departure jumps [%s] km; arrival jumps [%s] km\n', ...
            numel(H.t), H.tFirst*dayNd, H.tLast*dayNd, sprintf('%.1e ', H.jumpKm.dep), sprintf('%.1e ', H.jumpKm.arr));
    gates = addg(gates, 'H1', pass3(all(H.jumpKm.dep < 1)), max(H.jumpKm.dep), 'DRO -> transfer position jump < 1 km on all six');
    gates = addg(gates, 'H2', pass3(all(H.jumpKm.arr < 100)), max(H.jumpKm.arr), 'transfer -> slot position jump < 100 km on all six');
end
secWall(4) = toc(tSec);

%% ========================================================================
%  5. VALIDATION -- G0 the pole series vs the paper's constellation
%     construction; G0R the Table 1 floors (scoreRequirements only).
%     Anything but PASS (FAIL, UNRESOLVED, NOT RUN) stops before section 6.
%% ========================================================================
tSec = tic;
if isempty(stopped)
    if ~validate
        fprintf('\n5. VALIDATION: NOT RUN (validate = false)\n');
        gates = addg(gates, 'G0', 'NOT RUN', NaN, 'validate = false: the GDOP pipeline is not validated');
        stopped = 'STOPPED AT VALIDATION (G0 NOT RUN)';
    else
        fprintf('\n5. VALIDATION: on-slots pole series, one tulip period, vs tulipConstellation through the engine\n');
        Vpole = validate_pole_series(M, struct('dtauOverride', g0SlotOverride));
        if ~isempty(g0SlotOverride)
            fprintf('   (TEST SEAM g0SlotOverride: G0''s own series uses moved slot offsets)\n');
        end
        fprintf('   G0 %s | %s (%.1f s)\n', Vpole.status, Vpole.reason, Vpole.wallS);
        gates = addg(gates, 'G0', Vpole.status, Vpole.maxRel, Vpole.reason);
        if ~strcmp(Vpole.status, 'PASS')
            stopped = sprintf('STOPPED AT VALIDATION (G0 %s)', Vpole.status);
        elseif scoreRequirements
            vOpts = struct();
            if ~isempty(chunkSel), vOpts.chunkSel = chunkSel; end
            fprintf('   G0R: on-slots, one tulip period, paper grid (chunks %s) vs Table 1\n', chunk_text(chunkSel));
            V = mission_gdop('validate', M, vOpts);
            fprintf('   G0R %s | ours %.6f %.6f %.6f h | ref %.6f %.6f %.6f h (%.0f s)\n   %s\n', V.status, ...
                    V.ours.r1Floor, V.ours.r2Floor, V.ours.r3Floor, V.ref.r1Floor, V.ref.r2Floor, V.ref.r3Floor, ...
                    gfield(V, 'wallS'), V.reason);
            gates = addg(gates, 'G0R', V.status, V.ours.r3Floor, V.reason);
            if ~strcmp(V.status, 'PASS')
                stopped = sprintf('STOPPED AT VALIDATION (G0R %s)', V.status);
            end
        end
    end
end
secWall(5) = toc(tSec);

%% ========================================================================
%  6. THE SOUTH-POLE GDOP SERIES (one unscored call from clock 0), and the
%     R1-R3 scoring when scoreRequirements
%% ========================================================================
tSec = tic;
if isempty(stopped)
    fprintf('\n6. SOUTH-POLE GDOP SERIES\n');
    tSvcEnd = H.tLast + serviceDays*86400/M.tStar;
    assert(abs(tSvcEnd - H.tEnd) <= 1e-9, 'lunar_pnt_mission_study:window', ...
           'history ends at %.12g ND, the service window at %.12g ND', H.tEnd, tSvcEnd);
    gPole = mission_grid('pole');
    gPole.dtSec = dtSec;                               % the pole series shares the history's cadence
    Gp = mission_gdop(H, M, [0 tSvcEnd], gPole, ...
                      struct('score', false, 'keepSeries', true, 'quiet', true));
    pole = struct('tDays', Gp.t*dayNd, 'gdop', Gp.gdop(:, 1), 'nOn', Gp.nOn, ...
                  'nVis', double(Gp.vis(:, 1)), 'arrDays', [P.tArr]*dayNd, ...
                  'tLastDays', H.tLast*dayNd, 'latLon', Gp.poleLatLon, 'dtSec', Gp.grid.dtSec, ...
                  'wallS', Gp.wallS);
    gd  = pole.gdop;
    svc = pole.tDays >= pole.tLastDays;                % the service window (figures)
    st  = pole_series_stats(pole.tDays, gd, pole.tLastDays, pole.dtSec);   % NaN, Inf, >= 6: outages
    pole.stats = st;
    fprintf('   %d epochs at %g s from clock 0 to %.3f d (point %.1f, %g deg; %.2f s, %d engine calls)\n', ...
            numel(gd), pole.dtSec, tSvcEnd*dayNd, pole.latLon, Gp.wallS, Gp.nCalls);
    fprintf('   first finite GDOP %.3f d; first GDOP < 6 %.3f d; arrivals [%s] d\n', st.firstFiniteDay, ...
            st.firstLt6Day, sprintf('%.3f ', sort(pole.arrDays)));
    fprintf('   service window %.3f .. %.3f d (%d samples):\n', pole.tLastDays, tSvcEnd*dayNd, st.svcSamples);
    fprintf('     fraction GDOP < 6 %.4f, fraction finite %.4f, median %.3f, p95 %.3f (finite samples)\n', ...
            st.fracLt6, st.fracFinite, st.medianGdop, st.p95Gdop);
    fprintf('     longest outage (GDOP >= 6 or NaN) %.2f h; samples with 6 <= GDOP < 7: %d (GDOP < 6: %d)\n', ...
            st.longestOutageH, st.n6to7, st.nLt6);
    if scoreRequirements
        sGrid = mission_grid(scoreGrid);
        gOpts = struct();
        if ~isempty(chunkSel), gOpts.chunkSel = chunkSel; end
        Gservice = mission_gdop(H, M, [H.tLast tSvcEnd], sGrid, gOpts);
        on = V.G;
        fprintf('   R1-R3: service %.3f .. %.3f d (%s grid, chunks %s, %d points, %.0f s)\n', H.tLast*dayNd, tSvcEnd*dayNd, ...
                scoreGrid, chunk_text(chunkSel), numel(Gservice.lat), Gservice.wallS);
        fprintf('   %-24s %14s %14s %10s\n', 'requirement (floor)', sprintf('%g-d service', serviceDays), 'on-slots 1 per', 'threshold');
        fprintf('   %-24s %14.6f %14.6f %10s\n', 'R1 vis >= 2 [frac]', Gservice.r1Floor, on.r1Floor, '0.90');
        fprintf('   %-24s %14.6f %14.6f %10s\n', 'R2 GDOP < 6 [frac]', Gservice.r2Floor, on.r2Floor, '0.40');
        fprintf('   %-24s %14.3f %14.3f %10s\n', 'R3 GDOP < 6 [h]', Gservice.r3Floor, on.r3Floor, '9');
        fprintf('   %-24s %14.3f %14s %10s\n', 'R3 GDOP < 7 [h] (info)', Gservice.r3g7Floor, 'not scored', '9');
        gates = addg(gates, 'R1', pass3(Gservice.pass.r1), Gservice.r1Floor, sprintf('%d points below 0.90', Gservice.r1Nfail));
        gates = addg(gates, 'R2', pass3(Gservice.pass.r2), Gservice.r2Floor, sprintf('%d points below 0.40', Gservice.r2Nfail));
        gates = addg(gates, 'R3', pass3(Gservice.pass.r3), Gservice.r3Floor, sprintf('%d points below 9 h', Gservice.r3Nfail));
    end
end
secWall(6) = toc(tSec);

%% ========================================================================
%  7. FIGURES (slide-size fonts, a legend on every data axes, PNG)
%% ========================================================================
tSec = tic;
if isempty(stopped) && makeFigures
    fprintf('\n7. FIGURES -> %s\n', outDir);
    figs = gobjects(0);
    try
        FS = 14;
        cWait = [0.55 0.65 0.80];  cXfer = [0.90 0.55 0.15];  cOn = [0.20 0.60 0.30];
        cSc = lines(6);
        tEndD = H.tEnd*dayNd;

        % schedule.png
        figs(end+1) = figure('Theme', 'light', 'Color', 'w', 'Position', [60 60 1280 600], 'Visible', figVis);
        ax = axes(figs(end));  hold(ax, 'on');
        for kk = 1:6
            hv = 'off';  if kk == 1, hv = 'on'; end
            plot(ax, [0 P(kk).tDep*dayNd], [kk kk], '-', 'Color', cWait, 'LineWidth', 16, 'DisplayName', 'DRO wait', 'HandleVisibility', hv);
            plot(ax, [P(kk).tDep P(kk).tArr]*dayNd, [kk kk], '-', 'Color', cXfer, 'LineWidth', 16, 'DisplayName', 'transfer', 'HandleVisibility', hv);
            plot(ax, [P(kk).tArr*dayNd tEndD], [kk kk], '-', 'Color', cOn, 'LineWidth', 16, 'DisplayName', 'on station', 'HandleVisibility', hv);
            plot(ax, P(kk).tArr*dayNd, kk, 'kv', 'MarkerFaceColor', 'k', 'MarkerSize', 9, 'DisplayName', 'arrival', 'HandleVisibility', hv);
            text(ax, P(kk).tArr*dayNd + 0.3, kk + 0.33, sprintf('slot %d', S.assign(kk)), 'FontSize', FS - 2);
        end
        xline(ax, H.tLast*dayNd, '--k', 'LineWidth', 1.2, 'DisplayName', sprintf('last arrival %.2f d', H.tLast*dayNd));
        set(ax, 'YTick', 1:6, 'YTickLabel', compose('sc %d', 1:6), 'YDir', 'reverse', 'FontSize', FS, 'Color', 'w');
        xlim(ax, [0 tEndD]);  ylim(ax, [0.4 6.6]);  grid(ax, 'on');
        xlabel(ax, 'mission time [d]');
        legend(ax, 'Location', 'southoutside', 'Orientation', 'horizontal', 'FontSize', FS);
        title(ax, sprintf('Schedule: six spacecraft DRO -> tulip slots, last arrival %.2f d, then %g service days', ...
                          H.tLast*dayNd, serviceDays), 'FontSize', FS + 1);
        exportgraphics(figs(end), fullfile(outDir, 'schedule.png'), 'Resolution', 150);

        % assembly_3d.png
        figs(end+1) = figure('Theme', 'light', 'Color', 'w', 'Position', [60 60 1280 720], 'Visible', figVis);
        ax = axes(figs(end));  hold(ax, 'on');
        sPl = linspace(0, 1, 721);
        rD = M.depState(sPl);  rT = M.arrState(sPl);
        plot3(ax, rD(1, :), rD(2, :), rD(3, :), '-', 'Color', [0.1 0.3 0.8], 'LineWidth', 1.6, 'DisplayName', 'DRO (4.43 d)');
        plot3(ax, rT(1, :), rT(2, :), rT(3, :), '-', 'Color', [0.4 0.4 0.4], 'LineWidth', 1.0, 'DisplayName', 'tulip (23.2 d)');
        for kk = 1:6
            Y = P(kk).flight.Y;
            plot3(ax, Y(:, 1), Y(:, 2), Y(:, 3), '-', 'Color', cSc(kk, :), 'LineWidth', 1.2, ...
                  'DisplayName', sprintf('sc %d -> slot %d', kk, S.assign(kk)));
            plot3(ax, Y(end, 1), Y(end, 2), Y(end, 3), 'o', 'Color', cSc(kk, :), 'MarkerFaceColor', cSc(kk, :), ...
                  'MarkerSize', 7, 'HandleVisibility', 'off');
        end
        [sx, sy, sz] = sphere(32);  rMn = 1737.4/M.lStar;
        surf(ax, 1 - M.muStar + rMn*sx, rMn*sy, rMn*sz, 'FaceColor', [0.6 0.6 0.6], 'EdgeColor', 'none', 'DisplayName', 'Moon');
        axis(ax, 'equal');  grid(ax, 'on');  view(ax, [-35 22]);
        set(ax, 'FontSize', FS, 'Color', 'w');
        xlabel(ax, 'x [ND]');  ylabel(ax, 'y [ND]');  zlabel(ax, 'z [ND]');
        legend(ax, 'Location', 'eastoutside', 'FontSize', FS);
        title(ax, 'Assembly in the Earth-Moon rotating frame: six 70 mN min-time transfers', 'FontSize', FS + 1);
        exportgraphics(figs(end), fullfile(outDir, 'assembly_3d.png'), 'Resolution', 150);

        % pole_gdop.png: the full series from clock 0, and the service window
        figs(end+1) = figure('Theme', 'light', 'Color', 'w', 'Position', [60 60 1280 900], 'Visible', figVis);
        tlg = tiledlayout(figs(end), 2, 1, 'TileSpacing', 'compact', 'Padding', 'compact');
        gClip = pole.gdop;  gClip(gClip > 1e3) = 1e3;   % NaN (fewer than four in view) stays a gap
        for panel = 1:2
            ax = nexttile(tlg);  hold(ax, 'on');
            yyaxis(ax, 'left');
            if panel == 1
                plot(ax, pole.tDays, gClip, '-', 'LineWidth', 1.3, 'DisplayName', 'pole GDOP (clipped at 10^3; gap: < 4 in view)');
                set(ax, 'YScale', 'log');  ylim(ax, [1 1e3]);
                xlim(ax, [0 pole.tDays(end)]);
            else
                plot(ax, pole.tDays(svc), pole.gdop(svc), '-', 'LineWidth', 1.3, 'DisplayName', 'pole GDOP');
                ylim(ax, [0 10]);
                xlim(ax, [pole.tLastDays pole.tDays(end)]);
            end
            yline(ax, 6, '--r', 'LineWidth', 1.4, 'DisplayName', 'GDOP = 6');
            ylabel(ax, 'GDOP');
            yyaxis(ax, 'right');
            if panel == 1
                stairs(ax, pole.tDays, pole.nOn, '-', 'LineWidth', 1.6, 'DisplayName', 'satellites on station');
            end
            stairs(ax, pole.tDays, pole.nVis, ':', 'LineWidth', 1.2, 'DisplayName', 'satellites visible (0 deg mask)');
            ylim(ax, [0 6.5]);  ylabel(ax, 'satellites');
            if panel == 1
                xa = xline(ax, sort(pole.arrDays), ':k', 'LineWidth', 1.0);
                set(xa(2:end), 'HandleVisibility', 'off');  xa(1).DisplayName = 'arrivals';
                title(ax, sprintf('South-pole GDOP from clock 0: first finite %.2f d, first < 6 at %.2f d', ...
                                  st.firstFiniteDay, st.firstLt6Day), 'FontSize', FS);
            else
                title(ax, sprintf('Service window (%g d): GDOP < 6 %.1f%% of the time, longest outage %.2f h', ...
                                  serviceDays, 100*st.fracLt6, st.longestOutageH), 'FontSize', FS);
            end
            grid(ax, 'on');  set(ax, 'FontSize', FS, 'Color', 'w');
            xlabel(ax, 'mission time [d]');
            legend(ax, 'Location', 'northwest', 'FontSize', FS - 2);
        end
        exportgraphics(figs(end), fullfile(outDir, 'pole_gdop.png'), 'Resolution', 150);

        % service_map.png (scoreRequirements only)
        if scoreRequirements
            figs(end+1) = figure('Theme', 'light', 'Color', 'w', 'Position', [60 60 1000 900], 'Visible', figVis);
            ax = axes(figs(end));  hold(ax, 'on');
            rho = 90 + Gservice.lat;
            scatter(ax, rho.*cosd(Gservice.lon), rho.*sind(Gservice.lon), 10, Gservice.r2, 'filled', ...
                    'DisplayName', 'grid point (colour = R2 fraction)');
            wr = 90 + Gservice.worst.r2(1);
            plot(ax, wr*cosd(Gservice.worst.r2(2)), wr*sind(Gservice.worst.r2(2)), 'kp', 'MarkerSize', 16, ...
                 'MarkerFaceColor', 'y', 'DisplayName', sprintf('R2 floor %.3f at (%.1f, %g)', Gservice.r2Floor, Gservice.worst.r2));
            th = linspace(0, 360, 361);
            for latRing = [-85 -80 -75]
                plot(ax, (90 + latRing)*cosd(th), (90 + latRing)*sind(th), ':', 'Color', [0.5 0.5 0.5], 'HandleVisibility', 'off');
                text(ax, (90 + latRing)*cosd(45), (90 + latRing)*sind(45), sprintf('%d', latRing), 'FontSize', FS - 2);
            end
            axis(ax, 'equal');  axis(ax, [-16 16 -16 16]);  box(ax, 'on');
            set(ax, 'FontSize', FS, 'Color', 'w');
            xlabel(ax, 'colatitude x [deg] (longitude 0 to the right)');  ylabel(ax, 'colatitude y [deg]');
            cb = colorbar(ax);  cb.Label.String = 'R2: worst-day fraction with GDOP < 6';  cb.FontSize = FS;
            legend(ax, 'Location', 'southoutside', 'FontSize', FS);
            title(ax, sprintf('Service window (%g d): R2 per point, south pole, %s grid', serviceDays, scoreGrid), 'FontSize', FS + 1);
            exportgraphics(figs(end), fullfile(outDir, 'service_map.png'), 'Resolution', 150);
        end

        % propellant.png
        figs(end+1) = figure('Theme', 'light', 'Color', 'w', 'Position', [60 60 1280 600], 'Visible', figVis);
        tlp = tiledlayout(figs(end), 1, 2, 'TileSpacing', 'compact', 'Padding', 'compact');
        ax = nexttile(tlp);
        bar(ax, 1:6, [P.propKg], 0.6, 'FaceColor', cXfer, 'DisplayName', 'propellant');
        set(ax, 'XTick', 1:6, 'XTickLabel', compose('sc %d', 1:6), 'FontSize', FS, 'Color', 'w');
        ylabel(ax, 'propellant [kg]');  grid(ax, 'on');
        legend(ax, 'Location', 'northoutside', 'FontSize', FS);
        ax = nexttile(tlp);  axis(ax, 'off');
        rowsT = cell(8, 1);
        rowsT{1} = sprintf('%-3s %-4s %-7s %-4s %7s %7s %7s %7s', 'sc', 'slot', 'cell', 'wait', 't_f d', 'kg', 'dV', 'fly km');
        for kk = 1:6
            rowsT{kk + 1} = sprintf('%-3d %-4d %-7s %-4d %7.3f %7.3f %7.4f %7.1e', kk, S.assign(kk), ...
                     sprintf('%d,%d', S.iD(kk), S.iA(kk)), S.nRev(kk), P(kk).tf*dayNd, P(kk).propKg, ...
                     P(kk).dvKms, P(kk).gates.flyKm);
        end
        rowsT{8} = sprintf('total %.3f kg', sum([P.propKg]));
        text(ax, 0, 0.5, strjoin(rowsT, newline), 'FontName', 'Courier', 'FontSize', FS - 1, ...
             'Interpreter', 'none', 'VerticalAlignment', 'middle');
        title(tlp, 'Propellant and transfer summary (70 mN, Isp 900 s, 150 kg)', 'FontSize', FS + 1);
        exportgraphics(figs(end), fullfile(outDir, 'propellant.png'), 'Resolution', 150);
        if scoreRequirements, mapTxt = ', service_map.png'; else, mapTxt = ''; end
        fprintf('   schedule.png, assembly_3d.png, pole_gdop.png%s, propellant.png written\n', mapTxt);
    catch err
        if headless, close(figs(isgraphics(figs))); end
        rethrow(err);
    end
    if headless, close(figs(isgraphics(figs))); end
end
secWall(7) = toc(tSec);

%% ========================================================================
%  8. MOVIE
%% ========================================================================
tSec = tic;
movieDue = makeMovie && ~isempty(H) && ~isempty(pole);
if movieDue
    fprintf('\n8. MOVIE: rendered from the saved results after the save (end of section 9)\n');
elseif makeMovie
    fprintf('\n8. MOVIE: not rendered (the study stopped before the pole series)\n');
else
    fprintf('\n8. MOVIE: skipped (makeMovie = false)\n');
end
secWall(8) = toc(tSec);

%% ========================================================================
%  9. GATE TABLE, VERDICT, SAVE
%% ========================================================================
[~, hostName] = system('hostname');
run_meta = struct('date', char(datetime('now', 'Format', 'yyyy-MM-dd HH:mm')), 'matlab', version, ...
                  'host', strtrim(hostName), 'optimal_control', gitOC.head, ...
                  'optimal_controlDirty', gitOC.dirty, 'proj7', gitP7.head, 'proj7Dirty', gitP7.dirty, ...
                  'git', struct('optimal_control', gitOC, 'proj7', gitP7), ...
                  'knobs', knobs, 'reduced', {reduced}, 'sectionNames', {secName}, ...
                  'sectionWallS', secWall, 'studyWallS', toc(tStudy));
saveTarget = fullfile(outDir, 'lunar_pnt_mission_study.mat');
% fresh polishes are certified by THIS run: stamp their origin (the head
% carries '-dirty' when the working tree had modified tracked files)
here0 = struct('file', saveTarget, 'date', run_meta.date, 'head', gitOC.tag);
for kc = 1:numel(polishCache.origin)
    if isempty(polishCache.origin{kc}.file), polishCache.origin{kc} = here0; end
end
for kk = 1:6
    if strcmp(polishSource(kk).origin, 'fresh')
        polishSource(kk) = struct('sc', kk, 'origin', 'fresh', 'file', here0.file, 'date', here0.date, 'head', here0.head);
    end
end
run_meta.polishSource = polishSource;
% the VERDICT's disclosure suffix (reduced knobs, reused polishes, reused
% refusals), saved with the run so the movie and later readers see it
redTxt = '';
reu = polishSource(strcmp({polishSource.origin}, 'reused'));
if ~isempty(reu)
    srcs = unique(arrayfun(@(q) sprintf('%s, certified %s at %s', q.file, q.date, q.head), reu, ...
                           'UniformOutput', false), 'stable');
    redTxt = sprintf(' (REUSED polishes for sc %s from %s; not re-certified in this run)', ...
                     strtrim(sprintf('%d ', [reu.sc])), strjoin(srcs, ' and '));
end
if ~isempty(reusedRefusals)                           % said even when all six final polishes are fresh
    redTxt = sprintf('%s (exclusions reused: %d refused polish(es) from %s; not re-polished in this run)', ...
                     redTxt, size(reusedRefusals, 1), strjoin(unique(reusedRefusals(:, 2), 'stable').', ' and '));
end
if ~isempty(reduced), redTxt = sprintf('%s (reduced: %s)', redTxt, strjoin(reduced, ', ')); end
run_meta.verdictSuffix = redTxt;
run_meta.reusedRefusals = reusedRefusals;
save(saveTarget, 'S', 'S0', 'P', 'Pattempted', 'polishCache', 'H', 'Vpole', 'V', 'pole', 'Gservice', ...
     'gates', 'rounds', 'I2', 'fingerprint', 'run_meta', '-v7.3');

fprintf('\n== SUMMARY  lunar PNT mission: six spacecraft DRO -> tulip, south-pole GDOP, %g service days\n', serviceDays);
if ~isempty(P)
    fprintf('   last arrival: %.3f d unpolished (round 1) -> %.3f d realised after %d round(s)\n', ...
            S0.tLast*dayNd, max([P.tArr])*dayNd, numel(rounds));
end
fprintf('   wall per section [s]:');
for ks = 1:numel(secName), fprintf('  %s %.1f', secName{ks}, secWall(ks)); end
fprintf('  (total %.0f)\n', toc(tStudy));
fprintf('   %-4s %-11s %12s  %s\n', 'gate', 'status', 'value', 'detail');
for kg = 1:size(gates, 1)
    fprintf('   %-4s %-11s %12.5g  %s\n', gates{kg, 1}, gates{kg, 2}, gates{kg, 3}, gates{kg, 4});
end
failedIds = gates(~strcmp(gates(:, 2), 'PASS'), 1);
fprintf('   wrote %s\n', saveTarget);
if ~isempty(stopped)
    fprintf('   VERDICT: %s; not passed: %s%s\n', stopped, strjoin(failedIds.', ', '), redTxt);
elseif isempty(failedIds)
    fprintf('   VERDICT: PASS%s\n', redTxt);
else
    fprintf('   VERDICT: not passed: %s%s\n', strjoin(failedIds.', ', '), redTxt);
end
if movieDue
    movInfo = mission_movie(saveTarget, struct('outFile', fullfile(outDir, 'lunar_pnt_mission.mp4')));
    fprintf('   MOVIE: %s (%d frames, %d x %d, %.0f s wall)\n', movInfo.outFile, movInfo.nFrames, ...
            movInfo.frameSize, movInfo.wallS);
end


%% ------------------------------------------------------------------------
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

function G = setg(G, id, status, value, detail)
% SETG  Replace the row of gate id (append when absent).
%   INPUTS: as addg.  OUTPUTS: G.
row = find(strcmp(G(:, 1), id), 1);
if isempty(row), G = addg(G, id, status, value, detail); else, G(row, :) = {id, status, value, detail}; end
end

function s = setfield_(s, name, value)
% SETFIELD_  Copy of struct s with field name set.  INPUTS: s; name; value.  OUTPUTS: s.
s.(name) = value;
end

function v = gfield(s, name)
% GFIELD  s.(name), or NaN when s has no such field.  INPUTS: s struct; name.  OUTPUTS: v.
if isstruct(s) && isfield(s, name) && ~isempty(s.(name)), v = s.(name); else, v = NaN; end
end

function s = conj_status(Gs)
% CONJ_STATUS  PC over the six gate structs, three outcomes as in
%   transfer_study: PASS when every conjugate test passed and covered the
%   arc; FAIL when one was reproduced, covered and refuted; otherwise
%   UNRESOLVED (not reproduced, not covered, or not run).
%   INPUTS: Gs [1 x 6] polish gate structs.  OUTPUTS: s char.
if all(arrayfun(@(g) isequal(g.passConj, true), Gs))
    s = 'PASS';
elseif any(arrayfun(@(g) isequal(g.reproduced, true) && g.conjCovered == 1 && g.conjPass == 0, Gs))
    s = 'FAIL';
else
    s = 'UNRESOLVED';
end
end

function fp = catalog_fingerprint(M, thr, polishGates)
% CATALOG_FINGERPRINT  The identity a reused results file must match:
%   catalog name and created date, MD5 of sh.z8 and sh.tf_nd, the regime
%   constants, the slot phasing, the spacecraft DRO phases, MD5 of the DRO
%   and tulip orbit handles sampled at 64 phases (droMd5, tulipMd5), the
%   polish thresholds and polishGates.
%   INPUTS: M mission inputs; thr polish thresholds struct; polishGates
%   logical.  OUTPUTS: fp struct.
s64 = (0:63)/64;
fp = struct('catalogName', char(string(M.cat.name)), 'catalogCreated', char(string(M.cat.created)), ...
            'z8md5', md5_of(M.sh.z8), 'tfmd5', md5_of(M.sh.tf_nd), ...
            'thrustN', M.thrustN, 'ispS', M.ispS, 'm0kg', M.m0kg, 'muStar', M.muStar, ...
            'lStar', M.lStar, 'tStar', M.tStar, 'dtau6', M.dtau6, ...
            'droPhase0', M.droPhase0, 'droMd5', md5_of(M.depState(s64)), ...
            'tulipMd5', md5_of(M.arrState(s64)), ...
            'thresholds', thr, 'polishGates', logical(polishGates));
end

function g = git_state(gitBin, folder)
% GIT_STATE  Read-only git state of the repository holding a folder.
%   INPUTS:  gitBin char (git executable); folder char.
%   OUTPUTS: g struct .head (short HEAD, 'unknown' on failure) .root
%            (top-level folder) .dirty (logical: tracked files modified)
%            .nDirty (count of modified tracked paths) .tag (head, with
%            '-dirty' appended when dirty) .text (one line for the log).
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

function h = md5_of(x)
% MD5_OF  Hex MD5 of a numeric array's bytes (size included).
%   INPUTS: x numeric array.  OUTPUTS: h char [1 x 32].
md = java.security.MessageDigest.getInstance('MD5');
md.update(typecast([double(size(x)), double(x(:)).'], 'uint8'));
h = lower(reshape(dec2hex(typecast(md.digest(), 'uint8'), 2).', 1, []));
end

function L = load_reuse(file, fpNow)
% LOAD_REUSE  Load a results file for reuse, refusing (lunar_pnt_mission_
%   study:reuse) when it has no fingerprint or no polish origins, or when
%   any fingerprint field differs from this run's.
%   INPUTS: file char; fpNow struct.  OUTPUTS: L struct (.polishCache .run_meta).
if ~isfile(file)
    error('lunar_pnt_mission_study:reuse', 'reuse file not found: %s', file);
end
w = whos('-file', file);
if ~any(strcmp({w.name}, 'fingerprint'))
    error('lunar_pnt_mission_study:reuse', '%s carries no fingerprint (made before guarded reuse); refusing', file);
end
L = load(file, 'polishCache', 'run_meta', 'fingerprint');
if ~isfield(L.polishCache, 'origin')
    error('lunar_pnt_mission_study:reuse', '%s carries no polish origins; refusing', file);
end
fn = fieldnames(fpNow);
bad = fn(~cellfun(@(f) isfield(L.fingerprint, f) && isequal(L.fingerprint.(f), fpNow.(f)), fn));
if ~isempty(bad)
    error('lunar_pnt_mission_study:reuse', '%s fingerprint differs from this run in: %s; refusing', ...
          file, strjoin(bad.', ', '));
end
end

function why = reuse_entry_problem(Pk, tArr, rvf, rv0, thr)
% REUSE_ENTRY_PROBLEM  Why a cached polish cannot stand for the current
%   candidate ('' when it can): arrival epoch, epoch consistency, the
%   current slot state, the current DRO state at the cached departure
%   epoch, the current thresholds.
%   INPUTS: Pk polish struct; tArr scalar [ND]; rvf [6 x 1] slot state at
%   tArr; rv0 [6 x 1] the spacecraft's DRO state at Pk.tDep; thr struct.
%   OUTPUTS: why char.
why = '';
if ~(abs(Pk.tArr - tArr) < 1e-12)
    why = sprintf('cached tArr differs by %.3g ND', Pk.tArr - tArr);
elseif ~(abs(Pk.tDep + Pk.tf - Pk.tArr) < 1e-12)
    why = sprintf('cached tDep + tf - tArr = %.3g ND', Pk.tDep + Pk.tf - Pk.tArr);
elseif ~(max(abs(Pk.rvf(:) - rvf(:))) < 1e-12)
    why = sprintf('cached rvf differs from the current slot state by %.3g ND', max(abs(Pk.rvf(:) - rvf(:))));
elseif ~(max(abs(Pk.rv0(:) - rv0(:))) < 1e-12)
    why = sprintf('cached rv0 differs from the current DRO state at its tDep by %.3g ND', ...
                  max(abs(Pk.rv0(:) - rv0(:))));
elseif ~isequal(Pk.thresholds, thr)
    why = 'cached thresholds differ';
end
end

function s = chunk_text(sel)
% CHUNK_TEXT  'all' or the chunk list.  INPUTS: sel [1 x c].  OUTPUTS: s char.
if isempty(sel), s = 'all'; else, s = strtrim(sprintf('%d ', sel)); end
end

function s = knob_text(v)
% KNOB_TEXT  A knob value as short text.  INPUTS: v any.  OUTPUTS: s char.
if ischar(v) || isstring(v)
    [~, nm, ex] = fileparts(char(v));
    s = [nm ex];
elseif islogical(v)
    s = tern(v, 'true', 'false');
else
    s = mat2str(v, 6);
end
end
