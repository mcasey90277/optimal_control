function run_gto_catalog(sheetSel, maxCellsIn, batchSecIn, rungsIn)
%% Purpose:
%
%   FRONT DOOR for the GTO -> TULIP costate catalog: the proven three-gate
%   pipeline (direct solve -> ms_tfmin refinement -> pumpkyn tfMin
%   acceptance), walked down a thrust ladder, swept over (GTO departure
%   orientation x tulip petal count) sheets. THE SWATH DRIVER:
%   run_gto_catalog() solves all 16 sheets; run_gto_catalog(sheetSel)
%   solves the subset indexed by sheetSel. Resume/regeneration of unsolved
%   cells is automatic (thrust_ladder_library revisits no-OK cells on
%   every call; a batched shell driver is how an unattended run survives a
%   hung solver call -- see run_gto_batched.sh).
%
%   Endpoints are built through costate_common/get_family_orbit, so this
%   front door differs from the DRO/HALO catalogs' ONLY in the departure-
%   family parameters ('gto', keyed by orientDeg) -- the engine is shared,
%   unmodified (clone of HALO_tulip/run_halo_catalog.m's 'all' stage).
%
%   TWO-POINT TOP-RUNG MULTISTART (diagnostic-A, 2026-08-26; scoping fix
%   2026-08-26 review): a cold GTO cell has TWO failure basins a single
%   thrust_ladder_library call cannot escape (see diagnostic-A-report.md)
%   -- a mass-infeasible tf0=4.0 seed, and (even with thrLock pinning the
%   throttle) an occasional bad basin at tf0=0.30. Every sheet is
%   therefore walked TWICE per call: pass A at thrLock=true, tf0=0.30
%   (the higher-hit-rate recipe), then pass B at thrLock=true, cold tf0
%   (casadi_mintime_dro's own internal 4.0 ND default), SCOPED to the
%   mop-up set computed from the sheet file after pass A returns -- cells
%   with ATT>=1 and no OK rung, i.e. cells pass A itself just attempted
%   and could not clear (opts.cells, the engine's own index-pair subset
%   input). This scoping matters because thrust_ladder_library's todo
%   list is windowed by maxCells: without it, once pass A's window is
%   full, pass B's OWN resume-filtered todo would slide past pass A's
%   cells into FRESH ones that pass A never touched, burning their first
%   (and, at the default maxAtt=2, HALF their total) attempts on the cold
%   recipe instead of the warmer one.
%
%   Multistart here is FIRST-SUCCESS-WINS, NOT a true min-tf keep-best:
%   pass B never re-solves a cell pass A already cleared (thrust_ladder_
%   library's resume filter drops any cell with `any(OK(iD,iA,:))` before
%   pass B ever sees it), so there is no comparison of the two seeds'
%   t_f on a cell where pass A succeeds. This is justified by the
%   diagnostic's own measurement, not assumed: wherever both routes
%   converged on the same cell, they agreed to 5-6 significant digits in
%   t_f (diagnostic-A-report.md Section 3, "Independent cross-
%   validation") -- i.e. both seeds land the SAME basin whenever both are
%   reachable at all, so "first success" and "smaller t_f" pick the same
%   answer in every case measured. If a future family (different
%   endpoints, different thrust range) shows seed-dependent basins where
%   tf0=0.30 and cold tf0 both succeed but disagree, a true comparator
%   (solve both, keep the smaller t_f) is required here -- this driver
%   does not have one.
%
%   SHEET KEY RULE: opts.tauDRO = orientDeg for every sheet, NOT a Kepler
%   period -- the picker (and every consumer downstream) selects sheets by
%   nearest (tauDRO, Np); the GTO family's own Kepler period rides in
%   depParams.orientDeg / get_family_orbit's info, never in tauDRO.
%
%   DEPARTURE-PHASE GRID CAVEAT: thrust_ladder_library only accepts a
%   UNIFORM departure-phase grid (opts.sD0 + opts.nD -> mod(sD0+(0:nD-1)/
%   nD,1)); it has no explicit-fraction-vector input. A review requested a
%   grid uniform in TRUE ANOMALY (denser near perigee on this e~0.72
%   ellipse) mapped through Kepler's equation to non-uniform time
%   fractions -- that would require an engine change, which this task does
%   NOT make (house rule: house engines are not modified in place to serve
%   one campaign). sD0=0, nD=12 UNIFORM IN TIME FRACTION is used instead;
%   this still lands cleanly on perigee (index 1) and apogee (index 7 of
%   12, since 6/12=0.5) so both regimes are sampled, just not with extra
%   density at perigee. See task-3-report.md.
%
%% Inputs:
%
%  sheetSel                 [1 x k]                 Index vector into the
%                                                   16-sheet list (default
%                                                   1:16, i.e. all sheets)
%
%  maxCellsIn               double                  Optional batch
%                                                   override: fresh cells
%                                                   per sheet this call
%
%  batchSecIn                double                 Optional batch
%                                                   override: clean-exit
%                                                   wall budget (s), shared
%                                                   across the selected
%                                                   sheets this call
%
%  rungsIn                  [1 x r]                 Optional thrust-ladder
%                                                   override (N, high to
%                                                   low). Default (omitted
%                                                   or empty) is the full
%                                                   9-rung [15 12 10 7 5 3
%                                                   2 1.5 1] -- byte-
%                                                   identical to the prior
%                                                   no-arg behavior.
%                                                   Deliverable-7 v1 fleet
%                                                   calls this with the
%                                                   5-rung truncation
%                                                   [15 12 10 7 5] (see
%                                                   task-4-brief.md
%                                                   amendment, 2026-08-27
%                                                   re-census adjudication:
%                                                   pilot measured a
%                                                   closure wall below 5 N
%                                                   -- 0 entries at
%                                                   3/2/1.5/1 N by both
%                                                   routes).
%
%% Outputs:
%
%   none (files under results/catalog/, relative to this file's folder):
%     gto_o<ORIENT3>_Np<NP>.mat + _progress.txt   one pair per sheet
%     combined_progress.txt                       driver-level log +
%                                                  completion marker
%
%   Unattended: ./run_gto_batched.sh is the batched driver (OS-kill
%   backstop for a hung solver call, per the matlab-campaign pattern).
%
%% Revision History:
%  M. Casey                                                   (c) 08/26/2026
%  Copyright Coorbital Inc.
%% ------------------------ Begin Code Sequence ---------------------------

if ~exist('sheetSel','var') || isempty(sheetSel), sheetSel = 1:16; end
maxCells = inf;  batchSec = inf;
if exist('maxCellsIn','var') && ~isempty(maxCellsIn), maxCells = maxCellsIn; end
if exist('batchSecIn','var') && ~isempty(batchSecIn), batchSec = batchSecIn; end
if ~exist('rungsIn','var'), rungsIn = []; end

%% ===================== ADJUSTABLE PARAMETERS ============================

%% Departure: GTO pseudo-family, selected by orientation (sheet axis 1):
      orients = [0 90 180 270];               % deg

%% Arrival: tulip petal counts (sheet axis 2; periods locked by resonance):
       Nps     = [5 7 9 12];
       tulipPm = -1;

%% Phasing grid per sheet (departure x arrival; see the module-header
%% caveat above -- uniform-only, sD0=0/sA0=0 origins):
           nD = 12;
           nA = 6;
           sD0 = 0;
           sA0 = 0;

%% Thrust ladder (N, high to low -- the proven set; rungsIn overrides):
        rungs = [15 12 10 7 5 3 2 1.5 1];
if ~isempty(rungsIn), rungs = rungsIn; end
         ispS = 1710;                % specific impulse (s)
         m0kg = 150;                 % initial mass (kg)

%% Solver / gates (diagnostic-A recipe, 2026-08-26 -- the cold N=400/
%% maxIter=3000 defaults that work for DRO/HALO/DPO leave 67/72 GTO cells
%% unreached: see diagnostic-A-report.md):
            N = 400;                 % collocation intervals -- UNCHANGED,
                                      % 4.5 deg/interval on the true 0.27-
                                      % 0.30 ND solution is plenty
      floorKm = 500;                 % lunar altitude floor (km)
      maxIter = 6000;                % 3000 stopped a right-family cell one
                                      % iteration budget short (defect 1.3e-11)
    maxCpuSec = 600;                 % 300 truncates otherwise-healthy solves
       gateKm = 100;                 % flown-arrival gate (km)
       accTol = 1e-6;                % tfMin acceptance
      thrLock = true;                % THE decisive knob: pins thr==1,
                                      % deleting the mass-floor coast
                                      % manifold a free throttle falls into
                                      % on GTO's sub-revolution transfers
     tf0Fleet = 0.30;                % warm top-rung seed for pass A, ~0.72x
                                      % the thrLock hard cap 0.95*c/T @ 15N

%% ======================= END ADJUSTABLE PARAMETERS ======================

% Review fix (c): assert the departure-phase grid the engine will build
% internally has no duplicate phase (no sD=1 aliasing sD=0):
sDCheck = mod(sD0 + (0:nD-1)/nD, 1);
assert(numel(unique(sDCheck)) == nD, 'run_gto_catalog:sDgrid', ...
    'departure-phase grid has a duplicate phase');

setup_paths();
here = fileparts(mfilename('fullpath'));
catDir = fullfile(here, 'results', 'catalog');
if ~isfolder(catDir), mkdir(catDir); end

sheets = {};
for od = orients
    for np = Nps
        sheets{end+1} = struct('orientDeg', od, 'Np', np, ...
            'tag', sprintf('gto_o%03d_Np%d', od, np));
    end
end
nSheets = numel(sheets);
assert(nSheets == 16, 'run_gto_catalog:sheetCount', ...
    'expected 16 sheets, got %d', nSheets);

sheetSel = sheetSel(sheetSel >= 1 & sheetSel <= nSheets);
assert(~isempty(sheetSel), 'run_gto_catalog:sheetSel', ...
    'sheetSel selects no valid sheet (1..%d)', nSheets);

tAll = tic;
allDone = true;
logComb = fullfile(catDir, 'combined_progress.txt');
for ks = sheetSel
    if toc(tAll) > batchSec, allDone = false; break, end
    s = sheets{ks};
    outMat  = fullfile(catDir, [s.tag '.mat']);
    logFile = fullfile(catDir, [s.tag '_progress.txt']);
    fprintf('\n=== GTO sheet %d/%d: %s (orient=%d deg, Np=%d) ===\n', ...
            ks, nSheets, s.tag, s.orientDeg, s.Np);
    % GUARD (Important review finding, task-4 fix, 2026-08-27): a sheet's
    % on-disk rungs MUST match this call's rungs before the engine ever
    % touches it. thrust_ladder_library's resume path is a bare size
    % compare (isequal(size(Q.OK), [nD nA nR])); a caller that omits
    % rungsIn on a truncated sheet (silently reverting to the 9-rung
    % default) would fail that compare, the engine would treat every
    % cell as fresh (ATT reset to 0), and the very next save() would
    % overwrite the sheet's accumulated entries with all-NaN/all-false
    % arrays sized to the WRONG rung count -- the exact silent-wipe
    % mechanism the pilot trim (task 4) defused, standing re-armed for
    % every future call. This turns that into a loud stop instead.
    if isfile(outMat)
        Qchk = load(outMat, 'rungs');
        if ~isequal(Qchk.rungs(:)', rungs(:)')
            error('run_gto_catalog:rungMismatch', ...
                ['sheet %s on disk has rungs=[%s] but this call requests ' ...
                 'rungs=[%s] -- refusing to run (the engine''s resume ' ...
                 'check would silently re-solve from scratch and ' ...
                 'overwrite the sheet''s accumulated entries). Remedy: ' ...
                 'pass rungsIn matching the sheet, or move the sheet ' ...
                 'aside deliberately.'], s.tag, num2str(Qchk.rungs), ...
                num2str(rungs));
        end
    end
    optsCommon = struct('rungs',rungs, 'ispS',ispS, 'm0kg',m0kg, ...
        'N',N, 'floorKm',floorKm, 'maxIter',maxIter, 'maxCpuSec',maxCpuSec, ...
        'gateKm',gateKm, 'nD',nD, 'nA',nA, 'sD0',sD0, 'sA0',sA0, ...
        'depFamily','gto', 'depParams',struct('orientDeg',s.orientDeg), ...
        'arrFamily','tulip', 'arrParams',struct('Np',s.Np,'pm',tulipPm), ...
        'NpTulip',s.Np, 'tauTulip',2*pi*(s.Np-2)/(s.Np-1), 'pmTulip',tulipPm, ...
        'tauDRO',s.orientDeg, ...          % SHEET KEY RULE: orientDeg
        'accTol',accTol, 'maxCells',maxCells, ...
        'thrLock',thrLock, 'logFile',logFile);
    % TWO-POINT MULTISTART, FIRST-SUCCESS-WINS (diagnostic-A Step 2, scoping
    % fix 2026-08-26 review -- see the module header for the full
    % rationale). thrLock stays true on BOTH passes; only the top-rung tf0
    % seed varies.
    %
    % Pass A -- tf0=0.30, the higher-hit-rate seed (measured 6/9, 8.6-58s):
    optsA = optsCommon;
    optsA.tf0      = tf0Fleet;
    optsA.batchSec = max(0, batchSec - toc(tAll));
    % FAIL LOUD (fix round 2, 2026-08-27 review): no silent continue past a
    % sheet whose engine call errors -- log which sheet+pass failed to both
    % stdout and the combined log, then rethrow so the shell driver's own
    % non-zero-exit abort path fires instead of the loop quietly moving on
    % to the next sheet with this one left in an unknown state.
    try
        thrust_ladder_library(outMat, optsA);
    catch ME
        failMsg = sprintf('FATAL: sheet %s pass A (thrust_ladder_library) errored -- %s', ...
            s.tag, ME.message);
        fprintf(2, '%s\n', failMsg);
        fidF = fopen(logComb, 'a');
        if fidF > 0, fprintf(fidF, '%s\n', failMsg); fclose(fidF); end
        rethrow(ME);
    end
    % Pass B -- cold tf0 (casadi_mintime_dro's own internal default, 4.0
    % ND), SCOPED to the mop-up set: cells pass A itself just attempted
    % (ATT>=1) and could not clear (no OK rung). Computed fresh from the
    % sheet file pass A just wrote, NOT from P.ATT/P.OK in memory, so the
    % scope is exactly what pass A touched this call -- without this, a
    % maxCells-windowed pass B would resume-filter past pass A's cells
    % into fresh, never-attempted ones and burn their first attempt cold.
    % An empty mop-up set means pass A cleared (or exhausted) everything
    % it touched, so pass B is skipped for this sheet this call.
    if toc(tAll) <= batchSec
        mopCells = zeros(0,2);
        if isfile(outMat)
            S = load(outMat, 'ATT', 'OK');
            [mD, mA] = find(S.ATT >= 1 & ~any(S.OK,3));
            mopCells = [mD, mA];
        end
        if ~isempty(mopCells)
            optsB = optsCommon;
            optsB.cells    = mopCells;
            optsB.batchSec = max(0, batchSec - toc(tAll));
            try
                thrust_ladder_library(outMat, optsB);
            catch ME
                failMsg = sprintf('FATAL: sheet %s pass B (thrust_ladder_library) errored -- %s', ...
                    s.tag, ME.message);
                fprintf(2, '%s\n', failMsg);
                fidF = fopen(logComb, 'a');
                if fidF > 0, fprintf(fidF, '%s\n', failMsg); fclose(fidF); end
                rethrow(ME);
            end
        end
    end
    % OPEN-CELL SWEEP (fix round 2, 2026-08-27 review): reload THIS sheet's
    % just-saved data FRESH FROM DISK rather than trusting the in-memory
    % struct thrust_ladder_library returned. Both passes already save to
    % outMat before returning (on the happy path this is not a behavior
    % change) -- it is the review-required guarantee that the completion
    % marker below can NEVER be based on stale/in-memory state, only on
    % what is actually persisted for this sheet. Guard the degenerate case
    % where outMat still does not exist: a brand-new sheet whose pass A
    % ran out of this call's leftover batchSec budget before its very
    % first per-cell save (thrust_ladder_library's only save() sits INSIDE
    % its per-cell loop, gated by its own batchSec check, so zero cells
    % attempted means zero bytes written) -- that sheet is trivially open.
    if isfile(outMat)
        Qsweep = load(outMat, 'OK', 'ATT');
        open_ = false;
        for iD = 1:nD
            for iA = 1:nA
                if ~any(Qsweep.OK(iD,iA,:)) && Qsweep.ATT(iD,iA) < 2, open_ = true; end
            end
        end
    else
        open_ = true;
    end
    if open_, allDone = false; end
end

if allDone && numel(sheetSel) == nSheets
    fid = fopen(logComb, 'a');
    fprintf(fid, 'CATALOG ALL SHEETS COMPLETE\n');  fclose(fid);
elseif allDone
    fid = fopen(logComb, 'a');
    fprintf(fid, 'SELECTED SHEETS COMPLETE (%s)\n', mat2str(sheetSel));
    fclose(fid);
end
end
