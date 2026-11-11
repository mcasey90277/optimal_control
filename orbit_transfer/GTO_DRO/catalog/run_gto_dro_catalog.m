function plan = run_gto_dro_catalog(sheetSel, maxCellsIn, batchSecIn, rungsIn, optsIn)
%% Purpose:
%
%   FRONT DOOR for the GTO -> DRO costate catalog (stage 2 of the GTO -> DRO
%   min-time campaign, spec docs/superpowers/specs/2026-11-09-gto-dro-
%   mintime-design.md section 10): the shared three-gate ladder engine
%   (direct solve -> ms_tfmin refinement -> pumpkyn tfMin acceptance)
%   walked down the thrust ladder over (GTO departure orientation x DRO
%   period) sheets. Clone of GTO_tulip/catalog/run_gto_catalog.m with the
%   DRO arrival, the Earth-altitude floor (stage 1 measured a floored-free
%   solve passing 1722 km from Earth's centre) and a planOnly mode. The
%   two-point top-rung multistart (pass A tf0 warm, pass B cold mop-up,
%   first-success-wins) and the SHEET KEY RULE (opts.tauDRO = orientDeg)
%   are inherited unchanged; see that file's header for their rationale.
%
%   WHAT THE FLEET DOES NOT DO (stage 1 findings, GTO_DRO/README.md): it
%   does not select the top-rung basin (the floored cold start can land on
%   an extremal the conjugate test refutes) and it does not reach the
%   perigee column at 15 N. Both are repaired AFTER the fleet by
%   lib/repair_gto_dro_cells (reverse warm start, phase continuation)
%   through the engine's seedTop/redo options, on the sheets this front
%   door writes. The fleet removes ONLY the first flown gate (directGateKm =
%   inf: the open-loop flight of the N = 400 direct control is not a filter
%   on these transfers); the refined gate (gateKm = 100), tfMin acceptance
%   and the ms hard timeout stay.
%
%% Inputs:
%
%  sheetSel                 [1 x k]                 Index into the 16-sheet
%                                                   list (default 1:16);
%                                                   orientation outer
%                                                   [0 90 180 270], tau
%                                                   inner [0.5 1 2 3], so
%                                                   sheet 2 = (0 deg, tau 1)
%                                                   is the pilot
%  maxCellsIn               double                  Fresh cells per sheet
%                                                   this call [inf]
%  batchSecIn               double                  Clean-exit wall budget
%                                                   (s) shared by the
%                                                   selected sheets [inf]
%  rungsIn                  [1 x r]                 Thrust ladder override,
%                                                   N high to low
%                                                   [15 12 10 7 5]
%  optsIn                   struct (optional)       .planOnly [false]
%                                                   return the plan without
%                                                   solving; .tf0 [0.30]
%                                                   pass A top-rung seed;
%                                                   .catDir [results/
%                                                   catalog next to this
%                                                   folder's parent];
%                                                   .directGateKm [inf]
%                                                   first flown gate on the
%                                                   direct control (km);
%                                                   .msHardCapSec [300]
%                                                   hard wall-clock cap (s)
%                                                   on each ms_tfmin
%                                                   refinement attempt,
%                                                   passed to the engine
%                                                   (0 = off). Pilot
%                                                   2026-11-10: one cell
%                                                   hung 55 min inside
%                                                   ms_tfmin (pumpkyn
%                                                   tfMinProp restarting a
%                                                   failed ode45 forever)
%                                                   and the tfMin
%                                                   acceptance. The fence
%                                                   pool is opened ONCE per
%                                                   call (capped_pool(1))
%                                                   and reused by every
%                                                   engine call, deleted on
%                                                   return if this call
%                                                   opened it;
%                                                   .requireFence [true]
%                                                   when msHardCapSec > 0
%                                                   and no pool can be
%                                                   opened, REFUSE
%                                                   (run_gto_dro_catalog:
%                                                   noFence) rather than
%                                                   run a fleet unfenced;
%                                                   false = run UNFENCED
%
%% Outputs:
%
%  plan                     [k x 1] struct          .tag .orientDeg .tau
%                                                   .outMat .logFile .optsA
%                                                   .fence (fence status
%                                                   line, also written to
%                                                   combined_progress.txt)
%                                                   .optsB (pass B .cells
%                                                   empty); optsA.batchSec,
%                                                   optsB.batchSec and
%                                                   optsB.cells are set at
%                                                   run time, NOT in the
%                                                   plan; files under
%                                                   catDir: <tag>.mat +
%                                                   <tag>_progress.txt,
%                                                   combined_progress.txt
%
%% Revision History:
%  M. Casey                                                   (c) 11/10/2026
%  M. Casey  msHardCapSec knob + ode45 warning silenced        (c) 11/10/2026
%  Copyright Coorbital Inc.
%% ------------------------ Begin Code Sequence ---------------------------

if ~exist('sheetSel','var') || isempty(sheetSel), sheetSel = 1:16; end
maxCells = inf;  batchSec = inf;
if exist('maxCellsIn','var') && ~isempty(maxCellsIn), maxCells = maxCellsIn; end
if exist('batchSecIn','var') && ~isempty(batchSecIn), batchSec = batchSecIn; end
if ~exist('rungsIn','var'), rungsIn = []; end
if ~exist('optsIn','var') || isempty(optsIn), optsIn = struct(); end
planOnly = fdef(optsIn, 'planOnly', false);

%% ===================== ADJUSTABLE PARAMETERS ============================

%% Departure: GTO pseudo-family, selected by orientation (sheet axis 1):
      orients = [0 90 180 270];               % deg

%% Arrival: DRO periods (sheet axis 2; ND, tau 0.5/1/2/3 = 2.22/4.43/8.87/13.30 d):
         taus = [0.5 1 2 3];

%% Phasing grid per sheet (departure x arrival; uniform in time fraction):
           nD = 12;
           nA = 6;
          sD0 = 0;
          sA0 = 0;

%% Thrust ladder (N, high to low); rungsIn overrides, and must match any
%% sheet already on disk (guard below):
        rungs = [15 12 10 7 5];
if ~isempty(rungsIn), rungs = rungsIn; end
         ispS = 1710;                % specific impulse (s)
         m0kg = 150;                 % initial mass (kg)

%% Solver / gates (the stage 1 recipe, GTO_DRO/README.md):
            N = 400;                 % collocation intervals
      floorKm = 500;                 % lunar altitude floor (km)
 floorEarthKm = 200;                 % Earth altitude floor (km), nodes + midpoints
      maxIter = 6000;                % IPOPT iteration cap
    maxCpuSec = 600;                 % IPOPT CPU cap per rung (s)
       gateKm = 100;                 % refined flown-arrival gate (km)
 directGateKm = fdef(optsIn, 'directGateKm', inf);   % first flown gate on the DIRECT control (km); inf = none. The open-loop flight of the N = 400 direct control is not a filter on these transfers (pilot 2026-11-10: 54 ladders stopped at 1000 km; probe with no direct gate: 14 -> 26 of 30 rungs, misses up to 42,964 km closed by refinement); the refined gate, tfMin acceptance and the ms fence are the gates
       accTol = 1e-6;                % tfMin acceptance
 requireFence = fdef(optsIn, 'requireFence', true);  % msHardCapSec > 0 but no pool: refuse (true) or run UNFENCED (false)
 msHardCapSec = fdef(optsIn, 'msHardCapSec', 300);    % hard cap (s) per ms_tfmin attempt, run on a parallel worker (engine opts.msHardCapSec); 0 = off. msWallS (120 s) is soft and cannot stop a stuck residual evaluation
      thrLock = true;                % pin thr == 1 (classical min-time form)
     tf0Fleet = fdef(optsIn, 'tf0', 0.30);   % pass A top-rung seed (ND); the pilot re-bases it

%% ======================= END ADJUSTABLE PARAMETERS ======================

% Assert the departure-phase grid the engine will build internally has no
% duplicate phase (no sD=1 aliasing sD=0):
sDCheck = mod(sD0 + (0:nD-1)/nD, 1);
assert(numel(unique(sDCheck)) == nD, 'run_gto_dro_catalog:sDgrid', ...
    'departure-phase grid has a duplicate phase');

here = fileparts(mfilename('fullpath'));                 % GTO_DRO/catalog
run(fullfile(fileparts(here), 'setup_gto_dro_paths.m'));
catDir = fdef(optsIn, 'catDir', fullfile(fileparts(here), 'results', 'catalog'));
if ~planOnly && ~isfolder(catDir), mkdir(catDir); end

sheets = {};
for od = orients
    for ta = taus
        sheets{end+1} = struct('orientDeg', od, 'tau', ta, ...
            'tag', sprintf('gto_dro_o%03d_tau%g', od, ta));
    end
end
nSheets = numel(sheets);
assert(nSheets == 16, 'run_gto_dro_catalog:sheetCount', ...
    'expected 16 sheets, got %d', nSheets);

sheetSel = sheetSel(sheetSel >= 1 & sheetSel <= nSheets);
assert(~isempty(sheetSel), 'run_gto_dro_catalog:sheetSel', ...
    'sheetSel selects no valid sheet (1..%d)', nSheets);

logComb = fullfile(catDir, 'combined_progress.txt');
% THE FENCE (msHardCapSec > 0): one pool for the whole call, opened AFTER
% setup_gto_dro_paths so the workers inherit pumpkyn, costate_common,
% oclib and CasADi; every engine call finds it through current_pool and
% keeps it. A fleet must not run unfenced by accident: with requireFence
% (default) an unavailable pool is a refusal, not a fallback.
fenceLine = 'ms fence OFF (msHardCapSec = 0)';
if planOnly
    fenceLine = 'ms fence not opened (planOnly)';
elseif msHardCapSec > 0
    % WARNING FLOOD: a stuck flight emits one ode45 tolerance warning per
    % restart (2.5 million lines, 1.29 GB of log in 55 min on the pilot).
    % Off for this call when the cap is on, restored on exit; fenced
    % workers silence it themselves.
    odeWarn = warning('off', 'MATLAB:ode45:IntegrationTolNotMet');
    restoreOdeWarn = onCleanup(@() warning(odeWarn));
    hadPool = ~isempty(current_pool());
    fenceWhy = '';
    try
        fencePool = capped_pool(1);
    catch ME
        fencePool = [];  fenceWhy = ME.message;
    end
    if isempty(fencePool)
        if isempty(fenceWhy), [fenceWhy, ~] = lastwarn; end
        fenceLine = sprintf('ms fence UNAVAILABLE -- running UNFENCED: %s', fenceWhy);
    else
        fenceLine = sprintf('ms fence ACTIVE (pool of %d)', fencePool.NumWorkers);
        if ~hadPool, closeFence = onCleanup(@() delete_pool(fencePool)); end
    end
    fidF = fopen(logComb, 'a');
    if fidF > 0, fprintf(fidF, '=== %s  %s\n', char(datetime('now', 'Format', 'yyyy-MM-dd HH:mm:ss')), fenceLine); fclose(fidF); end
    fprintf('%s\n', fenceLine);
    if isempty(fencePool) && requireFence
        error('run_gto_dro_catalog:noFence', ...
              ['msHardCapSec = %g requested but no parallel pool could be opened (%s); refusing to run ' ...
               'unfenced (pass optsIn.requireFence = false to allow it)'], msHardCapSec, fenceWhy);
    end
end

plan = struct('tag',{},'orientDeg',{},'tau',{},'fence',{},'outMat',{},'logFile',{},'optsA',{},'optsB',{});
tAll = tic;
allDone = true;
for ks = sheetSel
    if ~planOnly && toc(tAll) > batchSec, allDone = false; break, end
    s = sheets{ks};
    outMat  = fullfile(catDir, [s.tag '.mat']);
    logFile = fullfile(catDir, [s.tag '_progress.txt']);
    optsCommon = struct('rungs',rungs, 'ispS',ispS, 'm0kg',m0kg, ...
        'N',N, 'floorKm',floorKm, 'floorEarthKm',floorEarthKm, ...
        'maxIter',maxIter, 'maxCpuSec',maxCpuSec, ...
        'gateKm',gateKm, 'directGateKm',directGateKm, 'msHardCapSec',msHardCapSec, 'nD',nD, 'nA',nA, 'sD0',sD0, 'sA0',sA0, ...
        'depFamily','gto', 'depParams',struct('orientDeg',s.orientDeg), ...
        'arrFamily','dro', 'arrParams',struct('tau',s.tau), ...
        'tauDRO',s.orientDeg, ...          % SHEET KEY RULE: orientDeg, never the DRO period
        'accTol',accTol, 'maxCells',maxCells, ...
        'thrLock',thrLock, 'logFile',logFile);
    optsA = optsCommon;  optsA.tf0 = tf0Fleet;
    optsB = optsCommon;  optsB.cells = zeros(0,2);
    plan(end+1,1) = struct('tag', s.tag, 'orientDeg', s.orientDeg, 'tau', s.tau, 'fence', fenceLine, ...
        'outMat', outMat, 'logFile', logFile, 'optsA', optsA, 'optsB', optsB);
    if planOnly, continue, end
    fprintf('\n=== GTO -> DRO sheet %d/%d: %s (orient=%d deg, tau=%g) ===\n', ...
            ks, nSheets, s.tag, s.orientDeg, s.tau);
    % GUARD: a sheet's on-disk rungs MUST match this call's rungs before the
    % engine ever touches it. thrust_ladder_library's resume path is a bare
    % size compare; a mismatched call would treat every cell as fresh (ATT
    % reset to 0) and the next save() would overwrite the accumulated
    % entries with arrays sized to the WRONG rung count. Loud stop instead.
    if isfile(outMat)
        Qchk = load(outMat, 'rungs');
        if ~isequal(Qchk.rungs(:)', rungs(:)')
            error('run_gto_dro_catalog:rungMismatch', ...
                ['sheet %s on disk has rungs=[%s] but this call requests ' ...
                 'rungs=[%s] -- refusing to run (the engine''s resume ' ...
                 'check would silently re-solve from scratch and ' ...
                 'overwrite the sheet''s accumulated entries). Remedy: ' ...
                 'pass rungsIn matching the sheet, or move the sheet ' ...
                 'aside deliberately.'], s.tag, num2str(Qchk.rungs), ...
                num2str(rungs));
        end
    end
    % TWO-POINT MULTISTART, FIRST-SUCCESS-WINS (see GTO_tulip/catalog/
    % run_gto_catalog.m for the rationale). thrLock stays true on BOTH
    % passes; only the top-rung tf0 seed varies.
    %
    % Pass A -- warm tf0 seed:
    optsA.batchSec = max(0, batchSec - toc(tAll));
    % FAIL LOUD: no silent continue past a sheet whose engine call errors --
    % log which sheet+pass failed to stdout and the combined log, then
    % rethrow so the shell driver's non-zero-exit abort path fires.
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
    % Pass B -- cold tf0, SCOPED to the mop-up set: cells pass A itself just
    % attempted (ATT>=1) and could not clear (no OK rung), computed fresh
    % from the sheet file pass A just wrote. Empty set -> pass B skipped.
    if toc(tAll) <= batchSec
        mopCells = zeros(0,2);
        if isfile(outMat)
            S = load(outMat, 'ATT', 'OK');
            [mD, mA] = find(S.ATT >= 1 & ~any(S.OK,3));
            mopCells = [mD, mA];
        end
        if ~isempty(mopCells)
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
    % OPEN-CELL SWEEP: reload THIS sheet FRESH FROM DISK so the completion
    % marker is never based on in-memory state. No sheet on disk (zero cells
    % attempted) is trivially open.
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

if planOnly, return, end

if allDone && numel(sheetSel) == nSheets
    fid = fopen(logComb, 'a');
    fprintf(fid, 'CATALOG ALL SHEETS COMPLETE\n');  fclose(fid);
elseif allDone
    fid = fopen(logComb, 'a');
    fprintf(fid, 'SELECTED SHEETS COMPLETE (%s)\n', mat2str(sheetSel));
    fclose(fid);
end
end

% ------------------------------------------------------------------------
function delete_pool(pool)
% DELETE_POOL  Delete the fence pool this call opened, if still valid.
% INPUTS: pool parallel.Pool.  OUTPUTS: none.
if ~isempty(pool) && isvalid(pool)
    try delete(pool); catch, end
end
end

% ------------------------------------------------------------------------
function v = fdef(s, f, v0)
% FDEF  s.(f) if present else v0.  INPUTS: s; f; v0.  OUTPUTS: v.
if isfield(s, f), v = s.(f); else, v = v0; end
end
