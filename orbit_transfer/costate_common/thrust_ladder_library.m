function P = thrust_ladder_library(outMat, opts)
%% Purpose:
%
%   Build a costate library with a THRUST axis.
%
% For every phase pair on the (departure x arrival) torus, anchor a direct
% minimum-time solve at the HIGH-thrust end -- where the transfer is nearly
% impulsive, sub-revolution, and converges cold in seconds -- then walk the
% thrust DOWN through a rung set, each rung warm-started from the rung above.
% Continuation down the ladder keeps every cell on ONE solution family;
% cold-solving each (cell, thrust) independently does not (measured: cold
% 15 N solves on two meshes land 50% apart in t_f, both feasible).
%
% Each rung runs the full pipeline: direct solve -> continuous verification
% (fly the control) -> ms_tfmin multiple-shooting refinement -> acceptance by
% pumpkyn.cr3bp.tfMin. Results save after EVERY rung.
%
%% Inputs:
%
%   outMat - output .mat path (rewritten after every rung)
%   opts   - (optional) struct:
%              .rungs     thrust rung set, N, high to low
%                         [15 12 10 7 5 3 2 1.5 1] -- 1.5 added
%                         because t_f(T) steepens below 2 N, so linear
%                         interpolation needs a tighter rung there
%              .ispS      specific impulse, s              [1710]
%              .m0kg      initial mass, kg                 [150]
%              .N         collocation intervals            [400]
%              .floorKm   lunar altitude floor, km         [500]
%              .floorEarthKm Earth altitude floor, km; [] = none,
%                         passed to casadi_mintime_dro as minAltEarthKm
%                         on every rung (Earth departures)   [[]]
%              .gateKm    flown-arrival gate, km           [100]
%              .directGateKm  gate on the DIRECT control's flown miss
%                         only (the first gate, before ms_tfmin refines).
%                         The refined-arrival gate, the msSeg loop and
%                         meta.gateKm keep using .gateKm. A repair call
%                         may raise it: on the GTO -> DRO demo cell the
%                         direct control of the reverse-warm-started 15 N
%                         basin flies 193.6 km off at N = 400 and
%                         refinement closes it to 7e-5 km. Not stored in
%                         meta                            [.gateKm]
%              .cells     [k x 2] (iD,iA) subset           [all 12x12]
%              .nD,.nA    torus resolution                 [12 12]
%              .sD0,.sA0  phase-grid origins; defaults MATCH the
%                         low-thrust library grid            [0, 0.075378]
%              .maxIter   NLP iteration cap                [3000]
%              .logFile   append-mode log path             [stdout]
%              .thrLock   pin thr==1 in casadi_mintime_dro (classical
%                         min-time form; deletes the mass-floor coast
%                         manifold a free throttle can fall into on a
%                         sub-revolution transfer)                [false]
%              .tf0       cold t_f guess (ND) for a cell's TOP rung only;
%                         [] leaves casadi_mintime_dro's own internal
%                         default (4.0 ND) in effect, same as today
%                                                                    [[]]
%              .maxCpuSec IPOPT max_cpu_time per rung (s), passed straight
%                         through to casadi_mintime_dro               [300]
%              .seedTop   struct .X [7 x N+1] [r;v;m], .U [4 x N+1]
%                         [alpha;thr], .tf (ND): an EXTERNAL warm start
%                         for the TOP rung of every cell this call
%                         solves (a campaign repair tool's reverse warm
%                         start or phase continuation); [] = the cold
%                         tf0 path, byte-identical to before. An
%                         external seed may land on a basin whose direct
%                         control misses the flown gate before refinement
%                         (measured on the GTO -> DRO demo cell: 193.6 km
%                         at N = 400), so a repair call may raise
%                         .directGateKm and must then enforce its own
%                         gate on the recorded FLYKM              [[]]
%              .msSeg     multiple-shooting segment counts tried in
%                         order; the first refinement that converges
%                         within the flown gate ends the loop     [[12 24]]
%              .msWallS   ms_tfmin SOFT wall budget per attempt (s). Soft:
%                         checked between fsolve iterations only, so it
%                         cannot stop a single residual evaluation that
%                         never returns                             [120]
%              .msHardCapSec  HARD wall-clock cap (s) on each refinement
%                         attempt (one ms_tfmin call plus the flight of
%                         its result), enforced by running the attempt on
%                         a parallel worker under run_capped, which
%                         cancels the worker at the cap. [] or 0 = off:
%                         the attempt runs in-process exactly as before
%                         (byte-identical). WHY: pumpkyn.cr3bp.tfMinProp
%                         restarts ode45 in a while loop with no progress
%                         check; an iterate whose segment flow hits an
%                         ode45 step-size collapse restarts at the same
%                         time forever, inside ONE residual evaluation,
%                         where .msWallS never fires (GTO -> DRO pilot,
%                         2026-11-10: 55 min, 2.5 million warnings). A
%                         capped or worker-failed attempt is logged
%                         ('ms capped after %g s' / 'ms worker error
%                         after %g s') and treated as a non-converged
%                         refinement for that msSeg entry. A CAPPED
%                         attempt also ends the msSeg loop for that rung
%                         (the rung is recorded as not refined and the
%                         ladder goes on): on the pilot's hung cell both
%                         K = 12 and K = 24 stuck (an iterate's junction
%                         mass runs to zero mid-segment), so trying the
%                         next K costs another full cap for nothing. A
%                         fast worker ERROR still tries the next K. Pool: the open
%                         pool if any (kept), else capped_pool(1) opened
%                         at the first attempt and deleted when this call
%                         returns. No pool obtainable (no PCT, licence
%                         held elsewhere): logged once, attempts run
%                         UNFENCED (the current_pool contract). Workers
%                         inherit the client path at pool creation. With
%                         the option on, MATLAB:ode45:IntegrationTolNotMet
%                         is silenced on the worker inside the fenced call
%                         AND on the client for the whole call (restored
%                         on return), which covers the no-pool fallback
%                         and the unfenced tfMin acceptance          [[]]
%              (.msFcn    TEST SEAM, not for campaigns: the refinement
%                         function, same signature as ms_tfmin
%                                                           [@ms_tfmin])
%              .redo      true = re-solve every cell in .cells even if
%                         it already holds OK rungs: its ATT is reset to
%                         0 and its TF/FLYKM/ACCDZ/RES/WALL/OK/Z8 slices
%                         are cleared before the ladder runs. Refuses
%                         without .cells (thrust_ladder_library:redoCells)
%                         and unless a size-matched sheet is loaded
%                         (needs resume = true and an existing sheet of
%                         the call's size: thrust_ladder_library:redoSheet)
%                         The caller owns the backup of what it clears
%                                                                 [false]
%              .keepDirect  true = keep every rung's DIRECT solution in a
%                         sidecar <outMat minus .mat>_direct.mat (the sheet
%                         keeps only refined costates, so a question about
%                         the direct solutions -- e.g. the GTO -> DRO pilot's
%                         46 refined entries flying below an Earth floor the
%                         direct solution holds -- needs a re-solve without
%                         it). Recorded for every rung whose direct solve
%                         passes the engine's direct gate (success and
%                         defect < 1e-9), whatever happens after: pre-flight
%                         rejection, flown-gate stop, refinement or
%                         acceptance failure, cap, or OK. Variables:
%                           DX  [7 x N+1 x nD x nA x nR] node states
%                           DU  [3 x N+1 x nD x nA x nR] node thrust
%                               directions (rows 1:3 of the solver's U)
%                           DTF [nD x nA x nR] direct t_f (ND)
%                           DEARTHKM, DMOONKM [nD x nA x nR] minimum Earth
%                               (centre [-muStar 0 0], 6378.137 km) and Moon
%                               (centre [1-muStar 0 0], 1737.4 km) altitude
%                               over the nodes AND the Hermite-Simpson
%                               midpoints o.Xm when present, km
%                           DDONE [nD x nA x nR] logical, slice recorded
%                           rungs, sD, sA, N, sheetFile (sheet file NAME)
%                         NaN / false where not recorded. Saved (temp file +
%                         movefile) after the sheet at every sheet save and
%                         at the end of every cell; a cell's slices are
%                         cleared and saved BEFORE its attempt-counter
%                         save, and redo cells are cleared up front, so a
%                         kill never leaves the sidecar claiming a solve the
%                         sheet does not hold. On resume the sidecar is
%                         loaded only when the sheet was loaded and its
%                         sizes, rungs, N, sD, sA and sheetFile match, else a
%                         fresh one starts (one log line). A failed sidecar
%                         write is logged and the run goes on (the sheet is
%                         never lost to it). About 12 MB per 12 x 6 x 5
%                         sheet. false = no sidecar, no extra work, sheet
%                         byte-identical (tests/test_ladder_keepdirect)
%                                                                 [false]
%
%% Outputs:
%
%   P - struct with [nD x nA x nRung] arrays .TF (ND), .FLYKM, .ACCDZ,
%       .RES, .WALL, .OK, plus .Z8 [8 x nD x nA x nRung] refined costates,
%       .rungs, .sD, .sA, .meta (orbit + thruster definition).
%
%% References:
%   [1] casadi_mintime_dro.m -- the direct transcription.
%   [2] ms_tfmin.m -- multiple-shooting refinement.
%   [3] process/COSTATE_LIBRARY_PIPELINE.md -- the three-step process.
%
%% Revision History:
%  M. Casey                                                   (c) 08/05/2026
%  M. Casey  moved into costate_common (cleanup step 12)      (c) 09/16/2026
%  M. Casey  floorEarthKm option (default off)                (c) 11/09/2026
%  M. Casey  seedTop and redo options (default off)           (c) 11/10/2026
%  M. Casey  msHardCapSec hard fence on ms_tfmin (default off) (c) 11/10/2026
%  M. Casey  keepDirect direct-solution sidecar (default off)  (c) 11/11/2026
%  Copyright Coorbital Inc.
%% ------------------------ Begin Code Sequence ---------------------------

if nargin < 2, opts = struct(); end
d = @(f,v) fdef(opts, f, v);
rungs   = d('rungs', [15 12 10 7 5 3 2 1.5 1]);
ispS    = d('ispS', 1710);
m0kg    = d('m0kg', 150);
N       = d('N', 400);
floorKm = d('floorKm', 500);
floorEarthKm = d('floorEarthKm', []);   % Earth floor, km; [] = none (2026-11-09)
gateKm  = d('gateKm', 100);
directGateKm = d('directGateKm', gateKm);   % first (direct-control) gate only
nD      = d('nD', 12);
nA      = d('nA', 12);
% Phase-grid ORIGINS. These must match the low-thrust library exactly or the
% two libraries sample different phase pairs and cannot be merged: that grid
% is offset in arrival phase by the demo's max-velocity-angle anchor.
sD0     = d('sD0', 0);
sA0     = d('sA0', 0.075378);
maxIter = d('maxIter', 3000);
% Diagnostic-A pass-throughs (2026-08-26): the two decisive GTO-recipe
% knobs already existed inside casadi_mintime_dro but were unreachable
% from here. Defaults reproduce today's behavior byte-for-byte: thrLock
% off (free throttle, the historical formulation test), tf0Top empty
% (casadi_mintime_dro's own cold default, 4.0 ND -- the DRO-campaign
% number), maxCpuSec 300 (the value that used to be hardcoded at the
% call site below).
thrLock   = d('thrLock',   false);
tf0Top    = d('tf0',       []);
maxCpuSec = d('maxCpuSec', 300);
% maxCells: process at most this many fresh cells then return. Lets a shell
% driver run the ladder in small batches under an OS-level timeout, so a
% hang inside an uninterruptible solver call costs minutes, not the run.
maxCells = d('maxCells', inf);
% batchSec: stop cleanly BETWEEN cells once this much wall time has passed.
% The driver's OS timeout is then only a backstop for genuine hangs, so a
% cell is never killed merely for being last in a batch (which would count an
% unfair failed attempt against it).
batchSec = d('batchSec', inf);
logFile = d('logFile', '');
lg = @(varargin) logmsg(logFile, sprintf(varargin{:}));

% --- constants and orbits -----------------------------------------------
% Orbit and solver settings are OPTIONS, not literals, so the front door
% (run_costate_library.m) genuinely controls them.
muStar = d('muStar', 0.012150585609624);
lStar  = d('lStar',  389703.264829278);
tStar  = d('tStar',  382981.289129055);
tauDRO = d('tauDRO', 1.0);
NpT    = d('NpTulip', 7);
tauT   = d('tauTulip', 5*2*pi/6);
pmT    = d('pmTulip', -1);
% FAMILY-AGNOSTIC ENDPOINTS: any pumpkyn family via the shared provider
% (costate_common/get_family_orbit). Defaults reproduce the DRO->tulip
% campaign exactly; a Halo campaign passes depFamily='halo' etc.
depFamily = d('depFamily', 'dro');
depParams = d('depParams', struct('tau', tauDRO));
arrFamily = d('arrFamily', 'tulip');
arrParams = d('arrParams', struct('Np', NpT, 'pm', pmT));
msSeg  = d('msSeg', [12 24]);      % multiple-shooting segment ladder
msWallS= d('msWallS', 120);        % ms_tfmin budget per attempt (s)
accTol = d('accTol', 1e-6);        % tfMin acceptance tolerance
msHardCapSec = d('msHardCapSec', []);   % hard fence per refinement (s); []/0 = off
msFcn  = d('msFcn', @ms_tfmin);    % test seam (header)
if ~isempty(msHardCapSec) && ~(isscalar(msHardCapSec) && isnumeric(msHardCapSec) && ...
        isreal(msHardCapSec) && isfinite(msHardCapSec) && msHardCapSec >= 0)
    error('thrust_ladder_library:msHardCapSec', ...
          'opts.msHardCapSec must be [] or a finite non-negative scalar (s)');
end
if ~isa(msFcn, 'function_handle')
    error('thrust_ladder_library:msFcn', 'opts.msFcn must be a function handle');
end
msCapOn = ~isempty(msHardCapSec) && msHardCapSec > 0;
capPool = [];  capNoPool = false;  capCleanup = [];   % fence pool state (msCapOn only)
nWorkerErr = 0;                    % CONSECUTIVE fenced worker errors (2 = fenceBroken)
if msCapOn
    % a stuck flight warns once per ode45 restart (2.5 million lines on the
    % pilot); option on only, so the default path's output is untouched
    odeWarn = warning('off', 'MATLAB:ode45:IntegrationTolNotMet');
    restoreOdeWarn = onCleanup(@() warning(odeWarn));
end
% THE shared propulsion conversion (costate_common/nd_propulsion)
ndp = nd_propulsion([], ispS, m0kg, lStar, tStar);
cnd = ndp.cnd;   ndT = ndp.ndT;

[tD, rvD] = get_family_orbit(depFamily, depParams);
[tT, rvT] = get_family_orbit(arrFamily, arrParams);
% THE shared endpoint rule (costate_common/phase_state, FINDINGS 44)
depAt = phase_state(tD, rvD);   arrAt = phase_state(tT, rvT);

sD = mod(sD0 + (0:nD-1)/nD, 1);  sA = mod(sA0 + (0:nA-1)/nA, 1);
todo = d('cells', []);
redo = d('redo', false);
if redo && isempty(todo)
    error('thrust_ladder_library:redoCells', ...
          'opts.redo requires opts.cells (the cells to re-solve)');
end
if isempty(todo)
    [gA, gD] = meshgrid(1:nA, 1:nD);
    todo = [gD(:), gA(:)];
end
seedTop = d('seedTop', []);
if ~isempty(seedTop)
    okSeed = isstruct(seedTop) && all(isfield(seedTop, {'X', 'U', 'tf'})) && ...
             isequal(size(seedTop.X), [7 N+1]) && isequal(size(seedTop.U), [4 N+1]) && ...
             isscalar(seedTop.tf) && seedTop.tf > 0 && ...
             all(isfinite(seedTop.X(:))) && all(isfinite(seedTop.U(:))) && isfinite(seedTop.tf);
    if ~okSeed
        error('thrust_ladder_library:seedTop', ...
              'opts.seedTop must carry .X [7 x %d], .U [4 x %d] and a positive .tf', N+1, N+1);
    end
end

keepDirect = d('keepDirect', false);
if ~(isscalar(keepDirect) && (islogical(keepDirect) || isnumeric(keepDirect)) && ...
        any(keepDirect == [0 1]))
    error('thrust_ladder_library:keepDirect', 'opts.keepDirect must be a logical scalar');
end
keepDirect = logical(keepDirect);
dsDirty = false;                       % sidecar holds a record not yet on disk

nR = numel(rungs);
TF = nan(nD,nA,nR); FLYKM = TF; ACCDZ = TF; RES = TF; WALL = TF;
OK = false(nD,nA,nR);  Z8 = nan(8,nD,nA,nR);
% ATT counts attempts per cell, written BEFORE the cell is solved so that a
% cell which HANGS (and is killed externally) still leaves a record. Without
% it a hung cell has no verified rung, resume would not skip it, and a
% batched driver would retry it forever.
ATT = zeros(nD,nA);
sheetLoaded = false;
maxAtt = d('maxAtt', 2);
% RESUME: a prior run's products are data -- reload them and skip cells that
% already produced verified rungs (opts.resume, default true). A redo call
% keeps every cell it names, resets its attempts and clears its slices.
if d('resume', true) && isfile(outMat)
    Q = load(outMat);
    if isequal(size(Q.OK), [nD nA nR])
        TF=Q.TF; FLYKM=Q.FLYKM; ACCDZ=Q.ACCDZ; RES=Q.RES; WALL=Q.WALL;
        OK=Q.OK; Z8=Q.Z8;
        if isfield(Q,'ATT'), ATT = Q.ATT; end
        sheetLoaded = true;
        if redo
            for kq = 1:size(todo,1)
                iq = todo(kq,1);  jq = todo(kq,2);
                ATT(iq,jq) = 0;
                TF(iq,jq,:) = NaN;  FLYKM(iq,jq,:) = NaN;  ACCDZ(iq,jq,:) = NaN;
                RES(iq,jq,:) = NaN; WALL(iq,jq,:) = NaN;   OK(iq,jq,:) = false;
                Z8(:,iq,jq,:) = NaN;
            end
            nSkip = 0;
        else
            keep = false(size(todo,1),1);
            for kq = 1:size(todo,1)
                iq = todo(kq,1);  jq = todo(kq,2);
                keep(kq) = ~any(OK(iq,jq,:)) && ATT(iq,jq) < maxAtt;
            end
            nSkip = nnz(~keep);  todo = todo(keep,:);
        end
    else
        nSkip = 0;
    end
else
    nSkip = 0;
end
if redo && ~sheetLoaded
    error('thrust_ladder_library:redoSheet', ...
          'opts.redo needs resume = true and an existing sheet of the same size');
end
% DIRECT SIDECAR (opts.keepDirect): loaded only alongside a loaded sheet,
% else fresh; redo cells cleared up front like their sheet slices
if keepDirect
    [~, sheetName, sheetExt] = fileparts(outMat);
    dirMat = [outMat(1:end-4) '_direct.mat'];
    [DS, dsReplaced] = direct_sidecar_load(dirMat, sheetLoaded, [sheetName sheetExt], ...
                                           rungs, sD, sA, N, nD, nA, lg);
    if dsReplaced || redo
        for kq = 1:size(todo,1)
            DS = direct_clear_cell(DS, todo(kq,1), todo(kq,2));
        end
        direct_sidecar_save(dirMat, DS, lg);
    end
end
meta = struct('muStar',muStar,'lStar',lStar,'tStar',tStar, ...
    'tauDRO',tauDRO,'NpTulip',NpT,'tauTulip',tauT,'pmTulip',pmT, ...
    'depFamily',depFamily,'depParams',depParams, ...
    'arrFamily',arrFamily,'arrParams',arrParams, ...
    'periodDRO',tD(end),'periodTulip',tT(end), ...
    'tauDep',tD(end),'tauArr',tT(end), ...     % schema-v2 canonical names
    'sheetSchema',2, 'env',env_pin(), ...
    'ispS',ispS,'m0kg',m0kg,'N',N,'floorKm',floorKm, ...
    'floorEarthKm',floorEarthKm,'gateKm',gateKm, ...
    'sD0',sD0,'sA0',sA0);

lg('THRUST LADDER: %d cells x %d rungs [%s] N, Isp=%d s, m0=%d kg, N=%d (%d cells resumed/skipped)', ...
   size(todo,1), nR, num2str(rungs), ispS, m0kg, N, nSkip);
if ~isempty(seedTop), lg('  top rung seeded externally (opts.seedTop, tf %.5f)', seedTop.tf); end
if redo, lg('  redo: %d cell(s) cleared and re-solved', size(todo,1)); end
tAll = tic;
for kc = 1:min(size(todo,1), maxCells)
    if toc(tAll) > batchSec
        lg('  [batch wall budget reached -- exiting cleanly after %d cells]', kc-1);
        break
    end
    iD = todo(kc,1);  iA = todo(kc,2);
    rv0 = depAt(sD(iD)).';          % row, as before
    rvf = arrAt(sA(iA)).';
    seedX = [];  seedU = [];  seedTf = tf0Top;  Tprev = [];
    if ~isempty(seedTop)
        seedX = seedTop.X;  seedU = seedTop.U;  seedTf = seedTop.tf;
    end
    ATT(iD,iA) = ATT(iD,iA) + 1;      % record BEFORE solving (hang-proof)
    if keepDirect                     % sidecar cleared BEFORE the attempt is recorded
        DS = direct_clear_cell(DS, iD, iA);
        direct_sidecar_save(dirMat, DS, lg);  dsDirty = false;
    end
    save(outMat, 'TF','FLYKM','ACCDZ','RES','WALL','OK','Z8','ATT', ...
         'rungs','sD','sA','meta');
    for kr = 1:nR
        TN = rungs(kr);  Tnd = ndT(TN);
        t0 = tic;
        if ~isempty(Tprev)
            % time-of-flight scales roughly as T^-0.6 in this band (pilot)
            seedTf = seedTf * (Tprev/TN)^0.6;
        end
        try
            o = casadi_mintime_dro(rv0(1:6), rvf(1:6), Tnd, cnd, muStar, N, ...
                    seedX, seedU, seedTf, ...
                    struct('maxIter',maxIter,'scheme','hermite-simpson', ...
                           'sundman',false,'returnModel',true, ...
                           'minAltKm',floorKm,'minAltEarthKm',floorEarthKm, ...
                           'maxCpuSec',maxCpuSec, ...
                           'thrLock',thrLock));
        catch ME
            lg('  (%2d,%2d) T=%5.1f  direct ERROR: %s', iD, iA, TN, ME.message);
            break
        end
        if ~(o.success && o.maxDefect < 1e-9)
            lg('  (%2d,%2d) T=%5.1f  direct fail (defect %.1e) -- ladder stops', ...
               iD, iA, TN, o.maxDefect);
            break
        end
        if keepDirect                 % every direct success, whatever follows
            try                       % a diagnostic must never cost a cell
                DS = direct_record(DS, o, iD, iA, kr, muStar, lStar);  dsDirty = true;
            catch MEds
                lg('  (%2d,%2d) T=%5.1f  direct sidecar record FAILED: %s', iD, iA, TN, MEds.message);
            end
        end
        % Cheap sanity pre-check before any integrator touches the solve
        % (preflight_screen, this folder).
        [pfOK, ~, pfAltKm] = preflight_screen(o, muStar, lStar, floorKm, seedTf);
        if ~pfOK
            lg('  (%2d,%2d) T=%5.1f  rejected pre-flight (minAlt %.0f km, tf %.4f) -- ladder stops', ...
               iD, iA, TN, pfAltKm, o.tf);
            break
        end
        C = certify_dro_mintime(o, struct('muStar',muStar,'lStar',lStar,'tStar',tStar), ...
                Tnd, cnd, struct('tfRef',[],'verbose',false,'posTolKm',inf));
        flyKm = C.globKm;
        TF(iD,iA,kr) = o.tf;  FLYKM(iD,iA,kr) = flyKm;  WALL(iD,iA,kr) = toc(t0);
        if flyKm > directGateKm
            lg('  (%2d,%2d) T=%5.1f  tf=%.5f flown %.1f km > gate -- ladder stops', ...
               iD, iA, TN, o.tf, flyKm);
            break
        end
        % --- indirect refinement + acceptance -----------------------------
        accDz = NaN;  okCell = false;
        if ~isempty(o.lamDef)
            % Covector mapping + seed assembly live in ONE home now
            % (costate_common: duals_to_costates via harvest_ms_seed);
            % the sign vote and the Hermite-Simpson midpoint association
            % are enforced there (harvest_ms_seed, this folder).
            for K = msSeg
                seed = harvest_ms_seed(o, K);
                if msCapOn
                    % HARD FENCE (opts.msHardCapSec): refinement + flight on
                    % a worker. Acquire (or re-acquire, after an idle or
                    % deleted pool) lazily, so a call that never refines
                    % never pays for a pool.
                    if ~capNoPool && (isempty(capPool) || ~isvalid(capPool))
                        capPool = current_pool();
                        if isempty(capPool) || ~isvalid(capPool)
                            capPool = capped_pool(1);
                            if ~isempty(capPool)
                                capCleanup = onCleanup(@() delete_pool(capPool));
                            end
                        end
                        if ~isempty(capPool)
                            lg('  ms hard cap %g s ACTIVE (pool of %d worker(s))', ...
                               msHardCapSec, capPool.NumWorkers);
                        else
                            capNoPool = true;
                            lg(['  ms hard cap %g s requested but no parallel pool could be ' ...
                                'obtained: ms_tfmin runs UNFENCED for the rest of this call'], msHardCapSec);
                        end
                    end
                end
                if msCapOn && ~capNoPool
                    [stMs, outMs, capPool] = fenced_call(capPool, msHardCapSec, ...
                        @ms_refine_fenced, 3, msFcn, rv0(1:6), rvf(1:6), seed, Tnd, cnd, ...
                        muStar, msWallS, lStar);
                    [nWorkerErr, stopK] = fence_outcome(stMs, nWorkerErr, lg, ...
                        sprintf('(%2d,%2d) T=%5.1f  K=%d  ms', iD, iA, TN, K));
                    if stopK, break, end
                    if ~strcmp(stMs.status, 'ok'), continue, end
                    [z, info, msKm] = outMs{:};
                else
                    [z, info] = msFcn(rv0(1:6), rvf(1:6), seed, Tnd, cnd, muStar, ...
                            struct('wallSec', msWallS));
                    [~, rvFly] = pumpkyn.cr3bp.tfMinProp(z(8), ...
                            [rv0(1:6)'; 1; z(1:7)], Tnd, cnd, muStar);
                    msKm = norm(rvFly(end,1:3) - rvf(1:3))*lStar;
                end
                if info.converged && msKm < gateKm
                    if msCapOn && ~capNoPool
                        % the acceptance flight runs tfMinProp too: same fence
                        [stA, outA, capPool] = fenced_call(capPool, msHardCapSec, ...
                            @tfmin_accept_fenced, 1, rv0(1:6), rvf(1:6), z, Tnd, cnd, muStar);
                        nWorkerErr = fence_outcome(stA, nWorkerErr, lg, ...
                            sprintf('(%2d,%2d) T=%5.1f  K=%d  tfMin acceptance', iD, iA, TN, K));
                        if ~strcmp(stA.status, 'ok'), break, end   % not accepted
                        zA = outA{1};
                    else
                        evalc('zA = pumpkyn.cr3bp.tfMin(rv0(1:6), rvf(1:6), z, Tnd, cnd, muStar);');
                    end
                    accDz = norm(zA - z);
                    if accDz < accTol
                        okCell = true;  Z8(:,iD,iA,kr) = z;
                        RES(iD,iA,kr) = info.normR;  TF(iD,iA,kr) = z(8);
                        FLYKM(iD,iA,kr) = msKm;
                    end
                    break
                end
            end
        end
        ACCDZ(iD,iA,kr) = accDz;  OK(iD,iA,kr) = okCell;
        WALL(iD,iA,kr) = toc(t0);
        lg('  (%2d,%2d) T=%5.1f  tf=%.5f (%.3f d)  fly=%.2fkm  acc=%.1e  %s (%.0fs)', ...
           iD, iA, TN, TF(iD,iA,kr), TF(iD,iA,kr)*tStar/86400, ...
           FLYKM(iD,iA,kr), accDz, okstr(okCell), toc(t0));
        % next rung warm-starts from THIS rung's direct solution. With the
        % uniform mesh the node grid already matches, so no resampling.
        seedX = o.X;
        seedU = o.U;
        seedU(1:3,:) = seedU(1:3,:) ./ max(vecnorm(seedU(1:3,:),2,1), eps);
        seedU(4,:)   = min(max(seedU(4,:),0),1);
        seedTf = o.tf;  Tprev = TN;
        save(outMat, 'TF','FLYKM','ACCDZ','RES','WALL','OK','Z8','ATT', ...
             'rungs','sD','sA','meta');                        % EVERY rung
        if keepDirect                 % after the sheet
            direct_sidecar_save(dirMat, DS, lg);  dsDirty = false;
        end
    end
    if keepDirect && dsDirty          % a rung recorded, then the ladder broke out
        direct_sidecar_save(dirMat, DS, lg);  dsDirty = false;
    end
    if mod(kc,10) == 0
        lg('--- %d/%d cells, %d (cell,rung) entries verified, %.1f min elapsed', ...
           kc, size(todo,1), nnz(OK), toc(tAll)/60);
    end
end
P = struct('TF',TF,'FLYKM',FLYKM,'ACCDZ',ACCDZ,'RES',RES,'WALL',WALL, ...
           'OK',OK,'Z8',Z8,'ATT',ATT,'rungs',rungs,'sD',sD,'sA',sA,'meta',meta);
lg('THRUST LADDER DONE: %d verified entries over %d cells, %.1f min', ...
   nnz(OK), size(todo,1), toc(tAll)/60);
end

% ------------------------------------------------------------------------
function [st, out, pool] = fenced_call(pool, capSec, fcn, nout, varargin)
% FENCED_CALL  One call under run_capped. fcn catches its own errors on
% the worker and returns them as two trailing outputs (errId, errMsg); a
% failure outside fcn (unresolvable on the worker) arrives through
% run_capped's opt-in cause output. Either way the CAUSE is logged.
% INPUTS:  pool parallel.Pool; capSec double; fcn fhandle returning
%          nout + 2 outputs; nout int; varargin fcn's arguments.
% OUTPUTS: st struct .status ('ok' | 'capped' | 'error' | 'poolLost'),
%          .sec wall time, .msg cause text; out {1 x nout} fcn's outputs
%          when ok; pool (emptied when run_capped deleted it).
out = cell(1, nout);
t0 = tic;
try
    res = cell(1, nout + 3);       % fcn's nout + errId + errMsg, then run_capped's cause
    [okRun, res{:}] = run_capped(pool, fcn, nout + 2, capSec, varargin{:});
catch ME
    % run_capped:poolUnusable -- the pool was deleted; the next fenced
    % call opens a fresh one
    pool = [];
    st = struct('status', 'poolLost', 'sec', toc(t0), 'msg', [ME.identifier ': ' ME.message]);
    return
end
sec = toc(t0);
if ~okRun
    if sec >= capSec
        st = struct('status', 'capped', 'sec', sec, 'msg', '');
    else
        % fcn catches everything inside itself, so this is a failure
        % OUTSIDE it -- typically the worker cannot resolve the function
        % (path): run_capped's opt-in cause carries fut.Error
        st = struct('status', 'error', 'sec', sec, 'msg', res{nout + 3});
    end
elseif ~isempty(res{nout + 1})
    st = struct('status', 'error', 'sec', sec, 'msg', [res{nout + 1} ': ' res{nout + 2}]);
else
    st = struct('status', 'ok', 'sec', sec, 'msg', '');
    out = res(1:nout);
end
end

function [nErr, stopK] = fence_outcome(st, nErr, lg, what)
% FENCE_OUTCOME  Log a fenced call's outcome and keep the CONSECUTIVE
% worker-error count. A capped call or a lost pool stops further K on the
% rung; a second consecutive worker error is infrastructure, not a
% numerical wall (CAMPAIGN_DISCIPLINE section 12), and throws
% thrust_ladder_library:fenceBroken.
% INPUTS:  st struct (fenced_call); nErr int; lg logger handle; what char
%          (cell / rung / K prefix).
% OUTPUTS: nErr updated count; stopK logical, end the msSeg loop.
stopK = false;
switch st.status
    case 'ok'
        nErr = 0;
    case 'capped'
        nErr = 0;  stopK = true;
        lg('  %s capped after %g s -- no further K on this rung', what, round(st.sec));
    case 'poolLost'
        nErr = 0;  stopK = true;
        lg('  %s fence pool lost after %g s (%s) -- no further K on this rung; next call opens a fresh pool', ...
           what, round(st.sec), st.msg);
    otherwise
        nErr = nErr + 1;
        lg('  %s worker error after %g s: %s', what, round(st.sec), st.msg);
        if nErr >= 2
            error('thrust_ladder_library:fenceBroken', ...
                  ['two consecutive fenced calls failed on the worker (last: %s) -- infrastructure, ' ...
                   'not a numerical wall; check the pool workers'' path (pumpkyn, costate_common, oclib)'], st.msg);
        end
end
end

function [z, info, msKm, errId, errMsg] = ms_refine_fenced(msFcn, rv0, rvf, seed, Tnd, cnd, muStar, msWallS, lStar)
% MS_REFINE_FENCED  One refinement attempt as run on a run_capped worker:
% the same ms_tfmin call and flight of its result as the in-process path,
% with the ode45 tolerance warning silenced on the worker (a stuck
% tfMinProp emits one per restart into the future's diary). Errors are
% caught here and returned, so the client logs their cause.
% INPUTS:  msFcn fhandle (ms_tfmin signature); rv0, rvf [1 x 6]; seed
%          struct (harvest_ms_seed); Tnd, cnd, muStar, msWallS, lStar double.
% OUTPUTS: z [8 x 1]; info struct (ms_tfmin); msKm double, flown arrival
%          miss of z (km); errId, errMsg char ('' when no error).
ws = warning('off', 'MATLAB:ode45:IntegrationTolNotMet');
restoreWarn = onCleanup(@() warning(ws));
z = [];  info = [];  msKm = NaN;  errId = '';  errMsg = '';
try
    [z, info] = msFcn(rv0, rvf, seed, Tnd, cnd, muStar, struct('wallSec', msWallS));
    [~, rvFly] = pumpkyn.cr3bp.tfMinProp(z(8), [rv0(:); 1; z(1:7)], Tnd, cnd, muStar);
    dr = rvFly(end,1:3) - rvf(1:3);
    msKm = sqrt(sum(dr.^2))*lStar;
catch ME
    errId = ME.identifier;  errMsg = ME.message;
    if isempty(errId), errId = 'unidentified'; end
end
clear restoreWarn
end

function [zA, errId, errMsg] = tfmin_accept_fenced(rv0, rvf, z, Tnd, cnd, muStar)
% TFMIN_ACCEPT_FENCED  The tfMin acceptance solve as run on a run_capped
% worker (same call as the in-process path), ode45 warning silenced,
% errors returned rather than thrown.
% INPUTS:  rv0, rvf [1 x 6]; z [8 x 1]; Tnd, cnd, muStar double.
% OUTPUTS: zA [8 x 1]; errId, errMsg char ('' when no error).
ws = warning('off', 'MATLAB:ode45:IntegrationTolNotMet');
restoreWarn = onCleanup(@() warning(ws));
zA = [];  errId = '';  errMsg = '';
try
    evalc('zA = pumpkyn.cr3bp.tfMin(rv0, rvf, z, Tnd, cnd, muStar);');
catch ME
    errId = ME.identifier;  errMsg = ME.message;
    if isempty(errId), errId = 'unidentified'; end
end
clear restoreWarn
end

function delete_pool(pool)
% DELETE_POOL  Delete a fence pool this call opened, if still valid.
% INPUTS: pool parallel.Pool or [].  OUTPUTS: none.
if ~isempty(pool) && isvalid(pool)
    try delete(pool); catch, end
end
end

% ------------------------------------------------------------------------
function [DS, replaced] = direct_sidecar_load(dirMat, sheetLoaded, sheetFile, rungs, sD, sA, N, nD, nA, lg)
% DIRECT_SIDECAR_LOAD  The keepDirect sidecar: the one on disk when the
% sheet was loaded and the sidecar matches the call, else a fresh one (a
% line is logged when a file on disk is set aside; the caller then writes
% the fresh one over it at once, so no stale record outlives the call).
% INPUTS:  dirMat char; sheetLoaded logical; sheetFile char (name only);
%          rungs [1 x nR]; sD [1 x nD]; sA [1 x nA]; N, nD, nA int; lg
%          logger handle.
% OUTPUTS: DS struct (sidecar variables, see the header's .keepDirect);
%          replaced logical, a file on disk was set aside.
replaced = false;
nR = numel(rungs);
DS = struct('DX', nan(7, N+1, nD, nA, nR), 'DU', nan(3, N+1, nD, nA, nR), ...
            'DTF', nan(nD, nA, nR), 'DEARTHKM', nan(nD, nA, nR), ...
            'DMOONKM', nan(nD, nA, nR), 'DDONE', false(nD, nA, nR), ...
            'rungs', rungs, 'sD', sD, 'sA', sA, 'N', N, 'sheetFile', sheetFile);
if ~isfile(dirMat), return, end
why = '';
if ~sheetLoaded
    why = 'no sheet was resumed';
else
    try
        Q = load(dirMat);
        fOK = all(isfield(Q, fieldnames(DS)));
        if ~fOK
            why = 'variables missing';
        elseif ~(isequal(size(Q.DX, 1:5), [7 N+1 nD nA nR]) && ...
                 isequal(size(Q.DU, 1:5), [3 N+1 nD nA nR]) && ndims(Q.DX) <= 5 && ...
                 isequal(size(Q.DTF, 1:3), [nD nA nR]) && isequal(size(Q.DEARTHKM, 1:3), [nD nA nR]) && ...
                 isequal(size(Q.DMOONKM, 1:3), [nD nA nR]) && isequal(size(Q.DDONE, 1:3), [nD nA nR]) && ...
                 ndims(Q.DTF) <= 3)
            why = 'sizes differ';
        elseif ~(isequal(Q.rungs(:).', rungs(:).') && isequal(Q.N, N) && ...
                 isequal(Q.sD(:).', sD(:).') && isequal(Q.sA(:).', sA(:).') && ...
                 strcmp(Q.sheetFile, sheetFile))
            why = 'rungs / N / grid / sheetFile differ';
        end
    catch ME
        why = ['unreadable: ' ME.message];
    end
end
if isempty(why)
    try
        Q.DDONE = logical(Q.DDONE);
        DS = Q;
    catch ME
        why = ['DDONE not logical: ' ME.message];
    end
end
if ~isempty(why)
    lg('  direct sidecar %s does not match this call (%s): starting a fresh one', dirMat, why);
    replaced = true;
    if isfile(dirMat)                 % keep the old data: set it aside, never overwrite
        aside = [dirMat(1:end-4) '_setaside_' char(datetime('now', 'Format', 'yyyyMMdd_HHmmss')) '.mat'];
        try
            movefile(dirMat, aside);
            lg('  old direct sidecar set aside as %s', aside);
        catch ME
            lg('  old direct sidecar could not be set aside (%s)', ME.message);
        end
    end
end
end

function DS = direct_clear_cell(DS, iD, iA)
% DIRECT_CLEAR_CELL  Mark one cell's sidecar slices not recorded.
% INPUTS: DS struct; iD, iA int.  OUTPUTS: DS.
DS.DX(:,:,iD,iA,:) = NaN;   DS.DU(:,:,iD,iA,:) = NaN;
DS.DTF(iD,iA,:) = NaN;      DS.DEARTHKM(iD,iA,:) = NaN;
DS.DMOONKM(iD,iA,:) = NaN;  DS.DDONE(iD,iA,:) = false;
end

function DS = direct_record(DS, o, iD, iA, kr, muStar, lStar)
% DIRECT_RECORD  Record one rung's direct solution into the sidecar.
% Altitudes are minima over the nodes and, when the solver returns them,
% the Hermite-Simpson midpoints o.Xm.
% INPUTS:  DS struct; o casadi_mintime_dro output (.X [7 x N+1], .U
%          [4 x N+1], .Xm [7 x N] or [], .tf); iD, iA, kr int; muStar,
%          lStar double.
% OUTPUTS: DS.
R = o.X(1:3,:);
if isfield(o, 'Xm') && ~isempty(o.Xm), R = [R, o.Xm(1:3,:)]; end
dE = sqrt(sum((R - [-muStar; 0; 0]).^2, 1));
dM = sqrt(sum((R - [1-muStar; 0; 0]).^2, 1));
DS.DX(:,:,iD,iA,kr) = o.X;
DS.DU(:,:,iD,iA,kr) = o.U(1:3,:);
DS.DTF(iD,iA,kr) = o.tf;
DS.DEARTHKM(iD,iA,kr) = min(dE)*lStar - 6378.137;
DS.DMOONKM(iD,iA,kr)  = min(dM)*lStar - 1737.4;
DS.DDONE(iD,iA,kr) = true;
end

function direct_sidecar_save(dirMat, DS, lg)
% DIRECT_SIDECAR_SAVE  Write the sidecar via a temp file in the same folder
% + movefile; a failure is logged, never thrown (the sheet must survive).
% INPUTS: dirMat char; DS struct; lg logger handle.  OUTPUTS: none.
fold = fileparts(dirMat);
if isempty(fold), fold = pwd; end
tmp = [tempname(fold) '.mat'];
try
    save(tmp, '-struct', 'DS');
    movefile(tmp, dirMat, 'f');
catch ME
    lg('  direct sidecar write FAILED (%s): %s -- the sheet is unaffected', dirMat, ME.message);
    if isfile(tmp)
        try delete(tmp); catch, end
    end
end
end

% ------------------------------------------------------------------------
function v = fdef(s, f, v0)
% FDEF  s.(f) if present else v0.  INPUTS: s;f;v0.  OUTPUTS: v.
if isfield(s, f), v = s.(f); else, v = v0; end
end

function logmsg(f, s)
% LOGMSG  Append a line to the log file (or stdout if none).
% INPUTS: f path ('' = stdout); s text.  OUTPUTS: none.
if isempty(f)
    fprintf('%s\n', s);
else
    fid = fopen(f, 'a');  fprintf(fid, '%s\n', s);  fclose(fid);
end
end

function s = okstr(t)
% OKSTR  OK/fail marker.  INPUTS: t logical.  OUTPUTS: s [char].
if t, s = 'OK'; else, s = 'fail'; end
end

% ------------------------------------------------------------------------
function e = env_pin()
% ENV_PIN  Environment pinning for sheet metadata (schema v2): MATLAB
% release plus the pumpkyn/pumpkynPie git revisions, best-effort.
% INPUTS: none.  OUTPUTS: e struct (.matlab .pumpkynRev .pumpkynPieRev).
e = struct('matlab', version, 'pumpkynRev', 'unknown', ...
           'pumpkynPieRev', 'unknown');
% git walks up to the nearest enclosing repo on its own -- no fragile
% parent-counting (both external reviewers flagged the counted version,
% and disagreed with each other on the right count)
try
    p = which('pumpkyn.cr3bp.prop');
    if ~isempty(p)
        [st, rev] = system(['git -C "' fileparts(p) '" rev-parse --short HEAD']);
        if st == 0, e.pumpkynRev = strtrim(rev); end
    end
    p = which('pumpkynPie.cr3bp.getDRO');
    if ~isempty(p)
        [st, rev] = system(['git -C "' fileparts(p) '" rev-parse --short HEAD']);
        if st == 0, e.pumpkynPieRev = strtrim(rev); end
    end
catch
end
end
