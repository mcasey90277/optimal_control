function T = pair_transfer(dep, arr, sD, sA, thrustN, opts)
%% Purpose:
%
%   ONE MINIMUM-TIME LOW-THRUST TRANSFER FOR ANY ORBIT PAIR, solved and
%   gated: run the shared ladder engine (thrust_ladder_library) on one
%   (departure phase x arrival phase) cell from the top rung DOWN to
%   thrustN, re-fly the accepted entry independently (fly_transfer), and
%   attach every verdict a consumer needs -- flown miss, tfMin acceptance,
%   flight admissibility, LOCATED minimum Moon and Earth altitude, the
%   conjugate test and the hypothesis gates. Returns the flight for
%   plotting and the entry in library format, always with its junction
%   states. Generalizes GTO_DRO/lib/gto_dro_entry (left as it is) to any
%   pair get_family_orbit knows; puts no campaign folder on the path.
%
%   LADDER. Default rungs: the standard set [15 12 10 7 5 3 2 1.5 1 0.75
%   0.5] entries above thrustN, then thrustN (unique, descending); thrustN
%   above 15 N gives [thrustN] alone. The higher rungs are warm-start
%   stepping stones; the entry AT thrustN is returned. A given opts.rungs is
%   used as is and must end at thrustN.
%
%   GATES. gates.pass requires every BLOCKING gate:
%     ok          the engine accepted the entry at thrustN
%     flown       re-flight miss flownKm < gateKm
%     accDz       tfMin acceptance |dz| accDz < accTol
%     admissible  validate_flight on the re-flight (reached t_f, mass law,
%                 sample clearances)
%     moon        located minimum Moon altitude >= floorKm, every local
%                 minimum located (located_min_altitude, Moon radius
%                 1737.4 km, propagated on pumpkyn's tfMinEoM)
%     earth       the same over the Earth (6378.137 km) against
%                 floorEarthKm, only when floorEarthKm is set
%   gates.failed lists the names of the failed ones. The engine's floorKm
%   constrains the DIRECT solve at its nodes; the gate uses the located
%   minimum of the re-flown refined trajectory. The second-order verdicts
%   (conjPass, minLamV, minQmt, dimS, secondStatus) are REPORTED, not
%   blocking: a re-solve ms_tfmin(seed_from_z8(z8, rv0, K), tolR 3e-11,
%   wallSec 600, conjTest) counts as reproduced when converged and
%   max|z - z8| < 1e-6; secondStatus 'pass' | 'refuted' | 'notReproduced'
%   | 'error' | 'skipped' as in leg1_gates, plus 'capped'.
%
%   SECOND-ORDER FENCE. With msHardCapSec > 0 the second-order block runs
%   on a parallel worker under run_capped with that hard cap (the same
%   fence the engine puts around ms_tfmin). A timeout gives conjPass NaN,
%   conjReason 'capped after <x> s', secondStatus 'capped'; a worker error
%   gives 'error: <cause>', secondStatus 'error'. msHardCapSec [] or 0, or
%   no obtainable pool, runs it in-process.
%
%   REASONS (T.reason, '' when ok):
%     'top rung failed'         the ladder broke at the first rung
%     'ladder stopped at <x> N' the ladder broke at a later rung
%     'not accepted'            the last rung was reached, not accepted
%     'engine error: <msg>'     thrust_ladder_library threw (caught here)
%     'mass exhausted'          the accepted t_f leaves no mass (guard)
%     'planOnly'                opts.planOnly: nothing was solved
%
%   REGATE (no solve). Re-apply the pass logic to a stored result with
%   other thresholds:
%
%     T2 = pair_transfer('regate', T, th)
%
%   th may carry .gateKm .accTol .floorKm .floorEarthKm (others refuse with
%   pair_transfer:regate); the measured values in T.gates are used as
%   stored (a floorKm change does not re-solve or re-locate anything).
%   T2.thresholds, T2.gates.pass / .failed and T2.entry.verdicts are
%   updated, and one record (.thresholds .overrides .at datetime) is
%   appended to the cell array T2.entry.provenance.regate, so provenance
%   and verdict thresholds cannot disagree. This makes every blocking
%   gate's failure testable.
%
%   OTHER MODES (no ladder solve):
%
%     T2 = pair_transfer('second', T, o)   re-run the second-order block
%                                          on a stored ok result; o takes
%                                          .msHardCapSec .K .secondDelaySec
%     M  = pair_transfer('mass', tf_nd, thrustN, o)   the all-burn mass
%                                          bookkeeping the solve uses (o
%                                          .ispS .m0kg): M.mfFrac .mf_kg
%                                          .dV_kms .ok .reason ('mass
%                                          exhausted' when mfFrac <= 0)
%
%   TEST SEAM: opts.secondDelaySec [0] pauses the second-order block that
%   long before it starts (on the worker when fenced), so a test can force
%   the cap without a hung computation.
%
%   Self-demo (nargin == 0): DRO tau 1 -> ELFO raanDeg 0, phases (0, 0),
%   15 N (stage 0 reference t_f 0.0900542599967 ND); about 4 minutes cold.
%
%  ASSUMPTIONS / NOTES:
%
% - Assumes costate_common, oclib, CasADi and pumpkyn are on the path; it
%   calls no setup script.
% - msHardCapSec > 0: the fence pool is the open pool if any (kept), else
%   a 1-worker capped_pool opened BEFORE the engine call (so the engine
%   and the second-order block share it) and deleted when this call
%   returns. Workers inherit the client path at pool creation, so a pool
%   opened before costate_common was on the path fails the fenced call on
%   the worker (reported as 'error: ...', not blocking).
% - T.flight.Y is fly_transfer's (pumpkyn tfMinProp) [M x 14]: columns 1:3
%   position, 4:6 velocity (rotating frame, ND), 7 mass fraction (1 at
%   t = 0), 8:10 lambda_r, 11:13 lambda_v, 14 lambda_m. alpha =
%   -lambda_v/|lambda_v| (tfMinEoM's thrust direction, all-burn).
% - Scratch: the engine writes a one-cell sheet and a log under
%   opts.outDir. The sheet is deleted on every path, a throw included
%   (onCleanup), unless keepScratch; the log is deleted on success and
%   kept on failure (its path in T.scratch.log); the folder is removed
%   when this call created it and it is empty.
%
%% Inputs:
%
%  dep, arr                 struct                  .family (char, a
%                                                   get_family_orbit family)
%                                                   .params (struct for it;
%                                                   [] or absent = defaults)
%
%  sD, sA                   double                  departure / arrival phase
%                                                   fraction (phase_state)
%
%  thrustN                  double                  thrust, N (positive): the
%                                                   returned rung
%
%  opts                     struct (optional)       .m0kg [150] .ispS [1710]
%                                                   .rungs [see LADDER]
%                                                   .tf0 [0.3] .thrLock [true]
%                                                   .maxCpuSec [600]
%                                                   .maxIter [6000] .N [400]
%                                                   .floorKm [500] Moon
%                                                   .floorEarthKm [[]] Earth
%                                                   .gateKm [100]
%                                                   .directGateKm [inf]
%                                                   .msHardCapSec [300]
%                                                   .accTol [1e-6] .K [24]
%                                                   junction segments
%                                                   .secondOrder [true]
%                                                   .outDir [tempname]
%                                                   .keepScratch [false]
%                                                   .planOnly [false] return
%                                                   the plan, no solve
%                                                   .quiet [false]
%                                                   .secondDelaySec [0] test
%                                                   seam (see TEST SEAM)
%                                                   unknown names throw
%                                                   pair_transfer:opts
%
%% Outputs:
%
%  T                        struct                  .ok .reason (see REASONS)
%                                                   .dep .arr .sD .sA
%                                                   .thrustN .m0kg .ispS
%                                                   .rv0 .rvf [1 x 6]
%                                                   .z8 [8 x 1] .tf_nd .tf_d
%                                                   .mf_kg .dV_kms
%                                                   .rungs .tfPerRung
%                                                   .okPerRung
%                                                   .flownKmPerRung [1 x nR]
%                                                   .flight (.t [M x 1]
%                                                   .Y [M x 14] .alpha
%                                                   [M x 3])
%                                                   .junctions (.tj
%                                                   [1 x K+1] .yj
%                                                   [14 x K+1])
%                                                   .gates (.msNormR
%                                                   .flownKm .flownMps
%                                                   .engineFlownKm .accDz
%                                                   .admissible
%                                                   .admissibleReason
%                                                   .moonMinKm .moonMinT
%                                                   .moonSampleMinKm
%                                                   .moonNFailed .moonReason
%                                                   .earthMinKm .earthMinT
%                                                   .earthSampleMinKm
%                                                   .earthNFailed
%                                                   .earthReason .conjPass
%                                                   (1 / 0 / NaN)
%                                                   .conjCovered .conjReason
%                                                   .minLamV .minQmt .dimS
%                                                   .secondStatus .pass
%                                                   .failed {cellstr})
%                                                   .thresholds (.gateKm
%                                                   .accTol .floorKm
%                                                   .floorEarthKm)
%                                                   .entry (library format:
%                                                   .sD .sA .z8 .tf_nd
%                                                   .thrustN .rv0 .rvf
%                                                   .junctions .verdicts
%                                                   .provenance with
%                                                   .engineOpts, .created
%                                                   datetime and .regate
%                                                   records)
%                                                   .engineOpts .scratch
%                                                   (.outDir .sheet .log)
%                                                   .wallS
%
%% References:
%   [1] GTO_DRO/lib/gto_dro_entry.m -- the single-entry driver generalized.
%   [2] lunar_pnt_mission/lib/leg1_gates.m -- re-flight gates and the
%       second-order block.
%   [3] thrust_ladder_library.m -- the run_capped fence pattern.
%
%% Revision History:
%  M. Casey                                                   (c) 11/19/2026
%  M. Casey  review round 1: fenced second order, engine-     (c) 11/19/2026
%            error and mass guards, regate provenance, modes
%  Copyright Coorbital Inc.
%% ------------------------ Begin Code Sequence ---------------------------

if nargin == 0
    dep = struct('family', 'dro',  'params', struct('tau', 1));
    arr = struct('family', 'elfo', 'params', struct('raanDeg', 0));
    T = pair_transfer(dep, arr, 0, 0, 15, struct());
    return
end
if ischar(dep) || isstring(dep)
    if nargin < 3, sD = struct(); end
    switch char(dep)
        case 'regate'
            T = regate(arr, sD);
        case 'second'
            T = second_mode(arr, sD);
        case 'mass'
            if nargin < 4, sA = struct(); end
            o = knobs(sA);
            T = mass_of(arr, sD, o.ispS, o.m0kg);
        otherwise
            error('pair_transfer:mode', 'unknown mode "%s" (regate, second, mass)', char(dep));
    end
    return
end
if nargin < 6, opts = struct(); end
o = knobs(opts);
tStart = tic;

%% Inputs:
dep = check_end(dep, 'dep');
arr = check_end(arr, 'arr');
if ~(isscalar(thrustN) && isnumeric(thrustN) && isreal(thrustN) && isfinite(thrustN) && thrustN > 0)
    error('pair_transfer:thrust', 'thrustN must be a positive finite scalar (N)');
end
if ~(isscalar(sD) && isscalar(sA) && isfinite(sD) && isfinite(sA))
    error('pair_transfer:phase', 'sD and sA must be finite scalars (phase fractions)');
end
rungs = rung_list(o, thrustN);
nR    = numel(rungs);

muStar = 0.012150585609624;
 lStar = 389703.264829278;
 tStar = 382981.289129055;

%% Scratch and engine options:
       outDir = o.outDir;
if isempty(outDir), outDir = tempname; end
          tag = sprintf('pt_%s_%s_sD%.6f_sA%.6f_%gN', dep.family, arr.family, sD, sA, thrustN);
       outMat = fullfile(outDir, [tag '.mat']);
      logFile = fullfile(outDir, [tag '.log']);
optsEngine = struct( ...
    'rungs',        rungs, ...
    'ispS',         o.ispS, ...
    'm0kg',         o.m0kg, ...
    'N',            o.N, ...
    'floorKm',      o.floorKm, ...
    'floorEarthKm', o.floorEarthKm, ...
    'maxIter',      o.maxIter, ...
    'gateKm',       o.gateKm, ...
    'directGateKm', o.directGateKm, ...
    'accTol',       o.accTol, ...
    'nD',           1, ...
    'nA',           1, ...
    'sD0',          sD, ...
    'sA0',          sA, ...
    'muStar',       muStar, 'lStar', lStar, 'tStar', tStar, ...
    'depFamily',    dep.family, ...
    'depParams',    dep.params, ...
    'arrFamily',    arr.family, ...
    'arrParams',    arr.params, ...
    'thrLock',      o.thrLock, ...
    'tf0',          o.tf0, ...
    'maxCpuSec',    o.maxCpuSec, ...
    'msHardCapSec', o.msHardCapSec, ...
    'resume',       false, ...
    'logFile',      logFile);

%% The result skeleton (what a failed or planned call returns):
T = skeleton(dep, arr, sD, sA, thrustN, o, rungs, optsEngine);
T.scratch = struct('outDir', outDir, 'sheet', outMat, 'log', '');
if o.planOnly
    T.reason = 'planOnly';
    T.wallS  = toc(tStart);
    return
end

%% Endpoints by THE shared endpoint rule (costate_common/phase_state):
[tD, rvD] = get_family_orbit(dep.family, dep.params);
[tA, rvA] = get_family_orbit(arr.family, arr.params);
depAt = phase_state(tD, rvD);   arrAt = phase_state(tA, rvA);
T.rv0 = depAt(sD).';
T.rvf = arrAt(sA).';

%% Solve (one cell, top rung down to thrustN):
madeDir = ~isfolder(outDir);
if madeDir, mkdir(outDir); end
if isfile(outMat),  delete(outMat);  end
if isfile(logFile), delete(logFile); end
% the sheet goes on EVERY path (a throw included), the folder when empty
cleanScratch = onCleanup(@() scratch_cleanup(outMat, outDir, madeDir, o.keepScratch));
% ONE fence pool for the engine and the second-order block (deleted on
% return only when opened here)
fenceHold = fence_pool(o.msHardCapSec);   % .pool, .dropper (onCleanup)
fence = fenceHold.pool;
try
    thrust_ladder_library(outMat, optsEngine);
catch errEngine
    T.reason = ['engine error: ' errEngine.message];
    T = apply_gates(T);
    T.entry = make_entry(T);
    if isfile(logFile), T.scratch.log = logFile; end
    T.wallS = toc(tStart);
    if ~o.quiet, print_summary(T); end
    clear fenceHold cleanScratch
    return
end
P = load(outMat);
if ~isequal(P.rungs, rungs)
    error('pair_transfer:rungIndex', 'engine rungs %s differ from the plan %s', ...
          mat2str(P.rungs), mat2str(rungs));
end
okR  = logical(reshape(P.OK(1, 1, :), 1, nR));
tfR  = reshape(P.TF(1, 1, :), 1, nR);
flyR = reshape(P.FLYKM(1, 1, :), 1, nR);
T.tfPerRung      = tfR;
T.okPerRung      = okR;
T.flownKmPerRung = flyR;
T.ok             = okR(nR);
T.reason         = stop_reason(okR, tfR, flyR, rungs, o.directGateKm);
T.gates.msNormR       = P.RES(1, 1, nR);
T.gates.engineFlownKm = P.FLYKM(1, 1, nR);
T.gates.accDz         = P.ACCDZ(1, 1, nR);

if T.ok
    T.z8    = P.Z8(:, 1, 1, nR);
    T.tf_nd = T.z8(8);
    T.tf_d  = T.tf_nd*tStar/86400;
    Ms      = mass_of(T.tf_nd, thrustN, o.ispS, o.m0kg);
    T.mf_kg  = Ms.mf_kg;
    T.dV_kms = Ms.dV_kms;
    if ~Ms.ok                                           % guard: no mass left
        T.ok = false;
        T.reason = Ms.reason;
    end
end
if T.ok
    phys = Ms.phys;

    %% Independent re-flight and its admissibility:
    F = fly_transfer(T.z8, T.rv0, T.rvf, phys);
    lamV = F.Y(:, 11:13);
    T.flight = struct('t', F.t, 'Y', F.Y, 'alpha', -lamV./vecnorm(lamV, 2, 2));
    T.gates.flownKm  = F.flyKm;
    T.gates.flownMps = F.flyVms;
    T.gates.admissible       = logical(F.admissibility.ok);
    T.gates.admissibleReason = char(string(F.admissibility.reason));

    %% Junction states (library entry):
    [yj, tj] = flight_to_junctions(F.t, F.Y, o.K);
    T.junctions = struct('tj', tj, 'yj', yj);

    %% Located clearance over the Moon and the Earth:
    pl = struct('lStar', lStar, 'phys', phys);
    Cm = located_min_altitude(F.t, F.Y, [1 - muStar 0 0], 1737.4, o.floorKm, pl);
    if isempty(o.floorEarthKm), fE = -inf; else, fE = o.floorEarthKm; end
    Ce = located_min_altitude(F.t, F.Y, [-muStar 0 0], 6378.137, fE, pl);
    T.gates.moonMinKm        = Cm.hMinKm;   T.gates.moonMinT  = Cm.tAtMin;
    T.gates.moonSampleMinKm  = Cm.sampleMinKm;
    T.gates.moonNFailed      = Cm.nFailed;  T.gates.moonReason = Cm.reason;
    T.gates.earthMinKm       = Ce.hMinKm;   T.gates.earthMinT  = Ce.tAtMin;
    T.gates.earthSampleMinKm = Ce.sampleMinKm;
    T.gates.earthNFailed     = Ce.nFailed;  T.gates.earthReason = Ce.reason;

    %% Second order (reported):
    if o.secondOrder
        T.gates = second_order(T.gates, T.z8, T.rv0, T.rvf, phys, o.K, ...
                               o.msHardCapSec, o.secondDelaySec, fence);
    else
        T.gates.secondStatus = 'skipped';
    end
end
clear fenceHold                                        % pool opened here goes now

%% Verdicts and the library entry:
T = apply_gates(T);
T.entry = make_entry(T);

%% Scratch (the sheet and an empty folder go in scratch_cleanup):
if T.ok && ~o.keepScratch
    drop_file(logFile);
elseif ~T.ok
    if isfile(logFile)
        T.scratch.log = logFile;
        if ~o.quiet
            fprintf('pair_transfer: %s -- engine log kept at %s\n', T.reason, logFile);
        end
    end
else
    T.scratch.log = logFile;
end
T.wallS = toc(tStart);
if ~o.quiet
    print_summary(T);
end
end

% ------------------------------------------------------------------------
function T = skeleton(dep, arr, sD, sA, thrustN, o, rungs, optsEngine)
% SKELETON  The result of a call that has not (yet) solved anything.
%   INPUTS: dep, arr; sD, sA; thrustN; o knobs; rungs [1 x nR];
%   optsEngine struct.  OUTPUTS: T struct (every field, NaN-filled).
nR = numel(rungs);
K  = o.K;
gates = struct('msNormR', NaN, 'flownKm', NaN, 'flownMps', NaN, 'engineFlownKm', NaN, ...
    'accDz', NaN, 'admissible', false, 'admissibleReason', 'not flown', ...
    'moonMinKm', NaN, 'moonMinT', NaN, 'moonSampleMinKm', NaN, 'moonNFailed', 0, 'moonReason', '', ...
    'earthMinKm', NaN, 'earthMinT', NaN, 'earthSampleMinKm', NaN, 'earthNFailed', 0, 'earthReason', '', ...
    'conjPass', NaN, 'conjCovered', NaN, 'conjReason', 'not run', ...
    'minLamV', NaN, 'minQmt', NaN, 'dimS', NaN, 'secondStatus', 'skipped', ...
    'pass', false, 'failed', {{}});
T = struct('ok', false, 'reason', '', 'dep', dep, 'arr', arr, 'sD', sD, 'sA', sA, ...
    'thrustN', thrustN, 'm0kg', o.m0kg, 'ispS', o.ispS, 'rv0', nan(1, 6), 'rvf', nan(1, 6), ...
    'z8', nan(8, 1), 'tf_nd', NaN, 'tf_d', NaN, 'mf_kg', NaN, 'dV_kms', NaN, ...
    'rungs', rungs, 'tfPerRung', nan(1, nR), 'okPerRung', false(1, nR), ...
    'flownKmPerRung', nan(1, nR), ...
    'flight', struct('t', zeros(0, 1), 'Y', zeros(0, 14), 'alpha', zeros(0, 3)), ...
    'junctions', struct('tj', nan(1, K + 1), 'yj', nan(14, K + 1)), ...
    'gates', gates, ...
    'thresholds', struct('gateKm', o.gateKm, 'accTol', o.accTol, 'floorKm', o.floorKm, ...
                         'floorEarthKm', o.floorEarthKm), ...
    'entry', struct(), 'engineOpts', optsEngine, 'scratch', struct(), 'wallS', NaN);
end

% ------------------------------------------------------------------------
function T = apply_gates(T)
% APPLY_GATES  THE pass logic (used by the solve and by 'regate'): every
%   blocking gate against T.thresholds.  INPUTS: T.  OUTPUTS: T with
%   .gates.pass and .gates.failed set.
g  = T.gates;
th = T.thresholds;
failed = {};
if ~T.ok,                     failed{end+1} = 'ok';         end
if ~(g.flownKm < th.gateKm),  failed{end+1} = 'flown';      end
if ~(g.accDz < th.accTol),    failed{end+1} = 'accDz';      end
if ~g.admissible,             failed{end+1} = 'admissible'; end
if ~(g.moonMinKm >= th.floorKm && g.moonNFailed == 0)
    failed{end+1} = 'moon';
end
if ~isempty(th.floorEarthKm) && ~(g.earthMinKm >= th.floorEarthKm && g.earthNFailed == 0)
    failed{end+1} = 'earth';
end
T.gates.failed = failed;
T.gates.pass   = isempty(failed);
end

% ------------------------------------------------------------------------
function T = regate(T, th)
% REGATE  Re-apply the pass logic with new thresholds, no solve; the
%   applied thresholds are recorded in T.entry.provenance.regate.
%   INPUTS: T a pair_transfer result; th struct (.gateKm .accTol .floorKm
%   .floorEarthKm).  OUTPUTS: T with .thresholds, .gates.pass/.failed,
%   .entry.verdicts and .entry.provenance.regate updated.
if ~(isstruct(T) && isfield(T, 'gates') && isfield(T, 'thresholds'))
    error('pair_transfer:regate', 'regate needs a pair_transfer result');
end
if isempty(th), th = struct(); end
fn = fieldnames(th);
for k = 1:numel(fn)
    if ~isfield(T.thresholds, fn{k})
        error('pair_transfer:regate', 'unknown threshold "%s" (gateKm, accTol, floorKm, floorEarthKm)', fn{k});
    end
    T.thresholds.(fn{k}) = th.(fn{k});
end
T = apply_gates(T);
if isfield(T, 'entry') && isfield(T.entry, 'verdicts')
    T.entry.verdicts = verdicts_of(T);
    rec = struct('thresholds', T.thresholds, 'overrides', th, ...
                 'at', datetime('now', 'TimeZone', 'local'));
    if isfield(T.entry.provenance, 'regate') && ~isempty(T.entry.provenance.regate)
        T.entry.provenance.regate = [T.entry.provenance.regate, {rec}];
    else
        T.entry.provenance.regate = {rec};
    end
end
end

% ------------------------------------------------------------------------
function E = make_entry(T)
% MAKE_ENTRY  The solved transfer in library-entry format (standing rule:
%   every solved transfer is captured).  INPUTS: T.  OUTPUTS: E struct.
eo = rmfield(T.engineOpts, 'logFile');
prov = struct('function', 'pair_transfer', 'dep', T.dep, 'arr', T.arr, ...
    'engine', 'thrust_ladder_library (direct -> ms_tfmin -> tfMin acceptance)', ...
    'engineOpts', eo, 'created', datetime('now', 'TimeZone', 'local'), ...
    'matlab', version, 'regate', {{}});
E = struct('sD', T.sD, 'sA', T.sA, 'z8', T.z8, 'tf_nd', T.tf_nd, 'thrustN', T.thrustN, ...
    'm0kg', T.m0kg, 'ispS', T.ispS, 'rv0', T.rv0, 'rvf', T.rvf, 'mf_kg', T.mf_kg, ...
    'dV_kms', T.dV_kms, 'junctions', T.junctions, 'verdicts', verdicts_of(T), ...
    'provenance', prov);
end

% ------------------------------------------------------------------------
function V = verdicts_of(T)
% VERDICTS_OF  What was checked on this entry, and what was not.
%   INPUTS: T.  OUTPUTS: V struct.
g = T.gates;
V = struct('ok', T.ok, 'reason', T.reason, 'pass', g.pass, 'failed', {g.failed}, ...
    'flownKm', g.flownKm, 'flownMps', g.flownMps, 'accDz', g.accDz, 'msNormR', g.msNormR, ...
    'admissible', g.admissible, 'moonMinKm', g.moonMinKm, 'earthMinKm', g.earthMinKm, ...
    'conjPass', g.conjPass, 'conjCovered', g.conjCovered, 'conjReason', g.conjReason, ...
    'minLamV', g.minLamV, 'minQmt', g.minQmt, 'dimS', g.dimS, ...
    'secondStatus', g.secondStatus, 'thresholds', T.thresholds, ...
    'notChecked', 'lift margin, H6, dense conjugate scan (the library status run)');
end

% ------------------------------------------------------------------------
function g = second_order(g, z8, rv0, rvf, phys, K, capSec, delaySec, pool)
% SECOND_ORDER  The second-order block, reported: fenced on a worker under
%   run_capped when capSec > 0 and a pool is at hand, else in-process.
%   INPUTS: g gates; z8 [8 x 1]; rv0, rvf [1 x 6]; phys; K segments;
%   capSec hard cap (s) or []/0; delaySec test seam (s); pool
%   parallel.Pool or [].  OUTPUTS: g with .conjPass .conjCovered
%   .conjReason .minLamV .minQmt .dimS .secondStatus.
fenced = ~isempty(capSec) && capSec > 0 && ~isempty(pool) && isvalid(pool);
if ~fenced
    S = second_core(z8, rv0, rvf, phys, K, delaySec);
else
    t0 = tic;
    try
        [okRun, S, cause] = run_capped(pool, @second_core, 1, capSec, z8, rv0, rvf, phys, K, delaySec);
    catch errCap                            % run_capped:poolUnusable
        okRun = false;  cause = [errCap.identifier ': ' errCap.message];
    end
    if ~okRun
        S = struct('conjPass', NaN, 'conjCovered', NaN, 'minLamV', NaN, 'minQmt', NaN, 'dimS', NaN);
        if strcmp(cause, 'capped')
            S.conjReason   = sprintf('capped after %.1f s', toc(t0));
            S.secondStatus = 'capped';
        else
            S.conjReason   = ['error: ' cause];
            S.secondStatus = 'error';
        end
    end
end
for f = {'conjPass', 'conjCovered', 'conjReason', 'minLamV', 'minQmt', 'dimS', 'secondStatus'}
    g.(f{1}) = S.(f{1});
end
end

% ------------------------------------------------------------------------
function S = second_core(z8, rv0, rvf, phys, K, delaySec)
% SECOND_CORE  Conjugate re-solve + mintime_hypothesis_gates (the
%   leg1_gates block), errors caught and returned. Runs in-process or on a
%   run_capped worker.  INPUTS: z8 [8 x 1]; rv0, rvf [1 x 6]; phys; K;
%   delaySec test seam (s).  OUTPUTS: S struct .conjPass .conjCovered
%   .conjReason .minLamV .minQmt .dimS .secondStatus.
S = struct('conjPass', NaN, 'conjCovered', NaN, 'conjReason', 'not run', ...
           'minLamV', NaN, 'minQmt', NaN, 'dimS', NaN, 'secondStatus', 'error');
if delaySec > 0, pause(delaySec); end
ws = warning('off', 'MATLAB:ode45:IntegrationTolNotMet');
restoreWarn = onCleanup(@() warning(ws));
Tnd = phys.Tnd;  cnd = phys.cnd;  mu = phys.muStar;
try
    seed = seed_from_z8(z8, rv0, K, Tnd, cnd, mu);
    [zc, info] = ms_tfmin(rv0, rvf, seed, Tnd, cnd, mu, ...
                          struct('tolR', 3e-11, 'wallSec', 600, 'conjTest', true));
    dz = max(abs(zc - z8));
    reproduced = info.converged && dz < 1e-6;
    if reproduced
        S.conjPass    = double(info.conj.pass);
        S.conjCovered = double(info.conj.covered);
        S.conjReason  = sprintf('reproduced; %s: %s', info.conj.verdict, info.conj.reason);
    else
        S.conjReason  = sprintf('not reproduced (converged %d, |dz| = %.1e)', info.converged, dz);
    end
    hy = mintime_hypothesis_gates(z8, rv0, Tnd, cnd, mu, struct('relTol', 1e-12));
    S.minLamV = hy.minLamV;  S.minQmt = hy.minQmt;  S.dimS = hy.dimS;
    passS = reproduced && S.conjPass == 1 && S.conjCovered == 1 && ...
            S.minLamV > 0 && S.minQmt > 0 && S.dimS == 1;
    if passS
        S.secondStatus = 'pass';
    elseif reproduced
        S.secondStatus = 'refuted';
    else
        S.secondStatus = 'notReproduced';
    end
catch err
    S.conjReason   = ['error: ' err.message];
    S.secondStatus = 'error';
end
clear restoreWarn
end

% ------------------------------------------------------------------------
function T = second_mode(T, opts)
% SECOND_MODE  Re-run the second-order block on a stored ok result (no
%   ladder solve).  INPUTS: T pair_transfer result; opts (.msHardCapSec
%   .K .secondDelaySec, knobs defaults otherwise).  OUTPUTS: T with the
%   second-order gates and .entry.verdicts updated.
if ~(isstruct(T) && isfield(T, 'ok') && T.ok)
    error('pair_transfer:second', 'the second mode needs an ok pair_transfer result');
end
if isempty(opts), opts = struct(); end
o  = knobs(opts);
Ms = mass_of(T.tf_nd, T.thrustN, T.ispS, T.m0kg);
fenceHold = fence_pool(o.msHardCapSec);
T.gates = second_order(T.gates, T.z8, T.rv0, T.rvf, Ms.phys, o.K, ...
                       o.msHardCapSec, o.secondDelaySec, fenceHold.pool);
clear fenceHold
if isfield(T, 'entry') && isfield(T.entry, 'verdicts')
    T.entry.verdicts = verdicts_of(T);
end
end

% ------------------------------------------------------------------------
function M = mass_of(tf_nd, thrustN, ispS, m0kg)
% MASS_OF  All-burn mass bookkeeping m(t) = 1 - T t / c, guarded.
%   INPUTS: tf_nd ND; thrustN N; ispS s; m0kg kg.  OUTPUTS: M struct
%   .mfFrac .mf_kg .dV_kms (NaN when exhausted) .ok .reason ('' or 'mass
%   exhausted') .phys (.Tnd .cnd .muStar .lStar .tStar .m0kg).
muStar = 0.012150585609624;
 lStar = 389703.264829278;
 tStar = 382981.289129055;
ndp    = nd_propulsion(thrustN, ispS, m0kg, lStar, tStar);
mfFrac = 1 - ndp.Tnd*tf_nd/ndp.cnd;
M = struct('mfFrac', mfFrac, 'mf_kg', NaN, 'dV_kms', NaN, 'ok', false, 'reason', 'mass exhausted', ...
           'phys', struct('Tnd', ndp.Tnd, 'cnd', ndp.cnd, 'muStar', muStar, 'lStar', lStar, ...
                          'tStar', tStar, 'm0kg', m0kg));
if isfinite(mfFrac) && mfFrac > 0
    M.mf_kg  = m0kg*mfFrac;
    M.dV_kms = ndp.cnd*log(1/mfFrac)*lStar/tStar;
    M.ok     = true;
    M.reason = '';
end
end

% ------------------------------------------------------------------------
function H = fence_pool(capSec)
% FENCE_POOL  The fence pool: the open pool (kept), else a capped_pool(1)
%   opened here and deleted when the caller clears H; pool [] when the
%   fence is off or no pool is obtainable.  INPUTS: capSec []/0 = off.
%   OUTPUTS: H struct .pool (parallel.Pool or []) .dropper (onCleanup
%   that deletes a pool opened here, or []).
H = struct('pool', [], 'dropper', []);
if isempty(capSec) || ~(capSec > 0)
    return
end
H.pool = current_pool();
if isempty(H.pool) || ~isvalid(H.pool)
    H.pool = capped_pool(1);
    if ~isempty(H.pool)
        pool = H.pool;
        H.dropper = onCleanup(@() delete_if_valid(pool));
    end
end
end

% ------------------------------------------------------------------------
function delete_if_valid(pool)
% DELETE_IF_VALID  Delete a pool opened here, if still valid.
%   INPUTS: pool parallel.Pool.  OUTPUTS: none.
if ~isempty(pool) && isvalid(pool)
    try delete(pool); catch, end
end
end

% ------------------------------------------------------------------------
function scratch_cleanup(outMat, outDir, madeDir, keep)
% SCRATCH_CLEANUP  Drop the engine sheet (unless keep) and the folder this
%   call made once it is empty.  INPUTS: outMat, outDir char; madeDir,
%   keep logical.  OUTPUTS: none.
if ~keep && isfile(outMat), delete(outMat); end
if madeDir && isfolder(outDir) && numel(dir(outDir)) <= 2   % only '.' and '..'
    rmdir(outDir);
end
end

% ------------------------------------------------------------------------
function r = stop_reason(okR, tfR, flyR, rungs, directGateKm)
% STOP_REASON  Why the returned rung has no accepted entry ('' if it has).
%   The ladder breaks on a direct failure (no t_f recorded) or on the
%   direct flown gate (t_f recorded, flown above directGateKm).
%   INPUTS: okR, tfR, flyR [1 x nR]; rungs; directGateKm.  OUTPUTS: r char.
nR = numel(rungs);
if okR(nR)
    r = '';
    return
end
kLast = find(isfinite(tfR), 1, 'last');
if isempty(kLast)
    kStop = 1;
elseif flyR(kLast) > directGateKm && ~okR(kLast)
    kStop = kLast;
else
    kStop = kLast + 1;
end
if kStop > nR
    r = 'not accepted';
elseif kStop == 1
    r = 'top rung failed';
else
    r = sprintf('ladder stopped at %g N', rungs(kStop));
end
end

% ------------------------------------------------------------------------
function rungs = rung_list(o, thrustN)
% RUNG_LIST  The ladder: opts.rungs as given (ending at thrustN), else the
%   standard rungs above thrustN, then thrustN.  INPUTS: o knobs; thrustN.
%   OUTPUTS: rungs [1 x nR], strictly descending.
if isempty(o.rungs)
    stdR  = [15 12 10 7 5 3 2 1.5 1 0.75 0.5];
    rungs = [stdR(stdR > thrustN), thrustN];
else
    rungs = o.rungs(:).';
    if ~(isnumeric(rungs) && all(isfinite(rungs)) && all(rungs > 0) && all(diff(rungs) < 0))
        error('pair_transfer:rungs', 'opts.rungs must be positive, finite and strictly descending');
    end
    if rungs(end) ~= thrustN
        error('pair_transfer:rungs', 'opts.rungs must end at thrustN = %g N (ends at %g N)', ...
              thrustN, rungs(end));
    end
end
end

% ------------------------------------------------------------------------
function e = check_end(e, name)
% CHECK_END  An orbit end: struct with .family (char) and .params (struct).
%   INPUTS: e; name ('dep' | 'arr').  OUTPUTS: e normalized.
if ~(isstruct(e) && isscalar(e) && isfield(e, 'family') && ...
        (ischar(e.family) || isstring(e.family)))
    error('pair_transfer:end', '%s must be a struct with .family (char) and .params (struct)', name);
end
e.family = char(e.family);
if ~isfield(e, 'params') || isempty(e.params), e.params = struct(); end
if ~isstruct(e.params)
    error('pair_transfer:end', '%s.params must be a struct', name);
end
end

% ------------------------------------------------------------------------
function o = knobs(opts)
% KNOBS  Defaults; unknown names throw pair_transfer:opts.
%   INPUTS: opts struct.  OUTPUTS: o struct.
o = struct( ...
    'm0kg',           150, ...
    'ispS',           1710, ...
    'rungs',          [], ...
    'tf0',            0.3, ...
    'thrLock',        true, ...
    'maxCpuSec',      600, ...
    'maxIter',        6000, ...
    'N',              400, ...
    'floorKm',        500, ...
    'floorEarthKm',   [], ...
    'gateKm',         100, ...
    'directGateKm',   inf, ...
    'msHardCapSec',   300, ...
    'accTol',         1e-6, ...
    'K',              24, ...
    'secondOrder',    true, ...
    'outDir',         '', ...
    'keepScratch',    false, ...
    'planOnly',       false, ...
    'quiet',          false, ...
    'secondDelaySec', 0);
fn = fieldnames(opts);
for k = 1:numel(fn)
    if ~isfield(o, fn{k})
        error('pair_transfer:opts', 'unknown option "%s"', fn{k});
    end
    o.(fn{k}) = opts.(fn{k});
end
end

% ------------------------------------------------------------------------
function print_summary(T)
% PRINT_SUMMARY  One block of the measured numbers.  INPUTS: T.
%   OUTPUTS: none (prints).
g = T.gates;
fprintf('pair_transfer: %s -> %s, phases (%.4f, %.4f), %g N: ok %d %s (%.0f s)\n', ...
        T.dep.family, T.arr.family, T.sD, T.sA, T.thrustN, T.ok, T.reason, T.wallS);
if ~T.ok, return, end
fprintf('  t_f %.10f ND = %.4f d, m_f %.3f kg, dV %.4f km/s\n', T.tf_nd, T.tf_d, T.mf_kg, T.dV_kms);
fprintf('  flown %.2e km %.2e m/s, accDz %.1e, |R| %.1e, admissible %d (%s)\n', g.flownKm, ...
        g.flownMps, g.accDz, g.msNormR, g.admissible, g.admissibleReason);
fprintf('  Moon located %.2f km (sample %.2f), Earth located %.1f km (sample %.1f)\n', ...
        g.moonMinKm, g.moonSampleMinKm, g.earthMinKm, g.earthSampleMinKm);
fprintf('  second %s: conjPass %g covered %g, minLamV %.3g minQmt %.3g dimS %g -- %s\n', ...
        g.secondStatus, g.conjPass, g.conjCovered, g.minLamV, g.minQmt, g.dimS, g.conjReason);
fprintf('  gates.pass %d%s\n', g.pass, failed_str(g.failed));
end

% ------------------------------------------------------------------------
function s = failed_str(f)
% FAILED_STR  ' (failed: a, b)' or ''.  INPUTS: f cellstr.  OUTPUTS: s.
if isempty(f), s = ''; else, s = [' (failed: ' strjoin(f, ', ') ')']; end
end

% ------------------------------------------------------------------------
function drop_file(f)
% DROP_FILE  Delete f if it exists.  INPUTS: f char.  OUTPUTS: none.
if isfile(f), delete(f); end
end
