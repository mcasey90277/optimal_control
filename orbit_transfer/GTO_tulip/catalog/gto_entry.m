function [entry, gates] = gto_entry(orientDeg, Np, depFrac, arrFrac, thrustN, opts)
%% Purpose:
%
%   THE SINGLE-ENTRY DRIVER for the GTO -> tulip costate campaign (Mike's
%   product request 2026-08-26): solve ONE (departure phase x arrival
%   phase) cell at a single requested thrust level, through the SAME
%   engine the sheet sweep uses (DRO_tulip/indirect/thrust_ladder_library,
%   family-agnostic via costate_common/get_family_orbit) -- never a
%   bespoke solver. Runs a 1x1 phasing grid.
%
%   LADDER DIRECTION (fixed 2026-08-27, review finding -- the original
%   construction was inverted): solves the standard ladder from 15 N DOWN
%   to the requested thrustN (higher rungs are warm-start stepping stones
%   only, never returned) and returns the entry AT thrustN. For thrustN =
%   15 this is a single solve, since there is no standard rung above 15 N.
%   A request for a deep rung (e.g. 1 N) is never cold-solved directly --
%   it is reached by continuation down from 15 N, exactly like the
%   campaign sweep, because a cold cold-start at a deep rung is precisely
%   the non-convergent case the pilot measured.
%
%   WARM RECIPE (matches run_gto_catalog's pass-A recipe, 2026-08-26
%   diagnostic-A): thrLock=true, tf0=0.30 (ND, seeds the ladder's TOP
%   rung only), maxCpuSec=600 -- the campaign abandoned the cold defaults
%   (free throttle, casadi_mintime_dro's own 4.0 ND seed) after measuring
%   67/72 GTO cells unreachable with them; this driver never falls back to
%   them.
%
%   Self-demo (nargin==0): orient 0 deg, Np 7, 15 N. Departure phase is
%   PERIGEE (0.0), not apogee (0.5), even though the original spec named
%   apogee: Task 2 of this campaign measured an apogee cold-start at this
%   exact (orient, Np, thrust) miss the flown-arrival gate by ~5 orders of
%   magnitude (11.95e6 km vs a 100 km gate) -- a bad local minimum for a
%   15 N cold start at that phase, not a scale/unit bug (the GTO locus
%   itself checks out physically). Perigee converges cleanly on the first
%   try. See task-2-report.md / task-3-report.md for the numbers.
%
%% Inputs:
%
%  orientDeg                double                  GTO departure-locus
%                                                   orientation (deg,
%                                                   Earth->Moon line to
%                                                   perigee direction)
%
%  Np                       double                  Tulip petal count
%
%  depFrac                  double                  Departure phase: TIME
%                                                   fraction since perigee
%                                                   over one Kepler period
%                                                   (0 = perigee, 0.5 =
%                                                   apogee)
%
%  arrFrac                  double                  Arrival phase: time
%                                                   fraction along the
%                                                   tulip's own period
%
%  thrustN                  double                  Thrust level (N): the
%                                                   BOTTOM (returned) rung
%                                                   of the ladder; standard
%                                                   rungs above it are
%                                                   warm-start stepping
%                                                   stones only
%
%  opts                     struct                  (optional):
%                                                   .conjTest [false] re-
%                                                     solve ms_tfmin seeded
%                                                     at the accepted z8
%                                                     with conjTest=true
%                                                     (costate_common
%                                                     conj_catalog_pass
%                                                     pattern); sets
%                                                     gates.conjPass
%                                                   .writeSheet [false]
%                                                     merge this entry into
%                                                     results/catalog/
%                                                     gto_o<..>_Np<..>.mat
%                                                     on the STANDARD
%                                                     nD=12/nA=6 campaign
%                                                     grid (nearest-phase
%                                                     snap); the rung AXIS
%                                                     is taken from the
%                                                     sheet already on
%                                                     disk if one exists,
%                                                     else a brand-new
%                                                     sheet is stamped with
%                                                     the campaign's actual
%                                                     5-rung fleet default
%                                                     [15 12 10 7 5] (fixed
%                                                     2026-08-29 review --
%                                                     a hardcoded 9-rung
%                                                     grid here made
%                                                     writeSheet THROW on
%                                                     every shipped 5-rung
%                                                     sheet), creating the
%                                                     sheet file if absent
%                                                     -- writeSheet also
%                                                     banks the
%                                                     stepping-stone rungs
%                                                     solved along the way
%                                                     to thrustN, not just
%                                                     thrustN itself
%                                                   .thrLock [true],
%                                                     .tf0 [0.30],
%                                                     .maxCpuSec [600]
%                                                     the warm-recipe knobs
%                                                     (see LADDER DIRECTION
%                                                     above); overridable,
%                                                     default = campaign
%                                                     pass-A recipe
%                                                   .force [false] required
%                                                     to let .writeSheet
%                                                     overwrite a cell the
%                                                     REAL campaign sweep
%                                                     already marked OK
%                                                     (refused otherwise --
%                                                     see merge_entry_
%                                                     into_sheet)
%                                                   .m0kg [150], .ispS
%                                                     [1710] thruster
%                                                     params (campaign
%                                                     standard)
%                                                   .N [400], .floorKm
%                                                     [500], .maxIter
%                                                     [6000], .gateKm
%                                                     [100], .accTol
%                                                     [1e-6] solver/gate
%                                                     knobs, campaign
%                                                     pass-A recipe
%                                                     (maxIter=6000, not
%                                                     the engine's own
%                                                     3000 default -- see
%                                                     WARM RECIPE above)
%
%% Outputs:
%
%  entry                    struct                  .ok, .z8 [8x1],
%                                                   .tf_nd, .dV_kms,
%                                                   .mf_kg, .coordinates
%                                                   (orientDeg, Np,
%                                                   depFrac, arrFrac,
%                                                   thrustN)
%
%  gates                    struct                  .msNormR, .flownKm,
%                                                   .acceptDz, .conjPass
%                                                   (NaN unless
%                                                   opts.conjTest)
%
%% Revision History:
%  M. Casey                                                   (c) 08/26/2026
%  Copyright Coorbital Inc.
%% ------------------------ Begin Code Sequence ---------------------------

if nargin == 0
    orientDeg = 0;  Np = 7;  depFrac = 0.0;  arrFrac = 0;  thrustN = 15;
end
if nargin < 6, opts = struct(); end
conjTest   = fdef(opts, 'conjTest',   false);
writeSheet = fdef(opts, 'writeSheet', false);
force      = fdef(opts, 'force',      false);
m0kg       = fdef(opts, 'm0kg',       150);
ispS       = fdef(opts, 'ispS',       1710);
thrLock    = fdef(opts, 'thrLock',    true);   % campaign pass-A recipe
tf0Warm    = fdef(opts, 'tf0',        0.30);   % campaign pass-A recipe
maxCpuSec  = fdef(opts, 'maxCpuSec',  600);    % campaign pass-A recipe
% N/floorKm/gateKm/accTol already equal thrust_ladder_library's own
% defaults; maxIter does NOT -- the engine default (3000) is the cold
% recipe run_gto_catalog abandoned (measured: 3000 stopped a right-family
% cell one iteration budget short, defect 1.3e-11). All five are wired
% here anyway so gto_entry tracks the campaign recipe by construction
% instead of by coincidence with the engine's defaults:
N          = fdef(opts, 'N',          400);    % campaign recipe (= engine default)
floorKm    = fdef(opts, 'floorKm',    500);    % campaign recipe (= engine default)
maxIterOpt = fdef(opts, 'maxIter',    6000);   % campaign pass-A recipe (NOT the
                                                 % engine's own 3000 default)
gateKm     = fdef(opts, 'gateKm',     100);    % campaign recipe (= engine default)
accTol     = fdef(opts, 'accTol',     1e-6);   % campaign recipe (= engine default)
tulipPm    = -1;

setup_paths();

here   = fileparts(mfilename('fullpath'));
resDir = fullfile(here, 'results');
if ~isfolder(resDir), mkdir(resDir); end

%% Ladder from 15 N DOWN to thrustN: standard rungs ABOVE the request are
%% warm-start stepping stones (never returned), thrustN is always last --
%% the accepted/returned entry is the BOTTOM rung. If thrustN matches a
%% standard rung exactly, the '>' filter drops it from the stepping-stone
%% set so it is not duplicated:
rungsStd = [15 12 10 7 5 3 2 1.5 1];
rungs = [rungsStd(rungsStd > thrustN), thrustN];
krEntry = numel(rungs);                % thrustN's position in Q.rungs

outMat  = fullfile(resDir, 'gto_entry_scratch.mat');
logFile = fullfile(resDir, 'gto_entry_scratch.log');
if isfile(logFile), delete(logFile); end

muStar = 0.012150585609624;
lStar  = 389703.264829278;
tStar  = 382981.289129055;

optsEngine = struct( ...
    'rungs',     rungs, ...
    'ispS',      ispS, ...
    'm0kg',      m0kg, ...
    'N',         N, ...
    'floorKm',   floorKm, ...
    'maxIter',   maxIterOpt, ...
    'gateKm',    gateKm, ...
    'accTol',    accTol, ...
    'nD',        1, ...
    'nA',        1, ...
    'sD0',       depFrac, ...
    'sA0',       arrFrac, ...
    'muStar',    muStar, 'lStar', lStar, 'tStar', tStar, ...
    'depFamily', 'gto', ...
    'depParams', struct('orientDeg', orientDeg), ...
    'arrFamily', 'tulip', ...
    'arrParams', struct('Np', Np, 'pm', tulipPm), ...
    'tauDRO',    orientDeg, ...        % sheet-key rule
    'thrLock',   thrLock, 'tf0', tf0Warm, 'maxCpuSec', maxCpuSec, ...
    'resume',    false, ...            % single-entry: never mixes with a
    ...                                % sheet file on disk
    'logFile',   logFile);

thrust_ladder_library(outMat, optsEngine);
Q = load(outMat);

% krEntry (computed above from the REQUESTED ladder, before the solve)
% must still index thrustN's slot in Q.rungs -- assert it does, since
% thrust_ladder_library saves opts.rungs verbatim and in order:
assert(isequal(Q.rungs(krEntry), thrustN), 'gto_entry:rungIndex', ...
    'Q.rungs(%d)=%.4g does not match requested thrustN=%.4g', ...
    krEntry, Q.rungs(krEntry), thrustN);
entry = struct('ok', Q.OK(1,1,krEntry), 'z8', Q.Z8(:,1,1,krEntry), ...
    'tf_nd', Q.TF(1,1,krEntry), 'dV_kms', NaN, 'mf_kg', NaN, ...
    'coordinates', struct('orientDeg', orientDeg, 'Np', Np, ...
        'depFrac', depFrac, 'arrFrac', arrFrac, 'thrustN', thrustN));
gates = struct('msNormR', Q.RES(1,1,krEntry), 'flownKm', Q.FLYKM(1,1,krEntry), ...
    'acceptDz', Q.ACCDZ(1,1,krEntry), 'conjPass', NaN);

if entry.ok
    %% Mass fraction + delta-V (all-burn min-time), formulas from
    %% costate_common/catalog_schema.m ('m_final' / 'deltav_kms'):
    g0  = 9.80665*tStar^2/(1000*lStar);
    cnd = (ispS/tStar)*g0;
    Tnd = (thrustN/m0kg)*tStar^2/(lStar*1000);
    mfFrac = 1 - Tnd*entry.tf_nd/cnd;
    entry.mf_kg = m0kg*mfFrac;
    entry.dV_kms = cnd*log(1/mfFrac)*lStar/tStar;

    if conjTest
        [tD, rvD] = get_family_orbit('gto', struct('orientDeg', orientDeg));
        [tA, rvA] = get_family_orbit('tulip', struct('Np', Np, 'pm', tulipPm));
        rv0 = interp1(tD, rvD, mod(depFrac,1)*tD(end), 'spline');
        rvf = interp1(tA, rvA, mod(arrFrac,1)*tA(end), 'spline');
        z8  = entry.z8;
        K   = fdef(opts, 'K', 24);
        seed = seed_from_z8(z8, rv0(1:6), K, Tnd, cnd, muStar);
        [z, info] = ms_tfmin(rv0(1:6), rvf(1:6), seed, Tnd, cnd, muStar, ...
            struct('conjTest', true, 'wallSec', 60));
        dz = norm(z - z8);
        if info.converged && dz < 1e-6
            gates.conjPass = info.conj.pass;
        else
            fprintf('gto_entry: conjTest re-solve did not reproduce z8 (dz=%.1e) -- conjPass left NaN\n', dz);
        end
    end

    if writeSheet
        sheetMat = fullfile(here, 'results', 'catalog', ...
            sprintf('gto_o%03d_Np%d.mat', orientDeg, Np));
        merge_entry_into_sheet(sheetMat, orientDeg, Np, depFrac, arrFrac, ...
            Q, force);
    end
end
end

% ------------------------------------------------------------------------
function v = fdef(s, f, v0)
% FDEF  s.(f) if present else v0.  INPUTS: s; f; v0.  OUTPUTS: v.
if isfield(s, f), v = s.(f); else, v = v0; end
end

% ------------------------------------------------------------------------
function merge_entry_into_sheet(sheetMat, orientDeg, Np, depFrac, arrFrac, ...
    Qlocal, force)
% MERGE_ENTRY_INTO_SHEET  Fold a gto_entry() 1x1 result into the STANDARD
% nD=12/nA=6 campaign grid (nearest-phase snap on depFrac/arrFrac),
% creating the sheet file if it does not already exist.
%
% RUNG-GRID SOURCE (Important review finding, 2026-08-29): the merge
% validates/writes against the SHEET'S OWN rung set, not a literal passed
% in from gto_entry's ladder-planning step. gto_entry's internal rungsStd
% (the full 9-rung [15 12 10 7 5 3 2 1.5 1] set it uses to plan warm-start
% stepping stones down to an arbitrary requested thrustN) is NOT the same
% thing as a campaign sheet's grid: the shipped deliverable-7 v1 fleet's
% sheets all carry the 5-rung [15 12 10 7 5] truncation (run_gto_catalog's
% own default, matching its own rung-mismatch guard). Validating a merge
% against the 9-rung literal made this function THROW on every shipped
% sheet the first time it ran post-solve (isequal(size(S.OK),[nD nA nR])
% failed at nR=9 vs the on-disk nR=5). Fix: when sheetMat already exists,
% its own S.rungs is the grid of record (nR and the match set both come
% from it); only a BRAND-NEW sheet is stamped with a rung set at all, and
% that stamp is the campaign's actual 5-rung fleet default [15 12 10 7 5]
% -- not gto_entry's internal ladder-planning literal.
%
% SAFETY (review finding, 2026-08-26): a real campaign sweep's OK cells are
% never clobbered by a single-entry write -- if the snapped (iD,iA,krFull)
% is already OK, the write is REFUSED unless force is true (a warning
% names the cell either way). The per-cell attempt counter S.ATT belongs
% to thrust_ladder_library's own resume/maxAtt bookkeeping alone; this
% merge path never touches it (gto_entry's scratch solve is accounted for
% in its own scratch .mat, not the campaign sheet).
%
% INPUTS: sheetMat path; orientDeg; Np; depFrac; arrFrac; Qlocal (the
% thrust_ladder_library output loaded from gto_entry's 1x1 scratch run);
% force (logical, overwrite-already-OK permission).
% OUTPUTS: none (writes sheetMat).
nD = 12;  nA = 6;
iD = mod(round(depFrac*nD), nD) + 1;
iA = mod(round(arrFrac*nA), nA) + 1;
sheetDir = fileparts(sheetMat);
if ~isfolder(sheetDir), mkdir(sheetDir); end
if isfile(sheetMat)
    S = load(sheetMat);
    rungsGrid = S.rungs;
    nR = numel(rungsGrid);
    if ~isequal(size(S.OK), [nD nA nR])
        error('gto_entry:writeSheet', ...
            '%s has an incompatible grid size for the merge', sheetMat);
    end
else
    rungsGrid = [15 12 10 7 5];    % campaign fleet default -- matches
                                    % run_gto_catalog.m's own rungsIn
                                    % default and rung-mismatch guard
    nR = numel(rungsGrid);
    S = struct('TF',nan(nD,nA,nR), 'FLYKM',nan(nD,nA,nR), ...
        'ACCDZ',nan(nD,nA,nR), 'RES',nan(nD,nA,nR), 'WALL',nan(nD,nA,nR), ...
        'OK',false(nD,nA,nR), 'Z8',nan(8,nD,nA,nR), 'ATT',zeros(nD,nA), ...
        'rungs',rungsGrid, 'sD',(0:nD-1)/nD, 'sA',(0:nA-1)/nA, ...
        'meta', struct('depFamily','gto', ...
            'depParams',struct('orientDeg',orientDeg), ...
            'arrFamily','tulip','arrParams',struct('Np',Np,'pm',-1), ...
            'tauDRO',orientDeg,'source','gto_entry writeSheet'));
end
for kr = 1:numel(Qlocal.rungs)
    if ~Qlocal.OK(1,1,kr), continue, end
    krFull = find(abs(rungsGrid - Qlocal.rungs(kr)) < 1e-9, 1);
    if isempty(krFull)
        fprintf(['gto_entry: writeSheet -- requested thrust %.6g N matches no ' ...
                 'rung in %s''s grid [%s]; not written\n'], ...
                 Qlocal.rungs(kr), sheetMat, num2str(rungsGrid));
        continue
    end
    if S.OK(iD,iA,krFull) && ~force
        fprintf(['gto_entry: writeSheet -- cell (iD=%d,iA=%d,kr=%d) in %s is ' ...
                 'already OK from a real campaign sweep; refusing to overwrite ' ...
                 '(pass opts.force=true to override)\n'], iD, iA, krFull, sheetMat);
        continue
    end
    S.TF(iD,iA,krFull)    = Qlocal.TF(1,1,kr);
    S.FLYKM(iD,iA,krFull) = Qlocal.FLYKM(1,1,kr);
    S.ACCDZ(iD,iA,krFull) = Qlocal.ACCDZ(1,1,kr);
    S.RES(iD,iA,krFull)   = Qlocal.RES(1,1,kr);
    S.WALL(iD,iA,krFull)  = Qlocal.WALL(1,1,kr);
    S.OK(iD,iA,krFull)    = true;
    S.Z8(:,iD,iA,krFull)  = Qlocal.Z8(:,1,1,kr);
end
save(sheetMat, '-struct', 'S');
end
