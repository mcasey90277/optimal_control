function [entry, gates, P] = gto_dro_entry(orientDeg, tauDRO, depFrac, arrFrac, thrustN, opts)
%% Purpose:
%
%   THE SINGLE-ENTRY DRIVER for GTO -> DRO minimum-time transfers: solve
%   ONE (departure phase x arrival phase) cell at a requested thrust level
%   through the SAME shared engine every costate catalog uses
%   (costate_common/thrust_ladder_library, family-agnostic via
%   get_family_orbit) -- never a bespoke solver. A copy of
%   GTO_tulip/catalog/gto_entry.m (2026-08-26) with the arrival family
%   switched from the tulip to a DRO of period tauDRO; the sheet-writing
%   path is dropped (stage 2's chain script owns sheets).
%
%   LADDER DIRECTION: solves from 15 N DOWN to thrustN; higher standard
%   rungs are warm-start stepping stones (never returned) and the entry
%   AT thrustN is returned. A deep rung is never cold-solved: the pilot
%   of the GTO -> tulip campaign measured a cold deep start as exactly
%   the non-convergent case.
%
%   WARM RECIPE (inherited from run_gto_catalog's pass A, 2026-08-26
%   diagnostic A; a DRO-specific recalibration is stage 2's pilot):
%   thrLock = true, tf0 = 0.30 ND (seeds the TOP rung only),
%   maxCpuSec = 600, maxIter = 6000. The engine's cold defaults (free
%   throttle, tf0 = 4.0) measured 67/72 GTO cells unreachable; this driver
%   never falls back to them.
%
%   Self-demo (nargin == 0): orient 0 deg, tau 1, PERIGEE departure
%   (phase 0), arrival phase 0, 15 N -- one solve. Perigee, not apogee:
%   the GTO campaign measured an apogee cold start at 15 N landing in a
%   bad local minimum (11.95e6 km flown miss).
%
%% Inputs:
%
%  orientDeg                double                  GTO locus orientation,
%                                                   deg: Earth->Moon line to
%                                                   the perigee direction,
%                                                   rotation sense
%
%  tauDRO                   double                  DRO period, ND (0.5, 1,
%                                                   2, 3 are the catalog
%                                                   periods; 1 = 4.43 d)
%
%  depFrac                  double                  departure phase: TIME
%                                                   fraction since perigee
%                                                   over one Kepler period
%                                                   (0 perigee, 0.5 apogee)
%
%  arrFrac                  double                  arrival phase: time
%                                                   fraction along the DRO
%                                                   from its +x crossing
%
%  thrustN                  double                  thrust, N, in [1, 15]:
%                                                   the BOTTOM (returned)
%                                                   rung
%
%  opts                     struct (optional)       .m0kg [150], .ispS [1710],
%                                                   .thrLock [true], .tf0
%                                                   [0.30], .maxCpuSec [600],
%                                                   .maxIter [6000], .N [400],
%                                                   .floorKm [500] lunar,
%                                                   .floorEarthKm [200] Earth
%                                                   altitude floor, km (engine
%                                                   minAltEarthKm; [] = none),
%                                                   .gateKm [100], .accTol
%                                                   [1e-6], .gtoParams
%                                                   (struct, .sma_km .ecc;
%                                                   default = get_family_
%                                                   orbit's GTO), .conjTest
%                                                   [false] re-solve ms_tfmin
%                                                   at the accepted z8 with
%                                                   the conjugate test, .K
%                                                   [24] segments for it,
%                                                   .outMat scratch path
%                                                   [results/scratch/...],
%                                                   .keepScratch [false] (the
%                                                   .log is kept whenever the
%                                                   solve fails),
%                                                   .planOnly [false] return
%                                                   the plan, no solve
%
%% Outputs:
%
%  entry                    struct                  .ok, .z8 [8 x 1], .tf_nd,
%                                                   .tf_d, .mf_kg, .dV_kms,
%                                                   .thrustN, .rungsWalked
%                                                   [1 x nR], .tfPerRung,
%                                                   .okPerRung, .rv0, .rvf
%                                                   [1 x 6], .depFrac,
%                                                   .arrFrac, .orientDeg,
%                                                   .tauDRO, .m0kg, .ispS,
%                                                   .engineOpts (the options
%                                                   handed to the engine),
%                                                   .scratchMat (default
%                                                   scratch path; carries
%                                                   orient, tau, both phases
%                                                   and thrust so parallel
%                                                   cells never collide)
%
%  gates                    struct                  .msNormR, .flownKm,
%                                                   .accDz, .wallS, .conjPass
%                                                   (NaN unless .conjTest
%                                                   reproduced z8),
%                                                   .conjReason (char: 'not run'
%                                                   | 'reproduced' | 'not
%                                                   reproduced, |dz| = abs and
%                                                   rel to max(1,max|z8|)' |
%                                                   'error: <msg>'),
%                                                   .ratioScreen [1 x nR-1],
%                                                   .ratioPass
%
%  P                        struct                  thrust_ladder_library's
%                                                   output for the cell
%                                                   (every rung); empty when
%                                                   .planOnly
%
%% Revision History:
%  M. Casey                                                   (c) 11/09/2026
%  M. Casey  floorEarthKm option (default 200 km)             (c) 11/09/2026
%  Copyright Coorbital Inc.
%% ------------------------ Begin Code Sequence ---------------------------

if nargin == 0
    orientDeg = 0;  tauDRO = 1;  depFrac = 0;  arrFrac = 0;  thrustN = 15;
end
if nargin < 6, opts = struct(); end
m0kg       = fdef(opts, 'm0kg',       150);
ispS       = fdef(opts, 'ispS',       1710);
thrLock    = fdef(opts, 'thrLock',    true);    % warm recipe
tf0Warm    = fdef(opts, 'tf0',        0.30);    % warm recipe (top rung only)
maxCpuSec  = fdef(opts, 'maxCpuSec',  600);     % warm recipe
maxIterOpt = fdef(opts, 'maxIter',    6000);    % warm recipe (engine default 3000 fell short)
N          = fdef(opts, 'N',          400);
floorKm    = fdef(opts, 'floorKm',    500);     % LUNAR floor
floorEarthKm = fdef(opts, 'floorEarthKm', 200);  % EARTH floor (engine option added 2026-11-09)
gateKm     = fdef(opts, 'gateKm',     100);
accTol     = fdef(opts, 'accTol',     1e-6);
gtoParams  = fdef(opts, 'gtoParams',  struct());
conjTest   = fdef(opts, 'conjTest',   false);
Kseg       = fdef(opts, 'K',          24);
keepScratch = fdef(opts, 'keepScratch', false);
planOnly   = fdef(opts, 'planOnly',   false);

if ~(isscalar(thrustN) && isreal(thrustN) && thrustN >= 1 && thrustN <= 15)
    error('gto_dro_entry:thrust', ...
          ['thrustN = %g N is outside [1, 15]: above 15 N there is no standard rung to anchor on; ' ...
           'below 1 N is the many-rev regime that ships junction states, not bare z8 (spec section 9)'], thrustN);
end

here   = fileparts(fileparts(mfilename('fullpath')));      % GTO_DRO/
run(fullfile(here, 'setup_gto_dro_paths.m'));

%% Ladder from 15 N DOWN to thrustN (stepping stones above, thrustN last)
rungsStd = [15 12 10 7 5 3 2 1.5 1];
rungs    = [rungsStd(rungsStd > thrustN), thrustN];
krEntry  = numel(rungs);

muStar = 0.012150585609624;
lStar  = 389703.264829278;
tStar  = 382981.289129055;

depParams = struct('orientDeg', orientDeg);
for f = {'sma_km', 'ecc'}
    if isfield(gtoParams, f{1}), depParams.(f{1}) = gtoParams.(f{1}); end
end

scratchDir = fullfile(here, 'results', 'scratch');
outMat  = fdef(opts, 'outMat', fullfile(scratchDir, sprintf('gto_dro_o%g_tau%g_sD%.4f_sA%.4f_%gN.mat', orientDeg, tauDRO, depFrac, arrFrac, thrustN)));
logFile = [outMat(1:end-4) '.log'];

optsEngine = struct( ...
    'rungs',     rungs, ...
    'ispS',      ispS, ...
    'm0kg',      m0kg, ...
    'N',         N, ...
    'floorKm',   floorKm, ...
    'floorEarthKm', floorEarthKm, ...
    'maxIter',   maxIterOpt, ...
    'gateKm',    gateKm, ...
    'accTol',    accTol, ...
    'nD',        1, ...
    'nA',        1, ...
    'sD0',       depFrac, ...
    'sA0',       arrFrac, ...
    'muStar',    muStar, 'lStar', lStar, 'tStar', tStar, ...
    'depFamily', 'gto', ...
    'depParams', depParams, ...
    'arrFamily', 'dro', ...
    'arrParams', struct('tau', tauDRO), ...
    'tauDRO',    orientDeg, ...        % SHEET KEY RULE (GTO catalog): the key is the orientation
    'thrLock',   thrLock, 'tf0', tf0Warm, 'maxCpuSec', maxCpuSec, ...
    'resume',    false, ...
    'logFile',   logFile);

entry = struct('ok', false, 'z8', nan(8, 1), 'tf_nd', NaN, 'tf_d', NaN, 'mf_kg', NaN, 'dV_kms', NaN, ...
    'thrustN', thrustN, 'rungsWalked', rungs, 'tfPerRung', nan(1, numel(rungs)), 'okPerRung', false(1, numel(rungs)), ...
    'rv0', nan(1, 6), 'rvf', nan(1, 6), 'depFrac', depFrac, 'arrFrac', arrFrac, ...
    'orientDeg', orientDeg, 'tauDRO', tauDRO, 'm0kg', m0kg, 'ispS', ispS, 'engineOpts', optsEngine, 'scratchMat', outMat);
gates = struct('msNormR', NaN, 'flownKm', NaN, 'accDz', NaN, 'wallS', NaN, 'conjPass', NaN, 'conjReason', 'not run', ...
    'ratioScreen', nan(1, max(numel(rungs) - 1, 0)), 'ratioPass', true);
P = struct();
if planOnly, return; end

%% Endpoints by THE shared endpoint rule (costate_common/phase_state)
[tD, rvD] = get_family_orbit('gto', depParams);
[tA, rvA] = get_family_orbit('dro', struct('tau', tauDRO));
depAt = phase_state(tD, rvD);   arrAt = phase_state(tA, rvA);
entry.rv0 = depAt(depFrac).';
entry.rvf = arrAt(arrFrac).';

%% Solve
if ~isfolder(fileparts(outMat)), mkdir(fileparts(outMat)); end
if isfile(outMat),  delete(outMat);  end
if isfile(logFile), delete(logFile); end
cleanMat = onCleanup(@() dropFile(outMat, keepScratch));   % .mat goes on every path; the .log only on success
thrust_ladder_library(outMat, optsEngine);
P = load(outMat);
assert(isequal(P.rungs(krEntry), thrustN), 'gto_dro_entry:rungIndex', ...
    'P.rungs(%d) = %.4g does not match the requested thrustN = %.4g', krEntry, P.rungs(krEntry), thrustN);

entry.ok        = logical(P.OK(1, 1, krEntry));
entry.z8        = P.Z8(:, 1, 1, krEntry);
entry.tf_nd     = P.TF(1, 1, krEntry);
entry.tf_d      = entry.tf_nd*tStar/86400;
entry.tfPerRung = squeeze(P.TF(1, 1, :)).';
entry.okPerRung = logical(squeeze(P.OK(1, 1, :)).');
gates.msNormR   = P.RES(1, 1, krEntry);
gates.flownKm   = P.FLYKM(1, 1, krEntry);
gates.accDz     = P.ACCDZ(1, 1, krEntry);
gates.wallS     = P.WALL(1, 1, krEntry);

% rung-to-rung screen on STEPS only (the top rung's guess has no branch behind it)
if numel(rungs) > 1
    gates.ratioScreen = entry.tfPerRung(2:end)./entry.tfPerRung(1:end-1);
    gates.ratioPass   = all(gates.ratioScreen > 1/3 & gates.ratioScreen < 3 | ~entry.okPerRung(2:end));
end

if entry.ok
    ndp = nd_propulsion(thrustN, ispS, m0kg, lStar, tStar);
    cnd = ndp.cnd;   Tnd = ndp.Tnd;
    mfFrac       = 1 - Tnd*entry.tf_nd/cnd;                 % all-burn: m(t) = 1 - T t / c
    entry.mf_kg  = m0kg*mfFrac;
    entry.dV_kms = cnd*log(1/mfFrac)*lStar/tStar;

    if conjTest
        try
            seed = seed_from_z8(entry.z8, entry.rv0(1:6), Kseg, Tnd, cnd, muStar);
            [z, info] = ms_tfmin(entry.rv0(1:6), entry.rvf(1:6), seed, Tnd, cnd, muStar, ...
                                 struct('conjTest', true, 'wallSec', 120));
            dzAbs = max(abs(z - entry.z8));
            dzRel = dzAbs/max(1, max(abs(entry.z8)));
            if info.converged && dzAbs < 1e-6
                gates.conjPass   = double(info.conj.pass);
                gates.conjReason = 'reproduced';
            else
                gates.conjReason = sprintf('not reproduced, |dz| = %.1e abs, %.1e rel', dzAbs, dzRel);
                fprintf('gto_dro_entry: conjTest %s -- conjPass left NaN\n', gates.conjReason);
            end
        catch cerr
            gates.conjReason = sprintf('error: %s', cerr.message);
            fprintf('gto_dro_entry: conjTest %s -- conjPass left NaN\n', gates.conjReason);
        end
    end
end

if entry.ok && ~keepScratch
    dropFile(logFile, false);
elseif ~entry.ok
    fprintf('gto_dro_entry: solve failed -- engine log kept at %s\n', logFile);
end

if nargin == 0 || nargout == 0
    fprintf('gto_dro_entry: orient %g deg, DRO tau %g, phases (%.3f, %.3f), %g N -> ok = %d\n', ...
            orientDeg, tauDRO, depFrac, arrFrac, thrustN, entry.ok);
    fprintf('   t_f = %.4f ND = %.3f d, m_f = %.2f kg, dV = %.3f km/s; ms |R| %.1e, flown %.3f km, accept |dz| %.1e\n', ...
            entry.tf_nd, entry.tf_d, entry.mf_kg, entry.dV_kms, gates.msNormR, gates.flownKm, gates.accDz);
end
end

% ------------------------------------------------------------------------
function v = fdef(s, f, v0)
% FDEF  s.(f) if present else v0.  INPUTS: s; f; v0.  OUTPUTS: v.
if isfield(s, f), v = s.(f); else, v = v0; end
end

% ------------------------------------------------------------------------
function dropFile(f, keep)
% DROPFILE  Delete file f unless keep.  INPUTS: f (char); keep (logical).  OUTPUTS: none.
if ~keep && isfile(f), delete(f); end
end
