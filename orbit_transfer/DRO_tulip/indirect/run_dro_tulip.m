function T = run_dro_tulip(sD, sA, opts)
%% Purpose:
%
%   FRONT DOOR for ONE DRO -> tulip minimum-time transfer. Give it a
%   departure phase and an arrival phase (fractions of the two orbits'
%   periods) and it returns a CERTIFIED transfer or a named reason why not:
%
%     >> T = run_dro_tulip;                    % the anchor, 70 mN
%     >> T = run_dro_tulip(0, 0.1587);         % another arrival phase
%     >> T = run_dro_tulip(0.75, 0.0754, struct('movie', true));
%
%   Three routes, tried in order:
%     LIBRARY   the pair is already certified on disk -> re-certify and
%               report (seconds).
%     WALK      the arrival phase differs -> pseudo-arclength continuation
%               in arrival phase from the nearest library solution, then a
%               departure walk. Minutes to hours; the arrival direction is
%               where the folds are.
%   Whatever the route, the answer goes through the SAME gate stack
%   (certify_root): normal-chart polish, flown position and velocity,
%   pumpkyn tfMin witness, conjugate test, hypothesis gates. A number is
%   never reported without its verdict.
%
%% Inputs:
%
%  sD                       double                  departure phase, fraction
%                                                   of the DRO period [0]
%  sA                       double                  arrival phase, fraction
%                                                   of the tulip period
%                                                   [0.0754, the anchor]
%  opts                     struct (optional)
%   .thrustN [0.070] .ispS [900] .m0kg [150]
%   .quiet [false] .movie [false] .wallSec [1800] walk budget
%   .tolPhase [1e-6] library match tolerance
%
%% Outputs:
%
%  T                        struct                  certify_root output plus
%                                                   .source ('library' |
%                                                   'walk-arrival' |
%                                                   'walk-both') .sD .sA
%                                                   .walkSec .movieStem
%
%% Revision History:
%  M. Casey                                                   (c) 09/09/2026
%  Copyright Coorbital Inc.
%% ------------------------ Begin Code Sequence ---------------------------

if nargin < 1 || isempty(sD), sD = 0; end
if nargin < 2 || isempty(sA), sA = 0.0754; end
if nargin < 3, opts = struct(); end
d = @(f,v) fieldd(opts, f, v);
quiet = d('quiet', false);  tolPhase = d('tolPhase', 1e-6);
wallSec = d('wallSec', 1800);
here = fileparts(mfilename('fullpath'));
addpath(fullfile(fileparts(here), '..', 'costate_common'));
say = @(varargin) sayif(quiet, varargin{:});
sD = mod(sD, 1);  sA = mod(sA, 1);
tStar = 382981.289129055;
t0 = tic;

say('DRO -> tulip min-time: departure phase %.4f, arrival phase %.4f, %.1f mN, Isp %g s', ...
    sD, sA, d('thrustN', 0.070)*1000, d('ispS', 900));

% The closure setup re-solves the DEFAULT anchor, so it must be asked for
% the anchor's OWN departure phase (0), not the target's -- at another sD
% that root does not exist and the setup's own guard rejects it. B is used
% here only for the closures (stateD, stateA, Tnd, cnd, mu); the target
% states are taken from those closures explicitly.
pool = capped_pool();          % the hard-timeout fence for external calls
setupOpts = opts;  setupOpts.sD = 0;  setupOpts.sA0 = 0.0754;
[B, anc] = arclength_arrival('setup', setupOpts);
lib = dro_tulip_library();

% ---- 1. the pair itself in the library ---------------------------------
hit = find(abs(wrapDiff([lib.sD], sD)) < tolPhase & abs(wrapDiff([lib.sA], sA)) < tolPhase, 1);
if ~isempty(hit)
    say('  found in the certified library (%s); re-certifying', lib(hit).src);
    seed = seed_of(lib(hit), B, sD);
    T = certify_root(seed, B.stateD(sD), B.stateA(sA), B, ...
                     struct('sA', sA, 'sD', sD, 'wallSec', wallSec, 'pool', pool));
    T.source = 'library';
else
    % ---- 2. walk: arrival phase first (folds live there), then departure
    near = nearest_point(lib, sD, sA);
    say('  not in the library; walking from (%.4f, %.4f) [%s]', near.sD, near.sA, near.src);
    T = walk_to(near, sD, sA, B, anc, opts, say, tolPhase, wallSec, pool);
end

% IDENTITY COMES FROM THE CERTIFICATE. Overwriting it with the request is
% how one phase gets certified and another reported. Assert instead.
T.walkSec = toc(t0);
if T.ok
    assert(abs(wrapDiff(T.sD, sD)) < 10*tolPhase && abs(wrapDiff(T.sA, sA)) < 10*tolPhase, ...
        ['certified (%.6f, %.6f) but (%.6f, %.6f) was requested -- refusing to ' ...
         'relabel'], T.sD, T.sA, sD, sA);
else
    T.sD = sD;  T.sA = sA;          % a refusal may name what was asked for
end
say('  %s', T.reason);
if T.ok
    say('  t_f = %.4f d   dV = %.4f km/s   fuel = %.2f kg', T.tfDays, T.dvKms, T.propellantKg);
    say('  flown miss %.3f km / %.3f m/s | tfMin witness |dz| = %.1e | conjugate %s', ...
        T.flyKm, T.flyVms, T.dz, tern(T.conj == 1, 'PASS', 'FAIL'));
    say('  gates: min|lam_v| = %.3e, min Q = %.3e, dim S = %d', T.g.minLamV, T.g.minQmt, T.g.dimS);
else
    say('  NOT certified: %s', T.reason);
end
say('  (%.1f s, source: %s)', T.walkSec, T.source);

T.movieStem = '';
if d('movie', false) && T.ok
    stem = fullfile(here, '..', 'direct', 'results', ...
                    sprintf('dro_tulip_sD%.4f_sA%.4f', sD, sA));
    T.movieStem = movie_70mN(stem, struct('z', T.z, 'rv0', B.stateD(sD), 'rvf', B.stateA(sA), ...
                                          'Tnd', B.Tnd, 'cnd', B.cnd, 'mu', B.mu));
end
end

% ------------------------------------------------------------------------
function T = walk_to(near, sD, sA, B, anc, opts, say, tolPhase, wallSec, pool)
% WALK_TO  Continuation from a library point to (sD, sA).  INPUTS: see
% caller.  OUTPUTS: T (certify_root output + .source).
T = struct('ok', false, 'reason', 'walk not attempted', 'source', 'walk', ...
           'z', nan(8,1), 'Y', [], 'tfDays', NaN, 'dvKms', NaN, ...
           'propellantKg', NaN, 'finalMassKg', NaN, ...
           'flyKm', NaN, 'flyVms', NaN, 'dz', NaN, 'conj', -1, 'g', [], 'rho', NaN, ...
           'normR', NaN, 'wallSec', NaN);
so = opts;  so.sD = near.sD;  so.sA0 = near.sA;  so.anchorMat = near.file;
so.root = struct('z', near.z, 'it', struct('Y', near.Y));   % layout-independent
if isfield(near, 'K') && ~isempty(near.K), so.K = near.K; end
[Bn, ancN] = arclength_arrival('setup', so);

if abs(wrapDiff(near.sA, sA)) > tolPhase
    dirn = sign(wrapDiff(near.sA, sA));
    lvl = near.sA + wrapDiff(near.sA, sA);              % unwrapped target
    say('  arrival walk %.4f -> %.4f (direction %+d)', near.sA, sA, dirn);
    ao = so;  ao.direction = dirn;  ao.levels = lvl;  ao.nStep = 2000;
    ao.sAStop = lvl + dirn*0.01;  ao.deadlineSec = wallSec;
    A = arclength_arrival(ancN, ao);
    hits = find([A.crossings.converged]);
    if isempty(hits)
        T.reason = sprintf('arrival walk reached %.4f but never crossed %.4f (%s)', A.q(end), lvl, A.stop);
        T.source = 'walk-arrival';  return
    end
    % TRY EVERY converged crossing, fastest first. The last one is not
    % necessarily the usable one: it can fail certification while an earlier
    % crossing at the same phase certifies. (Astra chain review 2026-09-10.)
    tfs = arrayfun(@(c) c.p(ancN.ctf), A.crossings(hits));
    [~, ord] = sort(tfs);
    T.source = 'walk-arrival';
    for kk = ord(:)'
        Tk = certify_crossing(A.crossings(hits(kk)).p, sA, Bn, ancN, ...
                              struct('wallSec', wallSec, 'pool', pool));
        if Tk.ok, T = Tk;  T.source = 'walk-arrival';  break, end
        T.reason = Tk.reason;  T.tfDays = Tk.tfDays;
    end
    if ~T.ok, return, end
else
    seedN = seed_of(near, Bn, near.sD);
    T = certify_root(seedN, Bn.stateD(near.sD), Bn.stateA(sA), Bn, ...
                     struct('sA', sA, 'sD', near.sD, 'wallSec', wallSec));
    T.source = 'library-seed';
    if ~T.ok, T.reason = ['seed point did not re-certify: ' T.reason];  return, end
end

if abs(wrapDiff(near.sD, sD)) > tolPhase
    % walk to the EXACT requested phase; never infer a step count by
    % rounding and then label the result with what was asked for
    say('  departure walk %.4f -> %.4f', near.sD, sD);
    C0 = T;  C0.sA = sA;
    R = rib_from_crossing(C0, Bn, ancN, struct('targets', wrapDiff(near.sD, sD), ...
        'nD', 24, 'wallSec', wallSec, 'logFile', ''));
    if isempty(R.pts)
        T.ok = false;  T.reason = ['departure walk: ' R.stop];  T.source = 'walk-both';  return
    end
    T = R.pts(end);  T.source = 'walk-both';
end
end

function seed = seed_of(P, B, sD)
% SEED_OF  ms seed from a library point.  INPUTS: P; B; sD.  OUTPUTS: seed.
K = size(P.Y, 2);
seed = struct('tf', P.z(8), 'tGrid', linspace(0, P.z(8), K+1), 'Y', [P.Y, P.Y(:,end)]);
rv0 = B.stateD(sD);
seed.Y(1:7, 1) = [rv0(1:6); 1];  seed.Y(8:14, 1) = P.z(1:7);
end

function P = nearest_point(lib, sD, sA)
% NEAREST_POINT  Library point closest in phase, ARRIVAL weighted heavily --
% the arrival direction carries the folds and costs the most to walk.
% INPUTS: lib; sD; sA.  OUTPUTS: P.
w = arrayfun(@(p) 10*abs(wrapDiff(p.sA, sA)) + abs(wrapDiff(p.sD, sD)), lib);
[~, k] = min(w);  P = lib(k);
end

function d_ = wrapDiff(a, b)
% WRAPDIFF  Signed difference b - a wrapped to [-0.5, 0.5).  INPUTS: a; b.
% OUTPUTS: d_.
d_ = mod(b - a + 0.5, 1) - 0.5;
end

function s = tern(c, a, b)
% TERN  Inline conditional.  INPUTS: c; a; b.  OUTPUTS: s.
if c, s = a; else, s = b; end
end

function sayif(quiet, varargin)
% SAYIF  printf unless quiet.  INPUTS: quiet; varargin.
if ~quiet, fprintf([varargin{1} '\n'], varargin{2:end}); end
end

function v = fieldd(s, f, d_)
% FIELDD  Field with default.  INPUTS: s; f; d_.  OUTPUTS: v.
if isfield(s, f) && ~isempty(s.(f)), v = s.(f); else, v = d_; end
end
