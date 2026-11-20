function C = located_min_altitude(t, X, center, radiusKm, floorKm, opts)
%% Purpose:
%
%   THE LOCATED MINIMUM ALTITUDE of a flight over one body: the closest
%   approach BETWEEN samples, not the smallest sample. A perilune pass that
%   falls between two propagator samples is invisible to a sample minimum;
%   this finds it on the true dynamics. Lifted from the local LOCATED_MIN
%   of lunar_pnt_mission/lib/leg1_gates (2026-11-11), which stays as it is,
%   into the shared library so every gate (pair_transfer, the clearance
%   sweep) measures clearance the same way.
%
%   METHOD. The samples (t, X) are refined densely by a cubic Hermite
%   interpolant of position with the sampled velocities, nSub points per
%   sample interval; every place where the range rate (r - r_b).v changes
%   sign from negative to non-negative brackets a local minimum. Its root
%   is found by fzero on the range rate of the state PROPAGATED from the
%   left sample of the interval that holds the query time, or on the
%   Hermite interpolant when no propagator is given. If the propagated
%   range rate does not change sign on the dense bracket, the bracket
%   widens to [t(k-1), t(k+2)] (clamped to the ends); if still no sign
%   change, fminbnd on the distance there. hMinKm is the smallest altitude
%   over all located minima and the two end points (with no interior
%   minimum the nearer end point is returned), and never above the sample
%   minimum: if a propagated minimum disagrees with the samples upward,
%   the sample minimum and its time are reported.
%
%   A located minimum that cannot be computed (NaN state: a re-integration
%   that stopped before its end time, or a throw) FAILS the gate with a
%   reason; it is never dropped.
%
%  ASSUMPTIONS / NOTES:
%
% - Duplicate sample times (an integrator restart can emit one) are
%   dropped, keeping the first sample of each time, before refinement.
% - The built-in min-time propagator (opts.phys) integrates pumpkyn's
%   tfMinEoM without events from one sample to the query time: one sample
%   interval of a tfMinProp all-burn flight holds no switch. It must REACH
%   the query time (relative 1e-12 of the span); otherwise the state is
%   NaN and the minimum counts as not located.
% - Body radii used by the shared library: Moon 1737.4 km, Earth
%   6378.137 km (casadi_mintime_dro, thrust_ladder_library, leg1_gates).
%
%% Inputs:
%
%  t                        [M x 1]                 sample times (any unit
%                                                   the propagator uses; ND
%                                                   for CR3BP flights),
%                                                   non-decreasing
%
%  X                        [M x n >= 6]            states, one ROW per
%                                                   time; columns 1:3
%                                                   position, 4:6 velocity
%                                                   (n >= 14 for the
%                                                   min-time propagator:
%                                                   7 mass, 8:14 costates)
%
%  center                   [1 x 3]                 body position, same
%                                                   frame and unit as X
%                                                   (CR3BP: Moon [1-mu 0 0],
%                                                   Earth [-mu 0 0])
%
%  radiusKm                 double                  body radius, km
%
%  floorKm                  double                  altitude floor, km
%                                                   (-inf = report only)
%
%  opts                     struct (optional)       .lStar [389703.264829278]
%                                                   length unit, km per
%                                                   state unit (1 for km)
%                                                   .nSub [16] dense points
%                                                   per sample interval
%                                                   .prop [] handle
%                                                   (tA, xA [n x 1], tB) ->
%                                                   [n x 1], the true flow
%                                                   .phys [] struct .Tnd
%                                                   .cnd .muStar: builds the
%                                                   min-time propagator when
%                                                   .prop is empty
%                                                   neither: Hermite only
%
%% Outputs:
%
%  C                        struct                  .hMinKm located minimum
%                                                   altitude, km
%                                                   .tAtMin its time
%                                                   .pass hMinKm >= floorKm
%                                                   and nFailed == 0
%                                                   .reason '' or why not
%                                                   .sampleMinKm
%                                                   .sampleTAtMin .samplePass
%                                                   (the sample minimum)
%                                                   .floorKm .radiusKm
%                                                   .minima [nMin x 2]
%                                                   [t hKm] (NaN hKm = not
%                                                   located) .nBrackets
%                                                   .nFailed .method
%                                                   ('hermite' |
%                                                   'propagated')
%
%% References:
%   [1] lunar_pnt_mission/lib/leg1_gates.m, local LOCATED_MIN (the source).
%
%% Revision History:
%  M. Casey                                                   (c) 11/19/2026
%  M. Casey  review round 1: clamp to the sample minimum      (c) 11/19/2026
%  Copyright Coorbital Inc.
%% ------------------------ Begin Code Sequence ---------------------------

%% Self-demo: straight-line flyby, closest approach 1 km between samples
if nargin == 0
    tD = (0:9).';
    XD = [[-4.5 1 0] + tD*[1 0 0], repmat([1 0 0], 10, 1)];
    C  = located_min_altitude(tD, XD, [0 0 0], 0, 0.5, struct('lStar', 1));
    fprintf('located_min_altitude: sample minimum %.4f km at t %.2f, located %.6f km at t %.6f (%s)\n', ...
            C.sampleMinKm, C.sampleTAtMin, C.hMinKm, C.tAtMin, C.method);
    return
end
if nargin < 6, opts = struct(); end
     lStar = fdef(opts, 'lStar', 389703.264829278);
      nSub = fdef(opts, 'nSub',  16);
      prop = fdef(opts, 'prop',  []);
      phys = fdef(opts, 'phys',  []);

t = t(:);
if size(X, 1) ~= numel(t) || size(X, 2) < 6 || numel(t) < 2
    error('located_min_altitude:shape', ...
          'X must be [numel(t) x >= 6] with at least 2 samples; got [%d x %d] for %d times', ...
          size(X, 1), size(X, 2), numel(t));
end
if ~(isnumeric(center) && numel(center) == 3)
    error('located_min_altitude:center', 'center must be [1 x 3]');
end
center = center(:).';
if isempty(prop) && ~isempty(phys)
    if size(X, 2) < 14
        error('located_min_altitude:shape', ...
              'the min-time propagator needs X with >= 14 columns (state, mass, costates); got %d', size(X, 2));
    end
    prop = @(tA, xA, tB) fly14(tA, xA, tB, phys);
end
keep = [true; diff(t) > 0];                       % drop duplicate times
t = t(keep);  X = X(keep, :);
M = numel(t);

%% The sample minimum:
hS = lStar*sqrt(sum((X(:, 1:3) - center).^2, 2)) - radiusKm;
[hSmin, iS] = min(hS);

%% Dense Hermite refinement of the range rate:
s   = (0:nSub - 1)/nSub;
h00 = 2*s.^3 - 3*s.^2 + 1;   h10 = s.^3 - 2*s.^2 + s;
h01 = -2*s.^3 + 3*s.^2;      h11 = s.^3 - s.^2;
d00 = 6*s.^2 - 6*s;          d10 = 3*s.^2 - 4*s + 1;
d01 = -6*s.^2 + 6*s;         d11 = 3*s.^2 - 2*s;
 nD = (M - 1)*nSub + 1;
 tt = zeros(nD, 1);  rr = zeros(nD, 1);  kInt = zeros(nD, 1);
for k = 1:M - 1
    dt  = t(k+1) - t(k);
    r0  = X(k, 1:3) - center;    v0 = X(k, 4:6);
    r1  = X(k+1, 1:3) - center;  v1 = X(k+1, 4:6);
    pos = h00.'*r0 + (dt*h10.')*v0 + h01.'*r1 + (dt*h11.')*v1;
    vel = (d00.'*r0 + (dt*d10.')*v0 + d01.'*r1 + (dt*d11.')*v1)/dt;
    idx = (k - 1)*nSub + (1:nSub);
    tt(idx)   = t(k) + s.'*dt;
    rr(idx)   = sum(pos.*vel, 2);
    kInt(idx) = k;
end
tt(nD)   = t(M);
rr(nD)   = (X(M, 1:3) - center)*X(M, 4:6).';
kInt(nD) = M - 1;

%% Every negative-to-non-negative range-rate change is a local minimum:
  br = find(rr(1:end-1) < 0 & rr(2:end) >= 0);
mins = nan(numel(br), 2);
 why = cell(1, numel(br));
for q = 1:numel(br)
    k = kInt(br(q));
    a = tt(br(q));  b = tt(br(q) + 1);
    try
        if isempty(prop)
            g  = @(tq) rr_of(hermite_state(tq, t(k), t(k+1), X(k, :), X(k+1, :)), center);
            tR = fzero(g, [a b], optimset('TolX', 1e-15));
            rv = hermite_state(tR, t(k), t(k+1), X(k, :), X(k+1, :));
        else
            g  = @(tq) rr_of(piece_state(tq, t, X, prop), center);
            ga = g(a);  gb = g(b);
            if ~(ga < 0 && gb >= 0) && all(isfinite([ga gb]))
                % Hermite bracket off: widen to [t(k-1), t(k+2)], clamped
                a = t(max(k - 1, 1));  b = t(min(k + 2, M));  ga = g(a);  gb = g(b);
            end
            if ~all(isfinite([ga gb]))
                error('located_min_altitude:prop', ...
                      're-integration stopped early or non-finite near t = %.6g', a);
            end
            if sign(ga) ~= sign(gb)
                tR = fzero(g, [a b], optimset('TolX', 1e-15));
            else                                  % no sign change in the exact flight
                tR = fminbnd(@(tq) dist_of(piece_state(tq, t, X, prop), center), a, b, ...
                             optimset('TolX', 1e-15));
            end
            rv = piece_state(tR, t, X, prop);
        end
        hq = lStar*sqrt(sum((rv(1:3).' - center).^2)) - radiusKm;
        if ~isfinite(hq)
            error('located_min_altitude:prop', ...
                  're-integration stopped early or non-finite at t = %.6g', tR);
        end
        mins(q, :) = [tR, hq];
    catch err
        mins(q, :) = [a, NaN];
        why{q} = sprintf('minimum near t = %.6g not located: %s', a, err.message);
    end
end

%% The verdict over the located minima and the two end points:
cand = [mins; t(1), hS(1); t(M), hS(M)];
[hMin, kMin] = min(cand(:, 2));                   % min ignores NaN
tMin = cand(kMin, 1);
if hSmin < hMin                                   % never report above a sample
    hMin = hSmin;  tMin = t(iS);
end
nFailed = nnz(isnan(mins(:, 2)));
why  = why(~cellfun(@isempty, why));
pass = hMin >= floorKm && nFailed == 0;
if hMin < floorKm
    why = [{sprintf('located minimum %.4g km below the %.4g km floor', hMin, floorKm)}, why];
end
if isempty(prop), method = 'hermite'; else, method = 'propagated'; end
C = struct('hMinKm', hMin, 'tAtMin', tMin, 'pass', pass, ...
           'reason', strjoin(why, '; '), ...
           'sampleMinKm', hSmin, 'sampleTAtMin', t(iS), 'samplePass', hSmin >= floorKm, ...
           'floorKm', floorKm, 'radiusKm', radiusKm, 'minima', mins, ...
           'nBrackets', numel(br), 'nFailed', nFailed, 'method', method);
end

% ------------------------------------------------------------------------
function x = piece_state(tq, t, X, prop)
% PIECE_STATE  State at tq propagated from the left sample of the sample
%   interval that holds tq (the last interval for tq = t(end)).
%   INPUTS: tq; t [M x 1]; X [M x n]; prop handle.  OUTPUTS: x [n x 1].
k = find(t <= tq, 1, 'last');
if isempty(k), k = 1; end
k = min(k, numel(t) - 1);
x = prop(t(k), X(k, :).', tq);
end

% ------------------------------------------------------------------------
function g = rr_of(x, center)
% RR_OF  Range rate (r - r_b).v of a state.  INPUTS: x [n >= 6 x 1];
%   center [1 x 3].  OUTPUTS: g scalar.
g = (x(1:3).' - center)*x(4:6);
end

% ------------------------------------------------------------------------
function d = dist_of(x, center)
% DIST_OF  Distance to the body.  INPUTS: x [n >= 6 x 1]; center [1 x 3].
%   OUTPUTS: d scalar.
d = sqrt(sum((x(1:3).' - center).^2));
end

% ------------------------------------------------------------------------
function x = hermite_state(tq, ta, tb, xa, xb)
% HERMITE_STATE  Cubic Hermite position and its derivative at tq from the
%   end states of one interval.  INPUTS: tq, ta, tb; xa, xb [1 x n >= 6].
%   OUTPUTS: x [6 x 1] (position, velocity of the interpolant).
dt = tb - ta;  s = (tq - ta)/dt;
pos = (2*s^3 - 3*s^2 + 1)*xa(1:3) + dt*(s^3 - 2*s^2 + s)*xa(4:6) + ...
      (-2*s^3 + 3*s^2)*xb(1:3) + dt*(s^3 - s^2)*xb(4:6);
vel = ((6*s^2 - 6*s)*xa(1:3) + dt*(3*s^2 - 4*s + 1)*xa(4:6) + ...
       (-6*s^2 + 6*s)*xb(1:3) + dt*(3*s^2 - 2*s)*xb(4:6))/dt;
x = [pos, vel].';
end

% ------------------------------------------------------------------------
function x = fly14(tA, xA, tB, phys)
% FLY14  Min-time transfer state (state + costates, 14) at tB from (tA, xA)
%   on pumpkyn's tfMinEoM, no events. NaN if tB is not reached.
%   INPUTS: tA, tB ND; xA [>= 14 x 1]; phys (.Tnd .cnd .muStar).
%   OUTPUTS: x [14 x 1].
if tB == tA
    x = xA(1:14);
    x = x(:);
    return
end
w  = warning('off', 'MATLAB:ode113:IntegrationTolNotMet');
cw = onCleanup(@() warning(w));
[tO, Y] = ode113(@pumpkyn.cr3bp.tfMinEoM, [tA, 0.5*(tA + tB), tB], xA(1:14), ...
                 odeset('RelTol', 1e-12, 'AbsTol', 1e-13), phys.Tnd, phys.cnd, phys.muStar);
if abs(tO(end) - tB) > 1e-12*max(abs(tB - tA), eps) || ~all(isfinite(Y(end, :)))
    x = NaN(14, 1);
else
    x = Y(end, :).';
end
end

% ------------------------------------------------------------------------
function v = fdef(s, f, v0)
% FDEF  s.(f) if present and non-empty, else v0.  INPUTS: s; f; v0.
%   OUTPUTS: v.
if isfield(s, f) && ~isempty(s.(f)), v = s.(f); else, v = v0; end
end
