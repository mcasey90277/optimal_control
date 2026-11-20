function ok = test_located_min_altitude()
%% Purpose:
%
%   Unit test for located_min_altitude, the minimum altitude over a body
%   BETWEEN flight samples (range-rate sign change refined by fzero on the
%   true dynamics). Every case has an analytic answer:
%
%     1. straight-line flyby, closest approach 1 km at t = 4.5 between the
%        samples t = 4 and t = 5 (sample minimum sqrt(1.25) = 1.1180 km):
%        found by the Hermite route AND by an exact straight-line
%        propagator, both below the sample minimum;
%     2. Kepler hyperbola about the Moon (mu 4902.800066 km^3/s^2, rp 2000
%        km, e 1.5), samples every 600 s with periapsis 300 s from the
%        nearest sample: the propagated route finds the analytic periapsis
%        altitude 262.6 km (radius 1737.4 km) to 1e-6 km, the sample
%        minimum is far above it;
%     3. no interior minimum: a receding flight returns its first sample,
%        an approaching flight its last, with no bracket;
%     4. deliberate floor failure: a floor between the located and the
%        sample minimum fails the located gate while the sample gate passes;
%     5. deliberate propagator failure: a propagator returning NaN counts
%        the minimum as not located and fails the gate (never dropped);
%     6. duplicate sample times (an integrator restart) are tolerated;
%        and a propagator that disagrees with the samples UPWARD never
%        lifts the reported minimum above the sample minimum;
%     7. bad shapes refuse by name.
%
%   Units: opts.lStar = 1 makes the states kilometres, so every number
%   above is in km.
%
%% Inputs:
%
%  (none)
%
%% Outputs:
%
%  ok                       logical                 true iff every check
%                                                   passed
%
%% Revision History:
%  M. Casey                                                   (c) 11/19/2026
%  Copyright Coorbital Inc.
%% ------------------------ Begin Code Sequence ---------------------------
here = fileparts(mfilename('fullpath'));
oldPath = path;
restorePath = onCleanup(@() path(oldPath));
addpath(fileparts(here));                          % costate_common
ok = true;

%% 1. Straight-line flyby, minimum between samples:
t   = (0:9).';
p0  = [-4.5 1 0];   v = [1 0 0];
X   = [p0 + t*v, repmat(v, numel(t), 1)];
lin = @(tA, xA, tB) [xA(1:3) + xA(4:6)*(tB - tA); xA(4:6)];
oL  = struct('lStar', 1);
C1h = located_min_altitude(t, X, [0 0 0], 0, 0.5, oL);
oL.prop = lin;
C1p = located_min_altitude(t, X, [0 0 0], 0, 0.5, oL);
fprintf('  line: sample min %.6f km at t %.3f; located (hermite) %.12f km at t %.12f; located (prop) %.12f km at t %.12f\n', ...
        C1h.sampleMinKm, C1h.sampleTAtMin, C1h.hMinKm, C1h.tAtMin, C1p.hMinKm, C1p.tAtMin);
ok = check('line, hermite: located 1 km at t 4.5 (1e-9)', ...
           abs(C1h.hMinKm - 1) < 1e-9 && abs(C1h.tAtMin - 4.5) < 1e-9 && strcmp(C1h.method, 'hermite')) && ok;
ok = check('line, propagated: located 1 km at t 4.5 (1e-9)', ...
           abs(C1p.hMinKm - 1) < 1e-9 && abs(C1p.tAtMin - 4.5) < 1e-9 && strcmp(C1p.method, 'propagated')) && ok;
ok = check(sprintf('line: located %.4f km below the sample minimum %.4f km', C1p.hMinKm, C1p.sampleMinKm), ...
           C1p.hMinKm < C1p.sampleMinKm - 0.1 && abs(C1p.sampleMinKm - sqrt(1.25)) < 1e-12) && ok;
ok = check('line: one bracket, none failed, pass at floor 0.5', ...
           C1p.nBrackets == 1 && C1p.nFailed == 0 && C1p.pass && isempty(C1p.reason)) && ok;

%% 2. Kepler hyperbola about the Moon, propagated route:
muM = 4902.800066;  rp = 2000;  ecc = 1.5;  rM = 1737.4;
vp  = sqrt(muM*(1 + ecc)/rp);
xp  = [rp 0 0 0 vp 0];                              % periapsis state, t = 0
dt  = 600;
tS  = ((-5:5) + 0.5).'*dt;                          % periapsis 300 s from samples
kep = @(tt, x) [x(4:6); -muM*x(1:3)/sqrt(sum(x(1:3).^2))^3];
ode = odeset('RelTol', 1e-13, 'AbsTol', 1e-12);
XS  = zeros(numel(tS), 6);
for k = 1:numel(tS)
    [~, Z] = ode113(kep, [0, tS(k)/2, tS(k)], xp, ode);
    XS(k, :) = Z(end, :);
end
prop2 = @(tA, xA, tB) kep_prop(kep, ode, tA, xA, tB);
C2 = located_min_altitude(tS, XS, [0 0 0], rM, 200, struct('lStar', 1, 'prop', prop2));
C2h = located_min_altitude(tS, XS, [0 0 0], rM, 200, struct('lStar', 1));
fprintf('  hyperbola: analytic %.6f km at t 0; located (prop) %.9f km at t %.6f s; hermite %.4f km; sample min %.4f km\n', ...
        rp - rM, C2.hMinKm, C2.tAtMin, C2h.hMinKm, C2.sampleMinKm);
ok = check('hyperbola: located periapsis altitude 262.6 km (1e-6 km) at t 0 (1e-3 s)', ...
           abs(C2.hMinKm - (rp - rM)) < 1e-6 && abs(C2.tAtMin) < 1e-3) && ok;
ok = check(sprintf('hyperbola: sample minimum %.1f km is above the located %.1f km', C2.sampleMinKm, C2.hMinKm), ...
           C2.sampleMinKm > C2.hMinKm + 50) && ok;               % measured 343.3 vs 262.6
% Hermite-only lands BELOW the true periapsis here (conservative on this
% case), but it is not the true minimum: the propagator is what finds it
ok = check(sprintf('hyperbola: Hermite-only %.4f km is off the true value by %+.2f km (needs the propagator)', ...
           C2h.hMinKm, C2h.hMinKm - (rp - rM)), abs(C2h.hMinKm - (rp - rM)) > 0.1) && ok;

%% 3. No interior minimum (end points returned):
tR  = (0:9).';
XR  = [[1 1 0] + tR*[1 0 0], repmat([1 0 0], 10, 1)];       % receding
C3a = located_min_altitude(tR, XR, [0 0 0], 0, 0, struct('lStar', 1, 'prop', lin));
XA  = [[-12 1 0] + tR*[1 0 0], repmat([1 0 0], 10, 1)];     % approaching
C3b = located_min_altitude(tR, XA, [0 0 0], 0, 0, struct('lStar', 1, 'prop', lin));
fprintf('  receding: %.6f km at t %g (%d brackets); approaching: %.6f km at t %g (%d brackets)\n', ...
        C3a.hMinKm, C3a.tAtMin, C3a.nBrackets, C3b.hMinKm, C3b.tAtMin, C3b.nBrackets);
ok = check('receding: first sample returned, no bracket', ...
           C3a.nBrackets == 0 && C3a.tAtMin == 0 && abs(C3a.hMinKm - sqrt(2)) < 1e-12) && ok;
ok = check('approaching: last sample returned, no bracket', ...
           C3b.nBrackets == 0 && C3b.tAtMin == 9 && abs(C3b.hMinKm - sqrt(10)) < 1e-12) && ok;

%% 4. Deliberate floor failure (between located and sample minimum):
C4 = located_min_altitude(t, X, [0 0 0], 0, 1.05, struct('lStar', 1, 'prop', lin));
ok = check(sprintf('floor 1.05 km: located fails (%s), sample passes', C4.reason), ...
           ~C4.pass && C4.samplePass && contains(C4.reason, 'below')) && ok;

%% 5. Deliberate propagator failure:
bad = @(tA, xA, tB) NaN(6, 1);
C5 = located_min_altitude(t, X, [0 0 0], 0, 0, struct('lStar', 1, 'prop', bad));
ok = check(sprintf('NaN propagator: nFailed %d, pass false (%s)', C5.nFailed, C5.reason), ...
           C5.nFailed == 1 && ~C5.pass && contains(C5.reason, 'not located')) && ok;

%% 6. Duplicate sample times:
tD = [t(1:5); t(5); t(6:end)];
XD = [X(1:5, :); X(5, :); X(6:end, :)];
C6 = located_min_altitude(tD, XD, [0 0 0], 0, 0.5, struct('lStar', 1, 'prop', lin));
ok = check('duplicate sample time tolerated, same answer', ...
           abs(C6.hMinKm - 1) < 1e-9 && C6.nFailed == 0) && ok;

%% 6b. A propagator that disagrees UPWARD never lifts the answer above the
%%     sample minimum (shifted 5 km off the samples: located would be 6 km):
up  = @(tA, xA, tB) [xA(1:3) + xA(4:6)*(tB - tA) + [0; 5; 0]; xA(4:6)];
C6b = located_min_altitude(t, X, [0 0 0], 0, 0.5, struct('lStar', 1, 'prop', up));
fprintf('  upward propagator: minima %s km, reported %.6f km at t %g (sample min %.6f at t %g)\n', ...
        mat2str(C6b.minima(:, 2).', 6), C6b.hMinKm, C6b.tAtMin, C6b.sampleMinKm, C6b.sampleTAtMin);
ok = check('upward propagator: reported = sample minimum 1.1180 km at its sample time 4', ...
           abs(C6b.hMinKm - sqrt(1.25)) < 1e-12 && C6b.tAtMin == 4 && C6b.minima(1, 2) > 5) && ok;

%% 7. Shapes refuse by name:
try
    located_min_altitude(t, X(:, 1:5), [0 0 0], 0, 0, struct('lStar', 1));
    ok = check('X with 5 columns refuses', false) && ok;
catch err
    ok = check(['X with 5 columns refuses: ' err.identifier], ...
               strcmp(err.identifier, 'located_min_altitude:shape')) && ok;
end

if ok, fprintf('test_located_min_altitude: ALL PASS\n');
else,  fprintf('test_located_min_altitude: FAILURES\n'); end
end

% ------------------------------------------------------------------------
function x = kep_prop(kep, ode, tA, xA, tB)
% KEP_PROP  Two-body state at tB from (tA, xA).  INPUTS: kep field handle;
%   ode odeset; tA, tB s; xA [6 x 1] km, km/s.  OUTPUTS: x [6 x 1].
if tB == tA
    x = xA(:);
    return
end
[~, Z] = ode113(kep, [tA, 0.5*(tA + tB), tB], xA(:), ode);
x = Z(end, :).';
end

% ------------------------------------------------------------------------
function ok = check(name, cond)
% CHECK  One gate line.  INPUTS: name; cond.  OUTPUTS: ok.
ok = logical(cond);
if ok, tag = 'PASS'; else, tag = 'FAIL'; end
fprintf('  [%s] %s\n', tag, name);
end
