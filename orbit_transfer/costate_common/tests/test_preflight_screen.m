function ok = test_preflight_screen()
%% Purpose:
%
%   Tests preflight_screen -- the node-only sanity check run before any
%   integrator touches a direct solution -- on synthetic node sets whose
%   lunar altitudes are known exactly:
%     1. minAltKm is the node minimum of |r - r_Moon| lStar - 1737.4.
%     2. ALTITUDE trips below HALF the floor, not above it.
%     3. TF band [0.3, 3] x seedTf: outside trips 'tf', inside passes.
%     4. seedTf = [] skips the band.
%     5. Altitude is reported first when both fail.
%
%% Inputs:
%
%  none
%
%% Outputs:
%
%  ok                       logical                 All checks passed
%
%% Revision History:
%  M. Casey                                                   (c) 09/16/2026
%  Copyright Coorbital Inc.
%% ------------------------ Begin Code Sequence ---------------------------

ok = true;
here = fileparts(fileparts(mfilename('fullpath')));
addpath(here);

mu      = 0.012150585609624;
lStar   = 389703.264829278;
rMoonKm = 1737.4;
floorKm = 500;

%% Nodes at chosen lunar altitudes along +x from the Moon:
alt = [5000 1200 800 3000];
X   = [1 - mu + (rMoonKm + alt)/lStar; zeros(2, 4)];
o   = struct('X', X, 'tf', 2);

%% 1. Node minimum:
[okS, why, minAlt] = preflight_screen(o, mu, lStar, floorKm, []);
ok = chk(ok, abs(minAlt - 800) < 1e-6, sprintf('minAltKm is the node minimum: %.9f km', minAlt));
ok = chk(ok, okS && isempty(why), 'a clean solution passes with an empty reason');

%% 2. Half-floor threshold:
oH = o;  oH.X(1,3) = 1 - mu + (rMoonKm + 251)/lStar;           % just above floor/2
[okH, whyH] = preflight_screen(oH, mu, lStar, floorKm, []);
oB = o;  oB.X(1,3) = 1 - mu + (rMoonKm + 249)/lStar;           % just below
[okB, whyB] = preflight_screen(oB, mu, lStar, floorKm, []);
ok = chk(ok, okH && isempty(whyH), 'an altitude just above half the floor passes');
ok = chk(ok, ~okB && strcmp(whyB, 'altitude'), 'an altitude below half the floor trips ''altitude''');

%% 3. Time-of-flight band:
[okHi, whyHi] = preflight_screen(o, mu, lStar, floorKm, 0.6);   % tf = 3.33 x seed
[okLo, whyLo] = preflight_screen(o, mu, lStar, floorKm, 7);     % tf = 0.29 x seed
okE1 = preflight_screen(o, mu, lStar, floorKm, 0.7);             % tf = 2.86 x seed
okE2 = preflight_screen(o, mu, lStar, floorKm, 6.5);             % tf = 0.31 x seed
ok = chk(ok, ~okHi && strcmp(whyHi, 'tf'), 'tf above 3 x seedTf trips ''tf''');
ok = chk(ok, ~okLo && strcmp(whyLo, 'tf'), 'tf below 0.3 x seedTf trips ''tf''');
ok = chk(ok, okE1 && okE2, 'just inside both edges passes');

%% 4. No seed:
ok = chk(ok, preflight_screen(o, mu, lStar, floorKm, []), 'seedTf = [] skips the band');

%% 5. Precedence:
[~, whyBoth] = preflight_screen(oB, mu, lStar, floorKm, 0.6);
ok = chk(ok, strcmp(whyBoth, 'altitude'), 'altitude is reported when both fail');

if ok, fprintf('TEST_PREFLIGHT_SCREEN: ALL PASS\n');
else,  fprintf('TEST_PREFLIGHT_SCREEN: FAILURE (see lines above)\n');
end
end

% ------------------------------------------------------------------------
function ok = chk(ok, cond, label)
%% Purpose:
%
%   Accumulate one labelled pass/fail.
%
if cond, fprintf('  PASS  %s\n', label);
else,    fprintf('  FAIL  %s\n', label);  ok = false;
end
end
