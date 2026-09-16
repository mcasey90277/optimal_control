function ok = test_true_min_altitude()
%% Purpose:
%
%   Tests true_min_altitude -- the minimum lunar altitude of the PROPAGATED
%   trajectory, not of the nodes -- on a coasting lunar flyby whose
%   periselene falls strictly BETWEEN two nodes:
%     1. It finds the between-node minimum: never below a dense (20,000
%        sample) propagation of the same dynamics, and within 2 km of it
%        (the resolution of its 65 samples per interval).
%     2. It is not the node minimum: preflight_screen's node-only altitude
%        is hundreds of km higher on the same solution.
%     3. rMoonKm enters linearly: +10 km radius, -10 km altitude.
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
Tmax    = 0.05;  c = 3;

%% A coasting flyby, 2000 km periselene, built backward and forward from it:
rP   = (rMoonKm + 2000) / lStar;
vP   = 1.25 * sqrt(mu / rP);                         % faster than circular: a flyby
zP   = [1 - mu + rP; 0; 0; 0; vP; 0; 1];
coast = @(t, z) cr3bp_thrust_rhs(z, [1; 0; 0; 0], mu, Tmax, c);
oo   = odeset('RelTol', 1e-13, 'AbsTol', 1e-15);
half = 0.015;
[~, Zb] = ode113(coast, [0 -half], zP, oo);          % back to the first node
z0   = Zb(end,:).';
tN   = linspace(0, 2*half, 4);                       % periselene near t = half, between nodes 2 and 3
[~, Z] = ode113(coast, tN, z0, oo);
o    = struct('X', Z.', 'U', repmat([1; 0; 0; 0], 1, 4), 'Um', [], 'tNodes', tN, 'tf', 2*half);

%% Dense reference:
[~, Zd] = ode113(coast, linspace(0, 2*half, 20001), z0, oo);
refKm = min(vecnorm(Zd(:,1:3) - [1-mu 0 0], 2, 2)) * lStar - rMoonKm;

%% 1. Between-node minimum:
amin = true_min_altitude(o, mu, Tmax, c, lStar, rMoonKm);
ok = chk(ok, amin >= refKm - 0.01 && amin - refKm < 2, ...
         sprintf('propagated minimum %.3f km vs dense %.3f km (65 samples per interval)', amin, refKm));

%% 2. Not the node minimum:
[~, ~, nodeKm] = preflight_screen(o, mu, lStar, 0, []);
ok = chk(ok, nodeKm > amin + 100, sprintf('node-only minimum %.0f km sits above it (periselene between nodes)', nodeKm));

%% 3. Radius enters linearly:
amin10 = true_min_altitude(o, mu, Tmax, c, lStar, rMoonKm + 10);
ok = chk(ok, abs((amin - amin10) - 10) < 1e-9, sprintf('+10 km radius -> %.9f km lower', amin - amin10));

if ok, fprintf('TEST_TRUE_MIN_ALTITUDE: ALL PASS\n');
else,  fprintf('TEST_TRUE_MIN_ALTITUDE: FAILURE (see lines above)\n');
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
