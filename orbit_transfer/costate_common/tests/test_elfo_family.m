function ok = test_elfo_family()
% TEST_ELFO_FAMILY  Unit test for the 'elfo' algebraic arrival locus.
%
% Checks the Kepler period, the Moon-centred radius extremes, the closed
% seam, strict monotonicity of tau, agreement with the proj7 oracle
% (im_elfo_states at t = 0 for all five IM satellites, interpolated through
% phase_state), RAAN wrapping, and the unknown-family error.
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
addpath(fileparts(here));                          % costate_common
ok = true;

muStar = 0.012150585609624;
lStar  = 389703.264829278;

% 1. Default locus: period, radii, seam, monotone tau
[tau, rv, info] = get_family_orbit('elfo', struct());
ok = check(sprintf('period %.5f ND ~ 0.3080 (5e-5)', info.periodND), ...
           abs(info.periodND - 0.3080) < 5e-5 && abs(tau(end) - info.periodND) < 1e-12) && ok;
rMoon = vecnorm(rv(:,1:3) - [1-muStar 0 0], 2, 2)*lStar;
ok = check(sprintf('radius min %.1f / max %.1f km ~ 3720 / 20280 (1 km)', ...
           min(rMoon), max(rMoon)), ...
           abs(min(rMoon) - 3720) < 1 && abs(max(rMoon) - 20280) < 1) && ok;
ok = check('rv(1,:) == rv(end,:) (1e-12)', max(abs(rv(1,:) - rv(end,:))) < 1e-12) && ok;
ok = check('tau(1) = 0, strictly increasing, [2001x1]', ...
           tau(1) == 0 && all(diff(tau) > 0) && isequal(size(tau), [2001 1]) ...
           && isequal(size(rv), [2001 6]) && strcmp(info.family, 'elfo')) && ok;

% 2. Oracle: the five IM satellites at t = 0
[raanOpt, M0opt, planeOf, oe] = im_elfo_optimum();
[~, rTgt] = im_elfo_states(raanOpt(planeOf), M0opt, [], 4);
errOracle = 0;
for k = 1:numel(M0opt)
    [tk, rvk] = get_family_orbit('elfo', struct('oe', oe, 'raanDeg', raanOpt(planeOf(k))));
    xA = phase_state(tk, rvk);
    xk = xA(M0opt(k)/360);
    errOracle = max(errOracle, max(abs(xk(1:3).' - rTgt(1,:,k))));
end
ok = check(sprintf('oracle position agreement (5 sats) %.1e < 1e-9', errOracle), ...
           errOracle < 1e-9) && ok;

% 2b. Oracle velocity: the oracle returns positions only, so its rotating-
%     frame velocity at t = 0 is the one-sided second-order difference of
%     its PROPAGATED positions (dt = 2*pi/(Nt-1)); a wrong frame term or a
%     wrong sense of motion in the locus is an O(1) error against 1e-4
NtV = 200001;   dtV = 2*pi/(NtV - 1);
[~, rFine] = im_elfo_states(raanOpt(planeOf), M0opt, [], NtV);
errVel = 0;
for k = 1:numel(M0opt)
    [tk, rvk] = get_family_orbit('elfo', struct('oe', oe, 'raanDeg', raanOpt(planeOf(k))));
    xA = phase_state(tk, rvk);
    xk = xA(M0opt(k)/360);
    vFd = (-3*rFine(1,:,k) + 4*rFine(2,:,k) - rFine(3,:,k))/(2*dtV);
    errVel = max(errVel, max(abs(xk(4:6).' - vFd)));
end
ok = check(sprintf('oracle velocity agreement (5 sats) %.1e < 1e-4', errVel), ...
           errVel < 1e-4) && ok;

% 3. RAAN wraps mod 360
[~, rv370] = get_family_orbit('elfo', struct('raanDeg', 370));
[~, rv10]  = get_family_orbit('elfo', struct('raanDeg', 10));
ok = check('raanDeg 370 == 10 (1e-12)', max(abs(rv370(:) - rv10(:))) < 1e-12) && ok;

% 4. Unknown family still errors by name
id = '';
try
    get_family_orbit('nosuch', struct());
catch e
    id = e.identifier;
end
ok = check('unknown family errors get_family_orbit:family', ...
           strcmp(id, 'get_family_orbit:family')) && ok;

fprintf('test_elfo_family: %s\n', string(ok));
end

% ------------------------------------------------------------------------
function ok = check(name, cond)
% CHECK  One gate line.  INPUTS: name;cond  OUTPUTS: ok
ok = logical(cond);
if ok, tag = 'PASS'; else, tag = 'FAIL'; end
fprintf('  [%s] %s\n', tag, name);
end
