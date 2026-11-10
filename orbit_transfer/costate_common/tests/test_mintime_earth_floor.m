function ok = test_mintime_earth_floor()
%% Purpose:
%
%   Tests the EARTH altitude floor of casadi_mintime_dro (opts.minAltEarthKm,
%   added 2026-11-09) on the GTO -> DRO cell that motivated it: orient 0,
%   DRO tau 1, departure phase 0.10, arrival phase 0, 15 N, 150 kg,
%   Isp 1710 s, N = 400, Hermite-Simpson, thrLock, cold t_f guess 0.30 ND.
%   Without the floor the cold solve lands on a defect-clean "solution" that
%   passes inside the Earth (t_f ~0.2937 ND); with a 200 km floor it lands on
%   a benign basin (t_f ~0.3040 ND) where the floor does not bind. A second
%   floor-satisfying basin exists at t_f ~0.2668 ND (warm-startable, not
%   reached by floored cold starts): the floor removes the Earth-grazing
%   basin, it does not pick among the remaining ones.
%     a. default off: .minAltEarthKm is NaN, .altEarthMinKm still computed.
%     b. the floored solve succeeds with a clean defect.
%     c. the floor is recorded and honoured at the nodes.
%     d. the floor changed the solution.
%     e. the floored solution is a DIFFERENT basin well clear of the surface
%        (node minimum > 1000 km, not an arc hugging the 200 km floor), the
%        floor is not binding, and t_f is a sane transfer time (0.20-0.35 ND).
%        It does NOT establish which benign basin is reached.
%     f. default-off identity: the unfloored t_f is the engine's known value.
%   Two direct solves: the floored one takes ~100 s (the whole test ~2
%   min). The bitwise default-path regression (floor off: X, U, Um, t_f and
%   lamDef identical, trapezoid and Hermite-Simpson) is
%   DRO_tulip/direct/certify/tests/test_minenergy_objective.
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
%  M. Casey                                                   (c) 11/09/2026
%  M. Casey  check (e) re-pinned to what the floor proves     (c) 11/09/2026
%  Copyright Coorbital Inc.
%% ------------------------ Begin Code Sequence ---------------------------

ok = true;
cc = fileparts(fileparts(mfilename('fullpath')));       % costate_common
addpath(cc);
if isempty(which('casadi.SX')), addpath(fullfile(getenv('HOME'), 'casadi-3.7.0')); end

%% The cell:
muStar = 0.012150585609624;  lStar = 389703.264829278;  tStar = 382981.289129055;
ndp = nd_propulsion(15, 1710, 150, lStar, tStar);
[tD, rvD] = get_family_orbit('gto', struct('orientDeg', 0, 'sma_km', 24446, 'ecc', 0.724781150290436));
[tA, rvA] = get_family_orbit('dro', struct('tau', 1));
depAt = phase_state(tD, rvD);  arrAt = phase_state(tA, rvA);
rv0 = depAt(0.10).';  rvf = arrAt(0).';
base = struct('maxIter', 6000, 'scheme', 'hermite-simpson', 'sundman', false, ...
              'returnModel', true, 'minAltKm', 500, 'maxCpuSec', 600, 'thrLock', true);
tfRef0 = 0.29373476167281287;           % the PRE-EDIT engine's value on this cell
                                        % (measured bitwise 2026-11-09, R2026a, CasADi 3.7.0)

%% The two solves:
t0 = tic;
o0 = casadi_mintime_dro(rv0(1:6), rvf(1:6), ndp.Tnd, ndp.cnd, muStar, 400, [], [], 0.30, base);
w0 = toc(t0);  t0 = tic;
bE = base;  bE.minAltEarthKm = 200;
oE = casadi_mintime_dro(rv0(1:6), rvf(1:6), ndp.Tnd, ndp.cnd, muStar, 400, [], [], 0.30, bE);
wE = toc(t0);
fprintf('  no floor : success %d  defect %.2e  tf %.10f ND  altEarthMin %.1f km  active %d  (%.1f s)\n', ...
        o0.success, o0.maxDefect, o0.tf, o0.altEarthMinKm, o0.altEarthActive, w0);
fprintf('  200 km   : success %d  defect %.2e  tf %.10f ND  altEarthMin %.1f km  active %d  (%.1f s)\n', ...
        oE.success, oE.maxDefect, oE.tf, oE.altEarthMinKm, oE.altEarthActive, wE);

%% a. Default off, diagnostic always on:
ok = chk(ok, isnan(o0.minAltEarthKm) && isfinite(o0.altEarthMinKm) && ~o0.altEarthActive, ...
         sprintf('default off: minAltEarthKm NaN, altEarthMinKm computed (%.1f km)', o0.altEarthMinKm));

%% b. The floored solve converges:
ok = chk(ok, oE.success && oE.maxDefect < 1e-9, ...
         sprintf('floored solve succeeds, defect %.2e', oE.maxDefect));

%% c. Recorded and honoured:
ok = chk(ok, oE.minAltEarthKm == 200 && oE.altEarthMinKm >= 200 - 1, ...
         sprintf('floor recorded (200 km) and honoured at the nodes (%.1f km)', oE.altEarthMinKm));

%% d. The floor changed the solution:
ok = chk(ok, oE.altEarthMinKm > o0.altEarthMinKm, ...
         sprintf('floor changed the solution: %.1f km > %.1f km', oE.altEarthMinKm, o0.altEarthMinKm));

%% e. A different basin, well clear, floor not binding, sane t_f:
ok = chk(ok, oE.altEarthMinKm > 1000, ...
         sprintf('e1 floored basin well clear of the surface: altEarthMin %.1f km (> 1000)', oE.altEarthMinKm));
ok = chk(ok, ~oE.altEarthActive, 'e2 the Earth floor is not binding');
ok = chk(ok, oE.tf > 0.20 && oE.tf < 0.35, ...
         sprintf('e3 sane transfer time: tf %.10f ND in (0.20, 0.35)', oE.tf));
fprintf('  NOTE  a second floor-satisfying basin exists at t_f ~0.2668 ND (warm-startable, not reached by floored cold starts, 2026-11-09)\n');

%% f. Default-off identity with the engine's known unconstrained value:
ok = chk(ok, abs(o0.tf - tfRef0) < 1e-9, ...
         sprintf('unfloored tf %.14f ND matches the pre-edit engine''s %.14f (|d| = %.1e)', ...
                 o0.tf, tfRef0, abs(o0.tf - tfRef0)));

if ok, fprintf('TEST_MINTIME_EARTH_FLOOR: ALL PASS\n');
else,  fprintf('TEST_MINTIME_EARTH_FLOOR: FAILURE (see lines above)\n');
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
