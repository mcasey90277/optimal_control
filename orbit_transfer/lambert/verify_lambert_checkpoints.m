% VERIFY_LAMBERT_CHECKPOINTS  Headless capture of tutorial checkpoint numbers.
%
%   Runs each planned checkpoint snippet against the reference solver and
%   prints true values, including element-by-element comparison against the
%   pyKep (Izzo) solutions produced by the standalone harness (the same
%   routine pumpkyn.pykep.lambert2Body wraps). pyKep numbers below were
%   generated on 2026-07-06 with lambert_harness.
%
% INPUTS: (none -- script)  OUTPUTS: (none -- prints)

here = fileparts(mfilename('fullpath'));
cd(here);

fprintf('===== CHECKPOINT A: Stumpff functions =====\n');
[C0, S0] = stumpff(0);
fprintf('C(0) = %.12f (expect 0.5)   S(0) = %.12f (expect 1/6 = %.12f)\n', C0, S0, 1/6);
[Cp, Sp] = stumpff(pi^2);
fprintf('C(pi^2) = %.12f (expect 2/pi^2 = %.12f)\n', Cp, 2/pi^2);
fprintf('S(pi^2) = %.12f (expect 1/pi^2 = %.12f)\n', Sp, 1/pi^2);
[Cn, Sn] = stumpff(-4);
fprintf('C(-4) = %.12f (expect (cosh2-1)/4 = %.12f)\n', Cn, (cosh(2)-1)/4);
fprintf('S(-4) = %.12f (expect (sinh2-2)/8 = %.12f)\n', Sn, (sinh(2)-2)/8);
% continuity across the series/closed-form switch at |z| = 1e-4
zL = 1e-4 - 1e-12;  zR = 1e-4 + 1e-12;
[CL, SL] = stumpff(zL);  [CR, SR] = stumpff(zR);
fprintf('switch continuity: |dC| = %.2e, |dS| = %.2e\n', abs(CL-CR), abs(SL-SR));

fprintf('\n===== CHECKPOINT B: the TOF curve =====\n');
r1 = [1; 0; 0];  r2 = [0; 1.2; 0];  mu = 1;
r1n = norm(r1);  r2n = norm(r2);
A = sin(pi/2) * sqrt(r1n*r2n / (1 - cos(pi/2)));
fprintf('A = %.12f (expect sqrt(1.2) = %.12f)\n', A, sqrt(1.2));
[t0, y0] = lambert_tof(0, r1n, r2n, A, mu);
fprintf('parabolic: y(0) = %.12f   t(0) = %.12f\n', y0, t0);
tm5 = lambert_tof(-5, r1n, r2n, A, mu);
t20 = lambert_tof(20, r1n, r2n, A, mu);
fprintf('t(-5) = %.12f   t(20) = %.12f  (monotone increasing)\n', tm5, t20);
zz = linspace(-20, (2*pi)^2 - 1e-6, 2000);
tt = lambert_tof(zz, r1n, r2n, A, mu);
fprintf('monotone on grid: %d violations\n', sum(diff(tt(~isnan(tt))) <= 0));

fprintf('\n===== CHECKPOINT C1: canonical elliptic single-rev =====\n');
dt = 2.0;
[v1, v2, info] = lambert_uv(r1, r2, dt, mu, +1);
% pyKep (harness, 2026-07-06):
pk_v1 = [ 1.7292152422437140e-01;  9.9659460660332588e-01; 0];
pk_v2 = [-8.3049550550277162e-01; -6.8224231238170905e-03; 0];
fprintf('v1 = [% .12f % .12f % .12f]\n', v1);
fprintf('v2 = [% .12f % .12f % .12f]\n', v2);
fprintf('vs pyKep: |dv1| = %.3e, |dv2| = %.3e   (z = %.10f, %d iters, resid %.1e s)\n', ...
        norm(v1-pk_v1), norm(v2-pk_v2), info.z, info.iters, info.resid);
% independent truth: propagate (r1,v1) for dt with tight-tol ode89
opts = odeset('RelTol', 3e-14, 'AbsTol', 1e-14);
f2b = @(t, x) [x(4:6); -mu*x(1:3)/norm(x(1:3))^3];
sol = ode89(f2b, [0 dt], [r1; v1], opts);
xf  = deval(sol, dt);
fprintf('round-trip: ||r(dt) - r2|| = %.3e   ||v(dt) - v2|| = %.3e\n', ...
        norm(xf(1:3) - r2), norm(xf(4:6) - v2));

fprintf('\n===== CHECKPOINT C2: hyperbolic =====\n');
[v1h, v2h, infoh] = lambert_uv(r1, r2, 0.5, mu, +1);
pk_v1h = [-1.7426333047972373e+00;  2.5599239004201517e+00; 0];
pk_v2h = [-2.1332699170167930e+00;  2.1692872882005960e+00; 0];
fprintf('z = %.10f (< 0: hyperbolic)   vs pyKep: |dv1| = %.3e, |dv2| = %.3e\n', ...
        infoh.z, norm(v1h-pk_v1h), norm(v2h-pk_v2h));
vinf_check = norm(v1h)^2/2 - mu/r1n;   % specific energy > 0
fprintf('specific energy = %.6f (> 0 confirms hyperbolic)\n', vinf_check);

fprintf('\n===== CHECKPOINT C3: retrograde =====\n');
[v1r, ~, infor] = lambert_uv(r1, r2, 2.0, mu, -1);
pk_v1r = [-7.3835272264009166e-01; -7.3862253462995409e-01; 0];
fprintf('dtheta = %.6f rad (expect 3*pi/2 = %.6f)   vs pyKep: |dv1| = %.3e\n', ...
        infor.dtheta, 3*pi/2, norm(v1r-pk_v1r));

fprintf('\n===== CHECKPOINT C3b: fast retrograde (regression: used to hang) =====\n');
[v1f, ~, infof] = lambert_uv(r1, r2, 0.1, mu, -1);
pk_v1f = [-2.1809220323093157e+01; -4.5772112892885693e-02; 0];
fprintf('z = %.4f   |v1| = %.4f   vs pyKep: |dv1| = %.3e\n', ...
        infof.z, norm(v1f), norm(v1f - pk_v1f));

fprintf('\n===== CHECKPOINT C4: Vallado Example 7-5 (dimensional) =====\n');
r1v = [15945.34; 0; 0];  r2v = [12214.83899; 10249.46731; 0];
[v1v, v2v] = lambert_uv(r1v, r2v, 4560, 398600.4418, +1);
fprintf('v1 = [% .6f % .6f % .6f] km/s  (Vallado: [2.058913 2.915965 0])\n', v1v);
fprintf('v2 = [% .6f % .6f % .6f] km/s  (Vallado: [-3.451565 0.910315 0])\n', v2v);
pk_v1v = [2.0589133537073092;  2.9159643516499401; 0];
pk_v2v = [-3.4515648446831912; 0.9103142481137412; 0];
fprintf('vs pyKep: |dv1| = %.3e, |dv2| = %.3e km/s\n', norm(v1v-pk_v1v), norm(v2v-pk_v2v));

fprintf('\n===== CHECKPOINT D: multi-rev =====\n');
dt = 40;
[V1, V2, Nrev, out] = lambert_uv_multirev(r1, r2, dt, mu, +1, 20);
fprintf('Nmax = %d   solutions = %d  (expect 6 and 13, matching pyKep)\n', ...
        out.Nmax, size(V1,2));
fprintf('t_min per band: '); fprintf('%.4f ', out.tmins); fprintf('\n');
% pyKep all 13 v1 solutions (harness, dt = 40, prograde, 2026-07-06):
PK1 = [ 1.1638842645798422e+00  6.0077127301169997e-01
        1.0757482493429424e+00  6.2600861812789099e-01
       -2.7802678771263317e-01  1.2748898978246603e+00
        9.9299877679860105e-01  6.5118780901843898e-01
       -1.9698872932539305e-01  1.2199961425710959e+00
        9.0919572923183256e-01  6.7824286360619634e-01
       -1.1795702571072726e-01  1.1685032300603886e+00
        8.1928049114088419e-01  7.0911458948558093e-01
       -3.3799176029372768e-02  1.1159123174428698e+00
        7.1530384696632887e-01  7.4733690273625841e-01
        6.4030691140523774e-02  1.0577001804443587e+00
        5.7264178790210662e-01  8.0447882017176531e-01
        2.0038760104803852e-01  9.8179097298348761e-01 ]';
% match each of my columns to nearest pyKep column
err = zeros(1, size(V1,2));  perm = zeros(1, size(V1,2));
for c = 1:size(V1,2)
    d = vecnorm(PK1 - V1(1:2, c));
    [err(c), perm(c)] = min(d);
end
fprintf('column matching to pyKep: perm = ['); fprintf('%d ', perm); fprintf(']\n');
fprintf('max |v1 - pyKep| across all %d solutions = %.3e\n', size(V1,2), max(err));
% round-trip every solution
rt = zeros(1, size(V1,2));
for c = 1:size(V1,2)
    solc = ode89(f2b, [0 dt], [r1; V1(:,c)], opts);
    xfc  = deval(solc, dt);
    rt(c) = norm(xfc(1:3) - r2);
end
fprintf('max round-trip ||r(dt) - r2|| over all solutions = %.3e\n', max(rt));

fprintf('\n===== CHECKPOINT D2: infeasible band =====\n');
[~, ~, ~, out7] = lambert_uv_multirev(r1, r2, 7.0, mu, +1, 20);
fprintf('dt = 7: Nmax = %d, t_min(1) = %.4f (dt below it -> 0-rev only? %d)\n', ...
        out7.Nmax, out7.tmins(1), out7.Nmax == 0);
