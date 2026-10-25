function ok = test_lagrange_dalembert()
%% Purpose:
%
%   mech.lagrange_dalembert on problems with exact answers: a free
%   particle (no constraint), a damped particle (generalized forces enter),
%   the planar pendulum written as a constrained particle (the multiplier
%   is the string tension), a knife-edge particle (v and w constant,
%   lambda = m v w), and the three refusals.
%
%% Inputs:
%
%  none
%
%% Outputs:
%
%  ok                       logical                 Every check passed
%
%% Revision History:
%  M. Casey                                                   (c) 10/25/2026
%  Copyright Coorbital Inc.
%% ------------------------ Begin Code Sequence ---------------------------

ok = true;
here = fileparts(mfilename('fullpath'));
addpath(fileparts(here));                      % mechlib (the +mech package root)

% 1. free particle: qdd = 0, lambda is [0 x 1]
syms x y xd yd real
L = (xd^2 + yd^2)/2;
sys = mech.lagrange_dalembert([x; y], [xd; yd], L, sym(zeros(0, 2)));
yd0 = sys.f(0, [1; 2; 3; 4]);
ok = check(ok, isequal(size(yd0), [4 1]) && max(abs(yd0 - [3; 4; 0; 0])) < 1e-15, ...
           'free particle: f = [qd; 0], column [4 x 1]');
ok = check(ok, isequal(size(sys.lambda([1; 2], [3; 4])), [0 1]) && sys.k == 0 && sys.n == 2, ...
           'free particle: lambda is [0 x 1], k = 0, n = 2');
ok = check(ok, abs(sys.energy([1; 2], [3; 4]) - 12.5) < 1e-15, 'free particle: energy = |qd|^2/2');

% 2. damped particle: Q = -c qd must reach the right-hand side
c = 0.7;
sysD = mech.lagrange_dalembert([x; y], [xd; yd], L, sym(zeros(0, 2)), -c*[xd; yd]);
ydD = sysD.f(0, [0; 0; 3; 4]);
ok = check(ok, max(abs(ydD(3:4) - (-c*[3; 4]))) < 1e-15, 'damped particle: qdd = -c qd (Q is used)');

% 3. pendulum as a constrained particle: A = [x, y], lambda = -T/l
m = 1; g = 9.81; l = 1.7;
Lp = m/2*(xd^2 + yd^2) - m*g*y;
sysP = mech.lagrange_dalembert([x; y], [xd; yd], Lp, [x, y]);
v = 0.7;
lam = sysP.lambda([0; -l], [v; 0]);
T = m*g + m*v^2/l;
ok = check(ok, abs(lam + T/l) < 1e-10, sprintf('pendulum at the bottom, l = %.1f: lambda = -T/l = %.10f (err %.2e)', l, -T/l, abs(lam + T/l)));
phi = 0.6;
qg = l*[sin(phi); -cos(phi)];  qdg = v*[cos(phi); sin(phi)];
lamg = sysP.lambda(qg, qdg);
Tg = m*g*cos(phi) + m*v^2/l;
ok = check(ok, abs(lamg + Tg/l) < 1e-10, sprintf('pendulum at phi = 0.6, l = %.1f: lambda = -T/l = %.10f (err %.2e)', l, -Tg/l, abs(lamg + Tg/l)));
Ep = sysP.energy([0; -l], [v; 0]);
ok = check(ok, abs(Ep - (m*v^2/2 + m*g*(-l))) < 1e-12, ...
           sprintf('pendulum energy = kinetic + potential (%.6f)', Ep));
ok = check(ok, isequal(size(sysP.qddSym), [2 1]) && isequal(size(sysP.lambdaSym), [1 1]), ...
           'symbolic reduced equations have the right shapes');

% 3b. integrate two periods: energy drift and the small-amplitude period
optsO = odeset('RelTol', 1e-12, 'AbsTol', 1e-12);
phi0 = 0.3;  y0 = [l*sin(phi0); -l*cos(phi0); 0; 0];
T0 = 2*pi*sqrt(l/g);
[tt, yy] = ode113(sysP.f, linspace(0, 2*T0, 801), y0, optsO);
E = arrayfun(@(kk) sysP.energy(yy(kk, 1:2).', yy(kk, 3:4).'), 1:numel(tt));
drift = max(abs(E - E(1)))/abs(E(1));
ok = check(ok, drift < 1e-9, sprintf('pendulum: energy drift over two periods %.2e (< 1e-9)', drift));
rad = max(abs(sqrt(sum(yy(:, 1:2).^2, 2)) - l));
ok = check(ok, rad < 1e-9, sprintf('pendulum: radius drift %.2e (< 1e-9)', rad));
optsE = odeset(optsO, 'Events', @xUpCrossing);
[~, ~, te] = ode113(sysP.f, [0 3*T0], [l*sin(1e-3); -l*cos(1e-3); 0; 0], optsE);
per = te(2) - te(1);
ok = check(ok, abs(per - T0)/T0 < 1e-6, sprintf('pendulum: small-amplitude period %.9f vs 2 pi sqrt(l/g) %.9f (rel %.1e < 1e-6)', per, T0, abs(per - T0)/T0));

% 4. knife-edge particle on q = (x, y, th): v, w constant, lambda = m v w
syms th thdot real
mk = 1.3; Jk = 0.4;
Lk = mk/2*(xd^2 + yd^2) + Jk/2*thdot^2;
sysK = mech.lagrange_dalembert([x; y; th], [xd; yd; thdot], Lk, [-sin(th), cos(th), 0]);
thk = 0.4; vk = 1.3; wk = 0.8;
qk = [0; 0; thk];  qdk = [vk*cos(thk); vk*sin(thk); wk];
lamk = sysK.lambda(qk, qdk);
ok = check(ok, abs(lamk - mk*vk*wk) < 1e-12, sprintf('knife edge: lambda = m v w (%.2e)', abs(lamk - mk*vk*wk)));
fk = sysK.f(0, [qk; qdk]);
vdot = fk(4)*cos(thk) + fk(5)*sin(thk) + wk*(-qdk(1)*sin(thk) + qdk(2)*cos(thk));
ok = check(ok, abs(vdot) < 1e-12 && abs(fk(6)) < 1e-12, sprintf('knife edge: vd = %.1e, wd = %.1e (both 0)', vdot, fk(6)));
ok = check(ok, abs(sysK.constraint(qk, qdk)) < 1e-15, 'knife edge: constraint handle returns A qd = 0 on an admissible velocity');

% 5. refusals
sysN = mech.lagrange_dalembert([x; y], [xd; yd], L, [1 0]);
fN = sysN.f(0, [0; 0; 0; 2]);
ok = check(ok, max(abs(fN - [0; 2; 0; 0])) < 1e-15 && abs(sysN.lambda([0; 0], [0; 2])) < 1e-15, ...
           'double A accepted (A = [1 0]: xd stays 0, lambda = 0)');
ok = check(ok, throws(@() mech.lagrange_dalembert([x; y], [xd; yd], 1.5, sym(zeros(0, 2))), 'mech:lagrange_dalembert:shape'), ...
           'double L refused with :shape');
ok = check(ok, throws(@() mech.lagrange_dalembert([x; y], xd, L, sym(zeros(0, 2))), 'mech:lagrange_dalembert:shape'), ...
           'mismatched qd refused with :shape');
ok = check(ok, throws(@() mech.lagrange_dalembert([x; y], [xd; yd], xd^2/2, sym(zeros(0, 2))), 'mech:lagrange_dalembert:M'), ...
           'singular mass matrix refused with :M');
ok = check(ok, throws(@() mech.lagrange_dalembert([x; y], [xd; yd], L, [x, y; 2*x, 2*y]), 'mech:lagrange_dalembert:A'), ...
           'rank-deficient A refused with :A');

if ok, fprintf('TEST_LAGRANGE_DALEMBERT: ALL PASS\n'); else, fprintf('TEST_LAGRANGE_DALEMBERT: FAILURES\n'); end
end

% ------------------------------------------------------------------------
function ok = check(ok, cond, msg)
%% Purpose:
%
%   Print one PASS/FAIL line and fold the condition into the verdict.
%
if cond, fprintf('  PASS  %s\n', msg); else, fprintf('  FAIL  %s\n', msg); end
ok = ok && logical(cond);
end

% ------------------------------------------------------------------------
function tf = throws(fun, id)
%% Purpose:
%
%   True when fun() throws an error with identifier id.
%
tf = false;
try
    fun();
catch err
    tf = strcmp(err.identifier, id);
end
end

% ------------------------------------------------------------------------
function [value, isterminal, direction] = xUpCrossing(~, y)
%% Purpose:
%
%   ode event: upward zero crossings of the first state, non-terminal.
%
value      = y(1);
isterminal = 0;
direction  = 1;
end
