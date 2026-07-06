% RUN_LAMBERT_DEMO  Demo: the universal-variables Lambert solver, multi-rev.
%
%   Canonical instance (mu = 1): r1 = [1;0;0], r2 = [0;1.2;0], 90-deg
%   prograde transfer, dt = 40 TU. Finds all 13 solutions (0-rev plus a
%   pair for each N = 1..6), prints the solution table, and draws:
%   (left)  the time-of-flight curve t(z) across the revolution bands with
%           the dt = 40 line and the 13 converged roots;
%   (right) the 13 transfer arcs, propagated with the returned v1.
%
% INPUTS:  (none -- script)   OUTPUTS: (none -- prints + one figure)

clear; clc;
r1 = [1; 0; 0];  r2 = [0; 1.2; 0];  mu = 1;  dt = 40;

[V1, V2, Nrev, out] = lambert_uv_multirev(r1, r2, dt, mu, +1, 20);
nsol = size(V1, 2);
fprintf('dt = %g: Nmax = %d, %d solutions\n', dt, out.Nmax, nsol);
fprintf('%3s %5s %12s %14s %14s\n', '#', 'N', 'z', 'v1_x', 'v1_y');
for c = 1:nsol
    fprintf('%3d %5d %12.6f %14.10f %14.10f\n', ...
            c, Nrev(c), out.zs(c), V1(1,c), V1(2,c));
end

% --- figure ---------------------------------------------------------------
fig = figure('Position', [60 60 1150 470], 'Color', 'w');
fig.Theme = 'light';

% left: the TOF curve
subplot(1,2,1); hold on; grid on;
r1n = norm(r1); r2n = norm(r2);
zz = linspace(-10, (2*pi*7)^2, 20000);
tt = lambert_tof(zz, r1n, r2n, out.A, mu);
plot(zz, tt, 'b-', 'LineWidth', 1.1);
yline(dt, 'r--', 'LineWidth', 1.2);
for N = 1:7
    xline((2*pi*N)^2, 'k:', 'Alpha', 0.4);
end
plot(out.zs, dt*ones(1, nsol), 'ko', 'MarkerFaceColor', 'y', 'MarkerSize', 6);
ylim([0 70]); xlabel('z'); ylabel('t(z)  [TU]');
title(sprintf('TOF curve: %d roots of t(z) = %g', nsol, dt));

% right: the transfer arcs
subplot(1,2,2); hold on; grid on; axis equal;
f2b = @(t, x) [x(4:6); -mu*x(1:3)/norm(x(1:3))^3];
opts = odeset('RelTol', 1e-12, 'AbsTol', 1e-12);
cmap = turbo(nsol);
for c = 1:nsol
    sol = ode89(f2b, [0 dt], [r1; V1(:,c)], opts);
    tq  = linspace(0, dt, 3000);
    xq  = deval(sol, tq);
    plot(xq(1,:), xq(2,:), '-', 'Color', [cmap(c,:) 0.75], 'LineWidth', 0.9);
end
plot(0, 0, 'k.', 'MarkerSize', 22);
plot(r1(1), r1(2), 'ks', 'MarkerFaceColor', 'g', 'MarkerSize', 8);
plot(r2(1), r2(2), 'kp', 'MarkerFaceColor', 'r', 'MarkerSize', 11);
xlabel('x'); ylabel('y');
title('All 13 transfer arcs (0 through 6 revolutions)');
