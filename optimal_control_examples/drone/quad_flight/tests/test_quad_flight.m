function ok = test_quad_flight()
%% Purpose:
%
%   The quadrotor flight script's pieces, pinned independently of each
%   other:
%
%     1. quad_frames round-trips (to_lib then from_lib is the identity) for
%        random states in both frames, and the NED convention is the
%        aerospace one: NED hover (zero angles, FRD body) is the library's
%        hover with x_b pointing north (yaw 90 deg), NED yaw 90 deg heads
%        east (library yaw 0), body rates flip q and r.
%     2. simulate_quad at hover trim holds the state; free fall from rest
%        drops g t^2 / 2 (mechanics, not the integrator's own output); the
%        zero-order hold honours a thrust step exactly at its breakpoint.
%     3. The same physical flight entered in ENU and in NED gives the same
%        world trajectory (the frame layer is transparent to the physics).
%     4. quad_georef returns the origin at the origin and moves latitude by
%        1/(M + h) rad per metre north (WGS-84 meridian radius).
%     5. plot_quad_flight and movie_quad write their files.
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
%  M. Casey                                                   (c) 10/23/2026
%  Copyright Coorbital Inc.
%% ------------------------ Begin Code Sequence ---------------------------

ok   = true;
here = fileparts(mfilename('fullpath'));
addpath(fileparts(here), fullfile(fileparts(fileparts(here)), 'drone_dynamics'));
p    = drone_params('rigid6dof');
rng(7);

%% 1. frames
P = [0 1 0; 1 0 0; 0 0 -1];                      % NED -> ENU (self-inverse)
worst = 0;
for frame = {'ENU', 'NED'}
    X = [randn(6, 5); 0.8*randn(3, 5); randn(3, 5)];
    Y = quad_frames('from_lib', quad_frames('to_lib', X, frame{1}), frame{1});
    worst = max(worst, max(max(abs(Y - X))));
end
ok = chk(ok, worst < 1e-13, sprintf('to_lib then from_lib is the identity in both frames (%.1e)', worst));

xE = quad_frames('to_lib', zeros(12, 1), 'ENU');
ok = chk(ok, isequal(xE, zeros(12, 1)), 'ENU is the library frame: identity');

xN = quad_frames('to_lib', zeros(12, 1), 'NED');
ok = chk(ok, max(abs(xN(1:6))) < 1e-15 && abs(xN(7)) < 1e-15 && abs(xN(8)) < 1e-15 ...
         && abs(xN(9) - pi/2) < 1e-15, 'NED hover (zero angles) = library hover with x_b north (yaw 90 deg)');
xN = quad_frames('to_lib', [zeros(8, 1); pi/2; zeros(3, 1)], 'NED');
ok = chk(ok, max(abs(xN(7:9))) < 1e-15, 'NED yaw 90 deg (heading east) = library yaw 0');
xN = quad_frames('to_lib', [1; 2; 3; 4; 5; 6; 0; 0; 0; 0.1; 0.2; 0.3], 'NED');
ok = chk(ok, isequal(xN(1:3), P*[1; 2; 3]) && isequal(xN(4:6), P*[4; 5; 6]) ...
         && isequal(xN(10:12), [0.1; -0.2; -0.3]), 'NED position/velocity through P, FRD body rates flip q and r');
% a general NED attitude: check the stated matrix identity directly
ThN = [20; -10; 35]*pi/180;
xN  = quad_frames('to_lib', [zeros(6, 1); ThN; zeros(3, 1)], 'NED');
Re  = euler_zyx(xN(7:9));
Rn  = euler_zyx(ThN);
ok = chk(ok, max(max(abs(Re - P*Rn*diag([1 -1 -1])))) < 1e-15, 'NED attitude: R_lib = P R_ned diag(1,-1,-1), a proper rotation');
ok = chk(ok, abs(det(Re) - 1) < 1e-14, 'converted attitude has det +1');

%% 2. simulate_quad
hover = p.m*p.g/4;
x0 = [0; 0; 1; zeros(9, 1)];
ctrl = struct('t', 0, 'f', hover*ones(4, 1));
[t, X, U] = simulate_quad(x0, ctrl, 2, p, struct('dt', 0.05));
ok = chk(ok, max(max(abs(X - x0))) < 1e-9, sprintf('hover trim holds the state for 2 s (%.1e)', max(max(abs(X - x0)))));
ok = chk(ok, abs(t(end) - 2) < 1e-15 && isequal(size(U), [4 numel(t)]), 'output grid reaches tEnd, U has one column per time');

ctrl = struct('t', 0, 'f', zeros(4, 1));
[t, X] = simulate_quad(x0, ctrl, 1, p, struct('dt', 0.1));
ok = chk(ok, abs(X(3, end) - (1 - 0.5*p.g*1^2)) < 1e-9, 'free fall from rest: z(1) = z0 - g/2 (mechanics)');
ok = chk(ok, abs(X(6, end) + p.g) < 1e-9, 'free fall: vz(1) = -g');

ctrl = struct('t', [0 0.5], 'f', [hover*ones(4, 1), 1.2*hover*ones(4, 1)]);
[t, X, U] = simulate_quad(x0, ctrl, 1, p, struct('dt', 0.1));
k = find(abs(t - 0.5) < 1e-12);
ok = chk(ok, numel(k) == 1 && max(abs(U(:, k) - 1.2*hover)) < 1e-15 && max(abs(U(:, k-1) - hover)) < 1e-15, ...
         'zero-order hold: the step is taken at its breakpoint, once');
% after the step the net acceleration is 0.2 g upward: z(1) = 1 + 0.5*0.2*g*0.5^2
ok = chk(ok, abs(X(3, end) - (1 + 0.1*p.g*0.25)) < 1e-9, 'thrust step at 0.5 s: climb matches a = 0.2 g for 0.5 s');

%% 3. ENU and NED entries of one physical flight
% a 2 % differential for 0.1 s each way: the Crazyflie's inertia is tiny
% (tau / I ~ 424 rad/s^2 per unit differential), so this is a 5-degree pitch
ctrl = struct('t', [0 0.4 0.5 0.6], 'f', hover*[1 1 1 1; 0.98 1 1.02 1; 1.02 1 0.98 1; 1 1 1 1].');
xUserE = [0; 0; 1; zeros(9, 1)];                             % ENU: at 1 m up, heading east
xUserN = [0; 0; -1; 0; 0; 0; 0; 0; pi/2; 0; 0; 0];           % NED: same point, FRD heading east
[tE, XE] = simulate_quad(quad_frames('to_lib', xUserE, 'ENU'), ctrl, 1.5, p, struct('dt', 0.05));
[tN, XN] = simulate_quad(quad_frames('to_lib', xUserN, 'NED'), ctrl, 1.5, p, struct('dt', 0.05));
worst = max(max(abs(XE - XN)));
ok = chk(ok, isequal(tE, tN) && worst < 1e-12, sprintf('ENU and NED entries fly the same world trajectory (%.1e)', worst));
XNu = quad_frames('from_lib', XN, 'NED');
ok = chk(ok, max(abs(XNu(3, :) + XE(3, :))) < 1e-12 && max(abs(XNu(1, :) - XE(2, :))) < 1e-12, ...
         'reported back in NED: down = -up, north = ENU y');

%% 4. georef
lla0 = [33.77, -118.39, 50];
[lat, lon, alt] = quad_georef([0 100 0; 0 0 100; 0 0 0], lla0);
ok = chk(ok, abs(lat(1) - lla0(1)) < 1e-12 && abs(lon(1) - lla0(2)) < 1e-12 && abs(alt(1) - lla0(3)) < 1e-6, ...
         'origin maps to lla0');
a = 6378137;  e2 = 0.00669437999014;  phi = lla0(1)*pi/180;
M = a*(1 - e2)/(1 - e2*sin(phi)^2)^1.5;                     % meridian radius of curvature
dlat = (lat(3) - lla0(1))*pi/180;
ok = chk(ok, abs(dlat*(M + lla0(3)) - 100) < 1e-3, sprintf('100 m north = 100/(M+h) rad of latitude (%.4f m)', dlat*(M + lla0(3))));
% 100 m along the tangent plane leaves the ellipsoid by d^2/(2R) ~ 0.8 mm
% and tilts latitude by ~ 5e-9 deg (second order); both are physics
ok = chk(ok, abs(lat(2) - lla0(1)) < 1e-7 && lon(2) > lla0(2) && abs(alt(2) - lla0(3)) < 2e-3, ...
         '100 m east: longitude grows, latitude and altitude hold to second order');

%% 5. plots and movie
outDir = fullfile(tempname, 'quad_flight_test');
mkdir(outDir);
info = struct('frame', 'ENU', 'lla0', lla0, 'label', 'test');
files = plot_quad_flight(tE, XE, hover*ones(4, numel(tE)), info, outDir);
ok = chk(ok, numel(files) == 2 && all(cellfun(@(f) exist(f, 'file') == 2, files)), 'plot_quad_flight writes its two PNGs');
mov = movie_quad(tE, XE, info, fullfile(outDir, 'test'), struct('fps', 10, 'seconds', 1));
d = dir(mov);
ok = chk(ok, ~isempty(d) && d.bytes > 1e4, sprintf('movie_quad writes an MP4 (%d bytes)', d.bytes));
close all;
rmdir(fileparts(outDir), 's');

if ok, fprintf('TEST_QUAD_FLIGHT: ALL PASS\n'); else, fprintf('TEST_QUAD_FLIGHT: FAILURE\n'); end
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
