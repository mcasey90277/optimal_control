function [v1, v2, info] = lambert_uv(r1, r2, dt, mu, dir)
% LAMBERT_UV  Single-revolution Lambert solver via universal variables.
%
%   Solves the two-body two-point BVP: given position vectors r1, r2 and
%   a time of flight dt, find the transfer orbit's terminal velocities.
%   Universal-variables formulation (one iteration variable z covers
%   hyperbolic / parabolic / elliptic seamlessly): bracketed bisection on
%   lambert_tof(z) - dt over z in (-4*pi^2, (2*pi)^2), then the f-and-g
%   functions give the velocities.
%
% INPUTS:
%   r1  - Initial position vector [3x1]
%   r2  - Final position vector [3x1]
%   dt  - Time of flight [scalar, > 0]
%   mu  - Gravitational parameter [scalar]
%   dir - Transfer sense: +1 = prograde (angular momentum along +z),
%         -1 = retrograde [scalar] (matches pumpkyn lambert2Body's flag)
%
% OUTPUTS:
%   v1   - Velocity at r1 [3x1]
%   v2   - Velocity at r2 [3x1]
%   info - Struct: z (converged), dtheta [rad], A, y, iters, resid [s]
%
% REFERENCES:
%   [1] Curtis, "Orbital Mechanics for Engineering Students", Alg. 5.2.
%   [2] Bate, Mueller, White, "Fundamentals of Astrodynamics", Ch. 5.

r1n = norm(r1);
r2n = norm(r2);

% --- transfer angle, resolved by the requested sense ---------------------
cross_z = r1(1)*r2(2) - r1(2)*r2(1);
dtheta  = acos(max(-1, min(1, dot(r1, r2)/(r1n*r2n))));
if (dir > 0 && cross_z < 0) || (dir < 0 && cross_z >= 0)
    dtheta = 2*pi - dtheta;
end
if abs(sin(dtheta)) < 1e-12
    error('lambert_uv:degenerate', ...
        'Transfer angle is 0 or pi: the transfer plane is undefined.');
end

A = sin(dtheta) * sqrt(r1n*r2n / (1 - cos(dtheta)));

% --- bracket the root of t(z) - dt by adaptive coarse scan ----------------
% On the valid (y > 0) part of the single-rev window t(z) crosses dt once,
% but the location of the valid domain's left edge depends on the geometry
% (sign of A), and very fast transfers need very negative z. Scan a window,
% expanding leftward until a faster-than-dt point is found. The hard floor
% z ~ -4.7e5 is where sinh(sqrt(-z)) overflows double precision -- the
% unbounded iteration variable is a real limitation of the z-formulation
% (Izzo's solver iterates on a bounded variable for exactly this reason).
zhi_edge = (2*pi)^2 - 1e-10;
window = 4*pi^2;
idx = [];
for attempt = 1:12
    zgrid = linspace(-window, zhi_edge, 6000);
    tgrid = lambert_tof(zgrid, r1n, r2n, A, mu);
    idx = find(tgrid < dt, 1, 'last');  % rightmost faster-than-dt point
    if ~isempty(idx) && idx < numel(zgrid)
        break
    end
    window = window * 8;
    if window > 4.7e5
        error('lambert_uv:tooFast', ...
            ['No single-rev bracket found: dt is faster than the scan ' ...
             'window supports for this geometry/sense.']);
    end
end
zlo = zgrid(idx);                       % t(zlo) < dt
zhi = zgrid(idx + 1);                   % t(zhi) > dt (or the window edge)

% --- bisection ------------------------------------------------------------
% The scan bracket is NaN-free by construction: t is finite at both ends
% and the valid (y > 0) set is an interval, so every midpoint is valid.
% (The infeasible y < 0 region, when it exists, lies strictly on the FAST
% side -- left of the achievable-time domain -- never inside a bracket.)
iters = 0;
for k = 1:200
    zm = 0.5*(zlo + zhi);
    tm = lambert_tof(zm, r1n, r2n, A, mu);
    iters = k;
    if tm > dt
        zhi = zm;                       % too slow: go left
    else
        zlo = zm;                       % too fast: go right
    end
    if zhi - zlo < 1e-14 * max(1, abs(zm))
        break
    end
end
z = 0.5*(zlo + zhi);

% --- f and g functions -> velocities --------------------------------------
[t, y] = lambert_tof(z, r1n, r2n, A, mu);
f    = 1 - y/r1n;
g    = A * sqrt(y/mu);
gdot = 1 - y/r2n;
v1 = (r2 - f*r1) / g;
v2 = (gdot*r2 - r1) / g;

info = struct('z', z, 'dtheta', dtheta, 'A', A, 'y', y, ...
              'iters', iters, 'resid', t - dt);
end
