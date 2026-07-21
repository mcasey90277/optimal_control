function [V1, V2, Ns, out] = lambert_uv_multirev(r1, r2, dt, mu, dir, Nmax_req)
% LAMBERT_UV_MULTIREV  All Lambert solutions up to Nmax revolutions.
%
%   For N full revolutions the universal variable lives in the band
%   z in ((2*pi*N)^2, (2*pi*(N+1))^2), where the time-of-flight curve
%   t(z) is U-shaped: it diverges at both band edges and has one interior
%   minimum t_min(N). If dt > t_min there are exactly TWO N-rev solutions
%   (one on each side of the minimum); if dt < t_min there are none, and
%   the largest N with t_min(N) <= dt is Nmax. Total solution count:
%   2*Nmax + 1 (the single 0-rev solution plus a pair per feasible N) --
%   the same convention as pyKep / pumpkyn.pykep.lambert2Body.
%
%   Method per band: golden-section search for the interior minimum of
%   t(z), then bracketed bisection on each monotone side.
%
% INPUTS:
%   r1       - Initial position vector [3x1]
%   r2       - Final position vector [3x1]
%   dt       - Time of flight [scalar]
%   mu       - Gravitational parameter [scalar]
%   dir      - +1 prograde, -1 retrograde [scalar]
%   Nmax_req - Highest revolution count to attempt [scalar int]
%
% OUTPUTS:
%   V1  - Initial velocities [3 x (2*Nmax+1)]; col 1 = 0-rev, then for
%         each N the left-branch (smaller z) and right-branch solutions
%   V2  - Final velocities [3 x (2*Nmax+1)], same column convention
%   Ns  - Revolution count per column [1 x (2*Nmax+1)]
%   out - Struct: Nmax, zs (converged z per column), tmins (per band)
%
% REFERENCES:
%   [1] Curtis, "Orbital Mechanics", Sec. 5.3 (universal Lambert).
%   [2] Izzo, "Revisiting Lambert's problem", CMDA 121:1-15, 2015.

r1n = norm(r1);  r2n = norm(r2);
cross_z = r1(1)*r2(2) - r1(2)*r2(1);
dtheta  = acos(max(-1, min(1, dot(r1, r2)/(r1n*r2n))));
if (dir > 0 && cross_z < 0) || (dir < 0 && cross_z >= 0)
    dtheta = 2*pi - dtheta;
end
A = sin(dtheta) * sqrt(r1n*r2n / (1 - cos(dtheta)));
tof = @(z) lambert_tof(z, r1n, r2n, A, mu);

% --- 0-rev solution --------------------------------------------------------
[v1, v2, info0] = lambert_uv(r1, r2, dt, mu, dir);
V1 = v1;  V2 = v2;  Ns = 0;  zs = info0.z;  tmins = [];

% --- N-rev bands -----------------------------------------------------------
for N = 1:Nmax_req
    zLo = (2*pi*N)^2   + 1e-9;
    zHi = (2*pi*(N+1))^2 - 1e-9;

    % golden-section search for the interior minimum of t(z); the golden
    % ratio makes the interior points reusable -- one new evaluation per
    % iteration. (t diverges at the band edges; y stays finite there, the
    % blow-up lives in chi = sqrt(y/C) as C -> 0.)
    gr = (sqrt(5) - 1)/2;
    a = zLo; b = zHi;
    c = b - gr*(b - a);  d = a + gr*(b - a);
    fc = tof(c);  fd = tof(d);
    for k = 1:200
        if fc < fd
            b = d;  d = c;  fd = fc;
            c = b - gr*(b - a);  fc = tof(c);
        else
            a = c;  c = d;  fc = fd;
            d = a + gr*(b - a);  fd = tof(d);
        end
        if b - a < 1e-12 * b, break; end
    end
    zmin = 0.5*(a + b);
    tmin = tof(zmin);
    tmins(end+1) = tmin; %#ok<AGROW>

    if ~(tmin <= dt)                    % band infeasible: no higher N either
        break
    end

    % left branch: t decreasing on [zLo, zmin]
    zL = bisect_branch(tof, zLo, zmin, dt, -1);
    % right branch: t increasing on [zmin, zHi]
    zR = bisect_branch(tof, zmin, zHi, dt, +1);

    for z = [zL, zR]
        [~, y] = tof(z);
        f    = 1 - y/r1n;
        g    = A * sqrt(y/mu);
        gdot = 1 - y/r2n;
        V1(:, end+1) = (r2 - f*r1) / g;       %#ok<AGROW>
        V2(:, end+1) = (gdot*r2 - r1) / g;    %#ok<AGROW>
        Ns(end+1) = N;                        %#ok<AGROW>
        zs(end+1) = z;                        %#ok<AGROW>
    end
end

out = struct('Nmax', max(Ns), 'zs', zs, 'tmins', tmins, ...
             'dtheta', dtheta, 'A', A);
end

% --------------------------------------------------------------------------
function z = bisect_branch(tof, za, zb, dt, slope)
% BISECT_BRANCH  Bisection for tof(z) = dt on a monotone branch.
%   slope = +1 if t increases from za to zb, -1 if it decreases.
lo = za;  hi = zb;
for k = 1:200
    zm = 0.5*(lo + hi);
    tm = tof(zm);
    high = isnan(tm) || tm > dt;        % NaN only at band edges: treat as long
    if (slope > 0) == high
        hi = zm;
    else
        lo = zm;
    end
    if hi - lo < 1e-13 * max(1, zm), break; end
end
z = 0.5*(lo + hi);
end
