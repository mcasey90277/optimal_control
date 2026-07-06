function [t, y] = lambert_tof(z, r1n, r2n, A, mu)
% LAMBERT_TOF  Time of flight as a function of the universal variable z.
%
%   The heart of the universal-variables Lambert formulation: for fixed
%   endpoint geometry (r1n, r2n, A), each z corresponds to one conic
%   through both endpoints, and its transfer time is
%
%     y(z)   = r1n + r2n + A*(z*S(z) - 1)/sqrt(C(z))
%     chi    = sqrt(y/C)
%     sqrt(mu)*t = chi^3 * S(z) + A*sqrt(y)
%
%   z < 0 hyperbolic, z = 0 parabolic, 0 < z < (2*pi)^2 elliptic
%   (single rev), ((2*pi*N)^2, (2*pi*(N+1))^2) the N-rev elliptic band.
%   Where y(z) < 0 the geometry is infeasible and t = NaN is returned.
%
% INPUTS:
%   z    - Universal-variable argument [scalar or vector]
%   r1n  - Norm of initial position [scalar]
%   r2n  - Norm of final position [scalar]
%   A    - Geometry constant sin(dth)*sqrt(r1n*r2n/(1-cos(dth))) [scalar]
%   mu   - Gravitational parameter [scalar]
%
% OUTPUTS:
%   t    - Time of flight for each z [same size as z]; NaN where invalid
%   y    - The y(z) intermediate [same size as z]
%
% REFERENCES:
%   [1] Curtis, "Orbital Mechanics for Engineering Students", Alg. 5.2.

[C, S] = stumpff(z);
y = r1n + r2n + A * (z .* S - 1) ./ sqrt(C);

t = nan(size(z));
ok = y >= 0 & C > 0;
chi = sqrt(y(ok) ./ C(ok));
t(ok) = (chi.^3 .* S(ok) + A * sqrt(y(ok))) / sqrt(mu);
end
