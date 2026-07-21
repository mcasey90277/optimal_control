function [C, S] = stumpff(z)
% STUMPFF  Stumpff functions C(z) and S(z) for universal variables.
%
%   C(z) = (1 - cos(sqrt(z)))/z                 (z > 0, elliptic)
%        = (cosh(sqrt(-z)) - 1)/(-z)            (z < 0, hyperbolic)
%   S(z) = (sqrt(z) - sin(sqrt(z)))/z^(3/2)     (z > 0)
%        = (sinh(sqrt(-z)) - sqrt(-z))/(-z)^(3/2) (z < 0)
%
%   Near z = 0 both closed forms are 0/0; a 4-term Taylor series is used
%   for |z| < 1e-4 (relative error < 1e-19 there):
%     C = 1/2 - z/24 + z^2/720 - z^3/40320
%     S = 1/6 - z/120 + z^2/5040 - z^3/362880
%
% INPUTS:
%   z - Universal-variable argument [scalar or any array]
%
% OUTPUTS:
%   C - Stumpff C(z) [same size as z]
%   S - Stumpff S(z) [same size as z]
%
% REFERENCES:
%   [1] Bate, Mueller, White, "Fundamentals of Astrodynamics", Ch. 4.
%   [2] Curtis, "Orbital Mechanics for Engineering Students", Sec. 3.5.

C = zeros(size(z));
S = zeros(size(z));

small = abs(z) < 1e-4;
zs = z(small);
C(small) = 1/2 - zs/24 + zs.^2/720 - zs.^3/40320;
S(small) = 1/6 - zs/120 + zs.^2/5040 - zs.^3/362880;

pos = z >= 1e-4;
zp  = z(pos);
sz  = sqrt(zp);
C(pos) = (1 - cos(sz)) ./ zp;
S(pos) = (sz - sin(sz)) ./ zp.^1.5;

neg = z <= -1e-4;
zn  = -z(neg);                       % zn > 0
sn  = sqrt(zn);
C(neg) = (cosh(sn) - 1) ./ zn;
S(neg) = (sinh(sn) - sn) ./ zn.^1.5;
end
