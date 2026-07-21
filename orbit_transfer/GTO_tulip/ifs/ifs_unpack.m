function [lam0, N, tau] = ifs_unpack(Z, k)
% IFS_UNPACK  Split the IFS unknown vector.
%
% INPUTS:
%   Z - unknown vector [(8+17k) x 1]
%   k - number of switches [scalar]
% OUTPUTS:
%   lam0 - initial costate [8x1]
%   N    - node states [16 x k]
%   tau  - switch times [k x 1]
%
% REFERENCES: docs/superpowers/specs/2026-07-11-ifs-design.md
lam0 = Z(1:8);
N    = reshape(Z(8 + (1:16*k)), 16, k);
tau  = Z(8 + 16*k + (1:k));
tau  = tau(:);
end
