function Z = ifs_pack(lam0, N, tau)
% IFS_PACK  Assemble the IFS unknown vector.
%
% INPUTS:
%   lam0 - initial costate [8x1]
%   N    - node augmented states at each switch [16 x k]
%   tau  - switch times [k x 1]
% OUTPUTS:
%   Z    - unknown vector [(8+17k) x 1]
%
% REFERENCES: docs/superpowers/specs/2026-07-11-ifs-design.md
Z = [lam0(:); N(:); tau(:)];
end
