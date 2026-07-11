function g = ifs_gseed(tau, tau0, tauf)
% IFS_GSEED  Inverse of IFS_TAUS: gap-params from switch times (seed-time, real).
%   g_i = logit((tau_i - tau0)/(tauf - tau0)); inputs clamped into (0,1).
%
% INPUTS:
%   tau  - switch times [kx1], in (tau0, tauf)
%   tau0 - initial time [scalar]
%   tauf - final time [scalar] (tau0 < tauf)
%
% OUTPUTS:
%   g    - unconstrained gap-params [kx1]
%
% REFERENCES: docs/superpowers/specs/2026-07-11-ifs-design.md
tau = tau(:);
p = (tau - tau0)/(tauf - tau0);
p = min(max(p, 1e-9), 1 - 1e-9);
g = log(p./(1 - p));
end
