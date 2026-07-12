function g = ifs_gseed(tau, tau0, tauf)
% IFS_GSEED  Inverse of IFS_TAUS (stick-breaking): gap-params from ordered switch
%   times. frac_i = (p_i - p_{i-1})/(1 - p_{i-1}); g_i = logit(frac_i). Seed-time
%   (real); clamps inputs strictly inside (0,1) and enforces strict increase.
%
% INPUTS:
%   tau  - switch times [kx1], strictly increasing in (tau0, tauf)
%   tau0 - initial time [scalar]
%   tauf - final time [scalar] (tau0 < tauf)
%
% OUTPUTS:
%   g    - unconstrained gap-params [kx1]
%
% REFERENCES: docs/superpowers/specs/2026-07-11-ifs-design.md
tau = tau(:);  k = numel(tau);
p = (tau - tau0)/(tauf - tau0);
p = min(max(p, 1e-9), 1 - 1e-9);
for ii = 2:k
    if p(ii) <= p(ii-1), p(ii) = p(ii-1) + 1e-9; end
end
g = zeros(k,1);  prev = 0;
for ii = 1:k
    frac = (p(ii) - prev)/(1 - prev);
    frac = min(max(frac, 1e-12), 1 - 1e-12);
    g(ii) = log(frac/(1 - frac));
    prev = p(ii);
end
end
