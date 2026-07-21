function g = ifs_gseed(tau, tau0, tauf, mode)
% IFS_GSEED  Inverse of IFS_TAUS: switch-time unknowns from ordered switch times.
%   mode='sigmoid' (default): stick-breaking inverse. frac_i = (p_i -
%     p_{i-1})/(1 - p_{i-1}); g_i = logit(frac_i). Clamps inside (0,1),
%     enforces strict increase.
%   mode='direct': identity -- the unknowns ARE the switch times.
%
% INPUTS:
%   tau  - switch times [kx1], strictly increasing in (tau0, tauf)
%   tau0 - initial time [scalar]
%   tauf - final time [scalar] (tau0 < tauf)
%   mode - 'sigmoid'(default)|'direct' [char]
%
% OUTPUTS:
%   g    - switch-time unknowns [kx1]
%
% REFERENCES: docs/superpowers/specs/2026-07-11-ifs-design.md
if nargin < 4 || isempty(mode), mode = 'sigmoid'; end
if strcmp(mode, 'direct')
    g = tau(:);
    return;
end
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
