function tau = ifs_taus(g, tau0, tauf, mode)
% IFS_TAUS  Map switch-time unknowns to bounded switch times.
%   mode='sigmoid' (default): stick-breaking p_1 = s_1; p_i = p_{i-1} +
%     s_i*(1-p_{i-1}); s_i = sigmoid(g_i); tau_i = tau0 + (tauf-tau0)*p_i
%     => tau0 < tau_1 < ... < tau_k < tauf by construction. CS-safe.
%   mode='direct': the unknowns g ARE the switch times (identity pass-through);
%     monotonicity/bounds are enforced by the solver (IFS_SOLVE2 projectTau),
%     not by this map. Retires the sigmoid whose dtau/dg->0 as a gap closes
%     (GPT-5.6-sol review) -- the map is now non-degenerate everywhere.
%
% INPUTS:
%   g    - switch-time unknowns [kx1] (may be complex under complex-step):
%          stick-breaking gap-params (sigmoid) or the times themselves (direct)
%   tau0 - initial time [scalar]
%   tauf - final time [scalar] (tau0 < tauf)
%   mode - 'sigmoid'(default)|'direct' [char]
%
% OUTPUTS:
%   tau  - switch times [kx1]
%
% REFERENCES: docs/superpowers/specs/2026-07-11-ifs-design.md
if nargin < 4 || isempty(mode), mode = 'sigmoid'; end
if strcmp(mode, 'direct')
    tau = g(:);
    return;
end
g = g(:);  k = numel(g);  p = zeros(k,1);  prev = 0;
for ii = 1:k
    if real(g(ii)) >= 0
        s = 1/(1 + exp(-g(ii)));
    else
        eg = exp(g(ii));  s = eg/(1 + eg);   % CS-safe, no overflow
    end
    p(ii) = prev + s*(1 - prev);
    prev = p(ii);
end
tau = tau0 + (tauf - tau0)*p;
end
