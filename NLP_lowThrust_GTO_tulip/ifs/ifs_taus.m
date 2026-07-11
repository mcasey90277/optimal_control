function tau = ifs_taus(g, tau0, tauf)
% IFS_TAUS  Map unconstrained gap-params to bounded switch times (CS-safe sigmoid).
%   tau_i = tau0 + (tauf-tau0)*sigmoid(g_i), each strictly in (tau0, tauf).
%   Independent per switch (block-sparsity preserving); ordering not enforced.
%
% INPUTS:
%   g    - unconstrained gap-params [kx1] (may be complex under complex-step)
%   tau0 - initial time [scalar]
%   tauf - final time [scalar] (tau0 < tauf)
%
% OUTPUTS:
%   tau  - switch times [kx1], each strictly in (tau0, tauf)
%
% REFERENCES: docs/superpowers/specs/2026-07-11-ifs-design.md
g = g(:);  sig = zeros(size(g));
for ii = 1:numel(g)
    if real(g(ii)) >= 0
        sig(ii) = 1/(1 + exp(-g(ii)));
    else
        eg = exp(g(ii));  sig(ii) = eg/(1 + eg);   % CS-safe, no overflow
    end
end
tau = tau0 + (tauf - tau0)*sig;
end
