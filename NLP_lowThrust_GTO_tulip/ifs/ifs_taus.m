function tau = ifs_taus(g, tau0, tauf)
% IFS_TAUS  Map unconstrained gap-params to MONOTONE bounded switch times.
%   Stick-breaking: p_1 = s_1; p_i = p_{i-1} + s_i*(1-p_{i-1}); s_i = sigmoid(g_i);
%   tau_i = tau0 + (tauf-tau0)*p_i  =>  tau0 < tau_1 < ... < tau_k < tauf.
%   CS-safe (branch on real(g), no overflow).
%
% INPUTS:
%   g    - unconstrained gap-params [kx1] (may be complex under complex-step)
%   tau0 - initial time [scalar]
%   tauf - final time [scalar] (tau0 < tauf)
%
% OUTPUTS:
%   tau  - switch times [kx1], strictly increasing in (tau0, tauf)
%
% REFERENCES: docs/superpowers/specs/2026-07-11-ifs-design.md
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
