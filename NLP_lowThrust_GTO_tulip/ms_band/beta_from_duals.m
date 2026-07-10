function [beta, info] = beta_from_duals(X, U, lamDef, c)
% BETA_FROM_DUALS  Empirical costate scale from KKT duals via the switching law.
%
% The discrete costates lamDef are the physical costates up to one positive
% scale beta and a global sign (-1, verified). Writing S = 1 - beta*W with
% W = ||lam_v||c/m - lam_m (sign -1 folded in), beta is pinned by least
% squares on S=0 at the throttle-switch intervals. The per-switch spread of
% 1/W doubles as a quality metric (~1% on certified solutions).
%
% INPUTS:
%   X      - direct-solution states [8x(N+1)] (row 7 mass, row 8 time)
%   U      - direct-solution controls [4x(N+1)] (row 4 throttle)
%   lamDef - KKT-dual discrete costates [8xN] (per interval)
%   c      - ND exhaust velocity [scalar]
%
% OUTPUTS:
%   beta - positive scale: physical costate = -beta*lamDef(1:7,:) [scalar]
%   info - struct: spreadPct, burnAgree, coastAgree, S [1xN], W [1xN]
%
% REFERENCES:
%   [1] sundman_minfuel/OPTIMALITY_VERIFICATION_PLAN.md, section D.

N     = size(lamDef, 2);
s     = U(4, :);
mMid  = 0.5*(X(7, 1:end-1) + X(7, 2:end));
nLamV = sqrt(sum(lamDef(4:6, :).^2, 1));
W     = nLamV.*c./mMid + (-1)*lamDef(7, :);
swI   = min(find(diff(double(s > 0.5)) ~= 0), N);
beta  = sum(W(swI))/sum(W(swI).^2);
S     = 1 - beta*W;
burnI  = (s(1:end-1) > 0.5) & (s(2:end) > 0.5);
coastI = (s(1:end-1) < 0.5) & (s(2:end) < 0.5);
info = struct('spreadPct', 100*std(1./W(swI))/abs(mean(1./W(swI))), ...
              'burnAgree', mean(S(burnI) < 0), ...
              'coastAgree', mean(S(coastI) > 0), 'S', S, 'W', W);
end
