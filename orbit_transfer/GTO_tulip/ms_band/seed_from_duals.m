function [Zseed, tJ, info] = seed_from_duals(matFile, factorTarget, M)
% SEED_FROM_DUALS  MS seed from a direct Sundman solution's KKT-dual costates.
%
% Loads a direct min-fuel solution (X, U, out.lamDef), converts duals to
% physical time-domain costates lam = -beta*lamDef (beta from the switching
% -law fit), linearly rescales the time axis source-tf -> target-tf, places
% M tau-uniform arc joints, and samples state+costate at the joints ON the
% source solution's own node grid (no collocation resampling). Verifies the
% costate scale by propagating one mid-trajectory arc and comparing end
% states (info.arcCheckErr; plan-doc section F.6 route).
%
% INPUTS:
%   matFile      - .mat with out.X [8xnN], out.U [4xnN], factor, out.lamDef
%                  [8x(nN-1)] (resolved on the path; NOTE: X/U/lamDef are
%                  nested under out in this campaign's saved files, not
%                  top-level -- verified directly against
%                  legacy_ms_f1120.mat)
%   factorTarget - target tf factor [scalar]
%   M            - number of arcs [scalar]
%
% OUTPUTS:
%   Zseed - MS unknown seed [(14M-7)x1]
%   tJ    - joint times at the TARGET tf [1x(M+1)]
%   info  - struct: beta, spreadPct, burnAgree, coastAgree, factorSrc,
%           arcCheckErr

S = load(matFile);
X = S.out.X;  U = S.out.U;  lamDef = S.out.lamDef;
p = cr3bp_lt_params(0.025, 15, 2100);

[beta, bInfo] = beta_from_duals(X, U, lamDef, p.c);
lamNode = -beta*[lamDef(1:7, :), lamDef(1:7, end)];   % pad final node [7x nN]

tSrc  = X(8, :);                                       % physical node times
scale = (factorTarget*6.290694)/tSrc(end);
tScl  = tSrc*scale;
yNode = [X(1:7, :); lamNode];                          % [14 x nN]

tJ = arc_boundaries_tau(tScl, X(1:3, :), M, p.muStar);
% sample at joints on the source grid (dedupe any repeated node times)
[tU, iu] = unique(tScl);
yJ = interp1(tU.', yNode(:, iu).', tJ.', 'linear').';  % [14 x (M+1)]
Zseed = ms_pack(yJ(8:14, 1), yJ(:, 2:M));

% one-arc scale verification at the SOURCE tf (mid-trajectory arc)
prob0 = ms_problem(S.factor, 1e-3);
kMid  = floor(M/2);
tJsrc = arc_boundaries_tau(tSrc, X(1:3, :), M, p.muStar);
ySrc  = interp1(unique(tSrc).', yNode(:, iu).', tJsrc.', 'linear').';
[~, Yarc] = ode113(@(t, y) lt_pmp_eom_minfuel(t, y, prob0.Tmax, prob0.c, ...
            prob0.muStar, prob0.epsSmooth), ...
            [tJsrc(kMid) tJsrc(kMid+1)], ySrc(:, kMid), prob0.odeOpts);
info = struct('beta', beta, 'spreadPct', bInfo.spreadPct, ...
    'burnAgree', bInfo.burnAgree, 'coastAgree', bInfo.coastAgree, ...
    'factorSrc', S.factor, ...
    'arcCheckErr', max(abs(Yarc(end, 1:7).' - ySrc(1:7, kMid+1))));
end
