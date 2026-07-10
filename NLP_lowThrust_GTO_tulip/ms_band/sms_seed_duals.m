function [Zseed, prob, info] = sms_seed_duals(matFile, M, epsSeed)
% SMS_SEED_DUALS  Sundman-MS seed from a direct solution's KKT-dual costates.
%
% NATIVE-domain seed: the direct solver (casadi_minfuel_sundman) already
% works in the Sundman variable tau (dt/dtau = kappa, sigma = tau/tauf0
% normalized), so its defect-constraint duals ARE the sigma-domain discrete
% costates for all 8 states [r; v; m; t] — including lamT — up to one
% positive scale beta and a global sign. No time-domain conversion (the
% (dkappa/dr)*L offset that killed the time-domain seed never appears).
% beta comes from BETA_FROM_DUALS (switching-law fit, unchanged); costates
% are lam = -beta*lamDef(1:8,:), final node padded with the last interval.
%
% The stored sigma grid is NOT uniform (measured: uniformity err 4.4e-4 vs
% mean spacing 2.5e-4 on legacy_ms_f1120), so joints are placed uniformly
% in tau = sigma*tauf0 and states/costates are sampled by linear
% interpolation on the source node grid (per-solution sampling, not
% collocation resampling; exact-node subsampling is used automatically
% wherever a joint lands on a node).
%
% INPUTS:
%   matFile - .mat with out.X [8x(N+1)], out.U [4x(N+1)], out.lamDef
%             [8xN], factor, tauf0, sigma [(N+1)x1] (direct_build_minfuel
%             save layout; tauf0/sigma reconstructed from X if absent)
%   M       - number of arcs [scalar]
%   epsSeed - smoothing value prob starts at (eps-march start) [scalar]
%
% OUTPUTS:
%   Zseed - MS unknown seed [(16M-8)x1]
%   prob  - problem struct from SMS_PROBLEM with tf, sigf, sJ set
%   info  - struct: beta, spreadPct, burnAgree, coastAgree, factorSrc,
%           lamTrelStd (rel. std of the lamT dual row — should be small:
%           lamT is a constant of the motion), node1Err (seed's node-1
%           state vs [rv0;1;0]), arcCheckErr (one mid-arc propagation at
%           eps=1e-3 vs the next joint's seeded state, rows 1:8)
%
% REFERENCES:
%   [1] sundman_minfuel/OPTIMALITY_VERIFICATION_PLAN.md, section D.
%   [2] .superpowers/sdd/task-S1-brief.md (native dual seeding).

S = load(matFile);
X = S.out.X;  U = S.out.U;  lamDef = S.out.lamDef;

prob = sms_problem(S.factor, epsSeed);

[beta, bInfo] = beta_from_duals(X, U, lamDef, prob.c);
lamNode = -beta*[lamDef(1:8, :), lamDef(1:8, end)];   % pad final node [8 x nN]

if isfield(S, 'tauf0')
    sigf = S.tauf0;
else
    r1   = sqrt(sum((X(1:3, :) - [-prob.muStar; 0; 0]).^2, 1));
    sigC = cumtrapz(X(8, :).', 1./r1(:).^prob.pSund);
    sigf = sigC(end);
end
if isfield(S, 'sigma')
    tauN = S.sigma(:).'*sigf;                          % native node taus
else
    r1   = sqrt(sum((X(1:3, :) - [-prob.muStar; 0; 0]).^2, 1));
    tauN = cumtrapz(X(8, :).', 1./r1(:).^prob.pSund).';
end

prob.sigf = sigf;
prob.sJ   = linspace(0, sigf, M+1);

Y16 = [X(1:8, :); lamNode];                            % [16 x nN]
[tU, iu] = unique(tauN);
yJ = interp1(tU.', Y16(:, iu).', prob.sJ.', 'linear').';   % [16 x (M+1)]

Zseed = sms_pack(yJ(9:16, 1), yJ(:, 2:M));

% one-arc propagation check (mid-trajectory, near-bang eps), rows 1:8
kMid  = floor(M/2);
probC = prob;  probC.epsSmooth = 1e-3;
[~, Yarc] = ode113(@(ss, y) sms_eom(ss, y, probC.Tmax, probC.c, ...
            probC.muStar, probC.epsSmooth, probC.pSund), ...
            [prob.sJ(kMid) prob.sJ(kMid+1)], yJ(:, kMid), probC.odeOpts);
info = struct('beta', beta, 'spreadPct', bInfo.spreadPct, ...
    'burnAgree', bInfo.burnAgree, 'coastAgree', bInfo.coastAgree, ...
    'factorSrc', S.factor, ...
    'lamTrelStd', std(lamDef(8, :))/abs(mean(lamDef(8, :))), ...
    'node1Err', max(abs(yJ(1:8, 1) - [prob.rv0; prob.m0; 0])), ...
    'arcCheckErr', max(abs(Yarc(end, 1:8).' - yJ(1:8, kMid+1))));
end
