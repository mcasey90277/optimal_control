function [Zseed, prob, info] = sms_seed_duals(matFile, M, epsSeed, mode)
% SMS_SEED_DUALS  Sundman-MS seed from a direct solution's KKT-dual costates.
%
% NATIVE-domain seed: the direct solver (casadi_minfuel_sundman) works in
% the Sundman variable tau (dt/dtau = kappa, sigma = tau/tauf0), so its
% defect-constraint duals are the sigma-domain discrete costates for all
% 8 states [r; v; m; t] — including lamT — up to one positive scale beta
% and a global sign. HOW the per-INTERVAL duals map to per-NODE costates
% is the delicate step (GPT-5.6 review, 2026-07-10): four candidate
% conversions are implemented, selected by `mode`:
%
%   'a' (baseline)  node k <- -beta*lamDef(:,k) (left interval), final
%                   node duplicated. The original Task-6 convention.
%   'b' (h-weight)  node k <- -beta_b*lamDef(:,k)/h_k, h_k = dsig(k)
%                   (nonuniform normalized interval widths), beta_b refit
%                   after weighting; node placement as 'a'. Tests the
%                   "duals carry local trapezoid weights" hypothesis.
%   'c' (adjacent-h) interior node k <- -beta*(h_{k-1}*lamDef(:,k-1)
%                   + h_k*lamDef(:,k))/(h_{k-1}+h_k); one-sided at the
%                   ends (trapezoid-adjoint averaging convention).
%   'd' (midpoint)  PRINCIPLED map for THIS transcription. KKT
%                   stationarity of casadi_minfuel_sundman w.r.t. X_k:
%                   defect D_k = X_{k+1}-X_k-tauf*(h_k/2)(F_k+F_{k+1})
%                   couples mu_{k-1}, mu_k to node k with weights
%                   (h_{k-1}, h_k)/2, and the objective quadrature carries
%                   the MATCHING node weight tauf*(h_{k-1}+h_k)/2 — the
%                   h's cancel, leaving mu_k ~ -beta*lambda(tau at the
%                   MIDPOINT of interval k) with NO weight. So: place
%                   -beta*lamDef at interval midpoints and linearly
%                   interpolate onto the nodes (linear extrapolation at
%                   the two boundary nodes). The baseline 'a' is this map
%                   with a half-interval shift — a trajectory-shaped
%                   O(h*dlam/dtau) error, largest where costates move
%                   fast (perigee).
%
% beta comes from BETA_FROM_DUALS (switching-law fit on the interval
% duals; those are midpoint quantities, consistent with the mMid the fit
% uses). Modes a,c,d share that beta (their maps are scale-preserving);
% mode b refits on the weighted duals.
%
% The stored sigma grid is validated (tau(1) = 0, strictly increasing,
% tau(end) = sigf) and is NOT uniform on the legacy files, so MS joints
% are placed uniformly in tau and states/costates are sampled by linear
% interpolation on the source node grid.
%
% INPUTS:
%   matFile - .mat with out.X [8x(N+1)], out.U [4x(N+1)], out.lamDef
%             [8xN], factor, tauf0, sigma [(N+1)x1] (direct_build_minfuel
%             save layout; tauf0/sigma reconstructed from X if absent)
%   M       - number of arcs [scalar]
%   epsSeed - smoothing value prob starts at [scalar]
%   mode    - dual->costate conversion: 'a'|'b'|'c'|'d' [char, default 'a']
%
% OUTPUTS:
%   Zseed - MS unknown seed [(16M-8)x1]
%   prob  - problem struct from SMS_PROBLEM with tf, sigf, sJ set
%   info  - struct: mode, beta, spreadPct, burnAgree, coastAgree,
%           factorSrc, lamTrelStd, node1Err, arcCheckErr, and the raw
%           per-node data for validation harnesses: tauN [1xnN], Y16
%           [16xnN] (direct states + candidate costates at nodes),
%           X [8xnN], U [4xnN]
%
% REFERENCES:
%   [1] sundman_minfuel/OPTIMALITY_VERIFICATION_PLAN.md, section D.
%   [2] .superpowers/sdd/task-S1-brief.md (native dual seeding).
%   [3] .superpowers/sdd/gpt56_review_S1.md (dual-map candidates).

if nargin < 4 || isempty(mode), mode = 'a'; end

S = load(matFile);
X = S.out.X;  U = S.out.U;  lamDef = S.out.lamDef;
nN = size(X, 2);  N = nN - 1;

prob = sms_problem(S.factor, epsSeed);

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
% grid validation (review robustness item)
if abs(tauN(1)) > 1e-12*sigf || any(diff(tauN) <= 0) ...
        || abs(tauN(end) - sigf) > 1e-9*sigf
    error('sms_seed_duals:grid', ...
          'tau grid invalid: tau(1)=%.3e, tau(end)-sigf=%.3e, monotone=%d', ...
          tauN(1), tauN(end) - sigf, all(diff(tauN) > 0));
end
h = diff(tauN)/sigf;                                   % normalized widths [1xN]

% ---- candidate dual -> node-costate maps -----------------------------------
switch mode
    case 'a'
        [beta, bInfo] = beta_from_duals(X, U, lamDef, prob.c);
        lamNode = -beta*[lamDef(1:8, :), lamDef(1:8, end)];
    case 'b'
        lamDefW = lamDef./h;                           % h broadcast over rows
        [beta, bInfo] = beta_from_duals(X, U, lamDefW, prob.c);
        lamNode = -beta*[lamDefW(1:8, :), lamDefW(1:8, end)];
    case 'c'
        [beta, bInfo] = beta_from_duals(X, U, lamDef, prob.c);
        lamNode = zeros(8, nN);
        lamNode(:, 1)  = lamDef(1:8, 1);
        lamNode(:, nN) = lamDef(1:8, N);
        for k = 2:N
            lamNode(:, k) = (h(k-1)*lamDef(1:8, k-1) + h(k)*lamDef(1:8, k)) ...
                            /(h(k-1) + h(k));
        end
        lamNode = -beta*lamNode;
    case 'd'
        [beta, bInfo] = beta_from_duals(X, U, lamDef, prob.c);
        tauMid  = (tauN(1:end-1) + tauN(2:end))/2;     % [1xN]
        lamNode = -beta*interp1(tauMid.', lamDef(1:8, :).', tauN.', ...
                                'linear', 'extrap').';
    otherwise
        error('sms_seed_duals:mode', 'unknown mode %s', mode);
end

prob.sigf = sigf;
prob.sJ   = linspace(0, sigf, M+1);

Y16 = [X(1:8, :); lamNode];                            % [16 x nN]
yJ  = interp1(tauN.', Y16.', prob.sJ.', 'linear').';   % [16 x (M+1)]

Zseed = sms_pack(yJ(9:16, 1), yJ(:, 2:M));

% one-arc propagation check (mid-trajectory, near-bang eps), rows 1:8
kMid  = floor(M/2);
probC = prob;  probC.epsSmooth = 1e-3;
[~, Yarc] = ode113(@(ss, y) sms_eom(ss, y, probC.Tmax, probC.c, ...
            probC.muStar, probC.epsSmooth, probC.pSund), ...
            [prob.sJ(kMid) prob.sJ(kMid+1)], yJ(:, kMid), probC.odeOpts);
info = struct('mode', mode, 'beta', beta, 'spreadPct', bInfo.spreadPct, ...
    'burnAgree', bInfo.burnAgree, 'coastAgree', bInfo.coastAgree, ...
    'factorSrc', S.factor, ...
    'lamTrelStd', std(lamDef(8, :))/abs(mean(lamDef(8, :))), ...
    'node1Err', max(abs(yJ(1:8, 1) - [prob.rv0; prob.m0; 0])), ...
    'arcCheckErr', max(abs(Yarc(end, 1:8).' - yJ(1:8, kMid+1))), ...
    'tauN', tauN, 'Y16', Y16, 'X', X, 'U', U);
end
