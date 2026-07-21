function [score, tauSwitch, diag] = pmp_refine_indicator(seedFile, opts)
% PMP_REFINE_INDICATOR  Per-interval mesh-refinement score from the PMP
%   switching function, plus passive Hamiltonian-residual diagnostics.
%
% Recovers the discrete costates from a direct solution's KKT defect duals
% (mode-'d' midpoint map, ms_band/sms_seed_duals), forms the min-fuel
% switching function S(tau) = 1 - ||lamV||*c/m - lamM on the node grid,
% and scores each collocation interval by how poorly the mesh localizes an
% S=0 crossing (a switch): wide intervals with a central crossing score
% highest. Also returns, as PASSIVE diagnostics only, the Sundman-domain
% Hamiltonian residual |kappa*(Ht+lamT)| per node and the count of
% switching-law sign violations outside a switch deadband.
%
% INPUTS:
%   seedFile - .mat with out.X [8x(N+1)], out.U [4x(N+1)], out.lamDef [8xN],
%              factor, tauf0, sigma [(N+1)x1] (sms_seed_duals input layout)
%   opts     - struct: M arcs for dual map [default 40], epsEval smoothing
%              [default 1e-4], mode dual->costate map [default 'd'], nbr
%              neighbor half-window for scoring [default 3] (also sets the
%              nViol switch-deadband half-width)
%
% OUTPUTS:
%   score     - per-interval refinement score [1xN], >= 0
%   tauSwitch - direct-throttle switch times in tau, sorted [1xnsw]
%   diag      - struct: Snode [1xnN], tauN [1xnN], tauCr [1xncr] (S=0
%               crossings), swI [1xnsw], Hres [1xnN], HresMax, nViol,
%               betaSpread, beta
%
% REFERENCES:
%   [1] docs/superpowers/specs/2026-07-11-pmp-mesh-refine-design.md
%   [2] ms_band/verify_direct_pmp.m (Snode + dual-map reuse)

if nargin < 2, opts = struct(); end
if ~isfield(opts, 'M'),       opts.M = 40;        end
if ~isfield(opts, 'epsEval'), opts.epsEval = 1e-4; end
if ~isfield(opts, 'mode'),    opts.mode = 'd';    end
if ~isfield(opts, 'nbr'),     opts.nbr = 3;       end

[~, prob, info] = sms_seed_duals(seedFile, opts.M, opts.epsEval, opts.mode);
tauN = info.tauN;   Y16 = info.Y16;   X = info.X;   U = info.U;
nN = size(X, 2);    N = nN - 1;

% switching function on the node grid (same construction as verify_direct_pmp)
lamV  = Y16(12:14, :);   lamM = Y16(15, :);   lamT = Y16(16, :);
Snode = 1 - sqrt(sum(lamV.^2, 1)).*prob.c./X(7, :) - lamM;

% S=0 crossings, localized by linear interpolation within each bracket
crossI = find(diff(sign(Snode)) ~= 0);        % node index just before crossing
tauCr  = zeros(1, numel(crossI));
for q = 1:numel(crossI)
    k = crossI(q);
    tauCr(q) = tauN(k) + (0 - Snode(k))*(tauN(k+1) - tauN(k))/(Snode(k+1) - Snode(k));
end

% per-interval score: for each crossing, spread weight over +-nbr intervals,
% weighted by (normalized local width) * (centrality of the crossing in [0,0.5])
score = zeros(1, N);
hInt  = diff(tauN);                            % [1xN] interval widths in tau
sigf  = tauN(end);
for q = 1:numel(crossI)
    kc  = crossI(q);
    off = min(tauCr(q) - tauN(kc), tauN(kc+1) - tauCr(q)) / (tauN(kc+1) - tauN(kc));
    for kk = max(1, kc-opts.nbr):min(N, kc+opts.nbr)
        score(kk) = score(kk) + (hInt(kk)/sigf)*off;
    end
end

% direct-throttle switch times
s     = U(4, :);
swI   = find(diff(double(s > 0.5)) ~= 0);      % switch interval indices [1xnsw]
tauSwitch = sort((tauN(swI) + tauN(swI+1))/2);

% ---- PASSIVE diagnostics ---------------------------------------------------
rE   = [-prob.muStar; 0; 0];
Hres = zeros(1, nN);
for k = 1:nN
    [~, Htk] = sms_eom(0, Y16(:, k), prob.Tmax, prob.c, prob.muStar, ...
                       opts.epsEval, prob.pSund);
    r1 = sqrt(sum((X(1:3, k) - rE).^2));
    Hres(k) = abs(r1^prob.pSund * (Htk + lamT(k)));
end
% switching-law sign violations outside a +-3-node deadband of a direct switch
viol = sign(Snode) ~= sign(0.5 - s);
dead = false(1, nN);
for w = -opts.nbr:opts.nbr
    idx = swI + w;  idx = idx(idx >= 1 & idx <= nN);  dead(idx) = true;
end
nViol = nnz(viol & ~dead);

diag = struct('Snode', Snode, 'tauN', tauN, 'tauCr', tauCr, 'swI', swI, ...
              'Hres', Hres, 'HresMax', max(Hres), 'nViol', nViol, ...
              'betaSpread', info.spreadPct, 'beta', info.beta);
end
