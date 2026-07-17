function res = run_transfer_mee(cfg)
% RUN_TRANSFER_MEE  One full MEE/L-domain pipeline: seed -> homotopy ->
% report -> save. The Task-4 cross-formulation VALIDATION GATE driver: this
% is the MEE analog of run_transfer.m, mirroring its stage structure (seed,
% guarded eps:1->0 homotopy, structure report, certified-only save) but
% built on the L-domain solver (casadi_lt_mee.m / lt_mee_rhs.m / mee_seed.m,
% Tasks 1-3) instead of the reviewed Cartesian/Sundman one. The whole point
% of this run is to confirm the two independently-built formulations land on
% the SAME physical transfer (m_f, revs, switch structure) before any MEE
% ladder work is trusted.
%
% Stages: (1) seed: two-pass mee_seed call -- a cheap low-N probe first to
% learn the achieved revolution count for cfg.seedThr (ode113 integration
% cost is independent of the node-sampling count N), then a second call at
% the paper's node density (cfg.nodesPerRev per revolution) for the actual
% collocation seed; (2) guarded eps:1->0 homotopy at fixed
% tf = cfg.ctf*cfg.tfMinAnchor (the Cartesian min-time anchor value -- valid
% because the physical transfer is identical regardless of which element set
% parameterizes the collocation); (3) structure report (revs, switches,
% m_f_kg, dV_kms, edge, apogee-burn ratio, terminal inclination); (4)
% independent cross-formulation check: reconstruct the Cartesian (r,v) path
% from the converged MEE solution via elements_to_cart and confirm it
% satisfies the 2-body-plus-thrust equation of motion by finite difference.
%
% INPUTS:  cfg - .thrustN (default 10), .ctf (default 1.5), .tfMinAnchor
%          (default 22.2248 ND, the Cartesian free-longitude min-time anchor
%          for 10 N -- see run_mintime.m / results/M2_manifold.mat), .tag
%          (results filename stem, default 'MEE_M2_10N'), .seedThr (mee_seed
%          constant throttle, default 0.4 -- see DEVIATION note below),
%          .betaMode (default 'tangential'), .nodesPerRev (default 25),
%          .maxIter (default 1500), .sched (homotopy eps schedule, optional),
%          .m0kg (default 1500), .ispS (default 2000)
% OUTPUTS: res - .cfg .seed .tf .fuel .tbl .report .recon (ALWAYS returned;
%          saved to results/<tag>.mat ONLY when best.certified -- same
%          "never cache uncertified" discipline as run_transfer.m); .report
%          = .revs .switches .m_f_kg .dV_kms .edge .apoBurnRatio .incDeg
%          .LdotMin .certified; .recon = .maxResidual .nInterior (the
%          cross-formulation EOM check, Step 3 of the Task-4 brief)
%
% REFERENCES: [1] earth_elliptic_to_geo/run_transfer.m (Cartesian template).
%             [2] earth_elliptic_to_geo/homotopy_mee.m (guarded eps sweep).
%             [3] .superpowers/sdd/task-4-brief.md (this task's spec).
here   = fileparts(mfilename('fullpath'));
resDir = fullfile(here, 'results');
if ~exist(resDir, 'dir'), mkdir(resDir); end
d = @(f,v) getdef6(cfg, f, v);

thrustN    = d('thrustN', 10);
ctf        = d('ctf', 1.5);
tfMinAnchor = d('tfMinAnchor', 22.2248);   % Cartesian M2 min-time anchor [ND]
tag        = d('tag', 'MEE_M2_10N');
seedThr    = d('seedThr', 0.4);
betaMode   = d('betaMode', 'tangential');
nodesPerRev = d('nodesPerRev', 25);
maxIter    = d('maxIter', 1500);
sched      = d('sched', []);
m0kg       = d('m0kg', 1500);
ispS       = d('ispS', 2000);

p  = kepler_lt_params(thrustN, m0kg, ispS);
tf = ctf * tfMinAnchor;

% Config fingerprint (cache drift-guard): every physics/solver-relevant cfg
% field the driver reads. Saved into every cache .mat written below; on load,
% compared field-by-field against the CURRENT config so a stale cache from a
% different configuration can never be silently reused just because it
% happens to share cfg.tag (the only cache key) -- see check_cache_fp below.
% N (actual node count) is appended once known, after the probe stage.
fpBase = struct('thrustN', thrustN, 'm0kg', m0kg, 'ispS', ispS, 'ctf', ctf, ...
    'tfMinAnchor', tfMinAnchor, 'seedThr', seedThr, 'betaMode', betaMode, ...
    'nodesPerRev', nodesPerRev, 'maxIter', maxIter, 'sched', sched);

% DEVIATION (seed throttle): seedThr defaults to ~0.4, not thr=1. A thr=1
% constant burn from the paper's e=0.75 apogee GTO start crosses GEO (P=1) at
% rev~3.06 and goes coordinate-singular past ~rev 4 (Task-2 finding -- the
% MEE seed ODE chases an e->1 singularity once P overshoots the target).
% Throttling to ~0.4 (tangential steering, stopP=1.0) stretches the P=1
% crossing to ~7.5 revs, matching the paper's optimal revolution count and
% giving a well-conditioned collocation seed.

% --- stage 1: seed (two-pass: cheap revs probe, then full-density sample) --
probeFile = fullfile(resDir, [tag '_seed_probe.mat']);
if exist(probeFile, 'file')
    S = load(probeFile);  infoP = S.infoP;
    check_cache_fp(S, fpBase, probeFile, tag);
else
    optsP = struct('thr', seedThr, 'betaMode', betaMode, 'N', 50, 'stopP', 1.0);
    [~, ~, ~, ~, infoP] = mee_seed(p, optsP);
    fp = fpBase;
    save(probeFile, 'infoP', 'fp');
end
assert(infoP.nRev >= 6.5 && infoP.nRev <= 9, 'run_transfer_mee:revsOutOfRange', ...
    'seedThr=%.3f gives nRev=%.3f, outside the required [6.5,9] window -- adjust cfg.seedThr', ...
    seedThr, infoP.nRev);
N = round(nodesPerRev * infoP.nRev);
fp = fpBase;  fp.N = N;

seedFile = fullfile(resDir, [tag '_seed.mat']);
if exist(seedFile, 'file')
    S = load(seedFile);
    sigma = S.sigma;  X0 = S.X0;  U0 = S.U0;  dL0 = S.dL0;  seedInfo = S.seedInfo;
    check_cache_fp(S, fp, seedFile, tag);
else
    opts = struct('thr', seedThr, 'betaMode', betaMode, 'N', N, 'stopP', 1.0);
    [sigma, X0, U0, dL0, seedInfo] = mee_seed(p, opts);
    save(seedFile, 'sigma', 'X0', 'U0', 'dL0', 'seedInfo', 'fp');
end
x0 = X0(:,1);
fprintf(['RUN_TRANSFER_MEE %s: T=%g N, ctf=%.2f, tf=%.4f ND (%.1f h), ' ...
         'seedThr=%.2f, N=%d nodes (%d nodes/rev), seed revs=%.4f\n'], ...
        tag, thrustN, ctf, tf, tf*p.TU_s/3600, seedThr, N, nodesPerRev, seedInfo.nRev);

% --- stage 2: guarded eps 1->0 homotopy at fixed tf ------------------------
ho = struct('par', p, 'x0', x0, 'tfTarget', tf, 'maxIter', maxIter, ...
            'resDir', resDir, 'tag', tag, 'printLevel', 0, 'fp', fp);
if ~isempty(sched), ho.sched = sched; end
[best, tbl] = homotopy_mee(sigma, X0, U0, dL0, ho);

% --- stage 3: structure report ---------------------------------------------
[rr, ~] = recon_cart(best.X, sigma, best.dL, p.mu);
revs   = best.dL / (2*pi);
ss     = best.U(4,:);
nEarly = round(0.8 * numel(ss));
bMask  = ss(1:nEarly) > 0.5;
apoBurnRatio = median(rr(bMask)) / median(rr(~bMask));
Lmod   = mod(pi + sigma(:).'*best.dL, 2*pi);
burnLmod = Lmod(ss > 0.5);
apoFrac  = mean(abs(angdiff(burnLmod, pi)) < pi/4);   % fraction of burn nodes near L=pi
report = struct('revs', revs, 'switches', best.switches, 'm_f_kg', best.m_f_kg, ...
    'dV_kms', best.dV_kms, 'edge', best.edge, 'apoBurnRatio', apoBurnRatio, ...
    'apoFrac', apoFrac, 'incDeg', best.incDeg, 'LdotMin', best.LdotMin, ...
    'defect', best.maxDefect, 'termErr', best.termErr, ...
    'ipoptStatus', best.ipoptStatus, 'certified', best.certified);

% --- stage 4: cross-formulation reconstruction check ------------------------
recon = check_reconstruction(best.X, best.U, best.dL, sigma, p);

res = struct('cfg', cfg, 'seed', seedInfo, 'tf', tf, 'fuel', best, ...
             'tbl', tbl, 'report', report, 'recon', recon, 'sigma', sigma, 'fp', fp);
if best.certified
    save(fullfile(resDir, [tag '.mat']), 'res');
else
    warning('run_transfer_mee:uncertified', ['%s: NOT saved (certified=0, ' ...
        'defect=%.2e) -- campaign rule: never cache uncertified results'], ...
        tag, best.maxDefect);
end
fprintf(['DONE %s: certified=%d revs=%.3f sw=%d edge=%.1f%% mf=%.2f kg ' ...
         'dV=%.3f km/s apoBurn=%.2f apoFrac=%.2f incDeg=%.4f LdotMin=%.3e ' ...
         'reconMaxRes=%.3e\n'], tag, report.certified, revs, best.switches, ...
         100*best.edge, best.m_f_kg, best.dV_kms, apoBurnRatio, apoFrac, ...
         best.incDeg, best.LdotMin, recon.maxResidual);
end

% ---------------------------------------------------------------------------
function v = getdef6(s, f, dflt)
% GETDEF6  Optional-field default (mirrors casadi_lt_2body's helper).
if isfield(s, f) && ~isempty(s.(f)), v = s.(f); else, v = dflt; end
end

% ---------------------------------------------------------------------------
function check_cache_fp(S, fp, file, tag)
% CHECK_CACHE_FP  Fail-loud cache-fingerprint guard. If loaded cache struct S
% carries a .fp field, compare it field-by-field against the current config
% fingerprint fp and error, naming the first mismatched field and the
% offending file, on any disagreement -- a stale cache built under a
% different thrustN/ctf/seedThr/... must never be silently reused just
% because it happens to share cfg.tag (the only cache key). BACKWARD COMPAT:
% a pre-fix cache file with NO .fp field (e.g. the MEE_M2_10N certified run
% predating this guard) only WARNs and is trusted as-is -- cfg.tag is already
% an exact match (that's how the file was found by name), and no per-field
% comparison is possible without a stored fingerprint.
if ~isfield(S, 'fp')
    warning('run_transfer_mee:noCachedFingerprint', ['%s has no cached ' ...
        'config fingerprint (pre-fix cache) -- trusting it because ' ...
        'tag=''%s'' matches; use a new cfg.tag to regain fingerprint ' ...
        'protection for this run'], file, tag);
    return;
end
flds = fieldnames(fp);
for k = 1:numel(flds)
    f = flds{k};
    if ~isfield(S.fp, f) || ~isequal(S.fp.(f), fp.(f))
        error('run_transfer_mee:fingerprintMismatch', ['cached config ' ...
            'fingerprint mismatch in %s: field ''%s'' differs between the ' ...
            'cache and the current cfg -- stale cache from a different ' ...
            'configuration under the same tag=''%s''; delete the file or ' ...
            'use a new cfg.tag'], file, f, tag);
    end
end
end

% ---------------------------------------------------------------------------
function [rr, vv] = recon_cart(X, sigma, dL, mu)
% RECON_CART  Reconstruct inertial (r,v) at every node from the converged
% MEE solution via elements_to_cart. rr is |r| [1xN]; vv is [3xN] velocity.
Nn = size(X, 2);
rr = zeros(1, Nn);  vv = zeros(3, Nn);
for k = 1:Nn
    Lk = pi + sigma(k)*dL;
    [rk, vk] = elements_to_cart(X(1,k), X(2,k), X(3,k), X(4,k), X(5,k), Lk, mu);
    rr(k) = norm(rk);  vv(:,k) = vk;
end
end

% ---------------------------------------------------------------------------
function d = angdiff(a, b)
% ANGDIFF  Signed shortest angular difference a-b, wrapped to (-pi,pi].
d = mod(a - b + pi, 2*pi) - pi;
end

% ---------------------------------------------------------------------------
function recon = check_reconstruction(X, U, dL, sigma, par)
% CHECK_RECONSTRUCTION  Independent cross-formulation validation (Task-4
% brief Step 3, item 4): map the converged MEE solution through
% elements_to_cart, then confirm the reconstructed inertial (r,v) path
% satisfies the 2-body-plus-thrust equation of motion,
%     vdot = -mu*r/|r|^3 + (Tmax/m)*thr*(RTN thrust direction in inertial),
% via central finite differencing of v against the physical time (X(7,:)).
% This is a genuinely independent check: elements_to_cart is algebraic
% (no ODE/collocation machinery shared with lt_mee_rhs.m or casadi_lt_mee.m),
% so a defect-free collocation solve that were nonetheless dynamically wrong
% (e.g. a sign error in the Gauss-to-Cartesian mapping) would NOT be caught
% by casadi_lt_mee's own internal defect re-check -- only by this.
%
% INPUTS:  X [7x(N+1)] converged MEE states; U [4x(N+1)] controls; dL scalar
%          converged span; sigma [(N+1)x1]; par - kepler_lt_params struct
% OUTPUTS: recon - struct .maxResidual (max over interior nodes of
%          |vdot_FD - vdot_model|, ND accel units) .nInterior (node count
%          used) .residuals [1x(N-2)] per-node values (FD-limited, NOT
%          gated at machine precision -- report, don't gate, per the brief)
Nn = size(X, 2);
mu = par.mu;  Tm = par.Tmax;
r = zeros(3, Nn);  v = zeros(3, Nn);
for k = 1:Nn
    Lk = pi + sigma(k)*dL;
    [r(:,k), v(:,k)] = elements_to_cart(X(1,k), X(2,k), X(3,k), X(4,k), X(5,k), Lk, mu);
end
t = X(7,:);
nInt = Nn - 2;
residuals = zeros(1, nInt);
for k = 2:Nn-1
    vdotFD = (v(:,k+1) - v(:,k-1)) / (t(k+1) - t(k-1));    % central FD, nonuniform-safe
    rk   = r(:,k);  rmag = norm(rk);
    m    = X(6,k);  thr = U(4,k);  beta = U(1:3,k);
    rhat = rk / rmag;
    hvec = cross(rk, v(:,k));  nhat = hvec / norm(hvec);
    that = cross(nhat, rhat);
    Rrtn2eci = [rhat, that, nhat];         % columns = radial, transverse, normal (inertial)
    accGrav   = -mu * rk / rmag^3;
    accThrust = (Tm/m) * thr * (Rrtn2eci * beta);
    residuals(k-1) = norm(vdotFD - (accGrav + accThrust));
end
recon = struct('maxResidual', max(residuals), 'nInterior', nInt, 'residuals', residuals);
end
