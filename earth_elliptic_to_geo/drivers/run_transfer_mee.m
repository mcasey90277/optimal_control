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
%          (results filename stem, default mee_fuel_tag(cfg.thrustN) --
%          equals 'MEE_M2_10N' at the default thrustN=10, final-review Fix 2:
%          was a hardcoded 'MEE_M2_10N' literal regardless of thrustN before
%          this fix), .seedThr (mee_seed
%          constant throttle, default 0.4 -- see DEVIATION note below),
%          .betaMode (default 'tangential'), .nodesPerRev (default 25),
%          .maxIter (default 1500), .sched (homotopy eps schedule, optional),
%          .m0kg (default 1500), .ispS (default 2000), .warmStart (Task 7
%          EXTERNAL WARM-START ENTRY POINT, optional -- struct .sigma
%          [(Mp+1)x1] .X [7x(Mp+1)] .U [4x(Mp+1)] .dL [scalar], the PREVIOUS
%          ladder rung's converged fuel solution, with .dL ALREADY C-law
%          rescaled by the caller (run_ladder.m: dL_guess =
%          dL_prev*(T_prev/T_new)) -- when supplied, Stage 1's cold
%          constant-throttle mee_seed pass is skipped entirely in favor of
%          mesh-refining this trajectory onto the new rung's own node grid
%          (interp_warmstart.m, factored out of what was formerly inline
%          interp1 code here -- linear for X and U's RTN thrust-direction
%          rows, nearest-plus-renormalize for U's throttle/beta rows; N
%          sized from warmStart.dL/(2*pi) at cfg.nodesPerRev), and
%          Stage 2 tries a SINGLE direct eps=0 solve first (same-basin
%          refinement precedent, mirroring solve_warm_node) before falling
%          back to a short eps continuation tail (cfg.warmFallbackSched,
%          default [0.3 0.2 0.12 0.07 0.04 0.025 0.015 0.008 0.004 0.002
%          0.001 0]) only if the direct attempt does not certify. .xf (Task
%          4, [5x1] terminal target [P;ex;ey;hx;hy], default [1;0;0;0;0] =
%          GEO -- forwarded into mee_seed's stopP (stopP = xf(1), so the seed
%          integrates to the CUSTOM final radius rather than the hardcoded
%          1.0), the homotopy_mee .xf opts field, and the warm-direct
%          casadi_lt_mee call), .initElems (Task 4, [7x1]
%          [P;ex;ey;hx;hy;m;t] initial-orbit override, default [] = the
%          paper's legacy literal GTO-at-apogee state -- forwarded into both
%          mee_seed calls; see mee_seed.m)
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
%             [4] earth_elliptic_to_geo/nodestudy_mee.m>solve_warm_node
%                 (the mesh-refine-and-single-solve pattern Task 7's
%                 cfg.warmStart direct-eps0 attempt mirrors).
%             [5] .superpowers/sdd/task-7-brief.md (warm-start entry point,
%                 this task's KNOWN GAP to close).
resDir = fullfile(module_root(), 'results');
if ~exist(resDir, 'dir'), mkdir(resDir); end
d = @(f,v) optdef(cfg, f, v);

thrustN    = d('thrustN', 10);
ctf        = d('ctf', 1.5);
tfMinAnchor = d('tfMinAnchor', 22.2248);   % Cartesian M2 min-time anchor [ND]
tag        = d('tag', mee_fuel_tag(thrustN));
seedThr    = d('seedThr', 0.4);
betaMode   = d('betaMode', 'tangential');
nodesPerRev = d('nodesPerRev', 25);
maxIter    = d('maxIter', 1500);
sched      = d('sched', []);
m0kg       = d('m0kg', 1500);
ispS       = d('ispS', 2000);
warmStart  = d('warmStart', []);
xf         = d('xf', [1;0;0;0;0]);   % Task 4: terminal target [P;ex;ey;hx;hy],
                                      % default GEO -- forwarded to mee_seed's
                                      % stopP, homotopy_mee, and the warm-
                                      % direct casadi_lt_mee call
initElems  = d('initElems', []);     % Task 4: initial-orbit override [7x1]
                                      % [P;ex;ey;hx;hy;m;t], default [] (paper
                                      % legacy literal) -- forwarded to mee_seed

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
    'nodesPerRev', nodesPerRev, 'maxIter', maxIter, 'sched', sched, ...
    'xf', xf, 'initElems_isset', ~isempty(initElems));

% DEVIATION (seed throttle): seedThr defaults to ~0.4, not thr=1. A thr=1
% constant burn from the paper's e=0.75 apogee GTO start crosses GEO (P=1) at
% rev~3.06 and goes coordinate-singular past ~rev 4 (Task-2 finding -- the
% MEE seed ODE chases an e->1 singularity once P overshoots the target).
% Throttling to ~0.4 (tangential steering, stopP=1.0) stretches the P=1
% crossing to ~7.5 revs, matching the paper's optimal revolution count and
% giving a well-conditioned collocation seed.

if isempty(warmStart)
% --- stage 1: seed (two-pass: cheap revs probe, then full-density sample) --
probeFile = fullfile(resDir, [tag '_seed_probe.mat']);
if exist(probeFile, 'file')
    S = load(probeFile);  infoP = S.infoP;
    check_cache_fp(S, fpBase, probeFile, tag);
else
    optsP = struct('thr', seedThr, 'betaMode', betaMode, 'N', 50, 'stopP', xf(1), ...
        'initElems', initElems);
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
    opts = struct('thr', seedThr, 'betaMode', betaMode, 'N', N, 'stopP', xf(1), ...
        'initElems', initElems);
    [sigma, X0, U0, dL0, seedInfo] = mee_seed(p, opts);
    save(seedFile, 'sigma', 'X0', 'U0', 'dL0', 'seedInfo', 'fp');
end
x0 = X0(:,1);
fprintf(['RUN_TRANSFER_MEE %s: T=%g N, ctf=%.2f, tf=%.4f ND (%.1f h), ' ...
         'seedThr=%.2f, N=%d nodes (%d nodes/rev), seed revs=%.4f\n'], ...
        tag, thrustN, ctf, tf, tf*p.TU_s/3600, seedThr, N, nodesPerRev, seedInfo.nRev);
else
% --- stage 1 (WARM-START PATH, Task 7): mesh-refine the PREVIOUS ladder
% rung's converged fuel trajectory onto THIS rung's own node grid instead of
% cold-seeding -- skips mee_seed entirely. warmStart.dL is the caller's
% (run_ladder.m) C-law rescale dL_guess = dL_prev*(T_prev/T_new); it both
% sizes the new node count (self-similar-shape assumption: revs_guess =
% dL_guess/2pi) and serves as the solver's dL0 initial guess. interp1
% pattern mirrors nodestudy_mee.m>solve_warm_node: LINEAR for the continuous
% state X and the RTN thrust-direction rows of U (rows 1-3), NEAREST for the
% throttle row (row 4, keeps bang-bang switch edges crisp instead of
% blurring them into intermediate throttle values).
revsGuess = warmStart.dL / (2*pi);
N = round(nodesPerRev * revsGuess);
assert(N >= 1, 'run_transfer_mee:warmStartTooFewNodes', ...
    'warmStart.dL=%.4f rad gives revsGuess=%.4f -- N would be < 1', warmStart.dL, revsGuess);
fp = fpBase;  fp.N = N;  fp.warmStartDL = warmStart.dL;

sigma = linspace(0, 1, N+1).';
W     = interp_warmstart(warmStart.X, warmStart.U, warmStart.dL, warmStart.sigma, sigma);
X0    = W.X;  U0 = W.U;  dL0 = W.dL;
% Rescale the interpolated physical-time row (X row 7) so the INITIAL GUESS
% lands near this rung's own tf target: the previous rung's absolute time
% values live on a different (shorter) tf scale. Every other row (orbital
% elements, mass, controls) is left as interpolated per-sigma-fraction
% values, UNSCALED -- the C-law's self-similar-shape argument (same
% physical spiral shape at fixed sigma-fraction, different total span) is
% exactly what makes this warm start meaningful, and is why m_f comes out
% nearly thrust-independent (paper Fig-23 / this task's gate).
if X0(7,end) > 0
    X0(7,:) = X0(7,:) * (tf / X0(7,end));
end
x0 = X0(:,1);
seedInfo = struct('nRev', revsGuess, 'tEnd', tf, 'mEnd', X0(6,end));
fprintf(['RUN_TRANSFER_MEE %s: T=%g N, ctf=%.2f, tf=%.4f ND (%.1f h), ' ...
         'WARM-STARTED from prior rung (dL_guess=%.4f rad -> revs_guess=%.4f, ' ...
         'N=%d nodes, %d nodes/rev)\n'], ...
        tag, thrustN, ctf, tf, tf*p.TU_s/3600, warmStart.dL, revsGuess, N, nodesPerRev);
end

% --- stage 2: guarded eps 1->0 homotopy at fixed tf -------------------------
% (or, warm-start path: a single direct eps=0 solve, falling back to a short
% continuation tail only if that does not certify -- see cfg.warmStart doc)
if isempty(warmStart)
    ho = struct('par', p, 'x0', x0, 'tfTarget', tf, 'maxIter', maxIter, ...
                'resDir', resDir, 'tag', tag, 'printLevel', 0, 'fp', fp, 'xf', xf);
    if ~isempty(sched), ho.sched = sched; end
    [best, tbl] = homotopy_mee(sigma, X0, U0, dL0, ho);
else
    % Direct eps=0 attempt (same-basin refinement precedent, mirrors
    % nodestudy_mee.m>solve_warm_node): warmTight=true, no sweep. Resume-safe
    % (cached to <tag>_warmdirect.mat, config-fingerprinted like every other
    % cache in this file).
    directFile = fullfile(resDir, [tag '_warmdirect.mat']);
    if exist(directFile, 'file')
        S = load(directFile);
        oD = S.oD;  okD = S.okD;
        check_cache_fp(S, fp, directFile, tag);
    else
        oD = casadi_lt_mee(sigma, X0, U0, dL0, struct('par', p, 'mode', 'fixedtf', ...
            'eps', 0, 'tfTarget', tf, 'x0', x0, 'maxIter', maxIter, ...
            'warmTight', true, 'printLevel', 0, 'xf', xf));
        okD = oD.success && oD.maxDefect < 1e-8;
        save(directFile, 'oD', 'okD', 'fp');
    end
    fprintf('  [warm direct eps=0] ok=%d defect=%.2e sw=%d edge=%.1f%% mf=%.2f kg\n', ...
            okD, oD.maxDefect, oD.switches, 100*oD.edge, oD.m_f_kg);
    if okD
        best = oD;  best.epsReached = 0;  best.certified = true;
        tbl  = [0, oD.maxDefect, oD.switches, oD.edge, oD.m_f_kg];
    else
        fprintf(['  [warm direct eps=0] did NOT certify (defect=%.2e, status=%s) -- ' ...
                 'falling back to a short eps continuation tail\n'], oD.maxDefect, oD.ipoptStatus);
        fallbackSched = d('warmFallbackSched', ...
            [0.3 0.2 0.12 0.07 0.04 0.025 0.015 0.008 0.004 0.002 0.001 0]);
        ho = struct('par', p, 'x0', x0, 'tfTarget', tf, 'maxIter', maxIter, ...
                    'resDir', resDir, 'tag', [tag '_warmtail'], 'printLevel', 0, ...
                    'fp', fp, 'sched', fallbackSched);
        [best, tbl] = homotopy_mee(sigma, X0, U0, dL0, ho);
    end
end

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
function check_cache_fp(S, fp, file, tag)
% CHECK_CACHE_FP  Fail-loud cache-fingerprint guard. If loaded cache struct S
% carries a .fp field, compare it field-by-field against the current config
% fingerprint fp and error, naming the first mismatched field and the
% offending file, on any disagreement -- a stale cache built under a
% different thrustN/ctf/seedThr/... must never be silently reused just
% because it happens to share cfg.tag (the only cache key). BACKWARD COMPAT,
% two distinct cases (Task 4, harmonized with homotopy_mee's check_cache_fp /
% run_mintime_mee.m's check_cache_fp_mt):
%   (1) NO .fp AT ALL (pre-fix cache file, e.g. the MEE_M2_10N certified run
%       predating this guard) -- WARN and trust as-is, no per-field
%       comparison possible.
%   (2) SCHEMA-OLDER .fp (Task 4): the cache HAS a .fp, but a field that
%       exists in the CURRENT fp (e.g. the newly added xf/initElems_isset)
%       is simply ABSENT from the cached one -- this is schema evolution,
%       not a configuration mismatch, and must not hard-error: the
%       documented default-cfg reuse of results/MEE_M2_10N*.mat (saved
%       before xf/initElems existed) is exactly this case. WARN (id
%       'run_transfer_mee:fpSchemaOlder') and treat as compatible. The hard
%       error is preserved for fields present on BOTH sides with different
%       values -- a genuine configuration drift under the same tag must
%       still fail loudly.
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
    if ~isfield(S.fp, f)
        warning('run_transfer_mee:fpSchemaOlder', ['%s: field ''%s'' ' ...
            'present in current fp but absent from cache (schema ' ...
            'evolution) -- trusting as compatible under tag=''%s'''], ...
            file, f, tag);
        continue;
    end
    if ~isequal(S.fp.(f), fp.(f))
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
    [rhat, that, nhat] = rtn_frame(rk, v(:,k));
    Rrtn2eci = [rhat, that, nhat];         % columns = radial, transverse, normal (inertial)
    accGrav   = -mu * rk / rmag^3;
    accThrust = (Tm/m) * thr * (Rrtn2eci * beta);
    residuals(k-1) = norm(vdotFD - (accGrav + accThrust));
end
recon = struct('maxResidual', max(residuals), 'nInterior', nInt, 'residuals', residuals);
end
