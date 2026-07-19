function [best, tbl] = homotopy_mee(sigma, X0, U0, dL0, opts)
% HOMOTOPY_MEE  Guarded energy->fuel sweep at fixed t_f (eps: 1 -> 0), MEE/
% L-domain analog of homotopy_2body.m.
%
% First step (eps=1) runs LOOSE (genuine move from a propagated seed); every
% later step warm-starts TIGHT from the previous converged iterate. GUARD: a
% step that fails to converge tight never advances the warm start and never
% overwrites best (same campaign lesson as homotopy_2body.m: a loose iterate
% must not poison the chain).
%
% RESUME-SAFE: each eps-step's solver output + the post-step warm-start state
% (Xk,Uk,dLk) is cached to opts.resDir/<opts.tag>_step<k>.mat; a step whose
% cache file already exists is loaded instead of re-solved. This lets a
% 10-minute Bash-tool call (or a sporadic MEX init crash) be resumed by
% simply re-invoking the caller -- already-completed eps-steps are skipped.
%
% INPUTS:  sigma/X0/U0/dL0 - seed (mee_seed layout); opts - .par .x0
%          .tfTarget .maxIter .sched .resDir .tag .printLevel .xf ([5x1]
%          terminal target [P;ex;ey;hx;hy], default [1;0;0;0;0] = GEO;
%          forwarded to every eps-step casadi_lt_mee call and recorded in
%          fp.xf -- Task 3)
% OUTPUTS: best - last tight solver out + .certified .epsReached;  tbl [Kx5]
%          = [eps, maxDefect, switches, edge, m_f_kg]
%
% REFERENCES: [1] earth_elliptic_to_geo/homotopy_2body.m (pattern this mirrors).
%             [2] earth_elliptic_to_geo/casadi_lt_mee.m (per-step solver).
d = @(f,v) optdef(opts, f, v);
sched   = d('sched', [1 0.6 0.35 0.2 0.12 0.07 0.04 0.025 0.015 0.008 0.004 0.002 0.001 0]);
maxIter = d('maxIter', 1500);
resDir  = d('resDir', pwd);
tag     = d('tag', 'mee_run');
printLvl = d('printLevel', 0);
xf       = d('xf', [1;0;0;0;0]);
fp       = d('fp', struct());
fp.sched = sched;   % overwrite with the ACTUALLY-resolved schedule (opts.fp
                     % may carry the caller's raw cfg.sched, which is [] when
                     % the caller means "use this function's own default")
fp.xf    = xf;      % Task 3: terminal target now user-definable (default
                     % GEO); recorded so a differing xf under the same tag
                     % is caught by check_cache_fp below

if ~exist(resDir, 'dir'), mkdir(resDir); end

Xk = X0;  Uk = U0;  dLk = dL0;  best = [];  tbl = zeros(numel(sched), 5);
for ke = 1:numel(sched)
    e = sched(ke);
    stepFile = fullfile(resDir, sprintf('%s_step%02d.mat', tag, ke));
    if exist(stepFile, 'file')
        S  = load(stepFile);
        o  = S.o;  ok = S.ok;
        Xk = S.Xk;  Uk = S.Uk;  dLk = S.dLk;
        check_cache_fp(S, fp, stepFile, tag);
        if ok, best = o;  best.epsReached = e; end
        tbl(ke,:) = [e, o.maxDefect, o.switches, o.edge, o.m_f_kg];
        fprintf('  [cached] eps=%6.4f ok=%d defect=%.2e sw=%3d edge=%5.1f%% mf=%.2f kg\n', ...
                e, ok, o.maxDefect, o.switches, 100*o.edge, o.m_f_kg);
        continue;
    end
    o = casadi_lt_mee(sigma, Xk, Uk, dLk, struct('par', opts.par, ...
        'mode', 'fixedtf', 'eps', e, 'tfTarget', opts.tfTarget, 'x0', opts.x0, ...
        'maxIter', maxIter, 'warmTight', ke > 1, 'printLevel', printLvl, 'xf', xf, ...
        'liftDL', optdef(opts, 'liftDL', false)));
    ok = o.success && o.maxDefect < 1e-8;
    tbl(ke,:) = [e, o.maxDefect, o.switches, o.edge, o.m_f_kg];
    fprintf('  eps=%6.4f ok=%d defect=%.2e sw=%3d edge=%5.1f%% mf=%.2f kg\n', ...
            e, ok, o.maxDefect, o.switches, 100*o.edge, o.m_f_kg);
    if ok
        Xk = o.X;  Uk = o.U;  dLk = o.dL;  best = o;  best.epsReached = e;
    end
    save(stepFile, 'o', 'ok', 'Xk', 'Uk', 'dLk', 'e', 'fp');
end
if isempty(best)
    best = o;  best.epsReached = NaN;  best.certified = false;
else
    best.certified = (best.epsReached == 0);
end
end

% ---------------------------------------------------------------------------
function check_cache_fp(S, fp, file, tag)
% CHECK_CACHE_FP  Fail-loud cache-fingerprint guard (mirrors
% run_transfer_mee.m's helper of the same name, and run_mintime_mee.m's
% check_cache_fp_mt). If loaded per-eps-step cache struct S carries a .fp
% field, compare it field-by-field against the current config fingerprint fp
% and error, naming the first mismatched field and the offending file, on
% any disagreement. BACKWARD COMPAT, two distinct cases:
%   (1) NO .fp AT ALL (pre-fingerprint-guard cache) -- WARN and trust as-is,
%       no per-field comparison possible.
%   (2) SCHEMA-OLDER .fp (Task 3): the cache HAS a .fp, but a field that
%       exists in the CURRENT fp (e.g. the newly added xf) is simply ABSENT
%       from the cached one -- this is schema evolution, not a
%       configuration mismatch, and must not hard-error. WARN (id
%       'homotopy_mee:fpSchemaOlder') and treat as compatible. The hard
%       error is preserved for fields present on BOTH sides with different
%       values -- a genuine configuration drift under the same tag must
%       still fail loudly.
if ~isfield(S, 'fp')
    warning('homotopy_mee:noCachedFingerprint', ['%s has no cached config ' ...
        'fingerprint (pre-fix cache) -- trusting it because tag=''%s'' ' ...
        'matches; use a new tag to regain fingerprint protection for this ' ...
        'run'], file, tag);
    return;
end
flds = fieldnames(fp);
for k = 1:numel(flds)
    f = flds{k};
    if ~isfield(S.fp, f)
        warning('homotopy_mee:fpSchemaOlder', ['%s: field ''%s'' present in ' ...
            'current fp but absent from cache (schema evolution) -- trusting ' ...
            'as compatible under tag=''%s'''], file, f, tag);
        continue;
    end
    if ~isequal(S.fp.(f), fp.(f))
        error('homotopy_mee:fingerprintMismatch', ['cached config ' ...
            'fingerprint mismatch in %s: field ''%s'' differs between the ' ...
            'cache and the current config -- stale cache from a different ' ...
            'configuration under the same tag=''%s''; delete the file or ' ...
            'use a new tag'], file, f, tag);
    end
end
end
