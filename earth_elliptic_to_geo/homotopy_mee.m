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
%          .tfTarget .maxIter .sched .resDir .tag .printLevel
% OUTPUTS: best - last tight solver out + .certified .epsReached;  tbl [Kx5]
%          = [eps, maxDefect, switches, edge, m_f_kg]
%
% REFERENCES: [1] earth_elliptic_to_geo/homotopy_2body.m (pattern this mirrors).
%             [2] earth_elliptic_to_geo/casadi_lt_mee.m (per-step solver).
d = @(f,v) getdef5(opts, f, v);
sched   = d('sched', [1 0.6 0.35 0.2 0.12 0.07 0.04 0.025 0.015 0.008 0.004 0.002 0.001 0]);
maxIter = d('maxIter', 1500);
resDir  = d('resDir', pwd);
tag     = d('tag', 'mee_run');
printLvl = d('printLevel', 0);

if ~exist(resDir, 'dir'), mkdir(resDir); end

Xk = X0;  Uk = U0;  dLk = dL0;  best = [];  tbl = zeros(numel(sched), 5);
for ke = 1:numel(sched)
    e = sched(ke);
    stepFile = fullfile(resDir, sprintf('%s_step%02d.mat', tag, ke));
    if exist(stepFile, 'file')
        S  = load(stepFile);
        o  = S.o;  ok = S.ok;
        Xk = S.Xk;  Uk = S.Uk;  dLk = S.dLk;
        if ok, best = o;  best.epsReached = e; end
        tbl(ke,:) = [e, o.maxDefect, o.switches, o.edge, o.m_f_kg];
        fprintf('  [cached] eps=%6.4f ok=%d defect=%.2e sw=%3d edge=%5.1f%% mf=%.2f kg\n', ...
                e, ok, o.maxDefect, o.switches, 100*o.edge, o.m_f_kg);
        continue;
    end
    o = casadi_lt_mee(sigma, Xk, Uk, dLk, struct('par', opts.par, ...
        'mode', 'fixedtf', 'eps', e, 'tfTarget', opts.tfTarget, 'x0', opts.x0, ...
        'maxIter', maxIter, 'warmTight', ke > 1, 'printLevel', printLvl));
    ok = o.success && o.maxDefect < 1e-8;
    tbl(ke,:) = [e, o.maxDefect, o.switches, o.edge, o.m_f_kg];
    fprintf('  eps=%6.4f ok=%d defect=%.2e sw=%3d edge=%5.1f%% mf=%.2f kg\n', ...
            e, ok, o.maxDefect, o.switches, 100*o.edge, o.m_f_kg);
    if ok
        Xk = o.X;  Uk = o.U;  dLk = o.dL;  best = o;  best.epsReached = e;
    end
    save(stepFile, 'o', 'ok', 'Xk', 'Uk', 'dLk', 'e');
end
if isempty(best)
    best = o;  best.epsReached = NaN;  best.certified = false;
else
    best.certified = (best.epsReached == 0);
end
end

% ---------------------------------------------------------------------------
function v = getdef5(s, f, dflt)
% GETDEF5  Optional-field default (mirrors casadi_lt_2body's helper).
if isfield(s, f) && ~isempty(s.(f)), v = s.(f); else, v = dflt; end
end
