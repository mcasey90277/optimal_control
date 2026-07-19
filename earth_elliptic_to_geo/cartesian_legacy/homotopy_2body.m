function [best, tbl] = homotopy_2body(sigma, X0, U0, tauf0, term, tf, opts)
% HOMOTOPY_2BODY  Guarded energy->fuel sweep at fixed t_f (eps: 1 -> 0).
%
% First step (eps=1) runs LOOSE (genuine move from a propagated seed); every
% later step warm-starts TIGHT from the previous converged iterate. GUARD: a
% step that fails to converge tight never advances the warm start and never
% overwrites best (campaign lesson: a loose iterate must not poison the chain).
%
% INPUTS:  sigma/X0/U0/tauf0 - seed (seed_2body layout);  term - geo_terminal
%          struct;  tf - fixed transfer time [ND];  opts - .par .rv0 .maxIter .sched
% OUTPUTS: best - last tight solver out + .certified .epsReached;  tbl [Kx5]
%          = [eps, maxDefect, switches, edge, m_f_kg]
%
% REFERENCES: [1] sundman_minfuel/sundman_homotopy.m (pattern). [2] DESIGN.md sec 4.
d = @(f,v) optdef(opts, f, v);
sched   = d('sched', [1 0.6 0.35 0.2 0.12 0.07 0.04 0.025 0.015 0.008 0.004 0.002 0.001 0]);
maxIter = d('maxIter', 1500);
Xk = X0;  Uk = U0;  best = [];  tbl = zeros(numel(sched), 5);
for ke = 1:numel(sched)
    e = sched(ke);
    o = casadi_lt_2body(sigma, Xk, Uk, tauf0, term, struct('par',opts.par, ...
        'mode','fixedtf', 'eps',e, 'tfTarget',tf, 'rv0',opts.rv0, ...
        'maxIter',maxIter, 'warmTight', ke > 1, 'printLevel',0));
    ok = o.success && o.maxDefect < 1e-8;
    tbl(ke,:) = [e, o.maxDefect, o.switches, o.edge, o.m_f_kg];
    fprintf('  eps=%6.4f ok=%d defect=%.2e sw=%3d edge=%5.1f%% mf=%.2f kg\n', ...
            e, ok, o.maxDefect, o.switches, 100*o.edge, o.m_f_kg);
    if ok
        Xk = o.X;  Uk = o.U;  best = o;  best.epsReached = e;
    end
end
if isempty(best)
    best = o;  best.epsReached = NaN;  best.certified = false;
else
    best.certified = (best.epsReached == 0);
end
end
