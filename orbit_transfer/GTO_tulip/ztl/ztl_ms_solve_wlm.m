function [z, out] = ztl_ms_solve_wlm(z0, prob, opts)
% ZTL_MS_SOLVE_WLM  Weighted Levenberg-Marquardt for the multiple-shooting
% system: FIXED affine-invariant row weighting + per-iter column scaling.
%
% Combines the three things the earlier variants each had only one of:
%   - two-sided conditioning (cond ~5e6, vs 3e11 raw) for accurate steps;
%   - a CONSISTENT objective so every step is genuine descent (no ascent bug);
%   - LM damping for globalization (the pure Newton step overshoots the seed).
% It minimizes the FIXED-weighted least squares ||Dr R||^2 with Dr =
% 1/rownorm(J at the seed) held CONSTANT (Deuflhard affine-invariant scaling
% -- makes every equation unit-sensitivity). Because Dr is fixed and positive,
% ||Dr R|| -> 0 iff ||R|| -> 0, so the weighted objective is legitimate and its
% LM step descends it. Within each iter the weighted Jacobian Dr*J is
% column-scaled (Dc, per-iter) and damped:
%       [ Dr J Dc ; sqrt(mu) I ] w = [ -Dr R ; 0 ],   dz = Dc w.
% Acceptance is on the weighted norm; the campaign GATE is the TRUE ||R||.
%
% INPUTS:
%   z0, prob  - as ztl_ms_residual
%   opts - .tolR [1e-9] (true ||R|| gate) .maxIter [200] .mu0 [1e-6]
%          .verbose [true]
%
% OUTPUTS:
%   z, out (.resNorm [true] .resNormW [weighted] .iters .flag
%          [1 conv|0 maxIter|-2 stall] .hist .termErr .maxCont .grazed)
%
% REFERENCES: Deuflhard, Newton Methods for Nonlinear Problems, 2004 (scaling);
%   Marquardt 1963.

if nargin < 3, opts = struct(); end
g = @(f,d) getdef(opts, f, d);
tolR = g('tolR', 1e-9);  maxIter = g('maxIter', 200);
mu = g('mu0', 1e-6);  verbose = g('verbose', true);

z = z0(:);
[R, J, info] = ztl_ms_residual(z, prob, true);

% FIXED affine-invariant row weighting from the seed Jacobian
rw = sqrt(sum(J.^2, 2));  rw(rw == 0) = 1;
dr = 1 ./ rw;
m = numel(dr);
Dr = spdiags(dr, 0, m, m);

rn = norm(R);  rnW = norm(Dr*R);
hist = nan(maxIter, 4);  flag = 0;  stall = 0;

for it = 1:maxIter
    hist(it, :) = [rn, rnW, mu, info.maxCont];
    if verbose
        fprintf('  it %3d: ||R||=%.4e  ||DrR||=%.4e  mu=%.1e  termErr=%.2e  maxCont=%.2e\n', ...
                it, rn, rnW, mu, info.termErr, info.maxCont);
    end
    if rn < tolR, flag = 1; break; end

    Jw = Dr * J;  Rw = Dr * R;
    cn = sqrt(sum(Jw.^2, 1)).';  cn(cn == 0) = 1;
    Dc = spdiags(1./cn, 0, numel(cn), numel(cn));
    Jwc = Jw * Dc;
    nz = numel(z);

    accepted = false;
    for tries = 1:14
        w  = [Jwc; sqrt(mu)*speye(nz)] \ [-Rw; zeros(nz,1)];
        dz = Dc * w;
        [Rt, Jt, it2] = ztl_ms_residual(z + dz, prob, true);
        if all(isfinite(Rt)) && norm(Dr*Rt) < rnW
            z = z + dz;  R = Rt;  J = Jt;  info = it2;
            rn = norm(R);  rnW = norm(Dr*R);  mu = max(mu*0.3, 1e-12);
            accepted = true;  break
        else
            mu = mu*4;
        end
    end
    if ~accepted
        stall = stall + 1;
        if mu > 1e12 || stall >= 4, flag = -2; break; end
    else
        stall = 0;
    end
end

out = struct('resNorm', rn, 'resNormW', rnW, 'iters', it, 'flag', flag, ...
             'hist', hist(1:it,:), 'termErr', info.termErr, ...
             'maxCont', info.maxCont, 'grazed', info.grazed);
end

% ---------------------------------------------------------------------------
function v = getdef(s, f, d)
if isfield(s, f) && ~isempty(s.(f)), v = s.(f); else, v = d; end
end
