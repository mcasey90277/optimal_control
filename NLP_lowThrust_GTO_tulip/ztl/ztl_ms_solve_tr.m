function [z, out] = ztl_ms_solve_tr(z0, prob, opts)
% ZTL_MS_SOLVE_TR  SVD-based Levenberg-Marquardt TRUST REGION for the
% multiple-shooting system (More-style), on the column-scaled Jacobian.
%
% Why this over the earlier variants: the residual valley is stiff AND
% nonlinear, so the full (equilibrated) Newton step overshoots far from the
% solution while a crude diagonal-LM under-moves the stiff directions. A
% proper trust region fixes both: it solves min ||Jc s + R||^2 s.t. ||s|| <=
% Delta EXACTLY (the LM step s(mu) with mu set to the boundary), takes a
% controlled-LENGTH step, and adapts Delta by the actual-vs-predicted
% reduction RATIO on the TRUE objective -- interpolating between a short
% steepest-descent-like step (globalizing) and the full Newton step (fast,
% quadratic near the root). The SVD of the column-scaled Jacobian
% Jc = J*diag(1/colnorm) (objective-preserving; cond ~8e8 vs 3e11 raw) makes
% the boundary mu-search a closed-form scalar solve, so it is free.
%
% Scaled coordinates s = D*dz, D = diag(column norms of J); the trust region
% ||s|| <= Delta is the standard MINPACK metric.
%
% INPUTS:
%   z0, prob - as ztl_ms_residual
%   opts - .tolR [1e-9] .maxIter [200] .Dmax [1e6] .eta [1e-4] .verbose [true]
%
% OUTPUTS:
%   z, out (.resNorm .iters .flag [1 conv|0 maxIter|-2 Delta collapse]
%          .hist [iters x 4: ||R|| Delta rho condJc] .termErr .maxCont .grazed)
%
% REFERENCES: More, "The Levenberg-Marquardt algorithm: implementation and
%   theory", 1978; Nocedal & Wright ch.4 (trust region), ch.10 (LSQ).

if nargin < 3, opts = struct(); end
g = @(f,d) getdef(opts, f, d);
tolR = g('tolR', 1e-9);  maxIter = g('maxIter', 200);
Dmax = g('Dmax', 1e6);   eta = g('eta', 1e-4);  verbose = g('verbose', true);

z = z0(:);
[R, J, info] = ztl_ms_residual(z, prob, true);
rn = norm(R);
[D, Jc, V, sig, UtR, sgn] = refactor(J, R);
% Start Delta modest: the full GN step (||sgn|| ~ 1e2) overshoots the stiff
% nonlinear valley and wastes iterations shrinking; a modest Delta lets the
% ratio test grow it into the quadratic phase near the solution. A warm
% restart passes the last accepted Delta via opts.Delta0.
Delta = g('Delta0', min(norm(sgn), 1));
hist = nan(maxIter, 4);  flag = 0;  it = 0;

while it < maxIter
    it = it + 1;
    condJc = sig(1)/max(sig(end), realmin);
    hist(it, :) = [rn, Delta, NaN, condJc];
    if verbose
        fprintf('  it %3d: ||R||=%.4e  Delta=%.2e  cond(Jc)=%.2e  termErr=%.2e  maxCont=%.2e\n', ...
                it, rn, Delta, condJc, info.termErr, info.maxCont);
    end
    if rn < tolR, flag = 1; break; end

    % --- trust-region subproblem (reuses the current SVD) -------------------
    if isfinite(norm(sgn)) && norm(sgn) <= Delta
        s = sgn;  atBoundary = false;
    else
        mu = lm_mu(sig, UtR, Delta);
        s  = -V * ((sig .* UtR) ./ (sig.^2 + mu));
        atBoundary = true;
    end
    dz = s ./ D;                                % dz = diag(1/D) s
    pred = 0.5*rn^2 - 0.5*norm(R + Jc*s)^2;      % model reduction (>=0)

    [Rt, Jt, it2] = ztl_ms_residual(z + dz, prob, true);
    rnt = norm(Rt);
    rho = (0.5*rn^2 - 0.5*rnt^2) / max(pred, realmin);
    hist(it, 3) = rho;

    if rho < 0.25
        Delta = 0.25*Delta;
    elseif rho > 0.75 && atBoundary
        Delta = min(2*Delta, Dmax);
    end

    if rho > eta && isfinite(rnt)               % accept
        z = z + dz;  R = Rt;  J = Jt;  info = it2;  rn = rnt;
        [D, Jc, V, sig, UtR, sgn] = refactor(J, R);
    end
    if Delta < 1e-13*max(norm(z),1), flag = -2; break; end
end

out = struct('resNorm', rn, 'iters', it, 'flag', flag, 'hist', hist(1:it,:), ...
             'termErr', info.termErr, 'maxCont', info.maxCont, 'grazed', info.grazed, ...
             'Delta', Delta);
end

% ---------------------------------------------------------------------------
function [D, Jc, V, sig, UtR, sgn] = refactor(J, R)
% Column scaling + SVD of the column-equilibrated Jacobian; scaled GN step.
D = sqrt(sum(J.^2, 1)).';  D(D == 0) = 1;
Jc = J ./ D.';                                  % J * diag(1/D)
[U, S, V] = svd(full(Jc));
sig = diag(S);
UtR = U.' * R;
tol = sig(1) * 1e-14;
if sig(end) > tol
    sgn = -V * (UtR ./ sig);                    % scaled Gauss-Newton step
else
    sgn = Inf(size(UtR));                        % rank-deficient -> force boundary
end
end

function mu = lm_mu(sig, UtR, Delta)
% Smallest mu >= 0 with ||s(mu)|| = Delta, s(mu) = -V (sig.*UtR)./(sig.^2+mu).
% ||s(mu)||^2 = sum( (sig.*UtR ./ (sig.^2+mu)).^2 ), monotone decreasing in mu.
a = sig .* UtR;
phi = @(mu) sqrt(sum((a ./ (sig.^2 + mu)).^2)) - Delta;
lo = 0;
hi = 1;
while phi(hi) > 0
    hi = hi*10;
    if hi > 1e25, mu = hi; return; end
end
% geometric bisection on [max(lo,tiny), hi]
lo = max(lo, hi*1e-25);
for k = 1:80
    mid = sqrt(lo*hi);
    if phi(mid) > 0, lo = mid; else, hi = mid; end
end
mu = sqrt(lo*hi);
end

function v = getdef(s, f, d)
if isfield(s, f) && ~isempty(s.(f)), v = s.(f); else, v = d; end
end
