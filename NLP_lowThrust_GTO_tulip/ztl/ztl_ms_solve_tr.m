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
geodesic = g('geodesic', true);   % 2nd-order (geodesic) acceleration
geoH = g('geoH', 0.1);  geoAlpha = g('geoAlpha', 0.75);

z = z0(:);
rf = @ztl_ms_residual;  if isfield(prob,'resFun'), rf = prob.resFun; end
[R, J, info] = rf(z, prob, true);
rn = norm(R);
[D, Jc, U, V, sig, UtR, sgn] = refactor(J, R);
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

    % --- trust-region subproblem -------------------------------------------
    % INTERIOR (near the solution): use the accurate TWO-SIDED-equilibrated
    % Newton step -- for a square system it equals -J\R but is computed at
    % cond ~5e6 (vs the column-only ~1e9), lowering the linear-solve NOISE
    % FLOOR from ~5e-7 to ~1e-9 (the last barrier at ||R||~1e-6). BOUNDARY
    % (far): the column-SVD LM mu-step globalizes.
    interior = isfinite(norm(sgn)) && norm(sgn) <= Delta;
    if interior
        atBoundary = false;
        dz = solve2s(J, -R);                     % accurate Newton (cond ~5e6)
        if geodesic
            Jdz = J*dz;
            Rh  = rf(z + geoH*dz, prob, false);
            fvv = (2/geoH^2)*(Rh - R - geoH*Jdz);
            da  = solve2s(J, -fvv);
            if norm(da) <= geoAlpha*norm(dz), dz = dz + 0.5*da; end
        end
        predResid = R + J*dz;
    else
        atBoundary = true;
        muUsed = lm_mu(sig, UtR, Delta);
        s  = -V * ((sig .* UtR) ./ (sig.^2 + muUsed));   % column-SVD damped step
        dz = s ./ D;
        if geodesic
            Jdz = Jc * s;
            Rh  = rf(z + geoH*dz, prob, false);
            fvv = (2/geoH^2)*(Rh - R - geoH*Jdz);
            sa  = -V * ((sig .* (U.'*fvv)) ./ (sig.^2 + muUsed));
            if norm(sa) <= geoAlpha*norm(s), s = s + 0.5*sa; dz = s ./ D; end
        end
        predResid = R + Jc*s;
    end
    pred = 0.5*rn^2 - 0.5*norm(predResid)^2;      % model reduction

    [Rt, Jt, it2] = rf(z + dz, prob, true);
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
        [D, Jc, U, V, sig, UtR, sgn] = refactor(J, R);
    end
    if Delta < 1e-13*max(norm(z),1), flag = -2; break; end
end

out = struct('resNorm', rn, 'iters', it, 'flag', flag, 'hist', hist(1:it,:), ...
             'termErr', info.termErr, 'maxCont', info.maxCont, 'grazed', info.grazed, ...
             'Delta', Delta);
end

% ---------------------------------------------------------------------------
function [D, Jc, U, V, sig, UtR, sgn] = refactor(J, R)
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

function s = solve_damped2s(Jc, R, mu)
% Damped LM step s = argmin ||Jc s + R||^2 + mu||s||^2, computed by two-sided
% (Ruiz) equilibration of the augmented system [Jc; sqrt(mu) I] s = [-R; 0].
% The equilibration lowers the augmented condition number so the step is
% accurate in the small-singular-value (stiff) directions where the
% column-SVD reconstruction floors at the conditioning noise.
n = size(Jc, 2);
A = [Jc; sqrt(mu)*speye(n)];
b = [-R; zeros(n, 1)];
[m2, n2] = size(A);  dr = ones(m2,1);  dc = ones(n2,1);
for it = 1:5
    M = (spdiags(dr,0,m2,m2)*A)*spdiags(dc,0,n2,n2);
    rr = sqrt(max(abs(M),[],2));  rr(rr==0)=1;
    cc = sqrt(max(abs(M),[],1)).'; cc(cc==0)=1;
    dr = dr./rr;  dc = dc./cc;
end
Ae = (spdiags(dr,0,m2,m2)*A)*spdiags(dc,0,n2,n2);
s = dc .* (Ae \ (dr .* b));
end

function x = solve2s(J, b)
% Solve J x = b via TWO-SIDED (Ruiz) equilibration: x = Dc*((Dr J Dc)\(Dr b)).
% For square nonsingular J this equals J\b but the solved system has cond ~5e6
% instead of ~1e9 -- lowering the linear-solve noise floor by ~2-3 orders.
[m, n] = size(J);  dr = ones(m,1);  dc = ones(n,1);
for it = 1:5
    A = (spdiags(dr,0,m,m)*J)*spdiags(dc,0,n,n);
    rr = sqrt(max(abs(A),[],2));  rr(rr==0)=1;
    cc = sqrt(max(abs(A),[],1)).'; cc(cc==0)=1;
    dr = dr./rr;  dc = dc./cc;
end
Je = (spdiags(dr,0,m,m)*J)*spdiags(dc,0,n,n);
x = dc .* (Je \ (dr .* b));
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
