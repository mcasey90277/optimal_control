function [z, out] = ztl_ms_solve_newton(z0, prob, opts)
% ZTL_MS_SOLVE_NEWTON  Newton on the square multiple-shooting system, with the
% Newton step computed via TWO-SIDED-EQUILIBRATED solve and globalized by a
% line search on the TRUE residual.
%
% Key idea (why this beats the column-scaled LM, which floored at ~3e-5 on the
% continuity rows): for a SQUARE nonsingular J the exact Newton step is
% p = -J\R, and two-sided equilibration is just a numerically STABLE way to
% compute it -- the row/col scalings cancel in exact arithmetic
% (Dc*((Dr J Dc)\(-Dr R)) = -J\R) but the solved system has cond ~5e6 instead
% of ~3e11, recovering ~10 accurate digits. The direction is therefore the
% true Newton direction (guaranteed descent for ||R||^2), so an Armijo
% backtracking line search on the TRUE ||R|| globalizes it -- no
% objective-changing row-weighting, no ascent steps. Near the root this is
% quadratic; far from it the line search shortens the (overshooting) step.
%
% INPUTS:
%   z0   - initial unknown [14M-7 x 1]
%   prob - problem struct for ztl_ms_residual
%   opts - .tolR [1e-9] .maxIter [120] .lamReg [1e-12] (Levenberg floor if J
%          is near-singular) .verbose [true]
%
% OUTPUTS:
%   z, out (.resNorm .iters .flag [1 conv|0 maxIter|-2 stall] .hist
%          [iters x 3: ||R|| alpha condEq] .termErr .maxCont .grazed)
%
% REFERENCES: Deuflhard, Newton Methods for Nonlinear Problems, 2004
%   (affine-invariant Newton + damping); Ruiz equilibration.

if nargin < 3, opts = struct(); end
g = @(f,d) getdef(opts, f, d);
tolR = g('tolR', 1e-9);  maxIter = g('maxIter', 120);
lamReg = g('lamReg', 1e-12);  verbose = g('verbose', true);

z = z0(:);
[R, J, info] = ztl_ms_residual(z, prob, true);
rn = norm(R);
hist = nan(maxIter, 3);  flag = 0;  stall = 0;

for it = 1:maxIter
    [Dr, Dc] = ruiz(J, 5);
    Je = (Dr * J) * Dc;
    condEq = condest(Je'*Je)^0.5;                 % cheap cond estimate
    hist(it, :) = [rn, NaN, condEq];
    if verbose
        fprintf('  it %3d: ||R||=%.4e  cond(eq)=%.2e  termErr=%.2e  maxCont=%.2e\n', ...
                it, rn, condEq, info.termErr, info.maxCont);
    end
    if rn < tolR, flag = 1; break; end

    % equilibrated Newton direction (= -J\R, stably): add a tiny Levenberg
    % floor in equilibrated space for near-singular J.
    nz = numel(z);
    w  = [Je; sqrt(lamReg)*speye(nz)] \ [-(Dr*R); zeros(nz,1)];
    p  = Dc * w;

    % Armijo backtracking on the TRUE residual norm
    alpha = 1;  accepted = false;
    for ls = 1:30
        [Rt, Jt, it2] = ztl_ms_residual(z + alpha*p, prob, true);
        if all(isfinite(Rt)) && norm(Rt) < (1 - 1e-4*alpha)*rn
            z = z + alpha*p;  R = Rt;  J = Jt;  info = it2;  rn = norm(R);
            accepted = true;  hist(it, 2) = alpha;  break
        end
        alpha = alpha/2;
    end
    if ~accepted
        stall = stall + 1;
        if stall >= 3, flag = -2; break; end
    else
        stall = 0;
    end
end

out = struct('resNorm', rn, 'iters', it, 'flag', flag, 'hist', hist(1:it,:), ...
             'termErr', info.termErr, 'maxCont', info.maxCont, 'grazed', info.grazed);
end

% ---------------------------------------------------------------------------
function [Dr, Dc] = ruiz(J, iters)
[m, n] = size(J);  dr = ones(m,1);  dc = ones(n,1);
for it = 1:iters
    A = (spdiags(dr,0,m,m)*J)*spdiags(dc,0,n,n);
    rr = sqrt(max(abs(A),[],2));  rr(rr==0)=1;
    cc = sqrt(max(abs(A),[],1)).'; cc(cc==0)=1;
    dr = dr./rr;  dc = dc./cc;
end
Dr = spdiags(dr,0,m,m);  Dc = spdiags(dc,0,n,n);
end

function v = getdef(s, f, d)
if isfield(s, f) && ~isempty(s.(f)), v = s.(f); else, v = d; end
end
