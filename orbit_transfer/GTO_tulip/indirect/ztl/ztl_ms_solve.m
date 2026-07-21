function [z, out] = ztl_ms_solve(z0, prob, opts)
% ZTL_MS_SOLVE  Levenberg-Marquardt on the square multiple-shooting system,
% via the AUGMENTED-QR form (never the normal equations).
%
% Minimizes ||R(z)||^2 with R from ztl_ms_residual and the exact block
% Jacobian. The MS Jacobian is better conditioned than the single-shooting
% 7x7 but is still ill-conditioned (cond ~1e11), so the normal equations
% J'J (cond ~1e22 -> numerically singular) MUST be avoided. Instead the LM
% step solves the STACKED least-squares system
%       [ J ; sqrt(mu) D ] dz = [ -R ; 0 ],   D = diag(column norms),
% by QR (MATLAB backslash on the 2n-by-n system). Its condition number is
% ~cond(J), not cond(J)^2 -- the difference between ~5 usable digits and
% none. mu is decreased on acceptance, increased on rejection.
%
% INPUTS:
%   z0   - initial unknown vector [14M-7 x 1]
%   prob - problem struct for ztl_ms_residual
%   opts - (optional) struct:
%          .tolR     convergence on ||R|| [1e-9]
%          .maxIter  [200]
%          .mu0      initial LM damping [1e-6]
%          .verbose  print per-iter line [true]
%
% OUTPUTS:
%   z   - solution [14M-7 x 1]
%   out - struct: .resNorm .iters .flag (1 converged | 0 maxIter |
%         -1 mu blowup | -2 stall) .hist [iters x 4: ||R|| mu condJ termErr]
%
% REFERENCES: Levenberg 1944; Marquardt 1963; Nocedal & Wright ch.10.

if nargin < 3, opts = struct(); end
g = @(f,d) getdef(opts, f, d);
tolR = g('tolR', 1e-9);  maxIter = g('maxIter', 200);
mu = g('mu0', 1e-6);  verbose = g('verbose', true);

z = z0(:);
[R, J, info] = ztl_ms_residual(z, prob, true);
rn = norm(R);
hist = nan(maxIter, 4);
flag = 0;  stall = 0;

for it = 1:maxIter
    condJ = cond(J);
    hist(it, :) = [rn, mu, condJ, info.termErr];
    if verbose
        fprintf('  it %3d: ||R||=%.4e  mu=%.1e  cond(J)=%.2e  termErr=%.2e  maxCont=%.2e\n', ...
                it, rn, mu, condJ, info.termErr, info.maxCont);
    end
    if rn < tolR, flag = 1; break; end

    % Column-scaled (Marquardt) LM step -- OBJECTIVE-PRESERVING. The MS
    % Jacobian is artificially ill-conditioned by the disparate scales of the
    % state vs costate columns (cond 3e11); column equilibration dz = Dc w
    % drops it to ~1e9 (col-only) WITHOUT changing the objective, so the step
    % is guaranteed descent. (Row-scaling the residual would lower cond
    % further but changes the objective and gives ascent steps far from the
    % solution -- do NOT do it.)
    cn = sqrt(sum(J.^2, 1)).';  cn(cn == 0) = 1;
    Dc = spdiags(1./cn, 0, numel(cn), numel(cn));
    Jc = J * Dc;                                  % column-equilibrated
    nz = numel(z);
    accepted = false;
    for tries = 1:12
        w  = [Jc; sqrt(mu)*speye(nz)] \ [-R; zeros(nz, 1)];
        dz = Dc * w;
        zt = z + dz;
        [Rt, Jt, it2] = ztl_ms_residual(zt, prob, true);
        if all(isfinite(Rt)) && norm(Rt) < rn
            z = zt;  R = Rt;  J = Jt;  info = it2;
            rn = norm(R);  mu = max(mu*0.3, 1e-12);
            accepted = true;  break
        else
            mu = mu*4;
        end
    end
    if ~accepted
        stall = stall + 1;
        if mu > 1e10, flag = -1; break; end
        if stall >= 4,  flag = -2; break; end
    else
        stall = 0;
    end
end

out = struct('resNorm', rn, 'iters', it, 'flag', flag, 'hist', hist(1:it, :), ...
             'termErr', info.termErr, 'maxCont', info.maxCont, 'grazed', info.grazed);
end

% ---------------------------------------------------------------------------
function [Dr, Dc] = ruiz(J, iters)
% Ruiz two-sided equilibration: diagonal Dr, Dc s.t. rows and columns of
% Dr*J*Dc have ~unit infinity-norm. Reduces the effective condition number.
[m, n] = size(J);
dr = ones(m, 1);  dc = ones(n, 1);
for it = 1:iters
    A = (spdiags(dr,0,m,m) * J) * spdiags(dc,0,n,n);
    rr = sqrt(max(abs(A), [], 2));  rr(rr == 0) = 1;
    cc = sqrt(max(abs(A), [], 1)).';  cc(cc == 0) = 1;
    dr = dr ./ rr;  dc = dc ./ cc;
end
Dr = spdiags(dr, 0, m, m);  Dc = spdiags(dc, 0, n, n);
end

function v = getdef(s, f, d)
if isfield(s, f) && ~isempty(s.(f)), v = s.(f); else, v = d; end
end
