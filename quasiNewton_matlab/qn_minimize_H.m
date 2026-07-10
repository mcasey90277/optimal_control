function [x, info] = qn_minimize_H(fun, x0, update_fn, opts)
% QN_MINIMIZE_H  Quasi-Newton minimization using an inverse-Hessian model.
%
% Minimizes a smooth f: R^n -> R by the quasi-Newton iteration
%   p_k = -H_k * grad f(x_k)            (search direction; matvec, no solve)
%   x_{k+1} = x_k + alpha_k p_k         (alpha_k by backtracking-Armijo)
%   H_{k+1} = update_fn(H_k, s_k, y_k)  (BFGS or DFP inverse-Hessian update)
% The inverse-Hessian form makes the step a matrix-vector product. The update
% is skipped when the curvature condition s'*y > 0 fails (keeps H_k PD).
%
% INPUTS:
%   fun       - handle returning [f, g]: f scalar, g gradient [n x 1]
%   x0        - initial point                                  [n x 1]
%   update_fn - inverse-Hessian update handle, e.g. @bfgs_update_H [fn]
%   opts      - (optional) struct with fields:
%                 H0      initial inverse-Hessian (default I)   [n x n]
%                 tol     stop when ||g|| <= tol (default 1e-8) [scalar]
%                 maxit   iteration cap (default 1000)          [scalar]
%                 c1      Armijo constant (default 1e-4)        [scalar]
%                 bt      backtracking factor (default 0.5)     [scalar]
%                 linesearch  handle alpha = ls(fun,x,f,g,p)    [fn]
%                             (default: backtracking Armijo below)
%
% OUTPUTS:
%   x    - approximate minimizer                               [n x 1]
%   info - struct: iters, gnorm (final), gnorm_hist [iters+1 x 1],
%          converged (logical), skips (# curvature-skipped updates)
%
% REFERENCES:
%   [1] Nocedal & Wright, Numerical Optimization, 2nd ed., Ch. 6.

    if nargin < 4, opts = struct(); end
    n     = numel(x0);
    H0    = getfield_default(opts, 'H0',    eye(n));
    tol   = getfield_default(opts, 'tol',   1e-8);
    maxit = getfield_default(opts, 'maxit', 1000);
    c1    = getfield_default(opts, 'c1',    1e-4);
    bt    = getfield_default(opts, 'bt',    0.5);
    ls    = getfield_default(opts, 'linesearch', ...
                             @(f_, x_, f0_, g_, p_) armijo(f_, x_, f0_, g_, p_, c1, bt));

    x = x0(:);
    [f, g] = fun(x);
    H = H0;
    gnorm_hist = zeros(maxit + 1, 1);
    gnorm_hist(1) = norm(g);
    skips = 0;  k = 0;  converged = false;

    while k < maxit
        if norm(g) <= tol, converged = true; break; end
        p = -H * g;                         % search direction (matvec)
        if g' * p >= 0                      % guard: H must give descent
            p = -g; H = eye(n);             % reset to steepest descent
        end
        alpha = ls(fun, x, f, g, p);        % step length (handle)
        xnew = x + alpha * p;
        [fnew, gnew] = fun(xnew);
        s = xnew - x;  y = gnew - g;
        if s' * y > 1e-12 * norm(s) * norm(y)   % curvature ok -> update
            H = update_fn(H, s, y);
        else
            skips = skips + 1;                  % skip to keep H PD
        end
        x = xnew; f = fnew; g = gnew; k = k + 1;
        gnorm_hist(k + 1) = norm(g);
    end

    info = struct('iters', k, 'gnorm', norm(g), ...
                  'gnorm_hist', gnorm_hist(1:k + 1), ...
                  'converged', converged, 'skips', skips);
end

function v = getfield_default(s, name, dflt)
    if isfield(s, name) && ~isempty(s.(name)), v = s.(name); else, v = dflt; end
end

function alpha = armijo(fun, x, f0, g, p, c1, bt)
% Backtracking line search enforcing the Armijo sufficient-decrease condition.
% (No curvature/Wolfe condition -- adequate for descent, but on a quadratic it
% forfeits the conjugate-direction finite termination; use an exact/Wolfe
% search there. See newton_quasinewton_bfgs.tex, globalization.)
    alpha = 1; gtp = g' * p;
    while true
        ft = fun(x + alpha * p);
        if ft <= f0 + c1 * alpha * gtp || alpha < 1e-16, break; end
        alpha = bt * alpha;
    end
end
