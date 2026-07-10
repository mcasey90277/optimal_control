% DEMO_PROB1  BFGS vs DFP (inverse-Hessian form) on the convex quadratic.
%
% Problem 1 of Kanamori & Ohara (arXiv:1010.2846, sec. 6):
%   min  f(x) = 1/2 x'A x - e'x,   A = tridiag(-1, 2, -1),  e = ones(n,1).
% Strongly convex; unique minimizer x* = A\e; gradient Ax - e; Hessian A.
% This is the minimal core check that bfgs_update_H / dfp_update_H drive a
% working quasi-Newton loop and both reach x* on the easy (quadratic) case.
%
% REFERENCES:
%   Kanamori & Ohara, arXiv:1010.2846 (2010), Problem 1.

clear; clc;
rng(0);                                   % reproducible x0

n = 20;
A = full(spdiags(ones(n,1)*[-1 2 -1], -1:1, n, n));   % SPD tridiagonal
e = ones(n, 1);
fun = @(x) prob1_obj(x, A, e);                        % [f] or [f, grad]

xstar = A \ e;                            % exact minimizer
x0 = sqrt(10) * randn(n, 1);              % x0 ~ N(0, 10 I)

% Exact line search for the quadratic: along p, alpha* = -(g'p)/(p'A p).
% (Exact search restores the conjugate-direction property -> both BFGS and DFP
%  terminate in <= n steps on a quadratic. The paper approximates this with
%  fminbnd; the inexact-search robustness gap is the next experiment.)
exact_ls = @(fun, x, f, g, p) -(g' * p) / (p' * A * p);
opts = struct('tol', n*1e-10, 'maxit', 2000, 'linesearch', exact_ls);

[xB, infoB] = qn_minimize_H(fun, x0, @bfgs_update_H, opts);
[xD, infoD] = qn_minimize_H(fun, x0, @dfp_update_H,  opts);

fprintf('\nProblem 1: convex quadratic, n = %d (exact line search)\n', n);
fprintf('%-6s  %6s  %12s  %14s  %6s\n', ...
        'method','iters','final |g|','||x-x*||','skips');
fprintf('%-6s  %6d  %12.3e  %14.3e  %6d\n', ...
        'BFGS', infoB.iters, infoB.gnorm, norm(xB-xstar), infoB.skips);
fprintf('%-6s  %6d  %12.3e  %14.3e  %6d\n', ...
        'DFP',  infoD.iters, infoD.gnorm, norm(xD-xstar), infoD.skips);
fprintf('both converged: BFGS=%d  DFP=%d\n', infoB.converged, infoD.converged);
