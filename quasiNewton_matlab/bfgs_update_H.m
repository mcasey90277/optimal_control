function H1 = bfgs_update_H(H, s, y)
% BFGS_UPDATE_H  BFGS update of the inverse-Hessian approximation.
%
% Implements the BFGS update of the inverse Hessian H (approximating
% (grad^2 f)^{-1}) so that the inverse secant equation H1*y = s holds while
% positive definiteness is preserved (provided s'*y > 0). This is the
% "V H V' + rho s s'" form, derived from the DFP Hessian update by the
% B<->H, s<->y duality (see newton_quasinewton_bfgs.tex, eq. bfgsH).
%
% INPUTS:
%   H  - current inverse-Hessian approximation, symmetric PD [n x n]
%   s  - step  s = x_{k+1} - x_k                             [n x 1]
%   y  - gradient change  y = grad f(x_{k+1}) - grad f(x_k)  [n x 1]
%
% OUTPUTS:
%   H1 - updated inverse-Hessian approximation, symmetric PD [n x n]
%        satisfies H1*y = s exactly; H1 > 0 when H > 0 and s'*y > 0.
%
% NOTES:
%   Curvature condition s'*y > 0 is required for PD preservation. The caller
%   should guarantee it (Wolfe line search) or skip/damp the update; here we
%   error on a non-positive denominator rather than silently produce a
%   non-PD matrix.
%
% REFERENCES:
%   [1] Nocedal & Wright, Numerical Optimization, 2nd ed., eq. (6.17).
%   [2] Kanamori & Ohara, "A Bregman Extension of Quasi-Newton Updates II,"
%       arXiv:1010.2846 (2010), eq. (3) (B-form) and its inverse.

    sy = s' * y;
    if sy <= 0
        error('bfgs_update_H:curvature', ...
              'Curvature condition s''*y > 0 violated (s''*y = %.3e).', sy);
    end
    n   = numel(s);
    rho = 1 / sy;
    V   = eye(n) - rho * (s * y');
    H1  = V * H * V' + rho * (s * s');
    H1  = (H1 + H1') / 2;   % symmetrize against round-off
end
