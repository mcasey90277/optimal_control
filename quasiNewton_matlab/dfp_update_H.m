function H1 = dfp_update_H(H, s, y)
% DFP_UPDATE_H  DFP update of the inverse-Hessian approximation.
%
% Implements the DFP update of the inverse Hessian H so that the inverse
% secant equation H1*y = s holds while positive definiteness is preserved
% (provided s'*y > 0). DFP is the dual of BFGS under the B<->H, s<->y swap;
% on the inverse Hessian it is the rank-two formula below.
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
%   Curvature condition s'*y > 0 is required (see bfgs_update_H). Empirically
%   DFP is more sensitive than BFGS to an inexact line search (Kanamori &
%   Ohara, Table 3) -- the motivation for the robustness comparison.
%
% REFERENCES:
%   [1] Nocedal & Wright, Numerical Optimization, 2nd ed., eq. (6.15) (B-form);
%       the inverse form is the rank-two update used here.
%   [2] Kanamori & Ohara, arXiv:1010.2846 (2010), eq. (2).

    sy = s' * y;
    if sy <= 0
        error('dfp_update_H:curvature', ...
              'Curvature condition s''*y > 0 violated (s''*y = %.3e).', sy);
    end
    Hy  = H * y;
    yHy = y' * Hy;
    H1  = H + (s * s') / sy - (Hy * Hy') / yHy;
    H1  = (H1 + H1') / 2;   % symmetrize against round-off
end
