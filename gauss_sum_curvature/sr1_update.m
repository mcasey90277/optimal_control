function [H1, denom] = sr1_update(H, s, y)
% SR1_UPDATE  Symmetric-rank-one update of a (possibly indefinite) Hessian estimate.
%
% The SR1 quasi-Newton update. Unlike BFGS/DFP it does NOT force positive
% definiteness, so it can represent and converge to an indefinite Hessian --
% exactly what is needed for measurement-model curvature (see demo_sr1_curvature).
%
%   H1 = H + (y - H s)(y - H s)' / ((y - H s)' s)
%
% INPUTS:
%   H - current symmetric Hessian estimate                 [n x n]
%   s - step in the operating point  (x_new - x_old)       [n x 1]
%   y - change in the gradient being matched               [n x 1]
%       (for a measurement component: change in that row of the Jacobian)
%
% OUTPUTS:
%   H1    - updated symmetric Hessian estimate             [n x n]
%   denom - the SR1 denominator (y - H s)' s; |denom| tiny => skip the update
%           (the SR1 breakdown), the only safeguard SR1 needs.
%
% REFERENCES:
%   [1] Nocedal & Wright, Numerical Optimization, 2nd ed., sec. 6.2 (SR1).

    v = y - H*s;
    denom = v' * s;
    H1 = H + (v * v') / denom;
    H1 = 0.5*(H1 + H1');   % symmetrize against round-off
end
