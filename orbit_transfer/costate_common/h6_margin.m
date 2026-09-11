function H = h6_margin(z8, Tmax, c, opts)
%% Purpose:
%
%   H6: exclude the REDUCED problem's spurious-zero mechanism.
%
%   The 6-state reduction is non-autonomous, so its reduced Hamiltonian
%   h = p.F is NOT conserved. Since p(t)'J(t) = 0 along the variational
%   equations, the conjugate determinant can vanish because h(t) = 0 rather
%   than because the rank drops -- a zero that is not a loss of optimality.
%   The autonomous form excludes this with h == -1; the reduction does not,
%   and none of H1-H4 covers it (Astra proof review, 2026-09-10).
%
%   Here it closes in one line. From H = 1 + lambda.f = 0,
%
%       lambda_rv . f_rv = -1 - lambda_m*(-T/c) = -1 + (T/c) lambda_m,
%
%   which vanishes exactly at lambda_m = c/T. And lambda_m decreases
%   MONOTONICALLY (lam_m_dot = -T|lam_v|/m^2 < 0) from lambda_m(0) to
%   lambda_m(t_f) = 0, so its largest value on the arc is lambda_m(0):
%
%       the mechanism cannot fire  <==>  lambda_m(0) < c/T.
%
%   By the identity integral(T Q_mt dt) = lambda_m(0), this is equivalently
%   "the total strict-bang margin stays below c/T".
%
%   Measured over the shipped 70 mN catalog: c/T = 49.38 and lambda_m(0) in
%   [2.11, 10.88] across all 53 entries -- 4.5x headroom, none can fire it.
%
%% Inputs:
%
%  z8                       [8 x 1]                 [lam0(7); tf]; only
%                                                   lam0(7) = lambda_m(0) is
%                                                   used
%  Tmax, c                  double                  ND thrust acceleration
%                                                   and exhaust speed
%  opts                     struct (optional)
%   .marginMin [1] required lambda_m(0) < c/T ratio (1 = the bare condition;
%   raise it to demand headroom), .Hresid [0] the numerical Hamiltonian
%   residual |lambda.f + 1| on the arc: the reduced Hamiltonian is known only
%   to that accuracy, so h_max must clear zero by MORE than it
%
%% Outputs:
%
%  H                        struct                  .lamM0 .threshold (c/T)
%                                                   .margin (threshold/lamM0)
%                                                   .hMax (the largest
%                                                   lambda_rv.f_rv on the arc,
%                                                   attained at t = 0; the
%                                                   mechanism fires iff it
%                                                   reaches 0) .clearance
%                                                   (-hMax) .ok .reason
%
%% Revision History:
%  M. Casey                                                   (c) 09/10/2026
%  Copyright Coorbital Inc.
%% ------------------------ Begin Code Sequence ---------------------------

if nargin < 4, opts = struct(); end
marginMin = fieldd(opts, 'marginMin', 1);  Hresid = fieldd(opts, 'Hresid', 0);
lamM0 = z8(7);
thresh = c/Tmax;
hMax = -1 + (Tmax/c)*lamM0;
H = struct('lamM0', lamM0, 'threshold', thresh, 'margin', thresh/max(lamM0, realmin), ...
           'hMax', hMax, 'clearance', -hMax, 'ok', false, 'reason', '');

if ~(isfinite(lamM0) && isreal(lamM0))
    H.reason = 'lambda_m(0) is not a real finite number';  return
end
if lamM0 < 0
    % lambda_m decreases to lambda_m(t_f) = 0, so a negative start means the
    % transversality condition cannot hold: this is not an admissible lift
    H.reason = sprintf('lambda_m(0) = %.4g is NEGATIVE, yet it decreases to 0 at t_f', lamM0);
    H.margin = NaN;  return
end
if ~(H.margin > marginMin)
    H.reason = sprintf(['H6 FAILS: lambda_m(0) = %.4f against c/T = %.4f (margin %.4f). ' ...
                        'The reduced Hamiltonian can vanish on the arc, so a determinant ' ...
                        'zero need not be a conjugate point.'], lamM0, thresh, H.margin);
    return
end
if ~(H.clearance > Hresid)
    H.reason = sprintf(['H6 NOT ESTABLISHED: h_max = %.3e clears zero by %.3e, but the ' ...
                        'Hamiltonian is only known to %.1e -- the clearance is inside the ' ...
                        'numerical error'], hMax, H.clearance, Hresid);
    return
end
H.ok = true;
H.reason = sprintf(['H6 holds: lambda_m(0) = %.4f is %.1fx below c/T = %.4f, so ' ...
                    'lambda_rv.f_rv stays at or below %.4f (clearance %.3f > residual %.1e)'], ...
                   lamM0, H.margin, thresh, hMax, H.clearance, Hresid);
end

function v = fieldd(s, f, d_)
% FIELDD  Field with default.  INPUTS: s; f; d_.  OUTPUTS: v.
if isfield(s, f) && ~isempty(s.(f)), v = s.(f); else, v = d_; end
end
