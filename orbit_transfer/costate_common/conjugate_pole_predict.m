function W = conjugate_pole_predict(s, condJ, opts)
%% Purpose:
%
%   EARLY WARNING that a continuation is walking into a conjugate point,
%   fitted from the conditioning the solver already reports at every step.
%
%   A conjugate point is a non-trivial solution of the linearised
%   boundary-value problem, so the multiple-shooting Jacobian goes singular
%   there and cond(J) has a POLE, not an exponential. That distinction is
%   the whole content of this function: on the measured DRO -> tulip
%   departure walk (FINDINGS 39), log-linear extrapolation from data ending
%   at sD = 0.0450 puts the singularity at 0.050, while the pole fit puts it
%   at 0.0469 against a measured 0.04665.
%
%   Model:  log(cond) = a - p*log(sPole - s),   p > 0
%   fitted on the last nFit points by eliminating (a, p): the ratio of
%   successive log-differences depends on sPole alone, so a bracketed
%   scalar solve gives it, and (a, p) follow by least squares.
%
%   The warning is deliberately conservative: it fires only when the
%   conditioning is both RISING and CONVEX in the log, when the fitted
%   exponent is positive, and when the pole lies ahead of the data but
%   within reach. Anything else returns .ok = false with a reason -- an
%   early warning that cries wolf on a benign stretch is worse than none.
%
%% Inputs:
%
%  s                        [1 x n]                 continuation parameter,
%                                                   monotone
%  condJ                    [1 x n]                 cond(J) at each point
%  opts                     struct (optional)
%   .nFit [4] points used (>= 3), .minDecades [1] rise required over them,
%   .maxAhead [inf] refuse a pole further than this beyond the last point
%
%% Outputs:
%
%  W                        struct                  .ok .reason .sPole
%                                                   .exponent .decades
%                                                   .residual (of the log
%                                                   fit) .sLast
%
%% Revision History:
%  M. Casey                                                   (c) 09/09/2026
%  Copyright Coorbital Inc.
%% ------------------------ Begin Code Sequence ---------------------------

if nargin < 3, opts = struct(); end
d = @(f,v) fieldd(opts, f, v);
nFit = max(3, d('nFit', 4));  minDec = d('minDecades', 1);  maxAhead = d('maxAhead', inf);
s = s(:).';  condJ = condJ(:).';

W = struct('ok', false, 'reason', '', 'sPole', NaN, 'exponent', NaN, ...
           'decades', NaN, 'residual', NaN, 'sLast', NaN);
if numel(s) < 3 || numel(condJ) ~= numel(s)
    W.reason = sprintf('need at least 3 points (got %d)', numel(s));  return
end
if any(~isfinite(condJ)) || any(condJ <= 0)
    W.reason = 'cond(J) must be finite and positive';  return
end
k = min(nFit, numel(s));
x = s(end-k+1:end);  y = log(condJ(end-k+1:end));
W.sLast = x(end);
dir = sign(x(end) - x(1));
if dir == 0 || any(dir*diff(x) <= 0), W.reason = 's must be monotone'; return, end
x = dir*x;                                     % walk left to right
W.decades = (y(end) - y(1))/log(10);
if W.decades < minDec
    W.reason = sprintf('conditioning rose only %.2f decades over the last %d points', W.decades, k);
    return
end
% Convexity must be measured with DIVIDED differences: the continuation's
% samples are not evenly spaced (the walk halves its step near trouble), and
% diff(diff(y)) on uneven x is not a curvature -- it rejected the very series
% this function was built from.
sl = diff(y(:))./diff(x(:));
if any(diff(sl) < 0)
    W.reason = 'log cond(J) is not convex: rising, but not toward a pole';  return
end

% sPole from the log-difference ratios alone
f = @(sp) poleResid(sp, x, y);
lo = x(end) + 1e-12*max(abs(x(end)), 1);
hi = x(end) + 10*(x(end) - x(1)) + eps;
try
    sp = fzero(f, [lo + (hi-lo)*1e-6, hi]);
catch
    W.reason = 'no pole bracketed ahead of the data';  return
end
A = [ones(k,1), -log(sp - x(:))];
c = A\y(:);
W.exponent = c(2);
W.residual = norm(A*c - y(:), inf);
W.sPole = dir*sp;
if ~(W.exponent > 0)
    W.reason = sprintf('fitted exponent %.2f is not positive', W.exponent);  W.ok = false;  return
end
ahead = abs(W.sPole - W.sLast);
if ahead > maxAhead
    W.reason = sprintf('pole predicted %.4g ahead, beyond maxAhead %.4g', ahead, maxAhead);  return
end
W.ok = true;
W.reason = sprintf(['cond(J) rose %.2f decades and is convex in the log; pole at %.5f ' ...
                    '(%.4g ahead), exponent %.2f, fit residual %.1e'], ...
                   W.decades, W.sPole, ahead, W.exponent, W.residual);
end

function r = poleResid(sp, x, y)
% POLERESID  Residual of the log-difference-ratio condition at a trial pole.
% With log(cond) = a - p log(sp - s), successive differences of y are
% p*log((sp-s_i)/(sp-s_{i+1})), so y's differences must be PROPORTIONAL to
% those logs -- a condition on sp alone. INPUTS: sp; x; y. OUTPUTS: r.
if sp <= x(end), r = 1e6; return, end
u = log(sp - x(:));
dy = diff(y(:));  du = -diff(u);
if any(du <= 0), r = 1e6; return, end
p = (du'*dy)/(du'*du);                 % best p in least squares
r = (dy(end)/du(end)) - (dy(1)/du(1)); % ends agree only at the right pole
r = r / max(p, eps);
end

function v = fieldd(s, f, d_)
% FIELDD  Field with default.  INPUTS: s; f; d_.  OUTPUTS: v.
if isfield(s, f) && ~isempty(s.(f)), v = s.(f); else, v = d_; end
end
