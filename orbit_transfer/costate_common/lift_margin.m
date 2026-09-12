function M = lift_margin(C1, C2, lam, opts)
%% Purpose:
%
%   The rank statement behind dim S, as a MEASURED MARGIN rather than a
%   threshold verdict.
%
%   dim S >= 1 is CONSTRUCTIVE: a lift is exhibited, and its residual is
%   reported -- as a BACKWARD error |C lam|/(|C||lam|) against liftTol (a
%   small residual constructs a nearby matrix with a null vector; it does
%   not construct an exact null vector of the continuous problem).
%   dim S <= 1 is a theorem plus a measurement. For an m x 7 matrix,
%   Eckart-Young says sigma_6 is the distance to the nearest matrix of rank
%   AT MOST FIVE (sigma_7 is the distance to the nearest rank-deficient
%   one), so
%
%       sigma_6(C_num) > ||C_exact - C_num||   ==>   rank(C_exact) >= 6,
%
%   and rank EXACTLY six -- dim S = 1 -- follows only with the separate
%   fact that an exact nonzero lift lies in the kernel, which the exhibited
%   normal extremal supplies up to its own residual (Astra review #2,
%   2026-09-11: the earlier header claimed "rank(C) = 6 exactly" from
%   sigma_6 alone). ||dC|| is the numerical error in C. It is MEASURED by
%   rebuilding C at a second numerical setting (a looser integration
%   tolerance, a different sample count) and taking the difference. The
%   output is the margin sigma_6/||dC||: how many times larger the smallest
%   retained singular value is than the uncertainty in the matrix itself.
%
%   This replaces a rule whose own header claimed a safety guarantee it did
%   not provide -- capping a tolerance at 1e-3*sigma_1 does not stop a noisy
%   lift from reporting a spurious nullity (Astra review, 2026-09-10).
%
%   NOTE ON RIGOUR. A two-setting difference is a SENSITIVITY estimate of
%   ONE error component (the frozen-control adjoint integration), not a
%   total error bound: the flown trajectory, its pchip interpolants, the
%   endpoint data and the matrix assembly are common to both builds and
%   cancel in the difference. Turning this into a proof needs a total
%   error enclosure; the margin is what makes that gap visible and
%   quantitative instead of hidden in a constant.
%
%% Inputs:
%
%  C1                       [m x 7]                 constraint matrix at the
%                                                   TIGHTER setting
%  C2                       [m x 7]                 the same at a looser one
%  lam                      [7 x 1]                 the exhibited lift
%  opts                     struct (optional)
%   .marginMin [10] margin needed to certify, .liftTol [1e-6] BACKWARD
%   error |C lam|/(sigma_1 |lam|) for `lam` to count as a lift at all (the
%   same normalisation mintime_hypothesis_gates reports as nullResidRel)
%
%% Outputs:
%
%  M                        struct                  .sigma6 .errEst .margin
%                                                   .dimS .certified .reason
%                                                   .nullResid .sv
%
%% Revision History:
%  M. Casey                                                   (c) 09/10/2026
%  Copyright Coorbital Inc.
%% ------------------------ Begin Code Sequence ---------------------------

if nargin < 4, opts = struct(); end
d = @(f,v) fieldd(opts, f, v);
marginMin = d('marginMin', 10);  liftTol = d('liftTol', 1e-6);
assert(isequal(size(C1), size(C2)), 'lift_margin:shape', 'the two builds must have the same shape');
assert(isnumeric(C1) && isreal(C1) && all(isfinite(C1(:))) && isnumeric(C2) && isreal(C2) && all(isfinite(C2(:))), ...
       'lift_margin:matrix', 'both constraint matrices must be real and finite');
assert(size(C1, 2) == 7 && size(C1, 1) >= 7, 'lift_margin:matrix', 'constraint matrices must be [m x 7] with m >= 7');
assert(isscalar(marginMin) && isfinite(marginMin) && marginMin > 0 && isscalar(liftTol) && isfinite(liftTol) && liftTol > 0, ...
       'lift_margin:opts', 'marginMin and liftTol must be positive finite scalars');
lam = lam(:);
assert(numel(lam) == 7, 'lift_margin:lift', 'the lift must have 7 components');

M = struct('sigma6', NaN, 'errEst', NaN, 'margin', NaN, 'dimS', NaN, ...
           'certified', false, 'reason', '', 'nullResid', NaN, 'sv', []);

sv = svd(C1);  M.sv = sv;
M.sigma6 = sv(6);
M.errEst = norm(C1 - C2);
% a lift must be a real, finite, NONZERO vector before its residual means
% anything: the zero vector has residual zero and proves nothing
if ~(isreal(lam) && all(isfinite(lam)))
    M.reason = 'the supplied lift is not real and finite';  return
end
if norm(lam) == 0
    M.reason = 'the supplied lift is the zero vector: no lift is exhibited';  return
end
M.nullResid = norm(C1*lam)/norm(lam);

% the CONSTRUCTIVE half: is the supplied vector actually a lift?
if ~(M.nullResid <= liftTol*max(sv(1), realmin))
    M.reason = sprintf(['the supplied vector is not a lift (backward error %.2e > ' ...
                        '%.0e): dim S >= 1 is not established'], ...
                       M.nullResid/max(sv(1), realmin), liftTol);
    M.dimS = nnz(sv <= max(M.errEst, eps*sv(1)));
    return
end

% the rank half, by Eckart-Young against the MEASURED error
floor_ = max(M.errEst, eps*sv(1));
M.dimS = nnz(sv <= floor_);
M.margin = M.sigma6/max(floor_, realmin);
if M.dimS ~= 1
    M.reason = sprintf(['dim S = %d at the measured error floor %.2e ' ...
                        '(sigma_6 = %.2e, sigma_7 = %.2e)'], M.dimS, floor_, sv(6), sv(7));
    return
end
if M.margin < marginMin
    M.reason = sprintf(['dim S = 1 but only by a margin of %.1f (sigma_6 = %.2e, ' ...
                        'error %.2e): not certified'], M.margin, M.sigma6, M.errEst);
    return
end
M.certified = true;
M.reason = sprintf(['nullity one NUMERICALLY SUPPORTED (rank >= 6 by Eckart-Young, dim S = 1 with the ' ...
                    'exhibited lift): sigma_6 = %.2e exceeds the measured error %.2e by %.0fx, and ' ...
                    'the lift''s backward error is %.1e'], ...
                   M.sigma6, M.errEst, M.margin, M.nullResid/max(sv(1), realmin));
end

function v = fieldd(s, f, d_)
% FIELDD  Field with default.  INPUTS: s; f; d_.  OUTPUTS: v.
if isfield(s, f) && ~isempty(s.(f)), v = s.(f); else, v = d_; end
end
