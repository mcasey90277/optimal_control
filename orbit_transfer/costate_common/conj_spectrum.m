function out = conj_spectrum(z8, rv0, Tmax, c, muStar, opts)
%% Purpose:
%
%   DENSE singular-spectrum scan of the free-time quotiented conjugate
%   matrix, closing the two blind spots of the sampled sign test:
%
%     * two conjugate times inside one segment -- invisible to a sign test
%       at the junctions, visible to a scan that samples inside them;
%     * an EVEN-MULTIPLICITY crossing -- no sign change at all, but k
%       singular values collapse together, which a spectrum sees and a
%       determinant cannot.
%
%   Built on the SAME variational integration the instrument uses
%   (mintime_prop_seg), sub-divided: a finite-difference reconstruction was
%   tried first on 2026-09-10 and is accurate to only ~2e-6 relative, which
%   on some entries is ABOVE the quantity being measured.
%
%   WHAT IT DOES NOT DO. It does not judge by the smallest singular value at
%   t_f. That value does not discriminate: a certified entry and a refuted
%   one were measured at 1.37e-7 and 1.39e-7 -- the same number. Near t_f the
%   whole spectrum is geometrically graded by the hyperbolic flow (2.4e0
%   down to 1.4e-7, each ~30x the next) and the determinant is their product,
%   ~1e-16, so its SIGN is meaningless there. The verdict here comes from
%   INTERIOR structure only, and the endpoint spectrum is reported for
%   information. (doc/conjugate_research_memo_2026-09-10.md.)
%
%% Inputs:
%
%  z8                       [8 x 1]                 [lam0(7); tf]
%  rv0                      [6 x 1]                 departure state
%  Tmax, c, muStar          double                  as tfMinProp
%  opts                     struct (optional)
%   .K [24] segments, .nSub [8] samples per segment, .tolCollapse [1e-3]
%   relative dip in sigma_min that counts as a candidate, .multTol [0.1]
%   sigma_5/sigma_6 below this at a candidate means MULTIPLICITY
%
%% Outputs:
%
%  out                      struct                  .t [1 x N] .sv [6 x N]
%                                                   .det [1 x N] .nInterior
%                                                   (sign changes strictly
%                                                   inside) .tFirst .tf
%                                                   .multiplicity (count of
%                                                   candidates where two or
%                                                   more values collapse)
%                                                   .svEnd [6 x 1] .minRel
%
%% Revision History:
%  M. Casey                                                   (c) 09/10/2026
%  Copyright Coorbital Inc.
%% ------------------------ Begin Code Sequence ---------------------------

if nargin < 6, opts = struct(); end
d = @(f,v) fieldd(opts, f, v);
K = d('K', 24);  nSub = max(2, d('nSub', 8));
tolCollapse = d('tolCollapse', 1e-3);  multTol = d('multTol', 0.1);
z8 = z8(:);  rv0 = rv0(:);
tf = z8(8);  lam0 = z8(1:7);

% the quotient: project out the scaling direction lambda(0)_{rv}
Pq = null(lam0(1:6).');                        % 6 x 5

dt = tf/(K*nSub);
y = [rv0; 1; lam0];
Phi = eye(14);
N = K*nSub;
out.t = zeros(1, N);  out.sv = zeros(6, N);  out.det = zeros(1, N);
for k = 1:N
    [y, P] = mintime_prop_seg(dt, y, true, Tmax, c, muStar);
    Phi = P*Phi;
    F = mintime_rhs_point(y, Tmax, c, muStar);
    M = [Phi(1:6, 8:13)*Pq, F(1:6)];
    nrm = max(vecnorm(M, 2, 1), realmin);
    M = M ./ nrm;                              % equilibrate: conjugacy is
    out.t(k) = k*dt;                           % dependence, not smallness
    out.sv(:, k) = svd(M);
    out.det(k) = det(M);
end
out.tf = tf;
out.svEnd = out.sv(:, end);

% INTERIOR sign changes only. The last sample is t_f, where the determinant
% is a product of graded values and its sign carries no information.
inner = 1:(N-1);
sg = sign(out.det(inner));
out.nInterior = nnz(diff(sg) ~= 0);
ix = find(diff(sg) ~= 0, 1);
if isempty(ix), out.tFirst = NaN; else, out.tFirst = out.t(ix); end

% MULTIPLICITY: at a dip, does more than one singular value collapse?
% sigma_5/sigma_6 near 1 means two directions went together.
rel = out.sv(6, inner) ./ max(median(out.sv(6, inner)), realmin);
out.minRel = min(rel);
cand = find(rel < tolCollapse);
out.multiplicity = 0;
for k = cand
    if out.sv(6,k)/max(out.sv(5,k), realmin) > multTol
        out.multiplicity = out.multiplicity + 1;
    end
end
end

function v = fieldd(s, f, d_)
% FIELDD  Field with default.  INPUTS: s; f; d_.  OUTPUTS: v.
if isfield(s, f) && ~isempty(s.(f)), v = s.(f); else, v = d_; end
end
