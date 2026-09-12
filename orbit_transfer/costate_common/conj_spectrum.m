function out = conj_spectrum(z8, rv0, Tmax, c, muStar, opts)
%% Purpose:
%
%   DENSE singular-spectrum scan of the free-time quotiented conjugate
%   matrix, adding sensitivity to the two blind spots of the sampled sign
%   test (it does not close them: no between-sample bound is proven):
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
%   This is a CANDIDATE-DETECTION scan. A dip of sigma_6 below tolCollapse
%   of its median, a TRUSTED determinant sign change between two samples,
%   or a vanishing raw column norm is a candidate; candidates are clustered
%   by contiguity (Astra reviews #2 and #3, 2026-09-11):
%
%     start     -- the cluster touching the first sample: the start-up
%                  transient, Phi_rv -> 0 at t = 0. Its extent is reported
%                  as UNCOVERED (.tUncovered); a local minimum INSIDE it is
%                  a candidate of its own and is resolved;
%     endpoint  -- the cluster reaching t_f (the last sample IS a candidate
%                  source and IS evaluated);
%     interior  -- anything else.
%
%   The resolution is conj_resolve (a pure function of M(t), so synthetic
%   matrices drive it in tests/test_conj_resolve). It KEEPS every
%   evaluation: shifted grids at 4x and 16x, then EVERY local minimum
%   bracketed and located by golden section with an a-posteriori
%   unimodality check; a trusted opposite-sign bracket anywhere in the
%   record is an established root; the smallest value ever seen is never
%   upgraded by a later search; a vanishing RAW column norm is a candidate
%   (column normalisation hides that rank loss). Per candidate:
%     zero        an established sign bracket, or the smallest value at or
%                    below zeroFloor [1e-7] x median (floor-level: possible
%                    rank loss). Both block;
%     unresolved  a vanishing column, a non-unimodal bracket, or a minimum
%                    inside the floor band: blocks;
%     near-miss   a positive minimum located, >= clearFactor [100] x floor.
%   .nSmall is the numerical corank at the located minimum, not a root
%   multiplicity. The floor is a POLICY value for the STM's accuracy, not a
%   proven bound. A scan that never leaves its clusters with a trusted sign
%   is NOT clear (.testable). (doc/conjugate_research_memo_2026-09-10.md.)
%
%% Inputs:
%
%  z8                       [8 x 1]                 [lam0(7); tf]
%  rv0                      [6 x 1]                 departure state
%  Tmax, c, muStar          double                  as tfMinProp
%  opts                     struct (optional)
%   .K [24] segments, .nSub [8] samples per segment, .tolCollapse [1e-3]
%   relative dip in sigma_6 that counts as a candidate, .refine [4]
%   sub-sampling factor (applied twice: 4x then 16x, shifted grids),
%   .zeroFloor [1e-7] smallest value / median at or below which a candidate
%   is a zero, .clearFactor [100] multiple of zeroFloor above which it is
%   cleared as a near-miss, .colTol [1e-10] raw column-norm ratio below
%   which a column is vanishing, .signTol [1e-10] sigma ratio above which
%   a determinant sign is trusted
%
%% Outputs:
%
%  out                      struct                  .t [1 x N] .sv [6 x N]
%                                                   .det [1 x N] .nInterior
%                                                   (sign changes strictly
%                                                   inside) .tFirst .tf
%                                                   .candidates (struct
%                                                   array, conj_resolve
%                                                   schema: .class .kind
%                                                   .tOverTf .tMinOverTf
%                                                   .sigMinRel .nSmall
%                                                   .nBrackets .unimodal
%                                                   .colBad .tWindow ...)
%                                                   .nInteriorCand .nEndCand
%                                                   .nStart .tUncovered
%                                                   .nZero .nZeroSign
%                                                   .nZeroFloor .nNearMiss
%                                                   .nUnresolved
%                                                   .multiplicity .testable
%                                                   .clear (the gate)
%                                                   .reason .nEval
%                                                   .sigRatio .colRatio
%                                                   .trusted .med
%                                                   .specConsistency
%                                                   .zeroFloor .clearFactor
%                                                   .svEnd .minRel
%
%% Revision History:
%  M. Casey                                                   (c) 09/10/2026
%  Copyright Coorbital Inc.
%% ------------------------ Begin Code Sequence ---------------------------

if nargin < 6, opts = struct(); end
d = @(f,v) fieldd(opts, f, v);
K = d('K', 24);  nSub = max(2, d('nSub', 8));
tolCollapse = d('tolCollapse', 1e-3);
refine = d('refine', 4);
z8 = z8(:);  rv0 = rv0(:);
tf = z8(8);  lam0 = z8(1:7);

% the quotient: project out the scaling direction lambda(0)_{rv}
Pq = null(lam0(1:6).');                        % 6 x 5

dt = tf/(K*nSub);
y = [rv0; 1; lam0];
Phi = eye(14);
N = K*nSub;
out.t = zeros(1, N);  out.sv = zeros(6, N);  out.det = zeros(1, N);
yS = zeros(14, N);  PhiS = zeros(14, 14, N);   % kept for the local refinement
for k = 1:N
    [y, P] = mintime_prop_seg(dt, y, true, Tmax, c, muStar);
    Phi = P*Phi;
    yS(:, k) = y;  PhiS(:, :, k) = Phi;
    out.t(k) = k*dt;
    [out.sv(:, k), out.det(k)] = specAt(y, Phi, Pq, Tmax, c, muStar);
end
out.tf = tf;
out.svEnd = out.sv(:, end);

% THE RESOLUTION is conj_resolve, a pure function of M(t) so that synthetic
% matrices can drive it in tests. M(t) here: propagate from the stored coarse
% sample at or before t (or from t = 0), then assemble the quotiented block
% with the flow column. A coarse time re-uses its stored sample exactly.
    function M = Mat(t)
        kk = floor(t/dt + 1e-9);  kk = min(max(kk, 0), N);   % (not the loop's k: nested workspace)
        if kk == 0, yk = [rv0; 1; lam0];  Pk = eye(14);  else, yk = yS(:, kk);  Pk = PhiS(:, :, kk); end
        tau = t - kk*dt;
        if tau > 1e-12*max(tf, 1)
            [yk, P] = mintime_prop_seg(tau, yk, true, Tmax, c, muStar);  Pk = P*Pk;
        end
        F = mintime_rhs_point(yk, Tmax, c, muStar);
        M = [Pk(1:6, 8:13)*Pq, F(1:6)];
    end
ro = struct('tolCollapse', tolCollapse, 'refine', refine, ...
            'zeroFloor', d('zeroFloor', 1e-7), 'clearFactor', d('clearFactor', 100), ...
            'colTol', d('colTol', 1e-10), 'signTol', d('signTol', 1e-10));
R = conj_resolve(@Mat, out.t, ro);
% the coarse spectrum from the resolver is the SAME evaluation as the loop
% above (stored samples, no propagation); keep the resolver's record as the
% one source and carry the loop's as a consistency check
out.specConsistency = max(abs(R.sv(:) - out.sv(:))) / max(max(abs(out.sv(:))), realmin);
for f = {'sv', 'det', 'sigRatio', 'colRatio', 'trusted', 'med', 'minRel', 'nInterior', 'tFirst', ...
         'candidates', 'nStart', 'tUncovered', 'nInteriorCand', 'nEndCand', 'nZero', 'nZeroSign', ...
         'nZeroFloor', 'nNearMiss', 'nUnresolved', 'multiplicity', 'testable', 'clear', 'reason', ...
         'nEval', 'zeroFloor', 'clearFactor'}
    out.(f{1}) = R.(f{1});
end
end

function [sv, dt_] = specAt(y, Phi, Pq, Tmax, c, muStar)
% SPECAT  Column-normalised spectrum and determinant of the quotiented
% conjugate matrix at one stored sample (the coarse pass; the resolver
% recomputes the same thing through Mat).  INPUTS: y; Phi; Pq; Tmax; c;
% muStar.  OUTPUTS: sv [6x1]; dt_.
F = mintime_rhs_point(y, Tmax, c, muStar);
M = [Phi(1:6, 8:13)*Pq, F(1:6)];
nrm = max(vecnorm(M, 2, 1), realmin);
M = M ./ nrm;
sv = svd(M);  dt_ = det(M);
end

function v = fieldd(s, f, d_)
% FIELDD  Field with default.  INPUTS: s; f; d_.  OUTPUTS: v.
if isfield(s, f) && ~isempty(s.(f)), v = s.(f); else, v = d_; end
end
