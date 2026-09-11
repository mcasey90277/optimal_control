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
%   This is a CANDIDATE-DETECTION scan. A dip of sigma_6 below tolCollapse
%   of its median, or a determinant sign change between two samples, is a
%   candidate, and every candidate is LOCATED and CLASSIFIED rather than
%   counted:
%
%     start     -- within the first segment, where Phi_rv -> 0 at t = 0 and
%                  the whole matrix is still growing from zero;
%     endpoint  -- within the last segment, where the hyperbolic flow grades
%                  the whole spectrum down together (a certified and a
%                  refuted entry measured sigma_min 1.37e-7 and 1.39e-7 at
%                  t_f: the endpoint value does not discriminate);
%     interior  -- the only class that can be a conjugate point. Each one is
%                  REFINED TWICE: the window around it is re-integrated at
%                  4x and then 16x the sampling and the sampled minima
%                  compared level by level. A zero keeps falling at BOTH
%                  levels (or carries a determinant sign change); a
%                  near-miss plateaus. One level is not enough: on the 70 mN
%                  sweep four entries read 0.49 at the first level against
%                  a 0.5 threshold and 0.8-1.0 at the second (FINDINGS 42),
%                  while the refuted entry's true crossing reads 0.06.
%
%   The determinant's last sample (t_f) is not used for the sign count
%   because that sample is the graded product above; the sign changes are
%   counted strictly inside. (doc/conjugate_research_memo_2026-09-10.md.)
%
%% Inputs:
%
%  z8                       [8 x 1]                 [lam0(7); tf]
%  rv0                      [6 x 1]                 departure state
%  Tmax, c, muStar          double                  as tfMinProp
%  opts                     struct (optional)
%   .K [24] segments, .nSub [8] samples per segment, .tolCollapse [1e-3]
%   relative dip in sigma_6 that counts as a candidate, .multTol [0.1]
%   sigma_6/sigma_5 ABOVE this at a candidate means two values collapsed
%   together (MULTIPLICITY), .refine [4] sub-sampling factor for interior
%   candidates (applied twice: 4x then 16x), .zeroRatio [0.5] level-to-
%   level minimum ratio; an interior candidate is a ZERO only if BOTH
%   levels fall below it
%
%% Outputs:
%
%  out                      struct                  .t [1 x N] .sv [6 x N]
%                                                   .det [1 x N] .nInterior
%                                                   (sign changes strictly
%                                                   inside) .tFirst .tf
%                                                   .candidates (struct
%                                                   array: .tOverTf .rel
%                                                   .ratio sigma_6/sigma_5
%                                                   .class .kind
%                                                   .refineRatio (4x / 1x)
%                                                   .refineRatio2 (16x / 4x)
%                                                   .signChange)
%                                                   .nInteriorCand .nNearMiss
%                                                   .nZero .multiplicity
%                                                   (INTERIOR candidates with
%                                                   two values collapsed)
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
refine = d('refine', 4);  zeroRatio = d('zeroRatio', 0.5);
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

% INTERIOR sign changes only. The last sample is t_f, where the determinant
% is a product of graded values and its sign carries no information.
inner = 1:(N-1);
sg = sign(out.det(inner));
out.nInterior = nnz(diff(sg) ~= 0);
ix = find(diff(sg) ~= 0, 1);
if isempty(ix), out.tFirst = NaN; else, out.tFirst = out.t(ix); end

% CANDIDATES: dips of sigma_6 below tolCollapse of its median, clustered
% into events, each located, classified and (if interior) refined.
med = max(median(out.sv(6, inner)), realmin);
rel = out.sv(6, inner) ./ med;
out.minRel = min(rel);
isCand = rel < tolCollapse;
% a determinant SIGN CHANGE between two inner samples is a candidate by
% definition, whether or not sigma_6 dipped below the depth threshold at
% either sample: a simple zero straddled by two samples need not
sg0 = sign(out.det(inner));
for ix = find(diff(sg0) ~= 0)
    isCand(ix:ix+1) = true;
end
edges = diff([false, isCand, false]);
starts = find(edges == 1);  ends = find(edges == -1) - 1;
cands = struct('tOverTf', {}, 'rel', {}, 'ratio', {}, 'class', {}, 'kind', {}, ...
               'refineRatio', {}, 'refineRatio2', {}, 'signChange', {}, 'kMin', {});
sg = sign(out.det(inner));
for e = 1:numel(starts)
    a = starts(e);  b = ends(e);
    [~, im] = min(rel(a:b));  km = a + im - 1;
    cd = struct('tOverTf', out.t(km)/tf, 'rel', rel(km), ...
                'ratio', out.sv(6,km)/max(out.sv(5,km), realmin), ...
                'class', 'interior', 'kind', 'n/a', 'refineRatio', NaN, ...
                'refineRatio2', NaN, 'signChange', false, 'kMin', km);
    if a <= nSub,              cd.class = 'start';
    elseif b >= (N-1) - nSub,  cd.class = 'endpoint';
    end
    if strcmp(cd.class, 'interior')
        % a determinant sign change inside the window is a zero outright
        lo = max(1, a-1);  hi = min(numel(inner), b+1);
        cd.signChange = any(diff(sg(lo:hi)) ~= 0);
        % re-integrate the window at `refine`x the sampling from the sample
        % before it (state and STM were kept), compare the sampled minimum
        k0 = max(1, a-2);
        if k0 == 1, yw = [rv0; 1; lam0];  Pw = eye(14);  kStart = 0;
        else,       yw = yS(:, k0-1);     Pw = PhiS(:, :, k0-1);  kStart = k0-1;
        end
        kEnd = min(N-1, b+2);
        m1 = windowMin(yw, Pw, kEnd - kStart, dt, refine,   Pq, Tmax, c, muStar);
        m2 = windowMin(yw, Pw, kEnd - kStart, dt, refine^2, Pq, Tmax, c, muStar);
        cd.refineRatio  = m1/max(out.sv(6, km), realmin);
        cd.refineRatio2 = m2/max(m1, realmin);
        if cd.signChange || (cd.refineRatio < zeroRatio && cd.refineRatio2 < zeroRatio) ...
                         || m2/med < 1e-8
            cd.kind = 'zero';
        else
            cd.kind = 'near-miss';
        end
    end
    cands(end+1) = cd; %#ok<AGROW>
end
out.candidates = cands;
isInt = strcmp({cands.class}, 'interior');
out.nInteriorCand = nnz(isInt);
out.nNearMiss = nnz(isInt & strcmp({cands.kind}, 'near-miss'));
out.nZero     = nnz(isInt & strcmp({cands.kind}, 'zero'));
% MULTIPLICITY: an INTERIOR candidate where sigma_5 collapsed with sigma_6
out.multiplicity = nnz(isInt & [cands.ratio] > multTol);
end

function m = windowMin(yw, Pw, nCoarse, dt, sub, Pq, Tmax, c, muStar)
% WINDOWMIN  Smallest sigma_6 over nCoarse coarse steps re-integrated at
% `sub` sub-steps each, from state yw / STM Pw.  INPUTS: yw; Pw; nCoarse;
% dt; sub; Pq; Tmax; c; muStar.  OUTPUTS: m.
m = inf;
for q = 1:nCoarse*sub
    [yw, P] = mintime_prop_seg(dt/sub, yw, true, Tmax, c, muStar);
    Pw = P*Pw;
    svq = specAt(yw, Pw, Pq, Tmax, c, muStar);
    m = min(m, svq(6));
end
end

function [sv, dt_] = specAt(y, Phi, Pq, Tmax, c, muStar)
% SPECAT  Equilibrated spectrum and determinant of the quotiented conjugate
% matrix at one sample.  INPUTS: y; Phi; Pq; Tmax; c; muStar.
% OUTPUTS: sv [6x1]; dt_.
F = mintime_rhs_point(y, Tmax, c, muStar);
M = [Phi(1:6, 8:13)*Pq, F(1:6)];
nrm = max(vecnorm(M, 2, 1), realmin);
M = M ./ nrm;                                  % equilibrate: conjugacy is
sv = svd(M);  dt_ = det(M);                    % dependence, not smallness
end

function v = fieldd(s, f, d_)
% FIELDD  Field with default.  INPUTS: s; f; d_.  OUTPUTS: v.
if isfield(s, f) && ~isempty(s.(f)), v = s.(f); else, v = d_; end
end
