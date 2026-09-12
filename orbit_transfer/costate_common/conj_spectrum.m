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
%   of its median, or a determinant sign change between two samples, is a
%   candidate, and every candidate is LOCATED and CLASSIFIED rather than
%   counted. Class is by CONTIGUITY (Astra review 2026-09-11):
%
%     start     -- the cluster touches the first sample: the start-up
%                  transient, Phi_rv -> 0 at t = 0. Not refined; the
%                  interval it covers is reported as UNCOVERED (.nStart),
%                  pending a short-time argument;
%     endpoint  -- the cluster reaches the last inner sample. REFINED like
%                  an interior one, through t_f;
%     interior  -- anything else.
%
%   Every refined candidate is RESOLVED, not merely re-sampled: the window
%   is re-integrated on SHIFTED grids at 4x and 16x (a nested grid keeps
%   the same nearest node when the root lies within h/32 of it, so a true
%   zero can plateau -- the old two-level ratio test was not a zero-
%   exclusion test), then the minimum of sigma_6(t) is bracketed one fine
%   step either side of the finest argmin and LOCATED by golden section.
%   The located minimum, relative to the median, is judged against a
%   numerical floor (.zeroFloor [1e-7], a POLICY value for the STM's
%   accuracy, not a proven bound):
%     zero        <= zeroFloor, or a determinant sign change in the window;
%     near-miss   >= clearFactor [100] x zeroFloor: a positive minimum was
%                    located and the candidate is cleared;
%     unresolved  in between: blocks a PASS.
%   Multiplicity is the number of singular values at the floor at a
%   located zero (resolved rank loss), not a ratio of the two smallest.
%   The determinant's last sample (t_f) is not used for the SIGN count
%   (its sign is judged by ms_conjugate_test with a resolution rule); its
%   spectrum IS scanned. (doc/conjugate_research_memo_2026-09-10.md.)
%
%% Inputs:
%
%  z8                       [8 x 1]                 [lam0(7); tf]
%  rv0                      [6 x 1]                 departure state
%  Tmax, c, muStar          double                  as tfMinProp
%  opts                     struct (optional)
%   .K [24] segments, .nSub [8] samples per segment, .tolCollapse [1e-3]
%   relative dip in sigma_6 that counts as a candidate (sigma_6/sigma_5 at
%   a candidate is reported as .ratio, diagnostic only), .refine [4]
%   sub-sampling factor for refined
%   candidates (applied twice: 4x then 16x, shifted grids), .zeroFloor
%   [1e-7] located minimum / median at or below which a candidate is a
%   zero, .clearFactor [100] multiple of zeroFloor above which it is
%   cleared as a near-miss
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
%                                                   .class .kind (zero |
%                                                   near-miss | unresolved |
%                                                   start) .refineRatio
%                                                   .refineRatio2 .signChange
%                                                   .tMinOverTf .sigMinRel
%                                                   .nSmall .svAtMin)
%                                                   .nInteriorCand .nEndCand
%                                                   .nStart .nNearMiss .nZero
%                                                   .nUnresolved .multiplicity
%                                                   .clear (the gate: no
%                                                   zero, nothing unresolved,
%                                                   no multiplicity, no coarse
%                                                   sign change) .zeroFloor
%                                                   .clearFactor .svEnd .minRel
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

% INTERIOR sign changes only. The last sample is t_f, where the determinant
% is a product of graded values and its sign carries no information.
inner = 1:(N-1);
sg = sign(out.det(inner));
out.nInterior = nnz(diff(sg) ~= 0);
ix = find(diff(sg) ~= 0, 1);
if isempty(ix), out.tFirst = NaN; else, out.tFirst = out.t(ix); end

% CANDIDATES: dips of sigma_6 below tolCollapse of its median, clustered
% into events, each LOCATED, classified and -- unless it is the structural
% start-up transient -- RESOLVED: the window is re-sampled on SHIFTED finer
% grids (a nested grid keeps the same nearest node when the root lies within
% h/32 of it, so the sampled minimum can plateau at a true zero -- Astra
% review 2026-09-11) and the minimum of sigma_6(t) is then bracketed and
% located by golden section. The located minimum is judged against a
% numerical floor: at or below it the candidate is a ZERO; a clear factor
% above it, a NEAR-MISS with a measured positive minimum; between, it is
% UNRESOLVED and must block a PASS. A plateau is never taken as proof.
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
               'refineRatio', {}, 'refineRatio2', {}, 'signChange', {}, 'kMin', {}, ...
               'tMinOverTf', {}, 'sigMinRel', {}, 'nSmall', {}, 'svAtMin', {});
sg = sign(out.det(inner));
zeroFloor = d('zeroFloor', 1e-7);  clearFactor = d('clearFactor', 100);
for e = 1:numel(starts)
    a = starts(e);  b = ends(e);
    [~, im] = min(rel(a:b));  km = a + im - 1;
    cd = struct('tOverTf', out.t(km)/tf, 'rel', rel(km), ...
                'ratio', out.sv(6,km)/max(out.sv(5,km), realmin), ...
                'class', 'interior', 'kind', 'n/a', 'refineRatio', NaN, ...
                'refineRatio2', NaN, 'signChange', false, 'kMin', km, ...
                'tMinOverTf', NaN, 'sigMinRel', NaN, 'nSmall', NaN, 'svAtMin', nan(6,1));
    % CLASS is by CONTIGUITY, not by band: a cluster that touches the first
    % sample is the start-up transient (Phi_rv -> 0 at t = 0, structural,
    % not refined -- the interval it covers is reported as uncovered); one
    % that reaches the last inner sample is an ENDPOINT candidate and IS
    % refined, through t_f; anything else is interior.
    if a == 1,             cd.class = 'start';
    elseif b == numel(inner), cd.class = 'endpoint';
    end
    if ~strcmp(cd.class, 'start')
        lo = max(1, a-1);  hi = min(numel(inner), b+1);
        cd.signChange = any(diff(sg(lo:hi)) ~= 0);
        % window [k0, kEnd] in coarse samples, re-integrated from the stored
        % state one sample before it (or from t = 0)
        k0 = max(1, a-2);
        if k0 == 1, yw = [rv0; 1; lam0];  Pw = eye(14);  kStart = 0;
        else,       yw = yS(:, k0-1);     Pw = PhiS(:, :, k0-1);  kStart = k0-1;
        end
        kEnd = min(N, b+2);                      % an endpoint window runs to t_f
        nC = kEnd - kStart;
        [m1, ~]  = windowMin(yw, Pw, nC, dt, refine,   0.5, Pq, Tmax, c, muStar);
        [m2, t2] = windowMin(yw, Pw, nC, dt, refine^2, 0.5, Pq, Tmax, c, muStar);
        cd.refineRatio  = m1/max(out.sv(6, km), realmin);
        cd.refineRatio2 = m2/max(m1, realmin);
        % LOCATE the minimum: golden section on the bracket one fine step
        % either side of the finest-grid argmin (sigma_6 is V-shaped at a
        % simple zero and smooth at a near-miss: unimodal on that bracket)
        hFine = dt/refine^2;
        tA = max(0, t2 - hFine);  tB = min(nC*dt, t2 + hFine);
        [cd.sigMinRel, tStar, cd.svAtMin] = locateMin(yw, Pw, kStart*dt, tA, tB, med, Pq, Tmax, c, muStar);
        cd.tMinOverTf = tStar/tf;
        cd.nSmall = nnz(cd.svAtMin/med <= zeroFloor);
        if cd.signChange || cd.sigMinRel <= zeroFloor
            cd.kind = 'zero';
        elseif cd.sigMinRel >= clearFactor*zeroFloor
            cd.kind = 'near-miss';                % a POSITIVE minimum was located
        else
            cd.kind = 'unresolved';               % inside the numerical floor band
        end
    else
        cd.kind = 'start';
    end
    cands(end+1) = cd; %#ok<AGROW>
end
out.candidates = cands;
cls = {cands.class};  knd = {cands.kind};
isInt = strcmp(cls, 'interior');  isEnd = strcmp(cls, 'endpoint');
refined = isInt | isEnd;
out.nInteriorCand = nnz(isInt);
out.nEndCand      = nnz(isEnd);
out.nStart        = nnz(strcmp(cls, 'start'));
out.nNearMiss     = nnz(refined & strcmp(knd, 'near-miss'));
out.nZero         = nnz(refined & strcmp(knd, 'zero'));
out.nUnresolved   = nnz(refined & strcmp(knd, 'unresolved'));
% MULTIPLICITY: a located zero at which TWO OR MORE singular values sit at
% the floor -- resolved rank loss, not a ratio between two small values
out.multiplicity  = nnz(refined & strcmp(knd, 'zero') & [cands.nSmall] >= 2);
out.zeroFloor = zeroFloor;  out.clearFactor = clearFactor;
% the verdict the caller should gate on: no located zero, nothing
% unresolved, no multiplicity; near-misses with a located positive minimum
% are cleared, start transients are reported as uncovered, not cleared
out.clear = out.nZero == 0 && out.nUnresolved == 0 && out.multiplicity == 0 && out.nInterior == 0;
end

function [m, tAt] = windowMin(yw, Pw, nCoarse, dt, sub, shift, Pq, Tmax, c, muStar)
% WINDOWMIN  Smallest sigma_6 over nCoarse coarse steps re-integrated at
% `sub` sub-steps each, on a grid SHIFTED by `shift` sub-steps so that no
% coarse node is re-sampled.  INPUTS: yw; Pw; nCoarse; dt; sub; shift; Pq;
% Tmax; c; muStar.  OUTPUTS: m; tAt (time of the minimum, from t = 0 of
% the window's start state, i.e. absolute time = tWindowStart + tAt).
m = inf;  tAt = NaN;  h = dt/sub;  tRun = 0;
first = shift*h;
if first > 0
    [yw, P] = mintime_prop_seg(first, yw, true, Tmax, c, muStar);  Pw = P*Pw;  tRun = first;
    svq = specAt(yw, Pw, Pq, Tmax, c, muStar);
    if svq(6) < m, m = svq(6);  tAt = tRun; end
end
for q = 1:nCoarse*sub - 1
    [yw, P] = mintime_prop_seg(h, yw, true, Tmax, c, muStar);
    Pw = P*Pw;  tRun = tRun + h;
    svq = specAt(yw, Pw, Pq, Tmax, c, muStar);
    if svq(6) < m, m = svq(6);  tAt = tRun; end
end
end

function [sigRel, tStar, svStar] = locateMin(yw, Pw, tW, tA, tB, med, Pq, Tmax, c, muStar)
% LOCATEMIN  Golden-section minimisation of sigma_6(t) on [tA, tB] (times
% relative to the window start tW), each evaluation one propagation from
% the window's stored state.  INPUTS: yw; Pw; tW; tA; tB; med; Pq; Tmax;
% c; muStar.  OUTPUTS: sigRel (min sigma_6 / med); tStar (absolute time);
% svStar [6x1].
    function [s6, sv] = f(t)
        if t <= 0, y = yw;  P = Pw;
        else,      [y, Pp] = mintime_prop_seg(t, yw, true, Tmax, c, muStar);  P = Pp*Pw;
        end
        sv = specAt(y, P, Pq, Tmax, c, muStar);  s6 = sv(6);
    end
gr = (sqrt(5) - 1)/2;
a = tA;  b = tB;
x1 = b - gr*(b - a);  x2 = a + gr*(b - a);
f1 = f(x1);  f2 = f(x2);
for it = 1:40
    if f1 < f2, b = x2;  x2 = x1;  f2 = f1;  x1 = b - gr*(b - a);  f1 = f(x1);
    else,       a = x1;  x1 = x2;  f1 = f2;  x2 = a + gr*(b - a);  f2 = f(x2);
    end
    if (b - a) < 1e-9*max(tB, 1), break, end
end
if f1 < f2, tBest = x1; else, tBest = x2; end
[s6, svStar] = f(tBest);
sigRel = s6/med;  tStar = tW + tBest;
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
