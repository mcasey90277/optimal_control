function out = ms_conjugate_test(info, spec)
%% Purpose:
%
%   CONJUGATE-POINT TEST on a converged multiple-shooting extremal -- the
%   first piece of second-order optimality checking this pipeline has had.
%   First-order (PMP) conditions admit maxima and saddle extremals too; a
%   conjugate point in (0, tf) means the extremal STOPS being locally
%   minimizing there (Jacobi's necessary condition).
%
%   Method (free-final-time form, after Caillau/Bonnard's hampath): state
%   variations due to initial-costate variations follow
%   dx(t) = Phi_xl(t,0) dl(0). Two EXACT degeneracies must be handled
%   before any determinant is meaningful (both measured as 1e-37..1e-13
%   noise dets with spurious sign flips when ignored):
%     1. COSTATE SCALING: the min-time control -lam_v/|lam_v| is invariant
%        under scaling of (lam_r, lam_v), so dl(0) along the converged
%        costate produces dx == 0 identically. Columns are restricted to
%        an orthonormal complement P [6 x 5] of that direction.
%     2. TIME REPARAMETRIZATION (free tf): variations along the flow f(t)
%        are trajectory-preserving. The flow completes the square:
%        the monitored quantity is det([Phi_xl(t,0)*P, f(t)]) [6 x 6].
%   A sign change of that det strictly inside (0, tf) is a conjugate
%   point. The segment STMs come free from the ms_bvp Jacobian machinery
%   (opts.keepSTMs); junction states and flows come from info.Y and
%   spec.flow.
%
%  ASSUMPTIONS / NOTES:
%
% • The mass costate lam_m is EXCLUDED (columns default 8:13): in constant-
%   Isp all-burn flight the control never depends on lam_m and mass is a
%   known function of time, so the lam_m column of Phi_xl is exactly zero
%   effect -- a third degeneracy if included.
% • Junction resolution: the det is sampled at junctions 2..K+1 (= tf),
%   so a conjugate PAIR closer together than one segment can hide. K =
%   12-48 in the catalogs; treat "no sign change" as strong but not
%   airtight.
% • det is scaled to det^(1/m) magnitude for readable reporting; the SIGN
%   is taken from the equilibrated block (robust to bad scaling).
% • INITIAL COAST: until the control has acted, the state block is
%   structurally zero (sigma_min/sigma_max <= spec.rankTol [1e-13]);
%   those samples are skipped, .firstFullRank reports where testing
%   starts. (Review 2026-09-05: the old code counted them as focal points.)
% • A sign change on the bracket ending at the last sample: if that sample
%   is t_f AND its determinant is RESOLVED (sigma_min/sigma_max above
%   spec.resolvedTol, so its sign is trustworthy), then opposite nonzero
%   signs put the zero STRICTLY inside (t_K, t_f) -- a conjugate point, and
%   the verdict is FAIL (Astra review 2026-09-11; the earlier rule called
%   every last-bracket crossing ENDPOINT, which was unnecessarily
%   inconclusive). A zero run reaching the last sample, an unresolved last
%   determinant, or a last sample short of t_f stays ENDPOINT: the root may
%   sit at t_f itself, which junction resolution cannot decide.
% • UNRESOLVED SAMPLES: a live sample whose block has sigma_min/sigma_max
%   at or below spec.resolvedTol has an untrustworthy sign. If any such
%   sample exists (including t_f) and no interior root was found, the
%   verdict is UNDETERMINED, never PASS: a same-sign tiny final determinant
%   is not evidence of no root (Astra review #2, 2026-09-11).
% • COVERAGE: without info.Yend a free-time test cannot sample the flow at
%   t_f and stops at t_K. That is reported as UNDETERMINED, never PASS
%   (Astra + Gemini reviews 2026-09-11): the final segment is unmonitored.
% • INITIAL INTERVAL: the sign test starts at the first full-rank junction;
%   the interval before it is UNCOVERED by this instrument (reported as
%   .tFirstFullRank). The block is structurally zero at t = 0, so a sign
%   established at the first junction is not compared against a known
%   short-time sign. conj_spectrum samples inside the first segment.
% • TWO EXACT IDENTITIES are measured as implementation diagnostics: the
%   scaling kernel J(t) p(0) = 0 (.kernelRight) and the left kernel
%   p(t)' J(t) = 0 (.kernelLeft, from symplecticity: the scaling direction
%   flows to (0, p(t)), and its symplectic product with any variation is
%   conserved). Both are relative residuals, max over the live samples.
% • FIXED-tf callers: use stateRows = 1:ny_state (ALL state rows, mass
%   included) -- the interior Jacobi field must vanish in the full state.
%   The mixed block [1:6 14] is the terminal shooting Jacobian and is valid
%   only AT tf (review 2026-09-05, P0.3).
% • This is a NECESSARY-condition check (Jacobi), not a full sufficiency
%   proof: passing means "no disqualifying focal point found at junction
%   resolution".
%
%% Inputs:
%
%  info                     struct                  ms_bvp info with .PHI
%                                                   {1 x K} segment STMs,
%                                                   .Y [ny x K] junction
%                                                   starts, .tGrid, and
%                                                   .Yend [ny x 1] y(tf)
%                                                   (needed to sample the
%                                                   flow column at tf)
%
%  spec                     struct
%   .flow                   fhandle                 f = flow(y): state rows
%                                                   of the dynamics at a
%                                                   junction [m x 1].
%                                                   REQUIRED unless
%                                                   .freeTime is false
%   .stateRows              [1 x m]                 State rows [default 1:6]
%   .costateCols            [1 x mc]                Initial-costate columns
%                                                   [default 8:13]
%   .quotientDir            [mc x 1]                Invariance direction to
%                                                   quotient out [default
%                                                   info.Y(costateCols,1)];
%                                                   pass [] to disable
%   .freeTime               logical                 Append the flow column
%                                                   [default true]
%   .rankTol                double                  sigma_min/sigma_max
%                                                   below which a sample is
%                                                   structurally rank-
%                                                   deficient [1e-13]
%   .zeroTol                double                  |detScaled| at or below
%                                                   which a sample counts as
%                                                   an exact zero [0]
%   .resolvedTol            double                  sigma_min/sigma_max of
%                                                   the LAST sample above
%                                                   which its sign is
%                                                   trusted, so a last-
%                                                   bracket crossing is an
%                                                   interior root [1e-10]
%
%% Outputs:
%
%  out                      struct
%   .t                      [1 x K]                 Junction times sampled
%                                                   (t_2 .. tf)
%   .detScaled              [1 x K]                 sign(det)*|det|^(1/m)
%   .sigRatio               [1 x K]                 sigma_min/sigma_max of
%                                                   the equilibrated block
%   .firstFullRank          double                  First sample index with
%                                                   full rank (K+1: never)
%   .tested                 logical                 At least one full-rank,
%                                                   finite sample existed
%   .sampledThrough         double                  Last sampled time (tf,
%                                                   or t_K without Yend)
%   .nCrossings             double                  All roots found
%                                                   (interior + endpoint)
%   .nInterior              double                  Roots strictly before
%                                                   the last bracket
%   .nTouch                 double                  Isolated zero runs
%                                                   (even-multiplicity
%                                                   touches), included in
%                                                   nInterior
%   .atFinal                logical                 A root on the bracket
%                                                   ending at the last sample
%                                                   that could NOT be placed
%                                                   strictly inside it
%   .nEndResolved           double                  Last-bracket crossings
%                                                   placed strictly inside
%                                                   (counted in nInterior)
%   .nUnresolved            double                  Live samples whose sign
%                                                   is untrustworthy
%                                                   (sigma ratio <= resolvedTol)
%   .covered                logical                 Sampled through t_f
%   .tFirstFullRank         double                  Time of the first full-
%                                                   rank sample (NaN: never);
%                                                   (0, tFirstFullRank) is
%                                                   uncovered
%   .kernelRight            double                  max |J p(0)| / (|J||p(0)|)
%   .kernelLeft             double                  max |p(t)' J| / (|p||J|)
%   .verdict                char                    'PASS' | 'FAIL' (interior
%                                                   root) | 'ENDPOINT' (an
%                                                   exact-zero run reaching
%                                                   the last sample: root at
%                                                   or before t_f, refine) |
%                                                   'UNDETERMINED' (nothing
%                                                   testable, the final
%                                                   segment not covered, or
%                                                   a live sample whose sign
%                                                   is not trustworthy)
%   .reason                 char                    One line on the verdict
%   .pass                   logical                 verdict == 'PASS'
%   .stateRows/.costateCols/.freeTime               Spec echo (provenance)
%
%% Revision History:
%  M. Casey                                                   (c) 08/08/2026
%  Copyright Coorbital Inc.
%% ------------------------ Begin Code Sequence ---------------------------

if nargin < 2, spec = struct(); end
rows = fieldd(spec, 'stateRows', 1:6);
cols = fieldd(spec, 'costateCols', 8:13);
freeT = fieldd(spec, 'freeTime', true);
m = numel(rows);
assert(isfield(info, 'PHI') && ~isempty(info.PHI), ...
    'ms_conjugate_test:noSTMs', ...
    'info.PHI missing: run ms_bvp with opts.keepSTMs = true');

% Column quotient: orthonormal complement of the invariance direction.
if isfield(spec, 'quotientDir')
    qDir = spec.quotientDir;
else
    qDir = info.Y(cols, 1);                    % converged initial costate
end
if isempty(qDir)
    P = eye(numel(cols));
else
    P = null(qDir(:)');                        % [mc x mc-1]
end
nc = size(P, 2) + double(freeT);
assert(nc == m, 'ms_conjugate_test:shape', ...
    'quotiented columns + flow = %d; need %d (= numel(stateRows))', nc, m);
if freeT
    assert(isfield(spec, 'flow') && ~isempty(spec.flow), ...
        'ms_conjugate_test:noFlow', ...
        'freeTime test needs spec.flow (state-rows dynamics at a junction)');
end

rankTol = fieldd(spec, 'rankTol', 1e-13);
resolvedTol = fieldd(spec, 'resolvedTol', 1e-10);

% Chain ALL K segment STMs: samples at t_2 .. t_{K+1} = t_f (review
% 2026-09-05: the old loop stopped at t_K and left the final segment
% unmonitored). For each sample the block is EQUILIBRATED (positive
% row/column scaling, which preserves sign(det)) before the sign test and
% the rank diagnostic sigma_min/sigma_max.
K = numel(info.PHI);
% The flow column at tf needs y(tf) (info.Yend, returned by ms_bvp's
% keepSTMs since 2026-09-05); a free-time caller without it is sampled to
% t_K as before.
nS = K;
if freeT && ~(isfield(info, 'Yend') && ~isempty(info.Yend)), nS = K - 1; end
covered = (nS == K);
PhiCum = eye(size(info.PHI{1}));
dets = zeros(1, nS);  sigR = zeros(1, nS);
kernR = zeros(1, nS);  kernL = zeros(1, nS);
for k = 1:nS
    PhiCum = info.PHI{k} * PhiCum;             % Phi(t_{k+1}, 0)
    J = PhiCum(rows, cols);                    % d state(t) / d costate(0)
    M = J * P;
    if k < K
        yk = info.Y(:,k+1);
    elseif isfield(info, 'Yend') && ~isempty(info.Yend)
        yk = info.Yend;
    else
        yk = [];                               % fixed-time caller without y(tf)
    end
    % the two exact identities of the HOMOGENEOUS (scaling-quotiented) form,
    % as relative residuals; both are skipped when no quotient is in use
    if ~isempty(qDir)
        nJ = max(norm(J), realmin);
        kernR(k) = norm(J*qDir(:)) / (nJ*max(norm(qDir), realmin));
        if ~isempty(yk) && max(rows) + 7 <= numel(yk)
            pk = yk(rows + 7);                 % costates conjugate to the state rows
            kernL(k) = norm(pk(:).'*J) / (nJ*max(norm(pk), realmin));
        end
    end
    if freeT
        M = [M, spec.flow(yk)];
    end
    % Sign from the pivoted LU of the EQUILIBRATED block (permutation
    % parity x diag(U) signs); magnitude as a log-determinant reassembled
    % with the scalings, so an underflowing raw det can never masquerade as
    % an exact zero root (Astra review #2, 2026-09-06).
    [Me, rs, cs] = equilibrate(M);
    if ~all(isfinite(Me(:)))
        dets(k) = NaN;  continue
    end
    sv = svd(Me);
    if sv(1) > 0, sigR(k) = sv(end)/sv(1); end
    [~, U, Pp] = lu(Me);
    du = diag(U);
    if any(du == 0)
        dets(k) = 0;
    else
        sgn = det(Pp) * prod(sign(du));
        logAbs = sum(log(abs(du))) + sum(log(rs)) + sum(log(cs));   % log|det M|
        dets(k) = sgn * exp(logAbs / m);                             % sign * |det|^(1/m)
    end
end

% STRUCTURAL rank deficiency: until the control has acted on the state
% (an initial coast, s == 0), the state block of Phi_xl is identically
% zero -- those samples are not focal points (Astra, review 2026-09-05).
% Samples before the first full-rank one are skipped; crossings are
% counted only from there on.
kFull = find(sigR > rankTol, 1);
if isempty(kFull), kFull = nS + 1; end         % never attained: nothing to test
live = kFull:nS;
tested = ~isempty(live) && ~any(isnan(dets(live)));
tFirstFullRank = NaN;
if kFull <= nS, tFirstFullRank = info.tGrid(kFull + 1); end

% ROOT COUNTING on the live samples (Astra review #2): classify each sample
% as +, - or 0 (exact zero, or |detScaled| <= spec.zeroTol); merge runs of
% zeros; a zero run between opposite signs is ONE root (the crossing), an
% isolated zero run (same sign both sides, or at an end) is one touch/root.
% A root whose bracket ends at the LAST sample (t = tf, or t_K when Yend is
% absent) is an ENDPOINT root: flagged, never counted as interior AND never
% a PASS. At junction resolution the root may sit at tf or strictly before
% it, so the verdict is INCONCLUSIVE -- refine the bracket before classifying
% (Astra doc review 2026-09-07; the earlier 'weak minimum' reading was
% unjustified).
zeroTol = fieldd(spec, 'zeroTol', 0);
nIn = 0;  nEnd = 0;  nTouch = 0;  nEndResolved = 0;
% the last sample's sign is trusted when it is t_f itself and its block is
% resolved; then a sign change on the last bracket is strictly interior
lastResolved = covered && tested && sigR(nS) > resolvedTol && abs(dets(nS)) > zeroTol;
if tested
    dl = dets(live);
    cls = sign(dl);  cls(abs(dl) <= zeroTol) = 0;
    n = numel(cls);
    k = 1;  lastSign = 0;  lastSignIdx = 0;
    while k <= n
        if cls(k) ~= 0
            if lastSign ~= 0 && cls(k) ~= lastSign        % sign change (zeros between merged)
                if k == n && ~lastResolved
                    nEnd = nEnd + 1;
                elseif k == n
                    nIn = nIn + 1;  nEndResolved = nEndResolved + 1;
                else
                    nIn = nIn + 1;
                end
            end
            lastSign = cls(k);  lastSignIdx = k;  k = k + 1;
        else
            k2 = k;  while k2 < n && cls(k2+1) == 0, k2 = k2 + 1; end   % zero run k..k2
            nextSign = 0;  if k2 < n, nextSign = cls(k2+1); end
            if lastSign ~= 0 && nextSign ~= 0 && nextSign ~= lastSign
                % crossing THROUGH the zero run: counted when nextSign is read
            elseif k2 == n
                nEnd = nEnd + 1;                                   % zero run reaching the last sample
            else
                nTouch = nTouch + 1;  nIn = nIn + 1;               % isolated touch: a root
            end
            k = k2 + 1;
        end
    end
end
atFinal = nEnd > 0;
nUnres = 0;
if tested, nUnres = nnz(sigR(live) <= resolvedTol); end
if ~tested
    verdict = 'UNDETERMINED';  reason = 'no full-rank finite sample to test';
elseif nIn > 0
    verdict = 'FAIL';
    reason = sprintf('%d interior root(s)', nIn);
    if nEndResolved > 0
        reason = sprintf('%s, %d placed strictly inside the last bracket (t_f resolved)', reason, nEndResolved);
    end
elseif ~covered
    verdict = 'UNDETERMINED';
    reason = sprintf('final segment (%.6g, %.6g] not covered: info.Yend missing', ...
                     info.tGrid(nS+1), info.tGrid(end));
elseif nUnres > 0
    verdict = 'UNDETERMINED';
    reason = sprintf('%d live sample(s) with sigma ratio <= %.0e: sign not trustworthy', nUnres, resolvedTol);
elseif atFinal
    verdict = 'ENDPOINT';  reason = 'root on the last bracket, not resolvable at junction resolution';
else
    verdict = 'PASS';  reason = 'no sign change on the live samples';
end
kMax = @(x) max([x, 0]);

out = struct('t', info.tGrid(2:nS+1), 'detScaled', dets, 'sigRatio', sigR, ...
             'firstFullRank', kFull, 'tFirstFullRank', tFirstFullRank, 'tested', tested, ...
             'sampledThrough', info.tGrid(nS+1), 'covered', covered, ...
             'nCrossings', nIn + nEnd, 'nInterior', nIn, 'nTouch', nTouch, ...
             'nEndResolved', nEndResolved, 'nUnresolved', nUnres, 'atFinal', atFinal, ...
             'kernelRight', kMax(kernR(live)), 'kernelLeft', kMax(kernL(live)), ...
             'verdict', verdict, 'reason', reason, 'pass', strcmp(verdict, 'PASS'), ...
             'stateRows', rows, 'costateCols', cols, 'freeTime', freeT);
end

% ------------------------------------------------------------------------
function [Me, r, c] = equilibrate(M)
% EQUILIBRATE  Positive row then column scaling by max-abs (sign(det)
% invariant); returns the scalings so log|det M| can be reassembled.
% INPUTS: M [m x m].  OUTPUTS: Me [m x m]; r [m x 1]; c [1 x m].
r = max(abs(M), [], 2);  r(r == 0) = 1;
Me = M ./ r;
c = max(abs(Me), [], 1);  c(c == 0) = 1;
Me = Me ./ c;
end

% ------------------------------------------------------------------------
function v = fieldd(s, f, v0)
% FIELDD  s.(f) if present else v0.  INPUTS: s;f;v0.  OUTPUTS: v.
if isfield(s, f), v = s.(f); else, v = v0; end
end
