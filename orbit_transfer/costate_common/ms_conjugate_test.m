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
% • A sign change on the bracket ending at tf IS counted (it lies in
%   (t_K, tf]) and is additionally flagged .atFinal.
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
%   .nCrossings             double                  Sign changes / zeros
%                                                   from firstFullRank on
%   .atFinal                logical                 The last crossing is on
%                                                   the bracket ending at tf
%   .pass                   logical                 nCrossings == 0
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
PhiCum = eye(size(info.PHI{1}));
dets = zeros(1, nS);  sigR = zeros(1, nS);
for k = 1:nS
    PhiCum = info.PHI{k} * PhiCum;             % Phi(t_{k+1}, 0)
    M = PhiCum(rows, cols) * P;
    if freeT
        if k < K, yk = info.Y(:,k+1); else, yk = info.Yend; end
        M = [M, spec.flow(yk)];
    end
    dk = det(M);
    dets(k) = sign(dk) * abs(dk)^(1/m);
    Me = equilibrate(M);
    sv = svd(Me);
    if all(isfinite(sv)) && sv(1) > 0, sigR(k) = sv(end)/sv(1); end
    if all(isfinite(Me(:))), dets(k) = sign(det(Me)) * abs(dets(k)); end
end

% STRUCTURAL rank deficiency: until the control has acted on the state
% (an initial coast, s == 0), the state block of Phi_xl is identically
% zero -- those samples are not focal points (Astra, review 2026-09-05).
% Samples before the first full-rank one are skipped; crossings are
% counted only from there on.
kFull = find(sigR > rankTol, 1);
if isempty(kFull), kFull = nS + 1; end         % never attained: nothing to test
live = kFull:nS;
dl = dets(live);
idxNZ = find(dl ~= 0);
sgn = sign(dl(idxNZ));
flips = find(diff(sgn) ~= 0);                  % between idxNZ(f), idxNZ(f+1)
% Every bracket is inside (0, t_f]; the one ending at t_f is COUNTED and
% additionally flagged .atFinal (a root exactly at t_f is a measure-zero
% boundary case a consumer may choose to examine; the old rule subtracted
% it, which demoted a strictly interior crossing -- sol/Astra, 2026-09-05).
atFinal = ~isempty(flips) && idxNZ(flips(end)+1) == numel(dl);
nIn = numel(flips) + nnz(dl == 0);

out = struct('t', info.tGrid(2:nS+1), 'detScaled', dets, 'sigRatio', sigR, ...
             'firstFullRank', kFull, ...
             'nCrossings', nIn, 'atFinal', atFinal, 'pass', nIn == 0, ...
             'stateRows', rows, 'costateCols', cols, 'freeTime', freeT);
end

% ------------------------------------------------------------------------
function Me = equilibrate(M)
% EQUILIBRATE  Positive row then column scaling by max-abs (sign(det)
% invariant).  INPUTS: M [m x m].  OUTPUTS: Me [m x m].
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
