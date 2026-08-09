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
% • Junction resolution: the det is sampled at junctions 2..K, so a
%   conjugate PAIR closer together than one segment can hide, and the
%   final endpoint itself is not sampled. K = 12-48 in the catalogs; treat
%   "no sign change" as strong but not airtight.
% • det is scaled to det^(1/m) magnitude for readable reporting.
% • A sign change on the LAST sampled interval is reported with
%   .atFinal = true; by itself it does not refute local minimality on
%   [0, tf).
% • This is a NECESSARY-condition check (Jacobi), not a full sufficiency
%   proof: passing means "no disqualifying focal point found at junction
%   resolution".
%
%% Inputs:
%
%  info                     struct                  ms_bvp info with .PHI
%                                                   {1 x K} segment STMs,
%                                                   .Y [ny x K] junction
%                                                   starts, .tGrid
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
%
%% Outputs:
%
%  out                      struct
%   .t                      [1 x K-1]               Junction times sampled
%   .detScaled              [1 x K-1]               sign(det)*|det|^(1/m)
%   .nCrossings             double                  Sign changes strictly
%                                                   inside
%   .atFinal                logical                 Sign change on the last
%                                                   sampled interval
%   .pass                   logical                 nCrossings == 0
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

K = numel(info.PHI);
PhiCum = eye(size(info.PHI{1}));
dets = zeros(1, K-1);
for k = 1:K-1
    PhiCum = info.PHI{k} * PhiCum;             % Phi(t_{k+1}, 0)
    M = PhiCum(rows, cols) * P;
    if freeT
        M = [M, spec.flow(info.Y(:,k+1))];
    end
    dk = det(M);
    dets(k) = sign(dk) * abs(dk)^(1/m);
end

% Crossing detection over nonzero samples, on ORIGINAL indices so that
% index compression cannot misclassify atFinal (review finding, GPT
% 2026-08-08). An exactly-zero sample is itself a sampled focal point.
idxNZ = find(dets ~= 0);
sgn = sign(dets(idxNZ));
flips = find(diff(sgn) ~= 0);                  % between idxNZ(f), idxNZ(f+1)
atFinal = ~isempty(flips) && idxNZ(flips(end)+1) == numel(dets);
nIn = numel(flips) - double(atFinal) + nnz(dets == 0);

out = struct('t', info.tGrid(2:K), 'detScaled', dets, ...
             'nCrossings', nIn, 'atFinal', atFinal, 'pass', nIn == 0);
end

% ------------------------------------------------------------------------
function v = fieldd(s, f, v0)
% FIELDD  s.(f) if present else v0.  INPUTS: s;f;v0.  OUTPUTS: v.
if isfield(s, f), v = s.(f); else, v = v0; end
end
