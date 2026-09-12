function R = conj_resolve(Mfun, tGrid, opts)
%% Purpose:
%
%   CANDIDATE DETECTION AND RESOLUTION for a dense conjugate-matrix scan,
%   as a pure function of a matrix-valued function of time, so that it can
%   be driven by synthetic matrices in tests (Astra review #3, 2026-09-11)
%   and by the propagated CR3BP conjugate matrix in conj_spectrum.
%
%   Mfun(t) returns the RAW n x n matrix M(t) whose rank loss is a conjugate
%   time. Each evaluation is column-normalised (rank and determinant sign
%   are unchanged by positive column scalings) and records the spectrum, the
%   determinant sign, the sign's trustworthiness (sigma_min/sigma_max above
%   signTol) and the RAW column-norm ratio: a column whose norm vanishes
%   loses rank without any dip in the normalised spectrum, so it is a
%   candidate in its own right.
%
%   Coarse pass on tGrid (the last sample is t_f and IS included). A sample
%   is a candidate if sigma_n / median < tolCollapse, or it is a LOCAL
%   MINIMUM of the coarse spectrum below tolLocal (a V-shaped zero between
%   two samples reads far above tolCollapse at both), or its raw column
%   ratio < colTol, or a TRUSTED sign differs from its trusted neighbour's.
%   Candidates are clustered by contiguity:
%
%     start     the cluster touching sample 1 -- the structural start-up
%               transient (the matrix grows from zero). Its record is
%               refined like any window; the INCREASING PREFIX from the
%               left edge is structural and is reported as UNCOVERED
%               (.tUncovered, up to the first local maximum); everything
%               after that maximum is resolved like an interior window,
%               so a zero merged into the transient shows as a dip after
%               some growth. A monotone transient stays kind 'start'.
%     endpoint  the cluster reaching sample N = t_f; resolved through t_f,
%               with t_f itself evaluated.
%     interior  anything else.
%
%   RESOLUTION of a cluster keeps EVERY evaluation it makes (nothing is
%   discarded): the window [a-2, b+2] is re-sampled on shifted grids at
%   refine x and refine^2 x, then EVERY local minimum of the assembled
%   evaluations is bracketed by its neighbours and located by golden
%   section, and the bracket's evaluations are checked for unimodality
%   afterwards. The cluster's kind is then decided from the whole record:
%     zero        a bracket of opposite TRUSTED signs (an established root,
%                 .nZeroSign), or the smallest value ever seen at or below
%                 zeroFloor x median (floor-level: possible rank loss,
%                 .nZeroFloor). Both block.
%     unresolved  a vanishing column, a non-unimodal bracket, or a smallest
%                 value inside the floor band [zeroFloor, clearFactor x
%                 zeroFloor). Blocks.
%     near-miss   every bracket unimodal, no bad column, and the smallest
%                 value ever seen at least clearFactor x zeroFloor above the
%                 floor: a positive minimum was located. Cleared.
%   A candidate once seen at the floor is never upgraded by a later search
%   that returns a larger value. .nSmall is the NUMERICAL CORANK at the
%   located minimum (singular values at the floor), not a root multiplicity.
%
%   CLEAR requires: the scan was testable (some sample outside every cluster
%   with a trusted sign -- a scan that never leaves the start transient is
%   not clear), no zero, nothing unresolved, no coarse trusted sign change.
%   The floor is a POLICY value for the matrix's numerical error; the
%   caller is responsible for setting it from a measurement.
%
%% Inputs:
%
%  Mfun                     fhandle                 M = Mfun(t), [n x n],
%                                                   t in (0, t_f]
%  tGrid                    [1 x N]                 uniform coarse times,
%                                                   tGrid(N) = t_f, > 0
%  opts                     struct (optional)
%   .tolCollapse [1e-3] absolute dip, .tolLocal [0.2] a LOCAL MINIMUM of
%   the coarse spectrum below this (relative to the median) is a candidate
%   whatever its depth -- a V-shaped zero between two coarse samples can
%   read far above tolCollapse at both of them (Astra review #3 synthetic
%   cases), .refine [4] .zeroFloor [1e-7] .clearFactor [100]
%   .colTol [1e-10] raw column-norm ratio below which a column is
%   vanishing, .signTol [1e-10] sigma_min/sigma_max above which a
%   determinant sign is trusted, .maxMinima [8] local minima resolved per
%   cluster, .goldenIter [40]
%
%% Outputs:
%
%  R                        struct                  .t .sv [n x N] .det
%                                                   .sigRatio .colRatio
%                                                   .trusted (coarse signs,
%                                                   0 = untrusted) .med
%                                                   .minRel .nInterior
%                                                   (coarse trusted sign
%                                                   changes) .tFirst .tf
%                                                   .svEnd .candidates
%                                                   (struct array) .nStart
%                                                   .tUncovered .nInteriorCand
%                                                   .nEndCand .nZero
%                                                   .nZeroSign .nZeroFloor
%                                                   .nNearMiss .nUnresolved
%                                                   .multiplicity .testable
%                                                   .clear .reason .nEval
%                                                   .zeroFloor .clearFactor
%
%% Revision History:
%  M. Casey                                                   (c) 09/11/2026
%  Copyright Coorbital Inc.
%% ------------------------ Begin Code Sequence ---------------------------

if nargin < 3, opts = struct(); end
d = @(f,v) fieldd(opts, f, v);
tolCollapse = d('tolCollapse', 1e-3);  tolLocal = d('tolLocal', 0.2);
refine = max(2, round(d('refine', 4)));
zeroFloor = d('zeroFloor', 1e-7);  clearFactor = d('clearFactor', 100);
colTol = d('colTol', 1e-10);  signTol = d('signTol', 1e-10);
maxMinima = d('maxMinima', 8);  goldenIter = d('goldenIter', 40);
assert(zeroFloor > 0 && clearFactor > 1 && tolCollapse > 0, 'conj_resolve:opts', 'floor, clear factor and collapse tolerance must be positive (clearFactor > 1)');

tGrid = tGrid(:).';  N = numel(tGrid);
assert(N >= 3 && all(diff(tGrid) > 0) && tGrid(1) > 0, 'conj_resolve:grid', 'tGrid must be increasing, positive, with >= 3 samples');
dt = tGrid(2) - tGrid(1);
assert(max(abs(diff(tGrid) - dt)) < 1e-9*dt, 'conj_resolve:grid', 'tGrid must be uniform');
tf = tGrid(end);

% ---- evaluation, with a persistent record ---------------------------------
nEval = 0;
    function e = evalAt(t)
        M = Mfun(t);
        assert(isnumeric(M) && ismatrix(M) && size(M,1) == size(M,2) && all(isfinite(M(:))), ...
               'conj_resolve:matrix', 'Mfun(%.6g) did not return a finite square matrix', t);
        nrm = vecnorm(M, 2, 1);
        colRatio = min(nrm) / max(max(nrm), realmin);
        Mn = M ./ max(nrm, realmin);
        sv = svd(Mn);
        sig = sv(end) / max(sv(1), realmin);
        dt_ = det(Mn);
        e = struct('t', t, 'sv', sv(:), 'det', dt_, 'sig', sig, 'colRatio', colRatio);
        nEval = nEval + 1;
    end

% ---- coarse pass ------------------------------------------------------------
n = [];
Ec = repmat(struct('t', [], 'sv', [], 'det', [], 'sig', [], 'colRatio', []), 1, N);
for k = 1:N
    Ec(k) = evalAt(tGrid(k));
    if isempty(n), n = numel(Ec(k).sv); end
end
sv6 = arrayfun(@(e) e.sv(end), Ec);
med = max(median(sv6(1:N-1)), realmin);
rel = sv6 / med;
trusted = arrayfun(@(e) sign(e.det) * (e.sig > signTol), Ec);
colR = [Ec.colRatio];

R.t = tGrid;  R.sv = [Ec.sv];  R.det = [Ec.det];  R.sigRatio = [Ec.sig];
R.colRatio = colR;  R.trusted = trusted;  R.med = med;  R.minRel = min(rel);
R.tf = tf;  R.svEnd = Ec(N).sv;
% coarse TRUSTED sign changes strictly inside (compat with the old field)
tr = trusted(1:N-1);
last = 0;  nInt = 0;  tFirst = NaN;
for k = 1:N-1
    if tr(k) ~= 0
        if last ~= 0 && tr(k) ~= last
            nInt = nInt + 1;  if isnan(tFirst), tFirst = tGrid(k); end
        end
        last = tr(k);
    end
end
R.nInterior = nInt;  R.tFirst = tFirst;

% ---- candidates and clusters ---------------------------------------------------
isCand = rel < tolCollapse | colR < colTol;
for k = 2:N-1                                  % coarse local minima, any depth below tolLocal
    if rel(k) < tolLocal && rel(k) < rel(k-1) && rel(k) <= rel(k+1), isCand(k-1:k+1) = true; end
end
last = 0;  lastK = 0;
for k = 1:N
    if trusted(k) ~= 0
        if last ~= 0 && trusted(k) ~= last, isCand(lastK:k) = true; end
        last = trusted(k);  lastK = k;
    end
end
edges = diff([false, isCand, false]);
starts = find(edges == 1);  ends = find(edges == -1) - 1;
testable = any(~isCand & trusted ~= 0);

cands = emptyCand();
nStart = 0;  tUncovered = 0;
for e = 1:numel(starts)
    a = starts(e);  b = ends(e);
    if a == 1
        % the start-up transient GROWS from zero: its record is increasing
        % from the left edge, and that increasing prefix is structural, not
        % evidence. Everything from its first local MAXIMUM onward is
        % resolved like any other window (a zero merged into the transient
        % shows as a dip after some growth); the prefix is UNCOVERED.
        nStart = nStart + 1;
        c = resolveWindow(1, min(b+2, N), 'start', 1);
        tUncovered = c.tUncovered;
        cands(end+1) = c; %#ok<AGROW>
        continue
    end
    cls = 'interior';  if b == N, cls = 'endpoint'; end
    [~, im] = min(rel(a:b));  km = a + im - 1;
    cands(end+1) = resolveWindow(max(a-2, 1), min(b+2, N), cls, km); %#ok<AGROW>
end

% ---- verdict ------------------------------------------------------------------
kinds = {cands.kind};  classes = {cands.class};
refined = ~strcmp(kinds, 'start');            % a start cluster resolved past its growth counts
R.candidates = cands;
R.nStart = nStart;  R.tUncovered = tUncovered;
R.nInteriorCand = nnz(strcmp(classes, 'interior') | (strcmp(classes, 'start') & refined));
R.nEndCand = nnz(strcmp(classes, 'endpoint'));
R.nZeroSign = nnz(refined & strcmp(kinds, 'zero') & [cands.nBrackets] > 0);
R.nZero = nnz(refined & strcmp(kinds, 'zero'));
R.nZeroFloor = R.nZero - R.nZeroSign;
R.nNearMiss = nnz(refined & strcmp(kinds, 'near-miss'));
R.nUnresolved = nnz(refined & strcmp(kinds, 'unresolved'));
R.multiplicity = nnz(refined & strcmp(kinds, 'zero') & [cands.nSmall] >= 2);
R.testable = testable;
R.clear = testable && R.nZero == 0 && R.nUnresolved == 0 && R.nInterior == 0;
if ~testable
    R.reason = 'not testable: no sample outside every candidate cluster with a trusted sign';
elseif R.nZero > 0
    R.reason = sprintf('%d zero(s): %d by trusted sign bracket, %d at the floor', R.nZero, R.nZeroSign, R.nZeroFloor);
elseif R.nUnresolved > 0
    R.reason = sprintf('%d unresolved candidate(s)', R.nUnresolved);
elseif R.nInterior > 0
    R.reason = sprintf('%d coarse trusted sign change(s)', R.nInterior);
else
    R.reason = sprintf('clear: %d near-miss(es) with a located positive minimum, uncovered [0, %.4g]', R.nNearMiss, tUncovered);
end
R.nEval = nEval;  R.zeroFloor = zeroFloor;  R.clearFactor = clearFactor;

% ---------------------------------------------------------------------------
    function c = resolveWindow(ia, ib, cls, km)
        % every evaluation in the window is KEPT: coarse samples, two shifted
        % grids, both golden-section searches
        c = blankCand();
        c.class = cls;  c.kMin = km;  c.tOverTf = tGrid(km)/tf;  c.rel = rel(km);
        c.ratio = Ec(km).sv(end)/max(Ec(km).sv(end-1), realmin);
        tA = tGrid(ia);  tB = tGrid(ib);
        c.tWindow = [tA, tB];
        E = Ec(ia:ib);
        coarseMin = min(arrayfun(@(x) x.sv(end), E));
        h1 = dt/refine;  h2 = dt/refine^2;
        t1 = (tA + h1/2):h1:(tB - h1/4);           % shifted: no coarse node is re-sampled
        t2 = (tA + h2/2):h2:(tB - h2/4);
        for t = t1, E(end+1) = evalAt(t); end %#ok<AGROW>
        m1 = min(arrayfun(@(x) x.sv(end), E(numel(Ec(ia:ib))+1:end)));
        n1 = numel(E);
        for t = t2, E(end+1) = evalAt(t); end %#ok<AGROW>
        m2 = min(arrayfun(@(x) x.sv(end), E(n1+1:end)));
        c.refineRatio  = m1/max(coarseMin, realmin);
        c.refineRatio2 = m2/max(m1, realmin);
        % sort, then locate EVERY local minimum by golden section on its
        % neighbour bracket, checking unimodality of what was evaluated
        [~, order] = sort([E.t]);  E = E(order);
        v = arrayfun(@(x) x.sv(end), E) / med;
        c.tUncovered = 0;
        if strcmp(cls, 'start')
            % strip the increasing prefix: uncovered up to the first local
            % maximum; if the whole record is non-decreasing there is nothing
            % to resolve and the cluster is the transient, full stop
            iMax = find(diff(v) < 0, 1);
            if isempty(iMax)
                c.tUncovered = E(end).t;  c.kind = 'start';
                [c.sigMinRel, iMin] = min(v);  c.tMinOverTf = E(iMin).t/tf;  c.svAtMin = E(iMin).sv;
                c.nSmall = NaN;  c.nEval = numel(E);
                return
            end
            c.tUncovered = E(iMax).t;
            E = E(iMax:end);  v = v(iMax:end);
        end
        locs = find(v(2:end-1) < v(1:end-2) & v(2:end-1) <= v(3:end)) + 1;
        if isempty(locs), [~, locs] = min(v); locs = locs(locs > 1 & locs < numel(v)); end
        [~, ord] = sort(v(locs));  locs = locs(ord(1:min(end, maxMinima)));
        allUnimodal = true;
        brackets = zeros(numel(locs), 2);
        for q = 1:numel(locs)
            L = locs(q);
            brackets(q, :) = [E(L-1).t, E(L+1).t];
            E = [E, golden(E(L-1).t, E(L+1).t)]; %#ok<AGROW>
        end
        [~, order] = sort([E.t]);  E = E(order);
        v = arrayfun(@(x) x.sv(end), E) / med;
        % unimodality of EVERYTHING evaluated inside each bracket, endpoints
        % included; values below the floor are noise and are not held to it
        for q = 1:size(brackets, 1)
            in = [E.t] >= brackets(q,1) - 1e-12 & [E.t] <= brackets(q,2) + 1e-12;
            w = v(in);  [~, iw] = min(w);
            uni = all(diff(w(1:iw)) <= zeroFloor) && all(diff(w(iw:end)) >= -zeroFloor);
            allUnimodal = allUnimodal && uni;
        end
        % trusted-sign brackets across the WHOLE record
        tsg = arrayfun(@(x) sign(x.det) * (x.sig > signTol), E);
        nB = 0;  lastS = 0;
        for q = 1:numel(tsg)
            if tsg(q) ~= 0
                if lastS ~= 0 && tsg(q) ~= lastS, nB = nB + 1; end
                lastS = tsg(q);
            end
        end
        c.nBrackets = nB;  c.signChange = nB > 0;
        [c.sigMinRel, iMin] = min(v);
        c.tMinOverTf = E(iMin).t / tf;
        c.svAtMin = E(iMin).sv;
        c.nSmall = nnz(E(iMin).sv / med <= zeroFloor);
        c.colBad = any([E.colRatio] < colTol);
        c.unimodal = allUnimodal;
        c.nEval = numel(E);
        if strcmp(cls, 'start') && numel(locs) == 0 && nB == 0
            c.kind = 'start';  return                % a dip that never became a local minimum
        end
        if nB > 0 || c.sigMinRel <= zeroFloor
            c.kind = 'zero';
        elseif c.colBad || ~allUnimodal
            c.kind = 'unresolved';
        elseif c.sigMinRel >= clearFactor*zeroFloor
            c.kind = 'near-miss';
        else
            c.kind = 'unresolved';
        end
    end

    function Eb = golden(a, b)
        % golden-section minimisation of sigma_n(t) on [a, b]; returns every
        % evaluation made (the caller judges unimodality on the full record)
        gr = (sqrt(5) - 1)/2;
        x1 = b - gr*(b - a);  x2 = a + gr*(b - a);
        e1 = evalAt(x1);  e2 = evalAt(x2);  Eb = [e1, e2];
        f1 = e1.sv(end);  f2 = e2.sv(end);
        for it = 1:goldenIter
            if f1 < f2, b = x2;  x2 = x1;  f2 = f1;  x1 = b - gr*(b - a);  e1 = evalAt(x1);  f1 = e1.sv(end);  Eb(end+1) = e1; %#ok<AGROW>
            else,       a = x1;  x1 = x2;  f1 = f2;  x2 = a + gr*(b - a);  e2 = evalAt(x2);  f2 = e2.sv(end);  Eb(end+1) = e2; %#ok<AGROW>
            end
            if (b - a) < 1e-9*max(abs(b), 1), break, end
        end
    end
end

function c = blankCand()
% BLANKCAND  One candidate record with every field present.
% INPUTS: none.  OUTPUTS: c.
c = struct('tOverTf', NaN, 'rel', NaN, 'ratio', NaN, 'class', '', 'kind', 'n/a', ...
           'refineRatio', NaN, 'refineRatio2', NaN, 'signChange', false, 'kMin', NaN, ...
           'tMinOverTf', NaN, 'sigMinRel', NaN, 'nSmall', NaN, 'svAtMin', [], ...
           'nBrackets', 0, 'colBad', false, 'unimodal', true, 'tWindow', [NaN NaN], 'nEval', 0, ...
           'tUncovered', 0);
end

function c = emptyCand()
% EMPTYCAND  An empty candidate array with the schema of blankCand.
% INPUTS: none.  OUTPUTS: c [1 x 0].
c = repmat(blankCand(), 1, 0);
end

function v = fieldd(s, f, v0)
% FIELDD  s.(f) if present else v0.  INPUTS: s; f; v0.  OUTPUTS: v.
if isfield(s, f) && ~isempty(s.(f)), v = s.(f); else, v = v0; end
end
