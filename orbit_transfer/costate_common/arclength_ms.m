function A = arclength_ms(resFactory, dRdq, p0, q0, opts)
%% Purpose:
%
%   GENERIC pseudo-arclength continuation of a root curve R(p, q) = 0 in
%   any parameter q, for any residual that comes with its Jacobian. This is
%   the engine behind branch tracing in thrust (arclength_thrust's job) and
%   in arrival phase (the sheet); it owns nothing about orbits.
%
%   Design (GPT-6 Astra reviews of 2026-09-08/09, all applied):
%
%   * SCALED coordinates throughout: p = Dx .* x, q = sq * t. An unscaled
%     arclength mixes O(50) costates with O(1) states and O(3) times.
%   * Tangent = right null vector of the FULL SVD of [R_x R_t]. Solving
%     R_x v = -R_t is exactly singular at a fold, which is where the tangent
%     matters; the null vector is not.
%   * Bordered Newton corrector F = [R; tau'(w - w_pred)] with BACKTRACKING
%     on the scaled residual, a cap on the correction relative to the
%     predictor step (a corrector that lands far from its predictor has
%     jumped branch), an admissibility hook, and convergence tested after
%     the last update too.
%   * Step control by Newton effort (target ~4 iterations), clamped; a
%     successful step never proposes a trial below dsMin.
%   * A sign change of the tangent's q-component is a FOLD CANDIDATE; it is
%     called a fold only when R_x loses rank (sigma_min small) while the
%     augmented [R_x R_t] stays regular.
%   * Every crossing of a requested q LEVEL is located: linear interpolation
%     between the bracketing roots, then Newton at FIXED q. Repeated
%     crossings after folds are all kept -- a grid point can hold several
%     candidates and it is the caller's job to certify and choose.
%   * Every stop is classified and logged. A budget stop keeps everything
%     accepted so far.
%
%% Inputs:
%
%  resFactory               fhandle                 h = resFactory(q), then
%                                                   [R, J] = h(p) with J the
%                                                   n x n Jacobian d R / d p
%
%  dRdq                     fhandle                 col = dRdq(p, q), the
%                                                   n x 1 derivative d R / d q
%                                                   (analytic or FD; the
%                                                   caller's responsibility)
%
%  p0, q0                   [n x 1], double         a ROOT to start from
%
%  opts                     struct (optional)
%   .direction              +1 | -1                 initial sense in q [+1]
%   .Dx                     [n x 1]                 unknown scales -- PASS
%                                                   THEM (per block); the
%                                                   default max(|p0|,
%                                                   0.1 max|p0|) is a
%                                                   placeholder
%   .sq                     double                  parameter scale [max(|q0|,1)]
%   .ds, .dsMin, .dsMax     double                  [0.05, 1e-4, 0.5]
%   .nStep                  int                     [200]
%   .qStop                  [lo hi]                 stop outside [-inf inf]
%   .levels                 [1 x m]                 q levels to record
%                                                   crossings of []
%   .newtonTol, .newtonMax  double, int             [1e-9, 12]
%   .newtonTarget           int                     [4]
%   .maxCorrFrac            double                  reject a corrector that
%                                                   moves more than this
%                                                   times ds from the
%                                                   predictor [2]
%   .foldRatio              double                  sminX/sminAug below this
%                                                   at a sign change = fold
%                                                   [1e-2]
%   .admissible             fhandle                 ok = admissible(p, q) []
%   .deadlineSec            double                  wall budget [inf]
%   .logFile                char                    '' = stdout
%
%% Outputs:
%
%  A                        struct                  .q [1xN] .p {1xN}
%                                                   .tauQ .sminX .sminAug
%                                                   .normR .dsUsed .nNewton
%                                                   .tanResid .folds (struct
%                                                   array: index, q, sminX,
%                                                   sminAug, classified)
%                                                   .crossings (struct array:
%                                                   level, q, p, converged,
%                                                   normR, afterIndex)
%                                                   .stop (reason) .nCalls
%
%% Revision History:
%  M. Casey                                                   (c) 09/09/2026
%  Copyright Coorbital Inc.
%% ------------------------ Begin Code Sequence ---------------------------

if nargin < 5, opts = struct(); end
d = @(f,v) fieldd(opts, f, v);
direction = d('direction', +1);
p0 = p0(:);  n = numel(p0);
% SCALING IS THE CALLER'S JOB. The default below is only a placeholder: a
% component that starts near zero but is O(1) along the curve gets a tiny
% scale and then dominates the arclength (measured: the unit-circle test
% crawled to q = 0.6 in 200 steps with x2(0) = 0 scaled by 1e-2). Pass
% per-BLOCK scales (states, costates, time, extras) for anything real.
Dx = d('Dx', max(abs(p0), 0.1*max(abs(p0))));  Dx = Dx(:);
sq = d('sq', max(abs(q0), 1));
ds = d('ds', 0.05);  dsMin = d('dsMin', 1e-4);  dsMax = d('dsMax', 0.5);
nStep = d('nStep', 200);
qStop = d('qStop', [-inf inf]);
levels = d('levels', []);
nTol = d('newtonTol', 1e-9);  nMax = d('newtonMax', 12);  nTarget = d('newtonTarget', 4);
maxCorrFrac = d('maxCorrFrac', 2);
foldRatio = d('foldRatio', 1e-2);
admissible = d('admissible', []);
deadline = d('deadlineSec', inf);
logFile = d('logFile', '');
lg = @(varargin) logmsg(logFile, sprintf(varargin{:}));
tStart = tic;  nCalls = 0;

pOf = @(x) x .* Dx;  qOf = @(t) t*sq;
x = p0 ./ Dx;  t = q0/sq;

A = struct('q', [], 'p', {{}}, 'tauQ', [], 'sminX', [], 'sminAug', [], ...
           'normR', [], 'dsUsed', [], 'nNewton', [], 'tanResid', [], ...
           'folds', struct('index', {}, 'q', {}, 'sminX', {}, 'sminAug', {}, 'classified', {}), ...
           'crossings', struct('level', {}, 'q', {}, 'p', {}, 'converged', {}, 'normR', {}, 'afterIndex', {}), ...
           'stop', '', 'nCalls', 0);

% the start must be a root
[R, J] = evalRJ(resFactory, qOf(t), pOf(x));  nCalls = nCalls + 1;
if norm(R, inf) > 1e-6
    A.stop = sprintf('start point is not a root (|R| = %.2e)', norm(R, inf));
    lg('arclength_ms: %s', A.stop);  return
end

tau = [];  dsUsed = 0;  nNw = 0;
for step = 0:nStep
    if toc(tStart) > deadline, A.stop = 'deadline'; break, end
    q = qOf(t);  p = pOf(x);
    [R, J] = evalRJ(resFactory, q, p);  nCalls = nCalls + 1;
    Jx = J .* Dx(:)';
    Rt = dRdq(p, q) * sq;
    if ~all(isfinite(Jx(:))) || ~all(isfinite(Rt)), A.stop = 'nonfinite Jacobian'; break, end
    sX = svd(Jx);  sminX = sX(end);
    [~, SA, VA] = svd([Jx Rt]);  sminAug = SA(n, n);
    tauNew = VA(:, end);
    tanResid = norm([Jx Rt]*tauNew) / max(norm([Jx Rt], 'fro'), realmin);
    if isempty(tau)
        if sign(tauNew(end)) ~= direction, tauNew = -tauNew; end
    elseif tauNew'*tau < 0
        tauNew = -tauNew;
    end
    tau = tauNew;

    % record this root
    k = numel(A.q) + 1;
    A.q(k) = q;  A.p{k} = p;  A.tauQ(k) = tau(end);
    A.sminX(k) = sminX;  A.sminAug(k) = sminAug;  A.normR(k) = norm(R, inf);
    A.dsUsed(k) = dsUsed;  A.nNewton(k) = nNw;  A.tanResid(k) = tanResid;
    lg('  step %3d: q = %+.6f  tau_q = %+.3e  sminX = %.2e  sminAug = %.2e  |R| = %.1e  ds = %.2e  nNw = %d', ...
       step, q, tau(end), sminX, sminAug, norm(R, inf), dsUsed, nNw);

    % fold candidate: sign change of the tangent's q-component. LOCALIZE it
    % before classifying: interpolate along the arc to tau_q = 0, correct
    % back onto the curve in the hyperplane through that point, and test
    % the rank structure THERE -- at the nearest accepted root, 0.1 in
    % arclength away, sigma_min(R_x) is still O(ds) and the test is blind
    % (measured on the unit circle: ratio 0.11 at the root, 1e-8 at the fold).
    if k > 1 && sign(A.tauQ(k)) ~= sign(A.tauQ(k-1)) && A.tauQ(k-1) ~= 0
        f = A.tauQ(k-1) / (A.tauQ(k-1) - A.tauQ(k));
        wa = [A.p{k-1} ./ Dx; A.q(k-1)/sq];  wb = [A.p{k} ./ Dx; A.q(k)/sq];
        wf = wa + f*(wb - wa);
        [wf, cvF, nc] = correctInPlane(resFactory, dRdq, wf, tau, Dx, sq, nTol, nMax, admissible);
        nCalls = nCalls + nc;
        qF = qOf(wf(end));  pF = pOf(wf(1:end-1));
        [~, JF_] = evalRJ(resFactory, qF, pF);  nCalls = nCalls + 1;
        sXf = svd(JF_ .* Dx(:)');  sAf = svd([JF_ .* Dx(:)', dRdq(pF, qF)*sq]);
        isFold = cvF && sXf(end) < foldRatio*sAf(end);
        A.folds(end+1) = struct('index', k, 'q', qF, 'sminX', sXf(end), ...
                                'sminAug', sAf(end), 'classified', isFold);
        lg('  *** tangent q-component changed sign; localized at q = %.6f: %s (sminX %.1e / sminAug %.1e, corrector conv %d)', ...
           qF, tern(isFold, 'FOLD', 'sign change, NOT classified as a fold'), sXf(end), sAf(end), cvF);
    end

    % level crossings between the previous root and this one
    if k > 1 && ~isempty(levels)
        qa = A.q(k-1);  qb = A.q(k);
        for L = levels(:)'
            if (L - qa)*(L - qb) < 0 || L == qb
                f = (L - qa)/(qb - qa);
                pL = A.p{k-1} + f*(A.p{k} - A.p{k-1});
                [pL, cv, nr, nc] = newtonFixedQ(resFactory, L, pL, Dx, nTol, nMax);
                nCalls = nCalls + nc;
                A.crossings(end+1) = struct('level', L, 'q', L, 'p', pL, 'converged', cv, ...
                                            'normR', nr, 'afterIndex', k-1);
                lg('    level %.6f crossed: re-solved at fixed q -> converged = %d, |R| = %.1e', L, cv, nr);
            end
        end
    end

    if q < qStop(1) || q > qStop(2), A.stop = 'qStop'; break, end
    if step == nStep, A.stop = 'nStep'; break, end

    % --- predictor / corrector -------------------------------------------
    ok = false;
    while ~ok
        if ds < dsMin, A.stop = 'stalled (ds < dsMin)'; break, end
        w0 = [x; t];  wp = w0 + ds*tau;  w = wp;
        nNw = 0;
        for it = 1:nMax
            qi = qOf(w(end));  pi_ = pOf(w(1:end-1));
            if ~isempty(admissible) && ~admissible(pi_, qi), break, end
            [Ri, Ji] = evalRJ(resFactory, qi, pi_);  nCalls = nCalls + 1;
            if ~all(isfinite(Ri)) || ~all(isfinite(Ji(:))), break, end
            F = [Ri; tau'*(w - wp)];
            if norm(F, inf) < nTol, ok = true; nNw = it - 1; break, end
            Rti = dRdq(pi_, qi) * sq;
            if ~all(isfinite(Rti)), break, end
            JF = [Ji .* Dx(:)', Rti; tau'];
            dw = -(JF \ F);
            if ~all(isfinite(dw)), break, end
            % backtracking on the scaled augmented residual
            alpha = 1;  accepted = false;
            for bt = 1:6
                wt = w + alpha*dw;
                qt = qOf(wt(end));  pt_ = pOf(wt(1:end-1));
                if ~isempty(admissible) && ~admissible(pt_, qt), alpha = alpha/2; continue, end
                Rtst = evalR(resFactory, qt, pt_);  nCalls = nCalls + 1;
                Ft = [Rtst; tau'*(wt - wp)];
                if all(isfinite(Ft)) && norm(Ft, inf) <= (1 - 1e-4*alpha)*norm(F, inf)
                    w = wt;  accepted = true;  break
                end
                alpha = alpha/2;
            end
            if ~accepted, break, end
            if it == nMax
                qi = qOf(w(end));  pi_ = pOf(w(1:end-1));
                Rl = evalR(resFactory, qi, pi_);  nCalls = nCalls + 1;
                if all(isfinite(Rl)) && norm([Rl; tau'*(w - wp)], inf) < nTol
                    ok = true;  nNw = nMax;
                end
            end
        end
        % branch-jump guard: the corrector must land near its predictor
        if ok && norm(w - wp) > maxCorrFrac*ds, ok = false; end
        if ok
            x = w(1:end-1);  t = w(end);  dsUsed = ds;
            fac = sqrt(nTarget/max(nNw, 1));
            ds = min(dsMax, max(dsMin, ds*min(max(fac, 0.5), 2.0)));
        else
            ds = ds/2;
        end
    end
    if ~ok, break, end
end
if isempty(A.stop), A.stop = 'unknown'; end
A.nCalls = nCalls;
lg('arclength_ms: %d roots, q %.6f -> %.6f, %d fold(s), %d level crossing(s), stop = %s, %d residual calls', ...
   numel(A.q), A.q(1), A.q(end), numel(A.folds), numel(A.crossings), A.stop, nCalls);
end

% ------------------------------------------------------------------------
function [R, J] = evalRJ(resFactory, q, p)
% EVALRJ  Residual and Jacobian at (p, q).  INPUTS: resFactory; q; p.
% OUTPUTS: R; J.
h = resFactory(q);  [R, J] = h(p);
end

function R = evalR(resFactory, q, p)
% EVALR  Residual only.  INPUTS: resFactory; q; p.  OUTPUTS: R.
h = resFactory(q);  R = h(p);
end

function [p, converged, normR, nCalls] = newtonFixedQ(resFactory, q, p, Dx, tol, nMax)
% NEWTONFIXEDQ  Plain scaled Newton at a FIXED parameter value, used to land
% exactly on a requested level.  INPUTS: resFactory; q; p (guess); Dx; tol;
% nMax.  OUTPUTS: p; converged; normR; nCalls.
converged = false;  nCalls = 0;  normR = inf;
h = resFactory(q);
for it = 1:nMax
    [R, J] = h(p);  nCalls = nCalls + 1;
    normR = norm(R, inf);
    if ~all(isfinite(R)), return, end
    if normR < tol, converged = true; return, end
    dx = -((J .* Dx(:)') \ R);
    if ~all(isfinite(dx)), return, end
    p = p + dx .* Dx;
end
R = h(p);  nCalls = nCalls + 1;  normR = norm(R, inf);
converged = all(isfinite(R)) && normR < tol;
end

function [w, converged, nCalls] = correctInPlane(resFactory, dRdq, wp, tau, Dx, sq, tol, nMax, admissible)
% CORRECTINPLANE  Bordered Newton from a predictor wp in the hyperplane
% through wp normal to tau -- the same corrector the main loop uses, exposed
% for fold localization.  INPUTS: resFactory; dRdq; wp; tau; Dx; sq; tol;
% nMax; admissible.  OUTPUTS: w; converged; nCalls.
w = wp;  converged = false;  nCalls = 0;
for it = 1:nMax
    q = w(end)*sq;  p = w(1:end-1) .* Dx;
    if ~isempty(admissible) && ~admissible(p, q), return, end
    h = resFactory(q);  [R, J] = h(p);  nCalls = nCalls + 1;
    if ~all(isfinite(R)), return, end
    F = [R; tau'*(w - wp)];
    if norm(F, inf) < tol, converged = true; return, end
    JF = [J .* Dx(:)', dRdq(p, q)*sq; tau'];
    dw = -(JF \ F);
    if ~all(isfinite(dw)), return, end
    w = w + dw;
end
q = w(end)*sq;  p = w(1:end-1) .* Dx;  h = resFactory(q);  R = h(p);  nCalls = nCalls + 1;
converged = all(isfinite(R)) && norm([R; tau'*(w - wp)], inf) < tol;
end

function s = tern(c, a, b)
% TERN  Inline conditional.  INPUTS: c; a; b.  OUTPUTS: s.
if c, s = a; else, s = b; end
end

function v = fieldd(s, f, d_)
% FIELDD  Field with default.  INPUTS: s; f; d_.  OUTPUTS: v.
if isfield(s, f) && ~isempty(s.(f)), v = s.(f); else, v = d_; end
end

function logmsg(f, s)
% LOGMSG  Append to a log file or stdout.  INPUTS: f; s.  OUTPUTS: none.
if isempty(f), fprintf('%s\n', s);
else, fid = fopen(f, 'a'); fprintf(fid, '%s\n', s); fclose(fid);
end
end
