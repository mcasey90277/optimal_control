function out = ifs_solve2(Z0, prob, opts)
% IFS_SOLVE2  Scaled, rank-revealing damped Gauss-Newton solve of the IFS system.
%
% Rung-1 replacement for IFS_SOLVE's lsqnonlin. Same square multiple-shooting
% residual (IFS_RESIDUAL) and its analytic complex-step Jacobian, but a custom
% step that (a) two-sided equilibrates the Jacobian (column then row scaling, so
% conditioning is judged on the balanced matrix per the GPT-5.6-sol review), and
% (b) takes a TRUNCATED / rank-revealing SVD Gauss-Newton step: singular
% directions below relTrunc*sigma_max are dropped, so the step moves fully along
% the well-determined directions instead of blowing up along the near-null
% lambda_r0 direction that stalled plain LM (RESULTS.md post-merge finding). A
% Levenberg (Tikhonov) filter with growing lambda is the fallback when the
% truncated step can't achieve a meaningful decrease -- and, crucially, a step
% is only accepted at alpha >= alphaFloor, so the solver damps instead of
% crawling along microscopic steps (the failure mode of the first cut). In
% direct-tau mode (prob.tauParam='direct') the switch times are projected
% monotone (min-gap) after each trial step.
%
% INPUTS:
%   Z0   - seed unknown vector [(8+17k)x1]
%   prob - problem struct (see IFS_RESIDUAL); optional prob.tauParam
%          'sigmoid'(default)|'direct'
%   opts - struct (all optional):
%          tolR       success threshold on ||R||_2         [1e-8]
%          maxIter    Gauss-Newton iterations              [200]
%          relTrunc   STARTING SVD truncation ratio         [1e-2]
%          relTruncFinal  floor the truncation relaxes to   [1e-11]
%          adaptTrunc relax relTrunc when the descent plateaus [true]
%          alphaFloor min line-search step to accept       [1e-2]
%          verbose    print per-iter line                  [true]
%
% TRUNCATION CONTINUATION: aggressive truncation (large relTrunc) enters the
% tiny 40-rev shooting basin fast but FLOORS ||R|| at the residual trapped in
% the dropped near-null subspace. So when the relative decrease plateaus with a
% healthy step, relTrunc is relaxed (x0.1, down to relTruncFinal), picking up
% the next directions and driving ||R|| further down -- an inner analogue of the
% continuation the campaign/papers use.
% OUTPUTS:
%   out - struct: Z, resNorm, iterations, flag, success (resNorm<=tolR),
%         seedResNorm, hist [struct array per iter], rankEq (final), condEq,
%         sigMinEq
%
% REFERENCES:
%   [1] ifs/RESULTS.md (post-merge scaled-SVD diagnosis).
%   [2] reviews/gpt56sol_2026-07-11.md (scale before judging conditioning;
%       rank-revealing Newton step).

if ~isfield(opts,'tolR'),          opts.tolR = 1e-8;       end
if ~isfield(opts,'maxIter'),       opts.maxIter = 200;     end
if ~isfield(opts,'relTrunc'),      opts.relTrunc = 1e-2;   end
if ~isfield(opts,'relTruncFinal'), opts.relTruncFinal = 1e-11; end
if ~isfield(opts,'adaptTrunc'),    opts.adaptTrunc = true; end
if ~isfield(opts,'plateauTol'),    opts.plateauTol = 2e-3; end
if ~isfield(opts,'alphaFloor'),    opts.alphaFloor = 1e-2; end
if ~isfield(opts,'verbose'),       opts.verbose = true;    end
if ~isfield(prob,'tauParam'),      prob.tauParam = 'sigmoid'; end

% ode113 chatters "IntegrationTolNotMet" on stiff perigee arcs; it still returns
% the last accepted state. Silence it here and restore on exit.
ws = warning('off', 'MATLAB:ode113:IntegrationTolNotMet');
cleaner = onCleanup(@() warning(ws));

Z = Z0(:);
k = prob.k;
tauIdx = 8 + 16*k + (1:k);
gapFloor = 1e-7*(prob.tauf - prob.tau0);

Z = projectTau(Z, prob, tauIdx, gapFloor);
[R, J] = ifs_residual(Z, prob);
rn = norm(R);
seedResNorm = rn;
hist = struct('rn',{},'rank',{},'sigMin',{},'sigMax',{},'alpha',{},'lam',{},'stepNorm',{});

if opts.verbose
    fprintf('ifs_solve2: k=%d tauParam=%s relTrunc=%.1e  seed ||R||=%.4e\n', ...
            k, prob.tauParam, opts.relTrunc, rn);
    fprintf('  it        ||R||     rankEq    sigMax    sigMin      cond      alpha     lam\n');
end

flag = 0;                                  % 0: maxIter, 1: converged, -2: stalled
rtc  = opts.relTrunc;                       % current truncation ratio (relaxes over the solve)
for it = 1:opts.maxIter
    if rn <= opts.tolR, flag = 1; break; end

    % --- two-sided equilibration + one SVD at this iterate --------------
    cN = sqrt(sum(J.^2, 1)).';
    dc = 1 ./ max(cN, 1e-300);
    Jc = J .* dc.';
    rNr = sqrt(sum(Jc.^2, 2));
    dr = 1 ./ max(rNr, 1e-300);
    Js = Jc .* dr;
    bs = -(dr .* R);
    [U, Smat, V] = svd(full(Js), 'econ');
    s = diag(Smat);  sMax = s(1);  sMin = s(end);
    tol = rtc * sMax;
    g = U.' * bs;

    % --- escalate damping until an alpha>=alphaFloor step is accepted ---
    lam = 0;  accepted = false;  alphaUsed = 0;  stepNorm = 0;  rnNew = rn;
    for dtry = 1:40
        [dZ, sN] = svdStep(s, g, V, dc, lam, tol);
        alpha = 1;  ok = false;
        for ls = 1:14
            Ztry = projectTau(Z + alpha*dZ, prob, tauIdx, gapFloor);
            rnt  = norm(ifs_residual(Ztry, prob));
            if rnt < (1 - 1e-4*alpha) * rn
                ok = alpha >= opts.alphaFloor;   % reject microsteps -> escalate damping
                break;
            end
            alpha = 0.5*alpha;
        end
        if ok
            accepted = true;  alphaUsed = alpha;  stepNorm = sN;  rnNew = rnt;
            Zacc = Ztry;
            break;
        end
        % no acceptable (large-enough) step: raise Levenberg damping and retry
        if lam == 0, lam = 1e-6*sMax^2; else, lam = lam*10; end
        if lam > 1e10*sMax^2, break; end
    end

    if ~accepted
        % Damping couldn't find a step at this truncation. If we can still
        % relax the truncation (add directions), do so and retry this iterate;
        % only declare a true stall once at the finest truncation.
        if opts.adaptTrunc && rtc > opts.relTruncFinal
            rtc = max(rtc*0.1, opts.relTruncFinal);
            if opts.verbose, fprintf('  %3d  relax relTrunc -> %.1e (stall recovery)\n', it, rtc); end
            continue;
        end
        flag = -2;
        hist(end+1) = mkrec(rn, nnz(s>tol), sMin, sMax, 0, lam, 0); %#ok<AGROW>
        if opts.verbose, fprintf('  %3d  STALLED (no alpha>=%.0e step; lam=%.1e; relTrunc=%.1e)\n', it, opts.alphaFloor, lam, rtc); end
        break;
    end

    prevRn = rn;
    Z = Zacc;  rn = rnNew;
    [R, J] = ifs_residual(Z, prob);
    rankEq = nnz(s > tol);
    hist(end+1) = mkrec(rn, rankEq, sMin, sMax, alphaUsed, lam, stepNorm); %#ok<AGROW>
    if opts.verbose
        fprintf('  %3d  %.4e   %4d/%-4d %.3e  %.3e  %.3e  %.2e  %.1e  rt=%.0e\n', ...
                it, rn, rankEq, numel(s), sMax, sMin, sMax/max(sMin,1e-300), alphaUsed, lam, rtc);
    end
    % Truncation continuation: if the accepted step barely moved ||R||, the
    % well-determined subspace at this truncation is exhausted -- relax to pick
    % up the next directions.
    if opts.adaptTrunc && rtc > opts.relTruncFinal && (prevRn - rn) < opts.plateauTol*prevRn
        rtc = max(rtc*0.1, opts.relTruncFinal);
        if opts.verbose, fprintf('       relax relTrunc -> %.1e (plateau)\n', rtc); end
    end
end

out = struct('Z',Z,'resNorm',rn,'iterations',numel(hist),'flag',flag, ...
             'success', rn <= opts.tolR, 'seedResNorm', seedResNorm, ...
             'hist', hist, 'rankEq', lastfield(hist,'rank'), ...
             'condEq', condFinal(hist), 'sigMinEq', lastfield(hist,'sigMin'));
if opts.verbose
    fprintf('ifs_solve2: DONE k=%d ||R0||=%.3e ||R||=%.3e iters=%d flag=%d success=%d\n', ...
            k, seedResNorm, out.resNorm, out.iterations, flag, out.success);
end
end

% ======================================================================
function [dZ, stepNorm] = svdStep(s, g, V, dc, lam, tol)
if lam <= 0
    f = zeros(size(s));  keep = s > tol;
    f(keep) = g(keep) ./ s(keep);                       % hard truncation
else
    f = (s ./ (s.^2 + lam)) .* g;                       % Tikhonov filter
end
dw = V * f;  dZ = dc .* dw;  stepNorm = norm(dw);
end

% ----------------------------------------------------------------------
function Z = projectTau(Z, prob, tauIdx, gapFloor)
if ~strcmp(prob.tauParam, 'direct') || isempty(tauIdx), return; end
tau = Z(tauIdx);
lo = prob.tau0 + gapFloor;
tau(1) = max(tau(1), lo);
for ii = 2:numel(tau)
    tau(ii) = max(tau(ii), tau(ii-1) + gapFloor);
end
hi = prob.tauf - gapFloor;
if tau(end) > hi
    tau(end) = hi;
    for ii = numel(tau)-1:-1:1
        if tau(ii) > tau(ii+1) - gapFloor, tau(ii) = tau(ii+1) - gapFloor; end
    end
end
Z(tauIdx) = tau;
end

% ----------------------------------------------------------------------
function r = mkrec(rn, rank, sMin, sMax, alpha, lam, stepNorm)
r = struct('rn',rn,'rank',rank,'sigMin',sMin,'sigMax',sMax, ...
           'alpha',alpha,'lam',lam,'stepNorm',stepNorm);
end
function v = lastfield(hist, f)
if isempty(hist), v = NaN; else, v = hist(end).(f); end
end
function v = condFinal(hist)
if isempty(hist), v = NaN; else, v = hist(end).sigMax/max(hist(end).sigMin,1e-300); end
end
