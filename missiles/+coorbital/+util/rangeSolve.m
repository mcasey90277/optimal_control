function [xSol,fAch,info] = rangeSolve(fTarget,fRange,xLo,xHi,tolF)
%% Purpose:
%
%  Solve fRange(x) = fTarget for a scalar parameter x by bisection on the
%  bracket [xLo,xHi], stopping when the ACHIEVED value is within tolF of the
%  target. Written for point-to-point targeting, where x is the
%  thrust-termination time (s) and fRange returns the surface range flown
%  (m), but nothing here knows that: fRange is any handle taking a scalar and
%  returning a scalar, and it is treated as an opaque, expensive black box.
%
%  Bisection rather than a secant or Newton method because each evaluation of
%  fRange is a full trajectory propagation of several seconds, the derivative
%  is unavailable, and the function is only piecewise smooth. Bisection never
%  leaves the bracket, never needs a derivative, and cannot be thrown off a
%  cliff by a local kink -- it just costs a fixed number of evaluations per
%  decade of tolerance, which is the right trade when the alternative is a
%  superlinear method that occasionally has to be restarted.
%
%% Inputs:
%
%  fTarget          [1 x 1]                     Value of fRange to solve for.
%                                               Units are whatever fRange
%                                               returns; metres of surface
%                                               range in the intended use
%
%  fRange           Function handle             y = fRange(x), scalar in,
%                                               finite real scalar out. Called
%                                               once per bracket endpoint and
%                                               once per iteration, never
%                                               twice at the same abscissa
%
%  xLo              [1 x 1]                     Lower bracket endpoint.
%                                               Units are the parameter's;
%                                               seconds in the intended use
%
%  xHi              [1 x 1]                     Upper bracket endpoint, same
%                                               units. Must exceed xLo
%
%  tolF             [1 x 1]                     Convergence tolerance on the
%                                               ACHIEVED value, in fRange's
%                                               units, not on x. Must be
%                                               finite and positive
%
%% Outputs:
%
%  xSol             [1 x 1]                     Solved parameter. On a failure
%                                               to converge, the bracket
%                                               endpoint whose achieved value
%                                               is nearest the target -- the
%                                               best the bracket can do
%
%  fAch             [1 x 1]                     Value fRange actually takes at
%                                               xSol
%
%  info             Struct                      Solve record:
%                                               xLo, xHi    [1 x 1] the
%                                                           original bracket
%                                               fLo, fHi    [1 x 1] achieved
%                                                           values there, in
%                                                           bracket order
%                                               fMin, fMax  [1 x 1] the same
%                                                           two ORDERED low to
%                                                           high: the
%                                                           achievable band
%                                               fTarget     [1 x 1] as passed
%                                               tolF        [1 x 1] as passed
%                                               residual    [1 x 1] fAch minus
%                                                           fTarget, signed
%                                               iterations  [1 x 1] bisection
%                                                           steps taken
%                                               nEval       [1 x 1] calls made
%                                                           to fRange, equal to
%                                                           iterations + 2
%                                               converged   logical
%                                               identifier  char, '' when
%                                                           converged, else a
%                                                           MATLAB error
%                                                           identifier the
%                                                           caller may rethrow
%                                               message     char, always
%                                                           populated, see
%                                                           Notes
%
%% Notes:
%
%  AN UNREACHABLE TARGET IS NOT AN ERROR. If both bracket endpoints fall on
%  the same side of the target -- for a monotonic fRange, exactly the case of
%  a target outside the achievable band -- this function does NOT throw. It
%  returns converged = false, sets info.identifier to
%  coorbital:rangeSolve:targetOutsideBracket, and writes into info.message a
%  sentence naming both bracket endpoints and both achieved values. The
%  refusal is deliberately left to the caller, because only the caller knows
%  how to phrase it: a targeting script wants to say "5100 km requested,
%  3200 to 4700 km reachable, lengthen the burn or pick a nearer target",
%  which needs the band and the units and the problem's vocabulary. A
%  library function that threw here would force every caller into a
%  try/catch merely to read numbers this struct already carries. What must
%  never happen is a silent nearest miss, and converged = false is what
%  prevents it -- a caller that ignores the flag and uses xSol anyway gets an
%  endpoint, not an interior point that looks like a solution.
%
%  A target sitting exactly ON an endpoint value is INSIDE the band, not
%  outside it, and converges immediately at that endpoint with zero
%  iterations. More generally, either endpoint already within tolF of the
%  target short-circuits the solve; if both are, the nearer wins.
%
%  Genuinely malformed input does throw:
%      coorbital:rangeSolve:badBracket      xLo >= xHi, or either non-finite
%      coorbital:rangeSolve:badTolerance    tolF non-finite or not positive
%      coorbital:rangeSolve:nonFiniteEval   fRange returned NaN or Inf
%  The last of these matters in the intended use: a propagation that fails
%  and hands back NaN would otherwise make every subsequent comparison false
%  and let the bracket wander to an arbitrary answer.
%
%  MONOTONICITY IS ASSUMED, NOT CHECKED. Bisection on a bracket whose
%  endpoints straddle the target converges to SOME crossing; with several
%  crossings it finds one of them and does not say which. This is why the
%  design chose thrust-termination time as the parameter -- range is
%  monotonic in it -- rather than the pitch-program loft angle, which has a
%  lofted and a depressed branch either side of a maximum-range hump, or the
%  glide handoff altitude, which bifurcates on phugoid troughs.
%
%  fRange IS NEVER CALLED TWICE AT THE SAME ABSCISSA. The two endpoint values
%  are evaluated once and carried through the loop; only the midpoint is new
%  each iteration. At several seconds a call this is the difference between a
%  usable script and an unusable one, and info.nEval is reported so the
%  caller can see what it paid.
%
%  The loop also stops if the bracket collapses to the floating-point spacing
%  of x before tolF is met, which is what happens if tolF is finer than
%  fRange's own noise floor. That exit reports converged = false with
%  identifier coorbital:rangeSolve:toleranceNotMet rather than looping to the
%  iteration cap.
%
%% References:
%   [1] Press, W.H., et al., "Numerical Recipes," 3rd ed., Cambridge, 2007,
%       Section 9.1. Bisection: guaranteed linear convergence, one bit of the
%       bracket per iteration, and the only bracketing method that cannot
%       fail on a continuous function whose endpoints straddle the root.
%
%% Revision History:
%  Michael Casey                                                08/07/2026
%  Copyright 2026 Coorbital, Inc.
%% ------------------------ Begin Code Sequence ---------------------------

%% Self-demo. A toy range law standing in for a trajectory: monotonic in the
%% cutoff time, saturating as the booster runs out. Both the solve and the
%% refusal are shown, because the refusal is the path a user meets first:
if nargin == 0
             fDemo = @(tc) 7.0e6.*(1 - exp(-tc./60));
    [tcD,rngD,iD] = coorbital.util.rangeSolve(4.0e6,fDemo,10,180,1.0e3);
    fprintf('want 4000.0 km: cutoff %.3f s, achieved %.3f km, %d evals\n', ...
            tcD,rngD./1000,iD.nEval);
       [~,~,iBad] = coorbital.util.rangeSolve(9.0e6,fDemo,10,180,1.0e3);
    fprintf('want 9000.0 km: converged = %d\n',iBad.converged);
    fprintf('  %s\n',iBad.message);
              xSol = [];
              fAch = [];
              info = [];
    return;
end

%% Validate the bracket. An inverted or degenerate interval is a caller
%% mistake, not a physically unreachable request, so it throws:
if ~isscalar(xLo) || ~isscalar(xHi) || ~isfinite(xLo) || ~isfinite(xHi)
    error('coorbital:rangeSolve:badBracket', ...
        'The bracket endpoints must be finite scalars; got %s and %s.', ...
        mat2str(xLo),mat2str(xHi));
end
if xLo >= xHi
    error('coorbital:rangeSolve:badBracket', ...
        ['The bracket must satisfy xLo < xHi; got xLo = %.10g and ', ...
         'xHi = %.10g.'],xLo,xHi);
end

%% Validate the tolerance. Zero or negative can never be met and would spin
%% the solve to its iteration cap for nothing:
if ~isscalar(tolF) || ~isfinite(tolF) || tolF <= 0
    error('coorbital:rangeSolve:badTolerance', ...
        'tolF must be a finite positive scalar; got %s.',mat2str(tolF));
end

%% Evaluate the two endpoints ONCE. These values are carried through the
%% whole solve; nothing below re-evaluates them:
             nEval = 0;
               fLo = evalOnce(fRange,xLo);
               fHi = evalOnce(fRange,xHi);

%% The achievable band, ordered low to high so that the caller does not have
%% to know the sign of the slope to read it:
              fMin = min(fLo,fHi);
              fMax = max(fLo,fHi);

%% Seed the record. Fields are set here in one place so that every exit path
%% below returns the same shape of struct:
          info.xLo = xLo;
          info.xHi = xHi;
          info.fLo = fLo;
          info.fHi = fHi;
         info.fMin = fMin;
         info.fMax = fMax;
      info.fTarget = fTarget;
         info.tolF = tolF;
     info.residual = NaN;
   info.iterations = 0;
        info.nEval = nEval;
    info.converged = false;
   info.identifier = '';
      info.message = '';

%% Either endpoint already within tolerance short-circuits the solve, the
%% nearer of the two winning if both are. This is also the exact-endpoint
%% case: a target sitting precisely on fLo or fHi is ON the band, not outside
%% it, and must converge rather than be refused:
            errLoE = abs(fLo - fTarget);
            errHiE = abs(fHi - fTarget);
if errLoE <= tolF || errHiE <= tolF
    if errLoE <= errHiE
              xSol = xLo;
              fAch = fLo;
    else
              xSol = xHi;
              fAch = fHi;
    end
              info = finishInfo(info,fAch,true,'', ...
                     sprintf(['Converged at the bracket endpoint x = %.10g, ', ...
                              'achieving %.10g against a target of %.10g.'], ...
                              xSol,fAch,fTarget));
    return;
end

%% Both endpoints on the same side of the target. For a monotonic fRange this
%% is precisely a target outside the achievable band. Refuse, but do it by
%% returning: the caller gets the band, the identifier and a ready-made
%% sentence, and decides for itself whether that is fatal:
if sign(fLo - fTarget) == sign(fHi - fTarget)
    if errLoE <= errHiE
              xSol = xLo;
              fAch = fLo;
    else
              xSol = xHi;
              fAch = fHi;
    end
               msg = sprintf(['A target of %.10g is outside the achievable ', ...
                              'band. The bracket x = %.6g and x = %.6g ', ...
                              'achieves %.6g and %.6g respectively, so only ', ...
                              '%.6g to %.6g is reachable. Widen the bracket ', ...
                              'or move the target into that band.'], ...
                              fTarget,xLo,xHi,fLo,fHi,fMin,fMax);
              info = finishInfo(info,fAch,false, ...
                     'coorbital:rangeSolve:targetOutsideBracket',msg);
    return;
end

%% Bisect. aLo and aHi are the live bracket; gLo is the achieved value at aLo,
%% kept only so the half-selection reads as a sign test rather than a
%% re-evaluation. maxIter is a safety net, not the normal exit: the tolerance
%% test or the bracket-collapse test below always fires first in practice,
%% and 200 halvings shrink any representable interval past its own spacing:
               aLo = xLo;
               aHi = xHi;
               gLo = fLo;
           maxIter = 200;
             kIter = 0;
              xSol = NaN;
              fAch = NaN;
         converged = false;
              stop = false;

while kIter < maxIter && ~stop
             kIter = kIter + 1;
              xMid = 0.5.*(aLo + aHi);
              fMid = evalOnce(fRange,xMid);
              xSol = xMid;
              fAch = fMid;

%% Converged on the ACHIEVED value, which is the quantity the user actually
%% cares about -- a kilometre of range, not a millisecond of cutoff time:
    if abs(fMid - fTarget) <= tolF
         converged = true;
              stop = true;

%% Otherwise keep the half whose endpoints still straddle the target. Written
%% as a sign comparison against the cached gLo so it is correct for a
%% decreasing fRange without a separate branch:
    elseif sign(fMid - fTarget) == sign(gLo - fTarget)
               aLo = xMid;
               gLo = fMid;
    else
               aHi = xMid;
    end

%% The bracket has collapsed to the spacing of the numbers themselves. No
%% further halving can change anything, so stop rather than burn evaluations
%% down to maxIter. Reached when tolF is finer than fRange's own noise floor:
    if ~stop && (aHi - aLo) <= 4*eps(max(abs(aLo),abs(aHi)))
              stop = true;
    end
end

%% Report:
   info.iterations = kIter;
if converged
              info = finishInfo(info,fAch,true,'', ...
                     sprintf(['Converged in %d iteration(s) at x = %.10g, ', ...
                              'achieving %.10g against a target of %.10g.'], ...
                              kIter,xSol,fAch,fTarget));
else
               msg = sprintf(['Failed to reach a tolerance of %.6g after %d ', ...
                              'iteration(s). The bracket x = %.6g and ', ...
                              'x = %.6g achieves %.6g and %.6g respectively; ', ...
                              'the closest reached was %.10g at x = %.10g. ', ...
                              'Loosen tolF or check fRange for noise.'], ...
                              tolF,kIter,xLo,xHi,fLo,fHi,fAch,xSol);
              info = finishInfo(info,fAch,false, ...
                     'coorbital:rangeSolve:toleranceNotMet',msg);
end

    function y = evalOnce(fh,x)
    %% Purpose:
    %
    %  Evaluate the caller's handle once, tally the call, and refuse a
    %  non-finite result. Nested so that the tally lives in one place and
    %  info.nEval cannot drift from the number of calls actually made.
    %
    %% Inputs:
    %
    %  fh               Function handle             The caller's fRange
    %
    %  x                [1 x 1]                     Abscissa, in the
    %                                               parameter's units
    %
    %% Outputs:
    %
    %  y                [1 x 1]                     fh(x), in fRange's units
    %
    %% Revision History:
    %  Michael Casey                                                08/07/2026
    %  Copyright 2026 Coorbital, Inc.
    %% ------------------------ Begin Code Sequence ---------------------------
             nEval = nEval + 1;
                 y = fh(x);
        if ~isscalar(y) || ~isfinite(y) || ~isreal(y)
            error('coorbital:rangeSolve:nonFiniteEval', ...
                ['fRange returned %s at x = %.10g. It must return a finite ', ...
                 'real scalar; a propagation that fails must be caught by ', ...
                 'the caller, not passed through as NaN.'],mat2str(y),x);
        end
    end

    function s = finishInfo(s,fOut,ok,ident,msg)
    %% Purpose:
    %
    %  Stamp the exit fields of the record on every return path, so that no
    %  path can forget one and hand the caller a struct missing a field.
    %
    %% Inputs:
    %
    %  s                Struct                      The partially filled record
    %
    %  fOut             [1 x 1]                     Returned achieved value
    %
    %  ok               [1 x 1] logical             Converged flag
    %
    %  ident            Char                        Error identifier, '' if ok
    %
    %  msg              Char                        Human-readable summary
    %
    %% Outputs:
    %
    %  s                Struct                      The completed record
    %
    %% Revision History:
    %  Michael Casey                                                08/07/2026
    %  Copyright 2026 Coorbital, Inc.
    %% ------------------------ Begin Code Sequence ---------------------------
        s.residual = fOut - s.fTarget;
           s.nEval = nEval;
       s.converged = ok;
      s.identifier = ident;
         s.message = msg;
    end
end
