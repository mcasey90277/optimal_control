function test_rangeSolve()
%% Purpose:
%
%  Verify the one-dimensional bisection solver against SYNTHETIC monotonic
%  functions whose roots are known in closed form, so that the solver is
%  exercised entirely independently of any trajectory propagation. Nothing
%  here integrates an equation of motion; the whole point is that a solver
%  bug and a dynamics bug must not be able to hide behind one another.
%
%  Two reference functions are used:
%
%      fCube(x) = x^3 + 2   on [0,3], strictly INCREASING, band [2,29]
%      fLin(x)  = 100 - 10x on [0,8], strictly DECREASING, band [20,100]
%
%  fCube(2) = 10 and fLin(6.5) = 35, both exactly, so the analytic roots are
%  x = 2 and x = 6.5 with no rounding to allow for. The decreasing case is
%  carried because a solver that assumes an increasing function updates the
%  wrong half of the bracket and diverges on it, and thrust-termination time
%  is not guaranteed to be the increasing sense in every future use.
%
%  Several assertions are written against a HAND-TRACED bisection sequence
%  rather than a tolerance. On [0,3] with target 10 the midpoints are
%  1.5, 2.25, 1.875, 2.0625 and the values there are 5.375, 13.390625,
%  8.591796875, 10.773681640625. Every one of those is a dyadic rational
%  and therefore exact in binary double precision, so they are asserted with
%  ==, not with a tolerance. The trace is worked out below where it is used.
%
%% Inputs:
%
%  none
%
%% Outputs:
%
%  none                                         Throws on any failed assertion
%
%% Revision History:
%  Michael Casey                                                08/07/2026
%  Copyright 2026 Coorbital, Inc.
%% ------------------------ Begin Code Sequence ---------------------------

%% The counted reference function shares this workspace through the nested
%% function at the bottom of the file. Preallocated, never grown; nCalled is
%% the fill level and the call count the TEST keeps for itself, independently
%% of whatever the solver reports in info.nEval:
             evalX = zeros(1,512);
           nCalled = 0;

%% The two synthetic references, uncounted:
             fCube = @(x) x.^3 + 2;
              fLin = @(x) 100 - 10.*x;

%% -------------------------------------------------------------------
%% Convergence on an INCREASING function to an analytic root
%% -------------------------------------------------------------------
%% fCube(2) = 8 + 2 = 10 exactly, so the root of fCube(x) = 10 is x = 2. The
%% derivative there is 3*2^2 = 12, so a residual of 1e-9 in f corresponds to
%% about 8e-11 in x; asserting |x-2| < 1e-9 is therefore comfortable but not
%% vacuous:
    [xSol,fAch,info] = coorbital.util.rangeSolve(10,fCube,0,3,1e-9);
    assert(info.converged,'increasing case failed to converge');
    assert(abs(fAch - 10) <= 1e-9, ...
        'achieved value outside the requested tolerance: %.17g',fAch);
    assert(abs(xSol - 2) < 1e-9, ...
        'solved parameter wrong: got %.17g, analytic root is 2',xSol);

%% fAch must be the value the CALLER'S function takes at xSol, not a
%% separately bookkept number that could have drifted from it:
    assert(fAch == fCube(xSol), ...
        'fAch %.17g is not fCube(xSol) = %.17g',fAch,fCube(xSol));

%% -------------------------------------------------------------------
%% Convergence on a DECREASING function
%% -------------------------------------------------------------------
%% fLin(6.5) = 100 - 65 = 35 exactly. Slope -10, so a 1e-9 residual is 1e-10
%% in x:
    [xSol,fAch,info] = coorbital.util.rangeSolve(35,fLin,0,8,1e-9);
    assert(info.converged,'decreasing case failed to converge');
    assert(abs(xSol - 6.5) < 1e-9, ...
        'decreasing case solved parameter wrong: got %.17g, analytic root is 6.5',xSol);
    assert(abs(fAch - 35) <= 1e-9, ...
        'decreasing case achieved value outside tolerance: %.17g',fAch);

%% -------------------------------------------------------------------
%% The hand-traced bisection, asserted exactly
%% -------------------------------------------------------------------
%% Target 10, bracket [0,3], tolerance 1. Pure bisection, stopping the moment
%% the midpoint value is within the tolerance of the target:
%%
%%   endpoints  fCube(0) = 2       below 10
%%              fCube(3) = 29      above 10        -> root is bracketed
%%   iter 1     m = 1.5            f = 5.375       below 10 -> [1.5 , 3   ]
%%   iter 2     m = 2.25           f = 13.390625   above 10 -> [1.5 , 2.25]
%%   iter 3     m = 1.875          f = 8.591796875 below 10 -> [1.875, 2.25]
%%   iter 4     m = 2.0625         f = 10.773681640625  (2.0625 = 33/16, and
%%                                 33^3/16^3 = 35937/4096 = 8.773681640625)
%%                                 |f - 10| = 0.773681640625 <= 1  -> stop
%%
%% Six evaluations in all, four of them iterations. Every midpoint is a
%% dyadic rational and every cube of one is exact in binary, so these are
%% equalities and not tolerances. This is the assertion that a solver
%% returning the bracket midpoint without checking the tolerance fails:
    [xSol,fAch,info] = coorbital.util.rangeSolve(10,fCube,0,3,1);
    assert(xSol == 2.0625, ...
        'hand-traced bisection expected x = 2.0625, got %.17g',xSol);
    assert(fAch == 10.773681640625, ...
        'hand-traced bisection expected f = 10.773681640625, got %.17g',fAch);
    assert(info.iterations == 4, ...
        'hand-traced bisection expected 4 iterations, got %d',info.iterations);
    assert(info.nEval == 6, ...
        'hand-traced bisection expected 6 evaluations, got %d',info.nEval);
    assert(info.converged,'the hand-traced case must report converged');

%% -------------------------------------------------------------------
%% A tighter tolerance gives a tighter result
%% -------------------------------------------------------------------
%% From the trace above, tolF = 1 stops with a residual of 0.773681640625.
%% tolF = 1e-9 must do materially better and must take more iterations to do
%% it. Both bounds are asserted as well as the comparison, so that a solver
%% which simply ignores tolF and always runs to machine precision would still
%% be caught by the iteration-count assertion in the hand trace above:
  [~,fLoose,infoLoose] = coorbital.util.rangeSolve(10,fCube,0,3,1);
  [~,fTight,infoTight] = coorbital.util.rangeSolve(10,fCube,0,3,1e-9);
          errLoose = abs(fLoose - 10);
          errTight = abs(fTight - 10);
    assert(errLoose <= 1, ...
        'loose solve violated its own tolerance: residual %.17g',errLoose);
    assert(errTight <= 1e-9, ...
        'tight solve violated its own tolerance: residual %.17g',errTight);
    assert(errTight < errLoose, ...
        'a tighter tolerance did not give a tighter result: %.17g vs %.17g', ...
        errTight,errLoose);
    assert(infoTight.iterations > infoLoose.iterations, ...
        'a tighter tolerance did not cost more iterations: %d vs %d', ...
        infoTight.iterations,infoLoose.iterations);

%% -------------------------------------------------------------------
%% Target ABOVE the reachable band: refuse, do not throw
%% -------------------------------------------------------------------
%% The band on [0,3] is [2,29]. A request for 1000 cannot be met. The
%% contract is that the solver returns rather than throws, so that the CALLER
%% can print a refusal naming what is achievable and what to change -- a
%% targeting script that silently returned the nearest miss would be worse
%% than one that refused, and the caller needs the band to say anything
%% useful. The call below is deliberately NOT wrapped in try/catch: if it
%% throws, the test fails on the throw itself, which is the point:
    [xSol,fAch,info] = coorbital.util.rangeSolve(1000,fCube,0,3,1e-9);
    assert(~info.converged, ...
        'an unreachable target was reported as converged');
    assert(info.fMin == 2 && info.fMax == 29, ...
        'the achievable band must be [2,29], got [%.17g,%.17g]',info.fMin,info.fMax);
    assert(info.xLo == 0 && info.xHi == 3, ...
        'info must carry the original bracket, got [%.17g,%.17g]',info.xLo,info.xHi);
    assert(info.fLo == 2 && info.fHi == 29, ...
        'info must carry the endpoint values, got %.17g and %.17g',info.fLo,info.fHi);

%% The best it could do is the top of the band, and it must say so rather
%% than hand back something in the interior:
    assert(xSol == 3 && fAch == 29, ...
        'an over-range request must return the nearest reachable endpoint, got x=%.17g f=%.17g', ...
        xSol,fAch);

%% A named identifier the caller can rethrow, and a message that states BOTH
%% bracket values and BOTH achieved values so a user sees the band at once:
    assert(strcmp(info.identifier,'coorbital:rangeSolve:targetOutsideBracket'), ...
        'wrong refusal identifier: %s',info.identifier);
    assert(contains(info.message,sprintf('%.6g',0)), ...
        'the refusal message does not state the lower bracket value: %s',info.message);
    assert(contains(info.message,sprintf('%.6g',3)), ...
        'the refusal message does not state the upper bracket value: %s',info.message);
    assert(contains(info.message,sprintf('%.6g',2)), ...
        'the refusal message does not state the lower achieved value: %s',info.message);
    assert(contains(info.message,sprintf('%.6g',29)), ...
        'the refusal message does not state the upper achieved value: %s',info.message);

%% Only the two endpoints were worth evaluating; there is nothing to bisect:
    assert(info.nEval == 2, ...
        'an out-of-band request should cost 2 evaluations, cost %d',info.nEval);
    assert(info.iterations == 0, ...
        'an out-of-band request should take 0 iterations, took %d',info.iterations);

%% -------------------------------------------------------------------
%% Target BELOW the reachable band, and on a decreasing function
%% -------------------------------------------------------------------
    [xSol,fAch,info] = coorbital.util.rangeSolve(-5,fCube,0,3,1e-9);
    assert(~info.converged,'an under-range target was reported as converged');
    assert(xSol == 0 && fAch == 2, ...
        'an under-range request must return the nearest reachable endpoint, got x=%.17g f=%.17g', ...
        xSol,fAch);
    assert(info.fMin == 2 && info.fMax == 29, ...
        'band wrong on the under-range case: [%.17g,%.17g]',info.fMin,info.fMax);

%% On fLin over [0,8] the band is [20,100] even though fLo = 100 exceeds
%% fHi = 20. fMin and fMax must be ORDERED, not merely copied across:
    [~,~,info] = coorbital.util.rangeSolve(500,fLin,0,8,1e-9);
    assert(~info.converged,'unreachable decreasing target reported as converged');
    assert(info.fMin == 20 && info.fMax == 100, ...
        'the band must be ordered low-to-high regardless of slope sign: [%.17g,%.17g]', ...
        info.fMin,info.fMax);
    assert(info.fLo == 100 && info.fHi == 20, ...
        'fLo and fHi must stay attached to xLo and xHi: %.17g, %.17g', ...
        info.fLo,info.fHi);

%% -------------------------------------------------------------------
%% Exact-endpoint targets
%% -------------------------------------------------------------------
%% fCube(0) = 2 exactly. The answer is the endpoint itself, reached with no
%% bisection at all, and it must be reported as converged rather than as an
%% out-of-band refusal -- the target is ON the boundary of the band, not
%% outside it:
    [xSol,fAch,info] = coorbital.util.rangeSolve(2,fCube,0,3,1e-9);
    assert(info.converged,'a target sitting exactly on fLo must converge');
    assert(xSol == 0 && fAch == 2, ...
        'exact lower endpoint wrong: x=%.17g f=%.17g',xSol,fAch);
    assert(info.iterations == 0, ...
        'an exact endpoint needs no iterations, took %d',info.iterations);
    assert(info.nEval == 2, ...
        'an exact endpoint costs the two endpoint evaluations, cost %d',info.nEval);

%% fCube(3) = 29 exactly, the upper end:
    [xSol,fAch,info] = coorbital.util.rangeSolve(29,fCube,0,3,1e-9);
    assert(info.converged,'a target sitting exactly on fHi must converge');
    assert(xSol == 3 && fAch == 29, ...
        'exact upper endpoint wrong: x=%.17g f=%.17g',xSol,fAch);
    assert(info.iterations == 0, ...
        'an exact endpoint needs no iterations, took %d',info.iterations);

%% And on the decreasing function, where the exact endpoint is the LOW one:
    [xSol,~,info] = coorbital.util.rangeSolve(20,fLin,0,8,1e-9);
    assert(info.converged && xSol == 8, ...
        'exact endpoint on the decreasing function wrong: x=%.17g',xSol);

%% -------------------------------------------------------------------
%% An inverted or degenerate bracket is an error, not a refusal
%% -------------------------------------------------------------------
%% This is a caller mistake rather than a physically unreachable request, so
%% it throws, and it throws with an identifier a caller can catch by name:
            threwE = false;
    try
        coorbital.util.rangeSolve(10,fCube,3,0,1e-9);
    catch err
            threwE = true;
        assert(strcmp(err.identifier,'coorbital:rangeSolve:badBracket'), ...
            'inverted bracket threw the wrong identifier: %s',err.identifier);
        assert(contains(err.message,'3') && contains(err.message,'0'), ...
            'the bad-bracket message must state both bracket values: %s',err.message);
    end
    assert(threwE,'an inverted bracket xLo > xHi did not throw');

%% Equal endpoints are just as bad -- there is no interval to bisect:
            threwE = false;
    try
        coorbital.util.rangeSolve(10,fCube,2,2,1e-9);
    catch err
            threwE = true;
        assert(strcmp(err.identifier,'coorbital:rangeSolve:badBracket'), ...
            'degenerate bracket threw the wrong identifier: %s',err.identifier);
    end
    assert(threwE,'a degenerate bracket xLo == xHi did not throw');

%% A non-positive or non-finite tolerance can never be satisfied and would
%% spin the solver to its iteration cap:
            threwE = false;
    try
        coorbital.util.rangeSolve(10,fCube,0,3,0);
    catch err
            threwE = true;
        assert(strcmp(err.identifier,'coorbital:rangeSolve:badTolerance'), ...
            'zero tolerance threw the wrong identifier: %s',err.identifier);
    end
    assert(threwE,'a zero tolerance did not throw');

%% A failed propagation handing back NaN must stop the solve rather than let
%% the bracket wander on a meaningless comparison:
            threwE = false;
    try
        coorbital.util.rangeSolve(10,@(x) NaN,0,3,1e-9);
    catch err
            threwE = true;
        assert(strcmp(err.identifier,'coorbital:rangeSolve:nonFiniteEval'), ...
            'a NaN evaluation threw the wrong identifier: %s',err.identifier);
    end
    assert(threwE,'a NaN-returning fRange did not throw');

%% -------------------------------------------------------------------
%% Evaluation accounting, and no endpoint evaluated twice
%% -------------------------------------------------------------------
%% Each evaluation of fRange costs a full trajectory propagation in the
%% intended use, several seconds apiece, so a solver that re-evaluates an
%% endpoint it has already paid for is a real defect and not a stylistic one.
%% The count kept here is the test's own, incremented inside the nested
%% counter, and is compared against what the solver claims:
             evalX = zeros(1,512);
           nCalled = 0;
    [~,~,info] = coorbital.util.rangeSolve(10,@countedCube,0,3,1);
    assert(nCalled == info.nEval, ...
        'info.nEval says %d but fRange was actually called %d times', ...
        info.nEval,nCalled);
    assert(nCalled == 6, ...
        'the hand-traced solve should call fRange 6 times, called %d',nCalled);
    assert(info.nEval == info.iterations + 2, ...
        'nEval must be the two endpoints plus one per iteration: %d vs %d + 2', ...
        info.nEval,info.iterations);

%% No abscissa was visited twice, which covers the endpoints and the
%% midpoints together:
             usedX = evalX(1:nCalled);
    assert(numel(unique(usedX)) == nCalled, ...
        'fRange was called twice at the same abscissa: %s',mat2str(usedX));

%% The endpoints specifically, once each and no more:
    assert(sum(usedX == 0) == 1, ...
        'the lower endpoint was evaluated %d times, not once',sum(usedX == 0));
    assert(sum(usedX == 3) == 1, ...
        'the upper endpoint was evaluated %d times, not once',sum(usedX == 3));

%% The first two calls are the bracket endpoints, in either order, and the
%% remaining four are the hand-traced midpoints in order:
    assert(isequal(sort(usedX(1:2)),[0 3]), ...
        'the first two evaluations must be the bracket endpoints, were %s', ...
        mat2str(usedX(1:2)));
    assert(isequal(usedX(3:6),[1.5 2.25 1.875 2.0625]), ...
        'midpoint sequence wrong: %s',mat2str(usedX(3:6)));

%% -------------------------------------------------------------------
%% The info struct carries everything the contract promises
%% -------------------------------------------------------------------
    [~,~,info] = coorbital.util.rangeSolve(10,fCube,0,3,1e-9);
            wanted = {'xLo','xHi','fLo','fHi','fMin','fMax','fTarget','tolF', ...
                    'residual','iterations','nEval','converged','identifier', ...
                    'message'};
    for k = 1:numel(wanted)
        assert(isfield(info,wanted{k}), ...
            'info is missing the field %s',wanted{k});
    end
    assert(info.fTarget == 10,'info.fTarget wrong: %.17g',info.fTarget);
    assert(info.tolF == 1e-9,'info.tolF wrong: %.17g',info.tolF);
    assert(abs(info.residual) <= 1e-9, ...
        'info.residual outside the tolerance: %.17g',info.residual);
    assert(islogical(info.converged), ...
        'info.converged must be logical, is %s',class(info.converged));
    assert(isempty(info.identifier), ...
        'a converged solve must leave identifier empty, got %s',info.identifier);

    function y = countedCube(x)
    %% Purpose:
    %
    %  The cube reference function, wrapped so that the test can count how
    %  many times the solver calls it and at which abscissae. Nested rather
    %  than local so that it shares evalX and nCalled with the parent
    %  workspace; the log is preallocated and filled by index, never grown.
    %
    %% Inputs:
    %
    %  x                [1 x 1]                     Abscissa (dimensionless)
    %
    %% Outputs:
    %
    %  y                [1 x 1]                     x^3 + 2 (dimensionless)
    %
    %% Revision History:
    %  Michael Casey                                                08/07/2026
    %  Copyright 2026 Coorbital, Inc.
    %% ------------------------ Begin Code Sequence ---------------------------
           nCalled = nCalled + 1;
    evalX(nCalled) = x;
                 y = x.^3 + 2;
    end
end
