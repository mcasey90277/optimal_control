function test_aimSolve()
%% Purpose:
%
%  Verify the two-axis damped-Newton solver against SYNTHETIC residuals whose
%  roots, Jacobians and failure modes are known in closed form, so that the
%  solver is exercised entirely independently of any trajectory propagation.
%  Nothing here integrates an equation of motion; the whole point is that a
%  solver bug and a dynamics bug must not be able to hide behind one another.
%
%  Seven reference residuals are used, each chosen to isolate one behaviour:
%
%      fLin     an exactly linear 2x2 system with root [0.42;0.31]. Newton
%               solves it in ONE iteration, so the evaluation accounting can
%               be asserted by hand: 1 + 3*1 = 4 calls
%      fMild    a quadratically coupled system with the same root, for
%               ordinary nonlinear convergence
%      fAtan    1e4*atan(x1 - 0.3) in the first component. The classic Newton
%               divergence: from x1 = 2.3 the FULL step lands at -3.24, where
%               the residual is LARGER than where it started, and an undamped
%               method never comes back
%      fNoRoot  1e4*(atan(x1) - 2). atan is bounded by pi/2 < 2, so there is
%               no root anywhere and the solve must refuse
%      fSing    two components that are the same linear combination of the
%               controls, so the Jacobian columns are parallel: the
%               ill-conditioning of a range maximum, in closed form
%      fGuard   fAtan behind a coordinate guard that THROWS for x1 < -3,
%               which is exactly where the undamped first step lands
%      fProbe   a residual that returns at ONE point and throws everywhere
%               else, so the guess flies but neither difference probe does
%
%  Several assertions are hand-traced rather than written as tolerances. On
%  fAtan from x1 = 2.3 the sequence of x1 - 0.3 is 2, -0.767872, 0.273082,
%  -0.013380, 1.6e-6, ~0, and the corresponding residual magnitudes are
%  11071.5, 6548.4, 2665.8, 133.8, 0.016, ~0 metres. Those numbers are worked
%  out below where they are used.
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
%  Michael Casey                                                08/08/2026
%  Copyright 2026 Coorbital, Inc.
%% ------------------------ Begin Code Sequence ---------------------------

%% The counted reference function shares this workspace through the nested
%% function at the bottom of the file. Preallocated, never grown; nCalled is
%% the fill level and the call count the TEST keeps for itself, independently
%% of whatever the solver reports in info.nEval:
             evalX = zeros(2,512);
           nCalled = 0;

%% The synthetic references. Residuals are scaled to metres so that the
%% tolerances below read the way a targeting caller's would:
              aLin = [ 3.0e6, 1.0e6 ; -0.5e6, 2.0e6 ];
             rRoot = [0.42;0.31];
              fLin = @(x) aLin*(x(:) - rRoot);
             fMild = @(x) [1.0e5.*((x(1) - 0.42) + 0.35.*(x(2) - 0.31).^2); ...
                           1.0e5.*(0.25.*(x(1) - 0.42).^2 + (x(2) - 0.31))];
             fAtan = @(x) [1.0e4.*atan(x(1) - 0.3);1.0e4.*(x(2) - 0.7)];
           fNoRoot = @(x) [1.0e4.*(atan(x(1)) - 2.0);1.0e4.*(x(2) - 0.5)];
             fSing = @(x) [1.0e4.*(x(1) + x(2) - 1.0); ...
                           2.0e4.*(x(1) + x(2) - 1.5)];
               hFD = [1.0e-5;1.0e-5];

%% -------------------------------------------------------------------
%% A well-conditioned LINEAR system, solved in one iteration
%% -------------------------------------------------------------------
%% aLin is nonsingular with determinant 6.5e12, so the root of
%% aLin*(x - rRoot) = 0 is x = rRoot exactly and the forward-difference
%% Jacobian IS aLin up to rounding. Newton on a linear system is exact, so
%% ONE iteration must do it and the evaluation count is hand-countable:
%% one call at the initial guess, two finite-difference probes, one trial
%% step, four in all. Anything more means a base point was re-evaluated:
    [xSol,fAch,info] = coorbital.util.aimSolve(fLin,[0.10;0.90],hFD,1.0e-3);
    assert(info.converged,'the linear case failed to converge: %s',info.message);
    assert(info.iterations == 1, ...
        'Newton on a linear system must take exactly 1 iteration, took %d', ...
        info.iterations);
    assert(info.nEval == 4, ...
        ['a one-iteration solve costs 1 initial + 2 Jacobian + 1 trial = 4 ', ...
         'evaluations, it cost %d'],info.nEval);
    assert(max(abs(xSol - rRoot)) < 1.0e-9, ...
        'linear root wrong: got [%.17g,%.17g], analytic root is [0.42,0.31]', ...
        xSol(1),xSol(2));
    assert(max(abs(fAch)) <= 1.0e-3, ...
        'linear case achieved [%.17g,%.17g], outside the 1e-3 tolerance', ...
        fAch(1),fAch(2));

%% The full Newton step was taken, so no damping was needed and none was
%% applied. A solver that damped unconditionally would fail here:
    assert(isscalar(info.lambda) && info.lambda == 1, ...
        'the linear case should take the full step, lambda = %s', ...
        mat2str(info.lambda));
    assert(info.nShrink == 0, ...
        'the linear case should need no step shrinkage, took %d',info.nShrink);

%% fAch must be the value the CALLER'S residual takes at xSol, not a
%% separately bookkept number that could have drifted from it:
    assert(isequal(fAch,fLin(xSol)), ...
        'fAch [%.17g,%.17g] is not fLin(xSol)',fAch(1),fAch(2));

%% The output is a COLUMN whatever the caller passed, and a row guess is
%% accepted. A targeting script that builds its guess with [psi, ctrl] must
%% not get back something that broadcasts into a 2x2 downstream:
    assert(iscolumn(xSol) && numel(xSol) == 2, ...
        'xSol must be a [2 x 1] column, is %s',mat2str(size(xSol)));
    assert(iscolumn(fAch) && numel(fAch) == 2, ...
        'fAch must be a [2 x 1] column, is %s',mat2str(size(fAch)));
              xRow = coorbital.util.aimSolve(fLin,[0.10, 0.90],hFD',1.0e-3);
    assert(iscolumn(xRow) && max(abs(xRow - rRoot)) < 1.0e-9, ...
        'a row-vector guess and step must give the same column answer: %s', ...
        mat2str(xRow));

%% -------------------------------------------------------------------
%% A guess already at the root costs one evaluation and no iterations
%% -------------------------------------------------------------------
    [xAt,fAt,infoAt] = coorbital.util.aimSolve(fLin,rRoot,hFD,1.0e-3);
    assert(infoAt.converged && infoAt.iterations == 0, ...
        'a guess at the root must converge with 0 iterations, took %d', ...
        infoAt.iterations);
    assert(infoAt.nEval == 1, ...
        'a guess at the root costs one evaluation, cost %d',infoAt.nEval);
    assert(isequal(xAt,rRoot) && all(fAt == 0), ...
        'a guess at the root must be returned untouched: %s',mat2str(xAt));

%% -------------------------------------------------------------------
%% A mildly NONLINEAR system converges
%% -------------------------------------------------------------------
%% fMild has the same root but a Jacobian that changes with x, so the first
%% Newton step overshoots and the method has to iterate. It is still
%% well-conditioned everywhere on the path, so no damping should be needed:
    [xSol,fAch,info] = coorbital.util.aimSolve(fMild,[0.90;-0.20],hFD,1.0e-3);
    assert(info.converged, ...
        'the mildly nonlinear case failed to converge: %s',info.message);
    assert(max(abs(xSol - rRoot)) < 1.0e-8, ...
        'nonlinear root wrong: got [%.17g,%.17g]',xSol(1),xSol(2));
    assert(max(abs(fAch)) <= 1.0e-3, ...
        'nonlinear case achieved [%.17g,%.17g], outside tolerance', ...
        fAch(1),fAch(2));
    assert(info.iterations > 1, ...
        ['a nonlinear system should need more than one Newton step; this ', ...
         'one took %d, which suggests the residual is not what the test ', ...
         'thinks it is'],info.iterations);
    assert(info.iterations <= 10, ...
        'quadratic convergence should not take %d iterations',info.iterations);

%% Every Jacobian on the path was reported, one per iteration, and each was
%% well conditioned. This is the field the caller reads to decide whether an
%% answer that DID converge is nonetheless sitting on a knife edge:
    assert(size(info.jacobians,1) == 2 && size(info.jacobians,2) == 2 && ...
            size(info.jacobians,3) == info.iterations, ...
        'the Jacobian history should be [2 2 %d], is %s', ...
        info.iterations,mat2str(size(info.jacobians)));
    assert(numel(info.condJ) == info.iterations, ...
        'one condition number per Jacobian: %d against %d iterations', ...
        numel(info.condJ),info.iterations);
    assert(all(isfinite(info.condJ)) && all(info.condJ >= 1), ...
        'condition numbers must be finite and at least one: %s', ...
        mat2str(info.condJ));

%% -------------------------------------------------------------------
%% DAMPING IS WHAT MAKES fAtan CONVERGE AT ALL
%% -------------------------------------------------------------------
%% Write u = x1 - 0.3. Newton on 1e4*atan(u) is the textbook divergent case:
%% the step is -atan(u)*(1 + u^2), which for |u| > 1.3917 overshoots past the
%% root and lands further out than it started. From u = 2:
%%
%%   full step  -atan(2)*5 = -5.5357   ->  u = -3.5357, |f| = 12951 > 11071
%%                                          REJECTED, halve
%%   half step  -2.7679                ->  u = -0.767872, |f| = 6548 < 11071
%%                                          accepted
%%   thereafter u = 0.273082, -0.013380, 1.6e-6, ~0 -- quadratic, with
%%              |f| = 2665.82, 133.79, 0.0160, ~0
%%
%% An undamped Newton from the same guess goes u = 2, -3.54, 13.9, -280, ...
%% and never returns, so this single case is the whole justification for the
%% line search. The second component starts AT its root, so f2 = 0 throughout
%% and max(abs(f)) is |f1| alone, which makes the trace above exact:
    [xSol,fAch,info] = coorbital.util.aimSolve(fAtan,[2.30;0.70],hFD,1.0e-3);
    assert(info.converged, ...
        'the damped case failed to converge: %s',info.message);
    assert(abs(xSol(1) - 0.3) < 1.0e-9 && abs(xSol(2) - 0.7) < 1.0e-12, ...
        'damped root wrong: got [%.17g,%.17g], analytic root is [0.3,0.7]', ...
        xSol(1),xSol(2));
    assert(max(abs(fAch)) <= 1.0e-3, ...
        'damped case achieved [%.17g,%.17g], outside tolerance', ...
        fAch(1),fAch(2));

%% The step really was shrunk. Without this assertion a solver that had
%% deleted its line search would still pass the convergence check above on
%% every OTHER case in this file, because every other case takes full steps:
    assert(info.nShrink >= 1, ...
        ['the divergent guess must have forced at least one step halving, ', ...
         'and none was recorded -- the line search is not running']);
    assert(any(info.lambda < 1), ...
        'no iteration reported a damped step: lambda = %s',mat2str(info.lambda));
    assert(all(info.lambda > 0) && all(info.lambda <= 1), ...
        'lambda must lie in (0,1]: %s',mat2str(info.lambda));

%% The first accepted iterate is the HALF step, not the full one. This is the
%% hand-traced number: 2.3 - 0.5*5.5357 = -0.4679 to four decimals, and the
%% full step would have put it at -3.2357 instead:
             xNext = info.xHist(1,2);
    assert(abs(xNext - (-0.467872)) < 1.0e-4, ...
        ['the first accepted iterate should be the HALVED step at ', ...
         '-0.467872, it was %.17g'],xNext);

%% -------------------------------------------------------------------
%% A tighter tolerance gives a tighter result, and costs more
%% -------------------------------------------------------------------
%% From the trace above the residual magnitudes at successive iterates are
%% 11071, 6548, 2665.8, 133.8, 0.016, ~0. A tolerance of 5000 m therefore
%% stops after the SECOND iteration at 2665.8 m; a tolerance of 1e-3 m has to
%% go further. Both bounds are asserted as well as the comparison, so a solver
%% that ignored tolF and always ran to machine precision would be caught by
%% the iteration-count assertion rather than slipping through:
    [~,fLoose,infoLoose] = coorbital.util.aimSolve(fAtan,[2.30;0.70],hFD,5.0e3);
    [~,fTight,infoTight] = coorbital.util.aimSolve(fAtan,[2.30;0.70],hFD,1.0e-3);
          errLoose = max(abs(fLoose));
          errTight = max(abs(fTight));
    assert(errLoose <= 5.0e3, ...
        'the loose solve violated its own tolerance: %.17g',errLoose);
    assert(errTight <= 1.0e-3, ...
        'the tight solve violated its own tolerance: %.17g',errTight);
    assert(errTight < errLoose, ...
        'a tighter tolerance did not give a tighter result: %.17g against %.17g', ...
        errTight,errLoose);
    assert(infoTight.iterations > infoLoose.iterations, ...
        'a tighter tolerance did not cost more iterations: %d against %d', ...
        infoTight.iterations,infoLoose.iterations);
    assert(infoLoose.iterations == 2, ...
        ['the hand trace stops the 5000 m solve after 2 iterations at ', ...
         '2665.8 m; it took %d'],infoLoose.iterations);
    assert(abs(errLoose - 2665.82) < 1.0, ...
        'the loose solve should end on 2665.82 m, it ended on %.17g',errLoose);

%% -------------------------------------------------------------------
%% NO ROOT IN REACH: refuse, do not throw, and carry the best point
%% -------------------------------------------------------------------
%% atan is bounded by pi/2 = 1.5708, so 1e4*(atan(x1) - 2) can never be
%% smaller than 4292 m in magnitude and there is no root anywhere on the
%% plane. Newton marches x1 off to the right forever, reducing the residual
%% the whole way but never reaching zero. opts.maxStep caps each march at
%% 0.5, which both keeps the Jacobian well conditioned -- it would otherwise
%% blow past the condMax ceiling and refuse for the wrong reason -- and makes
%% the answer exactly predictable: ten iterations of 0.5 from zero is x1 = 5.
%%
%% The call below is deliberately NOT wrapped in try/catch: if it throws, the
%% test fails on the throw itself, which is the point. A targeting script
%% that silently returned a near miss would be worse than one that refused,
%% and converged = false is the whole guard against that:
              oCap = struct('maxIter',10,'maxStep',[0.5;0.5]);
    [xSol,fAch,info] = coorbital.util.aimSolve(fNoRoot,[0.0;0.5],hFD,1.0,oCap);
    assert(~info.converged, ...
        'a residual with no root was reported as converged');
    assert(strcmp(info.identifier,'coorbital:aimSolve:toleranceNotMet'), ...
        'wrong refusal identifier on the no-root case: ''%s''',info.identifier);
    assert(info.iterations == 10, ...
        'the no-root case should run its 10-iteration cap, ran %d', ...
        info.iterations);

%% Best-so-far, and it really is better than where it started: 6266 m against
%% 20000 m at the guess. A solver that returned its last midpoint without
%% tracking the best would still pass here, but one that returned the initial
%% guess unchanged would not:
    assert(all(isfinite(xSol)) && all(isfinite(fAch)), ...
        'a refusal must still return finite numbers, got x = %s f = %s', ...
        mat2str(xSol),mat2str(fAch));
    assert(max(abs(fAch)) < max(abs(fNoRoot([0.0;0.5]))), ...
        ['the refusal returned %.17g, no better than the %.17g it started ', ...
         'with'],max(abs(fAch)),max(abs(fNoRoot([0.0;0.5]))));
    assert(abs(max(abs(fAch)) - 6265.99) < 1.0, ...
        'the no-root case should end on 6266 m, it ended on %.17g', ...
        max(abs(fAch)));
    assert(isequal(fAch,info.residual) && info.errMax == max(abs(fAch)), ...
        'info.residual and info.errMax must describe the returned point');

%% opts.maxStep was honoured on every single iteration. This is the guard
%% that stops a full Newton step from a poor guess throwing the propagation
%% into a coordinate guard, so it is asserted on the whole path, not just the
%% endpoints:
             dStep = abs(diff(info.xHist,1,2));
    assert(all(dStep(:) <= 0.5 + 1.0e-9), ...
        'a step of %.17g exceeded the 0.5 cap in opts.maxStep',max(dStep(:)));
    assert(abs(xSol(1) - 5.0) < 1.0e-9, ...
        'ten capped steps of 0.5 from zero must land on 5, landed on %.17g', ...
        xSol(1));

%% The message names what was wanted, what was reached and what to change:
    assert(contains(info.message,'Failed to reach'), ...
        'the refusal message must say so plainly: %s',info.message);
    assert(contains(info.message,'opts.maxIter'), ...
        'the refusal message must name a remedy: %s',info.message);

%% -------------------------------------------------------------------
%% A SINGULAR JACOBIAN is detected and reported, not turned into infinity
%% -------------------------------------------------------------------
%% Both components of fSing depend on x1 and x2 only through the sum x1 + x2,
%% so the Jacobian columns are identical up to scale: [1e4 1e4; 2e4 2e4],
%% exactly rank one. This is the closed-form version of what happens near a
%% maximum-range control, where the downrange sensitivity collapses and the
%% two controls stop moving the two residuals independently. The two
%% components also disagree -- one wants x1 + x2 = 1 and the other 1.5 -- so
%% there is no root to find even if the step could be computed.
%%
%% What must NOT happen is a backslash against a rank-one matrix producing an
%% infinite step, or a near-singularity warning escaping to the console:
    [xSol,fAch,info] = coorbital.util.aimSolve(fSing,[0.10;0.20],hFD,1.0);
    assert(~info.converged, ...
        'a singular system was reported as converged');
    assert(strcmp(info.identifier,'coorbital:aimSolve:singularJacobian'), ...
        'wrong identifier on the singular case: ''%s''',info.identifier);
    assert(info.iterations == 1, ...
        'the singular case should refuse on its first Jacobian, took %d', ...
        info.iterations);
    assert(all(isfinite(xSol)) && all(isfinite(fAch)), ...
        'the singular case produced non-finite output: x = %s f = %s', ...
        mat2str(xSol),mat2str(fAch));
    assert(isequal(xSol,[0.10;0.20]), ...
        'with no step taken the best point is the initial guess, got %s', ...
        mat2str(xSol));

%% The condition number is REPORTED, which is the whole point: the caller
%% wants to print it, not to guess at it:
    assert(isscalar(info.condJ) && info.condJ > 1.0e8, ...
        ['the rank-one Jacobian should report a condition number above the ', ...
         '1e8 ceiling, reported %s'],mat2str(info.condJ));
    assert(size(info.jacobians,3) == 1 && ...
            isequal(info.jacobians(:,:,1),info.jacobians), ...
        'the offending Jacobian must be in the history, size is %s', ...
        mat2str(size(info.jacobians)));
    assert(contains(info.message,'condition number'), ...
        'the singular message must quote the condition number: %s',info.message);

%% Two Jacobian evaluations on top of the initial one, and no trial step was
%% ever attempted because no trustworthy step existed:
    assert(info.nEval == 3, ...
        ['refusing on the first Jacobian costs 1 initial + 2 probes = 3 ', ...
         'evaluations, it cost %d'],info.nEval);

%% A ceiling loose enough to let the rank-one system through must still not
%% produce an infinity: the fallback guard on the step itself catches it:
            oLoose = struct('condMax',1.0e18,'maxIter',3);
    [xL,fL,infoL] = coorbital.util.aimSolve(fSing,[0.10;0.20],hFD,1.0,oLoose);
    assert(~infoL.converged, ...
        'the singular system converged with a loosened condMax');
    assert(all(isfinite(xL)) && all(isfinite(fL)), ...
        'a loosened condMax let a non-finite result out: x = %s f = %s', ...
        mat2str(xL),mat2str(fL));

%% -------------------------------------------------------------------
%% A RESIDUAL THAT THROWS is handled by shrinkage, not by aborting
%% -------------------------------------------------------------------
%% fGuard is fAtan behind a coordinate guard that throws for x1 < -3, which
%% is precisely where the undamped first step from x1 = 2.3 lands: -3.2357.
%% So the first trial evaluation does not merely fail to improve, it fails to
%% return at all, and the solve must treat that as an infeasible control pair
%% and halve. The halved step lands at -0.4679, well inside the guard.
%%
%% Not wrapped in try/catch, again deliberately: a propagation failure inside
%% the line search is an ordinary event in targeting and must not escape:
    [xSol,fAch,info] = coorbital.util.aimSolve(@guardResid,[2.30;0.70],hFD,1.0e-3);
    assert(info.converged, ...
        'the guarded case failed to converge: %s',info.message);
    assert(max(abs(xSol - [0.3;0.7])) < 1.0e-9, ...
        'guarded root wrong: got [%.17g,%.17g]',xSol(1),xSol(2));
    assert(max(abs(fAch)) <= 1.0e-3, ...
        'guarded case achieved [%.17g,%.17g], outside tolerance', ...
        fAch(1),fAch(2));
    assert(info.nInfeasible >= 1, ...
        ['the guard must have been hit at least once and counted; ', ...
         'nInfeasible is %d'],info.nInfeasible);
    assert(info.nShrink >= info.nInfeasible, ...
        'every infeasible trial must also have shrunk the step: %d against %d', ...
        info.nShrink,info.nInfeasible);

%% A guard that throws EVERYWHERE around the guess cannot be worked around by
%% shrinking, because the finite-difference probe fails too. That refuses,
%% with the underlying failure quoted so the caller can see what went wrong,
%% and it still does not throw:
    [~,~,infoG] = coorbital.util.aimSolve(@probeThrow,[1.0;1.0],hFD,1.0e-3);
    assert(~infoG.converged, ...
        'a residual that never evaluates was reported as converged');
    assert(strcmp(infoG.identifier,'coorbital:aimSolve:jacobianFailed'), ...
        'wrong identifier when the probe point fails: ''%s''',infoG.identifier);
    assert(contains(infoG.message,'no flight here'), ...
        'the message must quote the underlying failure: %s',infoG.message);

%% -------------------------------------------------------------------
%% THE EVALUATION COUNT IS TRUTHFUL, and no point is evaluated twice
%% -------------------------------------------------------------------
%% Each evaluation is a full trajectory propagation in the intended use, and
%% Task 3's caller nests a whole solve inside one of them, so a solver that
%% re-evaluates a point it has already paid for is a real defect and not a
%% stylistic one. The count kept here is the TEST'S own, incremented inside
%% the nested counter, and is compared against what the solver claims:
             evalX = zeros(2,512);
           nCalled = 0;
    [~,~,infoC] = coorbital.util.aimSolve(@countedMild,[0.90;-0.20],hFD,1.0e-3);
    assert(nCalled == infoC.nEval, ...
        'info.nEval says %d but fResid was actually called %d times', ...
        infoC.nEval,nCalled);
    assert(infoC.nEval == 1 + 3*infoC.iterations + infoC.nShrink, ...
        ['the accounting must be one initial call plus three per iteration ', ...
         'plus one per rejected trial: %d against 1 + 3*%d + %d'], ...
        infoC.nEval,infoC.iterations,infoC.nShrink);

%% No control pair was visited twice, which covers the base points and the
%% finite-difference probes together:
             usedX = evalX(:,1:nCalled);
    assert(size(unique(usedX','rows'),1) == nCalled, ...
        'fResid was called twice at the same control pair');

%% And specifically: every ACCEPTED iterate was evaluated exactly once. The
%% initial guess is evaluated at the top and never again, and each accepted
%% trial value is carried forward as the next base value rather than being
%% recomputed for the Jacobian. This is the assertion that a solver
%% re-evaluating its base point each iteration fails:
    for kIt = 1:size(infoC.xHist,2)
             nHits = sum(all(usedX == infoC.xHist(:,kIt),1));
        assert(nHits == 1, ...
            'iterate %d was evaluated %d times, not once',kIt,nHits);
    end
    assert(size(infoC.xHist,2) == infoC.iterations + 1, ...
        'xHist should hold the guess plus one point per iteration: %d against %d', ...
        size(infoC.xHist,2),infoC.iterations + 1);
    assert(isequal(size(infoC.fHist),size(infoC.xHist)), ...
        'xHist and fHist must match in size: %s against %s', ...
        mat2str(size(infoC.xHist)),mat2str(size(infoC.fHist)));

%% -------------------------------------------------------------------
%% The info struct carries everything the contract promises
%% -------------------------------------------------------------------
    [~,~,info] = coorbital.util.aimSolve(fLin,[0.10;0.90],hFD,1.0e-3);
            wanted = {'x0','dx0','tolF','opts','residual','errMax', ...
                      'iterations','condJ','jacobians','lambda','xHist', ...
                      'fHist','nEval','nShrink','nInfeasible','converged', ...
                      'identifier','message'};
    for k = 1:numel(wanted)
        assert(isfield(info,wanted{k}), ...
            'info is missing the field %s',wanted{k});
    end
    assert(isequal(info.x0,[0.10;0.90]), ...
        'info.x0 must carry the guess as a column: %s',mat2str(info.x0));
    assert(isequal(info.dx0,hFD), ...
        'info.dx0 must carry the step actually used: %s',mat2str(info.dx0));
    assert(info.tolF == 1.0e-3,'info.tolF wrong: %.17g',info.tolF);
    assert(islogical(info.converged), ...
        'info.converged must be logical, is %s',class(info.converged));
    assert(isempty(info.identifier), ...
        'a converged solve must leave identifier empty, got %s',info.identifier);
    assert(contains(info.message,'Converged'), ...
        'a converged solve must say so: %s',info.message);
    assert(info.opts.maxIter == 25 && info.opts.maxHalve == 12, ...
        'info.opts must report the defaults actually used');
    assert(isequal(info.opts.maxStep,[Inf;Inf]), ...
        'the default step cap should be unbounded: %s',mat2str(info.opts.maxStep));

%% -------------------------------------------------------------------
%% Malformed input is an error, not a refusal
%% -------------------------------------------------------------------
%% A NON-FINITE GUESS IS THE DEFECT THAT HIDES. Every test in the solver is
%% an inequality and every inequality against NaN is false, so a NaN guess
%% would neither converge nor trip a refusal: it would run to the iteration
%% cap and report a tolerance it could not meet, as though the caller's
%% residual were at fault:
    for badX = {[1;2;3],[NaN;1],[Inf;1],[1i;1],[],0.5,'ab'}
        expectThrow(@() coorbital.util.aimSolve(fLin,badX{1},hFD,1.0e-3), ...
            'coorbital:aimSolve:badGuess', ...
            sprintf('guess %s',mat2str(badX{1})));
    end

%% A SCALAR dx0 IS REFUSED, and this is a deliberate refusal rather than a
%% convenience broadcast. The two controls are on unrelated scales -- radians
%% of heading against radians of loft or a dimensionless burn fraction -- so
%% one difference step cannot be the right chord for both, and silently
%% accepting it would bury a wrong Jacobian column where nobody looks:
    expectThrow(@() coorbital.util.aimSolve(fLin,[0.10;0.90],1.0e-5,1.0e-3), ...
        'coorbital:aimSolve:badStep','a scalar difference step');
    for badH = {[0;1e-5],[-1e-5;1e-5],[NaN;1e-5],[1e-5;1e-5;1e-5]}
        expectThrow(@() coorbital.util.aimSolve(fLin,[0.10;0.90],badH{1},1.0e-3), ...
            'coorbital:aimSolve:badStep', ...
            sprintf('difference step %s',mat2str(badH{1})));
    end

%% A non-positive or non-finite tolerance can never be met and would spin the
%% solve to its iteration cap for nothing:
    for badT = {0,-1,NaN,Inf,[1 2]}
        expectThrow(@() coorbital.util.aimSolve(fLin,[0.10;0.90],hFD,badT{1}), ...
            'coorbital:aimSolve:badTolerance', ...
            sprintf('tolerance %s',mat2str(badT{1})));
    end

%% An unrecognised option is an ERROR and not a silent no-op. A caller who
%% writes opts.maxIterations and watches nothing change has no way to tell
%% that from a solver that honoured it:
    expectThrow(@() coorbital.util.aimSolve(fLin,[0.10;0.90],hFD,1.0e-3, ...
        struct('maxIterations',5)),'coorbital:aimSolve:badOptions', ...
        'a misspelt option');
    expectThrow(@() coorbital.util.aimSolve(fLin,[0.10;0.90],hFD,1.0e-3, ...
        struct('maxIter',0)),'coorbital:aimSolve:badOptions','maxIter of zero');
    expectThrow(@() coorbital.util.aimSolve(fLin,[0.10;0.90],hFD,1.0e-3, ...
        struct('maxStep',[0;1])),'coorbital:aimSolve:badOptions', ...
        'a zero step cap');
    expectThrow(@() coorbital.util.aimSolve(fLin,[0.10;0.90],hFD,1.0e-3, ...
        struct('condMax',0.5)),'coorbital:aimSolve:badOptions', ...
        'a condMax below one');

%% Not a handle at all:
    expectThrow(@() coorbital.util.aimSolve(42,[0.10;0.90],hFD,1.0e-3), ...
        'coorbital:aimSolve:badResidual','a non-handle residual');

%% A residual that fails AT THE INITIAL GUESS is an unusable guess rather
%% than a failure to converge: there is no step to shrink and no earlier
%% point to fall back on, so this one throws:
    expectThrow(@() coorbital.util.aimSolve(@(x) [NaN;0],[0.10;0.90],hFD,1.0e-3), ...
        'coorbital:aimSolve:nonFiniteEval','a NaN residual at the guess');
    expectThrow(@() coorbital.util.aimSolve(@(x) [1;2;3],[0.10;0.90],hFD,1.0e-3), ...
        'coorbital:aimSolve:nonFiniteEval','a three-vector residual');

%% -------------------------------------------------------------------
%% Defaults: dx0 and tolF may be omitted or left empty
%% -------------------------------------------------------------------
    [xDef,~,infoDef] = coorbital.util.aimSolve(fLin,[0.10;0.90]);
    assert(infoDef.converged && max(abs(xDef - rRoot)) < 1.0e-6, ...
        'the all-defaults call did not solve the linear case: %s', ...
        infoDef.message);
    assert(infoDef.tolF == 1.0, ...
        'the default tolerance should be one metre, is %.17g',infoDef.tolF);
    assert(isequal(infoDef.dx0,1.0e-4.*max(abs([0.10;0.90]),1)), ...
        ['the default difference step must be relative with an absolute ', ...
         'floor, per control: %s'],mat2str(infoDef.dx0));
    [~,~,infoEmpty] = coorbital.util.aimSolve(fLin,[0.10;0.90],[],[],struct());
    assert(isequal(infoEmpty.dx0,infoDef.dx0) && infoEmpty.tolF == infoDef.tolF, ...
        'an empty argument must mean the same as an omitted one');

    function f = guardResid(x)
    %% Purpose:
    %
    %  The atan reference behind a coordinate guard that THROWS below
    %  x1 = -3. The undamped first Newton step from x1 = 2.3 lands at
    %  -3.2357, inside the guard; the halved step lands at -0.4679, outside
    %  it. So the guard is hit exactly once, on the first trial of the first
    %  iteration, which is what makes the shrinkage path deterministic.
    %
    %% Inputs:
    %
    %  x                [2 x 1]                     Controls (rad)
    %
    %% Outputs:
    %
    %  f                [2 x 1]                     Residuals (m)
    %
    %% Revision History:
    %  Michael Casey                                                08/08/2026
    %  Copyright 2026 Coorbital, Inc.
    %% ------------------------ Begin Code Sequence ---------------------------
        if x(1) < -3.0
            error('test:aimSolve:guard', ...
                'the vehicle left the atmosphere model at x1 = %.6g',x(1));
        end
                 f = [1.0e4.*atan(x(1) - 0.3);1.0e4.*(x(2) - 0.7)];
    end

    function f = probeThrow(x)
    %% Purpose:
    %
    %  A residual that evaluates at ONE point and nowhere else, standing in
    %  for a guess whose whole neighbourhood fails to propagate. The initial
    %  guess flies, so the solve gets past the entry check, but both
    %  finite-difference probes fail and there is nothing to shrink toward.
    %  The solve must refuse on the Jacobian rather than grind through a line
    %  search it has no direction for.
    %
    %% Inputs:
    %
    %  x                [2 x 1]                     Controls (rad)
    %
    %% Outputs:
    %
    %  f                [2 x 1]                     Residuals (m), returned
    %                                               only at x = [1;1]
    %
    %% Revision History:
    %  Michael Casey                                                08/08/2026
    %  Copyright 2026 Coorbital, Inc.
    %% ------------------------ Begin Code Sequence ---------------------------
        if ~isequal(x(:),[1.0;1.0])
            error('test:aimSolve:probe','no flight here at all');
        end
                 f = [1.0e4;0];
    end

    function f = countedMild(x)
    %% Purpose:
    %
    %  The mildly nonlinear reference, wrapped so that the test can count how
    %  many times the solver calls it and at which control pairs. Nested
    %  rather than local so that it shares evalX and nCalled with the parent
    %  workspace; the log is preallocated and filled by column index, never
    %  grown.
    %
    %% Inputs:
    %
    %  x                [2 x 1]                     Controls (dimensionless)
    %
    %% Outputs:
    %
    %  f                [2 x 1]                     Residuals (m)
    %
    %% Revision History:
    %  Michael Casey                                                08/08/2026
    %  Copyright 2026 Coorbital, Inc.
    %% ------------------------ Begin Code Sequence ---------------------------
           nCalled = nCalled + 1;
  evalX(:,nCalled) = x(:);
                 f = [1.0e5.*((x(1) - 0.42) + 0.35.*(x(2) - 0.31).^2); ...
                      1.0e5.*(0.25.*(x(1) - 0.42).^2 + (x(2) - 0.31))];
    end
end

function expectThrow(fh,wantId,what)
%% Purpose:
%
%  Assert that a call throws, and throws with the identifier it promised.
%  Local rather than nested because it shares nothing with the test's
%  workspace, and factored out because the same three lines would otherwise
%  appear twenty times over the input-validation block.
%
%% Inputs:
%
%  fh               Function handle             Nullary handle wrapping the
%                                               call under test
%
%  wantId           Char                        Expected MATLAB error
%                                               identifier
%
%  what             Char                        Phrase naming the input, for
%                                               the failure message
%
%% Outputs:
%
%  none                                         Throws on any failed assertion
%
%% Revision History:
%  Michael Casey                                                08/08/2026
%  Copyright 2026 Coorbital, Inc.
%% ------------------------ Begin Code Sequence ---------------------------
            threwE = false;
    try
        fh();
    catch err
            threwE = true;
        assert(strcmp(err.identifier,wantId), ...
            '%s threw ''%s'' rather than ''%s''',what,err.identifier,wantId);
    end
    assert(threwE,'%s did not throw',what);
end
