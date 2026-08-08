function [xSol,fAch,info] = aimSolve(fResid,x0,dx0,tolF,opts)
%% Purpose:
%
%  Drive a TWO-component residual to zero in TWO controls by damped Newton
%  iteration with a forward-difference Jacobian. Written for two-axis
%  targeting, where the controls are the launch heading psi (rad) and a range
%  control -- a loft angle (rad) or a burn fraction (-) -- and the residuals
%  are the downrange miss and the cross-track miss in metres. Nothing here
%  knows that: fResid is any handle taking a two-vector and returning a
%  two-vector, and it is treated as an opaque, expensive black box.
%
%  Two-dimensional sibling of coorbital.util.rangeSolve, and deliberately the
%  same shape of contract: converge on the ACHIEVED residual rather than on
%  the step, never throw merely because the solve did not converge, and hand
%  the caller an info struct rich enough to write its own refusal. Where
%  rangeSolve could bisect -- one control, one monotone residual, a bracket
%  that cannot be left -- two coupled axes offer no bracket to preserve, so
%  the method must be Newton and the safeguards must be explicit instead:
%  a condition-number ceiling on the Jacobian, a cap on the step, and a
%  backtracking line search that treats a propagation failure as an
%  infeasible point rather than as an error.
%
%  Each residual evaluation is a full trajectory propagation in the intended
%  use, order 0.15 s apiece, and Task 3's closed-loop caller nests a whole
%  solve inside one evaluation. THE EVALUATION COUNT IS THE COST. A clean
%  iteration costs exactly three: two finite-difference probes and one trial
%  step. The base point is never re-evaluated -- the accepted trial value
%  becomes the next base value -- and info.nEval reports what was actually
%  spent.
%
%% Inputs:
%
%  fResid           Function handle             f = fResid(x), taking a
%                                               [2 x 1] control vector and
%                                               returning a [2 x 1] residual
%                                               vector. Units are whatever
%                                               fResid returns; metres of
%                                               downrange and cross-track
%                                               miss in the intended use.
%                                               MAY THROW for a control pair
%                                               that does not fly, see Notes
%
%  x0               [2 x 1]                     Initial guess for the two
%                                               controls. Row or column, in
%                                               the controls' own units:
%                                               rad for a heading, rad or a
%                                               dimensionless fraction for
%                                               the range control. Must be
%                                               finite and real
%
%  dx0              [2 x 1]                     Forward-difference step, ONE
%                                               PER CONTROL, in each
%                                               control's own units. Must be
%                                               finite, real and strictly
%                                               positive. A scalar is
%                                               REFUSED, see Notes. Pass []
%                                               or omit for the default
%                                               1e-4.*max(abs(x0),1),
%                                               elementwise
%
%  tolF             [1 x 1]                     Convergence tolerance on the
%                                               ACHIEVED residual, in
%                                               fResid's units, not on the
%                                               step. The test is
%                                               max(abs(f)) <= tolF, so both
%                                               axes must be inside it. Must
%                                               be finite and positive. Pass
%                                               [] or omit for 1.0, which is
%                                               one metre in the intended use
%
%  opts             Struct                      Optional safeguards. Any
%                                               field may be omitted or left
%                                               empty for its default; an
%                                               UNRECOGNISED field is an
%                                               error rather than a silent
%                                               no-op:
%                                               maxIter  [1 x 1] iteration
%                                                        cap, default 25
%                                               maxStep  [2 x 1] largest
%                                                        permitted change in
%                                                        each control per
%                                                        iteration, in the
%                                                        controls' units,
%                                                        default [Inf;Inf].
%                                                        A scalar is
%                                                        broadcast to both
%                                               maxHalve [1 x 1] step
%                                                        halvings allowed per
%                                                        iteration before the
%                                                        solve is declared
%                                                        stalled, default 12
%                                               condMax  [1 x 1] ceiling on
%                                                        the 2-norm condition
%                                                        number of the
%                                                        Jacobian, default
%                                                        1e8
%
%% Outputs:
%
%  xSol             [2 x 1]                     Solved controls, always a
%                                               COLUMN whatever the shape of
%                                               x0. On a failure to
%                                               converge, the BEST-SO-FAR
%                                               point: whichever control pair
%                                               on the search path achieved
%                                               the smallest max(abs(f)),
%                                               which is x0 itself when the
%                                               very first step could not be
%                                               taken
%
%  fAch             [2 x 1]                     The residual fResid actually
%                                               takes at xSol, signed, in
%                                               fResid's units
%
%  info             Struct                      Solve record:
%                                               x0          [2 x 1] as
%                                                           passed, as a
%                                                           column
%                                               dx0         [2 x 1] the
%                                                           difference step
%                                                           actually used
%                                               tolF        [1 x 1] as passed
%                                               opts        Struct, the
%                                                           safeguards
%                                                           actually used,
%                                                           defaults filled
%                                               residual    [2 x 1] equal to
%                                                           fAch, carried so
%                                                           the record stands
%                                                           alone
%                                               errMax      [1 x 1]
%                                                           max(abs(fAch)),
%                                                           the quantity
%                                                           compared against
%                                                           tolF
%                                               iterations  [1 x 1]
%                                                           iterations
%                                                           STARTED, which on
%                                                           a refusal
%                                                           includes the one
%                                                           that refused
%                                               nEval       [1 x 1] calls
%                                                           made to fResid,
%                                                           including the
%                                                           ones that threw
%                                               nShrink     [1 x 1] TRIAL
%                                                           evaluations
%                                                           rejected, each of
%                                                           which halved the
%                                                           step
%                                               nInfeasible [1 x 1] every
%                                                           evaluation that
%                                                           threw or returned
%                                                           a non-finite pair,
%                                                           at a TRIAL or a
%                                                           PROBE point. NOT a
%                                                           subset of nShrink:
%                                                           an infeasible
%                                                           trial shrinks the
%                                                           step and so counts
%                                                           in both, but an
%                                                           infeasible probe
%                                                           ends the solve
%                                                           without any
%                                                           shrinkage and
%                                                           counts only here
%                                               condJ       [1 x nJ] 2-norm
%                                                           condition number
%                                                           of every Jacobian
%                                                           built
%                                               jacobians   [2 x 2 x nJ] the
%                                                           Jacobians
%                                                           themselves, in
%                                                           order, units of
%                                                           residual per
%                                                           control
%                                               lambda      [1 x nL] damping
%                                                           factor of the last
%                                                           step ATTEMPTED in
%                                                           each line search.
%                                                           1 means the full
%                                                           Newton step was
%                                                           accepted at once;
%                                                           on a stall it is
%                                                           the shortest step
%                                                           tried, never the
%                                                           untried halving
%                                                           below it
%                                               xHist       [2 x (nA+1)] x0
%                                                           followed by every
%                                                           ACCEPTED iterate
%                                               fHist       [2 x (nA+1)] the
%                                                           residuals there
%                                               converged   logical
%                                               identifier  char, '' when
%                                                           converged, else a
%                                                           MATLAB error
%                                                           identifier the
%                                                           caller may rethrow
%                                               message     char, always
%                                                           populated
%
%% Notes:
%
%  A FAILURE TO CONVERGE IS NOT AN ERROR. Every way this solve can end badly
%  -- iterations exhausted, Jacobian ill-conditioned, step stalled, a probe
%  point that will not fly -- returns converged = false with an identifier,
%  a sentence, and the best point reached. The refusal is left to the caller
%  for the reason rangeSolve leaves it: only the caller knows the vocabulary.
%  A targeting script wants to say "aim point 4200 km downrange and 180 km
%  cross-track, best reachable miss 11 km, widen the loft bracket", which
%  needs the units and the problem's words. A library function that threw
%  here would force every caller into a try/catch to read numbers this struct
%  already carries. What must never happen is a silent near miss, and
%  converged = false is what prevents it.
%
%  The six refusals, all with converged = false:
%      coorbital:aimSolve:toleranceNotMet    opts.maxIter iterations went by
%                                            without max(abs(f)) reaching tolF
%      coorbital:aimSolve:singularJacobian   the Jacobian's condition number
%                                            exceeded opts.condMax
%      coorbital:aimSolve:stepStalled        opts.maxHalve halvings of the
%                                            Newton step all failed to reduce
%                                            max(abs(f))
%      coorbital:aimSolve:jacobianFailed     fResid would not evaluate at a
%                                            finite-difference probe point.
%                                            Remedy: SHRINK dx0
%      coorbital:aimSolve:stepUnderflow      a control has drifted to a
%                                            magnitude at which adding its
%                                            dx0 changes nothing, so no
%                                            derivative exists there.
%                                            Remedy: ENLARGE dx0, or bound
%                                            the controls with opts.maxStep
%      coorbital:aimSolve:jacobianNotFinite  both residuals were finite but
%                                            their difference quotient
%                                            overflowed. Remedy: rescale
%                                            fResid
%  The last three are kept apart rather than merged into one "the Jacobian
%  could not be built" because their remedies are not the same, and two of
%  them are exact opposites. A single message telling a caller to shrink dx0
%  when the fix is to enlarge it is worse than no message.
%
%  CONVERGENCE IS ON THE RESIDUAL, NEVER ON THE STEP. max(abs(f)) <= tolF is
%  the only success criterion. A small step is emphatically NOT evidence of a
%  solution here: this method damps, and a damped step is small precisely
%  when the line search is in trouble. A solver that also accepted a small
%  step would report its worst failures as successes.
%
%  THE JACOBIAN IS GUARDED BEFORE IT IS USED. Near a maximum-range control
%  the downrange sensitivity dR/dcontrol collapses, the two columns of the
%  Jacobian become parallel, and the Newton step is a wild extrapolation
%  along a direction the residuals barely respond to. The 2-norm condition
%  number is formed from the singular values, sigma_max/sigma_min, and
%  compared against opts.condMax BEFORE any step is computed, so the caller
%  is told the number rather than handed an infinity. The step itself is then
%  SOLVED THROUGH THE SAME THREE FACTORS -- rotate the residual into the
%  singular basis, divide by the singular values, rotate back -- which costs
%  no second factorisation and cannot emit the near-singularity warning a
%  backslash would print to the console. No inverse is formed at any point.
%  This is the single most common way a two-axis targeting solve goes wrong,
%  and it is a geometry statement about the trajectory, not a numerical
%  accident.
%
%  Because that test bounds the smallest singular value below by
%  sigma_max/condMax, it also bounds the step, and there is deliberately no
%  separate finiteness test on the step itself: such a test would be
%  unreachable, and unreachable code is worse than absent code because it
%  reads as protection that was never exercised. A Jacobian whose columns are
%  parallel to the last bit gives sigma_min exactly zero, an infinite
%  condition number, and the singularJacobian refusal.
%
%  THE STEP IS BOUNDED TWICE. opts.maxStep caps the change in each control
%  per iteration, applied as ONE scale factor on the whole step so the Newton
%  DIRECTION survives -- clipping the components independently would rotate
%  the step into a direction nobody computed. Underneath that, the
%  backtracking line search halves until max(abs(f)) strictly decreases.
%  Along a Newton direction the linearised residual is (1 - lambda)*f, so the
%  same infinity norm that defines convergence is a proper merit function for
%  the line search and no second measure is needed.
%
%  A FAILED PROPAGATION IS NOT A FAILED SOLVE. fResid is allowed to throw, or
%  to return a non-finite pair, for a control pair whose flight does not
%  complete -- a coordinate guard, a vehicle that never reaches the glide
%  phase, an integrator that gives up. At a TRIAL point that is caught,
%  counted in both nShrink and nInfeasible, and treated exactly like a step
%  that failed to improve: halve and try again. At a finite-difference PROBE
%  point there is nothing to halve, so the solve refuses with jacobianFailed,
%  counting the failure in nInfeasible only, and the message quotes the
%  underlying error text -- which is quoted on THAT path alone, so a caller
%  can never be shown an unrelated earlier failure as the cause of this one.
%
%  THE INITIAL POINT IS THE ONE EXCEPTION, AND IT THROWS. A non-finite
%  residual at x0 raises coorbital:aimSolve:nonFiniteEval, and an fResid that
%  throws at x0 has its own error propagate untouched. There is no step to
%  shrink and no previous point to fall back on, so this is not a failure to
%  converge -- it is an unusable starting guess, which is an input error, and
%  it is reported the way rangeSolve reports its nonFiniteEval.
%
%  dx0 IS A TWO-VECTOR AND A SCALAR IS REFUSED. The two controls do not share
%  a scale: a heading is radians of order one, while the range control may be
%  radians of order a tenth or a burn fraction in [0,1]. A single difference
%  step therefore cannot be right for both -- it is either lost in the
%  propagation's own noise on one axis or a chord across real curvature on
%  the other -- and accepting a scalar would let that mistake pass silently.
%  The default, 1e-4.*max(abs(x0),1) elementwise, is relative to each control
%  with an absolute floor, and it is deliberately far coarser than the
%  sqrt(eps) step appropriate to an analytic function: a residual produced by
%  integrating a trajectory is smooth only down to a few metres, so a tiny
%  difference step measures the integrator's noise instead of the derivative.
%  Callers who know their own scales should say so.
%
%  Genuinely malformed input throws:
%      coorbital:aimSolve:badResidual    fResid not a function handle
%      coorbital:aimSolve:badGuess       x0 not a finite real two-vector
%      coorbital:aimSolve:badStep        dx0 not a finite real positive
%                                        two-vector, a scalar included
%      coorbital:aimSolve:badTolerance   tolF non-finite or not positive
%      coorbital:aimSolve:badOptions     an unrecognised or invalid opts field
%      coorbital:aimSolve:nonFiniteEval  fResid returned a non-finite, complex
%                                        or wrongly sized value at x0
%  The tolerance and the residual checks matter for the same reason at the
%  two ends of the comparison: every test in this solver is an inequality,
%  and every inequality against NaN is false. A NaN loose in the loop would
%  neither converge nor trip a refusal.
%
%  The difference quotient divides by the REALISED step, xProbe(k) - x(k),
%  not by the requested dx0(k). For a control of order one and a step of
%  order 1e-4 the two differ in the last bits, and using the requested value
%  puts that discrepancy straight into the Jacobian.
%
%% References:
%   [1] Dennis, J.E. and Schnabel, R.B., "Numerical Methods for Unconstrained
%       Optimization and Nonlinear Equations," SIAM Classics in Applied
%       Mathematics 16, 1996, Chapters 5 and 6. Finite-difference Jacobians,
%       and the line-search damped Newton method for systems of nonlinear
%       equations: full step first, backtrack until the merit function falls.
%   [2] Press, W.H., et al., "Numerical Recipes," 3rd ed., Cambridge, 2007,
%       Section 9.7. Globally convergent Newton-Raphson for nonlinear systems.
%   [3] Golub, G.H. and Van Loan, C.F., "Matrix Computations," 4th ed., Johns
%       Hopkins, 2013, Section 2.6.2. The 2-norm condition number as the
%       ratio of extreme singular values.
%
%% Revision History:
%  Michael Casey                                                08/08/2026
%  Copyright 2026 Coorbital, Inc.
%% ------------------------ Begin Code Sequence ---------------------------

%% Self-demo. A toy two-axis targeting law standing in for a trajectory:
%% downrange saturates in the range control the way a real range-versus-loft
%% curve does, and cross-track is mostly heading error with a little coupling
%% into the loft. Both the solve and the refusal are shown, because the
%% refusal is the path a user meets first:
if nargin == 0
                cM = coorbital.util.missileConst();
             fDemo = @(x) [5.0e6.*sin(2.*x(2)).*cos(x(1) - 0.6) - 3.0e6; ...
                           cM.rE.*(x(1) - 0.6) + 3.0e5.*(x(2) - 0.35)];
    [xD,fD,infoD] = coorbital.util.aimSolve(fDemo,[0.55;0.30], ...
                    [1.0e-5;1.0e-5],1.0);
    fprintf('aim point reached: psi %.6f rad, control %.6f rad\n',xD(1),xD(2));
    fprintf('  miss %.4f m downrange, %.4f m cross-track, %d iteration(s), %d eval(s)\n', ...
            fD(1),fD(2),infoD.iterations,infoD.nEval);
              fBad = @(x) [5.0e6.*sin(2.*x(2)).*cos(x(1) - 0.6) - 9.0e6; ...
                          cM.rE.*(x(1) - 0.6) + 3.0e5.*(x(2) - 0.35)];
       [~,~,iBad] = coorbital.util.aimSolve(fBad,[0.55;0.30], ...
                    [1.0e-5;1.0e-5],1.0);
    fprintf('downrange beyond the envelope: converged = %d\n',iBad.converged);
    fprintf('  %s\n',iBad.message);
              xSol = [];
              fAch = [];
              info = [];
    return;
end

%% Fill the optional trailing arguments so that every check below sees a
%% value. An empty is treated as "not supplied", which lets a caller reach
%% opts without inventing a dx0:
if nargin < 3
               dx0 = [];
end
if nargin < 4
              tolF = [];
end
if nargin < 5
              opts = struct();
end

%% Validate the residual handle. Everything downstream calls it, and the
%% failure a non-handle produces otherwise is reported from inside a nested
%% function where the caller cannot see what went wrong:
if ~isa(fResid,'function_handle')
    error('coorbital:aimSolve:badResidual', ...
        ['fResid must be a function handle taking a [2 x 1] control ', ...
         'vector and returning a [2 x 1] residual; got a %s.'], ...
        class(fResid));
end

%% Validate the initial guess. A NaN here is invisible everywhere else:
%% every convergence and line-search test below is an inequality, and every
%% inequality against NaN is false, so the solve would neither converge nor
%% refuse -- it would run to the iteration cap and report a tolerance it
%% could not meet, as though the caller's fResid were at fault:
if ~isnumeric(x0) || ~isreal(x0) || numel(x0) ~= 2 || ~all(isfinite(x0(:)))
    error('coorbital:aimSolve:badGuess', ...
        ['x0 must be a finite real two-vector of controls; got %s. The ', ...
         'two entries are the heading and the range control, in that ', ...
         'order.'],mat2str(x0));
end
                x0 = x0(:);

%% Validate the difference step, and REFUSE A SCALAR outright. The two
%% controls are on unrelated scales, so one number cannot be the right chord
%% for both, and silently broadcasting it would hide that in a Jacobian
%% column nobody looks at:
if isempty(dx0)
               dx0 = 1.0e-4.*max(abs(x0),1);
end
if isscalar(dx0)
    error('coorbital:aimSolve:badStep', ...
        ['dx0 must be a two-vector, one difference step per control; got ', ...
         'the scalar %.10g. A heading in radians and a range control that ', ...
         'may be an angle or a fraction do not share a scale, so one step ', ...
         'cannot serve both.'],dx0);
end
if ~isnumeric(dx0) || ~isreal(dx0) || numel(dx0) ~= 2 || ...
        ~all(isfinite(dx0(:))) || any(dx0(:) <= 0)
    error('coorbital:aimSolve:badStep', ...
        ['dx0 must be a finite real two-vector of strictly positive ', ...
         'difference steps; got %s.'],mat2str(dx0));
end
               dx0 = dx0(:);

%% Validate the tolerance. Zero or negative can never be met and would spin
%% the solve to its iteration cap for nothing:
if isempty(tolF)
              tolF = 1.0;
end
if ~isscalar(tolF) || ~isreal(tolF) || ~isfinite(tolF) || tolF <= 0
    error('coorbital:aimSolve:badTolerance', ...
        ['tolF must be a finite positive scalar residual tolerance; got ', ...
         '%s.'],mat2str(tolF));
end

%% Merge the safeguards. An unrecognised field is an ERROR rather than a
%% silent no-op: a caller who writes opts.maxIterations and watches the
%% solver ignore it has no way to tell that from a solver that honoured it:
if ~isstruct(opts) || ~isscalar(opts)
    error('coorbital:aimSolve:badOptions', ...
        'opts must be a scalar struct; got a %s.',class(opts));
end
              dflt = struct('maxIter',25,'maxStep',[Inf;Inf], ...
                            'maxHalve',12,'condMax',1.0e8);
             known = fieldnames(dflt);
             given = fieldnames(opts);
for kOpt = 1:numel(given)
             fName = given{kOpt};
    if ~any(strcmp(fName,known))
        error('coorbital:aimSolve:badOptions', ...
            'Unrecognised option ''%s''. The options are %s.', ...
            fName,strjoin(known',', '));
    end
    if ~isempty(opts.(fName))
      dflt.(fName) = opts.(fName);
    end
end
              opts = dflt;

%% Each safeguard checked on its own, so the message names the one at fault:
if ~isscalar(opts.maxIter) || ~isfinite(opts.maxIter) || opts.maxIter < 1 || ...
        opts.maxIter ~= fix(opts.maxIter)
    error('coorbital:aimSolve:badOptions', ...
        'opts.maxIter must be a positive whole number; got %s.', ...
        mat2str(opts.maxIter));
end
if ~isscalar(opts.maxHalve) || ~isfinite(opts.maxHalve) || ...
        opts.maxHalve < 0 || opts.maxHalve ~= fix(opts.maxHalve)
    error('coorbital:aimSolve:badOptions', ...
        'opts.maxHalve must be a non-negative whole number; got %s.', ...
        mat2str(opts.maxHalve));
end
if ~isscalar(opts.condMax) || ~isreal(opts.condMax) || opts.condMax <= 1
    error('coorbital:aimSolve:badOptions', ...
        ['opts.condMax must be a real scalar greater than one; got %s. ', ...
         'The condition number of any matrix is at least one.'], ...
        mat2str(opts.condMax));
end
if isscalar(opts.maxStep)
      opts.maxStep = [opts.maxStep;opts.maxStep];
end
if ~isnumeric(opts.maxStep) || ~isreal(opts.maxStep) || ...
        numel(opts.maxStep) ~= 2 || any(opts.maxStep(:) <= 0) || ...
        any(isnan(opts.maxStep(:)))
    error('coorbital:aimSolve:badOptions', ...
        ['opts.maxStep must be one or two strictly positive real step ', ...
         'bounds, Inf permitted; got %s.'],mat2str(opts.maxStep));
end
      opts.maxStep = opts.maxStep(:);

%% Evaluate the initial guess. This is the ONE place a bad residual throws:
%% there is no step to shrink and no earlier point to fall back on, so an
%% unusable x0 is an input error and not a failure to converge:
             nEval = 0;
       nInfeasible = 0;
           nShrink = 0;
           msgLast = '';
              fCur = evalStrict(x0);
              xCur = x0;
            errCur = max(abs(fCur));

%% Preallocate the histories to the iteration cap and trim on the way out;
%% nothing here is grown inside the loop:
           maxIter = opts.maxIter;
             xHist = zeros(2,maxIter + 1);
             fHist = zeros(2,maxIter + 1);
           jacHist = zeros(2,2,maxIter);
          condHist = zeros(1,maxIter);
           lamHist = zeros(1,maxIter);
        xHist(:,1) = xCur;
        fHist(:,1) = fCur;

%% BEST-SO-FAR, seeded from the initial guess. The line search only accepts a
%% strict decrease in max(abs(f)), so the current iterate IS the best point
%% on the path -- but the output contract promises the best, so it is tracked
%% explicitly rather than inferred, and x0 stands as the answer if the very
%% first step cannot be taken:
              xSol = xCur;
              fAch = fCur;
           bestErr = errCur;

%% Counters kept apart because the exits differ in how far into an iteration
%% they get: nJac Jacobians were built, nLS line searches were run, nAcc
%% steps were accepted, and on a refusal these need not be equal:
             kIter = 0;
              nJac = 0;
               nLS = 0;
              nAcc = 0;
         converged = false;
             ident = '';

while true

%% Converged on the ACHIEVED residual, tested BEFORE any work is done, so a
%% guess that is already good costs one evaluation and nothing else. This is
%% the only success criterion: a small step proves nothing, because damping
%% makes the step small exactly when the solve is struggling:
    if errCur <= tolF
         converged = true;
        break;
    end
    if kIter >= maxIter
             ident = 'coorbital:aimSolve:toleranceNotMet';
        break;
    end
             kIter = kIter + 1;

%% Forward-difference Jacobian, two evaluations. The base value fCur is
%% reused, never recomputed -- that is the whole reason an iteration costs
%% three propagations and not four -- and the quotient divides by the
%% REALISED step so that the rounding in xCur + dx0 does not leak into the
%% derivative.
%%
%% THREE DISTINCT THINGS CAN GO WRONG HERE, and they are kept apart because
%% they call for opposite remedies. The realised step can vanish, which means
%% dx0 is too SMALL against the magnitude the control has drifted to and must
%% be enlarged. The probe point can fail to fly, which is fResid's business
%% and usually means dx0 is too LARGE or the guess is badly placed. The
%% quotient can overflow, which is a scaling problem in the residual itself.
%% Reporting all three as one refusal would hand the caller a remedy that is
%% sometimes exactly backwards:
              jMat = zeros(2,2);
             jFail = '';
              kBad = 0;
    for kCol = 1:2
            xProbe = xCur;
      xProbe(kCol) = xCur(kCol) + dx0(kCol);
             hReal = xProbe(kCol) - xCur(kCol);
              kBad = kCol;

%% Tested BEFORE the probe is evaluated, because a probe at a point that is
%% bit-for-bit the base point buys nothing but a propagation:
        if hReal == 0
             jFail = 'underflow';
            break;
        end
      [fProbe,okP] = evalSoft(xProbe);
        if ~okP
             jFail = 'probe';
            break;
        end
      jMat(:,kCol) = (fProbe - fCur)./hReal;
        if ~all(isfinite(jMat(:,kCol)))
             jFail = 'overflow';
            break;
        end
    end

%% None of the three leaves a Jacobian to guard or a step to halve, so each
%% refuses rather than shrinking -- but each with its own identifier, so the
%% caller is not told to shrink dx0 when the fix is to enlarge it:
    if strcmp(jFail,'underflow')
             ident = 'coorbital:aimSolve:stepUnderflow';
        break;
    elseif strcmp(jFail,'probe')
             ident = 'coorbital:aimSolve:jacobianFailed';
        break;
    elseif strcmp(jFail,'overflow')
             ident = 'coorbital:aimSolve:jacobianNotFinite';
        break;
    end

%% Factor the Jacobian ONCE, by singular value decomposition, and read both
%% the conditioning and the step out of the same factors. The 2-norm
%% condition number is the ratio of the singular values, formed and tested
%% BEFORE any step is computed, so an ill-conditioned Jacobian comes back to
%% the caller as a number it can print rather than as an infinity in the
%% step. Factoring here rather than leaving the step to a backslash also
%% keeps the solve silent: a Jacobian that LU decides is singular to working
%% precision would warn to the console, which a library function called
%% thousands of times inside a targeting loop must never do:
  [uJac,sJac,vJac] = svd(jMat);
              sVal = diag(sJac);
    if sVal(2) <= 0
             condJ = Inf;
    else
             condJ = sVal(1)./sVal(2);
    end
              nJac = nJac + 1;
 jacHist(:,:,nJac) = jMat;
    condHist(nJac) = condJ;
    if ~isfinite(condJ) || condJ > opts.condMax
             ident = 'coorbital:aimSolve:singularJacobian';
        break;
    end

%% The Newton direction, solved through the same three factors: rotate the
%% residual into the singular basis, divide by the singular values, rotate
%% back. No inverse is formed anywhere.
%%
%% There is deliberately no second finiteness test on the step. The two tests
%% above already bound it: the Jacobian entries are finite, and the condition
%% test bounds sVal(2) below by sVal(1)/condMax, so the division cannot
%% overflow unless the residual itself is within a factor of condMax of the
%% floating-point ceiling. A test for that would be unreachable code, and
%% unreachable code is worse than absent code because it reads as protection.
%% If it ever did happen the solve still refuses rather than misbehaves: an
%% infinite trial point is simply an infeasible one, and the line search
%% below shrinks and then reports stepStalled:
              sDir = -(vJac*((uJac'*fCur)./sVal));

%% Cap the step with ONE scale factor rather than clipping the components
%% independently: clipping would rotate the step off the Newton direction and
%% into a direction nobody computed, and the line search below assumes the
%% direction it is shortening is the Newton one:
              capF = min([1;opts.maxStep./abs(sDir)]);
              sDir = capF.*sDir;

%% Backtracking line search. Along a Newton direction the linearised residual
%% is (1 - lambda)*f, so max(abs(f)) falls for any small enough lambda and
%% serves as the merit function; using the same measure the tolerance uses
%% means the solve cannot accept a step that moves it away from converging.
%% A trial point that throws or returns a non-finite pair is not an error, it
%% is an infeasible control pair, and it shrinks the step exactly as a step
%% that merely failed to improve would:
               lam = 1;
          accepted = false;
            errTry = errCur;
              fTry = fCur;
              xTry = xCur;
    for kHalve = 0:opts.maxHalve
              xTry = xCur + lam.*sDir;
        [fTry,okT] = evalSoft(xTry);
        if okT
            errTry = max(abs(fTry));
            if errTry < bestErr
           bestErr = errTry;
              xSol = xTry;
              fAch = fTry;
            end
            if errTry < errCur
          accepted = true;
                break;
            end
        end
           nShrink = nShrink + 1;
               lam = lam./2;
    end

%% Undo the halving the rejected last attempt left behind, so that lambda
%% always names a factor that was actually TRIED. Reporting 2^-(maxHalve+1)
%% here would name a step the solve never took:
    if ~accepted
               lam = lam.*2;
    end
               nLS = nLS + 1;
      lamHist(nLS) = lam;
    if ~accepted
             ident = 'coorbital:aimSolve:stepStalled';
        break;
    end

%% Accept, and carry the trial residual forward as the next base value so
%% that the new iterate is never re-evaluated:
              xCur = xTry;
              fCur = fTry;
            errCur = errTry;
              nAcc = nAcc + 1;
 xHist(:,nAcc + 1) = xCur;
 fHist(:,nAcc + 1) = fCur;
end

%% Trim the histories to what was actually filled, so a caller reading
%% info.jacobians cannot pick up a preallocated zero and mistake it for a
%% Jacobian that was never built:
           info.x0 = x0;
          info.dx0 = dx0;
         info.tolF = tolF;
         info.opts = opts;
     info.residual = fAch;
       info.errMax = max(abs(fAch));
   info.iterations = kIter;
        info.condJ = condHist(1:nJac);
    info.jacobians = jacHist(:,:,1:nJac);
       info.lambda = lamHist(1:nLS);
        info.xHist = xHist(:,1:(nAcc + 1));
        info.fHist = fHist(:,1:(nAcc + 1));
        info.nEval = nEval;
      info.nShrink = nShrink;
  info.nInfeasible = nInfeasible;
    info.converged = converged;
   info.identifier = ident;

%% One sentence per exit, each naming the numbers a caller would otherwise
%% have to dig out of the struct to write its own refusal:
if converged
      info.message = sprintf(['Converged in %d iteration(s) and %d ', ...
                     'evaluation(s) at x = [%.10g, %.10g], achieving ', ...
                     'residuals of [%.6g, %.6g] against a tolerance of ', ...
                     '%.6g.'],kIter,nEval,xSol(1),xSol(2),fAch(1),fAch(2), ...
                     tolF);
elseif strcmp(ident,'coorbital:aimSolve:toleranceNotMet')
      info.message = sprintf(['Failed to reach a residual tolerance of ', ...
                     '%.6g in %d iteration(s) and %d evaluation(s). The ', ...
                     'closest point visited was x = [%.10g, %.10g], ', ...
                     'achieving residuals of [%.6g, %.6g], against ', ...
                     '[%.6g, %.6g] at the initial guess. Raise ', ...
                     'opts.maxIter, loosen tolF, or start nearer.'], ...
                     tolF,kIter,nEval,xSol(1),xSol(2),fAch(1),fAch(2), ...
                     fHist(1,1),fHist(2,1));
elseif strcmp(ident,'coorbital:aimSolve:singularJacobian')
      info.message = sprintf(['The finite-difference Jacobian is ', ...
                     'ill-conditioned at x = [%.10g, %.10g] after %d ', ...
                     'iteration(s): condition number %.6g against a ', ...
                     'ceiling of %.6g. The two controls have stopped ', ...
                     'moving the two residuals independently, which is ', ...
                     'what a range maximum looks like -- the downrange ', ...
                     'sensitivity collapses there -- so no Newton step is ', ...
                     'trustworthy. The best point visited was ', ...
                     'x = [%.10g, %.10g], achieving [%.6g, %.6g]. Move ', ...
                     'the aim point inside the reachable envelope or ', ...
                     'restart away from the maximum.'],xCur(1),xCur(2), ...
                     kIter,condHist(nJac),opts.condMax,xSol(1),xSol(2), ...
                     fAch(1),fAch(2));
elseif strcmp(ident,'coorbital:aimSolve:stepStalled')
      info.message = sprintf(['The damped Newton step stalled at ', ...
                     'x = [%.10g, %.10g] after %d iteration(s): %d trial ', ...
                     'step(s), the shortest a factor %.6g of the full ', ...
                     'Newton step, all failed to reduce the largest ', ...
                     'residual below %.6g. Over the whole solve %d trial ', ...
                     'point(s) would not fly at all. Achieved residuals ', ...
                     'there are [%.6g, %.6g] against a tolerance of %.6g. ', ...
                     'This is a near miss being REFUSED rather than ', ...
                     'returned: raise opts.maxHalve if the direction is ', ...
                     'right and only the length is wrong, otherwise the ', ...
                     'aim point is out of reach from this guess.'], ...
                     xSol(1),xSol(2),kIter,opts.maxHalve + 1,lam,errCur, ...
                     nInfeasible,fAch(1),fAch(2),tolF);
elseif strcmp(ident,'coorbital:aimSolve:stepUnderflow')
      info.message = sprintf(['The difference step underflowed at ', ...
                     'x = [%.10g, %.10g] after %d iteration(s): control %d ', ...
                     'has reached a magnitude at which adding its step of ', ...
                     '%.6g leaves it bit-for-bit unchanged, so no ', ...
                     'derivative can be measured there. ENLARGE dx0(%d) -- ', ...
                     'shrinking it makes this worse -- or bound the ', ...
                     'controls with opts.maxStep so they cannot run out to ', ...
                     'this magnitude in the first place. The best point ', ...
                     'visited was x = [%.10g, %.10g], achieving ', ...
                     '[%.6g, %.6g].'],xCur(1),xCur(2),kIter,kBad, ...
                     dx0(kBad),kBad,xSol(1),xSol(2),fAch(1),fAch(2));
elseif strcmp(ident,'coorbital:aimSolve:jacobianNotFinite')
      info.message = sprintf(['The difference quotient overflowed at ', ...
                     'x = [%.10g, %.10g] after %d iteration(s). fResid ', ...
                     'returned finite residuals at both the base point and ', ...
                     'the probe beside control %d, but they differ by more ', ...
                     'than the floating-point range spans once divided by ', ...
                     'the step of %.6g, so the Jacobian column is not ', ...
                     'finite. Rescale fResid, or shrink dx0(%d) so the ', ...
                     'probe does not straddle whatever the residual is ', ...
                     'doing there.'],xCur(1),xCur(2),kIter,kBad,dx0(kBad), ...
                     kBad);
else
      info.message = sprintf(['fResid would not evaluate at the ', ...
                     'finite-difference probe point beside control %d at ', ...
                     'x = [%.10g, %.10g] after %d iteration(s), so the ', ...
                     'Jacobian cannot be built there. The difference ', ...
                     'steps were [%.6g, %.6g]. Shrink dx0, or start from ', ...
                     'a guess whose neighbourhood flies. The underlying ', ...
                     'failure was: %s'],kBad,xCur(1),xCur(2),kIter,dx0(1), ...
                     dx0(2),msgLast);
end

    function fOut = evalStrict(xIn)
    %% Purpose:
    %
    %  Evaluate the caller's residual at the INITIAL point, tally the call,
    %  and refuse anything that is not a finite real pair. Nested so that the
    %  tally lives in one place and info.nEval cannot drift from the number
    %  of calls actually made. Strict rather than soft because there is no
    %  step to shrink at x0.
    %
    %% Inputs:
    %
    %  xIn              [2 x 1]                     Controls, in the
    %                                               controls' units
    %
    %% Outputs:
    %
    %  fOut             [2 x 1]                     fResid(xIn), in fResid's
    %                                               units
    %
    %% Revision History:
    %  Michael Casey                                                08/08/2026
    %  Copyright 2026 Coorbital, Inc.
    %% ------------------------ Begin Code Sequence ---------------------------
             nEval = nEval + 1;
              fRaw = fResid(xIn);
        if ~isnumeric(fRaw) || ~isreal(fRaw) || numel(fRaw) ~= 2 || ...
                ~all(isfinite(fRaw(:)))
            error('coorbital:aimSolve:nonFiniteEval', ...
                ['fResid returned %s at the initial guess x = [%.10g, ', ...
                 '%.10g]. It must return a finite real [2 x 1] residual ', ...
                 'there; a propagation that fails at the starting point ', ...
                 'is an unusable guess, not a failure to converge.'], ...
                mat2str(fRaw),xIn(1),xIn(2));
        end
              fOut = fRaw(:);
    end

    function [fOut,ok] = evalSoft(xIn)
    %% Purpose:
    %
    %  Evaluate the caller's residual at a TRIAL or PROBE point, tally the
    %  call, and report failure as a flag rather than an exception. A control
    %  pair whose flight does not complete is an ordinary event in targeting
    %  -- the propagation hits a coordinate guard, or the vehicle never
    %  reaches the phase the residual is measured in -- and the solve
    %  responds by shrinking the step, not by giving up. Both a thrown error
    %  and a non-finite return are treated the same way, because they mean
    %  the same thing.
    %
    %% Inputs:
    %
    %  xIn              [2 x 1]                     Controls, in the
    %                                               controls' units
    %
    %% Outputs:
    %
    %  fOut             [2 x 1]                     fResid(xIn) in fResid's
    %                                               units when ok, otherwise
    %                                               [NaN;NaN]
    %
    %  ok               [1 x 1] logical             True when the point flew
    %                                               and returned a finite
    %                                               real pair
    %
    %% Revision History:
    %  Michael Casey                                                08/08/2026
    %  Copyright 2026 Coorbital, Inc.
    %% ------------------------ Begin Code Sequence ---------------------------
                ok = false;
              fOut = [NaN;NaN];
             nEval = nEval + 1;
        try
              fRaw = fResid(xIn);
        catch errEval
       nInfeasible = nInfeasible + 1;
           msgLast = errEval.message;
            return;
        end
        if ~isnumeric(fRaw) || ~isreal(fRaw) || numel(fRaw) ~= 2 || ...
                ~all(isfinite(fRaw(:)))
       nInfeasible = nInfeasible + 1;
           msgLast = sprintf('fResid returned %s rather than a finite real pair', ...
                     mat2str(fRaw));
            return;
        end
              fOut = fRaw(:);
                ok = true;
    end
end
