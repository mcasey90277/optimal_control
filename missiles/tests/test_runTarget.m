function test_runTarget()
%% Purpose:
%
%  Pin HGV/run_target, the point-to-point targeting script: the closed-form
%  launch azimuth, the bisection on thrust-termination time, the separation
%  link that throws the whole booster away, the reachable-envelope refusal,
%  and the summary that reports all of it. Every other entry-script test in
%  this suite checks a trajectory that was flown FORWARD from a given azimuth;
%  this one checks a trajectory that was SOLVED to arrive somewhere, which is
%  a different claim and fails in different ways.
%
%  THE THREE THINGS THIS TEST EXISTS TO CATCH, and how each is caught:
%
%    A wrong azimuth.    The impact point walks off the launch-to-target great
%                        circle while the range solve still converges, so the
%                        RANGE residual stays small and the vehicle still
%                        misses. Caught by measuring the miss impact-to-target
%                        with coorbital.util.greatCircle, from the flown state,
%                        rather than by reading the solver's own residual.
%
%    A range solve that  Caught by the same measured miss: a cutoff returned
%    did not converge.   without iterating lands thousands of kilometres away.
%
%    A staging slip.     The separation link is what makes the post-cutoff
%                        vehicle the payload. Drop it and the carried mass is
%                        the whole stack while the equations of motion still
%                        divide by 900 kg. coorbital.eom.massConstant raises
%                        on that, so the mutation shows up as a thrown error
%                        rather than a wrong number -- which is the point of
%                        the guard.
%
%  WHY THE FLOWN CASE IS NEITHER EQUATORIAL NOR DUE EAST. From (0,0) the
%  central angle reduces to acos(cos(lat2) cos(lon2)), which is symmetric in
%  its two arguments, so a due-east equatorial case is provably BLIND to a
%  latitude/longitude transposition at every great-circle call site. The
%  shipped configuration is 20 N 155 W to 35 N 120 W on a 56.6 deg azimuth,
%  where all four coordinates are distinct; part 2 below proves that geometry
%  is discriminating rather than assuming it.
%
%  COST. A trajectory propagation here is about 0.15 s and the shipped solve
%  takes twelve of them, so the whole test runs in a few seconds and there was
%  no need to loosen the tolerance or shrink the bracket for it. The SHIPPED
%  configuration is therefore flown at its shipped settings, which is the
%  stronger test: what is pinned is exactly what a user gets. The second,
%  range-tracking case does use a narrower bracket and a 5 km tolerance,
%  because nothing there needs sub-kilometre accuracy -- only that the solve
%  moves when the target moves.
%
%  All reference numbers below were read out of a run of the committed code at
%  full printed precision, not rounded from a report.
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

                 c = coorbital.util.missileConst();

%% Figures open before anything runs. Nothing in this file asks for a plot, so
%% nothing in this file may leave one behind -- but the count is taken as a
%% DIFFERENCE, because another test in the suite is entitled to have left one
%% and that is not this file's failure to report:
            nFig0  = numel(findall(groot,'Type','figure'));

%% The shipped launch point and target, restated here in the test's own hand.
%% They are NOT read out of run_target: a test that took its reference from
%% the code under test could not see the code change them:
            latLnD = 20;
            lonLnD = -155;
            latTgD = 35;
            lonTgD = -120;
             latLn = deg2rad(latLnD);
             lonLn = deg2rad(lonLnD);
             latTg = deg2rad(latTgD);
             lonTg = deg2rad(lonTgD);

%% ---------------------------------------------------------------------
%% 1. The shipped configuration
%% ---------------------------------------------------------------------
%% Plots off so the suite creates no figures, and the summary captured rather
%% than printed so the suite output stays readable. The captured text is not
%% discarded -- part 4 asserts against it, because the summary is what a
%% reader of this library actually acts on:
            outDef = evalc(['[trDef,inDef] = ' ...
                            'run_target(struct(''showPlots'',false));']);

%% Termination was nominal, and no caution or cross-range warning fired:
    assert(contains(outDef,'(nominal)'), ...
        'the shipped run did not terminate nominally. Summary was:\n%s',outDef);
    assert(~contains(outDef,'CAUTION'), ...
        'the shipped run printed the truncation caution. Summary was:\n%s',outDef);
    assert(~contains(outDef,'WARNING'), ...
        'the shipped run printed a warning. Summary was:\n%s',outDef);
    assert(~contains(outDef,'REFUSED'), ...
        'the shipped target must be reachable. Summary was:\n%s',outDef);
    assert(~inDef.refused,'info.refused must be false on a converged solve');
    assert(~inDef.crossWarn,'info.crossWarn must be false at zero bank');
    assert(inDef.stopOK,'every phase must have ended as intended');

%% Terminal altitude is the configured stop altitude. eventAltitude is a
%% terminal one-sided event, so this is exact to the event solver:
             hEndM = trDef.x(end,1) - c.rE;
    assertAbs(hEndM,0,1e-6,'terminal altitude (m)');

%% THE MISS IS INSIDE THE REQUESTED TOLERANCE. Asserted FIRST, ahead of every
%% pinned literal, because it is the claim the whole script exists to make and
%% because a pinned literal that happens to sit near its budget would otherwise
%% fire first and report a 1e-4 drift where the real fault is a 200 km miss.
%% The shipped tolerance is 1 km, restated here rather than read from the code:
              tolM = 1000;
    assert(inDef.missM <= tolM, ...
        ['the solve converged but missed by %.2f m against a %.1f m ' ...
         'tolerance'],inDef.missM,tolM);

%% THE HEADLINE NUMBERS, to 1e-4 relative, read at full printed precision from
%% a run of the committed code:
    assertRel(inDef.tBurn    ,80.517757894736832  ,1e-4,'full burn time (s)');
    assertRel(inDef.rngReqM  ,3811239.9245379278  ,1e-4,'required range (m)');
    assertRel(inDef.rngAchM  ,3811751.1679983996  ,1e-4,'achieved range (m)');
    assertRel(inDef.missM    ,511.24346047081326  ,1e-4,'miss distance (m)');
    assertRel(inDef.residM   ,511.24346047174186  ,1e-4,'range residual (m)');
    assertRel(inDef.tCut     ,75.681974583675967  ,1e-4,'solved cutoff (s)');
    assertRel(inDef.cutFrac  ,0.93994140624999989 ,1e-4,'cutoff as a fraction of burn');
    assertRel(inDef.psiLaunch,0.98833115856379961 ,1e-4,'launch azimuth (rad)');
    assertRel(inDef.rngMinM  ,202861.23979221776  ,1e-4,'reachable minimum (m)');
    assertRel(inDef.rngMaxM  ,7737630.3673877586  ,1e-4,'reachable maximum (m)');
    assertRel(inDef.tFlight  ,1626.0865606304174  ,1e-4,'flight time (s)');
    assertRel(inDef.vImpact  ,229.2361700377173   ,1e-4,'impact speed (m/s)');
    assertRel(inDef.hGlideMax,118.0857834646022   ,1e-4,'peak glide altitude (km)');
    assertRel(inDef.nBstMax  ,22.748457688111941  ,1e-4,'peak boost sensed load (g)');
    assertRel(inDef.nAerMax  ,7.9220945536800578  ,1e-4,'peak aero load (g)');
    assertRel(inDef.mpBurned ,28198.242187499993  ,1e-4,'propellant burned (kg)');
    assertRel(inDef.mpWasted ,1801.7578125000073  ,1e-4,'propellant jettisoned (kg)');

%% The impact point itself, pinned separately from the range. A single scalar
%% range is built from both coordinates and cannot tell them apart:
    assertRel(rad2deg(inDef.latImpact), 35.00131692011076,1e-4,'impact latitude (deg)');
    assertRel(rad2deg(inDef.lonImpact),-119.9946288960424,1e-4,'impact longitude (deg)');

%% Bisection is deterministic: the same bracket and the same tolerance over
%% the same range function take the same number of halvings every time. A
%% change here means the range function moved, which the pinned ranges above
%% should already have caught -- or that the solver stopped iterating, which
%% they might not:
    assert(inDef.iterations == 10, ...
        'bisection took %d steps against the expected 10',inDef.iterations);
    assert(inDef.nEval == inDef.iterations + 2, ...
        ['%d evaluations for %d iterations; rangeSolve must cost exactly two ' ...
         'endpoint evaluations plus one per step'],inDef.nEval,inDef.iterations);

%% ---------------------------------------------------------------------
%% 2. The impact point, checked INDEPENDENTLY from the flown trajectory
%% ---------------------------------------------------------------------
%% Not from info, not from the summary: from traj.x(end,:), through
%% coorbital.util.greatCircle, in the documented (lat1,lon1,lat2,lon2) order.
%% This is the only assertion in the file that does not trust run_target's own
%% bookkeeping at all:
          missIndep = c.rE.*coorbital.util.greatCircle(trDef.x(end,3), ...
                                                       trDef.x(end,2), ...
                                                       latTg,lonTg);
    assert(missIndep <= tolM, ...
        ['the flown impact point is %.2f m from the target, outside the ' ...
         '%.1f m tolerance, whatever the summary reported'],missIndep,tolM);
    assertAbs(missIndep,inDef.missM,1e-6, ...
        'independently measured miss against the reported one (m)');

%% ...and the geometry is DISCRIMINATING, so the assertion above is shown to
%% be able to fail rather than assumed to be. Every way of transposing the
%% arguments must move the RANGE by far more than the miss assertion's budget;
%% if a future change to the shipped points ever brings one of them close, this
%% guard fails loudly instead of the check above going quietly blind.
%%
%% The margin is stated in metres against the 1 km tolerance, not in degrees
%% against a round number, because metres are what the assertions above are
%% written in. The worst of the three transpositions moves the range by 253 km,
%% a factor of 253 over the budget; a bare five-degree threshold would have
%% FAILED on this geometry while the range check it guards had a 253-to-1
%% margin, which would have been a threshold rejecting a sound test:
            angTrue = coorbital.util.greatCircle(latLn,lonLn,latTg,lonTg);
            angSwap = [coorbital.util.greatCircle(latLn,lonLn,lonTg,latTg); ...
                       coorbital.util.greatCircle(lonLn,latLn,latTg,lonTg); ...
                       coorbital.util.greatCircle(lonLn,latLn,lonTg,latTg)];
            dRngSwp = min(abs(c.rE.*(angSwap - angTrue)));
    assert(dRngSwp > 50.*tolM, ...
        ['the shipped geometry is not discriminating: a transposed ' ...
         'greatCircle call comes within %.1f m of the correct %.1f m, ' ...
         'against a %.1f m miss budget'],dRngSwp,c.rE.*angTrue,tolM);

%% ---------------------------------------------------------------------
%% 3. The azimuth, to machine precision, and the zero-cross-range property
%% ---------------------------------------------------------------------
%% The launch heading state must BE the closed-form bearing. Nothing
%% integrates the initial condition, so this is an equality and not a
%% tolerance -- 1e-13 rad is 6e-7 m of arc at Earth radius, far below any
%% physical scale in this problem:
            psiRef = coorbital.util.greatCircleBearing(latLn,lonLn,latTg,lonTg);
    assertAbs(trDef.x(1,6),psiRef,1e-13,'flown initial heading against the bearing (rad)');
    assertAbs(inDef.psiLaunch,psiRef,1e-13,'reported launch azimuth (rad)');

%% ...and the bearing must be discriminating too. A near-meridional pair is
%% what detects a latitude/longitude transposition inside greatCircleBearing;
%% here the check is that the shipped pair distinguishes the reverse course
%% and the transposed one from the true one:
            psiRev = coorbital.util.greatCircleBearing(latTg,lonTg,latLn,lonLn);
            psiTrn = coorbital.util.greatCircleBearing(lonLn,latLn,lonTg,latTg);
    assert(abs(rad2deg(wrapPi(psiRev - psiRef))) > 5, ...
        'the reverse bearing is only %.3f deg from the forward one', ...
        abs(rad2deg(wrapPi(psiRev - psiRef))));
    assert(abs(rad2deg(wrapPi(psiTrn - psiRef))) > 5, ...
        'a lat/lon transposed bearing is only %.3f deg from the true one', ...
        abs(rad2deg(wrapPi(psiTrn - psiRef))));

%% CROSS-RANGE IS ZERO IN THIS CONFIGURATION, and the whole miss is therefore
%% the range residual. That is the claim the summary makes in its limitations
%% block, and this is the measurement behind it. Not a general guarantee: it
%% holds because every commanded bank angle is zero:
    assertAbs(inDef.xTrackM,0,1e-3,'cross-track offset (m)');
    assertAbs(inDef.missM,abs(inDef.residM),1e-6, ...
        'the miss against the magnitude of the range residual (m)');

%% ---------------------------------------------------------------------
%% 4. The summary must REPORT what was computed
%% ---------------------------------------------------------------------
%% A number found correctly and then printed from the wrong variable is
%% invisible to every assertion above, and the printed summary is the only
%% part of this script most users ever read:
    assertAbs(summaryNumber(outDef,'required range +([-\d.]+) +km'), ...
        inDef.rngReqM./1000,1e-3,'reported required range (km)');
    assertAbs(summaryNumber(outDef,'achieved range +([-\d.]+) +km'), ...
        inDef.rngAchM./1000,1e-3,'reported achieved range (km)');
    assertAbs(summaryNumber(outDef,'MISS DISTANCE +([-\d.]+) +m'), ...
        inDef.missM,0.01,'reported miss distance (m)');
    assertAbs(summaryNumber(outDef,'solved cutoff +([-\d.]+) +s'), ...
        inDef.tCut,1e-3,'reported cutoff time (s)');
    assertAbs(summaryNumber(outDef,'reachable +([-\d.]+) to'), ...
        inDef.rngMinM./1000,1e-3,'reported envelope floor (km)');
    assertAbs(summaryNumber(outDef,'reachable +[-\d.]+ to ([-\d.]+) km'), ...
        inDef.rngMaxM./1000,1e-3,'reported envelope ceiling (km)');
    assertAbs(summaryNumber(outDef,'iterations +(\d+) '), ...
        inDef.iterations,0.5,'reported iteration count');
    assertAbs(summaryNumber(outDef,'launch azimuth +([-\d.]+) +deg'), ...
        rad2deg(psiRef),1e-5,'reported launch azimuth (deg)');

%% The cutoff must be reported as a FRACTION of the full burn as well as in
%% seconds -- the brief for this script asks for both, and seconds alone say
%% nothing without the burn they are a fraction of:
           fracTok = regexp(outDef,'\(([\d.]+) of the ([\d.]+) s full burn\)', ...
                            'tokens','once');
    assert(numel(fracTok) == 2, ...
        'the summary does not report the cutoff as a fraction of full burn:\n%s', ...
        outDef);
    assertAbs(str2double(fracTok{1}),inDef.cutFrac,1e-5,'reported cutoff fraction');
    assertAbs(str2double(fracTok{2}),inDef.tBurn  ,1e-3,'reported full burn time (s)');

%% THE TWO LIMITATIONS MUST BE STATED. Both are things a reader will otherwise
%% not notice, and both were required in writing; a future edit that quietly
%% drops either one is exactly what this checks:
    assert(contains(outDef,'NON-ROTATING EARTH'), ...
        'the summary must state that the azimuth is exact only for a non-rotating Earth');
    assert(contains(outDef,'OUTER ITERATION'), ...
        'the summary must say what rotation would require: an outer azimuth iteration');
    assert(contains(outDef,'NOT A GENERAL') && contains(outDef,'GUARANTEE'), ...
        'the summary must say that zero cross-range is a property of this configuration');

%% ---------------------------------------------------------------------
%% 5. A SECOND target, to prove the solve tracks it
%% ---------------------------------------------------------------------
%% Everything in part 1 is consistent with a script that ignores the target
%% and always returns the same cutoff. A second, nearer target on a different
%% azimuth is what rules that out. Narrower bracket and a 5 km tolerance,
%% because only the tracking matters here and each halving costs a
%% propagation:
            latTg2 = 30;
            lonTg2 = -140;
            tol2M  = 5000;
            out2   = evalc(sprintf( ...
                     ['[tr2,in2] = run_target(struct(''latTarget'',%.10g,' ...
                      '''lonTarget'',%.10g,''cutFracMin'',0.70,' ...
                      '''cutFracMax'',0.90,''tolRangeKm'',%.10g,' ...
                      '''showPlots'',false));'],latTg2,lonTg2,tol2M./1000));
    assert(contains(out2,'(nominal)') && ~contains(out2,'REFUSED'), ...
        'the second target did not solve. Summary was:\n%s',out2);
    assert(in2.missM <= tol2M, ...
        'the second solve missed by %.2f m against its %.0f m tolerance', ...
        in2.missM,tol2M);

%% Independently, from the flown state again:
            miss2  = c.rE.*coorbital.util.greatCircle(tr2.x(end,3),tr2.x(end,2), ...
                                                      deg2rad(latTg2),deg2rad(lonTg2));
    assert(miss2 <= tol2M, ...
        'the second flown impact point is %.2f m from its target',miss2);

%% ...and the solve genuinely MOVED. A nearer target must take a shorter burn
%% and a different azimuth, and the cutoff must not be the bracket midpoint,
%% which is what a solver that returned without iterating would hand back:
    assert(in2.rngReqM < inDef.rngReqM - 1e6, ...
        'the second target is not materially nearer; it cannot discriminate');
    assert(in2.tCut < inDef.tCut - 1, ...
        ['a %.0f km nearer target solved to a cutoff of %.4f s against the ' ...
         'shipped %.4f s; the solve is not responding to the target'], ...
        (inDef.rngReqM - in2.rngReqM)./1000,in2.tCut,inDef.tCut);
    assert(abs(rad2deg(wrapPi(in2.psiLaunch - inDef.psiLaunch))) > 1, ...
        'the two targets returned the same azimuth to within a degree');
            tMid2  = 0.5.*(in2.tCutLo + in2.tCutHi);
    assert(abs(in2.tCut - tMid2) > 1e-6, ...
        ['the solved cutoff %.9f s IS the bracket midpoint; the solver ' ...
         'returned without iterating'],in2.tCut);

%% ---------------------------------------------------------------------
%% 6. The reachable-envelope refusal, both ends
%% ---------------------------------------------------------------------
%% An unreachable target must be REFUSED, in words, with the band -- and must
%% NOT throw. Throwing would force every caller into a try/catch to read
%% numbers the returned struct already carries, and worse, a targeting script
%% that silently returned the nearest miss instead would hand back something
%% that LOOKS like a solution:
[outFar,infFar,threwFar,trFar] = tryRun(struct('latTarget',-30,'lonTarget',20, ...
                                               'showPlots',false));
    assert(~threwFar,'a target beyond the envelope must not throw:\n%s',outFar);
    assert(contains(outFar,'REFUSED'),'no refusal banner for a too-far target:\n%s',outFar);
    assert(contains(outFar,'TOO FAR'),'the refusal must say the target is too far:\n%s',outFar);
    assert(infFar.refused,'info.refused must be true on a refusal');
    assert(~infFar.solveInfo.converged,'the solve must report converged = false');
    assert(strcmp(infFar.solveInfo.identifier, ...
                  'coorbital:rangeSolve:targetOutsideBracket'), ...
        'expected the outside-bracket identifier, got "%s"', ...
        infFar.solveInfo.identifier);
    assert(infFar.rngReqM > infFar.rngMaxM, ...
        'the too-far case must have a required range above the band');
    assertRel(infFar.rngMaxM,7737630.3673877586,1e-4,'refused-case envelope ceiling (m)');

%% The band must be PRINTED, in kilometres, so a user can act on it without
%% reading the struct:
    assertAbs(summaryNumber(outFar,'reachable +([-\d.]+) to'), ...
        infFar.rngMinM./1000,0.01,'refusal message envelope floor (km)');
    assertAbs(summaryNumber(outFar,'reachable +[-\d.]+ to ([-\d.]+) km'), ...
        infFar.rngMaxM./1000,0.01,'refusal message envelope ceiling (km)');
    assertAbs(summaryNumber(outFar,'required range +([-\d.]+) +km'), ...
        infFar.rngReqM./1000,0.01,'refusal message required range (km)');

%% ...and the same at the near end, which is a different branch of the message
%% and a different thing to tell the user to change:
     [outNear,infNear,threwNear] = tryRun(struct('latTarget',20.4, ...
                                                 'lonTarget',-155, ...
                                                 'showPlots',false));
    assert(~threwNear,'a target inside the envelope floor must not throw:\n%s',outNear);
    assert(contains(outNear,'REFUSED'),'no refusal banner for a too-close target:\n%s',outNear);
    assert(contains(outNear,'TOO CLOSE'), ...
        'the refusal must say the target is too close:\n%s',outNear);
    assert(infNear.refused,'info.refused must be true on a refusal');
    assert(infNear.rngReqM < infNear.rngMinM, ...
        'the too-close case must have a required range below the band');
    assertAbs(summaryNumber(outNear,'reachable +([-\d.]+) to'), ...
        infNear.rngMinM./1000,0.01,'refusal message envelope floor (km)');

%% A refused run returns an EMPTY trajectory. Handing back a partial one would
%% invite a caller to plot the nearest miss as though it were the answer:
    assert(isempty(trFar),'a refused run must return an empty trajectory');

%% ---------------------------------------------------------------------
%% 7. The override hook itself
%% ---------------------------------------------------------------------
%% A misspelt parameter must raise rather than silently leave the shipped
%% value in place, or a future test could believe it flew a case it did not:
             threw = false;
    try
        run_target(struct('latTargetDeg',35));
    catch
             threw = true;
    end
    assert(threw,'an unrecognised override field must raise, not be ignored');

%% Nothing above asked for a figure, so nothing above may have made one:
            nFig1  = numel(findall(groot,'Type','figure'));
    assert(nFig1 == nFig0, ...
        ['run_target left %d new figure(s) open with showPlots false; the ' ...
         'suite must run headless'],nFig1 - nFig0);
end

function [out,info,threw,traj] = tryRun(optStruct)
%% Purpose:
%
%  Run run_target with the given overrides, capturing its printed output and
%  recording whether it threw. Exists because the refusal path must be shown
%  NOT to throw, and a bare call inside the test would abort the test file on
%  the very behaviour being denied.
%
%% Inputs:
%
%  optStruct        Struct                      USER PARAMETERS overrides, in
%                                               the block's own human units
%
%% Outputs:
%
%  out              Char [1 x n]                Captured summary text
%
%  info             Struct                      run_target's info output, or
%                                               an empty struct if it threw
%
%  threw            [1 x 1] logical             True if run_target raised
%
%  traj             Struct or []                run_target's trajectory output
%
%% Revision History:
%  Michael Casey                                                08/07/2026
%  Copyright 2026 Coorbital, Inc.
%% ------------------------ Begin Code Sequence ---------------------------

    assert(isstruct(optStruct),'tryRun takes a struct of USER PARAMETERS overrides.');
              info = struct();
              traj = [];
             threw = false;

%% evalc runs in THIS workspace, so optStruct is what run_target receives; the
%% assertion above is also what tells the code analyser that, since it cannot
%% see a variable used only inside an evaluated string:
    try
               out = evalc('[traj,info] = run_target(optStruct);');
    catch err
             threw = true;
               out = err.message;
    end
end

function a = wrapPi(a)
%% Purpose:
%
%  Wrap an angle difference into (-pi,pi], so that two bearings a hair either
%  side of north compare as nearly equal rather than as nearly a full turn
%  apart.
%
%% Inputs:
%
%  a                [1 x 1]                     Angle (rad)
%
%% Outputs:
%
%  a                [1 x 1]                     The same angle in (-pi,pi]
%                                               (rad)
%
%% Revision History:
%  Michael Casey                                                08/07/2026
%  Copyright 2026 Coorbital, Inc.
%% ------------------------ Begin Code Sequence ---------------------------

                 a = mod(a + pi,2.*pi) - pi;
end

function val = summaryNumber(txt,pat)
%% Purpose:
%
%  Pull one number out of a captured run_target summary. Parsing the printed
%  text rather than recomputing the quantity is deliberate: what a reader of
%  this library acts on is the summary, so the summary is what gets asserted.
%
%  The patterns use explicit space runs rather than \s, because in MATLAB \s
%  also matches a newline and would happily reach across into the next line of
%  the report.
%
%% Inputs:
%
%  txt              Char [1 x n]                Captured run_target summary
%
%  pat              Char [1 x m]                Regular expression with
%                                               exactly one capturing group
%                                               around the number
%
%% Outputs:
%
%  val              [1 x 1]                     The captured number
%
%% Revision History:
%  Michael Casey                                                08/07/2026
%  Copyright 2026 Coorbital, Inc.
%% ------------------------ Begin Code Sequence ---------------------------

               tok = regexp(txt,pat,'tokens','once');
    assert(~isempty(tok), ...
        'pattern "%s" matched nothing in the summary:\n%s',pat,txt);
               val = str2double(tok{1});
    assert(~isnan(val),'pattern "%s" captured "%s", which is not a number', ...
        pat,tok{1});
end

function assertRel(got,want,tol,what)
%% Purpose:
%
%  Assert a value against a reference to a relative tolerance, with a message
%  that reports both numbers and the actual relative error.
%
%% Inputs:
%
%  got              [1 x 1]                     Measured value
%
%  want             [1 x 1]                     Reference value (nonzero)
%
%  tol              [1 x 1]                     Relative tolerance (-)
%
%  what             Char [1 x n]                Name of the quantity, used in
%                                               the failure message
%
%% Outputs:
%
%  none                                         Throws on failure
%
%% Revision History:
%  Michael Casey                                                08/07/2026
%  Copyright 2026 Coorbital, Inc.
%% ------------------------ Begin Code Sequence ---------------------------

               rel = abs(got - want)./abs(want);
    assert(rel < tol, ...
        '%s: got %.12g against %.12g, %.3e relative, budget %.1e', ...
        what,got,want,rel,tol);
end

function assertAbs(got,want,tol,what)
%% Purpose:
%
%  Assert a value against a reference to an absolute tolerance. Used where the
%  reference is zero, or where the printed resolution rather than the value
%  sets the budget, so a relative test would be meaningless.
%
%% Inputs:
%
%  got              [1 x 1]                     Measured value
%
%  want             [1 x 1]                     Reference value
%
%  tol              [1 x 1]                     Absolute tolerance, in the
%                                               units of the quantity
%
%  what             Char [1 x n]                Name of the quantity, used in
%                                               the failure message
%
%% Outputs:
%
%  none                                         Throws on failure
%
%% Revision History:
%  Michael Casey                                                08/07/2026
%  Copyright 2026 Coorbital, Inc.
%% ------------------------ Begin Code Sequence ---------------------------

    assert(abs(got - want) < tol, ...
        '%s: got %.12g against %.12g, %.3e absolute, budget %.1e', ...
        what,got,want,abs(got - want),tol);
end
