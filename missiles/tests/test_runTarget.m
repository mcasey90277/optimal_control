function test_runTarget()
%% Purpose:
%
%  Pin HGV/run_target, the point-to-point targeting script: the TWO-AXIS solve
%  that drives the down-range and cross-range components of the miss to zero
%  in the launch azimuth and the thrust-termination time, the bisection that
%  seeds it, the separation link that throws the whole booster away, the
%  reachable-envelope refusal, and the summary that reports all of it. Every
%  other entry-script test in this suite checks a trajectory that was flown
%  FORWARD from a given azimuth; this one checks a trajectory that was SOLVED
%  to arrive somewhere, which is a different claim and fails in different ways.
%
%  THE FOUR THINGS THIS TEST EXISTS TO CATCH, and how each is caught:
%
%    A wrong azimuth.    The impact point walks off the launch-to-target great
%                        circle while the range residual stays small, so the
%                        solve looks converged and the vehicle still misses.
%                        Caught by measuring the miss impact-to-target with
%                        coorbital.util.greatCircle, from the flown state,
%                        rather than by reading any solver's own residual.
%
%    An azimuth that is  Caught by parts 6 and 8b, the banked and the rotating
%    not SOLVED at all.  cases, which are the two configurations where the
%                        closed-form great-circle bearing is NOT the answer.
%                        Force it back to the closed form and both of them
%                        miss -- by 21.5 km and by 231.6 km -- while the
%                        shipped zero-bank non-rotating case, where the seed
%                        IS the answer, carries on passing. That asymmetry is
%                        the whole point of those two sections.
%
%    A seed solve that   Caught by the same measured miss: a cutoff returned
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
%  takes fourteen of them -- twelve in the seed bisection, one in the Newton,
%  one to re-fly the answer -- so the whole test runs in a few seconds and
%  there was no need to loosen the tolerance or shrink the bracket for it. The
%  rotating case in part 8b is the dearest at twenty-three, because its Newton
%  has real work to do. The SHIPPED
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
%  Michael Casey  The display choices -- the adaptive altitude
%                 exaggeration, its cap and floor, and the
%                 hemisphere captions -- none of which was
%                 pinned by anything                            08/07/2026
%  Michael Casey  Two-axis targeting: the banked case now ARRIVES
%                 rather than warning, the rotating case is flown
%                 rather than refused, and both are pinned on the
%                 miss they used to have and the miss they have
%                 now                                           08/08/2026
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
    assert(~contains(outDef,'targeting: REFUSED'), ...
        'the shipped target must be reachable. Summary was:\n%s',outDef);
    assert(~inDef.refused,'info.refused must be false on a converged solve');
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
    assert(inDef.seedIter == 10, ...
        'the seed bisection took %d steps against the expected 10',inDef.seedIter);
    assert(inDef.seedEval == inDef.seedIter + 2, ...
        ['%d evaluations for %d iterations; rangeSolve must cost exactly two ' ...
         'endpoint evaluations plus one per step'],inDef.seedEval,inDef.seedIter);

%% THE SHIPPED CASE IS THE ONE WHERE THE SEED IS ALREADY THE ANSWER, and the
%% cheapest possible two-axis solve is the evidence: aimSolve tests its
%% convergence BEFORE it does any work, so a guess already inside tolerance
%% costs exactly one propagation and zero iterations, and the azimuth comes
%% back bit-for-bit unmoved. A non-zero iteration count here would mean the
%% closed form had stopped being exact on a non-rotating zero-bank run, which
%% is a physics change and not a solver change:
    assert(inDef.aimIter == 0, ...
        ['the two-axis solve took %d iteration(s) on the shipped case; the ' ...
         'seed is exact there and it must converge at its first evaluation'], ...
        inDef.aimIter);
    assert(inDef.aimEval == 1, ...
        'the two-axis solve cost %d propagation(s) where 1 is the whole bill', ...
        inDef.aimEval);
    assert(inDef.nProp == inDef.seedEval + inDef.aimEval + 1, ...
        ['%d propagations reported against %d + %d + 1; the total must be the ' ...
         'two solves plus the one re-flight'], ...
        inDef.nProp,inDef.seedEval,inDef.aimEval);
    assert(inDef.tCut == inDef.tSeed, ...
        ['the solved cutoff %.15f s differs from the seed %.15f s on a case ' ...
         'the seed already solved'],inDef.tCut,inDef.tSeed);

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
%% The SEED is the closed-form bearing, and on this configuration the SOLVED
%% azimuth is the same number to the bit -- not to a tolerance. That is the
%% claim step 5 of the brief asks for: a two-axis solve over a case whose
%% closed form is exact must return the closed form. Nothing integrates the
%% initial condition either, so the flown heading is the same equality again:
            psiRef = coorbital.util.greatCircleBearing(latLn,lonLn,latTg,lonTg);
    assertAbs(inDef.psiSeed,psiRef,1e-13,'reported seed azimuth (rad)');
    assert(inDef.psiLaunch == inDef.psiSeed, ...
        ['the solved azimuth %.17g rad is not the seed %.17g rad; on a ' ...
         'non-rotating zero-bank run the closed form IS the answer and the ' ...
         'solve must hand it back unchanged'],inDef.psiLaunch,inDef.psiSeed);
    assert(inDef.dPsiAim == 0, ...
        'the solve reported a %.3e rad azimuth correction where it made none', ...
        inDef.dPsiAim);
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
%% down-range. That is the claim the summary makes in its limitations block,
%% and this is the measurement behind it. Not a general guarantee: it holds
%% because every commanded bank angle is zero AND the Earth is not turning:
    assertAbs(inDef.xTrackM,0,1e-3,'cross-range component of the miss (m)');
    assertAbs(inDef.downM,inDef.missM,1e-6, ...
        'the down-range component against the whole miss (m)');
    assertAbs(inDef.missM,abs(inDef.residM),1e-6, ...
        'the miss against the magnitude of the range residual (m)');

%% ...and the decomposition must COMPOSE. hypot of the two components is the
%% miss by construction, so this catches an aim frame that is not a unit basis
%% -- a psiCourse taken at the wrong end of the arc, say, or a scale factor
%% applied to one component and not the other -- which every zero-cross-range
%% assertion above is blind to:
    assertAbs(hypot(inDef.downM,inDef.xTrackM),inDef.missM,1e-6, ...
        'the two components composed against the measured miss (m)');
    assertAbs(inDef.missHypM,inDef.missM,1e-6, ...
        'the reported composition against the measured miss (m)');

%% THE SEED MISS IS THE MISS, here and only here, because the solve had
%% nothing to do. It is pinned because parts 6 and 8b pin it too, and it is
%% the number that makes those two sections mean something: the same field
%% reads 21.5 km and 231.6 km there:
    assertAbs(inDef.seedMissM,inDef.missM,1e-6, ...
        'the seed miss against the solved miss on a case the seed solved (m)');

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
        inDef.seedIter + inDef.aimIter,0.5,'reported iteration count');
    assertAbs(summaryNumber(outDef,'propagations +(\d+) '), ...
        inDef.nProp,0.5,'reported propagation count');
    assertAbs(summaryNumber(outDef,'seed azimuth +([-\d.]+) +deg'), ...
        rad2deg(psiRef),1e-5,'reported seed azimuth (deg)');
    assertAbs(summaryNumber(outDef,'launch azimuth +([-\d.]+) +deg'), ...
        rad2deg(psiRef),1e-5,'reported launch azimuth (deg)');
    assertAbs(summaryNumber(outDef,'seed miss +([-\d.]+) +m'), ...
        inDef.seedMissM,0.01,'reported seed miss (m)');
    assertAbs(summaryNumber(outDef,'down-range miss +([-+\d.]+) +m'), ...
        inDef.downM,0.01,'reported down-range component (m)');
    assertAbs(summaryNumber(outDef,'cross-range miss +([-+\d.]+) +m'), ...
        inDef.xTrackM,0.01,'reported cross-range component (m)');

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

%% THE THREE STANDING STATEMENTS MUST BE MADE. Each is something a reader will
%% otherwise not notice, and a future edit that quietly drops one is exactly
%% what this checks. They are also the three the file used to get WRONG in the
%% opposite direction: it claimed a closed-form azimuth, claimed the solve saw
%% only downrange, and claimed neither could be helped:
    assert(contains(outDef,'THE AZIMUTH IS SOLVED, NOT CLOSED FORM'), ...
        'the summary must say the azimuth is solved and the closed form is only a seed');
    assert(contains(outDef,'THE SOLVE CONTROLS BOTH AXES'), ...
        'the summary must say that cross-range is driven to zero, not merely measured');
    assert(contains(outDef,'THE NEWTON IS LOCAL'), ...
        'the summary must say the solve is local and that a failure is refused');
    assert(~contains(outDef,'OUTER ITERATION'), ...
        ['the summary still asks for the outer azimuth iteration that this ' ...
         'script now has']);
%% Matched on CONTIGUOUS fragments rather than on whole sentences: the claims
%% are wrapped across printed lines, and a pattern spanning a wrap would break
%% on any rewording of the line breaks while saying nothing about the content:
    assert(contains(outDef,'PROPERTY OF THIS CONFIGURATION') && ...
           contains(outDef,'GENERAL GUARANTEE'), ...
        'the summary must say that zero cross-range is a property of this configuration');
    assert(contains(outDef,'IT IS NOT MODELLED ON THIS RUN'), ...
        'the shipped run must say plainly that it did not turn the Earth');

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
    assert(contains(out2,'(nominal)') && ~contains(out2,'targeting: REFUSED'), ...
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
%% 6. A BANKED case, which is one of the two cases where the closed-form
%%    azimuth is NOT the answer and the two-axis solve has to earn its keep
%% ---------------------------------------------------------------------
%% WHY THIS CASE EXISTS, and it is one of the two most important sections in
%% the file. Everything above flies at zero bank over a still Earth, where the
%% seed azimuth is exact, the Newton converges at its first evaluation without
%% moving anything, and the measured impact-to-target miss and the magnitude of
%% the range residual agree to 9.3e-10 m. Nothing up there can tell a SOLVED
%% azimuth from a closed-form one, and nothing up there can tell a measured
%% miss from a restated residual: two mutations that did exactly that -- forcing
%% crossWarn = false, and setting missM = abs(resM) -- survived an earlier
%% version of this suite.
%%
%% A 75 deg terminal bank breaks all of that. It walks the impact point 21.5 km
%% off the launch-to-target great circle while the range residual stays at
%% 639 m, which is what this configuration USED TO SHIP as a converged solve
%% beside a WARNING. It now arrives. What is pinned is therefore both numbers:
%% the miss the seed alone still produces, which is the old failure preserved
%% as a measurement, and the miss after the solve, which is the fix. Force the
%% azimuth back to the closed form and the second of those becomes the first:
            outBnk = evalc(['[trBnk,inBnk] = run_target(struct(' ...
                            '''descBank'',[75 75],''showPlots'',false));']);
    assert(contains(outBnk,'(nominal)') && ~contains(outBnk,'targeting: REFUSED'), ...
        'the banked case did not solve. Summary was:\n%s',outBnk);

%% THE SEED STILL MISSES BY 21.5 KM, and that is the number this script used to
%% report as its answer. Pinned so the section keeps its teeth: if a future
%% change made the bank harmless, the case would stop discriminating and this
%% would say so rather than passing quietly:
%% Both read at full precision out of a run of the COMMITTED code, like every
%% other literal in this file. The cross-range in particular is NOT the
%% 21515.222106821158 the previous version pinned: that figure came from the
%% old asin(sin(angFly)*sin(dPsi)) cross-track formula, and the tangent-plane
%% decomposition that replaced it gives 21515.222142760958. The two agree to
%% 3.6e-5 m, which is a useful independent check of the new construction
%% against the one it replaced -- but agreement is not provenance, and a
%% literal carried across from the previous formula would be a number this
%% suite never measured:
    assertRel(inBnk.seedMissM  ,21524.695285497211,1e-3,'banked-case SEED miss (m)');
    assertRel(inBnk.seedXTrackM,21515.222142760958,1e-3,'banked-case SEED cross-range (m)');
    assert(inBnk.seedMissM > 20.*tolM, ...
        ['the seed misses by only %.2f m at 75 deg of bank; this case no ' ...
         'longer discriminates a solved azimuth from a closed-form one'], ...
        inBnk.seedMissM);

%% ...AND THE SOLVE ARRIVES. This is the assertion the whole task turns on:
    assert(inBnk.missM <= tolM, ...
        ['the banked run missed by %.2f m against a %.1f m tolerance; the ' ...
         'two-axis solve did not close the cross-range'],inBnk.missM,tolM);
    assertRel(inBnk.missM  ,54.795200772654532,1e-3,'banked-case miss distance (m)');
    assertRel(inBnk.downM  ,-54.72406984419969,1e-3,'banked-case down-range miss (m)');
    assertRel(inBnk.xTrackM,-2.791094302003172,1e-2,'banked-case cross-range miss (m)');

%% THE AZIMUTH MOVED, and by the pinned amount. A solve that converged by
%% moving only the cutoff would leave the cross-range where it was, so the
%% correction is the mechanism and not a by-product:
    assertRel(rad2deg(inBnk.dPsiAim),-0.3436388265495105,1e-3, ...
        'banked-case azimuth correction (deg)');
    assert(abs(inBnk.dPsiAim) > 1e-4, ...
        ['the banked case solved to an azimuth %.3e rad from the closed form; ' ...
         'it cannot have closed 21.5 km of cross-range without aiming off'], ...
        inBnk.dPsiAim);
    assertRel(inBnk.tCut,75.860087783793503,1e-6,'banked-case solved cutoff (s)');
    assert(inBnk.aimIter >= 1, ...
        'the banked case took %d Newton iteration(s); the seed is not the answer there', ...
        inBnk.aimIter);

%% THE ASSERTION THAT KILLS THE RESIDUAL SUBSTITUTION, kept from the version
%% that shipped the warning. It has moved to the SEED, which is where the
%% separation between a measured miss and a restated range residual now lives:
%% 21525 m against 639 m, a factor of 33.7:
    assert(inBnk.seedMissM > 30.*abs(inBnk.seedDownM), ...
        ['the banked seed miss is %.2f m against a %.2f m down-range ' ...
         'component, a factor of only %.2f. The miss looks like it is being ' ...
         'read back out of a downrange-only solve rather than measured from ' ...
         'the flown impact point; at 75 deg of bank the two must differ by ' ...
         'more than thirty'],inBnk.seedMissM,abs(inBnk.seedDownM), ...
        inBnk.seedMissM./abs(inBnk.seedDownM));

%% ...and the arrival checked once more from the flown state alone, which at
%% nonzero bank is a genuinely independent number rather than the same one by
%% construction:
            missBnk = c.rE.*coorbital.util.greatCircle(trBnk.x(end,3), ...
                                                       trBnk.x(end,2), ...
                                                       latTg,lonTg);
    assertRel(missBnk,inBnk.missM,1e-9, ...
        'independently measured banked-case miss against the reported one (m)');
    assert(missBnk <= tolM, ...
        ['the INDEPENDENTLY measured banked impact point is %.2f m from the ' ...
         'target, outside the %.1f m tolerance, whatever the summary said'], ...
        missBnk,tolM);

%% NO WARNING, AND NO REFUSAL. The cross-range warning was the old script's way
%% of admitting it could not fly this case; a run that arrives must not print
%% it, and the absence is pinned so that reinstating the warning without
%% reinstating the failure would be caught:
    assert(~contains(outBnk,'WARNING'), ...
        ['the banked run printed a WARNING although it arrived within ' ...
         'tolerance. Summary was:\n%s'],outBnk);
    assert(~contains(outBnk,'solve CONVERGED and the vehicle still missed'), ...
        ['the banked run still prints the converged-and-missed warning, which ' ...
         'is no longer true of it. Summary was:\n%s'],outBnk);

%% The limitations paragraph must FORK on the bank, and each branch must be the
%% one that fits its run. The zero-bank wording -- "a zero-bank track over a
%% non-rotating sphere never leaves the great circle it departed on" -- printed
%% verbatim under a 75 deg descent bank would be self-refuting:
    assert(contains(outBnk,'NON-ZERO BANK'), ...
        ['the banked run must say plainly that it commands a non-zero bank. ' ...
         'Summary was:\n%s'],outBnk);
    assert(~contains(outBnk,'and rotation is off, so the seed left only'), ...
        'the banked run printed the zero-bank explanation. Summary was:\n%s', ...
        outBnk);
    assert(contains(outDef,'and rotation is off, so the seed left only'), ...
        'the zero-bank explanation vanished from the shipped run');
    assert(contains(outDef,'seed azimuth WAS the answer here'), ...
        'the shipped run must say that its seed was already the answer');
    assert(~contains(outDef,'NON-ZERO BANK'), ...
        'the shipped zero-bank run claimed a non-zero bank');

%% The handoff wording forks on the same fact, and for the same reason: at
%% unequal banks the handoff is a real control discontinuity, not bookkeeping:
    assert(contains(outDef,'banks are EQUAL'), ...
        'the shipped run must call the handoff pure bookkeeping');
    assert(contains(outBnk,'banks DIFFER'), ...
        'the banked run must not call the handoff pure bookkeeping');

%% ---------------------------------------------------------------------
%% 7. The reachable-envelope refusal, both ends
%% ---------------------------------------------------------------------
%% An unreachable target must be REFUSED, in words, with the band -- and must
%% NOT throw. Throwing would force every caller into a try/catch to read
%% numbers the returned struct already carries, and worse, a targeting script
%% that silently returned the nearest miss instead would hand back something
%% that LOOKS like a solution:
[outFar,infFar,threwFar,trFar] = tryRun(struct('latTarget',-30,'lonTarget',20, ...
                                               'showPlots',false));
    assert(~threwFar,'a target beyond the envelope must not throw:\n%s',outFar);
    assert(contains(outFar,'targeting: REFUSED'),'no refusal banner for a too-far target:\n%s',outFar);
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
    assert(contains(outNear,'targeting: REFUSED'),'no refusal banner for a too-close target:\n%s',outNear);
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
%% 8. The DISPLAY choices: the adaptive exaggeration and the captions
%% ---------------------------------------------------------------------
%% Neither moves a number in the summary, and both were shipped with nothing
%% asserting anything about them. That is not harmless. The exaggeration
%% decides whether a reader sees a glide or an orbit, and the captions are
%% where most readers meet these coordinates at all -- a caption that reports
%% 155 deg west as "155.00E" is a figure that lies about which ocean the
%% vehicle crossed.
%%
%% TRUE SCALE IS THE SHIPPED DEFAULT, altExag = 1. The movie's altitude inset is
%% true-scale already and carries the profile the globe used to exaggerate for,
%% and coorbital.viz drops the "(altitude exaggerated Nx)" caption clause at
%% unity, so a true-scale picture makes no claim it is not keeping:
    assert(inDef.altExag == 1, ...
        'the shipped run must draw at TRUE SCALE; it drew at %gx',inDef.altExag);

%% ...and the ADAPTIVE rule is one word away, altExag = 'auto'. It is still the
%% right answer when the skip phugoid is the point of the picture, and it still
%% gives 16x here and not the 30x this script used to hard-code:
%% 0.3 * 6378.137 / 118.0858 = 16.2, floored. The sentinel is FLOWN rather than
%% asserted about, because a sentinel nothing exercises is a sentinel that rots:
             hPkFlw = max(trDef.x(:,1)) - c.rE;
              outAu = evalc(['[trAu,inAu] = run_target(struct(''altExag'',' ...
                             '''auto'',''tolRangeKm'',25,''showPlots'',false));']);
    assert(contains(outAu,'(nominal)') && ~contains(outAu,'targeting: REFUSED'), ...
        'the auto-exaggeration run did not solve. Summary was:\n%s',outAu);
    assert(inAu.altExag == 16, ...
        ['altExag = ''auto'' drew at %gx against the 16x the adaptive rule ' ...
         'gives for its %.4f km apogee'],inAu.altExag,inAu.hGlideMax);
    assert(inAu.altExag ~= inDef.altExag, ...
        'the auto run drew at the shipped true scale; the sentinel did nothing');
    assert(inAu.altExag == inAu.altExagRule(max(trAu.x(:,1)) - c.rE), ...
        ['the auto run drew at %gx while its own rule gives %gx for the ' ...
         '%.4f km apogee it flew'],inAu.altExag, ...
        inAu.altExagRule(max(trAu.x(:,1)) - c.rE),(max(trAu.x(:,1)) - c.rE)./1000);
    assert(inDef.altExagRule(hPkFlw) == 16, ...
        ['the rule gives %gx for the shipped run''s %.4f km apogee; true ' ...
         'scale is a CHOICE made over it, not the absence of a rule'], ...
        inDef.altExagRule(hPkFlw),hPkFlw./1000);

%% THE CAP AND THE FLOOR, asked of the rule directly. Neither is reachable
%% from any configuration this script can fly -- the cap binds only under
%% 63.8 km of apogee and the floor only above 956.7 km, against a 118 km
%% shipped peak -- so a check that only ever flies trajectories cannot see
%% either of them, and both were unpinned. info.altExagRule exists for exactly
%% this: it is the rule itself, so the two ends can be asked about:
    assert(inDef.altExagRule(40e3) == 30, ...
        ['the rule returns %g at a 40 km apogee; it must CAP at 30, the ' ...
         'factor the rest of this library uses by hand. Uncapped it would ' ...
         'return %g and paint a low flight halfway off the frame.'], ...
        inDef.altExagRule(40e3),floor(0.3.*c.rE./40e3));
    assert(inDef.altExagRule(0) == 30, ...
        ['a zero apogee returns %g; the rule must survive a flight that never ' ...
         'left the pad and hand it the cap.'],inDef.altExagRule(0));
    assert(inDef.altExagRule(2e6) == 2, ...
        ['the rule returns %g at a 2000 km apogee; it must FLOOR at 2. Below ' ...
         '2x an exaggeration is indistinguishable from true scale, and a ' ...
         'floor of 0 or 1 would hand coorbital.viz an AltScale it refuses.'], ...
        inDef.altExagRule(2e6));

%% ...and between the two it holds its stated invariant, that the APPARENT
%% apogee stays at or under 0.3 rE. Checked across the whole band rather than
%% at one point, because a rule that happened to be right at 118 km and wrong
%% everywhere else would pass every assertion above:
    for hKmT = [70 100 150 250 400 700 950]
              eGot = inDef.altExagRule(hKmT.*1000);
        assert(eGot >= 2 && eGot <= 30 && eGot == fix(eGot), ...
            'the rule returns %g at a %g km apogee; it must be a whole number in [2 30].', ...
            eGot,hKmT);
        assert(eGot.*hKmT.*1000 <= 0.3.*c.rE, ...
            ['at a %g km apogee the rule returns %gx, which paints the arc ' ...
             '%.1f km off the surface against the 0.3 rE = %.1f km it holds ' ...
             'itself to'],hKmT,eGot,eGot.*hKmT,0.3.*c.rE./1000);
    end

%% ...and it is MONOTONE: a higher flight never gets a bigger exaggeration.
%% A rule that was not would paint two similar flights at wildly different
%% scales and give a reader no way to compare them:
            hSwp   = (20:20:1200).*1000;
            eSwp   = arrayfun(inDef.altExagRule,hSwp);
    assert(all(diff(eSwp) <= 0), ...
        'the exaggeration rule is not monotone decreasing in apogee.');

%% THE CAPTIONS CARRY HEMISPHERE LETTERS. The shipped launch and target are
%% both WEST of the meridian, so a formatter that never reached its west branch
%% would produce a caption that is wrong by 310 degrees of longitude and looks
%% perfectly reasonable:
    assert(strcmp(inDef.launchStr,'20.00N 155.00W'), ...
        'the launch caption reads "%s" against the expected "20.00N 155.00W"', ...
        inDef.launchStr);
    assert(strcmp(inDef.targetStr,'35.00N 120.00W'), ...
        'the target caption reads "%s" against the expected "35.00N 120.00W"', ...
        inDef.targetStr);
    assert(strcmp(in2.targetStr,'30.00N 140.00W'), ...
        ['the second target''s caption reads "%s"; the captions must follow ' ...
         'the overrides, not the shipped block'],in2.targetStr);

%% THE OVERRIDE THE HEADER PROMISES MUST WORK. It did not: with altExag
%% supplied, a guard skipped the assignment and the very next line read the
%% variable that was never assigned, so every documented use of this override
%% died with "Unrecognized function or variable 'altExag'". A bracket and a
%% loose tolerance keep the extra solve cheap; what is checked is the display
%% choice, not the trajectory:
             exgWnt = 7;
              outEx = evalc(sprintf( ...
                      ['[trEx,inEx] = run_target(struct(''altExag'',%d,' ...
                       '''tolRangeKm'',25,''showPlots'',false));'],exgWnt));
    assert(contains(outEx,'(nominal)') && ~contains(outEx,'targeting: REFUSED'), ...
        'the altExag override run did not solve. Summary was:\n%s',outEx);
    assert(inEx.altExag == exgWnt, ...
        'the run was given altExag = %d and drew at %g',exgWnt,inEx.altExag);
    assert(exgWnt ~= inDef.altExag, ...
        ['the override was set to the value the block already ships; this ' ...
         'check cannot tell an honoured override from an ignored one']);
    assert(~isempty(trEx),'the altExag override run returned no trajectory');

%% ...and a nonsense one must RAISE, before the range solve rather than after
%% it, or a mistyped sentinel would cost a whole solve before it was noticed:
            threwEx = false;
    try
        run_target(struct('altExag','atuo','showPlots',false));
    catch
            threwEx = true;
    end
    assert(threwEx, ...
        'a mistyped altExag sentinel must raise, not be silently accepted');

%% ---------------------------------------------------------------------
%% 8b. A ROTATING EARTH IS FLOWN, not refused
%% ---------------------------------------------------------------------
%% earthSpin true is the other case where the closed-form bearing is not the
%% answer, and it is the harder of the two: the deflection is 231.552 km rather
%% than 21.5 km, and it is almost entirely CROSSWISE, so a downrange-only solve
%% converges on range and lands a quarter of a thousand kilometres away. This
%% script first printed that as a converged solution, then refused the case
%% outright; it now flies it.
%%
%% Note what rotation does and does not do here, because the refusal this
%% replaces described it the other way round: the state is EARTH-FIXED with a
%% planet-relative speed, so the target does not slide east under the vehicle.
%% It sits still and the VEHICLE is deflected. The two descriptions differ in
%% sign as well as in mechanism, and only one of them matches the equations of
%% motion in coorbital.eom.glide3DOF:
             outSpn = evalc(['[trSpn,inSpn] = run_target(struct(''earthSpin'',' ...
                             'true,''showPlots'',false));']);
    assert(contains(outSpn,'(nominal)') && ~contains(outSpn,'targeting: REFUSED'), ...
        'the rotating case did not solve. Summary was:\n%s',outSpn);
    assert(~inSpn.refused,'info.refused must be false on a flown rotating case');
    assert(~isempty(trSpn),'a flown run must return a trajectory');
    assert(contains(outSpn,'Earth rotation ON'), ...
        'the rotating run must report Earth rotation ON:\n%s',outSpn);
    assert(contains(outSpn,'IT IS MODELLED ON THIS RUN'), ...
        'the rotating run must say the rotation was modelled:\n%s',outSpn);
    assert(contains(outSpn,'THE ROTATION and not from any commanded turn'), ...
        ['the rotating zero-bank run must attribute its cross-range to the ' ...
         'rotation rather than to a bank:\n%s'],outSpn);

%% THE SEED STILL MISSES BY 231.552 KM. This is the measurement that used to be
%% the answer, kept as evidence rather than as a memory in a comment. It is
%% almost all cross-range, which is exactly why the second control was needed:
    assertRel(inSpn.seedMissM  ,231551.6280300247,1e-4,'rotating-case SEED miss (m)');
    assertRel(inSpn.seedXTrackM,231477.49939778377,1e-4,'rotating-case SEED cross-range (m)');
    assert(abs(inSpn.seedXTrackM) > 0.99.*inSpn.seedMissM, ...
        ['only %.2f m of the %.2f m seed miss is crosswise; a downrange-only ' ...
         'solve would have caught the rest, and this case would not be the ' ...
         'discriminator it is claimed to be'],abs(inSpn.seedXTrackM), ...
        inSpn.seedMissM);

%% ...AND THE SOLVE ARRIVES:
    assert(inSpn.missM <= tolM, ...
        ['the rotating run missed by %.2f m against a %.1f m tolerance; the ' ...
         'two-axis solve did not remove the Coriolis deflection'], ...
        inSpn.missM,tolM);
    assertRel(inSpn.missM,4.3614169744792672,1e-2,'rotating-case miss distance (m)');
    assertRel(rad2deg(inSpn.dPsiAim),-3.5792824249972601,1e-4, ...
        'rotating-case azimuth correction (deg)');
    assertRel(inSpn.tCut,75.207196173000099,1e-6,'rotating-case solved cutoff (s)');
    assert(inSpn.aimIter >= 1, ...
        'the rotating case took %d Newton iteration(s); it cannot have solved itself', ...
        inSpn.aimIter);

%% ...measured independently from the flown terminal state, on the impact
%% sphere, exactly as part 2 does for the shipped case:
            missSpn = c.rE.*coorbital.util.greatCircle(trSpn.x(end,3), ...
                                                       trSpn.x(end,2), ...
                                                       latTg,lonTg);
    assertRel(missSpn,inSpn.missM,1e-9, ...
        'independently measured rotating-case miss against the reported one (m)');
    assert(missSpn <= tolM, ...
        ['the INDEPENDENTLY measured rotating impact point is %.2f m from the ' ...
         'target'],missSpn);

%% THE ROTATION MUST ACTUALLY BE ON in the run that was flown. Without this,
%% an override that silently failed to reach env.omegaE would make every
%% assertion above pass by flying the shipped case twice:
    assert(inSpn.env.omegaE > 0, ...
        'the rotating run flew with env.omegaE = %.6e; the override did nothing', ...
        inSpn.env.omegaE);
    assert(abs(rad2deg(wrapPi(inSpn.psiLaunch - inDef.psiLaunch))) > 1, ...
        ['the rotating case solved to the same azimuth as the shipped one to ' ...
         'within a degree; it cannot have removed a 231 km deflection']);

%% ---------------------------------------------------------------------
%% 9. The override hook itself
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
