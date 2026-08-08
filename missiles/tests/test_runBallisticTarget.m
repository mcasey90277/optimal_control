function test_runBallisticTarget()
%% Purpose:
%
%  Pin BM/run_ballistic_target, the ballistic point-to-point targeting script:
%  the closed-form launch azimuth, the bracketing of the max-range loft angle,
%  the TWO branch solves either side of it, the branch selector, the
%  after-the-fact measurement of which branch was actually flown, and the
%  refusals at both ends of the envelope and for an unreachable branch.
%
%  WHAT MAKES THIS DIFFERENT FROM test_runTarget, and it is the whole reason
%  this file is separate. HGV/run_target ranges on thrust-termination time,
%  which is monotonic in range, so there is ONE answer and pinning it is
%  enough. A ballistic trajectory has TWO answers for every range short of
%  maximum -- a lofted arc and a depressed one -- and almost everything that
%  can go wrong here goes wrong by returning the WRONG ONE while still
%  converging on range. So:
%
%    THE BRANCH IS CHECKED INDEPENDENTLY OF THE LABEL. run_ballistic_target
%    reports which branch it thinks it flew; this file re-derives it from the
%    flown apogee and the flown flight time against the max-range arc's own,
%    and asserts the two agree. A selector wired to a constant passes every
%    range assertion and fails this one.
%
%    THE TWO BRANCHES ARE ASSERTED TO BE DIFFERENT. Ten times the apogee and
%    nearly three times the flight time separate them on the shipped geometry.
%    A script that solved the same branch twice and labelled them differently
%    would satisfy every miss assertion in this file and fail that one.
%
%    THE SELECTOR IS SHOWN NOT TO BE A CONSTANT. minimum-energy picks the
%    LOFTED arc at the shipped 3175 km target and the DEPRESSED arc at 4218 km,
%    because the ranking of the two loft-angle distances flips between them.
%    Both are flown here. One case alone cannot tell a working rule from a
%    hard-wired answer.
%
%  WHY THE FLOWN CASE IS NEITHER EQUATORIAL NOR DUE EAST. From (0,0) the
%  central angle reduces to acos(cos(lat2) cos(lon2)), which is symmetric in
%  its two arguments, so a due-east equatorial case is provably BLIND to a
%  latitude/longitude transposition at every great-circle call site. The
%  shipped configuration is 45 N 100 W to 62 N 60 W on a 39.2 deg azimuth,
%  where all four coordinates are distinct; part 2 below PROVES that geometry
%  is discriminating rather than assuming it, in metres against the miss
%  budget the assertions are written in.
%
%  THE alphaMax GUARD IS PINNED TOO, in part 9, and it is not a formality. At
%  BM/run_ballistic's 6 deg clamp the vehicle cannot be pushed past the
%  max-range attitude at any commanded loft angle, range is monotone across the
%  whole bracket, and the depressed branch does not exist. The script refuses
%  with coorbital:runBallisticTarget:maximumNotBracketed rather than silently
%  reporting a one-branch problem as a two-branch one; that refusal is the only
%  thing standing between a tightened clamp and a "lofted" answer with nothing
%  on the other side of it.
%
%  COST. A trajectory propagation is about 0.05 s and a full solve takes 57 of
%  them, so a run is about 3 s and the whole file about 35 s. The SHIPPED
%  configuration is therefore flown at its shipped settings, which is the
%  stronger test: what is pinned is exactly what a user gets.
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
%% They are NOT read out of run_ballistic_target: a test that took its
%% reference from the code under test could not see the code change them:
            latLnD = 45;
            lonLnD = -100;
            latTgD = 62;
            lonTgD = -60;
             latLn = deg2rad(latLnD);
             lonLn = deg2rad(lonLnD);
             latTg = deg2rad(latTgD);
             lonTg = deg2rad(lonTgD);

%% The shipped range tolerance, restated rather than read from the code. Every
%% miss assertion in this file is written against it:
              tolM = 1000;

%% ---------------------------------------------------------------------
%% 1. The shipped configuration: branch = 'minimum-energy'
%% ---------------------------------------------------------------------
%% Plots off so the suite creates no figures, and the summary captured rather
%% than printed so the suite output stays readable. The captured text is not
%% discarded -- part 4 asserts against it, because the summary is what a reader
%% of this library actually acts on:
            outDef = evalc(['[trDef,inDef] = ' ...
                            'run_ballistic_target(struct(''showPlots'',false));']);

    assert(contains(outDef,'(nominal)'), ...
        'the shipped run did not terminate nominally. Summary was:\n%s',outDef);
    assert(~contains(outDef,'CAUTION'), ...
        'the shipped run printed a caution. Summary was:\n%s',outDef);
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
%% fire first and report a 1e-4 drift where the real fault is a 1000 km miss:
    assert(inDef.missM <= tolM, ...
        ['the solve converged but missed by %.2f m against a %.1f m ' ...
         'tolerance'],inDef.missM,tolM);

%% THE HUMP. These four numbers describe the max-range arc, and everything
%% about the two-branch structure hangs off them: the loft angle that splits
%% the bracket, the range that bounds the envelope, and the apogee and flight
%% time the flown branch is MEASURED against in part 5:
    assertRel(inDef.loftStarDeg,25.068448737863108 ,1e-4,'max-range loft angle (deg)');
    assertRel(inDef.rngMaxM    ,5055301.6893285597 ,1e-4,'maximum range (m)');
    assertRel(inDef.rngMinM    ,556603.26462887647 ,1e-4,'envelope floor (m)');
    assertRel(inDef.hApoStarM  ,1021485.8049903549 ,1e-4,'max-range arc apogee (m)');
    assertRel(inDef.tFlyStarS  ,1327.8118532153449 ,1e-4,'max-range arc flight time (s)');
    assert(inDef.maxRange.nTurn == 1, ...
        ['the coarse scan of range against loft changed direction %d times; a ' ...
         'single interior hump changes direction exactly once, and the ' ...
         'golden-section search below it is only valid on one'], ...
        inDef.maxRange.nTurn);
    assert(inDef.maxRange.nEval == inDef.maxRange.nScan + inDef.maxRange.nGolden + 3, ...
        ['the bracketing reported %d propagations for %d scan points and %d ' ...
         'golden steps; it must cost the scan, two initial interior points, ' ...
         'one per step and one final evaluation at the located maximum'], ...
        inDef.maxRange.nEval,inDef.maxRange.nScan,inDef.maxRange.nGolden);

%% THE HEADLINE NUMBERS of the flown arc, to 1e-4 relative, read at full
%% printed precision from a run of the committed code:
    assertRel(inDef.rngReqM  ,3174980.9089995562 ,1e-4,'required range (m)');
    assertRel(inDef.rngAchM  ,3174865.0504538282 ,1e-4,'achieved range (m)');
    assertRel(inDef.missM    ,115.85854572764042 ,1e-4,'miss distance (m)');
    assertRel(inDef.residM   ,-115.85854572802782,1e-4,'range residual (m)');
    assertRel(inDef.loftDeg  ,59.145739482374033 ,1e-4,'flown loft angle (deg)');
    assertRel(inDef.psiLaunch,0.68412945715650275,1e-4,'launch azimuth (rad)');
    assertRel(inDef.tFlight  ,1816.0099756049399 ,1e-4,'flight time (s)');
    assertRel(inDef.hApogee  ,2118448.9231126327 ,1e-4,'apogee (m)');
    assertRel(inDef.tApogee  ,936.83532260709865 ,1e-4,'time of apogee (s)');
    assertRel(inDef.vImpact  ,3502.7890627171873 ,1e-4,'impact speed (m/s)');
    assertRel(rad2deg(inDef.gamImpact),-62.741757597464222,1e-4,'impact angle (deg)');
    assertRel(inDef.tBurnout ,80.517757894736846 ,1e-4,'burnout time (s)');
    assertRel(inDef.hBurnout ,114579.20371766761 ,1e-4,'burnout altitude (m)');
    assertRel(inDef.vBurnout ,5741.8872569467831 ,1e-4,'burnout speed (m/s)');
    assertRel(rad2deg(inDef.gamBurnout),62.840462372477255,1e-4,'burnout gamma (deg)');
    assertRel(inDef.qBstMax  ,81489.445833564532 ,1e-4,'boost peak dynamic pressure (Pa)');
    assertRel(inDef.nBstMax  ,40.363728442775539 ,1e-4,'peak boost sensed load (g)');
    assertRel(inDef.qReMax   ,7523889.6014499972 ,1e-4,'re-entry peak dynamic pressure (Pa)');
    assertRel(inDef.nReMax   ,82.066669559692102 ,1e-4,'peak re-entry deceleration (g)');

%% The impact point itself, pinned separately from the range. A single scalar
%% range is built from both coordinates and cannot tell them apart:
    assertRel(rad2deg(inDef.latImpact), 61.999681084883072,1e-4,'impact latitude (deg)');
    assertRel(rad2deg(inDef.lonImpact),-60.002110254674157,1e-4,'impact longitude (deg)');

%% Cross-range is zero in this configuration, so the whole miss IS the range
%% residual. Not a general guarantee: it holds because the bank angle is zero:
    assertAbs(inDef.xTrackM,0,1e-3,'cross-track offset (m)');
    assertAbs(inDef.missM,abs(inDef.residM),1e-6, ...
        'the miss against the magnitude of the range residual (m)');

%% ---------------------------------------------------------------------
%% 2. The impact point, checked INDEPENDENTLY from the flown trajectory
%% ---------------------------------------------------------------------
%% Not from info, not from the summary: from traj.x(end,:), through
%% coorbital.util.greatCircle, in the documented (lat1,lon1,lat2,lon2) order.
%% This is the only assertion in the file that does not trust the script's own
%% bookkeeping at all:
          missIndep = c.rE.*coorbital.util.greatCircle(trDef.x(end,3), ...
                                                       trDef.x(end,2), ...
                                                       latTg,lonTg);
    assert(missIndep <= tolM, ...
        ['the flown impact point is %.2f m from the target, outside the ' ...
         '%.1f m tolerance, whatever the summary reported'],missIndep,tolM);
    assertAbs(missIndep,inDef.missM,1e-6, ...
        'independently measured miss against the reported one (m)');

%% ...and the geometry is DISCRIMINATING, so the assertion above is shown to be
%% able to fail rather than assumed to be. Every way of transposing the
%% arguments must move the RANGE by far more than the miss assertion's budget;
%% if a future change to the shipped points ever brings one of them close, this
%% guard fails loudly instead of the check above going quietly blind:
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
%% 3. The azimuth, to machine precision, and that it too is discriminating
%% ---------------------------------------------------------------------
%% The launch heading state must BE the closed-form bearing. Nothing integrates
%% the initial condition, so this is an equality and not a tolerance -- 1e-13
%% rad is 6e-7 m of arc at Earth radius, far below any physical scale here:
            psiRef = coorbital.util.greatCircleBearing(latLn,lonLn,latTg,lonTg);
    assertAbs(trDef.x(1,6),psiRef,1e-13, ...
        'flown initial heading against the bearing (rad)');
    assertAbs(inDef.psiLaunch,psiRef,1e-13,'reported launch azimuth (rad)');

%% ...and the shipped pair must distinguish the reverse course and the
%% transposed one from the true one, or an azimuth defect could hide in it:
            psiRev = coorbital.util.greatCircleBearing(latTg,lonTg,latLn,lonLn);
            psiTrn = coorbital.util.greatCircleBearing(lonLn,latLn,lonTg,latTg);
    assert(abs(rad2deg(wrapPi(psiRev - psiRef))) > 5, ...
        'the reverse bearing is only %.3f deg from the forward one', ...
        abs(rad2deg(wrapPi(psiRev - psiRef))));
    assert(abs(rad2deg(wrapPi(psiTrn - psiRef))) > 5, ...
        'a lat/lon transposed bearing is only %.3f deg from the true one', ...
        abs(rad2deg(wrapPi(psiTrn - psiRef))));

%% ---------------------------------------------------------------------
%% 4. The summary must REPORT what was computed, BOTH branches included
%% ---------------------------------------------------------------------
%% A number found correctly and then printed from the wrong variable is
%% invisible to every assertion above, and the printed summary is the only part
%% of this script most users ever read. The two-branch table is checked here
%% because reporting BOTH arcs is a requirement in its own right: only one is
%% flown, and the other is the trade the user is being asked to make:
    assertAbs(summaryNumber(outDef,'required range +([-\d.]+) +km'), ...
        inDef.rngReqM./1000,1e-3,'reported required range (km)');
    assertAbs(summaryNumber(outDef,'achieved range +([-\d.]+) +km'), ...
        inDef.rngAchM./1000,1e-3,'reported achieved range (km)');
    assertAbs(summaryNumber(outDef,'MISS DISTANCE +([-\d.]+) +m'), ...
        inDef.missM,0.01,'reported miss distance (m)');
    assertAbs(summaryNumber(outDef,'max-range loft +([-\d.]+) +deg'), ...
        inDef.loftStarDeg,1e-3,'reported max-range loft angle (deg)');
    assertAbs(summaryNumber(outDef,'MAXIMUM RANGE +([-\d.]+) +km'), ...
        inDef.rngMaxM./1000,1e-3,'reported maximum range (km)');
    assertAbs(summaryNumber(outDef,'reachable +([-\d.]+) to'), ...
        inDef.rngMinM./1000,1e-3,'reported envelope floor (km)');
    assertAbs(summaryNumber(outDef,'reachable +[-\d.]+ to ([-\d.]+) km'), ...
        inDef.rngMaxM./1000,1e-3,'reported envelope ceiling (km)');
    assertAbs(summaryNumber(outDef,'launch azimuth +([-\d.]+) +deg'), ...
        rad2deg(psiRef),1e-5,'reported launch azimuth (deg)');

%% The four quantities the brief asks for on BOTH arcs, parsed out of the
%% table as two columns each. Blank the tolerance argument's ambition here: what
%% is being checked is that the printed pair matches the computed pair, to the
%% resolution the table prints at:
    assertPair(outDef,'APOGEE \(km\) +([-\d.]+) +([-\d.]+)', ...
        inDef.depressed.hApoM./1000,inDef.lofted.hApoM./1000,1e-3,'apogee (km)');
    assertPair(outDef,'FLIGHT TIME \(s\) +([-\d.]+) +([-\d.]+)', ...
        inDef.depressed.tFlyS,inDef.lofted.tFlyS,1e-3,'flight time (s)');
    assertPair(outDef,'IMPACT SPEED \(m/s\) +([-\d.]+) +([-\d.]+)', ...
        inDef.depressed.vImpM,inDef.lofted.vImpM,1e-3,'impact speed (m/s)');
    assertPair(outDef,'IMPACT ANGLE \(deg\) +([-\d.]+) +([-\d.]+)', ...
        rad2deg(inDef.depressed.gamImR),rad2deg(inDef.lofted.gamImR),1e-3, ...
        'impact flight-path angle (deg)');

%% THE THINGS A READER MUST BE TOLD, and which a future edit could quietly
%% drop. The branch structure and the minimum-energy caveat are not decoration:
%% a user who reads "minimum-energy" and assumes the textbook meaning has been
%% misled, and the script says so at length precisely so it cannot be:
    assert(contains(outDef,'THE TWO ARCS THAT REACH THIS TARGET'), ...
        'the summary must report BOTH branch solutions, not only the flown one');
    assert(contains(outDef,'BRANCH SELECTION'), ...
        'the summary must say which branch was chosen and on what grounds');
    assert(contains(outDef,'MEASURED as'), ...
        'the summary must report the branch MEASURED from the flown arc');
    assert(contains(outDef,'WHAT MINIMUM-ENERGY MEANS HERE'), ...
        ['the minimum-energy default must state precisely what it implements; ' ...
         'the textbook meaning is not available to a fixed booster']);
    assert(contains(outDef,'NON-ROTATING EARTH'), ...
        'the summary must state that the azimuth is exact only for a non-rotating Earth');
    assert(contains(outDef,'COMMANDED ATTITUDE, NOT THE BURNOUT GAMMA'), ...
        ['the summary must warn that the loft angle is a command and not the ' ...
         'burnout flight-path angle; on this configuration they differ by ' ...
         'tens of degrees']);

%% ---------------------------------------------------------------------
%% 5. ALL THREE BRANCH SELECTIONS, and the branch measured rather than labelled
%% ---------------------------------------------------------------------
             outLof = evalc(['[trLof,inLof] = run_ballistic_target(struct(' ...
                             '''branch'',''lofted'',''showPlots'',false));']);
             outDep = evalc(['[trDep,inDep] = run_ballistic_target(struct(' ...
                             '''branch'',''depressed'',''showPlots'',false));']);
    assert(contains(outLof,'(nominal)') && ~contains(outLof,'REFUSED'), ...
        'the lofted request did not solve. Summary was:\n%s',outLof);
    assert(contains(outDep,'(nominal)') && ~contains(outDep,'REFUSED'), ...
        'the depressed request did not solve. Summary was:\n%s',outDep);

%% Each run must have FLOWN the branch it was asked for, and must MEASURE as
%% that branch too. The second is the one that matters: the first is a label
%% the script writes about itself:
    assert(strcmp(inLof.branchFlown,'lofted') && ...
           strcmp(inLof.branchMeasured,'lofted') && inLof.branchAgrees, ...
        ['the lofted request flew "%s" and measured as "%s"'], ...
        inLof.branchFlown,inLof.branchMeasured);
    assert(strcmp(inDep.branchFlown,'depressed') && ...
           strcmp(inDep.branchMeasured,'depressed') && inDep.branchAgrees, ...
        ['the depressed request flew "%s" and measured as "%s"'], ...
        inDep.branchFlown,inDep.branchMeasured);
    assert(inLof.branchTimeAgrees && inDep.branchTimeAgrees, ...
        'the apogee and the flight time disagreed about the branch');

%% ...and the measurement is REDONE HERE, from the flown trajectories and the
%% max-range arc's own apogee and flight time, so that a measureBranch wired to
%% a constant is caught rather than believed:
             hApoL  = maxAlt(trLof,c.rE);
             hApoD  = maxAlt(trDep,c.rE);
    assert(hApoL > inDef.hApoStarM && trLof.t(end) > inDef.tFlyStarS, ...
        ['the "lofted" arc apogees at %.3f km in %.1f s against the max-range ' ...
         'arc''s %.3f km in %.1f s; a lofted arc must exceed BOTH'], ...
        hApoL./1000,trLof.t(end),inDef.hApoStarM./1000,inDef.tFlyStarS);
    assert(hApoD < inDef.hApoStarM && trDep.t(end) < inDef.tFlyStarS, ...
        ['the "depressed" arc apogees at %.3f km in %.1f s against the ' ...
         'max-range arc''s %.3f km in %.1f s; a depressed arc must fall short ' ...
         'of BOTH'],hApoD./1000,trDep.t(end),inDef.hApoStarM./1000,inDef.tFlyStarS);

%% THE TWO BRANCHES MUST BE GENUINELY DIFFERENT TRAJECTORIES. A script that
%% solved one branch twice and labelled the copies differently would satisfy
%% every miss assertion in this file. Ten times the apogee and nearly three
%% times the flight time separate them here; the thresholds are deliberately
%% far below the measured margins so that a real change of geometry does not
%% trip them while a collapsed pair cannot pass:
    assert(hApoL > 5.*hApoD, ...
        ['the lofted and depressed arcs apogee at %.3f km and %.3f km, a ' ...
         'factor of only %.2f. They are not two different branches'], ...
        hApoL./1000,hApoD./1000,hApoL./hApoD);
    assert(trLof.t(end) > 2.*trDep.t(end), ...
        ['the lofted and depressed arcs fly for %.1f s and %.1f s, a factor of ' ...
         'only %.2f'],trLof.t(end),trDep.t(end),trLof.t(end)./trDep.t(end));
    assert(inLof.loftDeg > inDef.loftStarDeg && inDep.loftDeg < inDef.loftStarDeg, ...
        ['the two solved loft angles are %.4f and %.4f deg; they must lie on ' ...
         'OPPOSITE sides of the %.4f deg max-range angle'], ...
        inLof.loftDeg,inDep.loftDeg,inDef.loftStarDeg);

%% Both arcs hit the target, and both are verified independently from the flown
%% state rather than from the solver's bookkeeping:
    assert(inLof.missM <= tolM && inDep.missM <= tolM, ...
        'the lofted arc missed by %.2f m and the depressed by %.2f m against %.1f m', ...
        inLof.missM,inDep.missM,tolM);
             missIL = c.rE.*coorbital.util.greatCircle(trLof.x(end,3),trLof.x(end,2), ...
                                                       latTg,lonTg);
             missID = c.rE.*coorbital.util.greatCircle(trDep.x(end,3),trDep.x(end,2), ...
                                                       latTg,lonTg);
    assert(missIL <= tolM && missID <= tolM, ...
        ['independently measured: the lofted impact is %.2f m from the target ' ...
         'and the depressed %.2f m, against a %.1f m tolerance'], ...
        missIL,missID,tolM);
    assertAbs(missIL,inLof.missM,1e-6,'independent lofted miss (m)');
    assertAbs(missID,inDep.missM,1e-6,'independent depressed miss (m)');

%% The depressed arc's own headline numbers, pinned. It is the branch the
%% shipped run does NOT fly, so nothing in part 1 pins it:
    assertRel(inDep.loftDeg  ,-14.360235287375231,1e-4,'depressed loft angle (deg)');
    assertRel(inDep.rngAchM  ,3174314.8635772206 ,1e-4,'depressed achieved range (m)');
    assertRel(inDep.missM    ,666.04542233526706 ,1e-4,'depressed miss (m)');
    assertRel(inDep.hApogee  ,208357.41491637472 ,1e-4,'depressed apogee (m)');
    assertRel(inDep.tFlight  ,659.00850038788599 ,1e-4,'depressed flight time (s)');
    assertRel(inDep.vImpact  ,924.200491974653   ,1e-4,'depressed impact speed (m/s)');
    assertRel(rad2deg(inDep.gamImpact),-18.226928729261608,1e-4, ...
        'depressed impact angle (deg)');
    assertRel(rad2deg(inDep.latImpact), 61.998166186475828,1e-4, ...
        'depressed impact latitude (deg)');
    assertRel(rad2deg(inDep.lonImpact),-60.012130788730545,1e-4, ...
        'depressed impact longitude (deg)');

%% THE STEEPER ARC IS THE FASTER ONE HERE, which is the opposite of the vacuum
%% intuition and is a DRAG result: the shallow depressed descent spends far
%% longer in dense air. Pinned because the summary states it in words, and a
%% statement in words that no assertion checks is a statement that can go stale:
    assert(inLof.vImpact > 3.*inDep.vImpact, ...
        ['the lofted arc arrives at %.1f m/s and the depressed at %.1f m/s; ' ...
         'the summary claims the steep arc arrives much the faster of the two'], ...
        inLof.vImpact,inDep.vImpact);
    assert(inLof.gamImpact < inDep.gamImpact, ...
        'the lofted arc must arrive STEEPER: %.3f deg against %.3f deg', ...
        rad2deg(inLof.gamImpact),rad2deg(inDep.gamImpact));

%% ---------------------------------------------------------------------
%% 6. THE SELECTOR IS NOT A CONSTANT
%% ---------------------------------------------------------------------
%% minimum-energy took the LOFTED arc above. At a longer target the ranking of
%% the two loft-angle distances flips and it must take the DEPRESSED one. A
%% selector hard-wired either way passes exactly one of these two:
    assert(strcmp(inDef.branchFlown,'lofted') && ...
           inDef.dLoftLofDeg < inDef.dLoftDepDeg, ...
        ['at the shipped target minimum-energy must take the lofted arc: it ' ...
         'lies %.4f deg from the max-range angle against the depressed arc''s ' ...
         '%.4f deg, and it flew "%s"'], ...
        inDef.dLoftLofDeg,inDef.dLoftDepDeg,inDef.branchFlown);
             outFlp = evalc(['[trFlp,inFlp] = run_ballistic_target(struct(' ...
                             '''lonTarget'',-40,''showPlots'',false));']);
    assert(contains(outFlp,'(nominal)') && ~contains(outFlp,'REFUSED'), ...
        'the 4218 km minimum-energy case did not solve. Summary was:\n%s',outFlp);
    assert(inFlp.rngReqM > inDef.rngReqM + 1e6, ...
        'the second minimum-energy target is not materially further; it cannot discriminate');
    assert(strcmp(inFlp.branchFlown,'depressed') && ...
           strcmp(inFlp.branchMeasured,'depressed'), ...
        ['at 4218 km the depressed arc lies %.4f deg from the max-range angle ' ...
         'against the lofted arc''s %.4f deg, so minimum-energy must take it; ' ...
         'it flew "%s" and measured as "%s". A selector wired to a constant ' ...
         'passes the shipped case and fails here'], ...
        inFlp.dLoftDepDeg,inFlp.dLoftLofDeg,inFlp.branchFlown,inFlp.branchMeasured);
    assert(inFlp.dLoftDepDeg < inFlp.dLoftLofDeg, ...
        'the flip case does not actually flip: %.4f against %.4f deg', ...
        inFlp.dLoftDepDeg,inFlp.dLoftLofDeg);
    assert(inFlp.missM <= tolM, ...
        'the flip case missed by %.2f m',inFlp.missM);
             missIF = c.rE.*coorbital.util.greatCircle(trFlp.x(end,3),trFlp.x(end,2), ...
                                                       deg2rad(62),deg2rad(-40));
    assert(missIF <= tolM, ...
        'the flip case''s flown impact point is %.2f m from its target',missIF);
    assert(maxAlt(trFlp,c.rE) < inFlp.hApoStarM, ...
        'the flip case measured as depressed but apogees above the max-range arc');

%% ...and the SOLVE MOVED with the target. A longer target on a different
%% azimuth must give a different loft angle and a different bearing:
    assert(abs(inFlp.loftDeg - inDef.loftDeg) > 1, ...
        'the two minimum-energy targets solved to the same loft angle');
    assert(abs(rad2deg(wrapPi(inFlp.psiLaunch - inDef.psiLaunch))) > 1, ...
        'the two targets returned the same azimuth to within a degree');

%% ---------------------------------------------------------------------
%% 7. THE ILL-CONDITIONED BAND, where the selection is not worth trusting
%% ---------------------------------------------------------------------
%% Between those two targets the two loft distances pass through equality, and
%% the script CAUTIONS there rather than pretending the choice is meaningful.
%% Unpinned, that caution would be a paragraph nothing can fire:
             outCau = evalc(['[~,inCau] = run_ballistic_target(struct(' ...
                             '''lonTarget'',-48,''showPlots'',false));']);
    assert(inCau.meClose, ...
        ['at 3789 km the two loft distances are %.4f and %.4f deg, a ratio of ' ...
         '%.4f, and info.meClose must be true inside the 5 %% band'], ...
        inCau.dLoftDepDeg,inCau.dLoftLofDeg,inCau.meRatio);
    assertRel(inCau.meRatio,0.95856292689990052,1e-3,'ill-conditioning ratio');
    assert(contains(outCau,'turns on a rounding'), ...
        ['the near-degenerate case printed no caution. Summary was:\n%s'],outCau);
    assert(~inDef.meClose && ~contains(outDef,'turns on a rounding'), ...
        ['the shipped case is %.2f %% from degenerate and must NOT print the ' ...
         'caution, or the caution says nothing'],100.*(1 - inDef.meRatio));

%% ---------------------------------------------------------------------
%% 8. THE REFUSALS: beyond maximum range, too close, and a missing branch
%% ---------------------------------------------------------------------
%% An unreachable target must be REFUSED, in words, with the band and with the
%% maximum range and the loft angle that achieves it -- and must NOT throw.
%% Throwing would force every caller into a try/catch to read numbers the
%% returned struct already carries, and worse, a targeting script that silently
%% returned the nearest miss instead would hand back something that LOOKS like
%% a solution:
[outFar,infFar,threwFar,trFar] = tryRun(struct('latTarget',0,'lonTarget',20, ...
                                               'showPlots',false));
    assert(~threwFar,'a target beyond maximum range must not throw:\n%s',outFar);
    assert(contains(outFar,'REFUSED'),'no refusal banner for a too-far target:\n%s',outFar);
    assert(contains(outFar,'BEYOND MAXIMUM RANGE'), ...
        'the refusal must say the target is beyond maximum range:\n%s',outFar);
    assert(infFar.refused,'info.refused must be true on a refusal');
    assert(isempty(trFar),'a refused run must return an empty trajectory');
    assert(~infFar.depressed.exists && ~infFar.lofted.exists, ...
        'beyond maximum range NEITHER branch may report that it reaches the target');
    assert(strcmp(infFar.depressed.solveInfo.identifier, ...
                  'coorbital:rangeSolve:targetOutsideBracket'), ...
        'expected the outside-bracket identifier, got "%s"', ...
        infFar.depressed.solveInfo.identifier);
    assert(infFar.rngReqM > infFar.rngMaxM, ...
        'the too-far case must have a required range above the maximum');

%% THE MAXIMUM, THE LOFT ANGLE ACHIEVING IT AND THE SHORTFALL must all be
%% PRINTED, not merely returned, because a user acting on a refusal is reading
%% the page and not the struct:
    assertAbs(summaryNumber(outFar,'MAXIMUM RANGE +([-\d.]+) +km'), ...
        infFar.rngMaxM./1000,0.01,'refusal maximum range (km)');
    assertAbs(summaryNumber(outFar,'angle of ([-\d.]+) deg'), ...
        infFar.loftStarDeg,1e-3,'refusal max-range loft angle (deg)');
    assertAbs(summaryNumber(outFar,'shortfall +\+?([-\d.]+) +km'), ...
        infFar.shortfallM./1000,0.01,'refusal shortfall (km)');
    assertAbs(summaryNumber(outFar,'required range +([-\d.]+) +km'), ...
        infFar.rngReqM./1000,0.01,'refusal required range (km)');

%% THE BAND, both edges, parsed out of the sentence that states it:
            bandTok = regexp(outFar,'([\d.]+) to ([\d.]+) km reachable band', ...
                             'tokens','once');
    assert(numel(bandTok) == 2, ...
        'the refusal does not print the reachable band:\n%s',outFar);
    assertAbs(str2double(bandTok{1}),infFar.rngMinM./1000,0.01, ...
        'printed band floor (km)');
    assertAbs(str2double(bandTok{2}),infFar.rngMaxM./1000,0.01, ...
        'printed band ceiling (km)');
    assertRel(infFar.rngMaxM,5055301.6893285597,1e-4,'refused-case maximum range (m)');

%% ...and the same at the near end, which is a different branch of the message
%% and a different thing to tell the user to change:
   [outNear,infNear,threwNear] = tryRun(struct('latTarget',47,'lonTarget',-98, ...
                                               'showPlots',false));
    assert(~threwNear,'a target inside the envelope floor must not throw:\n%s',outNear);
    assert(contains(outNear,'REFUSED') && contains(outNear,'TOO CLOSE'), ...
        'the refusal must say the target is too close:\n%s',outNear);
    assert(infNear.refused && infNear.rngReqM < infNear.rngMinM, ...
        'the too-close case must have a required range below the band floor');

%% A THIRD REFUSAL WITH NO COUNTERPART IN HGV/run_target: the target is inside
%% the envelope, but only ONE branch reaches it, and the other one was asked
%% for. Silently flying the branch that does reach it would be the worst
%% outcome -- it converges, it hits, and it is not the trajectory that was
%% asked for:
  [outBr,infBr,threwBr,trBr] = tryRun(struct('latTarget',54,'lonTarget',-92, ...
                                             'branch','depressed','showPlots',false));
    assert(~threwBr,'an unreachable branch must not throw:\n%s',outBr);
    assert(contains(outBr,'REFUSED'),'no refusal banner for an unreachable branch:\n%s',outBr);
    assert(infBr.refused && isempty(trBr), ...
        'the unreachable-branch case must refuse and return an empty trajectory');
    assert(~infBr.depressed.exists && infBr.lofted.exists, ...
        ['this case exists to test the ONE-branch refusal: the depressed arc ' ...
         'must not reach %.3f km while the lofted one does'],infBr.rngReqM./1000);
    assert(infBr.rngReqM > infBr.rngMinM && infBr.rngReqM < infBr.rngMaxM, ...
        'the unreachable-branch target must sit INSIDE the overall envelope');
    assert(contains(outBr,'though the other one does'), ...
        'the refusal must say that the other branch reaches the target:\n%s',outBr);

%% ---------------------------------------------------------------------
%% 9. THE GUARDS: the alphaMax clamp, the branch word, the override hook
%% ---------------------------------------------------------------------
%% AT BM/run_ballistic's 6 deg CLAMP THERE IS NO DEPRESSED BRANCH AT ALL. The
%% vehicle cannot be pushed past the max-range attitude at any commanded loft
%% angle, so range is monotone across the whole bracket and the maximum sits on
%% an endpoint. That must be REFUSED with a named identifier, not reported as a
%% two-branch problem with one branch quietly missing:
             idClp = errorIdOf(@() run_ballistic_target( ...
                                    struct('alphaMax',6,'showPlots',false)));
    assert(strcmp(idClp,'coorbital:runBallisticTarget:maximumNotBracketed'), ...
        ['a 6 deg angle-of-attack clamp leaves the max-range loft angle ' ...
         'outside the bracket and there is then only one branch; expected ' ...
         'coorbital:runBallisticTarget:maximumNotBracketed, got "%s"'],idClp);

%% A misspelt branch word must raise rather than fall through to a default.
%% There is no sensible default: every range short of the maximum is reached by
%% two different trajectories:
             idBad = errorIdOf(@() run_ballistic_target( ...
                                    struct('branch','lofty','showPlots',false)));
    assert(~isempty(idBad),'an unrecognised branch word must raise, not be ignored');

%% ...and so must a misspelt override field, or a future test could believe it
%% flew a case it did not:
             idOvr = errorIdOf(@() run_ballistic_target(struct('latTargetDeg',35)));
    assert(~isempty(idOvr),'an unrecognised override field must raise, not be ignored');

%% ---------------------------------------------------------------------
%% 10. THE DISPLAY CHOICE, which the two branches make necessary
%% ---------------------------------------------------------------------
%% A fixed exaggeration cannot serve both arcs: they differ by a factor of ten
%% in apogee, and at the 3x BM/run_ballistic hard-codes the depressed arc is
%% invisible. Pinned as the two literals they are, with the apogees they come
%% from, and the rule's cap and floor asked of the rule directly because no
%% flyable configuration reaches either:
    assert(inDef.altExag == 2, ...
        'the lofted run drew its globe at %gx against the 2x its %.3f km apogee earns', ...
        inDef.altExag,inDef.hApogee./1000);
    assert(inDep.altExag == 9, ...
        'the depressed run drew its globe at %gx against the 9x its %.3f km apogee earns', ...
        inDep.altExag,inDep.hApogee./1000);
    assert(inDef.altExag ~= inDep.altExag, ...
        'the two branches drew at the same exaggeration; the rule is not adapting');
    assert(inDef.altExagRule(40e3) == 30 && inDef.altExagRule(0) == 30, ...
        'the exaggeration rule must CAP at 30, including at a zero apogee');
    assert(inDef.altExagRule(2e6) == 2, ...
        'the exaggeration rule must FLOOR at 2; coorbital.viz refuses 0 or 1');
    assert(inDef.altExag == inDef.altExagRule(maxAlt(trDef,c.rE)), ...
        'the reported exaggeration is not the one its own rule gives for the flown apogee');

%% The captions carry hemisphere letters. Both shipped coordinates are WEST of
%% the meridian, so a formatter that never reached its west branch would produce
%% a caption wrong by 200 degrees of longitude and looking perfectly reasonable:
    assert(strcmp(inDef.launchStr,'45.00N 100.00W'), ...
        'the launch caption reads "%s" against the expected "45.00N 100.00W"', ...
        inDef.launchStr);
    assert(strcmp(inDef.targetStr,'62.00N 60.00W'), ...
        'the target caption reads "%s" against the expected "62.00N 60.00W"', ...
        inDef.targetStr);

%% Nothing above asked for a figure, so nothing above may have made one:
            nFig1  = numel(findall(groot,'Type','figure'));
    assert(nFig1 == nFig0, ...
        ['run_ballistic_target left %d new figure(s) open with showPlots ' ...
         'false; the suite must run headless'],nFig1 - nFig0);
end

function h = maxAlt(traj,rE)
%% Purpose:
%
%  Peak altitude of a trajectory above the reference sphere, taken from the
%  state history itself rather than from the reported apogee. Used so that the
%  branch measurement in part 5 does not read back the very number it is
%  checking.
%
%% Inputs:
%
%  traj             Struct                      Trajectory from phaseRun
%
%  rE               [1 x 1]                     Reference sphere radius (m)
%
%% Outputs:
%
%  h                [1 x 1]                     Peak altitude (m)
%
%% Revision History:
%  Michael Casey                                                08/07/2026
%  Copyright 2026 Coorbital, Inc.
%% ------------------------ Begin Code Sequence ---------------------------

                 h = max(traj.x(:,1)) - rE;
end

function id = errorIdOf(fh)
%% Purpose:
%
%  Run a handle and return the identifier of whatever it raised, or an empty
%  char if it did not raise at all. Exists so that a guard can be checked for
%  raising the RIGHT error rather than merely for raising: a test that only
%  asks "did it throw" passes on a typo in the code it is guarding.
%
%% Inputs:
%
%  fh               Function handle             Called with no arguments
%
%% Outputs:
%
%  id               Char [1 x n]                Error identifier, '' if the
%                                               call returned normally
%
%% Revision History:
%  Michael Casey                                                08/07/2026
%  Copyright 2026 Coorbital, Inc.
%% ------------------------ Begin Code Sequence ---------------------------

                id = '';
    try
        fh();
    catch err
                id = err.identifier;
        if isempty(id)
                id = 'unidentified:error';
        end
    end
end

function [out,info,threw,traj] = tryRun(optStruct)
%% Purpose:
%
%  Run run_ballistic_target with the given overrides, capturing its printed
%  output and recording whether it threw. Exists because the refusal paths must
%  be shown NOT to throw, and a bare call inside the test would abort the test
%  file on the very behaviour being denied.
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
%  info             Struct                      The script's info output, or an
%                                               empty struct if it threw
%
%  threw            [1 x 1] logical             True if the script raised
%
%  traj             Struct or []                The script's trajectory output
%
%% Revision History:
%  Michael Casey                                                08/07/2026
%  Copyright 2026 Coorbital, Inc.
%% ------------------------ Begin Code Sequence ---------------------------

    assert(isstruct(optStruct),'tryRun takes a struct of USER PARAMETERS overrides.');
              info = struct();
              traj = [];
             threw = false;

%% evalc runs in THIS workspace, so optStruct is what the script receives; the
%% assertion above is also what tells the code analyser that, since it cannot
%% see a variable used only inside an evaluated string:
    try
               out = evalc('[traj,info] = run_ballistic_target(optStruct);');
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
%  a                [1 x 1]                     The same angle in (-pi,pi] (rad)
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
%  Pull one number out of a captured summary. Parsing the printed text rather
%  than recomputing the quantity is deliberate: what a reader of this library
%  acts on is the summary, so the summary is what gets asserted.
%
%  The patterns use explicit space runs rather than \s, because in MATLAB \s
%  also matches a newline and would happily reach across into the next line of
%  the report.
%
%% Inputs:
%
%  txt              Char [1 x n]                Captured summary
%
%  pat              Char [1 x m]                Regular expression with exactly
%                                               one capturing group around the
%                                               number
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

function assertPair(txt,pat,wantA,wantB,tol,what)
%% Purpose:
%
%  Assert one row of the printed two-branch table against the two values it is
%  supposed to be reporting. Exists because the branch table is the deliverable
%  of this script -- the user is being shown a trade -- and a table whose two
%  columns came from one branch would be invisible to every other assertion.
%
%% Inputs:
%
%  txt              Char [1 x n]                Captured summary
%
%  pat              Char [1 x m]                Regular expression with exactly
%                                               TWO capturing groups, depressed
%                                               column first
%
%  wantA            [1 x 1]                     Expected depressed-column value
%
%  wantB            [1 x 1]                     Expected lofted-column value
%
%  tol              [1 x 1]                     Absolute tolerance, in the
%                                               units of the quantity
%
%  what             Char [1 x k]                Name of the row, for the message
%
%% Outputs:
%
%  none                                         Throws on failure
%
%% Revision History:
%  Michael Casey                                                08/07/2026
%  Copyright 2026 Coorbital, Inc.
%% ------------------------ Begin Code Sequence ---------------------------

               tok = regexp(txt,pat,'tokens','once');
    assert(numel(tok) == 2, ...
        'the two-branch table has no "%s" row matching "%s":\n%s',what,pat,txt);
    assertAbs(str2double(tok{1}),wantA,tol,['printed depressed ' what]);
    assertAbs(str2double(tok{2}),wantB,tol,['printed lofted ' what]);
    assert(abs(wantA - wantB) > tol, ...
        ['the two columns of the "%s" row are %.6g and %.6g, closer than the ' ...
         '%.3g tolerance this row is checked to. It cannot tell one branch ' ...
         'printed twice from two branches printed once each'],what,wantA,wantB,tol);
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
%  reference is zero, or where the printed resolution rather than the value sets
%  the budget, so a relative test would be meaningless.
%
%% Inputs:
%
%  got              [1 x 1]                     Measured value
%
%  want             [1 x 1]                     Reference value
%
%  tol              [1 x 1]                     Absolute tolerance, in the units
%                                               of the quantity
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
