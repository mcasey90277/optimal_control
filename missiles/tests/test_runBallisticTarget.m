function test_runBallisticTarget()
%% Purpose:
%
%  Pin BM/run_ballistic_target, the ballistic point-to-point targeting script:
%  the closed-form launch azimuth, the bracketing of the max-range loft angle,
%  the TWO full-burn branch solves either side of it, the TEXTBOOK
%  MINIMUM-ENERGY solve on loft angle and cutoff fraction together, the
%  after-the-fact measurement of which branch was actually flown, and the
%  refusals at both ends of the envelope, for an unreachable branch and for a
%  minimum-energy request the vehicle cannot meet.
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
%    THE MINIMUM-ENERGY MODE IS GRADED AGAINST A CLOSED FORM THIS FILE COMPUTES
%    ITSELF, in parts 6, 7 and 13. It is the one mode whose answer is known
%    independently of the code: the classical minimum-energy trajectory between
%    two points at the same radius has V*^2 = (mu/rE) 2 sin(psi/2)/(1 +
%    sin(psi/2)) and gammaStar = 45 deg - psi/4, and classicalArc below writes
%    that out rather than calling the script's own. Mutate gammaStar in the
%    script to a flat 45 deg and this file still knows the right answer.
%
%    ...AND IT IS SHOWN NOT TO BE EITHER FULL-BURN ARC. Every range and miss
%    assertion here is satisfied by the lofted arc and by the depressed one, so
%    a minimum-energy mode that quietly fell back to SELECTING a branch -- which
%    is exactly what it used to do -- would pass all of them.
%    assertBetweenArcs is what sees the difference, on two separate targets.
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
%    THE MEASUREMENT IS EXERCISED WHERE IT CAN FAIL, in part 11. Every run in
%    part 5 is one where the measurement AGREES with the label, so all of them
%    also pass if measureBranch is never called and the flown name is copied
%    from the commanded one. Part 11 flies the one geometry on which the two
%    can legitimately part company -- a target within the range tolerance of
%    the MAXIMUM, where both brackets converge at their shared endpoint and a
%    "lofted" request is served by the max-range arc itself -- and asserts the
%    disagreement, the caution and the printed line.
%
%  THE alphaMax GUARD IS PINNED TOO, in part 9, and it is not a formality: at
%  BM/run_ballistic's 6 deg clamp the shipped bracket has no interior maximum,
%  so one branch is empty everywhere on it and the script must refuse with
%  coorbital:runBallisticTarget:maximumNotBracketed rather than report a
%  one-branch problem as a two-branch one. Part 9 also pins the CAUSE, which
%  this file used to get wrong: widen loftMin to -140 deg at the same 6 deg
%  clamp and the hump IS there, at -42.907 deg with 5211.525 km of range. The
%  refusal is about bracket width. What 6 deg costs is reach -- its depressed
%  branch bottoms out at 4708.463 km -- and what 12 deg costs is 156.224 km of
%  maximum range.
%
%  COST. A trajectory propagation is about 0.05 s. A 'lofted' or 'depressed' run
%  takes 59 of them, about 3 s; a 'minimum-energy' run takes 249, about 11 s,
%  because the nested solve is 190 of those on its own. The whole file is about
%  80 s. The SHIPPED configuration is flown at its shipped settings throughout,
%  which is the stronger test: what is pinned is exactly what a user gets.
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
    assertRel(inDef.rngAchM  ,3175569.3035038970 ,1e-4,'achieved range (m)');
    assertRel(inDef.missM    ,588.39450434150390 ,1e-4,'miss distance (m)');
    assertRel(inDef.residM   ,588.39450434083120 ,1e-4,'range residual (m)');
    assertRel(inDef.loftDeg  ,28.601997036545610 ,1e-4,'flown loft angle (deg)');
    assertRel(inDef.psiLaunch,0.68412945715650275,1e-4,'launch azimuth (rad)');
    assertRel(inDef.tFlight  ,1013.3601833344570 ,1e-4,'flight time (s)');
    assertRel(inDef.hApogee  ,735074.48724318016 ,1e-4,'apogee (m)');
    assertRel(inDef.tApogee  ,531.25092599567850 ,1e-4,'time of apogee (s)');
    assertRel(inDef.vImpact  ,2372.2975669065300 ,1e-4,'impact speed (m/s)');
    assertRel(rad2deg(inDef.gamImpact),-38.888486066477910,1e-4,'impact angle (deg)');
    assertRel(inDef.tBurnout ,77.559280706465860 ,1e-4,'cutoff time (s)');
    assertRel(inDef.hBurnout ,80560.646773226560 ,1e-4,'burnout altitude (m)');
    assertRel(inDef.vBurnout ,4833.3136917076080 ,1e-4,'burnout speed (m/s)');
    assertRel(rad2deg(inDef.gamBurnout),37.865407738253940,1e-4,'burnout gamma (deg)');
    assertRel(inDef.qBstMax  ,84875.700372397520 ,1e-4,'boost peak dynamic pressure (Pa)');
    assertRel(inDef.nBstMax  ,27.657742256519640 ,1e-4,'peak boost sensed load (g)');
    assertRel(inDef.qReMax   ,3761799.6834535259 ,1e-4,'re-entry peak dynamic pressure (Pa)');
    assertRel(inDef.nReMax   ,41.031751916221470 ,1e-4,'peak re-entry deceleration (g)');

%% The impact point itself, pinned separately from the range. A single scalar
%% range is built from both coordinates and cannot tell them apart:
    assertRel(rad2deg(inDef.latImpact), 62.001619132623560,1e-4,'impact latitude (deg)');
    assertRel(rad2deg(inDef.lonImpact),-59.989282264861000,1e-4,'impact longitude (deg)');

%% Cross-range is zero in this configuration, so the whole miss IS the range
%% residual. Not a general guarantee: it holds because the bank angle is zero:
    assertAbs(inDef.xTrackM,0,1e-3,'cross-track offset (m)');
    assertAbs(inDef.missM,abs(inDef.residM),1e-6, ...
        'the miss against the magnitude of the range residual (m)');

%% THE PROPAGATION COUNT MUST AGREE WITH ITSELF, and it did not: the summary
%% printed 59 while info.nProp returned 57. The two propagations between them
%% are the ones that re-create each branch's state history at its converged loft
%% angle -- br.nProp and not br.nEval -- and a user comparing the printed cost
%% of a run against the returned one had no way to tell which was wrong:
    assertAbs(summaryNumber(outDef,'propagations +(\d+) +\(whole'), ...
        inDef.nProp,0.5,'printed propagation count against info.nProp');
    assert(inDef.nProp == inDef.maxRange.nEval + 1 + ...
                          inDef.depressed.nProp + inDef.lofted.nProp + ...
                          inDef.minEnergy.nProp, ...
        ['info.nProp is %d against %d bracketing propagations, one re-flying ' ...
         'the maximum, %d and %d on the two full-burn branches and %d on the ' ...
         'minimum-energy solve'], ...
        inDef.nProp,inDef.maxRange.nEval,inDef.depressed.nProp, ...
        inDef.lofted.nProp,inDef.minEnergy.nProp);

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
    assertPair(outDef,'    miss \(m\) +([-\d.]+) +([-\d.]+)', ...
        inDef.depressed.missM,inDef.lofted.missM,1e-3,'branch miss (m)');
    assertPair(outDef,'APOGEE \(km\) +([-\d.]+) +([-\d.]+)', ...
        inDef.depressed.hApoM./1000,inDef.lofted.hApoM./1000,1e-3,'apogee (km)');
    assertPair(outDef,'FLIGHT TIME \(s\) +([-\d.]+) +([-\d.]+)', ...
        inDef.depressed.tFlyS,inDef.lofted.tFlyS,1e-3,'flight time (s)');
    assertPair(outDef,'IMPACT SPEED \(m/s\) +([-\d.]+) +([-\d.]+)', ...
        inDef.depressed.vImpM,inDef.lofted.vImpM,1e-3,'impact speed (m/s)');
    assertPair(outDef,'IMPACT ANGLE \(deg\) +([-\d.]+) +([-\d.]+)', ...
        rad2deg(inDef.depressed.gamImR),rad2deg(inDef.lofted.gamImR),1e-3, ...
        'impact flight-path angle (deg)');

%% ...AND THE PRINTED MISS ROW MUST BE A MISS. The row above only checks that
%% the printed pair matches info, and a transposed greatCircle inside
%% solveBranch moves the printed number and the returned one TOGETHER: the
%% table then reads 1354.24 and 235.58 m against the true 666.05 and 115.86,
%% still says "reaches the target: yes", and nothing objects under a stated
%% 1000 m tolerance. This is what objects. Every branch that claims to reach the
%% target is held to the tolerance, on every run this file makes:
    assertBranchMiss(inDef,tolM,'the shipped minimum-energy run');

%% THE THINGS A READER MUST BE TOLD, and which a future edit could quietly
%% drop. The branch structure and the minimum-energy caveat are not decoration:
%% a user who reads "minimum-energy" and assumes the textbook meaning has been
%% misled, and the script says so at length precisely so it cannot be:
    assert(contains(outDef,'THE TWO ARCS THAT REACH THIS TARGET'), ...
        'the summary must report BOTH branch solutions, not only the flown one');
    assert(contains(outDef,'WHAT WAS FLOWN'), ...
        'the summary must say which trajectory was flown and on what grounds');
    assert(contains(outDef,'THE MINIMUM-ENERGY SOLVE'), ...
        ['the minimum-energy default must show its solve: two unknowns, two ' ...
         'residuals, and the classical arc it was solved against']);
    assert(contains(outDef,'AGAINST THE CLASSICAL ARC'), ...
        ['the minimum-energy default must be reported against the closed-form ' ...
         'reference for its range angle, not only against itself']);
    assert(contains(outDef,'PHYSICS, NOT SOLVER ERROR'), ...
        ['the summary must attribute the gap to the classical arc rather than ' ...
         'leaving a reader to read it as convergence error']);
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
    assertBranchMiss(inLof,tolM,'the lofted request');
    assertBranchMiss(inDep,tolM,'the depressed request');

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
%% 6. THE MINIMUM-ENERGY SOLVE, against the CLASSICAL closed form
%% ---------------------------------------------------------------------
%% THIS IS THE PART THE MODE EXISTS FOR. The shipped default no longer selects
%% between the two full-burn arcs -- neither of them is a minimum-energy
%% trajectory, because at full burn the booster delivers a fixed delta-V and
%% both arcs leave burnout with essentially the same energy, 0.55 %% apart. It
%% solves the LOFT ANGLE and the CUTOFF FRACTION together against two residuals:
%% the range, and the burnout flight-path angle against gammaStar.
%%
%% THE REFERENCE IS COMPUTED HERE, IN THIS FILE'S OWN HAND, from the range angle
%% and coorbital.util.missileConst. It is NOT read out of info.classical: a test
%% that took its reference from the code under test could not see the code
%% change the formula, which is exactly the mutation this part has to catch:
%%
%%     V*^2      = (mu/rE) * 2 sin(psi/2) / (1 + sin(psi/2))
%%     gammaStar = 45 deg - psi/4
%%
            psiReq = coorbital.util.greatCircle(latLn,lonLn,latTg,lonTg);
            gamStr = pi./4 - psiReq./4;
             refKm = classicalArc(psiReq,c);

%% The closed form first, to machine precision, because everything below is
%% measured against it and a wrong reference makes every tolerance meaningless.
%% 1e-9 relative, not 1e-4: there is no integration anywhere in it:
    assertRel(rad2deg(gamStr),37.869665306645022,1e-9,'gammaStar (deg)');
    assertRel(refKm.V        ,4970.2903871732401,1e-9,'classical burnout speed (m/s)');
    assertRel(refKm.hApoM    ,687304.31264742557,1e-9,'classical apogee (m)');
    assertRel(refKm.tofS     ,927.80711586516302,1e-9,'classical flight time (s)');
    assertRel(inDef.classical.hApoM,refKm.hApoM,1e-9, ...
        'the script own classical apogee against this file (m)');
    assertRel(inDef.classical.tofS ,refKm.tofS ,1e-9, ...
        'the script own classical flight time against this file (s)');
    assertAbs(inDef.gamStarR,gamStr,1e-12,'the script own gammaStar (rad)');

%% THE BURN IS ACTUALLY CUT SHORT. This is the whole mechanism: gammaStar needs
%% a burnout speed the booster OVERSHOOTS at full burn, so the second control
%% has to be exercised. A solve that quietly ran to propellant exhaustion would
%% satisfy the range assertion and nothing else here:
    assert(inDef.minEnergy.solved, ...
        'the shipped minimum-energy run did not report a solved trajectory');
    assert(inDef.minEnergy.cutFrac < 0.999, ...
        ['the minimum-energy solve cut the burn at %.6f of full, which is no ' ...
         'cut at all. gammaStar is reachable only below full burn on this ' ...
         'vehicle'],inDef.minEnergy.cutFrac);
    assert(inDef.minEnergy.propLeftKg > 100, ...
        'a cut-short burn must leave propellant unburned; this one left %.3f kg', ...
        inDef.minEnergy.propLeftKg);
    assertRel(inDef.minEnergy.cutFrac   ,0.96325683593750000,1e-4,'cutoff fraction');
    assertRel(inDef.minEnergy.tCutS     ,77.559280706465860 ,1e-4,'cutoff time (s)');
    assertRel(inDef.minEnergy.propLeftKg,1102.2949218749980 ,1e-4, ...
        'propellant left unburned (kg)');
    assertAbs(inDef.minEnergy.tCutS,inDef.minEnergy.cutFrac.*inDef.tBurn,1e-6, ...
        'the cutoff time against the fraction times the full burn (s)');
    assert(inDef.minEnergy.tCutS < inDef.tBurn, ...
        'the cutoff at %.4f s must precede the %.4f s full burn', ...
        inDef.minEnergy.tCutS,inDef.tBurn);

%% RESIDUAL ONE: THE BURNOUT GAMMA IS gammaStar, to the outer tolerance. This is
%% the residual that makes the trajectory minimum-energy rather than merely
%% on-target, and it is asserted against the tolerance the script was ASKED for
%% and not against whatever it happened to achieve:
            tolGamR = deg2rad(0.01);
    assertAbs(inDef.gamBurnout,gamStr,tolGamR, ...
        'flown burnout gamma against gammaStar (rad)');
    assertAbs(inDef.meGamResR,inDef.gamBurnout - gamStr,1e-15, ...
        'the reported gamma residual against the flown one (rad)');

%% RESIDUAL TWO: THE RANGE, to the range tolerance. Part 1 asserts it through
%% the miss; it is restated here as the residual the solve drives, because the
%% PAIR is what makes the trajectory minimum-energy and half a pair proves
%% nothing -- either full-burn arc satisfies the range residual alone:
    assert(abs(inDef.residM) <= tolM, ...
        'the minimum-energy range residual is %.2f m against a %.1f m tolerance', ...
        inDef.residM,tolM);
    assert(abs(inDef.depressed.gamBoR - gamStr) > 10.*tolGamR && ...
           abs(inDef.lofted.gamBoR    - gamStr) > 10.*tolGamR, ...
        ['the two full-burn arcs sit %.4f and %.4f deg from gammaStar; if ' ...
         'either were within the outer tolerance this part could not tell a ' ...
         'genuine minimum-energy solve from a branch selection'], ...
        rad2deg(abs(inDef.depressed.gamBoR - gamStr)), ...
        rad2deg(abs(inDef.lofted.gamBoR - gamStr)));

%% THE FLOWN ARC IS THE CLASSICAL ARC, to a tolerance taken FROM THE PHYSICS
%% rather than chosen for comfort. The closed form assumes an IMPULSIVE burn AT
%% the impact radius in a VACUUM; this flight burns for 77.6 s and reaches
%% gammaStar 80.6 km up at 4833.3 m/s rather than 4970.3 m/s at zero altitude.
%% MEASURED on the committed code: the apogee comes out 6.95 %% high and the
%% flight time 9.22 %% long. The budgets below are set a little above those --
%% 10 %% and 12 %% -- loose enough that a re-tuned integrator does not trip
%% them, and tight enough that the arc cannot wander to the 2118 km, 30.3 min
%% LOFTED arc (+208 %% and +96 %%) or to the 208 km, 11.0 min DEPRESSED one
%% (-70 %% and -29 %%), which are the two answers the old selector gave:
    assert(abs(inDef.meApoRelE) < 0.10, ...
        ['the flown apogee is %.3f km against the classical %.3f km, %.2f %% ' ...
         'out. A finite boost is worth a few per cent; this is not a few per ' ...
         'cent'],inDef.hApogee./1000,refKm.hApoM./1000,100.*inDef.meApoRelE);
    assert(abs(inDef.meTofRelE) < 0.12, ...
        ['the flown flight time is %.4f min against the classical %.4f min, ' ...
         '%.2f %% out'],inDef.tFlight./60,refKm.tofS./60,100.*inDef.meTofRelE);
    assertAbs(inDef.meApoRelE,(inDef.hApogee - refKm.hApoM)./refKm.hApoM,1e-12, ...
        'the reported apogee residual against this file');
    assertAbs(inDef.meTofRelE,(inDef.tFlight - refKm.tofS)./refKm.tofS,1e-12, ...
        'the reported flight-time residual against this file');

%% ...AND THE TOLERANCES ARE JUSTIFIED, not merely stated. Almost the whole gap
%% is the FINITE BURNOUT ALTITUDE: a Keplerian arc from the flown burnout state
%% apogees within 20 m of what the vehicle flew, so what the classical arc gets
%% wrong is where the boost ENDED and not how the solve converged. If a future
%% change ever makes the coast lossy enough for this to fail, the percentage
%% budgets above stop being attributable and must be re-derived rather than
%% quietly widened:
    assertAbs(inDef.minEnergy.hApoKepM,inDef.hApogee,1000, ...
        'Keplerian apogee of the burnout state against the flown apogee (m)');
    assert(abs(inDef.minEnergy.hApoKepM - refKm.hApoM) > ...
           0.9.*abs(inDef.hApogee - refKm.hApoM), ...
        ['the finite burnout altitude accounts for only %.3f km of the %.3f ' ...
         'km gap to the classical arc; the attribution the summary prints is ' ...
         'wrong'],abs(inDef.minEnergy.hApoKepM - refKm.hApoM)./1000, ...
        abs(inDef.hApogee - refKm.hApoM)./1000);

%% ...and the arc really is the textbook one rather than either full-burn arc.
%% The apogee-to-range ratio of a minimum-energy ballistic trajectory sits near
%% a quarter; the lofted arc here is at 0.67 and the depressed at 0.066:
             ratME = inDef.hApogee./inDef.rngAchM;
    assert(ratME > 0.20 && ratME < 0.28, ...
        ['the flown apogee-to-range ratio is %.4f; a minimum-energy arc sits ' ...
         'near 0.25 and the two full-burn arcs here sit at %.4f and %.4f'], ...
        ratME,inDef.depressed.hApoM./inDef.rngAchM, ...
        inDef.lofted.hApoM./inDef.rngAchM);

%% ---------------------------------------------------------------------
%% 7. THE SOLVE IS NOT A CONSTANT, and it is neither full-burn arc
%% ---------------------------------------------------------------------
%% One case cannot tell a working two-parameter solve from a hard-wired answer.
%% A second target, 1031 km further out, must move BOTH unknowns and must meet
%% ITS OWN classical reference -- a different gammaStar, a different apogee and
%% a different flight time:
            outFar2 = evalc(['[trFar2,inFar2] = run_ballistic_target(struct(' ...
                             '''lonTarget'',-40,''showPlots'',false));']);
    assert(contains(outFar2,'(nominal)') && ~contains(outFar2,'REFUSED'), ...
        'the 4206 km minimum-energy case did not solve. Summary was:\n%s',outFar2);
    assert(inFar2.rngReqM > inDef.rngReqM + 1e6, ...
        'the second minimum-energy target is not materially further out');

%% BOTH unknowns moved. A solve wired to a constant loft angle, or one that
%% never touched the cutoff, passes the shipped case and fails here:
    assert(abs(inFar2.loftDeg - inDef.loftDeg) > 1, ...
        ['the two minimum-energy targets solved to loft angles %.4f and %.4f ' ...
         'deg, less than a degree apart'],inFar2.loftDeg,inDef.loftDeg);
    assert(abs(inFar2.minEnergy.cutFrac - inDef.minEnergy.cutFrac) > 0.01, ...
        ['the two targets solved to cutoff fractions %.6f and %.6f; the ' ...
         'cutoff is the range control and it must move with the range'], ...
        inFar2.minEnergy.cutFrac,inDef.minEnergy.cutFrac);
    assert(abs(rad2deg(wrapPi(inFar2.psiLaunch - inDef.psiLaunch))) > 1, ...
        'the two targets returned the same azimuth to within a degree');

%% ...and the SECOND case meets its OWN reference, recomputed here rather than
%% reused. gammaStar FALLS as the range angle grows -- 37.870 deg at 3175 km
%% against 35.554 deg at 4206 km -- so a gammaStar wired to a constant passes
%% the shipped case and fails this one. That is precisely the mutation this
%% assertion exists to catch:
           psiFar2 = coorbital.util.greatCircle(latLn,lonLn,latTg,deg2rad(-40));
           gamStr2 = pi./4 - psiFar2./4;
           refFar2 = classicalArc(psiFar2,c);
    assert(abs(rad2deg(gamStr2 - gamStr)) > 2, ...
        ['the two gammaStar values are %.4f and %.4f deg, only %.4f deg ' ...
         'apart; this case cannot see a gammaStar wired to a constant'], ...
        rad2deg(gamStr2),rad2deg(gamStr),abs(rad2deg(gamStr2 - gamStr)));
    assertAbs(inFar2.gamBurnout,gamStr2,tolGamR, ...
        'the 4206 km burnout gamma against ITS OWN gammaStar (rad)');
    assert(abs(inFar2.meApoRelE) < 0.10, ...
        ['the 4206 km apogee is %.3f km against the classical %.3f km, ' ...
         '%.2f %% out'],inFar2.hApogee./1000,refFar2.hApoM./1000, ...
        100.*inFar2.meApoRelE);
    assert(abs(inFar2.meTofRelE) < 0.12, ...
        ['the 4206 km flight time is %.4f min against the classical %.4f ' ...
         'min, %.2f %% out'],inFar2.tFlight./60,refFar2.tofS./60, ...
        100.*inFar2.meTofRelE);
    assert(inFar2.missM <= tolM, ...
        'the 4206 km case missed by %.2f m',inFar2.missM);
            missIF = c.rE.*coorbital.util.greatCircle(trFar2.x(end,3), ...
                                                      trFar2.x(end,2), ...
                                                      latTg,deg2rad(-40));
    assert(missIF <= tolM, ...
        'the 4206 km case flown impact point is %.2f m from its target',missIF);
    assertBranchMiss(inFar2,tolM,'the 4206 km minimum-energy case');

%% ...AND THE MINIMUM-ENERGY ARC IS NEITHER OF THE TWO FULL-BURN ARCS, on both
%% cases. It lies BETWEEN them in apogee and in flight time, which is what a
%% trajectory carrying less burnout energy than either of them must do, and it
%% is far from both -- so a mode that quietly fell back to selecting a branch
%% could not pass this:
    assertBetweenArcs(inDef ,'the shipped 3175 km case');
    assertBetweenArcs(inFar2,'the 4206 km case');

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
%% AT BM/run_ballistic's 6 deg CLAMP THE SHIPPED BRACKET HAS NO INTERIOR
%% MAXIMUM, so range is monotone across all of it, one of the two branches is
%% empty everywhere on it, and the two-branch structure this script is built
%% around does not exist. That must be REFUSED with a named identifier, not
%% reported as a two-branch problem with one branch quietly missing:
             idClp = errorIdOf(@() run_ballistic_target( ...
                                    struct('alphaMax',6,'showPlots',false)));
    assert(strcmp(idClp,'coorbital:runBallisticTarget:maximumNotBracketed'), ...
        ['a 6 deg angle-of-attack clamp leaves the max-range loft angle ' ...
         'outside the shipped bracket and there is then only one branch on ' ...
         'it; expected coorbital:runBallisticTarget:maximumNotBracketed, got ' ...
         '"%s"'],idClp);

%% ...AND THE CAUSE IS BRACKET WIDTH, NOT THE CLAMP. This is a correction: this
%% file, the script's header, its refusal message and docs/README.md all used to
%% say that at 6 deg "the depressed branch does not exist" and that range is
%% monotone in loft everywhere. Measured, it is not. Widen loftMin to -140 at
%% the SAME 6 deg clamp and the hump is there, interior, with both branches
%% around it and no refusal -- it simply sits BELOW the shipped loftMin of
%% -40 deg, which is why the shipped bracket cannot find it:
             outWid = evalc(['[~,inWid] = run_ballistic_target(struct(' ...
                             '''alphaMax'',6,''loftMin'',-140,''nScanLoft'',21,' ...
                             '''showPlots'',false));']);
    assert(~inWid.refused && ~contains(outWid,'REFUSED'), ...
        ['the widened 6 deg bracket must solve rather than refuse. Summary ' ...
         'was:\n%s'],outWid);
    assert(inWid.loftStarDeg < -40, ...
        ['the 6 deg max-range loft angle is %.4f deg and must lie BELOW the ' ...
         'shipped loftMin of -40 deg, or the shipped refusal above is not ' ...
         'about bracket width at all'],inWid.loftStarDeg);
    assert(inWid.loftStarDeg > -43.5 && inWid.loftStarDeg < -42.3, ...
        ['the 6 deg max-range loft angle came out at %.4f deg against the ' ...
         '-42.907 deg measured on the committed code'],inWid.loftStarDeg);
    assertRel(inWid.rngMaxM,5211525.2732533719,1e-4,'6 deg maximum range (m)');

%% WHAT THE CLAMP ACTUALLY COSTS IS REACH, and that is the real reason 12 deg is
%% shipped. At 6 deg the depressed branch exists but bottoms out at 4708.463 km,
%% so it cannot serve the 3175 km shipped target at any loft angle -- widening
%% the bracket clears the bracketing refusal only to meet the unreachable-branch
%% one:
    assertRel(inWid.depressed.bandLo,4708462.651675,1e-4, ...
        '6 deg depressed-branch floor (m)');
    assert(~inWid.depressed.exists && inWid.depressed.bandLo > inDef.rngReqM, ...
        ['at 6 deg the depressed branch spans %.3f to %.3f km and must not ' ...
         'reach the %.3f km shipped target; THAT, and not the absence of a ' ...
         'hump, is what the clamp costs'],inWid.depressed.bandLo./1000, ...
        inWid.depressed.bandHi./1000,inDef.rngReqM./1000);

%% ...and the 12 deg the script ships buys that reach with MAXIMUM RANGE. The
%% trade is a number, so it is pinned as one:
    assert(inWid.rngMaxM > inDef.rngMaxM, ...
        ['raising the clamp from 6 to 12 deg must COST maximum range: ' ...
         '%.3f km at 12 deg against %.3f km at 6 deg'], ...
        inDef.rngMaxM./1000,inWid.rngMaxM./1000);
    assertAbs((inWid.rngMaxM - inDef.rngMaxM)./1000,156.2236,0.01, ...
        'maximum range given up to raise the clamp from 6 to 12 deg (km)');

%% ...and a ONE-BRANCH minimum-energy run must still SOLVE, and must not print
%% arithmetic it does not have. Two separate claims:
%%
%%   THE SOLVE SURVIVES A MISSING BRANCH. The outer loft bracket is normally the
%%   two full-burn branch solutions; here the depressed one does not exist, so
%%   the bracket takes the user's own loftMin at that end, which is admissible
%%   because the full burn OVERSHOOTS the target there. The minimum-energy
%%   trajectory is a property of the GEOMETRY, so the answer must be very nearly
%%   the same one the shipped 12 deg clamp finds -- a different COMMANDED loft
%%   angle, because the clamp changes what a command delivers, but the same
%%   burnout gamma and very nearly the same arc:
    assert(inWid.minEnergy.solved, ...
        ['the one-branch case did not solve the minimum-energy trajectory. ' ...
         'Summary was:\n%s'],outWid);
    assert(~inWid.depressed.exists, ...
        'this case is meant to have exactly one full-burn branch');
    assertAbs(inWid.gamBurnout,inDef.gamStarR,deg2rad(0.01), ...
        'one-branch burnout gamma against gammaStar (rad)');
    assert(abs(inWid.hApogee - inDef.hApogee) < 0.05.*inDef.hApogee, ...
        ['the 6 deg clamp found a %.3f km apogee against the 12 deg clamp''s ' ...
         '%.3f km. The minimum-energy arc is set by the GEOMETRY and the two ' ...
         'must very nearly agree'],inWid.hApogee./1000,inDef.hApogee./1000);
    assert(abs(inWid.loftDeg - inDef.loftDeg) > 10, ...
        ['the two clamps solved to loft angles %.4f and %.4f deg. They must ' ...
         'DIFFER: the clamp decides what a commanded attitude delivers, and ' ...
         'if they agreed the loft angle would not be doing any work'], ...
        inWid.loftDeg,inDef.loftDeg);

%%   AND IT PRINTS NO NaN. The contrast paragraph compares the two full-burn
%%   burnout energies and the missing branch supplies NaN; unguarded it read
%%   "leave burnout at NaN and -44.9073 MJ/kg, NaN %% apart", which in a results
%%   paragraph reads as a measurement:
    assert(~contains(outWid,'NaN'), ...
        ['the one-branch run printed a NaN, which in a results paragraph ' ...
         'reads as a measurement. Summary was:\n%s'],outWid);
    assert(contains(outWid,'ONLY ONE FULL-BURN ARC REACHES THIS TARGET'), ...
        ['a minimum-energy run with only one full-burn arc must say so ' ...
         'instead of printing half a comparison. Summary was:\n%s'],outWid);

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
%% 10. THE DISPLAY CHOICE: true scale by default, adaptive on request
%% ---------------------------------------------------------------------
%% TRUE SCALE IS THE SHIPPED DEFAULT, altExag = 1, on BOTH branches. The movie's
%% altitude inset is true-scale already and carries the profile the globe used
%% to exaggerate for, and coorbital.viz drops the "(altitude exaggerated Nx)"
%% caption clause at unity, so a true-scale picture makes no claim it is not
%% keeping:
    assert(inDef.altExag == 1 && inDep.altExag == 1, ...
        ['the shipped runs must draw at TRUE SCALE; they drew at %gx and %gx'], ...
        inDef.altExag,inDep.altExag);

%% ...and the ADAPTIVE rule is one word away, which is what to ask for when the
%% two branches are being compared: they differ by a factor of ten in apogee and
%% no fixed factor suits both -- at the 3x BM/run_ballistic hard-codes the
%% depressed arc is invisible. The sentinel is flown rather than asserted about,
%% because a sentinel nothing exercises is a sentinel that can rot:
             outAut = evalc(['[trAut,inAut] = run_ballistic_target(struct(' ...
                             '''altExag'',''auto'',''showPlots'',false));']);
    assert(contains(outAut,'(nominal)') && ~contains(outAut,'REFUSED'), ...
        'the auto-exaggeration run did not solve. Summary was:\n%s',outAut);
    assert(inAut.altExag == 2, ...
        ['altExag = ''auto'' drew the lofted arc at %gx against the 2x its ' ...
         '%.3f km apogee earns'],inAut.altExag,inAut.hApogee./1000);
    assert(inAut.altExag ~= inDef.altExag, ...
        'the auto run drew at the shipped true scale; the sentinel did nothing');
    assert(inAut.altExag == inDef.altExagRule(maxAlt(trAut,c.rE)), ...
        ['the auto run drew at %gx while its own rule gives %gx for the ' ...
         '%.3f km apogee it flew'],inAut.altExag, ...
        inDef.altExagRule(maxAlt(trAut,c.rE)),maxAlt(trAut,c.rE)./1000);

%% ...and the rule really does adapt, which is the whole reason it exists here:
%% the depressed arc's apogee earns 9x where the lofted arc's earns 2x. Asked of
%% the rule with the flown apogee rather than by flying a second auto run:
    assert(inDef.altExagRule(maxAlt(trDep,c.rE)) == 9, ...
        ['the rule gives %gx for the depressed arc''s %.3f km apogee against ' ...
         '2x for the lofted arc''s; one fixed factor cannot serve both'], ...
        inDef.altExagRule(maxAlt(trDep,c.rE)),maxAlt(trDep,c.rE)./1000);
    assert(inDef.altExagRule(40e3) == 30 && inDef.altExagRule(0) == 30, ...
        'the exaggeration rule must CAP at 30, including at a zero apogee');
    assert(inDef.altExagRule(2e6) == 2, ...
        'the exaggeration rule must FLOOR at 2; coorbital.viz refuses 0 or 1');

%% ...and a nonsense exaggeration must RAISE, before the solve rather than after
%% it, or a mistyped sentinel would cost a full two-branch run before it was
%% noticed:
             idExg = errorIdOf(@() run_ballistic_target( ...
                                    struct('altExag','atuo','showPlots',false)));
    assert(~isempty(idExg), ...
        'a mistyped altExag sentinel must raise, not be silently accepted');

%% The captions carry hemisphere letters. Both shipped coordinates are WEST of
%% the meridian, so a formatter that never reached its west branch would produce
%% a caption wrong by 200 degrees of longitude and looking perfectly reasonable:
    assert(strcmp(inDef.launchStr,'45.00N 100.00W'), ...
        'the launch caption reads "%s" against the expected "45.00N 100.00W"', ...
        inDef.launchStr);
    assert(strcmp(inDef.targetStr,'62.00N 60.00W'), ...
        'the target caption reads "%s" against the expected "62.00N 60.00W"', ...
        inDef.targetStr);

%% ---------------------------------------------------------------------
%% 11. THE BRANCH MEASUREMENT ITSELF, exercised where it CAN disagree
%% ---------------------------------------------------------------------
%% Part 5 re-derives the branch from the flown apogee and flight time, but every
%% run it does that on is a run where the measurement AGREES with the label. So
%% every assertion in it also passes if measureBranch is never called: replace
%% the call in the "WHICH BRANCH DID IT ACTUALLY FLY" block with
%% flownName = pickName and flownAgree = true, and the whole of the rest of this
%% file still passes while info.branchMeasured, info.branchAgrees,
%% info.branchTimeAgrees and the summary's "MEASURED as" line all become
%% restatements of the commanded branch. This part is what closes that, and the
%% headline claim of the script -- THE BRANCH IS MEASURED, NOT ASSUMED -- is
%% what it is closing.
%%
%% FINDING A CASE WHERE THE TWO CAN LEGITIMATELY DIFFER TAKES SOME CARE, and the
%% reason is worth recording. Apogee rises monotonically with the loft angle,
%% and each branch is solved on a bracket that has the max-range angle at one
%% end, so a solution STRICTLY INSIDE its bracket always measures as the branch
%% it was solved on: a mislocated maximum moves the reference apogee and the
%% bracket end together and cannot separate them. Widening or narrowing loftMin
%% and loftMax does not break that. What does is the SHARED endpoint: when the
%% required range falls within the range tolerance of the maximum,
%% coorbital.util.rangeSolve short-circuits at that endpoint and BOTH branches
%% come back holding the max-range arc itself. A lofted request is then served
%% by an arc whose apogee EQUALS rather than exceeds the reference, the strict
%% test names it depressed, and the label and the measurement part company.
%%
%% The case below asks for 5030.000 km on a 50 km tolerance against a
%% 5055.302 km maximum, and asks for the LOFTED arc by name:
             outMx = evalc(['[trMx,inMx] = run_ballistic_target(struct(' ...
                            '''latTarget'',62.5130919700,' ...
                            '''lonTarget'',-23.7426042098,''tolRangeKm'',50,' ...
                            '''branch'',''lofted'',''showPlots'',false));']);
    assert(~inMx.refused && ~isempty(trMx), ...
        'the near-maximum-range case must solve. Summary was:\n%s',outMx);

%% THE PREMISE FIRST: the solve really did land on the shared endpoint, and the
%% flown apogee really is the max-range arc's. If a future change moves either,
%% this part stops discriminating and must say so rather than pass vacuously:
    assert(inMx.atMaxRange, ...
        ['the near-maximum case did not land ON the max-range angle: it flew ' ...
         '%.10f deg against a %.10f deg maximum, so it cannot exercise the ' ...
         'measurement at all'],inMx.loftDeg,inMx.loftStarDeg);
    assertAbs(inMx.hApogee,inMx.hApoStarM,1e-9, ...
        'flown apogee against the max-range arc apogee (m)');
    assertAbs(inMx.tFlight,inMx.tFlyStarS,1e-9, ...
        'flown flight time against the max-range arc flight time (s)');

%% THE LABEL AND THE MEASUREMENT DISAGREE, and the measurement is not a copy of
%% the label. This is the assertion the whole part exists for:
    assert(strcmp(inMx.branchFlown,'lofted'), ...
        'the run was asked for the lofted branch and flew "%s"',inMx.branchFlown);
    assert(strcmp(inMx.branchMeasured,'depressed'), ...
        ['the flown arc IS the max-range arc, so its apogee is not GREATER ' ...
         'than the max-range arc''s and it must measure as "depressed"; it ' ...
         'measured as "%s". A measurement wired to the commanded branch says ' ...
         '"lofted" here, and this is the only assertion in this file that ' ...
         'sees the difference'],inMx.branchMeasured);
    assert(~inMx.branchAgrees, ...
        'info.branchAgrees must be FALSE when the measured branch is not the commanded one');
    assert(inMx.branchTimeAgrees, ...
        ['the apogee and the flight time must still agree with each other ' ...
         'here -- both are EQUAL to the max-range arc''s, so both read ' ...
         '"depressed" -- and the disagreement is with the LABEL alone']);

%% ...AND THE CAUTION MUST BE PRINTED, not merely flagged in the struct. The
%% summary is what a reader acts on, and it must diagnose the degeneracy rather
%% than blaming a bracketing step that did nothing wrong:
    assert(contains(outMx,'*** CAUTION ***'), ...
        'the disagreement printed no caution. Summary was:\n%s',outMx);
    assert(contains(outMx,'MEASURES as the depressed branch'), ...
        'the caution must name the MEASURED branch. Summary was:\n%s',outMx);
    assert(contains(outMx,'converged at their SHARED endpoint'), ...
        ['the caution must diagnose the max-range degeneracy, not blame the ' ...
         'max-range angle. Summary was:\n%s'],outMx);
    assert(contains(outMx,'MEASURED as      depressed'), ...
        ['the summary''s MEASURED line must carry the measured branch and ' ...
         'not the commanded one. Summary was:\n%s'],outMx);
    assertAbs(summaryNumber(outMx, ...
        'MEASURED as +depressed +\(from the FLOWN apogee ([-\d.]+) km'), ...
        inMx.hApogee./1000,1e-3,'printed measured-branch apogee (km)');

%% ---------------------------------------------------------------------
%% 12. A ROTATING EARTH IS REFUSED, not flown behind a caution
%% ---------------------------------------------------------------------
%% earthSpin true used to run to completion and print three statements that
%% contradicted each other: a limitations paragraph reporting rotation ON and
%% then that the initial bearing was the WHOLE answer, a cross-track WARNING
%% telling the user to set a bank angle that was already zero, and a 315.690 km
%% miss printed in the same format as a converged solution. The closed-form
%% azimuth is simply not the answer on a turning Earth, so the run is refused --
%% and refused BEFORE the bracketing, because 57 propagations spent on an answer
%% that is then discarded are 57 propagations wasted:
  [outSpn,infSpn,threwSpn,trSpn] = tryRun(struct('earthSpin',true,'showPlots',false));
    assert(~threwSpn,'earthSpin true must refuse, not throw:\n%s',outSpn);
    assert(contains(outSpn,'REFUSED'), ...
        'no refusal banner for earthSpin true:\n%s',outSpn);
    assert(contains(outSpn,'OUTER AZIMUTH ITERATION'), ...
        'the refusal must name the iteration the script does not have:\n%s',outSpn);
    assert(infSpn.refused && strcmp(infSpn.refusedWhy,'earthSpin'), ...
        'info must carry refused with refusedWhy = "earthSpin"; got "%s"', ...
        infSpn.refusedWhy);
    assert(isempty(trSpn), ...
        'a refused run must return an empty trajectory');
    assert(infSpn.omegaE > 0, ...
        'the refusal must report the non-zero omegaE it refused on');
    assertRel(infSpn.tgSpeedM,c.omegaE.*c.rE.*cos(latTg),1e-9, ...
        'eastward ground speed under the target (m/s)');
    assert(~isfield(infSpn,'rngMaxM'), ...
        ['the earthSpin refusal must come BEFORE the bracketing; it handed ' ...
         'back a maximum range, so the whole two-branch solve was paid for ' ...
         'and then thrown away']);

%% NONE OF THE THREE CONTRADICTIONS MAY SURVIVE. Each is asserted as the string
%% it was, because each was printed and each was read:
    assert(~contains(outSpn,'MISS DISTANCE'), ...
        ['the refusal must not report a miss against a target that stood ' ...
         'still:\n%s'],outSpn);
    assert(~contains(outSpn,'WHOLE answer'), ...
        ['the refusal must not claim the closed-form bearing is the whole ' ...
         'answer on a rotating Earth:\n%s'],outSpn);
    assert(~contains(outSpn,'WARNING'), ...
        ['the refusal must not warn about a bank angle on a trajectory it ' ...
         'never flew:\n%s'],outSpn);

%% ...and the surviving prose in a NORMAL run is gated on a still Earth: it
%% names the refusal rather than interpolating a rotation state that can no
%% longer be anything but OFF:
    assert(contains(outDef,'Earth rotation OFF'), ...
        'the shipped summary must state that Earth rotation is off:\n%s',outDef);
    assert(contains(outDef,'earthSpin'), ...
        ['the limitations paragraph must say WHY the non-rotating case is the ' ...
         'only one described:\n%s'],outDef);

%% ---------------------------------------------------------------------
%% 13. PLESETSK TO NEW YORK: the published case, and what it costs
%% ---------------------------------------------------------------------
%% 7132.320 km on a 64.0707 deg range angle, which is the geometry published
%% ballistic-missile range figures are usually quoted on. The classical
%% minimum-energy arc for it is 6581.8 m/s at gammaStar = 28.9823 deg, apogeeing
%% at 1205.989 km after 26.012 min, and this file computes all four from the
%% closed form rather than quoting them:
            latPl = deg2rad(62.925);
            lonPl = deg2rad(40.577);
            latNy = deg2rad(40.7128);
            lonNy = deg2rad(-74.006);
            psiPny = coorbital.util.greatCircle(latPl,lonPl,latNy,lonNy);
            refPny = classicalArc(psiPny,c);
            gamPny = pi./4 - psiPny./4;
    assertRel(c.rE.*psiPny  ,7132319.5118074585,1e-9,'Plesetsk to New York (m)');
    assertRel(rad2deg(gamPny),28.982320928298690,1e-9,'its gammaStar (deg)');
    assertRel(refPny.V      ,6581.8444425394830,1e-9,'its classical burnout speed (m/s)');
    assertRel(refPny.hApoM  ,1205989.0518133850,1e-9,'its classical apogee (m)');
    assertRel(refPny.tofS   ,1560.7229921211660,1e-9,'its classical flight time (s)');

%% THE SHIPPED BOOSTER CANNOT FLY IT, and that is a measurement rather than an
%% opinion: 7132.320 km against a 5055.302 km maximum, 2077.018 km short. It
%% must REFUSE -- through the envelope path, because a burn CUT SHORTER than the
%% full one carries LESS energy and cannot fly further than the full burn does.
%% A minimum-energy request beyond maximum range is beyond it by a WIDER margin
%% than a lofted or depressed one, and the refusal has to say so rather than
%% silently flying something else:
  [outPny,infPny,threwPny,trPny] = tryRun(struct('latLaunch',62.925, ...
                                                 'lonLaunch',40.577, ...
                                                 'latTarget',40.7128, ...
                                                 'lonTarget',-74.006, ...
                                                 'showPlots',false));
    assert(~threwPny,'the Plesetsk case must refuse, not throw:\n%s',outPny);
    assert(infPny.refused && strcmp(infPny.refusedWhy,'envelope'), ...
        'expected an envelope refusal; got refused = %d, why = "%s"', ...
        infPny.refused,infPny.refusedWhy);
    assert(isempty(trPny),'a refused run must return an empty trajectory');
    assertRel(infPny.rngReqM,c.rE.*psiPny,1e-9,'the refused required range (m)');
    assert(infPny.rngReqM > infPny.rngMaxM + 2.0e6, ...
        ['the Plesetsk case is only %.3f km beyond the maximum; it is meant to ' ...
         'be far beyond it'],(infPny.rngReqM - infPny.rngMaxM)./1000);
    assert(contains(outPny,'MINIMUM-ENERGY IS THE HARDER ASK'), ...
        ['the refusal must explain that minimum-energy needs a burnout speed ' ...
         'BELOW full burn and is therefore further out of reach, not nearer. ' ...
         'Summary was:\n%s'],outPny);

%% ...and the refusal must PRINT the classical arc it declined to fly, because
%% that is the number a user came for:
    assertAbs(summaryNumber(outPny,'V\* = ([-\d.]+) m/s'),refPny.V,0.1, ...
        'printed classical burnout speed (m/s)');
    assertAbs(summaryNumber(outPny,'gamma\* = ([-\d.]+) deg'),rad2deg(gamPny),1e-3, ...
        'printed gammaStar (deg)');
    assertAbs(summaryNumber(outPny,'flown: ([-\d.]+) km of apogee'), ...
        refPny.hApoM./1000,1e-2,'printed classical apogee (km)');
    assertAbs(summaryNumber(outPny,'apogee in ([-\d.]+) min'),refPny.tofS./60,1e-2, ...
        'printed classical flight time (min)');

%% GIVE IT A BOOSTER SIZED TO THE RANGE AND THE MODE REPRODUCES THE REFERENCE,
%% which is the acceptance test for the whole mode. The classical result
%% PRESUMES the burnout energy is free to choose -- you size the booster to the
%% range -- and the shipped placeholder is sized to about 5000 km, so the case
%% is flown here on a booster with the delta-V for 7132 km. Nothing else
%% changes: same vehicle, same pitch programme, same clamp, same solver
%% settings. MEASURED: 1257.452 km of apogee against 1205.989 km (+4.27 %%) and
%% 27.5096 min against 26.0120 min (+5.76 %%), with the burnout gamma
%% 0.0058 deg from gammaStar. The budgets are 8 %% and 8 %%, set from those
%% measurements the same way part 6 sets its own:
  [outBig,infBig,threwBig,trBig] = tryRun(struct('latLaunch',62.925, ...
                                                 'lonLaunch',40.577, ...
                                                 'latTarget',40.7128, ...
                                                 'lonTarget',-74.006, ...
                                                 'boosterFn',@bigBooster, ...
                                                 'showPlots',false));
    assert(~threwBig,'the sized-booster Plesetsk case threw:\n%s',outBig);
    assert(~infBig.refused && ~isempty(trBig), ...
        ['the sized-booster Plesetsk case must SOLVE; a booster with the ' ...
         'delta-V for the range is what the classical result presumes. ' ...
         'Summary was:\n%s'],outBig);
    assert(infBig.minEnergy.solved && infBig.minEnergy.cutFrac < 0.999, ...
        'the sized-booster case must still CUT THE BURN SHORT; it cut at %.6f', ...
        infBig.minEnergy.cutFrac);
    assertAbs(infBig.gamBurnout,gamPny,deg2rad(0.01), ...
        'sized-booster burnout gamma against gammaStar (rad)');
    assert(abs(infBig.meApoRelE) < 0.08, ...
        ['the sized-booster apogee is %.3f km against the classical %.3f km, ' ...
         '%.2f %% out; the reference is 1205.989 km'], ...
        infBig.hApogee./1000,refPny.hApoM./1000,100.*infBig.meApoRelE);
    assert(abs(infBig.meTofRelE) < 0.08, ...
        ['the sized-booster flight time is %.4f min against the classical ' ...
         '%.4f min, %.2f %% out; the reference is 26.012 min'], ...
        infBig.tFlight./60,refPny.tofS./60,100.*infBig.meTofRelE);
    assert(infBig.missM <= tolM, ...
        'the sized-booster Plesetsk case missed by %.2f m',infBig.missM);
            missIB = c.rE.*coorbital.util.greatCircle(trBig.x(end,3), ...
                                                      trBig.x(end,2),latNy,lonNy);
    assert(missIB <= tolM, ...
        'the sized-booster flown impact point is %.2f m from New York',missIB);

%% ---------------------------------------------------------------------
%% 14. THE MINIMUM-ENERGY GUARDS
%% ---------------------------------------------------------------------
%% A CUT-SHORT BURN JETTISONS PROPELLANT, so the coast mass is only well defined
%% if the WHOLE booster goes overboard. separation = false cannot express that
%% and must RAISE rather than fly against a mass the equations of motion would
%% disagree with. It is a user-block contradiction, not an unreachable target,
%% so it throws where the envelope cases refuse:
            idSep = errorIdOf(@() run_ballistic_target( ...
                                   struct('separation',false,'showPlots',false)));
    assert(~isempty(idSep), ...
        'minimum-energy with separation = false must raise, not fly');

%% ...and the two FULL-BURN modes are unaffected by that restriction, because
%% they burn to exhaustion and leave nothing unburned to account for. Flown
%% rather than asserted about, or the exemption is only a claim:
           outNoSep = evalc(['[~,inNoSep] = run_ballistic_target(struct(' ...
                             '''separation'',false,''branch'',''lofted'',' ...
                             '''showPlots'',false));']);
    assert(~inNoSep.refused && ~contains(outNoSep,'REFUSED'), ...
        ['a LOFTED run with separation = false must still fly. Summary ' ...
         'was:\n%s'],outNoSep);

%% A cutFracMin at or above the full burn leaves the inner solve nothing to
%% bisect, and must be caught in the user-block checks rather than inside the
%% solver:
            idCut = errorIdOf(@() run_ballistic_target( ...
                                   struct('cutFracMin',1,'showPlots',false)));
    assert(~isempty(idCut),'cutFracMin = 1 must raise');
            idTol = errorIdOf(@() run_ballistic_target( ...
                                   struct('tolGamDeg',0,'showPlots',false)));
    assert(~isempty(idTol),'tolGamDeg = 0 must raise');

%% THE OUTER TOLERANCE IS ACTUALLY HONOURED, which nothing above can see: every
%% assertion so far is written against the SHIPPED 0.01 deg. Tighten it by a
%% factor of ten and the residual must fall below the tighter figure, which a
%% solve that ignored the tolerance -- or one that stopped after a fixed number
%% of steps -- could not do:
          outTight = evalc(['[~,inTight] = run_ballistic_target(struct(' ...
                            '''tolGamDeg'',0.001,''showPlots'',false));']);
    assert(~inTight.refused,'the tightened run did not solve:\n%s',outTight);
    assert(abs(rad2deg(inTight.meGamResR)) < 0.001, ...
        ['asked for 0.001 deg on the burnout gamma the solve delivered ' ...
         '%.6f deg'],abs(rad2deg(inTight.meGamResR)));
    assert(inTight.minEnergy.outerIters > inDef.minEnergy.outerIters, ...
        ['a ten times tighter outer tolerance took %d outer steps against ' ...
         'the shipped %d; it cannot cost the same'], ...
        inTight.minEnergy.outerIters,inDef.minEnergy.outerIters);
    assert(abs(rad2deg(inTight.meGamResR)) < abs(rad2deg(inDef.meGamResR)), ...
        'the tightened run is no closer to gammaStar than the shipped one');

%% Nothing above asked for a figure, so nothing above may have made one:
            nFig1  = numel(findall(groot,'Type','figure'));
    assert(nFig1 == nFig0, ...
        ['run_ballistic_target left %d new figure(s) open with showPlots ' ...
         'false; the suite must run headless'],nFig1 - nFig0);
end

function ref = classicalArc(psiR,c)
%% Purpose:
%
%  The CLASSICAL minimum-energy ballistic trajectory for a free-flight range
%  angle, in closed form, WRITTEN OUT HERE rather than called from
%  BM/run_ballistic_target. That duplication is the point: this is the reference
%  the minimum-energy mode is graded against, and a reference read out of the
%  code under test cannot see that code change the formula. Mutate gammaStar in
%  the script to a flat 45 deg and this function still returns the right answer,
%  so the assertions that compare the two fire.
%
%      V*^2      = (mu/rE) * 2 sin(psi/2) / (1 + sin(psi/2))
%      gammaStar = 45 deg - psi/4
%
%  and the flight time is twice the burnout-to-apogee time of the two-body arc
%  those two determine, by Kepler's equation.
%
%% Inputs:
%
%  psiR             [1 x 1]                     Free-flight range angle (rad)
%
%  c                Struct                      coorbital.util.missileConst
%
%% Outputs:
%
%  ref              Struct                      V (m/s), gamR (rad), hApoM (m)
%                                               apogee ALTITUDE, tofS (s)
%
%% References:
%   [1] Bate, R.R., Mueller, D.D., White, J.E., "Fundamentals of
%       Astrodynamics," Dover, 1971, Ch. 6.
%
%% Revision History:
%  Michael Casey                                                08/08/2026
%  Copyright 2026 Coorbital, Inc.
%% ------------------------ Begin Code Sequence ---------------------------

             sHalf = sin(psiR./2);
                V2 = (c.muE./c.rE).*2.*sHalf./(1 + sHalf);
             ref.V = sqrt(V2);
          ref.gamR = pi./4 - psiR./4;
               eps = V2./2 - c.muE./c.rE;
                aM = -c.muE./(2.*eps);
                hM = c.rE.*ref.V.*cos(ref.gamR);
                pM = hM.^2./c.muE;
               ecc = sqrt(max(0,1 - pM./aM));
         ref.hApoM = aM.*(1 + ecc) - c.rE;
                nu = acos((pM./c.rE - 1)./ecc);
                EA = 2.*atan(sqrt((1 - ecc)./(1 + ecc)).*tan(nu./2));
                MA = EA - ecc.*sin(EA);
             nMean = sqrt(c.muE./aM.^3);
          ref.tofS = 2.*(pi - MA)./nMean;
end

function assertBetweenArcs(info,what)
%% Purpose:
%
%  Assert that the flown minimum-energy trajectory is NEITHER of the two
%  full-burn arcs, and lies strictly between them in both apogee and flight
%  time. Exists because every range and miss assertion in this file is satisfied
%  by either full-burn arc: a minimum-energy mode that quietly fell back to
%  selecting a branch -- which is exactly what it used to do -- would pass all
%  of them. This is what sees the difference.
%
%  BOTH SEPARATIONS ARE CHECKED AS FACTORS rather than as differences, so the
%  thresholds do not have to be re-derived for every geometry: the trajectory
%  must be at least half again the depressed apogee and at most two thirds of
%  the lofted one, and its flight time strictly inside the pair.
%
%% Inputs:
%
%  info             Struct                      A non-refused run's info, with
%                                               depressed, lofted, hApogee and
%                                               tFlight
%
%  what             Char [1 x n]                Name of the case, for messages
%
%% Outputs:
%
%  none                                         Throws on any failed assertion
%
%% Revision History:
%  Michael Casey                                                08/08/2026
%  Copyright 2026 Coorbital, Inc.
%% ------------------------ Begin Code Sequence ---------------------------

    assert(info.depressed.exists && info.lofted.exists, ...
        '%s needs both full-burn arcs for this comparison',what);
    assert(info.hApogee > 1.5.*info.depressed.hApoM, ...
        ['%s apogees at %.3f km against the DEPRESSED arc''s %.3f km; the ' ...
         'minimum-energy trajectory must be clearly above it'], ...
        what,info.hApogee./1000,info.depressed.hApoM./1000);
    assert(info.hApogee < 0.67.*info.lofted.hApoM, ...
        ['%s apogees at %.3f km against the LOFTED arc''s %.3f km; the ' ...
         'minimum-energy trajectory must be clearly below it'], ...
        what,info.hApogee./1000,info.lofted.hApoM./1000);
    assert(info.tFlight > info.depressed.tFlyS && ...
           info.tFlight < info.lofted.tFlyS, ...
        ['%s flies for %.2f s against %.2f s depressed and %.2f s lofted; it ' ...
         'must lie strictly between them'],what,info.tFlight, ...
        info.depressed.tFlyS,info.lofted.tFlyS);
end

function bst = bigBooster()
%% Purpose:
%
%  A booster with the delta-V for the Plesetsk-to-New-York range, for part 13.
%  The classical minimum-energy result PRESUMES the burnout energy is free to
%  choose -- you size the booster to the range -- and the shipped placeholder is
%  sized to about 5000 km, so the 7132 km case cannot be flown on it at all.
%  This is the shipped booster with the propellant load, the dry structure and
%  the thrust scaled together, which leaves the burn time and the liftoff
%  thrust-to-weight essentially where they were and buys about 860 m/s of ideal
%  delta-V. PLACEHOLDER values, like everything else in this library.
%
%% Inputs:
%
%  none
%
%% Outputs:
%
%  bst              Struct                      Booster parameters; see
%                                               coorbital.util.boosterDefaults
%
%% Revision History:
%  Michael Casey                                                08/08/2026
%  Copyright 2026 Coorbital, Inc.
%% ------------------------ Begin Code Sequence ---------------------------

               bst = coorbital.util.boosterDefaults();
       bst.massDry = 2000;              %kg,  PLACEHOLDER, scaled with the load
      bst.massProp = 52000;             %kg,  PLACEHOLDER
     bst.thrustVac = 1615000;           %N,   PLACEHOLDER, holds T/W near 3
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

function assertBranchMiss(info,tolM,what)
%% Purpose:
%
%  Hold every branch that CLAIMS to reach the target to the range tolerance,
%  using the branch's own reported miss. Exists because the two-branch table is
%  the deliverable of run_ballistic_target and its miss row was unpinned: a
%  transposed coorbital.util.greatCircle inside solveBranch moves the printed
%  number and the returned one together, so comparing the printed pair against
%  info cannot see it. The table then reads 1354.24 and 235.58 m against the
%  true 666.05 and 115.86, still says "reaches the target: yes", and nothing
%  objects under a stated 1000 m tolerance. This objects.
%
%  A branch that does NOT reach the target is skipped rather than failed: its
%  miss is deliberately NaN, and NaN <= tolM is false, so an unguarded
%  assertion would fire on every legitimate one-branch geometry in this file.
%
%% Inputs:
%
%  info             Struct                      run_ballistic_target's info,
%                                               from a run that was not refused
%
%  tolM             [1 x 1]                     Range tolerance the run was
%                                               given (m)
%
%  what             Char [1 x n]                Name of the run, for the message
%
%% Outputs:
%
%  none                                         Throws on failure
%
%% Revision History:
%  Michael Casey                                                08/07/2026
%  Copyright 2026 Coorbital, Inc.
%% ------------------------ Begin Code Sequence ---------------------------

               brs = {info.depressed,info.lofted};
    for kb = 1:numel(brs)
                br = brs{kb};
        if br.exists
            assert(isfinite(br.missM) && br.missM <= tolM, ...
                ['on %s the %s arc reports "reaches the target" and then ' ...
                 'misses it by %.2f m against a %.1f m tolerance'], ...
                what,br.name,br.missM,tolM);
        end
    end
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
