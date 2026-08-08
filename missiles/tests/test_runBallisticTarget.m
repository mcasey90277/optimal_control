function test_runBallisticTarget()
%% Purpose:
%
%  Pin BM/run_ballistic_target, the ballistic point-to-point targeting script:
%  the closed-form launch azimuth, the CERTIFIED bracketing of the max-range
%  loft angle, the TWO full-burn branch solves on the certified one-sided
%  intervals either side of it, the MINIMUM-ENERGY CONSTRAINED MINIMISATION on
%  loft angle and cutoff fraction, the after-the-fact measurement of which
%  branch was actually flown, and the refusals at both ends of the envelope, for
%  an unreachable branch and for a minimum-energy request that cannot be solved.
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
%    flown loft angle against the certified maximiser interval, and asserts the
%    two agree. A selector wired to a constant passes every range assertion and
%    fails this one.
%
%    THE TWO BRANCHES ARE ASSERTED TO BE DIFFERENT. Measured on the shipped
%    geometry at the vehicle's 6 deg clamp: 2.55 times the apogee and 1.48
%    times the flight time separate them. A script that solved the same branch
%    twice and labelled the copies differently would satisfy every miss
%    assertion in this file and fail that one. The thresholds below are set from
%    those measurements; they were 5x and 2x when this script shipped a 12 deg
%    clamp, and the clamp -- not the assertion -- is what changed.
%
%    THE MINIMUM-ENERGY MODE IS GRADED ON ITS OWN OBJECTIVE, in parts 6, 6b, 7
%    and 13, and the objective is a MINIMUM of a stated quantity rather than a
%    residual against a formula:
%
%      part 6  asserts that the flown burnout energy is BELOW the energy at the
%              two neighbouring FEASIBLE points -- same range, different loft --
%              by a margin far larger than the measured noise of the inner
%              feasibility solve, and below both full-burn arcs' energies;
%
%      part 6b asserts the VACUUM EQUAL-RADIUS LIMIT: the same constrained
%              minimisation, applied to an impulsive burn at r = rE in a vacuum,
%              must return the classical V* and gammaStar = 45 deg - psi/4.
%              That is the check that says the objective is the right one. It is
%              written out in this file's own hand, in closed form, so a mutated
%              formula in the script cannot be mirrored by it.
%
%    The classical arc is still computed here and still compared against, but as
%    a DIAGNOSTIC: it is derived for two endpoints at the same radius and this
%    burnout is 82 km up, so agreement with it is evidence about the size of
%    that modelling gap and not evidence that the solve converged.
%
%    ...AND THE ANSWER IS SHOWN NOT TO BE EITHER FULL-BURN ARC. Every range and
%    miss assertion here is satisfied by the lofted arc and by the depressed
%    one, so a minimum-energy mode that quietly fell back to SELECTING a branch
%    -- which is what it once did -- would pass all of them. assertBetweenArcs
%    and the energy assertions are what see the difference.
%
%  WHY THE FLOWN CASE IS NEITHER EQUATORIAL NOR DUE EAST. From (0,0) the
%  central angle reduces to acos(cos(lat2) cos(lon2)), which is symmetric in
%  its two arguments, so a due-east equatorial case is provably BLIND to a
%  latitude/longitude transposition at every great-circle call site. The
%  shipped configuration is 45 N 100 W to 62 N 28 W on a 40.6 deg azimuth,
%  where all four coordinates are distinct; part 2 below PROVES that geometry
%  is discriminating rather than assuming it, in metres against the miss
%  budget the assertions are written in.
%
%    THE MEASUREMENT IS EXERCISED WHERE IT CAN FAIL, in part 11. Every run in
%    part 5 is one where the measurement AGREES with the label, so all of them
%    also pass if branchOfLoft is never called and the flown name is copied from
%    the commanded one. Part 11 flies the one geometry on which the two can
%    legitimately part company -- a target within the range tolerance of the
%    largest range the search can certify, where both branch solves converge
%    INSIDE the unresolved maximiser interval -- and asserts the disagreement,
%    the COALESCED classification, the caution and the printed line.
%
%  THE SHIPPED TARGET IS INSIDE THE 6 deg DEPRESSED BAND, which is why all three
%  modes fly on the shipped configuration. The clamp is now a VEHICLE limit,
%  BM/vehicle_bm's alphaMaxDeg, read by this script and by BM/run_ballistic
%  alike; part 9 pins that, pins the bracket-width refusal at the old -40 deg
%  loftMin, and pins what raising the clamp to 12 deg costs -- 156.224 km of
%  maximum range -- so the trade is a number rather than an argument.
%
%  COST. A trajectory propagation is about 0.05 s. A 'lofted' or 'depressed' run
%  takes 62 of them, about 3 s; a 'minimum-energy' run takes about 730, roughly
%  28 s, because the constrained minimisation is about 670 of those on its own
%  -- every outer evaluation is a whole inner feasibility solve. The SHIPPED
%  configuration is flown at its shipped settings throughout, which is the
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
            lonTgD = -28;
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
    assertRel(inDef.loftStarDeg,-42.907364725882300,1e-4,'max-range loft angle (deg)');
    assertRel(inDef.rngMaxM    ,5211525.2732604900 ,1e-4,'maximum range (m)');
    assertRel(inDef.rngMinM    ,556603.26462910500 ,1e-4,'envelope floor (m)');
    assertRel(inDef.hApoStarM  ,974102.46048605197 ,1e-4,'max-range arc apogee (m)');
    assertRel(inDef.tFlyStarS  ,1314.7768233764600 ,1e-4,'max-range arc flight time (s)');
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

%% THE MAXIMISER INTERVAL IS CERTIFIED, AND IT IS AN INTERVAL. Golden section
%% proves only that the maximiser lies inside it, and the two branch brackets
%% are built from its ENDS rather than from its midpoint -- which is the whole
%% of finding 3. So the interval must be returned, must be no wider than the
%% tolerance it was asked for, and must contain the reported max-range angle:
    assert(isfield(inDef.maxRange,'aL') && isfield(inDef.maxRange,'bL'), ...
        ['info.maxRange must carry the CERTIFIED maximiser interval aL and bL. ' ...
         'A midpoint alone cannot split the loft axis into two monotone ' ...
         'branches, because one of the two halves then straddles the maximum']);
    assertRel(rad2deg(inDef.maxRange.aL),-42.928957736544199,1e-4, ...
        'certified maximiser interval, lower end (deg)');
    assertRel(rad2deg(inDef.maxRange.bL),-42.885771715220397,1e-4, ...
        'certified maximiser interval, upper end (deg)');
    assert(rad2deg(inDef.maxRange.bL - inDef.maxRange.aL) <= 0.05 + 1e-12, ...
        ['the certified maximiser interval is %.6f deg wide against the ' ...
         '0.05 deg tolLoftDeg it was asked for'], ...
        rad2deg(inDef.maxRange.bL - inDef.maxRange.aL));
    assert(inDef.maxRange.aL <= deg2rad(inDef.loftStarDeg) && ...
           deg2rad(inDef.loftStarDeg) <= inDef.maxRange.bL, ...
        'the reported max-range angle is not inside the interval it came from');
    assertRel(inDef.maxRange.fAL,5211525.2583527200,1e-4, ...
        'range at the lower end of the maximiser interval (m)');
    assertRel(inDef.maxRange.fBL,5211525.2313836403,1e-4, ...
        'range at the upper end of the maximiser interval (m)');
    assertAbs(inDef.rngTopM,min(inDef.maxRange.fAL,inDef.maxRange.fBL),1e-9, ...
        'the top of the certifiable band against the smaller end range (m)');
    assert(~inDef.coalescedRq, ...
        ['the shipped target must sit clear of the unresolved band around the ' ...
         'maximum, or the branch labels carry no information']);

%% THE HEADLINE NUMBERS of the flown arc, to 1e-4 relative, read at full
%% printed precision from a run of the committed code:
    assertRel(inDef.rngReqM  ,4828045.2122497400 ,1e-4,'required range (m)');
    assertRel(inDef.rngAchM  ,4828006.2029809402 ,1e-4,'achieved range (m)');
    assertRel(inDef.missM    ,39.009268797573501 ,1e-4,'miss distance (m)');
    assertRel(inDef.residM   ,-39.00926879700270 ,1e-4,'range residual (m)');
    assertRel(inDef.loftDeg  ,-31.765561810950701,1e-4,'flown loft angle (deg)');
    assertRel(inDef.psiLaunch,0.70782521841714603,1e-4,'launch azimuth (rad)');
    assertRel(inDef.tFlight  ,1276.9020361444400 ,1e-4,'flight time (s)');
    assertRel(inDef.hApogee  ,965613.17589598403 ,1e-4,'apogee (m)');
    assertRel(inDef.tApogee  ,664.48120097246898 ,1e-4,'time of apogee (s)');
    assertRel(inDef.vImpact  ,2514.9710204880501 ,1e-4,'impact speed (m/s)');
    assertRel(rad2deg(inDef.gamImpact),-33.813601126209998,1e-4,'impact angle (deg)');
    assertRel(inDef.tBurnout ,80.113240197794095 ,1e-4,'cutoff time (s)');
    assertRel(inDef.hBurnout ,82181.682991809202 ,1e-4,'burnout altitude (m)');
    assertRel(inDef.vBurnout ,5682.3474321055897 ,1e-4,'burnout speed (m/s)');
    assertRel(rad2deg(inDef.gamBurnout),33.329037982044098,1e-4,'burnout gamma (deg)');
    assertRel(inDef.qBstMax  ,93438.812519846600 ,1e-4,'boost peak dynamic pressure (Pa)');
    assertRel(inDef.nBstMax  ,37.975521841799702 ,1e-4,'peak boost sensed load (g)');
    assertRel(inDef.qReMax   ,4539046.0693888399 ,1e-4,'re-entry peak dynamic pressure (Pa)');
    assertRel(inDef.nReMax   ,49.509550727719901 ,1e-4,'peak re-entry deceleration (g)');

%% THE CLAMP THE RUN FLEW IS THE VEHICLE'S, not a literal in the script. Both BM
%% entry scripts read BM/vehicle_bm's alphaMaxDeg, and this file restates the
%% 6 deg in its own hand rather than reading it back out of the vehicle file:
    assertAbs(inDef.alphaMaxDeg,6,1e-12, ...
        ['the run must fly the VEHICLE''S angle-of-attack limit (deg). It is ' ...
         'not a targeting degree of freedom, and the two BM scripts disagreeing ' ...
         'about it is what this pins']);
                vBM = vehicle_bm();
    assertAbs(vBM.alphaMaxDeg,inDef.alphaMaxDeg,1e-12, ...
        'the flown clamp against BM/vehicle_bm''s own alphaMaxDeg (deg)');

%% The impact point itself, pinned separately from the range. A single scalar
%% range is built from both coordinates and cannot tell them apart:
    assertRel(rad2deg(inDef.latImpact), 62.000070946438400,1e-4,'impact latitude (deg)');
    assertRel(rad2deg(inDef.lonImpact),-28.000730970067399,1e-4,'impact longitude (deg)');

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
        ['the minimum-energy default must show its solve: the objective, the ' ...
         'constraint, and what was minimised over what']);
    assert(contains(outDef,'A CONSTRAINED MINIMISATION, stated and then solved'), ...
        ['the summary must STATE the optimisation problem. It used to call ' ...
         'itself a two-residual solve while driving the burnout gamma to a ' ...
         'value derived for an equal-radius vacuum arc, which is not a ' ...
         'condition this trajectory satisfies']);
    assert(contains(outDef,'IT IS A MINIMUM, and here is the evidence'), ...
        ['a minimisation must show that it found a minimum, by reporting the ' ...
         'objective at the neighbouring feasible points']);
    assert(contains(outDef,'THE NOISE IS MEASURED, NOT ASSUMED SMALL'), ...
        ['the propagated effect of the inner feasibility tolerance on the outer ' ...
         'objective must be MEASURED and printed, not asserted to be small']);
    assert(contains(outDef,'AGAINST THE CLASSICAL VACUUM EQUAL-RADIUS ARC'), ...
        ['the minimum-energy default must be reported against the closed-form ' ...
         'reference for its range angle, and must name the limit that ' ...
         'reference is derived in']);
    assert(contains(outDef,'DIAGNOSTIC, not a residual'), ...
        ['the summary must say that the classical comparison is a diagnostic. ' ...
         'Reporting it as a residual is the critical finding this rewrite ' ...
         'answers']);
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

%% ...and the measurement is REDONE HERE, from the flown loft angles against the
%% CERTIFIED maximiser interval, so that a branchOfLoft wired to a constant is
%% caught rather than believed. Position against that interval is the only
%% branch invariant available: apogee and flight time are neither monotone in
%% the loft angle for a powered atmospheric arc nor separated near the maximum,
%% and they classified the branch until 08/08/2026:
    assert(deg2rad(inLof.loftDeg) > inDef.maxRange.bL, ...
        ['the "lofted" arc solved to %.6f deg, which is not above the ' ...
         '%.6f deg upper end of the certified maximiser interval'], ...
        inLof.loftDeg,rad2deg(inDef.maxRange.bL));
    assert(deg2rad(inDep.loftDeg) < inDef.maxRange.aL, ...
        ['the "depressed" arc solved to %.6f deg, which is not below the ' ...
         '%.6f deg lower end of the certified maximiser interval'], ...
        inDep.loftDeg,rad2deg(inDef.maxRange.aL));

%% The apogee and flight time still ORDER themselves the way a reader expects,
%% and that is asserted -- as a description of this geometry, which is all it
%% is. It is deliberately NOT the branch test any more:
             hApoL  = maxAlt(trLof,c.rE);
             hApoD  = maxAlt(trDep,c.rE);
    assert(hApoL > inDef.hApoStarM && trLof.t(end) > inDef.tFlyStarS, ...
        ['the "lofted" arc apogees at %.3f km in %.1f s against the max-range ' ...
         'arc''s %.3f km in %.1f s; on THIS geometry it exceeds both'], ...
        hApoL./1000,trLof.t(end),inDef.hApoStarM./1000,inDef.tFlyStarS);
    assert(hApoD < inDef.hApoStarM && trDep.t(end) < inDef.tFlyStarS, ...
        ['the "depressed" arc apogees at %.3f km in %.1f s against the ' ...
         'max-range arc''s %.3f km in %.1f s; on THIS geometry it falls short ' ...
         'of both'],hApoD./1000,trDep.t(end),inDef.hApoStarM./1000,inDef.tFlyStarS);
    assert(strcmp(inLof.branchByApogee,'lofted') && ...
           strcmp(inDep.branchByApogee,'depressed'), ...
        'the descriptive apogee reading disagrees with the certified branch here');

%% THE TWO BRANCHES MUST BE GENUINELY DIFFERENT TRAJECTORIES. A script that
%% solved one branch twice and labelled the copies differently would satisfy
%% every miss assertion in this file. MEASURED on the shipped geometry at the
%% vehicle's 6 deg clamp: 1453.661 km of apogee against 569.177 km, a factor of
%% 2.55, and 1550.500 s against 1047.949 s, a factor of 1.48. The thresholds are
%% set below those measurements -- and they are LOWER than the 5x and 2x this
%% file used to carry, because the shipped clamp fell from 12 deg to the
%% vehicle's 6 deg and a tighter clamp compresses the whole family:
    assert(hApoL > 2.0.*hApoD, ...
        ['the lofted and depressed arcs apogee at %.3f km and %.3f km, a ' ...
         'factor of only %.2f. They are not two different branches'], ...
        hApoL./1000,hApoD./1000,hApoL./hApoD);
    assert(trLof.t(end) > 1.3.*trDep.t(end), ...
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
    assertRel(inDep.loftDeg  ,-126.72856844054301,1e-4,'depressed loft angle (deg)');
    assertRel(inDep.rngAchM  ,4828502.4821094898 ,1e-4,'depressed achieved range (m)');
    assertRel(inDep.missM    ,457.26985976312702 ,1e-4,'depressed miss (m)');
    assertRel(inDep.hApogee  ,569176.71056946299 ,1e-4,'depressed apogee (m)');
    assertRel(inDep.tFlight  ,1047.9490771616999 ,1e-4,'depressed flight time (s)');
    assertRel(inDep.vImpact  ,1781.3245968933699 ,1e-4,'depressed impact speed (m/s)');
    assertRel(rad2deg(inDep.gamImpact),-22.967933808687101,1e-4, ...
        'depressed impact angle (deg)');
    assertRel(rad2deg(inDep.latImpact), 61.999168072388397,1e-4, ...
        'depressed impact latitude (deg)');
    assertRel(rad2deg(inDep.lonImpact),-27.991431762343499,1e-4, ...
        'depressed impact longitude (deg)');

%% THE DEPRESSED BRANCH AT A 6 deg CLAMP IS A DEGENERATE TAIL, and the shipped
%% loftMin of -140 deg is what it takes to hold it. The solved commanded
%% attitude is -126.729 deg, which is not a pitch programme anyone would fly; it
%% is where the branch goes once the clamp saturates the achievable pitch-over,
%% and it is pinned so that nobody quietly narrows loftMin and loses the branch:
    assert(inDep.loftDeg < -100, ...
        ['the depressed arc solved to %.4f deg. At the vehicle''s 6 deg clamp ' ...
         'this branch lives far below the horizon, and a loftMin above it makes ' ...
         'the branch unreachable rather than making the vehicle better behaved'], ...
        inDep.loftDeg);

%% THE STEEPER ARC IS THE FASTER ONE HERE, which is the opposite of the vacuum
%% intuition and is a DRAG result: the shallow depressed descent spends far
%% longer in dense air. Pinned because the summary states it in words, and a
%% statement in words that no assertion checks is a statement that can go stale:
    assert(inLof.vImpact > 1.5.*inDep.vImpact, ...
        ['the lofted arc arrives at %.1f m/s and the depressed at %.1f m/s; ' ...
         'the summary claims the steep arc arrives much the faster of the two'], ...
        inLof.vImpact,inDep.vImpact);
    assert(inLof.gamImpact < inDep.gamImpact, ...
        'the lofted arc must arrive STEEPER: %.3f deg against %.3f deg', ...
        rad2deg(inLof.gamImpact),rad2deg(inDep.gamImpact));

%% ---------------------------------------------------------------------
%% 6. THE MINIMUM-ENERGY SOLVE, graded on ITS OWN OBJECTIVE
%% ---------------------------------------------------------------------
%% THIS IS THE PART THE MODE EXISTS FOR, and it is not the part it used to be.
%% The mode solves an explicit constrained minimisation -- least burnout
%% specific energy subject to making the required range -- and until 08/08/2026
%% it instead drove the burnout flight-path angle to gammaStar = 45 deg - psi/4
%% with psi the PAD-TO-TARGET central angle. That closed form is derived for a
%% free-flight arc with both endpoints on the SAME radius; this burnout is 82 km
%% up and downrange of the pad, so the residual did not apply to the arc being
%% flown and the reported agreement verified the wrong condition. gammaStar is
%% still computed here and still compared against, as a DIAGNOSTIC.
%%
%% THE REFERENCE IS COMPUTED HERE, IN THIS FILE'S OWN HAND, from the range angle
%% and coorbital.util.missileConst. It is NOT read out of info.classical: a test
%% that took its reference from the code under test could not see the code
%% change the formula:
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
    assertRel(rad2deg(gamStr),34.157232983539901,1e-9,'gammaStar (deg)');
    assertRel(refKm.V        ,5807.2214117471003,1e-9,'classical burnout speed (m/s)');
    assertRel(refKm.hApoM    ,952696.11896354402,1e-9,'classical apogee (m)');
    assertRel(refKm.tofS     ,1210.1170457957700,1e-9,'classical flight time (s)');
    assertRel(inDef.classical.hApoM,refKm.hApoM,1e-9, ...
        'the script own classical apogee against this file (m)');
    assertRel(inDef.classical.tofS ,refKm.tofS ,1e-9, ...
        'the script own classical flight time against this file (s)');
    assertAbs(inDef.gamStarR,gamStr,1e-12,'the script own gammaStar (rad)');

%% THE BURN IS ACTUALLY CUT SHORT. This is the mechanism: a full burn has no
%% energy freedom left, so a solve that quietly ran to propellant exhaustion
%% would satisfy the range assertion and would be minimising nothing:
    assert(inDef.minEnergy.solved, ...
        'the shipped minimum-energy run did not report a solved trajectory');
    assert(inDef.minEnergy.cutFrac < 0.999, ...
        ['the minimum-energy solve cut the burn at %.6f of full, which is no ' ...
         'cut at all; there is then no second control and nothing to ' ...
         'minimise'],inDef.minEnergy.cutFrac);
    assert(inDef.minEnergy.propLeftKg > 50, ...
        'a cut-short burn must leave propellant unburned; this one left %.3f kg', ...
        inDef.minEnergy.propLeftKg);
    assertRel(inDef.minEnergy.cutFrac   ,0.99497604370117198,1e-4,'cutoff fraction');
    assertRel(inDef.minEnergy.tCutS     ,80.113240197794095 ,1e-4,'cutoff time (s)');
    assertRel(inDef.minEnergy.propLeftKg,150.71868896484500 ,1e-4, ...
        'propellant left unburned (kg)');
    assertAbs(inDef.minEnergy.tCutS,inDef.minEnergy.cutFrac.*inDef.tBurn,1e-6, ...
        'the cutoff time against the fraction times the full burn (s)');
    assert(inDef.minEnergy.tCutS < inDef.tBurn, ...
        'the cutoff at %.4f s must precede the %.4f s full burn', ...
        inDef.minEnergy.tCutS,inDef.tBurn);

%% THE CONSTRAINT IS MET, to the tolerance the INNER solve was given -- which is
%% tolRangeMEKm and not tolRangeKm, because here the range is a constraint on a
%% minimisation rather than the answer:
             tolMEM = 50;
    assert(abs(inDef.residM) <= tolMEM, ...
        ['the minimum-energy range residual is %.2f m against the %.1f m ' ...
         'constraint tolerance'],inDef.residM,tolMEM);
    assertRel(inDef.minEnergy.epsBo,-45555274.839652702,1e-6, ...
        'minimised burnout specific energy (J/kg)');

%% THE OBJECTIVE IS ACTUALLY MINIMISED, and this is the assertion the whole
%% rewrite turns on. Three claims, each of which a gamma-matched solve could
%% fail:
%%
%%   IT BEATS BOTH NEIGHBOURING FEASIBLE POINTS. Those are trajectories that
%%   make the SAME range at a different loft angle, so they are admissible
%%   competitors and not straw men. Measured: -45.4864 and -44.8189 MJ/kg
%%   against the solution's -45.5553.
    assert(inDef.minEnergy.epsBo < inDef.minEnergy.epsNeighLo && ...
           inDef.minEnergy.epsBo < inDef.minEnergy.epsNeighHi, ...
        ['the "minimum" burnout energy of %.6f MJ/kg is not below the ' ...
         '%.6f and %.6f MJ/kg of the two neighbouring FEASIBLE points. It is ' ...
         'not a minimum'],inDef.minEnergy.epsBo./1e6, ...
        inDef.minEnergy.epsNeighLo./1e6,inDef.minEnergy.epsNeighHi./1e6);

%%   AND IT BEATS BOTH FULL-BURN ARCS, which also make the required range and
%%   are the two answers the old selector chose between:
    assert(inDef.minEnergy.epsBo < inDef.depressed.epsBo && ...
           inDef.minEnergy.epsBo < inDef.lofted.epsBo, ...
        ['the minimum-energy arc leaves burnout at %.6f MJ/kg against the two ' ...
         'full-burn arcs'' %.6f and %.6f. A minimisation that cannot beat the ' ...
         'points it was seeded from has done nothing'], ...
        inDef.minEnergy.epsBo./1e6,inDef.depressed.epsBo./1e6, ...
        inDef.lofted.epsBo./1e6);

%%   AND THE MARGIN IS FAR LARGER THAN THE MEASURED NOISE. Every outer
%%   evaluation is an inner solve accurate only to the constraint tolerance, so
%%   the objective is sampled with an error; the script re-solves ten times
%%   tighter and reports the change. Measured: 101 J/kg of noise against the
%%   68911 J/kg shallower side of the valley, a ratio of 1.5e-3. The budget is
%%   1e-2, two orders above the propagation's own repeatability and still two
%%   orders below a valley that could be noise:
    assert(inDef.minEnergy.epsNoiseRel < 1e-2, ...
        ['the inner tolerance moves the burnout energy by %.4g J/kg against a ' ...
         '%.4g J/kg valley depth, a ratio of %.3e. At that ratio the minimum ' ...
         'is not resolved and tolRangeMEKm has to be tightened'], ...
        inDef.minEnergy.epsNoiseJkg, ...
        inDef.minEnergy.epsNoiseJkg./inDef.minEnergy.epsNoiseRel, ...
        inDef.minEnergy.epsNoiseRel);
    assertRel(inDef.minEnergy.epsNoiseJkg,101.00457912683500,0.5, ...
        'measured propagated noise of the constraint tolerance (J/kg)');

%% THE DIAGNOSTIC COMPARISON against the classical vacuum equal-radius arc, to a
%% tolerance taken FROM THE PHYSICS rather than chosen for comfort. The closed
%% form assumes an IMPULSIVE burn AT the impact radius in a VACUUM; this flight
%% burns for 80.1 s and finishes 82.2 km up at 5682.3 m/s rather than at zero
%% altitude at 5807.2 m/s. MEASURED: the apogee comes out 1.36 %% high, the
%% flight time 5.52 %% long and the burnout gamma 0.83 deg low. The budgets are
%% set above those and tight enough that the arc cannot wander to the 1453.7 km,
%% 25.8 min LOFTED arc (+53 %% and +28 %%) or to the 569.2 km, 17.5 min
%% DEPRESSED one (-40 %% and -13 %%), which are the two answers the old selector
%% gave. THESE ARE NOT RESIDUALS: nothing was driven to zero against them:
    assert(abs(inDef.meApoRelE) < 0.06, ...
        ['the flown apogee is %.3f km against the classical %.3f km, %.2f %% ' ...
         'out. A finite boost is worth a few per cent; this is not a few per ' ...
         'cent'],inDef.hApogee./1000,refKm.hApoM./1000,100.*inDef.meApoRelE);
    assert(abs(inDef.meTofRelE) < 0.09, ...
        ['the flown flight time is %.4f min against the classical %.4f min, ' ...
         '%.2f %% out'],inDef.tFlight./60,refKm.tofS./60,100.*inDef.meTofRelE);
    assert(abs(rad2deg(inDef.meGamResR)) < 2, ...
        ['the flown burnout gamma is %.4f deg from gammaStar. It is not driven ' ...
         'there, but a finite boost that ends 82 km up should still land ' ...
         'within a degree or two of the vacuum equal-radius value'], ...
        rad2deg(inDef.meGamResR));
    assertAbs(inDef.meApoRelE,(inDef.hApogee - refKm.hApoM)./refKm.hApoM,1e-12, ...
        'the reported apogee residual against this file');
    assertAbs(inDef.meTofRelE,(inDef.tFlight - refKm.tofS)./refKm.tofS,1e-12, ...
        'the reported flight-time residual against this file');

%% ...AND THE DIAGNOSTIC IS ATTRIBUTABLE, not merely small. Almost the whole gap
%% is the FINITE BURNOUT ALTITUDE: a Keplerian arc from the flown burnout state
%% apogees within 30 m of what the vehicle flew, so what the classical arc gets
%% wrong is where the boost ENDED and not how the solve converged:
    assertAbs(inDef.minEnergy.hApoKepM,inDef.hApogee,1000, ...
        'Keplerian apogee of the burnout state against the flown apogee (m)');
    assert(abs(inDef.minEnergy.hApoKepM - refKm.hApoM) > ...
           0.9.*abs(inDef.hApogee - refKm.hApoM), ...
        ['the finite burnout altitude accounts for only %.3f km of the %.3f ' ...
         'km gap to the classical arc; the attribution the summary prints is ' ...
         'wrong'],abs(inDef.minEnergy.hApoKepM - refKm.hApoM)./1000, ...
        abs(inDef.hApogee - refKm.hApoM)./1000);

%% ...and the arc really is a minimum-energy one rather than either full-burn
%% arc. The apogee-to-range ratio of a minimum-energy ballistic trajectory sits
%% near a quarter; measured here it is 0.200, against the lofted arc's 0.301 and
%% the depressed arc's 0.118:
             ratME = inDef.hApogee./inDef.rngAchM;
    assert(ratME > 0.15 && ratME < 0.26, ...
        ['the flown apogee-to-range ratio is %.4f; a minimum-energy arc sits ' ...
         'near a quarter and the two full-burn arcs here sit at %.4f and %.4f'], ...
        ratME,inDef.depressed.hApoM./inDef.rngAchM, ...
        inDef.lofted.hApoM./inDef.rngAchM);

%% ---------------------------------------------------------------------
%% 6b. THE VACUUM EQUAL-RADIUS LIMIT: the objective REPRODUCES the classical arc
%% ---------------------------------------------------------------------
%% THIS IS THE CHECK THAT SAYS THE OBJECTIVE IS THE RIGHT ONE, and it is the
%% only one in this file that grades the FORMULATION rather than the code.
%%
%% Take the boost to be IMPULSIVE at r = rE in a VACUUM, which is the limit the
%% classical result is derived in. The feasible family is then the
%% one-parameter family of Keplerian arcs of free-flight central angle psi from
%% rE back to rE: pick a burnout flight-path angle, and the speed that makes the
%% arc close on psi follows. Minimising the SAME objective the script minimises,
%% eps = V^2/2 - mu/rE, over that family must return exactly
%%
%%     gamma = 45 deg - psi/4   and   V = V*,
%%
%% because that IS the classical minimum-energy trajectory. Everything below is
%% written out in this file's own hand -- the free-flight range relation, the
%% golden-section minimisation, the closed-form answer -- so a mutated formula
%% in the script cannot be mirrored by it:
       [gamVac,vVac,nVac] = vacuumMinEnergy(psiReq,c);
    assert(abs(gamVac - gamStr) < 1e-7, ...
        ['minimising the burnout energy over the VACUUM EQUAL-RADIUS feasible ' ...
         'family gave a burnout flight-path angle of %.9f deg against the ' ...
         'classical %.9f deg, %.3e rad apart against a 1e-7 rad budget. If the ' ...
         'objective were the wrong one, this is where it would show'], ...
        rad2deg(gamVac),rad2deg(gamStr),abs(gamVac - gamStr));
    assertRel(vVac,refKm.V,1e-9, ...
        'vacuum-limit minimising burnout speed against V* (m/s)');
    assert(nVac > 20, ...
        ['the vacuum-limit minimisation converged in %d golden-section steps, ' ...
         'which is too few to have searched anything'],nVac);

%% ...and the limit is DISCRIMINATING: a minimiser that simply returned the
%% middle of its bracket, or the classical value read back from a formula, would
%% have to be told apart from one that searched. The bracket is deliberately
%% wide and its midpoint is far from the answer:
    assert(abs(rad2deg(gamVac - deg2rad(45))) > 5, ...
        ['the vacuum-limit answer is %.4f deg, within 5 deg of the 45 deg a ' ...
         'gammaStar mutated to a constant would give; this geometry cannot ' ...
         'discriminate'],rad2deg(gamVac));

%% ...and the SAME limit check on a second, very different range angle, so the
%% agreement above cannot be a coincidence of one geometry:
            psiTwo = deg2rad(100);
       [gamTwo,vTwo] = vacuumMinEnergy(psiTwo,c);
            refTwo = classicalArc(psiTwo,c);
    assertAbs(gamTwo,pi./4 - psiTwo./4,1e-7, ...
        'vacuum-limit burnout flight-path angle at a 100 deg range angle (rad)');
    assertRel(vTwo,refTwo.V,1e-9, ...
        'vacuum-limit burnout speed at a 100 deg range angle (m/s)');

%% ---------------------------------------------------------------------
%% 7. THE SOLVE IS NOT A CONSTANT, and it is neither full-burn arc
%% ---------------------------------------------------------------------
%% One case cannot tell a working minimisation from a hard-wired answer. A
%% second target, 1653 km CLOSER, must move both the loft angle and the cutoff
%% fraction, must minimise its own objective, and must meet its own diagnostic
%% reference -- a different gammaStar, a different apogee and a different flight
%% time.
%%
%% IT IS CLOSER RATHER THAN FURTHER BECAUSE THE DEPRESSED BAND IS ONLY 500 km
%% WIDE at the vehicle's 6 deg clamp -- 4708 to 5212 km -- so no second target
%% a thousand kilometres away can have both full-burn branches. This one has
%% ONLY the lofted arc, which is the point: it also exercises the one-branch
%% path, where the minimisation takes the user's own loftMin as the end of the
%% feasible family because the full burn OVERSHOOTS the target there:
            outFar2 = evalc(['[trFar2,inFar2] = run_ballistic_target(struct(' ...
                             '''lonTarget'',-60,''showPlots'',false));']);
    assert(contains(outFar2,'(nominal)') && ~contains(outFar2,'REFUSED'), ...
        'the 3175 km minimum-energy case did not solve. Summary was:\n%s',outFar2);
    assert(inDef.rngReqM > inFar2.rngReqM + 1e6, ...
        'the second minimum-energy target is not materially closer in');
    assert(~inFar2.depressed.exists && inFar2.lofted.exists, ...
        ['the second case is meant to have exactly ONE full-burn branch, so ' ...
         'that the minimisation is exercised where one end of the feasible ' ...
         'family is the user''s loft limit rather than a branch solution']);

%% BOTH unknowns moved. A solve wired to a constant loft angle, or one that
%% never touched the cutoff, passes the shipped case and fails here:
    assert(abs(inFar2.loftDeg - inDef.loftDeg) > 1, ...
        ['the two minimum-energy targets solved to loft angles %.4f and %.4f ' ...
         'deg, less than a degree apart'],inFar2.loftDeg,inDef.loftDeg);
    assert(abs(inFar2.minEnergy.cutFrac - inDef.minEnergy.cutFrac) > 0.01, ...
        ['the two targets solved to cutoff fractions %.6f and %.6f; the ' ...
         'cutoff enforces the range constraint and must move with the range'], ...
        inFar2.minEnergy.cutFrac,inDef.minEnergy.cutFrac);
    assert(abs(rad2deg(wrapPi(inFar2.psiLaunch - inDef.psiLaunch))) > 1, ...
        'the two targets returned the same azimuth to within a degree');

%% AND THE SECOND CASE IS ALSO A MINIMUM, on its own family. The objective is
%% what this mode claims to deliver, so it is asserted on every case that solves
%% and not only on the shipped one:
    assert(inFar2.minEnergy.epsBo < inFar2.minEnergy.epsNeighLo && ...
           inFar2.minEnergy.epsBo < inFar2.minEnergy.epsNeighHi, ...
        ['the 3175 km case''s burnout energy of %.6f MJ/kg is not below the ' ...
         '%.6f and %.6f MJ/kg of its neighbouring feasible points'], ...
        inFar2.minEnergy.epsBo./1e6,inFar2.minEnergy.epsNeighLo./1e6, ...
        inFar2.minEnergy.epsNeighHi./1e6);
    assert(inFar2.minEnergy.epsBo < inFar2.lofted.epsBo, ...
        ['the 3175 km case''s burnout energy of %.6f MJ/kg does not beat the ' ...
         'one full-burn arc that reaches it, at %.6f MJ/kg'], ...
        inFar2.minEnergy.epsBo./1e6,inFar2.lofted.epsBo./1e6);
    assert(inFar2.minEnergy.epsNoiseRel < 1e-2, ...
        'the 3175 km case''s valley is only %.3e above its own noise floor', ...
        inFar2.minEnergy.epsNoiseRel);

%% ...and the SECOND case meets its OWN diagnostic reference, recomputed here
%% rather than reused. gammaStar RISES as the range angle shrinks -- 34.157 deg
%% at 4828 km against 37.870 deg at 3175 km -- so a gammaStar wired to a
%% constant passes the shipped case and fails this one:
           psiFar2 = coorbital.util.greatCircle(latLn,lonLn,latTg,deg2rad(-60));
           gamStr2 = pi./4 - psiFar2./4;
           refFar2 = classicalArc(psiFar2,c);
    assert(abs(rad2deg(gamStr2 - gamStr)) > 2, ...
        ['the two gammaStar values are %.4f and %.4f deg, only %.4f deg ' ...
         'apart; this case cannot see a gammaStar wired to a constant'], ...
        rad2deg(gamStr2),rad2deg(gamStr),abs(rad2deg(gamStr2 - gamStr)));
    assertAbs(inFar2.gamStarR,gamStr2,1e-12, ...
        'the 3175 km case''s own gammaStar against this file (rad)');
    assertAbs(inFar2.gamBurnout,gamStr2,deg2rad(2), ...
        ['the 3175 km burnout gamma against ITS OWN gammaStar (rad). NOT a ' ...
         'residual: nothing drives it there, and the budget is the size of the ' ...
         'finite-boost modelling gap']);
    assert(abs(inFar2.meApoRelE) < 0.06, ...
        ['the 3175 km apogee is %.3f km against the classical %.3f km, ' ...
         '%.2f %% out'],inFar2.hApogee./1000,refFar2.hApoM./1000, ...
        100.*inFar2.meApoRelE);
    assert(abs(inFar2.meTofRelE) < 0.10, ...
        ['the 3175 km flight time is %.4f min against the classical %.4f ' ...
         'min, %.2f %% out'],inFar2.tFlight./60,refFar2.tofS./60, ...
        100.*inFar2.meTofRelE);
    assert(inFar2.missM <= tolM, ...
        'the 3175 km case missed by %.2f m',inFar2.missM);
            missIF = c.rE.*coorbital.util.greatCircle(trFar2.x(end,3), ...
                                                      trFar2.x(end,2), ...
                                                      latTg,deg2rad(-60));
    assert(missIF <= tolM, ...
        'the 3175 km case flown impact point is %.2f m from its target',missIF);
    assertBranchMiss(inFar2,tolM,'the 3175 km minimum-energy case');

%% ...AND THE MINIMUM-ENERGY ARC IS NEITHER OF THE TWO FULL-BURN ARCS on the
%% shipped case. It lies BETWEEN them in apogee and in flight time, which is
%% what a trajectory carrying less burnout energy than either of them must do,
%% so a mode that quietly fell back to selecting a branch could not pass this.
%% The 3175 km case has only ONE full-burn arc and is checked against that one
%% above instead:
    assertBetweenArcs(inDef,'the shipped 4828 km case');
    assert(inFar2.hApogee < 0.5.*inFar2.lofted.hApoM, ...
        ['the 3175 km minimum-energy arc apogees at %.3f km against the only ' ...
         'full-burn arc''s %.3f km; it must be clearly below it'], ...
        inFar2.hApogee./1000,inFar2.lofted.hApoM./1000);

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
    assertRel(infFar.rngMaxM,5211525.2734660003,1e-4,'refused-case maximum range (m)');

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
  [outBr,infBr,threwBr,trBr] = tryRun(struct('lonTarget',-60, ...
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
%% 9. THE GUARDS: the vehicle clamp, the loft bracket, the branch word
%% ---------------------------------------------------------------------
%% THE CLAMP IS A VEHICLE LIMIT AND BOTH BM SCRIPTS READ IT. Until 08/08/2026
%% BM/run_ballistic carried 6 deg and this script carried 12 deg for nominally
%% the same airframe, and the 12 had been chosen to bring the old demonstration
%% target inside the depressed branch. Two limits for one vehicle made the two
%% scripts' performance non-comparable, and a control-authority limit chosen for
%% reachability is not a limit at all:
                vB1 = vehicle_bm();
    assert(isfield(vB1,'alphaMaxDeg'), ...
        ['BM/vehicle_bm must define alphaMaxDeg. It is where a ' ...
         'control-authority limit belongs, and both BM entry scripts read it ' ...
         'from there']);
    assertAbs(vB1.alphaMaxDeg,6,1e-12, ...
        ['BM/vehicle_bm must carry the angle-of-attack clamp at the 6 deg ' ...
         'BM/run_ballistic has always flown (deg)']);
    assertAbs(inDef.alphaMaxDeg,vB1.alphaMaxDeg,1e-12, ...
        'the clamp this script flew against the vehicle''s own (deg)');

%% THE SHIPPED LOFT BRACKET HAS TO REACH BELOW -42.9 deg, because that is where
%% the max-range attitude sits at the vehicle's clamp. Put back the -40 deg
%% loftMin this script used to ship and the largest range on the bracket is
%% found AT an endpoint, so range is monotone across all of it, one branch is
%% empty everywhere and the two-branch structure does not exist. That must be
%% REFUSED with a named identifier, not reported as a two-branch problem with
%% one branch quietly missing:
             idClp = errorIdOf(@() run_ballistic_target( ...
                                    struct('loftMin',-40,'showPlots',false)));
    assert(strcmp(idClp,'coorbital:runBallisticTarget:maximumNotBracketed'), ...
        ['a -40 deg loftMin leaves the max-range loft angle outside the ' ...
         'bracket at the vehicle''s 6 deg clamp, and there is then only one ' ...
         'branch on it; expected ' ...
         'coorbital:runBallisticTarget:maximumNotBracketed, got "%s"'],idClp);

%% THE SHIPPED BRACKET DOES FIND IT, and the numbers are pinned so that nobody
%% narrows loftMin back without meeting this:
    assert(inDef.loftStarDeg < -40, ...
        ['the max-range loft angle is %.4f deg and must lie BELOW -40 deg, or ' ...
         'the refusal above is not about bracket width at all'], ...
        inDef.loftStarDeg);
    assert(inDef.loftStarDeg > -43.5 && inDef.loftStarDeg < -42.3, ...
        ['the max-range loft angle came out at %.4f deg against the ' ...
         '-42.907 deg measured on the committed code'],inDef.loftStarDeg);
    assertRel(inDef.depressed.bandLo,4708462.6516564805,1e-4, ...
        'depressed-branch floor at the vehicle clamp (m)');
    assert(inDef.rngReqM > inDef.depressed.bandLo && ...
           inDef.rngReqM < inDef.rngMaxM, ...
        ['the SHIPPED TARGET must lie inside the depressed band %.3f to ' ...
         '%.3f km, or the shipped configuration cannot fly all three modes. ' ...
         'It is %.3f km'],inDef.depressed.bandLo./1000,inDef.rngMaxM./1000, ...
        inDef.rngReqM./1000);

%% WHAT RAISING THE CLAMP WOULD BUY AND WHAT IT WOULD COST, measured rather than
%% argued. 12 deg gives the depressed branch reach down to 1684 km and pays
%% 156.224 km of maximum range for it. NEITHER number is a reason to move a
%% vehicle limit; they are here so that the trade is a number in the record and
%% so that raising the clamp cannot be done quietly:
             out12 = evalc(['[~,in12] = run_ballistic_target(struct(' ...
                            '''alphaMax'',12,''loftMin'',-40,''nScanLoft'',13,' ...
                            '''showPlots'',false));']);
    assert(~in12.refused && ~contains(out12,'REFUSED'), ...
        'the 12 deg sensitivity case must solve. Summary was:\n%s',out12);
    assertAbs(in12.alphaMaxDeg,12,1e-12, ...
        'the override must actually reach the equations of motion (deg)');
    assertRel(in12.rngMaxM,5055301.6893044896,1e-4,'12 deg maximum range (m)');
    assertRel(in12.depressed.bandLo,1684116.7693332899,1e-4, ...
        '12 deg depressed-branch floor (m)');
    assert(inDef.rngMaxM > in12.rngMaxM, ...
        ['raising the clamp from 6 to 12 deg must COST maximum range: ' ...
         '%.3f km at 6 deg against %.3f km at 12 deg'], ...
        inDef.rngMaxM./1000,in12.rngMaxM./1000);
    assertAbs((inDef.rngMaxM - in12.rngMaxM)./1000,156.2236,0.01, ...
        'maximum range given up to raise the clamp from 6 to 12 deg (km)');
    assert(in12.depressed.bandLo < 0.5.*inDef.depressed.bandLo, ...
        ['12 deg must buy real depressed-branch reach: %.3f km against the ' ...
         'vehicle clamp''s %.3f km'],in12.depressed.bandLo./1000, ...
        inDef.depressed.bandLo./1000);

%% ...AND THE MINIMUM-ENERGY ANSWER IS A PROPERTY OF THE GEOMETRY, not of the
%% clamp. The two clamps must reach very nearly the same ARC from very different
%% COMMANDED loft angles, because the clamp changes what a command delivers:
    assert(in12.minEnergy.solved, ...
        'the 12 deg case did not solve the minimum-energy trajectory');
    assert(abs(in12.hApogee - inDef.hApogee) < 0.05.*inDef.hApogee, ...
        ['the 12 deg clamp found a %.3f km apogee against the vehicle clamp''s ' ...
         '%.3f km. The minimum-energy arc is set by the GEOMETRY and the two ' ...
         'must very nearly agree'],in12.hApogee./1000,inDef.hApogee./1000);
    assert(abs(in12.loftDeg - inDef.loftDeg) > 10, ...
        ['the two clamps solved to loft angles %.4f and %.4f deg. They must ' ...
         'DIFFER: the clamp decides what a commanded attitude delivers, and ' ...
         'if they agreed the loft angle would not be doing any work'], ...
        in12.loftDeg,inDef.loftDeg);

%% A ONE-BRANCH minimum-energy run must still SOLVE and must not print
%% arithmetic it does not have. The contrast paragraph compares the two
%% full-burn burnout energies and a missing branch supplies NaN; unguarded it
%% read "leave burnout at NaN and -44.9073 MJ/kg, NaN %% apart", which in a
%% results paragraph reads as a measurement:
    assert(~contains(outFar2,'NaN'), ...
        ['the one-branch run printed a NaN, which in a results paragraph ' ...
         'reads as a measurement. Summary was:\n%s'],outFar2);
    assert(contains(outFar2,'ONLY ONE FULL-BURN ARC REACHES THIS TARGET'), ...
        ['a minimum-energy run with only one full-burn arc must say so ' ...
         'instead of printing half a comparison. Summary was:\n%s'],outFar2);

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
%% the depressed arc's 569.177 km apogee earns 3x where the minimum-energy arc's
%% 965.613 km earns 2x. Asked of the rule with the flown apogee rather than by
%% flying a second auto run:
    assert(inDef.altExagRule(maxAlt(trDep,c.rE)) == 3, ...
        ['the rule gives %gx for the depressed arc''s %.3f km apogee against ' ...
         '2x for the shipped arc''s; one fixed factor cannot serve both'], ...
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
    assert(strcmp(inDef.targetStr,'62.00N 28.00W'), ...
        'the target caption reads "%s" against the expected "62.00N 28.00W"', ...
        inDef.targetStr);

%% ---------------------------------------------------------------------
%% 11. THE BRANCH MEASUREMENT ITSELF, exercised where it CAN disagree
%% ---------------------------------------------------------------------
%% Part 5 re-derives the branch from the flown loft angle against the certified
%% maximiser interval, but every run it does that on is a run where the
%% measurement AGREES with the label. So every assertion in it also passes if
%% branchOfLoft is never called: replace the call in the "WHICH BRANCH DID IT
%% ACTUALLY FLY" block with flownName = pickName and branchOK = true, and the
%% whole of the rest of this file still passes while info.branchMeasured,
%% info.branchAgrees and the summary's "MEASURED as" line all become restatements
%% of the commanded branch. This part is what closes that, and the headline claim
%% of the script -- THE BRANCH IS MEASURED, NOT ASSUMED -- is what it is closing.
%%
%% FINDING A CASE WHERE THE TWO CAN LEGITIMATELY DIFFER TAKES SOME CARE, and the
%% reason is worth recording. Each branch is solved on an interval that lies
%% entirely on one CERTIFIED side of the maximiser interval, so a solution
%% strictly inside its own bracket always measures as the branch it was solved
%% on. What breaks that is the ENDPOINT: when the required range comes within the
%% range tolerance of the largest range the search can certify,
%% coorbital.util.rangeSolve short-circuits at the bracket endpoint, both branch
%% solves come back holding an arc INSIDE the unresolved maximiser interval, and
%% the position test can only answer COALESCED -- because at that range the two
%% arcs are the same arc and the distinction has no content.
%%
%% The case below asks for 5212.791 km on a 50 km tolerance against a
%% 5211.525 km certified maximum, and asks for the LOFTED arc by name:
             outMx = evalc(['[trMx,inMx] = run_ballistic_target(struct(' ...
                            '''latTarget'',55,''lonTarget'',-25,' ...
                            '''tolRangeKm'',50,' ...
                            '''branch'',''lofted'',''showPlots'',false));']);
    assert(~inMx.refused && ~isempty(trMx), ...
        'the near-maximum-range case must solve. Summary was:\n%s',outMx);

%% THE PREMISE FIRST: the solve really did land inside the unresolved interval,
%% and the required range really is inside the tolerance of the certifiable
%% maximum. If a future change moves either, this part stops discriminating and
%% must say so rather than pass vacuously:
    assert(inMx.atMaxRange, ...
        ['the near-maximum case did not land INSIDE the certified maximiser ' ...
         'interval: it flew %.10f deg against %.10f to %.10f deg, so it cannot ' ...
         'exercise the measurement at all'],inMx.loftDeg, ...
        rad2deg(inMx.maxRange.aL),rad2deg(inMx.maxRange.bL));
    assert(inMx.coalescedRq, ...
        ['the required %.3f km must fall in the unresolved band above ' ...
         '%.3f km for this case to mean anything'], ...
        inMx.rngReqM./1000,inMx.rngTopM./1000);
    assertAbs(deg2rad(inMx.loftDeg),inMx.maxRange.bL,1e-12, ...
        ['the lofted solve must have short-circuited AT the upper end of the ' ...
         'certified interval (rad)']);

%% THE LABEL AND THE MEASUREMENT DISAGREE, and the measurement is not a copy of
%% the label. This is the assertion the whole part exists for:
    assert(strcmp(inMx.branchFlown,'lofted'), ...
        'the run was asked for the lofted branch and flew "%s"',inMx.branchFlown);
    assert(strcmp(inMx.branchMeasured,'coalesced'), ...
        ['the flown arc sits INSIDE the certified maximiser interval, where no ' ...
         'evidence this search has can tell the two branches apart, so it must ' ...
         'measure as "coalesced"; it measured as "%s". A measurement wired to ' ...
         'the commanded branch says "lofted" here, and this is the only ' ...
         'assertion in this file that sees the difference'],inMx.branchMeasured);
    assert(~inMx.branchAgrees, ...
        'info.branchAgrees must be FALSE when the measured branch is not the commanded one');
    assert(strcmp(inMx.branchByApogee,'lofted') && ...
           strcmp(inMx.branchByTime,'lofted'), ...
        ['the DESCRIPTIVE apogee and time readings both say "lofted" here, ' ...
         'because the flown arc is a hair above the max-range one in both. ' ...
         'THAT IS EXACTLY WHY THEY CANNOT CLASSIFY: a difference of 129 m of ' ...
         'apogee and 0.07 s of flight time is search noise, not a branch, and ' ...
         'the old measurement read it as one']);

%% ...AND THE CAUTION MUST BE PRINTED, not merely flagged in the struct. The
%% summary is what a reader acts on, and it must diagnose the degeneracy rather
%% than blaming a bracketing step that did nothing wrong:
    assert(contains(outMx,'*** CAUTION ***'), ...
        'the disagreement printed no caution. Summary was:\n%s',outMx);
    assert(contains(outMx,'MEASURES as coalesced'), ...
        'the caution must name the MEASURED branch. Summary was:\n%s',outMx);
    assert(contains(outMx,'landed INSIDE the certified maximiser interval'), ...
        ['the caution must diagnose the degeneracy against the certified ' ...
         'interval, not blame the bracketing. Summary was:\n%s'],outMx);
    assert(contains(outMx,'bracketing is NOT suspect'), ...
        ['the caution must exonerate the bracketing explicitly; the old one ' ...
         'told the reader to suspect the max-range angle. Summary was:\n%s'],outMx);
    assert(contains(outMx,'MEASURED as      coalesced'), ...
        ['the summary''s MEASURED line must carry the measured branch and ' ...
         'not the commanded one. Summary was:\n%s'],outMx);
    assertAbs(summaryNumber(outMx, ...
        'MEASURED as +coalesced +\(from the flown ([-\d.]+) deg loft angle'), ...
        inMx.loftDeg,1e-3,'printed measured-branch loft angle (deg)');
    assert(contains(outMx,'UNRESOLVED band around the maximum'), ...
        ['a required range inside the unresolved band must be called out where ' ...
         'the hump is reported, not only in the branch caution. Summary ' ...
         'was:\n%s'],outMx);

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
%% opinion: 7132.320 km against a 5211.525 km maximum, 1920.794 km short. It
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
    assert(infPny.rngReqM > infPny.rngMaxM + 1.5e6, ...
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

%% GIVE IT A BOOSTER SIZED TO THE RANGE AND THE MODE LANDS ON THE REFERENCE,
%% which is the acceptance test for the whole mode against a published geometry.
%% The classical result PRESUMES the burnout energy is free to choose -- you size
%% the booster to the range -- and the shipped placeholder is sized to about
%% 5200 km, so the case is flown here on a booster with the delta-V for 7132 km.
%% Nothing else changes: same vehicle, same pitch programme, same clamp, same
%% solver settings. MEASURED: 1218.445 km of apogee against 1205.989 km
%% (+1.03 %%) and 27.161 min against 26.012 min (+4.42 %%), with the burnout
%% gamma 0.557 deg from gammaStar. The budgets are 5 %%, 8 %% and 2 deg, set
%% from those measurements the same way part 6 sets its own -- and the gamma one
%% is a DIAGNOSTIC budget, since nothing drives the burnout gamma there:
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
    assertAbs(infBig.gamBurnout,gamPny,deg2rad(2), ...
        ['sized-booster burnout gamma against gammaStar (rad). A DIAGNOSTIC ' ...
         'budget: nothing drives the burnout gamma there, and a finite boost ' ...
         'ending 83 km up cannot be expected to sit on a value derived for an ' ...
         'impulsive burn at the surface']);
    assert(infBig.minEnergy.epsBo < infBig.minEnergy.epsNeighLo && ...
           infBig.minEnergy.epsBo < infBig.minEnergy.epsNeighHi, ...
        ['the sized-booster case''s burnout energy of %.6f MJ/kg is not below ' ...
         'the %.6f and %.6f MJ/kg of its neighbouring feasible points'], ...
        infBig.minEnergy.epsBo./1e6,infBig.minEnergy.epsNeighLo./1e6, ...
        infBig.minEnergy.epsNeighHi./1e6);
    assert(abs(infBig.meApoRelE) < 0.05, ...
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
                                   struct('tolLoftMEDeg',0,'showPlots',false)));
    assert(~isempty(idTol),'tolLoftMEDeg = 0 must raise');

%% AND THE CONSTRAINT TOLERANCE MUST NOT BE LOOSER THAN THE BRANCH TOLERANCE.
%% It is the inner solve of a minimisation, so its error is NOISE on the
%% objective; a constraint looser than the range tolerance would make the
%% loosest thing in the run the one thing that has to hold:
            idMEt = errorIdOf(@() run_ballistic_target( ...
                                   struct('tolRangeMEKm',5,'tolRangeKm',1, ...
                                          'showPlots',false)));
    assert(~isempty(idMEt), ...
        'tolRangeMEKm above tolRangeKm must raise, not be quietly accepted');

%% THE MINIMISATION TOLERANCE IS ACTUALLY HONOURED, which nothing above can see:
%% every assertion so far is written against the SHIPPED 0.05 deg. Tighten it by
%% a factor of ten and the search must take more steps, close a narrower bracket
%% and find an objective no WORSE than the shipped one -- which a solve that
%% ignored the tolerance, or one that stopped after a fixed number of steps,
%% could not do:
          outTight = evalc(['[~,inTight] = run_ballistic_target(struct(' ...
                            '''tolLoftMEDeg'',0.005,''showPlots'',false));']);
    assert(~inTight.refused,'the tightened run did not solve:\n%s',outTight);
    assert(rad2deg(inTight.minEnergy.widthR) <= 0.005 + 1e-12, ...
        ['asked for 0.005 deg on the loft bracket the minimisation closed to ' ...
         '%.6f deg'],rad2deg(inTight.minEnergy.widthR));
    assert(inTight.minEnergy.nGolden > inDef.minEnergy.nGolden, ...
        ['a ten times tighter loft tolerance took %d golden-section steps ' ...
         'against the shipped %d; it cannot cost the same'], ...
        inTight.minEnergy.nGolden,inDef.minEnergy.nGolden);
    assert(inTight.minEnergy.epsBo <= inDef.minEnergy.epsBo + 1e-6, ...
        ['the tightened run found %.9f MJ/kg against the shipped run''s ' ...
         '%.9f MJ/kg. A finer search of the SAME objective cannot end up ' ...
         'higher'],inTight.minEnergy.epsBo./1e6,inDef.minEnergy.epsBo./1e6);

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

function [gamStarR,vStarM,nStep] = vacuumMinEnergy(psiR,c)
%% Purpose:
%
%  Solve the SAME constrained minimisation BM/run_ballistic_target solves, in
%  the VACUUM EQUAL-RADIUS LIMIT where the answer is known in closed form:
%
%      minimise    eps = V^2/2 - mu/rE
%      over        the burnout flight-path angle gamma
%      subject to  the free-flight arc from rE returning to rE after a central
%                  angle of psi.
%
%  This is the check that grades the OBJECTIVE rather than the code. If the
%  quantity being minimised is the right one, this must return the classical
%  gammaStar = 45 deg - psi/4 and V*; if it is not, no amount of solver
%  correctness downstream will make the mode mean what its name says.
%
%  THE CONSTRAINT IS CLOSED FORM HERE, which is what makes the check tight. At
%  r = rE, with q = p/rE the scaled semi-latus rectum,
%
%      e cos(nu) = q - 1,     e sin(nu) = q tan(gamma),
%
%  so the true anomaly at burnout is nu = atan2(q tan gamma, q - 1). The arc is
%  symmetric about apogee, so returning to rE after a central angle psi means
%  nu = pi - psi/2, and that inverts:
%
%      q = tan(nu) / (tan(nu) - tan(gamma)),   V = sqrt(mu q rE) / (rE cos gamma).
%
%  Nothing is integrated and nothing is bisected against a propagation, so the
%  only error is floating point.
%
%  THE MINIMISATION IS A GOLDEN SECTION over gamma, deliberately started on a
%  bracket whose midpoint is nowhere near the answer, so a routine that returned
%  its own midpoint could not pass.
%
%% Inputs:
%
%  psiR             [1 x 1]                     Free-flight range angle (rad),
%                                               strictly between 0 and pi
%
%  c                Struct                      coorbital.util.missileConst
%
%% Outputs:
%
%  gamStarR         [1 x 1]                     Minimising burnout flight-path
%                                               angle (rad)
%
%  vStarM           [1 x 1]                     Burnout speed there (m/s)
%
%  nStep            [1 x 1]                     Golden-section steps taken
%
%% References:
%   [1] Bate, R.R., Mueller, D.D., White, J.E., "Fundamentals of
%       Astrodynamics," Dover, 1971, Ch. 6.
%
%% Revision History:
%  Michael Casey                                                08/08/2026
%  Copyright 2026 Coorbital, Inc.
%% ------------------------ Begin Code Sequence ---------------------------

                nuR = pi - psiR./2;
              vOfGam = @(g) sqrt(c.muE.*qOf(g,nuR).*c.rE)./(c.rE.*cos(g));
              epsOf = @(g) vOfGam(g).^2./2 - c.muE./c.rE;
                 aG = deg2rad(1);
                 bG = deg2rad(85);
               gRat = (sqrt(5) - 1)./2;
                 x1 = bG - gRat.*(bG - aG);
                 x2 = aG + gRat.*(bG - aG);
                 e1 = epsOf(x1);
                 e2 = epsOf(x2);
              nStep = 0;
    while (bG - aG) > 1e-10 && nStep < 500
              nStep = nStep + 1;
        if e1 < e2
                 bG = x2;
                 x2 = x1;
                 e2 = e1;
                 x1 = bG - gRat.*(bG - aG);
                 e1 = epsOf(x1);
        else
                 aG = x1;
                 x1 = x2;
                 e1 = e2;
                 x2 = aG + gRat.*(bG - aG);
                 e2 = epsOf(x2);
        end
    end
           gamStarR = aG./2 + bG./2;
             vStarM = vOfGam(gamStarR);
end

function q = qOf(gamR,nuR)
%% Purpose:
%
%  The scaled semi-latus rectum p/rE of the free-flight arc that leaves r = rE
%  at flight-path angle gamR and reaches true anomaly nuR there. Inverts
%  tan(nu) = q tan(gamma) / (q - 1), which is the equal-radius free-flight
%  constraint written at the burnout point.
%
%% Inputs:
%
%  gamR             [1 x 1]                     Burnout flight-path angle (rad)
%
%  nuR              [1 x 1]                     True anomaly at burnout (rad),
%                                               in (pi/2, pi) for an ascending
%                                               burnout on a symmetric arc
%
%% Outputs:
%
%  q                [1 x 1]                     p/rE (-)
%
%% Revision History:
%  Michael Casey                                                08/08/2026
%  Copyright 2026 Coorbital, Inc.
%% ------------------------ Begin Code Sequence ---------------------------

                 q = tan(nuR)./(tan(nuR) - tan(gamR));
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
%  BOTH SEPARATIONS ARE CHECKED AS FACTORS rather than as differences, but the
%  factors DO have to be re-derived when the vehicle changes: at the 12 deg
%  clamp this script used to ship, the minimum-energy arc sat at 0.51 of the
%  lofted apogee, and at the vehicle's 6 deg it sits at 0.66, because a tighter
%  clamp compresses the whole family towards the max-range arc. Measured on the
%  shipped case: 965.613 km against 569.177 km depressed and 1453.661 km lofted.
%  The thresholds below are half again the depressed apogee and three quarters
%  of the lofted one, with the flight time strictly inside the pair.
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
    assert(info.hApogee < 0.75.*info.lofted.hApoM, ...
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
