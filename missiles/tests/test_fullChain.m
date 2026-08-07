function test_fullChain()
%% Purpose:
%
%  Pin the COMPOSITION of the full boost-glide-descent chain: the three phases
%  flying together on one seven-state vector, the two seams between them, the
%  three-vehicle binding, the mass bookkeeping across stage separation, and the
%  numbers the entry script prints. Every other test in this suite exercises
%  one library function, or one shorter chain; this one is the milestone's
%  headline deliverable -- boost, glide and descent in a single flight.
%
%  Without it the suite passes unchanged while run_boost_glide drops the
%  separation link, flies the glide on the booster's airframe, loses time
%  across a junction, runs the descent before the glide, transposes latitude
%  and longitude, or reports the wrong range.
%
%% Note -- how the seams are checked, and why a re-run would not do:
%
%  Part 4 re-integrates ACROSS each junction with ode89 at 1e-12, a different
%  solver at a tighter tolerance than the driver's ode45 at 1e-10, starting
%  from a driver sample two seconds before the junction and finishing at a
%  driver sample two seconds after it, applying the phase link in between.
%
%  THE CLAIM THAT MATTERS IS THE TOLERANCE ORDERING, not the solver name: a
%  reference no tighter than the driver measures its own error as much as the
%  driver's, and the residual then bounds nothing. That ordering is asserted in
%  crossJunction rather than left to a comment, so a later edit cannot
%  downgrade the reference -- or tighten the driver through env.odeRelTol --
%  and leave this check silently vacuous.
%
%  A stronger-sounding argument was tried and does NOT hold for this anchoring
%  scheme, so it is recorded here rather than repeated: re-running the driver's
%  own ode45 at its own 1e-10 does not reproduce the trajectory bitwise and
%  does not give a zero residual. It gives 1.337e-06 m at junction 1 against
%  the ode89 reference's 1.340e-06, agreeing to five digits. The bitwise trap
%  of the glide milestone (docs/LESSONS_LEARNED.md) needs the reference to
%  repeat the driver's STEP SEQUENCE, which happens when it restarts from the
%  same state at the same phase boundary; anchoring on an interior sample two
%  seconds inside the phase breaks that sequence, so the circularity is absent
%  here by accident of the anchoring. Relying on the accident would be
%  careless. The tolerance ordering is the load-bearing part.
%
%  The reference reaches the equations of motion through info.phases, which
%  run_boost_glide returns for this purpose. What is independent is the
%  integration, not the physics: nothing here claims to re-derive the equations
%  of motion, only to confirm that the chain's seams are seams and not steps.
%
%% Note -- the chain is NOT continuous in all seven states, and must not be:
%
%  The booster is jettisoned at burnout, so component 7 jumps DOWN by
%  bst.massDry at the first junction. That is physical. Part 3 asserts the
%  expected jump in component 7 and part 4 asserts continuity in components 1
%  through 6; a test demanding seven-state continuity would be asserting a
%  falsehood.
%
%% Note -- the shipped geometry is off-axis, and that is deliberate:
%
%  A due-east trajectory launched from (0,0) is provably BLIND to a
%  transposition at a greatCircle call site: the central angle reduces to
%  acos(cos(lat2) cos(lon2)), which is symmetric in its two arguments, so
%  swapping them leaves the reported range bit-identical. run_boost_glide ships
%  from 20 deg N, 155 deg W on a 60 deg azimuth and impacts at 34.0 deg N,
%  77.6 deg W, so the four coordinates are four genuinely distinct numbers, and
%  part 6 MEASURES that the geometry discriminates rather than assuming it.
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
               veh = vehicle_hgv();
               bst = coorbital.util.boosterDefaults();

%% ---------------------------------------------------------------------
%% 1. The shipped configuration, and that every phase ended as intended
%% ---------------------------------------------------------------------
%% Plots off so the suite creates no figures, and the summary captured rather
%% than printed so the suite output stays readable. The captured text is not
%% discarded -- it is what a user actually reads, so it is parsed and asserted
%% in part 7 alongside the returned struct:
            outDef = evalc(['[trDef,inDef] = run_boost_glide(' ...
                            'struct(''showPlots'',false));']);

    assert(numel(strfind(outDef,'(nominal)')) == 3, ...
        ['expected three nominal phase terminations, found %d. ' ...
         'Summary was:\n%s'],numel(strfind(outDef,'(nominal)')),outDef);
    assert(~contains(outDef,'CAUTION'), ...
        'the shipped run printed the truncation caution. Summary was:\n%s',outDef);
    assert(~contains(outDef,'horizon'), ...
        'a phase hit its integration horizon. Summary was:\n%s',outDef);
    assert(inDef.stopOK,'run_boost_glide reported a non-nominal termination');
    assert(contains(outDef,'propellant exhausted'), ...
        'boost did not end on the burnout event. Summary was:\n%s',outDef);
    assert(contains(outDef,'handoff altitude at'), ...
        'the glide did not end on the handoff event. Summary was:\n%s',outDef);
    assert(contains(outDef,'impact altitude at'), ...
        'the descent did not end on the impact event. Summary was:\n%s',outDef);
    assert(inDef.cfrOK, ...
        'the continued-glide counterfactual did not reach the stop altitude');

%% THE THREE PHASES, IN ORDER, EACH NON-EMPTY. phaseIdx must run 1, 2, 3 and
%% never go backwards; a chain that lost a phase, or ran them out of order, is
%% not the trajectory the rest of this test believes it is asserting:
    assert(isequal(unique(trDef.phaseIdx).',[1 2 3]), ...
        'phases present were %s, expected 1 2 3', ...
        mat2str(unique(trDef.phaseIdx).'));
    assert(all(diff(trDef.phaseIdx) >= 0), ...
        'phaseIdx is not monotonically non-decreasing; the phases are out of order');
    assert(trDef.phaseIdx(1) == 1 && trDef.phaseIdx(end) == 3, ...
        'the chain must begin in phase 1 and end in phase 3');
    for kp = 1:3
        assert(sum(trDef.phaseIdx == kp) >= 2, ...
            'phase %d contributed %d samples; a phase must be non-empty', ...
            kp,sum(trDef.phaseIdx == kp));
    end
    assert(numel(trDef.junction) == 2,'a three-phase chain has two junctions');

%% Seven states, two controls:
    assert(isequal(size(trDef.x),[numel(trDef.t) 7]), ...
        'state history must be [N x 7]; the mass state is component 7');
    assert(isequal(size(trDef.u),[numel(trDef.t) 2]), ...
        'control history must be [N x 2]');
    assert(numel(trDef.t) > 1000, ...
        'only %d samples returned; that is not a resolved trajectory',numel(trDef.t));

%% The three phases really are three DIFFERENT flight regimes -- powered,
%% shallow-banked and hard-banked. Asserted on the recorded control so that a
%% chain which quietly flew one schedule three times cannot pass:
               kBO = find(trDef.phaseIdx == 1,1,'last');
               kHO = find(trDef.phaseIdx == 2,1,'last');
                nS = numel(trDef.t);
    assert(max(abs(trDef.u(trDef.phaseIdx == 2,2))) < 1e-12, ...
        'the glide phase is not flying the zero-bank glide schedule');
    assertAbs(max(abs(trDef.u(trDef.phaseIdx == 3,2))),deg2rad(75),1e-12, ...
        'descent bank angle (rad)');

%% ---------------------------------------------------------------------
%% 2. Time: strictly monotonic across both junctions, and no lost interval
%% ---------------------------------------------------------------------
    assert(all(diff(trDef.t) > 0), ...
        ['the cumulative clock is not strictly increasing; smallest step ' ...
         '%.6e s at sample %d'],min(diff(trDef.t)),find(diff(trDef.t) <= 0,1));

%% phaseRun accumulates each phase's duration onto the clock and records that
%% running total as the junction time. If the two disagree the chain has lost
%% or duplicated an interval, and EXACT equality is the right budget: both are
%% the same floating-point sum, formed once:
    assert(trDef.junction(1).t == trDef.t(kBO), ...
        ['junction 1 is recorded at t = %.17g but the last phase-1 sample is ' ...
         'at %.17g; the cumulative clock has a gap'], ...
        trDef.junction(1).t,trDef.t(kBO));
    assert(trDef.junction(2).t == trDef.t(kHO), ...
        ['junction 2 is recorded at t = %.17g but the last phase-2 sample is ' ...
         'at %.17g; the cumulative clock has a gap'], ...
        trDef.junction(2).t,trDef.t(kHO));

%% ...and the step ACROSS each junction is an ordinary integration step, not a
%% jump. Measured against the largest step the following phase took on its own,
%% so this stays meaningful if the step pattern changes:
    checkJunctionStep(trDef,kBO,2,1);
    checkJunctionStep(trDef,kHO,3,2);

%% ---------------------------------------------------------------------
%% 3. Mass bookkeeping across stage separation
%% ---------------------------------------------------------------------
%% On the pad, EXACTLY the documented stack -- these are decimal literals added
%% together, so anything but exact equality is a bookkeeping change:
    assert(trDef.x(1,7) == veh.mass + bst.massDry + bst.massProp, ...
        'liftoff mass %.17g, expected %.17g', ...
        trDef.x(1,7),veh.mass + bst.massDry + bst.massProp);

%% Mass DECREASES STRICTLY while the motor runs. Not merely "ends lower":
%% a constant-flow motor cannot pause, and a sign error or a clamped thrust
%% would show up here as a flat or rising stretch:
             mBst  = trDef.x(trDef.phaseIdx == 1,7);
    assert(all(diff(mBst) < 0), ...
        ['the boost mass did not decrease strictly; %d of %d steps were flat ' ...
         'or rising'],sum(diff(mBst) >= 0),numel(mBst)-1);

%% At burnout, payload plus spent structure. This is the event target, so it is
%% exact only to the event solver, and it is the value BEFORE the link:
    assertAbs(trDef.x(kBO,7),veh.mass + bst.massDry,1e-6,'burnout mass (kg)');

%% AFTER THE LINK. junction(1).x is the state as handed to the glide, i.e. the
%% far side of the staging jump, so its mass component is the payload alone.
%% This is the assertion that removing the separation link cannot survive:
    assertAbs(trDef.junction(1).x(7),veh.mass,1e-6, ...
        'mass at the separation junction (kg), post-link');

%% ...and the jump itself is exactly the dry mass. Exact, because the link is a
%% subtraction of a literal and touches nothing the integrator produced:
    assert(trDef.x(kBO,7) - trDef.junction(1).x(7) == bst.massDry, ...
        ['the staging jump was %.17g kg against a %.17g kg dry mass; the link ' ...
         'is not jettisoning the booster'], ...
        trDef.x(kBO,7) - trDef.junction(1).x(7),bst.massDry);

%% Components 1 through 6 pass through the link untouched -- the chain is
%% discontinuous in mass at that junction and in NOTHING else:
    assert(isequal(trDef.junction(1).x(1:6),trDef.x(kBO,1:6).'), ...
        'the separation link disturbed a state component other than the mass');

%% The SECOND junction jettisons nothing, so it is continuous in all seven:
    assert(isequal(trDef.junction(2).x,trDef.x(kHO,:).'), ...
        ['the glide-to-descent junction changed the state; nothing is ' ...
         'separated there and it must be an identity handoff']);

%% MASS IS EXACTLY CONSTANT AFTER BURNOUT. dm/dt is identically zero in the
%% unpowered phases, so every Runge-Kutta stage adds exactly zero and the mass
%% must be bit-identical across every sample of BOTH of them, not merely close:
             mUnp  = trDef.x(trDef.phaseIdx >= 2,7);
    assert(max(mUnp) - min(mUnp) == 0, ...
        ['the unpowered mass varied by %.3e kg across %d samples spanning two ' ...
         'phases; with dm/dt = 0 it must be bit-exactly constant'], ...
        max(mUnp) - min(mUnp),numel(mUnp));
    assertAbs(mUnp(1),veh.mass,1e-6,'constant unpowered mass (kg)');
    assertAbs(trDef.junction(2).x(7),veh.mass,1e-6, ...
        'mass at the handoff junction (kg)');

%% ---------------------------------------------------------------------
%% 4. Continuity across BOTH junctions, against an independent reference
%% ---------------------------------------------------------------------
%% See the note in the header for why the reference is ode89 at 1e-12 and not a
%% re-run of the driver. Budgets are absolute and in the natural unit of each
%% state; the measured residuals sit 300 to 3000 times inside them, while any
%% real discontinuity at a seam -- a dropped link, a state handed on unlinked,
%% a phase started from the wrong row -- is metres to kilometres:
              tolX = [1e-3; 1e-9; 1e-9; 1e-3; 1e-8; 1e-9];
    crossJunction(trDef,inDef,1,tolX);
    crossJunction(trDef,inDef,2,tolX);

%% ---------------------------------------------------------------------
%% 5. Headline numbers, pinned to 1e-4 relative
%% ---------------------------------------------------------------------
    assertRel(inDef.tBurnout,80.517757894736846,1e-4,'burnout time (s)');
    assertRel(inDef.hBurnout,46997.974129904062,1e-4,'burnout altitude (m)');
    assertRel(inDef.vBurnout,5840.7450141833751,1e-4,'burnout speed (m/s)');
    assertRel(rad2deg(inDef.gamBurnout),9.7577955576856557,1e-4, ...
        'burnout flight path angle (deg)');
    assertRel(inDef.downBOKm,113.67617200129315,1e-4,'burnout downrange (km)');
    assertAbs(inDef.mGlide,veh.mass,1e-12,'mass into the glide (kg)');

    assertRel(inDef.tHandoff,2120.5260710551343,1e-4,'handoff time (s)');
    assertRel(inDef.vHandoff,694.88963225448549,1e-4,'handoff speed (m/s)');
    assertRel(rad2deg(inDef.gamHandoff),-4.2125785961442999,1e-4, ...
        'handoff flight path angle (deg)');
    assertRel(inDef.hGlideMax,168.27897961791419,1e-4,'peak glide altitude (km)');
    assertRel(inDef.glideKm,7527.5175161512025,1e-4,'glide ground leg (km)');
    assertRel(inDef.descKm,30.13417876516138,1e-4,'descent ground leg (km)');

    assertRel(inDef.tFlight,2194.7735422243827,1e-4,'flight time (s)');
    assertRel(inDef.rangeKm,7663.0518386937174,1e-4,'ground range (km)');
    assertRel(rad2deg(inDef.angTot),68.838365896987682,1e-4,'central angle (deg)');
    assertRel(inDef.vImpact,401.5707782584081,1e-4,'impact speed (m/s)');
    assertRel(rad2deg(inDef.gamImpact),-32.22774392588861,1e-4, ...
        'impact flight path angle (deg)');

    assertRel(inDef.qBstMax,113155.71498409571,1e-4,'boost peak dynamic pressure (Pa)');
    assertRel(inDef.nBstMax,39.886332191981168,1e-4,'boost peak sensed load (g)');
    assertRel(inDef.qGldMax,260188.99657384725,1e-4,'glide peak dynamic pressure (Pa)');
    assertRel(inDef.nGldMax,8.334586218596451,1e-4,'glide peak deceleration (g)');
    assertRel(inDef.qDscMax,98954.597663692897,1e-4,'descent peak dynamic pressure (Pa)');
    assertRel(inDef.nDscMax,3.1697944064305994,1e-4,'descent peak aero load (g)');
    assertRel(inDef.nAerMax,8.334586218596451,1e-4,'peak aero load (g)');
    assertRel(inDef.tNAerMax,569.7723246998296,1e-4,'time of the peak aero load (s)');

%% The peak aero LOAD is not a deceleration, and the split is what proves it.
%% With constLD the ratio is fixed by the drag polar, so the load can be
%% reassembled from the vehicle's own CL and LD -- an identity the summary's
%% two "of which" lines must satisfy, and one that a future edit relabelling
%% the load as a deceleration again would break:
    assertRel(inDef.nAerLift,7.73846903214032,1e-4,'lift part of the peak load (g)');
    assertRel(inDef.nAerDrag,3.0953876128561273,1e-4,'drag part of the peak load (g)');
    assertRel(inDef.nAerLift./inDef.nAerDrag,veh.LD,1e-12, ...
        'lift-to-drag ratio at the peak load');
    assertRel(sqrt(inDef.nAerLift.^2 + inDef.nAerDrag.^2),inDef.nAerMax,1e-12, ...
        'peak load reassembled from its lift and drag parts');
    assert(inDef.nAerDrag < inDef.nAerMax./2, ...
        ['the drag part is %.3f g of a %.3f g load; if these ever come close ' ...
         'the "not a deceleration" caveat has stopped being the point'], ...
        inDef.nAerDrag,inDef.nAerMax);

%% The boost peak sensed load is the end-of-burn thrust over the burnout mass,
%% LESS the drag the stack is still feeling at 47 km -- unlike the ballistic
%% ascent, this one burns out low enough for that term to matter. Checking the
%% peak against that closed form catches a peak found correctly from the wrong
%% acceleration, and the one-sided form catches a sign error on the drag:
          nBstIdeal = bst.thrustVac./((veh.mass + bst.massDry).*c.g0);
    assert(inDef.nBstMax < nBstIdeal, ...
        ['the boost peak sensed load %.6f g exceeds thrust over burnout ' ...
         'weight, %.6f g. Drag can only reduce it.'],inDef.nBstMax,nBstIdeal);
    assertRel(inDef.nBstMax,nBstIdeal,0.03, ...
        'boost peak sensed load against thrust over burnout weight');

%% Terminal and handoff altitudes are the configured stop altitudes.
%% eventAltitude is a terminal one-sided event, so both are exact to the event
%% solver, and the handoff one in particular is what a mis-ordered phase chain
%% cannot reproduce:
    assertAbs(inDef.hHandoff,15000,1e-3,'handoff altitude (m)');
    assertAbs(trDef.x(nS,1) - c.rE,0,1e-3,'impact altitude (m)');

%% ---------------------------------------------------------------------
%% 5b. What the descent phase buys, against the continued-glide counterfactual
%% ---------------------------------------------------------------------
%% The SIZE of each difference is a model observation and is pinned only as a
%% regression. The DIRECTION is fixed by theory: the descent keeps cos(75 deg)
%% of its lift supporting the vehicle against the glide's cos(0), so from the
%% same state it must fall steeper, arrive sooner, arrive faster, and meet more
%% dynamic pressure. Those four are asserted as inequalities:
             tDesc = inDef.tFlight - inDef.tHandoff;
    assert(tDesc < inDef.tCfr, ...
        ['the descent took %.2f s against the continued glide''s %.2f s. With ' ...
         'less lift holding it up it cannot take longer.'],tDesc,inDef.tCfr);
    assert(abs(inDef.gamDscAvg) > abs(inDef.gamCfrAvg), ...
        ['the descent averaged %.3f deg of flight path against the continued ' ...
         'glide''s %.3f deg; it must be the steeper segment'], ...
        inDef.gamDscAvg,inDef.gamCfrAvg);
    assert(inDef.vImpact > inDef.vCfr, ...
        ['the descent arrived at %.2f m/s against the continued glide''s ' ...
         '%.2f m/s'],inDef.vImpact,inDef.vCfr);
    assert(inDef.qDscMax > inDef.qCfrMax, ...
        ['the descent peaked at %.2f kPa against the continued glide''s ' ...
         '%.2f kPa'],inDef.qDscMax./1000,inDef.qCfrMax./1000);
    assertRel(inDef.tCfr,235.3715908462047,1e-4,'continued-glide time to impact (s)');
    assertRel(inDef.cfrKm,96.436684498197963,1e-4,'continued-glide ground leg (km)');
    assertRel(inDef.vCfr,228.88751714419135,1e-4,'continued-glide impact speed (m/s)');
    assertRel(inDef.qCfrMax,36826.257582635859,1e-4, ...
        'continued-glide peak dynamic pressure (Pa)');

%% ---------------------------------------------------------------------
%% 6. The greatCircle argument order, which the shipped geometry does pin
%% ---------------------------------------------------------------------
%% The angle the launch-to-impact call site SHOULD have produced, formed here
%% in the documented (lat1,lon1,lat2,lon2) order:
              lat1 = deg2rad(20);
              lon1 = deg2rad(-155);
              lat2 = trDef.x(nS,3);
              lon2 = trDef.x(nS,2);
            angRef = coorbital.util.greatCircle(lat1,lon1,lat2,lon2);

%% Every way of getting the order wrong, so the assertion that follows is SHOWN
%% to be discriminating rather than assumed to be. If a future change to the
%% launch site ever makes one of these coincide with the truth, this guard
%% fails loudly instead of the order check going quietly blind:
           angSwap = [coorbital.util.greatCircle(lat1,lon1,lon2,lat2); ...
                      coorbital.util.greatCircle(lon1,lat1,lat2,lon2); ...
                      coorbital.util.greatCircle(lon1,lat1,lon2,lat2)];
    assert(min(abs(rad2deg(angSwap - angRef))) > 5, ...
        ['the shipped geometry is not discriminating: a transposed greatCircle ' ...
         'call comes within %.3f deg of the correct %.4f deg'], ...
        min(abs(rad2deg(angSwap - angRef))),rad2deg(angRef));

%% ...and now the check itself, against both the returned angle and the printed
%% one:
    assertAbs(rad2deg(inDef.angTot),rad2deg(angRef),1e-9, ...
        'returned central angle against the terminal state in (lat,lon) order');
            angPrn = summaryNumber(outDef,'central angle +([-\d.]+) +deg');
    assert(abs(angPrn - rad2deg(angRef)) < 1e-3, ...
        ['reported central angle %.4f deg against %.4f deg from the terminal ' ...
         'state in (lat,lon) order. The arguments at the greatCircle call site ' ...
         'in run_boost_glide look transposed.'],angPrn,rad2deg(angRef));

%% Terminal latitude and longitude separately, which is what catches a
%% transposition in the STATE indexing, x(:,2) against x(:,3):
    assertRel(rad2deg(inDef.latImpact),33.985342136064382,1e-4, ...
        'impact latitude (deg)');
    assertRel(rad2deg(inDef.lonImpact),-77.58837898057115,1e-4, ...
        'impact longitude (deg)');
    assert(abs(rad2deg(lat2) - rad2deg(lon2)) > 10, ...
        'impact latitude and longitude are too close together to tell apart');

%% The flight really did leave the launch meridian and the launch parallel, so
%% the four coordinates in the call above are four distinct numbers rather than
%% an accident of the endpoint:
    assert(abs(rad2deg(lat2) - 20) > 10 && abs(rad2deg(lon2) - (-155)) > 10, ...
        'the trajectory did not move far enough in both coordinates');

%% ---------------------------------------------------------------------
%% 7. The summary text reports what the struct returned
%% ---------------------------------------------------------------------
%% A number computed correctly and then printed from the wrong variable is a
%% failure the struct alone cannot see:
    assertAbs(summaryNumber(outDef,'ground range +([-\d.]+) +km'), ...
        inDef.rangeKm,0.005 + 1e-9,'printed ground range (km)');
    assertAbs(summaryNumber(outDef,'flight time +([-\d.]+) +s'), ...
        inDef.tFlight,0.005 + 1e-9,'printed flight time (s)');
    assertAbs(summaryNumber(outDef,'peak aero load +([-\d.]+) +g'), ...
        inDef.nAerMax,0.005 + 1e-9,'printed peak aero load (g)');
    assertAbs(summaryNumber(outDef,'of which drag +([-\d.]+) +g'), ...
        inDef.nAerDrag,0.005 + 1e-9,'printed drag part of the peak load (g)');

%% The label itself is asserted. The quantity is a load factor and calling it a
%% deceleration overstates the braking 2.69-fold, so the word must not come
%% back -- and the honest wording must be present:
    assert(~contains(outDef,'decel  ') && ~contains(outDef,'peak decel'), ...
        ['the summary is labelling the aerodynamic load factor as a ' ...
         'deceleration again. It is %.1f %% lift.'], ...
        100.*inDef.nAerLift./inDef.nAerMax);
    assert(contains(outDef,'lift and drag, thrust excluded'), ...
        'the peak-load lines must say what the load factor is made of');
    assert(contains(outDef,'ONLY THE DRAG PART'), ...
        'the summary must say which part of the load actually brakes');
    assertAbs(summaryNumber(outDef,'handoff speed +([-\d.]+) +m/s'), ...
        inDef.vHandoff,0.005 + 1e-9,'printed handoff speed (m/s)');
    assertAbs(summaryNumber(outDef,'mass into glide +([-\d.]+) +kg'), ...
        veh.mass,0.05 + 1e-9,'printed mass into the glide (kg)');
    assert(contains(outDef,'booster JETTISONED at burnout'), ...
        'the summary did not report that the booster was jettisoned');
    assert(contains(outDef,'MUST fall steeper and arrive sooner'), ...
        'the summary did not reach the expected verdict on the descent bank');
    assert(contains(outDef,'PLACEHOLDER'), ...
        'the summary must carry the PLACEHOLDER caveat');
    assert(contains(outDef,'Validity:'), ...
        'the summary must carry a validity note');

%% ---------------------------------------------------------------------
%% 8. Separation switched OFF, which is a supported configuration
%% ---------------------------------------------------------------------
%% Keeping the dead booster is legitimate to model, and the user block exposes
%% it. Flying it here proves the link = [] path works AND, by differing from the
%% shipped run, proves the shipped run really is jettisoning something. It is
%% also the honest way to price a retained booster: with the separation link
%% simply deleted the chain does not fly at all, because massConstant refuses a
%% carried mass that disagrees with the vehicle it is flying:
            outNos = evalc(['[trNos,inNos] = run_boost_glide(struct(' ...
                            '''showPlots'',false,''separation'',false));']);
    assert(inNos.stopOK && ~contains(outNos,'CAUTION'), ...
        'the no-separation run did not terminate nominally. Summary was:\n%s',outNos);
    assert(contains(outNos,'booster RETAINED'), ...
        'the no-separation run did not report a retained booster');
    assertAbs(inNos.mGlide,veh.mass + bst.massDry,1e-12,'retained glide mass (kg)');
    assertAbs(trNos.junction(1).x(7),veh.mass + bst.massDry,1e-6, ...
        'retained mass at the separation junction (kg)');
    assertRel(inNos.rangeKm,2853.708329419961,1e-4,'no-separation ground range (km)');
    assertRel(inNos.tFlight,638.64034384572642,1e-4,'no-separation flight time (s)');
    assertRel(inNos.vImpact,408.63318460877383,1e-4,'no-separation impact speed (m/s)');

%% ...and the two runs must be materially different flights. Dragging a dead
%% booster through the glide costs range, and it costs an enormous amount of it
%% here: the stack is 2.7 times the mass but a far worse lifting shape, so it
%% does not glide at all and arrives 4800 km short:
    assert(inDef.rangeKm - inNos.rangeKm > 1000, ...
        ['keeping the booster changed the range by only %.1f km; the two ' ...
         'configurations are not distinguishable, so nothing here would ' ...
         'notice a dropped separation link'],inDef.rangeKm - inNos.rangeKm);

%% ---------------------------------------------------------------------
%% 9. The phugoid-trough warning on hHandoff
%% ---------------------------------------------------------------------
%% hHandoff is advertised as freely settable, and it is the one entry in the
%% user block that can change the answer by a quarter while every phase still
%% reports NOMINAL: the glide descends in a damped skip phugoid, so a handoff
%% above a trough ends the glide a skip early. The script detects it from the
%% counterfactual it already flies, and this part asserts that the detector is
%% quiet on the shipped case and loud on a case that really is truncated:
    assert(~inDef.troughWarn && ~contains(outDef,'PHUGOID TROUGH'), ...
        ['the shipped 15 km handoff raised the trough warning; rebound was ' ...
         '%.1f m'],inDef.reboundM);
    assertAbs(inDef.reboundM,0,1e-6,'shipped continued-glide rebound (m)');

            out30 = evalc(['[~,in30] = run_boost_glide(struct(' ...
                           '''showPlots'',false,''hHandoff'',30));']);
    assert(in30.stopOK && ~contains(out30,'CAUTION'), ...
        ['the 30 km handoff case must still terminate nominally -- that it ' ...
         'does is the whole hazard']);
    assert(in30.troughWarn && contains(out30,'PHUGOID TROUGH'), ...
        ['a 30 km handoff cuts the glide short by a whole skip and the ' ...
         'warning did not fire; rebound was %.1f m'],in30.reboundM);
    assertRel(in30.reboundM,31325.514742694795,1e-4,'30 km continued-glide rebound (m)');
    assertRel(in30.rangeKm,5781.5704976898978,1e-4,'30 km handoff ground range (km)');
    assert(inDef.rangeKm - in30.rangeKm > 1500, ...
        ['the 30 km handoff cost only %.1f km; if that ever falls the warning ' ...
         'is warning about nothing'],inDef.rangeKm - in30.rangeKm);

%% The blind spot, asserted so it stays documented rather than discovered. The
%% check is one-sided: a handoff that truncates a shallow skip WITHOUT the
%% continued glide climbing back out of it is invisible to it. 25 km is exactly
%% that case -- it costs 173 km of range and raises no warning:
            out25 = evalc(['[~,in25] = run_boost_glide(struct(' ...
                           '''showPlots'',false,''hHandoff'',25));']);
    assert(~in25.troughWarn && ~contains(out25,'PHUGOID TROUGH'), ...
        ['the 25 km case is the documented blind spot of the rebound check. ' ...
         'If it now fires, the check got better and this assertion, the ' ...
         'header note and the user-block comment all need rewriting.']);
    assert(inDef.rangeKm - in25.rangeKm > 100, ...
        ['the 25 km blind spot cost only %.1f km of range; the caveat may no ' ...
         'longer be worth its words'],inDef.rangeKm - in25.rangeKm);

%% ---------------------------------------------------------------------
%% 10. The unpowered guidance runs on the PHASE-LOCAL clock
%% ---------------------------------------------------------------------
%% phaseRun evaluates each phase's guide at that phase's own tspan-relative
%% time, so a glide schedule written against [0 .. tMaxGlide] starts at the
%% handoff from the BOOST, not at liftoff. Both shipped unpowered schedules are
%% constant, which makes the whole convention invisible: an 80 s offset in the
%% clock handed to the guide would change nothing at all and pass every
%% assertion above. A time-varying bank makes it visible:
            gTime = [0 400 800 6000];
            gBank = [0  70   0    0];
            outTv = evalc(['[~,inTv] = run_boost_glide(struct(''showPlots'',false,' ...
                           '''glideTime'',gTime,''glideBank'',gBank,' ...
                           '''glideAlpha'',[0 0 0 0]));']);
    assert(inTv.stopOK && ~contains(outTv,'CAUTION') && ~contains(outTv,'PHUGOID'), ...
        'the time-varying-bank case did not fly cleanly. Summary was:\n%s',outTv);
    assertRel(inTv.rangeKm,6101.6034301223772,1e-4, ...
        'time-varying glide bank, ground range (km)');
    assertRel(inTv.tFlight,1951.9689001876234,1e-4, ...
        'time-varying glide bank, flight time (s)');

%% ...and the discrimination is MEASURED, not assumed. The same schedule shifted
%% by the boost duration is what the flight would see if the guide were being
%% handed cumulative time instead of phase-local time. It must be a different
%% flight, and it is a wildly different one -- 1839 km, 30 percent:
            gShft = gTime + inDef.tBurnout;
            outSh = evalc(['[~,inSh] = run_boost_glide(struct(''showPlots'',false,' ...
                           '''glideTime'',gShft,''glideBank'',gBank,' ...
                           '''glideAlpha'',[0 0 0 0]));']);
    assert(abs(inTv.rangeKm - inSh.rangeKm) > 100, ...
        ['shifting the glide schedule by the %.1f s boost duration moved the ' ...
         'range only %.3f km. This case exists to pin the phase-local clock ' ...
         'convention, and at that sensitivity it does not pin it.'], ...
        inDef.tBurnout,abs(inTv.rangeKm - inSh.rangeKm));

%% ---------------------------------------------------------------------
%% 11. The override hook itself
%% ---------------------------------------------------------------------
%% A misspelt parameter must raise rather than silently leave the shipped value
%% in place, or a future test could believe it flew a case it did not:
             threw = false;
    try
        run_boost_glide(struct('hHandoffKm',15));
    catch
             threw = true;
    end
    assert(threw,'an unrecognised override field must raise, not be ignored');
end

function checkJunctionStep(traj,kJ,kNext,kJun)
%% Purpose:
%
%  Confirm that the step across a phase junction is an ordinary integration
%  step and not a lost interval. A chain that dropped time at a seam -- by
%  referencing a phase to the wrong origin, or by discarding samples -- shows
%  up as a single step at the junction far larger than any step the following
%  phase takes on its own.
%
%% Inputs:
%
%  traj             Struct                      Trajectory from phaseRun
%
%  kJ               [1 x 1]                     Index of the LAST sample of the
%                                               phase that is ending
%
%  kNext            [1 x 1]                     Index of the phase that follows
%
%  kJun             [1 x 1]                     Which junction this is, for the
%                                               failure message
%
%% Outputs:
%
%  none                                         Throws on failure
%
%% Revision History:
%  Michael Casey                                                08/07/2026
%  Copyright 2026 Coorbital, Inc.
%% ------------------------ Begin Code Sequence ---------------------------

               dtJ = traj.t(kJ+1) - traj.t(kJ);
            dtNext = diff(traj.t(traj.phaseIdx == kNext));
    assert(dtJ > 0,'junction %d did not advance the clock',kJun);
    assert(dtJ <= max(dtNext), ...
        ['the step across junction %d is %.6g s against a largest in-phase ' ...
         'step of %.6g s in phase %d; that is a gap, not a step'], ...
        kJun,dtJ,max(dtNext),kNext);
end

function crossJunction(traj,info,kJun,tolX)
%% Purpose:
%
%  Re-integrate ACROSS one phase junction with an independent solver and
%  confirm that the driver's trajectory is continuous through it in states 1
%  through 6. The reference starts at a driver sample two seconds before the
%  junction, runs the OUTGOING phase's dynamics up to the junction time,
%  applies that phase's link, runs the INCOMING phase's dynamics on to a driver
%  sample two seconds after it, and compares.
%
%  ode89 at 1e-12 against the driver's ode45 at 1e-10: a different method at a
%  tighter tolerance. The tolerance ordering is what the check rests on and is
%  asserted below; see the note in the file header for why the more familiar
%  "a same-solver re-run goes bitwise-zero" argument does not apply to this
%  anchoring scheme and must not be relied on.
%
%  Mass is checked too, but as a JUMP, not as a continuity: at the separation
%  junction component 7 must drop by the booster dry mass, and demanding
%  continuity there would be asserting a falsehood.
%
%% Inputs:
%
%  traj             Struct                      Trajectory from phaseRun
%
%  info             Struct                      run_boost_glide summary; uses
%                                               the phases and env fields
%
%  kJun             [1 x 1]                     Junction index, 1 or 2; the
%                                               outgoing phase is kJun
%
%  tolX             [6 x 1]                     Absolute budgets on states 1
%                                               through 6, in m, rad, rad, m/s,
%                                               rad, rad
%
%% Outputs:
%
%  none                                         Throws on failure
%
%% Revision History:
%  Michael Casey                                                08/07/2026
%  Copyright 2026 Coorbital, Inc.
%% ------------------------ Begin Code Sequence ---------------------------

                ph = info.phases;
               env = info.env;
                kJ = find(traj.phaseIdx == kJun,1,'last');
                tJ = traj.t(kJ);

%% Phase origins on the cumulative clock. Each phase's guide is written against
%% its OWN tspan, which starts at zero, so the reference must subtract the
%% phase's start time before calling it:
              tOff = [0; traj.junction(1).t; traj.junction(2).t];

%% Anchor on ACTUAL driver samples either side, so nothing is interpolated:
                ka = find(traj.t <= tJ - 2 & traj.phaseIdx == kJun,1,'last');
                kb = find(traj.t >= tJ + 2 & traj.phaseIdx == kJun+1,1,'first');
    assert(~isempty(ka) && ~isempty(kb), ...
        ['junction %d has no driver sample two seconds either side; the ' ...
         'cross-junction check has nothing to anchor on'],kJun);

                o1 = tOff(kJun);
                o2 = tOff(kJun+1);
              odeA = @(t,x) ph(kJun).eom(t-o1,x,ph(kJun).guide(t-o1,x),[],env);
              odeB = @(t,x) ph(kJun+1).eom(t-o2,x,ph(kJun+1).guide(t-o2,x),[],env);

%% The reference must be TIGHTER than the driver, or the residual is measuring
%% its own error as much as the driver's and bounds nothing. phaseRun defaults
%% to 1e-10 and lets env override it, so the driver's figure is read rather
%% than assumed, and the ordering is asserted rather than trusted to a comment:
            refTol = 1e-12;
            drvTol = 1e-10;
    if isfield(env,'odeRelTol')
            drvTol = env.odeRelTol;
    end
    assert(refTol <= drvTol./100, ...
        ['the cross-junction reference runs at %.1e against the driver''s ' ...
         '%.1e. A reference that is not at least a hundred times tighter ' ...
         'cannot bound the driver''s error, and this check would be measuring ' ...
         'nothing.'],refTol,drvTol);
              opts = odeset('RelTol',refTol,'AbsTol',refTol);

%% Leg one, up to the junction:
              solA = ode89(odeA,[traj.t(ka) tJ],traj.x(ka,:).',opts);
               xJm = solA.y(:,end);
    for k = 1:6
        assert(abs(xJm(k) - traj.x(kJ,k)) < tolX(k), ...
            ['junction %d, state %d on the NEAR side: ode89 says %.15g, the ' ...
             'driver says %.15g, %.3e apart against a %.1e budget'], ...
            kJun,k,xJm(k),traj.x(kJ,k),abs(xJm(k) - traj.x(kJ,k)),tolX(k));
    end
    assertAbs(xJm(7),traj.x(kJ,7),1e-6, ...
        sprintf('junction %d near-side mass (kg)',kJun));

%% The link, which is the whole point of the junction when there is one:
               xJp = xJm;
    if ~isempty(ph(kJun).link)
               xJp = ph(kJun).link(xJm);
    end
    assertAbs(xJp(7),traj.junction(kJun).x(7),1e-6, ...
        sprintf('junction %d post-link mass (kg)',kJun));

%% Leg two, on into the next phase:
              solB = ode89(odeB,[tJ traj.t(kb)],xJp,opts);
                xb = solB.y(:,end);
    for k = 1:6
        assert(abs(xb(k) - traj.x(kb,k)) < tolX(k), ...
            ['junction %d, state %d two seconds PAST the seam: ode89 says ' ...
             '%.15g, the driver says %.15g, %.3e apart against a %.1e budget. ' ...
             'The chain is not continuous through that junction.'], ...
            kJun,k,xb(k),traj.x(kb,k),abs(xb(k) - traj.x(kb,k)),tolX(k));
    end
    assertAbs(xb(7),traj.x(kb,7),1e-6, ...
        sprintf('junction %d far-side mass (kg)',kJun));
end

function val = summaryNumber(txt,pat)
%% Purpose:
%
%  Pull one number out of a captured run_boost_glide summary. Parsing the
%  printed text rather than recomputing the quantity is deliberate: what a
%  reader of this library acts on is the summary, so the summary is what gets
%  asserted.
%
%  The patterns use explicit space runs rather than \s, because in MATLAB \s
%  also matches a newline and would happily reach across into the next line of
%  the report.
%
%% Inputs:
%
%  txt              Char [1 x n]                Captured run_boost_glide summary
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
