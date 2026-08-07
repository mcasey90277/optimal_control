function test_runBallistic()
%% Purpose:
%
%  Pin the COMPOSITION of the three-phase ballistic chain. Every other test in
%  this suite exercises one library function in isolation; this one exercises
%  the entry script that wires them together -- unit conversion, the two-vehicle
%  binding, the mass bookkeeping across stage separation, the boost/apogee/
%  impact event chain, the derived-quantity block, the great-circle call sites,
%  the Keplerian cross-check, the termination diagnosis and the printed
%  summary. Without it the suite passes unchanged while run_ballistic drops
%  the separation link, flies the boost on the wrong airframe, transposes
%  latitude and longitude, or reports the wrong range.
%
%% Note -- the shipped geometry is ALREADY off-axis, and that is deliberate:
%
%  A due-east trajectory launched from (0,0) is provably BLIND to a
%  transposition at a greatCircle call site: the central angle reduces to
%  acos(cos(lat2) cos(lon2)), which is symmetric in its two arguments, so
%  swapping them leaves the reported range bit-identical. Asserting the
%  terminal latitude and longitude separately does not close the hole either,
%  because those come from traj, which such a mutation never touches. run_glide
%  ships due east and has to fly a SECOND, off-axis case to catch it.
%  run_ballistic ships from 45 deg N, 100 deg W on a 35 deg azimuth precisely
%  so that no second run is needed: the four coordinates are four genuinely
%  distinct numbers, and part 4 below MEASURES that the geometry discriminates
%  rather than assuming it -- the nearest transposed reading is 25.8 deg away
%  from the true 40.75 deg central angle.
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
               veh = vehicle_bm();
               bst = coorbital.util.boosterDefaults();

%% ---------------------------------------------------------------------
%% 1. vehicle_bm is a re-entry body, not a second copy of the glide vehicle
%% ---------------------------------------------------------------------
%% The whole point of the file is that it differs from the library defaults.
%% If a future tidy-up ever collapses it back onto them, the "ballistic"
%% trajectory silently becomes a glide and the Keplerian cross-check below
%% starts comparing two different flight regimes:
               vdf = coorbital.util.vehicleDefaults();
    assert(veh.Sref ~= vdf.Sref && veh.CL ~= vdf.CL && veh.LD ~= vdf.LD, ...
        ['vehicle_bm has collapsed onto vehicleDefaults. It exists BECAUSE ' ...
         'it must differ: a waverider flown ballistic glides, and the ' ...
         'trajectory stops being the Keplerian arc this test checks.']);

%% The mass, by contrast, MUST match, because the mass bookkeeping quoted in
%% coorbital.util.boosterDefaults -- 32400 kg on the pad, 2400 kg at burnout,
%% 900 kg after separation -- is written against that payload:
    assert(veh.mass == vdf.mass, ...
        ['vehicle_bm.mass = %.17g against vehicleDefaults.mass = %.17g. The ' ...
         'mass bookkeeping in boosterDefaults is quoted against the default ' ...
         'payload; changing it here invalidates every number in that header.'], ...
        veh.mass,vdf.mass);
    assert(veh.noseRadius > 0,'vehicle_bm must carry a positive noseRadius');

%% CD is derived, not stored, so the pair below is what fixes the drag. Pinned
%% because the ballistic coefficient is what makes the trajectory ballistic:
             CDref = veh.CL./veh.LD;
    assertRel(CDref,0.25,1e-12,'vehicle_bm derived CD');
    assertRel(veh.mass./(CDref.*veh.Sref),9350.649350649351,1e-9, ...
        'vehicle_bm ballistic coefficient (kg/m^2)');

%% ...and CL is pinned ABSOLUTELY as well, which the two joint pins above do
%% NOT cover. They constrain only the RATIO: doubling CL to 0.010 and LD to
%% 0.040 together holds CD at 0.25 and the ballistic coefficient at
%% 9350.649350649351 exactly, so both assertions above stay green. The
%% trajectory does not. coorbital.eom.glide3DOF reads CL directly for lift --
%% aLift = qbar*Sref*CL/mass, alongside the aDrag term that reads CD -- so the
%% split of a fixed CD between CL and LD is OBSERVABLE, through lift and not
%% through drag. Measured on that mutation: impact speed moves 3156.49 to
%% 3144.86 m/s and the pinned literal further down this file breaks at
%% 3.683e-03 with nothing saying why. This pin makes the cause fire first.
%% With CL and the derived CD both pinned, LD = CL/CD is fixed too:
    assertRel(veh.CL,0.005,1e-12,'vehicle_bm CL');

%% ---------------------------------------------------------------------
%% 2. The shipped configuration
%% ---------------------------------------------------------------------
%% Plots off so the suite creates no figures, and the summary captured rather
%% than printed so the suite output stays readable. The captured text is not
%% discarded -- it is what a user actually reads, so it is parsed and asserted
%% alongside the returned struct:
            outDef = evalc(['[trDef,inDef] = run_ballistic(' ...
                            'struct(''showPlots'',false));']);

%% Every phase terminated as intended. All three halves matter: three nominal
%% phase lines must be present AND the truncation banner absent AND the
%% script's own verdict true, so that a future change to the diagnosis cannot
%% pass by printing contradictory things:
    assert(numel(strfind(outDef,'(nominal)')) == 3, ...
        ['expected three nominal phase terminations, found %d. ' ...
         'Summary was:\n%s'],numel(strfind(outDef,'(nominal)')),outDef);
    assert(~contains(outDef,'CAUTION'), ...
        'the shipped run printed the truncation caution. Summary was:\n%s',outDef);
    assert(~contains(outDef,'horizon'), ...
        'a phase hit its integration horizon. Summary was:\n%s',outDef);
    assert(inDef.stopOK,'run_ballistic reported a non-nominal termination');
    assert(contains(outDef,'propellant exhausted'), ...
        'boost did not end on the burnout event. Summary was:\n%s',outDef);
    assert(contains(outDef,'apogee (nominal)'), ...
        'the coast did not end on the apogee event. Summary was:\n%s',outDef);

%% THE THREE PHASES, IN ORDER. phaseIdx must run 1, 2, 3 and never go
%% backwards; a chain that lost a phase, or ran them out of order, is not the
%% trajectory the rest of this test believes it is asserting:
    assert(isequal(unique(trDef.phaseIdx).',[1 2 3]), ...
        'phases present were %s, expected 1 2 3', ...
        mat2str(unique(trDef.phaseIdx).'));
    assert(all(diff(trDef.phaseIdx) >= 0), ...
        'phaseIdx is not monotonically non-decreasing; the phases are out of order');
    assert(trDef.phaseIdx(1) == 1 && trDef.phaseIdx(end) == 3, ...
        'the chain must begin in phase 1 and end in phase 3');
    assert(numel(trDef.junction) == 2,'a three-phase chain has two junctions');

%% Seven states, two controls:
    assert(isequal(size(trDef.x),[numel(trDef.t) 7]), ...
        'state history must be [N x 7]; the mass state is component 7');
    assert(isequal(size(trDef.u),[numel(trDef.t) 2]), ...
        'control history must be [N x 2]');
    assert(numel(trDef.t) > 500, ...
        'only %d samples returned; that is not a resolved trajectory',numel(trDef.t));

%% ---------------------------------------------------------------------
%% 3. Mass bookkeeping across stage separation
%% ---------------------------------------------------------------------
%% On the pad, EXACTLY the documented stack -- these are decimal literals
%% added together, so anything but exact equality is a bookkeeping change:
    assert(trDef.x(1,7) == veh.mass + bst.massDry + bst.massProp, ...
        'liftoff mass %.17g, expected %.17g', ...
        trDef.x(1,7),veh.mass + bst.massDry + bst.massProp);

%% At burnout, payload plus spent structure. This is the event target, so it
%% is exact only to the event solver:
               kBO = find(trDef.phaseIdx == 1,1,'last');
    assertAbs(trDef.x(kBO,7),veh.mass + bst.massDry,1e-6,'burnout mass (kg)');

%% AFTER THE LINK. junction(1).x is the state as handed to the coast, i.e. the
%% far side of the staging jump, so its mass component is the payload alone.
%% This is the assertion that removing the separation link cannot survive:
    assertAbs(trDef.junction(1).x(7),veh.mass,1e-6, ...
        'mass at the separation junction (kg), post-link');

%% ...and the jump itself is exactly the dry mass. Exact, because the link is
%% a subtraction of a literal and touches nothing the integrator produced:
    assert(trDef.x(kBO,7) - trDef.junction(1).x(7) == bst.massDry, ...
        ['the staging jump was %.17g kg against a %.17g kg dry mass; the ' ...
         'link is not jettisoning the booster'], ...
        trDef.x(kBO,7) - trDef.junction(1).x(7),bst.massDry);

%% Components 1 through 6 pass through the link untouched -- the chain is
%% discontinuous in mass at that junction and in NOTHING else:
    assert(isequal(trDef.junction(1).x(1:6),trDef.x(kBO,1:6).'), ...
        'the separation link disturbed a state component other than the mass');

%% MASS IS EXACTLY CONSTANT AFTER BURNOUT. dm/dt is identically zero in the
%% unpowered phases, so every Runge-Kutta stage adds exactly zero and the mass
%% must be bit-identical across every sample, not merely close:
             mUnp  = trDef.x(trDef.phaseIdx >= 2,7);
    assert(max(mUnp) - min(mUnp) == 0, ...
        ['the unpowered mass varied by %.3e kg across %d samples; with ' ...
         'dm/dt = 0 it must be bit-exactly constant'], ...
        max(mUnp) - min(mUnp),numel(mUnp));
    assertAbs(mUnp(1),veh.mass,1e-6,'constant coast mass (kg)');
    assertAbs(trDef.junction(2).x(7),veh.mass,1e-6, ...
        'mass at the apogee junction (kg)');

%% ---------------------------------------------------------------------
%% 4. Headline numbers, pinned to 1e-4 relative
%% ---------------------------------------------------------------------
    assertRel(inDef.tBurnout,80.517757894737741,1e-4,'burnout time (s)');
    assertRel(inDef.hBurnout,102000.27979884855,1e-4,'burnout altitude (m)');
    assertRel(inDef.vBurnout,5768.1977161798968,1e-4,'burnout speed (m/s)');
    assertRel(rad2deg(inDef.gamBurnout),47.599708657243582,1e-4, ...
        'burnout flight path angle (deg)');
    assertRel(inDef.downBOKm,68.640084819172984,1e-4,'burnout downrange (km)');

    assertRel(inDef.tApogee,838.20774215245854,1e-4,'apogee time (s)');
    assertRel(inDef.hApogee,1619234.597011527,1e-4,'apogee altitude (m)');
    assertRel(inDef.vApogee,3151.6209945263226,1e-4,'apogee speed (m/s)');

    assertRel(inDef.tFlight,1620.6138710598609,1e-4,'flight time (s)');
    assertRel(inDef.rangeKm,4536.361953009874,1e-4,'ground range (km)');
    assertRel(rad2deg(inDef.angTot),40.750832766870516,1e-4,'central angle (deg)');

    assertRel(inDef.vImpact,3156.4873607478771,1e-4,'impact speed (m/s)');
    assertRel(rad2deg(inDef.gamImpact),-47.724959925777831,1e-4, ...
        'impact flight path angle (deg)');

%% Peak loads. The re-entry peak is the "peak deceleration" of a ballistic
%% missile and is the sharpest quantity here, so it is worth saying that it is
%% NOT fragile: the samples either side of it differ from it by 5e-5 relative,
%% five times inside the 1e-4 budget, so a modest change in where ode45 places
%% its steps cannot move it out of tolerance:
    assertRel(inDef.qBstMax,84127.762151770759,1e-4,'boost peak dynamic pressure (Pa)');
    assertRel(inDef.nBstMax,40.363543279094941,1e-4,'boost peak sensed load (g)');
    assertRel(inDef.qReMax,6281459.0322939027,1e-4,'re-entry peak dynamic pressure (Pa)');
    assertRel(inDef.nReMax,68.514883931424833,1e-4,'re-entry peak deceleration (g)');

%% The boost peak sensed load is the end-of-burn thrust over the burnout mass,
%% because the placeholder motor does not throttle. Checking it against that
%% closed form catches a peak found correctly from the wrong acceleration:
    assertRel(inDef.nBstMax, ...
        bst.thrustVac./((veh.mass + bst.massDry).*c.g0),2e-3, ...
        'boost peak sensed load against thrust over burnout weight');

%% Terminal altitude is the configured stop altitude. eventAltitude is a
%% terminal one-sided event, so this is exact to the event solver:
    assertAbs(trDef.x(end,1) - c.rE,0,1e-3,'impact altitude (m)');

%% ---------------------------------------------------------------------
%% 5. The Keplerian cross-check, recomputed here at 1e-6
%% ---------------------------------------------------------------------
%% run_ballistic asserts this internally. It is recomputed here from the
%% burnout state so that the check is stated independently of the file under
%% test: if the formula in run_ballistic were quietly edited to agree with a
%% broken propagation, this copy would still disagree with it.
                xBO = trDef.junction(1).x;
                 r0 = xBO(1);
                 V0 = xBO(4);
               gam0 = xBO(5);
                 rI = c.rE;
                 hK = r0.*V0.*cos(gam0);
                 pK = hK.^2./c.muE;
               epsK = V0.^2./2 - c.muE./r0;
                 eK = sqrt(1 + 2.*epsK.*hK.^2./c.muE.^2);
                 f0 = atan2( hK.*V0.*sin(gam0)./c.muE , pK./r0 - 1 );
              cosfI = (pK./rI - 1)./eK;
    assert(abs(cosfI) <= 1, ...
        ['the Keplerian arc through the burnout state does not reach the ' ...
         'impact radius: (p/rI - 1)/e = %.6f'],cosfI);
                 fI = 2.*pi - acos(cosfI);
               sVac = rI.*mod(fI - f0,2.*pi);

%% Both sides are exact descriptions of one two-body arc, so the budget is
%% tight on purpose. A loose tolerance here would hide a defect rather than
%% absorb a modelling difference:
    assertRel(inDef.sVacProp,sVac,1e-6, ...
        'propagated vacuum range against the closed-form Keplerian range');
    assertRel(inDef.sVacClosed,sVac,1e-12, ...
        'the closed form as run_ballistic computes it');
    assert(inDef.relVac < 1e-6, ...
        'run_ballistic reported %.3e relative disagreement, budget 1e-6',inDef.relVac);
    assertRel(sVac,4467597.5075559439,1e-6,'closed-form vacuum range (m)');

%% The drag and lift range effects are REGRESSION pins on the two extra
%% propagations and their great-circle call sites, not physics tolerances --
%% nothing predicts their size, which is exactly why run_ballistic prints them
%% rather than asserting them. What theory DOES fix is the sign of each, and
%% that is asserted:
    assert(inDef.dragDiff < 0, ...
        ['drag lengthened the arc by %.3f m. Drag removes energy; it cannot ' ...
         'do that.'],inDef.dragDiff);
    assert(inDef.liftDiff > 0, ...
        ['upward lift shortened the arc by %.3f m. With CL > 0 at zero bank ' ...
         'the lift vector is up and must extend the arc.'],-inDef.liftDiff);
    assertRel(inDef.sDrgB,4467557.9029780682,1e-4,'drag-only range (m)');
    assertRel(inDef.sAeroB,4467721.8681907002,1e-4,'full-aerodynamics range (m)');
    assertRel(inDef.vImpVac,5936.2865265314022,1e-4,'vacuum impact speed (m/s)');

%% ---------------------------------------------------------------------
%% 6. The greatCircle argument order, which the shipped geometry does pin
%% ---------------------------------------------------------------------
%% The angle the launch-to-impact call site SHOULD have produced, formed here
%% in the documented (lat1,lon1,lat2,lon2) order:
               lat1 = deg2rad(45);
               lon1 = deg2rad(-100);
               lat2 = trDef.x(end,3);
               lon2 = trDef.x(end,2);
             angRef = coorbital.util.greatCircle(lat1,lon1,lat2,lon2);

%% Every way of getting the order wrong, so the assertion that follows is
%% SHOWN to be discriminating rather than assumed to be. If a future change to
%% the launch site ever makes one of these coincide with the truth, this guard
%% fails loudly instead of the order check going quietly blind:
            angSwap = [coorbital.util.greatCircle(lat1,lon1,lon2,lat2); ...
                       coorbital.util.greatCircle(lon1,lat1,lat2,lon2); ...
                       coorbital.util.greatCircle(lon1,lat1,lon2,lat2)];
    assert(min(abs(rad2deg(angSwap - angRef))) > 5, ...
        ['the shipped geometry is not discriminating: a transposed ' ...
         'greatCircle call comes within %.3f deg of the correct %.4f deg'], ...
        min(abs(rad2deg(angSwap - angRef))),rad2deg(angRef));

%% ...and now the check itself, against both the returned angle and the
%% printed one:
    assertAbs(rad2deg(inDef.angTot),rad2deg(angRef),1e-9, ...
        'returned central angle against the terminal state in (lat,lon) order');
             angPrn = summaryNumber(outDef,'central angle +([-\d.]+) +deg');
    assert(abs(angPrn - rad2deg(angRef)) < 1e-3, ...
        ['reported central angle %.4f deg against %.4f deg from the terminal ' ...
         'state in (lat,lon) order. The arguments at the greatCircle call ' ...
         'site in run_ballistic look transposed.'],angPrn,rad2deg(angRef));

%% Terminal latitude and longitude separately, which is what catches a
%% transposition in the STATE indexing, x(:,2) against x(:,3):
    assertRel(rad2deg(inDef.latImpact),66.032455598584903,1e-4, ...
        'impact latitude (deg)');
    assertRel(rad2deg(inDef.lonImpact),-32.823459847022562,1e-4, ...
        'impact longitude (deg)');
    assert(abs(rad2deg(lat2) - rad2deg(lon2)) > 10, ...
        'impact latitude and longitude are too close together to tell apart');

%% The flight really did leave the launch meridian and the launch parallel, so
%% the four coordinates in the call above are four distinct numbers rather
%% than an accident of the endpoint:
    assert(abs(rad2deg(lat2) - 45) > 10 && abs(rad2deg(lon2) - (-100)) > 10, ...
        'the trajectory did not move far enough in both coordinates');

%% ---------------------------------------------------------------------
%% 7. The summary text reports what the struct returned
%% ---------------------------------------------------------------------
%% A number computed correctly and then printed from the wrong variable is a
%% failure the struct alone cannot see:
    assertAbs(summaryNumber(outDef,'ground range +([-\d.]+) +km'), ...
        inDef.rangeKm,0.005 + 1e-9,'printed ground range (km)');
    assertAbs(summaryNumber(outDef,'apogee +([-\d.]+) +km'), ...
        inDef.hApogee./1000,0.005 + 1e-9,'printed apogee (km)');
    assertAbs(summaryNumber(outDef,'flight time +([-\d.]+) +s'), ...
        inDef.tFlight,0.005 + 1e-9,'printed flight time (s)');
    assertAbs(summaryNumber(outDef,'re-entry decel +([-\d.]+) +g'), ...
        inDef.nReMax,0.005 + 1e-9,'printed re-entry deceleration (g)');
    assertAbs(summaryNumber(outDef,'mass into coast +([-\d.]+) +kg'), ...
        veh.mass,0.05 + 1e-9,'printed coast mass (kg)');
    assert(contains(outDef,'booster JETTISONED at burnout'), ...
        'the summary did not report that the booster was jettisoned');
    assert(contains(outDef,'PLACEHOLDER'), ...
        'the summary must carry the PLACEHOLDER caveat');
    assert(contains(outDef,'Validity:'), ...
        'the summary must carry a validity note');

%% ---------------------------------------------------------------------
%% 8. Separation switched OFF, which is a supported configuration
%% ---------------------------------------------------------------------
%% Keeping the dead booster is legitimate to model, and the user block exposes
%% it. Flying it here proves the link = [] path works AND, by differing from
%% the shipped run, proves the shipped run really is jettisoning something:
            outNos = evalc(['[trNos,inNos] = run_ballistic(struct(' ...
                            '''showPlots'',false,''separation'',false));']);
    assert(inNos.stopOK && ~contains(outNos,'CAUTION'), ...
        'the no-separation run did not terminate nominally. Summary was:\n%s',outNos);
    assert(contains(outNos,'booster RETAINED'), ...
        'the no-separation run did not report a retained booster');
    assertAbs(inNos.mCoast,veh.mass + bst.massDry,1e-12,'retained coast mass (kg)');
    assertAbs(trNos.junction(1).x(7),veh.mass + bst.massDry,1e-6, ...
        'retained mass at the separation junction (kg)');
    assertRel(inNos.rangeKm,4539.5318036615363,1e-4,'no-separation ground range (km)');
    assertRel(inNos.tFlight,1622.0127303584172,1e-4,'no-separation flight time (s)');
    assertRel(inNos.vImpact,2226.7226722199403,1e-4,'no-separation impact speed (m/s)');

%% ...and the two runs must be materially different flights. This is the
%% statement that a missing separation link cannot satisfy: the heavier,
%% draggier stack arrives 930 m/s slower:
    assert(abs(inNos.vImpact - inDef.vImpact) > 500, ...
        ['keeping the booster changed the impact speed by only %.1f m/s; ' ...
         'the two configurations are not distinguishable, so nothing here ' ...
         'would notice a dropped separation link'], ...
        abs(inNos.vImpact - inDef.vImpact));
    assert(abs(inNos.rangeKm - inDef.rangeKm) > 1, ...
        'keeping the booster changed the range by only %.3f km', ...
        abs(inNos.rangeKm - inDef.rangeKm));

%% The Keplerian cross-check must hold in that configuration too -- it is a
%% property of the two-body arc, not of the mass on it:
    assert(inNos.relVac < 1e-6, ...
        'the no-separation run reported %.3e relative disagreement',inNos.relVac);

%% ---------------------------------------------------------------------
%% 9. The override hook itself
%% ---------------------------------------------------------------------
%% A misspelt parameter must raise rather than silently leave the shipped
%% value in place, or a future test could believe it flew a case it did not:
             threw = false;
    try
        run_ballistic(struct('latLaunchDeg',45));
    catch
             threw = true;
    end
    assert(threw,'an unrecognised override field must raise, not be ignored');
end

function val = summaryNumber(txt,pat)
%% Purpose:
%
%  Pull one number out of a captured run_ballistic summary. Parsing the
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
%  txt              Char [1 x n]                Captured run_ballistic summary
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
