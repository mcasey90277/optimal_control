function test_runGlide()
%% Purpose:
%
%  Pin the COMPOSITION. Every other test in this suite exercises one library
%  function in isolation; this one exercises the entry script that wires them
%  together -- unit conversion, vehicle selection, environment assembly, the
%  propagation, the derived-quantity block, the great-circle call site, the
%  peak search, the termination diagnosis, and the printed summary. Without
%  it the suite passes unchanged while run_glide transposes latitude and
%  longitude, mislabels a peak, or reports the wrong range.
%
%  WHY THERE ARE TWO RUNS, AND WHY THE SECOND ONE IS NOT DUE EAST. The
%  shipped configuration is a due-east glide launched from (0 deg, 0 deg).
%  That geometry is provably BLIND to a transposition at the greatCircle call
%  site: with the entry point at the origin the central angle reduces to
%  acos(cos(lat2) cos(lon2)), which is symmetric in lat2 and lon2, so swapping
%  the two terminal arguments leaves the reported range bit-identical.
%  Measured, not assumed -- swapping them changes the printed summary by
%  nothing at all. Asserting terminal latitude and longitude separately does
%  not close the hole either: those come from traj, which the mutation never
%  touches. The only thing that catches it is flying a geometry in which the
%  four coordinates are genuinely distinct, which is what the second run does.
%  Its heading is south-east rather than east, which as a bonus exercises the
%  cos(psi) turning terms that a due-east arc leaves at zero.
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
%  Michael Casey                                                08/06/2026
%  Copyright 2026 Coorbital, Inc.
%% ------------------------ Begin Code Sequence ---------------------------

                 c = coorbital.util.missileConst();

%% ---------------------------------------------------------------------
%% 1. vehicle_hgv must not silently diverge from the library defaults
%% ---------------------------------------------------------------------
%% vehicle_hgv deliberately re-states the four generic placeholders so that
%% real numbers have exactly one home when they arrive. The duplication is by
%% design, but until now nothing guarded it: editing one copy and not the
%% other would leave two different vehicles in the tree with no complaint.
%% Exact equality is right here -- both sides are the same decimal literals,
%% and any difference at all is the divergence this check exists to catch:
               veh = vehicle_hgv();
               vdf = coorbital.util.vehicleDefaults();
         duplicate = {'mass','Sref','CL','LD'};
    for kf = 1:numel(duplicate)
              name = duplicate{kf};
        assert(isfield(veh,name) && isfield(vdf,name), ...
            'both vehicle_hgv and vehicleDefaults must define %s',name);
        assert(veh.(name) == vdf.(name), ...
            ['vehicle_hgv.%s = %.17g has diverged from ' ...
             'vehicleDefaults.%s = %.17g. The two are duplicated ON PURPOSE ' ...
             'and must be changed together; see the note in vehicle_hgv.m'], ...
            name,veh.(name),name,vdf.(name));
    end

%% Fields vehicle_hgv does not override must still arrive from the defaults:
    assert(veh.noseRadius == vdf.noseRadius, ...
        'noseRadius must pass through vehicle_hgv unchanged');

%% ...but agreeing with each other is not the same as being RIGHT. The check
%% above pins the two files together and pins neither to a number, so a
%% retune applied consistently to both passes it in silence. The three
%% aerodynamic placeholders below are the library's only physics-bearing
%% constants with no absolute anchor anywhere, and the analytic suite cannot
%% supply one: test_equilibriumGlide builds its closed-form Veq(r) from
%% veh.mass, veh.Sref and veh.CL, so all three cancel between the model and
%% its own reference and the check passes at any values. That is exactly the
%% Hscale blindness of docs/LESSONS_LEARNED.md, and the cure is the same one
%% -- pin the constant where it is DEFINED, so a retune fails here and says
%% why, rather than surfacing as an unexplained shift in a pinned trajectory
%% literal further down this file. veh.mass needs no pin here; it is already
%% pinned to 900 kg in test_constThrust, where the booster mass bookkeeping
%% is quoted against it:
    assert(vdf.Sref == 0.75,'expected vehicleDefaults Sref = 0.75 m^2');
    assert(vdf.CL   == 0.35,'expected vehicleDefaults CL   = 0.35');
    assert(vdf.LD   == 2.5, 'expected vehicleDefaults LD   = 2.5');

%% ---------------------------------------------------------------------
%% 2. The shipped configuration
%% ---------------------------------------------------------------------
%% Plots off so the suite creates no figures, and the summary captured rather
%% than printed so the suite output stays readable. The captured text is not
%% discarded -- it is the only place the reported range, angle and peaks
%% appear, so it is parsed and asserted below:
            outDef = evalc('trDef = run_glide(struct(''showPlots'',false));');

%% Termination was nominal, not a horizon timeout. Both halves matter: the
%% nominal wording must be present AND the truncation banner absent, so that a
%% future change to the diagnosis cannot pass by printing both:
    assert(contains(outDef,'(nominal)'), ...
        'the shipped run did not terminate nominally. Summary was:\n%s',outDef);
    assert(~contains(outDef,'CAUTION'), ...
        'the shipped run printed the truncation caution. Summary was:\n%s',outDef);
    assert(~contains(outDef,'integration horizon'), ...
        'the shipped run hit the integration horizon. Summary was:\n%s',outDef);

%% Terminal altitude is the configured stop altitude. eventAltitude is a
%% terminal one-sided event, so this is exact to the event solver, not to a
%% sampling interval:
             hEndM = trDef.x(end,1) - c.rE;
    assert(abs(hEndM - 5000) < 1e-6, ...
        'terminal altitude %.6f m, expected the configured 5000 m stop',hEndM);

%% Headline results, to 1e-4 relative:
    assertRel(trDef.t(end),2073.7667711260506,1e-4,'flight time (s)');

              rngD = summaryNumber(outDef,'ground range +([-\d.]+) +km');
    assertRel(rngD,6986.8174506370033,1e-4,'reported ground range (km)');

              angD = summaryNumber(outDef,'central angle +([-\d.]+) +deg');
    assertRel(angD,62.763649032602110,1e-4,'reported central angle (deg)');

%% The reported range must be the reported angle on the reported sphere. This
%% catches a units slip between the two lines, which agreeing with the pinned
%% literals separately would not:
    assertRel(rngD,c.rE.*deg2rad(angD)./1000,1e-5,'range against its own angle');

%% Peak sensed load, recomputed from the trajectory so that the pinned literal
%% is checked at full precision rather than at the summary's two decimals:
    [nMax,tAtNMax] = peakLoad(trDef,veh,c);
    assertRel(nMax,1.1079843546543626,1e-4,'peak load factor (g)');
    assertRel(nMax.*c.g0,10.865614771571204,1e-4,'peak sensed load (m/s^2)');
    assertRel(tAtNMax,1178.648029900691,1e-4,'time of peak load (s)');

%% ...and the summary must REPORT that peak, at both the printed precision of
%% the g value and of the m/s^2 value. This is the half the recomputation
%% cannot see: a peak found correctly and then printed from the wrong variable:
            decTok = regexp(outDef, ...
                     ['deceleration +([-\d.]+) +g +\(([-\d.]+) m/s\^2 sensed ' ...
                      'aero load, at t = ([-\d.]+) s\)'],'tokens','once');
    assert(numel(decTok) == 3, ...
        'could not parse the deceleration line. Summary was:\n%s',outDef);
    assert(abs(str2double(decTok{1}) - nMax) < 0.005 + 1e-9, ...
        'summary reported %s g against the trajectory''s %.6f g',decTok{1},nMax);
    assert(abs(str2double(decTok{2}) - nMax.*c.g0) < 0.005 + 1e-9, ...
        'summary reported %s m/s^2 against the trajectory''s %.6f', ...
        decTok{2},nMax.*c.g0);
    assert(abs(str2double(decTok{3}) - tAtNMax) < 0.05 + 1e-9, ...
        'summary put the peak at t = %s s against the trajectory''s %.4f s', ...
        decTok{3},tAtNMax);

%% The terminal block of the summary, parsed from the tail so the entry block's
%% identically-labelled lines cannot be matched by mistake:
           tailDef = summaryTail(outDef);
    assertAbs(summaryNumber(tailDef,'altitude +([-\d.]+) +km'),5, ...
        1e-3,'reported terminal altitude (km)');
    assertRel(summaryNumber(tailDef,'speed +([-\d.]+) +m/s'), ...
        321.00896892076537,1e-4,'reported terminal speed (m/s)');

%% TERMINAL LATITUDE AND LONGITUDE, SEPARATELY. A due-east glide launched from
%% the equator must end at latitude zero with the whole of the motion in
%% longitude. Asserting the range alone cannot tell those two apart -- it is a
%% single scalar built from both -- so they are pinned individually. This is
%% what catches a transposition in the STATE indexing, x(:,2) against x(:,3);
%% the call-site order is a different failure and is covered in part 3:
    assertAbs(rad2deg(trDef.x(end,3)),0,1e-9,'terminal latitude (deg)');
    assertRel(rad2deg(trDef.x(end,2)),62.76364903260211,1e-4, ...
        'terminal longitude (deg)');
    assertAbs(summaryNumber(tailDef,'latitude +([-\d.]+) +deg'),0,1e-3, ...
        'reported terminal latitude (deg)');
    assertRel(summaryNumber(tailDef,'longitude +([-\d.]+) +deg'), ...
        62.76364903260211,1e-4,'reported terminal longitude (deg)');

%% The arc never leaves the equator at all, which is the premise of the two
%% assertions above rather than an accident of the endpoint:
    assert(max(abs(trDef.x(:,3))) < 1e-12, ...
        'the due-east equatorial arc wandered to %.3e rad of latitude', ...
        max(abs(trDef.x(:,3))));

%% Terminal flight path angle, and a monotone descent -- the peak altitude is
%% the entry altitude, reached at t = 0:
    assertRel(rad2deg(trDef.x(end,5)),-12.59797093289877,1e-4, ...
        'terminal flight path angle (deg)');
    assertAbs(max(trDef.x(:,1) - c.rE)./1000,60,1e-9,'peak altitude (km)');

%% Sample count is left deliberately loose. It is an ode45 adaptive-step
%% artefact, not a physical result, and pinning it would make the suite
%% hostage to a solver revision. It is checked only for being a real
%% trajectory rather than a two-point degenerate one:
    assert(numel(trDef.t) > 1000, ...
        'only %d samples returned; that is not a resolved trajectory', ...
        numel(trDef.t));
    assert(isequal(size(trDef.x),[numel(trDef.t) 6]), ...
        'state history must be [N x 6]');
    assert(isequal(size(trDef.u),[numel(trDef.t) 2]), ...
        'control history must be [N x 2]');

%% ---------------------------------------------------------------------
%% 3. An off-axis run, which is what pins the greatCircle argument order
%% ---------------------------------------------------------------------
%% Same vehicle, same entry altitude, speed and flight path angle; only the
%% entry point and heading move, to a place where latitude and longitude are
%% four genuinely distinct numbers:
            latOff = 25;
            lonOff = -100;
            psiOff = 135;
            outOff = evalc(['trOff = run_glide(struct(''latEntry'',latOff,' ...
                            '''lonEntry'',lonOff,''psiEntry'',psiOff,' ...
                            '''showPlots'',false));']);

    assert(contains(outOff,'(nominal)') && ~contains(outOff,'CAUTION'), ...
        'the off-axis run did not terminate nominally. Summary was:\n%s',outOff);

%% The angle the call site SHOULD have produced, formed here in the documented
%% (lat1,lon1,lat2,lon2) order:
              lat1 = deg2rad(latOff);
              lon1 = deg2rad(lonOff);
              lat2 = trOff.x(end,3);
              lon2 = trOff.x(end,2);
            angRef = coorbital.util.greatCircle(lat1,lon1,lat2,lon2);

%% Every way of getting the order wrong, so the assertion that follows is
%% shown to be discriminating rather than assumed to be. If a future change to
%% the entry point ever makes one of these coincide with the truth again, this
%% guard fails loudly instead of the order check going quietly blind -- which
%% is exactly what happens in the shipped due-east geometry:
           angSwap = [coorbital.util.greatCircle(lat1,lon1,lon2,lat2); ...
                      coorbital.util.greatCircle(lon1,lat1,lat2,lon2); ...
                      coorbital.util.greatCircle(lon1,lat1,lon2,lat2)];
    assert(min(abs(rad2deg(angSwap - angRef))) > 5, ...
        ['the off-axis geometry is not discriminating: a transposed ' ...
         'greatCircle call comes within %.3f deg of the correct %.4f deg'], ...
        min(abs(rad2deg(angSwap - angRef))),rad2deg(angRef));

%% ...and now the check itself. The summary prints the angle to four decimals,
%% so 1e-3 deg is comfortably inside the printed resolution while a
%% transposition is off by tens of degrees:
            angOff = summaryNumber(outOff,'central angle +([-\d.]+) +deg');
    assert(abs(angOff - rad2deg(angRef)) < 1e-3, ...
        ['reported central angle %.4f deg against %.4f deg from the ' ...
         'terminal state in (lat,lon) order. The arguments at the ' ...
         'greatCircle call site in run_glide look transposed.'], ...
        angOff,rad2deg(angRef));

%% Range and angle are invariant under moving the entry point and heading on a
%% non-rotating sphere: the same vehicle flying the same entry state must fly
%% the same arc length wherever it is put. An independent statement from the
%% pinned literal, and one no single-run test can make:
              rngO = summaryNumber(outOff,'ground range +([-\d.]+) +km');
    assertRel(rngO,rngD,1e-6,'off-axis ground range against the shipped run');
    assertRel(angOff,angD,1e-6,'off-axis central angle against the shipped run');

%% The arc really did fly off the equator and really did turn, so the cos(psi)
%% terms in latdot, Vdot and gammadot were exercised rather than multiplied by
%% zero as they are on the due-east arc:
    assert(abs(rad2deg(lat2) - (-22.109991177599976)) < 1e-3, ...
        'off-axis terminal latitude %.6f deg, expected -22.1100 deg',rad2deg(lat2));
    assert(abs(rad2deg(lon2) - (-57.264873193306954)) < 1e-3, ...
        'off-axis terminal longitude %.6f deg, expected -57.2649 deg',rad2deg(lon2));
    assert(abs(rad2deg(lat2) - rad2deg(lon2)) > 10, ...
        'terminal latitude and longitude are too close together to tell apart');
    assert(abs(rad2deg(trOff.x(end,6)) - psiOff) > 1, ...
        'heading never changed; the tan(lat) turning term was not exercised');

%% ---------------------------------------------------------------------
%% 4. The override hook itself
%% ---------------------------------------------------------------------
%% A misspelt parameter must raise rather than silently leave the shipped
%% value in place, or a future test could believe it flew a case it did not:
             threw = false;
    try
        run_glide(struct('latEntryDeg',25));
    catch
             threw = true;
    end
    assert(threw,'an unrecognised override field must raise, not be ignored');
end

function [nMax,tAtNMax] = peakLoad(traj,veh,c)
%% Purpose:
%
%  Recompute the peak sensed aerodynamic load factor and the time at which it
%  occurs, from the returned trajectory and the vehicle's declared parameters.
%  Gravity is not sensed and is excluded, matching what run_glide reports.
%
%  This duplicates the summary's arithmetic on purpose: the duplicate is
%  compared against a pinned absolute literal in the caller, so a fault
%  present in both copies still fails. The duplicate exists only to recover
%  full precision, which the summary's two printed decimals discard.
%
%% Inputs:
%
%  traj             Struct                      Trajectory from run_glide;
%                                               uses t [N x 1] (s), x [N x 6],
%                                               u [N x 2] (rad)
%
%  veh              Struct                      Vehicle; uses mass (kg),
%                                               Sref (m^2)
%
%  c                Struct                      Constants from missileConst;
%                                               uses rE (m), g0 (m/s^2)
%
%% Outputs:
%
%  nMax             [1 x 1]                     Peak sensed aerodynamic load
%                                               factor (g)
%
%  tAtNMax          [1 x 1]                     Time of that peak (s)
%
%% Revision History:
%  Michael Casey                                                08/06/2026
%  Copyright 2026 Coorbital, Inc.
%% ------------------------ Begin Code Sequence ---------------------------

                nS = numel(traj.t);
    [rho,~,~,aSnd] = coorbital.atmos.expAtmos(traj.x(:,1) - c.rE);
                 V = traj.x(:,4);
              qbar = 0.5.*rho.*V.^2;
              mach = V./aSnd;
               CLv = zeros(nS,1);
               CDv = zeros(nS,1);
    for ks = 1:nS
        [CLv(ks),CDv(ks)] = coorbital.aero.constLD(traj.u(ks,1),mach(ks),veh);
    end
             aLift = qbar.*veh.Sref.*CLv./veh.mass;
             aDrag = qbar.*veh.Sref.*CDv./veh.mass;
             nLoad = sqrt(aLift.^2 + aDrag.^2)./c.g0;
        [nMax,kPk] = max(nLoad);
    assert(kPk > 1 && kPk < nS,'the load peak is not interior to the arc');
           tAtNMax = traj.t(kPk);
end

function tail = summaryTail(txt)
%% Purpose:
%
%  Return the part of a run_glide summary from the "Terminal conditions"
%  heading onwards. The entry block carries identically labelled altitude,
%  speed, latitude and longitude lines, so a pattern matched against the whole
%  summary would silently read the entry values instead.
%
%% Inputs:
%
%  txt              Char [1 x n]                Captured run_glide summary
%
%% Outputs:
%
%  tail             Char [1 x m]                Summary from the terminal
%                                               heading to the end
%
%% Revision History:
%  Michael Casey                                                08/06/2026
%  Copyright 2026 Coorbital, Inc.
%% ------------------------ Begin Code Sequence ---------------------------

                k0 = strfind(txt,'Terminal conditions');
    assert(~isempty(k0), ...
        'the summary has no "Terminal conditions" block:\n%s',txt);
              tail = txt(k0(1):end);
end

function val = summaryNumber(txt,pat)
%% Purpose:
%
%  Pull one number out of a captured run_glide summary. Parsing the printed
%  text rather than recomputing the quantity is deliberate: what a reader of
%  this library acts on is the summary, so the summary is what gets asserted.
%
%  The patterns use explicit space runs rather than \s, because in MATLAB \s
%  also matches a newline and would happily reach across into the next line
%  of the report.
%
%% Inputs:
%
%  txt              Char [1 x n]                Captured run_glide summary
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
%  Michael Casey                                                08/06/2026
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
%  Michael Casey                                                08/06/2026
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
%  Michael Casey                                                08/06/2026
%  Copyright 2026 Coorbital, Inc.
%% ------------------------ Begin Code Sequence ---------------------------

    assert(abs(got - want) < tol, ...
        '%s: got %.12g against %.12g, %.3e absolute, budget %.1e', ...
        what,got,want,abs(got - want),tol);
end
