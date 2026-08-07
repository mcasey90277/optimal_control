function test_boostEvents()
%% Purpose:
%
%  Verify the burnout and apogee ODE events, both as the three-output
%  interface ode45 requires and, more importantly, behaviourally: propagate
%  a real boost phase through coorbital.prop.phaseRun and confirm it
%  terminates at the TOTAL burnout mass (payload plus dry booster structure,
%  never dry structure alone) at the time the propellant budget itself
%  predicts; and propagate a lofted unpowered arc and confirm apogee is
%  found where the flight-path angle actually crosses zero descending, with
%  a still-climbing vehicle left untouched.
%
%  Every literal below is worked out by hand from
%  coorbital.util.boosterDefaults and coorbital.util.missileConst, or is a
%  quantity recomputed independently of the event under test (an integrated
%  trajectory's own recorded time, an independently re-evaluated dr/dt).
%  None is obtained by re-running an event's own formula on its own outputs.
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
               veh = coorbital.util.vehicleDefaults();
               bst = coorbital.util.boosterDefaults();

%% =========================================================================
%% eventBurnout -- interface. mBurnout is an arbitrary literal, deliberately
%  NOT equal to any single field of bst or veh, so a mutant that reaches back
%  into bst.massDry instead of using the passed argument cannot pass by
%  accident:
%% =========================================================================
            xProbeB = [c.rE + 5e3; 0; 0; 300; deg2rad(60); deg2rad(90); 5000];
    [vB,itB,dirB] = coorbital.prop.eventBurnout(0,xProbeB,3000);
    assert(abs(vB - 2000) < 1e-9, ...
        'eventBurnout value = %.6f, expected 5000 - 3000 = 2000 kg',vB);
    assert(itB == 1,'eventBurnout must be terminal');
    assert(dirB == -1,'eventBurnout direction must be -1 (mass only decreases)');

%% The value must vanish exactly at mBurnout, and mBurnout here is the
%  documented TOTAL (payload + dry structure), not bst.massDry alone:
              mBurnout = veh.mass + bst.massDry;
    assert(mBurnout == 2400,'sanity: veh.mass + bst.massDry should be 2400 kg');
              xAtBO = [c.rE + 5e3; 0; 0; 300; deg2rad(60); deg2rad(90); mBurnout];
    [vZeroB,~,~] = coorbital.prop.eventBurnout(0,xAtBO,mBurnout);
    assert(abs(vZeroB) < 1e-9,'eventBurnout must vanish at mBurnout, got %.6e',vZeroB);

%% Using bst.massDry ALONE (the payload-omitting form the header warns
%  against) would put the crossing at 1500 kg instead of 2400 kg -- confirm
%  the two are different numbers, or this whole test would be vacuous:
    assert(bst.massDry ~= mBurnout, ...
        'bst.massDry must differ from the total burnout mass for this test to bite');

%% =========================================================================
%% eventApogee -- interface. A nonzero gamma probes the value formula
%  directly; a nonzero-value 6-state vector confirms the event reads x(5)
%  regardless of state width:
%% =========================================================================
            xProbeA = [c.rE + 40e3; 0; 0; 4000; deg2rad(12); deg2rad(90)];
    [vA,itA,dirA] = coorbital.prop.eventApogee(0,xProbeA);
    assert(abs(vA - deg2rad(12)) < 1e-12, ...
        'eventApogee value = %.9f, expected gamma = 12 deg in radians',vA);
    assert(itA == 1,'eventApogee must be terminal');
    assert(dirA == -1,'eventApogee direction must be -1 (descending crossings only)');

              xAtApo = [c.rE + 40e3; 0; 0; 4000; 0; deg2rad(90)];
    [vZeroA,~,~] = coorbital.prop.eventApogee(0,xAtApo);
    assert(abs(vZeroA) < 1e-12,'eventApogee must vanish at gamma = 0, got %.6e',vZeroA);

%% A seven-state (boost) vector must work identically -- only x(5) is read:
            xProbeA7 = [c.rE + 40e3; 0; 0; 4000; deg2rad(-7); deg2rad(90); 12345];
    [vA7,~,~] = coorbital.prop.eventApogee(0,xProbeA7);
    assert(abs(vA7 - deg2rad(-7)) < 1e-12, ...
        'eventApogee on a 7-state vector must still read x(5), got %.9f',vA7);

%% =========================================================================
%% BEHAVIOURAL: a real boost phase, propagated through phaseRun with the
%  actual constant-mass-flow engine, must terminate at x(7) = mBurnout, and
%  the elapsed time must match the propellant budget's own prediction --
%  the assertion that ties the event to the physics rather than to itself.
%% =========================================================================
          envB.atmos = @coorbital.atmos.expAtmos;
           envB.grav = @coorbital.grav.sphereGrav;
           envB.aero = @coorbital.aero.constLD;
           envB.prop = @coorbital.prop.constThrust;
         envB.omegaE = 0;
            schedB = struct('tGrid',[0 200],'alpha',[0 0],'sigma',[0 0]);
             ph.eom = @coorbital.eom.boost3DOF;
           ph.guide = @(t,x) coorbital.guide.prescribed(t,x,schedB);
       ph.terminate = @(t,x) coorbital.prop.eventBurnout(t,x,mBurnout);
           ph.tspan = [0 200];

%% Off-vertical, off-zero-alpha-and-sigma-free lofted liftoff state, matching
%  boost3DOF's own self-demo convention: gamma = 90 deg is singular (heading
%  is undefined in vertical flight), so the state starts already pitched over:
                m0 = veh.mass + bst.massDry + bst.massProp;
    assert(m0 == 32400,'sanity: liftoff mass should be 32400 kg');
               x0B = [c.rE + 5e3; 0; 0; 300; deg2rad(60); deg2rad(90); m0];
              traj = coorbital.prop.phaseRun(ph,x0B,bst,envB);

%% It stopped on the event, not the 200 s horizon:
    assert(traj.t(end) < 200, ...
        'boost phase hit the 200 s horizon instead of the burnout event');

%% Mass at termination is the TOTAL burnout mass to within 1e-6 kg:
    assert(abs(traj.x(end,7) - mBurnout) < 1e-6, ...
        'terminal mass = %.9f kg, expected mBurnout = %.9f kg', ...
        traj.x(end,7),mBurnout);

%% Mass decreased monotonically throughout -- the propellant only ever burns:
    assert(all(diff(traj.x(:,7)) <= 0), ...
        'mass must never increase during the boost phase');
    assert(traj.x(1,7) == m0,'the boost phase must start at the full liftoff mass');

%% The burn duration, read off the trajectory's OWN recorded clock, must
%  match bst.massProp/mdot -- computed here directly from boosterDefaults and
%  missileConst, never by rearranging eventBurnout's or phaseRun's output.
%  Hand arithmetic (Sutton & Biblarz Ch. 2-3): exhaust speed = Isp*g0 =
%  260*9.80665 = 2549.729 m/s; mdot = 950000/2549.729 = 372.5886162803969
%  kg/s; burn time = 30000/372.5886162803969 = 80.51775789473683 s:
            mdotVal = bst.thrustVac./(bst.Isp.*c.g0);
    assert(abs(mdotVal - 372.5886162803969) < 1e-9, ...
        'sanity: mdot should be 372.5886162803969 kg/s, got %.10f',mdotVal);
        burnTimeExpected = bst.massProp./mdotVal;
    assert(abs(burnTimeExpected - 80.51775789473683) < 1e-9, ...
        'sanity: burn time should be 80.51775789473683 s, got %.11f',burnTimeExpected);

    assert(abs(traj.t(end) - burnTimeExpected) < 1e-6, ...
        ['boost phase terminated at %.9f s against the propellant budget''s ' ...
         '%.9f s -- the event is not pinned to the physics'], ...
        traj.t(end),burnTimeExpected);

%% =========================================================================
%% BEHAVIOURAL: a lofted, unpowered arc must terminate at apogee -- gamma
%  crossing zero descending, with dr/dt changing sign there -- and a
%  still-climbing truncation of the SAME arc must not trigger at all.
%% =========================================================================
         envA.atmos = @coorbital.atmos.expAtmos;
          envA.grav = @coorbital.grav.sphereGrav;
          envA.aero = @coorbital.aero.constLD;
        envA.omegaE = 0;
             schedA = struct('tGrid',[0 5000],'alpha',[0 0],'sigma',[0 0]);
               odeA = @(t,x) coorbital.eom.glide3DOF(t,x, ...
                             coorbital.guide.prescribed(t,x,schedA),veh,envA);
              xLoft = [c.rE + 30e3; 0; 0; 4000; deg2rad(15); deg2rad(90)];

           optsApo = odeset('RelTol',1e-10,'AbsTol',1e-10, ...
                            'Events',@(t,x) coorbital.prop.eventApogee(t,x));
    [tApo,xApo,teApo] = ode45(odeA,[0 4000],xLoft,optsApo);

%% Exactly one termination, before the horizon:
    assert(numel(teApo) == 1,'expected exactly one apogee event, got %d',numel(teApo));
    assert(tApo(end) < 4000,'the lofted arc hit the horizon instead of apogee');

%% Gamma at termination is (numerically) exactly zero:
    assert(abs(xApo(end,5)) < 1e-6, ...
        'gamma at the apogee event = %.6e rad, expected ~0',xApo(end,5));

%% dr/dt changes sign AT the event, checked against an INDEPENDENT
%  re-integration that carries no event at all, so the crossing time is not
%  taken from the event's own output. rdot = V sin(gamma) is evaluated one
%  second either side of the recorded event time by interpolating that
%  independent solution:
              solApo = ode45(odeA,[0 4000],xLoft, ...
                             odeset('RelTol',1e-10,'AbsTol',1e-10));
               rdotF = @(xs) xs(4).*sin(xs(5));
              xPreEv = deval(solApo,tApo(end) - 1);
             xPostEv = deval(solApo,tApo(end) + 1);
          rdotPreEv = rdotF(xPreEv);
         rdotPostEv = rdotF(xPostEv);
    assert(rdotPreEv > 0, ...
        'dr/dt one second before apogee should be positive (still climbing), got %.6f', ...
        rdotPreEv);
    assert(rdotPostEv < 0, ...
        'dr/dt one second after apogee should be negative (descending), got %.6f', ...
        rdotPostEv);

%% A still-climbing truncation of the SAME arc, cut off at half the time to
%  apogee, must NOT trigger the event at all:
             tHalf = tApo(end)/2;
    [tClimb,xClimb,teClimb,~,~] = ode45(odeA,[0 tHalf],xLoft,optsApo);
    assert(isempty(teClimb), ...
        'a still-climbing vehicle triggered eventApogee at t = %.3f s',teClimb);
    assert(abs(tClimb(end) - tHalf) < 1e-9, ...
        'the still-climbing run must reach its horizon, not stop early');
    assert(xClimb(end,5) > 0, ...
        'the still-climbing run must end with a positive flight-path angle, got %.6f rad', ...
        xClimb(end,5));
end
