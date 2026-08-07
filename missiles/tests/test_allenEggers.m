function test_allenEggers()
%% Purpose:
%
%  Validate a zero-lift ballistic entry against the Allen-Eggers closed-form
%  solution -- again, a reference derived independently of the propagator. For
%  a steep entry into an isothermal exponential atmosphere, neglecting gravity
%  in the along-track equation and holding the flight path angle fixed, the
%  peak deceleration is
%
%      aMax = V_e^2 |sin(gamma_e)| / (2 e H)
%
%  with e Euler's number and H the density scale height. The striking feature
%  is what is absent: no mass, no reference area, no drag coefficient. The
%  peak deceleration is independent of the ballistic coefficient. That
%  independence is the most distinctive property this test has, so it is
%  exercised directly -- the SAME entry is flown at two ballistic coefficients
%  eight times apart, and both must reproduce the same peak.
%
%  MEASURING THE RIGHT QUANTITY. Allen-Eggers predicts the AERODYNAMIC
%  deceleration alone. Differencing the propagated speed history gives the NET
%  rate, which also carries gravity's along-track component and does not match.
%  The drag term 0.5 rho V^2 S CD / m is formed explicitly instead, from the
%  propagated density and speed and the vehicle's declared parameters.
%
%  TOLERANCE, AND WHY IT IS 10 PERCENT AND NOT 2. The closed form drops
%  gravity from the speed equation and freezes the flight path angle, and both
%  assumptions leak at a 30 deg entry. Their sizes are known, not guessed. At
%  the higher ballistic coefficient the measured peak occurs at V = 3736 m/s
%  against the closed form's V_e/sqrt(e) = 3639 m/s -- gravity has added
%  2.7 percent to the speed on the way down, which is 5.4 percent on V^2 --
%  and at gamma = -31.4 deg against the assumed -30 deg, which is a further
%  4.2 percent on |sin(gamma)|. Together those predict +9.8 percent; the
%  measured excess is +10.5 percent. The assertion is therefore set at
%  12 percent, and it is made ONE-SIDED BELOW as well: both neglected effects
%  can only ADD speed and steepen the path for a descending ballistic entry,
%  so the simulated peak must exceed the closed form. A simulated peak that
%  came in under the analytic value would mean too much drag, not a better
%  approximation.
%
%  PEAK LOCATION, the sharper of the two closed-form statements. Allen-Eggers
%  also places the peak, at the density
%
%      rho_pk = beta |sin(gamma_e)| / H,     beta = m / (S CD)
%
%  and this is compared in DENSITY, not in altitude. Altitude is the logarithm
%  of the quantity actually predicted, so a relative altitude error is not a
%  consistent measure of the same physical discrepancy: the two cases below
%  disagree by 4.84 and 4.36 percent in density, but the SAME disagreements
%  read as 6.55 and 1.52 percent in altitude, purely because the log offset is
%  being divided by two different absolute altitudes. A tolerance set in
%  altitude would mean different things at different beta. In density it means
%  one thing, and 6 percent is the budget.
%
%  At 4.8 percent this is roughly twice as sharp as the 10.5 percent
%  peak-VALUE comparison. The two are not redundant: the peak value is
%  beta-independent, whereas this is a statement about where the peak sits for
%  a GIVEN beta, so it is the one that would notice a drag area or a mass that
%  the beta-independent check is deliberately built to ignore.
%
%  WHAT THIS TEST CANNOT SEE: THE SCALE HEIGHT. Every comparison in this file
%  is structurally blind to H. It appears in the closed form and in the
%  atmosphere, so changing it moves prediction and simulation together --
%  measured, inflating c.Hscale by 15 percent moves the peak-value agreement
%  only from +10.50 to +10.44 percent, and the peak-location check shifts by
%  a comparable trickle. This is not a defect that a tighter tolerance fixes;
%  Allen-Eggers IS the exact solution for that atmosphere, so H cancels. The
%  guard lives where it belongs, in test_missileConst, which PINS Hscale to
%  7200 m rather than range-checking it. An earlier draft of this file tried
%  to anchor H here against the hydrostatic value R*T0/g0; that was wrong on
%  two counts, and is recorded in LESSONS_LEARNED.md -- Allen-Eggers needs
%  only an exponential DENSITY profile, so temperature and pressure are not
%  preconditions of it at all, and the library's 7200 m sits 1.64 percent from
%  R*T0/g0 = 7317.8 m anyway, so the anchor would have been asserting little
%  while false-alarming on any legitimate change to T0.
%
%% Inputs:
%
%  none
%
%% Outputs:
%
%  none                                         Throws on any failed assertion
%
%% References:
%   [1] Allen, H.J., Eggers, A.J., "A Study of the Motion and Aerodynamic
%       Heating of Ballistic Missiles Entering the Earth's Atmosphere at High
%       Supersonic Speeds," NACA TR-1381, 1958.
%   [2] Vinh, N.X., Busemann, A., Culp, R.D., "Hypersonic and Planetary Entry
%       Flight Mechanics," Univ. Michigan Press, 1980, Ch. 6.
%
%% Revision History:
%  Michael Casey                                                08/06/2026
%  Copyright 2026 Coorbital, Inc.
%% ------------------------ Begin Code Sequence ---------------------------

                 c = coorbital.util.missileConst();
               veh = coorbital.util.vehicleDefaults();

%% Zero lift so the vehicle falls ballistically. CD is left at the vehicle's
%% own CL/LD; the ballistic coefficient is varied below through the reference
%% area, which is the one knob that changes beta and nothing else:
         env.atmos = @coorbital.atmos.expAtmos;
          env.grav = @coorbital.grav.sphereGrav;
          env.aero = @(al,ma,vv) deal(0,vv.CL/vv.LD);
        env.omegaE = 0;

%% Steep entry, where the frozen-gamma and no-gravity assumptions are least
%% bad:
            gammaE = deg2rad(-30);
                vE = 6000;
                h0 = 120e3;

%% Allen-Eggers closed form. Nothing in it depends on the vehicle:
           aMaxRef = vE.^2.*abs(sin(gammaE))./(2.*exp(1).*c.Hscale);

%% Two ballistic coefficients eight times apart. Each run is carried well past
%% its OWN peak: the heavy case peaks near 5 km so it must reach the ground,
%% while the light case peaks near 20 km and is stopped at 5 km, short of the
%% near-vertical terminal fall where cos(gamma) approaches the singularity the
%% 3DOF equations guard against:
         srefScale = [1; 8];
             hStop = [0; 5e3];
                nB = numel(srefScale);
           aMaxNum = zeros(nB,1);
             hPeak = zeros(nB,1);
             vPeak = zeros(nB,1);
             gPeak = zeros(nB,1);
           rhoPeak = zeros(nB,1);
             betaB = zeros(nB,1);

    for kb = 1:nB
                vv = veh;
           vv.Sref = veh.Sref.*srefScale(kb);
           [~,CDb] = env.aero(0,0,vv);
         betaB(kb) = vv.mass./(vv.Sref.*CDb);

             sched = struct('tGrid',[0 600],'alpha',[0 0],'sigma',[0 0]);
            ph.eom = @coorbital.eom.glide3DOF;
          ph.guide = @(t,x) coorbital.guide.prescribed(t,x,sched);
      ph.terminate = @(t,x) coorbital.prop.eventAltitude(t,x,hStop(kb));
          ph.tspan = [0 600];
                x0 = [c.rE + h0; 0; 0; vE; gammaE; deg2rad(90)];
              traj = coorbital.prop.phaseRun(ph,x0,vv,env);

        assert(traj.t(end) < 600, ...
            'case %d hit the time horizon instead of the altitude event',kb);

%% Peak AERODYNAMIC deceleration, formed explicitly rather than differenced
%% out of the speed history:
              rArc = traj.x(:,1);
              vArc = traj.x(:,4);
    [rhoArc,~,~,~] = coorbital.atmos.expAtmos(rArc - c.rE);
             accel = 0.5.*rhoArc.*vArc.^2.*vv.Sref.*CDb./vv.mass;
         [aPk,kPk] = max(accel);

%% The peak must be interior, or the run was stopped before the vehicle
%% reached it and the maximum is an artefact of where integration ended:
        assert(kPk > 1 && kPk < numel(accel), ...
            'case %d peaked at sample %d of %d; the peak is not interior', ...
            kb,kPk,numel(accel));

       aMaxNum(kb) = aPk;
         hPeak(kb) = rArc(kPk) - c.rE;
         vPeak(kb) = vArc(kPk);
         gPeak(kb) = traj.x(kPk,5);
       rhoPeak(kb) = rhoArc(kPk);
    end

%% Each case against the closed form, one-sided both ways:
    for kb = 1:nB
               rel = (aMaxNum(kb) - aMaxRef)./aMaxRef;
        assert(rel > 0, ...
            ['beta = %.0f kg/m^2 peaked at %.1f m/s^2, BELOW the Allen-Eggers ' ...
             '%.1f m/s^2 (%.2f percent). Gravity and flight path steepening can ' ...
             'only raise the peak, so this means too much drag, not a better ' ...
             'approximation'],betaB(kb),aMaxNum(kb),aMaxRef,100*rel);
        assert(rel < 0.12, ...
            ['beta = %.0f kg/m^2 peaked at %.1f m/s^2 against Allen-Eggers ' ...
             '%.1f m/s^2 (%.2f percent high, budget 12); peak at %.0f m, ' ...
             'V = %.0f m/s, gamma = %.2f deg'], ...
            betaB(kb),aMaxNum(kb),aMaxRef,100*rel,hPeak(kb),vPeak(kb), ...
            rad2deg(gPeak(kb)));
    end

%% PEAK LOCATION. Allen-Eggers puts the peak where the density reaches
%% beta |sin(gamma_e)| / H. Compared in density, not altitude -- the same two
%% discrepancies read as 4.84 and 4.36 percent here but as 6.55 and 1.52
%% percent in altitude, so an altitude budget would mean different things at
%% different beta. This is a per-beta statement, unlike the beta-independent
%% peak VALUE above, so it is a genuinely separate claim, and about twice as
%% sharp:
    for kb = 1:nB
            rhoRef = betaB(kb).*abs(sin(gammaE))./c.Hscale;
            relLoc = (rhoPeak(kb) - rhoRef)./rhoRef;
        assert(abs(relLoc) < 0.06, ...
            ['beta = %.0f kg/m^2 peaked at rho = %.5f kg/m^3 against the ' ...
             'Allen-Eggers location rho = beta|sin(gammaE)|/H = %.5f kg/m^3 ' ...
             '(%.2f percent, budget 6); peak at %.0f m'], ...
            betaB(kb),rhoPeak(kb),rhoRef,100*relLoc,hPeak(kb));
    end

%% BALLISTIC COEFFICIENT INDEPENDENCE. An
%% eightfold change in beta must leave the peak deceleration essentially
%% untouched while moving the altitude at which it occurs by kilometres. If
%% the peaks tracked beta, the agreement above would be a coincidence of the
%% chosen vehicle rather than a property of the drag model and atmosphere:
    assert(abs(betaB(1)./betaB(2) - 8) < 1e-9, ...
        'the two cases differ in beta by %.3f, expected 8',betaB(1)./betaB(2));

           aSpread = abs(aMaxNum(1) - aMaxNum(2))./mean(aMaxNum);
    assert(aSpread < 0.03, ...
        ['peak deceleration moved %.2f percent (%.1f vs %.1f m/s^2) for an ' ...
         'eightfold change in ballistic coefficient; Allen-Eggers says it ' ...
         'must not move'],100*aSpread,aMaxNum(1),aMaxNum(2));

    assert(abs(hPeak(2) - hPeak(1)) > 10e3, ...
        ['the peak altitude moved only %.0f m between beta = %.0f and ' ...
         'beta = %.0f kg/m^2; the two cases are not in different regimes, ' ...
         'so the independence check above is vacuous'], ...
        abs(hPeak(2) - hPeak(1)),betaB(1),betaB(2));

%% ...and the lower ballistic coefficient must be the one that peaks higher,
%% which is the direction Allen-Eggers gives for the peak location:
    assert(hPeak(2) > hPeak(1), ...
        ['the lighter-loaded case peaked at %.0f m, below the heavier case at ' ...
         '%.0f m; the peak altitude depends on beta the wrong way'], ...
        hPeak(2),hPeak(1));
end
