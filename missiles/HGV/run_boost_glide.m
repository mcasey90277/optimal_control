function [traj,info] = run_boost_glide(opts)
%% Purpose:
%
%  Fly a boost-glide vehicle from the pad to impact and report what it did:
%  a powered boost on a prescribed pitch program, stage separation at burnout,
%  an unpowered hypersonic glide down to a handoff altitude, and a steepened
%  terminal descent to impact. Everything a routine run needs to change lives
%  in the USER PARAMETERS block; nothing below it should require editing.
%
%  Angles and distances in the USER PARAMETERS block are in degrees and
%  kilometres because that is how a user thinks about them. They are converted
%  to the library's SI units (m, m/s, rad, s) immediately after the block, in
%  one place, and nothing below that point works in any other unit.
%
%% Note -- what makes the descent a SEPARATE phase and not an arbitrary cut:
%
%  Phases 2 and 3 run the same equations of motion on the same airframe. What
%  differs is the CONTROL and what the control does to the flight:
%
%    Phase 2, glide.   The bank angle is small, so nearly the whole lift vector
%                      points up. Lift very nearly balances the weight less the
%                      centrifugal relief, the flight path stays within a few
%                      degrees of horizontal, and the vehicle converts speed
%                      into range at the shallowest sink rate the airframe can
%                      hold. This is the segment that makes the vehicle a
%                      GLIDER: it flies far further than the ballistic arc its
%                      burnout state would otherwise put it on.
%
%    Phase 3, descent. The vehicle rolls to a large bank angle. The lift
%                      magnitude is unchanged, but only cos(sigma) of it is
%                      left holding the vehicle up, so the support collapses,
%                      the flight path steepens sharply, and the vehicle dives
%                      into the dense lower atmosphere carrying far more speed
%                      than a continued glide would have left it. Dynamic
%                      pressure and drag deceleration both rise steeply. This
%                      is the terminal segment, and it trades the range the
%                      glide was buying for a steep, fast arrival.
%
%  The summary below MEASURES that contrast -- mean flight-path angle, mean
%  sink rate, peak dynamic pressure and mean drag deceleration, phase by phase
%  -- rather than asserting it, so a reader can see whether the split earned
%  itself on this particular run.
%
%  The handoff is an altitude, and the operator sets it. It must sit BELOW
%  every trough of the glide phugoid, or the glide is cut short at its first
%  dip instead of at its end; the shipped 15 km clears the lowest earlier
%  trough of the shipped case by about 6 km. Raise it and check the glide
%  duration in the summary before believing the result.
%
%% Note -- three vehicles, one chain:
%
%  coorbital.prop.phaseRun carries ONE vehicle struct for the whole chain, but
%  the boosted stack and the separated glide vehicle have different reference
%  areas, different aerodynamic coefficients and different masses. Each phase
%  therefore BINDS its own vehicle in the equation-of-motion closure below, and
%  the struct handed to phaseRun is never read by the equations of motion. The
%  glide and the descent share one closure because they share one airframe --
%  nothing is jettisoned at the handoff, only the bank angle changes.
%
%  This is not optional bookkeeping. coorbital.eom.glide3DOF divides by
%  veh.mass and never reads the mass STATE x(7), so a chain that dropped the
%  booster from x(7) while still handing the stack to the equations of motion
%  would fly the whole glide at the wrong weight. coorbital.eom.massConstant
%  refuses to run when those two disagree, which is what turns that silent
%  error into an immediate one.
%
%% Inputs:
%
%  opts             Struct, optional            Overrides for named USER
%                                               PARAMETERS entries, given in
%                                               the SAME human units as the
%                                               block itself (deg, km). Omit
%                                               it -- or pass [] -- for the
%                                               shipped configuration. An
%                                               unrecognised field name is an
%                                               error, not a silent no-op.
%                                               Exists so an automated test or
%                                               a batch sweep can drive this
%                                               script at more than one
%                                               operating point without
%                                               editing it.
%
%% Outputs:
%
%  traj             Struct                      Trajectory from
%                                               coorbital.prop.phaseRun:
%                                               t        [N x 1] (s)
%                                               x        [N x 7] state, the
%                                                        seventh component
%                                                        being mass (kg)
%                                               u        [N x 2] control (rad)
%                                               phaseIdx [N x 1]
%                                               junction [2 x 1] struct; the
%                                                        first is the state
%                                                        AFTER separation
%
%  info             Struct                      The summary's numbers at full
%                                               precision, so a test does not
%                                               have to read them back out of
%                                               printed text; all SI except the
%                                               *Km distances. It also carries
%                                               the assembled phase array, the
%                                               environment and the liftoff
%                                               state (fields phases, env, x0,
%                                               boostVeh, glideVeh), which is
%                                               everything an independent
%                                               checker needs to re-integrate
%                                               the same chain with a different
%                                               solver.
%
%% References:
%   [1] Vinh, N.X., Busemann, A., Culp, R.D., "Hypersonic and Planetary Entry
%       Flight Mechanics," Univ. Michigan Press, 1980, Ch. 5 (gliding entry,
%       the equilibrium glide, and the skipping phugoid this trajectory flies).
%   [2] Eggers, A.J., Allen, H.J., Neice, S.E., "A Comparative Analysis of the
%       Performance of Long-Range Hypervelocity Vehicles," NACA TR-1382, 1958.
%       (Boost-glide against ballistic and skip trajectories.)
%
%% Revision History:
%  Michael Casey                                                08/07/2026
%  Copyright 2026 Coorbital, Inc.
%% ------------------------ Begin Code Sequence ---------------------------

%% Resolve paths so the script runs from anywhere:
              here = fileparts(mfilename('fullpath'));
    addpath(here);
    addpath(fullfile(here,'..'));

%% ========================= USER PARAMETERS ==============================
%% Launch site -- where the vehicle leaves the pad:
         latLaunch = 20;               %deg, geocentric latitude [-89 .. 89]; the poles are singular
         lonLaunch = -155;             %deg, longitude [-180 .. 180]
           hLaunch = 0;                %km,  pad altitude above the sphere [0 .. 5]
          azLaunch = 60;               %deg, launch azimuth clockwise from north [0 .. 360]
           vLaunch = 10;               %m/s, speed at the first integrated point [2 .. 50].
                                       %     NOT a physical launch speed: the 3DOF equations
                                       %     are singular at V = 0, so the integration starts
                                       %     a moment after first motion

%% Phase 1, boost -- commanded pitch ATTITUDE against time since liftoff. The
%% angle of attack follows as alpha = theta - gamma, so a schedule that tracks
%% the natural gravity turn keeps alpha small. The launch flight-path angle is
%% taken from the FIRST node, so the vehicle leaves the pad at zero incidence.
%% This program is DEPRESSED on purpose -- it pitches over hard and burns out
%% shallow and low, which is what turns a booster into a boost-glide first
%% stage rather than the lofted ballistic arc BM/run_ballistic flies:
         pitchTime = [0  6  15  30  50  70  82];   %s,   nodes, strictly increasing [0 .. burn time]
        pitchAngle = [89 80  55  28  12   5   2];  %deg, commanded attitude [0 .. 90], descending
          alphaMax = 10;               %deg, clamp on |angle of attack| [1 .. 15]; this clamp,
                                       %     not the schedule, is what limits how fast the
                                       %     flight path can be pushed over
         bankBoost = 0;                %deg, bank during boost [-90 .. 90]; 0 keeps the ascent
                                       %     in the launch plane
         tMaxBoost = 200;              %s,  boost horizon; must exceed the burn time

%% Staging -- whether the spent booster is thrown away at burnout:
        separation = true;             %true jettisons bst.massDry and flies the glide vehicle
                                       %     alone; false keeps the dead booster attached, which
                                       %     is a legitimate configuration and flies the whole
                                       %     unpowered flight on the STACK mass and stack
                                       %     aerodynamics -- and, that stack being a poor lifting
                                       %     shape, barely glides at all

%% Phase 2, glide -- prescribed control against time SINCE THE START OF THE
%% GLIDE, not since liftoff. Small bank keeps the lift vector up and buys
%% range; see the note in the header for what separates this from phase 3:
         glideTime = [0 6000];         %s,   nodes, strictly increasing, phase-local
        glideAlpha = [0 0];            %deg, angle of attack. NO EFFECT with the default constLD
                                       %     aero model, which ignores it; becomes live only when
                                       %     aeroFn below is swapped for an alpha-dependent model
         glideBank = [0 0];            %deg, bank [-90 .. 90]; 0 = all lift up, +/-90 = none up
          hHandoff = 15;               %km, glide-to-descent handoff, DESCENDING crossing only
                                       %    (0 .. hLaunch); must sit below every trough of the
                                       %    glide phugoid, see the header note
         tMaxGlide = 6000;             %s, glide horizon; raise it if the glide is cut short

%% Phase 3, descent -- the terminal dive. Same airframe, same equations; the
%% large bank angle is the whole difference, and it is what steepens the flight
%% path and drives dynamic pressure up on the way in:
          descTime = [0 1000];         %s,   nodes, strictly increasing, phase-local
         descAlpha = [0 0];            %deg, angle of attack; see glideAlpha
          descBank = [75 75];          %deg, bank [-90 .. 90]; only cos(bank) of the lift is
                                       %     left supporting the vehicle, so 75 deg keeps 26 %
                                       %     of it and the rest of the weight is unopposed
             hStop = 0;                %km, impact altitude, DESCENDING crossing only [0 .. 10]
          tMaxDesc = 1000;             %s, descent horizon

%% Vehicle and booster -- point these at any file returning the right struct:
         vehicleFn = @vehicle_hgv;                       %handle, see vehicle_hgv.m
         boosterFn = @coorbital.util.boosterDefaults;    %handle, returns the booster struct

%% Model selection -- swap a handle here to raise fidelity; each replacement
%% must keep the same signature as the one it replaces:
           atmosFn = @coorbital.atmos.expAtmos;   %[rho,P,T,a] = atmosFn(h)
            gravFn = @coorbital.grav.sphereGrav;  %[gr,gLat]   = gravFn(r,lat)
            aeroFn = @coorbital.aero.constLD;     %[CL,CD]     = aeroFn(alpha,mach,veh)
            propFn = @coorbital.prop.constThrust; %[T,mdot]    = propFn(t,P,veh)
         earthSpin = false;            %true enables the Coriolis and centrifugal terms

%% Output:
         showPlots = true;             %false to skip the figures, e.g. under matlab -batch
%% ======================= END USER PARAMETERS ============================

%% Apply the caller's overrides, if any. This is the ONE supported way to
%% change a run without editing the block above. Overrides are applied BEFORE
%% the unit conversion below, so a caller writes degrees and kilometres exactly
%% as the block does. Every entry in the block is overridable and nothing else
%% is, so a misspelt name raises rather than quietly doing nothing:
    if nargin < 1 || isempty(opts)
              opts = struct();
    end
    assert(isstruct(opts),'opts must be a struct of USER PARAMETERS overrides.');
       overridable = {'latLaunch','lonLaunch','hLaunch','azLaunch','vLaunch', ...
                      'pitchTime','pitchAngle','alphaMax','bankBoost','tMaxBoost', ...
                      'separation','glideTime','glideAlpha','glideBank', ...
                      'hHandoff','tMaxGlide','descTime','descAlpha','descBank', ...
                      'hStop','tMaxDesc','vehicleFn','boosterFn','atmosFn', ...
                      'gravFn','aeroFn','propFn','earthSpin','showPlots'};
             given = fieldnames(opts);
    for ko = 1:numel(given)
        assert(any(strcmp(given{ko},overridable)), ...
            '"%s" is not a USER PARAMETERS entry; there is nothing to override.', ...
            given{ko});
    end
         latLaunch = overrideOf(opts,'latLaunch',latLaunch);
         lonLaunch = overrideOf(opts,'lonLaunch',lonLaunch);
           hLaunch = overrideOf(opts,'hLaunch',hLaunch);
          azLaunch = overrideOf(opts,'azLaunch',azLaunch);
           vLaunch = overrideOf(opts,'vLaunch',vLaunch);
         pitchTime = overrideOf(opts,'pitchTime',pitchTime);
        pitchAngle = overrideOf(opts,'pitchAngle',pitchAngle);
          alphaMax = overrideOf(opts,'alphaMax',alphaMax);
         bankBoost = overrideOf(opts,'bankBoost',bankBoost);
         tMaxBoost = overrideOf(opts,'tMaxBoost',tMaxBoost);
        separation = overrideOf(opts,'separation',separation);
         glideTime = overrideOf(opts,'glideTime',glideTime);
        glideAlpha = overrideOf(opts,'glideAlpha',glideAlpha);
         glideBank = overrideOf(opts,'glideBank',glideBank);
          hHandoff = overrideOf(opts,'hHandoff',hHandoff);
         tMaxGlide = overrideOf(opts,'tMaxGlide',tMaxGlide);
          descTime = overrideOf(opts,'descTime',descTime);
         descAlpha = overrideOf(opts,'descAlpha',descAlpha);
          descBank = overrideOf(opts,'descBank',descBank);
             hStop = overrideOf(opts,'hStop',hStop);
          tMaxDesc = overrideOf(opts,'tMaxDesc',tMaxDesc);
         vehicleFn = overrideOf(opts,'vehicleFn',vehicleFn);
         boosterFn = overrideOf(opts,'boosterFn',boosterFn);
           atmosFn = overrideOf(opts,'atmosFn',atmosFn);
            gravFn = overrideOf(opts,'gravFn',gravFn);
            aeroFn = overrideOf(opts,'aeroFn',aeroFn);
            propFn = overrideOf(opts,'propFn',propFn);
         earthSpin = overrideOf(opts,'earthSpin',earthSpin);
         showPlots = overrideOf(opts,'showPlots',showPlots);

%% Convert the user block to library SI units. This is the ONLY unit
%% conversion in the file; everything past this point is m, m/s, rad and s:
          hLaunchM = hLaunch.*1000;
         hHandoffM = hHandoff.*1000;
            hStopM = hStop.*1000;
        latLaunchR = deg2rad(latLaunch);
        lonLaunchR = deg2rad(lonLaunch);
         azLaunchR = deg2rad(azLaunch);
         pitchTimS = pitchTime(:).';
         pitchAngR = deg2rad(pitchAngle(:).');
         alphaMaxR = deg2rad(alphaMax);
        bankBoostR = deg2rad(bankBoost);
         glideTimS = glideTime(:).';
         glideAlfR = deg2rad(glideAlpha(:).');
         glideBnkR = deg2rad(glideBank(:).');
          descTimS = descTime(:).';
          descAlfR = deg2rad(descAlpha(:).');
          descBnkR = deg2rad(descBank(:).');

%% Sanity-check the user block before spending time in the integrator:
    checkSchedule(pitchTimS,pitchAngR,'pitchTime','pitchAngle');
    checkSchedule(glideTimS,glideAlfR,'glideTime','glideAlpha');
    checkSchedule(glideTimS,glideBnkR,'glideTime','glideBank');
    checkSchedule(descTimS ,descAlfR ,'descTime' ,'descAlpha');
    checkSchedule(descTimS ,descBnkR ,'descTime' ,'descBank');
    assert(hHandoffM > hStopM, ...
        'hHandoff (%.1f km) must be above hStop (%.1f km).',hHandoff,hStop);
    assert(hStopM < hLaunchM + 1000, ...
        'hStop (%.1f km) is at or above the pad; the run would end at t = 0.',hStop);
    assert(vLaunch > 1, ...
        'vLaunch must exceed 1 m/s; the equations of motion are singular below that.');
    assert(all([tMaxBoost tMaxGlide tMaxDesc] > 0),'every phase horizon must be positive.');
    assert(islogical(separation) || isnumeric(separation), ...
        'separation must be true or false.');

                 c = coorbital.util.missileConst();
               veh = vehicleFn();
               bst = boosterFn();

%% Mass bookkeeping -- the state mass is ALWAYS the total mass carried; see
%% coorbital.util.boosterDefaults:
          mLiftoff = veh.mass + bst.massDry + bst.massProp;
         mBurnoutT = veh.mass + bst.massDry;
    if separation
            mGlide = veh.mass;
    else
            mGlide = mBurnoutT;
    end

%% Assemble the environment from the handles chosen above:
         env.atmos = atmosFn;
          env.grav = gravFn;
          env.aero = aeroFn;
          env.prop = propFn;
        env.omegaE = 0;
    if earthSpin
        env.omegaE = c.omegaE;
    end

%% The unpowered vehicle, flown by BOTH the glide and the descent. Its
%% aerodynamics are the glide vehicle's when the booster is jettisoned and the
%% STACK's when it is not, and its mass field must equal the mass actually
%% carried. Rebuilding the mass alone would leave Sref, CL and LD describing
%% the jettisoned stack, which is the worse bug because it looks repaired:
    if separation
          glideVeh = veh;
    else
          glideVeh = bst;
    end
     glideVeh.mass = mGlide;

%% Guidance, one schedule per phase, each written against ITS OWN phase clock:
        schedBoost = struct('tGrid',pitchTimS, ...
                            'theta',pitchAngR, ...
                            'sigma',bankBoostR.*ones(1,numel(pitchTimS)), ...
                            'alphaMax',alphaMaxR);
        schedGlide = struct('tGrid',glideTimS,'alpha',glideAlfR,'sigma',glideBnkR);
         schedDesc = struct('tGrid',descTimS ,'alpha',descAlfR ,'sigma',descBnkR);

%% Equations of motion, each with ITS OWN vehicle bound in. The vehArg argument
%% is the struct phaseRun forwards; it is deliberately not used, because a
%% single chain-wide vehicle cannot describe both the stack and the separated
%% glide vehicle:
           eomFree = coorbital.eom.massConstant(@coorbital.eom.glide3DOF);
            eomBst = @(t,x,u,vehArg,envArg) coorbital.eom.boost3DOF(t,x,u,bst,envArg);
            eomGld = @(t,x,u,vehArg,envArg) eomFree(t,x,u,glideVeh,envArg);

%% Phase 1, boost: burns to propellant exhaustion, then jettisons the spent
%% booster. The link is the staging discontinuity -- components 1 through 6
%% pass through untouched and component 7 drops by the dry mass. With
%% separation switched off the link is [], which phaseRun reads as identity:
         ph(1).eom = eomBst;
       ph(1).guide = @(t,x) coorbital.guide.pitchProgram(t,x,schedBoost);
   ph(1).terminate = @(t,x) coorbital.prop.eventBurnout(t,x,mBurnoutT);
       ph(1).tspan = [0 tMaxBoost];
        ph(1).link = [];
    if separation
        ph(1).link = @(x) [x(1:6); x(7) - bst.massDry];
    end

%% Phase 2, glide: unpowered, small bank, ends on the descending crossing of
%% the handoff altitude:
         ph(2).eom = eomGld;
       ph(2).guide = @(t,x) coorbital.guide.prescribed(t,x,schedGlide);
   ph(2).terminate = @(t,x) coorbital.prop.eventAltitude(t,x,hHandoffM);
       ph(2).tspan = [0 tMaxGlide];
        ph(2).link = [];

%% Phase 3, descent: same airframe and same equations as the glide, hard bank,
%% ends at impact. Nothing is jettisoned here, so there is no link:
         ph(3).eom = eomGld;
       ph(3).guide = @(t,x) coorbital.guide.prescribed(t,x,schedDesc);
   ph(3).terminate = @(t,x) coorbital.prop.eventAltitude(t,x,hStopM);
       ph(3).tspan = [0 tMaxDesc];
        ph(3).link = [];

%% Liftoff state. The pad flight-path angle is the first commanded pitch
%% attitude, so the angle of attack starts at zero:
                x0 = [c.rE + hLaunchM; ...
                      lonLaunchR; ...
                      latLaunchR; ...
                      vLaunch; ...
                      pitchAngR(1); ...
                      azLaunchR; ...
                      mLiftoff];

%% Propagate. The vehicle argument is a formality here -- every phase binds its
%% own, see the note in the header -- so the stack is passed to keep the call
%% honest about what is on the pad:
              traj = coorbital.prop.phaseRun(ph,x0,bst,env);

%% Phase boundaries. phaseRun records the boundary sample ONCE, carrying the
%% outgoing phase's control, so the last row of phase k is that phase's
%% terminal state on the NEAR side of any staging jump:
               kBO = find(traj.phaseIdx == 1,1,'last');
               kHO = find(traj.phaseIdx == 2,1,'last');
                nS = numel(traj.t);
    assert(~isempty(kBO) && ~isempty(kHO), ...
        'the chain did not produce all three phases; phases present: %s', ...
        mat2str(unique(traj.phaseIdx).'));

%% Derived quantities for the summary:
               hKm = (traj.x(:,1) - c.rE)./1000;
                 V = traj.x(:,4);
             massS = traj.x(:,7);
 [rho,Pamb,~,aSnd] = atmosFn(traj.x(:,1) - c.rE);
              qbar = 0.5.*rho.*V.^2;
              mach = V./aSnd;

%% Aerodynamic and thrust accelerations sample by sample, with the phase's own
%% vehicle and the mass actually carried at that instant, so a table-driven
%% aero model or a throttled motor reports correctly here and not just the
%% constant ones:
             aLift = zeros(nS,1);
             aDrag = zeros(nS,1);
             aThrV = zeros(nS,1);
             aThrN = zeros(nS,1);
    for k = 1:nS
        if traj.phaseIdx(k) == 1
              vehK = bst;
        else
              vehK = glideVeh;
        end
       [CLk,CDk]   = aeroFn(traj.u(k,1),mach(k),vehK);
          aLift(k) = qbar(k).*vehK.Sref.*CLk./massS(k);
          aDrag(k) = qbar(k).*vehK.Sref.*CDk./massS(k);
        if traj.phaseIdx(k) == 1
            [Tk,~] = propFn(traj.t(k),Pamb(k),bst);
          aThrV(k) = Tk.*cos(traj.u(k,1))./massS(k);
          aThrN(k) = Tk.*sin(traj.u(k,1))./massS(k);
        end
    end

%% Sensed specific force -- what an accelerometer on board would read. Gravity
%% is not sensed and is excluded. nAero drops the thrust, which is the quantity
%% meant by "deceleration" on an unpowered entry:
             nSens = sqrt((aThrV - aDrag).^2 + (aThrN + aLift).^2)./c.g0;
             nAero = sqrt(aLift.^2 + aDrag.^2)./c.g0;

%% Mass must be constant once the motor is out, and it must be the mass the
%% unpowered vehicle was built around, because the six-state glide equations
%% divide by glideVeh.mass rather than by the mass state. Two separate claims,
%% asserted separately because they can fail separately and have different
%% budgets. Constancy is EXACT -- massConstant returns dm/dt = 0 identically,
%% so every Runge-Kutta stage adds exactly zero. Agreement with mGlide is only
%% as good as the burnout event solve:
             isUnp = traj.phaseIdx >= 2;
             mSpan = max(massS(isUnp)) - min(massS(isUnp));
    assert(mSpan == 0, ...
        ['the unpowered mass varied by %.3e kg; dm/dt is identically zero ' ...
         'after burnout, so it must be bit-exactly constant'],mSpan);
    assert(abs(massS(end) - mGlide) < 1e-6, ...
        ['the unpowered flight is carrying %.9f kg while its vehicle was ' ...
         'built around %.9f kg; the aerodynamic accelerations are being ' ...
         'divided by the wrong mass'],massS(end),mGlide);

%% Masks and peaks, taken over the phase where each is meaningful:
             isBst = traj.phaseIdx == 1;
             isGld = traj.phaseIdx == 2;
             isDsc = traj.phaseIdx == 3;
     [qBstMax,kQB] = maxOver(qbar ,isBst);
     [nBstMax,kNB] = maxOver(nSens,isBst);
     [qGldMax,kQG] = maxOver(qbar ,isGld);
     [nGldMax,kNG] = maxOver(nAero,isGld);
     [qDscMax,kQD] = maxOver(qbar ,isDsc);
     [nDscMax,kND] = maxOver(nAero,isDsc);
     [nAerMax,kNA] = maxOver(nAero,isUnp);
       [hGlMax,kG] = maxOver(hKm  ,isGld);

%% Phase means, time-weighted rather than sample-weighted so that an adaptive
%% step clustering samples in the interesting part of a phase cannot bias them:
    [gamGldAvg,sinkGld,drgGldAvg] = phaseMeans(traj.t,traj.x,aDrag,isGld,c.rE);
    [gamDscAvg,sinkDsc,drgDscAvg] = phaseMeans(traj.t,traj.x,aDrag,isDsc,c.rE);

%% Great-circle ranges. Everything is measured on the impact sphere
%% r = rE + hStop, so the launch-to-impact range and the leg ranges share one
%% baseline:
                rI = c.rE + hStopM;
            angTot = coorbital.util.greatCircle(latLaunchR,lonLaunchR, ...
                                                traj.x(end,3),traj.x(end,2));
           rangeKm = rI.*angTot./1000;
             angBO = coorbital.util.greatCircle(latLaunchR,lonLaunchR, ...
                                                traj.x(kBO,3),traj.x(kBO,2));
          downBOKm = rI.*angBO./1000;
             angGl = coorbital.util.greatCircle(traj.x(kBO,3),traj.x(kBO,2), ...
                                                traj.x(kHO,3),traj.x(kHO,2));
          glideKm  = rI.*angGl./1000;
             angDs = coorbital.util.greatCircle(traj.x(kHO,3),traj.x(kHO,2), ...
                                                traj.x(end,3),traj.x(end,2));
           descKm  = rI.*angDs./1000;

%% ---------------------------------------------------------------------
%% What the descent phase actually buys, measured against the right baseline
%% ---------------------------------------------------------------------
%% Comparing the descent's peak dynamic pressure against the GLIDE's would be
%% meaningless, and the numbers say so: the glide's peak is set by the entry
%% pull-up, half an hour and five kilometres per second earlier, and no
%% terminal segment at Mach 2 is going to beat it. The comparison that means
%% something holds the state fixed and changes only the control -- the phase-3
%% initial state flown twice to the same stop altitude, once on the descent
%% schedule (the flown trajectory) and once on the GLIDE schedule continued,
%% which is exactly the trajectory the phase split replaces.
%%
%% These differences are MODEL OBSERVATIONS, reported and NOT asserted:
%% nothing in this file predicts their size. What theory does fix is their
%% direction, and the judgement printed below is passed on that:
              xHOP = traj.junction(2).x;
         phG       = ph(3);
         phG.guide = @(t,x) coorbital.guide.prescribed(t,x,schedGlide);
         phG.tspan = [0 tMaxGlide];
             trajG = coorbital.prop.phaseRun(phG,xHOP,glideVeh,env);
             cfrOK = abs((trajG.x(end,1) - c.rE) - hStopM) < 1e-3;
               nSG = numel(trajG.t);
   [rhoG,~,~,aSnG] = atmosFn(trajG.x(:,1) - c.rE);
              qbrG = 0.5.*rhoG.*trajG.x(:,4).^2;
             aDrgG = zeros(nSG,1);
    for k = 1:nSG
            [~,CDk] = aeroFn(trajG.u(k,1),trajG.x(k,4)./aSnG(k),glideVeh);
          aDrgG(k) = qbrG(k).*glideVeh.Sref.*CDk./glideVeh.mass;
    end
            angCfr = coorbital.util.greatCircle(xHOP(3),xHOP(2), ...
                                                trajG.x(end,3),trajG.x(end,2));
            cfrKm  = rI.*angCfr./1000;
              allG = true(nSG,1);
   [gamCfrAvg,sinkCfr,drgCfrAvg] = phaseMeans(trajG.t,trajG.x,aDrgG,allG,c.rE);

%% ---------------------------------------------------------------------
%% Termination diagnosis, one line per phase
%% ---------------------------------------------------------------------
%% An unexpected termination must SAY so. A boost that ran out of horizon
%% before it ran out of propellant, or a glide that never came down to the
%% handoff, both produce a short trajectory that otherwise reads as a completed
%% flight:
     [why1,ok1] = whyBurnout(traj,kBO,mBurnoutT,tMaxBoost);
     [why2,ok2] = whyAltitude(traj.t(kHO) - traj.t(kBO),traj.t(kHO), ...
                              traj.x(kHO,1) - c.rE,hHandoffM,hHandoff, ...
                              tMaxGlide,'handoff');
     [why3,ok3] = whyAltitude(traj.t(nS) - traj.t(kHO),traj.t(nS), ...
                              traj.x(nS,1) - c.rE,hStopM,hStop, ...
                              tMaxDesc,'impact');
          allOK = ok1 && ok2 && ok3;

%% Labels for the switches in the model line:
           spinTxt = 'OFF';
    if earthSpin
           spinTxt = 'ON';
    end
            sepTxt = 'booster JETTISONED at burnout';
    if ~separation
            sepTxt = 'booster RETAINED through impact (no separation)';
    end

%% Conventional lower edge of the hypersonic regime. A modelling convention,
%% not a physical constant of the Earth or the air, so it does not belong in
%% missileConst; it is the threshold below which holding CL and L/D constant
%% stops being defensible:
         machHyper = 5;

%% Report:
    fprintf('\n');
    fprintf('===== Boost-glide trajectory summary ====================\n');
    fprintf('  Phase termination\n');
    fprintf('    1 boost          %s\n',why1);
    fprintf('    2 glide          %s\n',why2);
    fprintf('    3 descent        %s\n',why3);
    if ~allOK
        fprintf('  *** CAUTION ***  at least one phase did NOT end as intended; the\n');
        fprintf('                   trajectory below is TRUNCATED, not a completed flight\n');
    end
    fprintf('\n');
    fprintf('  Launch\n');
    fprintf('    latitude         %10.4f deg\n',latLaunch);
    fprintf('    longitude        %10.4f deg\n',lonLaunch);
    fprintf('    pad altitude     %10.3f km\n',hLaunch);
    fprintf('    azimuth          %10.2f %-5s  (clockwise from north)\n',azLaunch,'deg');
    fprintf('    initial speed    %10.2f %-5s  (integration start, not a physical launch speed)\n', ...
            vLaunch,'m/s');
    fprintf('    liftoff mass     %10.1f %-5s  (payload %.1f + booster dry %.1f + propellant %.1f)\n', ...
            mLiftoff,'kg',veh.mass,bst.massDry,bst.massProp);
    fprintf('    liftoff T/W      %10.3f %-5s  (vacuum thrust over liftoff weight)\n', ...
            bst.thrustVac./(mLiftoff.*c.g0),'');
    fprintf('\n');
    fprintf('  Burnout and separation, end of phase 1\n');
    fprintf('    time             %10.3f s\n',traj.t(kBO));
    fprintf('    altitude         %10.3f km\n',hKm(kBO));
    fprintf('    speed            %10.2f %-5s  (Mach %.2f, planet-relative)\n', ...
            V(kBO),'m/s',mach(kBO));
    fprintf('    flight path      %10.3f %-5s  (positive is climbing; a DEPRESSED burnout,\n', ...
            rad2deg(traj.x(kBO,5)),'deg');
    fprintf('                                        which is what makes the glide possible)\n');
    fprintf('    heading          %10.3f deg\n',rad2deg(traj.x(kBO,6)));
    fprintf('    downrange        %10.2f %-5s  (great circle from the pad)\n',downBOKm,'km');
    fprintf('    mass at burnout  %10.1f %-5s  (payload + spent booster, before staging)\n', ...
            massS(kBO),'kg');
    fprintf('    mass into glide  %10.1f %-5s  (%s)\n',mGlide,'kg',sepTxt);
    fprintf('\n');
    fprintf('  Glide, phase 2, from separation to the handoff\n');
    fprintf('    duration         %10.2f %-5s  (%.2f min)\n', ...
            traj.t(kHO) - traj.t(kBO),'s',(traj.t(kHO) - traj.t(kBO))./60);
    fprintf('    ground covered   %10.2f %-5s  (great circle, burnout point to handoff point)\n', ...
            glideKm,'km');
    fprintf('    peak altitude    %10.2f %-5s  (at t = %.1f s; the first skip, thrown by the\n', ...
            hGlMax,'km',traj.t(kG));
    fprintf('                                        pull-up as the vehicle meets the air)\n');
    fprintf('    handoff time     %10.3f s\n',traj.t(kHO));
    fprintf('    handoff altitude %10.4f %-5s  (%+.2e m from the %.3f km target, event residual)\n', ...
            hKm(kHO),'km',traj.x(kHO,1) - c.rE - hHandoffM,hHandoff);
    fprintf('    handoff speed    %10.2f %-5s  (Mach %.2f)\n',V(kHO),'m/s',mach(kHO));
    fprintf('    handoff gamma    %10.3f %-5s  (negative is descending)\n', ...
            rad2deg(traj.x(kHO,5)),'deg');
    fprintf('\n');
    fprintf('  Descent, phase 3, from the handoff to impact\n');
    fprintf('    duration         %10.2f s\n',traj.t(nS) - traj.t(kHO));
    fprintf('    ground covered   %10.2f %-5s  (great circle, handoff point to impact)\n', ...
            descKm,'km');
    fprintf('    commanded bank   %10.2f %-5s  (cos = %.4f, the fraction of lift still\n', ...
            rad2deg(descBnkR(1)),'deg',cos(descBnkR(1)));
    fprintf('                                        supporting the vehicle)\n');
    fprintf('\n');
    fprintf('  What the descent phase buys, from the SAME handoff state\n');
    fprintf('    Glide and descent share an airframe and share the equations of motion; the\n');
    fprintf('    bank angle is the whole difference. So the phase-3 initial state is flown\n');
    fprintf('    twice to the %.1f km stop -- once on the %.0f deg descent bank, which is the\n', ...
            hStop,rad2deg(descBnkR(1)));
    fprintf('    trajectory above, and once on the %.0f deg GLIDE bank continued, which is the\n', ...
            rad2deg(glideBnkR(end)));
    fprintf('    trajectory the phase split replaces. Comparing the descent against the\n');
    fprintf('    glide PHASE instead would be meaningless: that phase peaked at the entry\n');
    fprintf('    pull-up, %.0f s and %.0f m/s earlier, and no Mach %.1f terminal segment beats it.\n', ...
            traj.t(kHO) - traj.t(kQG),V(kQG) - V(kHO),mach(kHO));
    if ~cfrOK
        fprintf('    *** the continued glide did NOT reach the stop altitude within its %.0f s\n',tMaxGlide);
        fprintf('        horizon; the column below is truncated and must not be read as a\n');
        fprintf('        completed comparison ***\n');
    end
    fprintf('                                 flown descent   glide continued\n');
    fprintf('    time to the stop      %14.2f  %14.2f   s\n', ...
            traj.t(nS) - traj.t(kHO),trajG.t(end));
    fprintf('    ground covered        %14.2f  %14.2f   km\n',descKm,cfrKm);
    fprintf('    mean flight path      %14.3f  %14.3f   deg (time-weighted)\n', ...
            gamDscAvg,gamCfrAvg);
    fprintf('    mean sink rate        %14.2f  %14.2f   m/s\n',sinkDsc,sinkCfr);
    fprintf('    speed at the stop     %14.2f  %14.2f   m/s\n',V(end),trajG.x(end,4));
    fprintf('    peak dynamic pressure %14.2f  %14.2f   kPa\n',qDscMax./1000,max(qbrG)./1000);
    fprintf('    mean drag deceleration%14.3f  %14.3f   g\n', ...
            drgDscAvg./c.g0,drgCfrAvg./c.g0);
    fprintf('    These differences are MODEL OBSERVATIONS, reported and NOT asserted; nothing\n');
    fprintf('    in this file predicts their size. Their DIRECTION is fixed by cos(sigma):\n');
    if cos(descBnkR(1)) < cos(glideBnkR(end))
        fprintf('    the descent keeps %.4f of the lift vector supporting the vehicle against\n', ...
                cos(descBnkR(1)));
        fprintf('    the glide''s %.4f, so it MUST fall steeper and arrive sooner, and it does.\n', ...
                cos(glideBnkR(end)));
    else
        fprintf('    the descent keeps %.4f of the lift vector supporting the vehicle against\n', ...
                cos(descBnkR(1)));
        fprintf('    the glide''s %.4f, so it is NOT the steeper segment. The phase split is not\n', ...
                cos(glideBnkR(end)));
        fprintf('    earning itself on this configuration; check descBank against glideBank.\n');
    end
    fprintf('\n');
    fprintf('  Phase means over the flown trajectory, time-weighted\n');
    fprintf('                                     glide         descent\n');
    fprintf('    mean flight path      %14.3f  %14.3f   deg\n',gamGldAvg,gamDscAvg);
    fprintf('    mean sink rate        %14.2f  %14.2f   m/s\n',sinkGld,sinkDsc);
    fprintf('    mean drag deceleration%14.3f  %14.3f   g\n', ...
            drgGldAvg./c.g0,drgDscAvg./c.g0);
    fprintf('\n');
    fprintf('  Overall\n');
    fprintf('    flight time      %10.2f %-5s  (%.2f min)\n',traj.t(end),'s',traj.t(end)./60);
    fprintf('    ground range     %10.2f %-5s  (great circle on the r = %.1f km impact sphere)\n', ...
            rangeKm,'km',rI./1000);
    fprintf('    central angle    %10.4f deg\n',rad2deg(angTot));
    fprintf('    leg sum          %10.2f %-5s  = %.2f (boost) + %.2f (glide) + %.2f (descent).\n', ...
            downBOKm + glideKm + descKm,'km',downBOKm,glideKm,descKm);
    fprintf('                                        The legs need NOT add to the total: each is\n');
    fprintf('                                        its own great circle, and a banked flight\n');
    fprintf('                                        does not stay on one.\n');
    fprintf('    samples          %10d %-5s  (ode45 adaptive steps over 3 phases)\n',nS,'');
    fprintf('\n');
    fprintf('  Impact\n');
    fprintf('    altitude         %10.6f %-5s  (%+.2e m from the %.3f km stop, event residual)\n', ...
            hKm(end),'km',traj.x(end,1) - c.rE - hStopM,hStop);
    fprintf('    speed            %10.2f %-5s  (Mach %.2f)\n',V(end),'m/s',mach(end));
    fprintf('    flight path      %10.3f %-5s  (negative is descending)\n', ...
            rad2deg(traj.x(end,5)),'deg');
    fprintf('    heading          %10.3f deg\n',rad2deg(traj.x(end,6)));
    fprintf('    latitude         %10.4f deg\n',rad2deg(traj.x(end,3)));
    fprintf('    longitude        %10.4f deg\n',rad2deg(traj.x(end,2)));
    fprintf('\n');
    fprintf('  Peak loads\n');
    fprintf('    boost max q      %10.2f %-5s  (at t = %.1f s, h = %.2f km)\n', ...
            qBstMax./1000,'kPa',traj.t(kQB),hKm(kQB));
    fprintf('    boost sensed     %10.2f %-5s  (%.2f m/s^2, thrust included, at t = %.1f s)\n', ...
            nBstMax,'g',nBstMax.*c.g0,traj.t(kNB));
    fprintf('    glide max q      %10.2f %-5s  (at t = %.1f s, h = %.2f km)\n', ...
            qGldMax./1000,'kPa',traj.t(kQG),hKm(kQG));
    fprintf('    glide decel      %10.2f %-5s  (%.2f m/s^2 aero only, at t = %.1f s, h = %.2f km)\n', ...
            nGldMax,'g',nGldMax.*c.g0,traj.t(kNG),hKm(kNG));
    fprintf('    descent max q    %10.2f %-5s  (at t = %.1f s, h = %.2f km)\n', ...
            qDscMax./1000,'kPa',traj.t(kQD),hKm(kQD));
    fprintf('    descent decel    %10.2f %-5s  (%.2f m/s^2 aero only, at t = %.1f s, h = %.2f km)\n', ...
            nDscMax,'g',nDscMax.*c.g0,traj.t(kND),hKm(kND));
    fprintf('    peak decel       %10.2f %-5s  (%.2f m/s^2, over the whole unpowered flight,\n', ...
            nAerMax,'g',nAerMax.*c.g0);
    fprintf('                                        at t = %.1f s, h = %.2f km, in phase %d)\n', ...
            traj.t(kNA),hKm(kNA),traj.phaseIdx(kNA));
    fprintf('\n');
    fprintf('  Models: atmos %s | grav %s | aero %s\n', ...
            func2str(atmosFn),func2str(gravFn),func2str(aeroFn));
    fprintf('          prop %s | Earth rotation %s\n',func2str(propFn),spinTxt);
    fprintf('  Vehicle: %s, mass %.1f kg, Sref %.3f m^2, CL %.3f, L/D %.2f (PLACEHOLDER values)\n', ...
            func2str(vehicleFn),veh.mass,veh.Sref,veh.CL,veh.LD);
    fprintf('  Booster: %s, dry %.1f kg, prop %.1f kg, Tvac %.1f kN, Isp %.1f s (PLACEHOLDER values)\n', ...
            func2str(boosterFn),bst.massDry,bst.massProp,bst.thrustVac./1000,bst.Isp);
    fprintf('  Validity: %s holds CL and L/D constant, a hypersonic approximation,\n', ...
            func2str(aeroFn));
    fprintf('            and this flight runs from a subsonic liftoff through Mach %.0f and\n', ...
            max(mach));
    fprintf('            back down again, so the transonic ends of the boost and of the\n');
    fprintf('            terminal descent lie outside it. The isothermal %s and\n',func2str(atmosFn));
    fprintf('            the spherical %s carry no J2 and no winds.\n',func2str(gravFn));
    fprintf('            Treat this as an indicative trajectory, not a performance prediction.\n');
    if mach(kHO) < machHyper
        fprintf('            The handoff is at Mach %.2f, below Mach %.0f, so the ENTIRE\n', ...
                mach(kHO),machHyper);
        fprintf('            descent phase is outside the constant-coefficient regime. It is\n');
        fprintf('            reported to show what the phase split does, not as a terminal\n');
        fprintf('            performance figure.\n');
    end
    if nBstMax > 15
        fprintf('            The %.1f g sensed load at end of burn is an artefact of the\n',nBstMax);
        fprintf('            constant-thrust placeholder motor, which does not throttle as the\n');
        fprintf('            stack empties. A real stage would tail off or stage before that.\n');
    end
    fprintf('=========================================================\n\n');

%% Hand the summary's numbers back at full precision, so a test does not have
%% to read them out of printed text:
      info.tBurnout  = traj.t(kBO);
      info.hBurnout  = traj.x(kBO,1) - c.rE;
      info.vBurnout  = V(kBO);
    info.gamBurnout  = traj.x(kBO,5);
    info.psiBurnout  = traj.x(kBO,6);
      info.mBurnout  = massS(kBO);
        info.mGlide  = mGlide;
      info.mLiftoff  = mLiftoff;
      info.downBOKm  = downBOKm;
      info.tHandoff  = traj.t(kHO);
      info.hHandoff  = traj.x(kHO,1) - c.rE;
      info.vHandoff  = V(kHO);
    info.gamHandoff  = traj.x(kHO,5);
      info.hGlideMax = hGlMax;
       info.glideKm  = glideKm;
        info.descKm  = descKm;
       info.tFlight  = traj.t(end);
       info.rangeKm  = rangeKm;
        info.angTot  = angTot;
       info.vImpact  = V(end);
     info.gamImpact  = traj.x(end,5);
     info.latImpact  = traj.x(end,3);
     info.lonImpact  = traj.x(end,2);
      info.machImp   = mach(end);
       info.qBstMax  = qBstMax;
       info.nBstMax  = nBstMax;
       info.qGldMax  = qGldMax;
       info.nGldMax  = nGldMax;
       info.qDscMax  = qDscMax;
       info.nDscMax  = nDscMax;
       info.nAerMax  = nAerMax;
      info.tNAerMax  = traj.t(kNA);
    info.gamGldAvg   = gamGldAvg;
    info.gamDscAvg   = gamDscAvg;
       info.sinkGld  = sinkGld;
       info.sinkDsc  = sinkDsc;
     info.drgGldAvg  = drgGldAvg;
     info.drgDscAvg  = drgDscAvg;
        info.cfrOK   = cfrOK;
        info.tCfr    = trajG.t(end);
       info.cfrKm    = cfrKm;
        info.vCfr    = trajG.x(end,4);
      info.qCfrMax   = max(qbrG);
    info.gamCfrAvg   = gamCfrAvg;
       info.sinkCfr  = sinkCfr;
     info.drgCfrAvg  = drgCfrAvg;
        info.stopOK  = allOK;
       info.stopWhy  = {why1;why2;why3};

%% ...and the machinery itself, so an independent checker can re-integrate the
%% very same chain with a different solver instead of trusting this one:
        info.phases  = ph;
           info.env  = env;
            info.x0  = x0;
      info.boostVeh  = bst;
      info.glideVeh  = glideVeh;

%% Plots:
    if showPlots
        figure('color',[1 1 1],'name','Boost-glide time histories');
        subplot(3,1,1); plot(traj.t,hKm,'linewidth',1.5); grid on;
        ylabel('altitude (km)'); title('Boost, glide and terminal descent');
        subplot(3,1,2); plot(traj.t,V,'linewidth',1.5); grid on;
        ylabel('speed (m/s)');
        subplot(3,1,3); plot(traj.t,massS,'linewidth',1.5); grid on;
        ylabel('mass (kg)'); xlabel('time (s)');

        figure('color',[1 1 1],'name','Ground track');
        plot(rad2deg(traj.x(:,2)),rad2deg(traj.x(:,3)),'linewidth',1.5); grid on; hold on;
        plot(lonLaunch,latLaunch,'o','markersize',8,'linewidth',1.5);
        plot(rad2deg(traj.x(kHO,2)),rad2deg(traj.x(kHO,3)),'d','markersize',8,'linewidth',1.5);
        plot(rad2deg(traj.x(end,2)),rad2deg(traj.x(end,3)),'s','markersize',8,'linewidth',1.5);
        xlabel('longitude (deg)'); ylabel('latitude (deg)');
        legend('ground track','launch','handoff','impact','location','best');
        title(sprintf('Ground track, %.0f km great-circle range',rangeKm));

        figure('color',[1 1 1],'name','Loads and flight path');
        subplot(3,1,1); plot(traj.t,qbar./1000,'linewidth',1.5); grid on;
        ylabel('q (kPa)'); title('Dynamic pressure, load and flight path angle');
        subplot(3,1,2); plot(traj.t,nAero,'linewidth',1.5); grid on;
        ylabel('aero load (g)');
        subplot(3,1,3); plot(traj.t,rad2deg(traj.x(:,5)),'linewidth',1.5); grid on;
        ylabel('gamma (deg)'); xlabel('time (s)');
    end

%% Hand the trajectory back only when the caller asked for it. Typing
%% "run_boost_glide" at the prompt should leave the summary on screen, not bury
%% it under a dump of the whole struct as ans:
    if nargout == 0
        clear traj info;
    end
end

function checkSchedule(tGrid,vGrid,tName,vName)
%% Purpose:
%
%  Reject a malformed control schedule before it reaches the integrator, where
%  interp1 would either error deep inside a derivative evaluation or, worse,
%  silently extrapolate off a non-monotonic grid.
%
%% Inputs:
%
%  tGrid            [1 x K]                     Schedule node times (s)
%
%  vGrid            [1 x K]                     Scheduled values (rad)
%
%  tName            Char [1 x n]                Name of the time entry in the
%                                               USER PARAMETERS block, used in
%                                               the failure message
%
%  vName            Char [1 x m]                Name of the value entry
%
%% Outputs:
%
%  none                                         Throws on a malformed schedule
%
%% Revision History:
%  Michael Casey                                                08/07/2026
%  Copyright 2026 Coorbital, Inc.
%% ------------------------ Begin Code Sequence ---------------------------

    assert(numel(tGrid) == numel(vGrid), ...
        '%s has %d nodes and %s has %d; they must match.', ...
        tName,numel(tGrid),vName,numel(vGrid));
    assert(numel(tGrid) >= 2,'%s needs at least two nodes.',tName);
    assert(all(diff(tGrid) > 0),'%s must be strictly increasing.',tName);
end

function [gamAvg,sinkAvg,drgAvg] = phaseMeans(t,x,aDrag,mask,rE)
%% Purpose:
%
%  Time-weighted means of the quantities that separate a glide from a terminal
%  descent, over one phase. Time-weighted rather than sample-weighted because
%  an adaptive integrator clusters samples where the dynamics are stiff, and a
%  plain mean over samples would report the interesting part of a phase as
%  though it lasted longer than it did.
%
%  The sink rate is the phase's altitude drop over its duration, which is the
%  average vertical speed with the sign flipped -- positive means descending --
%  and is the right average even when the phase climbs part way, as a skipping
%  glide does.
%
%% Inputs:
%
%  t                [N x 1]                     Cumulative time (s)
%
%  x                [N x 7]                     State history; uses x(:,1) = r
%                                               (m) and x(:,5) = gamma (rad)
%
%  aDrag            [N x 1]                     Drag acceleration (m/s^2)
%
%  mask             [N x 1] logical              Samples belonging to the phase;
%                                               at least two must be true
%
%  rE               [1 x 1]                     Planet radius (m)
%
%% Outputs:
%
%  gamAvg           [1 x 1]                     Mean flight path angle (deg)
%
%  sinkAvg          [1 x 1]                     Mean sink rate (m/s), positive
%                                               descending
%
%  drgAvg           [1 x 1]                     Mean drag acceleration (m/s^2)
%
%% Revision History:
%  Michael Casey                                                08/07/2026
%  Copyright 2026 Coorbital, Inc.
%% ------------------------ Begin Code Sequence ---------------------------

               idx = find(mask);
    assert(numel(idx) >= 2,'a phase mean needs at least two samples.');
                tp = t(idx);
              span = tp(end) - tp(1);
    assert(span > 0,'a phase of zero duration has no time-weighted mean.');
            gamAvg = rad2deg(trapz(tp,x(idx,5))./span);
           sinkAvg = ((x(idx(1),1) - rE) - (x(idx(end),1) - rE))./span;
            drgAvg = trapz(tp,aDrag(idx))./span;
end

function [pk,kAt] = maxOver(v,mask)
%% Purpose:
%
%  Largest value of a history within a masked subset of its samples, and the
%  index into the FULL history at which it occurs. Used so that a boost peak, a
%  glide peak and a descent peak can be reported separately without any search
%  wandering into another phase.
%
%% Inputs:
%
%  v                [N x 1]                     History to search
%
%  mask             [N x 1] logical             Samples to consider; must have
%                                               at least one true entry
%
%% Outputs:
%
%  pk               [1 x 1]                     Largest masked value
%
%  kAt              [1 x 1]                     Index into v, not into the
%                                               masked subset
%
%% Revision History:
%  Michael Casey                                                08/07/2026
%  Copyright 2026 Coorbital, Inc.
%% ------------------------ Begin Code Sequence ---------------------------

               idx = find(mask);
    assert(~isempty(idx),'the mask selected no samples to search.');
        [pk,kLoc]  = max(v(idx));
               kAt = idx(kLoc);
end

function [why,ok] = whyBurnout(traj,kBO,mBurnoutT,tMaxBoost)
%% Purpose:
%
%  Diagnose why the boost phase stopped. Nominal termination is propellant
%  exhaustion, which coorbital.prop.eventBurnout detects as the mass state
%  descending through the total burnout mass. Anything else -- a horizon
%  timeout, an integrator failure -- must be named rather than left to look
%  like a completed burn.
%
%% Inputs:
%
%  traj             Struct                      Trajectory from phaseRun
%
%  kBO              [1 x 1]                     Index of the last phase-1 sample
%
%  mBurnoutT        [1 x 1]                     Total burnout mass (kg)
%
%  tMaxBoost        [1 x 1]                     Boost horizon (s)
%
%% Outputs:
%
%  why              Char [1 x n]                Human-readable reason
%
%  ok               [1 x 1] logical             True only for a nominal burnout
%
%% Revision History:
%  Michael Casey                                                08/07/2026
%  Copyright 2026 Coorbital, Inc.
%% ------------------------ Begin Code Sequence ---------------------------

              tEnd = traj.t(kBO);
              mEnd = traj.x(kBO,7);
                ok = false;
    if abs(mEnd - mBurnoutT) < 1e-6
               why = sprintf('propellant exhausted at t = %.3f s, m = %.1f kg (nominal)', ...
                             tEnd,mEnd);
                ok = true;
    elseif tEnd >= tMaxBoost - 1e-6
               why = sprintf(['hit the %.0f s boost horizon with %.1f kg of ' ...
                              'propellant still aboard'],tMaxBoost,mEnd - mBurnoutT);
    else
               why = sprintf(['stopped early at t = %.3f s with m = %.1f kg against a ' ...
                              '%.1f kg burnout mass (integrator failure or an ' ...
                              'unmodelled event)'],tEnd,mEnd,mBurnoutT);
    end
end

function [why,ok] = whyAltitude(tPhase,tEnd,hEnd,hTargetM,hTargetKm,tMax,what)
%% Purpose:
%
%  Diagnose why an altitude-terminated phase stopped -- the glide at its
%  handoff, the descent at impact. Nominal termination is the descending
%  crossing of the target altitude. A phase that ran out of horizon in mid-air
%  is the failure this exists to name, and for the glide it is the likely one:
%  a handoff set above a phugoid trough is reached early, and a handoff set
%  below the terminal glide never arrives at all.
%
%% Inputs:
%
%  tPhase           [1 x 1]                     Duration of the phase (s)
%
%  tEnd             [1 x 1]                     Cumulative time at the end of
%                                               the phase (s)
%
%  hEnd             [1 x 1]                     Altitude at the end of the
%                                               phase (m)
%
%  hTargetM         [1 x 1]                     Target altitude (m)
%
%  hTargetKm        [1 x 1]                     Target altitude (km), for the
%                                               message
%
%  tMax             [1 x 1]                     Phase horizon (s)
%
%  what             Char [1 x n]                Name of the event, e.g.
%                                               'handoff' or 'impact'
%
%% Outputs:
%
%  why              Char [1 x n]                Human-readable reason
%
%  ok               [1 x 1] logical             True only for a nominal stop
%
%% Revision History:
%  Michael Casey                                                08/07/2026
%  Copyright 2026 Coorbital, Inc.
%% ------------------------ Begin Code Sequence ---------------------------

                ok = false;
    if abs(hEnd - hTargetM) < 1e-3
               why = sprintf(['reached the %.1f km %s altitude at t = %.3f s ' ...
                              '(nominal)'],hTargetKm,what,tEnd);
                ok = true;
    elseif tPhase >= tMax - 1e-6
               why = sprintf(['hit the %.0f s horizon at %.3f km, still flying; the ' ...
                              '%s altitude was never reached'],tMax,hEnd./1000,what);
    else
               why = sprintf(['stopped early at t = %.3f s, h = %.3f km ' ...
                              '(integrator failure or an unmodelled event)'], ...
                             tEnd,hEnd./1000);
    end
end

function v = overrideOf(opts,name,v)
%% Purpose:
%
%  Return the caller's override for one USER PARAMETERS entry when it was
%  supplied, and the script's own value otherwise. Factored out so the override
%  block above reads as one unbranched line per parameter.
%
%% Inputs:
%
%  opts             Struct                      Caller overrides, already
%                                               checked for unknown fields
%
%  name             Char [1 x n]                Name of the USER PARAMETERS
%                                               entry being resolved
%
%  v                Scalar / vector / handle    The script's own value, in the
%                                               user block's own units
%
%% Outputs:
%
%  v                Scalar / vector / handle    opts.(name) when that field is
%                                               present, otherwise the input v
%                                               unchanged
%
%% Revision History:
%  Michael Casey                                                08/07/2026
%  Copyright 2026 Coorbital, Inc.
%% ------------------------ Begin Code Sequence ---------------------------

    if isfield(opts,name)
                 v = opts.(name);
    end
end
