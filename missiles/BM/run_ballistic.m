function [traj,info] = run_ballistic(opts)
%% Purpose:
%
%  Fly a ballistic missile from the pad to impact and report what it did:
%  a powered boost on a prescribed pitch program, stage separation at
%  burnout, an unpowered coast to apogee, and a descent to impact. Everything
%  a routine run needs to change lives in the USER PARAMETERS block; nothing
%  below it should require editing.
%
%  Angles and distances in the USER PARAMETERS block are in degrees and
%  kilometres because that is how a user thinks about them. They are converted
%  to the library's SI units (m, m/s, rad, s) immediately after the block, in
%  one place, and nothing below that point works in any other unit.
%
%% Note -- why THREE phases, and why apogee splits the coast:
%
%  The flight is powered then unpowered, so two phases would do. It is run as
%  three because coorbital.prop.eventApogee terminates on the flight-path
%  angle crossing zero DESCENDING, which makes apogee an exactly solved event
%  rather than the largest sample on an adaptive grid -- apogee altitude and
%  time are then good to the event solver's tolerance, not to a step size.
%  The split also leaves phase 3 monotonically descending, so the one-sided
%  coorbital.prop.eventAltitude that ends the flight cannot be reached from
%  the wrong side.
%
%  eventApogee is NOT used to terminate a coast that has to continue past
%  apogee to impact; that is precisely the misuse it invites. It is used here
%  only because the coast is deliberately cut in two at that point.
%
%% Note -- two vehicles, one chain:
%
%  coorbital.prop.phaseRun carries ONE vehicle struct for the whole chain, but
%  the boosted stack and the separated re-entry body have different reference
%  areas, different aerodynamic coefficients and different masses. Each phase
%  therefore BINDS its own vehicle in the equation-of-motion closure below,
%  and the struct handed to phaseRun is never read by the equations of motion.
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
%                                               printed text. Fields are named
%                                               in the assembly block near the
%                                               end of this file; all SI
%                                               except the *Km distances.
%
%% References:
%   [1] Vinh, N.X., Busemann, A., Culp, R.D., "Hypersonic and Planetary Entry
%       Flight Mechanics," Univ. Michigan Press, 1980.
%   [2] Bate, R.R., Mueller, D.D., White, J.E., "Fundamentals of
%       Astrodynamics," Dover, 1971, Ch. 6 (ballistic missile trajectories:
%       the free-flight range angle of a Keplerian arc).
%
%% Revision History:
%  Michael Casey                                                08/07/2026
%  Michael Casey  alphaMax read from the vehicle; drag and     08/08/2026
%                 lift attributed on arcs that share one
%                 gravity model and one rotation rate
%  Copyright 2026 Coorbital, Inc.
%% ------------------------ Begin Code Sequence ---------------------------

%% Resolve paths so the script runs from anywhere:
              here = fileparts(mfilename('fullpath'));
    addpath(here);
    addpath(fullfile(here,'..'));

%% ========================= USER PARAMETERS ==============================
%% Launch site -- where the missile leaves the pad:
          latLaunch = 45;              %deg, geocentric latitude [-89 .. 89]; the poles are singular
          lonLaunch = -100;            %deg, longitude [-180 .. 180]
            hLaunch = 0;               %km,  pad altitude above the sphere [0 .. 5]
           azLaunch = 35;              %deg, launch azimuth clockwise from north [0 .. 360]
            vLaunch = 10;              %m/s, speed at the first integrated point [2 .. 50].
                                       %     NOT a physical launch speed: the 3DOF equations
                                       %     are singular at V = 0, so the integration starts
                                       %     a moment after first motion

%% Pitch program -- commanded pitch ATTITUDE against time since liftoff. The
%% angle of attack follows as alpha = theta - gamma, so a schedule that tracks
%% the natural gravity turn keeps alpha small. The launch flight-path angle is
%% taken from the FIRST node, so the vehicle leaves the pad at zero incidence:
          pitchTime = [0  6  15  30  50  70  82];   %s,   nodes, strictly increasing [0 .. burn time]
         pitchAngle = [89 86  76  60  46  38  34];  %deg, commanded attitude [5 .. 90], descending
           alphaMax = [];              %deg, clamp on |angle of attack|. EMPTY -- the
                                       %     shipped value -- reads the VEHICLE'S OWN
                                       %     control-authority limit, veh.alphaMaxDeg,
                                       %     which is 6 deg in BM/vehicle_bm.m. That is
                                       %     where a structural, actuator and
                                       %     aerothermal limit belongs, and it is why
                                       %     this script and BM/run_ballistic_target
                                       %     can no longer disagree about it. Put a
                                       %     number here [1 .. 15] to OVERRIDE the
                                       %     vehicle for a deliberate sensitivity
                                       %     study; the 6 deg is a PLACEHOLDER awaiting
                                       %     a qualification basis, not a cleared value
          bankAngle = 0;               %deg, bank [-90 .. 90]; 0 keeps the flight in the launch plane

%% Staging -- whether the spent booster is thrown away at burnout:
         separation = true;            %true jettisons bst.massDry and flies the re-entry body
                                       %     alone; false keeps the dead booster attached, which
                                       %     is a legitimate configuration and flies the whole
                                       %     coast on the STACK mass and stack aerodynamics

%% Termination -- when to stop integrating each phase:
              hStop = 0;               %km, impact altitude, DESCENDING crossing only [0 .. 20]
          tMaxBoost = 200;             %s,  boost horizon; must exceed the burn time
          tMaxCoast = 3000;            %s,  ascending-coast horizon, liftoff to apogee
           tMaxDesc = 3000;            %s,  descending horizon, apogee to impact

%% Vehicle and booster -- point these at any file returning the right struct:
          vehicleFn = @vehicle_bm;                        %handle, see vehicle_bm.m
          boosterFn = @coorbital.util.boosterDefaults;    %handle, returns the booster struct

%% Model selection -- swap a handle here to raise fidelity; each replacement
%% must keep the same signature as the one it replaces.
%%
%% ONE EXCEPTION, and it applies to gravFn only: the Keplerian cross-check
%% reported below always re-propagates on coorbital.grav.sphereGrav, whatever
%% is selected here, because the closed form it is compared against is
%% two-body by construction. Select an oblate model and the FLOWN trajectory
%% is oblate while that single cross-check stays Keplerian -- so it keeps
%% validating the propagator instead of aborting on the J2 signal:
            atmosFn = @coorbital.atmos.expAtmos;   %[rho,P,T,a] = atmosFn(h)
             gravFn = @coorbital.grav.sphereGrav;  %[gr,gLat]   = gravFn(r,lat)
             aeroFn = @coorbital.aero.constLD;     %[CL,CD]     = aeroFn(alpha,mach,veh)
             propFn = @coorbital.prop.constThrust; %[T,mdot]    = propFn(t,P,veh)
          earthSpin = false;           %true enables the Coriolis and centrifugal terms

%% Output:
          showPlots = true;            %false to skip the figures, e.g. under matlab -batch
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
                      'pitchTime','pitchAngle','alphaMax','bankAngle', ...
                      'separation','hStop','tMaxBoost','tMaxCoast','tMaxDesc', ...
                      'vehicleFn','boosterFn','atmosFn','gravFn','aeroFn', ...
                      'propFn','earthSpin','showPlots'};
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
         bankAngle = overrideOf(opts,'bankAngle',bankAngle);
        separation = overrideOf(opts,'separation',separation);
             hStop = overrideOf(opts,'hStop',hStop);
         tMaxBoost = overrideOf(opts,'tMaxBoost',tMaxBoost);
         tMaxCoast = overrideOf(opts,'tMaxCoast',tMaxCoast);
          tMaxDesc = overrideOf(opts,'tMaxDesc',tMaxDesc);
         vehicleFn = overrideOf(opts,'vehicleFn',vehicleFn);
         boosterFn = overrideOf(opts,'boosterFn',boosterFn);
           atmosFn = overrideOf(opts,'atmosFn',atmosFn);
            gravFn = overrideOf(opts,'gravFn',gravFn);
            aeroFn = overrideOf(opts,'aeroFn',aeroFn);
            propFn = overrideOf(opts,'propFn',propFn);
         earthSpin = overrideOf(opts,'earthSpin',earthSpin);
         showPlots = overrideOf(opts,'showPlots',showPlots);

%% The constants, the vehicle and the booster, read HERE rather than after the
%% sanity checks, because the ANGLE-OF-ATTACK CLAMP is a property of the
%% vehicle and the unit conversion below has to convert it:
                 c = coorbital.util.missileConst();
               veh = vehicleFn();
               bst = boosterFn();

%% Resolve the angle-of-attack clamp. It is a CONTROL-AUTHORITY LIMIT of the
%% airframe, not a targeting knob, so its home is BM/vehicle_bm.m and an empty
%% user-block entry reads it from there. This script and
%% BM/run_ballistic_target used to carry two different literals -- 6 deg and
%% 12 deg -- for nominally the same vehicle, which made their performance
%% non-comparable; a number in the block still overrides it, deliberately and
%% visibly, for a sensitivity study:
    assert(isfield(veh,'alphaMaxDeg'), ...
        ['%s returned a vehicle with no alphaMaxDeg field. The clamp on the ' ...
         'angle of attack is a vehicle limit and this script reads it from ' ...
         'the vehicle; set alphaMax in the USER PARAMETERS block to fly a ' ...
         'vehicle file that does not carry one.'],func2str(vehicleFn));
    if isempty(alphaMax)
          alphaMax = veh.alphaMaxDeg;
    end
    assert(isscalar(alphaMax) && isnumeric(alphaMax) && isfinite(alphaMax) && ...
           alphaMax > 0, ...
        'alphaMax must be a positive finite scalar in degrees; got %s.', ...
        mat2str(alphaMax));

%% Convert the user block to library SI units. This is the ONLY unit
%% conversion in the file; everything past this point is m, m/s, rad and s:
          hLaunchM = hLaunch.*1000;
            hStopM = hStop.*1000;
        latLaunchR = deg2rad(latLaunch);
        lonLaunchR = deg2rad(lonLaunch);
         azLaunchR = deg2rad(azLaunch);
        pitchAngR  = deg2rad(pitchAngle(:).');
         alphaMaxR = deg2rad(alphaMax);
           bankRad = deg2rad(bankAngle);
         pitchTimS = pitchTime(:).';

%% Sanity-check the user block before spending time in the integrator:
    assert(numel(pitchTimS) == numel(pitchAngR), ...
        'pitchTime has %d nodes and pitchAngle has %d; they must match.', ...
        numel(pitchTimS),numel(pitchAngR));
    assert(numel(pitchTimS) >= 2,'the pitch program needs at least two nodes.');
    assert(all(diff(pitchTimS) > 0),'pitchTime must be strictly increasing.');
    assert(hStopM < hLaunchM + 1000, ...
        'hStop (%.1f km) is at or above the pad; the run would end at t = 0.',hStop);
    assert(vLaunch > 1, ...
        'vLaunch must exceed 1 m/s; the equations of motion are singular below that.');
    assert(all([tMaxBoost tMaxCoast tMaxDesc] > 0),'every phase horizon must be positive.');
    assert(islogical(separation) || isnumeric(separation), ...
        'separation must be true or false.');

%% Mass bookkeeping -- the state mass is ALWAYS the total mass carried; see
%% coorbital.util.boosterDefaults:
           mLiftoff = veh.mass + bst.massDry + bst.massProp;
          mBurnoutT = veh.mass + bst.massDry;
    if separation
             mCoast = veh.mass;
    else
             mCoast = mBurnoutT;
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

%% The coast vehicle. Its aerodynamics are the re-entry body's when the
%% booster is jettisoned and the STACK's when it is not, and its mass field
%% must equal the mass actually carried.
%%
%% READ THIS BEFORE CHANGING ANYTHING ABOUT SEPARATION. The unpowered phases
%% run coorbital.eom.glide3DOF under coorbital.eom.massConstant, and
%% glide3DOF divides by veh.mass, not by the mass STATE. So x(7) is inert
%% after burnout: it is carried for bookkeeping and read by nothing. Deleting
%% the separation link therefore leaves the flown trajectory BIT-IDENTICAL and
%% changes only the recorded mass -- measured, not assumed. What actually
%% makes separation physical is the branch below, which rebuilds the coast
%% vehicle's mass AND its aerodynamics from the flag.
%%
%% The two representations of the same mass are tied together by the assertion
%% after the propagation, and by nothing else. That assertion is not
%% decoration; it is the only thing standing between a dropped link and a
%% coast flown at the wrong weight:
    if separation
           coastVeh = veh;
    else
           coastVeh = bst;
    end
   coastVeh.mass    = mCoast;

%% Guidance. Boost flies the commanded pitch attitude; the unpowered phases
%% fly at zero incidence, which for the near-symmetric re-entry body is what
%% "ballistic" means:
         schedBoost = struct('tGrid',pitchTimS, ...
                             'theta',pitchAngR, ...
                             'sigma',bankRad.*ones(1,numel(pitchTimS)), ...
                             'alphaMax',alphaMaxR);
         schedCoast = struct('tGrid',[0 max(tMaxCoast,tMaxDesc)], ...
                             'alpha',[0 0], ...
                             'sigma',[bankRad bankRad]);

%% Equations of motion, each with ITS OWN vehicle bound in. The vehArg
%% argument is the struct phaseRun forwards; it is deliberately not used,
%% because a single chain-wide vehicle cannot describe both the stack and the
%% separated body:
            eomFree = coorbital.eom.massConstant(@coorbital.eom.glide3DOF);
            eomBst  = @(t,x,u,vehArg,envArg) coorbital.eom.boost3DOF(t,x,u,bst,envArg);
            eomCst  = @(t,x,u,vehArg,envArg) eomFree(t,x,u,coastVeh,envArg);

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

%% Phase 2, ascending coast: unpowered, ends exactly at apogee:
          ph(2).eom = eomCst;
        ph(2).guide = @(t,x) coorbital.guide.prescribed(t,x,schedCoast);
    ph(2).terminate = @(t,x) coorbital.prop.eventApogee(t,x);
        ph(2).tspan = [0 tMaxCoast];
         ph(2).link = [];

%% Phase 3, descent to impact: unpowered, ends at the impact altitude:
          ph(3).eom = eomCst;
        ph(3).guide = @(t,x) coorbital.guide.prescribed(t,x,schedCoast);
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

%% Propagate. The vehicle argument is a formality here -- every phase binds
%% its own, see the note in the header -- so the stack is passed to keep the
%% call honest about what is on the pad:
              traj = coorbital.prop.phaseRun(ph,x0,bst,env);

%% Phase boundaries. phaseRun records the boundary sample ONCE, carrying the
%% outgoing phase's control, so the last row of phase k is that phase's
%% terminal state on the NEAR side of any staging jump:
               kBO = find(traj.phaseIdx == 1,1,'last');
              kApo = find(traj.phaseIdx == 2,1,'last');
                nS = numel(traj.t);
    assert(~isempty(kBO) && ~isempty(kApo), ...
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
              vehK = coastVeh;
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

%% Sensed specific force -- what an accelerometer on board would read.
%% Gravity is not sensed and is excluded. nAero drops the thrust, which is the
%% quantity meant by "deceleration" during re-entry:
             nSens = sqrt((aThrV - aDrag).^2 + (aThrN + aLift).^2)./c.g0;
             nAero = sqrt(aLift.^2 + aDrag.^2)./c.g0;

%% Mass must be constant once the motor is out, and it must be the mass the
%% coast vehicle was built around, because the six-state glide equations
%% divide by coastVeh.mass rather than by the mass state. Two separate claims,
%% asserted separately because they can fail separately and have different
%% budgets. Constancy is EXACT -- massConstant returns dm/dt = 0 identically,
%% so every Runge-Kutta stage adds exactly zero. Agreement with mCoast is only
%% as good as the burnout event solve, which lands a fraction of a nanogram
%% off the target mass:
             isUnp = traj.phaseIdx >= 2;
             mSpan = max(massS(isUnp)) - min(massS(isUnp));
    assert(mSpan == 0, ...
        ['the unpowered mass varied by %.3e kg; dm/dt is identically zero ' ...
         'after burnout, so it must be bit-exactly constant'],mSpan);
    assert(abs(massS(end) - mCoast) < 1e-6, ...
        ['the coast is flying %.9f kg while the coast vehicle was built ' ...
         'around %.9f kg; the aerodynamic accelerations are being divided ' ...
         'by the wrong mass'],massS(end),mCoast);

%% Peaks, taken over the phases where each is meaningful:
             isBst = traj.phaseIdx == 1;
             isRe  = traj.phaseIdx == 3;
     [qBstMax,kQB] = maxOver(qbar,isBst);
     [qReMax ,kQR] = maxOver(qbar,isRe);
     [nBstMax,kNB] = maxOver(nSens,isBst);
     [nReMax ,kNR] = maxOver(nAero,isRe);

%% Great-circle ranges. Everything is measured on the impact sphere r = rE +
%% hStop, so the launch-to-impact range and the burnout-to-impact range in the
%% cross-check below share one baseline:
                rI = c.rE + hStopM;
            angTot = coorbital.util.greatCircle(latLaunchR,lonLaunchR, ...
                                                traj.x(end,3),traj.x(end,2));
           rangeKm = rI.*angTot./1000;
             angBO = coorbital.util.greatCircle(latLaunchR,lonLaunchR, ...
                                                traj.x(kBO,3),traj.x(kBO,2));
          downBOKm = rI.*angBO./1000;

%% ---------------------------------------------------------------------
%% Keplerian cross-check, posed in two SEPARATE steps
%% ---------------------------------------------------------------------
%% Step one is a validation and is asserted. The SAME burnout state is
%% re-propagated with lift and drag switched off and rotation off, and
%% compared against the closed-form free-flight range of the Keplerian arc
%% through that state. Both sides are then exact descriptions of the same
%% two-body problem, so they must agree to integration tolerance -- and a
%% loose budget here would hide a real defect rather than absorb a modelling
%% difference.
%%
%% Step two is NOT a validation. Differencing the vacuum range against the
%% drag-included range measures a model discrepancy whose size nothing in this
%% file predicts, so it is printed and judged, never asserted.
%%
%% Step two is run THREE TIMES, and it has to be. Switching env.aero off
%% removes LIFT as well as drag, so a single vacuum-against-full difference is
%% the total aerodynamic effect, not the drag effect -- and the two have
%% opposite signs. Drag can only shorten a ballistic arc; upward lift, which is
%% what a positive CL at zero bank produces, lengthens it. Reporting their sum
%% as "the drag effect" would invite exactly the wrong verdict on the sign
%% whenever lift wins, which for a high-ballistic-coefficient re-entry body on
%% a steep arc it can. A propagation with lift alone suppressed therefore
%% separates them, and the sign judgement below is passed on the drag-only
%% figure, which theory does predict.
%%
%% AND ALL THREE COMPARISON ARCS FLY THE SAME GRAVITY AND THE SAME ROTATION.
%% This is the third of them, and it exists because the FLOWN trajectory does
%% not. The flown arc carries whatever gravFn and earthSpin the user selected;
%% the vacuum arc is pinned to sphereGrav with rotation off, because only a
%% two-body arc can be compared against a two-body closed form. Differencing
%% one against the other therefore folded the GRAVITY-MODEL difference into a
%% number the summary called drag, and the same contamination reached lift
%% through the flown arc. Every difference below is now taken between arcs that
%% differ in ONE model and nothing else; the gap between the pinned
%% full-aerodynamics arc and the flown one is reported separately, as what it
%% is.

%% Burnout state as HANDED TO THE COAST -- the far side of the staging jump:
               xBO = traj.junction(1).x;
                r0 = xBO(1);
                V0 = xBO(4);
              gam0 = xBO(5);

%% Free propagation from that state: same equations, zero aerodynamics, no
%% rotation, POINT-MASS gravity, one phase all the way down to the impact
%% altitude.
%%
%% All three of those pins are deliberate, and env.grav is the one a reader
%% will not expect. The closed form below is Keplerian BY CONSTRUCTION -- it
%% is built from c.muE and nothing else, so it describes the arc of a point
%% mass. Inheriting a user-selected oblate gravity here would propagate a J2
%% arc and compare it against a two-body reference, firing the 1e-6 assertion
%% on a MODELLING DIFFERENCE, which is precisely what that assertion's message
%% says a failure is not. Measured with a scratch J2 model before this pin was
%% added: the script aborted at 3.2e-03 relative. The MAIN run still flies
%% whatever gravity the user selected in the block above; only this
%% cross-check is forced back to sphereGrav, because only this cross-check has
%% a two-body closed form to answer to:
        envVac      = env;
        envVac.aero = @(alphaArg,machArg,vehArg) deal(0,0);
        envVac.grav = @coorbital.grav.sphereGrav;
      envVac.omegaE = 0;
          phV.eom   = eomCst;
        phV.guide   = @(t,x) coorbital.guide.prescribed(t,x,schedCoast);
    phV.terminate   = @(t,x) coorbital.prop.eventAltitude(t,x,hStopM);
        phV.tspan   = [0 tMaxCoast + tMaxDesc];
             trajV  = coorbital.prop.phaseRun(phV,xBO,coastVeh,envVac);
    assert(abs((trajV.x(end,1) - c.rE) - hStopM) < 1e-3, ...
        ['the vacuum re-propagation did not reach the %.1f km impact ' ...
         'altitude; it stopped at %.3f km after %.1f s'], ...
        hStop,(trajV.x(end,1) - c.rE)./1000,trajV.t(end));
            angVacP = coorbital.util.greatCircle(xBO(3),xBO(2), ...
                                                 trajV.x(end,3),trajV.x(end,2));
           sVacProp = rI.*angVacP;

%% Closed-form free-flight range of the same arc:
                 hM = r0.*V0.*cos(gam0);
                 pM = hM.^2./c.muE;
               epsM = V0.^2./2 - c.muE./r0;
                eOr = sqrt(1 + 2.*epsM.*hM.^2./c.muE.^2);
    assert(eOr > 1e-12, ...
        'the burnout arc is numerically circular (e = %.3e); the range formula is degenerate.',eOr);
                 f0 = atan2( hM.*V0.*sin(gam0)./c.muE , pM./r0 - 1 );
              cosfI = (pM./rI - 1)./eOr;
    assert(abs(cosfI) <= 1, ...
        ['the Keplerian arc through the burnout state never reaches the ' ...
         'impact radius: (p/rI - 1)/e = %.6f is outside [-1,1]. Burnout was ' ...
         'r = %.1f km, V = %.1f m/s, gamma = %.3f deg, giving e = %.6f and ' ...
         'p = %.1f km against rI = %.1f km.'], ...
        cosfI,r0./1000,V0,rad2deg(gam0),eOr,pM./1000,rI./1000);
                 fI = 2.*pi - acos(cosfI);
               sVac = rI.*mod(fI - f0,2.*pi);

%% Agreement between the two, which is the assertable half:
             relVac = abs(sVacProp - sVac)./sVac;
             tolVac = 1e-6;
    assert(relVac < tolVac, ...
        ['the vacuum propagation and the closed-form Keplerian range ' ...
         'disagree by %.3e relative (%.3f km on %.3f km), budget %.1e. ' ...
         'These are two exact descriptions of the same two-body arc, so a ' ...
         'disagreement is a defect, not a modelling difference.'], ...
        relVac,(sVacProp - sVac)./1000,sVac./1000,tolVac);

%% The full-aerodynamics range, on the same baseline as the vacuum figures:
%% burnout point to impact point, on the impact sphere:
            angAeroB = coorbital.util.greatCircle(xBO(3),xBO(2), ...
                                                  traj.x(end,3),traj.x(end,2));
              sAeroB = rI.*angAeroB;

%% ...and the same arc again with LIFT alone suppressed, so drag's own
%% contribution can be separated from lift's. Same burnout state, same
%% equations, same drag coefficient, CL forced to zero.
%%
%% IT IS BUILT FROM envVac AND NOT FROM env, and that is the whole point of
%% this line. envVac pins the gravity model to sphereGrav and the rotation
%% rate to zero; inheriting env instead would leave a user-selected oblate
%% gravity, or Earth rotation, INSIDE the difference that the summary then
%% calls the drag effect. Measured before the pin was added: with a scratch J2
%% model the "drag effect" moved by kilometres that had nothing to do with
%% drag, and the sign judgement below is passed on that number:
        envDrg      = envVac;
        envDrg.aero = @(alphaArg,machArg,vehArg) dragOnly(aeroFn,alphaArg,machArg,vehArg);
          phD       = phV;
             trajD  = coorbital.prop.phaseRun(phD,xBO,coastVeh,envDrg);
    assert(abs((trajD.x(end,1) - c.rE) - hStopM) < 1e-3, ...
        ['the drag-only re-propagation did not reach the %.1f km impact ' ...
         'altitude; it stopped at %.3f km after %.1f s'], ...
        hStop,(trajD.x(end,1) - c.rE)./1000,trajD.t(end));
            angDrgB = coorbital.util.greatCircle(xBO(3),xBO(2), ...
                                                 trajD.x(end,3),trajD.x(end,2));
              sDrgB = rI.*angDrgB;

%% ...and a THIRD re-propagation, with the full aerodynamic model but the SAME
%% pinned gravity and rotation. The flown trajectory cannot serve here: it
%% carries whatever gravity and rotation the user selected, so differencing it
%% against the drag-only arc would attribute a gravity-model change to LIFT.
%% With the shipped handles this arc and the flown one are the same physics and
%% agree to integration tolerance, which the summary reports as a check rather
%% than assumes:
        envFul      = envVac;
        envFul.aero = aeroFn;
          phF       = phV;
             trajF  = coorbital.prop.phaseRun(phF,xBO,coastVeh,envFul);
    assert(abs((trajF.x(end,1) - c.rE) - hStopM) < 1e-3, ...
        ['the full-aerodynamics re-propagation did not reach the %.1f km ' ...
         'impact altitude; it stopped at %.3f km after %.1f s'], ...
        hStop,(trajF.x(end,1) - c.rE)./1000,trajF.t(end));
            angFulB = coorbital.util.greatCircle(xBO(3),xBO(2), ...
                                                 trajF.x(end,3),trajF.x(end,2));
              sFulB = rI.*angFulB;

%% The three differences. ALL FOUR ARCS -- vacuum, drag-only, full-aero and the
%% closed form -- now share one gravity model and one rotation rate, so each
%% difference isolates exactly the one model that was switched:
            aeroDiff = sFulB - sVacProp;
             aeroPct = 100.*aeroDiff./sVacProp;
            dragDiff = sDrgB - sVacProp;
             dragPct = 100.*dragDiff./sVacProp;
            liftDiff = sFulB - sDrgB;
             liftPct = 100.*liftDiff./sVacProp;

%% How far the pinned full-aerodynamics arc sits from the FLOWN one. Zero when
%% the shipped handles are in place, and a measure of how much of the flown
%% range comes from a non-default gravity or from rotation when they are not:
            flownGap = sAeroB - sFulB;

%% ---------------------------------------------------------------------
%% Termination diagnosis, one line per phase
%% ---------------------------------------------------------------------
%% An unexpected termination must SAY so. A boost that ran out of horizon
%% before it ran out of propellant, or a coast that never turned over, both
%% produce a short trajectory that otherwise reads as a completed flight:
     [why1,ok1] = whyBoost(traj,kBO,mBurnoutT,tMaxBoost);
     [why2,ok2] = whyApogee(traj,kBO,kApo,tMaxCoast);
     [why3,ok3] = whyImpact(traj,kApo,nS,c.rE,hStopM,hStop,tMaxDesc);
        allOK   = ok1 && ok2 && ok3;

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

%% Which of the swappable models are still the library defaults. The validity
%% paragraph below describes the DEFAULT models' limitations by name -- an
%% isothermal windless atmosphere, a spherical no-J2 gravity field, a motor
%% that cannot throttle -- and every one of those statements becomes a
%% FALSEHOOD the moment the corresponding handle is replaced. A summary that
%% prints "the spherical j2Grav carries no J2" is worse than one that prints
%% nothing, so each caveat is gated on its own handle:
       defaultAtmos = isequal(atmosFn,@coorbital.atmos.expAtmos);
        defaultGrav = isequal(gravFn,@coorbital.grav.sphereGrav);
        defaultProp = isequal(propFn,@coorbital.prop.constThrust);

%% Report:
    fprintf('\n');
    fprintf('===== Ballistic trajectory summary ======================\n');
    fprintf('  Phase termination\n');
    fprintf('    1 boost          %s\n',why1);
    fprintf('    2 coast          %s\n',why2);
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
    fprintf('  Burnout, end of phase 1\n');
    fprintf('    time             %10.3f s\n',traj.t(kBO));
    fprintf('    altitude         %10.3f km\n',hKm(kBO));
    fprintf('    speed            %10.2f %-5s  (Mach %.2f, planet-relative)\n', ...
            V(kBO),'m/s',mach(kBO));
    fprintf('    flight path      %10.3f %-5s  (positive is climbing)\n', ...
            rad2deg(traj.x(kBO,5)),'deg');
    fprintf('    heading          %10.3f deg\n',rad2deg(traj.x(kBO,6)));
    fprintf('    downrange        %10.2f %-5s  (great circle from the pad)\n',downBOKm,'km');
    fprintf('    mass at burnout  %10.1f %-5s  (payload + spent booster)\n',massS(kBO),'kg');
    fprintf('    mass into coast  %10.1f %-5s  (%s)\n',mCoast,'kg',sepTxt);
    fprintf('\n');
    fprintf('  Apogee, end of phase 2\n');
    fprintf('    time             %10.3f s\n',traj.t(kApo));
    fprintf('    altitude         %10.3f %-5s  (exact event, not a grid maximum)\n',hKm(kApo),'km');
    fprintf('    speed            %10.2f m/s\n',V(kApo));
    fprintf('    flight path      %10.3e %-5s  (zero by construction at apogee)\n', ...
            rad2deg(traj.x(kApo,5)),'deg');
    fprintf('\n');
    fprintf('  Overall\n');
    fprintf('    flight time      %10.2f %-5s  (%.2f min)\n',traj.t(end),'s',traj.t(end)./60);
    fprintf('    ground range     %10.2f %-5s  (great circle on the r = %.1f km impact sphere)\n', ...
            rangeKm,'km',rI./1000);
    fprintf('    central angle    %10.4f deg\n',rad2deg(angTot));
    fprintf('    range check      %10.3f %-5s  = %.3f (pad to burnout) + %.3f (burnout to\n', ...
            downBOKm + sAeroB./1000,'km',downBOKm,sAeroB./1000);
    fprintf('                                         impact). The legs add exactly because\n');
    fprintf('                                         the unbanked non-rotating arc is a great circle.\n');
    fprintf('    apogee           %10.2f %-5s  (%.4f of the range; a minimum-energy arc sits\n', ...
            hKm(kApo),'km',hKm(kApo)./rangeKm);
    fprintf('                                        near 0.25, so this trajectory is LOFTED)\n');
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
    fprintf('    re-entry max q   %10.2f %-5s  (at t = %.1f s, h = %.2f km)\n', ...
            qReMax./1000,'kPa',traj.t(kQR),hKm(kQR));
    fprintf('    re-entry decel   %10.2f %-5s  (%.2f m/s^2 aero only, at t = %.1f s, h = %.2f km)\n', ...
            nReMax,'g',nReMax.*c.g0,traj.t(kNR),hKm(kNR));
    fprintf('\n');
    fprintf('  Keplerian cross-check, from the burnout state, on the r = %.1f km sphere\n',rI./1000);
    fprintf('    EVERY range in this block is BURNOUT to impact, not pad to impact. Add the\n');
    fprintf('    %.3f km flown under power to compare against the %.2f km ground range\n', ...
            downBOKm,rangeKm);
    fprintf('    above: %.3f + %.3f = %.3f km. Do not difference the two directly.\n', ...
            downBOKm,sAeroB./1000,downBOKm + sAeroB./1000);
    fprintf('    closed-form vacuum range %10.4f km\n',sVac./1000);
    fprintf('    propagated vacuum range  %10.4f %-5s  (lift and drag OFF, rotation OFF)\n', ...
            sVacProp./1000,'km');
    fprintf('    agreement                %10.3e %-5s  (relative; budget %.1e, asserted)\n', ...
            relVac,'',tolVac);
    fprintf('    ---\n');
    fprintf('    ALL THREE ARCS BELOW share the propagated vacuum arc''s gravity model and\n');
    fprintf('    its zero rotation rate, and differ from it in ONE model each, so every\n');
    fprintf('    difference isolates the model that was switched and nothing else.\n');
    fprintf('    drag-only range          %10.4f %-5s  (CL forced to zero, CD kept)\n', ...
            sDrgB./1000,'km');
    fprintf('    full-aerodynamics range  %10.4f %-5s  (re-propagated, NOT the flown arc)\n', ...
            sFulB./1000,'km');
    fprintf('    drag effect              %10.4f %-5s  (%+.4f %% of the vacuum range)\n', ...
            dragDiff./1000,'km',dragPct);
    fprintf('    lift effect              %10.4f %-5s  (%+.4f %%)\n', ...
            liftDiff./1000,'km',liftPct);
    fprintf('    net aerodynamic effect   %10.4f %-5s  (%+.4f %%)\n', ...
            aeroDiff./1000,'km',aeroPct);
    fprintf('    flown minus re-propagated %9.4f %-5s  (what the SELECTED gravity model and\n', ...
            flownGap./1000,'km');
    fprintf('                                         rotation add to the flown arc; it is\n');
    fprintf('                                         NOT part of any figure above)\n');
    if ~defaultGrav || earthSpin
        fprintf('    THE FLOWN ARC IS NOT THE ARC ABOVE. Gravity is %s and Earth rotation\n', ...
                func2str(gravFn));
        fprintf('    is %s, so the flown burnout-to-impact range of %.4f km is %.4f km from\n', ...
                spinTxt,sAeroB./1000,flownGap./1000);
        fprintf('    the re-propagated one. That gap is a GRAVITY and ROTATION effect and is\n');
        fprintf('    deliberately kept out of the drag and lift figures.\n');
    end
    fprintf('    These differences are MODEL DISCREPANCIES, reported and NOT asserted:\n');
    fprintf('    nothing in this file predicts their size, so no tolerance is attached.\n');
    if dragDiff < 0
        fprintf('    Drag sign: NEGATIVE, which is what drag must do. It removes energy on the\n');
        fprintf('    way up and again on re-entry, and can only shorten a ballistic arc.\n');
    else
        fprintf('    Drag sign: POSITIVE, which drag CANNOT do. RED FLAG: check the sign of CD\n');
        fprintf('    in the aero model and whether the Earth-rotation switch is contaminating\n');
        fprintf('    the comparison.\n');
    end
    if liftDiff > 0
        fprintf('    Lift sign: POSITIVE, as expected FOR THIS CONFIGURATION -- CL = %.4f at\n', ...
                coastVeh.CL);
        fprintf('    %.1f deg of bank puts the lift vector up, and at a %.1f km burnout altitude\n', ...
                bankAngle,(r0 - c.rE)./1000);
        fprintf('    the air is thin enough that the arc meets it only on the terminal dive,\n');
        fprintf('    where raising gamma flattens the descent and extends the range. That\n');
        fprintf('    premise is what fixes the sign; upward lift does NOT lengthen every\n');
        fprintf('    ballistic arc -- see the note in the other branch.\n');
    else
        fprintf('    Lift sign: NEGATIVE. That is not automatically a defect. Two ways it is\n');
        fprintf('    legitimate: (1) if the ascending leg is lofted enough for lift to raise\n');
        fprintf('    the burnout-equivalent gamma ABOVE the max-range value, extra loft costs\n');
        fprintf('    range rather than buying it -- at this burnout state the max-range gamma\n');
        fprintf('    is near 33 deg against the %.1f deg flown; (2) a nonzero bank, or a\n', ...
                rad2deg(gam0));
        fprintf('    negative CL, points the lift vector sideways or down. Check both before\n');
        fprintf('    treating it as a bug.\n');
    end
    fprintf('    Magnitude: small in RANGE because the atmosphere is met only on a steep\n');
    fprintf('    terminal dive, where the aerodynamic force acts nearly along the descent\n');
    fprintf('    direction and so moves the impact point very little. It is not small in\n');
    fprintf('    SPEED, which is where the same drag shows up:\n');
    fprintf('    impact speed, vacuum     %10.2f m/s\n',trajV.x(end,4));
    fprintf('    impact speed, flown      %10.2f %-5s  (%+.2f m/s)\n', ...
            V(end),'m/s',V(end) - trajV.x(end,4));
    if earthSpin
        fprintf('    NOTE: Earth rotation is ON for the full run and OFF for the vacuum\n');
        fprintf('    propagation, so the difference above is drag PLUS rotation, not drag.\n');
    end
    fprintf('\n');
    fprintf('  Models: atmos %s | grav %s | aero %s\n', ...
            func2str(atmosFn),func2str(gravFn),func2str(aeroFn));
    fprintf('          prop %s | Earth rotation %s\n',func2str(propFn),spinTxt);
    fprintf('  Vehicle: %s, mass %.1f kg, Sref %.3f m^2, CL %.3f, L/D %.2f (PLACEHOLDER values)\n', ...
            func2str(vehicleFn),veh.mass,veh.Sref,veh.CL,veh.LD);
    fprintf('  Booster: %s, dry %.1f kg, prop %.1f kg, Tvac %.1f kN, Isp %.1f s (PLACEHOLDER values)\n', ...
            func2str(boosterFn),bst.massDry,bst.massProp,bst.thrustVac./1000,bst.Isp);
    fprintf('  Validity: %s holds CL and L/D constant, a hypersonic approximation, and\n', ...
            func2str(aeroFn));
    fprintf('            this flight runs from a subsonic liftoff through Mach %.0f and back\n',machMaxOf(mach));
    fprintf('            down again, so the transonic ends of the boost and of the re-entry\n');
    fprintf('            lie outside it.\n');
    if defaultAtmos
        fprintf('            The atmosphere %s is isothermal and carries no winds.\n', ...
                func2str(atmosFn));
    else
        fprintf('            The atmosphere is %s, NOT the library default, so the\n', ...
                func2str(atmosFn));
        fprintf('            isothermal-and-windless caveat this script prints by default does\n');
        fprintf('            not describe it. Its validity range is its own to state.\n');
    end
    if defaultGrav
        fprintf('            The gravity model %s is spherical and carries no J2.\n', ...
                func2str(gravFn));
    else
        fprintf('            Gravity is %s, NOT the default sphereGrav, so the no-J2\n', ...
                func2str(gravFn));
        fprintf('            caveat this script prints by default does not apply. NOTE that the\n');
        fprintf('            Keplerian cross-check above is still flown two-body on purpose.\n');
    end
    fprintf('            Treat this as an indicative trajectory, not a performance\n');
    fprintf('            prediction.\n');
    if mach(end) < machHyper
        fprintf('            Impact is at Mach %.2f, below Mach %.0f, so the terminal state in\n', ...
                mach(end),machHyper);
        fprintf('            particular is outside the constant-coefficient regime.\n');
    end
    if nBstMax > 15 && defaultProp
        fprintf('            The %.1f g sensed load at end of burn is an artefact of the\n',nBstMax);
        fprintf('            constant-thrust placeholder motor, which does not throttle as the\n');
        fprintf('            stack empties. A real stage would tail off or stage before that.\n');
    elseif nBstMax > 15
        fprintf('            The %.1f g sensed load at end of burn is high. This run uses\n',nBstMax);
        fprintf('            %s, NOT the default constant-thrust motor, so\n',func2str(propFn));
        fprintf('            whether that load is a throttling artefact depends on that model.\n');
    end
    fprintf('=========================================================\n\n');

%% Hand the summary's numbers back at full precision, so a test does not have
%% to read them out of printed text:
       info.tBurnout = traj.t(kBO);
       info.hBurnout = traj.x(kBO,1) - c.rE;
       info.vBurnout = V(kBO);
     info.gamBurnout = traj.x(kBO,5);
       info.mBurnout = massS(kBO);
         info.mCoast = mCoast;
       info.mLiftoff = mLiftoff;
     info.downBOKm   = downBOKm;
        info.tApogee = traj.t(kApo);
        info.hApogee = traj.x(kApo,1) - c.rE;
        info.vApogee = V(kApo);
       info.tFlight  = traj.t(end);
       info.rangeKm  = rangeKm;
        info.angTot  = angTot;
       info.vImpact  = V(end);
     info.gamImpact  = traj.x(end,5);
     info.latImpact  = traj.x(end,3);
     info.lonImpact  = traj.x(end,2);
       info.qBstMax  = qBstMax;
       info.qReMax   = qReMax;
       info.nBstMax  = nBstMax;
       info.nReMax   = nReMax;
       info.tNReMax  = traj.t(kNR);
     info.sVacClosed = sVac;
       info.sVacProp = sVacProp;
      info.vImpVac   = trajV.x(end,4);
        info.relVac  = relVac;
        info.tolVac  = tolVac;
         info.sAeroB = sAeroB;
          info.sFulB = sFulB;
      info.flownGap  = flownGap;
          info.sDrgB = sDrgB;
      info.dragDiff  = dragDiff;
       info.dragPct  = dragPct;
      info.liftDiff  = liftDiff;
       info.liftPct  = liftPct;
      info.aeroDiff  = aeroDiff;
       info.aeroPct  = aeroPct;
      info.stopOK    = allOK;
      info.stopWhy   = {why1;why2;why3};

%% Vertical exaggeration for the globe. A DISPLAY choice, not a physical
%% constant, so it does not belong in missileConst. This flight apogees near
%% 1619 km on a 6378 km sphere, which is already visible at true scale, so the
%% factor is small and only lifts the arc clear of the surface at its ends;
%% coorbital.viz.globe3D writes it into the title either way:
           altExag = 3;

%% Plots. Every figure comes from coorbital.viz, which reads the trajectory
%% and never writes it -- see the header of each. Nothing below this line can
%% move a number in the summary above.
%%
%% ONE VEHICLE PER PHASE is handed over, because this chain does not fly one:
%% phase 1 is the stack on the pad and phases 2 and 3 are whatever survived
%% staging. The load-factor panel is the only one that reads it, and given
%% nothing it would draw the whole flight on the booster's reference area.
%% See the "two vehicles, one chain" note in the header of this file.
%%
%% SIX CHANNELS AND A SEVENTH PANEL. The 'nAero' channel is the aerodynamic
%% load with thrust EXCLUDED, and on a boosted flight that is near nothing
%% exactly where the interesting load is -- the sensed load peaks near 40 g at
%% end of burn on this booster while the air is still thin. So the sensed
%% load, which this script already computes for its own summary, is handed to
%% the plot as an extra panel rather than left to a single printed peak: a
%% peak is a number, and what a reader needs to see is the build-up and when
%% it happened. coorbital.viz.profilePlot does not compute it, and should not
%% -- the thrust term needs a propulsion model, an ambient pressure and a list
%% of which phases are burning, none of which a trajectory carries:
    if showPlots
        coorbital.viz.profilePlot(traj,bst,env, ...
            struct('Name','Ballistic profile', ...
                   'Channels',{{'altitude','speed','mach','q','nAero','mass'}}, ...
                   'Extra',{{nSens,'sensed load factor (g)', ...
                             'Sensed load factor, thrust included'}}, ...
                   'VehPhase',{{bst,coastVeh,coastVeh}}));
        coorbital.viz.groundTrack(traj,bst,env, ...
            struct('PhaseName',{{'boost','coast','descent'}}, ...
                   'Title',sprintf('Ground track, %.0f km great-circle range', ...
                                   rangeKm)));
        coorbital.viz.globe3D(traj,bst,env, ...
            struct('AltScale',altExag, ...
                   'Title','Boost, coast and descent to impact'));
    end

%% Hand the trajectory back only when the caller asked for it. Typing
%% "run_ballistic" at the prompt should leave the summary on screen, not bury
%% it under a dump of the whole struct as ans:
    if nargout == 0
        clear traj info;
    end
end

function [CL,CD] = dragOnly(aeroFn,alpha,mach,veh)
%% Purpose:
%
%  Wrap an aerodynamic model so that it returns its own drag coefficient and
%  ZERO lift. Used by the cross-check to fly the same arc a third time with
%  lift alone suppressed, which is what separates drag's effect on range from
%  lift's. Written against the model handle rather than against constLD, so
%  the separation survives a swap to a table-driven or Mach-dependent model.
%
%% Inputs:
%
%  aeroFn           Function handle             [CL,CD] = aeroFn(alpha,mach,veh)
%
%  alpha            [N x 1]                     Angle of attack (rad)
%
%  mach             [N x 1]                     Mach number (-)
%
%  veh              Struct                      Vehicle parameters
%
%% Outputs:
%
%  CL               scalar                      Zero, always
%
%  CD               scalar                      Drag coefficient from aeroFn (-)
%
%% Revision History:
%  Michael Casey                                                08/07/2026
%  Copyright 2026 Coorbital, Inc.
%% ------------------------ Begin Code Sequence ---------------------------

             [~,CD] = aeroFn(alpha,mach,veh);
                 CL = 0;
end

function m = machMaxOf(mach)
%% Purpose:
%
%  Largest Mach number reached anywhere on the trajectory, rounded nowhere.
%  Factored out only so the validity paragraph reads as one line.
%
%% Inputs:
%
%  mach             [N x 1]                     Mach number history (-)
%
%% Outputs:
%
%  m                [1 x 1]                     Peak Mach number (-)
%
%% Revision History:
%  Michael Casey                                                08/07/2026
%  Copyright 2026 Coorbital, Inc.
%% ------------------------ Begin Code Sequence ---------------------------

                 m = max(mach);
end

function [pk,kAt] = maxOver(v,mask)
%% Purpose:
%
%  Largest value of a history within a masked subset of its samples, and the
%  index into the FULL history at which it occurs. Used so that a boost peak
%  and a re-entry peak can be reported separately without either search
%  wandering into the other's phase.
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

function [why,ok] = whyBoost(traj,kBO,mBurnoutT,tMaxBoost)
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

function [why,ok] = whyApogee(traj,kBO,kApo,tMaxCoast)
%% Purpose:
%
%  Diagnose why the ascending coast stopped. Nominal termination is apogee,
%  which coorbital.prop.eventApogee detects as the flight-path angle crossing
%  zero descending. A coast that simply ran out of horizon while still
%  climbing is the failure this exists to name.
%
%% Inputs:
%
%  traj             Struct                      Trajectory from phaseRun
%
%  kBO              [1 x 1]                     Index of the last phase-1 sample
%
%  kApo             [1 x 1]                     Index of the last phase-2 sample
%
%  tMaxCoast        [1 x 1]                     Ascending-coast horizon (s)
%
%% Outputs:
%
%  why              Char [1 x n]                Human-readable reason
%
%  ok               [1 x 1] logical             True only for a nominal apogee
%
%% Revision History:
%  Michael Casey                                                08/07/2026
%  Copyright 2026 Coorbital, Inc.
%% ------------------------ Begin Code Sequence ---------------------------

              tPh2 = traj.t(kApo) - traj.t(kBO);
              gEnd = traj.x(kApo,5);
                ok = false;
    if abs(gEnd) < 1e-9
               why = sprintf(['flight path angle reached zero at t = %.3f s, ' ...
                              'apogee (nominal)'],traj.t(kApo));
                ok = true;
    elseif tPh2 >= tMaxCoast - 1e-6
               why = sprintf(['hit the %.0f s coast horizon still at gamma = %.4f deg; ' ...
                              'apogee was never reached'],tMaxCoast,rad2deg(gEnd));
    else
               why = sprintf(['stopped early at t = %.3f s with gamma = %.4f deg ' ...
                              '(integrator failure or an unmodelled event)'], ...
                             traj.t(kApo),rad2deg(gEnd));
    end
end

function [why,ok] = whyImpact(traj,kApo,nS,rE,hStopM,hStop,tMaxDesc)
%% Purpose:
%
%  Diagnose why the descent stopped. Nominal termination is the descending
%  crossing of the impact altitude. A descent that ran out of horizon in
%  mid-air is the failure this exists to name.
%
%% Inputs:
%
%  traj             Struct                      Trajectory from phaseRun
%
%  kApo             [1 x 1]                     Index of the last phase-2 sample
%
%  nS               [1 x 1]                     Number of samples in traj
%
%  rE               [1 x 1]                     Planet radius (m)
%
%  hStopM           [1 x 1]                     Impact altitude (m)
%
%  hStop            [1 x 1]                     Impact altitude (km), for the message
%
%  tMaxDesc         [1 x 1]                     Descent horizon (s)
%
%% Outputs:
%
%  why              Char [1 x n]                Human-readable reason
%
%  ok               [1 x 1] logical             True only for a nominal impact
%
%% Revision History:
%  Michael Casey                                                08/07/2026
%  Copyright 2026 Coorbital, Inc.
%% ------------------------ Begin Code Sequence ---------------------------

              tPh3 = traj.t(nS) - traj.t(kApo);
              hEnd = traj.x(nS,1) - rE;
                ok = false;
    if abs(hEnd - hStopM) < 1e-3
               why = sprintf(['reached the %.1f km impact altitude at t = %.3f s ' ...
                              '(nominal)'],hStop,traj.t(nS));
                ok = true;
    elseif tPh3 >= tMaxDesc - 1e-6
               why = sprintf(['hit the %.0f s descent horizon at %.3f km, still ' ...
                              'flying'],tMaxDesc,hEnd./1000);
    else
               why = sprintf(['stopped early at t = %.3f s, h = %.3f km ' ...
                              '(integrator failure or an unmodelled event)'], ...
                             traj.t(nS),hEnd./1000);
    end
end

function v = overrideOf(opts,name,v)
%% Purpose:
%
%  Return the caller's override for one USER PARAMETERS entry when it was
%  supplied, and the script's own value otherwise. Factored out so the
%  override block above reads as one unbranched line per parameter.
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
