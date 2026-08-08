function [traj,info] = run_target(opts)
%% Purpose:
%
%  Fly a boost-glide vehicle from a LAUNCH POINT to a DESTINATION POINT. Give
%  it two latitude/longitude pairs and it solves the trajectory that connects
%  them, flies the full boost-glide-descent chain, reports what it did, draws
%  the figures, and optionally renders a movie of the Earth with the
%  trajectory developing over it. Everything a routine run needs to change
%  lives in the USER PARAMETERS block; nothing below it should require
%  editing.
%
%  Angles and distances in the USER PARAMETERS block are in degrees and
%  kilometres because that is how a user thinks about them. They are converted
%  to the library's SI units (m, m/s, rad, s) immediately after the block, in
%  one place, and nothing below that point works in any other unit.
%
%% Note -- the targeting problem is TWO separate solves:
%
%  HGV/run_boost_glide flies launch site + azimuth + pitch program and lands
%  wherever the physics puts it. This script inverts that. Two unknowns, two
%  independent solves, and they do not interact on a non-rotating Earth:
%
%    Azimuth, CLOSED FORM.  coorbital.util.greatCircleBearing returns the
%                           initial bearing of the great-circle arc from the
%                           launch point to the target, clockwise from north,
%                           which is exactly this library's heading convention
%                           psi. It drops straight into the launch state. No
%                           iteration, no propagation, no tolerance.
%
%    Range, ONE-DIMENSIONAL ROOT SOLVE ON THRUST-TERMINATION TIME.  The boost
%                           phase already ends at tspan(2) when the burnout
%                           event has not fired first, so a shortened tspan IS
%                           a thrust cutoff and no new machinery is needed.
%                           Less burn gives less energy gives less range:
%                           monotonic and single-valued, which is what makes
%                           bisection safe. coorbital.util.rangeSolve does the
%                           bisection; this script supplies the propagation.
%
%  WHY CUTOFF TIME AND NOT SOMETHING ELSE. The pitch program's loft angle has
%  a lofted and a depressed branch either side of a maximum-range hump, so a
%  root finder lands on whichever branch it started nearest -- a bad default
%  for an example script. The glide handoff altitude bifurcates on the glide
%  phugoid's troughs; 30 km costs 1882 km of range in run_boost_glide's own
%  configuration while every phase still reports nominal. Thrust termination
%  is also what real systems do for energy management.
%
%  EARLY CUTOFF LEAVES UNBURNED PROPELLANT, and it is thrown away. At
%  separation the ENTIRE booster goes -- dry structure and the propellant
%  still in it -- so the post-separation vehicle is exactly the payload and
%  the link is phases(1).link = @(x) [x(1:6); veh.mass]. The summary prints
%  how much propellant went overboard unburned, because at a deep cutoff that
%  is most of the load.
%
%% Note -- TWO LIMITATIONS, both printed in the summary rather than buried:
%
%  THE AZIMUTH IS EXACT ONLY BECAUSE THE EARTH DOES NOT ROTATE HERE. With
%  env.omegaE = 0 the ground track of a zero-bank trajectory is a great
%  circle, so the initial bearing of the launch-to-target arc is the whole
%  answer. Switch earthSpin on and the target is carried east under the
%  vehicle for the whole flight -- of order 600 km at mid latitude over a
%  27-minute flight -- so the initial bearing is no longer the answer, the
%  azimuth needs an OUTER ITERATION around the range solve, and this script
%  does not have one. It refuses to pretend: with earthSpin true it prints a
%  caution and the reported miss is against a target that stood still.
%
%  THE SOLVE CONTROLS DOWNRANGE ONLY. Bisection matches the great-circle
%  distance from launch to impact; nothing in it steers the track sideways.
%  Cross-range comes out at zero in the SHIPPED configuration because the
%  bank angle is zero throughout, and a zero-bank trajectory over a
%  non-rotating sphere stays on the great circle it left on. That is a
%  property of THIS CONFIGURATION, not a general guarantee, and this script
%  MEASURES it rather than asserting it: the cross-track offset of the impact
%  point from the launch-to-target great circle is computed from the flown
%  state and printed every run, with a warning when it exceeds the range
%  tolerance. Set descBank to run_boost_glide's 75 deg and it does: measured
%  21.5 km of miss against a 0.6 km range residual, on the shipped geometry.
%
%% Note -- three vehicles, one chain:
%
%  coorbital.prop.phaseRun carries ONE vehicle struct for the whole chain, but
%  the boosted stack and the separated payload have different reference areas,
%  different aerodynamic coefficients and different masses. Each phase
%  therefore BINDS its own vehicle in the equation-of-motion closure below,
%  and the struct handed to phaseRun is never read by the equations of motion.
%  coorbital.eom.massConstant refuses to run when the carried mass state
%  disagrees with the vehicle the equations divide by, which is what turns a
%  staging slip from a silent 33-fold weight error into an immediate stop.
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
%                                               editing it
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
%                                               EMPTY when the target is
%                                               outside the reachable envelope
%                                               -- see info.refused
%
%  info             Struct                      The summary's numbers at full
%                                               precision, so a test does not
%                                               have to read them back out of
%                                               printed text; all SI except the
%                                               *Km distances. Always carries
%                                               refused (logical), rngReqM,
%                                               psiLaunch, rngMinM, rngMaxM and
%                                               the rangeSolve record in
%                                               solveInfo, so a refused run
%                                               still hands back the envelope
%                                               it refused against. A run that
%                                               was not refused additionally
%                                               carries the DISPLAY choices the
%                                               figures were drawn with:
%                                               altExag     [1 x 1] the
%                                                           exaggeration used
%                                                           (-)
%                                               altExagRule handle, the rule
%                                                           itself, so a
%                                                           checker can ask
%                                                           what it would do at
%                                                           an apogee this run
%                                                           did not fly -- its
%                                                           cap and floor are
%                                                           otherwise
%                                                           unreachable
%                                               launchStr,  char, the caption
%                                               targetStr   text, hemisphere
%                                                           letters and all
%
%% References:
%   [1] Bowditch, N., "The American Practical Navigator," Pub. No. 9, NGA,
%       chapter on Great-Circle Sailing. The initial course used for the
%       launch azimuth; see coorbital.util.greatCircleBearing.
%   [2] Press, W.H., et al., "Numerical Recipes," 3rd ed., Cambridge, 2007,
%       Section 9.1. Bisection, used by coorbital.util.rangeSolve.
%   [3] Vinh, N.X., Busemann, A., Culp, R.D., "Hypersonic and Planetary Entry
%       Flight Mechanics," Univ. Michigan Press, 1980, Ch. 5.
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
%% Launch point -- where the vehicle leaves the pad:
         latLaunch = 20;               %deg, geocentric latitude [-89 .. 89]; the poles are singular
         lonLaunch = -155;             %deg, longitude [-180 .. 180]
           hLaunch = 0;                %km,  pad altitude above the sphere [0 .. 5]
           vLaunch = 10;               %m/s, speed at the first integrated point [2 .. 50].
                                       %     NOT a physical launch speed: the 3DOF equations
                                       %     are singular at V = 0, so the integration starts
                                       %     a moment after first motion

%% Target point -- where the vehicle is to arrive. There is NO azimuth entry:
%% the launch azimuth is SOLVED from these two points, in closed form, and any
%% azimuth given here could only contradict it:
         latTarget = 35;               %deg, geocentric latitude [-89 .. 89]
         lonTarget = -120;             %deg, longitude [-180 .. 180]
                                       %     The shipped pair is a mid-Pacific pad to the
                                       %     California coast: 3811 km on a 56.6 deg azimuth.
                                       %     Deliberately NOT equatorial and NOT due east --
                                       %     that geometry makes latitude and longitude
                                       %     interchangeable in the range formula and hides a
                                       %     transposition at every call site that uses them

%% Range solve -- bisection on thrust-termination time, expressed as a
%% FRACTION of the booster's full burn so the bracket survives a change of
%% booster. The endpoints are propagated before the solve starts and their
%% ranges ARE the reachable envelope, which the summary prints:
        cutFracMin = 0.50;             %-, shortest burn tried [0.3 .. 0.95]. Too small and the
                                       %   vehicle never climbs to hHandoff, the glide event
                                       %   never fires and the propagation is rejected; the
                                       %   shipped 0.50 reaches 203 km of range
        cutFracMax = 1.00;             %-, longest burn tried (cutFracMin .. 1.0]. 1.0 is the
                                       %   full burn and is also the maximum range this
                                       %   configuration can fly, 7738 km
        tolRangeKm = 1.0;              %km, convergence tolerance on the ACHIEVED range
                                       %    [0.05 .. 50]. Each halving of this costs one more
                                       %    trajectory propagation, and a propagation is about
                                       %    0.2 s, so tightening it is cheap

%% Phase 1, boost -- commanded pitch ATTITUDE against time since liftoff. The
%% angle of attack follows as alpha = theta - gamma, so a schedule that tracks
%% the natural gravity turn keeps alpha small. This program is DEPRESSED on
%% purpose -- it pitches over hard and burns out shallow and low, which is what
%% turns a booster into a boost-glide first stage rather than a lofted arc.
%% The schedule is followed only as far as the SOLVED cutoff:
         pitchTime = [0  6  15  30  50  70  82];   %s,   nodes, strictly increasing [0 .. burn time]
        pitchAngle = [89 80  55  28  12   5   2];  %deg, commanded attitude [0 .. 90], descending
          alphaMax = 10;               %deg, clamp on |angle of attack| [1 .. 15]; this clamp,
                                       %     not the schedule, is what limits how fast the
                                       %     flight path can be pushed over
         bankBoost = 0;                %deg, bank during boost [-90 .. 90]. MUST BE 0 for the
                                       %     closed-form azimuth to be the answer -- see the
                                       %     cross-range note in the header
         tMaxBoost = 200;              %s,  boost horizon; must exceed the burn time

%% Phase 2, glide -- prescribed control against time SINCE THE START OF THE
%% GLIDE, not since liftoff:
         glideTime = [0 6000];         %s,   nodes, strictly increasing, phase-local
        glideAlpha = [0 0];            %deg, angle of attack. NO EFFECT with the default constLD
                                       %     aero model, which ignores it; becomes live only when
                                       %     aeroFn below is swapped for an alpha-dependent model
         glideBank = [0 0];            %deg, bank [-90 .. 90]. 0 keeps the track on the
                                       %     launch-to-target great circle; anything else steers
                                       %     the vehicle off it and the solve cannot see that
          hHandoff = 15;               %km, glide-to-descent handoff, DESCENDING crossing only.
                                       %    Hard bound: hStop < hHandoff < the burnout altitude.
                                       %    With glideBank and descBank equal, as shipped, the
                                       %    handoff is BOOKKEEPING ONLY -- both phases fly the
                                       %    same control -- and its placement cannot change the
                                       %    range. Give descBank a different value and it can:
                                       %    see the phugoid-trough warning in run_boost_glide
         tMaxGlide = 6000;             %s, glide horizon; raise it if the glide is cut short

%% Phase 3, descent -- the terminal segment. Same airframe, same equations;
%% only the bank angle can differ from the glide:
          descTime = [0 1000];         %s,   nodes, strictly increasing, phase-local
         descAlpha = [0 0];            %deg, angle of attack; see glideAlpha
          descBank = [0 0];            %deg, bank [-90 .. 90]. SHIPPED AT ZERO, and that is a
                                       %     targeting decision, not an aerodynamic one.
                                       %     run_boost_glide dives at 75 deg for a steep, fast
                                       %     arrival, but a banked segment turns the heading --
                                       %     166.4 deg of it between handoff and impact on the
                                       %     SHIPPED 20N 155W to 35N 120W geometry, measured --
                                       %     and walks the impact point off the
                                       %     launch-to-target great circle. On that same
                                       %     geometry it is 21.52 km of miss against a 0.59 km
                                       %     range residual: the range solve converges and the
                                       %     vehicle still misses. Raise it only if you
                                       %     want the terminal dive more than you want the
                                       %     target, and read the cross-range line when you do
             hStop = 0;                %km, impact altitude, DESCENDING crossing only
                                       %    [0 .. hHandoff); a target on high terrain is a
                                       %    legitimate setting, not just sea level
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
         earthSpin = false;            %MUST BE false for the closed-form azimuth to hold.
                                       %true enables the Coriolis and centrifugal terms and
                                       %INVALIDATES the targeting solve; the run then prints a
                                       %caution and reports the miss against a target that the
                                       %rotating Earth has moved out from under it

%% Output:
         showPlots = true;             %false to skip the figures, e.g. under matlab -batch
          movieOn  = false;            %true renders the globe movie; it is the expensive part
                                       %of a run, roughly 0.3 s a frame, so it is off by default
       movieFrames = 120;              %-, frames to render [30 .. 600]; ignored unless movieOn
         movieFile = fullfile(tempdir,'run_target.mp4');   %char, MP4 output path. Defaults
                                       %into tempdir on purpose: a movie is a build artefact
                                       %and does not belong in the source tree
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
       overridable = {'latLaunch','lonLaunch','hLaunch','vLaunch', ...
                      'latTarget','lonTarget','cutFracMin','cutFracMax', ...
                      'tolRangeKm','pitchTime','pitchAngle','alphaMax', ...
                      'bankBoost','tMaxBoost','glideTime','glideAlpha', ...
                      'glideBank','hHandoff','tMaxGlide','descTime', ...
                      'descAlpha','descBank','hStop','tMaxDesc','vehicleFn', ...
                      'boosterFn','atmosFn','gravFn','aeroFn','propFn', ...
                      'earthSpin','showPlots','movieOn','movieFrames', ...
                      'movieFile','altExag'};
             given = fieldnames(opts);
    for ko = 1:numel(given)
        assert(any(strcmp(given{ko},overridable)), ...
            '"%s" is not a USER PARAMETERS entry; there is nothing to override.', ...
            given{ko});
    end
         latLaunch = overrideOf(opts,'latLaunch',latLaunch);
         lonLaunch = overrideOf(opts,'lonLaunch',lonLaunch);
           hLaunch = overrideOf(opts,'hLaunch',hLaunch);
           vLaunch = overrideOf(opts,'vLaunch',vLaunch);
         latTarget = overrideOf(opts,'latTarget',latTarget);
         lonTarget = overrideOf(opts,'lonTarget',lonTarget);
        cutFracMin = overrideOf(opts,'cutFracMin',cutFracMin);
        cutFracMax = overrideOf(opts,'cutFracMax',cutFracMax);
        tolRangeKm = overrideOf(opts,'tolRangeKm',tolRangeKm);
         pitchTime = overrideOf(opts,'pitchTime',pitchTime);
        pitchAngle = overrideOf(opts,'pitchAngle',pitchAngle);
          alphaMax = overrideOf(opts,'alphaMax',alphaMax);
         bankBoost = overrideOf(opts,'bankBoost',bankBoost);
         tMaxBoost = overrideOf(opts,'tMaxBoost',tMaxBoost);
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
          movieOn  = overrideOf(opts,'movieOn',movieOn);
       movieFrames = overrideOf(opts,'movieFrames',movieFrames);
         movieFile = overrideOf(opts,'movieFile',movieFile);

%% Convert the user block to library SI units. This is the ONLY unit
%% conversion in the file; everything past this point is m, m/s, rad and s:
          hLaunchM = hLaunch.*1000;
         hHandoffM = hHandoff.*1000;
            hStopM = hStop.*1000;
           tolRngM = tolRangeKm.*1000;
        latLaunchR = deg2rad(latLaunch);
        lonLaunchR = deg2rad(lonLaunch);
        latTargetR = deg2rad(latTarget);
        lonTargetR = deg2rad(lonTarget);
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
    assert(hStopM >= 0, ...
        ['hStop (%.1f km) is below the spherical datum; the descent ends on a ' ...
         'descending crossing of that altitude and a negative one has no ' ...
         'meaning against this Earth model.'],hStop);
    assert(vLaunch > 1, ...
        'vLaunch must exceed 1 m/s; the equations of motion are singular below that.');
    assert(all([tMaxBoost tMaxGlide tMaxDesc] > 0),'every phase horizon must be positive.');
    assert(cutFracMin > 0 && cutFracMin < cutFracMax && cutFracMax <= 1, ...
        ['the cutoff bracket must satisfy 0 < cutFracMin < cutFracMax <= 1; ' ...
         'got %.4f and %.4f. A fraction above 1 cannot be flown -- the burnout ' ...
         'event stops the boost at propellant exhaustion whatever tspan says.'], ...
        cutFracMin,cutFracMax);
    assert(tolRngM > 0,'tolRangeKm must be positive.');
    assert(isscalar(latLaunch) && isscalar(lonLaunch) && ...
           isscalar(latTarget) && isscalar(lonTarget), ...
        ['the launch point and the target must each be ONE latitude and ONE ' ...
         'longitude; this script solves a single trajectory, and a vector ' ...
         'here would silently broadcast through every great-circle call.']);
    assert(abs(latLaunch - latTarget) + abs(lonLaunch - lonTarget) > 0, ...
        'the launch point and the target are the same point; there is nothing to solve.');

                 c = coorbital.util.missileConst();
               veh = vehicleFn();
               bst = boosterFn();

%% Mass bookkeeping. The state mass is ALWAYS the total mass carried. At
%% separation the WHOLE booster goes -- dry structure and whatever propellant
%% an early cutoff left in it -- so the unpowered vehicle is exactly the
%% payload and its mass does not depend on where the cutoff landed:
          mLiftoff = veh.mass + bst.massDry + bst.massProp;
         mBurnoutT = veh.mass + bst.massDry;
            mGlide = veh.mass;

%% Full burn time, from the motor model itself rather than a literal, so a
%% different booster or a different propFn moves it. Vacuum flow is used
%% because coorbital.prop.constThrust chokes the flow and does not throttle it
%% with back pressure; a propFn that DOES vary its flow makes this a nominal
%% figure and the cutoff fractions nominal with it:
         [~,mdot0] = propFn(0,0,bst);
    assert(mdot0 > 0, ...
        'the propulsion model reports zero mass flow; there is no burn to cut short.');
             tBurn = bst.massProp./mdot0;
    assert(tMaxBoost > tBurn, ...
        ['tMaxBoost (%.1f s) must exceed the %.3f s full burn, or the horizon ' ...
         'and not the cutoff sets the top of the bracket.'],tMaxBoost,tBurn);

%% Assemble the environment from the handles chosen above:
         env.atmos = atmosFn;
          env.grav = gravFn;
          env.aero = aeroFn;
          env.prop = propFn;
        env.omegaE = 0;
    if earthSpin
        env.omegaE = c.omegaE;
    end

%% The unpowered vehicle, flown by BOTH the glide and the descent. Mass AND
%% aerodynamics come from the payload, because the booster is gone. Rebuilding
%% the mass alone would leave Sref, CL and LD describing a jettisoned stack,
%% which is the worse bug because it looks repaired:
          glideVeh = veh;
     glideVeh.mass = mGlide;

%% Guidance, one schedule per phase, each written against ITS OWN phase clock:
        schedBoost = struct('tGrid',pitchTimS, ...
                            'theta',pitchAngR, ...
                            'sigma',bankBoostR.*ones(1,numel(pitchTimS)), ...
                            'alphaMax',alphaMaxR);
        schedGlide = struct('tGrid',glideTimS,'alpha',glideAlfR,'sigma',glideBnkR);
         schedDesc = struct('tGrid',descTimS ,'alpha',descAlfR ,'sigma',descBnkR);

%% Equations of motion, each with ITS OWN vehicle bound in. The vehArg
%% argument is the struct phaseRun forwards; it is deliberately not used,
%% because a single chain-wide vehicle cannot describe both the stack and the
%% separated payload -- and coorbital.eom.massConstant raises
%% coorbital:massConstant:massMismatch the moment it is got wrong:
           eomFree = coorbital.eom.massConstant(@coorbital.eom.glide3DOF);
            eomBst = @(t,x,u,vehArg,envArg) coorbital.eom.boost3DOF(t,x,u,bst,envArg);
            eomGld = @(t,x,u,vehArg,envArg) eomFree(t,x,u,glideVeh,envArg);

%% -----------------------------------------------------------------
%% The azimuth solve: closed form, one call, no iteration
%% -----------------------------------------------------------------
%% The initial bearing of the launch-to-target great circle IS the library's
%% heading state psi -- both are measured clockwise from north -- so it goes
%% straight into the launch state below with no conversion:
         psiLaunch = coorbital.util.greatCircleBearing(latLaunchR,lonLaunchR, ...
                                                       latTargetR,lonTargetR);

%% Liftoff state. Everything but the mass is fixed; the pad flight-path angle
%% is the first commanded pitch attitude, so the angle of attack starts at
%% zero:
                x0 = [c.rE + hLaunchM; ...
                      lonLaunchR; ...
                      latLaunchR; ...
                      vLaunch; ...
                      pitchAngR(1); ...
                      psiLaunch; ...
                      mLiftoff];

%% Everything the propagation needs, gathered once so the range function is a
%% single-argument handle over an immutable configuration:
       cfg.eomBst  = eomBst;
       cfg.eomGld  = eomGld;
      cfg.schedBst = schedBoost;
      cfg.schedGld = schedGlide;
      cfg.schedDsc = schedDesc;
      cfg.mBurnout = mBurnoutT;
        cfg.mGlide = mGlide;
      cfg.hHandoff = hHandoffM;
         cfg.hStop = hStopM;
     cfg.tMaxGlide = tMaxGlide;
      cfg.tMaxDesc = tMaxDesc;
            cfg.x0 = x0;
           cfg.bst = bst;
           cfg.env = env;
           cfg.rE  = c.rE;
          cfg.lat0 = latLaunchR;
          cfg.lon0 = lonLaunchR;

%% Every surface range in this file, required and achieved alike, is measured
%% on the IMPACT sphere r = rE + hStop, so the target distance and the flown
%% distance share one baseline and their difference is a real miss and not a
%% radius mismatch:
                rI = c.rE + hStopM;
           angReqR = coorbital.util.greatCircle(latLaunchR,lonLaunchR, ...
                                                latTargetR,lonTargetR);
            rngReq = rI.*angReqR;

%% -----------------------------------------------------------------
%% The range solve: bracket the envelope, then bisect
%% -----------------------------------------------------------------
%% The bracket endpoints in seconds. rangeSolve evaluates both before it
%% iterates and returns their achieved ranges in info.fMin and info.fMax,
%% which IS the reachable envelope -- so bracketing and reporting the envelope
%% are the same two propagations, not two pairs:
              tLo  = cutFracMin.*tBurn;
              tHi  = cutFracMax.*tBurn;
            fRange = @(tCut) flyRange(tCut,cfg,rI);
  [tCut,rngAch,sv] = coorbital.util.rangeSolve(rngReq,fRange,tLo,tHi,tolRngM);

%% Fractions of full burn, for the report. The cutoff as a bare number of
%% seconds says nothing without the burn it is a fraction of:
          fracSol  = tCut./tBurn;
          fracLo   = tLo./tBurn;
          fracHi   = tHi./tBurn;

%% -----------------------------------------------------------------
%% Refuse, loudly, when the target is outside the reachable envelope
%% -----------------------------------------------------------------
%% coorbital.util.rangeSolve returns converged = false with the achievable
%% band rather than throwing, precisely so the refusal can be phrased here in
%% the problem's own vocabulary -- kilometres, a launch point, a target. A
%% targeting script that silently returned the nearest miss would be worse
%% than one that refuses, because the nearest miss LOOKS like a solution:
    if ~sv.converged

%% Signed distance to the nearer edge of the band, computed ONCE: it is both
%% the number printed and the test that chooses which half of the advice below
%% applies, and computing it twice would let the two disagree:
            shortM = bandShortfall(rngReq,sv.fMin,sv.fMax);
        fprintf('\n');
        fprintf('===== Point-to-point targeting: REFUSED =================\n');
        fprintf('  No cutoff time in the bracket puts the vehicle on the target.\n');
        fprintf('\n');
        fprintf('    launch           %10.4f %-5s %10.4f deg  (lat, lon)\n', ...
                latLaunch,'deg',lonLaunch);
        fprintf('    target           %10.4f %-5s %10.4f deg  (lat, lon)\n', ...
                latTarget,'deg',lonTarget);
        fprintf('    required range   %10.2f %-5s (great circle on the r = %.3f km sphere)\n', ...
                rngReq./1000,'km',rI./1000);
        fprintf('    reachable        %10.2f to %.2f km\n',sv.fMin./1000,sv.fMax./1000);
        fprintf('    cutoff bracket   %10.3f to %.3f s  (%.3f to %.3f of the %.3f s full burn)\n', ...
                tLo,tHi,fracLo,fracHi,tBurn);
        fprintf('    shortfall        %10.2f %-5s (required minus the nearer edge of the band)\n', ...
                shortM./1000,'km');
        fprintf('    solver said      %s\n',sv.identifier);
        fprintf('\n');
        if shortM > 0
        fprintf('  The target is TOO FAR for this vehicle and this bracket. Either raise\n');
        fprintf('  cutFracMax towards 1.0 if it is below it (it is %.3f), give the booster\n',cutFracMax);
        fprintf('  more propellant, loft the pitch program, or choose a nearer target. At\n');
        fprintf('  cutFracMax = %.3f the vehicle flies %.2f km and no cutoff can add to that:\n', ...
                cutFracMax,sv.fMax./1000);
        fprintf('  cutting the burn SHORT is the only control this script has, and it can\n');
        fprintf('  only ever subtract range.\n');
        else
        fprintf('  The target is TOO CLOSE. Even the shortest burn in the bracket overflies\n');
        fprintf('  it by %.2f km. Lower cutFracMin below the shipped %.3f -- but NOT BELOW THE\n', ...
                (sv.fMin - rngReq)./1000,cutFracMin);
        fprintf('  DOCUMENTED 0.3 FLOOR, which is not a style preference: at 0.20 the boost\n');
        fprintf('  leaves the vehicle too low and too slow to reach the handoff at all, the\n');
        fprintf('  glide runs its whole horizon descending, and the propagation dies inside\n');
        fprintf('  the equations of motion with coorbital:glide3DOF:zeroSpeed at a negative\n');
        fprintf('  speed rather than reaching the propagationIncomplete refusal this script\n');
        fprintf('  would otherwise raise. 0.30, 0.35 and 0.45 are all sound.\n');
        fprintf('  Or accept that a boost-glide vehicle has a minimum range: below roughly\n');
        fprintf('  %.0f km the boost cannot lift it high enough for the glide phase to start,\n', ...
                sv.fMin./1000);
        fprintf('  the handoff event never fires, and the propagation is rejected rather\n');
        fprintf('  than believed.\n');
        end
        fprintf('=========================================================\n\n');
              traj = [];
       info.refused  = true;
       info.rngReqM  = rngReq;
     info.psiLaunch  = psiLaunch;
       info.rngMinM  = sv.fMin;
       info.rngMaxM  = sv.fMax;
         info.tBurn  = tBurn;
        info.tCutLo  = tLo;
        info.tCutHi  = tHi;
     info.solveInfo  = sv;
        if nargout == 0
            clear traj info;
        end
        return;
    end

%% -----------------------------------------------------------------
%% Fly the solved trajectory
%% -----------------------------------------------------------------
%% One more propagation at the converged cutoff. rangeSolve already flew this
%% abscissa but returns only the scalar range, so the state history is
%% re-created here rather than carried out through a persistent variable:
     [rngFly,traj] = flyRange(tCut,cfg,rI);
    assert(abs(rngFly - rngAch) < 1e-6, ...
        ['re-flying the solved cutoff gave %.6f m of range against the ' ...
         'solver''s %.6f m; the propagation is not repeatable and nothing ' ...
         'below can be trusted'],rngFly,rngAch);

%% Phase boundaries. phaseRun records the boundary sample ONCE, carrying the
%% outgoing phase's control, so the last row of phase k is that phase's
%% terminal state on the NEAR side of any staging jump:
               kBO = find(traj.phaseIdx == 1,1,'last');
               kHO = find(traj.phaseIdx == 2,1,'last');
                nS = numel(traj.t);

%% -----------------------------------------------------------------
%% Independent check of the two solves, from the FLOWN state
%% -----------------------------------------------------------------
%% The miss is measured impact-to-target, not read back out of the solver's
%% own bookkeeping. With zero bank the two agree, and the agreement is the
%% evidence for the zero-cross-range claim rather than a restatement of it:
             latIm = traj.x(end,3);
             lonIm = traj.x(end,2);
             missM = rI.*coorbital.util.greatCircle(latIm,lonIm, ...
                                                    latTargetR,lonTargetR);
              resM = rngAch - rngReq;

%% Cross-track offset of the impact point from the launch-to-target great
%% circle, positive to the right of the intended course. The standard
%% cross-track formula: the sine of the offset angle is the sine of the flown
%% central angle times the sine of the bearing error at the launch point:
            angFly = coorbital.util.greatCircle(latLaunchR,lonLaunchR,latIm,lonIm);
            psiFly = coorbital.util.greatCircleBearing(latLaunchR,lonLaunchR, ...
                                                       latIm,lonIm);
            dPsiR  = wrapPi(psiFly - psiLaunch);
           xTrackM = rI.*asin(sin(angFly).*sin(dPsiR));

%% The heading the vehicle actually left the pad on must be the bearing the
%% closed form asked for, to machine precision -- x0(6) is written from
%% psiLaunch above and ode45 never touches the initial condition, so any
%% disagreement here means the state was assembled wrong:
           psiSeed = traj.x(1,6);
    assert(abs(wrapPi(psiSeed - psiLaunch)) < 1e-12, ...
        ['the flown initial heading %.15f rad is not the solved launch ' ...
         'azimuth %.15f rad; the launch state was assembled wrong'], ...
        psiSeed,psiLaunch);

%% Down-track and cross-track should account for the whole miss when both are
%% small compared with the Earth. Reported, not asserted: it is a spherical
%% right-triangle approximation and it degrades at long offsets, so it is
%% shown as a consistency figure the reader can judge:
          missHypM = hypot(resM,xTrackM);
        crossWarn  = abs(xTrackM) > tolRngM;

%% Propellant actually burned, and the propellant thrown away still in the
%% booster. Measured from the flown mass state rather than from mdot times the
%% cutoff, so a propFn whose flow is not constant is still reported correctly.
%% At a deep cutoff the second number is the bigger one, and a summary that
%% reported only the first would hide the whole cost of solving the range this
%% way:
          mpBurned = mLiftoff - traj.x(kBO,7);
          mpWasted = bst.massProp - mpBurned;

%% -----------------------------------------------------------------
%% Derived quantities for the trajectory summary
%% -----------------------------------------------------------------
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
%% is not sensed and is excluded. nSens carries the thrust; nAero drops it.
%% nAero IS THE AERODYNAMIC LOAD FACTOR, lift and drag together, and it is NOT
%% a deceleration: with constLD the drag is only CD/CL = 1/LD of it:
             nSens = sqrt((aThrV - aDrag).^2 + (aThrN + aLift).^2)./c.g0;
             nAero = sqrt(aLift.^2 + aDrag.^2)./c.g0;

%% The mass carried through the unpowered flight must be exactly the payload,
%% because the six-state glide equations divide by glideVeh.mass and never
%% read the mass state. massConstant already refuses any disagreement inside
%% the derivative; this is the same claim checked once at the end, where the
%% failure message can name the staging link that caused it:
             isUnp = traj.phaseIdx >= 2;
             mSpan = max(massS(isUnp)) - min(massS(isUnp));
    assert(mSpan == 0, ...
        ['the unpowered mass varied by %.3e kg; dm/dt is identically zero ' ...
         'after separation, so it must be bit-exactly constant'],mSpan);
    assert(abs(massS(end) - mGlide) < 1e-6, ...
        ['the unpowered flight is carrying %.9f kg while its vehicle was ' ...
         'built around %.9f kg; the separation link is wrong'],massS(end),mGlide);

%% Masks and peaks, taken over the phase where each is meaningful:
             isBst = traj.phaseIdx == 1;
             isGld = traj.phaseIdx == 2;
             isDsc = traj.phaseIdx == 3;
     [qBstMax,kQB] = maxOver(qbar ,isBst);
     [nBstMax,kNB] = maxOver(nSens,isBst);
     [qGldMax,kQG] = maxOver(qbar ,isGld);
     [qDscMax,kQD] = maxOver(qbar ,isDsc);
     [nAerMax,kNA] = maxOver(nAero,isUnp);
       [hGlMax,kG] = maxOver(hKm  ,isGld);

%% Leg ranges, each its own great circle on the impact sphere:
             angBO = coorbital.util.greatCircle(latLaunchR,lonLaunchR, ...
                                                traj.x(kBO,3),traj.x(kBO,2));
          downBOKm = rI.*angBO./1000;
             angGl = coorbital.util.greatCircle(traj.x(kBO,3),traj.x(kBO,2), ...
                                                traj.x(kHO,3),traj.x(kHO,2));
           glideKm = rI.*angGl./1000;
             angDs = coorbital.util.greatCircle(traj.x(kHO,3),traj.x(kHO,2), ...
                                                latIm,lonIm);
            descKm = rI.*angDs./1000;

%% Termination diagnosis, one line per phase. The boost line is different in
%% THIS script from every other: a boost that stops at the horizon is a
%% failure elsewhere and is the whole POINT here, so the nominal case is the
%% commanded cutoff and propellant exhaustion is the special case:
     [why1,ok1] = whyCutoff(traj.t(kBO),massS(kBO),mBurnoutT,tCut);
     [why2,ok2] = whyAltitude(traj.t(kHO) - traj.t(kBO),traj.t(kHO), ...
                              traj.x(kHO,1) - c.rE,hHandoffM,hHandoff, ...
                              tMaxGlide,'handoff');
     [why3,ok3] = whyAltitude(traj.t(nS) - traj.t(kHO),traj.t(nS), ...
                              traj.x(nS,1) - c.rE,hStopM,hStop, ...
                              tMaxDesc,'impact');
          allOK = ok1 && ok2 && ok3;

%% How far the ground under the target would have moved during this flight had
%% the Earth been turning. Quantifies the non-rotating caveat instead of
%% merely stating it: this is the size of the error the outer azimuth
%% iteration that does not exist here would have to remove:
            driftM = c.omegaE.*rI.*cos(latTargetR).*traj.t(end);

%% Whether ANY phase commands a bank. The zero-bank case and the banked case
%% need two different paragraphs in the limitations block below, not one
%% paragraph with a number substituted into it: at 75 deg of descent bank the
%% sentence "cross-range came out at 21515.22 m because the bank angle is zero
%% throughout" is self-refuting, and it sat directly under a WARNING saying the
%% opposite. Gated on the RADIAN schedules, so a user block written in any
%% units that convert to zero is still recognised as zero:
           anyBank = any([bankBoostR glideBnkR descBnkR] ~= 0);

%% Labels for the switches in the model line:
           spinTxt = 'OFF';
    if earthSpin
           spinTxt = 'ON';
    end

%% Conventional lower edge of the hypersonic regime. A modelling convention,
%% not a physical constant of the Earth or the air, so it does not belong in
%% missileConst; it is the threshold below which holding CL and L/D constant
%% stops being defensible:
         machHyper = 5;

%% Which of the swappable models are still the library defaults. The validity
%% paragraph below describes the DEFAULT models by name, and every one of
%% those statements becomes a FALSEHOOD the moment the handle is replaced:
      defaultAtmos = isequal(atmosFn,@coorbital.atmos.expAtmos);
       defaultGrav = isequal(gravFn,@coorbital.grav.sphereGrav);
       defaultProp = isequal(propFn,@coorbital.prop.constThrust);

%% Report:
    fprintf('\n');
    fprintf('===== Point-to-point targeting summary ==================\n');
    fprintf('  Phase termination\n');
    fprintf('    1 boost          %s\n',why1);
    fprintf('    2 glide          %s\n',why2);
    fprintf('    3 descent        %s\n',why3);
    if ~allOK
        fprintf('  *** CAUTION ***  at least one phase did NOT end as intended; the\n');
        fprintf('                   trajectory below is TRUNCATED, not a completed flight\n');
    end
    fprintf('\n');
    fprintf('  Targeting solve\n');
    fprintf('    launch           %10.4f %-5s %10.4f deg   (latitude, longitude)\n', ...
            latLaunch,'deg',lonLaunch);
    fprintf('    target           %10.4f %-5s %10.4f deg   (latitude, longitude)\n', ...
            latTarget,'deg',lonTarget);
    fprintf('    required range   %10.3f %-5s  (great circle, launch to target, on the\n', ...
            rngReq./1000,'km');
    fprintf('                                        r = %.3f km impact sphere)\n',rI./1000);
    fprintf('    achieved range   %10.3f %-5s  (great circle, launch to the flown impact)\n', ...
            rngAch./1000,'km');
    fprintf('    range residual   %+10.2f %-5s  (achieved minus required, signed)\n',resM,'m');
    fprintf('    MISS DISTANCE    %10.2f %-5s  (great circle, impact point to target,\n', ...
            missM,'m');
    fprintf('                                        computed from the FLOWN state and not\n');
    fprintf('                                        from the solver''s own bookkeeping)\n');
    fprintf('    tolerance        %10.2f %-5s  (requested on the achieved range)\n', ...
            tolRngM,'m');
    fprintf('\n');
    fprintf('    launch azimuth   %10.6f %-5s  (clockwise from north, CLOSED FORM: the\n', ...
            rad2deg(psiLaunch),'deg');
    fprintf('                                        initial bearing of the launch-to-target\n');
    fprintf('                                        great circle. No iteration)\n');
    fprintf('    flown azimuth    %10.6f %-5s  (bearing launch to impact, from the flown\n', ...
            rad2deg(psiFly),'deg');
    fprintf('                                        state; %+.3e deg from the commanded one)\n', ...
            rad2deg(dPsiR));
    fprintf('    cross-track      %+10.2f %-5s  (%+.3e m; offset of the impact point from the\n', ...
            xTrackM,'m',xTrackM);
    fprintf('                                        launch-to-target great circle, positive\n');
    fprintf('                                        right of course. THE SOLVE DOES NOT\n');
    fprintf('                                        CONTROL THIS -- see the limitations)\n');
    fprintf('    residual + cross %10.2f %-5s  (hypotenuse of the two above, against the\n', ...
            missHypM,'m');
    fprintf('                                        %.2f m miss measured directly. They\n',missM);
    fprintf('                                        agree to %.2e m)\n',abs(missHypM - missM));
    if crossWarn
        fprintf('  *** WARNING ***  the impact point is %.2f km off the launch-to-target great\n', ...
                abs(xTrackM)./1000);
        fprintf('                   circle, more than the %.2f km range tolerance. The range\n', ...
                tolRngM./1000);
        fprintf('                   solve CONVERGED and the vehicle still missed, because a\n');
        fprintf('                   banked segment turns the heading and walks the track off\n');
        fprintf('                   the arc it left on. Set glideBank and descBank to zero, or\n');
        fprintf('                   accept the cross-range: closing it needs an outer azimuth\n');
        fprintf('                   iteration, which this script does not have.\n');
    end
    fprintf('\n');
    fprintf('    solved cutoff    %10.4f %-5s  (%.6f of the %.4f s full burn)\n', ...
            tCut,'s',fracSol,tBurn);
    fprintf('    propellant burned%10.1f %-5s  (of %.1f kg loaded)\n', ...
            mpBurned,'kg',bst.massProp);
    fprintf('    thrown away      %10.1f %-5s  (unburned, jettisoned inside the booster at\n', ...
            mpWasted,'kg');
    fprintf('                                        separation; %.1f %% of the load)\n', ...
            100.*mpWasted./bst.massProp);
    fprintf('    reachable        %10.3f to %.3f km  (the bracket endpoints, flown)\n', ...
            sv.fMin./1000,sv.fMax./1000);
    fprintf('    cutoff bracket   %10.4f to %.4f s   (%.3f to %.3f of full burn)\n', ...
            tLo,tHi,fracLo,fracHi);
    fprintf('    iterations       %10d %-5s  (bisection steps; %d trajectory propagations\n', ...
            sv.iterations,'',sv.nEval);
    fprintf('                                        in the solve, plus one to re-fly it)\n');
    fprintf('\n');
    fprintf('  Launch\n');
    fprintf('    pad altitude     %10.3f km\n',hLaunch);
    fprintf('    initial speed    %10.2f %-5s  (integration start, not a physical launch speed)\n', ...
            vLaunch,'m/s');
    fprintf('    liftoff mass     %10.1f %-5s  (payload %.1f + booster dry %.1f + propellant %.1f)\n', ...
            mLiftoff,'kg',veh.mass,bst.massDry,bst.massProp);
    fprintf('    liftoff T/W      %10.3f %-5s  (vacuum thrust over liftoff weight)\n', ...
            bst.thrustVac./(mLiftoff.*c.g0),'');
    fprintf('\n');
    fprintf('  Cutoff and separation, end of phase 1\n');
    fprintf('    time             %10.4f %-5s  (the SOLVED cutoff, not propellant exhaustion)\n', ...
            traj.t(kBO),'s');
    fprintf('    altitude         %10.3f km\n',hKm(kBO));
    fprintf('    speed            %10.2f %-5s  (Mach %.2f, planet-relative)\n', ...
            V(kBO),'m/s',mach(kBO));
    fprintf('    flight path      %10.3f %-5s  (positive is climbing)\n', ...
            rad2deg(traj.x(kBO,5)),'deg');
    fprintf('    heading          %10.3f deg\n',rad2deg(traj.x(kBO,6)));
    fprintf('    downrange        %10.2f %-5s  (great circle from the pad)\n',downBOKm,'km');
    fprintf('    mass at cutoff   %10.1f %-5s  (payload + booster + unburned propellant)\n', ...
            massS(kBO),'kg');
    fprintf('    mass into glide  %10.1f %-5s  (the payload alone; the WHOLE booster is\n', ...
            mGlide,'kg');
    fprintf('                                        jettisoned, unburned propellant included)\n');
    fprintf('\n');
    fprintf('  Glide, phase 2, from separation to the handoff\n');
    fprintf('    duration         %10.2f %-5s  (%.2f min)\n', ...
            traj.t(kHO) - traj.t(kBO),'s',(traj.t(kHO) - traj.t(kBO))./60);
    fprintf('    ground covered   %10.2f %-5s  (great circle, cutoff point to handoff point)\n', ...
            glideKm,'km');
    fprintf('    peak altitude    %10.2f %-5s  (at t = %.1f s; the first skip, thrown by the\n', ...
            hGlMax,'km',traj.t(kG));
    fprintf('                                        pull-up as the vehicle meets the air)\n');
    fprintf('    handoff time     %10.3f s\n',traj.t(kHO));
    fprintf('    handoff altitude %10.4f %-5s  (%+.2e m from the %.3f km target, event residual)\n', ...
            hKm(kHO),'km',traj.x(kHO,1) - c.rE - hHandoffM,hHandoff);
    fprintf('    handoff speed    %10.2f %-5s  (Mach %.2f)\n',V(kHO),'m/s',mach(kHO));
    fprintf('\n');
    fprintf('  Descent, phase 3, from the handoff to impact\n');
    fprintf('    duration         %10.2f s\n',traj.t(nS) - traj.t(kHO));
    fprintf('    ground covered   %10.2f %-5s  (great circle, handoff point to impact)\n', ...
            descKm,'km');
    fprintf('    commanded bank   %10.2f %-5s  (glide bank %.2f deg)\n', ...
            rad2deg(descBnkR(1)),'deg',rad2deg(glideBnkR(end)));
    if descBnkR(1) == glideBnkR(end)
    fprintf('                                        The two banks are EQUAL, so the handoff is\n');
    fprintf('                                        pure bookkeeping: both phases fly one\n');
    fprintf('                                        control and hHandoff cannot move the range.\n');
    else
    fprintf('                                        The two banks DIFFER, so the handoff is a\n');
    fprintf('                                        real control discontinuity and hHandoff can\n');
    fprintf('                                        move where the vehicle lands.\n');
    end
    fprintf('\n');
    fprintf('  Overall\n');
    fprintf('    flight time      %10.2f %-5s  (%.2f min)\n',traj.t(end),'s',traj.t(end)./60);
    fprintf('    ground range     %10.2f %-5s  (great circle on the r = %.1f km impact sphere)\n', ...
            rngAch./1000,'km',rI./1000);
    fprintf('    central angle    %10.4f deg\n',rad2deg(angFly));
    fprintf('    leg sum          %10.2f %-5s  = %.2f (boost) + %.2f (glide) + %.2f (descent).\n', ...
            downBOKm + glideKm + descKm,'km',downBOKm,glideKm,descKm);
    fprintf('                                        The legs need NOT add to the total: each is\n');
    fprintf('                                        its own great circle.\n');
    fprintf('    samples          %10d %-5s  (ode45 adaptive steps over 3 phases)\n',nS,'');
    fprintf('\n');
    fprintf('  Impact\n');
    fprintf('    altitude         %10.6f %-5s  (%+.2e m from the %.3f km stop, event residual)\n', ...
            hKm(end),'km',traj.x(end,1) - c.rE - hStopM,hStop);
    fprintf('    speed            %10.2f %-5s  (Mach %.2f)\n',V(end),'m/s',mach(end));
    fprintf('    flight path      %10.3f %-5s  (negative is descending)\n', ...
            rad2deg(traj.x(end,5)),'deg');
    fprintf('    heading          %10.3f deg\n',rad2deg(traj.x(end,6)));
    fprintf('    latitude         %11.6f %-4s  (target %.6f deg)\n', ...
            rad2deg(latIm),'deg',latTarget);
    fprintf('    longitude        %11.6f %-4s  (target %.6f deg)\n', ...
            rad2deg(lonIm),'deg',lonTarget);
    fprintf('\n');
    fprintf('  Peak loads\n');
    fprintf('    boost max q      %10.2f %-5s  (at t = %.1f s, h = %.2f km)\n', ...
            qBstMax./1000,'kPa',traj.t(kQB),hKm(kQB));
    fprintf('    boost sensed     %10.2f %-5s  (%.2f m/s^2, thrust included, at t = %.1f s)\n', ...
            nBstMax,'g',nBstMax.*c.g0,traj.t(kNB));
    fprintf('    glide max q      %10.2f %-5s  (at t = %.1f s, h = %.2f km)\n', ...
            qGldMax./1000,'kPa',traj.t(kQG),hKm(kQG));
    fprintf('    descent max q    %10.2f %-5s  (at t = %.1f s, h = %.2f km)\n', ...
            qDscMax./1000,'kPa',traj.t(kQD),hKm(kQD));
    fprintf('    peak aero load   %10.2f %-5s  (%.2f m/s^2 lift and drag, thrust excluded, at\n', ...
            nAerMax,'g',nAerMax.*c.g0);
    fprintf('                                        t = %.1f s, h = %.2f km, in phase %d).\n', ...
            traj.t(kNA),hKm(kNA),traj.phaseIdx(kNA));
    fprintf('                                        A LOAD FACTOR, what the structure feels,\n');
    fprintf('                                        and NOT a deceleration: %s\n',func2str(aeroFn));
    fprintf('                                        fixes CD/CL = 1/(L/D) = %.3f, so only\n',1./veh.LD);
    fprintf('                                        that fraction of it brakes the vehicle.\n');
    fprintf('\n');
    fprintf('  LIMITATIONS OF THIS TARGETING SOLUTION\n');
    fprintf('    1. THE AZIMUTH IS EXACT ONLY FOR A NON-ROTATING EARTH. This run has Earth\n');
    fprintf('       rotation %s, env.omegaE = %.6e rad/s, so the ground track of a zero-bank\n', ...
            spinTxt,env.omegaE);
    fprintf('       trajectory is a great circle and the initial bearing of the\n');
    fprintf('       launch-to-target arc is the WHOLE answer -- one closed-form call, no\n');
    fprintf('       iteration. Turn the Earth on and it stops being the answer: over this\n');
    fprintf('       %.0f s flight the ground beneath the target sweeps %.0f km east at\n', ...
            traj.t(end),driftM./1000);
    fprintf('       latitude %.1f deg, so the vehicle would have to be aimed where the target\n', ...
            latTarget);
    fprintf('       is GOING TO BE. That makes the azimuth depend on the flight time, the\n');
    fprintf('       flight time depend on the cutoff, and the cutoff depend on the azimuth:\n');
    fprintf('       an OUTER ITERATION around the range solve. This script has none, and\n');
    fprintf('       rotating-Earth targeting is out of scope for it.\n');
    if earthSpin
    fprintf('       *** earthSpin IS ON. The azimuth above was computed as though it were\n');
    fprintf('       not, so the miss reported above is against a target that stood still and\n');
    fprintf('       is NOT this vehicle''s miss. Do not read it as one. ***\n');
    end
    fprintf('    2. THE MISS IS THE RESIDUAL OF THE RANGE SOLVE, ALONG THE GREAT CIRCLE.\n');
    fprintf('       Bisection matches a DISTANCE; nothing in it steers sideways.\n');
    if anyBank
    fprintf('       THIS RUN COMMANDS A NON-ZERO BANK -- boost %.1f deg, glide %.1f deg,\n', ...
            rad2deg(bankBoostR),rad2deg(glideBnkR(end)));
    fprintf('       descent %.1f deg -- so the ground track is NOT the launch-to-target great\n', ...
            rad2deg(descBnkR(1)));
    fprintf('       circle, and the range solve cannot see the difference. Measured on this\n');
    fprintf('       run: %.2f m of cross-track offset, a %.2f m total miss, against a range\n', ...
            abs(xTrackM),missM);
    fprintf('       residual of only %.2f m. The solve converged on DOWNRANGE and the vehicle\n', ...
            abs(resM));
    fprintf('       still missed, by a factor of %.1f over the residual. Closing that needs an\n', ...
            missM./max(abs(resM),eps));
    fprintf('       OUTER AZIMUTH ITERATION, which this script does not have. Set every bank\n');
    fprintf('       to zero, or read the range residual as a downrange-only figure and the\n');
    fprintf('       MISS DISTANCE line above as the answer.\n');
    else
    fprintf('       Cross-range came out at %.2f m here because the bank angle is zero\n', ...
            abs(xTrackM));
    fprintf('       throughout -- boost %.1f deg, glide %.1f deg, descent %.1f deg -- and a\n', ...
            rad2deg(bankBoostR),rad2deg(glideBnkR(end)),rad2deg(descBnkR(1)));
    fprintf('       zero-bank trajectory over a non-rotating sphere never leaves the great\n');
    fprintf('       circle it departed on. THAT IS A PROPERTY OF THIS CONFIGURATION, NOT A\n');
    fprintf('       GENERAL GUARANTEE. It is measured above, not assumed. Bank any phase and\n');
    fprintf('       the track curves off the arc while the range solve still converges: at\n');
    fprintf('       run_boost_glide''s 75 deg terminal bank this geometry misses by 21.52 km\n');
    fprintf('       on a 0.59 km range residual.\n');
    end
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
    fprintf('            terminal descent lie outside it.\n');
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
        fprintf('            Gravity is %s, NOT the default sphereGrav, so the\n', ...
                func2str(gravFn));
        fprintf('            no-J2 caveat this script prints by default does not apply. Its\n');
        fprintf('            validity range is its own to state.\n');
    end
    if ~defaultProp
        fprintf('            The motor is %s, NOT the default constant-thrust\n', ...
                func2str(propFn));
        fprintf('            model, so the %.4f s full burn this script solves fractions of\n',tBurn);
        fprintf('            was taken from its VACUUM mass flow and is nominal only.\n');
    end
    fprintf('            The impact point above is where a PLACEHOLDER vehicle lands under\n');
    fprintf('            PLACEHOLDER aerodynamics. Treat it as an indicative trajectory and a\n');
    fprintf('            demonstration of the targeting method, NOT as an accuracy claim.\n');
    if mach(end) < machHyper
        fprintf('            Impact is at Mach %.2f, below Mach %.0f, so the terminal segment\n', ...
                mach(end),machHyper);
        fprintf('            is outside the constant-coefficient regime that sets where it\n');
        fprintf('            lands.\n');
    end
    if nBstMax > 15 && defaultProp
        fprintf('            The %.1f g sensed load at cutoff is an artefact of the\n',nBstMax);
        fprintf('            constant-thrust placeholder motor, which does not throttle as the\n');
        fprintf('            stack empties. A real stage would tail off before that.\n');
    end
    fprintf('=========================================================\n\n');

%% Hand the summary's numbers back at full precision, so a test does not have
%% to read them out of printed text:
       info.refused  = false;
       info.rngReqM  = rngReq;
       info.rngAchM  = rngAch;
         info.missM  = missM;
        info.residM  = resM;
       info.xTrackM  = xTrackM;
      info.missHypM  = missHypM;
     info.crossWarn  = crossWarn;
     info.psiLaunch  = psiLaunch;
        info.psiFly  = psiFly;
         info.dPsiR  = dPsiR;
          info.tCut  = tCut;
       info.cutFrac  = fracSol;
         info.tBurn  = tBurn;
        info.tCutLo  = tLo;
        info.tCutHi  = tHi;
       info.rngMinM  = sv.fMin;
       info.rngMaxM  = sv.fMax;
    info.iterations  = sv.iterations;
         info.nEval  = sv.nEval;
     info.solveInfo  = sv;
      info.mpBurned  = mpBurned;
      info.mpWasted  = mpWasted;
      info.tBurnout  = traj.t(kBO);
      info.hBurnout  = traj.x(kBO,1) - c.rE;
      info.vBurnout  = V(kBO);
    info.gamBurnout  = traj.x(kBO,5);
      info.mBurnout  = massS(kBO);
        info.mGlide  = mGlide;
      info.mLiftoff  = mLiftoff;
      info.downBOKm  = downBOKm;
      info.tHandoff  = traj.t(kHO);
      info.hHandoff  = traj.x(kHO,1) - c.rE;
      info.vHandoff  = V(kHO);
     info.hGlideMax  = hGlMax;
       info.glideKm  = glideKm;
        info.descKm  = descKm;
       info.tFlight  = traj.t(end);
       info.rangeKm  = rngAch./1000;
        info.angTot  = angFly;
       info.vImpact  = V(end);
     info.gamImpact  = traj.x(end,5);
     info.latImpact  = latIm;
     info.lonImpact  = lonIm;
       info.machImp  = mach(end);
       info.qBstMax  = qBstMax;
       info.nBstMax  = nBstMax;
       info.qGldMax  = qGldMax;
       info.qDscMax  = qDscMax;
       info.nAerMax  = nAerMax;
        info.driftM  = driftM;
        info.stopOK  = allOK;
       info.stopWhy  = {why1;why2;why3};

%% ...and the machinery itself, so an independent checker can re-integrate the
%% very same chain with a different solver instead of trusting this one:
        info.phases  = buildPhases(tCut,cfg);
           info.env  = env;
            info.x0  = x0;
      info.boostVeh  = bst;
      info.glideVeh  = glideVeh;

%% Vertical exaggeration for the globe. A DISPLAY choice, not a physical
%% constant, so it does not belong in missileConst. A true-scale arc would lie
%% on the surface and show nothing of the skip phugoid -- a 50 km glide on a
%% 6378 km sphere is one part in 128.
%%
%% But a FIXED exaggeration does not travel. The shipped 30x suits a glide that
%% stays under about 60 km; on a lofted intercontinental shot that peaks near
%% 260 km it paints the arc 7800 km off the surface, further out than the Earth
%% is wide, and the track leaves the frame entirely. So it is scaled to the
%% flown apogee by exagFor below, and remains overridable.
%%
%% THE SHIPPED CASE DID CHANGE, and the comment this replaces claimed it did
%% not. This script's own default peaks at 118.09 km, so the rule returns 16x,
%% not the 30x it used to hard-code -- the cap never comes near binding here.
%% Sixteen is the better picture and not a regression: at 30x a 118 km apogee
%% is painted 3543 km off the surface, more than half an Earth radius, and the
%% arc reads as an orbit rather than as a glide; at 16x it stands 1889 km off,
%% which is the 0.3 rE the rule is aiming at. What IS unchanged is the other
%% three entry scripts, and only because each hard-codes its own factor and
%% never reaches this rule at all -- HGV/run_glide and HGV/run_boost_glide at
%% 30x, BM/run_ballistic at 3x:
           hPeakM = max(traj.x(:,1)) - c.rE;
          altExag = overrideOf(opts,'altExag',exagFor(hPeakM,c.rE));

%% The display choices go into info beside the flight's own numbers. The RULE
%% goes with them, as a handle: its cap and its floor sit at apogees no
%% flyable configuration of this script reaches -- under 63.8 km and over
%% 956.7 km against a 118 km shipped peak -- so a checker that can only observe
%% flown cases cannot see either of them work, and both were unpinned:
      info.altExag = altExag;
  info.altExagRule = @(hM) exagFor(hM,c.rE);

%% Plots. Every figure comes from coorbital.viz, which reads the trajectory and
%% never writes it. Nothing below this line can move a number in the summary
%% above. ONE VEHICLE PER PHASE is handed to the profile plot, because this
%% chain does not fly one: phase 1 is the stack under power and phases 2 and 3
%% are the separated payload. The ground track and the globe are given the
%% TARGET, so the aim point and the impact point are both on the picture and
%% the miss is something a reader can see rather than only read.
%%
%% The two captions are formatted ONCE, with hemisphere letters, and reported
%% in info: a caption is the only place most readers meet these coordinates,
%% and a caption that no test can read is a caption nothing checks:
     info.launchStr = llStr(latLaunch,lonLaunch);
     info.targetStr = llStr(latTarget,lonTarget);
    if showPlots
        coorbital.viz.profilePlot(traj,bst,env, ...
            struct('Name','Point-to-point targeting profile', ...
                   'Channels',{{'altitude','speed','mach','q','nAero', ...
                                'mass','gamma'}}, ...
                   'VehPhase',{{bst,glideVeh,glideVeh}}));
        coorbital.viz.groundTrack(traj,bst,env, ...
            struct('Target',[latTargetR; lonTargetR], ...
                   'PhaseName',{{'boost','glide','descent'}}, ...
                   'Title',sprintf(['Ground track, %.0f km required, %.0f m miss ' ...
                                    '(cutoff %.2f s)'],rngReq./1000,missM,tCut)));
        coorbital.viz.globe3D(traj,bst,env, ...
            struct('Target',[latTargetR; lonTargetR], ...
                   'AltScale',altExag, ...
                   'Title',sprintf('Solved transfer, %s to %s', ...
                                   info.launchStr,info.targetStr)));
    end

%% The movie, off by default because it is the expensive part of a run. The
%% four-argument (traj,veh,env,opts) form is the coorbital.viz contract; the
%% vehicle and the environment are not read by it and are passed for that
%% reason alone:
    if movieOn
                mv = coorbital.viz.globeMovie(traj,bst,env, ...
            struct('File',movieFile, ...
                   'NFrame',movieFrames, ...
                   'AltScale',altExag, ...
                   'PhaseName',{{'boost','glide','descent'}}, ...
                   'Title',sprintf('Launch %s  to  target %s', ...
                                   info.launchStr,info.targetStr)));
        fprintf('  Movie: %d frames at %g fps written to\n    %s\n', ...
                mv.nFrame,mv.frameRate,mv.file);
        fprintf('         Earth texture %s, background %s\n\n',mv.texture,mv.background);
         info.movie  = mv;
    end

%% Hand the trajectory back only when the caller asked for it. Typing
%% "run_target" at the prompt should leave the summary on screen, not bury it
%% under a dump of the whole struct as ans:
    if nargout == 0
        clear traj info;
    end
end

function ph = buildPhases(tCut,cfg)
%% Purpose:
%
%  Assemble the three-phase chain for one thrust-termination time. Factored
%  out because the range solve builds it once per bisection step and the
%  reported info struct hands the same array back to an independent checker;
%  two copies of this would be two chains that could drift apart.
%
%  THE CUTOFF IS THE tspan, NOT AN EVENT. Phase 1 keeps the burnout event, so
%  a cutoff at or beyond propellant exhaustion still stops at exhaustion; a
%  cutoff before it simply runs out of tspan, which is exactly a commanded
%  thrust termination and needs no new machinery.
%
%  THE SEPARATION LINK THROWS THE WHOLE BOOSTER AWAY -- dry structure and the
%  propellant an early cutoff left in it -- so the state mass handed to the
%  unpowered phases is exactly cfg.mGlide, whatever the cutoff was. That is
%  what lets one glide vehicle struct serve every cutoff, and it is what
%  coorbital.eom.massConstant checks on every derivative evaluation.
%
%% Inputs:
%
%  tCut             [1 x 1]                     Thrust-termination time (s)
%
%  cfg              Struct                      Immutable configuration built
%                                               by the caller: eomBst, eomGld
%                                               (function handles), schedBst,
%                                               schedGld, schedDsc (guidance
%                                               schedules), mBurnout (kg),
%                                               mGlide (kg), hHandoff (m),
%                                               hStop (m), tMaxGlide (s),
%                                               tMaxDesc (s)
%
%% Outputs:
%
%  ph               [1 x 3] struct              Phase array for
%                                               coorbital.prop.phaseRun
%
%% Revision History:
%  Michael Casey                                                08/07/2026
%  Copyright 2026 Coorbital, Inc.
%% ------------------------ Begin Code Sequence ---------------------------

%% Phase 1, boost: burns until the commanded cutoff or propellant exhaustion,
%% whichever comes first, then jettisons the ENTIRE booster:
         ph(1).eom = cfg.eomBst;
       ph(1).guide = @(t,x) coorbital.guide.pitchProgram(t,x,cfg.schedBst);
   ph(1).terminate = @(t,x) coorbital.prop.eventBurnout(t,x,cfg.mBurnout);
       ph(1).tspan = [0 tCut];
        ph(1).link = @(x) [x(1:6); cfg.mGlide];

%% Phase 2, glide: unpowered, ends on the descending crossing of the handoff
%% altitude:
         ph(2).eom = cfg.eomGld;
       ph(2).guide = @(t,x) coorbital.guide.prescribed(t,x,cfg.schedGld);
   ph(2).terminate = @(t,x) coorbital.prop.eventAltitude(t,x,cfg.hHandoff);
       ph(2).tspan = [0 cfg.tMaxGlide];
        ph(2).link = [];

%% Phase 3, descent: same airframe and same equations as the glide, ends at
%% impact. Nothing is jettisoned here, so there is no link:
         ph(3).eom = cfg.eomGld;
       ph(3).guide = @(t,x) coorbital.guide.prescribed(t,x,cfg.schedDsc);
   ph(3).terminate = @(t,x) coorbital.prop.eventAltitude(t,x,cfg.hStop);
       ph(3).tspan = [0 cfg.tMaxDesc];
        ph(3).link = [];
end

function [rngM,traj] = flyRange(tCut,cfg,rI)
%% Purpose:
%
%  Fly the whole boost-glide-descent chain for one thrust-termination time and
%  return the great-circle surface range from the launch point to the impact
%  point. This is the function coorbital.util.rangeSolve bisects on, and it is
%  the only place the range solve touches the physics.
%
%  A PROPAGATION THAT DID NOT COMPLETE IS AN ERROR, NOT A SHORT RANGE. If a
%  phase runs out of horizon -- the usual cause being a cutoff so early that
%  the vehicle never climbs to the handoff altitude, so the descending-crossing
%  event never fires -- the trajectory still ends SOMEWHERE and still has a
%  great-circle range. Returning it would feed the bisection a number from a
%  flight that never happened and quietly break the monotonicity the whole
%  method rests on. It throws instead, naming the bracket entry to change.
%
%% Inputs:
%
%  tCut             [1 x 1]                     Thrust-termination time (s)
%
%  cfg              Struct                      Immutable configuration; see
%                                               buildPhases, plus x0 [7 x 1],
%                                               bst (struct), env (struct),
%                                               rE (m), lat0 (rad), lon0 (rad)
%
%  rI               [1 x 1]                     Radius of the impact sphere
%                                               (m); every range in the script
%                                               is measured on it
%
%% Outputs:
%
%  rngM             [1 x 1]                     Surface range from the launch
%                                               point to impact (m)
%
%  traj             Struct                      The trajectory itself, from
%                                               coorbital.prop.phaseRun
%
%% Revision History:
%  Michael Casey                                                08/07/2026
%  Copyright 2026 Coorbital, Inc.
%% ------------------------ Begin Code Sequence ---------------------------

                ph = buildPhases(tCut,cfg);
              traj = coorbital.prop.phaseRun(ph,cfg.x0,cfg.bst,cfg.env);

%% All three phases must have run, and the last must have stopped ON the
%% impact altitude rather than at its horizon:
             hEndM = traj.x(end,1) - cfg.rE;
             nPhRn = numel(unique(traj.phaseIdx));
    if nPhRn < 3 || abs(hEndM - cfg.hStop) > 1e-3
        error('coorbital:runTarget:propagationIncomplete', ...
            ['A cutoff at t = %.4f s produced %d of 3 phases and ended at ' ...
             'h = %.3f km against a %.3f km stop altitude. The flight did ' ...
             'not complete, so its range is meaningless to the bisection. ' ...
             'The usual cause is a cutoff too early for the vehicle to reach ' ...
             'the %.3f km handoff at all, so the descending-crossing event ' ...
             'never fires: raise cutFracMin.'], ...
            tCut,nPhRn,hEndM./1000,cfg.hStop./1000,cfg.hHandoff./1000);
    end

              rngM = rI.*coorbital.util.greatCircle(cfg.lat0,cfg.lon0, ...
                                                    traj.x(end,3),traj.x(end,2));
end

function d = bandShortfall(fTarget,fMin,fMax)
%% Purpose:
%
%  Signed distance from a target value to the nearer edge of an achievable
%  band, for the refusal message. Positive when the target is beyond the top
%  of the band, negative when it is below the bottom, zero when it is inside.
%
%% Inputs:
%
%  fTarget          [1 x 1]                     Required value (m)
%
%  fMin             [1 x 1]                     Bottom of the band (m)
%
%  fMax             [1 x 1]                     Top of the band (m)
%
%% Outputs:
%
%  d                [1 x 1]                     Signed shortfall (m)
%
%% Revision History:
%  Michael Casey                                                08/07/2026
%  Copyright 2026 Coorbital, Inc.
%% ------------------------ Begin Code Sequence ---------------------------

                 d = 0;
    if fTarget > fMax
                 d = fTarget - fMax;
    elseif fTarget < fMin
                 d = fTarget - fMin;
    end
end

function a = wrapPi(a)
%% Purpose:
%
%  Wrap an angle difference into (-pi,pi], so that a bearing error of a
%  billionth of a degree either side of north reads as a billionth of a degree
%  and not as very nearly a full turn.
%
%% Inputs:
%
%  a                [1 x 1]                     Angle (rad)
%
%% Outputs:
%
%  a                [1 x 1]                     The same angle in (-pi,pi]
%                                               (rad)
%
%% Revision History:
%  Michael Casey                                                08/07/2026
%  Copyright 2026 Coorbital, Inc.
%% ------------------------ Begin Code Sequence ---------------------------

                 a = mod(a + pi,2.*pi) - pi;
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

function [why,ok] = whyCutoff(tEnd,mEnd,mBurnoutT,tCut)
%% Purpose:
%
%  Diagnose why the boost phase stopped. In THIS script the nominal case is
%  the SOLVED thrust termination -- the phase running out of its shortened
%  tspan -- which in every other entry script would be a failure. Propellant
%  exhaustion is the other legitimate outcome and happens only when the solve
%  wanted the full burn; anything else is named as a fault.
%
%% Inputs:
%
%  tEnd             [1 x 1]                     Time at the end of phase 1 (s)
%
%  mEnd             [1 x 1]                     Mass state there (kg)
%
%  mBurnoutT        [1 x 1]                     Total burnout mass (kg)
%
%  tCut             [1 x 1]                     Commanded cutoff time (s)
%
%% Outputs:
%
%  why              Char [1 x n]                Human-readable reason
%
%  ok               [1 x 1] logical             True for either legitimate
%                                               outcome
%
%% Revision History:
%  Michael Casey                                                08/07/2026
%  Copyright 2026 Coorbital, Inc.
%% ------------------------ Begin Code Sequence ---------------------------

                ok = false;
    if abs(mEnd - mBurnoutT) < 1e-6
               why = sprintf(['propellant exhausted at t = %.4f s, m = %.1f kg; ' ...
                              'the solve wanted the FULL burn (nominal)'],tEnd,mEnd);
                ok = true;
    elseif abs(tEnd - tCut) < 1e-6
               why = sprintf(['thrust terminated on command at t = %.4f s with ' ...
                              '%.1f kg of propellant still aboard (nominal)'], ...
                             tEnd,mEnd - mBurnoutT);
                ok = true;
    else
               why = sprintf(['stopped at t = %.4f s, neither the %.4f s commanded ' ...
                              'cutoff nor propellant exhaustion (integrator failure ' ...
                              'or an unmodelled event)'],tEnd,tCut);
    end
end

function [why,ok] = whyAltitude(tPhase,tEnd,hEnd,hTargetM,hTargetKm,tMax,what)
%% Purpose:
%
%  Diagnose why an altitude-terminated phase stopped -- the glide at its
%  handoff, the descent at impact. Nominal termination is the descending
%  crossing of the target altitude; a phase that ran out of horizon in mid-air
%  is the failure this exists to name.
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

function e = exagFor(hPeakM,rE)
%% Purpose:
%
%  The vertical exaggeration the globe figures are drawn at, scaled to the
%  apogee the run actually flew. A DISPLAY rule and nothing else: it moves no
%  number in the summary and enters no equation of motion.
%
%  THE INVARIANT, and exactly where it stops holding. The factor is chosen to
%  hold the APPARENT apogee at or under 0.3 rE -- 1913 km, a little under a
%  third of the disc -- which is high enough that the skip phugoid is plainly
%  visible and low enough that the arc still reads as a flight over a planet
%  rather than as an orbit around one. Between them the two limits mean the
%  invariant holds only over a band of apogees:
%
%    Under 63.8 km of apogee the CAP wins. The raw ratio would exceed 30, the
%    factor the rest of this library uses by hand, and the apparent apogee
%    comes out under 0.3 rE rather than at it. Nothing is lost: a low flight
%    exaggerated 30x is already legible, and letting the factor run to 300x on
%    a pad abort would paint a stationary vehicle halfway to the Moon. The
%    max(hPeakM,1) guard is what keeps a zero-apogee case out of a division by
%    zero and hands it the cap.
%
%    Over 956.7 km of apogee the FLOOR wins and the invariant is ABANDONED,
%    not merely approached. At 2x the apparent apogee is 2 hPeak, which passes
%    0.3 rE at 956.7 km and grows without bound after it: a 3000 km apogee is
%    drawn 6000 km off the surface, well outside the invariant. That is
%    deliberate. Below 2x an exaggeration is indistinguishable from true scale,
%    so the alternative is not a better picture but an honest-looking one that
%    does nothing; and a flight that apogees above 957 km needs no help being
%    seen -- BM/run_ballistic apogees at 1620 km and ships 3x by hand. The rule
%    exists for flights the sphere would otherwise swallow.
%
%  The factor is floored to an integer so the caption reads "exaggerated 16x"
%  rather than "16.2038x", which is a precision the picture does not have.
%
%% Inputs:
%
%  hPeakM           [1 x 1]                     Peak altitude above the
%                                               reference sphere, over the
%                                               whole flight (m)
%
%  rE               [1 x 1]                     Reference sphere radius (m),
%                                               from coorbital.util.missileConst
%
%% Outputs:
%
%  e                [1 x 1]                     Altitude exaggeration to hand
%                                               coorbital.viz as AltScale (-),
%                                               a whole number in [2 .. 30]
%
%% Revision History:
%  Michael Casey                                                08/07/2026
%  Copyright 2026 Coorbital, Inc.
%% ------------------------ Begin Code Sequence ---------------------------

    assert(isscalar(hPeakM) && isfinite(hPeakM), ...
        'exagFor needs a finite scalar apogee in metres; %s given.', ...
        mat2str(hPeakM));
                 e = min(30,max(2,floor(0.3.*rE./max(hPeakM,1))));
end

function txt = llStr(latDeg,lonDeg)
%% Purpose:
%
%  Format a latitude and longitude for a human-readable figure title, with
%  hemisphere letters rather than signed decimals. A caption reading
%  "40.96N 100.29E" is unambiguous where "40.96, 100.29" leaves the reader to
%  remember which way the signs run.
%
%  TWO EDGE CASES, both handled here rather than left to surprise a reader.
%
%  A COORDINATE THAT ROUNDS TO ZERO takes its hemisphere from the ROUNDED
%  value, not from the raw sign. Longitude -0.001 deg used to print "0.00W", a
%  west longitude of zero, which is a contradiction manufactured entirely by
%  the rounding. Exact zero is reported N and E: the equator and the prime
%  meridian are the lines the letters are measured from and have to fall on one
%  side in a two-letter scheme.
%
%  A NON-FINITE COORDINATE IS REFUSED. Formatted, NaN produced "NaNN NaNE" --
%  a caption that reads like a hemisphere letter stuck to a number and is
%  neither, printed over a figure of a flight that cannot have happened. A
%  title is the last place a silent NaN should be allowed to reach.
%
%% Inputs:
%
%  latDeg           scalar                      Latitude (deg), positive north,
%                                               finite
%
%  lonDeg           scalar                      Longitude (deg), positive east,
%                                               finite
%
%% Outputs:
%
%  txt              char                        e.g. '40.96N 100.29E'
%
%% Revision History:
%  Michael Casey                                                08/07/2026
%  Copyright 2026 Coorbital, Inc.
%% ------------------------ Begin Code Sequence ---------------------------

    assert(isscalar(latDeg) && isfinite(latDeg) && ...
           isscalar(lonDeg) && isfinite(lonDeg), ...
        'llStr needs finite scalar coordinates; got %s and %s.', ...
        mat2str(latDeg),mat2str(lonDeg));

%% Round FIRST, then take the hemisphere letters from the rounded values, so a
%% coordinate a thousandth of a degree west of the meridian does not print as
%% a west longitude of zero:
              latR = round(latDeg,2);
              lonR = round(lonDeg,2);
              hLat = 'N';
    if latR < 0
              hLat = 'S';
    end
              hLon = 'E';
    if lonR < 0
              hLon = 'W';
    end
               txt = sprintf('%.2f%c %.2f%c', ...
                             abs(latR),hLat,abs(lonR),hLon);
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
