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
%% Note -- the targeting problem is TWO controls and TWO residuals:
%
%  HGV/run_boost_glide flies launch site + azimuth + pitch program and lands
%  wherever the physics puts it. This script inverts that. Two unknowns and
%  two misses, and on anything but a non-rotating zero-bank run they COUPLE:
%
%    LAUNCH AZIMUTH psi.    Rotating the launch heading rotates the whole
%                           ground track about the pad, so it moves the impact
%                           point mostly CROSSWISE to the intended course.
%
%    THRUST-TERMINATION     The boost phase already ends at tspan(2) when the
%    TIME tCut.             burnout event has not fired first, so a shortened
%                           tspan IS a thrust cutoff and no new machinery is
%                           needed. Less burn gives less energy gives less
%                           range, so it moves the impact point mostly ALONG
%                           the course.
%
%  THE TWO RESIDUALS ARE THE COMPONENTS OF THE MISS VECTOR, target to impact,
%  resolved at the TARGET along the intended course: down-range, positive
%  long, and cross-range, positive right of course. Both zero is exactly "the
%  impact point IS the target", and it refers to no great circle at all --
%  which matters, because once the Earth turns or a bank is commanded the
%  flown track is not one. Writing the pair instead as a great-circle range
%  residual beside a cross-track offset would be the same information in a
%  worse basis: it would keep quoting an arc the vehicle no longer flies, and
%  its Jacobian would be further off diagonal.
%
%  SOLVED IN TWO STAGES, and the first is what makes the second safe:
%
%    STAGE 1, A SEED.       coorbital.util.rangeSolve bisects the cutoff, at
%                           the closed-form great-circle bearing, until the
%                           flown range matches the required one. Bisection
%                           cannot leave its bracket, its two endpoint
%                           propagations ARE the reachable envelope this
%                           script refuses against, and it puts the second
%                           stage inside the basin a Newton step needs.
%
%    STAGE 2, THE SOLVE.    coorbital.util.aimSolve runs a damped Newton on
%                           BOTH controls at once, over a finite-difference
%                           Jacobian, until both components of the miss are
%                           inside tolerance. On a non-rotating zero-bank run
%                           the seed is already the answer, so this stage
%                           converges at its first evaluation and hands back
%                           the closed-form bearing bit for bit.
%
%  WHY CUTOFF TIME AND NOT SOMETHING ELSE for the range control. The pitch
%  program's loft angle has a lofted and a depressed branch either side of a
%  maximum-range hump, so a root finder lands on whichever branch it started
%  nearest -- a bad default for an example script. The glide handoff altitude
%  bifurcates on the glide phugoid's troughs; 30 km costs 1882 km of range in
%  run_boost_glide's own configuration while every phase still reports
%  nominal. Thrust termination is also what real systems do for energy
%  management.
%
%  EARLY CUTOFF LEAVES UNBURNED PROPELLANT, and it is thrown away. At
%  separation the ENTIRE booster goes -- dry structure and the propellant
%  still in it -- so the post-separation vehicle is exactly the payload and
%  the link is phases(1).link = @(x) [x(1:6); veh.mass]. The summary prints
%  how much propellant went overboard unburned, because at a deep cutoff that
%  is most of the load.
%
%% Note -- WHAT THIS SOLUTION IS AND IS NOT, printed in the summary rather
%% than buried:
%
%  THE AZIMUTH IS SOLVED, NOT CLOSED FORM, and that is what lets earthSpin be
%  true. Note what rotation actually does here: the state longitude and
%  latitude are EARTH-FIXED and the speed is planet-relative, so the target
%  does NOT slide east under the vehicle -- it sits still, and the VEHICLE is
%  deflected off the arc it left on by the Coriolis and centrifugal terms, to
%  the right in the northern hemisphere. Measured on the shipped geometry, the
%  closed-form bearing alone left the impact point 231.552 km from the target.
%  The second stage removes that, so the case is flown instead of refused.
%
%  THE SOLVE CONTROLS BOTH AXES, so a banked trajectory is targetable.
%  Cross-range is driven to zero rather than measured and warned about. At
%  run_boost_glide's 75 deg terminal bank this geometry used to miss by
%  21.5 km on a 0.6 km range residual; it now arrives, and the cross-range
%  line in the summary is the evidence that it did rather than a caveat.
%
%  WHAT REMAINS TRUE, and all three are stated in the summary. The Newton is
%  LOCAL: seeded, damped and capped, and when it cannot converge the run is
%  REFUSED with the best point it reached rather than returned as a near miss.
%  The reachable band the envelope refusal quotes is flown at the SEED
%  azimuth, so it is indicative for the solved one and not exact. And the
%  tolerance is on the MISS ITSELF: each component is held to
%  tolRangeKm/sqrt(2), so the two together cannot exceed the distance the user
%  asked for.
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
%                                               EMPTY on EITHER of the two
%                                               refusals -- a target outside the
%                                               reachable envelope, or a
%                                               two-axis solve that did not
%                                               converge -- see info.refused and
%                                               info.refusedWhy
%
%  info             Struct                      The summary's numbers at full
%                                               precision, so a test does not
%                                               have to read them back out of
%                                               printed text; all SI except the
%                                               *Km distances.
%
%                                               ALWAYS carries refused
%                                               (logical), rngReqM, psiSeed and
%                                               psiLaunch. refusedWhy is present
%                                               ONLY when refused is true, so
%                                               read it behind that flag and not
%                                               on its own -- on a converged run
%                                               the field does not exist and
%                                               touching it raises.
%
%                                               A REFUSAL carries refusedWhy,
%                                               either 'envelope' or 'aimSolve',
%                                               and in both cases rngMinM,
%                                               rngMaxM and the rangeSolve
%                                               record in solveInfo, so it hands
%                                               back the envelope it was
%                                               measured against. An aimSolve
%                                               refusal adds the aimSolve record
%                                               in aimInfo as well, whose xHist
%                                               and message say where the Newton
%                                               got to and why it stopped.
%
%                                               A RUN THAT WAS NOT REFUSED
%                                               additionally
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
%       chapter on Great-Circle Sailing. The initial course used to SEED the
%       launch azimuth; see coorbital.util.greatCircleBearing.
%   [2] Press, W.H., et al., "Numerical Recipes," 3rd ed., Cambridge, 2007,
%       Section 9.1. Bisection, used by coorbital.util.rangeSolve.
%   [3] Dennis, J.E. and Schnabel, R.B., "Numerical Methods for Unconstrained
%       Optimization and Nonlinear Equations," SIAM, 1996, Chapters 5 and 6.
%       The finite-difference Jacobian and line-search damped Newton method
%       used by coorbital.util.aimSolve for the two-axis solve.
%   [4] Vinh, N.X., Busemann, A., Culp, R.D., "Hypersonic and Planetary Entry
%       Flight Mechanics," Univ. Michigan Press, 1980, Ch. 5.
%
%% Revision History:
%  Michael Casey                                                08/07/2026
%  Michael Casey  Two-axis targeting: the azimuth is SOLVED beside
%                 the cutoff by coorbital.util.aimSolve, which
%                 retires the rotating-Earth refusal and the
%                 downrange-only cross-track warning              08/08/2026
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
%% the launch azimuth is SOLVED from these two points, seeded on the
%% great-circle bearing between them, and any azimuth given here could only
%% contradict it:
         latTarget = 35;               %deg, geocentric latitude [-89 .. 89]
         lonTarget = -120;             %deg, longitude [-180 .. 180]
                                       %     The shipped pair is a mid-Pacific pad to the
                                       %     California coast: 3811 km on a 56.6 deg azimuth.
                                       %     Deliberately NOT equatorial and NOT due east --
                                       %     that geometry makes latitude and longitude
                                       %     interchangeable in the range formula and hides a
                                       %     transposition at every call site that uses them

%% Seed bracket -- bisection on thrust-termination time, expressed as a
%% FRACTION of the booster's full burn so the bracket survives a change of
%% booster. The endpoints are propagated before the seed solve starts and
%% their ranges ARE the reachable envelope, which the summary prints:
        cutFracMin = 0.50;             %-, shortest burn tried [0.3 .. 0.95]. Too small and the
                                       %   vehicle never climbs to hHandoff, the glide event
                                       %   never fires and the propagation is rejected; the
                                       %   shipped 0.50 reaches 203 km of range
        cutFracMax = 1.00;             %-, longest burn tried (cutFracMin .. 1.0]. 1.0 is the
                                       %   full burn and is also the maximum range this
                                       %   configuration can fly, 7738 km
        tolRangeKm = 1.0;              %km, convergence tolerance on the MISS DISTANCE, impact
                                       %    point to target [0.05 .. 50]. The two-axis solve is
                                       %    held to this over sqrt(2) on EACH component, so the
                                       %    down-range and cross-range misses together cannot
                                       %    exceed the figure asked for here. The seed bisection
                                       %    uses the bare figure on range. Each halving costs
                                       %    about one more trajectory propagation, and a
                                       %    propagation is about 0.2 s, so tightening it is cheap

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
         bankBoost = 0;                %deg, bank during boost [-90 .. 90]. A non-zero value is
                                       %     no longer a targeting fault -- the solve absorbs
                                       %     the cross-range it produces -- but it is still a
                                       %     large control authority applied where the vehicle
                                       %     is fastest, so it moves the whole envelope
         tMaxBoost = 200;              %s,  boost horizon; must exceed the burn time

%% Phase 2, glide -- prescribed control against time SINCE THE START OF THE
%% GLIDE, not since liftoff:
         glideTime = [0 6000];         %s,   nodes, strictly increasing, phase-local
        glideAlpha = [0 0];            %deg, angle of attack. NO EFFECT with the default constLD
                                       %     aero model, which ignores it; becomes live only when
                                       %     aeroFn below is swapped for an alpha-dependent model
         glideBank = [0 0];            %deg, bank [-90 .. 90]. 0 keeps the track on the
                                       %     launch-to-target great circle ONLY while the Earth
                                       %     is also still -- with earthSpin true the Coriolis
                                       %     deflection walks a zero-bank track off that arc by
                                       %     231.552 km on the shipped geometry. Either way the
                                       %     SOLVE NOW SEES IT: a bank, a rotation, or both are
                                       %     one cross-range residual to it, and the launch
                                       %     azimuth is aimed off to cancel whatever produced it
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
          descBank = [0 0];            %deg, bank [-90 .. 90]. SHIPPED AT ZERO because the
                                       %     simplest case is the one an example script should
                                       %     ship, NOT because the solve cannot handle a bank.
                                       %     run_boost_glide dives at 75 deg for a steep, fast
                                       %     arrival, and a banked segment does turn the heading
                                       %     -- 166.4 deg of it between handoff and impact on
                                       %     the SHIPPED 20N 155W to 35N 120W geometry --
                                       %     walking the impact point off the launch-to-target
                                       %     great circle. That used to be 21.52 km of miss
                                       %     against a 0.59 km range residual; the two-axis
                                       %     solve now aims off to cancel it and arrives. Set it
                                       %     to 75 when you want the terminal dive, and read the
                                       %     cross-range line to see what it cost the azimuth
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
         earthSpin = false;            %true enables the Coriolis and centrifugal terms. It is
                                       %SUPPORTED: the state is Earth-fixed, so the target does
                                       %not move, but the vehicle is deflected off the arc it
                                       %left on -- 231.552 km of miss on the shipped geometry
                                       %with the closed-form bearing alone -- and the two-axis
                                       %solve aims off to remove it. Shipped false only so the
                                       %default run is the simplest one; turning it on costs a
                                       %few more propagations and changes no other setting

%% Output:
         showPlots = true;             %false to skip the figures, e.g. under matlab -batch
           altExag = 1;                %-, VERTICAL EXAGGERATION of the globe still and of the
                                       %   movie. 1 is TRUE SCALE and is the shipped default:
                                       %   the movie's altitude inset is true-scale already and
                                       %   carries the profile, so the globe need not lie about
                                       %   it, and coorbital.viz drops the "(altitude
                                       %   exaggerated Nx)" caption clause at unity. Any
                                       %   positive number is used as given. The char 'auto'
                                       %   selects the ADAPTIVE rule instead -- exagFor below,
                                       %   which holds the apparent apogee at or under 0.3 rE,
                                       %   capped at 30 and floored at 2, and returns 16 on the
                                       %   shipped case. Ask for it when the skip phugoid is
                                       %   the point of the picture: at true scale a 118 km
                                       %   glide is one part in 54 of an Earth radius
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
           altExag = overrideOf(opts,'altExag',altExag);
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

%% The display exaggeration, checked here rather than at the figures: a
%% mistyped one would otherwise cost a whole range solve before it was noticed.
%% Either a positive number, used as given, or the char 'auto', which defers to
%% the adaptive rule once the flown apogee is known:
    if ischar(altExag) || isstring(altExag)
        assert(strcmpi(char(altExag),'auto'), ...
            ['altExag = "%s" is not understood; it must be a positive number ' ...
             'or the word ''auto'', which selects the adaptive rule.'], ...
            char(altExag));
    else
        assert(isnumeric(altExag) && isscalar(altExag) && ...
               isfinite(altExag) && altExag > 0, ...
            ['altExag must be a positive finite scalar or the word ''auto''; ' ...
             'got %s. coorbital.viz refuses a zero or negative AltScale.'], ...
            mat2str(altExag));
    end

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
%% The seed azimuth: the closed-form bearing, which is a GUESS
%% -----------------------------------------------------------------
%% The initial bearing of the launch-to-target great circle IS the library's
%% heading state psi -- both are measured clockwise from north -- so it goes
%% straight into the launch state below with no conversion. It is the whole
%% answer only on a non-rotating zero-bank run; everywhere else it is the
%% starting point the two-axis solve improves on, and it is the right starting
%% point precisely because it is exact in the easy case and close in the hard
%% ones:
            psiSeed = coorbital.util.greatCircleBearing(latLaunchR,lonLaunchR, ...
                                                        latTargetR,lonTargetR);

%% Liftoff state TEMPLATE. Everything but the mass is fixed; the pad
%% flight-path angle is the first commanded pitch attitude, so the angle of
%% attack starts at zero. Entry 6 is the heading, and it is the ONE entry the
%% solve rewrites -- flyChain overwrites it per trial rather than rebuilding
%% the state, so a trial flight can differ from the seed in the azimuth and in
%% nothing else:
                x0 = [c.rE + hLaunchM; ...
                      lonLaunchR; ...
                      latLaunchR; ...
                      vLaunch; ...
                      pitchAngR(1); ...
                      psiSeed; ...
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
%% The aim frame: the basis both residuals are resolved in
%% -----------------------------------------------------------------
%% One direction, fixed for the whole solve: the course of the launch-to-target
%% great circle AT THE TARGET, which is the arriving bearing and therefore the
%% direction "long" points in. greatCircleBearing is not symmetric, so the
%% arriving course is the DEPARTING course of the reversed arc turned through
%% pi; taking the launch-end bearing instead would tilt the basis by the
%% convergence of the meridians, several degrees at this range.
%%
%% Fixed rather than recomputed per trial on purpose. A basis that moved with
%% the impact point would make the Jacobian describe the basis as much as the
%% trajectory, and the finite-difference columns would stop meaning what the
%% Newton step assumes they mean:
          psiCourse = mod(coorbital.util.greatCircleBearing(latTargetR,lonTargetR, ...
                                                            latLaunchR,lonLaunchR) + pi, ...
                          2.*pi);
           aim.latT = latTargetR;
           aim.lonT = lonTargetR;
           aim.psiC = psiCourse;
           aim.rI   = rI;

%% -----------------------------------------------------------------
%% Stage 1: bracket the envelope and bisect the cutoff at the seed azimuth
%% -----------------------------------------------------------------
%% The bracket endpoints in seconds. rangeSolve evaluates both before it
%% iterates and returns their achieved ranges in info.fMin and info.fMax,
%% which IS the reachable envelope -- so bracketing and reporting the envelope
%% are the same two propagations, not two pairs. Everything here flies the SEED
%% azimuth, which is what makes the envelope indicative rather than exact for
%% the solved one; it is quoted that way in the summary:
              tLo  = cutFracMin.*tBurn;
              tHi  = cutFracMax.*tBurn;
            fRange = @(tCut) flyRange(psiSeed,tCut,cfg,rI);
[tSeed,rngSeed,sv] = coorbital.util.rangeSolve(rngReq,fRange,tLo,tHi,tolRngM);

%% Fractions of full burn, for the report. The cutoff as a bare number of
%% seconds says nothing without the burn it is a fraction of:
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
    info.refusedWhy  = 'envelope';
       info.rngReqM  = rngReq;
      info.psiSeed   = psiSeed;
     info.psiLaunch  = psiSeed;
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
%% Stage 2: drive BOTH components of the miss to zero
%% -----------------------------------------------------------------
%% Two controls, two residuals, one damped Newton. The safeguards are the
%% caller's job because only the caller knows the scales:
%%
%%    dx0        1e-4 rad of heading swings the impact point about 380 m
%%               crosswise at this range, and 1e-4 of the burn time is about
%%               8 ms of cutoff, worth roughly 1.5 km of range. Both are far
%%               above the few metres an adaptive integrator's own step
%%               selection moves the impact point by, which is the floor a
%%               difference step has to clear, and both are small enough to be
%%               a chord and not a secant across real curvature. The cutoff
%%               step is written against tBurn rather than as a literal so a
%%               different booster carries it.
%%
%%    maxStep    10 deg of heading and 5 % of the burn per iteration. aimSolve
%%               defaults to unbounded because a library function cannot know
%%               these; this one can. The seed is already range-converged, so
%%               a step anywhere near either cap means the Newton direction is
%%               wrong, and capping it costs an iteration where an uncapped
%%               first step could swing the aim to the far side of the planet
%%               and refuse from somewhere nothing flies.
%%
%%    tolF       the user's tolerance over sqrt(2), so that BOTH components
%%               inside it puts the total miss inside the figure asked for.
%%               aimSolve tests max(abs(f)), which is the infinity norm, and
%%               the worst case for the hypotenuse is the two components equal.
%%
%% A trial pair that will not fly -- a cutoff too early to reach the handoff,
%% say -- throws out of flyChain, and aimSolve treats that as an infeasible
%% point and halves the step rather than giving up:
           tolAimM = tolRngM./sqrt(2);
            residF = @(xAim) flyMiss(xAim,cfg,aim);
             xSeed = [psiSeed; tSeed];
             dxAim = [1.0e-4; 1.0e-4.*tBurn];
           optsAim = struct('maxStep',[deg2rad(10); 0.05.*tBurn]);
    [xAim,fAch,av] = coorbital.util.aimSolve(residF,xSeed,dxAim,tolAimM,optsAim);
         psiLaunch = xAim(1);
              tCut = xAim(2);
          fracSol  = tCut./tBurn;

%% The residual at the SEED, which aimSolve records as the first column of its
%% history and which therefore costs nothing to read. It is the miss the
%% closed-form bearing and the bisected cutoff produce on their own -- exactly
%% what this script used to report as its answer -- so the summary can print
%% the before and the after side by side:
          seedDwnM = av.fHist(1,1);
          seedXtkM = av.fHist(2,1);
         seedMissM = hypot(seedDwnM,seedXtkM);

%% -----------------------------------------------------------------
%% Refuse, loudly, when the two-axis solve did not converge
%% -----------------------------------------------------------------
%% aimSolve returns converged = false with the best point it reached rather
%% than throwing, for the same reason rangeSolve does: only this script has the
%% vocabulary. A NEAR MISS RETURNED AS AN ANSWER IS THE FAILURE MODE THIS
%% WHOLE FILE EXISTS TO PREVENT, so it is refused on exactly the terms the
%% envelope refusal uses -- nothing thrown, empty traj, info.refused true:
    if ~av.converged
        fprintf('\n');
        fprintf('===== Point-to-point targeting: REFUSED =================\n');
        fprintf('  The two-axis aim solve did not converge on the target.\n');
        fprintf('\n');
        fprintf('    launch           %10.4f %-5s %10.4f deg  (lat, lon)\n', ...
                latLaunch,'deg',lonLaunch);
        fprintf('    target           %10.4f %-5s %10.4f deg  (lat, lon)\n', ...
                latTarget,'deg',lonTarget);
        fprintf('    required range   %10.2f %-5s (great circle on the r = %.3f km sphere)\n', ...
                rngReq./1000,'km',rI./1000);
        fprintf('    reachable        %10.2f to %.2f km  (flown at the SEED azimuth)\n', ...
                sv.fMin./1000,sv.fMax./1000);
        fprintf('    seed azimuth     %10.6f %-5s (the closed-form great-circle bearing)\n', ...
                rad2deg(psiSeed),'deg');
        fprintf('    best azimuth     %10.6f %-5s (%+.6f deg from the seed)\n', ...
                rad2deg(psiLaunch),'deg',rad2deg(wrapPi(psiLaunch - psiSeed)));
        fprintf('    best cutoff      %10.4f %-5s (%.6f of the %.4f s full burn)\n', ...
                tCut,'s',fracSol,tBurn);
        fprintf('    best miss        %10.2f %-5s (down-range %+.2f, cross-range %+.2f)\n', ...
                hypot(fAch(1),fAch(2)),'m',fAch(1),fAch(2));
        fprintf('    tolerance        %10.2f %-5s (per component; %.2f m on the total miss)\n', ...
                tolAimM,'m',tolRngM);
        fprintf('    propagations     %10d %-5s (%d in the seed bisection, %d in the Newton)\n', ...
                sv.nEval + av.nEval,'',sv.nEval,av.nEval);
        fprintf('    solver said      %s\n',av.identifier);
        fprintf('\n');
        fprintf('  %s\n',av.message);
        fprintf('\n');
        fprintf('  THE BEST POINT REACHED IS NOT AN ANSWER, and it is not being returned as\n');
        fprintf('  one. A targeting script that handed back its nearest miss would hand back\n');
        fprintf('  something that LOOKS like a solution: a cutoff, an azimuth, a trajectory\n');
        fprintf('  and a summary, all of them describing a flight that does not arrive. The\n');
        fprintf('  numbers above are in info.aimInfo for anyone who wants to restart from\n');
        fprintf('  them.\n');
        fprintf('=========================================================\n\n');
              traj = [];
       info.refused  = true;
    info.refusedWhy  = 'aimSolve';
       info.rngReqM  = rngReq;
      info.psiSeed   = psiSeed;
     info.psiLaunch  = psiLaunch;
          info.tCut  = tCut;
         info.missM  = hypot(fAch(1),fAch(2));
         info.downM  = fAch(1);
       info.xTrackM  = fAch(2);
       info.rngMinM  = sv.fMin;
       info.rngMaxM  = sv.fMax;
         info.tBurn  = tBurn;
        info.tCutLo  = tLo;
        info.tCutHi  = tHi;
     info.solveInfo  = sv;
       info.aimInfo  = av;
        if nargout == 0
            clear traj info;
        end
        return;
    end

%% -----------------------------------------------------------------
%% Fly the solved trajectory
%% -----------------------------------------------------------------
%% One more propagation at the converged control pair. aimSolve already flew
%% it but returns only the two-component residual, so the state history is
%% re-created here rather than carried out through a persistent variable. The
%% assert is the repeatability check that licenses everything below: if the
%% same controls give a different impact point, no number in the summary means
%% anything:
[traj,latIm,lonIm] = flyChain(psiLaunch,tCut,cfg);
            misVec = missVector(latIm,lonIm,aim);
    assert(max(abs(misVec - fAch)) < 1e-6, ...
        ['re-flying the solved controls gave a miss of [%.6f, %.6f] m against ' ...
         'the solver''s [%.6f, %.6f] m; the propagation is not repeatable and ' ...
         'nothing below can be trusted'],misVec(1),misVec(2),fAch(1),fAch(2));

%% Phase boundaries. phaseRun records the boundary sample ONCE, carrying the
%% outgoing phase's control, so the last row of phase k is that phase's
%% terminal state on the NEAR side of any staging jump:
               kBO = find(traj.phaseIdx == 1,1,'last');
               kHO = find(traj.phaseIdx == 2,1,'last');
                nS = numel(traj.t);

%% -----------------------------------------------------------------
%% The miss, measured from the FLOWN state
%% -----------------------------------------------------------------
%% missM is the great-circle distance impact-to-target, measured directly, and
%% the two components are the same vector resolved in the aim frame. They agree
%% BY CONSTRUCTION -- the components are the arc times a unit direction -- so
%% the assert below is checking the arithmetic and not the physics, which is
%% why it is an equality and not a tolerance a reader is invited to judge:
             missM = rI.*coorbital.util.greatCircle(latIm,lonIm, ...
                                                    latTargetR,lonTargetR);
             downM = misVec(1);
           xTrackM = misVec(2);
          missHypM = hypot(downM,xTrackM);
    assert(abs(missHypM - missM) < 1e-6, ...
        ['the miss components [%.6f, %.6f] m compose to %.6f m against the ' ...
         '%.6f m measured directly; the aim frame is not a unit basis'], ...
        downM,xTrackM,missHypM,missM);

%% The achieved great-circle range, launch to impact, and its residual against
%% the required one. NOT what the solve drives any more -- that is the miss
%% above -- but it is the number a reader has always used to judge whether the
%% cutoff came out sensible, and on a rotating or banked run the gap between it
%% and the down-range component is itself informative: it is the part of the
%% miss the old downrange-only solve could not see:
            rngAch = rI.*coorbital.util.greatCircle(latLaunchR,lonLaunchR, ...
                                                    latIm,lonIm);
              resM = rngAch - rngReq;

%% How far the flown track ended up from the arc it left on. With rotation off
%% and every bank zero this is zero to rounding; with either on it is the
%% deflection the azimuth solve had to aim off against, so it is a measurement
%% of the coupling rather than a caveat about it:
            angFly = coorbital.util.greatCircle(latLaunchR,lonLaunchR,latIm,lonIm);
            psiFly = coorbital.util.greatCircleBearing(latLaunchR,lonLaunchR, ...
                                                       latIm,lonIm);
            dPsiR  = wrapPi(psiFly - psiLaunch);

%% The heading the vehicle actually left the pad on must be the one the solve
%% asked for, to machine precision -- flyChain writes entry 6 of the state from
%% it and ode45 never touches the initial condition, so any disagreement here
%% means the launch state was assembled wrong:
             psiX0 = traj.x(1,6);
    assert(abs(wrapPi(psiX0 - psiLaunch)) < 1e-12, ...
        ['the flown initial heading %.15f rad is not the solved launch ' ...
         'azimuth %.15f rad; the launch state was assembled wrong'], ...
        psiX0,psiLaunch);

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

%% THERE IS DELIBERATELY NO EASTWARD-DRIFT FIGURE HERE ANY MORE, and the
%% deletion is the point. This script used to compute
%% c.omegaE*rI*cos(latTarget)*tFlight -- how far the ground at the target's
%% latitude sweeps east over the flight, 628 km on the rotating shipped case --
%% and print it as "the scale of the rotation effect". It is not. It is the
%% scale of a mechanism this file no longer claims: the state is Earth-fixed,
%% so the target does not travel that 628 km relative to anything the solve can
%% see. The quantity that IS the scale of the rotation effect is the deflection
%% of the VEHICLE, and it is measured rather than derived -- seedMissM, 231.55
%% km on that same run, a factor of 2.7 smaller. A retracted mechanism's number
%% left standing under a corrected sentence is the harder error to notice of
%% the two, because the prose reads right.

%% Whether ANY phase commands a bank, and whether the Earth turned. Together
%% they decide which paragraph the limitations block prints: a run with
%% neither is the one case where the seed azimuth was already the answer, and
%% saying so is a different statement from saying the solve removed a
%% deflection. One paragraph with numbers substituted into it cannot do both --
%% at 75 deg of descent bank the sentence "cross-range came out at 21515.22 m
%% because the bank angle is zero throughout" is self-refuting, and it used to
%% print directly under a WARNING saying the opposite. Gated on the RADIAN
%% schedules, so a user block written in any units that convert to zero is
%% still recognised as zero:
           anyBank = any([bankBoostR glideBnkR descBnkR] ~= 0);
           spinsUp = env.omegaE ~= 0;
           spinTxt = 'OFF';
    if spinsUp
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
    fprintf('    range residual   %+10.2f %-5s  (achieved minus required, signed. NOT what the\n',resM,'m');
    fprintf('                                        solve drives; see the two components below)\n');
    fprintf('    MISS DISTANCE    %10.2f %-5s  (great circle, impact point to target,\n', ...
            missM,'m');
    fprintf('                                        computed from the FLOWN state and not\n');
    fprintf('                                        from the solver''s own bookkeeping)\n');
    fprintf('    down-range miss  %+10.2f %-5s  (component along the intended course at the\n', ...
            downM,'m');
    fprintf('                                        target, positive LONG)\n');
    fprintf('    cross-range miss %+10.2f %-5s  (component across it, positive RIGHT of\n', ...
            xTrackM,'m');
    fprintf('                                        course. BOTH are what the solve drives to\n');
    fprintf('                                        zero, and they compose to the %.2f m miss)\n', ...
            missHypM);
    fprintf('    tolerance        %10.2f %-5s  (on EACH component; %.2f m on the total miss,\n', ...
            tolAimM,'m',tolRngM);
    fprintf('                                        which is what tolRangeKm asks for)\n');
    fprintf('\n');
    fprintf('    seed azimuth     %10.6f %-5s  (clockwise from north; the closed-form initial\n', ...
            rad2deg(psiSeed),'deg');
    fprintf('                                        bearing of the launch-to-target great\n');
    fprintf('                                        circle, which is the GUESS)\n');
    fprintf('    launch azimuth   %10.6f %-5s  (SOLVED, %+.6f deg from the seed)\n', ...
            rad2deg(psiLaunch),'deg',rad2deg(wrapPi(psiLaunch - psiSeed)));
    fprintf('    seed miss        %10.2f %-5s  (down-range %+.2f, cross-range %+.2f: where the\n', ...
            seedMissM,'m',seedDwnM,seedXtkM);
    fprintf('                                        seed azimuth and the bisected cutoff put\n');
    fprintf('                                        the vehicle on their own, before stage 2)\n');
    fprintf('    flown azimuth    %10.6f %-5s  (bearing launch to impact, from the flown\n', ...
            rad2deg(psiFly),'deg');
    fprintf('                                        state; %+.3e deg from the commanded one.\n', ...
            rad2deg(dPsiR));
    fprintf('                                        Non-zero means the track was deflected off\n');
    fprintf('                                        the arc it left on -- by rotation, by bank,\n');
    fprintf('                                        or by both)\n');
    fprintf('\n');
    fprintf('    solved cutoff    %10.4f %-5s  (%.6f of the %.4f s full burn)\n', ...
            tCut,'s',fracSol,tBurn);
    fprintf('    seed cutoff      %10.4f %-5s  (%.3f km of range at the seed azimuth; the\n', ...
            tSeed,'s',rngSeed./1000);
    fprintf('                                        bisection''s answer, which stage 2 moved by\n');
    fprintf('                                        %+.4f s)\n',tCut - tSeed);
    fprintf('    propellant burned%10.1f %-5s  (of %.1f kg loaded)\n', ...
            mpBurned,'kg',bst.massProp);
    fprintf('    thrown away      %10.1f %-5s  (unburned, jettisoned inside the booster at\n', ...
            mpWasted,'kg');
    fprintf('                                        separation; %.1f %% of the load)\n', ...
            100.*mpWasted./bst.massProp);
    fprintf('    reachable        %10.3f to %.3f km  (the bracket endpoints, flown at the\n', ...
            sv.fMin./1000,sv.fMax./1000);
    fprintf('                                        SEED azimuth, so indicative for the solved\n');
    fprintf('                                        one rather than exact)\n');
    fprintf('    cutoff bracket   %10.4f to %.4f s   (%.3f to %.3f of full burn)\n', ...
            tLo,tHi,fracLo,fracHi);
    fprintf('    iterations       %10d %-5s  (%d bisection steps seeding %d Newton step(s))\n', ...
            sv.iterations + av.iterations,'',sv.iterations,av.iterations);
    fprintf('    propagations     %10d %-5s  (%d in the seed bisection, %d in the Newton, plus\n', ...
            sv.nEval + av.nEval + 1,'',sv.nEval,av.nEval);
    fprintf('                                        one to re-fly the answer)\n');
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
    fprintf('  WHAT THIS TARGETING SOLUTION IS, AND WHAT IT IS NOT\n');
    fprintf('    1. THE AZIMUTH IS SOLVED, NOT CLOSED FORM. The closed-form great-circle\n');
    fprintf('       bearing is only the SEED. Earth rotation is %s on this run, env.omegaE =\n',spinTxt);
    fprintf('       %.6e rad/s, and the state is EARTH-FIXED with a planet-relative\n',env.omegaE);
    fprintf('       speed, so the target never slides east under the vehicle -- it sits still,\n');
    fprintf('       and rotation instead DEFLECTS THE VEHICLE off the arc it left on.\n');
    if spinsUp
    fprintf('       IT IS MODELLED ON THIS RUN, and the size of it is MEASURED and not\n');
    fprintf('       derived: the seed azimuth alone left the vehicle %.2f km from the target,\n', ...
            seedMissM./1000);
    fprintf('       %.2f km of that crosswise. Aiming off by %+.6f deg brought it to %.2f m.\n', ...
            abs(seedXtkM)./1000,rad2deg(wrapPi(psiLaunch - psiSeed)),missM);
    else
    fprintf('       IT IS NOT MODELLED ON THIS RUN, because earthSpin is false. Set it true and\n');
    fprintf('       the case is flown, not refused: on the shipped geometry the deflection is\n');
    fprintf('       231.552 km of miss, almost all of it crosswise, and the solve removes it.\n');
    end
    fprintf('    2. THE SOLVE CONTROLS BOTH AXES. Down-range and cross-range are driven to zero\n');
    fprintf('       together by a damped Newton on the launch azimuth and the cutoff time, so\n');
    fprintf('       cross-range is a MEASUREMENT of the solve here and not a caveat about it.\n');
    if anyBank
    fprintf('       THIS RUN COMMANDS A NON-ZERO BANK -- boost %.1f deg, glide %.1f deg,\n', ...
            rad2deg(bankBoostR),rad2deg(glideBnkR(end)));
    fprintf('       descent %.1f deg -- so the ground track is NOT the launch-to-target great\n', ...
            rad2deg(descBnkR(1)));
    fprintf('       circle, and a downrange-only solve could not see the difference. Measured\n');
    fprintf('       on this run: the seed put the vehicle %.2f m off, %.2f m of it crosswise,\n', ...
            seedMissM,abs(seedXtkM));
    fprintf('       and the second control closed it to %.2f m. That is the case that used to\n', ...
            missM);
    fprintf('       converge on range and miss by 21.5 km.\n');
    elseif spinsUp
    fprintf('       EVERY BANK IS ZERO ON THIS RUN -- boost %.1f deg, glide %.1f deg, descent\n', ...
            rad2deg(bankBoostR),rad2deg(glideBnkR(end)));
    fprintf('       %.1f deg -- so the %.2f m of cross-range the seed left came ENTIRELY FROM\n', ...
            rad2deg(descBnkR(1)),abs(seedXtkM));
    fprintf('       THE ROTATION and not from any commanded turn. That is the whole of the\n');
    fprintf('       231.552 km this script used to refuse rather than fly, and it is now\n');
    fprintf('       %.2f m. Bank a phase as well and the two effects simply add into the same\n', ...
            abs(xTrackM));
    fprintf('       residual; the solve does not need to tell them apart.\n');
    else
    fprintf('       EVERY BANK IS ZERO ON THIS RUN -- boost %.1f deg, glide %.1f deg, descent\n', ...
            rad2deg(bankBoostR),rad2deg(glideBnkR(end)));
    fprintf('       %.1f deg -- and rotation is off, so the seed left only %.2f m of\n', ...
            rad2deg(descBnkR(1)),abs(seedXtkM));
    fprintf('       cross-range: a zero-bank track over a non-rotating sphere never leaves the\n');
    fprintf('       great circle it departed on, so the seed azimuth WAS the answer here and\n');
    fprintf('       the Newton returned it unchanged.\n');
    fprintf('       THAT IS A PROPERTY OF THIS CONFIGURATION, NOT A GENERAL GUARANTEE, and it\n');
    fprintf('       is measured above rather than assumed. Bank a phase, or turn the Earth on,\n');
    fprintf('       and the track curves off the arc: at run_boost_glide''s 75 deg terminal\n');
    fprintf('       bank this geometry puts the seed 21.5 km wide, and the azimuth is aimed\n');
    fprintf('       off to cancel it.\n');
    end
    fprintf('    3. THE NEWTON IS LOCAL, AND A FAILURE IS REFUSED RATHER THAN RETURNED. It is\n');
    fprintf('       seeded on the bisection, damped, and capped at 10 deg of heading and 5 %% of\n');
    fprintf('       the burn per iteration. When it cannot reach the tolerance it prints a\n');
    fprintf('       REFUSED banner with the best point it found; it never hands back a near\n');
    fprintf('       miss dressed as a solution. Two consequences worth reading: the reachable\n');
    fprintf('       band above was flown at the SEED azimuth, so it is indicative for the\n');
    fprintf('       solved one and not exact, and the tolerance is %.2f m PER COMPONENT so that\n', ...
            tolAimM);
    fprintf('       the two together stay inside the %.2f m tolRangeKm asks for.\n',tolRngM);
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
         info.downM  = downM;
       info.xTrackM  = xTrackM;
      info.missHypM  = missHypM;
        info.residM  = resM;
      info.psiSeed   = psiSeed;
     info.psiLaunch  = psiLaunch;
        info.dPsiAim = wrapPi(psiLaunch - psiSeed);
        info.psiFly  = psiFly;
         info.dPsiR  = dPsiR;
     info.seedMissM  = seedMissM;
      info.seedDownM = seedDwnM;
    info.seedXTrackM = seedXtkM;
         info.tSeed  = tSeed;
      info.rngSeedM  = rngSeed;
          info.tCut  = tCut;
       info.cutFrac  = fracSol;
         info.tBurn  = tBurn;
        info.tCutLo  = tLo;
        info.tCutHi  = tHi;
       info.rngMinM  = sv.fMin;
       info.rngMaxM  = sv.fMax;
       info.tolAimM  = tolAimM;
      info.seedIter  = sv.iterations;
      info.seedEval  = sv.nEval;
       info.aimIter  = av.iterations;
       info.aimEval  = av.nEval;
         info.nProp  = sv.nEval + av.nEval + 1;
     info.solveInfo  = sv;
       info.aimInfo  = av;
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
        info.stopOK  = allOK;
       info.stopWhy  = {why1;why2;why3};

%% ...and the machinery itself, so an independent checker can re-integrate the
%% very same chain with a different solver instead of trusting this one. The
%% launch state handed out is the SOLVED one, taken from the flown trajectory
%% rather than rebuilt from the template, so it carries the solved heading and
%% cannot drift from what was actually integrated:
        info.phases  = buildPhases(tCut,cfg);
           info.env  = env;
            info.x0  = traj.x(1,:).';
      info.boostVeh  = bst;
      info.glideVeh  = glideVeh;

%% Vertical exaggeration for the globe. A DISPLAY choice, not a physical
%% constant, so it does not belong in missileConst. A true-scale arc would lie
%% on the surface and show nothing of the skip phugoid -- a 50 km glide on a
%% 6378 km sphere is one part in 128.
%%
%% TRUE SCALE IS NOW THE SHIPPED DEFAULT, altExag = 1. The movie carries a
%% true-scale altitude inset, which is precisely what the globe used to have to
%% exaggerate for, and coorbital.viz.globeMovie drops the "(altitude exaggerated
%% Nx)" caption clause at unity, so a true-scale picture makes no claim it is
%% not keeping.
%%
%% The ADAPTIVE rule is one word away -- altExag = 'auto' -- and it is what to
%% ask for when the skip phugoid is the point of the picture. A FIXED
%% exaggeration does not travel: 30x suits a glide that stays under about 60 km,
%% but on a lofted intercontinental shot peaking near 260 km it paints the arc
%% 7800 km off the surface, further out than the Earth is wide, and the track
%% leaves the frame. exagFor scales to the flown apogee instead. On this
%% script's own default, which peaks at 118.09 km, it returns 16x -- 1889 km of
%% apparent altitude, the 0.3 rE the rule aims at -- where 30x would paint the
%% same arc 3543 km off the surface and make a glide read as an orbit:
           hPeakM = max(traj.x(:,1)) - c.rE;
    if ischar(altExag) || isstring(altExag)
          altExag = exagFor(hPeakM,c.rE);
    end

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

function [traj,latIm,lonIm] = flyChain(psiL,tCut,cfg)
%% Purpose:
%
%  Fly the whole boost-glide-descent chain for ONE control pair -- a launch
%  azimuth and a thrust-termination time -- and hand back the trajectory and
%  the impact point. This is the only place in the file that touches the
%  physics: both solvers reach it, the seed bisection through flyRange and the
%  two-axis Newton through flyMiss, and the final answer is re-flown through it
%  as well. One propagation, one completeness check, one definition of where
%  the vehicle landed.
%
%  THE AZIMUTH IS WRITTEN INTO THE STATE, NOT REBUILT AROUND IT. cfg.x0 is a
%  template and entry 6 is the heading, so a trial flight differs from the seed
%  in the azimuth and in nothing else. Assembling a fresh launch state here
%  would be a second copy of the six lines that build it, free to drift from
%  the first.
%
%  A PROPAGATION THAT DID NOT COMPLETE IS AN ERROR, NOT A SHORT RANGE. If a
%  phase runs out of horizon -- the usual cause being a cutoff so early that
%  the vehicle never climbs to the handoff altitude, so the descending-crossing
%  event never fires -- the trajectory still ends SOMEWHERE and still has an
%  impact point. Returning it would feed the bisection a number from a flight
%  that never happened and quietly break the monotonicity that method rests on.
%  It throws instead. The two solvers read that throw differently and both
%  readings are right: rangeSolve has a bracket to protect and lets it
%  propagate, while aimSolve treats an infeasible control pair as a trial to
%  shorten and halves its step.
%
%% Inputs:
%
%  psiL             [1 x 1]                     Launch azimuth (rad),
%                                               clockwise from north
%
%  tCut             [1 x 1]                     Thrust-termination time (s)
%
%  cfg              Struct                      Immutable configuration; see
%                                               buildPhases, plus x0 [7 x 1]
%                                               (the launch-state template),
%                                               bst (struct), env (struct),
%                                               rE (m), lat0 (rad), lon0 (rad)
%
%% Outputs:
%
%  traj             Struct                      The trajectory itself, from
%                                               coorbital.prop.phaseRun
%
%  latIm            [1 x 1]                     Impact latitude (rad)
%
%  lonIm            [1 x 1]                     Impact longitude (rad)
%
%% Revision History:
%  Michael Casey                                                08/08/2026
%  Copyright 2026 Coorbital, Inc.
%% ------------------------ Begin Code Sequence ---------------------------

              xLau = cfg.x0;
           xLau(6) = psiL;
                ph = buildPhases(tCut,cfg);
              traj = coorbital.prop.phaseRun(ph,xLau,cfg.bst,cfg.env);

%% All three phases must have run, and the last must have stopped ON the
%% impact altitude rather than at its horizon:
             hEndM = traj.x(end,1) - cfg.rE;
             nPhRn = numel(unique(traj.phaseIdx));
    if nPhRn < 3 || abs(hEndM - cfg.hStop) > 1e-3
        error('coorbital:runTarget:propagationIncomplete', ...
            ['A cutoff at t = %.4f s on a %.6f deg azimuth produced %d of 3 ' ...
             'phases and ended at h = %.3f km against a %.3f km stop ' ...
             'altitude. The flight did not complete, so where it stopped is ' ...
             'not an impact point. The usual cause is a cutoff too early for ' ...
             'the vehicle to reach the %.3f km handoff at all, so the ' ...
             'descending-crossing event never fires: raise cutFracMin.'], ...
            tCut,rad2deg(psiL),nPhRn,hEndM./1000,cfg.hStop./1000, ...
            cfg.hHandoff./1000);
    end

             latIm = traj.x(end,3);
             lonIm = traj.x(end,2);
end

function rngM = flyRange(psiL,tCut,cfg,rI)
%% Purpose:
%
%  Great-circle surface range from the launch point to the impact point, for
%  one control pair. This is the scalar function coorbital.util.rangeSolve
%  bisects on during stage 1, with the azimuth held at the seed.
%
%  ONE SCALAR OUT AND NOTHING ELSE. rangeSolve's contract is scalar in, scalar
%  out, and a range function that also handed back a trajectory would invite a
%  caller to keep the last one and skip the re-flight the summary is built
%  from. The trajectory is flyChain's to return, and the answer is re-flown
%  once at the end where its repeatability is asserted.
%
%% Inputs:
%
%  psiL             [1 x 1]                     Launch azimuth (rad),
%                                               clockwise from north
%
%  tCut             [1 x 1]                     Thrust-termination time (s)
%
%  cfg              Struct                      Immutable configuration; see
%                                               flyChain
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
%% Revision History:
%  Michael Casey                                                08/07/2026
%  Michael Casey  Azimuth is now an argument, not a fixed
%                 property of cfg.x0                            08/08/2026
%  Copyright 2026 Coorbital, Inc.
%% ------------------------ Begin Code Sequence ---------------------------

   [~,latIm,lonIm] = flyChain(psiL,tCut,cfg);
              rngM = rI.*coorbital.util.greatCircle(cfg.lat0,cfg.lon0, ...
                                                    latIm,lonIm);
end

function fMiss = flyMiss(xAim,cfg,aim)
%% Purpose:
%
%  The two-component residual coorbital.util.aimSolve drives to zero in stage
%  2: fly the chain for one control pair and return the miss vector from the
%  target to the impact point, resolved in the aim frame. Sibling of flyRange
%  and built on the same flyChain, differing only in what it measures at the
%  end -- a distance from the LAUNCH point there, a two-component offset from
%  the TARGET here.
%
%  WHY A VECTOR MISS AND NOT A RANGE RESIDUAL BESIDE A CROSS-TRACK OFFSET. The
%  two carry the same information and one of them is better. Both components of
%  this one zero is exactly "the impact point IS the target", with no reference
%  to a great circle the flown track no longer follows once the Earth turns or
%  a bank is commanded; and the Jacobian is nearly diagonal, because the
%  azimuth moves the impact point almost purely crosswise and the cutoff almost
%  purely along the course. A range-residual pairing would keep quoting an arc
%  the vehicle is not on and would mix both controls into both residuals.
%
%% Inputs:
%
%  xAim             [2 x 1]                     Controls: the launch azimuth
%                                               (rad) and the
%                                               thrust-termination time (s),
%                                               in that order
%
%  cfg              Struct                      Immutable configuration; see
%                                               flyChain
%
%  aim              Struct                      The aim frame; see missVector
%
%% Outputs:
%
%  fMiss            [2 x 1]                     Down-range and cross-range
%                                               components of the miss (m);
%                                               see missVector. Throws
%                                               coorbital:runTarget:propagationIncomplete
%                                               for a control pair that does
%                                               not fly, which aimSolve treats
%                                               as an infeasible trial
%
%% Revision History:
%  Michael Casey                                                08/08/2026
%  Copyright 2026 Coorbital, Inc.
%% ------------------------ Begin Code Sequence ---------------------------

   [~,latIm,lonIm] = flyChain(xAim(1),xAim(2),cfg);
             fMiss = missVector(latIm,lonIm,aim);
end

function fMiss = missVector(latIm,lonIm,aim)
%% Purpose:
%
%  The miss vector from the TARGET to the IMPACT POINT, as surface distances
%  along and across the intended course. Zero in both components is exactly
%  coincidence, and hypot of the two is exactly the great-circle miss distance,
%  so the decomposition loses nothing and adds a direction.
%
%  THE COMPONENTS ARE FORMED DIRECTLY, NOT FROM A BEARING AND A DISTANCE. The
%  pair (wNrth,wEast) below is the same pair coorbital.util.greatCircleBearing
%  forms before it takes their arctangent: the northward and eastward parts of
%  the direction from the target to the impact point, both scaled by the sine
%  of the central angle between them. Taken as COMPONENTS they are smooth
%  straight through coincidence, which is where a bearing does not exist at all
%  -- greatCircleBearing rightly REFUSES there, and a residual that threw as it
%  converged would have the Newton halving its way into an infeasible point
%  every time it succeeded.
%
%  THE MAGNITUDE COMES FROM THE HAVERSINE, NOT FROM ASIN. hypot(wNrth,wEast) is
%  sin(Delta), and asin of it would fold every central angle past 90 degrees
%  back on itself -- reachable early in a solve, when the trial pair can land
%  most of a hemisphere away. coorbital.util.greatCircle gives Delta over the
%  whole range, so the components are rescaled by Delta/sin(Delta) instead: the
%  direction is taken from the pair, the length from the haversine. The ratio
%  tends to one at coincidence, which is why the guard below hands back unity
%  rather than dividing by zero.
%
%  SIGNS. Down-range is positive LONG -- the impact point is beyond the target
%  along the course -- and cross-range is positive RIGHT of that course.
%
%% Inputs:
%
%  latIm            [1 x 1]                     Impact latitude (rad)
%
%  lonIm            [1 x 1]                     Impact longitude (rad)
%
%  aim              Struct                      The aim frame, built once by
%                                               the caller and constant for the
%                                               whole solve:
%                                               latT [1 x 1] target latitude
%                                                    (rad)
%                                               lonT [1 x 1] target longitude
%                                                    (rad)
%                                               psiC [1 x 1] course of the
%                                                    launch-to-target great
%                                                    circle AT THE TARGET
%                                                    (rad), clockwise from
%                                                    north
%                                               rI   [1 x 1] radius of the
%                                                    impact sphere (m)
%
%% Outputs:
%
%  fMiss            [2 x 1]                     [down-range; cross-range]
%                                               components of the miss (m)
%
%% References:
%   [1] Bowditch, N., "The American Practical Navigator," Pub. No. 9, NGA,
%       chapter on Great-Circle Sailing. The four-parts course components are
%       the same ones coorbital.util.greatCircleBearing builds.
%
%% Revision History:
%  Michael Casey                                                08/08/2026
%  Copyright 2026 Coorbital, Inc.
%% ------------------------ Begin Code Sequence ---------------------------

             dLonI = lonIm - aim.lonT;
             wEast = sin(dLonI).*cos(latIm);
             wNrth = cos(aim.latT).*sin(latIm) - ...
                     sin(aim.latT).*cos(latIm).*cos(dLonI);

%% Direction from the pair, length from the haversine, joined by the ratio
%% Delta/sin(Delta). The guard is on the DIVISOR and not on the answer, and it
%% catches two geometries rather than one -- sin(Delta) vanishes at both ends
%% of the arc:
%%
%%    AT COINCIDENCE, Delta = 0, the ratio tends to unity and the components
%%    are zero anyway, so unity is the right value and the branch is exact.
%%
%%    AT THE ANTIPODE, Delta = pi, sin(Delta) is zero again but the true miss
%%    is rI*pi -- a hemisphere -- and this branch would report an exact hit.
%%    It is left unguarded deliberately, and the reason is worth stating
%%    because the obvious fix is to add a second test. The zero is ANALYTIC,
%%    not numerical: at the exact antipode sin(pi) evaluates to 1.22e-16 rather
%%    than to zero, and the two components come back finite, so sMag > 0 holds
%%    and the ratio pi/1.22e-16 is taken. Measured on that geometry the
%%    magnitude is right to 0.19 m in 20 000 km, with a direction that is
%%    meaningless because every great circle joins an antipodal pair -- the
%%    same ill-conditioning coorbital.util.greatCircleBearing refuses on a
%%    1e-12 test. An elseif here would therefore be a branch no double can
%%    reach, and an unreachable guard reads as protection that was never
%%    exercised. What actually keeps this problem away from the antipode is
%%    the envelope refusal: a target that far is unreachable and never gets
%%    as far as the Newton.
              sMag = hypot(wNrth,wEast);
              dAng = coorbital.util.greatCircle(aim.latT,aim.lonT,latIm,lonIm);
               kSc = 1;
    if sMag > 0
               kSc = dAng./sMag;
    end
              arcN = aim.rI.*kSc.*wNrth;
              arcE = aim.rI.*kSc.*wEast;

%% Rotate north-east into along-across. A course psiC clockwise from north
%% points along cos(psiC) north + sin(psiC) east, and right of it is the same
%% turned a further quarter turn:
             fMiss = [ arcN.*cos(aim.psiC) + arcE.*sin(aim.psiC); ...
                      -arcN.*sin(aim.psiC) + arcE.*cos(aim.psiC)];
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
