function [traj,info] = run_ballistic_target(opts)
%% Purpose:
%
%  Fly a BALLISTIC missile from a LAUNCH POINT to a DESTINATION POINT. Give it
%  two latitude/longitude pairs and it solves the trajectory that connects
%  them, flies the full boost-coast-descent chain, reports what it did and
%  draws the figures. Everything a routine run needs to change lives in the
%  USER PARAMETERS block; nothing below it should require editing.
%
%  Angles and distances in the USER PARAMETERS block are in degrees and
%  kilometres because that is how a user thinks about them. They are converted
%  to the library's SI units (m, m/s, rad, s) immediately after the block, in
%  one place, and nothing below that point works in any other unit.
%
%% Note -- the targeting problem is TWO separate solves:
%
%  BM/run_ballistic flies launch site + azimuth + pitch program and lands
%  wherever the physics puts it. This script inverts that. Two unknowns, two
%  independent solves, and they do not interact on a non-rotating Earth:
%
%    Azimuth, CLOSED FORM.  coorbital.util.greatCircleBearing returns the
%                           initial bearing of the great-circle arc from the
%                           launch point to the target, clockwise from north,
%                           which is exactly this library's heading convention
%                           psi. It drops straight into the launch state. No
%                           iteration, no propagation, no tolerance. It REFUSES
%                           coincident, antipodal and polar geometries rather
%                           than returning a plausible azimuth for them, and
%                           the refusal is caught below and re-phrased against
%                           the USER PARAMETERS entries that caused it.
%
%    Range, ONE-DIMENSIONAL ROOT SOLVE ON THE LOFT ANGLE.  See the next note.
%                           This is where a ballistic missile differs from the
%                           boost-glide vehicle of HGV/run_target, and it is
%                           the whole reason this script exists separately.
%
%% Note -- the ranging control is the LOFT ANGLE, and it has TWO BRANCHES:
%
%  HGV/run_target ranges on THRUST-TERMINATION TIME, because range is monotonic
%  and single-valued in it: less burn, less energy, less range. Bisection is
%  safe on a monotonic function and that is the whole reason it was chosen
%  there.
%
%  A ballistic missile is not flown that way. The classical control is the LOFT
%  ANGLE -- the terminal attitude of the pitch program -- and range is NOT
%  monotonic in it. Range rises to a MAXIMUM at some max-range loft angle and
%  falls away on both sides, so every range short of that maximum is reached by
%  TWO trajectories:
%
%      the LOFTED arc,    above the max-range angle: a high, slow, long-flight
%                         trajectory that arrives steeply;
%      the DEPRESSED arc, below it: a low, fast, short-flight trajectory that
%                         arrives shallow.
%
%  A ballistic user expects to CHOOSE. So this script always brackets the
%  max-range loft angle FIRST, then solves BOTH branches, reports both, and
%  flies the one the branch selector asks for. The two rangeSolve calls are
%  each made on a bracket that lies entirely on ONE side of the maximum, which
%  is what restores the monotonicity bisection needs -- see the MONOTONICITY IS
%  ASSUMED, NOT CHECKED note in coorbital.util.rangeSolve.
%
%  THE BRANCH IS MEASURED, NOT ASSUMED. A bracket that lies on one side of the
%  maximum SHOULD hold the solution on that side, but "should" is not evidence:
%  a max-range angle located slightly wrong, or a range function with a kink,
%  puts the answer on the other branch while the solve still converges. So the
%  branch of the FLOWN trajectory is measured after the fact, from its apogee
%  and its flight time against the max-range arc's own apogee and flight time,
%  and the summary reports the measured branch beside the requested one. They
%  disagreeing is a CAUTION, printed, not an assertion swallowed.
%
%% Note -- what "minimum-energy" can and cannot mean for a FIXED booster:
%
%  READ THIS BEFORE TRUSTING THE branch = 'minimum-energy' DEFAULT. The
%  textbook minimum-energy ballistic trajectory is the one that reaches a
%  required range with the LEAST burnout energy, and its burnout flight-path
%  angle is gammaStar = 45 deg - Lambda/4 for a free-flight range angle Lambda.
%  That statement presumes the burnout energy is FREE TO CHOOSE -- you size the
%  booster to the range.
%
%  This script cannot choose it. The booster carries a fixed propellant load
%  and burns to exhaustion on every branch, so BOTH arcs leave burnout with
%  essentially the same energy: measured on the shipped case, the two differ by
%  0.55 %%, and that difference is boost-phase drag and gravity loss, not a
%  design variable. No member of this family is the minimum-energy trajectory
%  unless the target happens to sit at the maximum range.
%
%  WHAT IS IMPLEMENTED, precisely: minimum-energy selects THE BRANCH WHOSE
%  SOLVED LOFT ANGLE LIES NEARER THE MAX-RANGE LOFT ANGLE. The reasoning is
%  that the max-range arc IS the minimum-energy arc for its own range -- that
%  is the same theorem read the other way round -- so the arc nearest it in the
%  control variable is the arc that comes nearest the minimum-energy geometry,
%  and the selection becomes exact as the target approaches maximum range.
%
%  THREE THINGS THIS RULE IS NOT, all of them PRINTED every run so the reader
%  can judge rather than take it on trust:
%
%    It is not "least burnout energy" read literally. That reading selects the
%    arc that WASTED more energy climbing out, because with the propellant load
%    fixed a lower burnout energy means a costlier ascent, not a cheaper
%    trajectory. On the shipped case it would select the DEPRESSED arc while
%    the rule above selects the lofted one. Both energies are reported.
%
%    It is not the gammaStar test. The burnout flight-path angle nearest
%    45 deg - Lambda/4 is the textbook criterion and it is reported too, but it
%    is measured on a burnout state this script does not control -- the alphaMax
%    clamp does, see the next note -- so it is a diagnostic here, not the rule.
%
%    It is not always well conditioned. The range-versus-loft hump is roughly
%    symmetric, so the two branch solutions can fall very nearly equidistant
%    from the max-range angle and the selection then turns on a rounding.
%    Measured on this vehicle, along the shipped azimuth: at 3175 km the two
%    distances are 39.43 deg and 34.08 deg, 14 %% apart, and the lofted arc
%    wins; by 3789 km they agree to 4 %% and the caution below fires; and by
%    4218 km the ranking has FLIPPED and the depressed arc wins by 10 %%. The
%    summary prints both distances and CAUTIONS when they fall within 5 %% of
%    each other, and the advice it gives then is the right one: ask for
%    'lofted' or 'depressed' explicitly.
%
%% Note -- the loft angle is a COMMANDED ATTITUDE, not the burnout gamma:
%
%  The loft angle is the last node of the commanded pitch-attitude schedule.
%  What the vehicle actually achieves at burnout is a different number, and on
%  this configuration it is a long way different: the shipped alphaMax clamp
%  limits how fast the flight path can be pushed over, so a commanded terminal
%  attitude of -14.3 deg produces a burnout flight-path angle of +11.8 deg. The
%  summary prints BOTH for every branch. Read the loft angle as the control it
%  is, not as a prediction of the burnout state.
%
%  THE CLAMP IS ALSO WHY alphaMax IS 12 deg HERE AND 6 deg IN BM/run_ballistic,
%  and that is a targeting decision rather than an aerodynamic one. At 6 deg the
%  clamp binds over the WHOLE loft bracket: a commanded -30 deg still burns out
%  at 33.4 deg of flight-path angle, above the max-range value, so the max-range
%  angle is never reached, range is monotone decreasing in loft across the
%  entire bracket, and THE DEPRESSED BRANCH DOES NOT EXIST. Measured, by
%  scanning the shipped bracket at 6, 8, 10 and 12 deg: the same commanded
%  -30 deg gives 23.6, 15.0 and 7.6 deg at the other three. The two-branch structure this script is built around
%  needs enough pitch-over authority to fly past the max-range angle, and the
%  bracketing step below refuses with that message if it is not there.
%
%% Note -- TWO LIMITATIONS, both printed in the summary rather than buried:
%
%  THE AZIMUTH IS EXACT ONLY BECAUSE THE EARTH DOES NOT ROTATE HERE. With
%  env.omegaE = 0 the ground track of an unbanked trajectory is a great circle,
%  so the initial bearing of the launch-to-target arc is the whole answer.
%  Switch earthSpin on and the target is carried east under the vehicle for the
%  whole flight, so the azimuth needs an OUTER ITERATION around the range solve
%  and this script does not have one. It refuses to pretend: with earthSpin
%  true it prints a caution and the reported miss is against a target that
%  stood still. The two branches make that error VISIBLE rather than abstract,
%  because they fly for very different times -- 659 s and 1816 s on the shipped
%  case -- so the same geometry drifts by very different amounts.
%
%  THE SOLVE CONTROLS DOWNRANGE ONLY. Bisection matches the great-circle
%  distance from launch to impact; nothing in it steers the track sideways.
%  Cross-range comes out at zero in the SHIPPED configuration because the bank
%  angle is zero throughout, and an unbanked trajectory over a non-rotating
%  sphere stays on the great circle it left on. That is a property of THIS
%  CONFIGURATION, not a general guarantee, and this script MEASURES it rather
%  than asserting it.
%
%% Note -- two vehicles, one chain, carried on ph.veh:
%
%  The boosted stack and the separated re-entry body have different reference
%  areas, different aerodynamic coefficients and different masses, and a phase
%  link can move the STATE across separation but can never replace the
%  PARAMETERS the equations divide by. coorbital.prop.phaseRun therefore accepts
%  an optional per-phase ph.veh, and this script uses it: phase 1 carries the
%  booster, phases 2 and 3 carry the re-entry body, and the equations of motion
%  are the library handles themselves rather than closures that ignore the
%  forwarded argument. The older entry scripts bind their vehicle inside the
%  closure instead; both work, and coorbital.eom.massConstant still guards the
%  mass half of the divergence either way.
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
%  traj             Struct                      Trajectory of the FLOWN branch,
%                                               from coorbital.prop.phaseRun:
%                                               t        [N x 1] (s)
%                                               x        [N x 7] state, the
%                                                        seventh component
%                                                        being mass (kg)
%                                               u        [N x 2] control (rad)
%                                               phaseIdx [N x 1]
%                                               junction [2 x 1] struct; the
%                                                        first is the state
%                                                        AFTER separation
%                                               EMPTY when the request is
%                                               refused -- see info.refused
%
%  info             Struct                      The summary's numbers at full
%                                               precision, so a test does not
%                                               have to read them back out of
%                                               printed text; all SI except the
%                                               *Km distances and the *Deg
%                                               angles. Always carries refused
%                                               (logical), rngReqM, psiLaunch,
%                                               loftStarDeg, rngMaxM, rngMinM
%                                               and the per-branch records in
%                                               depressed and lofted, so a
%                                               refused run still hands back the
%                                               envelope it refused against. A
%                                               run that was not refused
%                                               additionally carries the flown
%                                               branch's own numbers and the
%                                               DISPLAY choices the figures were
%                                               drawn with -- altExag, its rule
%                                               as a handle, and the two
%                                               hemisphere captions
%
%% References:
%   [1] Bate, R.R., Mueller, D.D., White, J.E., "Fundamentals of
%       Astrodynamics," Dover, 1971, Ch. 6. The free-flight range angle of a
%       Keplerian arc, the two burnout flight-path angles that reach a given
%       range, and the minimum-energy value gammaStar = 45 deg - Lambda/4 that
%       separates them.
%   [2] Bowditch, N., "The American Practical Navigator," Pub. No. 9, NGA,
%       chapter on Great-Circle Sailing. The initial course used for the launch
%       azimuth; see coorbital.util.greatCircleBearing.
%   [3] Press, W.H., et al., "Numerical Recipes," 3rd ed., Cambridge, 2007,
%       Sections 9.1 and 10.2. Bisection, used by coorbital.util.rangeSolve,
%       and golden-section search, used below to bracket the maximum.
%   [4] Vinh, N.X., Busemann, A., Culp, R.D., "Hypersonic and Planetary Entry
%       Flight Mechanics," Univ. Michigan Press, 1980.
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
%% Launch point -- where the missile leaves the pad:
         latLaunch = 45;               %deg, geocentric latitude [-89 .. 89]; the poles are singular
         lonLaunch = -100;             %deg, longitude [-180 .. 180]
           hLaunch = 0;                %km,  pad altitude above the sphere [0 .. 5]
           vLaunch = 10;               %m/s, speed at the first integrated point [2 .. 50].
                                       %     NOT a physical launch speed: the 3DOF equations
                                       %     are singular at V = 0, so the integration starts
                                       %     a moment after first motion

%% Target point -- where the missile is to arrive. There is NO azimuth entry:
%% the launch azimuth is SOLVED from these two points, in closed form, and any
%% azimuth given here could only contradict it:
         latTarget = 62;               %deg, geocentric latitude [-89 .. 89]
         lonTarget = -60;              %deg, longitude [-180 .. 180]
                                       %     The shipped pair is the central United States to
                                       %     Baffin Island: 3175 km on a 39.2 deg azimuth.
                                       %     Deliberately NOT equatorial and NOT due east --
                                       %     that geometry makes latitude and longitude
                                       %     interchangeable in the range formula and hides a
                                       %     transposition at every call site that uses them

%% WHICH OF THE TWO ARCS TO FLY. Both are solved and both are reported whatever
%% is asked for here; this chooses only which one is flown, plotted and returned:
            branch = 'minimum-energy'; %'minimum-energy' | 'lofted' | 'depressed'
                                       %  minimum-energy: the branch whose loft angle lies
                                       %                  NEARER the max-range loft angle. Read
                                       %                  the "what minimum-energy can and
                                       %                  cannot mean" note in the header before
                                       %                  relying on it -- with a fixed
                                       %                  propellant load burned to exhaustion
                                       %                  neither branch is minimum-energy in
                                       %                  the textbook sense, and the summary
                                       %                  prints every number needed to see that
                                       %  lofted:         the high, slow, long-flight arc
                                       %  depressed:      the low, fast, short-flight arc

%% Loft-angle bracket -- the search interval for BOTH the max-range angle and
%% the two branch solves. The loft angle is the LAST node of the commanded pitch
%% schedule; see the header for why it is not the burnout flight-path angle. The
%% max-range angle must lie strictly INSIDE this bracket or there is no
%% two-branch structure to solve, and the bracketing step below says so:
           loftMin = -40;              %deg, most depressed commanded terminal attitude
                                       %     [-60 .. 0]. Negative means the schedule commands
                                       %     the nose BELOW the horizon at end of burn, which
                                       %     the alphaMax clamp then declines to deliver in
                                       %     full: at -40 the achieved burnout gamma is +5.2 deg.
                                       %     Below about -50 the vehicle burns out already
                                       %     descending, the apogee event never fires and the
                                       %     propagation is rejected
           loftMax = 85;               %deg, most lofted commanded terminal attitude
                                       %     (loftMin .. 89). 89 is the pad attitude itself,
                                       %     i.e. no pitch-over at all and no downrange travel
        tolLoftDeg = 0.05;             %deg, bracket width the golden-section search closes the
                                       %     max-range angle to [0.005 .. 1]. The hump is flat
                                       %     at its top -- about 8 km of range per square
                                       %     degree -- so 0.05 deg is 0.01 km of range, far
                                       %     inside the range tolerance below
         nScanLoft = 13;               %-, coarse scan points across the bracket before the
                                       %   golden section starts [7 .. 41]. The scan is what
                                       %   makes the search safe: golden section needs a
                                       %   unimodal bracket, and the scan both finds one and
                                       %   provides the evidence that the range curve has a
                                       %   single interior hump
        tolRangeKm = 1.0;              %km, convergence tolerance on the ACHIEVED range of each
                                       %    branch [0.05 .. 50]. Each halving of this costs one
                                       %    more trajectory propagation per branch, and a
                                       %    propagation is about 0.05 s, so tightening it is cheap

%% Pitch program -- commanded pitch ATTITUDE against time since liftoff. ONLY
%% ITS SHAPE IS FLOWN. The profile is rescaled every propagation so that its
%% first node stays at pitchRef(1) and its LAST node becomes the solved loft
%% angle, with the intermediate nodes moved in proportion; at a loft angle equal
%% to pitchRef(end) it reproduces this reference exactly, and the shipped
%% reference IS BM/run_ballistic's pitch program, so the two scripts fly the
%% same schedule at loft = 34 deg:
         pitchTime = [0  6  15  30  50  70  82];   %s,   nodes, strictly increasing [0 .. burn time]
          pitchRef = [89 86 76 60 46 38 34];       %deg, reference attitudes [0 .. 90], strictly
                                       %     descending. pitchRef(1) is the pad attitude and IS
                                       %     flown as written; pitchRef(end) sets the shape only
          alphaMax = 12;               %deg, clamp on |angle of attack| [1 .. 15]. THIS CLAMP,
                                       %     NOT THE SCHEDULE, is what limits how fast the
                                       %     flight path can be pushed over, and at
                                       %     BM/run_ballistic's 6 deg the depressed branch does
                                       %     not exist at all -- see the header note
         bankAngle = 0;                %deg, bank [-90 .. 90], all phases. MUST BE 0 for the
                                       %     closed-form azimuth to be the answer; anything else
                                       %     turns the heading and walks the impact point off
                                       %     the launch-to-target great circle, which the solve
                                       %     cannot see. Measured and warned about, not assumed

%% Staging -- whether the spent booster is thrown away at burnout:
        separation = true;             %true jettisons bst.massDry and flies the re-entry body
                                       %     alone; false keeps the dead booster attached, which
                                       %     is a legitimate configuration and flies the whole
                                       %     coast on the STACK mass and stack aerodynamics

%% Termination -- when to stop integrating each phase:
             hStop = 0;                %km, impact altitude, DESCENDING crossing only [0 .. 20];
                                       %    a target on high terrain is a legitimate setting
         tMaxBoost = 200;              %s,  boost horizon; must exceed the burn time
         tMaxCoast = 3000;             %s,  ascending-coast horizon, burnout to apogee
          tMaxDesc = 3000;             %s,  descending horizon, apogee to impact

%% Vehicle and booster -- point these at any file returning the right struct:
         vehicleFn = @vehicle_bm;                        %handle, see vehicle_bm.m
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
         movieFile = fullfile(tempdir,'run_ballistic_target.mp4');   %char, MP4 output path.
                                       %Defaults into tempdir on purpose: a movie is a build
                                       %artefact and does not belong in the source tree
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
                      'latTarget','lonTarget','branch','loftMin','loftMax', ...
                      'tolLoftDeg','nScanLoft','tolRangeKm','pitchTime', ...
                      'pitchRef','alphaMax','bankAngle','separation','hStop', ...
                      'tMaxBoost','tMaxCoast','tMaxDesc','vehicleFn', ...
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
            branch = overrideOf(opts,'branch',branch);
           loftMin = overrideOf(opts,'loftMin',loftMin);
           loftMax = overrideOf(opts,'loftMax',loftMax);
        tolLoftDeg = overrideOf(opts,'tolLoftDeg',tolLoftDeg);
         nScanLoft = overrideOf(opts,'nScanLoft',nScanLoft);
        tolRangeKm = overrideOf(opts,'tolRangeKm',tolRangeKm);
         pitchTime = overrideOf(opts,'pitchTime',pitchTime);
          pitchRef = overrideOf(opts,'pitchRef',pitchRef);
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
          movieOn  = overrideOf(opts,'movieOn',movieOn);
       movieFrames = overrideOf(opts,'movieFrames',movieFrames);
         movieFile = overrideOf(opts,'movieFile',movieFile);

%% Convert the user block to library SI units. This is the ONLY unit
%% conversion in the file; everything past this point is m, m/s, rad and s:
          hLaunchM = hLaunch.*1000;
            hStopM = hStop.*1000;
           tolRngM = tolRangeKm.*1000;
        latLaunchR = deg2rad(latLaunch);
        lonLaunchR = deg2rad(lonLaunch);
        latTargetR = deg2rad(latTarget);
        lonTargetR = deg2rad(lonTarget);
         pitchTimS = pitchTime(:).';
         pitchRefR = deg2rad(pitchRef(:).');
         alphaMaxR = deg2rad(alphaMax);
           bankRad = deg2rad(bankAngle);
           loftLoR = deg2rad(loftMin);
           loftHiR = deg2rad(loftMax);
          tolLoftR = deg2rad(tolLoftDeg);

%% Sanity-check the user block before spending time in the integrator:
    assert(numel(pitchTimS) == numel(pitchRefR), ...
        'pitchTime has %d nodes and pitchRef has %d; they must match.', ...
        numel(pitchTimS),numel(pitchRefR));
    assert(numel(pitchTimS) >= 2,'the pitch program needs at least two nodes.');
    assert(all(diff(pitchTimS) > 0),'pitchTime must be strictly increasing.');
    assert(all(diff(pitchRefR) < 0), ...
        ['pitchRef must be strictly DESCENDING; its shape is what is flown and ' ...
         'a non-monotone reference makes the rescaling below meaningless.']);
    assert(pitchRefR(1) > pitchRefR(end), ...
        ['pitchRef(1) must exceed pitchRef(end); their difference is the ' ...
         'denominator of the shape normalisation.']);
    assert(loftLoR < loftHiR, ...
        'the loft bracket must satisfy loftMin < loftMax; got %.4f and %.4f deg.', ...
        loftMin,loftMax);
    assert(loftHiR < pitchRefR(1), ...
        ['loftMax (%.4f deg) must stay below the pad attitude pitchRef(1) ' ...
         '(%.4f deg); at or above it the schedule stops descending and there ' ...
         'is no pitch-over left to command.'],loftMax,pitchRef(1));
    assert(tolLoftR > 0,'tolLoftDeg must be positive.');
    assert(nScanLoft >= 5 && nScanLoft == fix(nScanLoft), ...
        'nScanLoft must be a whole number of at least 5; got %s.',mat2str(nScanLoft));
    assert(tolRngM > 0,'tolRangeKm must be positive.');
    assert(hStopM >= 0, ...
        ['hStop (%.1f km) is below the spherical datum; the descent ends on a ' ...
         'descending crossing of that altitude and a negative one has no ' ...
         'meaning against this Earth model.'],hStop);
    assert(hStopM < hLaunchM + 1000, ...
        'hStop (%.1f km) is at or above the pad; the run would end at t = 0.',hStop);
    assert(vLaunch > 1, ...
        'vLaunch must exceed 1 m/s; the equations of motion are singular below that.');
    assert(all([tMaxBoost tMaxCoast tMaxDesc] > 0),'every phase horizon must be positive.');
    assert(islogical(separation) || isnumeric(separation), ...
        'separation must be true or false.');
    assert(isscalar(latLaunch) && isscalar(lonLaunch) && ...
           isscalar(latTarget) && isscalar(lonTarget), ...
        ['the launch point and the target must each be ONE latitude and ONE ' ...
         'longitude; this script solves a single trajectory, and a vector ' ...
         'here would silently broadcast through every great-circle call.']);
    assert(abs(latLaunch - latTarget) + abs(lonLaunch - lonTarget) > 0, ...
        ['the launch point and the target are the same point; there is nothing ' ...
         'to solve. coorbital.util.greatCircleBearing would refuse this a few ' ...
         'lines further down, but it can only speak of coincident points, not ' ...
         'of latTarget and lonTarget.']);

%% The branch selector, checked against the three legal words rather than left
%% to fall through a switch's otherwise clause: a misspelling here would
%% otherwise choose a branch the user did not ask for and never say so:
         branchLeg = {'minimum-energy','lofted','depressed'};
    assert(ischar(branch) || isstring(branch), ...
        'branch must be one of ''minimum-energy'', ''lofted'' or ''depressed''.');
            branch = char(branch);
    assert(any(strcmp(branch,branchLeg)), ...
        ['branch = "%s" is not one of ''minimum-energy'', ''lofted'' or ' ...
         '''depressed''. There is no default: every range short of the maximum ' ...
         'is reached by two different trajectories and guessing which one was ' ...
         'wanted is exactly what this parameter exists to prevent.'],branch);

                 c = coorbital.util.missileConst();
               veh = vehicleFn();
               bst = boosterFn();

%% Mass bookkeeping -- the state mass is ALWAYS the total mass carried; see
%% coorbital.util.boosterDefaults:
          mLiftoff = veh.mass + bst.massDry + bst.massProp;
         mBurnoutT = veh.mass + bst.massDry;
    if separation
            mCoast = veh.mass;
    else
            mCoast = mBurnoutT;
    end

%% Full burn time, from the motor model itself rather than a literal, so a
%% different booster or a different propFn moves it. It bounds nothing here --
%% the boost always runs to propellant exhaustion, this script ranging on the
%% loft angle instead -- but the horizon must clear it or the horizon and not
%% the motor would end the burn:
         [~,mdot0] = propFn(0,0,bst);
    assert(mdot0 > 0, ...
        'the propulsion model reports zero mass flow; there is no burn to fly.');
             tBurn = bst.massProp./mdot0;
    assert(tMaxBoost > tBurn, ...
        ['tMaxBoost (%.1f s) must exceed the %.3f s full burn, or the horizon ' ...
         'and not propellant exhaustion ends the boost.'],tMaxBoost,tBurn);

%% Assemble the environment from the handles chosen above:
         env.atmos = atmosFn;
          env.grav = gravFn;
          env.aero = aeroFn;
          env.prop = propFn;
        env.omegaE = 0;
    if earthSpin
        env.omegaE = c.omegaE;
    end

%% The coast vehicle. Its aerodynamics are the re-entry body's when the booster
%% is jettisoned and the STACK's when it is not, and its mass field must equal
%% the mass actually carried. Rebuilding the mass alone would leave Sref, CL
%% and LD describing a jettisoned stack, which is the worse bug because it
%% looks repaired:
    if separation
          coastVeh = veh;
    else
          coastVeh = bst;
    end
     coastVeh.mass = mCoast;

%% The pitch-program SHAPE, normalised once. shape(1) is exactly 0 and
%% shape(end) exactly 1, so a schedule rebuilt from it starts at pitchRef(1)
%% and ends at whatever loft angle is asked for:
             shape = (pitchRefR(1) - pitchRefR)./(pitchRefR(1) - pitchRefR(end));

%% Guidance for the unpowered phases. The coast and the descent fly at zero
%% incidence, which for the near-symmetric re-entry body is what "ballistic"
%% means; the boost schedule is rebuilt per loft angle inside buildPhases:
        schedCoast = struct('tGrid',[0 max(tMaxCoast,tMaxDesc)], ...
                            'alpha',[0 0], ...
                            'sigma',[bankRad bankRad]);

%% -----------------------------------------------------------------
%% The azimuth solve: closed form, one call, no iteration
%% -----------------------------------------------------------------
%% The initial bearing of the launch-to-target great circle IS the library's
%% heading state psi -- both are measured clockwise from north -- so it goes
%% straight into the launch state below with no conversion.
%%
%% coorbital.util.greatCircleBearing REFUSES the geometries that have no
%% azimuth rather than returning the rounding artefact atan2 makes of them, and
%% the refusal arrives here as an error identifier. It is re-phrased against the
%% USER PARAMETERS entries that caused it, because the library function cannot
%% know that the two points it was handed came from latLaunch and latTarget:
    try
         psiLaunch = coorbital.util.greatCircleBearing(latLaunchR,lonLaunchR, ...
                                                       latTargetR,lonTargetR);
    catch errAz
        if strcmp(errAz.identifier,'coorbital:greatCircleBearing:degenerateArc')
            error('coorbital:runBallisticTarget:degenerateGeometry', ...
                ['There is no launch azimuth from (%.6f, %.6f) to (%.6f, %.6f): ' ...
                 'the two points are coincident or antipodal to within rounding. ' ...
                 'A point has no arc to itself, and every great circle joins an ' ...
                 'antipodal pair, so no azimuth exists to fly. Move latTarget or ' ...
                 'lonTarget. The library said: %s'], ...
                latLaunch,lonLaunch,latTarget,lonTarget,errAz.message);
        elseif strcmp(errAz.identifier,'coorbital:greatCircleBearing:polarOrigin')
            error('coorbital:runBallisticTarget:polarLaunch', ...
                ['The launch point (%.6f, %.6f) is at a pole, where the local ' ...
                 'north a bearing is measured from does not exist. Move ' ...
                 'latLaunch off the pole; the USER PARAMETERS block documents ' ...
                 'the usable band as -89 to 89 deg for exactly this reason. The ' ...
                 'library said: %s'],latLaunch,lonLaunch,errAz.message);
        else
            rethrow(errAz);
        end
    end

%% Liftoff state. Everything but the loft angle is fixed across every
%% propagation this script makes; the pad flight-path angle is the first
%% commanded pitch attitude, so the angle of attack starts at zero:
                x0 = [c.rE + hLaunchM; ...
                      lonLaunchR; ...
                      latLaunchR; ...
                      vLaunch; ...
                      pitchRefR(1); ...
                      psiLaunch; ...
                      mLiftoff];

%% Everything the propagation needs, gathered once so the range function is a
%% single-argument handle over an immutable configuration:
      cfg.eomCoast = coorbital.eom.massConstant(@coorbital.eom.glide3DOF);
      cfg.schedCst = schedCoast;
        cfg.pitchT = pitchTimS;
         cfg.shape = shape;
      cfg.thetaTop = pitchRefR(1);
      cfg.alphaMax = alphaMaxR;
      cfg.bankRad  = bankRad;
      cfg.mBurnout = mBurnoutT;
        cfg.mCoast = mCoast;
    cfg.separation = logical(separation);
      cfg.massDry  = bst.massDry;
         cfg.hStop = hStopM;
     cfg.tMaxBoost = tMaxBoost;
     cfg.tMaxCoast = tMaxCoast;
      cfg.tMaxDesc = tMaxDesc;
            cfg.x0 = x0;
           cfg.bst = bst;
      cfg.coastVeh = coastVeh;
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
            fRange = @(loftR) flyLoft(loftR,cfg,rI);

%% -----------------------------------------------------------------
%% Bracket the MAXIMUM-RANGE loft angle, before either branch is solved
%% -----------------------------------------------------------------
%% This is the step that has no counterpart in HGV/run_target, and everything
%% below depends on it: it splits the loft bracket into two intervals on each of
%% which range IS monotonic, which is the precondition bisection needs, and it
%% supplies the max-range arc whose apogee and flight time the flown branch is
%% MEASURED against afterwards:
[loftStarR,rngMax,mr] = maxRangeLoft(fRange,loftLoR,loftHiR,nScanLoft,tolLoftR);
       [~,trajStar]   = fRange(loftStarR);
              kApStar = find(trajStar.phaseIdx == 2,1,'last');
              hApoStar = trajStar.x(kApStar,1) - c.rE;
              tFlyStar = trajStar.t(end);

%% -----------------------------------------------------------------
%% Solve BOTH branches, each on its own side of the maximum
%% -----------------------------------------------------------------
%% Both are solved whatever the branch selector asks for, because the trade is
%% the point: a user choosing between a 659 s arrival at 924 m/s and an 1816 s
%% arrival at 3503 m/s needs both numbers in front of them, and only one of the
%% two is flown. coorbital.util.rangeSolve returns converged = false with the
%% achievable band rather than throwing, so a branch that cannot reach this
%% target simply reports that it does not exist:
               dep = solveBranch('depressed',rngReq,fRange,loftLoR,loftStarR, ...
                                 tolRngM,cfg,rI,c,latTargetR,lonTargetR, ...
                                 hApoStar,tFlyStar);
               lof = solveBranch('lofted'   ,rngReq,fRange,loftStarR,loftHiR, ...
                                 tolRngM,cfg,rI,c,latTargetR,lonTargetR, ...
                                 hApoStar,tFlyStar);

%% The reachable envelope, as the union of the two branch bands. Its ceiling is
%% the maximum range by construction -- both branches meet there -- and its
%% floor is whichever branch endpoint reaches least far:
            rngMin = min(dep.bandLo,lof.bandLo);

%% -----------------------------------------------------------------
%% Choose the branch, and say on what grounds
%% -----------------------------------------------------------------
%% dLoftDep and dLoftLof are the two distances the minimum-energy rule compares.
%% They are computed here, once, whatever the selector says, because the summary
%% prints them on every run: a reader deciding whether to trust the default
%% needs to see how close the call was:
          dLoftDep = NaN;
          dLoftLof = NaN;
    if dep.exists
          dLoftDep = abs(dep.loftR - loftStarR);
    end
    if lof.exists
          dLoftLof = abs(lof.loftR - loftStarR);
    end
          bothHere = dep.exists && lof.exists;
          meRatio  = NaN;
          meClose  = false;
    if bothHere && max(dLoftDep,dLoftLof) > 0
          meRatio  = min(dLoftDep,dLoftLof)./max(dLoftDep,dLoftLof);
          meClose  = meRatio > 0.95;
    end

%% The selection itself. minimum-energy takes the nearer branch when both
%% exist and the surviving one when only one does; a named branch takes that
%% branch or nothing:
          pickName = branch;
          pickWhy  = '';
         pickShort = 'the branch selector';
    switch branch
        case 'minimum-energy'
         pickShort = 'nearest the max-range angle';
            if bothHere && dLoftLof <= dLoftDep
          pickName = 'lofted';
          pickWhy  = sprintf(['nearer the %.4f deg max-range loft angle: %.4f ' ...
                              'deg against the depressed arc''s %.4f deg'], ...
                             rad2deg(loftStarR),rad2deg(dLoftLof),rad2deg(dLoftDep));
            elseif bothHere
          pickName = 'depressed';
          pickWhy  = sprintf(['nearer the %.4f deg max-range loft angle: %.4f ' ...
                              'deg against the lofted arc''s %.4f deg'], ...
                             rad2deg(loftStarR),rad2deg(dLoftDep),rad2deg(dLoftLof));
            elseif lof.exists
          pickName = 'lofted';
          pickWhy  = 'the only branch that reaches this target';
            elseif dep.exists
          pickName = 'depressed';
          pickWhy  = 'the only branch that reaches this target';
            else
          pickWhy  = 'neither branch reaches this target';
            end
        case 'lofted'
         pickShort = 'named in the user block';
          pickWhy  = 'the arc the user asked for by name, not a solver choice';
        case 'depressed'
         pickShort = 'named in the user block';
          pickWhy  = 'the arc the user asked for by name, not a solver choice';
    end
    if strcmp(pickName,'lofted')
              pick = lof;
    elseif strcmp(pickName,'depressed')
              pick = dep;
    else
              pick = dep;
    end

%% -----------------------------------------------------------------
%% Refuse, loudly, when the request cannot be flown
%% -----------------------------------------------------------------
%% Two different refusals, and they are told apart because they ask the user to
%% change two different things. NEITHER of them throws: a targeting script that
%% silently returned the nearest miss would be worse than one that refuses,
%% because the nearest miss LOOKS like a solution, and a script that threw
%% would force every caller into a try/catch to read numbers the returned
%% struct already carries:
          noBranch = ~dep.exists && ~lof.exists;
          badPick  = ~noBranch && ~pick.exists;
    if noBranch || badPick
            shortM = bandShortfall(rngReq,rngMin,rngMax);
        fprintf('\n');
        fprintf('===== Ballistic point-to-point targeting: REFUSED =======\n');
        if noBranch
        fprintf('  NEITHER arc reaches this target. No loft angle in the bracket does.\n');
        else
        fprintf('  The %s arc does not reach this target, though the other one does.\n', ...
                branch);
        end
        fprintf('\n');
        fprintf('    launch           %10.4f %-5s %10.4f deg  (lat, lon)\n', ...
                latLaunch,'deg',lonLaunch);
        fprintf('    target           %10.4f %-5s %10.4f deg  (lat, lon)\n', ...
                latTarget,'deg',lonTarget);
        fprintf('    required range   %10.3f %-5s (great circle on the r = %.3f km sphere)\n', ...
                rngReq./1000,'km',rI./1000);
        fprintf('    MAXIMUM RANGE    %10.3f %-5s (the best this vehicle can do, at a loft\n', ...
                rngMax./1000,'km');
        fprintf('                                       angle of %.4f deg)\n',rad2deg(loftStarR));
        fprintf('    shortfall        %+10.3f %-5s (required minus the nearer edge of the\n', ...
                shortM./1000,'km');
        fprintf('                                       %.3f to %.3f km reachable band)\n', ...
                rngMin./1000,rngMax./1000);
        fprintf('\n');
        fprintf('    depressed arc    %10.3f to %.3f km over loft %.4f to %.4f deg  [%s]\n', ...
                dep.bandLo./1000,dep.bandHi./1000,loftMin,rad2deg(loftStarR),dep.why);
        fprintf('    lofted arc       %10.3f to %.3f km over loft %.4f to %.4f deg  [%s]\n', ...
                lof.bandLo./1000,lof.bandHi./1000,rad2deg(loftStarR),loftMax,lof.why);
        fprintf('\n');
        if shortM > 0
        fprintf('  The target is BEYOND MAXIMUM RANGE by %.3f km. No loft angle can add to\n', ...
                shortM./1000);
        fprintf('  %.3f km: that IS the maximum, found by bracketing the hump rather than\n', ...
                rngMax./1000);
        fprintf('  assumed. Give the booster more propellant, raise the specific impulse,\n');
        fprintf('  stage it, or choose a nearer target. Widening loftMin or loftMax will not\n');
        fprintf('  help -- the maximum is INTERIOR to the bracket and both ends fly shorter.\n');
        elseif shortM < 0
        fprintf('  The target is TOO CLOSE. Even the shortest arc in the bracket overflies it\n');
        fprintf('  by %.3f km. Raise loftMax towards the %.1f deg pad attitude, which flies\n', ...
                -shortM./1000,pitchRef(1));
        fprintf('  almost straight up and almost nowhere downrange, or lower loftMin -- but\n');
        fprintf('  below about -50 deg the vehicle burns out already descending, the apogee\n');
        fprintf('  event never fires and the propagation is rejected rather than believed.\n');
        else
        fprintf('  The required range is inside the reachable band, but not on the arc that\n');
        fprintf('  was asked for. Ask for the other one, or for ''minimum-energy'', which\n');
        fprintf('  takes whichever branch exists.\n');
        end
        fprintf('=========================================================\n\n');
              traj = [];
       info.refused = true;
      info.branchAsked = branch;
       info.rngReqM = rngReq;
     info.psiLaunch = psiLaunch;
   info.loftStarDeg = rad2deg(loftStarR);
       info.rngMaxM = rngMax;
       info.rngMinM = rngMin;
     info.hApoStarM = hApoStar;
     info.tFlyStarS = tFlyStar;
    info.shortfallM = shortM;
     info.depressed = dep;
        info.lofted = lof;
      info.maxRange = mr;
        if nargout == 0
            clear traj info;
        end
        return;
    end

%% -----------------------------------------------------------------
%% Fly the chosen branch
%% -----------------------------------------------------------------
%% solveBranch already flew this loft angle and kept the trajectory, so nothing
%% is re-integrated here. The check is that the kept trajectory is the one the
%% solver reported, which a struct copied from the wrong branch would fail:
              traj = pick.traj;
    assert(abs(pick.rngM - rI.*coorbital.util.greatCircle(latLaunchR,lonLaunchR, ...
                                                          traj.x(end,3),traj.x(end,2))) < 1e-6, ...
        ['the trajectory carried by the %s branch does not land where that ' ...
         'branch reported; the branch records have been crossed'],pickName);

%% Phase boundaries. phaseRun records the boundary sample ONCE, carrying the
%% outgoing phase's control, so the last row of phase k is that phase's
%% terminal state on the NEAR side of any staging jump:
               kBO = find(traj.phaseIdx == 1,1,'last');
              kApo = find(traj.phaseIdx == 2,1,'last');
                nS = numel(traj.t);
    assert(~isempty(kBO) && ~isempty(kApo), ...
        'the chain did not produce all three phases; phases present: %s', ...
        mat2str(unique(traj.phaseIdx).'));

%% -----------------------------------------------------------------
%% Independent check of the two solves, from the FLOWN state
%% -----------------------------------------------------------------
%% The miss is measured impact-to-target, not read back out of the solver's own
%% bookkeeping. With zero bank the two agree, and the agreement is the evidence
%% for the zero-cross-range claim rather than a restatement of it:
             latIm = traj.x(end,3);
             lonIm = traj.x(end,2);
             missM = rI.*coorbital.util.greatCircle(latIm,lonIm, ...
                                                    latTargetR,lonTargetR);
              resM = pick.rngM - rngReq;

%% Cross-track offset of the impact point from the launch-to-target great
%% circle, positive to the right of the intended course. The standard
%% cross-track formula: the sine of the offset angle is the sine of the flown
%% central angle times the sine of the bearing error at the launch point:
            angFly = coorbital.util.greatCircle(latLaunchR,lonLaunchR,latIm,lonIm);
            psiFly = coorbital.util.greatCircleBearing(latLaunchR,lonLaunchR, ...
                                                       latIm,lonIm);
            dPsiR  = wrapPi(psiFly - psiLaunch);
           xTrackM = rI.*asin(sin(angFly).*sin(dPsiR));
          missHypM = hypot(resM,xTrackM);
        crossWarn  = abs(xTrackM) > tolRngM;

%% The heading the vehicle actually left the pad on must be the bearing the
%% closed form asked for, to machine precision -- x0(6) is written from
%% psiLaunch above and ode45 never touches the initial condition, so any
%% disagreement here means the state was assembled wrong:
           psiSeed = traj.x(1,6);
    assert(abs(wrapPi(psiSeed - psiLaunch)) < 1e-12, ...
        ['the flown initial heading %.15f rad is not the solved launch ' ...
         'azimuth %.15f rad; the launch state was assembled wrong'], ...
        psiSeed,psiLaunch);

%% -----------------------------------------------------------------
%% WHICH BRANCH DID IT ACTUALLY FLY
%% -----------------------------------------------------------------
%% Measured, not assumed. The max-range arc sits BETWEEN the two branches in
%% both apogee and flight time, so either quantity alone identifies the branch
%% and the two together cross-check each other. A bracket lying on one side of
%% the maximum should have held the answer there, but a max-range angle located
%% slightly wrong would break that silently, and the branch is the whole
%% product of this script:
        [flownName,flownAgree] = measureBranch(pick.hApoM,pick.tFlyS, ...
                                               hApoStar,tFlyStar);
         branchOK = strcmp(flownName,pickName) && flownAgree;

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
%% vehicle and the mass actually carried at that instant, so a table-driven aero
%% model or a throttled motor reports correctly here and not just the constant
%% ones:
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

%% Sensed specific force -- what an accelerometer on board would read. Gravity
%% is not sensed and is excluded. nAero drops the thrust, which is the quantity
%% meant by "deceleration" during re-entry:
             nSens = sqrt((aThrV - aDrag).^2 + (aThrN + aLift).^2)./c.g0;
             nAero = sqrt(aLift.^2 + aDrag.^2)./c.g0;

%% Mass must be constant once the motor is out, and it must be the mass the
%% coast vehicle was built around, because the six-state glide equations divide
%% by coastVeh.mass rather than by the mass state. Two separate claims, asserted
%% separately because they can fail separately and have different budgets:
             isUnp = traj.phaseIdx >= 2;
             mSpan = max(massS(isUnp)) - min(massS(isUnp));
    assert(mSpan == 0, ...
        ['the unpowered mass varied by %.3e kg; dm/dt is identically zero ' ...
         'after burnout, so it must be bit-exactly constant'],mSpan);
    assert(abs(massS(end) - mCoast) < 1e-6, ...
        ['the coast is flying %.9f kg while the coast vehicle was built ' ...
         'around %.9f kg; the separation link is wrong'],massS(end),mCoast);

%% Peaks, taken over the phases where each is meaningful:
             isBst = traj.phaseIdx == 1;
             isRe  = traj.phaseIdx == 3;
     [qBstMax,kQB] = maxOver(qbar ,isBst);
     [qReMax ,kQR] = maxOver(qbar ,isRe);
     [nBstMax,kNB] = maxOver(nSens,isBst);
     [nReMax ,kNR] = maxOver(nAero,isRe);

%% Leg ranges, each its own great circle on the impact sphere:
             angBO = coorbital.util.greatCircle(latLaunchR,lonLaunchR, ...
                                                traj.x(kBO,3),traj.x(kBO,2));
          downBOKm = rI.*angBO./1000;

%% Termination diagnosis, one line per phase. An unexpected termination must SAY
%% so: a boost that ran out of horizon before it ran out of propellant, or a
%% coast that never turned over, both produce a short trajectory that otherwise
%% reads as a completed flight:
        [why1,ok1] = whyBurnout(traj.t(kBO),massS(kBO),mBurnoutT,tMaxBoost);
        [why2,ok2] = whyApogee(traj.t(kApo) - traj.t(kBO),traj.t(kApo), ...
                               traj.x(kApo,5),tMaxCoast);
        [why3,ok3] = whyImpact(traj.t(nS) - traj.t(kApo),traj.t(nS), ...
                               traj.x(nS,1) - c.rE,hStopM,hStop,tMaxDesc);
             allOK = ok1 && ok2 && ok3;

%% The textbook minimum-energy burnout flight-path angle for the REQUIRED range
%% angle, gammaStar = 45 deg - Lambda/4. Reported as a diagnostic and used by
%% nothing: see the "what minimum-energy can and cannot mean" note. The 45 and
%% the 4 are constants of the closed-form result in reference [1], not physical
%% constants of the Earth, so they do not belong in missileConst:
          gamStarR = pi./4 - angReqR./4;

%% How far the ground under the target would have moved during this flight had
%% the Earth been turning. Quantifies the non-rotating caveat instead of merely
%% stating it, and it is quantified for BOTH branches because their flight times
%% differ by a factor of nearly three:
          driftDep = c.omegaE.*rI.*cos(latTargetR).*dep.tFlyS;
          driftLof = c.omegaE.*rI.*cos(latTargetR).*lof.tFlyS;

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
%% paragraph below describes the DEFAULT models by name, and every one of those
%% statements becomes a FALSEHOOD the moment the handle is replaced:
      defaultAtmos = isequal(atmosFn,@coorbital.atmos.expAtmos);
       defaultGrav = isequal(gravFn,@coorbital.grav.sphereGrav);
       defaultProp = isequal(propFn,@coorbital.prop.constThrust);

%% Report:
    fprintf('\n');
    fprintf('===== Ballistic point-to-point targeting summary ========\n');
    fprintf('  Phase termination\n');
    fprintf('    1 boost          %s\n',why1);
    fprintf('    2 coast          %s\n',why2);
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
            pick.rngM./1000,'km');
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
    fprintf('    cross-track      %+10.2f %-5s  (offset of the impact point from the\n', ...
            xTrackM,'m');
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
        fprintf('                   the arc it left on. Set bankAngle to zero, or accept the\n');
        fprintf('                   cross-range: closing it needs an outer azimuth iteration,\n');
        fprintf('                   which this script does not have.\n');
    end
    fprintf('\n');
    fprintf('  The range-versus-loft hump, bracketed before either branch was solved\n');
    fprintf('    max-range loft   %10.4f %-5s  (commanded terminal pitch attitude)\n', ...
            rad2deg(loftStarR),'deg');
    fprintf('    MAXIMUM RANGE    %10.3f %-5s  (nothing this vehicle can do exceeds it)\n', ...
            rngMax./1000,'km');
    fprintf('    its apogee       %10.3f %-5s  and flight time %.3f s. These two are the\n', ...
            hApoStar./1000,'km',tFlyStar);
    fprintf('                                        REFERENCE the flown branch is measured\n');
    fprintf('                                        against; the max-range arc sits between\n');
    fprintf('                                        the lofted and depressed ones in both.\n');
    fprintf('    loft bracket     %10.4f to %.4f deg   (searched; the maximum is INTERIOR,\n', ...
            loftMin,loftMax);
    fprintf('                                        which is what gives two branches at all)\n');
    fprintf('    reachable        %10.3f to %.3f km   (union of the two branch bands)\n', ...
            rngMin./1000,rngMax./1000);
    fprintf('    search cost      %10d %-5s  (%d coarse scan points, %d golden-section\n', ...
            mr.nEval,'',mr.nScan,mr.nGolden);
    fprintf('                                        steps, closed to %.4f deg of loft)\n', ...
            rad2deg(mr.widthR));
    if mr.nTurn > 1
        fprintf('  *** CAUTION ***  the coarse scan changed direction %d times, so the range\n', ...
                mr.nTurn);
        fprintf('                   curve is NOT the single hump this method assumes. The\n');
        fprintf('                   maximum reported above may be a local one. Raise\n');
        fprintf('                   nScanLoft and look at info.maxRange.scanRngM.\n');
    end
    fprintf('\n');
    fprintf('  THE TWO ARCS THAT REACH THIS TARGET\n');
    fprintf('    %-22s %16s %16s\n','','depressed','lofted');
    fprintf('    %-22s %16s %16s\n','reaches the target', ...
            yesNo(dep.exists),yesNo(lof.exists));
    printPair('loft angle (deg)'      ,rad2deg(dep.loftR)   ,rad2deg(lof.loftR)   ,dep,lof,'%16.4f');
    printPair('burnout gamma (deg)'   ,rad2deg(dep.gamBoR)  ,rad2deg(lof.gamBoR)  ,dep,lof,'%16.4f');
    printPair('achieved range (km)'   ,dep.rngM./1000       ,lof.rngM./1000       ,dep,lof,'%16.4f');
    printPair('miss (m)'              ,dep.missM            ,lof.missM            ,dep,lof,'%16.4f');
    printPair('APOGEE (km)'           ,dep.hApoM./1000      ,lof.hApoM./1000      ,dep,lof,'%16.4f');
    printPair('FLIGHT TIME (s)'       ,dep.tFlyS            ,lof.tFlyS            ,dep,lof,'%16.4f');
    printPair('IMPACT SPEED (m/s)'    ,dep.vImpM            ,lof.vImpM            ,dep,lof,'%16.4f');
    printPair('IMPACT ANGLE (deg)'    ,rad2deg(dep.gamImR)  ,rad2deg(lof.gamImR)  ,dep,lof,'%16.4f');
    printPair('burnout energy (MJ/kg)',dep.epsBo./1e6       ,lof.epsBo./1e6       ,dep,lof,'%16.4f');
    printPair('loft from max-range'   ,rad2deg(dLoftDep)    ,rad2deg(dLoftLof)    ,dep,lof,'%16.4f');
    printPair('bisection steps'       ,dep.iterations       ,lof.iterations       ,dep,lof,'%16d');
    printPair('propagations'          ,dep.nProp            ,lof.nProp            ,dep,lof,'%16d');
    if bothHere
    fprintf('    The lofted arc arrives STEEPER: %.3f deg against %.3f deg. It does NOT\n', ...
            rad2deg(lof.gamImR),rad2deg(dep.gamImR));
    if lof.vImpM > dep.vImpM
    fprintf('    arrive slower -- it arrives at %.1f m/s against the depressed arc''s %.1f,\n', ...
            lof.vImpM,dep.vImpM);
    fprintf('    which is the OPPOSITE of the vacuum intuition and is a DRAG result. In\n');
    fprintf('    vacuum the two arcs carry almost the same burnout energy (they differ by\n');
    fprintf('    %.2f %% here) and so would arrive at almost the same speed; with an\n', ...
            100.*abs(lof.epsBo - dep.epsBo)./abs(dep.epsBo));
    fprintf('    atmosphere the shallow %.1f deg descent of the depressed arc spends far\n', ...
            rad2deg(dep.gamImR));
    fprintf('    longer in dense air than the steep %.1f deg plunge of the lofted one, and\n', ...
            rad2deg(lof.gamImR));
    fprintf('    is braked to %.1f %% of the lofted arrival speed.\n', ...
            100.*dep.vImpM./lof.vImpM);
    else
    fprintf('    and arrives SLOWER: %.1f m/s against %.1f m/s, which is the vacuum\n', ...
            lof.vImpM,dep.vImpM);
    fprintf('    intuition. Note that it need not hold with an atmosphere -- a shallow\n');
    fprintf('    depressed descent is braked far harder than a steep lofted plunge.\n');
    end
    fprintf('    Only ONE of these two was flown. The other is reported so the trade is\n');
    fprintf('    visible: %.0f s and %.0f km of apogee separate them on this geometry.\n', ...
            abs(lof.tFlyS - dep.tFlyS),abs(lof.hApoM - dep.hApoM)./1000);
    else
    fprintf('    Only one arc reaches this target, so there is no trade to make. The other\n');
    fprintf('    branch''s band is printed above as [%.3f, %.3f] km against a required\n', ...
            rngMin./1000,rngMax./1000);
    fprintf('    %.3f km.\n',rngReq./1000);
    end
    fprintf('\n');
    fprintf('  BRANCH SELECTION\n');
    fprintf('    asked for        %-16s (%s)\n',branch,pickShort);
    fprintf('    on the grounds that it is %s\n',pickWhy);
    fprintf('    flew             %-16s (the branch this run actually solved for)\n',pickName);
    fprintf('    MEASURED as      %-16s (from the FLOWN apogee %.3f km against the\n', ...
            flownName,pick.hApoM./1000);
    fprintf('                                      max-range arc''s %.3f km, and the flown\n', ...
            hApoStar./1000);
    fprintf('                                      %.3f s flight time against its %.3f s)\n', ...
            pick.tFlyS,tFlyStar);
    if ~flownAgree
        fprintf('  *** CAUTION ***  the apogee and the flight time DISAGREE about which branch\n');
        fprintf('                   this is. One of them puts the flight above the max-range\n');
        fprintf('                   arc and the other below it, which a genuine lofted or\n');
        fprintf('                   depressed arc cannot do. Suspect the max-range angle.\n');
    end
    if ~branchOK
        fprintf('  *** CAUTION ***  the flown trajectory MEASURES as the %s branch while the\n', ...
                flownName);
        fprintf('                   solve was run on the %s bracket. The bracket did not hold\n', ...
                pickName);
        fprintf('                   the solution on its own side of the maximum, so the\n');
        fprintf('                   max-range angle above is suspect. Do not read the branch\n');
        fprintf('                   label; read the apogee and the flight time.\n');
    end
    if strcmp(branch,'minimum-energy')
    fprintf('    WHAT MINIMUM-ENERGY MEANS HERE, precisely, because it does not mean what\n');
    fprintf('    the textbook means. The textbook trajectory reaches the required range on\n');
    fprintf('    the LEAST burnout energy, at a burnout flight-path angle of\n');
    fprintf('    gamma* = 45 - Lambda/4 = %.4f deg for this %.4f deg range angle. That\n', ...
            rad2deg(gamStarR),rad2deg(angReqR));
    fprintf('    presumes the energy is free to choose. It is not: this booster carries a\n');
    fprintf('    fixed propellant load and burns to exhaustion on BOTH branches, so both\n');
    fprintf('    leave burnout with essentially the same energy -- %.4f and %.4f MJ/kg,\n', ...
            dep.epsBo./1e6,lof.epsBo./1e6);
    fprintf('    %.2f %% apart, and that difference is boost drag and gravity loss rather\n', ...
            100.*abs(lof.epsBo - dep.epsBo)./abs(dep.epsBo));
    fprintf('    than a design variable.\n');
    fprintf('    SO WHAT IS IMPLEMENTED IS: the branch whose loft angle lies NEARER the\n');
    fprintf('    %.4f deg max-range angle, on the grounds that the max-range arc IS the\n', ...
            rad2deg(loftStarR));
    fprintf('    minimum-energy arc for its own range. Here that is %.4f deg against\n', ...
            rad2deg(min(dLoftDep,dLoftLof)));
    fprintf('    %.4f deg, so the %s arc wins.\n', ...
            rad2deg(max(dLoftDep,dLoftLof)),pickName);
    fprintf('    TWO CHECKS ON THAT CHOICE, both reported and neither used:\n');
    fprintf('      least burnout energy read literally would pick the %s arc, the\n', ...
            lowerOf('depressed','lofted',dep.epsBo,lof.epsBo));
    fprintf('      opposite reading -- lower burnout energy with the propellant fixed means\n');
    fprintf('      a COSTLIER ascent, not a cheaper trajectory, so it is not used;\n');
    fprintf('      nearest gamma* would pick the %s arc (%.3f deg from %.3f deg against\n', ...
            nearerOf('depressed','lofted',abs(dep.gamBoR - gamStarR), ...
                     abs(lof.gamBoR - gamStarR)), ...
            rad2deg(min(abs(dep.gamBoR - gamStarR),abs(lof.gamBoR - gamStarR))), ...
            rad2deg(gamStarR));
    fprintf('      %.3f deg), but it is measured on a burnout state the alphaMax clamp and\n', ...
            rad2deg(max(abs(dep.gamBoR - gamStarR),abs(lof.gamBoR - gamStarR))));
    fprintf('      not this script controls, so it is a diagnostic.\n');
    if meClose
    fprintf('  *** CAUTION ***  the two loft distances agree to %.2f %%, inside the 5 %% band\n', ...
            100.*(1 - meRatio));
    fprintf('                   where this selection turns on a rounding rather than on the\n');
    fprintf('                   physics. The range-versus-loft hump is roughly symmetric and\n');
    fprintf('                   this happens over a wide band of targets. ASK FOR ''lofted''\n');
    fprintf('                   OR ''depressed'' EXPLICITLY at this range.\n');
    end
    end
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
    fprintf('  Burnout, end of phase 1\n');
    fprintf('    time             %10.3f %-5s  (propellant exhaustion; the boost is never cut\n', ...
            traj.t(kBO),'s');
    fprintf('                                        short here, the LOFT ANGLE being the\n');
    fprintf('                                        ranging control instead of the cutoff)\n');
    fprintf('    altitude         %10.3f km\n',hKm(kBO));
    fprintf('    speed            %10.2f %-5s  (Mach %.2f, planet-relative)\n', ...
            V(kBO),'m/s',mach(kBO));
    fprintf('    flight path      %10.3f %-5s  (positive is climbing; the COMMANDED terminal\n', ...
            rad2deg(traj.x(kBO,5)),'deg');
    fprintf('                                        attitude was %.4f deg. The %.1f deg clamp\n', ...
            rad2deg(pick.loftR),alphaMax);
    fprintf('                                        on angle of attack is the difference)\n');
    fprintf('    heading          %10.3f deg\n',rad2deg(traj.x(kBO,6)));
    fprintf('    downrange        %10.2f %-5s  (great circle from the pad)\n',downBOKm,'km');
    fprintf('    mass at burnout  %10.1f %-5s  (payload + spent booster)\n',massS(kBO),'kg');
    fprintf('    mass into coast  %10.1f %-5s  (%s)\n',mCoast,'kg',sepTxt);
    fprintf('\n');
    fprintf('  Apogee, end of phase 2\n');
    fprintf('    time             %10.3f s\n',traj.t(kApo));
    fprintf('    altitude         %10.3f %-5s  (exact event, not a grid maximum)\n', ...
            hKm(kApo),'km');
    fprintf('    speed            %10.2f m/s\n',V(kApo));
    fprintf('    flight path      %10.3e %-5s  (zero by construction at apogee)\n', ...
            rad2deg(traj.x(kApo,5)),'deg');
    fprintf('\n');
    fprintf('  Overall\n');
    fprintf('    flight time      %10.2f %-5s  (%.2f min)\n',traj.t(end),'s',traj.t(end)./60);
    fprintf('    ground range     %10.2f %-5s  (great circle on the r = %.1f km impact sphere)\n', ...
            pick.rngM./1000,'km',rI./1000);
    fprintf('    central angle    %10.4f deg\n',rad2deg(angFly));
    fprintf('    apogee ratio     %10.4f %-5s  (apogee over range; a minimum-energy arc sits\n', ...
            (traj.x(kApo,1) - c.rE)./pick.rngM,'');
    fprintf('                                        near 0.25, so this arc is %s)\n', ...
            loftedOrNot((traj.x(kApo,1) - c.rE)./pick.rngM));
    fprintf('    samples          %10d %-5s  (ode45 adaptive steps over 3 phases)\n',nS,'');
    fprintf('    propagations     %10d %-5s  (whole trajectories flown to produce this run:\n', ...
            mr.nEval + 1 + dep.nProp + lof.nProp,'');
    fprintf('                                        %d bracketing the maximum, 1 re-flying it,\n', ...
            mr.nEval);
    fprintf('                                        %d on the depressed branch, %d on the lofted)\n', ...
            dep.nProp,lof.nProp);
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
    fprintf('    re-entry max q   %10.2f %-5s  (at t = %.1f s, h = %.2f km)\n', ...
            qReMax./1000,'kPa',traj.t(kQR),hKm(kQR));
    fprintf('    re-entry decel   %10.2f %-5s  (%.2f m/s^2 aero only, at t = %.1f s, h = %.2f km)\n', ...
            nReMax,'g',nReMax.*c.g0,traj.t(kNR),hKm(kNR));
    fprintf('\n');
    fprintf('  LIMITATIONS OF THIS TARGETING SOLUTION\n');
    fprintf('    1. THE AZIMUTH IS EXACT ONLY FOR A NON-ROTATING EARTH. This run has Earth\n');
    fprintf('       rotation %s, env.omegaE = %.6e rad/s, so the ground track of an\n', ...
            spinTxt,env.omegaE);
    fprintf('       unbanked trajectory is a great circle and the initial bearing of the\n');
    fprintf('       launch-to-target arc is the WHOLE answer -- one closed-form call, no\n');
    fprintf('       iteration. Turn the Earth on and it stops being the answer, and the two\n');
    fprintf('       branches show how much it matters: the ground under the target sweeps\n');
    fprintf('       %.0f km east during the depressed arc''s %.0f s and %.0f km during the\n', ...
            driftDep./1000,dep.tFlyS,driftLof./1000);
    fprintf('       lofted arc''s %.0f s. The azimuth would then depend on WHICH BRANCH IS\n', ...
            lof.tFlyS);
    fprintf('       FLOWN, so the outer iteration this script does not have would have to\n');
    fprintf('       run per branch. Rotating-Earth targeting is out of scope for it.\n');
    if earthSpin
    fprintf('       *** earthSpin IS ON. The azimuth above was computed as though it were\n');
    fprintf('       not, so the miss reported above is against a target that stood still and\n');
    fprintf('       is NOT this vehicle''s miss. Do not read it as one. ***\n');
    end
    fprintf('    2. THE MISS IS THE RESIDUAL OF THE RANGE SOLVE, ALONG THE GREAT CIRCLE.\n');
    fprintf('       Bisection matches a DISTANCE; nothing in it steers sideways.\n');
    if bankRad ~= 0
    fprintf('       THIS RUN COMMANDS A NON-ZERO BANK of %.1f deg, so the ground track is NOT\n', ...
            bankAngle);
    fprintf('       the launch-to-target great circle and the range solve cannot see the\n');
    fprintf('       difference. Measured on this run: %.2f m of cross-track offset against a\n', ...
            abs(xTrackM));
    fprintf('       %.2f m range residual, for a %.2f m total miss.\n',abs(resM),missM);
    else
    fprintf('       Cross-range came out at %.2e m here because the bank angle is zero\n', ...
            abs(xTrackM));
    fprintf('       throughout, and an unbanked trajectory over a non-rotating sphere never\n');
    fprintf('       leaves the great circle it departed on. THAT IS A PROPERTY OF THIS\n');
    fprintf('       CONFIGURATION, NOT A GENERAL GUARANTEE. It is measured above, not assumed.\n');
    end
    fprintf('    3. THE LOFT ANGLE IS A COMMANDED ATTITUDE, NOT THE BURNOUT GAMMA. The %.1f deg\n', ...
            alphaMax);
    fprintf('       clamp on angle of attack limits how fast the flight path can be pushed\n');
    fprintf('       over, so the %.4f deg commanded here produced a %.4f deg burnout flight\n', ...
            rad2deg(pick.loftR),rad2deg(pick.gamBoR));
    fprintf('       path. At BM/run_ballistic''s 6 deg clamp the depressed branch does not\n');
    fprintf('       exist at all: the clamp holds burnout gamma above the max-range value\n');
    fprintf('       over the whole bracket and range becomes monotone in loft.\n');
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
    fprintf('            re-entry lie outside it.\n');
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
        fprintf('            caveat this script prints by default does not apply.\n');
    end
    if ~defaultProp
        fprintf('            The motor is %s, NOT the default constant-thrust\n', ...
                func2str(propFn));
        fprintf('            model, so the %.4f s burn quoted above is nominal only.\n',tBurn);
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
        fprintf('            The %.1f g sensed load at end of burn is an artefact of the\n',nBstMax);
        fprintf('            constant-thrust placeholder motor, which does not throttle as the\n');
        fprintf('            stack empties. A real stage would tail off or stage before that.\n');
    end
    fprintf('=========================================================\n\n');

%% Hand the summary's numbers back at full precision, so a test does not have to
%% read them out of printed text:
       info.refused = false;
   info.branchAsked = branch;
   info.branchFlown = pickName;
info.branchMeasured = flownName;
  info.branchAgrees = branchOK;
   info.branchWhy   = pickWhy;
info.branchTimeAgrees = flownAgree;
       info.rngReqM = rngReq;
       info.rngAchM = pick.rngM;
         info.missM = missM;
        info.residM = resM;
       info.xTrackM = xTrackM;
      info.missHypM = missHypM;
     info.crossWarn = crossWarn;
     info.psiLaunch = psiLaunch;
        info.psiFly = psiFly;
         info.dPsiR = dPsiR;
   info.loftStarDeg = rad2deg(loftStarR);
       info.rngMaxM = rngMax;
       info.rngMinM = rngMin;
     info.hApoStarM = hApoStar;
     info.tFlyStarS = tFlyStar;
      info.maxRange = mr;
     info.depressed = dep;
        info.lofted = lof;
    info.dLoftDepDeg = rad2deg(dLoftDep);
    info.dLoftLofDeg = rad2deg(dLoftLof);
       info.meRatio = meRatio;
       info.meClose = meClose;
     info.gamStarR  = gamStarR;
      info.loftDeg  = rad2deg(pick.loftR);
        info.tBurn  = tBurn;
     info.tBurnout  = traj.t(kBO);
     info.hBurnout  = traj.x(kBO,1) - c.rE;
     info.vBurnout  = V(kBO);
   info.gamBurnout  = traj.x(kBO,5);
     info.mBurnout  = massS(kBO);
       info.mCoast  = mCoast;
     info.mLiftoff  = mLiftoff;
     info.downBOKm  = downBOKm;
      info.tApogee  = traj.t(kApo);
      info.hApogee  = traj.x(kApo,1) - c.rE;
      info.vApogee  = V(kApo);
      info.tFlight  = traj.t(end);
      info.rangeKm  = pick.rngM./1000;
       info.angTot  = angFly;
      info.vImpact  = V(end);
    info.gamImpact  = traj.x(end,5);
    info.latImpact  = latIm;
    info.lonImpact  = lonIm;
      info.machImp  = mach(end);
      info.qBstMax  = qBstMax;
      info.qReMax   = qReMax;
      info.nBstMax  = nBstMax;
      info.nReMax   = nReMax;
       info.driftDepM = driftDep;
       info.driftLofM = driftLof;
        info.stopOK = allOK;
       info.stopWhy = {why1;why2;why3};
        info.nProp  = mr.nEval + 1 + dep.nEval + lof.nEval;

%% Vertical exaggeration for the globe, scaled to the apogee the run actually
%% flew rather than fixed, because the two branches differ by a factor of ten in
%% apogee and no single factor suits both. Same rule and same cap and floor as
%% HGV/run_target, whose exagFor header carries the full rationale:
            hPeakM = max(traj.x(:,1)) - c.rE;
           altExag = overrideOf(opts,'altExag',exagFor(hPeakM,c.rE));
       info.altExag = altExag;
   info.altExagRule = @(hM) exagFor(hM,c.rE);

%% Plots. Every figure comes from coorbital.viz, which reads the trajectory and
%% never writes it. Nothing below this line can move a number in the summary
%% above. ONE VEHICLE PER PHASE is handed to the profile plot, because this
%% chain does not fly one; the ground track and the globe are given the TARGET,
%% so the aim point and the impact point are both on the picture:
     info.launchStr = llStr(latLaunch,lonLaunch);
     info.targetStr = llStr(latTarget,lonTarget);
    if showPlots
        coorbital.viz.profilePlot(traj,bst,env, ...
            struct('Name','Ballistic targeting profile', ...
                   'Channels',{{'altitude','speed','mach','q','nAero','gamma'}}, ...
                   'Extra',{{nSens,'sensed load factor (g)', ...
                             'Sensed load factor, thrust included'}}, ...
                   'VehPhase',{{bst,coastVeh,coastVeh}}));
        coorbital.viz.groundTrack(traj,bst,env, ...
            struct('Target',[latTargetR; lonTargetR], ...
                   'PhaseName',{{'boost','coast','descent'}}, ...
                   'Title',sprintf(['Ground track, %s arc, %.0f km required, ' ...
                                    '%.0f m miss (loft %.2f deg)'], ...
                                   pickName,rngReq./1000,missM,rad2deg(pick.loftR))));
        coorbital.viz.globe3D(traj,bst,env, ...
            struct('Target',[latTargetR; lonTargetR], ...
                   'AltScale',altExag, ...
                   'Title',sprintf('Solved %s arc, %s to %s', ...
                                   pickName,info.launchStr,info.targetStr)));
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
                   'PhaseName',{{'boost','coast','descent'}}, ...
                   'Title',sprintf('%s arc, launch %s to target %s', ...
                                   pickName,info.launchStr,info.targetStr)));
        fprintf('  Movie: %d frames at %g fps written to\n    %s\n', ...
                mv.nFrame,mv.frameRate,mv.file);
        fprintf('         Earth texture %s, background %s\n\n',mv.texture,mv.background);
         info.movie  = mv;
    end

%% Hand the trajectory back only when the caller asked for it. Typing
%% "run_ballistic_target" at the prompt should leave the summary on screen, not
%% bury it under a dump of the whole struct as ans:
    if nargout == 0
        clear traj info;
    end
end

function ph = buildPhases(loftR,cfg)
%% Purpose:
%
%  Assemble the three-phase boost-coast-descent chain for ONE loft angle.
%  Factored out because the max-range search, both branch solves and the final
%  re-fly all need the same chain, and two copies of it would be two chains
%  that could drift apart.
%
%  THE LOFT ANGLE ENTERS THROUGH THE PITCH SCHEDULE AND NOWHERE ELSE. The
%  schedule is the reference profile's SHAPE rescaled so that its first node
%  stays at cfg.thetaTop and its last node becomes loftR; nothing else about the
%  chain depends on it. The boost still ends at propellant exhaustion, so a
%  change of loft angle changes WHERE the vehicle is going when the motor quits,
%  never how much propellant it burned.
%
%  EACH PHASE CARRIES ITS OWN VEHICLE on ph.veh, which
%  coorbital.prop.phaseRun forwards to the equations of motion. That is why the
%  eom fields below are the library handles themselves rather than closures
%  binding a vehicle and ignoring the forwarded argument.
%
%% Inputs:
%
%  loftR            [1 x 1]                     Commanded terminal pitch
%                                               attitude (rad)
%
%  cfg              Struct                      Immutable configuration built
%                                               by the caller: eomCoast
%                                               (handle), schedCst (struct),
%                                               pitchT [1 x K] (s), shape
%                                               [1 x K] (-), thetaTop (rad),
%                                               alphaMax (rad), bankRad (rad),
%                                               mBurnout (kg), mCoast (kg),
%                                               separation (logical), massDry
%                                               (kg), hStop (m), tMaxBoost,
%                                               tMaxCoast, tMaxDesc (s), bst and
%                                               coastVeh (structs)
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

            thetaR = cfg.thetaTop - cfg.shape.*(cfg.thetaTop - loftR);
            schedB = struct('tGrid',cfg.pitchT, ...
                            'theta',thetaR, ...
                            'sigma',cfg.bankRad.*ones(1,numel(cfg.pitchT)), ...
                            'alphaMax',cfg.alphaMax);

%% Phase 1, boost: burns to propellant exhaustion, then jettisons the spent
%% booster. With separation switched off the link is [], which phaseRun reads
%% as identity:
         ph(1).eom = @coorbital.eom.boost3DOF;
       ph(1).guide = @(t,x) coorbital.guide.pitchProgram(t,x,schedB);
   ph(1).terminate = @(t,x) coorbital.prop.eventBurnout(t,x,cfg.mBurnout);
       ph(1).tspan = [0 cfg.tMaxBoost];
        ph(1).link = [];
         ph(1).veh = cfg.bst;
    if cfg.separation
        ph(1).link = @(x) [x(1:6); x(7) - cfg.massDry];
    end

%% Phase 2, ascending coast: unpowered, ends exactly at apogee. Splitting the
%% coast there makes apogee an exactly solved event rather than the largest
%% sample on an adaptive grid, and leaves phase 3 monotonically descending so
%% the one-sided impact event cannot be reached from the wrong side:
         ph(2).eom = cfg.eomCoast;
       ph(2).guide = @(t,x) coorbital.guide.prescribed(t,x,cfg.schedCst);
   ph(2).terminate = @(t,x) coorbital.prop.eventApogee(t,x);
       ph(2).tspan = [0 cfg.tMaxCoast];
        ph(2).link = [];
         ph(2).veh = cfg.coastVeh;

%% Phase 3, descent to impact: unpowered, ends at the impact altitude:
         ph(3).eom = cfg.eomCoast;
       ph(3).guide = @(t,x) coorbital.guide.prescribed(t,x,cfg.schedCst);
   ph(3).terminate = @(t,x) coorbital.prop.eventAltitude(t,x,cfg.hStop);
       ph(3).tspan = [0 cfg.tMaxDesc];
        ph(3).link = [];
         ph(3).veh = cfg.coastVeh;
end

function [rngM,traj] = flyLoft(loftR,cfg,rI)
%% Purpose:
%
%  Fly the whole boost-coast-descent chain at one loft angle and return the
%  great-circle surface range from the launch point to the impact point. This
%  is the function the max-range search maximises and the function
%  coorbital.util.rangeSolve bisects on, and it is the only place either of them
%  touches the physics.
%
%  A PROPAGATION THAT DID NOT COMPLETE IS AN ERROR, NOT A SHORT RANGE. If a
%  phase runs out of horizon -- the usual cause being a loft angle so depressed
%  that the vehicle burns out already descending, so the apogee event never
%  fires -- the trajectory still ends SOMEWHERE and still has a great-circle
%  range. Returning it would feed the search a number from a flight that never
%  happened and quietly break both the unimodality the bracketing rests on and
%  the monotonicity the bisection rests on. It throws instead, naming the
%  bracket entry to change.
%
%% Inputs:
%
%  loftR            [1 x 1]                     Commanded terminal pitch
%                                               attitude (rad)
%
%  cfg              Struct                      Immutable configuration; see
%                                               buildPhases, plus x0 [7 x 1],
%                                               env (struct), rE (m), lat0
%                                               (rad), lon0 (rad)
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

                ph = buildPhases(loftR,cfg);
              traj = coorbital.prop.phaseRun(ph,cfg.x0,cfg.bst,cfg.env);

%% All three phases must have run, and the last must have stopped ON the impact
%% altitude rather than at its horizon:
             hEndM = traj.x(end,1) - cfg.rE;
             nPhRn = numel(unique(traj.phaseIdx));
    if nPhRn < 3 || abs(hEndM - cfg.hStop) > 1e-3
        error('coorbital:runBallisticTarget:propagationIncomplete', ...
            ['A loft angle of %.6f deg produced %d of 3 phases and ended at ' ...
             'h = %.3f km against a %.3f km stop altitude. The flight did not ' ...
             'complete, so its range is meaningless to the search. The usual ' ...
             'cause is a loft angle so depressed that the vehicle burns out ' ...
             'already descending and the apogee event never fires: raise ' ...
             'loftMin.'],rad2deg(loftR),nPhRn,hEndM./1000,cfg.hStop./1000);
    end

              rngM = rI.*coorbital.util.greatCircle(cfg.lat0,cfg.lon0, ...
                                                    traj.x(end,3),traj.x(end,2));
end

function [loftStarR,rngStar,mr] = maxRangeLoft(fRange,loftLoR,loftHiR,nScan,tolLoftR)
%% Purpose:
%
%  Locate the loft angle of MAXIMUM range on the bracket, which is what splits
%  the loft axis into a lofted and a depressed branch. Two stages, and the first
%  is what makes the second safe:
%
%    A COARSE SCAN of nScan equally spaced loft angles. Golden-section search
%    requires a UNIMODAL bracket and cannot detect that it was not given one, so
%    the scan both finds a three-point bracket around the best node and provides
%    the evidence: the number of direction changes in the scanned range curve is
%    counted and returned, and one is what a single hump looks like.
%
%    A GOLDEN-SECTION SEARCH on that three-point bracket, which needs no
%    derivative -- there is none available, each evaluation being a trajectory
%    propagation -- and shrinks the bracket by a fixed factor per evaluation
%    while reusing one of the two interior points every step.
%
%  THE MAXIMUM MUST BE INTERIOR. A coarse-scan maximum sitting on either
%  endpoint means the hump is outside the bracket, and then there is no
%  two-branch structure to solve at all: range is monotonic across the whole
%  bracket and one of the two branches is empty everywhere. That is not a
%  numerical near-miss to be nudged, it is a configuration error, and the
%  message below names its usual cause -- an alphaMax clamp too tight for the
%  vehicle to be pushed past the max-range attitude.
%
%% Inputs:
%
%  fRange           Function handle             rng = fRange(loftR), metres of
%                                               surface range for a commanded
%                                               terminal attitude in radians
%
%  loftLoR          [1 x 1]                     Bottom of the loft bracket (rad)
%
%  loftHiR          [1 x 1]                     Top of the loft bracket (rad),
%                                               greater than loftLoR
%
%  nScan            [1 x 1]                     Coarse scan points, at least 5
%
%  tolLoftR         [1 x 1]                     Bracket width to close the
%                                               golden section to (rad)
%
%% Outputs:
%
%  loftStarR        [1 x 1]                     Max-range loft angle (rad)
%
%  rngStar          [1 x 1]                     Range achieved there (m)
%
%  mr               Struct                      Search record:
%                                               scanLoftR [1 x nScan] (rad)
%                                               scanRngM  [1 x nScan] (m)
%                                               nScan     [1 x 1]
%                                               nGolden   [1 x 1] golden steps
%                                               nEval     [1 x 1] total calls
%                                                         to fRange
%                                               widthR    [1 x 1] final bracket
%                                                         width (rad)
%                                               nTurn     [1 x 1] direction
%                                                         changes in the scan;
%                                                         1 is a single hump
%
%% References:
%   [1] Press, W.H., et al., "Numerical Recipes," 3rd ed., Cambridge, 2007,
%       Section 10.2. Golden-section search in one dimension.
%
%% Revision History:
%  Michael Casey                                                08/07/2026
%  Copyright 2026 Coorbital, Inc.
%% ------------------------ Begin Code Sequence ---------------------------

%% The coarse scan:
          scanLoftR = linspace(loftLoR,loftHiR,nScan);
           scanRngM = zeros(1,nScan);
    for ks = 1:nScan
       scanRngM(ks) = fRange(scanLoftR(ks));
    end
             [~,kB] = max(scanRngM);

%% Direction changes in the scanned curve. A single hump changes direction
%% exactly once; more than that and the golden section below is being run on a
%% bracket it cannot be trusted on, which the caller reports rather than
%% suppresses:
              dSign = sign(diff(scanRngM));
              dSign = dSign(dSign ~= 0);
              nTurn = sum(diff(dSign) ~= 0);

%% The maximum must be interior, or there are not two branches:
    if kB == 1 || kB == nScan
        error('coorbital:runBallisticTarget:maximumNotBracketed', ...
            ['The largest range over the loft bracket %.4f to %.4f deg was ' ...
             'found AT the %.4f deg endpoint, so the max-range loft angle is ' ...
             'outside the bracket and range is monotonic across all of it. ' ...
             'There is then only ONE branch, and the lofted/depressed structure ' ...
             'this script is built around does not exist. Two usual causes, in ' ...
             'order of likelihood: (1) the alphaMax clamp is too tight for the ' ...
             'vehicle to be pushed past the max-range attitude -- at 6 deg on ' ...
             'the shipped booster a commanded -30 deg still burns out at ' ...
             '33.4 deg of flight-path angle, above the max-range value, and ' ...
             'the depressed branch is unreachable at every commanded loft; ' ...
             '(2) the bracket itself is too narrow. Widen loftMin and ' ...
             'loftMax, or raise alphaMax.'], ...
            rad2deg(loftLoR),rad2deg(loftHiR),rad2deg(scanLoftR(kB)));
    end

%% Golden-section search on the three-point bracket around the best scan node.
%% gRat is the reciprocal golden ratio; the two interior points are placed so
%% that whichever half survives already contains one of them:
                 aL = scanLoftR(kB-1);
                 bL = scanLoftR(kB+1);
               gRat = (sqrt(5) - 1)./2;
                 x1 = bL - gRat.*(bL - aL);
                 x2 = aL + gRat.*(bL - aL);
                 f1 = fRange(x1);
                 f2 = fRange(x2);
            nGolden = 0;
            maxStep = 200;
    while (bL - aL) > tolLoftR && nGolden < maxStep
            nGolden = nGolden + 1;
        if f1 < f2
                 aL = x1;
                 x1 = x2;
                 f1 = f2;
                 x2 = aL + gRat.*(bL - aL);
                 f2 = fRange(x2);
        else
                 bL = x2;
                 x2 = x1;
                 f2 = f1;
                 x1 = bL - gRat.*(bL - aL);
                 f1 = fRange(x1);
        end
    end
          loftStarR = 0.5.*(aL + bL);
            rngStar = fRange(loftStarR);

%% The located maximum must beat both bracket endpoints, or the search has
%% wandered off a curve that is not the hump the scan showed:
    if rngStar < max(scanRngM(1),scanRngM(nScan))
        error('coorbital:runBallisticTarget:maximumNotFound', ...
            ['The golden section settled at %.6f deg with %.4f km of range, ' ...
             'which is less than the %.4f km reached at a bracket endpoint. ' ...
             'The range curve is not the single hump the coarse scan showed.'], ...
            rad2deg(loftStarR),rngStar./1000, ...
            max(scanRngM(1),scanRngM(nScan))./1000);
    end

      mr.scanLoftR  = scanLoftR;
       mr.scanRngM  = scanRngM;
          mr.nScan  = nScan;
        mr.nGolden  = nGolden;
         mr.nEval   = nScan + 2 + nGolden + 1;
         mr.widthR  = bL - aL;
          mr.nTurn  = nTurn;
end

function br = solveBranch(name,rngReq,fRange,loftAR,loftBR,tolRngM,cfg,rI,c, ...
                          latTgR,lonTgR,hApoStar,tFlyStar)
%% Purpose:
%
%  Solve ONE branch of the two-branch loft problem and gather everything the
%  summary reports about it. The bracket handed in lies entirely on one side of
%  the max-range loft angle, which is what makes range monotonic on it and
%  bisection safe -- see the MONOTONICITY IS ASSUMED, NOT CHECKED note in
%  coorbital.util.rangeSolve.
%
%  A BRANCH THAT CANNOT REACH THE TARGET IS NOT AN ERROR. rangeSolve returns
%  converged = false with the achievable band rather than throwing, and that is
%  reported here as exists = false with the band intact, so the caller can
%  refuse in the problem's own vocabulary or fall back to the other branch.
%
%  The trajectory of the converged solution is KEPT rather than re-flown by the
%  caller: rangeSolve returns only the scalar range, so the state history is
%  re-created once here and handed back, which is one propagation rather than
%  two.
%
%% Inputs:
%
%  name             Char [1 x n]                'lofted' or 'depressed', for
%                                               the messages
%
%  rngReq           [1 x 1]                     Required surface range (m)
%
%  fRange           Function handle             rng = fRange(loftR)
%
%  loftAR           [1 x 1]                     One end of this branch's loft
%                                               bracket (rad)
%
%  loftBR           [1 x 1]                     The other end (rad), greater
%                                               than loftAR
%
%  tolRngM          [1 x 1]                     Convergence tolerance on the
%                                               achieved range (m)
%
%  cfg              Struct                      Chain configuration; see
%                                               buildPhases
%
%  rI               [1 x 1]                     Impact sphere radius (m)
%
%  c                Struct                      coorbital.util.missileConst
%
%  latTgR, lonTgR   [1 x 1]                     Target coordinates (rad)
%
%  hApoStar         [1 x 1]                     Apogee of the max-range arc (m)
%
%  tFlyStar         [1 x 1]                     Flight time of the max-range
%                                               arc (s)
%
%% Outputs:
%
%  br               Struct                      Branch record:
%                                               name, why   char
%                                               exists      logical
%                                               loftR       (rad)
%                                               rngM, missM (m)
%                                               hApoM       (m)
%                                               tFlyS       (s)
%                                               vImpM       (m/s)
%                                               gamImR      (rad)
%                                               gamBoR      (rad)
%                                               vBoM        (m/s)
%                                               hBoM        (m)
%                                               epsBo       (J/kg) specific
%                                                           orbital energy at
%                                                           burnout
%                                               measured    char, the branch
%                                                           this arc MEASURES
%                                                           as
%                                               bandLo,bandHi (m)
%                                               iterations, nEval
%                                               solveInfo   the rangeSolve record
%                                               traj        struct or []
%
%% Revision History:
%  Michael Casey                                                08/07/2026
%  Copyright 2026 Coorbital, Inc.
%% ------------------------ Begin Code Sequence ---------------------------

    [loftR,rngAch,sv] = coorbital.util.rangeSolve(rngReq,fRange,loftAR,loftBR,tolRngM);
           br.name    = name;
           br.exists  = sv.converged;
           br.loftR   = loftR;
           br.rngM    = rngAch;
           br.bandLo  = sv.fMin;
           br.bandHi  = sv.fMax;
      br.iterations   = sv.iterations;
           br.nEval   = sv.nEval;
           br.nProp   = sv.nEval;
       br.solveInfo   = sv;
           br.traj    = [];

%% A branch that did not converge has no trajectory to describe, and every
%% derived field is left NaN rather than filled from the nearest miss: a nearest
%% miss printed in a results table is indistinguishable from a solution:
    if ~sv.converged
           br.why     = sprintf('does not reach %.3f km',rngReq./1000);
           br.missM   = NaN;
           br.hApoM   = NaN;
           br.tFlyS   = NaN;
           br.vImpM   = NaN;
           br.gamImR  = NaN;
           br.gamBoR  = NaN;
           br.vBoM    = NaN;
           br.hBoM    = NaN;
           br.epsBo   = NaN;
        br.measured   = 'none';
        return;
    end

%% One propagation at the converged loft angle, kept:
      [rngFly,traj]   = fRange(loftR);
    assert(abs(rngFly - rngAch) < 1e-6, ...
        ['re-flying the solved %s loft angle gave %.6f m of range against the ' ...
         'solver''s %.6f m; the propagation is not repeatable and nothing ' ...
         'below can be trusted'],name,rngFly,rngAch);
           br.traj    = traj;
           br.nProp   = sv.nEval + 1;
           br.why     = 'reaches the target';
              kBO     = find(traj.phaseIdx == 1,1,'last');
              kApo    = find(traj.phaseIdx == 2,1,'last');
           br.missM   = rI.*coorbital.util.greatCircle(traj.x(end,3),traj.x(end,2), ...
                                                       latTgR,lonTgR);
           br.hApoM   = traj.x(kApo,1) - c.rE;
           br.tFlyS   = traj.t(end);
           br.vImpM   = traj.x(end,4);
           br.gamImR  = traj.x(end,5);
           br.gamBoR  = traj.x(kBO,5);
           br.vBoM    = traj.x(kBO,4);
           br.hBoM    = traj.x(kBO,1) - c.rE;

%% Specific orbital energy at burnout, the quantity the minimum-energy
%% discussion in the header turns on. Measured from the flown burnout state, so
%% it carries whatever the boost actually lost to drag and to gravity:
           br.epsBo   = traj.x(kBO,4).^2./2 - c.muE./traj.x(kBO,1);
        br.measured   = measureBranch(br.hApoM,br.tFlyS,hApoStar,tFlyStar);
end

function [name,agree] = measureBranch(hApoM,tFlyS,hApoStar,tFlyStar)
%% Purpose:
%
%  Name the branch a trajectory actually flew, from ITS OWN apogee and flight
%  time against the max-range arc's. The max-range arc lies between the lofted
%  and depressed arcs in both quantities, so either one identifies the branch
%  and the two together cross-check each other.
%
%  This exists so that a lofted or depressed request is VERIFIED rather than
%  taken on trust. A bracket lying on one side of the maximum should keep the
%  solution on that side, but a max-range angle located slightly wrong breaks
%  that silently and the branch is the whole product of the script.
%
%% Inputs:
%
%  hApoM            [1 x 1]                     Apogee of the arc in question
%                                               (m)
%
%  tFlyS            [1 x 1]                     Its flight time (s)
%
%  hApoStar         [1 x 1]                     Apogee of the max-range arc (m)
%
%  tFlyStar         [1 x 1]                     Flight time of the max-range
%                                               arc (s)
%
%% Outputs:
%
%  name             Char [1 x n]                'lofted' or 'depressed', taken
%                                               from the apogee, which is the
%                                               more direct of the two
%
%  agree            [1 x 1] logical             True when the flight time says
%                                               the same thing as the apogee
%
%% Revision History:
%  Michael Casey                                                08/07/2026
%  Copyright 2026 Coorbital, Inc.
%% ------------------------ Begin Code Sequence ---------------------------

              name = 'depressed';
    if hApoM > hApoStar
              name = 'lofted';
    end
             byTim = 'depressed';
    if tFlyS > tFlyStar
             byTim = 'lofted';
    end
             agree = strcmp(name,byTim);
end

function printPair(label,vDep,vLof,dep,lof,fmt)
%% Purpose:
%
%  Print one row of the two-branch comparison table, blanking the entry of any
%  branch that does not reach the target. A NaN printed as a number in a results
%  table reads as a measurement; a dash reads as what it is.
%
%% Inputs:
%
%  label            Char [1 x n]                Row label with its unit
%
%  vDep             [1 x 1]                     Depressed-branch value
%
%  vLof             [1 x 1]                     Lofted-branch value
%
%  dep              Struct                      Depressed record; only exists
%                                               is read
%
%  lof              Struct                      Lofted record; only exists is
%                                               read
%
%  fmt              Char [1 x m]                Field format, sized to 16
%                                               columns, e.g. '%16.4f' or
%                                               '%16d'. A count printed as
%                                               11.0000 reads as a measurement
%                                               that was not one
%
%% Outputs:
%
%  none                                         Prints one line
%
%% Revision History:
%  Michael Casey                                                08/07/2026
%  Copyright 2026 Coorbital, Inc.
%% ------------------------ Begin Code Sequence ---------------------------

              txtD = '               -';
              txtL = '               -';
    if dep.exists
              txtD = sprintf(fmt,vDep);
    end
    if lof.exists
              txtL = sprintf(fmt,vLof);
    end
    fprintf('    %-22s %s %s\n',label,txtD,txtL);
end

function s = yesNo(tf)
%% Purpose:
%
%  Render a logical as the word a reader of a results table wants to see.
%
%% Inputs:
%
%  tf               [1 x 1] logical             The flag
%
%% Outputs:
%
%  s                Char [1 x n]                'yes' or 'NO'
%
%% Revision History:
%  Michael Casey                                                08/07/2026
%  Copyright 2026 Coorbital, Inc.
%% ------------------------ Begin Code Sequence ---------------------------

                 s = 'NO';
    if tf
                 s = 'yes';
    end
end

function s = lowerOf(nameA,nameB,vA,vB)
%% Purpose:
%
%  Name whichever of two branches carries the smaller value. Used by the
%  minimum-energy paragraph to report what the literal least-burnout-energy
%  reading WOULD have selected, which is not what the script selects.
%
%% Inputs:
%
%  nameA, nameB     Char                        Branch names
%
%  vA, vB           [1 x 1]                     Their values
%
%% Outputs:
%
%  s                Char [1 x n]                nameA if vA <= vB, else nameB
%
%% Revision History:
%  Michael Casey                                                08/07/2026
%  Copyright 2026 Coorbital, Inc.
%% ------------------------ Begin Code Sequence ---------------------------

                 s = nameB;
    if vA <= vB
                 s = nameA;
    end
end

function s = nearerOf(nameA,nameB,dA,dB)
%% Purpose:
%
%  Name whichever of two branches lies nearer a reference, given the two
%  distances. Used by the gammaStar diagnostic.
%
%% Inputs:
%
%  nameA, nameB     Char                        Branch names
%
%  dA, dB           [1 x 1]                     Their distances from the
%                                               reference, non-negative
%
%% Outputs:
%
%  s                Char [1 x n]                nameA if dA <= dB, else nameB
%
%% Revision History:
%  Michael Casey                                                08/07/2026
%  Copyright 2026 Coorbital, Inc.
%% ------------------------ Begin Code Sequence ---------------------------

                 s = nameB;
    if dA <= dB
                 s = nameA;
    end
end

function s = loftedOrNot(ratio)
%% Purpose:
%
%  The one-word verdict the overall block prints on the apogee-to-range ratio. A
%  minimum-energy ballistic arc sits near 0.25, so a ratio well above that is a
%  lofted arc and well below it a depressed one.
%
%% Inputs:
%
%  ratio            [1 x 1]                     Apogee over surface range (-)
%
%% Outputs:
%
%  s                Char [1 x n]                'LOFTED', 'DEPRESSED' or
%                                               'near minimum-energy'
%
%% Revision History:
%  Michael Casey                                                08/07/2026
%  Copyright 2026 Coorbital, Inc.
%% ------------------------ Begin Code Sequence ---------------------------

                 s = 'near minimum-energy';
    if ratio > 0.28
                 s = 'LOFTED';
    elseif ratio < 0.22
                 s = 'DEPRESSED';
    end
end

function d = bandShortfall(fTarget,fMin,fMax)
%% Purpose:
%
%  Signed distance from a target value to the nearer edge of an achievable band,
%  for the refusal message. Positive when the target is beyond the top of the
%  band, negative when it is below the bottom, zero when it is inside.
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
%  a                [1 x 1]                     The same angle in (-pi,pi] (rad)
%
%% Revision History:
%  Michael Casey                                                08/07/2026
%  Copyright 2026 Coorbital, Inc.
%% ------------------------ Begin Code Sequence ---------------------------

                 a = mod(a + pi,2.*pi) - pi;
end

function [pk,kAt] = maxOver(v,mask)
%% Purpose:
%
%  Largest value of a history within a masked subset of its samples, and the
%  index into the FULL history at which it occurs. Used so that a boost peak and
%  a re-entry peak can be reported separately without either search wandering
%  into the other's phase.
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

function [why,ok] = whyBurnout(tEnd,mEnd,mBurnoutT,tMaxBoost)
%% Purpose:
%
%  Diagnose why the boost phase stopped. Nominal termination is propellant
%  exhaustion, which coorbital.prop.eventBurnout detects as the mass state
%  descending through the total burnout mass. THE BOOST IS NEVER CUT SHORT IN
%  THIS SCRIPT -- the loft angle is the ranging control, not the cutoff time --
%  so a horizon timeout or anything else is a fault and must be named rather
%  than left to look like a completed burn.
%
%% Inputs:
%
%  tEnd             [1 x 1]                     Time at the end of phase 1 (s)
%
%  mEnd             [1 x 1]                     Mass state there (kg)
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

function [why,ok] = whyApogee(tPhase,tEnd,gamEnd,tMaxCoast)
%% Purpose:
%
%  Diagnose why the ascending coast stopped. Nominal termination is apogee,
%  which coorbital.prop.eventApogee detects as the flight-path angle crossing
%  zero descending. A coast that simply ran out of horizon while still climbing
%  is the failure this exists to name, and on a heavily lofted arc it is a real
%  possibility rather than a formality.
%
%% Inputs:
%
%  tPhase           [1 x 1]                     Duration of the coast (s)
%
%  tEnd             [1 x 1]                     Cumulative time at its end (s)
%
%  gamEnd           [1 x 1]                     Flight-path angle there (rad)
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

                ok = false;
    if abs(gamEnd) < 1e-9
               why = sprintf(['flight path angle reached zero at t = %.3f s, ' ...
                              'apogee (nominal)'],tEnd);
                ok = true;
    elseif tPhase >= tMaxCoast - 1e-6
               why = sprintf(['hit the %.0f s coast horizon still at gamma = %.4f deg; ' ...
                              'apogee was never reached'],tMaxCoast,rad2deg(gamEnd));
    else
               why = sprintf(['stopped early at t = %.3f s with gamma = %.4f deg ' ...
                              '(integrator failure or an unmodelled event)'], ...
                             tEnd,rad2deg(gamEnd));
    end
end

function [why,ok] = whyImpact(tPhase,tEnd,hEnd,hStopM,hStop,tMaxDesc)
%% Purpose:
%
%  Diagnose why the descent stopped. Nominal termination is the descending
%  crossing of the impact altitude. A descent that ran out of horizon in mid-air
%  is the failure this exists to name.
%
%% Inputs:
%
%  tPhase           [1 x 1]                     Duration of the descent (s)
%
%  tEnd             [1 x 1]                     Cumulative time at its end (s)
%
%  hEnd             [1 x 1]                     Altitude there (m)
%
%  hStopM           [1 x 1]                     Impact altitude (m)
%
%  hStop            [1 x 1]                     Impact altitude (km), for the
%                                               message
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

                ok = false;
    if abs(hEnd - hStopM) < 1e-3
               why = sprintf(['reached the %.1f km impact altitude at t = %.3f s ' ...
                              '(nominal)'],hStop,tEnd);
                ok = true;
    elseif tPhase >= tMaxDesc - 1e-6
               why = sprintf(['hit the %.0f s descent horizon at %.3f km, still ' ...
                              'flying'],tMaxDesc,hEnd./1000);
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
%  The factor holds the APPARENT apogee at or under 0.3 rE, capped at 30 and
%  floored at 2. THIS SCRIPT NEEDS THE RULE RATHER THAN A FIXED FACTOR because
%  its two branches differ by a factor of ten in apogee -- 208 km depressed
%  against 2118 km lofted on the shipped case -- and no single number suits
%  both: at the 3x BM/run_ballistic hard-codes, the depressed arc is invisible.
%  The full rationale for the cap, the floor and the 0.3 rE target, and for
%  where the invariant stops holding, is in the exagFor header of
%  HGV/run_target.m; this is the same rule.
%
%% Inputs:
%
%  hPeakM           [1 x 1]                     Peak altitude above the
%                                               reference sphere, over the whole
%                                               flight (m)
%
%  rE               [1 x 1]                     Reference sphere radius (m),
%                                               from coorbital.util.missileConst
%
%% Outputs:
%
%  e                [1 x 1]                     Altitude exaggeration to hand
%                                               coorbital.viz as AltScale (-), a
%                                               whole number in [2 .. 30]
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
%  "45.00N 100.00W" is unambiguous where "45.00, -100.00" leaves the reader to
%  remember which way the signs run.
%
%  A coordinate that ROUNDS to zero takes its hemisphere from the rounded value,
%  so a longitude a thousandth of a degree west of the meridian does not print
%  as a west longitude of zero; and a non-finite coordinate is refused rather
%  than formatted into "NaNN NaNE".
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
%  txt              char                        e.g. '45.00N 100.00W'
%
%% Revision History:
%  Michael Casey                                                08/07/2026
%  Copyright 2026 Coorbital, Inc.
%% ------------------------ Begin Code Sequence ---------------------------

    assert(isscalar(latDeg) && isfinite(latDeg) && ...
           isscalar(lonDeg) && isfinite(lonDeg), ...
        'llStr needs finite scalar coordinates; got %s and %s.', ...
        mat2str(latDeg),mat2str(lonDeg));
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
               txt = sprintf('%.2f%c %.2f%c',abs(latR),hLat,abs(lonR),hLon);
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
