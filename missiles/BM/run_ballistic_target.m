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
%% Note -- THE THREE MODES ARE NOT THREE POINTS ON ONE CURVE:
%
%  TWO of the three fly the FULL BURN and range on the loft angle alone. They
%  are the two arcs described above, one either side of the maximum-range loft
%  angle, and they are what a fixed booster with no thrust termination can do:
%
%      'lofted'      the high, slow, long-flight arc above the max-range angle;
%      'depressed'   the low, fast, short-flight arc below it.
%
%  THE THIRD IS THE TEXTBOOK TRAJECTORY, and it is a different problem. The
%  minimum-energy ballistic trajectory between two points at the same radius is
%  the one that reaches the required free-flight range angle Lambda on the LEAST
%  burnout energy. It is a closed-form Keplerian result, reference [1]:
%
%      V*^2      = (mu/rE) * 2 sin(Lambda/2) / (1 + sin(Lambda/2))
%      gammaStar = 45 deg - Lambda/4
%
%  and it is the trajectory published ballistic-missile range tables are quoted
%  on. NEITHER FULL-BURN ARC IS IT. At full burn the booster delivers a fixed
%  delta-V, so the lofted and depressed arcs leave burnout with essentially the
%  same energy -- measured on the shipped case, 0.55 %% apart, and that
%  difference is boost drag and gravity loss rather than a design variable.
%  V* is a speed the booster OVERSHOOTS, so neither arc can carry it.
%
%  WHICH IS WHY 'minimum-energy' HAS A SECOND CONTROL. The burnout speed has to
%  be free, and the only way this vehicle can lower it is to CUT THE BURN SHORT
%  -- the same thrust-termination control HGV/run_target ranges on. So this mode
%  solves TWO unknowns against TWO residuals:
%
%      unknowns    the LOFT ANGLE and the CUTOFF FRACTION of the full burn;
%      residual 1  achieved range      minus the required range, to tolRangeKm;
%      residual 2  burnout gamma       minus gammaStar, to tolGamDeg.
%
%  THE SOLVE IS NESTED, NOT A 2x2 NEWTON, and the reason is that each level is
%  then MONOTONIC and can reuse coorbital.util.rangeSolve exactly as designed:
%
%      INNER, on the CUTOFF FRACTION at a held loft angle. Range rises
%             monotonically with cutoff -- more burn, more energy, more range --
%             which is the property rangeSolve's bisection needs and the same
%             one HGV/run_target rests on. Solved to tolRangeKm.
%
%      OUTER, on the LOFT ANGLE, driving the burnout gamma of that inner
%             solution to gammaStar. Along the constant-range locus the burnout
%             gamma rises monotonically with loft angle: measured on the shipped
%             case, 11.815 deg at the depressed end through 62.840 deg at the
%             lofted end, over nine equally spaced samples with no reversal.
%             Solved to tolGamDeg, again by rangeSolve.
%
%  A 2x2 damped Newton on the same residuals would also work -- the Jacobian is
%  dominated by its anti-diagonal, range answering mostly to cutoff and gamma
%  mostly to loft, so it is not ill-conditioned -- but it needs a finite
%  difference Jacobian, four propagations a step, and a damping rule, and it can
%  leave the region where the chain propagates at all. The nested form cannot:
%  bisection never leaves its bracket.
%
%  THE OUTER BRACKET IS THE TWO FULL-BURN BRANCH SOLUTIONS, which are solved
%  before it and cost nothing extra. Between them, and only between them, the
%  full burn reaches AT LEAST the required range, so the inner solve is
%  guaranteed a cutoff fraction at or below 1 that lands on it. Where a branch
%  does not exist the corresponding end of the loft bracket is used instead,
%  which is admissible for the same reason: that branch is missing precisely
%  because the full burn OVERSHOOTS the target there.
%
%  WHAT IT COSTS AND WHAT IT DELIVERS. Measured on the shipped 3174.981 km case:
%  11 outer bisection steps and 190 propagations, about 8 s of an 11 s run --
%  the two full-burn arcs are still solved and still reported. Against the classical
%  reference for that range angle -- 687.304 km of apogee, 15.4635 min of flight
%  and gammaStar = 37.8697 deg -- the solved trajectory flies 735.074 km,
%  16.8893 min and 37.8654 deg: +6.95 %%, +9.22 %% and -0.0043 deg.
%
%  THE TWO PERCENTAGES ARE PHYSICS, NOT SOLVER ERROR, and the summary says so.
%  The classical result assumes an IMPULSIVE burn at the impact radius in a
%  VACUUM. This flight burns for 77.6 s and reaches gammaStar 80.56 km up at
%  4833.3 m/s rather than 4970.3 m/s at zero altitude, and a Keplerian arc from
%  THAT burnout state apogees at 735.093 km -- the flown figure to 0.019 km, so
%  the whole 47.8 km gap to the classical arc is where the boost ENDED and not
%  how the solve converged. Coast drag supplies what is left. Expect a few per
%  cent from any finite boost; expect machine precision from neither.
%
%  MINIMUM-ENERGY REQUIRES SEPARATION. A cut-short burn leaves propellant in the
%  booster, and the coast vehicle's mass is fixed before the cutoff is known, so
%  the whole booster -- dry structure and the unburned propellant with it -- is
%  jettisoned at cutoff, exactly as HGV/run_target does. separation = false
%  cannot express that and is refused for this mode rather than flown with a
%  mass the equations of motion would disagree with.
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
%  and that is a targeting decision rather than an aerodynamic one. IT IS NOT
%  THAT 6 deg ABOLISHES THE DEPRESSED BRANCH. Measured, by scanning loft from
%  -140 to 85 deg at a 6 deg clamp: the max-range angle sits at a commanded
%  -42.902 deg with 5211.525 km of range, so the hump is there and so are both
%  branches. What is wrong at 6 deg is that the hump sits 2.902 deg BELOW the
%  shipped loftMin of -40 deg, so the shipped bracket finds its largest range ON
%  an endpoint and the bracketing step below refuses. THE CAUSE OF THAT REFUSAL
%  IS BRACKET WIDTH, NOT THE CLAMP.
%
%  What the clamp decides is WHERE the branches lie and what they can reach. At
%  6 deg the depressed branch spans only 4708.463 to 5211.525 km -- its floor is
%  the range still flown at a commanded -140 deg, by which point the clamp has
%  long since saturated -- so it cannot reach the 3175 km shipped target at any
%  loft angle, and widening loftMin instead of raising alphaMax would clear the
%  bracketing refusal only to meet the unreachable-branch one. At 12 deg the
%  depressed branch reaches down to 1684.117 km and the shipped target is inside
%  it. THE CLAMP IS PAID FOR IN RANGE: the maximum falls from 5211.525 km at
%  6 deg to 5055.302 km at 12 deg, 156.224 km given up to buy a depressed branch
%  that reaches targets a user would actually ask for.
%
%% Note -- TWO LIMITATIONS, one REFUSED outright and one printed in the summary:
%
%  THE AZIMUTH IS EXACT ONLY BECAUSE THE EARTH DOES NOT ROTATE HERE, and
%  earthSpin true is therefore REFUSED rather than flown. With env.omegaE = 0
%  the ground track of an unbanked trajectory is a great circle, so the initial
%  bearing of the launch-to-target arc is the whole answer. Switch earthSpin on
%  and the target is carried east under the vehicle for the whole flight, the
%  vehicle has to be aimed where the target is GOING TO BE, and closing that
%  needs an OUTER AZIMUTH ITERATION around the two branch solves -- PER BRANCH,
%  the two arcs flying for 659 s and 1816 s on the shipped case. This script
%  does not have one. It used to run anyway and print a caution; measured, that
%  produced a 315.690 km miss printed as a converged solution, beside a
%  cross-track WARNING blaming a banked segment when the bank was already zero.
%  It now refuses on the same terms as the envelope refusals -- empty traj,
%  info.refused true, nothing thrown -- and names the iteration it lacks.
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
%                                               (logical), refusedWhy
%                                               ('envelope', 'earthSpin' or
%                                               'minimumEnergy'), rngReqM and
%                                               psiLaunch. An ENVELOPE refusal
%                                               adds loftStarDeg, rngMaxM,
%                                               rngMinM, classical, gamStarR
%                                               and the per-branch records in
%                                               depressed and lofted, so it
%                                               hands back the envelope it
%                                               refused against; a
%                                               MINIMUMENERGY refusal adds all
%                                               of those plus minEnergy, whose
%                                               why and gamma band say what
%                                               could not be reached; an
%                                               earthSpin refusal happens
%                                               before any propagation and adds
%                                               omegaE and tgSpeedM instead. A
%                                               run that was not refused
%                                               additionally carries the flown
%                                               trajectory's own numbers, the
%                                               classical minimum-energy
%                                               reference in classical, the
%                                               minimum-energy record in
%                                               minEnergy (solved false on the
%                                               two full-burn modes), the three
%                                               residuals meGamResR, meApoRelE
%                                               and meTofRelE, and the DISPLAY
%                                               choices the figures were drawn
%                                               with -- altExag, its rule as a
%                                               handle, and the two hemisphere
%                                               captions
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

%% WHICH TRAJECTORY TO FLY. The two FULL-BURN arcs are solved and reported
%% whatever is asked for here, because the trade between them is the point; this
%% chooses what is flown, plotted and returned. Read the THREE MODES note in the
%% header: two of these are full-burn arcs either side of maximum range and the
%% third is the textbook trajectory, which is a different solve on a second
%% control and is NOT a member of the full-burn family:
            branch = 'minimum-energy'; %'minimum-energy' | 'lofted' | 'depressed'
                                       %  minimum-energy: the TEXTBOOK minimum-energy ballistic
                                       %                  trajectory -- the one published range
                                       %                  tables are quoted on. Solves the LOFT
                                       %                  ANGLE and the CUTOFF FRACTION together
                                       %                  so that the range is met AND the
                                       %                  burnout flight-path angle equals
                                       %                  gammaStar = 45 deg - Lambda/4. The
                                       %                  burn is CUT SHORT: it is the only way
                                       %                  to reach the burnout speed the
                                       %                  textbook result asks for. Requires
                                       %                  separation = true
                                       %  lofted:         the FULL-BURN high, slow, long-flight
                                       %                  arc, above the max-range loft angle
                                       %  depressed:      the FULL-BURN low, fast, short-flight
                                       %                  arc, below it

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

%% MINIMUM-ENERGY MODE ONLY. Ignored by 'lofted' and 'depressed', which fly the
%% full burn and have no cutoff to solve for. See the THREE MODES note:
        cutFracMin = 0.5;              %-, floor of the thrust-termination bracket the INNER
                                       %   solve bisects, as a fraction of the full burn
                                       %   [0.2 .. 0.95]. Both ends of that bracket are
                                       %   propagated before the first bisection step, so the
                                       %   floor must be a cutoff the chain can still FLY: too
                                       %   early and the vehicle never reaches apogee, the event
                                       %   never fires, and flyLoft refuses the propagation
                                       %   rather than returning a range from a flight that did
                                       %   not happen. Measured on this vehicle, 0.5 of the
                                       %   80.518 s burn flies between 13 and 127 km depending
                                       %   on loft angle, which is well clear of that failure
         tolGamDeg = 0.01;             %deg, convergence tolerance of the OUTER solve on the
                                       %     BURNOUT FLIGHT-PATH ANGLE against gammaStar
                                       %     [0.001 .. 1]. This is the tolerance on the residual
                                       %     that makes the trajectory minimum-energy at all, so
                                       %     it is the one to tighten if the classical
                                       %     comparison printed in the summary matters. Each
                                       %     halving costs one more OUTER step, and an outer
                                       %     step is a whole inner range solve -- about 16
                                       %     propagations -- so it is NOT as cheap as tolRangeKm

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
                                       %     flight path can be pushed over. It decides WHERE
                                       %     the two branches lie: at BM/run_ballistic's 6 deg
                                       %     the depressed branch still exists but bottoms out
                                       %     at 4708.463 km and cannot reach the shipped target,
                                       %     and its hump sits below loftMin so the bracket
                                       %     below refuses on WIDTH. 12 deg buys a depressed
                                       %     branch reaching 1684.117 km and pays 156.224 km of
                                       %     maximum range for it -- see the header note
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
         earthSpin = false;            %MUST BE false. true enables the Coriolis and centrifugal
                                       %terms, which INVALIDATE the closed-form azimuth, and the
                                       %run is then REFUSED rather than flown: aiming at where
                                       %the target is going to be needs an outer azimuth
                                       %iteration this script does not have. Measured before the
                                       %refusal existed, on the shipped geometry: a 315.690 km
                                       %miss printed as though it were a converged solution

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
                                       %   capped at 30 and floored at 2. 'auto' is worth asking
                                       %   for when the two branches are being compared: they
                                       %   differ by ten times in apogee and it draws them at 2x
                                       %   and 9x, where one fixed factor suits neither
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
                      'tolLoftDeg','nScanLoft','tolRangeKm','cutFracMin', ...
                      'tolGamDeg','pitchTime', ...
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
        cutFracMin = overrideOf(opts,'cutFracMin',cutFracMin);
         tolGamDeg = overrideOf(opts,'tolGamDeg',tolGamDeg);
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
           altExag = overrideOf(opts,'altExag',altExag);
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
           tolGamR = deg2rad(tolGamDeg);

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
    assert(isscalar(cutFracMin) && isfinite(cutFracMin) && ...
           cutFracMin > 0 && cutFracMin < 1, ...
        ['cutFracMin must lie strictly between 0 and 1; got %s. It is the ' ...
         'FLOOR of a bracket whose ceiling is the full burn, and a floor at ' ...
         'or above 1 leaves the minimum-energy solve nothing to bisect.'], ...
        mat2str(cutFracMin));
    assert(tolGamR > 0,'tolGamDeg must be positive.');
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
        'coorbital:runBallisticTarget:coincidentPoints', ...
        ['the launch point and the target are the same point; there is nothing ' ...
         'to solve. coorbital.util.greatCircleBearing would refuse this a few ' ...
         'lines further down, but it can only speak of coincident points, not ' ...
         'of latTarget and lonTarget.']);

%% The display exaggeration, checked here rather than at the figures: a
%% mistyped one would otherwise cost a full two-branch solve before it was
%% noticed. Either a positive number, used as given, or the char 'auto', which
%% defers to the adaptive rule once the flown apogee is known:
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
             isMinE = strcmp(branch,'minimum-energy');

%% MINIMUM-ENERGY CUTS THE BURN SHORT, so the booster it throws away still has
%% propellant in it. The coast vehicle's mass is built once, before the cutoff
%% is known, and the six-state glide equations divide by THAT mass rather than
%% by the mass state -- so a retained booster of unknown mass is not something
%% this chain can express. Refused here rather than flown against a mass the
%% equations disagree with:
    assert(~isMinE || logical(separation), ...
        ['branch = ''minimum-energy'' needs separation = true. It solves a ' ...
         'THRUST-TERMINATION FRACTION as well as a loft angle, so the booster ' ...
         'it drops carries unburned propellant, and the mass flying the coast ' ...
         'is therefore only well defined if the whole booster goes overboard. ' ...
         'Ask for ''lofted'' or ''depressed'' to keep a dead booster attached; ' ...
         'both burn to exhaustion and leave nothing unburned to account for.']);

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
         cfg.tBurn = tBurn;
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

%% The FULL-BURN range function, the one the max-range bracketing and both
%% branch solves ride on. The second argument of flyLoft is the cutoff fraction
%% and 1 means "no commanded cutoff": phase 1 keeps its horizon tspan and ends
%% on propellant exhaustion, exactly as it always has. Only the minimum-energy
%% solve ever passes anything else:
            fRange = @(loftR) flyLoft(loftR,1,cfg,rI);

%% -----------------------------------------------------------------
%% Refuse a ROTATING Earth, before a single trajectory is propagated
%% -----------------------------------------------------------------
%% The closed-form azimuth above is the whole answer only while the ground does
%% not move under the vehicle. With env.omegaE non-zero it is not an answer at
%% all, and everything downstream of it -- the branch solves, the miss, the
%% cross-track, the limitations paragraph -- would be describing a target that
%% stood still. This used to run anyway behind a caution, and what it produced
%% was three mutually contradictory statements around a 315.690 km miss. It
%% refuses instead, on the same terms as the envelope refusals below: nothing is
%% thrown, traj comes back empty and info.refused is true. It refuses HERE,
%% before the bracketing, because 57 propagations spent on an answer that will
%% be discarded are 57 propagations wasted:
    if env.omegaE ~= 0
            tgSpdM = c.omegaE.*rI.*cos(latTargetR);
        fprintf('\n');
        fprintf('===== Ballistic point-to-point targeting: REFUSED =======\n');
        fprintf('  earthSpin is TRUE. This script cannot target a ROTATING Earth.\n');
        fprintf('\n');
        fprintf('    launch           %10.4f %-5s %10.4f deg  (lat, lon)\n', ...
                latLaunch,'deg',lonLaunch);
        fprintf('    target           %10.4f %-5s %10.4f deg  (lat, lon)\n', ...
                latTarget,'deg',lonTarget);
        fprintf('    required range   %10.3f %-5s (great circle on the r = %.3f km sphere,\n', ...
                rngReq./1000,'km',rI./1000);
        fprintf('                                      measured to where the target is NOW)\n');
        fprintf('    launch azimuth   %10.6f %-5s (the closed form, computed and then NOT\n', ...
                rad2deg(psiLaunch),'deg');
        fprintf('                                      used: it is not the answer here)\n');
        fprintf('    env.omegaE     %.6e rad/s (Coriolis and centrifugal terms ON)\n', ...
                env.omegaE);
        fprintf('    target moves     %10.2f %-5s (eastward ground speed under the target at\n', ...
                tgSpdM,'m/s');
        fprintf('                                      %.4f deg latitude; every second of\n',latTarget);
        fprintf('                                      flight carries the aim point this far)\n');
        fprintf('\n');
        fprintf('  WHAT IS MISSING IS AN OUTER AZIMUTH ITERATION, and this script does not\n');
        fprintf('  have one. On a turning Earth the vehicle must be aimed where the target is\n');
        fprintf('  GOING TO BE, which makes the azimuth depend on the flight time, the flight\n');
        fprintf('  time depend on the loft angle, and the loft angle depend on the azimuth.\n');
        fprintf('  Worse here than in HGV/run_target: the loop would have to run PER BRANCH,\n');
        fprintf('  the lofted and depressed arcs flying for very different times over the\n');
        fprintf('  same geometry, so the two branches do not even share an aim point.\n');
        fprintf('  Running anyway is what this refusal replaces. It converged, reported a\n');
        fprintf('  315.690 km miss on the shipped geometry as though it were a solution, and\n');
        fprintf('  blamed the cross-range on a banked segment while the bank was zero. Set\n');
        fprintf('  earthSpin false, or bring an azimuth iteration.\n');
        fprintf('=========================================================\n\n');
              traj = [];
       info.refused = true;
    info.refusedWhy = 'earthSpin';
      info.branchAsked = branch;
       info.rngReqM = rngReq;
     info.psiLaunch = psiLaunch;
       info.omegaE  = env.omegaE;
      info.tgSpeedM = tgSpdM;
        if nargout == 0
            clear traj info;
        end
        return;
    end

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
%% dLoftDep and dLoftLof are how far each full-burn branch sits from the
%% max-range loft angle. They are computed here, once, whatever the selector
%% says, because the two-arc table prints them on every run: they are the
%% clearest single measure of how far apart the two arcs are as CONTROLS rather
%% than as outcomes. Nothing selects on them any more -- see the THREE MODES
%% note for why "nearest the max-range angle" was the wrong reading of
%% minimum-energy and what replaced it:
          dLoftDep = NaN;
          dLoftLof = NaN;
    if dep.exists
          dLoftDep = abs(dep.loftR - loftStarR);
    end
    if lof.exists
          dLoftLof = abs(lof.loftR - loftStarR);
    end
          bothHere = dep.exists && lof.exists;

%% The CLASSICAL reference for this range angle, computed from the geometry
%% alone and from nothing this script flew. It is what the minimum-energy mode
%% is solved against and what its residuals are reported against; the two
%% full-burn modes print it as a yardstick and are not expected to match it:
             refME = minEnergyRef(angReqR,c);
          gamStarR = refME.gamR;

%% The selection itself. A named branch takes that branch or nothing;
%% minimum-energy is not a branch at all and is solved below, after the envelope
%% refusal has established that the loft bracket it needs exists:
          pickName = branch;
          pickWhy  = '';
         pickShort = 'the branch selector';
    switch branch
        case 'minimum-energy'
         pickShort = 'solved on loft AND cutoff, not chosen between the two arcs';
          pickWhy  = 'the textbook minimum-energy trajectory for this range angle';
        case 'lofted'
         pickShort = 'named in the user block';
          pickWhy  = 'the arc the user asked for by name, not a solver choice';
        case 'depressed'
         pickShort = 'named in the user block';
          pickWhy  = 'the arc the user asked for by name, not a solver choice';
    end
    if strcmp(pickName,'lofted')
              pick = lof;
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
%% minimum-energy is exempt from the badPick test: it does not fly either
%% branch, so a branch that does not exist costs it only ONE END of the loft
%% bracket it searches, and the other end of the user's bracket serves instead.
%% It is NOT exempt from noBranch -- if the full burn cannot reach this target
%% at any loft angle then a burn CUT SHORTER certainly cannot:
          noBranch = ~dep.exists && ~lof.exists;
          badPick  = ~noBranch && ~isMinE && ~pick.exists;
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
        fprintf('  was asked for. Ask for the other one.\n');
        end
        if isMinE
        fprintf('\n');
        fprintf('  AND MINIMUM-ENERGY IS THE HARDER ASK, NOT THE EASIER ONE. The textbook\n');
        fprintf('  trajectory for this %.4f deg range angle wants a burnout speed of\n', ...
                rad2deg(angReqR));
        fprintf('  V* = %.1f m/s at gamma* = %.4f deg, which is BELOW what a full burn\n', ...
                refME.V,rad2deg(refME.gamR));
        fprintf('  delivers -- it is reached by CUTTING THE BURN SHORT. A trajectory that\n');
        fprintf('  carries less energy than the full burn cannot fly further than the full\n');
        fprintf('  burn does, so a target the full burn cannot reach is beyond this mode by\n');
        fprintf('  a wider margin than it is beyond ''lofted'' or ''depressed''. The\n');
        fprintf('  classical arc it would have flown: %.3f km of apogee in %.3f min.\n', ...
                refME.hApoM./1000,refME.tofS./60);
        end
        fprintf('=========================================================\n\n');
              traj = [];
       info.refused = true;
    info.refusedWhy = 'envelope';
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
     info.classical = refME;
      info.gamStarR = gamStarR;
        if nargout == 0
            clear traj info;
        end
        return;
    end

%% -----------------------------------------------------------------
%% The MINIMUM-ENERGY solve: two unknowns, two residuals, nested
%% -----------------------------------------------------------------
%% Only now, with the envelope refusal past, is the loft bracket the outer solve
%% needs known to exist. See the THREE MODES note for the method and for why the
%% two levels are nested rather than Newton-solved together:
                me = struct('solved',false);
    if isMinE
                me = minEnergySolve(rngReq,gamStarR,cfg,rI,c, ...
                                    latTargetR,lonTargetR,dep,lof, ...
                                    loftLoR,loftHiR,tolRngM,tolGamR,cutFracMin);
        if ~me.solved
        fprintf('\n');
        fprintf('===== Ballistic point-to-point targeting: REFUSED =======\n');
        fprintf('  The MINIMUM-ENERGY trajectory for this geometry is not in reach.\n');
        fprintf('\n');
        fprintf('    launch           %10.4f %-5s %10.4f deg  (lat, lon)\n', ...
                latLaunch,'deg',lonLaunch);
        fprintf('    target           %10.4f %-5s %10.4f deg  (lat, lon)\n', ...
                latTarget,'deg',lonTarget);
        fprintf('    required range   %10.3f %-5s (great circle on the r = %.3f km sphere,\n', ...
                rngReq./1000,'km',rI./1000);
        fprintf('                                       a %.4f deg range angle)\n',rad2deg(angReqR));
        fprintf('    gamma* wanted    %10.4f %-5s (45 - Lambda/4, the burnout flight-path\n', ...
                rad2deg(gamStarR),'deg');
        fprintf('                                       angle that makes the arc minimum-energy)\n');
        fprintf('    V* wanted        %10.1f %-5s (the burnout speed that goes with it)\n', ...
                refME.V,'m/s');
        if isfinite(me.gamLoR) && isfinite(me.gamHiR)
        fprintf('    gamma reachable  %10.4f to %.4f deg  (over the loft bracket %.4f to\n', ...
                rad2deg(me.gamLoR),rad2deg(me.gamHiR),rad2deg(me.loftAR));
        fprintf('                                       %.4f deg, each end already RANGE-SOLVED\n', ...
                rad2deg(me.loftBR));
        fprintf('                                       on its cutoff fraction)\n');
        else
        fprintf('    gamma reachable      NOT MEASURED  (the outer bisection never ran; there\n');
        fprintf('                                       was no loft interval to run it on)\n');
        end
        fprintf('\n');
        fprintf('  %s\n',me.why);
        fprintf('  The two FULL-BURN arcs below DO reach this target and are unaffected; ask\n');
        fprintf('  for ''lofted'' or ''depressed'' to fly one of them. What is refused is the\n');
        fprintf('  claim that either of them is the minimum-energy trajectory, which is the\n');
        fprintf('  claim this mode exists to stop being made silently.\n');
        fprintf('\n');
        fprintf('    depressed arc    %10.4f deg loft, burnout gamma %.4f deg  [%s]\n', ...
                rad2deg(dep.loftR),rad2deg(dep.gamBoR),dep.why);
        fprintf('    lofted arc       %10.4f deg loft, burnout gamma %.4f deg  [%s]\n', ...
                rad2deg(lof.loftR),rad2deg(lof.gamBoR),lof.why);
        fprintf('=========================================================\n\n');
              traj = [];
       info.refused = true;
    info.refusedWhy = 'minimumEnergy';
  info.branchAsked  = branch;
       info.rngReqM = rngReq;
     info.psiLaunch = psiLaunch;
   info.loftStarDeg = rad2deg(loftStarR);
       info.rngMaxM = rngMax;
       info.rngMinM = rngMin;
     info.hApoStarM = hApoStar;
     info.tFlyStarS = tFlyStar;
     info.depressed = dep;
        info.lofted = lof;
      info.maxRange = mr;
     info.classical = refME;
      info.gamStarR = gamStarR;
    info.minEnergy  = me;
            if nargout == 0
                clear traj info;
            end
            return;
        end
              pick = me;
          pickName = 'minimum-energy';
          pickWhy  = sprintf(['SOLVED, not selected: loft %.4f deg, burn cut at %.6f ' ...
                              'of full'],rad2deg(me.loftR),me.cutFrac);
    end

%% -----------------------------------------------------------------
%% Fly the chosen trajectory
%% -----------------------------------------------------------------
%% Whichever solve produced it already flew this trajectory and kept it, so
%% nothing is re-integrated here. The check is that the kept trajectory is the
%% one the solver reported, which a struct copied from the wrong record would
%% fail:
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
%% MINIMUM-ENERGY IS NOT MEASURED THIS WAY, and pretending otherwise would print
%% a falsehood. The test above compares against the FULL-BURN max-range arc, and
%% the minimum-energy trajectory does not burn to exhaustion: it is lower and
%% quicker than that arc for the same reason a smaller booster would be, so it
%% would always "measure as depressed" while being nothing of the sort. Its own
%% verification is the residual pair against gamma* and the required range, which
%% the minimum-energy paragraph below prints:
         flownName = 'minimum-energy';
        flownAgree = true;
          branchOK = true;
    if ~isMinE
        [flownName,flownAgree] = measureBranch(pick.hApoM,pick.tFlyS, ...
                                               hApoStar,tFlyStar);
          branchOK = strcmp(flownName,pickName) && flownAgree;
    end

%% THE ONE CASE IN WHICH THE MEASUREMENT LEGITIMATELY DISAGREES, and it is worth
%% naming because otherwise the caution below misdiagnoses it. Apogee rises
%% monotonically with loft angle and each branch is solved on a bracket with the
%% max-range angle at one end, so a solution strictly inside its bracket always
%% measures as the branch it was solved on. The exception is the shared endpoint:
%% when the required range comes within the range tolerance of the MAXIMUM,
%% coorbital.util.rangeSolve short-circuits at that endpoint and BOTH branches
%% return the max-range arc itself. The lofted request then measures as
%% depressed, because the apogee test is strict and the apogee is equal rather
%% than greater. That is a degeneracy of the request, not a mislocated maximum,
%% and the summary says so instead of blaming the bracketing:
        atMaxRange = pick.loftR == loftStarR;

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
%% The commanded cutoff time, or NaN when the burn was not cut short. It is what
%% tells whyBurnout apart a boost that ended on a COMMANDED thrust termination
%% -- nominal for minimum-energy -- from one that hit a horizon, which is a
%% fault. Passing the full-burn case a NaN is deliberate: every comparison
%% against NaN is false, so the cutoff clause simply cannot fire there:
              tCut = NaN;
            meProp = 0;
    if isMinE
              tCut = me.tCutS;
            meProp = me.nProp;
    end
        [why1,ok1] = whyBurnout(traj.t(kBO),massS(kBO),mBurnoutT,tMaxBoost,tCut);
        [why2,ok2] = whyApogee(traj.t(kApo) - traj.t(kBO),traj.t(kApo), ...
                               traj.x(kApo,5),tMaxCoast);
        [why3,ok3] = whyImpact(traj.t(nS) - traj.t(kApo),traj.t(nS), ...
                               traj.x(nS,1) - c.rE,hStopM,hStop,tMaxDesc);
             allOK = ok1 && ok2 && ok3;

%% How far the ground under the target would have moved during this flight had
%% the Earth been turning. Quantifies the non-rotating caveat instead of merely
%% stating it, and it is quantified for BOTH branches because their flight times
%% differ by a factor of nearly three. A branch that does not reach the target
%% has no flight time, so the pair is only reported when both exist -- a NaN
%% printed in that sentence would read as a measurement:
          driftDep = c.omegaE.*rI.*cos(latTargetR).*dep.tFlyS;
          driftLof = c.omegaE.*rI.*cos(latTargetR).*lof.tFlyS;

%% Earth rotation is necessarily OFF from here down: earthSpin true is refused
%% above, before the bracketing, and the assert at the head of the report
%% re-checks it so that no paragraph below can describe a run it did not
%% describe:
           spinTxt = 'OFF';
            sepTxt = 'booster JETTISONED at burnout';
    if isMinE
            sepTxt = 'WHOLE booster JETTISONED at cutoff, unburned propellant with it';
    elseif ~separation
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

%% EVERY STATEMENT THE SUMMARY MAKES ABOUT THE AZIMUTH, the cross-range and the
%% great-circle ground track is conditional on a non-rotating Earth. earthSpin
%% true is refused above precisely so that the prose below cannot contradict the
%% run it is describing; this re-checks the condition at the point of use rather
%% than trusting a refusal two hundred lines away:
    assert(env.omegaE == 0, ...
        ['the summary is about to describe a closed-form azimuth on a ' ...
         'non-rotating Earth while env.omegaE = %.6e rad/s; the earthSpin ' ...
         'refusal above should have made this unreachable'],env.omegaE);

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
        fprintf('                   solve CONVERGED and the vehicle still missed.\n');
        if bankRad ~= 0
        fprintf('                   THE CAUSE IS THE %.1f deg BANK: a banked segment turns the\n', ...
                bankAngle);
        fprintf('                   heading and walks the track off the arc it left on. Set\n');
        fprintf('                   bankAngle to zero, or accept the cross-range -- closing it\n');
        fprintf('                   needs an outer azimuth iteration this script does not have.\n');
        else
        fprintf('                   THE BANK ANGLE IS ALREADY ZERO, so a banked segment is NOT\n');
        fprintf('                   the cause and setting bankAngle to zero would change\n');
        fprintf('                   nothing. An unbanked track over a NON-ROTATING sphere\n');
        fprintf('                   cannot leave the great circle it departed on, and a\n');
        fprintf('                   rotating one is refused above, so the sideways motion came\n');
        fprintf('                   from neither. Look at the guidance schedules: something is\n');
        fprintf('                   commanding a sigma other than bankAngle.\n');
        end
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
    printPair('loft from max (deg)'   ,rad2deg(dLoftDep)    ,rad2deg(dLoftLof)    ,dep,lof,'%16.4f');
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
    fprintf('  WHAT WAS FLOWN\n');
    fprintf('    asked for        %-16s (%s)\n',branch,pickShort);
    fprintf('    on the grounds that it is %s\n',pickWhy);
    fprintf('    flew             %-16s (the trajectory this run actually solved for)\n',pickName);
    if isMinE
    fprintf('    NOT MEASURED against the max-range arc, unlike the two full-burn modes:\n');
    fprintf('    that arc burns to exhaustion and this one does not, so it would always\n');
    fprintf('    read as "depressed" while being nothing of the sort. The minimum-energy\n');
    fprintf('    paragraph below carries the residuals that DO verify it.\n');
    else
    fprintf('    MEASURED as      %-16s (from the FLOWN apogee %.3f km against the\n', ...
            flownName,pick.hApoM./1000);
    fprintf('                                      max-range arc''s %.3f km, and the flown\n', ...
            hApoStar./1000);
    fprintf('                                      %.3f s flight time against its %.3f s)\n', ...
            pick.tFlyS,tFlyStar);
    end
    if ~flownAgree
        fprintf('  *** CAUTION ***  the apogee and the flight time DISAGREE about which branch\n');
        fprintf('                   this is. One of them puts the flight above the max-range\n');
        fprintf('                   arc and the other below it, which a genuine lofted or\n');
        fprintf('                   depressed arc cannot do. Suspect the max-range angle.\n');
    end
    if ~branchOK && atMaxRange
        fprintf('  *** CAUTION ***  the flown trajectory MEASURES as the %s branch while the\n', ...
                flownName);
        fprintf('                   solve was run on the %s bracket, because BOTH brackets\n', ...
                pickName);
        fprintf('                   converged at their SHARED endpoint, the %.4f deg\n', ...
                rad2deg(loftStarR));
        fprintf('                   max-range angle: the required %.3f km is inside the\n', ...
                rngReq./1000);
        fprintf('                   %.3f km tolerance of the %.3f km maximum. The two arcs\n', ...
                tolRngM./1000,rngMax./1000);
        fprintf('                   ARE the same arc at this range and the lofted/depressed\n');
        fprintf('                   distinction has no meaning here. The max-range angle is\n');
        fprintf('                   NOT suspect. Tighten tolRangeKm, or take the max-range arc\n');
        fprintf('                   and stop asking which side of itself it is on.\n');
    elseif ~branchOK
        fprintf('  *** CAUTION ***  the flown trajectory MEASURES as the %s branch while the\n', ...
                flownName);
        fprintf('                   solve was run on the %s bracket. The bracket did not hold\n', ...
                pickName);
        fprintf('                   the solution on its own side of the maximum, so the\n');
        fprintf('                   max-range angle above is suspect. Do not read the branch\n');
        fprintf('                   label; read the apogee and the flight time.\n');
    end
%% THE MINIMUM-ENERGY PARAGRAPH, and it is the only place in this summary that
%% compares a flown trajectory against a result computed from the geometry
%% alone. It prints the two solved residuals -- what makes the arc
%% minimum-energy -- and then the two MODELLING residuals against the classical
%% arc, which are physics rather than solver error and are attributed as such:
    if isMinE
    fprintf('  THE MINIMUM-ENERGY SOLVE\n');
    fprintf('    TWO UNKNOWNS, TWO RESIDUALS. The textbook minimum-energy trajectory for a\n');
    fprintf('    %.4f deg range angle wants a burnout flight-path angle of\n',rad2deg(angReqR));
    fprintf('    gamma* = 45 - Lambda/4 = %.4f deg carrying V* = %.1f m/s, and a fixed\n', ...
            rad2deg(gamStarR),refME.V);
    fprintf('    booster burnt to exhaustion cannot deliver that speed -- it OVERSHOOTS it.\n');
    fprintf('    So the burn is CUT SHORT and the loft angle and the cutoff are solved\n');
    fprintf('    together, nested, each level monotonic and each bisected by\n');
    fprintf('    coorbital.util.rangeSolve.\n');
    fprintf('    %-22s %16s %16s\n','','solved','residual');
    fprintf('    %-22s %16.4f %16s\n','loft angle (deg)',rad2deg(me.loftR),'-');
    fprintf('    %-22s %16.6f %16s\n','cutoff fraction (-)',me.cutFrac,'-');
    fprintf('    %-22s %16.4f %16s\n','cutoff time (s)',me.tCutS,'-');
    fprintf('    %-22s %16.4f %16.2e\n','burnout gamma (deg)', ...
            rad2deg(me.gamBoR),rad2deg(me.gamBoR - gamStarR));
    fprintf('    %-22s %16.4f %16.2f\n','achieved range (km)', ...
            me.rngM./1000,me.rngM - rngReq);
    fprintf('      the gamma residual is in DEGREES against a %.4f deg tolerance; the range\n', ...
            tolGamDeg);
    fprintf('      residual is in METRES against a %.1f m one. %d outer steps, %d\n', ...
            tolRngM,me.outerIters,me.nProp);
    fprintf('      propagations.\n');
    fprintf('    AGAINST THE CLASSICAL ARC for this range angle, which is computed from the\n');
    fprintf('    geometry alone and from nothing this script flew:\n');
    fprintf('    %-22s %16s %16s %14s\n','','flown','classical','difference');
    fprintf('    %-22s %16.3f %16.3f %13.2f %%\n','APOGEE (km)', ...
            me.hApoM./1000,refME.hApoM./1000, ...
            100.*(me.hApoM - refME.hApoM)./refME.hApoM);
    fprintf('    %-22s %16.4f %16.4f %13.2f %%\n','FLIGHT TIME (min)', ...
            me.tFlyS./60,refME.tofS./60, ...
            100.*(me.tFlyS - refME.tofS)./refME.tofS);
    fprintf('    %-22s %16.4f %16.4f %13.2e\n','burnout gamma (deg)', ...
            rad2deg(me.gamBoR),rad2deg(refME.gamR),rad2deg(me.gamBoR - refME.gamR));
    fprintf('    %-22s %16.1f %16.1f %13.2f %%\n','burnout speed (m/s)', ...
            me.vBoM,refME.V,100.*(me.vBoM - refME.V)./refME.V);
    fprintf('    THE APOGEE AND TIME DIFFERENCES ARE PHYSICS, NOT SOLVER ERROR. The\n');
    fprintf('    classical result assumes an IMPULSIVE burn at the impact radius in a\n');
    fprintf('    VACUUM. This flight burns for %.3f s and arrives at gamma* %.2f km up at\n', ...
            me.tCutS,me.hBoM./1000);
    fprintf('    %.1f m/s, not at zero altitude at %.1f m/s, and a Keplerian arc from THAT\n', ...
            me.vBoM,refME.V);
    fprintf('    burnout state apogees at %.3f km -- which is the flown figure to %.3f km.\n', ...
            me.hApoKepM./1000,abs(me.hApoM - me.hApoKepM)./1000);
    fprintf('    Coast drag supplies the rest. Expect a few per cent from any finite boost.\n');
%% The contrast against the two full-burn arcs needs BOTH of them, and a branch
%% that does not reach the target supplies NaN for its energy. Printed unguarded
%% that reads as a measurement -- "NaN and -44.9073 MJ/kg, NaN %% apart" -- so
%% the one-branch case gets a paragraph that only claims what it has:
    if bothHere
    fprintf('    WHAT THE TWO FULL-BURN ARCS WOULD HAVE DONE INSTEAD, for contrast: they\n');
    fprintf('    leave burnout at %.4f and %.4f MJ/kg, %.2f %% apart, because the\n', ...
            dep.epsBo./1e6,lof.epsBo./1e6, ...
            100.*abs(lof.epsBo - dep.epsBo)./abs(dep.epsBo));
    fprintf('    propellant load is fixed and both burn all of it; this one leaves at\n');
    fprintf('    %.4f MJ/kg, having thrown %.1f kg of propellant away unburned. That gap\n', ...
            me.epsBo./1e6,me.propLeftKg);
    fprintf('    is the whole reason neither full-burn arc is the trajectory above.\n');
    else
    fprintf('    ONLY ONE FULL-BURN ARC REACHES THIS TARGET, so there is no pair of burnout\n');
    fprintf('    energies to contrast this one against; the %s arc alone leaves at\n', ...
            onlyBranch(dep.exists));
    fprintf('    %.4f MJ/kg against this trajectory''s %.4f MJ/kg, which threw %.1f kg of\n', ...
            onlyEps(dep,lof),me.epsBo./1e6,me.propLeftKg);
    fprintf('    propellant away unburned. THE MINIMUM-ENERGY SOLVE IS UNAFFECTED by the\n');
    fprintf('    missing branch: it uses the loft bracket''s own end there, and the arc it\n');
    fprintf('    found is the same one the geometry asks for either way.\n');
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
    if isMinE
    fprintf('    time             %10.3f %-5s  (COMMANDED THRUST TERMINATION at %.6f of the\n', ...
            traj.t(kBO),'s',me.cutFrac);
    fprintf('                                        %.4f s full burn. It is the SECOND control\n',tBurn);
    fprintf('                                        of this mode: the loft angle sets gamma*\n');
    fprintf('                                        and the cutoff sets the range)\n');
    else
    fprintf('    time             %10.3f %-5s  (propellant exhaustion; the boost is never cut\n', ...
            traj.t(kBO),'s');
    fprintf('                                        short on this branch, the LOFT ANGLE being\n');
    fprintf('                                        the ranging control instead of the cutoff)\n');
    end
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
    if isMinE
    fprintf('    mass at cutoff   %10.1f %-5s  (payload + booster + %.1f kg of UNBURNED\n', ...
            massS(kBO),'kg',me.propLeftKg);
    fprintf('                                        propellant, %.2f %% of the load, thrown\n', ...
            100.*me.propLeftKg./bst.massProp);
    fprintf('                                        away with the booster)\n');
    else
    fprintf('    mass at burnout  %10.1f %-5s  (payload + spent booster)\n',massS(kBO),'kg');
    end
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
            mr.nEval + 1 + dep.nProp + lof.nProp + meProp,'');
    fprintf('                                        %d bracketing the maximum, 1 re-flying it,\n', ...
            mr.nEval);
    if isMinE
    fprintf('                                        %d on the depressed branch, %d on the lofted,\n', ...
            dep.nProp,lof.nProp);
    fprintf('                                        %d on the minimum-energy solve)\n',meProp);
    else
    fprintf('                                        %d on the depressed branch, %d on the lofted)\n', ...
            dep.nProp,lof.nProp);
    end
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
    fprintf('    1. THE AZIMUTH IS EXACT ONLY FOR A NON-ROTATING EARTH, which is why earthSpin\n');
    fprintf('       true is declined at the top of this script. This run has Earth rotation %s,\n',spinTxt);
    fprintf('       env.omegaE = %.6e rad/s, so the ground track of an unbanked\n',env.omegaE);
    fprintf('       trajectory is a great circle and the initial bearing of the\n');
    fprintf('       launch-to-target arc is the WHOLE answer -- one closed-form call, no\n');
    fprintf('       iteration. Turn the Earth on and it stops being the answer, which is what\n');
    fprintf('       the refusal is for.\n');
    if bothHere
    fprintf('       The two branches show how much it would matter: the ground under the\n');
    fprintf('       target sweeps %.0f km east during the depressed arc''s %.0f s and %.0f km\n', ...
            driftDep./1000,dep.tFlyS,driftLof./1000);
    fprintf('       during the lofted arc''s %.0f s, so the aim point is not even the same for\n', ...
            lof.tFlyS);
    fprintf('       the two arcs and the outer azimuth iteration this script lacks would have\n');
    fprintf('       to run PER BRANCH. Rotating-Earth targeting is out of scope for it.\n');
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
    elseif crossWarn
    fprintf('       Cross-range came out at %.2e m WITH THE BANK ANGLE AT ZERO, which the\n', ...
            abs(xTrackM));
    fprintf('       WARNING above already declines to blame on a bank. The sentence normally\n');
    fprintf('       printed here -- that an unbanked track over a non-rotating sphere never\n');
    fprintf('       leaves the great circle it departed on -- is CONTRADICTED by that\n');
    fprintf('       measurement, so it is not printed. Read the warning instead.\n');
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
    fprintf('       path. The clamp also decides WHERE the two branches lie. At\n');
    fprintf('       BM/run_ballistic''s 6 deg the hump still exists, but at a commanded\n');
    fprintf('       -42.902 deg, 2.902 deg BELOW the shipped loftMin of -40 deg, so the\n');
    fprintf('       shipped bracket finds its maximum on an endpoint and is refused for\n');
    fprintf('       BRACKET WIDTH rather than for the clamp. Widening it recovers the hump\n');
    fprintf('       but not a useful second branch: at 6 deg the depressed arc bottoms out at\n');
    fprintf('       4708.463 km and cannot reach this %.3f km target at any loft angle.\n', ...
            rngReq./1000);
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
    info.atMaxRange = atMaxRange;
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

%% The minimum-energy record and the classical arc it was solved against, always
%% present so a caller does not have to test the mode before reading them. On a
%% 'lofted' or 'depressed' run minEnergy.solved is false and every other field
%% of it is absent -- there was no solve -- while classical is populated
%% regardless, because it depends on the GEOMETRY and not on what was flown:
    info.minEnergy  = me;
    info.classical  = refME;
     info.gamStarR  = gamStarR;
    info.meGamResR  = NaN;
    info.meApoRelE  = NaN;
    info.meTofRelE  = NaN;
    if isMinE
    info.meGamResR  = me.gamBoR - gamStarR;
    info.meApoRelE  = (me.hApoM - refME.hApoM)./refME.hApoM;
    info.meTofRelE  = (me.tFlyS - refME.tofS)./refME.tofS;
    end
%% The SAME sum the summary prints above, from the SAME fields. nProp counts
%% whole trajectories flown, and a branch costs its solver evaluations plus the
%% one propagation that re-creates the state history at the converged loft
%% angle, which is br.nProp and not br.nEval:
        info.nProp  = mr.nEval + 1 + dep.nProp + lof.nProp + meProp;

%% Vertical exaggeration for the globe. TRUE SCALE IS THE SHIPPED DEFAULT: the
%% movie carries a true-scale altitude inset, which is what the globe used to
%% have to lie for, and coorbital.viz drops the "(altitude exaggerated Nx)"
%% caption clause at unity so a true-scale picture makes no claim it is not
%% keeping. The ADAPTIVE rule is still one word away -- altExag = 'auto' --
%% because it is the right answer when the two branches are being compared:
%% they differ by a factor of ten in apogee, 208 km depressed against 2118 km
%% lofted, and 'auto' draws them at 9x and 2x where one fixed factor suits
%% neither. Same rule, cap and floor as HGV/run_target, whose exagFor header
%% carries the full rationale:
            hPeakM = max(traj.x(:,1)) - c.rE;
    if ischar(altExag) || isstring(altExag)
           altExag = exagFor(hPeakM,c.rE);
    end
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

function ph = buildPhases(loftR,cutFrac,cfg)
%% Purpose:
%
%  Assemble the three-phase boost-coast-descent chain for ONE loft angle and ONE
%  cutoff fraction. Factored out because the max-range search, both branch
%  solves, the minimum-energy solve and the final re-fly all need the same
%  chain, and two copies of it would be two chains that could drift apart.
%
%  THE LOFT ANGLE ENTERS THROUGH THE PITCH SCHEDULE AND NOWHERE ELSE. The
%  schedule is the reference profile's SHAPE rescaled so that its first node
%  stays at cfg.thetaTop and its last node becomes loftR; nothing else about the
%  chain depends on it. A change of loft angle changes WHERE the vehicle is
%  going when the motor quits, never how much propellant it burned.
%
%  THE CUTOFF FRACTION IS THE tspan, NOT AN EVENT, which is HGV/run_target's
%  rule and needs no new machinery: phase 1 keeps its burnout event, so a burn
%  that is not cut short still ends on propellant exhaustion. cutFrac = 1 is the
%  SENTINEL for "not cut short" and reproduces the full-burn chain exactly --
%  the horizon tspan and the jettison-the-dry-structure link -- which is what
%  keeps the 'lofted' and 'depressed' modes bit-for-bit what they were.
%
%  BELOW 1 THE WHOLE BOOSTER GOES OVERBOARD, unburned propellant included, and
%  the link WRITES cfg.mCoast rather than subtracting the dry mass. It has to:
%  the coast vehicle was built around cfg.mCoast before the cutoff was known,
%  and coorbital.eom.massConstant checks the state mass against it on every
%  derivative evaluation. That is also why the caller refuses
%  minimum-energy with separation = false.
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
%  cutFrac          [1 x 1]                     Thrust-termination time as a
%                                               fraction of the full burn (-),
%                                               in (0 .. 1]. 1 commands no
%                                               cutoff at all
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

%% Phase 1, boost. Not cut short: burns to propellant exhaustion, then jettisons
%% the spent structure. With separation switched off the link is [], which
%% phaseRun reads as identity. CUT SHORT: the tspan ends it at the commanded
%% fraction of the burn and the link writes the coast mass outright, so the
%% propellant still aboard leaves with the booster:
         ph(1).eom = @coorbital.eom.boost3DOF;
       ph(1).guide = @(t,x) coorbital.guide.pitchProgram(t,x,schedB);
   ph(1).terminate = @(t,x) coorbital.prop.eventBurnout(t,x,cfg.mBurnout);
       ph(1).tspan = [0 cfg.tMaxBoost];
        ph(1).link = [];
         ph(1).veh = cfg.bst;
    if cutFrac < 1
       ph(1).tspan = [0 cutFrac.*cfg.tBurn];
        ph(1).link = @(x) [x(1:6); cfg.mCoast];
    elseif cfg.separation
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

function [rngM,traj] = flyLoft(loftR,cutFrac,cfg,rI)
%% Purpose:
%
%  Fly the whole boost-coast-descent chain at one loft angle and one cutoff
%  fraction and return the great-circle surface range from the launch point to
%  the impact point. This is the function the max-range search maximises, the
%  function coorbital.util.rangeSolve bisects on at BOTH levels of the
%  minimum-energy solve, and the only place any of them touches the physics.
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
%  cutFrac          [1 x 1]                     Thrust-termination time as a
%                                               fraction of the full burn (-);
%                                               1 commands no cutoff
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

                ph = buildPhases(loftR,cutFrac,cfg);
              traj = coorbital.prop.phaseRun(ph,cfg.x0,cfg.bst,cfg.env);

%% All three phases must have run, and the last must have stopped ON the impact
%% altitude rather than at its horizon:
             hEndM = traj.x(end,1) - cfg.rE;
             nPhRn = numel(unique(traj.phaseIdx));
    if nPhRn < 3 || abs(hEndM - cfg.hStop) > 1e-3
        error('coorbital:runBallisticTarget:propagationIncomplete', ...
            ['A loft angle of %.6f deg at a cutoff fraction of %.6f produced ' ...
             '%d of 3 phases and ended at h = %.3f km against a %.3f km stop ' ...
             'altitude. The flight did not complete, so its range is ' ...
             'meaningless to the search. Two usual causes: a loft angle so ' ...
             'depressed that the vehicle burns out already descending and the ' ...
             'apogee event never fires -- raise loftMin; or a cutoff so early ' ...
             'that it never leaves the atmosphere -- raise cutFracMin.'], ...
            rad2deg(loftR),cutFrac,nPhRn,hEndM./1000,cfg.hStop./1000);
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
             'outside THIS BRACKET and range is monotonic across all of it. ' ...
             'There is then only ONE branch here, and the lofted/depressed ' ...
             'structure this script is built around does not exist. Two usual ' ...
             'causes, in order of likelihood: (1) THE BRACKET IS TOO NARROW ' ...
             'and the hump lies outside it -- widen loftMin and loftMax. This ' ...
             'is what happens at BM/run_ballistic''s 6 deg clamp on the ' ...
             'shipped booster: the maximum is at a commanded -42.902 deg, ' ...
             '2.902 deg below the shipped loftMin of -40 deg, and a loftMin of ' ...
             '-140 finds it. (2) The alphaMax clamp is too tight for the ' ...
             'vehicle to be pushed past the max-range attitude at any loft ' ...
             'angle the bracket allows -- raise alphaMax. Note that widening ' ...
             'the bracket recovers the hump but not necessarily a USEFUL ' ...
             'second branch: at 6 deg the depressed arc spans only 4708.463 ' ...
             'to 5211.525 km.'], ...
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

function me = minEnergySolve(rngReq,gamStarR,cfg,rI,c,latTgR,lonTgR,dep,lof, ...
                             loftLoR,loftHiR,tolRngM,tolGamR,cutFracMin)
%% Purpose:
%
%  Solve the TEXTBOOK minimum-energy ballistic trajectory: the loft angle AND
%  the thrust-termination fraction that together put the vehicle on the required
%  range with a burnout flight-path angle of gammaStar = 45 deg - Lambda/4.
%
%  TWO UNKNOWNS AGAINST TWO RESIDUALS, and they are solved NESTED rather than as
%  a 2x2 system, because nesting makes each level monotonic and lets both reuse
%  coorbital.util.rangeSolve exactly as it was designed to be used:
%
%    INNER, on the CUTOFF FRACTION at a held loft angle. Range rises
%           monotonically with the cutoff -- more burn, more energy, more range
%           -- which is the property bisection needs. Bracket [cutFracMin, 1].
%
%    OUTER, on the LOFT ANGLE, driving the burnout gamma of the inner solution
%           to gammaStar. Along the constant-range locus the burnout gamma rises
%           monotonically with the loft angle.
%
%  THE OUTER BRACKET IS THE TWO FULL-BURN BRANCH SOLUTIONS, already solved by
%  the caller and therefore free. Between them and only between them the FULL
%  BURN reaches at least the required range, so the inner solve is guaranteed a
%  cutoff fraction at or below 1 that lands on it. Where a branch does not exist
%  the caller's own loft limit is used at that end instead, which is admissible
%  for precisely the reason the branch is missing: the full burn OVERSHOOTS the
%  target there, so a cut-short burn can still be made to land on it.
%
%  A GAMMA STAR OUTSIDE THE REACHABLE BAND IS NOT AN ERROR. rangeSolve reports
%  it as converged = false with the band, and this returns solved = false with
%  the two end values so the caller can refuse in the problem's own vocabulary.
%  Every OTHER failure -- an inner solve that cannot make its range on a bracket
%  that guarantees it can -- is a broken invariant and throws.
%
%% Inputs:
%
%  rngReq           [1 x 1]                     Required surface range (m)
%
%  gamStarR         [1 x 1]                     Minimum-energy burnout
%                                               flight-path angle (rad)
%
%  cfg              Struct                      Chain configuration; see
%                                               buildPhases and flyLoft
%
%  rI               [1 x 1]                     Impact sphere radius (m)
%
%  c                Struct                      coorbital.util.missileConst
%
%  latTgR, lonTgR   [1 x 1]                     Target coordinates (rad)
%
%  dep, lof         Struct                      The two full-burn branch
%                                               records from solveBranch; only
%                                               exists and loftR are read
%
%  loftLoR, loftHiR [1 x 1]                     The user's loft bracket (rad),
%                                               used at whichever end has no
%                                               branch solution
%
%  tolRngM          [1 x 1]                     Inner tolerance, on the achieved
%                                               range (m)
%
%  tolGamR          [1 x 1]                     Outer tolerance, on the burnout
%                                               flight-path angle (rad)
%
%  cutFracMin       [1 x 1]                     Floor of the inner bracket, as a
%                                               fraction of the full burn (-)
%
%% Outputs:
%
%  me               Struct                      Solve record:
%                                               solved      logical
%                                               why         char, populated on
%                                                           a refusal
%                                               loftAR,loftBR (rad) the outer
%                                                           bracket actually
%                                                           searched
%                                               gamLoR,gamHiR (rad) burnout
%                                                           gamma at its ends
%                                               ALL BELOW ONLY WHEN solved:
%                                               exists      logical, always true
%                                               loftR       (rad)
%                                               cutFrac     (-)
%                                               tCutS       (s)
%                                               rngM, missM (m)
%                                               gamBoR      (rad)
%                                               vBoM        (m/s)
%                                               hBoM        (m)
%                                               hApoM       (m)
%                                               hApoKepM    (m) the KEPLERIAN
%                                                           apogee of the
%                                                           burnout state, for
%                                                           attributing the gap
%                                                           to the classical arc
%                                               tFlyS       (s)
%                                               vImpM       (m/s)
%                                               gamImR      (rad)
%                                               epsBo       (J/kg)
%                                               mCutKg,propLeftKg (kg)
%                                               outerIters, innerIters
%                                               nProp       propagations spent
%                                               solveInfo   the outer rangeSolve
%                                                           record
%                                               traj        struct
%
%% References:
%   [1] Bate, R.R., Mueller, D.D., White, J.E., "Fundamentals of
%       Astrodynamics," Dover, 1971, Ch. 6.
%   [2] Press, W.H., et al., "Numerical Recipes," 3rd ed., Cambridge, 2007,
%       Section 9.1. Bisection, used at both levels.
%
%% Revision History:
%  Michael Casey                                                08/08/2026
%  Copyright 2026 Coorbital, Inc.
%% ------------------------ Begin Code Sequence ---------------------------

%% The outer bracket, from the full-burn branch solutions where they exist and
%% from the user's own loft limits where they do not:
            loftAR = loftLoR;
    if dep.exists
            loftAR = dep.loftR;
    end
            loftBR = loftHiR;
    if lof.exists
            loftBR = lof.loftR;
    end
         me.solved = false;
            me.why = '';
         me.loftAR = loftAR;
         me.loftBR = loftBR;
         me.gamLoR = NaN;
         me.gamHiR = NaN;
      me.innerIters = 0;
      me.outerIters = 0;
          me.nProp = 0;

%% A collapsed bracket means the required range IS the maximum range to within
%% the range tolerance: the two branches have met, there is one loft angle left
%% and no freedom to trade gamma against it. Refused rather than bisected on an
%% interval of zero width, which rangeSolve would throw on:
    if ~(loftAR < loftBR)
            me.why = sprintf(['The two full-burn arcs have MERGED at a loft angle of ' ...
                              '%.4f deg, so there is no loft interval left to search. ' ...
                              'The required range is within the range tolerance of the ' ...
                              'MAXIMUM range, where the lofted and depressed arcs are the ' ...
                              'same arc; tighten tolRangeKm or move the target inside the ' ...
                              'envelope.'],rad2deg(loftAR));
        return;
    end

             nProp = 0;
           innerIt = 0;

%% The OUTER objective: hold a loft angle, solve the cutoff fraction that makes
%% the range, and report the burnout gamma that came out. Nested so the
%% propagation tally cannot drift from the number of flights actually made:
    function gamR = gammaAtRange(loftR)
        [cutF,rngA,svIn] = coorbital.util.rangeSolve(rngReq, ...
                               @(cf) flyCount(loftR,cf),cutFracMin,1,tolRngM);
        if ~svIn.converged
            error('coorbital:runBallisticTarget:minEnergyInnerFailed', ...
                ['The minimum-energy INNER solve could not make %.3f km of range ' ...
                 'at a loft angle of %.6f deg over a cutoff bracket of %.4f to ' ...
                 '1.0 of the full burn, achieving %.3f km. The outer bracket is ' ...
                 'supposed to guarantee this cannot happen -- between the two ' ...
                 'full-burn branch solutions the FULL burn already reaches this ' ...
                 'target -- so either the range-versus-loft curve is not the ' ...
                 'single hump the bracketing assumed, or cutFracMin is above the ' ...
                 'cutoff this loft angle needs. The library said: %s'], ...
                rngReq./1000,rad2deg(loftR),cutFracMin,rngA./1000,svIn.message);
        end
           innerIt = innerIt + svIn.iterations;
              oInn = flyMeasure(loftR,cutF);
              gamR = oInn.gamBoR;
    end

%% The INNER objective, and the only place either level touches the physics:
    function rngM = flyCount(loftR,cutF)
             nProp = nProp + 1;
              rngM = flyLoft(loftR,cutF,cfg,rI);
    end

%% One propagation, kept, and reduced to the quantities the record needs:
    function o = flyMeasure(loftR,cutF)
             nProp = nProp + 1;
        [rngM,trj] = flyLoft(loftR,cutF,cfg,rI);
               kBO = find(trj.phaseIdx == 1,1,'last');
              kApo = find(trj.phaseIdx == 2,1,'last');
            o.traj = trj;
           o.loftR = loftR;
         o.cutFrac = cutF;
           o.tCutS = trj.t(kBO);
            o.rngM = rngM;
          o.gamBoR = trj.x(kBO,5);
            o.vBoM = trj.x(kBO,4);
            o.hBoM = trj.x(kBO,1) - c.rE;
           o.mCutKg = trj.x(kBO,7);
           o.hApoM = trj.x(kApo,1) - c.rE;
           o.tFlyS = trj.t(end);
           o.vImpM = trj.x(end,4);
          o.gamImR = trj.x(end,5);
           o.epsBo = trj.x(kBO,4).^2./2 - c.muE./trj.x(kBO,1);
           o.missM = rI.*coorbital.util.greatCircle(trj.x(end,3),trj.x(end,2), ...
                                                    latTgR,lonTgR);
    end

%% The OUTER solve. rangeSolve is generic in its abscissa and its ordinate --
%% nothing in it knows that the parameter is an angle and the value another
%% angle -- so the same bisection that solves range against loft solves gamma
%% against loft here, with no second implementation to keep in step:
    [loftS,gamAch,svOut] = coorbital.util.rangeSolve(gamStarR,@gammaAtRange, ...
                                                     loftAR,loftBR,tolGamR);
          me.gamLoR = svOut.fLo;
          me.gamHiR = svOut.fHi;
    me.outerIters   = svOut.iterations;
      me.innerIters = innerIt;
       me.solveInfo = svOut;
          me.nProp  = nProp;

%% gammaStar outside the reachable band is a refusal, not an error. It means the
%% minimum-energy geometry for this range angle is not something this vehicle
%% can be flown into at this range, whatever the cutoff:
    if ~svOut.converged
            me.why = sprintf(['gamma* = %.4f deg is not reachable. Over the loft ' ...
                              'bracket %.4f to %.4f deg -- each end already solved ' ...
                              'for the cutoff fraction that makes the range -- the ' ...
                              'burnout flight-path angle spans only %.4f to %.4f ' ...
                              'deg, and the minimum-energy condition lies outside ' ...
                              'it. The library said: %s'], ...
                             rad2deg(gamStarR),rad2deg(loftAR),rad2deg(loftBR), ...
                             rad2deg(min(svOut.fMin,svOut.fMax)), ...
                             rad2deg(max(svOut.fMin,svOut.fMax)),svOut.message);
        return;
    end

%% Converged. Re-solve the cutoff at the settled loft angle and keep that
%% trajectory, which is one propagation rather than carrying a cache through the
%% bisection; the assert is that the kept flight is the one the outer solve was
%% told about, which a stale cache would fail:
    [cutS,~,svFin] = coorbital.util.rangeSolve(rngReq,@(cf) flyCount(loftS,cf), ...
                                               cutFracMin,1,tolRngM);
    assert(svFin.converged, ...
        ['the minimum-energy inner solve converged during the outer bisection ' ...
         'and then failed when repeated at the settled %.6f deg loft angle; ' ...
         'the propagation is not repeatable'],rad2deg(loftS));
                me = mergeInto(me,flyMeasure(loftS,cutS));
    assert(abs(me.gamBoR - gamAch) < 1e-9, ...
        ['the re-flown minimum-energy trajectory has a burnout gamma of %.12f ' ...
         'rad against the %.12f rad the outer solve reported; the two levels ' ...
         'have come out of step'],me.gamBoR,gamAch);
      me.innerIters = innerIt + svFin.iterations;
          me.nProp  = nProp;
         me.solved  = true;
         me.exists  = true;
      me.propLeftKg = me.mCutKg - cfg.mBurnout;
       me.hApoKepM  = keplerApogee(c.rE + me.hBoM,me.vBoM,me.gamBoR,c) - c.rE;
end

function s = mergeInto(s,add)
%% Purpose:
%
%  Copy every field of one struct onto another, overwriting collisions. Exists
%  so the minimum-energy record can be filled from the measured-flight struct in
%  one line instead of twenty, and so a field added to the measurement cannot be
%  forgotten in the copy.
%
%% Inputs:
%
%  s                Struct                      Destination
%
%  add              Struct                      Source; its fields win
%
%% Outputs:
%
%  s                Struct                      Destination with the source's
%                                               fields written onto it
%
%% Revision History:
%  Michael Casey                                                08/08/2026
%  Copyright 2026 Coorbital, Inc.
%% ------------------------ Begin Code Sequence ---------------------------

                fn = fieldnames(add);
    for kf = 1:numel(fn)
          s.(fn{kf}) = add.(fn{kf});
    end
end

function rApoM = keplerApogee(rM,vM,gamR,c)
%% Purpose:
%
%  Apogee radius of the Keplerian arc through one state. Used to attribute the
%  gap between the flown minimum-energy apogee and the classical one: the
%  classical result assumes an impulsive burn at the impact radius, and most of
%  the difference is simply that the real burnout happens tens of kilometres up.
%  What this leaves unexplained is drag over the coast, and that is the point of
%  computing it.
%
%% Inputs:
%
%  rM               [1 x 1]                     Radius at the state (m)
%
%  vM               [1 x 1]                     Speed there (m/s)
%
%  gamR             [1 x 1]                     Flight-path angle there (rad)
%
%  c                Struct                      coorbital.util.missileConst
%
%% Outputs:
%
%  rApoM            [1 x 1]                     Apogee radius (m), NaN if the
%                                               state is not bound
%
%% References:
%   [1] Bate, R.R., Mueller, D.D., White, J.E., "Fundamentals of
%       Astrodynamics," Dover, 1971, Ch. 1.
%
%% Revision History:
%  Michael Casey                                                08/08/2026
%  Copyright 2026 Coorbital, Inc.
%% ------------------------ Begin Code Sequence ---------------------------

               eps = vM.^2./2 - c.muE./rM;
             rApoM = NaN;
    if eps < 0
                aM = -c.muE./(2.*eps);
                hM = rM.*vM.*cos(gamR);
                pM = hM.^2./c.muE;
               ecc = sqrt(max(0,1 - pM./aM));
             rApoM = aM.*(1 + ecc);
    end
end

function ref = minEnergyRef(psiR,c)
%% Purpose:
%
%  The CLASSICAL minimum-energy ballistic trajectory between two points at the
%  same radius, in closed form, for a free-flight range angle psi. It is
%  computed from the geometry alone -- no propagation, no vehicle, no
%  atmosphere -- and it is what the minimum-energy mode is solved against and
%  reported against.
%
%  THE ASSUMPTIONS ARE STRONG AND THEY MATTER when the numbers are compared:
%  an IMPULSIVE burn AT the impact radius, a VACUUM, a SPHERICAL non-rotating
%  Earth. A real boost takes tens of seconds and finishes tens of kilometres up,
%  which is worth a few per cent of apogee and of flight time. The caller
%  attributes that gap rather than absorbing it.
%
%      V*^2      = (mu/rE) * 2 sin(psi/2) / (1 + sin(psi/2))
%      gammaStar = 45 deg - psi/4
%
%  and the rest is the two-body orbit those two determine. The 45 and the 4 are
%  constants of the closed-form result in reference [1], not physical constants
%  of the Earth, so they do not belong in coorbital.util.missileConst.
%
%% Inputs:
%
%  psiR             [1 x 1]                     Free-flight range angle (rad),
%                                               strictly between 0 and pi
%
%  c                Struct                      coorbital.util.missileConst;
%                                               muE (m^3/s^2) and rE (m) are
%                                               read
%
%% Outputs:
%
%  ref              Struct                      Closed-form reference:
%                                               psiR   [1 x 1] (rad) as passed
%                                               V      [1 x 1] (m/s) burnout
%                                                      speed
%                                               gamR   [1 x 1] (rad) burnout
%                                                      flight-path angle
%                                               aM     [1 x 1] (m) semi-major
%                                                      axis
%                                               ecc    [1 x 1] (-)
%                                               hApoM  [1 x 1] (m) apogee
%                                                      ALTITUDE above rE
%                                               tofS   [1 x 1] (s) burnout to
%                                                      impact
%
%% References:
%   [1] Bate, R.R., Mueller, D.D., White, J.E., "Fundamentals of
%       Astrodynamics," Dover, 1971, Ch. 6. The minimum-energy free-flight
%       trajectory and gammaStar = 45 deg - psi/4.
%
%% Revision History:
%  Michael Casey                                                08/08/2026
%  Copyright 2026 Coorbital, Inc.
%% ------------------------ Begin Code Sequence ---------------------------

    assert(isscalar(psiR) && isfinite(psiR) && psiR > 0 && psiR < pi, ...
        ['the free-flight range angle must lie strictly between 0 and pi rad; ' ...
         'got %s. A zero angle has no trajectory and pi is the antipodal ' ...
         'degeneracy every great circle satisfies.'],mat2str(psiR));

%% Burnout speed and flight-path angle, the two closed-form results:
              sHalf = sin(psiR./2);
                 V2 = (c.muE./c.rE).*2.*sHalf./(1 + sHalf);
             ref.V  = sqrt(V2);
           ref.gamR = pi./4 - psiR./4;
           ref.psiR = psiR;

%% The orbit those two determine, taken at r = rE because that is the radius the
%% closed form is stated at:
                eps = V2./2 - c.muE./c.rE;
    assert(eps < 0, ...
        ['the minimum-energy burnout speed for a %.6f rad range angle is not ' ...
         'bound to the Earth; the closed form has been applied outside its ' ...
         'range.'],psiR);
             ref.aM = -c.muE./(2.*eps);
                 hM = c.rE.*ref.V.*cos(ref.gamR);
                 pM = hM.^2./c.muE;
            ref.ecc = sqrt(max(0,1 - pM./ref.aM));
          ref.hApoM = ref.aM.*(1 + ref.ecc) - c.rE;

%% Time of flight, burnout to impact, by Kepler's equation. The arc is symmetric
%% about apogee, so it is twice the time from burnout up to apogee, and the
%% burnout point is on the ASCENDING side -- true anomaly in (0,pi):
                 nu = acos((pM./c.rE - 1)./ref.ecc);
                 EA = 2.*atan(sqrt((1 - ref.ecc)./(1 + ref.ecc)).*tan(nu./2));
                 MA = EA - ref.ecc.*sin(EA);
              nMean = sqrt(c.muE./ref.aM.^3);
           ref.tofS = 2.*(pi - MA)./nMean;
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

function s = onlyBranch(depExists)
%% Purpose:
%
%  Name whichever of the two full-burn arcs is the one that exists, for the
%  one-branch paragraph of the minimum-energy summary. Exactly one of them does
%  when that paragraph is reached: neither existing is refused above.
%
%% Inputs:
%
%  depExists        [1 x 1] logical             Whether the depressed arc
%                                               reaches the target
%
%% Outputs:
%
%  s                Char [1 x n]                'depressed' or 'lofted'
%
%% Revision History:
%  Michael Casey                                                08/08/2026
%  Copyright 2026 Coorbital, Inc.
%% ------------------------ Begin Code Sequence ---------------------------

                 s = 'lofted';
    if depExists
                 s = 'depressed';
    end
end

function e = onlyEps(dep,lof)
%% Purpose:
%
%  Burnout specific energy of whichever full-burn arc exists, in MJ/kg, for the
%  same paragraph. Reading the missing branch's field would print the NaN this
%  exists to keep out of a results sentence.
%
%% Inputs:
%
%  dep, lof         Struct                      The two branch records; exists
%                                               and epsBo are read
%
%% Outputs:
%
%  e                [1 x 1]                     Burnout specific energy of the
%                                               surviving arc (MJ/kg)
%
%% Revision History:
%  Michael Casey                                                08/08/2026
%  Copyright 2026 Coorbital, Inc.
%% ------------------------ Begin Code Sequence ---------------------------

                 e = lof.epsBo./1e6;
    if dep.exists
                 e = dep.epsBo./1e6;
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

function [why,ok] = whyBurnout(tEnd,mEnd,mBurnoutT,tMaxBoost,tCut)
%% Purpose:
%
%  Diagnose why the boost phase stopped. There are TWO nominal terminations and
%  which of them applies depends on the mode:
%
%    PROPELLANT EXHAUSTION, which coorbital.prop.eventBurnout detects as the
%    mass state descending through the total burnout mass. This is the only
%    nominal end for 'lofted' and 'depressed', which range on the loft angle and
%    never cut the burn short.
%
%    A COMMANDED THRUST TERMINATION at tCut, which is how 'minimum-energy'
%    reaches a burnout speed below what the full burn delivers. Recognised only
%    when a cutoff was actually commanded; tCut is NaN otherwise and the clause
%    cannot fire, every comparison against NaN being false.
%
%  Anything else -- a horizon timeout above all -- is a fault, and must be named
%  rather than left to look like a completed burn.
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
%  tCut             [1 x 1]                     Commanded thrust-termination
%                                               time (s), or NaN when the burn
%                                               was not cut short
%
%% Outputs:
%
%  why              Char [1 x n]                Human-readable reason
%
%  ok               [1 x 1] logical             True for either nominal end
%
%% Revision History:
%  Michael Casey                                                08/07/2026
%  Michael Casey  Recognise a commanded cutoff                   08/08/2026
%  Copyright 2026 Coorbital, Inc.
%% ------------------------ Begin Code Sequence ---------------------------

                ok = false;
    if abs(mEnd - mBurnoutT) < 1e-6
               why = sprintf('propellant exhausted at t = %.3f s, m = %.1f kg (nominal)', ...
                             tEnd,mEnd);
                ok = true;
    elseif abs(tEnd - tCut) < 1e-6
               why = sprintf(['thrust terminated on command at t = %.4f s with %.1f kg ' ...
                              'of propellant still aboard (nominal for ' ...
                              'minimum-energy)'],tEnd,mEnd - mBurnoutT);
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
%  floored at 2. IT IS NOT THE SHIPPED DEFAULT ANY MORE -- altExag ships at 1,
%  true scale -- but it is reached by altExag = 'auto', and that is worth asking
%  for when the two branches are compared: they differ by a factor of ten in
%  apogee, 208 km depressed against 2118 km lofted on the shipped case, so no
%  single fixed factor suits both. At the 3x BM/run_ballistic hard-codes, the
%  depressed arc is invisible. The full rationale for the cap, the floor and the
%  0.3 rE target, and for where the invariant stops holding, is in the exagFor
%  header of HGV/run_target.m; this is the same rule.
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
