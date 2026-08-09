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
%% Note -- the targeting problem is TWO controls and TWO residuals:
%
%  BM/run_ballistic flies launch site + azimuth + pitch program and lands
%  wherever the physics puts it. This script inverts that. Two unknowns and two
%  misses, and on anything but a non-rotating unbanked run they COUPLE:
%
%    LAUNCH AZIMUTH psi.    Rotating the launch heading rotates the whole ground
%                           track about the pad, so it moves the impact point
%                           mostly CROSSWISE to the intended course.
%
%    A RANGE CONTROL.       Which one depends on the mode, and the next two
%                           notes are about that: the LOFT ANGLE at the full
%                           burn for 'lofted' and 'depressed', the CUTOFF
%                           FRACTION for 'minimum-energy'. Either way it moves
%                           the impact point mostly ALONG the course. This is
%                           where a ballistic missile differs from the
%                           boost-glide vehicle of HGV/run_target, which ranges
%                           on thrust termination in every mode, and it is the
%                           whole reason this script exists separately.
%
%  THE TWO RESIDUALS ARE THE COMPONENTS OF THE MISS VECTOR, target to impact,
%  resolved at the TARGET along the intended course: down-range, positive long,
%  and cross-range, positive right of course. Both zero is exactly "the impact
%  point IS the target", and it refers to no great circle at all -- which
%  matters, because once the Earth turns or a bank is commanded the flown track
%  is not one. Writing the pair instead as a great-circle range residual beside
%  a cross-track offset would be the same information in a worse basis: it would
%  keep quoting an arc the vehicle no longer flies, and its Jacobian would be
%  further off diagonal.
%
%  SOLVED IN TWO STAGES, and the first is what makes the second safe:
%
%    STAGE 1, A SEED, AT     coorbital.util.greatCircleBearing returns the
%    ONE FIXED AZIMUTH.      initial bearing of the launch-to-target great
%                            circle, clockwise from north, which is exactly this
%                            library's heading convention psi. It REFUSES
%                            coincident, antipodal and polar geometries rather
%                            than returning a plausible azimuth for them, and
%                            the refusal is caught below and re-phrased against
%                            the USER PARAMETERS entries that caused it -- which
%                            still belongs at the seed, because a geometry with
%                            no bearing has no starting point and a Newton is no
%                            remedy for that. EVERYTHING ELSE IN STAGE 1 IS
%                            FLOWN AT THAT ONE AZIMUTH: the max-range hump, the
%                            certified maximiser interval, the two monotone
%                            branch brackets, both branch solves and the
%                            minimum-energy minimisation. That is deliberate and
%                            it is what preserves them. Range against loft is a
%                            CURVE at a fixed azimuth and a SURFACE otherwise,
%                            and a golden section on a surface certifies
%                            nothing.
%
%    STAGE 2, THE SOLVE.     coorbital.util.aimSolve runs a damped Newton on
%                            BOTH controls at once, over a finite-difference
%                            Jacobian, until both components of the miss are
%                            inside tolerance. On a non-rotating unbanked run
%                            the seed is usually already inside it, so this
%                            stage converges at its first evaluation and hands
%                            the stage-1 answer back unmoved.
%
%  WHAT STAGE 1 MEASURES IS THEREFORE APPROXIMATE FOR THE AZIMUTH STAGE 2
%  SETTLES ON, and this file says so in four places rather than one: here, in
%  the reachable line of the summary, in the aim-solve refusal banner and in
%  item 3 of the summary's closing block. The reachable band, the maximum range
%  and the too-far/too-close refusals are all quantities of the SEED azimuth.
%  The alternative -- re-bracketing the hump inside every Newton evaluation --
%  would multiply an already 700-propagation mode by the cost of a whole
%  bracketing, and would still not certify anything, because the certification
%  is a one-dimensional argument.
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
%  flies the one the branch selector asks for.
%
%  THE BRANCH BRACKETS ARE CERTIFIED, NOT SPLIT AT A MIDPOINT. A golden-section
%  search does not return a maximiser; it returns an interval [aL,bL] that is
%  PROVED to contain one. Splitting the loft axis at the midpoint of that
%  interval gives two brackets of which one crosses the true maximum, so range
%  is not monotone on it and bisection may pick the wrong root or reject a
%  reachable target. The two branch solves therefore run on the CERTIFIED
%  one-sided intervals
%
%      depressed   [loftMin, aL]        lofted   [bL, loftMax]
%
%  on which range provably is monotone -- see the MONOTONICITY IS ASSUMED, NOT
%  CHECKED note in coorbital.util.rangeSolve. The unresolved band between aL and
%  bL is not assigned to either branch: a required range within tolRangeKm of
%  min(R(aL),R(bL)) is reported as COALESCED, because at that range the two arcs
%  are the same arc and the distinction has no content.
%
%  THE BRANCH IS MEASURED, NOT ASSUMED, AND IT IS MEASURED FROM THE ROOT'S
%  POSITION. The flown loft angle is compared against the certified maximiser
%  interval: below aL is depressed, above bL is lofted, inside is coalesced.
%  That is the only invariant available. It used to be taken from the flown
%  apogee and flight time against the max-range arc's, and those are NOT branch
%  invariants for a finite powered atmospheric arc -- drag, lift, burnout
%  altitude and boost duration can make either non-monotone in the loft angle,
%  and near the maximum both differences vanish quadratically and their signs
%  are set by search and integration error rather than by the branch. Apogee and
%  flight time are still reported, as the descriptive quantities they are, and a
%  disagreement between them is a printed CAUTION rather than a classification.
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
%  THE THIRD IS A CONSTRAINED MINIMISATION, and it is a different problem. It
%  solves, explicitly,
%
%      minimise    eps_BO = V_BO^2/2 - mu/r_BO
%      over        the LOFT ANGLE and the CUTOFF FRACTION
%      subject to  R(loft,cutFrac) = the required range.
%
%  Two parameters and one constraint leave a ONE-DIMENSIONAL feasible family, so
%  the problem is well posed: parameterise the family by the loft angle, let the
%  cutoff fraction be whatever makes the range, and minimise the burnout
%  specific energy along it. That is a real objective with a real minimum, which
%  is what makes the name 'minimum-energy' accurate.
%
%  WHERE THE AZIMUTH GOES IN THIS MODE, and why it is not a third search
%  variable. The minimisation above runs at the SEED azimuth and settles a loft
%  angle. Stage 2 then HOLDS that loft angle and solves the azimuth beside the
%  CUTOFF FRACTION, so the two-component miss goes to zero while the range
%  control that enforces it stays the one this mode was built around. The
%  constraint is therefore met exactly at the solved azimuth; what is
%  approximate is the MINIMISATION, whose valley was mapped at the seed one.
%
%  THAT IS A CHOICE, AND HERE IS THE ONE IT WAS MADE AGAINST. Adding the azimuth
%  properly gives three parameters and two constraints -- still a
%  one-dimensional family, still well posed -- solved by putting a two-axis
%  aimSolve on (psi, cutFrac) INSIDE every outer evaluation, where the
%  one-dimensional cutoff bisection sits now. It is the more honest problem and
%  it was rejected on cost and on evidence: every outer evaluation is already
%  about 25 propagations, there are about 27 of them, and wrapping each in a
%  Newton multiplies the mode's 700 propagations by something between one and
%  two. What is bought is a loft angle re-minimised at the solved azimuth, and
%  the valley is FLAT AT ITS BOTTOM -- that is why tolLoftMEDeg is 0.05 deg
%  against an energy resolved to parts in a thousand -- so the answer moves by
%  far less than the cost. The approximation is stated in the summary rather
%  than hidden, and the aim solve reports the azimuth correction it applied so a
%  reader can judge how far the family was moved.
%
%  WHY IT IS NOT SOLVED AS A GAMMA MATCH ANY MORE. The classical closed form,
%  reference [1],
%
%      V*^2      = (mu/rE) * 2 sin(Lambda/2) / (1 + sin(Lambda/2))
%      gammaStar = 45 deg - Lambda/4
%
%  is derived for a free-flight arc whose two endpoints lie at the SAME radius,
%  with Lambda the free-flight central angle between them. This script's burnout
%  is neither: it happens downrange of the pad and tens of kilometres above the
%  impact sphere. Driving the burnout flight-path angle to 45 deg - Lambda/4
%  with Lambda taken as the PAD-TO-TARGET angle therefore solves a residual that
%  does not apply to the arc being flown, and reporting agreement with it
%  verifies the wrong condition. gammaStar is still computed and still printed --
%  it is the right yardstick for how far the flight is from the vacuum
%  equal-radius idealisation -- but it is a DIAGNOSTIC beside the achieved
%  burnout gamma, not the thing being solved.
%
%  NEITHER FULL-BURN ARC IS THE ANSWER. At full burn the booster delivers a fixed
%  delta-V, so the lofted and depressed arcs leave burnout with essentially the
%  same energy and there is no minimisation to do at all. The burnout energy has
%  to be free, and the only way this vehicle can lower it is to CUT THE BURN
%  SHORT -- the same thrust-termination control HGV/run_target ranges on.
%
%  THE SOLVE IS A ONE-DIMENSIONAL MINIMISATION OVER A NESTED FEASIBILITY SOLVE:
%
%      INNER, on the CUTOFF FRACTION at a held loft angle, enforcing the
%             CONSTRAINT. The bracket is found by SAMPLING the cutoff axis and
%             taking a sign change of R(loft,cutFrac) - R_req, and monotonicity
%             is then VERIFIED on the interval selected rather than assumed
%             across [cutFracMin, 1]: more burn does not have to mean more
%             range, because added burn also moves the burnout position,
%             altitude, flight-path angle and losses. Where several cutoff roots
%             exist, the one with the lowest burnout energy is taken. Solved to
%             tolRangeMEKm, which is deliberately TIGHTER than tolRangeKm --
%             see the noise note below.
%
%      OUTER, on the LOFT ANGLE, minimising eps_BO along that feasible family.
%             A coarse scan of nScanME points certifies a single valley -- the
%             two ends of the family are the two full-burn arcs, where the burn
%             is not cut at all and the energy is therefore highest, so a valley
%             between them is what the physics predicts -- and a golden-section
%             search then closes the bracket to tolLoftMEDeg.
%
%  THE OUTER BRACKET IS THE TWO FULL-BURN BRANCH SOLUTIONS, which are solved
%  before it and cost nothing extra. Between them, and only between them, the
%  full burn reaches AT LEAST the required range, so the inner solve is
%  guaranteed a cutoff fraction at or below 1 that lands on it. Where a branch
%  does not exist the corresponding end of the loft bracket is used instead,
%  which is admissible for the same reason: that branch is missing precisely
%  because the full burn OVERSHOOTS the target there.
%
%  THE INNER TOLERANCE IS NOISE ON THE OUTER OBJECTIVE, AND IT IS MEASURED
%  RATHER THAN WAVED AWAY. Every outer evaluation is an inner solve accurate
%  only to tolRangeMEKm, so eps_BO(loft) is sampled with an error whose size
%  nothing about the method predicts. The solve therefore re-runs the inner
%  problem at the settled loft angle with a tolerance ten times tighter and
%  reports the change in eps_BO as me.epsNoiseJkg, beside the DEPTH of the
%  energy valley over the neighbouring feasible points. The claim that the
%  minimum is real is the ratio of those two numbers, printed on every run.
%
%  IT IS SHOWN TO BE A MINIMUM, not merely a converged point. The coarse scan is
%  kept, and the summary prints the burnout energy at the two feasible points
%  either side of the solution: a minimum has to be below both.
%
%  IT REDUCES TO THE CLASSICAL ARC IN THE VACUUM EQUAL-RADIUS LIMIT, which is
%  the sanity check the objective earns. Take the boost to be impulsive at
%  r = rE in a vacuum; the feasible family is then the one-parameter family of
%  Keplerian arcs of central angle Lambda from rE back to rE, and minimising
%  V^2/2 - mu/rE over it returns exactly V* and gammaStar above. That limit is
%  asserted in tests/test_runBallisticTarget.
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
%  this configuration it is a long way different: the alphaMax clamp limits how
%  fast the flight path can be pushed over, so a commanded terminal attitude far
%  below the horizon still produces a climbing burnout. The summary prints BOTH
%  for every branch. Read the loft angle as the control it is, not as a
%  prediction of the burnout state.
%
%% Note -- alphaMax IS A VEHICLE LIMIT AND IS NOT SET HERE:
%
%  The clamp on the magnitude of the angle of attack lives in BM/vehicle_bm.m,
%  as veh.alphaMaxDeg, and BOTH this script and BM/run_ballistic read it from
%  there. It used to be a USER PARAMETERS entry in each of them, 12 deg here and
%  6 deg there, for nominally the same airframe -- and the 12 was chosen because
%  it brought the old demonstration target inside the depressed branch. A
%  control-authority limit is not a targeting degree of freedom, and two
%  different limits made the two scripts' performance non-comparable.
%
%  6 deg IS A PLACEHOLDER AWAITING A QUALIFICATION BASIS. It is the value
%  BM/run_ballistic has always flown and the one that was NOT chosen to make a
%  feature work. What it costs, measured: maximum range 5211.525 km at a
%  max-range commanded attitude of about -42.9 deg, and a depressed branch
%  spanning roughly 4708 to 5212 km. Raising the clamp to 12 deg gives the
%  depressed branch reach down to about 1684 km and PAYS about 156 km of maximum
%  range for it. Neither figure is a reason to move the limit; the qualified
%  number is whatever the airframe is cleared for, and it does not exist yet.
%  Put a number in alphaMax below to run a deliberate sensitivity study.
%
%  THE SHIPPED DEMONSTRATION TARGET IS INSIDE THE 6 deg DEPRESSED BAND, so all
%  three modes fly on the shipped configuration. That constrains it to roughly
%  4708 to 5212 km, and the loft bracket has to reach below -42.9 deg to find
%  the hump at all -- which is why loftMin ships at -140 deg. A commanded
%  terminal attitude that far below the horizon is not a pitch programme anyone
%  would fly; it is what the depressed branch DEGENERATES INTO once the clamp
%  saturates the achievable pitch-over, and it is reported as such rather than
%  dressed up.
%
%% Note -- WHAT THIS SOLUTION IS AND IS NOT, printed in the summary rather than
%% buried:
%
%  THE AZIMUTH IS SOLVED, NOT CLOSED FORM, and that is what lets earthSpin be
%  true. Note what rotation actually does here: the state longitude and latitude
%  are EARTH-FIXED and the speed is planet-relative, so the target does NOT
%  slide east under the vehicle -- it sits still, and the VEHICLE is deflected
%  off the arc it left on by the Coriolis and centrifugal terms, to the right in
%  the northern hemisphere. Measured on the shipped geometry with earthSpin
%  true, the closed-form bearing alone left the impact point 463.211 km from the
%  target, of which 462.870 km is CROSSWISE and only 17.764 km down-range -- a
%  ratio of 26 to 1, which is the signature of a deflection and not of a moving
%  aim point. The second stage removes it, so the case is flown instead of
%  refused.
%
%  THE SOLVE CONTROLS BOTH AXES, so a banked trajectory is targetable too.
%  Cross-range is driven to zero rather than measured and warned about, and a
%  bank and a rotation are one cross-range residual to it rather than two
%  separate faults.
%
%  A BANK NEEDS loftMin RAISED WITH IT, though, and that is a property of the
%  VEHICLE and not of the solve. Banking tilts the lift vector out of the
%  vertical plane, so less of it is left to push the nose over, and the shipped
%  loftMin of -140 deg -- the degenerate tail the depressed branch collapses
%  into at zero bank -- stops being a commandable attitude at all. Even 0.5 deg
%  of bank puts the first scan point outside what the alphaMax clamp can
%  deliver. The run is REFUSED for that rather than left to throw a library
%  identifier, and the refusal names both entries. At loftMin = -60 deg on the
%  shipped geometry a 2 deg bank flies and misses by 1.59 m with the azimuth
%  aimed off +3.271 deg, and 5 deg by 173.75 m aimed off +8.205 deg -- which is
%  the capability claim above, measured. It follows that the two bank-specific
%  paragraphs in the summary are UNREACHABLE at shipped settings, and reaching
%  them takes the companion change.
%
%  WHAT REMAINS TRUE, and all of it is stated in the summary. Everything stage 1
%  certifies -- the hump, the maximiser interval, the branch brackets, the
%  reachable band, and therefore the too-far and too-close refusals -- is
%  measured at the SEED azimuth and is indicative rather than exact for the
%  solved one. On a rotating run that is not a small caveat: at earthSpin true
%  the shipped geometry's maximum range rises from 5211.5 km to 5439.9 km and
%  the depressed band floor from 4708.5 km to 5085.8 km, so the shipped 4828 km
%  target leaves the depressed branch altogether and 'depressed' is REFUSED
%  there. That is a real physical result and not an artefact -- an easterly
%  launch gains from the planet's rotation -- but it is measured on one azimuth.
%  The Newton is LOCAL: seeded, damped and capped, and when it cannot converge
%  the run is REFUSED with the best point it reached rather than returned as a
%  near miss. And the tolerance is on the MISS ITSELF: each component is held to
%  tolRangeKm/sqrt(2), so the two together cannot exceed the distance the user
%  asked for.
%
%  THE NEWTON IS NOT BOXED TO THE CERTIFIED BRANCH SIDE, and that is a decision
%  rather than an omission. Stage 1's root came from a bisection that cannot
%  leave its bracket; stage 2's Newton has no bracket, so on the two full-burn
%  modes it CAN walk the loft angle towards the maximiser interval and past it.
%  Observed: rotating, 'lofted', a target at 5420 km -- the loft marched 1.33 deg
%  past the interval before the solve refused. A box is the obvious fix and it
%  is the wrong one, because the only interval available to draw it on is the
%  SEED-AZIMUTH one, and the maximum moves with the aim: about 0.28 deg of loft
%  per degree of azimuth on this geometry, against an interval 0.043 deg wide.
%  A box drawn there would exclude genuine roots. What stands instead is
%  DIAGNOSIS rather than prohibition, in three layers: opts.maxStep caps the
%  travel per iteration, branchOfLoft MEASURES the flown angle against the
%  interval and the summary prints a CAUTION naming stage 2 when they disagree,
%  and a Newton that cannot converge refuses outright rather than returning the
%  arc it wandered onto.
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
%                                               angles.
%
%                                               ALWAYS carries refused
%                                               (logical), branchAsked, rngReqM,
%                                               psiSeed and psiLaunch.
%                                               refusedWhy is present ONLY when
%                                               refused is true, so read it
%                                               behind that flag and not on its
%                                               own.
%
%                                               A REFUSAL carries refusedWhy,
%                                               one of 'envelope',
%                                               'minimumEnergy' or 'aimSolve'.
%                                               ALL THREE add loftStarDeg,
%                                               rngMaxM, rngMinM, hApoStarM,
%                                               tFlyStarS, maxRange, classical,
%                                               gamStarR and the per-branch
%                                               records in depressed and lofted,
%                                               so every one of them hands back
%                                               the envelope it was measured
%                                               against -- AT THE SEED AZIMUTH.
%                                               An ENVELOPE refusal adds
%                                               shortfallM; a MINIMUMENERGY
%                                               refusal adds minEnergy, whose
%                                               why says what could not be
%                                               solved; an AIMSOLVE refusal adds
%                                               minEnergy as well plus the
%                                               aimSolve record in aimInfo,
%                                               whose xHist and message say
%                                               where the Newton got to and why
%                                               it stopped, beside the seed and
%                                               best miss components.
%
%                                               A RUN THAT WAS NOT REFUSED
%                                               additionally carries the
%                                               flown trajectory's own numbers,
%                                               the CERTIFIED maximiser interval
%                                               in maxRange.aL and maxRange.bL
%                                               with the ranges maxRange.fAL and
%                                               maxRange.fBL achieved there, the
%                                               classical vacuum equal-radius
%                                               reference in classical, the
%                                               minimum-energy record in
%                                               minEnergy (solved false on the
%                                               two full-burn modes), the three
%                                               DIAGNOSTIC residuals meGamResR,
%                                               meApoRelE and meTofRelE against
%                                               that reference, and the DISPLAY
%                                               choices the figures were drawn
%                                               with -- altExag, its rule as a
%                                               handle, and the two hemisphere
%                                               captions. plotFiles, a cellstr
%                                               of the PNG paths written, is
%                                               present only when the figures
%                                               were drawn AND saved
%
%% References:
%   [1] Bate, R.R., Mueller, D.D., White, J.E., "Fundamentals of
%       Astrodynamics," Dover, 1971, Ch. 6. The free-flight range angle of a
%       Keplerian arc, the two burnout flight-path angles that reach a given
%       range, and the VACUUM EQUAL-RADIUS minimum-energy values V* and
%       gammaStar = 45 deg - Lambda/4, used here as a diagnostic reference.
%   [2] Bowditch, N., "The American Practical Navigator," Pub. No. 9, NGA,
%       chapter on Great-Circle Sailing. The initial course used to SEED the
%       launch azimuth, and the four-parts components the miss vector is built
%       from; see coorbital.util.greatCircleBearing and missVector below.
%   [3] Press, W.H., et al., "Numerical Recipes," 3rd ed., Cambridge, 2007,
%       Sections 9.1 and 10.2. Bisection, used by coorbital.util.rangeSolve,
%       and golden-section search, used below both to bracket the maximum-range
%       loft angle and to minimise the burnout energy along the feasible family.
%   [4] Dennis, J.E. and Schnabel, R.B., "Numerical Methods for Unconstrained
%       Optimization and Nonlinear Equations," SIAM, 1996, Chapters 5 and 6.
%       The finite-difference Jacobian and line-search damped Newton method
%       used by coorbital.util.aimSolve for the two-axis solve.
%   [5] Vinh, N.X., Busemann, A., Culp, R.D., "Hypersonic and Planetary Entry
%       Flight Mechanics," Univ. Michigan Press, 1980.
%
%% Revision History:
%  Michael Casey                                                08/07/2026
%  Michael Casey  minimum-energy solved as an explicit          08/08/2026
%                 constrained minimisation; certified branch
%                 brackets; alphaMax read from the vehicle
%  Michael Casey  Two-axis targeting: the azimuth is SOLVED     08/08/2026
%                 beside the mode's own range control by
%                 coorbital.util.aimSolve, which retires the
%                 rotating-Earth refusal and the downrange-only
%                 cross-track warning
%  Michael Casey  plotFile: the three figures are written to    08/09/2026
%                 disk through coorbital.viz.saveFigure and
%                 their paths printed, gated by showPlots and
%                 an empty stem
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
%% the launch azimuth is SOLVED from these two points, seeded on the
%% great-circle bearing between them, and any azimuth given here could only
%% contradict it:
         latTarget = 62;               %deg, geocentric latitude [-89 .. 89]
         lonTarget = -28;              %deg, longitude [-180 .. 180]
                                       %     The shipped pair is the central United States to
                                       %     the North Atlantic east of Iceland: about 4828 km
                                       %     on a 41 deg azimuth. Deliberately NOT equatorial
                                       %     and NOT due east -- that geometry makes latitude
                                       %     and longitude interchangeable in the range formula
                                       %     and hides a transposition at every call site that
                                       %     uses them.
                                       %     IT IS ALSO INSIDE THE DEPRESSED BAND, which at the
                                       %     vehicle's 6 deg angle-of-attack clamp spans only
                                       %     about 4708 to 5212 km. That is what lets all THREE
                                       %     modes fly on the shipped configuration; a target
                                       %     below 4708 km is reached by the lofted arc and by
                                       %     minimum-energy, and REFUSED for 'depressed'

%% WHICH TRAJECTORY TO FLY. The two FULL-BURN arcs are solved and reported
%% whatever is asked for here, because the trade between them is the point; this
%% chooses what is flown, plotted and returned. Read the THREE MODES note in the
%% header: two of these are full-burn arcs either side of maximum range and the
%% third is the textbook trajectory, which is a different solve on a second
%% control and is NOT a member of the full-burn family:
            branch = 'minimum-energy'; %'minimum-energy' | 'lofted' | 'depressed'
                                       %  minimum-energy: MINIMISES the burnout specific energy
                                       %                  V^2/2 - mu/r over the LOFT ANGLE and
                                       %                  the CUTOFF FRACTION, SUBJECT TO the
                                       %                  achieved range being the required one.
                                       %                  Two parameters and one constraint, so
                                       %                  a one-dimensional feasible family to
                                       %                  minimise along. The burn is CUT SHORT:
                                       %                  a full burn has no energy freedom
                                       %                  left. The classical gammaStar is
                                       %                  REPORTED beside the achieved burnout
                                       %                  gamma as a vacuum equal-radius
                                       %                  yardstick, not solved for. Requires
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
           loftMin = -140;             %deg, most depressed commanded terminal attitude
                                       %     [-200 .. 0]. Negative means the schedule commands
                                       %     the nose BELOW the horizon at end of burn, which
                                       %     the alphaMax clamp then declines to deliver in
                                       %     full. AT THE VEHICLE'S 6 deg CLAMP THE MAX-RANGE
                                       %     ATTITUDE SITS NEAR -42.9 deg, so a bracket that
                                       %     stops at -40 finds its largest range ON an endpoint
                                       %     and the bracketing below refuses for WIDTH. -140 is
                                       %     what it takes to hold the whole depressed branch;
                                       %     read it as the degenerate tail the branch collapses
                                       %     into once the clamp saturates, not as a pitch
                                       %     programme anyone would fly
           loftMax = 85;               %deg, most lofted commanded terminal attitude
                                       %     (loftMin .. 89). 89 is the pad attitude itself,
                                       %     i.e. no pitch-over at all and no downrange travel
        tolLoftDeg = 0.05;             %deg, bracket width the golden-section search closes the
                                       %     max-range angle to [0.005 .. 1]. The hump is flat
                                       %     at its top, so 0.05 deg is far inside the range
                                       %     tolerance below. The search ERRORS rather than
                                       %     returning an unconverged maximiser if 200 steps
                                       %     cannot close it
         nScanLoft = 21;               %-, coarse scan points across the bracket before the
                                       %   golden section starts [7 .. 81]. The scan is what
                                       %   makes the search safe: golden section needs a
                                       %   unimodal bracket and cannot detect that it was not
                                       %   given one, so the scan has to CERTIFY the single hump
                                       %   -- strictly rising to the best node and strictly
                                       %   falling after it. Where it cannot, the scan is
                                       %   REFINED by interleaving midpoints, up to four times,
                                       %   and only a refinement that still cannot certify is an
                                       %   error
        tolRangeKm = 1.0;              %km, convergence tolerance on the MISS DISTANCE, impact
                                       %    point to target [0.05 .. 50]. The two-axis solve is
                                       %    held to this over sqrt(2) on EACH component, so the
                                       %    down-range and cross-range misses together cannot
                                       %    exceed the figure asked for here. The stage-1 branch
                                       %    bisections use the bare figure on range. Each halving
                                       %    costs about one more trajectory propagation per
                                       %    branch, and a propagation is about 0.05 s, so
                                       %    tightening it is cheap

%% MINIMUM-ENERGY MODE ONLY. Ignored by 'lofted' and 'depressed', which fly the
%% full burn and have no cutoff to solve for. See the THREE MODES note:
        cutFracMin = 0.5;              %-, floor of the thrust-termination interval the INNER
                                       %   feasibility solve samples, as a fraction of the full
                                       %   burn [0.2 .. 0.95]. Every sample is propagated, so
                                       %   the floor must be a cutoff the chain can still FLY:
                                       %   too early and the vehicle never reaches apogee, the
                                       %   event never fires, and flyLoft refuses the
                                       %   propagation rather than returning a range from a
                                       %   flight that did not happen
      tolRangeMEKm = 0.05;             %km, convergence tolerance of the INNER feasibility solve
                                       %    [0.001 .. tolRangeKm]. Deliberately TIGHTER than
                                       %    tolRangeKm, because here the inner residual is NOISE
                                       %    on the outer objective rather than the answer: the
                                       %    solve measures the propagated effect on the burnout
                                       %    energy and prints it beside the depth of the energy
                                       %    valley, so the margin is reported and not assumed
           nScanME = 5;                %-, coarse scan points along the feasible family before
                                       %   the energy minimisation starts [3 .. 15]. Each one is
                                       %   a whole inner solve -- about 25 propagations -- so
                                       %   this is the expensive parameter of the mode. It also
                                       %   supplies the evidence that the energy has a SINGLE
                                       %   valley, and the two neighbouring feasible energies
                                       %   the solution is reported against
      tolLoftMEDeg = 0.05;             %deg, bracket width the energy minimisation closes the
                                       %     loft angle to [0.005 .. 2]. The valley is flat at
                                       %     its bottom, so the energy is far better resolved
                                       %     than the angle; each halving costs about 1.5 more
                                       %     inner solves

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
          alphaMax = [];               %deg, clamp on |angle of attack|. EMPTY -- the shipped
                                       %     value -- reads the VEHICLE'S OWN control-authority
                                       %     limit, veh.alphaMaxDeg, which is 6 deg in
                                       %     BM/vehicle_bm.m and is the SAME limit
                                       %     BM/run_ballistic flies. THIS CLAMP, NOT THE
                                       %     SCHEDULE, is what limits how fast the flight path
                                       %     can be pushed over, and it decides where the two
                                       %     branches lie and what they reach. Put a number here
                                       %     [1 .. 15] to OVERRIDE the vehicle for a deliberate
                                       %     sensitivity study -- 12 deg, for instance, buys a
                                       %     depressed branch reaching about 1684 km and pays
                                       %     about 156 km of maximum range for it. The 6 deg is
                                       %     a PLACEHOLDER awaiting a qualification basis, not a
                                       %     cleared value, and it is not a targeting degree of
                                       %     freedom -- see the alphaMax note in the header
         bankAngle = 0;                %deg, bank, all phases. SHIPPED AT ZERO because the
                                       %     simplest case is the one an example script should
                                       %     ship, NOT because the solve cannot handle a bank: a
                                       %     banked segment turns the heading and walks the
                                       %     impact point off the launch-to-target great circle,
                                       %     and the SOLVE NOW SEES THAT -- the launch azimuth is
                                       %     aimed off to compensate.
                                       %     A NON-ZERO VALUE NEEDS loftMin RAISED WITH IT, and
                                       %     the run is REFUSED rather than left to throw a
                                       %     library error if it is not. The bank tilts the lift
                                       %     vector out of the vertical plane, so less of it is
                                       %     left to push the nose over and the shipped
                                       %     loftMin = -140 deg stops being a commandable
                                       %     attitude at all -- even 0.5 deg of bank is enough.
                                       %     Measured, at loftMin = -60 deg on the shipped
                                       %     geometry: 2 deg of bank misses by 1.59 m with the
                                       %     azimuth aimed off +3.271 deg, and 5 deg by 173.75 m
                                       %     aimed off +8.205 deg. Useful range is therefore
                                       %     roughly [-10 .. 10] with a raised loftMin, not the
                                       %     [-90 .. 90] the attitude schedule would accept:
                                       %     this is a near-symmetric re-entry body with almost
                                       %     no lift, and a bank is a large authority on it

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
         earthSpin = false;            %true enables the Coriolis and centrifugal terms. It is
                                       %SUPPORTED: the state is Earth-fixed, so the target does
                                       %not move, but the vehicle is deflected off the arc it
                                       %left on -- 463.211 km of miss on the shipped geometry
                                       %with the closed-form bearing alone, 462.870 km of it
                                       %crosswise -- and the two-axis solve aims off to remove
                                       %it. Shipped false only so the default run is the simplest
                                       %one. IT ALSO MOVES THE ENVELOPE, which stage 1 measures
                                       %at the SEED azimuth: on this north-easterly geometry it
                                       %raises maximum range from 5211.5 to 5439.9 km and the
                                       %depressed band floor from 4708.5 to 5085.8 km, so
                                       %branch = 'depressed' is then REFUSED for the shipped
                                       %4828 km target while the other two modes fly

%% Output:
         showPlots = true;             %false to skip the figures, e.g. under matlab -batch
          plotFile = fullfile(tempdir,'run_ballistic_target'); %char, figure output path STEM,
                                       %WITHOUT an extension: each figure appends its own suffix
                                       %and '.png' -- _profile, _ground_track and _globe. Saving
                                       %happens whenever showPlots is true and this is non-empty,
                                       %so there is no second on/off flag to keep in step; set it
                                       %to '' to draw the figures without writing them. Defaults
                                       %into tempdir for the same reason movieFile does: a
                                       %picture is a build artefact and does not belong in the
                                       %source tree. The names are FIXED, so a re-run overwrites
                                       %its own pictures rather than piling up new ones
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
                      'tolRangeMEKm','nScanME','tolLoftMEDeg','pitchTime', ...
                      'pitchRef','alphaMax','bankAngle','separation','hStop', ...
                      'tMaxBoost','tMaxCoast','tMaxDesc','vehicleFn', ...
                      'boosterFn','atmosFn','gravFn','aeroFn','propFn', ...
                      'earthSpin','showPlots','plotFile','movieOn', ...
                      'movieFrames','movieFile','altExag'};
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
      tolRangeMEKm = overrideOf(opts,'tolRangeMEKm',tolRangeMEKm);
           nScanME = overrideOf(opts,'nScanME',nScanME);
      tolLoftMEDeg = overrideOf(opts,'tolLoftMEDeg',tolLoftMEDeg);
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
          plotFile = overrideOf(opts,'plotFile',plotFile);
           altExag = overrideOf(opts,'altExag',altExag);
          movieOn  = overrideOf(opts,'movieOn',movieOn);
       movieFrames = overrideOf(opts,'movieFrames',movieFrames);
         movieFile = overrideOf(opts,'movieFile',movieFile);

%% The constants, the vehicle and the booster, read HERE rather than after the
%% sanity checks, because the ANGLE-OF-ATTACK CLAMP is a property of the vehicle
%% and the unit conversion below has to convert it:
                 c = coorbital.util.missileConst();
               veh = vehicleFn();
               bst = boosterFn();

%% Resolve the angle-of-attack clamp. It is a CONTROL-AUTHORITY LIMIT of the
%% airframe, not a targeting knob, so its home is BM/vehicle_bm.m and an empty
%% user-block entry reads it from there -- the same limit BM/run_ballistic
%% flies. The two scripts used to carry 12 deg and 6 deg for one airframe, and
%% the 12 had been chosen to bring a demonstration target inside the depressed
%% branch. A number in the block still overrides it, deliberately and visibly:
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
           tolRngM = tolRangeKm.*1000;
         tolRngMEM = tolRangeMEKm.*1000;
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
        tolLoftMER = deg2rad(tolLoftMEDeg);

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
    assert(tolRngMEM > 0,'tolRangeMEKm must be positive.');
    assert(tolRngMEM <= tolRngM, ...
        ['tolRangeMEKm (%.6f km) must not exceed tolRangeKm (%.6f km). It is ' ...
         'the tolerance of the INNER feasibility solve, whose error is NOISE ' ...
         'on the energy the outer minimisation is comparing; a looser inner ' ...
         'tolerance than the outer range tolerance would make the constraint ' ...
         'the loosest thing in the solve.'],tolRangeMEKm,tolRangeKm);
    assert(tolLoftMER > 0,'tolLoftMEDeg must be positive.');
    assert(nScanME >= 3 && nScanME == fix(nScanME), ...
        'nScanME must be a whole number of at least 3; got %s.',mat2str(nScanME));
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
%% The seed azimuth: the closed-form bearing, which is a GUESS
%% -----------------------------------------------------------------
%% The initial bearing of the launch-to-target great circle IS the library's
%% heading state psi -- both are measured clockwise from north -- so it goes
%% straight into the launch state below with no conversion. It is the whole
%% answer only on a non-rotating unbanked run; everywhere else it is the
%% starting point the two-axis solve improves on, and it is the right starting
%% point precisely because it is exact in the easy case and close in the hard
%% ones.
%%
%% coorbital.util.greatCircleBearing REFUSES the geometries that have no
%% azimuth rather than returning the rounding artefact atan2 makes of them, and
%% the refusal arrives here as an error identifier. It is re-phrased against the
%% USER PARAMETERS entries that caused it, because the library function cannot
%% know that the two points it was handed came from latLaunch and latTarget.
%% The translation still belongs here and not in the solve: a geometry with no
%% bearing has no SEED, and a Newton with no starting point is not a remedy for
%% one:
    try
           psiSeed = coorbital.util.greatCircleBearing(latLaunchR,lonLaunchR, ...
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

%% Liftoff state TEMPLATE. Everything but the loft angle and the HEADING is
%% fixed across every propagation this script makes; the pad flight-path angle
%% is the first commanded pitch attitude, so the angle of attack starts at zero.
%% Entry 6 is the heading, and flyLoft overwrites it per trial rather than
%% rebuilding the state, so a trial flight can differ from the seed in the
%% azimuth and in nothing else:
                x0 = [c.rE + hLaunchM; ...
                      lonLaunchR; ...
                      latLaunchR; ...
                      vLaunch; ...
                      pitchRefR(1); ...
                      psiSeed; ...
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
%% branch solves ride on, bound to the SEED AZIMUTH. The cutoff-fraction
%% argument of flyLoft is 1 here, meaning "no commanded cutoff": phase 1 keeps
%% its horizon tspan and ends on propellant exhaustion, exactly as it always
%% has. Only the minimum-energy solve ever passes anything else.
%%
%% BINDING THE AZIMUTH IS WHAT KEEPS STAGE 1 ONE-DIMENSIONAL, and that is not a
%% shortcut, it is what preserves everything the bracketing certifies. Range
%% against loft has one interior hump; range against loft AND azimuth is a
%% surface, and a golden section on a surface certifies nothing. The hump, the
%% certified maximiser interval, the two monotone one-sided brackets and the
%% two branch solves are all statements about one azimuth, and they are made
%% about the seed one:
            fRange = @(loftR) flyLoft(psiSeed,loftR,1,cfg,rI);

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
%% Bracket the MAXIMUM-RANGE loft angle, before either branch is solved
%% -----------------------------------------------------------------
%% This is the step that has no counterpart in HGV/run_target, and everything
%% below depends on it. It returns a CERTIFIED INTERVAL [aL,bL] containing the
%% maximiser -- which is all a golden-section search can honestly claim -- and
%% that interval, not its midpoint, is what splits the loft axis into the two
%% one-sided brackets on which range provably IS monotonic. It also supplies the
%% max-range arc, whose apogee and flight time are reported beside the flown
%% ones as DESCRIPTIVE quantities.
%%
%% THE ONE LIBRARY REFUSAL THAT CAN SURFACE HERE IS TRANSLATED, on the same
%% terms as greatCircleBearing's above. coorbital.guide.pitchProgram refuses an
%% attitude the alphaMax clamp cannot deliver, and a NON-ZERO BANK makes the
%% deepest end of the loft bracket exactly that: the bank tilts the lift vector
%% out of the vertical plane, so less of it is available to push the nose over,
%% and a commanded terminal attitude of loftMin = -140 deg becomes unreachable.
%% The scan evaluates loftLoR FIRST, so this fires on the very first propagation
%% of the run and it fires for a reason the user can act on. Left raw it is an
%% unhandled library identifier out of a script whose bank entry now invites a
%% non-zero value; caught, it is a refusal that names the two entries that have
%% to move together:
    try
[loftStarR,rngMax,mr] = maxRangeLoft(fRange,loftLoR,loftHiR,nScanLoft,tolLoftR);
    catch errLoft
        if ~strcmp(errLoft.identifier,'coorbital:pitchProgram:unreachableAttitude')
            rethrow(errLoft);
        end
        fprintf('\n');
        fprintf('===== Ballistic point-to-point targeting: REFUSED =======\n');
        fprintf('  The loft bracket reaches attitudes this vehicle cannot be commanded to\n');
        fprintf('  at a %.4f deg bank.\n',bankAngle);
        fprintf('\n');
        fprintf('    bank angle       %10.4f %-5s (all phases)\n',bankAngle,'deg');
        fprintf('    loft bracket     %10.4f to %.4f deg   (loftMin is the end that fails:\n', ...
                loftMin,loftMax);
        fprintf('                                      the scan evaluates it first)\n');
        fprintf('    alphaMax         %10.4f %-5s (the vehicle''s own clamp, from %s)\n', ...
                alphaMax,'deg',func2str(vehicleFn));
        fprintf('\n');
        fprintf('  %s\n',errLoft.message);
        fprintf('\n');
        fprintf('  WHY THE TWO ENTRIES HAVE TO MOVE TOGETHER. A bank tilts the lift vector out\n');
        fprintf('  of the vertical plane, so less of it is left to push the nose over, and the\n');
        fprintf('  deepest commanded attitudes in the bracket stop being reachable at all. The\n');
        fprintf('  shipped loftMin of %.1f deg is not a pitch programme anyone would fly; it is\n', ...
                loftMin);
        fprintf('  the degenerate tail the DEPRESSED branch collapses into at zero bank once\n');
        fprintf('  the clamp saturates, and a banked run simply does not have that tail.\n');
        fprintf('  RAISE loftMin. Measured on the shipped geometry at loftMin = -60 deg: a\n');
        fprintf('  2 deg bank flies and misses by 1.59 m, aimed off +3.271 deg; a 5 deg bank\n');
        fprintf('  flies and misses by 173.75 m, aimed off +8.205 deg. THE TWO-AXIS SOLVE\n');
        fprintf('  HANDLES THE CROSS-RANGE A BANK PRODUCES -- that is what this refusal is NOT\n');
        fprintf('  about. What it is about is a bracket the vehicle cannot be flown across.\n');
        fprintf('=========================================================\n\n');
              traj = [];
       info.refused = true;
    info.refusedWhy = 'loftBracket';
  info.branchAsked  = branch;
       info.rngReqM = rngReq;
      info.psiSeed  = psiSeed;
     info.psiLaunch = psiSeed;
    info.bankAngleDeg = bankAngle;
     info.loftMinDeg = loftMin;
     info.loftMaxDeg = loftMax;
   info.alphaMaxDeg  = alphaMax;
    info.libraryErr  = errLoft.identifier;
        if nargout == 0
            clear traj info;
        end
        return;
    end
       [~,trajStar]   = fRange(loftStarR);
              kApStar = find(trajStar.phaseIdx == 2,1,'last');
              hApoStar = trajStar.x(kApStar,1) - c.rE;
              tFlyStar = trajStar.t(end);

%% The certified maximiser interval and the ranges achieved at its two ends. The
%% smaller of those two ranges is the highest range this search can place on a
%% NAMED branch: above it the required range falls inside the unresolved band
%% where the two arcs have merged, and the answer is one arc rather than two:
              loftAML = mr.aL;
              loftBML = mr.bL;
            rngTopM   = min(mr.fAL,mr.fBL);
          coalescedRq = rngReq >= rngTopM - tolRngM;

%% -----------------------------------------------------------------
%% Solve BOTH branches, each on its own CERTIFIED side of the maximum
%% -----------------------------------------------------------------
%% Both are solved whatever the branch selector asks for, because the trade is
%% the point: only one of the two is flown, and a user choosing between them
%% needs both sets of numbers in front of them. coorbital.util.rangeSolve
%% returns converged = false with the achievable band rather than throwing, so a
%% branch that cannot reach this target simply reports that it does not exist:
               dep = solveBranch('depressed',psiSeed,rngReq,fRange,loftLoR, ...
                                 loftAML,tolRngM,cfg,rI,c,latTargetR, ...
                                 lonTargetR,loftAML,loftBML,tolLoftR);
               lof = solveBranch('lofted'   ,psiSeed,rngReq,fRange,loftBML, ...
                                 loftHiR,tolRngM,cfg,rI,c,latTargetR, ...
                                 lonTargetR,loftAML,loftBML,tolLoftR);

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

%% The CLASSICAL VACUUM EQUAL-RADIUS reference for this range angle, computed
%% from the geometry alone and from nothing this script flew. IT IS A DIAGNOSTIC
%% AND NOT A RESIDUAL: it is derived for an impulsive burn at the impact radius
%% with both endpoints of the free-flight arc on the same sphere, and this
%% script's burnout is downrange and tens of kilometres up, so matching its
%% gammaStar would enforce a condition the flown arc does not satisfy. Every
%% mode prints it as a yardstick; nothing is solved against it:
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
         pickShort = 'minimised on loft AND cutoff, not chosen between the two arcs';
          pickWhy  = ['the least-burnout-energy trajectory that meets this ' ...
                      'required range'];
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
        fprintf('  AND MINIMUM-ENERGY IS THE HARDER ASK, NOT THE EASIER ONE. It MINIMISES the\n');
        fprintf('  burnout energy subject to making the range, so every trajectory it will\n');
        fprintf('  consider carries LESS energy than the full burn -- and a trajectory that\n');
        fprintf('  carries less energy than the full burn cannot fly further than the full\n');
        fprintf('  burn does. A target the full burn cannot reach is therefore beyond this\n');
        fprintf('  mode by a wider margin than it is beyond ''lofted'' or ''depressed''. For\n');
        fprintf('  scale, the VACUUM EQUAL-RADIUS classical arc for this %.4f deg range\n', ...
                rad2deg(angReqR));
        fprintf('  angle -- a diagnostic, not the thing this mode solves -- wants\n');
        fprintf('  V* = %.1f m/s at gamma* = %.4f deg. The classical arc it would have\n', ...
                refME.V,rad2deg(refME.gamR));
        fprintf('  flown: %.3f km of apogee in %.3f min.\n', ...
                refME.hApoM./1000,refME.tofS./60);
        end
        fprintf('=========================================================\n\n');
              traj = [];
       info.refused = true;
    info.refusedWhy = 'envelope';
      info.branchAsked = branch;
       info.rngReqM = rngReq;
      info.psiSeed  = psiSeed;
     info.psiLaunch = psiSeed;
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
%% The MINIMUM-ENERGY solve: minimise eps_BO subject to the range
%% -----------------------------------------------------------------
%% Only now, with the envelope refusal past, is the loft bracket the minimisation
%% needs known to exist. See the THREE MODES note for the method:
                me = struct('solved',false);
    if isMinE
                me = minEnergySolve(psiSeed,rngReq,cfg,rI,c, ...
                                    latTargetR,lonTargetR,dep,lof, ...
                                    loftLoR,loftHiR,tolRngMEM,tolLoftMER, ...
                                    cutFracMin,nScanME);
        if ~me.solved
        fprintf('\n');
        fprintf('===== Ballistic point-to-point targeting: REFUSED =======\n');
        fprintf('  The MINIMUM-ENERGY trajectory for this geometry could not be solved.\n');
        fprintf('\n');
        fprintf('    launch           %10.4f %-5s %10.4f deg  (lat, lon)\n', ...
                latLaunch,'deg',lonLaunch);
        fprintf('    target           %10.4f %-5s %10.4f deg  (lat, lon)\n', ...
                latTarget,'deg',lonTarget);
        fprintf('    required range   %10.3f %-5s (great circle on the r = %.3f km sphere,\n', ...
                rngReq./1000,'km',rI./1000);
        fprintf('                                       a %.4f deg range angle)\n',rad2deg(angReqR));
        fprintf('    loft interval    %10.4f to %.4f deg  (the feasible family the burnout\n', ...
                rad2deg(me.loftAR),rad2deg(me.loftBR));
        fprintf('                                       energy was to be minimised along)\n');
        fprintf('    gamma* (diag.)   %10.4f %-5s (45 - Lambda/4, the VACUUM EQUAL-RADIUS\n', ...
                rad2deg(gamStarR),'deg');
        fprintf('                                       reference; not a residual of this mode)\n');
        fprintf('    V* (diagnostic)  %10.1f %-5s (the burnout speed that goes with it)\n', ...
                refME.V,'m/s');
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
      info.psiSeed  = psiSeed;
     info.psiLaunch = psiSeed;
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
          pickWhy  = sprintf(['MINIMISED, not selected: loft %.4f deg, burn cut at ' ...
                              '%.6f of full, burnout energy %.6f MJ/kg'], ...
                             rad2deg(me.loftR),me.cutFrac,me.epsBo./1e6);
    end

%% What stage 1 cost, gathered before stage 2 adds to it. A mode that ran no
%% minimisation spent nothing on one, and NaN would be wrong here rather than
%% merely unhelpful -- this figure is summed, not printed on its own:
            meProp = 0;
    if isMinE
            meProp = me.nProp;
    end
           propOne = mr.nEval + 1 + dep.nProp + lof.nProp + meProp;

%% -----------------------------------------------------------------
%% Stage 2: drive BOTH components of the miss to zero
%% -----------------------------------------------------------------
%% Everything above ran at the SEED azimuth and closed the DOWN-RANGE half of
%% the problem. This closes the other half, on the mode's own two controls:
%%
%%    'lofted' and 'depressed'   the azimuth and the LOFT ANGLE, at the full
%%                               burn. Their seed is the branch root, which
%%                               stage 1 already flew.
%%
%%    'minimum-energy'           the azimuth and the CUTOFF FRACTION, with the
%%                               loft angle HELD at the value the minimisation
%%                               settled on. See the THREE MODES note for why
%%                               the loft angle is held rather than solved a
%%                               third time.
%%
%% The safeguards are the caller's job because only the caller knows the scales:
%%
%%    dx0        1e-4 rad of heading swings the impact point about 640 m
%%               crosswise at this range, and 1e-4 of the range control's own
%%               span is about 0.0225 deg of loft or 5e-5 of the burn -- both
%%               worth hundreds of metres of range here. All three are far above
%%               the few metres an adaptive integrator's own step selection
%%               moves the impact point by, which is the floor a difference step
%%               has to clear, and all three are small enough to be a chord and
%%               not a secant across real curvature. Written against the SPANS
%%               rather than as literals so a different bracket carries them.
%%
%%    maxStep    10 deg of heading and 5 % of the range control's span per
%%               iteration. aimSolve defaults to unbounded because a library
%%               function cannot know these; this one can. The seed is already
%%               range-converged, so a step anywhere near either cap means the
%%               Newton direction is wrong -- and on THIS script the loft cap
%%               does a second job the HGV one does not have to: an unbounded
%%               loft step can cross the max-range hump and land the solve on
%%               the other branch, which the branch measurement below would then
%%               report as a disagreement rather than as a step nobody wanted.
%%
%%    tolF       the user's tolerance over sqrt(2), so that BOTH components
%%               inside it puts the total miss inside the figure asked for.
%%               aimSolve tests max(abs(f)), which is the infinity norm, and the
%%               worst case for the hypotenuse is the two components equal.
%%
%% THE SEED IS A CONTROL PAIR THIS RUN HAS ALREADY FLOWN TO COMPLETION, which is
%% the sanity check aimSolve's contract asks of a caller: it throws at x0 rather
%% than shrinking a step it has not got, so a seed that will not fly is an input
%% error and not a failure to converge. Every path into this block came through
%% a solve that propagated exactly these controls and kept the trajectory. The
%% try/catch is belt and braces for the one thing that would defeat that -- a
%% propagation that is not repeatable -- and it turns the throw into this
%% script's own refusal rather than letting a library identifier out:
           spanCtl = loftHiR - loftLoR;
    if isMinE
           spanCtl = 1 - cutFracMin;
    end
           tolAimM = tolRngM./sqrt(2);
            loftHd = pick.loftR;
    if isMinE
            residF = @(xAim) flyMiss(xAim(1),loftHd,xAim(2),cfg,aim);
             xSeed = [psiSeed; pick.cutFrac];
    else
            residF = @(xAim) flyMiss(xAim(1),xAim(2),1,cfg,aim);
             xSeed = [psiSeed; pick.loftR];
    end
             dxAim = [1.0e-4; 1.0e-4.*spanCtl];
           optsAim = struct('maxStep',[deg2rad(10); 0.05.*spanCtl]);
    try
    [xAim,fAch,av] = coorbital.util.aimSolve(residF,xSeed,dxAim,tolAimM,optsAim);
    catch errAim
        error('coorbital:runBallisticTarget:seedWillNotFly', ...
            ['The two-axis solve could not evaluate its own SEED, which is a ' ...
             'control pair this run had already propagated to completion: ' ...
             'azimuth %.6f deg with %s. That means the propagation is not ' ...
             'repeatable, and nothing downstream of it could be trusted. The ' ...
             'solver said: %s'],rad2deg(xSeed(1)), ...
            ctlWord(isMinE,loftHd,xSeed(2)),errAim.message);
    end
         psiLaunch = xAim(1);
    if isMinE
           loftSol = loftHd;
            cutSol = xAim(2);
    else
           loftSol = xAim(2);
            cutSol = 1;
    end

%% The residual at the SEED, which aimSolve records as the first column of its
%% history and which therefore costs nothing to read. It is the miss the
%% closed-form bearing and the stage-1 range control produce on their own --
%% exactly what this script used to report as its answer -- so the summary can
%% print the before and the after side by side:
          seedDwnM = av.fHist(1,1);
          seedXtkM = av.fHist(2,1);
         seedMissM = hypot(seedDwnM,seedXtkM);

%% -----------------------------------------------------------------
%% Refuse, loudly, when the two-axis solve did not converge
%% -----------------------------------------------------------------
%% aimSolve returns converged = false with the best point it reached rather than
%% throwing, for the same reason rangeSolve does: only this script has the
%% vocabulary. A NEAR MISS RETURNED AS AN ANSWER IS THE FAILURE MODE THIS WHOLE
%% FILE EXISTS TO PREVENT, so it is refused on exactly the terms the envelope
%% refusal uses -- nothing thrown, empty traj, info.refused true:
    if ~av.converged
        fprintf('\n');
        fprintf('===== Ballistic point-to-point targeting: REFUSED =======\n');
        fprintf('  The two-axis aim solve did not converge on the target.\n');
        fprintf('\n');
        fprintf('    launch           %10.4f %-5s %10.4f deg  (lat, lon)\n', ...
                latLaunch,'deg',lonLaunch);
        fprintf('    target           %10.4f %-5s %10.4f deg  (lat, lon)\n', ...
                latTarget,'deg',lonTarget);
        fprintf('    required range   %10.3f %-5s (great circle on the r = %.3f km sphere)\n', ...
                rngReq./1000,'km',rI./1000);
        fprintf('    branch asked     %-16s (%s)\n',branch,pickShort);
        fprintf('    reachable        %10.3f to %.3f km  (flown at the SEED azimuth)\n', ...
                rngMin./1000,rngMax./1000);
        fprintf('    seed azimuth     %10.6f %-5s (the closed-form great-circle bearing)\n', ...
                rad2deg(psiSeed),'deg');
        fprintf('    best azimuth     %10.6f %-5s (%+.6f deg from the seed)\n', ...
                rad2deg(psiLaunch),'deg',rad2deg(wrapPi(psiLaunch - psiSeed)));
        fprintf('    best %s\n',ctlWord(isMinE,loftSol,cutSol));
        fprintf('    seed miss        %10.2f %-5s (down-range %+.2f, cross-range %+.2f)\n', ...
                seedMissM,'m',seedDwnM,seedXtkM);
        fprintf('    best miss        %10.2f %-5s (down-range %+.2f, cross-range %+.2f)\n', ...
                hypot(fAch(1),fAch(2)),'m',fAch(1),fAch(2));
        fprintf('    tolerance        %10.2f %-5s (per component; %.2f m on the total miss)\n', ...
                tolAimM,'m',tolRngM);
        fprintf('    propagations     %10d %-5s (%d before the Newton, %d in it)\n', ...
                propOne + av.nEval,'',propOne,av.nEval);
        fprintf('    solver said      %s\n',av.identifier);
        fprintf('\n');
        fprintf('  %s\n',av.message);
        fprintf('\n');
        fprintf('  THE BEST POINT REACHED IS NOT AN ANSWER, and it is not being returned as\n');
        fprintf('  one. A targeting script that handed back its nearest miss would hand back\n');
        fprintf('  something that LOOKS like a solution: a loft angle, an azimuth, a\n');
        fprintf('  trajectory and a summary, all of them describing a flight that does not\n');
        fprintf('  arrive. The numbers above are in info.aimInfo for anyone who wants to\n');
        fprintf('  restart from them.\n');
        fprintf('=========================================================\n\n');
              traj = [];
       info.refused = true;
    info.refusedWhy = 'aimSolve';
  info.branchAsked  = branch;
       info.rngReqM = rngReq;
      info.psiSeed  = psiSeed;
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
      info.loftDeg  = rad2deg(loftSol);
      info.cutFrac  = cutSol;
        info.missM  = hypot(fAch(1),fAch(2));
        info.downM  = fAch(1);
      info.xTrackM  = fAch(2);
     info.seedMissM = seedMissM;
     info.seedDownM = seedDwnM;
   info.seedXTrackM = seedXtkM;
      info.tolAimM  = tolAimM;
      info.aimIter  = av.iterations;
      info.aimEval  = av.nEval;
      info.aimInfo  = av;
      info.propOne  = propOne;
        info.nProp  = propOne + av.nEval;
        if nargout == 0
            clear traj info;
        end
        return;
    end

%% -----------------------------------------------------------------
%% Fly the solved trajectory
%% -----------------------------------------------------------------
%% One more propagation at the converged controls. aimSolve already flew them
%% but returns only the two-component residual, so the state history is
%% re-created here rather than carried out through a persistent variable. The
%% assert is the repeatability check that licenses everything below: if the same
%% controls give a different impact point, no number in the summary means
%% anything:
     [misVec,traj]  = flyMiss(psiLaunch,loftSol,cutSol,cfg,aim);
    assert(max(abs(misVec - fAch)) < 1e-6, ...
        ['re-flying the solved controls gave a miss of [%.6f, %.6f] m against ' ...
         'the solver''s [%.6f, %.6f] m; the propagation is not repeatable and ' ...
         'nothing below can be trusted'],misVec(1),misVec(2),fAch(1),fAch(2));

%% The flown arc, re-measured from the trajectory that was actually integrated
%% rather than carried forward from the stage-1 record it started as. Written
%% back onto pick -- and onto me, which pick is a copy of in that mode -- so
%% every line of the summary below describes the FLOWN controls and not the
%% seed ones. The stage-1 branch records dep and lof are deliberately left
%% alone: they are the two-arc TRADE TABLE, measured at the seed azimuth, and
%% only one of the two was ever refined:
              flwn = arcRecord(traj,psiLaunch,loftSol,cutSol,cfg,rI,c, ...
                               latTargetR,lonTargetR);
              pick = mergeInto(pick,flwn);
    if isMinE

%% THE MINIMISATION'S OWN ANSWER IS SAVED BEFORE THE FLOWN ONE OVERWRITES IT,
%% and this is the whole of the fix to a report that used to mix two operating
%% points inside one paragraph. The merge below is right for everything that
%% DESCRIBES THE FLIGHT -- apogee, flight time, burnout state, mass, miss -- and
%% wrong for the four numbers that SUBSTANTIATE THE MINIMUM, because those were
%% all measured at the seed azimuth and one of them is the minimum itself.
%%
%% What went wrong when it did not: on the rotating shipped run the minimisation
%% settled -45.977595 MJ/kg at a cutoff of 0.992294, stage 2 moved the cutoff to
%% 0.992621 to close the cross-range and the vehicle flew -45.924475 MJ/kg, and
%% the summary printed the FLOWN energy in a column headed MINIMISED and then
%% differenced it against neighbours measured at the OTHER azimuth. The valley
%% depths came out 284203.3 and 59729.3 J/kg against the true seed-azimuth
%% 337323.6 and 112849.6 -- the far side understated by 47 per cent -- and the
%% printed noise ratio did not even close against the printed depths, because it
%% was formed from the unprinted one. Every one of those numbers was individually
%% computed correctly; what was wrong was putting them in the same sentence.
%%
%% The still-Earth shipped run is honest either way, because there the aim solve
%% is a no-op and all four numbers share one basis. That is exactly why nothing
%% caught it, and why the assertions that now do are in part 12 rather than
%% part 6:
       me.seedEpsBo   = me.epsBo;
      me.seedCutFrac  = me.cutFrac;
        me.seedTCutS  = me.tCutS;
         me.seedRngM  = me.rngM;
                me = mergeInto(me,flwn);
        me.propLeftKg = me.mCutKg - cfg.mBurnout;
        me.hApoKepM   = keplerApogee(c.rE + me.hBoM,me.vBoM,me.gamBoR,c) - c.rE;

%% How far stage 2 moved the range control and what it cost in burnout energy.
%% Zero on any run where the aim solve is a no-op, which is the case the two
%% columns of the table below then agree in:
        me.dCutFrac   = me.cutFrac - me.seedCutFrac;
        me.dEpsBoJkg  = me.epsBo   - me.seedEpsBo;

%% ...and the one-line grounds this mode reports itself on. It names BOTH
%% points, because naming only the minimised one describes an arc the vehicle
%% did not fly and naming only the flown one calls a number MINIMISED that
%% nothing minimised:
          pickWhy  = sprintf(['MINIMISED at the seed azimuth to %.6f MJ/kg at loft ' ...
                              '%.4f deg and cutoff %.6f, FLOWN at %.6f MJ/kg at ' ...
                              'cutoff %.6f after the aim solve'], ...
                             me.seedEpsBo./1e6,rad2deg(me.loftR),me.seedCutFrac, ...
                             me.epsBo./1e6,me.cutFrac);
    end

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
%% The miss, measured from the FLOWN state
%% -----------------------------------------------------------------
%% missM is the great-circle distance impact-to-target, measured directly, and
%% the two components are the same vector resolved in the aim frame. They agree
%% BY CONSTRUCTION -- the components are the arc times a unit direction -- so the
%% assert below is checking the arithmetic and not the physics, which is why it
%% is an equality and not a tolerance a reader is invited to judge:
             latIm = traj.x(end,3);
             lonIm = traj.x(end,2);
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
%% the required one. NOT what the solve drives any more -- that is the miss above
%% -- but it is the number a reader has always used to judge whether the loft
%% angle came out sensible, and on a rotating or banked run the gap between it
%% and the down-range component is itself informative: it is the part of the
%% miss the old downrange-only solve could not see:
              resM = pick.rngM - rngReq;

%% How far the flown track ended up from the arc it left on. With rotation off
%% and the bank zero this is zero to rounding; with either on it is the
%% deflection the azimuth solve had to aim off against, so it is a measurement
%% of the coupling rather than a caveat about it:
            angFly = coorbital.util.greatCircle(latLaunchR,lonLaunchR,latIm,lonIm);
            psiFly = coorbital.util.greatCircleBearing(latLaunchR,lonLaunchR, ...
                                                       latIm,lonIm);
            dPsiR  = wrapPi(psiFly - psiLaunch);

%% The heading the vehicle actually left the pad on must be the one the solve
%% asked for, to machine precision -- flyLoft writes entry 6 of the state from it
%% and ode45 never touches the initial condition, so any disagreement here means
%% the launch state was assembled wrong:
             psiX0 = traj.x(1,6);
    assert(abs(wrapPi(psiX0 - psiLaunch)) < 1e-12, ...
        ['the flown initial heading %.15f rad is not the solved launch ' ...
         'azimuth %.15f rad; the launch state was assembled wrong'], ...
        psiX0,psiLaunch);

%% -----------------------------------------------------------------
%% WHICH BRANCH DID IT ACTUALLY FLY
%% -----------------------------------------------------------------
%% Measured from the ROOT'S POSITION relative to the certified maximiser
%% interval [aL,bL], which is the only invariant available: below aL the loft
%% angle is on the depressed side of every point the maximum could occupy, above
%% bL it is on the lofted side, and inside it the two arcs cannot be told apart
%% by any evidence this search has.
%%
%% IT USED TO BE MEASURED FROM THE FLOWN APOGEE AND FLIGHT TIME against the
%% max-range arc's, and that was wrong. Neither is a branch invariant for a
%% finite powered atmospheric arc -- drag, lift, burnout altitude and boost
%% duration can make either non-monotone in the loft angle -- and near the
%% maximum both differences vanish quadratically, so their signs are set by
%% search and integration error rather than by the branch. They are still
%% computed, below, and reported as the descriptive quantities they are.
%%
%% MINIMUM-ENERGY IS NOT ON A BRANCH AT ALL. Its loft angle lies inside the
%% interval the two full-burn arcs span, and it does not burn to exhaustion, so
%% the lofted/depressed vocabulary simply does not apply to it. Its verification
%% is the minimisation record the paragraph below prints:
         flownName = 'minimum-energy';
          branchOK = true;
    if ~isMinE
         flownName = branchOfLoft(pick.loftR,loftAML,loftBML,tolLoftR);
          branchOK = strcmp(flownName,pickName);
    end

%% The DESCRIPTIVE apogee and flight-time readings, kept because they are what a
%% reader recognises a lofted arc by, and because their disagreeing with each
%% other is worth printing. They classify nothing:
   [byApoName,flownAgree] = measureBranch(pick.hApoM,pick.tFlyS, ...
                                          hApoStar,tFlyStar);
         byTimName = 'depressed';
    if pick.tFlyS > tFlyStar
         byTimName = 'lofted';
    end

%% THE CASE IN WHICH THE LABEL AND THE MEASUREMENT LEGITIMATELY PART COMPANY:
%% the required range is within the range tolerance of the largest range this
%% search can certify, so coorbital.util.rangeSolve short-circuits at the
%% bracket endpoint and BOTH branches come back holding an arc inside the
%% unresolved maximiser interval. The two arcs ARE the same arc at that range
%% and the lofted/depressed distinction has no content, which the summary says
%% instead of blaming the bracketing. Tested on the ACHIEVED range against the
%% maximum and on the flown loft angle against the certified interval -- never
%% on floating-point equality with a midpoint, which depended on whether
%% rangeSolve happened to return an endpoint's exact bit pattern:
        nearMaxRng = pick.rngM >= rngMax - tolRngM;
        inMaxIntvl = pick.loftR >= loftAML - tolLoftR && ...
                     pick.loftR <= loftBML + tolLoftR;
        atMaxRange = nearMaxRng && inMaxIntvl;

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
    if isMinE
              tCut = me.tCutS;
    end
        [why1,ok1] = whyBurnout(traj.t(kBO),massS(kBO),mBurnoutT,tMaxBoost,tCut);
        [why2,ok2] = whyApogee(traj.t(kApo) - traj.t(kBO),traj.t(kApo), ...
                               traj.x(kApo,5),tMaxCoast);
        [why3,ok3] = whyImpact(traj.t(nS) - traj.t(kApo),traj.t(nS), ...
                               traj.x(nS,1) - c.rE,hStopM,hStop,tMaxDesc);
             allOK = ok1 && ok2 && ok3;

%% THERE ARE DELIBERATELY NO EASTWARD-DRIFT FIGURES HERE ANY MORE, and the
%% deletion is the point. This script used to compute
%% c.omegaE*rI*cos(latTarget)*tFly for each branch -- how far the ground at the
%% target's latitude sweeps east over that arc's flight, hundreds of kilometres
%% -- and print the pair as "how much it would matter" that the Earth turns. It
%% is not. It is the scale of a mechanism this file no longer claims: the state
%% is Earth-fixed and the speed planet-relative, so the target does not travel
%% relative to anything the solve can see. The quantity that IS the scale of the
%% rotation effect is the deflection of the VEHICLE, and it is measured rather
%% than derived -- seedMissM, and its cross-range component beside it. A
%% retracted mechanism's number left standing under a corrected sentence is the
%% harder error to notice of the two, because the prose reads right.

%% Whether the Earth turned on this run. It decides which paragraph the
%% limitations block prints, and it is read from env.omegaE rather than from the
%% user block so that a run driven through opts cannot be described wrongly:
           spinsUp = env.omegaE ~= 0;
           spinTxt = 'OFF';
    if spinsUp
           spinTxt = 'ON';
    end
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
    fprintf('                                        seed azimuth and the stage-1 range control\n');
    fprintf('                                        put the vehicle on their own)\n');
    fprintf('    flown azimuth    %10.6f %-5s  (bearing launch to impact, from the flown\n', ...
            rad2deg(psiFly),'deg');
    fprintf('                                        state; %+.3e deg from the commanded one.\n', ...
            rad2deg(dPsiR));
    fprintf('                                        Non-zero means the track was deflected off\n');
    fprintf('                                        the arc it left on -- by rotation, by bank,\n');
    fprintf('                                        or by both)\n');
    fprintf('\n');
    fprintf('  The range-versus-loft hump, bracketed before either branch was solved\n');
    fprintf('    max-range loft   %10.4f %-5s  (commanded terminal pitch attitude)\n', ...
            rad2deg(loftStarR),'deg');
    fprintf('    MAXIMUM RANGE    %10.3f %-5s  (nothing this vehicle can do exceeds it)\n', ...
            rngMax./1000,'km');
    fprintf('    CERTIFIED interval %8.4f to %.4f deg   (all a golden section can claim: the\n', ...
            rad2deg(loftAML),rad2deg(loftBML));
    fprintf('                                        maximiser is IN here. The angle above is\n');
    fprintf('                                        its midpoint, reported, not relied on)\n');
    fprintf('    ranges at its ends %8.3f and %.3f km   (the smaller of the two is the\n', ...
            mr.fAL./1000,mr.fBL./1000);
    fprintf('                                        highest range that can be put on a NAMED\n');
    fprintf('                                        branch; above it the two arcs have merged)\n');
    fprintf('    its apogee       %10.3f %-5s  and flight time %.3f s. DESCRIPTIVE only --\n', ...
            hApoStar./1000,'km',tFlyStar);
    fprintf('                                        neither is a branch invariant for a\n');
    fprintf('                                        powered atmospheric arc, so neither\n');
    fprintf('                                        classifies anything below.\n');
    fprintf('    loft bracket     %10.4f to %.4f deg   (searched; the maximum is INTERIOR,\n', ...
            loftMin,loftMax);
    fprintf('                                        which is what gives two branches at all)\n');
    fprintf('    reachable        %10.3f to %.3f km   (union of the two branch bands, flown\n', ...
            rngMin./1000,rngMax./1000);
    fprintf('                                        at the SEED azimuth, so indicative for the\n');
    fprintf('                                        solved one rather than exact)\n');
    fprintf('    search cost      %10d %-5s  (%d coarse scan points after %d refinement(s),\n', ...
            mr.nEval,'',mr.nScan,mr.nRefine);
    fprintf('                                        %d golden-section steps, closed to %.4f deg)\n', ...
            mr.nGolden,rad2deg(mr.widthR));
    if coalescedRq
        fprintf('  *** CAUTION ***  the required %.3f km is inside the %.3f km tolerance of the\n', ...
                rngReq./1000,tolRngM./1000);
        fprintf('                   %.3f km the search can certify, so it falls in the\n', ...
                rngTopM./1000);
        fprintf('                   UNRESOLVED band around the maximum. The lofted and depressed\n');
        fprintf('                   arcs are the same arc at this range and the branch label\n');
        fprintf('                   below carries no information. Tighten tolRangeKm and\n');
        fprintf('                   tolLoftDeg, or move the target inside the envelope.\n');
    end
    fprintf('\n');
    fprintf('  THE TWO ARCS THAT REACH THIS TARGET, both at the SEED azimuth\n');
    fprintf('    Only the flown one was refined by the two-axis solve; this table is the\n');
    fprintf('    TRADE, and the flown arc''s own solved numbers are above and below it.\n');
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
    fprintf('    NOT ON A BRANCH. Its loft angle lies inside the interval the two full-burn\n');
    fprintf('    arcs span and it does not burn to exhaustion, so the lofted/depressed\n');
    fprintf('    vocabulary does not apply. The minimisation record below is what verifies\n');
    fprintf('    it. For orientation only: %.3f km of apogee in %.2f s.\n', ...
            pick.hApoM./1000,pick.tFlyS);
    else
    fprintf('    MEASURED as      %-16s (from the flown %.4f deg loft angle against\n', ...
            flownName,rad2deg(pick.loftR));
    fprintf('                                      the CERTIFIED maximiser interval %.4f to\n', ...
            rad2deg(loftAML));
    fprintf('                                      %.4f deg. Position is the only branch\n', ...
            rad2deg(loftBML));
    fprintf('                                      invariant this search has)\n');
    fprintf('    described by     apogee %.3f km says %s, flight time %.3f s says %s\n', ...
            pick.hApoM./1000,byApoName,pick.tFlyS,byTimName);
    end
    if ~isMinE && ~flownAgree
        fprintf('  *** CAUTION ***  the apogee and the flight time DESCRIBE the flight\n');
        fprintf('                   differently -- one puts it above the max-range arc and the\n');
        fprintf('                   other below. Neither classifies the branch, so this changes\n');
        fprintf('                   nothing above; it says the flight is close enough to the\n');
        fprintf('                   maximum that those two readings have stopped agreeing.\n');
    end
    if ~branchOK && atMaxRange
        fprintf('  *** CAUTION ***  the flown trajectory MEASURES as %s while the solve was run\n', ...
                flownName);
        fprintf('                   on the %s bracket, because the required %.3f km is inside\n', ...
                pickName,rngReq./1000);
        fprintf('                   the %.3f km tolerance of the %.3f km maximum and the flown\n', ...
                tolRngM./1000,rngMax./1000);
        fprintf('                   loft angle landed INSIDE the certified maximiser interval.\n');
        fprintf('                   The two arcs ARE the same arc at this range and the\n');
        fprintf('                   lofted/depressed distinction has no meaning here. The\n');
        fprintf('                   bracketing is NOT suspect. Tighten tolRangeKm and\n');
        fprintf('                   tolLoftDeg, or take the max-range arc and stop asking which\n');
        fprintf('                   side of itself it is on.\n');
    elseif ~branchOK
        fprintf('  *** CAUTION ***  the flown trajectory MEASURES as %s while the solve was run\n', ...
                flownName);
        fprintf('                   on the %s bracket. STAGE 1''S ROOT WAS ON THE RIGHT SIDE:\n', ...
                pickName);
        fprintf('                   that bracket lies entirely on one certified side of the\n');
        fprintf('                   maximiser interval and bisection cannot leave a bracket.\n');
        fprintf('                   STAGE 2''S NEWTON IS NOT BRACKETED, and it is what moved the\n');
        fprintf('                   loft angle here: %.4f deg at the seed against %.4f deg\n', ...
                rad2deg(xSeed(2)),rad2deg(loftSol));
        fprintf('                   flown, %+.4f deg of travel while it aimed off %+.4f deg to\n', ...
                rad2deg(loftSol - xSeed(2)),rad2deg(wrapPi(psiLaunch - psiSeed)));
        fprintf('                   close the cross-range. THE BRACKETING IS NOT SUSPECT. Two\n');
        fprintf('                   things to know before reading this as a fault: the certified\n');
        fprintf('                   interval printed above is a SEED-AZIMUTH quantity, and the\n');
        fprintf('                   maximum itself moves with the aim -- about 0.28 deg of loft\n');
        fprintf('                   per degree of azimuth on this geometry, against an interval\n');
        fprintf('                   only %.4f deg wide -- so a flown angle just outside it may\n', ...
                rad2deg(loftBML - loftAML));
        fprintf('                   still be on the named branch AT THE AZIMUTH IT WAS FLOWN.\n');
        fprintf('                   That is why the Newton is not boxed to the interval: a box\n');
        fprintf('                   drawn at the seed azimuth could exclude the genuine root.\n');
        fprintf('                   Do not read the branch label; read the loft angle against\n');
        fprintf('                   the interval, and treat both as approximate by that much.\n');
    end
%% THE MINIMUM-ENERGY PARAGRAPH, and it is the only place in this summary that
%% compares a flown trajectory against a result computed from the geometry
%% alone. It prints the two solved residuals -- what makes the arc
%% minimum-energy -- and then the two MODELLING residuals against the classical
%% arc, which are physics rather than solver error and are attributed as such:
    if isMinE
    fprintf('  THE MINIMUM-ENERGY SOLVE\n');
    fprintf('    A CONSTRAINED MINIMISATION, stated and then solved: minimise the burnout\n');
    fprintf('    specific energy V^2/2 - mu/r over the LOFT ANGLE and the CUTOFF FRACTION,\n');
    fprintf('    subject to the achieved range being the required %.3f km. Two parameters\n', ...
            rngReq./1000);
    fprintf('    and one constraint leave a ONE-DIMENSIONAL feasible family; the cutoff is\n');
    fprintf('    solved for the constraint at each loft angle and the energy is minimised\n');
    fprintf('    along what is left.\n');
%% TWO COLUMNS, AND THE COLUMN HEADINGS ARE THE POINT. The minimisation ran at
%% the SEED azimuth; the two-axis solve then moved the cutoff to close the
%% cross-range, and what the vehicle flew is the second column. They are equal
%% on any run where the aim solve is a no-op -- which is the shipped one -- and
%% they must be shown apart rather than silently merged wherever it is not,
%% because only the first column is the MINIMISED quantity and only the second
%% describes the flight:
    fprintf('    %-22s %18s %18s\n','','MINIMISED (seed)','FLOWN (solved)');
    fprintf('    %-22s %18.4f %18s\n','loft angle (deg)',rad2deg(me.loftR),'held, same');
    fprintf('    %-22s %18.6f %18.6f\n','cutoff fraction (-)',me.seedCutFrac,me.cutFrac);
    fprintf('    %-22s %18.4f %18.4f\n','cutoff time (s)',me.seedTCutS,me.tCutS);
    fprintf('    %-22s %18.6f %18.6f\n','BURNOUT ENERGY (MJ/kg)', ...
            me.seedEpsBo./1e6,me.epsBo./1e6);
    fprintf('    %-22s %18.4f %18.4f\n','achieved range (km)', ...
            me.seedRngM./1000,me.rngM./1000);
    if me.dCutFrac == 0
    fprintf('      THE TWO COLUMNS ARE THE SAME RUN. The seed azimuth was already inside the\n');
    fprintf('      aim tolerance, so the two-axis solve converged at iteration 0 and returned\n');
    fprintf('      the minimisation''s own answer unmoved. The minimised energy IS the flown\n');
    fprintf('      energy here, and every number below shares one basis with both.\n');
    else
    fprintf('      THE TWO COLUMNS ARE DIFFERENT OPERATING POINTS, and the difference is the\n');
    fprintf('      price of the aim: stage 2 moved the cutoff by %+.6f of the burn to close\n', ...
            me.dCutFrac);
    fprintf('      %.2f km of cross-range, and that cost %+.1f J/kg of burnout energy. THE\n', ...
            abs(seedXtkM)./1000,me.dEpsBoJkg);
    fprintf('      FLOWN COLUMN IS NOT A MINIMUM OF ANYTHING -- the minimisation is the first\n');
    fprintf('      column, at the seed azimuth, and the valley evidence below is on THAT\n');
    fprintf('      basis. Re-minimising the loft angle at the solved azimuth is the work this\n');
    fprintf('      mode does not do; see the closing block.\n');
    end
    fprintf('      The range residual is %+.2f m against the flown column. The %.1f m\n', ...
            me.rngM - rngReq,tolRngMEM);
    fprintf('      tolerance quoted for the inner constraint governed STAGE 1, at the seed\n');
    fprintf('      azimuth; what governs the flown column is the two-axis tolerance of\n');
    fprintf('      %.1f m PER COMPONENT.\n',tolAimM);
    fprintf('      %d scan points and %d golden-section steps closed the loft angle to\n', ...
            me.nScan,me.nGolden);
    fprintf('      %.4f deg, over %d propagations.\n',rad2deg(me.widthR),me.nProp);
    fprintf('    IT IS A MINIMUM, and here is the evidence rather than the assertion. ALL OF\n');
    fprintf('    IT IS AT THE SEED AZIMUTH, against the MINIMISED column above and not the\n');
    fprintf('    flown one. The two neighbouring FEASIBLE points -- same range, different\n');
    fprintf('    loft -- leave burnout at %.6f and %.6f MJ/kg against that one''s\n', ...
            me.epsNeighLo./1e6,me.epsNeighHi./1e6);
    fprintf('    %.6f, so the valley is %.1f J/kg deep on the near side and %.1f J/kg\n', ...
            me.seedEpsBo./1e6,me.epsNeighLo - me.seedEpsBo, ...
            me.epsNeighHi - me.seedEpsBo);
    fprintf('    on the far side.\n');
    fprintf('    AND THE NOISE IS MEASURED, NOT ASSUMED SMALL. Re-solving the constraint at\n');
    fprintf('    the settled loft angle with a tolerance ten times tighter moves the burnout\n');
    fprintf('    energy by %.3e J/kg, which is %.2e of the shallower side of that valley\n', ...
            me.epsNoiseJkg,me.epsNoiseRel);
    fprintf('    -- %.1f J/kg, the shallower depth, which is the one the ratio divides by\n', ...
            me.epsDepthJkg);
    fprintf('    so that the arithmetic closes against the numbers printed here.\n');
    if me.epsNoiseRel > 0.1
    fprintf('    *** CAUTION *** that ratio is not small. The minimum is not resolved against\n');
    fprintf('    the inner tolerance; tighten tolRangeMEKm before believing the loft angle.\n');
    end
    fprintf('    AGAINST THE CLASSICAL VACUUM EQUAL-RADIUS ARC for this range angle, which is\n');
    fprintf('    computed from the geometry alone and is a DIAGNOSTIC, not a residual: it\n');
    fprintf('    assumes an impulsive burn with both ends of the free-flight arc on the same\n');
    fprintf('    sphere, and this burnout is %.2f km up and downrange of the pad. gamma* is\n', ...
            me.hBoM./1000);
    fprintf('    printed BESIDE the achieved burnout gamma, not driven to it:\n');
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
    fprintf('    THE DIFFERENCES ARE PHYSICS, NOT SOLVER ERROR, AND THEY ARE NOT RESIDUALS.\n');
    fprintf('    The classical result assumes an IMPULSIVE burn at the impact radius in a\n');
    fprintf('    VACUUM. This flight burns for %.3f s and finishes %.2f km up at %.1f m/s,\n', ...
            me.tCutS,me.hBoM./1000,me.vBoM);
    fprintf('    not at zero altitude at %.1f m/s, and a Keplerian arc from THAT burnout\n', ...
            refME.V);
    fprintf('    state apogees at %.3f km -- which is the flown figure to %.3f km. Coast\n', ...
            me.hApoKepM./1000,abs(me.hApoM - me.hApoKepM)./1000);
    fprintf('    drag supplies the rest. Expect a few per cent from any finite boost, and do\n');
    fprintf('    NOT read the gamma line as a convergence check: nothing above was solved\n');
    fprintf('    against gamma*, and driving the burnout gamma to a value derived for an\n');
    fprintf('    equal-radius vacuum arc is exactly the wrong condition to impose here.\n');
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
            me.seedEpsBo./1e6,me.propLeftKg);
    fprintf('    is the whole reason neither full-burn arc is the trajectory above. All\n');
    fprintf('    three energies here are at the SEED azimuth, which is the only basis the\n');
    fprintf('    two full-burn arcs were ever measured on.\n');
    else
    fprintf('    ONLY ONE FULL-BURN ARC REACHES THIS TARGET, so there is no pair of burnout\n');
    fprintf('    energies to contrast this one against; the %s arc alone leaves at\n', ...
            onlyBranch(dep.exists));
    fprintf('    %.4f MJ/kg against this trajectory''s %.4f MJ/kg -- both at the SEED\n', ...
            onlyEps(dep,lof),me.seedEpsBo./1e6);
    fprintf('    azimuth, which is the only basis the full-burn arc was measured on. This\n');
    fprintf('    one threw %.1f kg of propellant away unburned.\n',me.propLeftKg);
    fprintf('    THE MINIMUM-ENERGY SOLVE IS UNAFFECTED by the\n');
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
    fprintf('                                        of this mode: the cutoff enforces the range\n');
    fprintf('                                        constraint and the loft angle is what the\n');
    fprintf('                                        burnout energy is minimised over)\n');
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
            propOne + av.nEval + 1,'');
    fprintf('                                        %d bracketing the maximum, 1 re-flying it,\n', ...
            mr.nEval);
    if isMinE
    fprintf('                                        %d on the depressed branch, %d on the lofted,\n', ...
            dep.nProp,lof.nProp);
    fprintf('                                        %d on the minimum-energy solve,\n',meProp);
    else
    fprintf('                                        %d on the depressed branch, %d on the lofted,\n', ...
            dep.nProp,lof.nProp);
    end
    fprintf('                                        %d in the two-axis aim solve and 1 to re-fly\n', ...
            av.nEval);
    fprintf('                                        the answer)\n');
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
    fprintf('       463.211 km of miss, 462.870 km of it crosswise, and the solve removes it.\n');
    end
    fprintf('    2. THE SOLVE CONTROLS BOTH AXES. Down-range and cross-range are driven to zero\n');
    fprintf('       together by a damped Newton on the launch azimuth and the mode''s own range\n');
    fprintf('       control, so cross-range is a MEASUREMENT of the solve here and not a caveat\n');
    fprintf('       about it.\n');
    if bankRad ~= 0
    fprintf('       THIS RUN COMMANDS A NON-ZERO BANK of %.1f deg, so the ground track is NOT\n', ...
            bankAngle);
    fprintf('       the launch-to-target great circle and a downrange-only solve could not see\n');
    fprintf('       the difference. Measured on this run: the seed put the vehicle %.2f m off,\n', ...
            seedMissM);
    fprintf('       %.2f m of it crosswise, and the second control closed it to %.2f m.\n', ...
            abs(seedXtkM),missM);
    elseif spinsUp
    fprintf('       THE BANK ANGLE IS ZERO ON THIS RUN, so the %.2f m of cross-range the seed\n', ...
            abs(seedXtkM));
    fprintf('       left came ENTIRELY FROM THE ROTATION and not from any commanded turn. It is\n');
    fprintf('       now %.2f m. Bank as well and the two effects simply add into the same\n', ...
            abs(xTrackM));
    fprintf('       residual; the solve does not need to tell them apart.\n');
    else
    fprintf('       THE BANK ANGLE IS ZERO ON THIS RUN and rotation is off, so the seed left\n');
    fprintf('       only %.2e m of cross-range: an unbanked track over a non-rotating sphere\n', ...
            abs(seedXtkM));
    fprintf('       never leaves the great circle it departed on, so the seed azimuth WAS the\n');
    fprintf('       answer here and the Newton returned it unchanged.\n');
    fprintf('       THAT IS A PROPERTY OF THIS CONFIGURATION, NOT A GENERAL GUARANTEE, and it\n');
    fprintf('       is measured above rather than assumed. Bank a phase, or turn the Earth on,\n');
    fprintf('       and the track curves off the arc: at earthSpin true this geometry puts the\n');
    fprintf('       seed 463.2 km wide, and the azimuth is aimed off to cancel it.\n');
    end
    fprintf('    3. THE NEWTON IS LOCAL, AND WHAT IT WAS SEEDED FROM WAS MEASURED AT ONE\n');
    fprintf('       AZIMUTH. The hump, the CERTIFIED maximiser interval, the two monotone\n');
    fprintf('       branch brackets and the reachable band above were all flown at the SEED\n');
    fprintf('       azimuth -- range against loft is a curve there and a surface otherwise, and\n');
    fprintf('       a golden section on a surface certifies nothing. They are therefore\n');
    fprintf('       indicative for the %.6f deg azimuth the solve settled on rather than\n', ...
            rad2deg(psiLaunch));
    fprintf('       exact -- it moved %+.6f deg off the seed -- and\n', ...
            rad2deg(wrapPi(psiLaunch - psiSeed)));
    fprintf('       the too-far and too-close refusals are approximate by the same amount. The\n');
    fprintf('       Newton itself is seeded, damped and capped at 10 deg of heading per\n');
    fprintf('       iteration; when it cannot reach the tolerance the run is REFUSED with the\n');
    fprintf('       best point it found, never returned as a near miss. The tolerance is\n');
    fprintf('       %.2f m PER COMPONENT so the two together stay inside the %.2f m tolRangeKm\n', ...
            tolAimM,tolRngM);
    fprintf('       asks for.\n');
    if isMinE
    fprintf('       THE LOFT ANGLE WAS HELD, NOT SOLVED, in stage 2 of this mode: the two-axis\n');
    fprintf('       Newton moved the azimuth and the CUTOFF FRACTION and left the loft angle at\n');
    fprintf('       the %.4f deg the energy minimisation settled on at the seed azimuth. The\n', ...
            rad2deg(loftSol));
    fprintf('       constraint is therefore met EXACTLY at the solved azimuth -- the vehicle\n');
    fprintf('       arrives -- and the MINIMISATION is the part that is indicative. SO THE\n');
    fprintf('       BURNOUT ENERGY IS REPORTED TWICE ABOVE, in two columns: the MINIMISED one,\n');
    fprintf('       at the seed azimuth, which is the number the valley evidence is\n');
    fprintf('       differenced against, and the FLOWN one, which is what the vehicle left\n');
    fprintf('       burnout with after the cutoff moved by %+.6f of the burn. Neither column\n', ...
            me.dCutFrac);
    fprintf('       may stand in for the other: the flown energy is not a minimum of anything,\n');
    fprintf('       and the minimised one is not what was flown.\n');
    fprintf('       WHAT IS NOT DONE is re-minimising the loft angle at the solved azimuth. The\n');
    fprintf('       valley is flat at its bottom -- tolLoftMEDeg is %.3f deg against an energy\n', ...
            tolLoftMEDeg);
    fprintf('       resolved to parts in a thousand -- so the answer moves by far less than the\n');
    fprintf('       cost, which is a whole inner solve inside every outer evaluation. The size\n');
    fprintf('       of the approximation is the %+.1f J/kg between the two columns.\n', ...
            me.dEpsBoJkg);
    end
    fprintf('    4. THE LOFT ANGLE IS A COMMANDED ATTITUDE, NOT THE BURNOUT GAMMA. The %.1f deg\n', ...
            alphaMax);
    fprintf('       clamp on angle of attack limits how fast the flight path can be pushed\n');
    fprintf('       over, so the %.4f deg commanded here produced a %.4f deg burnout flight\n', ...
            rad2deg(pick.loftR),rad2deg(pick.gamBoR));
    fprintf('       path. The clamp also decides WHERE the two branches lie and how far the\n');
    fprintf('       depressed one reaches: on this run the depressed arc spans %.3f to\n', ...
            dep.bandLo./1000);
    fprintf('       %.3f km, and the required %.3f km has to be inside that to be flown on\n', ...
            dep.bandHi./1000,rngReq./1000);
    fprintf('       it. THE CLAMP IS A VEHICLE LIMIT AND IS NOT SET HERE: %s carries it as\n', ...
            func2str(vehicleFn));
    fprintf('       alphaMaxDeg and BM/run_ballistic reads the same one. It is a PLACEHOLDER\n');
    fprintf('       awaiting a qualification basis, and it must not be raised to bring a\n');
    fprintf('       target inside a branch -- which is what the 12 deg this script used to\n');
    fprintf('       ship had been chosen to do, at a cost of about 156 km of maximum range.\n');
    if ~isempty(overrideOf(opts,'alphaMax',[]))
    fprintf('       THIS RUN OVERRODE THE VEHICLE with alphaMax = %.4f deg, so nothing above\n', ...
            alphaMax);
    fprintf('       describes %s as it is defined.\n',func2str(vehicleFn));
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
 info.branchByApogee = byApoName;
   info.branchByTime = byTimName;
    info.atMaxRange = atMaxRange;
   info.coalescedRq = coalescedRq;
      info.rngTopM  = rngTopM;
     info.alphaMaxDeg = alphaMax;
       info.rngReqM = rngReq;
       info.rngAchM = pick.rngM;
         info.missM = missM;
        info.residM = resM;
         info.downM = downM;
       info.xTrackM = xTrackM;
      info.missHypM = missHypM;
      info.psiSeed  = psiSeed;
     info.psiLaunch = psiLaunch;
      info.dPsiAim  = wrapPi(psiLaunch - psiSeed);
        info.psiFly = psiFly;
         info.dPsiR = dPsiR;
    info.seedMissM  = seedMissM;
    info.seedDownM  = seedDwnM;
   info.seedXTrackM = seedXtkM;
      info.tolAimM  = tolAimM;
      info.aimIter  = av.iterations;
      info.aimEval  = av.nEval;
      info.aimInfo  = av;
      info.cutFrac  = cutSol;
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
        info.stopOK = allOK;
       info.stopWhy = {why1;why2;why3};

%% The minimum-energy record and the classical vacuum equal-radius arc it is
%% reported AGAINST, always present so a caller does not have to test the mode
%% before reading them. On a 'lofted' or 'depressed' run minEnergy.solved is
%% false and every other field of it is absent -- there was no solve -- while
%% classical is populated regardless, because it depends on the GEOMETRY and not
%% on what was flown. The three me*RelE / me*ResR fields are DIAGNOSTICS against
%% that reference and not residuals of anything this script drove to zero.
%%
%% minEnergy CARRIES TWO OPERATING POINTS AND THE PREFIX SAYS WHICH. The seed*
%% fields -- seedEpsBo, seedCutFrac, seedTCutS, seedRngM -- are the
%% MINIMISATION'S OWN answer, at the seed azimuth, and they are the basis the
%% valley evidence epsNeighLo, epsNeighHi, epsDepthJkg, epsNoiseJkg and
%% epsNoiseRel is differenced against. The unprefixed epsBo, cutFrac, tCutS and
%% rngM are what the vehicle FLEW, after the two-axis solve moved the cutoff to
%% close the cross-range. dCutFrac and dEpsBoJkg are the gap, and they are
%% exactly zero whenever the aim solve is a no-op. A caller comparing epsBo
%% against epsNeighLo is comparing two azimuths and will be wrong by the size of
%% dEpsBoJkg; that is what the seed prefix exists to prevent:
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
%% angle, which is br.nProp and not br.nEval. propOne is everything stage 1
%% spent; the two-axis solve and the single re-flight of its answer are what
%% stage 2 adds:
      info.propOne  = propOne;
        info.nProp  = propOne + av.nEval + 1;

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
             hProf = coorbital.viz.profilePlot(traj,bst,env, ...
            struct('Name','Ballistic targeting profile', ...
                   'Channels',{{'altitude','speed','mach','q','nAero','gamma'}}, ...
                   'Extra',{{nSens,'sensed load factor (g)', ...
                             'Sensed load factor, thrust included'}}, ...
                   'VehPhase',{{bst,coastVeh,coastVeh}}));
             hTrak = coorbital.viz.groundTrack(traj,bst,env, ...
            struct('Target',[latTargetR; lonTargetR], ...
                   'PhaseName',{{'boost','coast','descent'}}, ...
                   'Title',sprintf(['Ground track, %s arc, %.0f km required, ' ...
                                    '%.0f m miss (loft %.2f deg)'], ...
                                   pickName,rngReq./1000,missM,rad2deg(pick.loftR))));
             hGlob = coorbital.viz.globe3D(traj,bst,env, ...
            struct('Target',[latTargetR; lonTargetR], ...
                   'AltScale',altExag, ...
                   'Title',sprintf('Solved %s arc, %s to %s', ...
                                   pickName,info.launchStr,info.targetStr)));

%% Saved to disk when a stem is set, which it is by default. The suffixes name
%% what each figure IS rather than paraphrasing its title: a title carries a
%% miss distance, a minus sign and a degree symbol, and sanitising that into a
%% file name gives a different name on every run and a truncated one at that.
%% Fixed suffixes mean a re-run overwrites. An empty stem is the off switch;
%% there is deliberately no second flag beside showPlots:
        if ~isempty(plotFile)
           figFile = {coorbital.viz.saveFigure(hProf,[plotFile '_profile']), ...
                      coorbital.viz.saveFigure(hTrak,[plotFile '_ground_track']), ...
                      coorbital.viz.saveFigure(hGlob,[plotFile '_globe'])};
            fprintf('  Figures: %d PNG files written to\n',numel(figFile));
            for kf = 1:numel(figFile)
                fprintf('    %s\n',figFile{kf});
            end
            fprintf('\n');
     info.plotFiles = figFile;
        end
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

function [rngM,traj] = flyLoft(psiL,loftR,cutFrac,cfg,rI)
%% Purpose:
%
%  Fly the whole boost-coast-descent chain at one launch azimuth, one loft angle
%  and one cutoff fraction and return the great-circle surface range from the
%  launch point to the impact point. This is the function the max-range search
%  maximises, the function coorbital.util.rangeSolve bisects on at BOTH levels of
%  the minimum-energy solve, the function flyMiss builds the two-axis residual
%  from, and the only place any of them touches the physics.
%
%  THE AZIMUTH IS WRITTEN INTO THE STATE, NOT REBUILT AROUND IT. cfg.x0 is a
%  launch-state TEMPLATE and entry 6 is the heading, so a trial flight differs
%  from the seed in the azimuth and in nothing else. Assembling a fresh launch
%  state here would be a second copy of the seven lines that build it, free to
%  drift from the first.
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
%  EACH PHASE IS CHECKED AGAINST THE EVENT IT WAS SUPPOSED TO END ON, which
%  three phase labels plus a final altitude do not establish. A coast truncated
%  by its horizon still starts phase 3 and phase 3 still hits the altitude
%  event, so the old test -- three phases present, final altitude on the stop --
%  passed a physically malformed flight straight into every root solve, and the
%  later termination diagnostics run only on the SELECTED trajectory and cannot
%  protect the optimiser. The three claims now made are:
%
%    boost   ended at the commanded cutoff time when one was commanded, and at
%            the burnout MASS when the burn ran to exhaustion;
%    coast   ended at apogee, gamma through zero, with horizon to spare;
%    impact  ended on the stop altitude, DESCENDING.
%
%% Inputs:
%
%  psiL             [1 x 1]                     Launch azimuth (rad), clockwise
%                                               from north
%
%  loftR            [1 x 1]                     Commanded terminal pitch
%                                               attitude (rad)
%
%  cutFrac          [1 x 1]                     Thrust-termination time as a
%                                               fraction of the full burn (-);
%                                               1 commands no cutoff
%
%  cfg              Struct                      Immutable configuration; see
%                                               buildPhases, plus x0 [7 x 1]
%                                               (the launch-state TEMPLATE),
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
%  Michael Casey  Azimuth is now an argument, not a fixed
%                 property of cfg.x0                            08/08/2026
%  Copyright 2026 Coorbital, Inc.
%% ------------------------ Begin Code Sequence ---------------------------

              xLau = cfg.x0;
           xLau(6) = psiL;
                ph = buildPhases(loftR,cutFrac,cfg);
              traj = coorbital.prop.phaseRun(ph,xLau,cfg.bst,cfg.env);

%% All three phases must have run at all, before any of their endpoints can be
%% asked about:
             hEndM = traj.x(end,1) - cfg.rE;
             nPhRn = numel(unique(traj.phaseIdx));
    if nPhRn < 3
        error('coorbital:runBallisticTarget:propagationIncomplete', ...
            ['A loft angle of %.6f deg at a cutoff fraction of %.6f on a ' ...
             '%.6f deg azimuth produced %d of 3 phases and ended at ' ...
             'h = %.3f km against a %.3f km stop altitude. The flight did not ' ...
             'complete, so its range is meaningless to the search. Two usual ' ...
             'causes: a loft angle so depressed that the vehicle burns out ' ...
             'already descending and the apogee event never fires -- raise ' ...
             'loftMin; or a cutoff so early that it never leaves the ' ...
             'atmosphere -- raise cutFracMin.'], ...
            rad2deg(loftR),cutFrac,rad2deg(psiL),nPhRn,hEndM./1000, ...
            cfg.hStop./1000);
    end

%% Each phase against the event it was supposed to end on. The budgets are the
%% event solver's, not physical margins: the burnout event lands within a
%% nanogram of the target mass, the apogee event within 1e-14 rad of zero, and
%% the altitude event within a micrometre of the stop:
            tolTimS = 1e-6;
            tolMasK = 1e-6;
            tolGamZ = 1e-9;
            tolAltM = 1e-3;
               kBOf = find(traj.phaseIdx == 1,1,'last');
              kApoF = find(traj.phaseIdx == 2,1,'last');
    if cutFrac < 1
            boostOK = abs(traj.t(kBOf) - cutFrac.*cfg.tBurn) <= tolTimS;
            boostTx = sprintf(['boost was to end on the commanded %.6f s cutoff ' ...
                               'and ended at %.6f s'],cutFrac.*cfg.tBurn,traj.t(kBOf));
    else
            boostOK = abs(traj.x(kBOf,7) - cfg.mBurnout) <= tolMasK;
            boostTx = sprintf(['boost was to end on propellant exhaustion at ' ...
                               '%.6f kg and ended at %.6f kg'], ...
                              cfg.mBurnout,traj.x(kBOf,7));
    end
            coastOK = abs(traj.x(kApoF,5)) <= tolGamZ && ...
                      (traj.t(kApoF) - traj.t(kBOf)) < cfg.tMaxCoast - tolTimS;
            coastTx = sprintf(['coast was to end at apogee and ended at gamma = ' ...
                               '%.6e deg after %.3f s of a %.1f s horizon'], ...
                              rad2deg(traj.x(kApoF,5)), ...
                              traj.t(kApoF) - traj.t(kBOf),cfg.tMaxCoast);
           impactOK = abs(hEndM - cfg.hStop) <= tolAltM && traj.x(end,5) < 0;
           impactTx = sprintf(['descent was to end DESCENDING on the %.3f km stop ' ...
                               'and ended at %.6f km with gamma = %.4f deg'], ...
                              cfg.hStop./1000,hEndM./1000,rad2deg(traj.x(end,5)));
    if ~(boostOK && coastOK && impactOK)
        error('coorbital:runBallisticTarget:propagationIncomplete', ...
            ['A loft angle of %.6f deg at a cutoff fraction of %.6f on a ' ...
             '%.6f deg azimuth produced a trajectory that did not end as ' ...
             'intended, so its range is meaningless to the search. boost %s: ' ...
             '%s. coast %s: %s. impact %s: %s. Three phase labels and a final ' ...
             'altitude do not prove a completed flight, which is why each is ' ...
             'checked against its own event.'], ...
            rad2deg(loftR),cutFrac,rad2deg(psiL),okWord(boostOK),boostTx, ...
            okWord(coastOK),coastTx,okWord(impactOK),impactTx);
    end

              rngM = rI.*coorbital.util.greatCircle(cfg.lat0,cfg.lon0, ...
                                                    traj.x(end,3),traj.x(end,2));
end

function [fMiss,traj] = flyMiss(psiL,loftR,cutFrac,cfg,aim)
%% Purpose:
%
%  The two-component residual coorbital.util.aimSolve drives to zero in stage 2:
%  fly the chain for one control triple and return the miss vector from the
%  TARGET to the impact point, resolved in the aim frame. Sibling of flyLoft and
%  built on it, differing only in what it measures at the end -- a distance from
%  the LAUNCH point there, a two-component offset from the TARGET here.
%
%  WHY A VECTOR MISS AND NOT A RANGE RESIDUAL BESIDE A CROSS-TRACK OFFSET. The
%  two carry the same information and one of them is better. Both components of
%  this one zero is exactly "the impact point IS the target", with no reference
%  to a great circle the flown track no longer follows once the Earth turns or a
%  bank is commanded; and the Jacobian is nearly diagonal, because the azimuth
%  moves the impact point almost purely crosswise and the range control almost
%  purely along the course. A range-residual pairing would keep quoting an arc
%  the vehicle is not on and would mix both controls into both residuals.
%
%  ONLY TWO OF THE THREE ARGUMENTS ARE EVER SOLVED FOR AT ONCE. The caller binds
%  the third: the two full-burn branches hold cutFrac at 1 and solve
%  (psiL,loftR), while minimum-energy holds loftR at the angle its minimisation
%  settled on and solves (psiL,cutFrac). Which two is a property of the MODE, not
%  of this function, so this one takes all three and closes over none of them.
%
%  A CUTOFF FRACTION ABOVE 1 IS REFUSED HERE AND NOT CLAMPED, and the reason is
%  that clamping would be silent. buildPhases reads any value at or above 1 as
%  "no commanded cutoff", so a Newton step that overshot past a full burn would
%  be handed back the FULL-BURN trajectory for a control it did not ask for --
%  the residual would go flat in that direction, the Jacobian column with it,
%  and the solve would stall or step wildly with nothing in the record saying
%  why. Refusing turns it into an infeasible trial, which aimSolve already knows
%  how to handle: it halves the step and tries again. The unbounded search is
%  otherwise unavoidable -- aimSolve is a root finder and has no box.
%
%% Inputs:
%
%  psiL             [1 x 1]                     Launch azimuth (rad), clockwise
%                                               from north
%
%  loftR            [1 x 1]                     Commanded terminal pitch
%                                               attitude (rad)
%
%  cutFrac          [1 x 1]                     Thrust-termination time as a
%                                               fraction of the full burn (-);
%                                               1 commands no cutoff
%
%  cfg              Struct                      Immutable configuration; see
%                                               flyLoft
%
%  aim              Struct                      The aim frame; see missVector
%
%% Outputs:
%
%  fMiss            [2 x 1]                     Down-range and cross-range
%                                               components of the miss (m); see
%                                               missVector. Throws
%                                               coorbital:runBallisticTarget:propagationIncomplete
%                                               for a control set that does not
%                                               fly, which aimSolve treats as an
%                                               infeasible trial
%
%  traj             Struct                      The trajectory itself, from
%                                               coorbital.prop.phaseRun
%
%% Revision History:
%  Michael Casey                                                08/08/2026
%  Copyright 2026 Coorbital, Inc.
%% ------------------------ Begin Code Sequence ---------------------------

    if cutFrac > 1
        error('coorbital:runBallisticTarget:cutoffAboveFullBurn', ...
            ['A cutoff fraction of %.6f asks for more burn than the booster ' ...
             'carries. It is refused rather than clamped to 1, because ' ...
             'buildPhases reads 1 as "no commanded cutoff" and a clamp would ' ...
             'silently return the full-burn trajectory for a control the ' ...
             'solver did not ask for.'],cutFrac);
    end
            [~,traj] = flyLoft(psiL,loftR,cutFrac,cfg,aim.rI);
               fMiss = missVector(traj.x(end,3),traj.x(end,2),aim);
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
%  the direction from the target to the impact point, both scaled by the sine of
%  the central angle between them. Taken as COMPONENTS they are smooth straight
%  through coincidence, which is where a bearing does not exist at all --
%  greatCircleBearing rightly REFUSES there, and a residual that threw as it
%  converged would have the Newton halving its way into an infeasible point
%  every time it succeeded.
%
%  THE MAGNITUDE COMES FROM THE HAVERSINE, NOT FROM ASIN. hypot(wNrth,wEast) is
%  sin(Delta), and asin of it would fold every central angle past 90 degrees
%  back on itself -- reachable early in a solve, when the trial controls can land
%  most of a hemisphere away. coorbital.util.greatCircle gives Delta over the
%  whole range, so the components are rescaled by Delta/sin(Delta) instead: the
%  direction is taken from the pair, the length from the haversine.
%
%  SIGNS. Down-range is positive LONG -- the impact point is beyond the target
%  along the course -- and cross-range is positive RIGHT of that course.
%
%  This is the same construction HGV/run_target uses, deliberately duplicated
%  rather than shared: the two scripts are standalone worked examples that a
%  reader is meant to be able to follow end to end, and neither imports local
%  machinery from the other. If a third caller ever wants it, it belongs in
%  +coorbital/+util and all three should take it from there.
%
%% Inputs:
%
%  latIm            [1 x 1]                     Impact latitude (rad)
%
%  lonIm            [1 x 1]                     Impact longitude (rad)
%
%  aim              Struct                      The aim frame, built once by the
%                                               caller and constant for the
%                                               whole solve:
%                                               latT [1 x 1] target latitude
%                                                    (rad)
%                                               lonT [1 x 1] target longitude
%                                                    (rad)
%                                               psiC [1 x 1] course of the
%                                                    launch-to-target great
%                                                    circle AT THE TARGET (rad),
%                                                    clockwise from north
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
%% Delta/sin(Delta). The guard is on the DIVISOR and not on the answer, and
%% sin(Delta) vanishes at both ends of the arc rather than at one:
%%
%%    AT COINCIDENCE, Delta = 0, the ratio tends to unity and the components are
%%    zero anyway, so unity is the right value and the branch is exact.
%%
%%    AT THE ANTIPODE, Delta = pi, sin(Delta) is analytically zero again and
%%    this branch would report an exact hit on a hemisphere of miss. It is left
%%    unguarded deliberately: the zero is ANALYTIC, not numerical -- sin(pi)
%%    evaluates to 1.22e-16 rather than to zero -- so sMag > 0 holds, the ratio
%%    is taken and the magnitude comes back right, with a direction that is
%%    meaningless because every great circle joins an antipodal pair. An elseif
%%    here would be a branch no double can reach, and an unreachable guard reads
%%    as protection that was never exercised. What actually keeps this problem
%%    away from the antipode is the envelope refusal above: a target that far is
%%    unreachable and never gets as far as the Newton:
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

function o = arcRecord(traj,psiL,loftR,cutFrac,cfg,rI,c,latTgR,lonTgR)
%% Purpose:
%
%  Reduce ONE flown trajectory to the scalar quantities every record in this
%  file is built from -- burnout state, apogee, impact, burnout specific energy,
%  miss -- so that the branch solves, the minimum-energy solve and the final
%  re-flight all describe an arc in the same words and cannot drift apart.
%
%  IT MEASURES; IT DOES NOT FLY. The caller supplies the trajectory, which is
%  what lets the same reduction serve a propagation made for a bisection, one
%  made for a minimisation, and the single re-flight at the solved controls.
%
%% Inputs:
%
%  traj             Struct                      Trajectory from
%                                               coorbital.prop.phaseRun
%
%  psiL             [1 x 1]                     Launch azimuth flown (rad)
%
%  loftR            [1 x 1]                     Loft angle flown (rad)
%
%  cutFrac          [1 x 1]                     Cutoff fraction flown (-)
%
%  cfg              Struct                      Immutable configuration; lat0
%                                               and lon0 are read
%
%  rI               [1 x 1]                     Impact sphere radius (m)
%
%  c                Struct                      coorbital.util.missileConst
%
%  latTgR, lonTgR   [1 x 1]                     Target coordinates (rad)
%
%% Outputs:
%
%  o                Struct                      psiL, loftR, cutFrac (controls);
%                                               tCutS (s); rngM, missM (m);
%                                               gamBoR (rad); vBoM (m/s); hBoM
%                                               (m); mCutKg (kg); hApoM (m);
%                                               tFlyS (s); vImpM (m/s); gamImR
%                                               (rad); epsBo (J/kg); traj; ok
%                                               (logical, always true)
%
%% Revision History:
%  Michael Casey                                                08/08/2026
%  Copyright 2026 Coorbital, Inc.
%% ------------------------ Begin Code Sequence ---------------------------

               kBO = find(traj.phaseIdx == 1,1,'last');
              kApo = find(traj.phaseIdx == 2,1,'last');
            o.traj = traj;
            o.psiL = psiL;
           o.loftR = loftR;
         o.cutFrac = cutFrac;
           o.tCutS = traj.t(kBO);
            o.rngM = rI.*coorbital.util.greatCircle(cfg.lat0,cfg.lon0, ...
                                                    traj.x(end,3),traj.x(end,2));
          o.gamBoR = traj.x(kBO,5);
            o.vBoM = traj.x(kBO,4);
            o.hBoM = traj.x(kBO,1) - c.rE;
          o.mCutKg = traj.x(kBO,7);
           o.hApoM = traj.x(kApo,1) - c.rE;
           o.tFlyS = traj.t(end);
           o.vImpM = traj.x(end,4);
          o.gamImR = traj.x(end,5);

%% Specific orbital energy at burnout, the quantity the minimum-energy
%% discussion in the header turns on. Measured from the flown burnout state, so
%% it carries whatever the boost actually lost to drag and to gravity:
           o.epsBo = traj.x(kBO,4).^2./2 - c.muE./traj.x(kBO,1);
           o.missM = rI.*coorbital.util.greatCircle(traj.x(end,3),traj.x(end,2), ...
                                                    latTgR,lonTgR);
              o.ok = true;
end

function [loftStarR,rngStar,mr] = maxRangeLoft(fRange,loftLoR,loftHiR,nScan,tolLoftR)
%% Purpose:
%
%  Locate the loft angle of MAXIMUM range on the bracket, which is what splits
%  the loft axis into a lofted and a depressed branch. Three stages, and the
%  first two are what make the third safe:
%
%    A COARSE SCAN of nScan equally spaced loft angles. Golden-section search
%    requires a UNIMODAL bracket and cannot detect that it was not given one, so
%    the scan has to supply the evidence rather than the hope.
%
%    ADAPTIVE REFINEMENT until the evidence is conclusive. Unimodality is
%    CERTIFIED, not counted: the sampled range must rise STRICTLY to the best
%    node and fall STRICTLY after it. One direction change on a coarse grid does
%    not establish that -- flat runs and extrema hidden between samples both
%    survive it -- so where the scan cannot certify, midpoints are interleaved
%    into every interval and the test is repeated, up to maxRefine times. Every
%    earlier sample is reused exactly, so a refinement costs only its new
%    points. Only a refinement that STILL cannot certify is an error; erroring
%    on the first coarse scan would be brittle.
%
%    A GOLDEN-SECTION SEARCH on the three-point bracket around the best node,
%    which needs no derivative -- there is none available, each evaluation being
%    a trajectory propagation -- and shrinks the bracket by a fixed factor per
%    evaluation while reusing one of the two interior points every step.
%
%  WHAT IT RETURNS IS AN INTERVAL, AND THAT IS THE POINT. Golden section proves
%  only that the maximiser lies in the final bracket [aL,bL]; it does not
%  produce a maximiser. The midpoint is returned as loftStarR for reporting, but
%  mr.aL and mr.bL are what the caller must split the loft axis on, because
%  [loftMin,aL] and [bL,loftMax] are the intervals on which monotonicity is
%  CERTIFIED. Splitting at the midpoint instead gives one bracket that straddles
%  the true maximum, and bisection on it can converge to the wrong root or
%  reject a reachable target. mr.fAL and mr.fBL are the ranges achieved at the
%  two ends, which the caller needs to say where the branches stop being
%  distinguishable.
%
%  IT ERRORS RATHER THAN RETURNING AN UNCONVERGED MAXIMISER. If the step cap is
%  reached with the bracket still wider than tolLoftR, both branch brackets
%  would be wrong by an unstated amount, so that is a failure and not a result.
%
%  THE MAXIMUM MUST BE INTERIOR. A coarse-scan maximum sitting on either
%  endpoint means the hump is outside the bracket, and then there is no
%  two-branch structure to solve at all: range is monotonic across the whole
%  bracket and one of the two branches is empty everywhere. That is not a
%  numerical near-miss to be nudged and refinement cannot cure it, so it errors
%  at once.
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
%                                               scanLoftR [1 x nScan] (rad),
%                                                         AFTER refinement
%                                               scanRngM  [1 x nScan] (m)
%                                               nScan     [1 x 1] final count
%                                               nRefine   [1 x 1] refinement
%                                                         rounds taken
%                                               nGolden   [1 x 1] golden steps
%                                               nEval     [1 x 1] total calls
%                                                         to fRange
%                                               aL, bL    [1 x 1] the CERTIFIED
%                                                         maximiser interval
%                                                         (rad)
%                                               fAL, fBL  [1 x 1] ranges
%                                                         achieved at aL and
%                                                         bL (m)
%                                               widthR    [1 x 1] bL - aL (rad)
%                                               nTurn     [1 x 1] direction
%                                                         changes in the final
%                                                         scan; 1 is a single
%                                                         hump
%
%% References:
%   [1] Press, W.H., et al., "Numerical Recipes," 3rd ed., Cambridge, 2007,
%       Section 10.2. Golden-section search in one dimension.
%
%% Revision History:
%  Michael Casey                                                08/07/2026
%  Michael Casey  Certified interval, adaptive refinement,      08/08/2026
%                 and a non-convergence error
%  Copyright 2026 Coorbital, Inc.
%% ------------------------ Begin Code Sequence ---------------------------

%% The coarse scan:
          scanLoftR = linspace(loftLoR,loftHiR,nScan);
           scanRngM = zeros(1,nScan);
    for ks = 1:nScan
       scanRngM(ks) = fRange(scanLoftR(ks));
    end

%% Certify a single interior hump, refining until the evidence resolves. The
%% test is STRICT monotonicity either side of the best node, which a mere count
%% of direction changes does not give: a flat run reads as no turn at all, and
%% an extremum between two samples reads as none either. A refinement round
%% interleaves the midpoint of every interval, so every earlier evaluation is
%% reused exactly and the round costs only its new points:
          maxRefine = 4;
            nRefine = 0;
              [~,kB] = max(scanRngM);
    checkInterior(kB,numel(scanLoftR),scanLoftR,loftLoR,loftHiR);
    while ~unimodal(scanRngM,kB) && nRefine < maxRefine
            nRefine = nRefine + 1;
               nOld = numel(scanLoftR);
              newLo = zeros(1,2.*nOld - 1);
              newRn = zeros(1,2.*nOld - 1);
        newLo(1:2:end) = scanLoftR;
        newRn(1:2:end) = scanRngM;
        for ks = 1:(nOld - 1)
              xMidR = scanLoftR(ks)./2 + scanLoftR(ks+1)./2;
          newLo(2.*ks) = xMidR;
          newRn(2.*ks) = fRange(xMidR);
        end
          scanLoftR = newLo;
           scanRngM = newRn;
             [~,kB] = max(scanRngM);
        checkInterior(kB,numel(scanLoftR),scanLoftR,loftLoR,loftHiR);
    end
              nScanF = numel(scanLoftR);
    if ~unimodal(scanRngM,kB)
        error('coorbital:runBallisticTarget:rangeNotUnimodal', ...
            ['The sampled range curve does not certify two monotone branches, ' ...
             'and %d refinement(s) to %d points across the %.4f to %.4f deg ' ...
             'loft bracket did not resolve it. A golden-section search needs a ' ...
             'unimodal bracket and cannot detect that it was not given one, ' ...
             'and the two branch bisections below need range to be monotone on ' ...
             'their own side of the maximum -- so continuing would return a ' ...
             'converged answer on an unstated branch. Raise nScanLoft, narrow ' ...
             'the loft bracket around the hump, or look at ' ...
             'info.maxRange.scanRngM to see what the curve is actually ' ...
             'doing.'],nRefine,nScanF,rad2deg(loftLoR),rad2deg(loftHiR));
    end

%% Direction changes in the certified curve, reported for the record:
              dSign = sign(diff(scanRngM));
              dSign = dSign(dSign ~= 0);
              nTurn = sum(diff(dSign) ~= 0);

%% Golden-section search on the three-point bracket around the best scan node.
%% gRat is the reciprocal golden ratio; the two interior points are placed so
%% that whichever half survives already contains one of them. fA and fB track
%% the ranges at the LIVE bracket ends, which is free -- every value the bracket
%% ever takes has already been evaluated -- and is what lets the caller say
%% where the two branches stop being distinguishable:
                 aL = scanLoftR(kB-1);
                 bL = scanLoftR(kB+1);
                 fA = scanRngM(kB-1);
                 fB = scanRngM(kB+1);
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
                 fA = f1;
                 x1 = x2;
                 f1 = f2;
                 x2 = aL + gRat.*(bL - aL);
                 f2 = fRange(x2);
        else
                 bL = x2;
                 fB = f2;
                 x2 = x1;
                 f2 = f1;
                 x1 = bL - gRat.*(bL - aL);
                 f1 = fRange(x1);
        end
    end

%% AN UNCONVERGED MAXIMISER IS A FAILURE, NOT A RESULT. Both branch brackets are
%% built from aL and bL, so a bracket wider than asked for makes both of them
%% wrong by an amount nothing downstream knows about:
    if (bL - aL) > tolLoftR
        error('coorbital:runBallisticTarget:maxRangeNoConvergence', ...
            ['Golden-section search reached %d steps with the maximiser ' ...
             'interval still %.16g rad wide against a tolerance of %.16g rad. ' ...
             'The two branch brackets are built from that interval, so ' ...
             'returning it would make both of them wrong by an unstated ' ...
             'amount. Loosen tolLoftDeg.'],maxStep,bL - aL,tolLoftR);
    end
          loftStarR = aL./2 + bL./2;
            rngStar = fRange(loftStarR);

%% The located maximum must beat both bracket endpoints, or the search has
%% wandered off a curve that is not the hump the scan showed:
    if rngStar < max(scanRngM(1),scanRngM(nScanF))
        error('coorbital:runBallisticTarget:maximumNotFound', ...
            ['The golden section settled at %.6f deg with %.4f km of range, ' ...
             'which is less than the %.4f km reached at a bracket endpoint. ' ...
             'The range curve is not the single hump the coarse scan showed.'], ...
            rad2deg(loftStarR),rngStar./1000, ...
            max(scanRngM(1),scanRngM(nScanF))./1000);
    end

      mr.scanLoftR  = scanLoftR;
       mr.scanRngM  = scanRngM;
          mr.nScan  = nScanF;
       mr.nRefine   = nRefine;
        mr.nGolden  = nGolden;
         mr.nEval   = nScanF + 2 + nGolden + 1;
             mr.aL  = aL;
             mr.bL  = bL;
            mr.fAL  = fA;
            mr.fBL  = fB;
         mr.widthR  = bL - aL;
          mr.nTurn  = nTurn;
end

function tf = unimodal(fVals,kB)
%% Purpose:
%
%  Certify that a sampled curve rises STRICTLY to its best node and falls
%  STRICTLY after it, which is what "unimodal on this grid" actually means. A
%  count of direction changes is not the same test: a flat run contributes no
%  change and hides a plateau, and a pair of extrema between two samples
%  contributes none either.
%
%% Inputs:
%
%  fVals            [1 x n]                     Sampled values
%
%  kB               [1 x 1]                     Index of the largest, interior
%
%% Outputs:
%
%  tf               [1 x 1] logical             True when the samples certify a
%                                               single interior maximum
%
%% Revision History:
%  Michael Casey                                                08/08/2026
%  Copyright 2026 Coorbital, Inc.
%% ------------------------ Begin Code Sequence ---------------------------

                tf = all(diff(fVals(1:kB)) > 0) && all(diff(fVals(kB:end)) < 0);
end

function checkInterior(kB,nPts,scanLoftR,loftLoR,loftHiR)
%% Purpose:
%
%  Refuse a coarse-scan maximum that sits ON an endpoint of the loft bracket.
%  Refinement cannot cure that -- the hump is outside the interval being
%  searched, range is monotone across all of it, and one of the two branches is
%  empty everywhere on it -- so it is raised at once rather than after four
%  rounds of pointless propagation.
%
%% Inputs:
%
%  kB               [1 x 1]                     Index of the largest sample
%
%  nPts             [1 x 1]                     Number of samples
%
%  scanLoftR        [1 x nPts]                  The sampled loft angles (rad)
%
%  loftLoR, loftHiR [1 x 1]                     The bracket ends (rad)
%
%% Outputs:
%
%  none                                         Throws when the maximum is not
%                                               interior
%
%% Revision History:
%  Michael Casey                                                08/08/2026
%  Copyright 2026 Coorbital, Inc.
%% ------------------------ Begin Code Sequence ---------------------------

    if kB == 1 || kB == nPts
        error('coorbital:runBallisticTarget:maximumNotBracketed', ...
            ['The largest range over the loft bracket %.4f to %.4f deg was ' ...
             'found AT the %.4f deg endpoint, so the max-range loft angle is ' ...
             'outside THIS BRACKET and range is monotonic across all of it. ' ...
             'There is then only ONE branch here, and the lofted/depressed ' ...
             'structure this script is built around does not exist. Two usual ' ...
             'causes, in order of likelihood: (1) THE BRACKET IS TOO NARROW ' ...
             'and the hump lies outside it -- lower loftMin or raise loftMax. ' ...
             'At the vehicle''s 6 deg clamp the maximum is at a commanded ' ...
             '-42.9 deg, so a loftMin of -40 misses it and the shipped -140 ' ...
             'finds it. (2) The angle-of-attack clamp is too tight for the ' ...
             'vehicle to be pushed past the max-range attitude at any loft ' ...
             'angle the bracket allows. That clamp is a VEHICLE limit, in ' ...
             'BM/vehicle_bm.m, and widening the bracket is the fix to try ' ...
             'first: at 6 deg the depressed arc spans only about 4708 to ' ...
             '5212 km even once the hump is found.'], ...
            rad2deg(loftLoR),rad2deg(loftHiR),rad2deg(scanLoftR(kB)));
    end
end

function s = okWord(tf)
%% Purpose:
%
%  Render a logical as the word a failure message wants, so flyLoft's
%  completion diagnostic reads as a sentence rather than as three ones and
%  zeros.
%
%% Inputs:
%
%  tf               [1 x 1] logical             The flag
%
%% Outputs:
%
%  s                Char [1 x n]                'OK' or 'FAILED'
%
%% Revision History:
%  Michael Casey                                                08/08/2026
%  Copyright 2026 Coorbital, Inc.
%% ------------------------ Begin Code Sequence ---------------------------

                 s = 'FAILED';
    if tf
                 s = 'OK';
    end
end

function name = branchOfLoft(loftR,aL,bL,tolLoftR)
%% Purpose:
%
%  Name the branch a loft angle lies on, from its POSITION relative to the
%  certified maximiser interval [aL,bL]. That is the only branch invariant this
%  search has: below aL the angle is on the depressed side of every point the
%  maximum could occupy, above bL it is on the lofted side, and inside it no
%  evidence available here can tell the two arcs apart.
%
%  It replaces a test on the flown apogee and flight time against the max-range
%  arc's. Neither of those is a branch invariant for a finite powered
%  atmospheric arc, and near the maximum both differences vanish quadratically,
%  so their signs are set by search and integration error rather than by the
%  branch.
%
%% Inputs:
%
%  loftR            [1 x 1]                     Loft angle in question (rad)
%
%  aL, bL           [1 x 1]                     Certified maximiser interval
%                                               (rad)
%
%  tolLoftR         [1 x 1]                     Slack either side, the same
%                                               tolerance the interval was
%                                               closed to (rad)
%
%% Outputs:
%
%  name             Char [1 x n]                'depressed', 'lofted' or
%                                               'coalesced'
%
%% Revision History:
%  Michael Casey                                                08/08/2026
%  Copyright 2026 Coorbital, Inc.
%% ------------------------ Begin Code Sequence ---------------------------

              name = 'coalesced';
    if loftR < aL - tolLoftR
              name = 'depressed';
    elseif loftR > bL + tolLoftR
              name = 'lofted';
    end
end

function br = solveBranch(name,psiL,rngReq,fRange,loftAR,loftBR,tolRngM,cfg, ...
                          rI,c,latTgR,lonTgR,aMaxR,bMaxR,tolLoftR)
%% Purpose:
%
%  Solve ONE branch of the two-branch loft problem and gather everything the
%  summary reports about it. The bracket handed in lies entirely on one side of
%  the CERTIFIED MAXIMISER INTERVAL -- [loftMin,aL] or [bL,loftMax], never a
%  split at the interval's midpoint -- which is what makes range monotonic on it
%  and bisection safe. See the MONOTONICITY IS ASSUMED, NOT CHECKED note in
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
%  psiL             [1 x 1]                     Launch azimuth (rad) fRange is
%                                               bound to, carried so the kept
%                                               record says which one it was
%                                               flown at. It is the SEED
%                                               azimuth; the cross-range half of
%                                               the aim is closed afterwards
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
%  aMaxR, bMaxR     [1 x 1]                     The CERTIFIED maximiser
%                                               interval (rad); the branch this
%                                               arc measures as is its loft
%                                               angle's position against it
%
%  tolLoftR         [1 x 1]                     Slack either side of that
%                                               interval (rad)
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

%% One propagation at the converged loft angle, kept, and reduced by the same
%% arcRecord the minimum-energy solve and the final re-flight use:
      [rngFly,traj]   = fRange(loftR);
    assert(abs(rngFly - rngAch) < 1e-6, ...
        ['re-flying the solved %s loft angle gave %.6f m of range against the ' ...
         'solver''s %.6f m; the propagation is not repeatable and nothing ' ...
         'below can be trusted'],name,rngFly,rngAch);
           br.nProp   = sv.nEval + 1;
           br.why     = 'reaches the target';
                br    = mergeInto(br,arcRecord(traj,psiL,loftR,1,cfg,rI,c, ...
                                               latTgR,lonTgR));
        br.measured   = branchOfLoft(loftR,aMaxR,bMaxR,tolLoftR);
end

function me = minEnergySolve(psiL,rngReq,cfg,rI,c,latTgR,lonTgR,dep,lof, ...
                             loftLoR,loftHiR,tolRngMEM,tolLoftMER, ...
                             cutFracMin,nScanME)
%% Purpose:
%
%  Solve the minimum-energy ballistic trajectory as the CONSTRAINED
%  MINIMISATION it is:
%
%      minimise    eps_BO = V_BO^2/2 - mu/r_BO
%      over        the loft angle and the cutoff fraction
%      subject to  R(loft,cutFrac) = rngReq.
%
%  Two parameters and one constraint leave a ONE-DIMENSIONAL feasible family.
%  Parameterise it by the loft angle, let the cutoff fraction be whatever
%  enforces the constraint there, and minimise the burnout specific energy along
%  what is left. That is a real objective with a real minimum, and it is what
%  makes the mode's name true.
%
%  IT IS NOT A GAMMA MATCH. The previous implementation drove the burnout
%  flight-path angle to gammaStar = 45 deg - Lambda/4 with Lambda the
%  PAD-TO-TARGET central angle. That closed form is derived for a free-flight
%  arc with both endpoints on the SAME radius, and this burnout is downrange and
%  tens of kilometres up, so the residual did not apply to the arc being flown
%  and driving it to zero verified the wrong condition. gammaStar is still
%  computed by the caller and printed beside the achieved burnout gamma, as a
%  vacuum equal-radius yardstick.
%
%  THE INNER PROBLEM IS THE CONSTRAINT, and it is bracketed honestly:
%
%    The cutoff axis [cutFracMin, 1] is SAMPLED and a SIGN CHANGE of
%    R(loft,cutFrac) - rngReq is looked for. More burn does not have to mean
%    more range at fixed loft -- added speed can carry the trajectory through
%    its own range maximum, and added burn time also moves the burnout position,
%    altitude, flight-path angle, gravity loss and drag loss -- so the interval
%    is not assumed to hold one monotone root. Monotonicity is VERIFIED on the
%    interval actually selected, and where several roots exist the one with the
%    LOWEST burnout energy is taken, which is the declared objective.
%
%  THE OUTER PROBLEM IS A MINIMISATION, and its bracket is physical: the two
%  ends of the feasible family are the two FULL-BURN branch solutions, where the
%  burn is not cut at all and the burnout energy is therefore at its highest, so
%  a valley in between is what the physics predicts. A coarse scan of nScanME
%  points certifies a single valley -- refining once by interleaving midpoints
%  if it cannot -- and a golden-section minimisation then closes the loft
%  bracket to tolLoftMER.
%
%  WHERE A BRANCH DOES NOT EXIST the caller's own loft limit is used at that
%  end, which is admissible for precisely the reason the branch is missing: the
%  full burn OVERSHOOTS the target there, so a cut-short burn can still be made
%  to land on it. The energy at that end is then already below full-burn, and
%  the minimum may sit AT the end rather than inside -- which the golden section
%  handles, a monotone function being unimodal.
%
%  THE INNER TOLERANCE IS NOISE ON THE OUTER OBJECTIVE, and it is measured. At
%  the settled loft angle the constraint is re-solved with a tolerance ten times
%  tighter and the change in eps_BO is reported as epsNoiseJkg, against the
%  shallower side of the energy valley as epsNoiseRel. Nothing here asserts that
%  the noise is small; the ratio is handed back and printed.
%
%% Inputs:
%
%  psiL             [1 x 1]                     Launch azimuth (rad) every
%                                               propagation this minimisation
%                                               makes is flown on. It is the
%                                               SEED azimuth: the cross-range
%                                               half of the aim is closed
%                                               afterwards, by the two-axis
%                                               solve, over the cutoff fraction
%                                               and the azimuth at the loft
%                                               angle this returns
%
%  rngReq           [1 x 1]                     Required surface range (m)
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
%  tolRngMEM        [1 x 1]                     Tolerance of the inner
%                                               feasibility solve, on the
%                                               achieved range (m)
%
%  tolLoftMER       [1 x 1]                     Bracket width the energy
%                                               minimisation closes the loft
%                                               angle to (rad)
%
%  cutFracMin       [1 x 1]                     Floor of the sampled cutoff
%                                               interval, as a fraction of the
%                                               full burn (-)
%
%  nScanME          [1 x 1]                     Coarse scan points along the
%                                               feasible family, at least 3
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
%                                               scanLoftR [1 x n] (rad) and
%                                               scanEps   [1 x n] (J/kg), the
%                                                           coarse scan of the
%                                                           feasible family
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
%                                               epsBo       (J/kg) THE MINIMISED
%                                                           OBJECTIVE
%                                               epsNeighLo, epsNeighHi (J/kg)
%                                                           burnout energy at
%                                                           the two neighbouring
%                                                           FEASIBLE points
%                                               epsNoiseJkg (J/kg) change in
%                                                           epsBo when the
%                                                           constraint is
%                                                           re-solved ten times
%                                                           tighter
%                                               epsNoiseRel (-) that noise over
%                                                           the shallower side
%                                                           of the valley
%                                               mCutKg,propLeftKg (kg)
%                                               nScan, nGolden, nRefine
%                                               widthR      (rad) final bracket
%                                               nProp       propagations spent
%                                               traj        struct
%
%% References:
%   [1] Bate, R.R., Mueller, D.D., White, J.E., "Fundamentals of
%       Astrodynamics," Dover, 1971, Ch. 6.
%   [2] Press, W.H., et al., "Numerical Recipes," 3rd ed., Cambridge, 2007,
%       Sections 9.1 and 10.2. Bisection for the constraint, golden-section
%       minimisation for the objective.
%
%% Revision History:
%  Michael Casey                                                08/08/2026
%  Michael Casey  Rewritten as an explicit constrained          08/08/2026
%                 minimisation of the burnout energy
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
      me.scanLoftR = [];
        me.scanEps = [];
          me.nProp = 0;
          me.nScan = 0;
        me.nGolden = 0;
        me.nRefine = 0;

%% A collapsed bracket means the required range IS the maximum range to within
%% the range tolerance: the two branches have met, there is one loft angle left
%% and no family to minimise along. Refused rather than searched on an interval
%% of zero width:
    if ~(loftAR < loftBR)
            me.why = sprintf(['The two full-burn arcs have MERGED at a loft angle of ' ...
                              '%.4f deg, so there is no feasible family left to ' ...
                              'minimise along. The required range is within the range ' ...
                              'tolerance of the MAXIMUM range, where the lofted and ' ...
                              'depressed arcs are the same arc; tighten tolRangeKm or ' ...
                              'move the target inside the envelope.'],rad2deg(loftAR));
        return;
    end

             nProp = 0;

%% The INNER problem, and the only place this solve touches the physics. It
%% returns the feasible cutoff fraction at a held loft angle, or an empty
%% record when no cutoff on the sampled interval makes the range:
    function o = feasibleAt(loftR)
                 o = solveCutoff(loftR,rngReq,cutFracMin,tolRngMEM, ...
                                 @flyMeasure,@flyCount);
    end

%% The OUTER objective. Infeasibility is +Inf rather than an error, so a scan
%% point that cannot be made to reach the target simply loses:
    function e = epsAt(loftR)
                 o = feasibleAt(loftR);
                 e = Inf;
        if o.ok
                 e = o.epsBo;
        end
    end

%% Range only, for the sign-change sampling and the bisection. Every propagation
%% in this minimisation is flown at the SEED azimuth psiL:
    function rngM = flyCount(loftR,cutF)
             nProp = nProp + 1;
              rngM = flyLoft(psiL,loftR,cutF,cfg,rI);
    end

%% One propagation, kept, and reduced to the quantities the record needs:
    function o = flyMeasure(loftR,cutF)
             nProp = nProp + 1;
             [~,trj] = flyLoft(psiL,loftR,cutF,cfg,rI);
                 o = arcRecord(trj,psiL,loftR,cutF,cfg,rI,c,latTgR,lonTgR);
    end

%% THE COARSE SCAN of the feasible family. It does three jobs: it finds a
%% bracket for the minimisation, it CERTIFIES that the energy has a single
%% valley rather than several, and it supplies the two neighbouring feasible
%% energies the answer is reported against:
         scanLoftR = linspace(loftAR,loftBR,nScanME);
           scanEps = zeros(1,nScanME);
    for ks = 1:nScanME
        scanEps(ks) = epsAt(scanLoftR(ks));
    end
%% Certify a CONTIGUOUS feasible run carrying a SINGLE VALLEY, refining once by
%% interleaving midpoints if the coarse samples cannot.
%%
%% THE ENDS OF THE INTERVAL ARE ALLOWED TO BE INFEASIBLE, and that is not
%% slackness. The outer bracket ends are the two full-burn branch solutions,
%% which coorbital.util.rangeSolve places within tolRangeKm of the range --
%% EITHER SIDE of it. Where a branch settled a few hundred metres SHORT, the
%% full burn at that exact loft angle does not quite reach the target and no
%% cutoff at or below 1 can make it, so that one endpoint is infeasible by an
%% amount that is an artefact of the branch tolerance rather than a property of
%% the family. What must hold is that the feasible samples form ONE
%% uninterrupted run of at least three points -- a HOLE in the middle would mean
%% the family is not connected and the minimisation would be searching across
%% something it cannot see:
           maxRefME = 1;
            nRefine = 0;
       [okS,iLo,iHi,kB] = scanOK(scanEps);
    while ~okS && nRefine < maxRefME
            nRefine = nRefine + 1;
               nOld = numel(scanLoftR);
              newLo = zeros(1,2.*nOld - 1);
              newEp = zeros(1,2.*nOld - 1);
        newLo(1:2:end) = scanLoftR;
        newEp(1:2:end) = scanEps;
        for ks = 1:(nOld - 1)
              xMidR = scanLoftR(ks)./2 + scanLoftR(ks+1)./2;
          newLo(2.*ks) = xMidR;
          newEp(2.*ks) = epsAt(xMidR);
        end
          scanLoftR = newLo;
            scanEps = newEp;
       [okS,iLo,iHi,kB] = scanOK(scanEps);
    end
      me.scanLoftR = scanLoftR;
        me.scanEps = scanEps;
        me.nRefine = nRefine;
        me.nScan   = numel(scanLoftR);
       me.nFeasLo  = iLo;
       me.nFeasHi  = iHi;
          me.nProp = nProp;
    if ~okS
            me.why = sprintf(['The burnout energy sampled along the loft interval ' ...
                              '%.4f to %.4f deg does not certify a connected feasible ' ...
                              'family with a SINGLE valley, and one refinement to %d ' ...
                              'points did not resolve it: %d of them could not be made ' ...
                              'to reach %.3f km at any cutoff fraction at or above ' ...
                              '%.4f of the full burn. Either those points are not a ' ...
                              'run at the ends of the interval -- a hole in the middle ' ...
                              'means the family is not connected -- or the energy has ' ...
                              'more than one valley, and a golden-section minimisation ' ...
                              'would return a local minimum without saying so. Raise ' ...
                              'nScanME and look at info.minEnergy.scanEps.'], ...
                             rad2deg(loftAR),rad2deg(loftBR),numel(scanLoftR), ...
                             sum(~isfinite(scanEps)),rngReq./1000,cutFracMin);
        return;
    end

%% The GOLDEN-SECTION MINIMISATION, on the three-point bracket around the best
%% node -- or on the one-sided interval at the end of the feasible run, where
%% the valley bottom sits at an endpoint because that end has no full-burn
%% branch and is therefore already cut short. gRat is the reciprocal golden
%% ratio:
                 kLo = max(iLo,kB - 1);
                 kHi = min(iHi,kB + 1);
                 aE = scanLoftR(kLo);
                 bE = scanLoftR(kHi);
          epsNeighA = scanEps(kLo);
          epsNeighB = scanEps(kHi);
               gRat = (sqrt(5) - 1)./2;
                 x1 = bE - gRat.*(bE - aE);
                 x2 = aE + gRat.*(bE - aE);
                 e1 = epsAt(x1);
                 e2 = epsAt(x2);
            nGolden = 0;
            maxStep = 200;
    while (bE - aE) > tolLoftMER && nGolden < maxStep
            nGolden = nGolden + 1;
        if e1 < e2
                 bE = x2;
                 x2 = x1;
                 e2 = e1;
                 x1 = bE - gRat.*(bE - aE);
                 e1 = epsAt(x1);
        else
                 aE = x1;
                 x1 = x2;
                 e1 = e2;
                 x2 = aE + gRat.*(bE - aE);
                 e2 = epsAt(x2);
        end
    end
        me.nGolden = nGolden;
         me.widthR = bE - aE;
          me.nProp = nProp;
    if (bE - aE) > tolLoftMER
            me.why = sprintf(['The energy minimisation reached %d golden-section ' ...
                              'steps with the loft bracket still %.6g deg wide against ' ...
                              'a tolerance of %.6g deg. Loosen tolLoftMEDeg.'], ...
                             maxStep,rad2deg(bE - aE),rad2deg(tolLoftMER));
        return;
    end

%% The settled point, re-solved and KEPT. The minimiser is the midpoint of the
%% final bracket, which is what the bracket certifies and nothing finer:
             loftS = aE./2 + bE./2;
              oFin = feasibleAt(loftS);
    if ~oFin.ok
            me.why = sprintf(['The settled loft angle of %.6f deg turned out to be ' ...
                              'infeasible when the constraint was re-solved there, ' ...
                              'though every scan point around it was feasible. The ' ...
                              'feasible family is not connected across this interval; ' ...
                              'raise nScanME.'],rad2deg(loftS));
        return;
    end
                me = mergeInto(me,oFin);

%% THE EVIDENCE THAT IT IS A MINIMUM. The two neighbouring FEASIBLE points --
%% same range, different loft -- must both carry more burnout energy than the
%% answer, and how much more is the depth of the valley, printed rather than
%% asserted:
     me.epsNeighLo = epsNeighA;
     me.epsNeighHi = epsNeighB;
         me.nProp  = nProp;

%% THE PROPAGATED NOISE OF THE CONSTRAINT TOLERANCE, measured rather than waved
%% away. Every outer evaluation above was an inner solve accurate only to
%% tolRngMEM, so the sampled energies carry an error whose size nothing about
%% the method predicts. Re-solving at the settled loft angle ten times tighter
%% and differencing the burnout energy is that error, and the ratio against the
%% shallower side of the valley is what says whether the minimum is resolved:
             oFine = solveCutoff(loftS,rngReq,cutFracMin,tolRngMEM./10, ...
                                 @flyMeasure,@flyCount);
    assert(oFine.ok, ...
        ['the minimum-energy constraint solved at the settled %.6f deg loft ' ...
         'angle and then failed when repeated ten times tighter; the ' ...
         'propagation is not repeatable'],rad2deg(loftS));
   me.epsNoiseJkg = abs(oFine.epsBo - me.epsBo);

%% The depth of the valley, taken over the neighbours that EXIST. A minimum at
%% the end of the family has one neighbour, not two, and differencing the best
%% node against itself would report a zero-deep valley and a caution that was
%% only an artefact of counting:
            depths = [Inf Inf];
    if kLo < kB
         depths(1) = epsNeighA - me.epsBo;
    end
    if kHi > kB
         depths(2) = epsNeighB - me.epsBo;
    end
            depthM = min(depths);
   me.epsDepthJkg = depthM;
    me.epsNoiseRel = Inf;
    if isfinite(depthM) && depthM > 0
    me.epsNoiseRel = me.epsNoiseJkg./depthM;
    end
          me.nProp = nProp;
         me.solved = true;
         me.exists = true;
      me.propLeftKg = me.mCutKg - cfg.mBurnout;
       me.hApoKepM  = keplerApogee(c.rE + me.hBoM,me.vBoM,me.gamBoR,c) - c.rE;
end

function [tf,iLo,iHi,kB] = scanOK(eVals)
%% Purpose:
%
%  Certify that the sampled objective has ONE uninterrupted run of feasible
%  points, at least three of them, carrying a SINGLE valley -- and report where
%  that run is.
%
%  INFEASIBLE POINTS AT THE ENDS ARE ALLOWED. The interval being scanned is
%  bounded by the two full-burn branch solutions, which are placed within the
%  branch range tolerance EITHER SIDE of the required range, so an endpoint can
%  be a few hundred metres short of feasible for a reason that is an artefact of
%  that tolerance. A hole in the MIDDLE is a different thing entirely: the
%  feasible family is then not connected and a bracketing minimisation across it
%  is searching something it cannot see.
%
%  The valley test tolerates a FLAT step, unlike the strict test the range hump
%  gets, because two loft angles either side of a flat minimum can return
%  burnout energies that differ by less than the propagation resolves. What it
%  must not tolerate is a RISE before the minimum or a FALL after it, which is a
%  second valley.
%
%% Inputs:
%
%  eVals            [1 x n]                     Sampled objective; +Inf marks
%                                               an infeasible point
%
%% Outputs:
%
%  tf               [1 x 1] logical             True when the samples certify a
%                                               connected feasible run of at
%                                               least three points with a
%                                               single valley
%
%  iLo, iHi         [1 x 1]                     First and last index of the
%                                               feasible run; 1 and numel(eVals)
%                                               when there is none
%
%  kB               [1 x 1]                     Index of the smallest value
%
%% Revision History:
%  Michael Casey                                                08/08/2026
%  Copyright 2026 Coorbital, Inc.
%% ------------------------ Begin Code Sequence ---------------------------

              feas = isfinite(eVals);
            [~,kB] = min(eVals);
               iLo = 1;
               iHi = numel(eVals);
                tf = false;
    if ~any(feas)
        return;
    end
               iLo = find(feas,1,'first');
               iHi = find(feas,1,'last');
                tf = all(feas(iLo:iHi)) && (iHi - iLo) >= 2 && ...
                     all(diff(eVals(iLo:kB)) <= 0) && ...
                     all(diff(eVals(kB:iHi)) >= 0);
end

function o = solveCutoff(loftR,rngReq,cutFracMin,tolRngM,fMeasure,fCount)
%% Purpose:
%
%  Enforce the range constraint at ONE held loft angle: find the cutoff fraction
%  at which the chain lands on the required range, and return the measured
%  burnout state there.
%
%  IT DOES NOT ASSUME [cutFracMin, 1] HOLDS ONE MONOTONE ROOT. More burn does
%  not have to mean more range at a fixed loft angle: added speed can carry the
%  trajectory through its own range maximum, and added burn time also moves the
%  burnout position, altitude, flight-path angle, gravity loss and drag loss.
%  The endpoint at cutFrac = 1 being at or beyond the target proves only that
%  ONE end is on the far side; it proves neither that cutFracMin undershoots nor
%  that the interval holds a unique root. So:
%
%    the interval is SAMPLED and every SIGN CHANGE of R - rngReq is collected;
%    the interval carrying each sign change is CHECKED for monotonicity on a
%      finer grid before it is bisected, and one that is not monotone is
%      skipped rather than bisected on an assumption it fails;
%    where several roots survive, the one with the LOWEST BURNOUT ENERGY is
%      returned, which is the objective the caller is minimising.
%
%  AN INFEASIBLE LOFT ANGLE IS NOT AN ERROR. It returns ok = false, because the
%  caller is scanning a family and a point it cannot reach simply loses.
%
%% Inputs:
%
%  loftR            [1 x 1]                     Held loft angle (rad)
%
%  rngReq           [1 x 1]                     Required surface range (m)
%
%  cutFracMin       [1 x 1]                     Floor of the sampled cutoff
%                                               interval (-)
%
%  tolRngM          [1 x 1]                     Tolerance on the achieved range
%                                               (m)
%
%  fMeasure         Function handle             o = fMeasure(loftR,cutFrac),
%                                               one propagation reduced to the
%                                               burnout and impact record
%
%  fCount           Function handle             rng = fCount(loftR,cutFrac),
%                                               range only
%
%% Outputs:
%
%  o                Struct                      The fMeasure record at the
%                                               chosen root, with ok = true; or
%                                               a struct with ok = false and
%                                               why populated
%
%% References:
%   [1] Press, W.H., et al., "Numerical Recipes," 3rd ed., Cambridge, 2007,
%       Section 9.1. Bracketing by sampled sign change, then bisection.
%
%% Revision History:
%  Michael Casey                                                08/08/2026
%  Copyright 2026 Coorbital, Inc.
%% ------------------------ Begin Code Sequence ---------------------------

%% The sampled cutoff axis. Five points across the interval, which is enough to
%% separate roots that a physically meaningful cutoff schedule could produce
%% while costing five propagations; the monotonicity check below adds three more
%% on whichever sub-interval is selected:
             nSamp = 5;
             nFine = 3;
           cutGrid = linspace(cutFracMin,1,nSamp);
           rngGrid = zeros(1,nSamp);
    for ks = 1:nSamp
       rngGrid(ks) = fCount(loftR,cutGrid(ks));
    end
             resid = rngGrid - rngReq;

%% An endpoint already on the target is a root, and the cheapest one:
              o.ok = false;
             o.why = '';
    for ks = 1:nSamp
        if abs(resid(ks)) <= tolRngM
                 o = fMeasure(loftR,cutGrid(ks));
            return;
        end
    end

%% Every sign change, collected. Each is a candidate root:
              cand = zeros(1,nSamp - 1);
              nCan = 0;
    for ks = 1:(nSamp - 1)
        if sign(resid(ks)) ~= sign(resid(ks+1))
              nCan = nCan + 1;
        cand(nCan) = ks;
        end
    end
              cand = cand(1:nCan);
    if isempty(cand)
             o.why = sprintf(['no cutoff fraction between %.4f and 1 of the full ' ...
                              'burn changes the sign of the range residual at a loft ' ...
                              'angle of %.6f deg; the sampled ranges span %.3f to ' ...
                              '%.3f km against a required %.3f km'], ...
                             cutFracMin,rad2deg(loftR),min(rngGrid)./1000, ...
                             max(rngGrid)./1000,rngReq./1000);
        return;
    end

%% Solve each candidate, verifying monotonicity on its own interval first. A
%% sub-interval that is not monotone is skipped: bisection there would converge
%% to one of several crossings and not say which:
              best = [];
    for kc = 1:numel(cand)
                ks = cand(kc);
              cLoF = cutGrid(ks);
              cHiF = cutGrid(ks+1);
             fineC = linspace(cLoF,cHiF,nFine + 2);
             fineR = zeros(1,nFine + 2);
        fineR(1)   = rngGrid(ks);
        fineR(end) = rngGrid(ks+1);
        for kf = 2:(nFine + 1)
          fineR(kf) = fCount(loftR,fineC(kf));
        end
             dFine = diff(fineR);
        if ~(all(dFine > 0) || all(dFine < 0))
            continue;
        end
      [cutF,~,svIn] = coorbital.util.rangeSolve(rngReq, ...
                          @(cf) fCount(loftR,cf),cLoF,cHiF,tolRngM);
        if ~svIn.converged
            continue;
        end
              oCan = fMeasure(loftR,cutF);
        if isempty(best) || oCan.epsBo < best.epsBo
              best = oCan;
        end
    end
    if isempty(best)
             o.why = sprintf(['every sign change on the cutoff interval at a loft ' ...
                              'angle of %.6f deg either failed its monotonicity check ' ...
                              'or failed to converge, so no cutoff fraction can be ' ...
                              'certified to make %.3f km there'], ...
                             rad2deg(loftR),rngReq./1000);
        return;
    end
                 o = best;
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
%  DESCRIBE a trajectory by its apogee and flight time against the max-range
%  arc's, in the vocabulary a reader recognises a lofted arc by.
%
%  IT DOES NOT CLASSIFY THE BRANCH ANY MORE, and it must not be made to. Neither
%  apogee nor flight time is a branch invariant for a finite powered
%  atmospheric arc: drag, lift, burnout altitude and boost duration can make
%  either non-monotone in the loft angle, and near the maximum both differences
%  vanish quadratically, so their signs are set by search and integration error
%  rather than by the branch. The classification comes from branchOfLoft, on the
%  root's position against the certified maximiser interval. What this is for is
%  the printed description, and the agree flag it returns is worth a caution
%  because two readings that stop agreeing say the flight is close to the
%  maximum.
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
%                                               more direct of the two.
%                                               DESCRIPTIVE
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

function w = ctlWord(isMinE,loftR,cutFrac)
%% Purpose:
%
%  Name the mode's SECOND control and give its value, in one phrase, so a
%  refusal banner can report what the two-axis solve was moving without the
%  banner itself having to know which mode it is in. The two full-burn branches
%  solve the loft angle at a burn that is never cut; minimum-energy holds the
%  loft angle and solves the cutoff fraction instead.
%
%  Written as a phrase rather than as a label-and-number pair because the two
%  cases want different units, different precisions and different companion
%  facts, and a printf format flexible enough for both would be less readable
%  than this.
%
%% Inputs:
%
%  isMinE           Logical [1 x 1]             True on the minimum-energy mode
%
%  loftR            [1 x 1]                     Loft angle (rad)
%
%  cutFrac          [1 x 1]                     Cutoff fraction (-)
%
%% Outputs:
%
%  w                Char [1 x n]                The phrase
%
%% Revision History:
%  Michael Casey                                                08/08/2026
%  Copyright 2026 Coorbital, Inc.
%% ------------------------ Begin Code Sequence ---------------------------

                 w = sprintf(['loft angle   %10.6f deg (at the full burn; the ' ...
                              'cutoff is not a control here)'],rad2deg(loftR));
    if isMinE
                 w = sprintf(['cutoff       %10.6f of full burn (at the held ' ...
                              '%.4f deg loft angle)'],cutFrac,rad2deg(loftR));
    end
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
