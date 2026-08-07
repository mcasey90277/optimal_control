function test_boost3DOF()
%% Purpose:
%
%  Analytic validation of coorbital.eom.boost3DOF, the powered equations of
%  motion. Four checks, each covering something the others cannot:
%
%    1. REDUCTION TO GLIDE. With env.prop returning a dead engine and the
%       state mass set to veh.mass, the first six components must equal
%       coorbital.eom.glide3DOF at the same state, over several states and
%       with the Earth rate both off and on.
%
%    2. THRUST FORCE INCREMENT. The same state evaluated twice, once with the
%       real engine and once with a dead one, must differ by exactly the
%       thrust projection
%
%           d(dV/dt)     = T cos(alpha)/m
%           d(dgamma/dt) = T sin(alpha) cos(sigma)/(m V)
%           d(dpsi/dt)   = T sin(alpha) sin(sigma)/(m V cos(gamma))
%
%       evaluated at a NONZERO angle of attack and a NONZERO bank angle, with
%       a state mass deliberately unequal to any mass field in either
%       parameter struct.
%
%    3. TSIOLKOVSKY. In vacuum, with gravity switched off and thrust along
%       the velocity vector, the speed gained over the whole burn must equal
%       Isp*g0*log(mLiftoff/mBurnout).
%
%    4. MASS LINEARITY. Constant mass flow makes m(t) = mLiftoff - mdot*t an
%       exact line, and the propagated mass history must lie on it.
%
%% Note -- what check 1 proves, and what it emphatically does not:
%
%  The reduction holds BIT-EXACTLY, measured residual 0.000e+00. That is not
%  luck. Every shared term in boost3DOF is a character-identical copy of the
%  glide3DOF expression evaluated in the same order, and with T = 0 the
%  thrust perturbations (0 - aDrag) and (0 + aLift) are IEEE-exact
%  identities. So the check PROVES the six shared equations ARE glide3DOF --
%  no transcription drift, no dropped term, no altered sign is possible --
%  and it is simultaneously ZERO-BIT EVIDENCE OF CORRECTNESS. Any error in
%  glide3DOF is inherited by boost3DOF and cancels invisibly here. It
%  validates IDENTITY, not PHYSICS; the physics of the shared terms is
%  validated in test_glide3DOF, not here.
%
%  It also catches NONE of the thrust mutations, because every term it could
%  catch them in vanishes under the reduction. Swapping cos(sigma) and
%  sin(sigma) on the thrust normal terms, dropping cos(alpha) from the
%  tangential term, or using the wrong mass in a thrust denominator all leave
%  the reduction bit-exact. Check 2 is the one that sees them. The general
%  form of this, recorded in docs/LESSONS_LEARNED.md: a reduction test
%  validates only the terms that SURVIVE the reduction.
%
%% Note -- the trap check 3 does not escape:
%
%  Isp and g0 appear on both sides of the Tsiolkovsky comparison -- inside
%  constThrust's mdot = thrustVac/(Isp*g0) on the simulation side, and inside
%  the closed form Isp*g0*log(mass ratio) on the reference side -- so they
%  CANCEL. A wrong Isp is invisible to this test at any tolerance. This is
%  exactly the shared-constant blindness recorded in docs/LESSONS_LEARNED.md
%  for Hscale, and the only defence is the same one: pin the constant where
%  it is defined. Isp is pinned to 260 s in test_constThrust, and g0 is
%  pinned in test_missileConst. Those pins are load-bearing for this file.
%
%% Note -- provenance of the reference values:
%
%  Every reference below is computed from coorbital.util.missileConst,
%  coorbital.util.boosterDefaults, coorbital.util.vehicleDefaults and the
%  state itself, or from a direct call to the environment's own propulsion
%  model. None is obtained by rearranging boost3DOF's output.
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
%   [1] Vinh, N.X., Busemann, A., Culp, R.D., "Hypersonic and Planetary Entry
%       Flight Mechanics," Univ. Michigan Press, 1980, Eqs. (2.28)-(2.33).
%       (The flight-mechanics equations and the thrust projection onto the
%       velocity, flight-path and heading channels.)
%   [2] Sutton, G.P., Biblarz, O., "Rocket Propulsion Elements," 9th ed.,
%       Wiley, 2017, Sec. 4.1. (The ideal rocket equation.)
%
%% Revision History:
%  Michael Casey                                                08/07/2026
%  Copyright 2026 Coorbital, Inc.
%% ------------------------ Begin Code Sequence ---------------------------

                 c = coorbital.util.missileConst();
               veh = coorbital.util.vehicleDefaults();
               bst = coorbital.util.boosterDefaults();

%% =========================================================================
%% CHECK 1 -- reduction to glide.
%%
%% ONE environment struct serves both routines, so the non-thrust models are
%% identical by construction rather than by inspection; glide3DOF simply
%% ignores the prop field. The engine is dead, and the state mass is
%% veh.mass, which is the mass glide3DOF divides by.
%% =========================================================================
      envRed.atmos = @coorbital.atmos.expAtmos;
       envRed.grav = @coorbital.grav.sphereGrav;
       envRed.aero = @coorbital.aero.constLD;
       envRed.prop = @(tp,Pp,vp) deal(0,0);
     envRed.omegaE = 0;

%% Columns: h (m), lat (deg), V (m/s), gamma (deg), psi (deg), alpha (deg),
%% sigma (deg), and a flag scaling env.omegaE, so 0 is the non-rotating
%% branch and 1 the real Earth rate. Case 1 is equatorial, eastward and
%% unbanked; case 2 is off-equator on a south-westerly heading with a nonzero
%% bank; case 3 is case 2 with the rotation terms live; case 4 is a steeply
%% climbing southern-hemisphere state, also rotating:
          redCases = [30e3,   0, 4000, -2,  90,  0,   0, 0; ...
                      45e3,  38, 5200,  3, 215,  6, -40, 0; ...
                      45e3,  38, 5200,  3, 215,  6, -40, 1; ...
                      12e3, -25, 1200, 55, 140, 10,  25, 1];

          redWorst = 0;
for kc = 1:size(redCases,1)
     envRed.omegaE = redCases(kc,8).*c.omegaE;
            xGlide = [c.rE + redCases(kc,1); 0; deg2rad(redCases(kc,2)); ...
                      redCases(kc,3); deg2rad(redCases(kc,4)); ...
                      deg2rad(redCases(kc,5))];
              uRed = deg2rad([redCases(kc,6); redCases(kc,7)]);
            xBoost = [xGlide; veh.mass];

            dGlide = coorbital.eom.glide3DOF(0,xGlide,uRed,veh,envRed);
            dBoost = coorbital.eom.boost3DOF(0,xBoost,uRed,veh,envRed);

%% A trivially zero glide derivative would make the comparison vacuous:
    assert(max(abs(dGlide)) > 0, ...
        'reduction case %d produced an identically zero glide derivative',kc);

%% Relative residual, floored by eps so an exactly zero component cannot
%% divide by zero. This is measured at 0 for every case -- see the header:
             resid = max(abs(dBoost(1:6) - dGlide)./max(abs(dGlide),eps));
    assert(resid < 1e-12, ...
        ['boost3DOF does not reduce to glide3DOF on case %d: relative ' ...
         'residual %.6e, expected < 1e-12'],kc,resid);
          redWorst = max(redWorst,resid);

%% A dead engine burns no propellant:
    assert(dBoost(7) == 0, ...
        'a dead engine must give dm/dt = 0, got %.6e on case %d',dBoost(7),kc);
end
    assert(redWorst < 1e-12, ...
        'worst reduction residual %.6e over all cases',redWorst);

%% =========================================================================
%% CHECK 2 -- thrust force increment. THE check that validates the new terms.
%%
%% The two environments differ in exactly one field, env.prop, so every
%% non-thrust model, the atmosphere, gravity, aerodynamics and the Earth
%% rate, is bit-identical between the two evaluations and cancels in the
%% difference. What is left is the thrust projection alone.
%% =========================================================================
       envOn.atmos = @coorbital.atmos.expAtmos;
        envOn.grav = @coorbital.grav.sphereGrav;
        envOn.aero = @coorbital.aero.constLD;
        envOn.prop = @coorbital.prop.constThrust;
      envOn.omegaE = c.omegaE;
            envOff = envOn;
       envOff.prop = @(tp,Pp,vp) deal(0,0);

%% Columns: t (s), h (m), lat (deg), V (m/s), gamma (deg), psi (deg), m (kg),
%% alpha (deg), sigma (deg). Both angles of attack and both bank angles are
%% nonzero, and neither bank angle has |cos| equal to |sin|, so the gamma and
%% psi channels carry different coefficients and cannot be swapped unseen:
          thrCases = [25, 20e3,  32, 1500, 48, 115, 21000,   8,  35; ...
                      60, 48e3, -17, 3800, 22, 250,  7000, -12, -70];

for kc = 1:size(thrCases,1)
               tFi = thrCases(kc,1);
               hFi = thrCases(kc,2);
               Vfi = thrCases(kc,4);
             gamFi = deg2rad(thrCases(kc,5));
               mFi = thrCases(kc,7);
             alpFi = deg2rad(thrCases(kc,8));
             sigFi = deg2rad(thrCases(kc,9));
               xFi = [c.rE + hFi; 0; deg2rad(thrCases(kc,3)); Vfi; gamFi; ...
                      deg2rad(thrCases(kc,6)); mFi];
               uFi = [alpFi; sigFi];

%% The state mass must differ from every mass the routine could reach for
%% instead, or a mutant substituting one of them would pass by accident:
    assert(mFi ~= veh.mass && mFi ~= bst.massDry && mFi ~= bst.massProp, ...
        'thrust case %d state mass collides with a parameter mass',kc);
    assert(mFi ~= veh.mass + bst.massDry + bst.massProp && ...
           mFi ~= veh.mass + bst.massDry, ...
        'thrust case %d state mass collides with a stack total',kc);

%% The reference thrust and flow come from the environment's OWN propulsion
%% model, called directly at the ambient pressure the state implies. Nothing
%% here is read back out of boost3DOF:
      [~,Pamb,~,~] = envOn.atmos(hFi);
    [Tref,mdotRef] = envOn.prop(tFi,Pamb,bst);
    assert(Tref > 0,'thrust case %d must run the engine, got T = %.3f N',kc,Tref);
    assert(abs(sin(alpFi)) > 0.1, ...
        'thrust case %d needs a genuinely nonzero angle of attack',kc);
    assert(abs(abs(cos(sigFi)) - abs(sin(sigFi))) > 0.1, ...
        'thrust case %d bank angle cannot separate the gamma and psi channels',kc);

%% Vinh Eqs. (2.28)-(2.33) with a thrust of magnitude T at incidence alpha,
%% banked through sigma, over the STATE mass:
             dVref = Tref.*cos(alpFi)./mFi;
           dGamRef = Tref.*sin(alpFi).*cos(sigFi)./(mFi.*Vfi);
           dPsiRef = Tref.*sin(alpFi).*sin(sigFi)./(mFi.*Vfi.*cos(gamFi));

              dxOn = coorbital.eom.boost3DOF(tFi,xFi,uFi,bst,envOn);
             dxOff = coorbital.eom.boost3DOF(tFi,xFi,uFi,bst,envOff);
              dInc = dxOn - dxOff;

%% Thrust is a force, so it cannot move the kinematic channels at all:
    assert(all(dInc(1:3) == 0), ...
        'thrust case %d moved the kinematic channels; increment %.6e', ...
        kc,max(abs(dInc(1:3))));

%% Tangential channel. This is the assertion that sees a missing cos(alpha):
              errV = abs(dInc(4) - dVref)./abs(dVref);
    assert(errV < 1e-12, ...
        ['thrust case %d dV/dt increment %.15e against T cos(alpha)/m = ' ...
         '%.15e, relative error %.6e'],kc,dInc(4),dVref,errV);

%% Flight-path channel. This is the assertion that sees a cos/sin(sigma) swap:
            errGam = abs(dInc(5) - dGamRef)./abs(dGamRef);
    assert(errGam < 1e-12, ...
        ['thrust case %d dgamma/dt increment %.15e against ' ...
         'T sin(alpha) cos(sigma)/(m V) = %.15e, relative error %.6e'], ...
        kc,dInc(5),dGamRef,errGam);

%% Heading channel:
            errPsi = abs(dInc(6) - dPsiRef)./abs(dPsiRef);
    assert(errPsi < 1e-12, ...
        ['thrust case %d dpsi/dt increment %.15e against ' ...
         'T sin(alpha) sin(sigma)/(m V cos(gamma)) = %.15e, relative ' ...
         'error %.6e'],kc,dInc(6),dPsiRef,errPsi);

%% Mass channel. mdot is positive by convention and the equations of motion
%% apply the minus, so a sign flip in dm/dt fails right here:
              errM = abs(dInc(7) + mdotRef)./abs(mdotRef);
    assert(errM < 1e-12, ...
        ['thrust case %d dm/dt increment %.15e against -mdot = %.15e, ' ...
         'relative error %.6e'],kc,dInc(7),-mdotRef,errM);
end

%% =========================================================================
%% CHECK 3 -- Tsiolkovsky. Vacuum, no gravity, thrust along the velocity
%% vector, so dV/dt = T/m is the whole speed equation and the integral is the
%% rocket equation in closed form. See the header for what this cannot see.
%% =========================================================================
%% Vacuum: zero density kills lift and drag whatever the coefficients are,
%% and zero ambient pressure gives constThrust its full vacuum thrust. The
%% sound speed is held at one so the Mach argument stays finite:
      envTsi.atmos = @(hh) deal(zeros(size(hh)),zeros(size(hh)), ...
                                c.T0.*ones(size(hh)),ones(size(hh)));
       envTsi.grav = @(rr,latr) deal(0,0);
       envTsi.aero = @coorbital.aero.constLD;
       envTsi.prop = @coorbital.prop.constThrust;
     envTsi.omegaE = 0;

%% Mass bookkeeping straight out of boosterDefaults; the propagated state
%% mass is always the TOTAL mass carried:
          mLiftoff = veh.mass + bst.massDry + bst.massProp;
          mBurnout = veh.mass + bst.massDry;
    assert(mLiftoff == 32400,'sanity: liftoff mass should be 32400 kg');
    assert(mBurnout == 2400, 'sanity: burnout mass should be 2400 kg');

%% Hand arithmetic (Sutton & Biblarz Sec. 4.1): exhaust speed = Isp*g0 =
%% 260*9.80665 = 2549.729 m/s; mass ratio = 32400/2400 = 13.5;
%% deltaV = 2549.729*log(13.5) = 6636.153368978423 m/s:
           dVideal = bst.Isp.*c.g0.*log(mLiftoff./mBurnout);
    assert(abs(dVideal - 6636.153368978423) < 1e-9, ...
        'sanity: ideal delta-V should be 6636.153368978423 m/s, got %.9f',dVideal);

%% Thrust along the velocity vector, so alpha is zero for the whole burn. The
%% bank angle is irrelevant with no normal force and is held at zero too:
          schedTsi = struct('tGrid',[0 200],'alpha',[0 0],'sigma',[0 0]);
         phTsi.eom = @coorbital.eom.boost3DOF;
       phTsi.guide = @(t,x) coorbital.guide.prescribed(t,x,schedTsi);
   phTsi.terminate = @(t,x) coorbital.prop.eventBurnout(t,x,mBurnout);
       phTsi.tspan = [0 200];

%% A pitched-over liftoff state: gamma = 90 deg is the vertical-flight
%% singularity, so the ascent starts already turned, exactly as boost3DOF's
%% own self-demo does. With gravity off, the flight-path angle drifts by only
%% about 1.5 deg over the burn, nowhere near the guard:
             x0Tsi = [c.rE + 5e3; 0; 0; 300; deg2rad(60); deg2rad(90); mLiftoff];
           trajTsi = coorbital.prop.phaseRun(phTsi,x0Tsi,bst,envTsi);

%% It stopped on burnout, not on the horizon:
    assert(trajTsi.t(end) < 200, ...
        'the vacuum burn hit the 200 s horizon instead of the burnout event');
    assert(abs(trajTsi.x(end,7) - mBurnout) < 1e-6, ...
        'the vacuum burn ended at %.9f kg, expected %.9f kg', ...
        trajTsi.x(end,7),mBurnout);

             dVsim = trajTsi.x(end,4) - trajTsi.x(1,4);
            errTsi = abs(dVsim - dVideal)./dVideal;
    assert(errTsi < 1e-8, ...
        ['propagated vacuum delta-V %.9f m/s against Tsiolkovsky ' ...
         '%.9f m/s, relative error %.6e, expected < 1e-8'], ...
        dVsim,dVideal,errTsi);

%% =========================================================================
%% CHECK 4 -- mass linearity. Constant mass flow makes the mass history an
%% exact line, so any error in dm/dt shows up as curvature or slope.
%% =========================================================================
%% mdot from boosterDefaults and missileConst, never from the trajectory:
           mdotTsi = bst.thrustVac./(bst.Isp.*c.g0);
    assert(abs(mdotTsi - 372.5886162803969) < 1e-9, ...
        'sanity: mdot should be 372.5886162803969 kg/s, got %.10f',mdotTsi);

             mLine = mLiftoff - mdotTsi.*trajTsi.t;
           errMass = max(abs(trajTsi.x(:,7) - mLine)./mLine);
    assert(errMass < 1e-9, ...
        ['propagated mass departs the analytic line m = %.1f - %.9f t by ' ...
         '%.6e relative, expected < 1e-9'],mLiftoff,mdotTsi,errMass);

%% The line must actually descend across the full propellant budget, or the
%% comparison above would be satisfied by a trajectory that barely moved:
    assert(abs((mLine(1) - mLine(end)) - bst.massProp) < 1e-3, ...
        ['the analytic line spans %.6f kg over the recorded burn against a ' ...
         'propellant budget of %.1f kg'],mLine(1) - mLine(end),bst.massProp);
end
