function test_equilibriumGlide()
%% Purpose:
%
%  Validate the propagated glide against the closed-form equilibrium glide
%  condition -- a solution derived independently of the propagator, not
%  another run of the propagator. Setting dgamma/dt = 0 in the flight path
%  angle equation with the bank angle at zero gives
%
%      0.5 rho V^2 S CL = (m g - m V^2 / r) cos(gamma)
%
%  and therefore a closed-form speed at every radius,
%
%      V_eq^2(r) = m g(r) cos(gamma) / (0.5 rho(r) S CL + m cos(gamma) / r).
%
%  The cos(gamma) factor is the exact consequence of dgamma/dt = 0, and it is
%  carried so that the reference is the equilibrium condition itself rather
%  than an approximation to it. On this arc it is not decorative: the glide
%  starts at gamma = -0.16 deg but steepens to -3.13 deg by the slow end,
%  where cos(gamma) is 1.5e-3 below unity. Dropping the factor would bias
%  V_eq^2 by that amount, about 0.07 percent in V_eq -- small against the
%  1.5 percent quasi-steady departure measured below, but a bias rather than
%  a scatter, and there is no reason to accept it.
%
%  Two points of discipline, both learned the hard way on this task:
%
%  (1) The reference is built from the VEHICLE's declared parameters --
%      veh.mass, veh.Sref, veh.CL -- and not from whatever the aerodynamic
%      model happens to hand back. The analytic solution is a function of the
%      physical vehicle. If the aero model misreports the lift coefficient
%      then the propagator is flying a different vehicle from the one the
%      closed form describes, which is a real defect and must show up here.
%      Reading CL back out of coorbital.aero.constLD would make this test
%      blind to precisely that error, because the mistake would then enter
%      both sides of the comparison equally and cancel. Measured: scaling CL
%      by 1.2 inside constLD LOWERS the mismatch when the reference reads CL
%      from the model (14.2 down to 9.9 percent on the arc that formulation
%      flew), and raises it from 1.51 to 11.52 percent when the reference
%      reads it from the vehicle.
%
%  (2) An open-loop constant-alpha arc is not obliged to sit on the
%      equilibrium curve. It phugoids about it, and at these altitudes the
%      oscillation damps over many hundreds of seconds -- far too slowly to
%      be waited out by skipping a fixed fraction of the samples. So rather
%      than launch off the curve and hope, the initial state is placed ON it.
%      The entry speed is V_eq at the entry radius, and the entry flight path
%      angle is the one that makes the trajectory TANGENT to the curve:
%
%          d/dt (V^2 - V_eq^2(r)) = 0,   rdot = V sin(gamma),
%          Vdot = -D/m - g sin(gamma)
%
%          ==>  sin(gamma) = -(2 D / m) / (2 g + dV_eq^2/dr).
%
%      With no transient to wait out, the WHOLE arc is checked rather than a
%      tail of it, and the tolerance can be set by the physics instead of by
%      how long the transient happens to take.
%
%  Scope of the claim. The reference carries no drag term, so V_eq(r) is
%  independent of CD: what is tested here is the lift / weight / centrifugal
%  balance in the gamma equation, and the trajectory's continued adherence to
%  that balance across a 25 km altitude band over which the density changes by
%  a factor of 30 and the speed by a factor of nearly 5. Drag is validated
%  separately, in test_allenEggers.
%
%  Tolerance. The measured worst-case departure over the whole arc is about
%  1.5 percent, reached at the low, slow end where the quasi-steady assumption
%  starts to give way (drag is no longer small against the terms retained).
%  The assertion is set at 3 percent: roughly a factor of two of headroom over
%  the measured quasi-steady departure, and well below the 11.52 percent a
%  20 percent lift error produces.
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
%   [1] Saenger, E., Bredt, I., "A Rocket Drive for Long Range Bombers,"
%       Deutsche Luftfahrtforschung UM-3538, 1944. (equilibrium glide)
%   [2] Eggers, A.J., Allen, H.J., Neice, S.E., "A Comparative Analysis of the
%       Performance of Long-Range Hypervelocity Vehicles," NACA TR-1382, 1958.
%   [3] Vinh, N.X., Busemann, A., Culp, R.D., "Hypersonic and Planetary Entry
%       Flight Mechanics," Univ. Michigan Press, 1980, Sec. 5.3.
%
%% Revision History:
%  Michael Casey                                                08/06/2026
%  Copyright 2026 Coorbital, Inc.
%% ------------------------ Begin Code Sequence ---------------------------

                 c = coorbital.util.missileConst();
               veh = coorbital.util.vehicleDefaults();
         env.atmos = @coorbital.atmos.expAtmos;
          env.grav = @coorbital.grav.sphereGrav;
          env.aero = @coorbital.aero.constLD;
        env.omegaE = 0;

%% Entry radius and the equilibrium speed there. cos(gamma) is taken as unity
%% for the initial condition; the angle solved for below is under 0.2 deg, so
%% the neglected feedback on V_eq is at the 5e-6 level:
                h0 = 45e3;
                r0 = c.rE + h0;
              vEq0 = sqrt(eqGlideSpeedSq(r0,1,veh,env));

%% Tangency condition: the entry flight path angle for which the arc starts on
%% the equilibrium manifold and moves along it, exciting no phugoid. dV_eq^2/dr
%% comes from a central difference on the closed form itself, 1 m step:
                dr = 1;
             dv2dr = (eqGlideSpeedSq(r0+dr,1,veh,env) ...
                    - eqGlideSpeedSq(r0-dr,1,veh,env))./(2.*dr);
      [rho0,~,~,~] = coorbital.atmos.expAtmos(h0);
           [gr0,~] = coorbital.grav.sphereGrav(r0,0);
            [~,CD] = coorbital.aero.constLD(0,0,veh);
              aDr0 = 0.5.*rho0.*vEq0.^2.*veh.Sref.*CD./veh.mass;
              gam0 = asin(-2.*aDr0./(2.*gr0 + dv2dr));

%% The tangency angle must come out shallow and descending, or the premise of
%% the whole comparison is already broken:
    assert(gam0 < 0 && gam0 > deg2rad(-2), ...
        'tangency gave gamma0 = %.5f deg, expected a shallow descent', ...
        rad2deg(gam0));

%% Propagate unbanked at constant alpha, due east along the equator, from
%% 45 km down to 20 km:
             sched = struct('tGrid',[0 20000],'alpha',[0 0],'sigma',[0 0]);
            ph.eom = @coorbital.eom.glide3DOF;
          ph.guide = @(t,x) coorbital.guide.prescribed(t,x,sched);
      ph.terminate = @(t,x) coorbital.prop.eventAltitude(t,x,20e3);
          ph.tspan = [0 20000];
                x0 = [r0; 0; 0; vEq0; gam0; deg2rad(90)];
              traj = coorbital.prop.phaseRun(ph,x0,veh,env);

    assert(traj.t(end) < 20000, ...
        'the glide hit the time horizon instead of the 20 km event');
    assert(abs(traj.x(end,1) - c.rE - 20e3) < 1, ...
        'the glide terminated at %.1f m, expected 20000 m',traj.x(end,1) - c.rE);

%% An unbanked due-east arc launched from the equator must stay on it, which is
%% what licenses evaluating gravity at zero latitude in the reference below:
              rArc = traj.x(:,1);
              vArc = traj.x(:,4);
              gArc = traj.x(:,5);
    assert(max(abs(traj.x(:,3))) < 1e-12, ...
        'the arc left the equator by %.3e rad; the reference assumes it does not', ...
        max(abs(traj.x(:,3))));

%% The arc must actually sweep a wide range, or the comparison is being made
%% at effectively one operating point. These guard against anyone rescuing a
%% future failure by shortening the arc rather than fixing the physics:
    [rhoArc,~,~,~] = coorbital.atmos.expAtmos(rArc - c.rE);
            rhoSpan = max(rhoArc)./min(rhoArc);
              vSpan = max(vArc)./min(vArc);
    assert(rhoSpan > 20, ...
        'density spanned only a factor of %.1f; the arc is too short to validate', ...
        rhoSpan);
    assert(vSpan > 3, ...
        'speed spanned only a factor of %.1f; the arc is too short to validate', ...
        vSpan);

%% The shallow-glide premise, checked rather than assumed. The measured worst
%% case is -3.13 deg at the slow end, where cos(gamma) is 1.5e-3 below unity.
%% These are GUARDS, not the physics assertion, and they are deliberately loose:
%% every field of vehicleDefaults is marked PLACEHOLDER, and a routine change to
%% one of them moves this angle. Doubling Sref, for instance, steepens the glide
%% enough to trip a 5 deg bound -- and a run that failed HERE rather than on the
%% relative-speed assertion below would report a broken premise when what had
%% actually happened is that the physics check was never reached. 8 deg is still
%% unambiguously a glide and leaves that room:
    assert(max(abs(gArc)) < deg2rad(8), ...
        'flight path angle reached %.4f deg; the glide is no longer shallow', ...
        rad2deg(max(abs(gArc))));
    assert(max(abs(1 - cos(gArc))) < 1e-2, ...
        'cos(gamma) departed unity by %.3e; the small-angle premise has gone', ...
        max(abs(1 - cos(gArc))));

%% The analytic reference, evaluated at each propagated radius from the
%% VEHICLE's declared parameters. Nothing on this line came out of the
%% propagator except the radius and the flight path angle at which to
%% evaluate it:
               vEq = sqrt(eqGlideSpeedSq(rArc,cos(gArc),veh,env));
             relEr = abs(vArc - vEq)./vEq;

    assert(max(relEr) < 0.03, ...
        ['equilibrium glide mismatch: max relative speed error %.4f over the ' ...
         'whole arc, expected < 0.03'],max(relEr));
end

function v2 = eqGlideSpeedSq(r,cgam,veh,env)
%% Purpose:
%
%  Closed-form equilibrium glide speed squared at a given radius. This is the
%  analytic reference the test measures against, so it is written in terms of
%  the vehicle's declared physical parameters and the environment models, and
%  never in terms of anything the propagator returned. In particular the lift
%  coefficient is read from veh.CL, NOT from env.aero -- see the discipline
%  note in the parent function.
%
%  Latitude is passed to the gravity model as zero, which the caller asserts
%  is true for the equatorial arc it flies.
%
%% Inputs:
%
%  r                [N x 1]                     Geocentric radius (m)
%
%  cgam             [N x 1]                     Cosine of the flight path
%                                               angle (-), or scalar 1
%
%  veh              Struct                      Vehicle parameters; uses
%                                               mass (kg), Sref (m^2), CL (-)
%
%  env              Struct                      Environment model handles;
%                                               uses atmos and grav
%
%% Outputs:
%
%  v2               [N x 1]                     Equilibrium glide speed
%                                               squared (m^2/s^2)
%
%% References:
%   [1] Vinh, N.X., Busemann, A., Culp, R.D., "Hypersonic and Planetary Entry
%       Flight Mechanics," Univ. Michigan Press, 1980, Sec. 5.3.
%
%% Revision History:
%  Michael Casey                                                08/06/2026
%  Copyright 2026 Coorbital, Inc.
%% ------------------------ Begin Code Sequence ---------------------------

                 c = coorbital.util.missileConst();
       [rho,~,~,~] = env.atmos(r - c.rE);
            [gr,~] = env.grav(r,zeros(size(r)));

%% Lift balances weight less centrifugal relief, both projected by cos(gamma):
                v2 = (veh.mass.*gr.*cgam)./ ...
                     (0.5.*rho.*veh.Sref.*veh.CL + veh.mass.*cgam./r);
end
