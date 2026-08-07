function xdot = boost3DOF(t,x,u,veh,env)
%% Purpose:
%
%  Three-degree-of-freedom point-mass equations of motion for powered flight
%  over a rotating spherical Earth: the powered sibling of
%  coorbital.eom.glide3DOF, with a seventh state for the total mass currently
%  carried. Thrust of magnitude T acts at angle of attack alpha from the
%  velocity vector, in the plane set by bank angle sigma -- the same
%  convention lift already uses, so thrust and lift add directly. The mass in
%  every acceleration denominator is the STATE x(7), so propellant burn feels
%  the vehicle getting lighter; dm/dt = -mdot closes the rocket equation. The
%  rotation terms are gated on env.omegaE, so setting it to zero recovers the
%  non-rotating case exactly; the latitudinal gravity component gLat is
%  carried through so that an oblate gravity model drops in without editing
%  this file. With env.prop returning zeros the first six components reduce
%  exactly to glide3DOF.
%
%% Inputs:
%
%  t                scalar                      Time since phase start (s).
%                                               Passed to env.prop, whose
%                                               thrust may be time-varying.
%
%  x                [7 x 1]                     State:
%                                               r     (m) geocentric radius
%                                               lon   (rad) longitude
%                                               lat   (rad) geocentric latitude
%                                               V     (m/s) planet-relative speed
%                                               gamma (rad) flight path angle
%                                               psi   (rad) heading, clockwise
%                                                     from north
%                                               m     (kg) TOTAL mass currently
%                                                     carried
%
%  u                [2 x 1]                     Control:
%                                               alpha (rad) angle of attack
%                                               sigma (rad) bank angle
%
%  veh              Struct                      Boosted-stack parameters from
%                                               coorbital.util.boosterDefaults
%
%  env              Struct                      Environment model handles:
%                                               atmos, grav, aero, prop, omegaE
%
%% Outputs:
%
%  xdot             [7 x 1]                     State derivative
%
%% References:
%   [1] Vinh, N.X., Busemann, A., Culp, R.D., "Hypersonic and Planetary Entry
%       Flight Mechanics," Univ. Michigan Press, 1980, Eqs. (2.28)-(2.33).
%   [2] Sutton, G.P., Biblarz, O., "Rocket Propulsion Elements," 9th ed.,
%       Wiley, 2017, Sec. 4.1 (the ideal rocket equation, dm/dt = -mdot).
%
%% Revision History:
%  Michael Casey                                                08/07/2026
%  Copyright 2026 Coorbital, Inc.
%% ------------------------ Begin Code Sequence ---------------------------

%% Self-demo:
if nargin == 0
                 c = coorbital.util.missileConst();
               bst = coorbital.util.boosterDefaults();
               pay = coorbital.util.vehicleDefaults();
         env.atmos = @coorbital.atmos.expAtmos;
          env.grav = @coorbital.grav.sphereGrav;
          env.aero = @coorbital.aero.constLD;
          env.prop = @coorbital.prop.constThrust;
        env.omegaE = 0;
                m0 = pay.mass + bst.massDry + bst.massProp;
                x0 = [c.rE + 5e3; 0; 0; 300; deg2rad(60); deg2rad(90); m0];
              xdot = coorbital.eom.boost3DOF(10,x0,[0;0],bst,env);
    fprintf('xdot = [%.4g %.4g %.4g %.4g %.4g %.4g %.4g]''\n',xdot);
    return;
end

%% Unpack the state and control:
                 r = x(1);
               lat = x(3);
                 V = x(4);
             gamma = x(5);
               psi = x(6);
                 m = x(7);
             alpha = u(1);
             sigma = u(2);

                 c = coorbital.util.missileConst();
                 h = r - c.rE;

%% Guard the coordinate singularities so they fail loudly, not as silent NaN:
if abs(cos(lat)) < 1e-8
    error('coorbital:boost3DOF:polarSingularity', ...
        'Longitude rate is singular at the pole (lat = %.6f rad).',lat);
end
if abs(cos(gamma)) < 1e-8
    error('coorbital:boost3DOF:verticalFlight', ...
        'Heading rate is singular in vertical flight (gamma = %.6f rad).',gamma);
end
if V < 1
    error('coorbital:boost3DOF:zeroSpeed', ...
        'Equations are singular as V approaches zero (V = %.3f m/s).',V);
end
if m <= 0
    error('coorbital:boost3DOF:massDepleted', ...
        'Mass state is non-positive (m = %.3f kg); accelerations would invert.',m);
end

%% Environment: density, pressure, gravity, aerodynamic coefficients:
   [rho,P,~,aSnd] = env.atmos(h);
        [gr,gLat] = env.grav(r,lat);
          [CL,CD] = env.aero(alpha,V./aSnd,veh);

%% Propulsion: thrust and (positive) mass flow at ambient pressure P:
         [T,mdot] = env.prop(t,P,veh);

%% Aerodynamic and thrust accelerations, all over the STATE mass m = x(7):
              qbar = 0.5.*rho.*V.^2;
             aLift = qbar.*veh.Sref.*CL./m;
             aDrag = qbar.*veh.Sref.*CD./m;
             aThrV = T.*cos(alpha)./m;
             aThrN = T.*sin(alpha)./m;

%% Kinematics:
              rdot = V.*sin(gamma);
            londot = V.*cos(gamma).*sin(psi)./(r.*cos(lat));
            latdot = V.*cos(gamma).*cos(psi)./r;

%% Dynamics over a non-rotating sphere:
              Vdot = aThrV - aDrag - gr.*sin(gamma) ...
                     + gLat.*cos(gamma).*cos(psi);
          gammadot = (aThrN + aLift).*cos(sigma)./V - (gr./V - V./r).*cos(gamma) ...
                     - gLat.*sin(gamma).*cos(psi)./V;
            psidot = (aThrN + aLift).*sin(sigma)./(V.*cos(gamma)) ...
                     + V.*cos(gamma).*sin(psi).*tan(lat)./r ...
                     - gLat.*sin(psi)./(V.*cos(gamma));

%% Rotating-Earth Coriolis and centrifugal terms, gated on omegaE:
                om = env.omegaE;
if om ~= 0
              Vdot = Vdot + om.^2.*r.*cos(lat).* ...
                     (sin(gamma).*cos(lat) - cos(gamma).*sin(lat).*cos(psi));
          gammadot = gammadot + 2.*om.*cos(lat).*sin(psi) ...
                     + om.^2.*r.*cos(lat).* ...
                       (cos(gamma).*cos(lat) + sin(gamma).*cos(psi).*sin(lat))./V;
            psidot = psidot + 2.*om.*(sin(lat) - cos(lat).*cos(psi).*tan(gamma)) ...
                     + om.^2.*r.*sin(lat).*cos(lat).*sin(psi)./(V.*cos(gamma));
end

%% Assemble; the equations of motion apply the minus sign to the positive mdot:
              xdot = [rdot; londot; latdot; Vdot; gammadot; psidot; -mdot];
end
