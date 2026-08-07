function xdot = glide3DOF(t,x,u,veh,env)
%% Purpose:
%
%  Three-degree-of-freedom point-mass equations of motion for atmospheric
%  glide over a rotating spherical Earth. The rotation terms are gated on
%  env.omegaE, so setting it to zero recovers the non-rotating case exactly;
%  the latitudinal gravity component gLat is carried through so that an
%  oblate gravity model drops in without editing this file.
%
%% Inputs:
%
%  t                scalar                      Time since phase start (s).
%                                               Unused; present for ode45.
%
%  x                [6 x 1]                     State:
%                                               r     (m) geocentric radius
%                                               lon   (rad) longitude
%                                               lat   (rad) geocentric latitude
%                                               V     (m/s) planet-relative speed
%                                               gamma (rad) flight path angle
%                                               psi   (rad) heading, clockwise
%                                                     from north
%
%  u                [2 x 1]                     Control:
%                                               alpha (rad) angle of attack
%                                               sigma (rad) bank angle
%
%  veh              Struct                      Vehicle parameters from
%                                               coorbital.util.vehicleDefaults
%
%  env              Struct                      Environment model handles:
%                                               atmos, grav, aero, omegaE
%
%% Outputs:
%
%  xdot             [6 x 1]                     State derivative
%
%% References:
%   [1] Vinh, N.X., Busemann, A., Culp, R.D., "Hypersonic and Planetary Entry
%       Flight Mechanics," Univ. Michigan Press, 1980, Eqs. (2.28)-(2.33).
%
%% Revision History:
%  Michael Casey                                                08/06/2026
%  Copyright 2026 Coorbital, Inc.
%% ------------------------ Begin Code Sequence ---------------------------

%% Self-demo:
if nargin == 0
                 c = coorbital.util.missileConst();
               veh = coorbital.util.vehicleDefaults();
         env.atmos = @coorbital.atmos.expAtmos;
          env.grav = @coorbital.grav.sphereGrav;
          env.aero = @coorbital.aero.constLD;
        env.omegaE = 0;
                x0 = [c.rE + 60e3; 0; 0; 6000; deg2rad(-1); deg2rad(90)];
              xdot = coorbital.eom.glide3DOF(0,x0,[0;0],veh,env);
    fprintf('xdot = [%.4g %.4g %.4g %.4g %.4g %.4g]''\n',xdot);
    return;
end

%% Unpack the state and control:
                 r = x(1);
               lat = x(3);
                 V = x(4);
             gamma = x(5);
               psi = x(6);
             alpha = u(1);
             sigma = u(2);

                 c = coorbital.util.missileConst();
                 h = r - c.rE;

%% Guard the coordinate singularities so they fail loudly, not as silent NaN:
if abs(cos(lat)) < 1e-8
    error('coorbital:glide3DOF:polarSingularity', ...
        'Longitude rate is singular at the pole (lat = %.6f rad).',lat);
end
if abs(cos(gamma)) < 1e-8
    error('coorbital:glide3DOF:verticalFlight', ...
        'Heading rate is singular in vertical flight (gamma = %.6f rad).',gamma);
end
if V < 1
    error('coorbital:glide3DOF:zeroSpeed', ...
        'Equations are singular as V approaches zero (V = %.3f m/s).',V);
end

%% Environment: density, gravity, aerodynamic coefficients:
   [rho,~,~,aSnd] = env.atmos(h);
        [gr,gLat] = env.grav(r,lat);
          [CL,CD] = env.aero(alpha,V./aSnd,veh);

%% Aerodynamic accelerations:
              qbar = 0.5.*rho.*V.^2;
             aLift = qbar.*veh.Sref.*CL./veh.mass;
             aDrag = qbar.*veh.Sref.*CD./veh.mass;

%% Kinematics:
              rdot = V.*sin(gamma);
            londot = V.*cos(gamma).*sin(psi)./(r.*cos(lat));
            latdot = V.*cos(gamma).*cos(psi)./r;

%% Dynamics over a non-rotating sphere:
              Vdot = -aDrag - gr.*sin(gamma) ...
                     + gLat.*cos(gamma).*cos(psi);
          gammadot = aLift.*cos(sigma)./V - (gr./V - V./r).*cos(gamma) ...
                     - gLat.*sin(gamma).*cos(psi)./V;
            psidot = aLift.*sin(sigma)./(V.*cos(gamma)) ...
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

%% Assemble:
              xdot = [rdot; londot; latdot; Vdot; gammadot; psidot];
end
