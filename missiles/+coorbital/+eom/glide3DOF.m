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
%  t                scalar                      Phase clock (s): the phase's
%                                               OWN tspan value, not rebased
%                                               to zero. A phase given
%                                               tspan = [10 50] is evaluated
%                                               at 10 through 50. See TWO
%                                               CLOCKS in
%                                               coorbital.prop.phaseRun.
%                                               Unused; present for ode45
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
%  xdot             [6 x 1]                     State derivative. Raises
%                                               coorbital:glide3DOF:nonFiniteState,
%                                               :nonFiniteControl, :badRadius,
%                                               :badVehicleMass,
%                                               :badAtmosphere, :badGravity,
%                                               :badAero or
%                                               :nonFiniteDerivative rather
%                                               than return a NaN, in addition
%                                               to the three coordinate
%                                               singularity identifiers
%
%% References:
%   [1] Vinh, N.X., Busemann, A., Culp, R.D., "Hypersonic and Planetary Entry
%       Flight Mechanics," Univ. Michigan Press, 1980, Eqs. (2.28)-(2.33).
%
%% Revision History:
%  Michael Casey                                                08/06/2026
%  Michael Casey  Finiteness guards; NaN passes every threshold  08/07/2026
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

%% Nothing below this line can be believed unless the inputs are numbers, and
%% A COMPARISON IS NOT A VALIDATION. Every threshold guard in this file is
%% blind to NaN -- NaN < 1 is false, abs(cos(NaN)) < 1e-8 is false, and so is
%% every other inequality -- so a NaN arriving in the state or the control
%% walks straight past all of them, poisons xdot, and lets the integrator
%% march a trajectory of NaNs to tspan(end) without a word. The finiteness
%% tests therefore come FIRST, and each carries its own identifier so the
%% failure names where the bad number came from:
if numel(x) < 6 || ~isreal(x(1:6)) || ~all(isfinite(x(1:6)))
    error('coorbital:glide3DOF:nonFiniteState', ...
        ['The state must be six finite real components; got %s. NaN passes ' ...
         'every threshold guard in this file, so it is refused here.'], ...
        mat2str(x(:).'));
end
if numel(u) < 2 || ~isreal(u(1:2)) || ~all(isfinite(u(1:2)))
    error('coorbital:glide3DOF:nonFiniteControl', ...
        'The control must be two finite real components; got %s.', ...
        mat2str(u(:).'));
end
if r <= 0
    error('coorbital:glide3DOF:badRadius', ...
        ['The geocentric radius must be positive, the kinematics dividing ' ...
         'by it; got r = %.6g m.'],r);
end
if ~isscalar(veh.mass) || ~isreal(veh.mass) || ~isfinite(veh.mass) || ...
        veh.mass <= 0
    error('coorbital:glide3DOF:badVehicleMass', ...
        ['veh.mass must be a finite positive real scalar, both aerodynamic ' ...
         'accelerations dividing by it; got %s.'],mat2str(veh.mass));
end

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

%% The models are handles the caller installed, so their outputs are inputs
%% to this file and get the same treatment. aSnd is checked before it is used
%% as the Mach divisor, which is why the aero call sits below rather than
%% beside its two siblings:
if ~isreal(rho) || ~isfinite(rho) || rho < 0 || ...
   ~isreal(aSnd) || ~isfinite(aSnd) || aSnd <= 0
    error('coorbital:glide3DOF:badAtmosphere', ...
        ['env.atmos returned rho = %s and aSnd = %s at h = %.6g m. Density ' ...
         'must be finite and non-negative; the speed of sound must be ' ...
         'finite and positive, the Mach number dividing by it.'], ...
        mat2str(rho),mat2str(aSnd),h);
end
if ~isreal(gr) || ~isfinite(gr) || ~isreal(gLat) || ~isfinite(gLat)
    error('coorbital:glide3DOF:badGravity', ...
        'env.grav returned gr = %s and gLat = %s at r = %.6g m.', ...
        mat2str(gr),mat2str(gLat),r);
end

          [CL,CD] = env.aero(alpha,V./aSnd,veh);
if ~isreal(CL) || ~isfinite(CL) || ~isreal(CD) || ~isfinite(CD)
    error('coorbital:glide3DOF:badAero', ...
        'env.aero returned CL = %s and CD = %s at Mach %.6g.', ...
        mat2str(CL),mat2str(CD),V./aSnd);
end

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

%% Assemble, and check the result. The guards above cover every input and
%% every model output; this covers the arithmetic between them, and it is the
%% only place a cancellation such as Inf - Inf can be caught:
              xdot = [rdot; londot; latdot; Vdot; gammadot; psidot];
if ~isreal(xdot) || ~all(isfinite(xdot))
    error('coorbital:glide3DOF:nonFiniteDerivative', ...
        'The equations of motion produced a non-finite derivative, %s.', ...
        mat2str(xdot.'));
end
end
