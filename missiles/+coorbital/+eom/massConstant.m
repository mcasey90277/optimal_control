function eomFn = massConstant(baseEom)
%% Purpose:
%
%  Adapt a six-state atmospheric-flight EOM to the seven-state vector used by
%  a chain that includes a powered phase, by appending dm/dt = 0. Unpowered
%  phases carry mass so that every phase in a chain shares one state vector,
%  which removes the need for a state mapping at the junctions.
%
%% Inputs:
%
%  baseEom          Function handle             Six-state EOM with the
%                                               signature xdot = f(t,x,u,veh,env)
%
%% Outputs:
%
%  eomFn            Function handle             Seven-state EOM with the same
%                                               signature; component 7 is mass
%                                               (kg) and its derivative is zero
%
%% Revision History:
%  Michael Casey                                                08/07/2026
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
                fD = coorbital.eom.massConstant(@coorbital.eom.glide3DOF);
                xD = [c.rE + 40e3; 0; 0; 4000; deg2rad(-2); deg2rad(90); veh.mass];
             xdotD = fD(0,xD,[0;0],veh,env);
    fprintf('dm/dt = %.1f kg/s over a %d-state vector\n',xdotD(7),numel(xdotD));
             eomFn = [];
    return;
end

             eomFn = @(t,x,u,veh,env) [baseEom(t,x(1:6),u,veh,env); 0];
end
