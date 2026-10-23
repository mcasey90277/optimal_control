function p = drone_params(model)
%% Purpose:
%
%   THE physical constants of the three drone models, one home. Every
%   solver and every test reads them from here; nothing else holds a copy.
%
%   kinematic   p = position, u = velocity              (the paper's drone)
%   pointmass   [p; v], u = thrust acceleration         (the planner's drone)
%   rigid6dof   [p; v; Theta; omega], u = rotor thrusts (the vehicle)
%
%   rigid6dof constants are Crazyflie-class: Luis & Le Ny, "Design of a
%   Trajectory Tracking Controller for a Nanoquadcopter", 2016, Table 1, and
%   Forster, "System Identification of the Crazyflie 2.0 Nano Quadrocopter",
%   ETH 2015. thetaMax is the pitch/roll box every direct solve imposes: the
%   ZYX Euler kinematics are singular at 90 degrees.
%
%% Inputs:
%
%  model                    char                    'kinematic' |
%                                                   'pointmass' |
%                                                   'rigid6dof'
%
%% Outputs:
%
%  p                        struct                  .model .nx .nu always;
%                                                   pointmass: .g .aMax;
%                                                   rigid6dof: .m .I [3x1]
%                                                   .l .kappa .g .fMax
%                                                   .thetaMax
%
%% Revision History:
%  M. Casey                                                   (c) 10/23/2026
%  Copyright Coorbital Inc.
%% ------------------------ Begin Code Sequence ---------------------------

if nargin == 0
   %Demo: the three parameter sets:
     disp(drone_params('kinematic'));
     disp(drone_params('pointmass'));
     disp(drone_params('rigid6dof'));
     return
end

m      = 0.027;
fMax   = 0.15;
g      = 9.81;
switch lower(model)
    case 'kinematic'
        p = struct('model', 'kinematic', 'nx', 3, 'nu', 3);
    case 'pointmass'
        p = struct('model', 'pointmass', 'nx', 6, 'nu', 3, 'g', g, ...
                   'aMax', 4*fMax/m - g);
    case 'rigid6dof'
        p = struct('model', 'rigid6dof', 'nx', 12, 'nu', 4, 'm', m, ...
                   'I', [1.395e-5; 1.436e-5; 2.173e-5], 'l', 0.046, ...
                   'kappa', 0.005964, 'g', g, 'fMax', fMax, 'thetaMax', pi/3);
    otherwise
        error('drone_params:model', 'unknown model ''%s''', model);
end
end
