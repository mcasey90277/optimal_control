function [CL,CD] = constLD(alpha,mach,veh)
%% Purpose:
%
%  Constant lift coefficient and constant lift-to-drag ratio. The crudest
%  usable aerodynamic model, and the one for which equilibrium glide has a
%  closed-form solution. Signature matches table-driven models so the
%  environment struct can swap them.
%
%% Inputs:
%
%  alpha            [N x 1]                     Angle of attack (rad). Accepted
%                                               and ignored by this model.
%
%  mach             [N x 1]                     Mach number. Accepted and
%                                               ignored by this model.
%
%  veh              Struct                      Vehicle parameters; uses the
%                                               CL and LD fields
%
%% Outputs:
%
%  CL               scalar                      Lift coefficient (-)
%
%  CD               scalar                      Drag coefficient (-)
%
%% Revision History:
%  Michael Casey                                                08/06/2026
%  Copyright 2026 Coorbital, Inc.
%% ------------------------ Begin Code Sequence ---------------------------

%% Self-demo:
if nargin == 0
               veh = coorbital.util.vehicleDefaults();
       [CLd,CDd] = coorbital.aero.constLD(0,20,veh);
    fprintf('CL = %.4f, CD = %.4f, L/D = %.2f\n',CLd,CDd,CLd/CDd);
    [CL,CD] = deal([]);
    return;
end

%% Constant by construction; alpha and mach are deliberately unused:
                CL = veh.CL;
                CD = veh.CL/veh.LD;
end
