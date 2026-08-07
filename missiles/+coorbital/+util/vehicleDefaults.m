function veh = vehicleDefaults()
%% Purpose:
%
%  Default hypersonic glide vehicle parameters. PLACEHOLDER VALUES drawn from
%  the open literature for a generic lifting entry body -- they are of the
%  right order but are not any specific vehicle. Replace before any result
%  leaves this machine.
%
%% Inputs:
%
%  none
%
%% Outputs:
%
%  veh              Struct                      Vehicle parameters:
%                                               mass       (kg)
%                                               Sref       (m^2) reference area
%                                               CL         (-) lift coefficient
%                                               LD         (-) lift-to-drag ratio
%                                               noseRadius (m) for heating
%
%% Revision History:
%  Michael Casey                                                08/06/2026
%  Copyright 2026 Coorbital, Inc.
%% ------------------------ Begin Code Sequence ---------------------------

          veh.mass = 900;              %kg, PLACEHOLDER
          veh.Sref = 0.75;             %m^2, PLACEHOLDER
            veh.CL = 0.35;             %PLACEHOLDER
            veh.LD = 2.5;              %PLACEHOLDER, typical waverider range 2-4
    veh.noseRadius = 0.05;             %m, PLACEHOLDER
end
