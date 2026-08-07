function veh = vehicle_hgv()
%% Purpose:
%
%  Hypersonic glide vehicle parameters for the run_glide entry script.
%  Starts from the library defaults so there is one place to change a value
%  that should apply everywhere, and overrides only what is vehicle-specific.
%
%% Inputs:
%
%  none
%
%% Outputs:
%
%  veh              Struct                      Vehicle parameters; see
%                                               coorbital.util.vehicleDefaults
%
%% Revision History:
%  Michael Casey                                                08/06/2026
%  Copyright 2026 Coorbital, Inc.
%% ------------------------ Begin Code Sequence ---------------------------

               veh = coorbital.util.vehicleDefaults();

%% Vehicle-specific overrides (PLACEHOLDER values):
          veh.mass = 900;              %kg
          veh.Sref = 0.75;             %m^2
            veh.CL = 0.35;
            veh.LD = 2.5;
end
