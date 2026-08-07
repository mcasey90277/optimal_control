function veh = vehicle_hgv()
%% Purpose:
%
%  Hypersonic glide vehicle parameters for the run_glide entry script.
%  Starts from the library defaults so there is one place to change a value
%  that should apply everywhere, and overrides only what is vehicle-specific.
%
%% Note -- this file currently changes nothing:
%
%  Every value set below deliberately COINCIDES with the generic placeholders
%  in coorbital.util.vehicleDefaults, because no airframe-specific numbers
%  have been supplied for this vehicle yet. The file is not redundant: it
%  exists so that real numbers have exactly one home when they arrive, and so
%  that run_glide reaches for a named vehicle rather than for the library
%  defaults directly. Read the assignments below as a schedule of what must be
%  replaced, not as a statement that anything has been.
%
%  There is no CD field to set. coorbital.aero.constLD derives CD = CL/LD, so
%  drag is fixed by the two values below and cannot fall out of sync with
%  them. A real airframe with a Mach-dependent drag polar needs a different
%  aero model, not an extra field here.
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
