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
%                                               noseRadius (m) stagnation-point
%                                                  nose radius; see the note
%                                                  below -- nothing reads it yet
%
%% Note -- noseRadius is a placeholder for deferred work, not dead code:
%
%  No routine in the library consumes noseRadius today. It is here because the
%  design spec calls for a Sutton-Graves stagnation-point heating rate among
%  the derived quantities a trajectory reports, and heating was deliberately
%  cut from this milestone (see docs/DESIGN.md, "As built"). The field is
%  carried now, and asserted present and positive in test_constLD, so that the
%  vehicle interface does not have to change when the heating model lands.
%  Delete it only together with that plan, not as an unused-value cleanup.
%
%% Revision History:
%  Michael Casey                                                08/06/2026
%  Copyright 2026 Coorbital, Inc.
%% ------------------------ Begin Code Sequence ---------------------------

          veh.mass = 900;              %kg, PLACEHOLDER
          veh.Sref = 0.75;             %m^2, PLACEHOLDER
            veh.CL = 0.35;             %PLACEHOLDER
            veh.LD = 2.5;              %PLACEHOLDER, typical waverider range 2-4
    veh.noseRadius = 0.05;             %m, PLACEHOLDER; unused until Sutton-Graves
                                       %   heating arrives -- see the header note
end
