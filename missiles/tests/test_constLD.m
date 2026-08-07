function test_constLD()
%% Purpose:
%
%  Verify the constant-L/D aerodynamic model returns the vehicle's CL, a CD
%  consistent with the requested lift-to-drag ratio, and ignores alpha and
%  Mach as documented.
%
%% Revision History:
%  Michael Casey                                                08/06/2026
%  Copyright 2026 Coorbital, Inc.
%% ------------------------ Begin Code Sequence ---------------------------

               veh = coorbital.util.vehicleDefaults();

%% Required vehicle fields exist and are positive:
          required = {'mass','Sref','CL','LD','noseRadius'};
    for k = 1:numel(required)
        assert(isfield(veh,required{k}),'vehicleDefaults missing %s',required{k});
        assert(veh.(required{k}) > 0,'%s must be positive',required{k});
    end

%% Hand-computed expected values from vehicleDefaults' documented placeholders
%% (CL = 0.35, LD = 2.5, so CD = CL/LD = 0.35/2.5 = 0.14) -- literals, not the
%% implementation's own formula fed back on itself -- evaluated at a positive
%% angle of attack:
          [CL,CD] = coorbital.aero.constLD(0.1,20,veh);
    assert(abs(CL - 0.35) < 1e-12,'CL should be the vehicle CL, 0.35');
    assert(abs(CD - 0.14) < 1e-12,'CD should be CL/LD = 0.35/2.5 = 0.14');
    assert(CL > 0,'CL must be positive for a positive angle of attack');
    assert(CD > 0,'CD must be positive for a positive angle of attack');

%% CL/CD reproduces the vehicle's specified L/D to tight tolerance, from the
%% returned values (this is the one property the glide trajectory hinges on):
    assert(abs(CL/CD - veh.LD) < 1e-12,'CL/CD should equal veh.LD');

%% alpha and Mach are ignored by this model -- a very different alpha and
%% Mach must return bit-identical outputs:
        [CL2,CD2] = coorbital.aero.constLD(0.4,6,veh);
    assert(CL2 == CL && CD2 == CD,'constLD must not depend on alpha or Mach');
end
