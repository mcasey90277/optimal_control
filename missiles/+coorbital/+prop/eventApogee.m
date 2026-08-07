function [value,isterminal,direction] = eventApogee(t,x)
%% Purpose:
%
%  ODE event that fires at apogee: the flight-path angle x(5) crossing zero
%  DESCENDING. Terminal and one-sided, so a vehicle that is still climbing --
%  gamma positive and moving away from zero, not toward it -- cannot trigger
%  it. Works on any state whose fifth component is gamma, six-state glide or
%  seven-state boost alike, since only x(5) is read.
%
%% Inputs:
%
%  t                scalar                      Time (s). Unused.
%
%  x                [>=5 x 1]                   State; uses x(5) = gamma
%
%% Outputs:
%
%  value            scalar                      Flight-path angle gamma (rad)
%
%  isterminal       scalar                      1, always terminal
%
%  direction        scalar                      -1, descending crossings only
%
%% Revision History:
%  Michael Casey                                                08/07/2026
%  Copyright 2026 Coorbital, Inc.
%% ------------------------ Begin Code Sequence ---------------------------

%% Self-demo:
if nargin == 0
                cD = coorbital.util.missileConst();
                xD = [cD.rE + 40e3; 0; 0; 4000; deg2rad(-2); deg2rad(90)];
    [vD,itD,dirD] = coorbital.prop.eventApogee(0,xD);
    fprintf('value = %.4f rad, isterminal = %d, direction = %+d\n',vD,itD,dirD);
    [value,isterminal,direction] = deal([]);
    return;
end

             value = x(5);
        isterminal = 1;
         direction = -1;
end
