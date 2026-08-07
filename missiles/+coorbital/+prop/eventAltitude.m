function [value,isterminal,direction] = eventAltitude(t,x,hStop)
%% Purpose:
%
%  ODE event that fires when the vehicle descends through a given altitude.
%  Terminal and one-sided, so a climbing trajectory passing the same altitude
%  does not stop the integration.
%
%% Inputs:
%
%  t                scalar                      Time (s). Unused.
%
%  x                [6 x 1]                     Glide state; uses x(1) = r
%
%  hStop            scalar                      Altitude to stop at (m)
%
%% Outputs:
%
%  value            scalar                      Altitude above hStop (m)
%
%  isterminal       scalar                      1, always terminal
%
%  direction        scalar                      -1, descending crossings only
%
%% Revision History:
%  Michael Casey                                                08/06/2026
%  Copyright 2026 Coorbital, Inc.
%% ------------------------ Begin Code Sequence ---------------------------

%% Self-demo:
if nargin == 0
                cD = coorbital.util.missileConst();
                xD = [cD.rE + 32e3; 0; 0; 4000; 0; deg2rad(90)];
    [vD,itD,dirD] = coorbital.prop.eventAltitude(0,xD,20e3);
    fprintf('value = %.1f m, isterminal = %d, direction = %+d\n',vD,itD,dirD);
    [value,isterminal,direction] = deal([]);
    return;
end

                 c = coorbital.util.missileConst();
             value = (x(1) - c.rE) - hStop;
        isterminal = 1;
         direction = -1;
end
