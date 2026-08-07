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

                 c = coorbital.util.missileConst();
             value = (x(1) - c.rE) - hStop;
        isterminal = 1;
         direction = -1;
end
