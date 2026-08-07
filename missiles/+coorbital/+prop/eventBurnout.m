function [value,isterminal,direction] = eventBurnout(t,x,mBurnout)
%% Purpose:
%
%  ODE event that fires when the boost phase's propellant is exhausted:
%  the state mass x(7) descends through the TOTAL burnout mass, payload plus
%  spent booster structure. Terminal and one-sided in the descending
%  direction. Mass only ever decreases during a single burn, so direction
%  could not fire the other way regardless, but the one-sided form documents
%  that intent and guards a future engine-restart phase where mass could hold
%  flat before resuming its decrease.
%
%  The burnout mass is taken as an explicit scalar argument rather than
%  reached for inside the event via veh and bst structs. Comparing x(7)
%  against bst.massDry alone would treat the payload as part of the
%  propellant budget and burn straight through it before stopping; the
%  caller must pass veh.mass + bst.massDry (see
%  docs/plan_2026-08-07_boost_descent_chain.md, "Mass bookkeeping").
%
%% Inputs:
%
%  t                scalar                      Time (s). Unused.
%
%  x                [7 x 1]                     Boost state; uses x(7) = m
%
%  mBurnout         scalar                      Total burnout mass (kg),
%                                               veh.mass + bst.massDry
%
%% Outputs:
%
%  value            scalar                      Mass above mBurnout (kg)
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
                veh = coorbital.util.vehicleDefaults();
                bst = coorbital.util.boosterDefaults();
                 cD = coorbital.util.missileConst();
                 xD = [cD.rE + 5e3; 0; 0; 300; deg2rad(60); deg2rad(90); ...
                      veh.mass + bst.massDry + 100];
           mBurnoutD = veh.mass + bst.massDry;
      [vD,itD,dirD] = coorbital.prop.eventBurnout(0,xD,mBurnoutD);
    fprintf('value = %.1f kg, isterminal = %d, direction = %+d\n',vD,itD,dirD);
    [value,isterminal,direction] = deal([]);
    return;
end

             value = x(7) - mBurnout;
        isterminal = 1;
         direction = -1;
end
