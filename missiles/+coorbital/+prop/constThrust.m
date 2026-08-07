function [T,mdot] = constThrust(t,P,veh)
%% Purpose:
%
%  Constant-mass-flow rocket motor with an ambient back-pressure correction on
%  thrust. The crudest usable propulsion model, and the one for which the
%  boost phase reduces exactly to the Tsiolkovsky rocket equation in vacuum
%  with gravity and drag switched off. Signature matches richer models
%  (throttled, tabulated, pressure-fed blowdown) so the environment struct can
%  swap them without touching the equations of motion.
%
%     T    = max(veh.thrustVac - veh.Aexit*P, 0)
%     mdot = veh.thrustVac/(veh.Isp*g0)          while T > 0, and 0 otherwise
%
%% Note -- the choked-nozzle assumption, which is why this is self-consistent:
%
%  The nozzle throat is choked, and the chamber pressure is held fixed, so the
%  propellant mass flow through the throat is FIXED: it is set upstream of the
%  throat and cannot be signalled by conditions downstream of it. Ambient back
%  pressure therefore changes the DELIVERED THRUST -- through the pressure
%  term Aexit*(Pexit - P) acting on the exit plane -- and not the flow. That is
%  exactly why deriving a constant mdot from the VACUUM thrust and the VACUUM
%  Isp is consistent: those two are quoted at the same ambient condition
%  (P = 0), their ratio is the flow, and the flow so obtained is the same flow
%  at every altitude. Pairing a vacuum thrust with a sea-level Isp, or letting
%  mdot follow the corrected thrust, would both break that.
%
%% Note -- the clamp zeroes the flow too, deliberately:
%
%  The linear back-pressure debit is only valid while it is small. Driven far
%  enough it would return a negative thrust, which the max() clamps to zero.
%  When that clamp fires the model has left its domain of validity, so the
%  motor is reported as not running at all: mdot goes to zero WITH the thrust.
%  A clamped engine that keeps consuming propellant is not a physical model,
%  and the inconsistency would silently corrupt the integrated mass history --
%  propellant would disappear with nothing to show for it, and the burnout
%  event would fire early. The two outputs are always zero together.
%
%% Inputs:
%
%  t                [N x 1]                     Time since the start of the
%                                               phase (s). Accepted and
%                                               ignored by this model, whose
%                                               thrust is constant in time.
%
%  P                [N x 1]                     Ambient static pressure (Pa).
%                                               Pressure rather than altitude
%                                               is deliberate: it is what the
%                                               physics needs, the equations
%                                               of motion already hold it from
%                                               env.atmos, and it keeps this
%                                               model independent of which
%                                               atmosphere is installed.
%
%  veh              Struct                      Vehicle/stack parameters; uses
%                                               the thrustVac (N), Isp (s,
%                                               VACUUM) and Aexit (m^2) fields
%                                               of coorbital.util.boosterDefaults
%
%% Outputs:
%
%  T                [N x 1]                     Delivered thrust (N),
%                                               non-negative
%
%  mdot             [N x 1]                     Propellant mass flow (kg/s),
%                                               POSITIVE by convention. The
%                                               equations of motion apply the
%                                               minus sign in dm/dt = -mdot.
%
%% References:
%   [1] Sutton, G.P., Biblarz, O., "Rocket Propulsion Elements," 9th ed.,
%       Wiley, 2017, Sec. 3.3 (thrust and the pressure term) and Sec. 3.2
%       (choked flow at the throat).
%
%% Revision History:
%  Michael Casey                                                08/07/2026
%  Copyright 2026 Coorbital, Inc.
%% ------------------------ Begin Code Sequence ---------------------------

%% Self-demo:
if nargin == 0
               bst = coorbital.util.boosterDefaults();
            [~,P0] = coorbital.atmos.expAtmos(0);
      [Tsl,mdotSl] = coorbital.prop.constThrust(0,P0,bst);
      [Tvc,mdotVc] = coorbital.prop.constThrust(0,0,bst);
    fprintf('sea level: T = %.1f N, mdot = %.4f kg/s\n',Tsl,mdotSl);
    fprintf('vacuum   : T = %.1f N, mdot = %.4f kg/s\n',Tvc,mdotVc);
    fprintf('burn time: %.3f s\n',bst.massProp./mdotVc);
    [T,mdot] = deal([]);
    return;
end

                 c = coorbital.util.missileConst();

%% Delivered thrust: vacuum thrust less the back-pressure debit on the exit
%% plane, clamped at zero. t is deliberately unused -- the thrust is constant:
                 T = max(veh.thrustVac - veh.Aexit.*P,0);

%% Choked flow, so the mass flow is fixed by the vacuum thrust and vacuum Isp
%% and does not respond to P at all -- except that a clamped motor is not
%% running, and then the flow is zero along with the thrust:
              mdot = (veh.thrustVac./(veh.Isp.*c.g0)).*(T > 0);
end
