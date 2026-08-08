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
%     T    = veh.thrustVac - veh.Aexit*P,   required to be POSITIVE
%     mdot = veh.thrustVac/(veh.Isp*g0)
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
%% Note -- a nonpositive net thrust is an ERROR, not a shutdown:
%
%  The linear back-pressure debit is only valid while it is small. Driven far
%  enough the expression returns a nonpositive thrust, and this function then
%  REFUSES with coorbital:constThrust:invalidBackPressure.
%
%  It used to clamp instead -- max(...,0) on the thrust, and mdot zeroed with
%  it -- and that was wrong in the way that matters most, because it was
%  quiet. A choked motor does not stop burning because ambient back pressure
%  has made a SIMPLIFIED net-thrust expression go negative; the flow is set
%  upstream of the throat and knows nothing about the exit plane. So the
%  condition says the NOZZLE MODEL is out of its domain, not that combustion
%  ceased, and either way of resolving it silently is a lie: zeroing both
%  outputs freezes the mass state above burnout, after which the burnout event
%  never fires and coorbital.prop.phaseRun runs the phase out to tspan(end)
%  with the motor "running" and the propellant untouched; zeroing thrust alone
%  would drain propellant for nothing. Raise, and let the caller install a
%  separated-flow model if it wants to fly there.
%
%  THIS BRANCH IS UNREACHABLE IN EVERY SHIPPED CONFIGURATION and is expected
%  to stay that way. coorbital.util.boosterDefaults gives Aexit = 1.25 m^2
%  against thrustVac = 950 kN, so the largest possible debit is
%  1.25 * 101325 = 126.7 kN at sea level -- 13.3 percent of the vacuum thrust,
%  with a margin of better than 7x. Nothing in the library changes
%  numerically because of this edit.
%
%% Inputs:
%
%  t                [N x 1]                     Phase clock (s): the phase's
%                                               OWN tspan value, not rebased
%                                               to zero. See TWO CLOCKS in
%                                               coorbital.prop.phaseRun.
%                                               Accepted and ignored by this
%                                               model, whose thrust is
%                                               constant in time.
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
%  T                [N x 1]                     Delivered thrust (N), strictly
%                                               positive. Raises
%                                               coorbital:constThrust:invalidBackPressure
%                                               rather than clamp, see the
%                                               Note above
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
%  Michael Casey  Refuse a nonpositive net thrust, not clamp     08/07/2026
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
%% plane. t is deliberately unused -- the thrust is constant:
              Traw = veh.thrustVac - veh.Aexit.*P;
if any(Traw(:) <= 0) || any(~isfinite(Traw(:))) || ~isreal(Traw)
    error('coorbital:constThrust:invalidBackPressure', ...
        ['The simplified net thrust veh.thrustVac - veh.Aexit*P came out ' ...
         'nonpositive: thrustVac = %.6g N, Aexit = %.6g m^2, worst debit ' ...
         '%.6g N at P = %.6g Pa. The nozzle model is outside its domain of ' ...
         'validity -- a choked motor does not stop burning because the ' ...
         'back-pressure term overran -- so this is refused rather than ' ...
         'clamped. Install a separated-flow model, or reduce Aexit.'], ...
        veh.thrustVac,veh.Aexit,max(veh.Aexit.*P(:)),max(P(:)));
end
                 T = Traw;

%% Choked flow, so the mass flow is fixed by the vacuum thrust and vacuum Isp
%% and does not respond to P at all. Broadcast to the shape of T so a column
%% of pressures still returns a matching column of flows:
              mdot = (veh.thrustVac./(veh.Isp.*c.g0)).*ones(size(T));
end
