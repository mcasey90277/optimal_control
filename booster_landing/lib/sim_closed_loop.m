function out = sim_closed_loop(sol, ctrl, P, dsp)
% SIM_CLOSED_LOOP  Truth-model landing sim under TVLQR tracking.
%
% Plant = pdg_dynamics with dispersion multipliers (thrust/Isp bias, wind);
% controller = T*(t) - K(t) dx, magnitude-saturated to [Tmin, Tmax] with
% direction preserved, ALWAYS evaluated against the nominal model P (the
% flight computer does not see the true dispersed plant -- see the
% control_law note below). Integrates to the first of three terminal
% events: a real z=0 touchdown, a vertical arrest just above the pad
% followed by climb-away, or (safety net) propellant depletion to
% P.mdry -- see the ADAPTATION note.
%
% INPUTS:
%   sol,ctrl - Task 3 / Task 6 interfaces
%   P        - booster_params
%   dsp      - dispersions: .dr0 .dv0 [3x1], .thrust_scale .isp_scale
%              [scalars, def 1], .wind [3x1 m/s, needs P.drag.on] -- all
%              optional, defaults = nominal. thrust_scale and isp_scale
%              model two DIFFERENT physical error sources: thrust_scale
%              biases the delivered force vector itself (T is scaled
%              before pdg_dynamics sees it, so the mass-flow computed
%              there, |T|/(Isp*g0), scales WITH it -- a throttle/thrust
%              calibration error, force and propellant consumption move
%              together); isp_scale instead biases Pp.Isp, changing
%              mass-flow for the SAME delivered thrust (an efficiency
%              error). A real engine can have either or both.
% OUTPUTS:
%   out - .t(Mx1) .X(Mx7) .Tcmd(Mx3) .sat_frac (time-weighted duty cycle
%         on a thrust bound, trapz over out.t -- NOT a sample-count
%         fraction, which would be biased toward whatever region ode45
%         happens to step densely through)
%         .td struct: .r .v .m .miss .vtd .alt(=r(3)) .landed .stop
%           .landed - true only when the terminating event was a REAL
%                     z=0 touchdown (see ADAPTATION: an arrest reports
%                     vtd~0 by construction, since vz=0 IS the arrest
%                     event -- .landed is what distinguishes a genuine
%                     safe landing from that structurally optimistic
%                     number).
%           .stop   - 'touchdown' | 'arrest' | 'horizon' (no terminal
%                     event fired before the 1.5*sol.tf integration
%                     horizon; .td is then just the last integrated
%                     state, not any kind of landing)
%
% ADAPTATION FROM BRIEF (documented, forced by a genuine runtime issue):
% P.Tmin (338 kN) exceeds the vehicle's weight at every mass on this
% trajectory (m0*g0=294 kN, mdry*g0=251 kN), so a near-vertical, >=Tmin
% thrust vector always nets a small upward (decelerating) acceleration.
% Near tf the terminal-arc gain (~6.5e5) drives the saturated thrust
% direction close to vertical to null a (tiny) horizontal error, which can
% leave the closed loop unable to complete the last stretch of descent:
% instead of crossing z=0 it undergoes a vertical arrest followed by a
% slow climb-away, asymptoting to a stable near-ground hover. The brief's
% single z=0 touchdown event never fires in that regime -- confirmed by a
% diagnostic probe (ode45 OutputFcn) that showed z climbing 0.879->0.880 m
% over hundreds of steps of ~1e-6 s -- and ode45 stalls chasing that
% asymptote at tight tolerances. Two more terminal events are added:
% [2] the first zero-crossing of vz from negative to positive (closest
% approach / local min altitude -- fires at the arrest, before the climb
% asymptotes, giving a well-defined stop point instead of a stall) and
% [3] mass falling through P.mdry (a stalled integration at Tmax burns
% ~305.6 kg/s = Tmax/(Isp*g0); the ~7.85 s of slack between sol.tf and the
% 1.5*sol.tf horizon can burn far more than the touchdown mass margin
% above mdry, which would drive m<=0 in pdg_dynamics -- this event caps
% that before it happens). In a trajectory that truly lands, z reaches 0
% while still descending (vz<0), so event [1] fires first exactly as
% specified and out.td.landed=true. Event [2] or [3] firing instead means
% the closed loop did NOT complete the landing -- out.td.stop records
% which, and out.td.landed=false, so a caller can never mistake an arrest
% for a touchdown by reading .miss/.vtd alone. No Task-6 weight or gain
% was touched by this file.
%
% REFERENCES:
%   [1] docs/superpowers/specs/2026-08-08-booster-landing-design.md
%   [2] certify/certify_pdg.m G2 gate -- proves the OPEN-LOOP guidance
%       solution (the trajectory this controller tracks) lands for real,
%       to a ~1.6e-4 m continuous-time position residual; touchdown
%       failure in THIS file is a closed-loop tracking/saturation effect,
%       not evidence the underlying guidance trajectory doesn't land.
d = struct('dr0',zeros(3,1), 'dv0',zeros(3,1), 'thrust_scale',1, ...
           'isp_scale',1, 'wind',zeros(3,1));
if nargin >= 4
    fn = fieldnames(dsp);
    for k = 1:numel(fn), d.(fn{k}) = dsp.(fn{k}); end
end

Pp = P;  Pp.Isp = P.Isp * d.isp_scale;      % plant params differ from model
x0 = [P.r0 + d.dr0; P.v0 + d.dv0; P.m0];

% Tolerances loosened from 1e-8/1e-8 to 1e-6/1e-6 (ADAPTATION, task-7 fix
% report): the direction-preserving magnitude clamp in control_law makes
% plant_rhs non-smooth (a kink whenever |Traw| crosses a bound), and
% 1e-8 invites min-step stalls in exactly the near-saturated terminal arc
% this campaign cares about (and will matter more once this runs inside a
% Monte Carlo sweep over many dispersions).
oo = odeset('RelTol',1e-6, 'AbsTol',1e-6, 'Events', @(t,x) touchdown_event(t,x,P), ...
            'MaxStep', 0.25);
[tt, XX, te, xe, ie] = ode45(@(t,x) plant_rhs(t, x, ctrl, P, Pp, d), ...
                 [0, 1.5*sol.tf], x0, oo); %#ok<ASGLU>

out.t = tt;  out.X = XX;
Tc = zeros(numel(tt), 3);
for k = 1:numel(tt)
    Tc(k,:) = control_law(tt(k), XX(k,:).', ctrl, P).';
end
out.Tcmd = Tc;
Tmag = sqrt(sum(Tc.^2, 2));
onBound = double(Tmag > P.Tmax*0.999 | Tmag < P.Tmin*1.001);
out.sat_frac = trapz(tt, onBound) / max(tt(end) - tt(1), eps);

xend = XX(end,:).';
if isempty(ie)
    stop = 'horizon';                  % ran to 1.5*sol.tf, no event fired
elseif any(ie == 1)
    stop = 'touchdown';                % real z=0 crossing, still descending
else
    stop = 'arrest';                   % vz=0 proxy [2] or mdry cap [3]
end
out.td = struct('r', xend(1:3), 'v', xend(4:6), 'm', xend(7), ...
                'miss', sqrt(sum(xend(1:2).^2)), 'vtd', sqrt(sum(xend(4:6).^2)), ...
                'alt', xend(3), 'landed', strcmp(stop,'touchdown'), 'stop', stop);
end

function T = control_law(t, x, ctrl, P)
% TVLQR + magnitude saturation, direction preserved. K interp per column.
% Takes the NOMINAL model P (not the dispersed plant) deliberately: the
% flight computer's gain schedule and feedforward are built from the
% guidance model, not from ground truth it cannot observe -- callers must
% pass the model, dispersions only ever enter through the plant/dynamics
% side (plant_rhs below).
tq   = min(max(t, ctrl.tgrid(1)), ctrl.tgrid(end));
Kt   = zeros(3,7);
for r = 1:3
    Kt(r,:) = interp1(ctrl.tgrid.', squeeze(ctrl.K(r,:,:)).', tq, 'linear');
end
Traw = ctrl.Tnom(tq) - Kt * (x - ctrl.xnom(tq));
Tm   = sqrt(sum(Traw.^2));
T    = Traw * min(max(Tm, P.Tmin), P.Tmax) / max(Tm, 1e-9);
end

function xdot = plant_rhs(t, x, ctrl, P, Pp, d)
T    = control_law(t, x, ctrl, P) * d.thrust_scale;  % model P for the law
if Pp.drag.on                     % wind enters as airspeed shift
    xw = x;  xw(4:6) = x(4:6) - d.wind;
    xdot = pdg_dynamics(xw, T, Pp);   % dispersed plant Pp for TRUTH dynamics
    xdot(1:3) = x(4:6);           % kinematics use ground velocity
else
    xdot = pdg_dynamics(x, T, Pp);
end
end

function [val, isterm, dir_] = touchdown_event(~, x, P)
% Three terminal conditions (see ADAPTATION note above):
%  [1] z falling through 0        -- true touchdown, fires first in a
%                                     normal landing.
%  [2] vz rising through 0        -- closest-approach / arrest proxy;
%                                     fires only if the closed loop
%                                     arrests vertical velocity before z
%                                     reaches 0.
%  [3] m falling through P.mdry   -- propellant-depletion safety net for
%                                     a stalled integration that neither
%                                     [1] nor [2] catches quickly.
val    = [x(3); x(6); x(7) - P.mdry];
isterm = [1; 1; 1];
dir_   = [-1; 1; -1];
end
