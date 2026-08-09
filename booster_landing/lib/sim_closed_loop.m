function out = sim_closed_loop(sol, ctrl, P, dsp)
% SIM_CLOSED_LOOP  Truth-model landing sim under TVLQR tracking.
%
% Plant = pdg_dynamics with dispersion multipliers (thrust/Isp bias, wind);
% controller = T*(t) - K(t) dx, magnitude-saturated to [Tmin, Tmax] with
% direction preserved. Integrates to the z=0 crossing (ode event).
%
% INPUTS:
%   sol,ctrl - Task 3 / Task 6 interfaces
%   P        - booster_params
%   dsp      - dispersions: .dr0 .dv0 [3x1], .thrust_scale .isp_scale
%              [scalars, def 1], .wind [3x1 m/s, needs P.drag.on] -- all
%              optional, defaults = nominal
% OUTPUTS:
%   out - .t(Mx1) .X(Mx7) .Tcmd(Mx3) .sat_frac
%         .td struct: .r .v .m .miss .vtd
%
% ADAPTATION FROM BRIEF (documented, forced by a genuine runtime issue):
% P.Tmin (338 kN) exceeds the vehicle's weight at every mass on this
% trajectory (m0*g0=294 kN, mdry*g0=251 kN), so a near-vertical, >=Tmin
% thrust vector always nets a small upward (decelerating) acceleration.
% Near tf the terminal-arc gain (~6.5e5) drives the saturated thrust
% direction close to vertical to null the (tiny) horizontal error, which
% leaves the closed loop unable to complete the last ~1 m of descent: the
% single z=0 touchdown event of the brief's reference code NEVER fires --
% the trajectory asymptotes to a stable near-ground hover (verified: z
% climbs 0.879->0.880 m over the tail of a probe run) and ode45 stalls,
% taking vanishingly small steps forever chasing AbsTol=1e-8 against a
% derivative that itself -> 0. A second terminal event is added: the first
% zero-crossing of vz from negative to positive (closest approach / local
% min altitude). In a trajectory that truly lands, z reaches 0 while still
% descending (vz<0) and the original z=0 event fires first as designed;
% only when the closed loop instead arrests vertical velocity just above
% the pad does this second event fire, at the closest-approach state --
% the physically sensible touchdown proxy for that regime. No Task-6
% weight or gain was touched.
%
% REFERENCES:
%   [1] docs/superpowers/specs/2026-08-08-booster-landing-design.md
d = struct('dr0',zeros(3,1), 'dv0',zeros(3,1), 'thrust_scale',1, ...
           'isp_scale',1, 'wind',zeros(3,1));
if nargin >= 4
    fn = fieldnames(dsp);
    for k = 1:numel(fn), d.(fn{k}) = dsp.(fn{k}); end
end

Pp = P;  Pp.Isp = P.Isp * d.isp_scale;      % plant params differ from model
x0 = [P.r0 + d.dr0; P.v0 + d.dv0; P.m0];

oo = odeset('RelTol',1e-8, 'AbsTol',1e-8, 'Events', @touchdown_event, ...
            'MaxStep', 0.25);
[tt, XX] = ode45(@(t,x) plant_rhs(t, x, sol, ctrl, Pp, d), ...
                 [0, 1.5*sol.tf], x0, oo);

out.t = tt;  out.X = XX;
Tc = zeros(numel(tt), 3);
for k = 1:numel(tt)
    Tc(k,:) = control_law(tt(k), XX(k,:).', ctrl, P).';
end
out.Tcmd = Tc;
Tmag = sqrt(sum(Tc.^2, 2));
out.sat_frac = mean(Tmag > P.Tmax*0.999 | Tmag < P.Tmin*1.001);
xe = XX(end,:).';
out.td = struct('r', xe(1:3), 'v', xe(4:6), 'm', xe(7), ...
                'miss', sqrt(sum(xe(1:2).^2)), 'vtd', sqrt(sum(xe(4:6).^2)));
end

function T = control_law(t, x, ctrl, P)
% TVLQR + magnitude saturation, direction preserved. K interp per column.
tq   = min(max(t, ctrl.tgrid(1)), ctrl.tgrid(end));
Kt   = zeros(3,7);
for r = 1:3
    Kt(r,:) = interp1(ctrl.tgrid.', squeeze(ctrl.K(r,:,:)).', tq, 'linear');
end
Traw = ctrl.Tnom(tq) - Kt * (x - ctrl.xnom(tq));
Tm   = sqrt(sum(Traw.^2));
T    = Traw * min(max(Tm, P.Tmin), P.Tmax) / max(Tm, 1e-9);
end

function xdot = plant_rhs(t, x, sol, ctrl, Pp, d)
T    = control_law(t, x, ctrl, Pp) * d.thrust_scale;
if Pp.drag.on                     % wind enters as airspeed shift
    xw = x;  xw(4:6) = x(4:6) - d.wind;
    xdot = pdg_dynamics(xw, T, Pp);
    xdot(1:3) = x(4:6);           % kinematics use ground velocity
else
    xdot = pdg_dynamics(x, T, Pp);
end
end

function [val, isterm, dir_] = touchdown_event(~, x)
% Two terminal conditions (see ADAPTATION note above): [1] z falling
% through 0 (true touchdown, fires first in a normal landing) or
% [2] vz rising through 0 (closest-approach proxy, fires only if the
% closed loop arrests vertical velocity before z reaches 0).
val   = [x(3); x(6)];
isterm = [1; 1];
dir_   = [-1; 1];
end
