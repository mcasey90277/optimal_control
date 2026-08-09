function out = sim_closed_loop(sol, ctrl, P, dsp)
% SIM_CLOSED_LOOP  Truth-model landing sim under TVLQR tracking.
%
% Plant = pdg_dynamics with dispersion multipliers (thrust/Isp bias, wind);
% controller = T*(t) - K(t) dx, saturated into the [Tmin, Tmax] annulus by
% allocate_thrust (direction-preserving, except that the upper bound is
% served vertical-first -- see that function), ALWAYS evaluated against
% the nominal model P (the flight computer does not see the true dispersed
% plant -- see the control_law note below). The reference trajectory and
% gains are indexed by the vehicle's ALTITUDE, not wall-clock time -- see
% the ALTITUDE-INDEXED GUIDANCE note. Integrates to the first
% of three terminal events: a real z=0 touchdown, a vertical arrest just
% above the pad followed by climb-away, or (safety net) propellant
% depletion to P.mdry -- see the ADAPTATION note.
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
% ALTITUDE-INDEXED GUIDANCE (task-7b, 2026-08-08 -- supersedes the round-4
% "terminal-phase altitude schedule", which was a partial version of this
% same idea applied only below P.zTermBand):
%
% The tracker indexes the reference trajectory by the vehicle's own
% ALTITUDE, not by wall-clock time. At altitude z it flies toward
% (r*_xy(z), v*(z)) with feedforward T*(t*(z)) and gains K(t*(z)), where
% t*(z) is the nominal time at that altitude.
%
% WHY (measured, task-7b): a min-fuel landing trajectory has NO time
% margin -- it ends exactly at z=0 with the tanks near dry -- so a pure
% altitude dispersion is unrecoverable under time-indexed tracking. With
% dr0=[0;0;-50] and a well-tuned time-indexed vertical loop the vehicle
% tracked vznom(t) essentially perfectly (vz=-169.38 vs -169.84 at t=6.76),
% which is exactly the problem: tracking time perfectly means STAYING 50 m
% low all the way down and hitting the ground 50 m early, at -31.8 m/s,
% while the guidance still had 2 s of braking left to do. Tightening the
% loop made it WORSE (9.0 -> 43.4 m/s touchdown), because the better it
% tracked the wrong schedule the harder it drove into the ground. The
% round-4/5 code only escaped this by having a vertical loop so weak
% (bandwidth ~0.12 rad/s) that velocity SAGGED toward the altitude-
% appropriate profile on its own -- accidentally doing altitude indexing,
% badly, in the last 150 m.
%
% Under altitude indexing the same dispersion is nearly a non-event: at
% z=1950 the reference is v*(1950) = -179.83 m/s against the vehicle's
% -180.00, an error of 0.17 m/s. The vehicle simply flies the profile it
% is actually on and brakes where the guidance says to brake IN ALTITUDE,
% which is the invariant the fuel-optimal solution actually encodes.
%
% Three structural simplifications fall out, all of which delete
% special-case code rather than add it:
%  * the ALTITUDE-ERROR CHANNEL VANISHES identically (xref(3) := x(3)), so
%    round 4's explicit Kt(:,3)=0 terminal zeroing is no longer needed --
%    it is now true by construction, everywhere, not below a band.
%  * the "PAST-tf FREEZE" cannot occur. K was previously indexed by
%    wall-clock t and clamped at tgrid(end)=tf, where P(tf)=Qf is diagonal
%    and B's position rows are zero, making K(tf)'s position columns
%    structurally zero -- so any run lasting past nominal tf silently flew
%    with NO position feedback at all. Indexing by t*(z) means the gain
%    reaches its terminal value only as the vehicle reaches the ground.
%  * P.zTermBand and its hand-picked 150 m are no longer used by this file.
%
% Limitation, asserted rather than assumed: altitude indexing needs a
% strictly monotone descent. That holds for this campaign (vz<0 throughout
% the guidance solution) and is checked below; a trajectory with a hover or
% ascent segment would need a different schedule variable (e.g. arclength).
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
% for a touchdown by reading .miss/.vtd alone.
%
% REFERENCES:
%   [1] docs/superpowers/specs/2026-08-08-booster-landing-design.md
%   [2] certify/certify_pdg.m G2 gate -- proves the OPEN-LOOP guidance
%       solution (the trajectory this controller tracks) lands for real,
%       to a ~4e-4 m continuous-time position residual; touchdown failure
%       in THIS file is a closed-loop tracking/saturation effect, not
%       evidence the underlying guidance trajectory doesn't land.
d = struct('dr0',zeros(3,1), 'dv0',zeros(3,1), 'thrust_scale',1, ...
           'isp_scale',1, 'wind',zeros(3,1));
if nargin >= 4
    fn = fieldnames(dsp);
    for k = 1:numel(fn), d.(fn{k}) = dsp.(fn{k}); end
end

Pp = P;  Pp.Isp = P.Isp * d.isp_scale;      % plant params differ from model
x0 = [P.r0 + d.dr0; P.v0 + d.dv0; P.m0];

%% ALTITUDE-INDEXED reference schedule, built ONCE from sol (see the
%% ALTITUDE-INDEXED GUIDANCE note above). Columns of REF, per altitude:
%% [rx*, ry*, vx*, vy*, vz*, m*, t*]. z* needs no column -- it IS the
%% index, which is exactly why the altitude-error channel vanishes.
[zSorted, sortIdx] = sort(sol.X(3,:));          % ascending for interp1
REF = [sol.X(1:2,sortIdx); sol.X(4:7,sortIdx); sol.t(sortIdx)].';
assert(all(diff(zSorted) > 0), 'sim_closed_loop:zNotMonotone', ...
    ['guidance altitude is not strictly monotone -- altitude indexing ' ...
     'needs a monotone descent (true for this campaign: vz<0 throughout)']);
ct.zLo = zSorted(1);  ct.zHi = zSorted(end);
ct.ref = @(zq) interp1(zSorted, REF, min(max(zq, zSorted(1)), zSorted(end)), 'pchip').';
ct.tf  = sol.tf;

% Tolerances loosened from 1e-8/1e-8 to 1e-6/1e-6 (ADAPTATION, task-7 fix
% report): the direction-preserving magnitude clamp in control_law makes
% plant_rhs non-smooth (a kink whenever |Traw| crosses a bound), and
% 1e-8 invites min-step stalls in exactly the near-saturated terminal arc
% this campaign cares about (and will matter more once this runs inside a
% Monte Carlo sweep over many dispersions).
oo = odeset('RelTol',1e-6, 'AbsTol',1e-6, 'Events', @(t,x) touchdown_event(t,x,P), ...
            'MaxStep', 0.25);
[tt, XX, te, xe, ie] = ode45(@(t,x) plant_rhs(t, x, ctrl, P, Pp, d, ct), ...
                 [0, 1.5*sol.tf], x0, oo); %#ok<ASGLU>

out.t = tt;  out.X = XX;
Tc = zeros(numel(tt), 3);
for k = 1:numel(tt)
    Tc(k,:) = control_law(tt(k), XX(k,:).', ctrl, P, ct).';
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

function T = control_law(t, x, ctrl, P, ct)
% TVLQR + magnitude saturation, direction preserved. K interp per column.
% Takes the NOMINAL model P (not the dispersed plant) deliberately: the
% flight computer's gain schedule and feedforward are built from the
% guidance model, not from ground truth it cannot observe -- callers must
% pass the model, dispersions only ever enter through the plant/dynamics
% side (plant_rhs below).
%
% ct - ALTITUDE-INDEXED reference schedule (see the note in
%      sim_closed_loop's header): .ref(z) -> [rx*;ry*;vx*;vy*;vz*;m*;t*],
%      .zLo/.zHi clamp range, .tf.
%
% NOTE the wall-clock t argument is deliberately UNUSED for the reference:
% the schedule variable is the vehicle's own altitude x(3). It is kept in
% the signature because ode45 passes it and because a future
% time-referenced variant would need it.
ref  = ct.ref(x(3));                 % reference at the ACTUAL altitude
tqK  = min(max(ref(7), ctrl.tgrid(1)), ctrl.tgrid(end));   % nominal t at z
Kt   = zeros(3,7);
for r = 1:3
    Kt(r,:) = interp1(ctrl.tgrid.', squeeze(ctrl.K(r,:,:)).', tqK, 'linear');
end
xref = [ref(1:2); x(3); ref(3:6)];   % z-reference IS the actual altitude,
                                     % so the altitude-error channel is
                                     % identically zero by construction --
                                     % no gain zeroing special case needed.

Traw = ctrl.Tnom(tqK) - Kt * (x - xref);
T    = allocate_thrust(Traw, P);
end

function T = allocate_thrust(Traw, P)
% ALLOCATE_THRUST  Annulus saturation with VERTICAL PRIORITY on the upper
% bound.
%
% The brief's law is "magnitude clamped to [Tmin,Tmax], direction kept".
% That is exactly what this returns whenever the raw command already fits
% under Tmax. That is the COMMON case but NOT the whole acceptance case --
% an earlier draft of this comment claimed the over-Tmax branch was a
% no-op after round 5's gain redesign, and measurement says otherwise:
% the raw command exceeds Tmax on 1.9% of the dispersed run's flight time,
% peaking at 968 kN (+14.6% over the bound), and on 2.1% of the nominal
% run's -- though the nominal excursions are a knife-edge artifact of the
% guidance riding exactly on Tmax (peak 850 kN, only +0.6% over), not a
% real demand. The dispersed 14.6% overshoot is real, and it is the whole
% reason this branch's policy matters.
% It differs ONLY in the over-Tmax branch, where preserving direction
% means scaling the VERTICAL component down in lockstep with a
% (necessarily unaffordable) lateral demand -- i.e. paying for a lateral
% correction with braking authority, at the one moment the vehicle can
% least afford it. Large lateral offsets are precisely where that matters:
% a 212 m offset commands >2 MN of lateral thrust, the direction-preserving
% clamp then points the (Tmax-limited) vector nearly horizontal, and the
% vehicle free-falls. Vertical priority instead keeps the commanded
% vertical component and spends only the REMAINING annulus radius,
% sqrt(Tmax^2 - Tz^2), on lateral -- the correct mission priority (miss the
% pad before you crash into it), and a no-op wherever the command fits.
%
% The lower bound (|T| >= Tmin) stays direction-preserving: Tmin is a
% throttle floor, not a budget to allocate.
%
% NEGATIVE-Tz GUARD (round-5 re-review): vertical priority is only
% meaningful when the vertical demand is a BRAKING demand. A commanded
% Traw(3) < 0 means "accelerate the descent" -- reachable whenever the
% vehicle sits ABOVE the guidance's v(z) schedule, a regime task 8's Monte
% Carlo will draw routinely -- and there is no braking authority to
% protect there. Serving it first would be actively perverse: the
% unguarded form turned Traw=[1e6;1e6;-1e6] into T=[0;0;-845 kN], i.e.
% spend the ENTIRE annulus thrusting downward and zero the lateral
% correction outright. So the priority branch is taken only for
% Traw(3) >= 0; a negative commanded Tz falls back to plain
% direction-preserving scaling, which keeps the lateral share intact.
%
% |T| == Tmax identically on BOTH over-budget branches, so no Tmin restore
% is needed there (an earlier draft carried one; it was dead code).
% Priority branch: if Traw(3) <= Tmax then Tz=Traw(3) and |Traw|>Tmax
% forces |Txy_raw| > lat, so Txy is always scaled and |T| = sqrt(lat^2 +
% Tz^2) = Tmax; if Traw(3) > Tmax then Tz=Tmax, lat=0 and |T| = Tmax.
% Fallback branch: direction-preserving scaling of a vector with
% |Traw| > Tmax lands exactly on Tmax.
%
% INPUTS:  Traw - raw commanded thrust [3x1 N];  P - booster_params
% OUTPUTS: T    - annulus-feasible thrust [3x1 N]
Tm = sqrt(sum(Traw.^2));
if Tm > P.Tmax && Traw(3) >= 0
    Tz  = min(Traw(3), P.Tmax);                    % vertical served first
    lat = sqrt(max(P.Tmax^2 - Tz^2, 0));           % annulus radius left over
    Txy = Traw(1:2);  n = sqrt(sum(Txy.^2));
    if n > lat, Txy = Txy * (lat / n); end
    T   = [Txy; Tz];
else
    T = Traw * min(max(Tm, P.Tmin), P.Tmax) / max(Tm, 1e-9);
end
end

function xdot = plant_rhs(t, x, ctrl, P, Pp, d, ct)
T    = control_law(t, x, ctrl, P, ct) * d.thrust_scale;  % model P for the law
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
