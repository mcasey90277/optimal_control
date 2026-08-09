### Task 7: Closed-loop simulation (`sim_closed_loop`)

**Files:**
- Create: `lib/sim_closed_loop.m`
- Create: `tests/test_closed_loop_nominal.m`

**Interfaces:**
- Consumes: `sol`, `ctrl`, `pdg_dynamics`.
- Produces: `out = sim_closed_loop(sol, ctrl, P, dsp)`:
  - `dsp` dispersion struct (all optional, default 0): `.dr0(3x1) .dv0(3x1) .thrust_scale (1=nominal) .isp_scale .wind(3x1 const m/s, drag-frame — only felt when P.drag.on)`
  - `out.t`, `out.X` (Mx7), `out.Tcmd` (Mx3, post-saturation), `out.td` struct: `.r(3x1) .v(3x1) .m .miss` (horizontal distance to pad) `.vtd` (touchdown speed) — from the z=0 ode event
  - `out.sat_frac` — fraction of time on a thrust bound (diagnostic for the expected terminal-arc saturation)
- Controller law inside: `Traw = ctrl.Tnom(t) - K(t)*(x - ctrl.xnom(t))`, magnitude clamped to `[P.Tmin, P.Tmax]` (direction kept), `K(t)` linear-interpolated from `ctrl.K`. Truth dynamics apply `dsp.thrust_scale` and `dsp.isp_scale` (the plant differs from the model — that's the point).

- [ ] **Step 1: Write the failing test**

`tests/test_closed_loop_nominal.m`:

```matlab
% TEST_CLOSED_LOOP_NOMINAL  Zero-dispersion closed loop must reproduce the
% guidance: touchdown within 1 m and 0.1 m/s of the pad-at-rest target.
% Also: with a 50 m initial position offset, the tracker must still land
% inside the pad radius (the whole point of feedback).
%
% INPUTS: none   OUTPUTS: none (throws on failure)
here_ = fileparts(mfilename('fullpath'));
addpath(fullfile(here_, '..'));  setup_paths;

P    = booster_params();
sol  = solve_pdg_colloc(P, struct('N', 30));
ctrl = tvlqr_design(sol, P);

out0 = sim_closed_loop(sol, ctrl, P, struct());
assert(out0.td.miss < 1.0,  'nominal miss %.2f m', out0.td.miss);
assert(out0.td.vtd  < 0.1 + 1e-9, 'nominal touchdown speed %.3f', out0.td.vtd);

dsp  = struct('dr0', [50; -30; 0]);
out1 = sim_closed_loop(sol, ctrl, P, dsp);
assert(out1.td.miss < P.pad_radius, 'dispersed miss %.1f m', out1.td.miss);
assert(out1.td.vtd  < P.vtd_max,    'dispersed vtd %.2f m/s', out1.td.vtd);
fprintf('test_closed_loop_nominal PASS  (nom miss %.3f m, disp miss %.2f m)\n', ...
        out0.td.miss, out1.td.miss);
```

- [ ] **Step 2: Run test to verify it fails**

Run: `/Applications/MATLAB_R2025b.app/bin/matlab -batch "cd('/Users/msc/Desktop/optimal_control/booster_landing'); setup_paths; test_closed_loop_nominal"`
Expected: FAIL — `sim_closed_loop` undefined.

- [ ] **Step 3: Implement `lib/sim_closed_loop.m`**

```matlab
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
val = x(3);  isterm = 1;  dir_ = -1;   % z falling through 0
end
```

- [ ] **Step 4: Run test to verify it passes**

Same command. Expected: PASS with nominal miss well under 1 m. If nominal miss is large: the guidance interpolants and the plant integrate slightly differently — check `ctrl.Tnom` uses nodes+midpoints (it does), then check the event isn't firing early (initial z must be positive). If the dispersed case misses: gains too soft — revisit `opts.Q/R` in Task 6, document the change there.

- [ ] **Step 5: Commit**

```bash
cd /Users/msc/Desktop/optimal_control
git add booster_landing/lib/sim_closed_loop.m booster_landing/tests/test_closed_loop_nominal.m
git commit -m "booster_landing: closed-loop TVLQR landing sim with touchdown event

Co-Authored-By: Claude Fable 5 <noreply@anthropic.com>"
```

---

