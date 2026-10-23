# Drone Dynamics Library Implementation Plan (Stage 0)

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** One tested home for three 3D drone models (kinematic, point mass, 6-DOF rigid-body quadrotor) in control-affine form, plus the shared geometry (layouts, pairs, clearance, the paper's log barrier with analytic derivatives), with an independent physics oracle, so Stages 1-3 solve on a plant nothing can quietly get wrong.

**Architecture:** `optimal_control_examples/drone/drone_dynamics/` mirrors `cartpole_common`: `drone_params` / `drone_field` / `drone_state_jac` (generated from the field itself by the Symbolic Math Toolbox) / `euler_zyx`, with `tests/test_drone_physics.m` as the oracle derived from mechanics, not from the implemented rows. Geometry helpers (`drone_layouts`, `pair_list`, `clearance_terms`, `barrier_terms`, `lower_bound`, `hermite_resample`) are objective-independent and complex-step safe. A folder runner `drone/run_tests.m` is registered in `optimal_control_examples/run_all_tests.m`.

**Tech Stack:** MATLAB R2026a (`/Applications/MATLAB_R2026a.app/bin/matlab`), Symbolic Math Toolbox (`matlabFunction`), `ode113`. No CasADi, no oclib, no orbit_transfer in this stage.

**Spec:** `docs/superpowers/specs/2026-10-23-drone-multiagent-design.md` (Stage 0 section and the "What the paper's problem is" table).

## Global Constraints

- House MATLAB style: every function has the header quartet (`%% Purpose:` / `%% Inputs:` / `%% Outputs:` / `%% Revision History:` with `M. Casey  (c) 10/23/2026` and `Copyright Coorbital Inc.`), `if nargin == 0` self-demo in library functions, never `i`/`j` as variables, never `norm` (write `sqrt(sum(v.^2))`), no `%#ok`, no emojis.
- Nothing in `drone_dynamics/` knows a cost, horizon, boundary condition or solver. Nothing references `orbit_transfer` or pumpkyn.
- `drone_field` and `barrier_terms` are complex-step safe: no `abs`, `norm`, `max`, `min`, and branches only on `real(.)`.
- Every function adds paths relative to `mfilename('fullpath')`; no absolute paths in shipped files.
- Constants (verbatim): `m = 0.027`, `I = [1.395e-5; 1.436e-5; 2.173e-5]`, `l = 0.046`, `kappa = 0.005964`, `g = 9.81`, `fMax = 0.15`, `thetaMax = pi/3`; barrier `eps = 1e-5`, `l = 1e-4`.
- Tests return a logical, print `PASS`/`FAIL` lines, and are run through `drone/run_tests.m`, which THROWS on any failure.
- Run MATLAB tests through the MCP `run_matlab_file`/`evaluate_matlab_code` tools or `matlab -batch`; never commit from a subagent; commit messages end with `Co-Authored-By: Claude Fable 5.1 <noreply@anthropic.com>`; do not stage the unrelated `missiles/` deletions.
- Generated files (`drone_state_jac_rigid6dof_gen.m`) are not hand-edited; the generator is re-run instead.

## Review Focus

1. `rigid6dof` near pitch 90°: `W(Θ)` divides by `cos θ`; the field must stay finite at `thetaMax = 60°` and `drone_params` must expose `thetaMax` for the solvers' box. (Task 3 test: finite field at `θ = thetaMax`, non-finite only at `θ = π/2`.)
2. Barrier at `h ≤ 0` (interpenetration): the quadratic branch must keep `U`, `U'`, `U''` finite, with no `log` of a non-positive number. (Task 8 test at `h = −0.01` and `h = 0`.)
3. `M = 1`: `pair_list(1)` is `[2 x 0]`, `barrier_terms` returns `U = 0`, zero gradient, empty-sparse Hessian, and `clearance_terms` returns `[0 x N]` without erroring. (Tasks 6-8.)
4. Ring layouts with `M` not divisible by the ring count: `drone_layouts` must refuse by identifier, not produce a lopsided layout. (Task 6.)
5. Complex inputs to `clearance_terms`/`barrier_terms` (a caller complex-stepping through them): the branch is on `real(h)` and `sqrt` is never applied to a magnitude taken with `abs`. (Task 8 test: the complex step of `U` matches `gU`.)

---

### Task 1: `drone_params`

**Files:**
- Create: `optimal_control_examples/drone/drone_dynamics/drone_params.m`
- Test: `optimal_control_examples/drone/drone_dynamics/tests/test_drone_params.m`

**Interfaces:**
- Produces: `p = drone_params(model)` with `model ∈ {'kinematic','pointmass','rigid6dof'}`; always `.model .nx .nu`; `pointmass` adds `.g .aMax`; `rigid6dof` adds `.m .I [3x1] .l .kappa .g .fMax .thetaMax`.

- [ ] **Step 1: Write the failing test**

```matlab
function ok = test_drone_params()
%% Purpose:
%
%   drone_params returns the specified constants for each model, with the
%   dimensions the field will read, and refuses an unknown model by
%   identifier. The constants are the spec's verbatim (Luis & Le Ny 2016,
%   Table 1; Forster 2015); aMax is derived from them.
%
%% Inputs:
%
%  none
%
%% Outputs:
%
%  ok                       logical                 All checks passed
%
%% Revision History:
%  M. Casey                                                   (c) 10/23/2026
%  Copyright Coorbital Inc.
%% ------------------------ Begin Code Sequence ---------------------------

ok = true;
addpath(fileparts(fileparts(mfilename('fullpath'))));

k = drone_params('kinematic');
ok = chk(ok, strcmp(k.model, 'kinematic') && k.nx == 3 && k.nu == 3, 'kinematic: nx = 3, nu = 3');
ok = chk(ok, isequal(sort(fieldnames(k)), sort({'model'; 'nx'; 'nu'})), 'kinematic: exactly model, nx, nu');

q = drone_params('pointmass');
ok = chk(ok, q.nx == 6 && q.nu == 3 && q.g == 9.81, 'pointmass: nx = 6, nu = 3, g = 9.81');
ok = chk(ok, abs(q.aMax - (4*0.15/0.027 - 9.81)) < 1e-12, 'pointmass: aMax = 4 fMax / m - g');

r = drone_params('rigid6dof');
want = struct('m', 0.027, 'l', 0.046, 'kappa', 0.005964, 'g', 9.81, 'fMax', 0.15, 'thetaMax', pi/3);
f = fieldnames(want);
for n = 1:numel(f)
    ok = chk(ok, isfield(r, f{n}) && r.(f{n}) == want.(f{n}), sprintf('rigid6dof: %s = %g', f{n}, want.(f{n})));
end
ok = chk(ok, isequal(size(r.I), [3 1]) && isequal(r.I, [1.395e-5; 1.436e-5; 2.173e-5]), 'rigid6dof: I is the 3x1 Crazyflie diagonal');
ok = chk(ok, r.nx == 12 && r.nu == 4, 'rigid6dof: nx = 12, nu = 4');
ok = chk(ok, strcmp(r.model, 'rigid6dof') && strcmp(drone_params('RIGID6DOF').model, 'rigid6dof'), 'model name is case-insensitive on input, lower-case on output');

try
    drone_params('hexacopter');
    ok = chk(ok, false, 'unknown model refused');
catch err
    ok = chk(ok, strcmp(err.identifier, 'drone_params:model'), 'unknown model refused by identifier drone_params:model');
end

if ok, fprintf('TEST_DRONE_PARAMS: ALL PASS\n'); else, fprintf('TEST_DRONE_PARAMS: FAILURE\n'); end
end

% ------------------------------------------------------------------------
function ok = chk(ok, cond, label)
%% Purpose:
%
%   Accumulate one labelled pass/fail.
%
if cond, fprintf('  PASS  %s\n', label);
else,    fprintf('  FAIL  %s\n', label);  ok = false;
end
end
```

- [ ] **Step 2: Run it to verify it fails**

Run (MCP `run_matlab_file` on the test, or): `matlab -batch "cd('optimal_control_examples/drone/drone_dynamics/tests'); test_drone_params"`
Expected: error `Unrecognized function or variable 'drone_params'`.

- [ ] **Step 3: Implement**

```matlab
function p = drone_params(model)
%% Purpose:
%
%   THE physical constants of the three drone models, one home. Every
%   solver and every test reads them from here; nothing else holds a copy.
%
%   kinematic   p = position, u = velocity              (the paper's drone)
%   pointmass   [p; v], u = thrust acceleration         (the planner's drone)
%   rigid6dof   [p; v; Theta; omega], u = rotor thrusts (the vehicle)
%
%   rigid6dof constants are Crazyflie-class: Luis & Le Ny, "Design of a
%   Trajectory Tracking Controller for a Nanoquadcopter", 2016, Table 1, and
%   Forster, "System Identification of the Crazyflie 2.0 Nano Quadrocopter",
%   ETH 2015. thetaMax is the pitch/roll box every direct solve imposes: the
%   ZYX Euler kinematics are singular at 90 degrees.
%
%% Inputs:
%
%  model                    char                    'kinematic' |
%                                                   'pointmass' |
%                                                   'rigid6dof'
%
%% Outputs:
%
%  p                        struct                  .model .nx .nu always;
%                                                   pointmass: .g .aMax;
%                                                   rigid6dof: .m .I [3x1]
%                                                   .l .kappa .g .fMax
%                                                   .thetaMax
%
%% Revision History:
%  M. Casey                                                   (c) 10/23/2026
%  Copyright Coorbital Inc.
%% ------------------------ Begin Code Sequence ---------------------------

if nargin == 0
   %Demo: the three parameter sets:
     disp(drone_params('kinematic'));
     disp(drone_params('pointmass'));
     disp(drone_params('rigid6dof'));
     return
end

m      = 0.027;
fMax   = 0.15;
g      = 9.81;
switch lower(model)
    case 'kinematic'
        p = struct('model', 'kinematic', 'nx', 3, 'nu', 3);
    case 'pointmass'
        p = struct('model', 'pointmass', 'nx', 6, 'nu', 3, 'g', g, ...
                   'aMax', 4*fMax/m - g);
    case 'rigid6dof'
        p = struct('model', 'rigid6dof', 'nx', 12, 'nu', 4, 'm', m, ...
                   'I', [1.395e-5; 1.436e-5; 2.173e-5], 'l', 0.046, ...
                   'kappa', 0.005964, 'g', g, 'fMax', fMax, 'thetaMax', pi/3);
    otherwise
        error('drone_params:model', 'unknown model ''%s''', model);
end
end
```

- [ ] **Step 4: Run the test**

Expected: every line `PASS`, `TEST_DRONE_PARAMS: ALL PASS`.

- [ ] **Step 5: Commit**

```bash
cd ~/Desktop/optimal_control
/Library/Developer/CommandLineTools/usr/bin/git add optimal_control_examples/drone/drone_dynamics/drone_params.m optimal_control_examples/drone/drone_dynamics/tests/test_drone_params.m
/Library/Developer/CommandLineTools/usr/bin/git commit -m "drone_dynamics: drone_params, the three models' constants in one home, with its test

Co-Authored-By: Claude Fable 5.1 <noreply@anthropic.com>"
```

---

### Task 2: `euler_zyx`

**Files:**
- Create: `optimal_control_examples/drone/drone_dynamics/euler_zyx.m`
- Test: `optimal_control_examples/drone/drone_dynamics/tests/test_euler_zyx.m`

**Interfaces:**
- Produces: `[R, W, dRdTh] = euler_zyx(Theta)`, `Theta = [phi; theta; psi]`; `R = Rz(psi) Ry(theta) Rx(phi)` [3x3] body-to-inertial; `W` [3x3] with `Thetadot = W omega`; `dRdTh` {1x3} cell of `dR/dTheta_k`. Complex-step safe.

- [ ] **Step 1: Write the failing test**

```matlab
function ok = test_euler_zyx()
%% Purpose:
%
%   euler_zyx against things that do not come from euler_zyx:
%     1. R is a rotation: R'R = I, det R = +1, R(0) = I;
%     2. R equals the three elementary rotations written out longhand here;
%     3. dRdTh matches a complex step through euler_zyx's own R (this one
%        IS self-referential in the field, but the complex step is an
%        independent differentiation of it);
%     4. W is the inverse of the body-rate map E(Theta) built here from the
%        kinematics of successive rotations:
%            omega = [phidot; 0; 0] + Rx' [0; thetadot; 0] + Rx' Ry' [0; 0; psidot]
%        so E = [1, 0, -sin th; 0, cos ph, sin ph cos th; 0, -sin ph, cos ph cos th]
%        and W E = I at every attitude away from the singularity;
%     5. W(0) = I.
%
%% Inputs:
%
%  none
%
%% Outputs:
%
%  ok                       logical                 All checks passed
%
%% Revision History:
%  M. Casey                                                   (c) 10/23/2026
%  Copyright Coorbital Inc.
%% ------------------------ Begin Code Sequence ---------------------------

ok = true;
addpath(fileparts(fileparts(mfilename('fullpath'))));
rngWas = rng(23);  restore = onCleanup(@() rng(rngWas));

worstOrth = 0;  worstDet = 0;  worstLong = 0;  worstD = 0;  worstW = 0;
for trial = 1:50
    Th = (rand(3,1) - 0.5)*2;  Th(2) = Th(2)*1.2;        % |theta| < 1.2 rad, away from 90 deg
    [R, W, dR] = euler_zyx(Th);
    worstOrth = max(worstOrth, max(max(abs(R.'*R - eye(3)))));
    worstDet  = max(worstDet, abs(det(R) - 1));
    ph = Th(1); th = Th(2); ps = Th(3);
    Rx = [1 0 0; 0 cos(ph) -sin(ph); 0 sin(ph) cos(ph)];
    Ry = [cos(th) 0 sin(th); 0 1 0; -sin(th) 0 cos(th)];
    Rz = [cos(ps) -sin(ps) 0; sin(ps) cos(ps) 0; 0 0 1];
    worstLong = max(worstLong, max(max(abs(R - Rz*Ry*Rx))));
    hs = 1e-20;
    for k = 1:3
        Tc = Th;  Tc(k) = Tc(k) + 1i*hs;
        Rc = euler_zyx(Tc);
        worstD = max(worstD, max(max(abs(imag(Rc)/hs - dR{k}))));
    end
    E = [1, 0, -sin(th); 0, cos(ph), sin(ph)*cos(th); 0, -sin(ph), cos(ph)*cos(th)];
    worstW = max(worstW, max(max(abs(W*E - eye(3)))));
end
ok = chk(ok, worstOrth < 1e-14, sprintf('R''R = I to %.1e', worstOrth));
ok = chk(ok, worstDet < 1e-14, sprintf('det R = 1 to %.1e', worstDet));
ok = chk(ok, worstLong < 1e-15, sprintf('R = Rz Ry Rx longhand to %.1e', worstLong));
ok = chk(ok, worstD < 1e-13, sprintf('dR/dTheta vs complex step to %.1e', worstD));
ok = chk(ok, worstW < 1e-13, sprintf('W E = I (E from successive-rotation kinematics) to %.1e', worstW));
[R0, W0] = euler_zyx([0; 0; 0]);
ok = chk(ok, isequal(R0, eye(3)) && isequal(W0, eye(3)), 'R(0) = I and W(0) = I exactly');

if ok, fprintf('TEST_EULER_ZYX: ALL PASS\n'); else, fprintf('TEST_EULER_ZYX: FAILURE\n'); end
end

% ------------------------------------------------------------------------
function ok = chk(ok, cond, label)
%% Purpose:
%
%   Accumulate one labelled pass/fail.
%
if cond, fprintf('  PASS  %s\n', label);
else,    fprintf('  FAIL  %s\n', label);  ok = false;
end
end
```

- [ ] **Step 2: Run it to verify it fails**

Expected: `Unrecognized function or variable 'euler_zyx'`.

- [ ] **Step 3: Implement**

```matlab
function [R, W, dRdTh] = euler_zyx(Theta)
%% Purpose:
%
%   ZYX (yaw-pitch-roll) Euler attitude: the body-to-inertial rotation
%   R = Rz(psi) Ry(theta) Rx(phi), the map W from body rates to Euler rates
%   (Thetadot = W omega), and dR/dTheta. Complex-step safe (sin, cos, tan
%   and division only; no abs, no branch on a value).
%
%  ASSUMPTIONS / NOTES:
%
% • W is singular at theta = +/- pi/2 (division by cos theta). The models
%   that use it carry a thetaMax box; this function does not guard, so a
%   caller sees Inf/NaN rather than a silent clip.
%
%% Inputs:
%
%  Theta                    [3 x 1]                 [phi; theta; psi]
%                                                   roll, pitch, yaw (rad)
%
%% Outputs:
%
%  R                        [3 x 3]                 Body-to-inertial
%
%  W                        [3 x 3]                 Thetadot = W omega
%
%  dRdTh                    {1 x 3}                 dR/dphi, dR/dtheta,
%                                                   dR/dpsi
%
%% Revision History:
%  M. Casey                                                   (c) 10/23/2026
%  Copyright Coorbital Inc.
%% ------------------------ Begin Code Sequence ---------------------------

if nargin == 0
   %Demo: a 30-degree roll, 20-degree pitch, 45-degree yaw:
     [R, W] = euler_zyx([30; 20; 45]*pi/180);
     fprintf('R =\n'); disp(R); fprintf('W =\n'); disp(W);
     return
end

ph = Theta(1);  th = Theta(2);  ps = Theta(3);
cph = cos(ph);  sph = sin(ph);
cth = cos(th);  sth = sin(th);
cps = cos(ps);  sps = sin(ps);

Rx = [1 0 0; 0 cph -sph; 0 sph cph];
Ry = [cth 0 sth; 0 1 0; -sth 0 cth];
Rz = [cps -sps 0; sps cps 0; 0 0 1];
R  = Rz*Ry*Rx;

if nargout > 1
    tth = sth/cth;
    W = [1, sph*tth, cph*tth;
         0, cph,     -sph;
         0, sph/cth, cph/cth];
end
if nargout > 2
    dRx = [0 0 0; 0 -sph -cph; 0 cph -sph];
    dRy = [-sth 0 cth; 0 0 0; -cth 0 -sth];
    dRz = [-sps -cps 0; cps -sps 0; 0 0 0];
    dRdTh = {Rz*Ry*dRx, Rz*dRy*Rx, dRz*Ry*Rx};
end
end
```

- [ ] **Step 4: Run the test**

Expected: all `PASS`.

- [ ] **Step 5: Commit**

```bash
/Library/Developer/CommandLineTools/usr/bin/git add optimal_control_examples/drone/drone_dynamics/euler_zyx.m optimal_control_examples/drone/drone_dynamics/tests/test_euler_zyx.m
/Library/Developer/CommandLineTools/usr/bin/git commit -m "drone_dynamics: euler_zyx (R, W, dR/dTheta), tested against longhand rotations, complex step and the successive-rotation rate map

Co-Authored-By: Claude Fable 5.1 <noreply@anthropic.com>"
```

---

### Task 3: `drone_field` and the physics oracle

**Files:**
- Create: `optimal_control_examples/drone/drone_dynamics/drone_field.m`
- Test: `optimal_control_examples/drone/drone_dynamics/tests/test_drone_physics.m` (THE oracle, written first)
- Test: `optimal_control_examples/drone/drone_dynamics/tests/test_drone_field.m`

**Interfaces:**
- Consumes: `drone_params`, `euler_zyx`.
- Produces: `[F, B] = drone_field(x, p)`; `xdot = F + B*u`. Shapes `[p.nx x 1]`, `[p.nx x p.nu]`. State layouts: kinematic `x = p`; pointmass `[p; v]`; rigid6dof `[p; v; Theta; omega]`. Rotor mixer (rotor 1 on +x_b, 2 on +y_b, 3 on -x_b, 4 on -y_b): `T = sum f`, `tau = [l (f2 - f4); l (f3 - f1); kappa (f1 - f2 + f3 - f4)]`.

- [ ] **Step 1: Write the oracle test (fails: no field yet)**

```matlab
function ok = test_drone_physics()
%% Purpose:
%
%   THE INDEPENDENT PHYSICAL ORACLE for the drone models. Every other test
%   in this folder compares the field against itself (a complex step of it,
%   a generated Jacobian of it); this one derives its expectations from
%   mechanics and geometry and never from the implemented rows. The house
%   rule behind it: on 2026-09-17 a cart-pole sign error survived nine
%   self-comparing reviews, and a power-balance identity found it.
%
%   Checks:
%     1. kinematic: F = 0 and B = I exactly (the model is its own oracle).
%     2. pointmass: free fall vdot = -g e3; hover a = g e3 holds v; power
%        balance d/dt(1/2 |v|^2 + g pz) = a . v pointwise.
%     3. rigid6dof hover: f_k = m g / 4, Theta = 0, omega = 0, v = 0 is an
%        equilibrium.
%     4. rigid6dof free fall: f = 0 gives vdot = -g e3 at every attitude.
%     5. rigid6dof POWER BALANCE, pointwise at random states and controls:
%            d/dt (1/2 m |v|^2 + m g pz + 1/2 omega' I omega)
%              = T (R e3) . v + tau . omega
%        The gyroscopic term omega x I omega does no work, so it cancels;
%        a sign error in any row breaks the identity.
%     6. rigid6dof KINEMATIC CONSISTENCY, pointwise: Rdot = R [omega]x with
%        Rdot = sum_k dR/dTheta_k Thetadot_k -- pins W independently of the
%        rates row.
%     7. rigid6dof torque-free flight (equal rotor thrusts) conserves the
%        INERTIAL angular momentum R I omega; a hover-trim flight stays put.
%
%% Inputs:
%
%  none
%
%% Outputs:
%
%  ok                       logical                 All checks passed
%
%% Revision History:
%  M. Casey                                                   (c) 10/23/2026
%  Copyright Coorbital Inc.
%% ------------------------ Begin Code Sequence ---------------------------

ok = true;
addpath(fileparts(fileparts(mfilename('fullpath'))));
rngWas = rng(31);  restore = onCleanup(@() rng(rngWas));
e3 = [0; 0; 1];

%% 1. kinematic
pk = drone_params('kinematic');
[F, B] = drone_field(randn(3,1), pk);
ok = chk(ok, isequal(F, zeros(3,1)) && isequal(B, eye(3)), 'kinematic: F = 0, B = I exactly');

%% 2. pointmass
pq = drone_params('pointmass');
worstFall = 0;  worstHover = 0;  worstPow = 0;
for trial = 1:20
    x = randn(6,1);  a = randn(3,1);
    [F, B] = drone_field(x, pq);
    xd = F + B*a;
    worstFall  = max(worstFall, max(abs(F(4:6) + pq.g*e3)));
    xdH = F + B*(pq.g*e3);
    worstHover = max(worstHover, max(abs(xdH(4:6))));
    v = x(4:6);
    lhs = v.'*xd(4:6) + pq.g*xd(3);                 % d/dt (1/2|v|^2 + g pz) along the field
    rhs = a.'*v;
    worstPow = max(worstPow, abs(lhs - rhs)/max(1, abs(rhs)));
end
ok = chk(ok, worstFall < 1e-15, sprintf('pointmass free fall: vdot = -g e3 (%.1e)', worstFall));
ok = chk(ok, worstHover < 1e-15, sprintf('pointmass hover a = g e3 holds v (%.1e)', worstHover));
ok = chk(ok, worstPow < 1e-13, sprintf('pointmass power balance d/dt(1/2|v|^2 + g pz) = a.v (%.1e)', worstPow));

%% 3. rigid6dof hover equilibrium
pr = drone_params('rigid6dof');
xH = [randn(3,1); zeros(9,1)];
fH = (pr.m*pr.g/4)*ones(4,1);
[F, B] = drone_field(xH, pr);
ok = chk(ok, max(abs(F + B*fH)) < 1e-15, 'rigid6dof hover (f_k = m g/4, level, at rest) is an equilibrium');

%% 4-6. rigid6dof pointwise identities at random states and controls
worstFall = 0;  worstPow = 0;  worstKin = 0;
Imat = diag(pr.I);
for trial = 1:50
    x = randn(12,1);  x(7:8) = x(7:8)*0.5;           % roll/pitch within +/- ~1.2 rad
    f = pr.fMax*rand(4,1);
    [F, B] = drone_field(x, pr);
    xd = F + B*f;
    worstFall = max(worstFall, max(abs(F(4:6) + pr.g*e3)));
    v = x(4:6);  Th = x(7:9);  om = x(10:12);
    [R, ~, dR] = euler_zyx(Th);
    T   = sum(f);
    tau = [pr.l*(f(2) - f(4)); pr.l*(f(3) - f(1)); pr.kappa*(f(1) - f(2) + f(3) - f(4))];
    % 5. power balance: E = 1/2 m|v|^2 + m g pz + 1/2 om' I om, dE/dt along the field
    lhs = pr.m*(v.'*xd(4:6)) + pr.m*pr.g*xd(3) + om.'*(Imat*xd(10:12));
    rhs = T*((R*e3).'*v) + tau.'*om;
    worstPow = max(worstPow, abs(lhs - rhs)/max(1e-3, abs(rhs)));
    % 6. kinematic consistency: Rdot from Thetadot vs R [omega]x
    Thd = xd(7:9);
    Rdot = dR{1}*Thd(1) + dR{2}*Thd(2) + dR{3}*Thd(3);
    omx = [0 -om(3) om(2); om(3) 0 -om(1); -om(2) om(1) 0];
    worstKin = max(worstKin, max(max(abs(Rdot - R*omx))));
end
ok = chk(ok, worstFall < 1e-15, sprintf('rigid6dof free fall: vdot = -g e3 at every attitude (%.1e)', worstFall));
ok = chk(ok, worstPow < 1e-10, sprintf('rigid6dof POWER BALANCE dE/dt = T (R e3).v + tau.omega (%.1e rel)', worstPow));
ok = chk(ok, worstKin < 1e-12, sprintf('rigid6dof KINEMATICS Rdot = R [omega]x (%.1e)', worstKin));

%% 7. torque-free flight conserves R I omega; hover trim stays put
x0 = [zeros(6,1); 0.3; -0.2; 0.5; 1.0; -0.7; 0.4];
fEq = 0.05*ones(4,1);                                % equal thrusts: tau = 0, T = 0.2 N
oo = odeset('RelTol', 1e-12, 'AbsTol', 1e-14);
[~, X] = ode113(@(t, x) fieldAt(x, fEq, pr), [0 2], x0, oo);
L0 = euler_zyx(x0(7:9))*(Imat*x0(10:12));
Lend = euler_zyx(X(end,7:9).')*(Imat*X(end,10:12).');
ok = chk(ok, max(abs(Lend - L0))/max(abs(L0)) < 1e-8, ...
         sprintf('torque-free flight conserves inertial angular momentum R I omega (%.1e rel over 2 s)', ...
                 max(abs(Lend - L0))/max(abs(L0))));
[~, XH] = ode113(@(t, x) fieldAt(x, fH, pr), [0 2], xH, oo);
ok = chk(ok, max(abs(XH(end,:).' - xH)) < 1e-12, 'hover trim holds for 2 s');

if ok, fprintf('TEST_DRONE_PHYSICS: ALL PASS\n'); else, fprintf('TEST_DRONE_PHYSICS: FAILURE\n'); end
end

% ------------------------------------------------------------------------
function xd = fieldAt(x, u, p)
%% Purpose:
%
%   xdot = F(x) + B(x) u for the integrator.
%
[F, B] = drone_field(x, p);
xd = F + B*u;
end

function ok = chk(ok, cond, label)
%% Purpose:
%
%   Accumulate one labelled pass/fail.
%
if cond, fprintf('  PASS  %s\n', label);
else,    fprintf('  FAIL  %s\n', label);  ok = false;
end
end
```

- [ ] **Step 2: Write the shape/consistency test**

```matlab
function ok = test_drone_field()
%% Purpose:
%
%   Shape and structure of drone_field for each model, and complex-step
%   safety: a complex perturbation of one state must produce a finite
%   complex result whose imaginary part is a derivative (compared against
%   central differences of F here; the full Jacobian test is
%   test_drone_state_jac). These are CONSISTENCY checks on the field's
%   form; the physics is pinned by test_drone_physics.
%
%% Inputs:
%
%  none
%
%% Outputs:
%
%  ok                       logical                 All checks passed
%
%% Revision History:
%  M. Casey                                                   (c) 10/23/2026
%  Copyright Coorbital Inc.
%% ------------------------ Begin Code Sequence ---------------------------

ok = true;
addpath(fileparts(fileparts(mfilename('fullpath'))));
rngWas = rng(7);  restore = onCleanup(@() rng(rngWas));

models = {'kinematic', 'pointmass', 'rigid6dof'};
for mdl = 1:3
    p = drone_params(models{mdl});
    x = randn(p.nx, 1);  if p.nx == 12, x(7:8) = 0.3*x(7:8); end
    [F, B] = drone_field(x, p);
    ok = chk(ok, isequal(size(F), [p.nx 1]) && isequal(size(B), [p.nx p.nu]), ...
             sprintf('%s: F is %dx1, B is %dx%d', p.model, p.nx, p.nx, p.nu));
    ok = chk(ok, all(isfinite(F(:))) && all(isfinite(B(:))), sprintf('%s: finite', p.model));
    % complex-step safety: imag part of F(x + i h e_k)/h vs central FD
    worst = 0;  hc = 1e-20;  hd = 1e-6;
    for k = 1:p.nx
        xc = x;  xc(k) = xc(k) + 1i*hc;
        Fc = drone_field(xc, p);
        xp = x;  xp(k) = xp(k) + hd;  xm = x;  xm(k) = xm(k) - hd;
        Fd = (drone_field(xp, p) - drone_field(xm, p))/(2*hd);
        worst = max(worst, max(abs(imag(Fc)/hc - Fd)));
        ok = ok && all(isfinite(Fc(:)));
    end
    ok = chk(ok, worst < 1e-7, sprintf('%s: complex step of F matches central differences (%.1e)', p.model, worst));
end

p = drone_params('rigid6dof');
x = randn(12,1);  x(7:8) = 0.3*x(7:8);
[~, B] = drone_field(x, p);
ok = chk(ok, max(max(abs(B(4:6,:) - B(4:6,1)*ones(1,4)))) < 1e-15, ...
         'rigid6dof: every rotor pushes along the same body z axis (B(4:6,:) columns equal)');
ok = chk(ok, isequal(B(1:3,:), zeros(3,4)) && isequal(B(7:9,:), zeros(3,4)), ...
         'rigid6dof: thrusts do not enter the position or Euler-rate rows');
Mix = [0 p.l 0 -p.l; -p.l 0 p.l 0; p.kappa -p.kappa p.kappa -p.kappa];
ok = chk(ok, max(max(abs(B(10:12,:) - diag(1./p.I)*Mix))) < 1e-12, 'rigid6dof: B(10:12,:) = I^-1 Mix (the "+" mixer)');
xT = x;  xT(8) = p.thetaMax;
F = drone_field(xT, p);
ok = chk(ok, all(isfinite(F)), 'rigid6dof: finite at theta = thetaMax');
xS = x;  xS(8) = pi/2;
F = drone_field(xS, p);
ok = chk(ok, ~all(isfinite(F)), 'rigid6dof: NOT finite at theta = pi/2 (the singularity is not hidden)');

if ok, fprintf('TEST_DRONE_FIELD: ALL PASS\n'); else, fprintf('TEST_DRONE_FIELD: FAILURE\n'); end
end

% ------------------------------------------------------------------------
function ok = chk(ok, cond, label)
%% Purpose:
%
%   Accumulate one labelled pass/fail.
%
if cond, fprintf('  PASS  %s\n', label);
else,    fprintf('  FAIL  %s\n', label);  ok = false;
end
end
```

- [ ] **Step 3: Run both to verify they fail**

Expected: `Unrecognized function or variable 'drone_field'` from each.

- [ ] **Step 4: Implement**

```matlab
function [F, B] = drone_field(x, p)
%% Purpose:
%
%   The drone dynamics in control-affine form, xdot = F(x) + B(x) u, for
%   the three models drone_params knows. One home for the physics: the
%   generated Jacobians are differentiated FROM this function, and the
%   oracle (tests/test_drone_physics) pins it by mechanics, not by
%   agreement with a copy.
%
%   kinematic   x = p (3)                     u = v (3):  pdot = v
%   pointmass   x = [p; v] (6)                u = a (3):  vdot = a - g e3
%   rigid6dof   x = [p; v; Theta; omega] (12) u = f (4):
%                 pdot     = v
%                 vdot     = -g e3 + (T/m) R(Theta) e3
%                 Thetadot = W(Theta) omega
%                 omegadot = I^-1 (tau - omega x I omega)
%               T = sum f,  tau = [l (f2 - f4); l (f3 - f1); kappa (f1 - f2 + f3 - f4)]
%               rotor 1 on +x_b, 2 on +y_b, 3 on -x_b, 4 on -y_b ("+" layout)
%
%  ASSUMPTIONS / NOTES:
%
% • Complex-step safe: sin/cos/tan, products and divisions only; the cross
%   product is written out; no abs, norm, max, min, no branch on a value.
% • W is singular at theta = +/- pi/2; p.thetaMax is the box the solvers
%   impose. Nothing is clipped here.
%
%% Inputs:
%
%  x                        [nx x 1]                State, layout per model
%
%  p                        struct                  drone_params(model)
%
%% Outputs:
%
%  F                        [nx x 1]                Drift
%
%  B                        [nx x nu]               Control columns
%
%% Revision History:
%  M. Casey                                                   (c) 10/23/2026
%  Copyright Coorbital Inc.
%% ------------------------ Begin Code Sequence ---------------------------

if nargin == 0
   %Demo: the rigid body at a 10-degree roll with hover thrust:
     pr = drone_params('rigid6dof');
     x0 = [zeros(6,1); 10*pi/180; 0; 0; zeros(3,1)];
     [F, B] = drone_field(x0, pr);
     xd = F + B*(pr.m*pr.g/4)*ones(4,1);
     fprintf('xdot at 10 deg roll, hover thrust: [%s]\n', sprintf('%.4f ', xd));
     return
end

switch p.model
    case 'kinematic'
        F = zeros(3, 1);
        B = eye(3);
    case 'pointmass'
        v = x(4:6);
        F = [v; 0; 0; -p.g];
        B = [zeros(3); eye(3)];
    case 'rigid6dof'
        v  = x(4:6);
        Th = x(7:9);
        om = x(10:12);
        [R, W] = euler_zyx(Th);
        Iw = p.I.*om;                                  % I omega, I diagonal
        omxIw = [om(2)*Iw(3) - om(3)*Iw(2);
                 om(3)*Iw(1) - om(1)*Iw(3);
                 om(1)*Iw(2) - om(2)*Iw(1)];
        F = [v;
             0; 0; -p.g;
             W*om;
             -omxIw./p.I];
        Mix = [0, p.l, 0, -p.l;
               -p.l, 0, p.l, 0;
               p.kappa, -p.kappa, p.kappa, -p.kappa];
        B = [zeros(3, 4);
             (R(:,3)/p.m)*ones(1, 4);
             zeros(3, 4);
             Mix./p.I];
    otherwise
        error('drone_field:model', 'unknown model ''%s''', p.model);
end
end
```

- [ ] **Step 5: Run both tests**

Expected: `TEST_DRONE_PHYSICS: ALL PASS` and `TEST_DRONE_FIELD: ALL PASS`. If the power balance fails, the sign is in a row of `F` or `B`; do not loosen the bar.

- [ ] **Step 6: Commit**

```bash
/Library/Developer/CommandLineTools/usr/bin/git add optimal_control_examples/drone/drone_dynamics/drone_field.m optimal_control_examples/drone/drone_dynamics/tests/test_drone_physics.m optimal_control_examples/drone/drone_dynamics/tests/test_drone_field.m
/Library/Developer/CommandLineTools/usr/bin/git commit -m "drone_dynamics: drone_field (kinematic / pointmass / rigid6dof, control-affine) with the mechanics oracle (power balance, kinematic consistency, angular-momentum conservation) and the form test

Co-Authored-By: Claude Fable 5.1 <noreply@anthropic.com>"
```

---

### Task 4: Generated state Jacobian

**Files:**
- Create: `optimal_control_examples/drone/drone_dynamics/gen_drone_jacs.m` (generator)
- Create: `optimal_control_examples/drone/drone_dynamics/drone_state_jac_rigid6dof_gen.m` (generated, committed)
- Create: `optimal_control_examples/drone/drone_dynamics/drone_state_jac.m`
- Test: `optimal_control_examples/drone/drone_dynamics/tests/test_drone_state_jac.m`

**Interfaces:**
- Produces: `A = drone_state_jac(x, u, p)` = `d(F + B u)/dx` [nx x nx]. Kinematic and pointmass are constant and written in closed form in the dispatcher (`zeros(3)` and `[0 I; 0 0]`); rigid6dof dispatches to the generated `drone_state_jac_rigid6dof_gen(x, u, pv)` with `pv = [m; I1; I2; I3; l; kappa; g]`.

- [ ] **Step 1: Write the failing test**

```matlab
function ok = test_drone_state_jac()
%% Purpose:
%
%   drone_state_jac against two independent differentiations of
%   drone_field: a complex step (exact to round-off, possible because the
%   field is complex-step safe) and central finite differences. Also: the
%   generated file exists, and regenerating it reproduces it byte for byte
%   (so nobody has hand-edited the generated body).
%
%% Inputs:
%
%  none
%
%% Outputs:
%
%  ok                       logical                 All checks passed
%
%% Revision History:
%  M. Casey                                                   (c) 10/23/2026
%  Copyright Coorbital Inc.
%% ------------------------ Begin Code Sequence ---------------------------

ok = true;
here = fileparts(fileparts(mfilename('fullpath')));
addpath(here);
rngWas = rng(11);  restore = onCleanup(@() rng(rngWas));

models = {'kinematic', 'pointmass', 'rigid6dof'};
for mdl = 1:3
    p = drone_params(models{mdl});
    worstCS = 0;  worstFD = 0;
    for trial = 1:10
        x = randn(p.nx, 1);  if p.nx == 12, x(7:8) = 0.4*x(7:8); end
        u = randn(p.nu, 1);
        A = drone_state_jac(x, u, p);
        Acs = zeros(p.nx);  Afd = zeros(p.nx);
        hc = 1e-20;  hd = 1e-6;
        for k = 1:p.nx
            xc = x;  xc(k) = xc(k) + 1i*hc;
            [Fc, Bc] = drone_field(xc, p);
            Acs(:,k) = imag(Fc + Bc*u)/hc;
            xp = x;  xp(k) = xp(k) + hd;  xm = x;  xm(k) = xm(k) - hd;
            [Fp, Bp] = drone_field(xp, p);  [Fm, Bm] = drone_field(xm, p);
            Afd(:,k) = ((Fp + Bp*u) - (Fm + Bm*u))/(2*hd);
        end
        worstCS = max(worstCS, max(max(abs(A - Acs))));
        worstFD = max(worstFD, max(max(abs(A - Afd))));
    end
    ok = chk(ok, isequal(size(A), [p.nx p.nx]), sprintf('%s: A is %dx%d', p.model, p.nx, p.nx));
    ok = chk(ok, worstCS < 1e-12, sprintf('%s: A vs complex step (%.1e)', p.model, worstCS));
    ok = chk(ok, worstFD < 1e-6, sprintf('%s: A vs central differences (%.1e)', p.model, worstFD));
end

genFile = fullfile(here, 'drone_state_jac_rigid6dof_gen.m');
ok = chk(ok, isfile(genFile), 'the generated rigid6dof Jacobian file exists');
before = fileread(genFile);
gen_drone_jacs();
after = fileread(genFile);
ok = chk(ok, strcmp(before, after), 'regenerating reproduces the committed file byte for byte');

if ok, fprintf('TEST_DRONE_STATE_JAC: ALL PASS\n'); else, fprintf('TEST_DRONE_STATE_JAC: FAILURE\n'); end
end

% ------------------------------------------------------------------------
function ok = chk(ok, cond, label)
%% Purpose:
%
%   Accumulate one labelled pass/fail.
%
if cond, fprintf('  PASS  %s\n', label);
else,    fprintf('  FAIL  %s\n', label);  ok = false;
end
end
```

- [ ] **Step 2: Run it to verify it fails**

Expected: `Unrecognized function or variable 'drone_state_jac'`.

- [ ] **Step 3: Write the generator and run it once**

```matlab
function gen_drone_jacs()
%% Purpose:
%
%   Generate drone_state_jac_rigid6dof_gen.m: the exact 12x12 state
%   Jacobian d(F + B u)/dx of the rigid-body model, by symbolic
%   differentiation of drone_field ITSELF on symbolic inputs, so this file
%   holds no second copy of the physics. Re-run whenever drone_field
%   changes; the generated file is committed and never hand-edited
%   (test_drone_state_jac regenerates and compares byte for byte).
%
%   The kinematic and pointmass Jacobians are constant matrices written in
%   drone_state_jac directly; generating them would add machinery for
%   nothing.
%
%% Inputs:
%
%  none
%
%% Outputs:
%
%  none (writes drone_state_jac_rigid6dof_gen.m beside this file)
%
%% Revision History:
%  M. Casey                                                   (c) 10/23/2026
%  Copyright Coorbital Inc.
%% ------------------------ Begin Code Sequence ---------------------------

here = fileparts(mfilename('fullpath'));
addpath(here);

xs = sym('x', [12 1], 'real');
us = sym('u', [4 1], 'real');
pv = sym('pv', [7 1], 'real');                  % [m; I1; I2; I3; l; kappa; g]
ps = drone_params('rigid6dof');
ps.m = pv(1);  ps.I = pv(2:4);  ps.l = pv(5);  ps.kappa = pv(6);  ps.g = pv(7);

[F, B] = drone_field(xs, ps);
A = jacobian(F + B*us, xs);
outFile = fullfile(here, 'drone_state_jac_rigid6dof_gen.m');
matlabFunction(A, 'File', outFile, 'Vars', {xs, us, pv}, 'Outputs', {'A'}, 'Optimize', true);

% matlabFunction stamps the file with the generation date, which would make
% every regeneration differ; strip that one line so the committed file is
% reproducible byte for byte (test_drone_state_jac checks exactly that), and
% say where the file came from.
txt = fileread(outFile);
lines = regexp(txt, '\r?\n', 'split');
keep = cellfun(@isempty, regexp(lines, '^%\s+\d{2}-\w{3}-\d{4}\s+\d{2}:\d{2}:\d{2}\s*$', 'once'));
lines = lines(keep);
lines = [lines(1), {'%    GENERATED by gen_drone_jacs.m from drone_field.m -- do not edit.'}, lines(2:end)];
fid = fopen(outFile, 'w');
fprintf(fid, '%s\n', lines{:});
fclose(fid);
fprintf('gen_drone_jacs: wrote drone_state_jac_rigid6dof_gen.m\n');
end
```

Run: `matlab -batch "cd('optimal_control_examples/drone/drone_dynamics'); gen_drone_jacs"`
Expected: the generated file appears, its second line is the GENERATED comment, and no date stamp remains (`grep -c '20[0-9][0-9] [0-9][0-9]:' drone_state_jac_rigid6dof_gen.m` prints 0). Run the generator twice and `diff` the results: no output.

- [ ] **Step 4: Write the dispatcher**

```matlab
function A = drone_state_jac(x, u, p)
%% Purpose:
%
%   The exact state Jacobian A = d(F(x) + B(x) u)/dx at fixed u, for each
%   model. kinematic and pointmass are constant and written here;
%   rigid6dof is the generated drone_state_jac_rigid6dof_gen (from
%   gen_drone_jacs, differentiated from drone_field itself).
%
%% Inputs:
%
%  x                        [nx x 1]                State
%
%  u                        [nu x 1]                Control
%
%  p                        struct                  drone_params(model)
%
%% Outputs:
%
%  A                        [nx x nx]               d(F + B u)/dx
%
%% Revision History:
%  M. Casey                                                   (c) 10/23/2026
%  Copyright Coorbital Inc.
%% ------------------------ Begin Code Sequence ---------------------------

if nargin == 0
   %Demo: the rigid body's Jacobian at a hover state:
     pr = drone_params('rigid6dof');
     A = drone_state_jac(zeros(12,1), (pr.m*pr.g/4)*ones(4,1), pr);
     fprintf('A at hover: nnz = %d, eigenvalue magnitudes max %.3e\n', nnz(A), max(abs(eig(A))));
     return
end

switch p.model
    case 'kinematic'
        A = zeros(3);
    case 'pointmass'
        A = [zeros(3), eye(3); zeros(3, 6)];
    case 'rigid6dof'
        A = drone_state_jac_rigid6dof_gen(x(:), u(:), [p.m; p.I(:); p.l; p.kappa; p.g]);
    otherwise
        error('drone_state_jac:model', 'unknown model ''%s''', p.model);
end
end
```

- [ ] **Step 5: Run the test**

Expected: all `PASS`, including the byte-for-byte regeneration line (the generator strips the date stamp, so two runs are identical).

- [ ] **Step 6: Commit**

```bash
/Library/Developer/CommandLineTools/usr/bin/git add optimal_control_examples/drone/drone_dynamics/gen_drone_jacs.m optimal_control_examples/drone/drone_dynamics/drone_state_jac_rigid6dof_gen.m optimal_control_examples/drone/drone_dynamics/drone_state_jac.m optimal_control_examples/drone/drone_dynamics/tests/test_drone_state_jac.m
/Library/Developer/CommandLineTools/usr/bin/git commit -m "drone_dynamics: drone_state_jac (generated rigid6dof Jacobian from drone_field itself; constant kinematic/pointmass), tested against complex step and central differences, regeneration byte-identical

Co-Authored-By: Claude Fable 5.1 <noreply@anthropic.com>"
```

---

### Task 5: `pair_list`, `lower_bound`, `drone_layouts`

**Files:**
- Create: `optimal_control_examples/drone/drone_dynamics/pair_list.m`
- Create: `optimal_control_examples/drone/drone_dynamics/lower_bound.m`
- Create: `optimal_control_examples/drone/drone_dynamics/drone_layouts.m`
- Test: `optimal_control_examples/drone/drone_dynamics/tests/test_layouts.m`

**Interfaces:**
- Produces: `pairs = pair_list(M)` [2 x M(M-1)/2], first row < second row, lexicographic; `JLB = lower_bound(x0, xT, T)`; `[x0, xT, geom] = drone_layouts(which, M, r)` with `which ∈ {'paper','circle','rings','sphere','headon'}`, `x0, xT` [3 x M], `geom = struct('M','r','d','which','eps',1e-5,'l',1e-4,'ball',[])`, `geom.ball = struct('c',[0;0;0],'R',0.1)` for `'sphere'`. `'paper'` selects by M: M ≤ 16 circle radius 0.5; M = 32 rings `[1/3, 0.5]` 16 each; M = 64 rings `[0.33, 0.5]` 32 each; M = 128 rings `[0.35 0.40 0.45 0.50]` 32 each; `r` defaults to 0.02 for M ≤ 32 and 0.01 above. (The inner radii are inferred from the paper's own `J_LB`: 1.156 at M = 32 implies 1/3, 2.297 at M = 64 implies 0.33.) `'headon'`: M = 2, `x_1: (-0.4,0,0) -> (0.4,0,0)`, `x_2` the reverse. `'sphere'`: M points on a Fibonacci sphere of radius 0.4, antipodal targets, `d = 3`, ball obstacle.

- [ ] **Step 1: Write the failing test**

```matlab
function ok = test_layouts()
%% Purpose:
%
%   The layouts reproduce the paper's lower bounds (its Table 1 J_LB
%   column: 0.200, 0.400, 0.800, 1.156, 2.297, 4.704 at T = 10), targets
%   are antipodal, the 2D layouts have z = 0, pair_list is complete and
%   ordered, and malformed requests are refused by identifier.
%
%% Inputs:
%
%  none
%
%% Outputs:
%
%  ok                       logical                 All checks passed
%
%% Revision History:
%  M. Casey                                                   (c) 10/23/2026
%  Copyright Coorbital Inc.
%% ------------------------ Begin Code Sequence ---------------------------

ok = true;
addpath(fileparts(fileparts(mfilename('fullpath'))));
T = 10;

Ms   = [4 8 16 32 64 128];
JLBp = [0.200 0.400 0.800 1.156 2.297 4.704];
rs   = [0.02 0.02 0.02 0.02 0.01 0.01];
for k = 1:numel(Ms)
    [x0, xT, geom] = drone_layouts('paper', Ms(k));
    J = lower_bound(x0, xT, T);
    ok = chk(ok, abs(J - JLBp(k)) < 6e-4, sprintf('M = %3d: J_LB = %.4f (paper %.3f)', Ms(k), J, JLBp(k)));
    ok = chk(ok, max(max(abs(xT + x0))) < 1e-15, sprintf('M = %3d: targets antipodal', Ms(k)));
    ok = chk(ok, all(x0(3,:) == 0) && all(xT(3,:) == 0) && geom.d == 2, sprintf('M = %3d: planar, d = 2', Ms(k)));
    ok = chk(ok, geom.r == rs(k) && geom.M == Ms(k) && geom.eps == 1e-5 && geom.l == 1e-4 && isempty(geom.ball), ...
             sprintf('M = %3d: geom r = %.2f, eps, l, no ball', Ms(k), rs(k)));
end
[x0, ~, geom] = drone_layouts('paper', 8, 0.05);
ok = chk(ok, geom.r == 0.05 && abs(sqrt(sum(x0(:,1).^2)) - 0.5) < 1e-15, 'explicit r is honoured; circle radius 0.5');

[x0, xT, geom] = drone_layouts('headon');
ok = chk(ok, isequal(x0, [-0.4 0.4; 0 0; 0 0]) && isequal(xT, [0.4 -0.4; 0 0; 0 0]) && geom.M == 2 && geom.r == 0.02, ...
         'headon: the two-agent oracle layout');

[x0, xT, geom] = drone_layouts('sphere', 16);
rad = sqrt(sum(x0.^2, 1));
ok = chk(ok, max(abs(rad - 0.4)) < 1e-12 && max(max(abs(xT + x0))) < 1e-15 && geom.d == 3, 'sphere: radius 0.4, antipodal, d = 3');
ok = chk(ok, isstruct(geom.ball) && isequal(geom.ball.c, [0; 0; 0]) && geom.ball.R == 0.1, 'sphere: ball obstacle c = 0, R = 0.1');
dmin = inf;
for a = 1:15
    for b = a+1:16
        dmin = min(dmin, sqrt(sum((x0(:,a) - x0(:,b)).^2)));
    end
end
ok = chk(ok, dmin > 2*geom.r, sprintf('sphere: initial points clear of each other (min %.3f > 2r)', dmin));

P = pair_list(5);
ok = chk(ok, isequal(P, [1 1 1 1 2 2 2 3 3 4; 2 3 4 5 3 4 5 4 5 5]), 'pair_list(5): the 10 pairs, lexicographic');
ok = chk(ok, isequal(size(pair_list(1)), [2 0]), 'pair_list(1) is [2 x 0]');

try
    drone_layouts('rings', 30);
    ok = chk(ok, false, 'rings with M not divisible refused');
catch err
    ok = chk(ok, strcmp(err.identifier, 'drone_layouts:M'), 'rings with M not divisible refused by drone_layouts:M');
end
try
    drone_layouts('hexagon', 6);
    ok = chk(ok, false, 'unknown layout refused');
catch err
    ok = chk(ok, strcmp(err.identifier, 'drone_layouts:which'), 'unknown layout refused by drone_layouts:which');
end

if ok, fprintf('TEST_LAYOUTS: ALL PASS\n'); else, fprintf('TEST_LAYOUTS: FAILURE\n'); end
end

% ------------------------------------------------------------------------
function ok = chk(ok, cond, label)
%% Purpose:
%
%   Accumulate one labelled pass/fail.
%
if cond, fprintf('  PASS  %s\n', label);
else,    fprintf('  FAIL  %s\n', label);  ok = false;
end
end
```

- [ ] **Step 2: Run it to verify it fails**

Expected: `Unrecognized function or variable 'drone_layouts'`.

- [ ] **Step 3: Implement the three functions**

```matlab
function pairs = pair_list(M)
%% Purpose:
%
%   Every unordered pair of M agents, as a 2 x P index array with the
%   smaller index first, lexicographic order. P = M (M - 1) / 2.
%
%% Inputs:
%
%  M                        double                  Number of agents
%
%% Outputs:
%
%  pairs                    [2 x P]                 [a; b] with a < b
%
%% Revision History:
%  M. Casey                                                   (c) 10/23/2026
%  Copyright Coorbital Inc.
%% ------------------------ Begin Code Sequence ---------------------------

if nargin == 0
   %Demo:
     disp(pair_list(4));
     return
end

P = M*(M - 1)/2;
pairs = zeros(2, P);
n = 0;
for a = 1:M-1
    for b = a+1:M
        n = n + 1;
        pairs(:, n) = [a; b];
    end
end
end
```

```matlab
function JLB = lower_bound(x0, xT, T)
%% Purpose:
%
%   The unconstrained lower bound on the multi-agent minimum-energy cost
%   J = integral sum_i 1/2 |v_i|^2 dt over fixed T: straight lines at
%   constant speed, J_LB = sum_i |xT_i - x0_i|^2 / (2 T). Every clearance
%   constraint can only raise the cost above this.
%
%% Inputs:
%
%  x0                       [3 x M]                 Initial positions
%
%  xT                       [3 x M]                 Final positions
%
%  T                        double                  Horizon
%
%% Outputs:
%
%  JLB                      double                  The bound
%
%% Revision History:
%  M. Casey                                                   (c) 10/23/2026
%  Copyright Coorbital Inc.
%% ------------------------ Begin Code Sequence ---------------------------

if nargin == 0
   %Demo: the paper's M = 4 circle, J_LB = 0.200:
     [a, b] = drone_layouts('paper', 4);
     fprintf('J_LB(M = 4) = %.4f\n', lower_bound(a, b, 10));
     return
end

JLB = sum(sum((xT - x0).^2))/(2*T);
end
```

```matlab
function [x0, xT, geom] = drone_layouts(which, M, r)
%% Purpose:
%
%   Initial and final positions for the multi-agent problems, as 3 x M
%   arrays (2D layouts carry z = 0), with the geometry struct the clearance
%   and barrier functions read.
%
%   'paper'   Zhang, Liu, Darbon & Karniadakis (arXiv:2610.30388) Table 1
%             by M: M <= 16 one circle of radius 0.5; M = 32 two rings of
%             16 at radii [1/3, 0.5]; M = 64 two rings of 32 at
%             [0.33, 0.5]; M = 128 four rings of 32 at [0.35 0.40 0.45
%             0.50]. Targets antipodal. r = 0.02 for M <= 32, else 0.01.
%             The inner radii are INFERRED from the paper's own J_LB
%             column (1.156 at M = 32 needs 1/3; 2.297 at M = 64 needs
%             0.33); its text says "0.33" for both.
%   'circle'  one circle of radius 0.5, any M
%   'rings'   M/16 rings ... see the switch (M must divide evenly)
%   'sphere'  OUR 3D layout: M points on a Fibonacci sphere of radius 0.4,
%             antipodal targets, ball obstacle at the origin, R = 0.1
%   'headon'  the two-agent closed-form oracle: (-0.4,0,0) <-> (0.4,0,0)
%
%% Inputs:
%
%  which                    char                    Layout name
%
%  M                        double                  Agents (ignored for
%                                                   'headon')
%
%  r                        double (optional)       Agent radius; default
%                                                   per layout
%
%% Outputs:
%
%  x0                       [3 x M]                 Initial positions
%
%  xT                       [3 x M]                 Final positions
%
%  geom                     struct                  .M .r .d (2|3) .which
%                                                   .eps [1e-5] .l [1e-4]
%                                                   .ball ([] or .c .R)
%
%% Revision History:
%  M. Casey                                                   (c) 10/23/2026
%  Copyright Coorbital Inc.
%% ------------------------ Begin Code Sequence ---------------------------

if nargin == 0
   %Demo: the M = 8 circle:
     [a, b] = drone_layouts('paper', 8);
     figure('color', [1 1 1]);  plot(a(1,:), a(2,:), 'o', b(1,:), b(2,:), 'x');  axis equal;  grid on
     return
end
if nargin < 3, r = []; end
if nargin < 2, M = 2; end

d = 2;  ball = [];
switch lower(which)
    case 'paper'
        if M <= 16,      radii = 0.5;                 per = M;
        elseif M == 32,  radii = [1/3, 0.5];          per = 16;
        elseif M == 64,  radii = [0.33, 0.5];         per = 32;
        elseif M == 128, radii = [0.35 0.40 0.45 0.50]; per = 32;
        else, error('drone_layouts:M', 'the paper has no layout for M = %d', M);
        end
        if isempty(r), if M <= 32, r = 0.02; else, r = 0.01; end, end
        x0 = rings(radii, per);
    case 'circle'
        if isempty(r), r = 0.02; end
        x0 = rings(0.5, M);
    case 'rings'
        nRing = 2;
        if mod(M, nRing) ~= 0
            error('drone_layouts:M', 'M = %d is not divisible by %d rings', M, nRing);
        end
        if isempty(r), r = 0.02; end
        x0 = rings([1/3, 0.5], M/nRing);
    case 'sphere'
        if isempty(r), r = 0.02; end
        d = 3;
        ball = struct('c', [0; 0; 0], 'R', 0.1);
        x0 = fibonacci_sphere(M, 0.4);
    case 'headon'
        M = 2;  if isempty(r), r = 0.02; end
        x0 = [-0.4 0.4; 0 0; 0 0];
    otherwise
        error('drone_layouts:which', 'unknown layout ''%s''', which);
end
xT = -x0;
geom = struct('M', M, 'r', r, 'd', d, 'which', lower(which), 'eps', 1e-5, 'l', 1e-4, 'ball', ball);
end

% ------------------------------------------------------------------------
function X = rings(radii, per)
%% Purpose:
%
%   per points on each ring, equally spaced from angle 0, z = 0.
%
X = zeros(3, per*numel(radii));
ang = 2*pi*(0:per-1)/per;
for k = 1:numel(radii)
    cols = (k-1)*per + (1:per);
    X(:, cols) = [radii(k)*cos(ang); radii(k)*sin(ang); zeros(1, per)];
end
end

function X = fibonacci_sphere(M, radius)
%% Purpose:
%
%   M approximately equally spaced points on a sphere (golden-angle
%   spiral).
%
X = zeros(3, M);
gold = pi*(3 - sqrt(5));
for k = 1:M
    z = 1 - 2*(k - 0.5)/M;
    rho = sqrt(1 - z^2);
    X(:, k) = radius*[rho*cos(gold*(k-1)); rho*sin(gold*(k-1)); z];
end
end
```

- [ ] **Step 4: Run the test**

Expected: all `PASS`. The J_LB lines must land within 6e-4 of the paper's rounded values (1.1556 vs 1.156, 2.2970 vs 2.297).

- [ ] **Step 5: Commit**

```bash
/Library/Developer/CommandLineTools/usr/bin/git add optimal_control_examples/drone/drone_dynamics/pair_list.m optimal_control_examples/drone/drone_dynamics/lower_bound.m optimal_control_examples/drone/drone_dynamics/drone_layouts.m optimal_control_examples/drone/drone_dynamics/tests/test_layouts.m
/Library/Developer/CommandLineTools/usr/bin/git commit -m "drone_dynamics: the paper's layouts (inner ring radii inferred from its J_LB column), our sphere + head-on layouts, pair_list, lower_bound, with the Table 1 J_LB test

Co-Authored-By: Claude Fable 5.1 <noreply@anthropic.com>"
```

---

### Task 6: `clearance_terms` and `hermite_resample`

**Files:**
- Create: `optimal_control_examples/drone/drone_dynamics/clearance_terms.m`
- Create: `optimal_control_examples/drone/drone_dynamics/hermite_resample.m`
- Test: `optimal_control_examples/drone/drone_dynamics/tests/test_clearance_terms.m`

**Interfaces:**
- Produces: `[h, rho, n, hBall] = clearance_terms(X, pairs, geom)` for `X` [3M x N] (agent-major stacked columns, one configuration per column): `h = rho - 2r` [P x N], `rho` [P x N], `n` [3 x P x N] unit vectors `x_a - x_b` over `rho`, `hBall` [M x N] (`|x_i - c| - R`, or `[]` without a ball). Complex-step safe. `[Xf, Vf, tf] = hermite_resample(X, V, tN, nSub)`: cubic Hermite from node values and derivatives, `nSub` sub-points per interval, `Xf` [nx x N nSub + 1].

- [ ] **Step 1: Write the failing test**

```matlab
function ok = test_clearance_terms()
%% Purpose:
%
%   clearance_terms against a plain loop written here; M = 1 and N > 1
%   shapes; complex-step safety (the imaginary part of h under a complex
%   perturbation equals the analytic directional derivative n . dx);
%   the ball term. hermite_resample reproduces a cubic exactly and returns
%   the node values at the node positions.
%
%% Inputs:
%
%  none
%
%% Outputs:
%
%  ok                       logical                 All checks passed
%
%% Revision History:
%  M. Casey                                                   (c) 10/23/2026
%  Copyright Coorbital Inc.
%% ------------------------ Begin Code Sequence ---------------------------

ok = true;
addpath(fileparts(fileparts(mfilename('fullpath'))));
rngWas = rng(5);  restore = onCleanup(@() rng(rngWas));

[~, ~, geom] = drone_layouts('sphere', 6);
M = 6;  N = 4;  pairs = pair_list(M);
X = randn(3*M, N);
[h, rho, n, hB] = clearance_terms(X, pairs, geom);
ok = chk(ok, isequal(size(h), [15 N]) && isequal(size(rho), [15 N]) && isequal(size(n), [3 15 N]) && isequal(size(hB), [M N]), ...
         'shapes: h, rho [P x N], n [3 x P x N], hBall [M x N]');
worst = 0;  worstB = 0;
for c = 1:N
    P = reshape(X(:,c), 3, M);
    for q = 1:15
        a = pairs(1,q);  b = pairs(2,q);
        dd = sqrt(sum((P(:,a) - P(:,b)).^2));
        worst = max(worst, abs(h(q,c) - (dd - 2*geom.r)) + abs(rho(q,c) - dd) + max(abs(n(:,q,c) - (P(:,a) - P(:,b))/dd)));
    end
    for a = 1:M
        worstB = max(worstB, abs(hB(a,c) - (sqrt(sum((P(:,a) - geom.ball.c).^2)) - geom.ball.R)));
    end
end
ok = chk(ok, worst < 1e-14, sprintf('pair terms vs a plain loop (%.1e)', worst));
ok = chk(ok, worstB < 1e-14, sprintf('ball terms vs a plain loop (%.1e)', worstB));

% complex step: d h_q / d x along a random direction equals n . (dx_a - dx_b)
dx = randn(3*M, 1);  hc = 1e-20;
hcs = clearance_terms(X(:,1) + 1i*hc*dx, pairs, geom);
worst = 0;
for q = 1:15
    a = pairs(1,q);  b = pairs(2,q);
    dxa = dx(3*(a-1)+(1:3));  dxb = dx(3*(b-1)+(1:3));
    worst = max(worst, abs(imag(hcs(q))/hc - n(:,q,1).'*(dxa - dxb)));
end
ok = chk(ok, worst < 1e-13, sprintf('complex step of h matches n . dx (%.1e)', worst));

[~, ~, g1] = drone_layouts('circle', 1);
[h1, ~, ~, hB1] = clearance_terms(randn(3, 3), pair_list(1), g1);
ok = chk(ok, isequal(size(h1), [0 3]) && isempty(hB1), 'M = 1: h is [0 x N], no ball term');

% hermite_resample: a cubic is reproduced exactly
tN = linspace(0, 2, 6);
c3 = randn(2, 4);
poly  = @(t) c3(:,1) + c3(:,2)*t + c3(:,3)*t.^2 + c3(:,4)*t.^3;
dpoly = @(t) c3(:,2) + 2*c3(:,3)*t + 3*c3(:,4)*t.^2;
[Xf, Vf, tf] = hermite_resample(poly(tN), dpoly(tN), tN, 10);
ok = chk(ok, numel(tf) == 51 && isequal(size(Xf), [2 51]) && isequal(size(Vf), [2 51]), 'hermite_resample: N nSub + 1 points');
ok = chk(ok, max(max(abs(Xf - poly(tf)))) < 1e-12 && max(max(abs(Vf - dpoly(tf)))) < 1e-11, 'hermite_resample reproduces a cubic and its derivative');
ok = chk(ok, max(abs(tf(1:10:end) - tN)) < 1e-15 && max(max(abs(Xf(:,1:10:end) - poly(tN)))) < 1e-15, 'nodes are returned exactly');

if ok, fprintf('TEST_CLEARANCE_TERMS: ALL PASS\n'); else, fprintf('TEST_CLEARANCE_TERMS: FAILURE\n'); end
end

% ------------------------------------------------------------------------
function ok = chk(ok, cond, label)
%% Purpose:
%
%   Accumulate one labelled pass/fail.
%
if cond, fprintf('  PASS  %s\n', label);
else,    fprintf('  FAIL  %s\n', label);  ok = false;
end
end
```

- [ ] **Step 2: Run it to verify it fails**

Expected: `Unrecognized function or variable 'clearance_terms'`.

- [ ] **Step 3: Implement**

```matlab
function [h, rho, n, hBall] = clearance_terms(X, pairs, geom)
%% Purpose:
%
%   Pairwise clearance of a set of agent configurations: for every pair
%   (a, b) and every column of X, the gap h = |x_a - x_b| - 2 r, the
%   distance rho, and the unit vector n = (x_a - x_b)/rho. With geom.ball
%   set, also the obstacle gap |x_i - c| - R per agent. Complex-step safe:
%   distances are sqrt of sums of squares, never a magnitude.
%
%  ASSUMPTIONS / NOTES:
%
% • rho = 0 (coincident agents) gives n = NaN; the barrier and the solver
%   never evaluate there (the NLP uses the squared form, the barrier's
%   quadratic branch is only reached at h <= l > -2r).
%
%% Inputs:
%
%  X                        [3M x N]                Agent-major stacked
%                                                   positions, one
%                                                   configuration per
%                                                   column
%
%  pairs                    [2 x P]                 From pair_list
%
%  geom                     struct                  .r, .ball ([] or .c .R)
%
%% Outputs:
%
%  h                        [P x N]                 rho - 2r
%
%  rho                      [P x N]                 Pair distances
%
%  n                        [3 x P x N]             Unit vectors a -> b
%
%  hBall                    [M x N] or []           |x_i - c| - R
%
%% Revision History:
%  M. Casey                                                   (c) 10/23/2026
%  Copyright Coorbital Inc.
%% ------------------------ Begin Code Sequence ---------------------------

if nargin == 0
   %Demo: the M = 4 circle's initial clearances:
     [a, ~, g] = drone_layouts('paper', 4);
     hh = clearance_terms(a(:), pair_list(4), g);
     fprintf('h at t = 0: [%s]\n', sprintf('%.3f ', hh));
     return
end

M = size(X, 1)/3;
N = size(X, 2);
P = size(pairs, 2);
h   = zeros(P, N, 'like', X);
rho = zeros(P, N, 'like', X);
n   = zeros(3, P, N, 'like', X);
for c = 1:N
    Pc = reshape(X(:, c), 3, M);
    for q = 1:P
        dv = Pc(:, pairs(1, q)) - Pc(:, pairs(2, q));
        dd = sqrt(sum(dv.^2));
        rho(q, c) = dd;
        h(q, c) = dd - 2*geom.r;
        n(:, q, c) = dv/dd;
    end
end
hBall = [];
if isfield(geom, 'ball') && ~isempty(geom.ball)
    hBall = zeros(M, N, 'like', X);
    for c = 1:N
        Pc = reshape(X(:, c), 3, M);
        for a = 1:M
            dv = Pc(:, a) - geom.ball.c;
            hBall(a, c) = sqrt(sum(dv.^2)) - geom.ball.R;
        end
    end
end
end
```

```matlab
function [Xf, Vf, tf] = hermite_resample(X, V, tN, nSub)
%% Purpose:
%
%   Cubic-Hermite reconstruction of a trajectory from node values X and
%   node derivatives V (the trapezoid transcription's own interpolant: v
%   linear on each interval makes x cubic), sampled nSub times per
%   interval. Used by the clearance gate to look BETWEEN the nodes, where
%   a nodal constraint says nothing.
%
%% Inputs:
%
%  X                        [nx x N+1]              Node values
%
%  V                        [nx x N+1]              Node derivatives
%
%  tN                       [1 x N+1]               Node times
%
%  nSub                     double                  Sub-intervals per
%                                                   interval
%
%% Outputs:
%
%  Xf                       [nx x N nSub + 1]       Values on the fine grid
%
%  Vf                       [nx x N nSub + 1]       Derivatives there
%
%  tf                       [1 x N nSub + 1]        The fine grid
%
%% Revision History:
%  M. Casey                                                   (c) 10/23/2026
%  Copyright Coorbital Inc.
%% ------------------------ Begin Code Sequence ---------------------------

if nargin == 0
   %Demo: a sine from 5 nodes:
     t = linspace(0, 2*pi, 6);
     [xf, ~, tt] = hermite_resample(sin(t), cos(t), t, 20);
     fprintf('max error vs sin on the fine grid: %.2e\n', max(abs(xf - sin(tt))));
     return
end

N  = numel(tN) - 1;
nx = size(X, 1);
nF = N*nSub + 1;
Xf = zeros(nx, nF);  Vf = zeros(nx, nF);  tf = zeros(1, nF);
s  = (0:nSub-1)/nSub;
for k = 1:N
    hk = tN(k+1) - tN(k);
    h00 = 2*s.^3 - 3*s.^2 + 1;   h10 = s.^3 - 2*s.^2 + s;
    h01 = -2*s.^3 + 3*s.^2;      h11 = s.^3 - s.^2;
    d00 = 6*s.^2 - 6*s;          d10 = 3*s.^2 - 4*s + 1;
    d01 = -6*s.^2 + 6*s;         d11 = 3*s.^2 - 2*s;
    cols = (k-1)*nSub + (1:nSub);
    Xf(:, cols) = X(:,k)*h00 + hk*V(:,k)*h10 + X(:,k+1)*h01 + hk*V(:,k+1)*h11;
    Vf(:, cols) = (X(:,k)*d00 + hk*V(:,k)*d10 + X(:,k+1)*d01 + hk*V(:,k+1)*d11)/hk;
    tf(cols) = tN(k) + hk*s;
end
Xf(:, end) = X(:, end);  Vf(:, end) = V(:, end);  tf(end) = tN(end);
end
```

- [ ] **Step 4: Run the test**

Expected: all `PASS`.

- [ ] **Step 5: Commit**

```bash
/Library/Developer/CommandLineTools/usr/bin/git add optimal_control_examples/drone/drone_dynamics/clearance_terms.m optimal_control_examples/drone/drone_dynamics/hermite_resample.m optimal_control_examples/drone/drone_dynamics/tests/test_clearance_terms.m
/Library/Developer/CommandLineTools/usr/bin/git commit -m "drone_dynamics: clearance_terms (pair + ball gaps, unit vectors, complex-step safe) and hermite_resample (the trapezoid interpolant on a fine grid), tested against loops and a cubic

Co-Authored-By: Claude Fable 5.1 <noreply@anthropic.com>"
```

---

### Task 7: `barrier_terms` with analytic gradient and Hessian

**Files:**
- Create: `optimal_control_examples/drone/drone_dynamics/barrier_shape.m`
- Create: `optimal_control_examples/drone/drone_dynamics/barrier_terms.m`
- Test: `optimal_control_examples/drone/drone_dynamics/tests/test_barrier_terms.m`

**Interfaces:**
- Produces: `[U, U1, U2] = barrier_shape(h, epsB, l)` elementwise, the paper's C¹ log barrier and its first two derivatives (branch on `real(h) > l`). `[U, gU, HU] = barrier_terms(x, pairs, geom)` for ONE configuration `x` [3M x 1]: `U` scalar sum over pairs (and ball agents), `gU` [3M x 1], `HU` sparse [3M x 3M]. Route 1's CasADi barrier expression is a second copy of `barrier_shape`'s two formulas (CasADi needs `if_else`); its test pins that copy to `barrier_terms` at random nodes.

- [ ] **Step 1: Write the failing test**

```matlab
function ok = test_barrier_terms()
%% Purpose:
%
%   The paper's C1 log barrier and its assembled derivatives:
%     1. barrier_shape is C1 at h = l: U and U' continuous, U'' continuous
%        too (both branches give eps/l^2 there), finite at h = 0 and h < 0;
%     2. U' and U'' vs complex step and central differences of U;
%     3. barrier_terms' gradient vs a complex step of U (complex-step safe);
%     4. barrier_terms' Hessian vs a complex step of the gradient, and
%        symmetric;
%     5. the ball term enters the i-block only;
%     6. M = 1: U = 0, zero gradient, empty-sparse Hessian;
%     7. U equals a plain sum over barrier_shape(h) with h from
%        clearance_terms (the assembly agrees with the pieces).
%
%% Inputs:
%
%  none
%
%% Outputs:
%
%  ok                       logical                 All checks passed
%
%% Revision History:
%  M. Casey                                                   (c) 10/23/2026
%  Copyright Coorbital Inc.
%% ------------------------ Begin Code Sequence ---------------------------

ok = true;
addpath(fileparts(fileparts(mfilename('fullpath'))));
rngWas = rng(3);  restore = onCleanup(@() rng(rngWas));
epsB = 1e-5;  l = 1e-4;

%% 1. the shape
[Ua, U1a, U2a] = barrier_shape(l*(1 + 1e-12), epsB, l);
[Ub, U1b, U2b] = barrier_shape(l*(1 - 1e-12), epsB, l);
ok = chk(ok, abs(Ua - Ub) < 1e-12*epsB && abs(U1a - U1b) < 1e-9*epsB/l && abs(U2a - U2b) < 1e-9*epsB/l^2, ...
         'barrier_shape is C1 at h = l (U, U'' and U'''' all continuous)');
[U0, U10, U20] = barrier_shape([0, -0.01], epsB, l);
ok = chk(ok, all(isfinite([U0 U10 U20])), 'finite at h = 0 and h < 0 (quadratic branch)');
ok = chk(ok, abs(barrier_shape(1, epsB, l)) < 1e-20 && abs(barrier_shape(0.5, epsB, l) + epsB*log(0.5)) < 1e-20, ...
         'U = -eps log h on the log branch');

%% 2. shape derivatives
hs = [2e-4, 1e-3, 0.05, 0.4, 5e-5, -2e-3];
hc = 1e-20;  hd = 1e-7;
worst1 = 0;  worst2 = 0;
for k = 1:numel(hs)
    [~, U1, U2] = barrier_shape(hs(k), epsB, l);
    Uc = barrier_shape(hs(k) + 1i*hc, epsB, l);
    worst1 = max(worst1, abs(imag(Uc)/hc - U1)/max(1, abs(U1)));
    [~, U1p] = barrier_shape(hs(k) + hd, epsB, l);
    [~, U1m] = barrier_shape(hs(k) - hd, epsB, l);
    worst2 = max(worst2, abs((U1p - U1m)/(2*hd) - U2)/max(1, abs(U2)));
end
ok = chk(ok, worst1 < 1e-13, sprintf('U'' vs complex step (%.1e rel)', worst1));
ok = chk(ok, worst2 < 1e-6, sprintf('U'''' vs central differences of U'' (%.1e rel)', worst2));

%% 3-4. assembled gradient and Hessian, with a ball
[x0, ~, geom] = drone_layouts('sphere', 5);
M = 5;  pairs = pair_list(M);
x = x0(:) + 0.05*randn(3*M, 1);
[U, gU, HU] = barrier_terms(x, pairs, geom);
ok = chk(ok, isscalar(U) && isequal(size(gU), [3*M 1]) && isequal(size(HU), [3*M 3*M]) && issparse(HU), 'shapes; HU sparse');
worstG = 0;  worstH = 0;
for k = 1:3*M
    xc = x;  xc(k) = xc(k) + 1i*hc;
    [Uc, gc] = barrier_terms(xc, pairs, geom);
    worstG = max(worstG, abs(imag(Uc)/hc - gU(k))/max(1, abs(gU(k))));
    worstH = max(worstH, max(abs(imag(gc)/hc - HU(:,k)))/max(1, max(abs(HU(:,k)))));
end
ok = chk(ok, worstG < 1e-12, sprintf('gradient vs complex step of U (%.1e rel)', worstG));
ok = chk(ok, worstH < 1e-10, sprintf('Hessian vs complex step of the gradient (%.1e rel)', worstH));
ok = chk(ok, max(max(abs(HU - HU.'))) < 1e-14*max(max(abs(HU))), 'Hessian symmetric');

%% 5. the ball enters only the diagonal agent block
geomNB = geom;  geomNB.ball = [];
[~, gNB, HNB] = barrier_terms(x, pairs, geomNB);
D = HU - HNB;
offBlock = 0;
for a = 1:M
    for b = 1:M
        if a ~= b, offBlock = max(offBlock, max(max(abs(D(3*(a-1)+(1:3), 3*(b-1)+(1:3)))))); end
    end
end
ok = chk(ok, offBlock == 0 && max(abs(gU - gNB)) > 0, 'ball term: i-block only, and it is not zero');

%% 6. M = 1
[~, ~, g1] = drone_layouts('circle', 1);
[U1, g1v, H1] = barrier_terms(randn(3,1), pair_list(1), g1);
ok = chk(ok, U1 == 0 && isequal(g1v, zeros(3,1)) && issparse(H1) && nnz(H1) == 0, 'M = 1: U = 0, zero gradient, empty Hessian');

%% 7. assembly agrees with the pieces
[h, ~, ~, hB] = clearance_terms(x, pairs, geom);
Uplain = sum(barrier_shape(h, geom.eps, geom.l)) + sum(barrier_shape(hB, geom.eps, geom.l));
ok = chk(ok, abs(U - Uplain) < 1e-15*max(1, abs(U)), 'U = plain sum of barrier_shape over clearance_terms');

if ok, fprintf('TEST_BARRIER_TERMS: ALL PASS\n'); else, fprintf('TEST_BARRIER_TERMS: FAILURE\n'); end
end

% ------------------------------------------------------------------------
function ok = chk(ok, cond, label)
%% Purpose:
%
%   Accumulate one labelled pass/fail.
%
if cond, fprintf('  PASS  %s\n', label);
else,    fprintf('  FAIL  %s\n', label);  ok = false;
end
end
```

- [ ] **Step 2: Run it to verify it fails**

Expected: `Unrecognized function or variable 'barrier_shape'`.

- [ ] **Step 3: Implement**

```matlab
function [U, U1, U2] = barrier_shape(h, epsB, l)
%% Purpose:
%
%   The C1 log barrier of Zhang, Liu, Darbon & Karniadakis (arXiv:
%   2610.30388, eq. 3b), elementwise in h, with its first two derivatives:
%
%       U  = -eps log h                                    h > l
%          = -eps log l + (eps/2) (((h - 2l)/l)^2 - 1)     h <= l
%       U' = -eps/h                | (eps/l^2)(h - 2l)
%       U''=  eps/h^2              |  eps/l^2
%
%   Both branches meet at h = l in U, U' AND U''. The quadratic branch keeps
%   everything finite through h = 0 (interpenetration), so an optimizer or
%   an integrator never evaluates log of a non-positive gap.
%
%  ASSUMPTIONS / NOTES:
%
% • Complex-step safe: the branch is on real(h). Route 1's CasADi barrier
%   expression is a SECOND COPY of these two formulas (CasADi needs
%   if_else); its test pins that copy to barrier_terms.
%
%% Inputs:
%
%  h                        [any]                   Gaps (may be complex)
%
%  epsB                     double                  Barrier weight eps
%
%  l                        double                  Switch gap
%
%% Outputs:
%
%  U, U1, U2                [size(h)]               U, dU/dh, d2U/dh2
%
%% Revision History:
%  M. Casey                                                   (c) 10/23/2026
%  Copyright Coorbital Inc.
%% ------------------------ Begin Code Sequence ---------------------------

if nargin == 0
   %Demo: the shape across the switch:
     hh = linspace(-1e-4, 1e-3, 400);
     figure('color', [1 1 1]);  plot(hh, barrier_shape(hh, 1e-5, 1e-4));  grid on
     xlabel('h');  ylabel('U_{\epsilon,l}(h)');
     return
end

logBranch = real(h) > l;
U  = zeros(size(h), 'like', h);
U1 = zeros(size(h), 'like', h);
U2 = zeros(size(h), 'like', h);
hl = h(logBranch);
U(logBranch)  = -epsB*log(hl);
U1(logBranch) = -epsB./hl;
U2(logBranch) =  epsB./hl.^2;
hq = h(~logBranch);
U(~logBranch)  = -epsB*log(l) + (epsB/2)*(((hq - 2*l)/l).^2 - 1);
U1(~logBranch) = (epsB/l^2)*(hq - 2*l);
U2(~logBranch) = (epsB/l^2)*ones(size(hq), 'like', hq);
end
```

```matlab
function [U, gU, HU] = barrier_terms(x, pairs, geom)
%% Purpose:
%
%   The summed barrier of one agent configuration, with its gradient and
%   sparse Hessian: U(x) = sum_pairs U(h_ab) + sum_i U(h_ball,i). For one
%   pair with rho = |x_a - x_b|, n = (x_a - x_b)/rho, h = rho - 2r:
%
%       grad_a =  U'(h) n,   grad_b = -U'(h) n
%       H_aa   =  U''(h) n n' + U'(h) (I - n n')/rho = H_bb,   H_ab = H_ba = -H_aa
%
%   and the ball term has the a-block only (c fixed). This is the
%   objective-independent penalty SHAPE on geometry; the objective that
%   integrates it, and the PMP field A = [0 -I; -HU 0] that reads it, live
%   in the routes. Complex-step safe (branch on real(h) inside
%   barrier_shape; no magnitudes).
%
%% Inputs:
%
%  x                        [3M x 1]                One configuration,
%                                                   agent-major
%
%  pairs                    [2 x P]                 From pair_list
%
%  geom                     struct                  .r .eps .l .ball
%
%% Outputs:
%
%  U                        double                  The summed barrier
%
%  gU                       [3M x 1]                dU/dx
%
%  HU                       [3M x 3M] sparse        d2U/dx2
%
%% Revision History:
%  M. Casey                                                   (c) 10/23/2026
%  Copyright Coorbital Inc.
%% ------------------------ Begin Code Sequence ---------------------------

if nargin == 0
   %Demo: the M = 4 circle at t = 0 (gaps 0.67 and 0.96):
     [a, ~, g] = drone_layouts('paper', 4);
     [U, gU] = barrier_terms(a(:), pair_list(4), g);
     fprintf('U = %.3e, |gU| = %.3e\n', U, sqrt(sum(gU.^2)));
     return
end

M = numel(x)/3;
P = size(pairs, 2);
[h, rho, n, hB] = clearance_terms(x(:), pairs, geom);
[Uq, U1, U2] = barrier_shape(h, geom.eps, geom.l);
U  = sum(Uq);
gU = zeros(3*M, 1, 'like', x);
needH = nargout > 2;
if needH
    rows = zeros(36*P + 9*M, 1);  cols = rows;  vals = zeros(36*P + 9*M, 1, 'like', x);
    nz = 0;
end
E3 = eye(3);
for q = 1:P
    a = pairs(1, q);  b = pairs(2, q);
    ia = 3*(a-1) + (1:3);  ib = 3*(b-1) + (1:3);
    nq = n(:, q);
    gq = U1(q)*nq;
    gU(ia) = gU(ia) + gq;
    gU(ib) = gU(ib) - gq;
    if needH
        Hq = U2(q)*(nq*nq.') + (U1(q)/rho(q))*(E3 - nq*nq.');
        [rows, cols, vals, nz] = addBlock(rows, cols, vals, nz, ia, ia,  Hq);
        [rows, cols, vals, nz] = addBlock(rows, cols, vals, nz, ib, ib,  Hq);
        [rows, cols, vals, nz] = addBlock(rows, cols, vals, nz, ia, ib, -Hq);
        [rows, cols, vals, nz] = addBlock(rows, cols, vals, nz, ib, ia, -Hq);
    end
end
if ~isempty(hB)
    [Ub, U1b, U2b] = barrier_shape(hB, geom.eps, geom.l);
    U = U + sum(Ub);
    for a = 1:M
        ia = 3*(a-1) + (1:3);
        dv = x(ia) - geom.ball.c;
        rb = hB(a) + geom.ball.R;
        nb = dv/rb;
        gU(ia) = gU(ia) + U1b(a)*nb;
        if needH
            Hb = U2b(a)*(nb*nb.') + (U1b(a)/rb)*(E3 - nb*nb.');
            [rows, cols, vals, nz] = addBlock(rows, cols, vals, nz, ia, ia, Hb);
        end
    end
end
if needH
    HU = sparse(rows(1:nz), cols(1:nz), vals(1:nz), 3*M, 3*M);
end
end

% ------------------------------------------------------------------------
function [rows, cols, vals, nz] = addBlock(rows, cols, vals, nz, ia, ib, Hq)
%% Purpose:
%
%   Append one 3x3 block to the sparse triplet list.
%
[C, R] = meshgrid(ib, ia);
idx = nz + (1:9);
rows(idx) = R(:);  cols(idx) = C(:);  vals(idx) = Hq(:);
nz = nz + 9;
end
```

- [ ] **Step 4: Run the test**

Expected: all `PASS`. The Hessian-vs-complex-step line is the one that catches a wrong `(I - n n')/rho` term.

- [ ] **Step 5: Commit**

```bash
/Library/Developer/CommandLineTools/usr/bin/git add optimal_control_examples/drone/drone_dynamics/barrier_shape.m optimal_control_examples/drone/drone_dynamics/barrier_terms.m optimal_control_examples/drone/drone_dynamics/tests/test_barrier_terms.m
/Library/Developer/CommandLineTools/usr/bin/git commit -m "drone_dynamics: the paper's C1 log barrier (barrier_shape) and its summed gradient / sparse Hessian over pairs and a ball (barrier_terms), tested against complex steps of U and of the gradient

Co-Authored-By: Claude Fable 5.1 <noreply@anthropic.com>"
```

---

### Task 8: Folder runner, README, registration in `run_all_tests`

**Files:**
- Create: `optimal_control_examples/drone/run_tests.m`
- Create: `optimal_control_examples/drone/drone_dynamics/README.md`
- Create: `optimal_control_examples/drone/README.md`
- Modify: `optimal_control_examples/run_all_tests.m` (the `suites` cell, after `conjugate_points`)

**Interfaces:**
- Produces: `ok = run_tests()` in `drone/` (throws on failure), running the eight Stage 0 tests in dependency order; later stages append their tests to its `names` list.

- [ ] **Step 1: Write the runner**

```matlab
function ok = run_tests()
%% Purpose:
%
%   THE runner for optimal_control_examples/drone: every test in the
%   folder, cheapest and most fundamental first (the physics oracle before
%   anything that trusts the field), one logical verdict, and a THROW on
%   any failure so `matlab -batch "exit(~run_tests())"` returns a real exit
%   code (the ex3 runner's contract). Later stages append their tests to
%   names below, after the library's.
%
%% Inputs:
%
%  none
%
%% Outputs:
%
%  ok                       logical                 Every test passed (also
%                                                   THROWS on failure)
%
%% Revision History:
%  M. Casey                                                   (c) 10/23/2026
%  Copyright Coorbital Inc.
%% ------------------------ Begin Code Sequence ---------------------------

here = fileparts(mfilename('fullpath'));
root = fileparts(fileparts(here));
addpath(fullfile(here, 'drone_dynamics'), fullfile(here, 'drone_dynamics', 'tests'), ...
        fullfile(root, 'oclib'));

names = {'test_drone_params'          % the constants
         'test_euler_zyx'             % attitude kinematics vs longhand + successive rotations
         'test_drone_physics'         % THE ORACLE: mechanics, not a mirror
         'test_drone_field'           % form and complex-step safety
         'test_drone_state_jac'       % generated Jacobian vs complex step + FD, regeneration
         'test_layouts'               % the paper's J_LB column, pair_list
         'test_clearance_terms'       % gaps, unit vectors, hermite_resample
         'test_barrier_terms'};       % the barrier and its derivatives

res = false(1, numel(names));
secs = zeros(1, numel(names));
for k = 1:numel(names)
    fprintf('\n==== %s\n', names{k});
    t0 = tic;
    try
        v = feval(names{k});
        if ~isscalar(v)
            error('run_tests:badVerdict', '%s did not return a scalar verdict', names{k});
        end
        res(k) = logical(v);
    catch err
        fprintf('  THREW  %s (%s)\n', err.message, err.identifier);
        res(k) = false;
    end
    secs(k) = toc(t0);
end

fprintf('\n==== SUMMARY (drone) ====\n');
for k = 1:numel(names)
    if res(k), verdict = 'PASS'; else, verdict = 'FAIL'; end
    fprintf('  %-28s %-4s %5.0f s\n', names{k}, verdict, secs(k));
end
ok = all(res);
fprintf('  %d of %d passed, %.0f s total\n', nnz(res), numel(res), sum(secs));
if ~ok
    error('run_tests:failed', '%d of %d tests FAILED: %s', nnz(~res), numel(res), strjoin(names(~res), ', '));
end
end
```

- [ ] **Step 2: Register the suite in `run_all_tests.m`**

In `optimal_control_examples/run_all_tests.m`, change the `suites` cell to:

```matlab
suites = { 'cartpole_common',   {'test_cartpole_params', 'test_cartpole_physics', ...
                                 'test_cartpole_field',  'test_cartpole_state_jac'}
           'ex3_cart_pole_pmp', {'run_tests'}
           'conjugate_points',  {'run_conjugate_tests'}
           'drone',             {'run_tests'} };
```

Read the name-collision comment block already in that file (it foresaw a second `run_tests`): the loop resolves by `addpath` of the suite folder immediately before `feval`, and MATLAB's current-folder precedence can override it. Make the dispatch unambiguous by replacing the `feval(names{k})`-style call for entries named `run_tests` with a `cd` into the suite folder guarded by `onCleanup`:

```matlab
        if strcmp(suites{k,2}{m}, 'run_tests')
            back = cd(fullfile(here, folder));  restoreDir = onCleanup(@() cd(back));
            v = run_tests();
            clear restoreDir
        else
            v = feval(suites{k,2}{m});
        end
```

(keep the surrounding scalar-verdict guard and try/catch exactly as they are).

- [ ] **Step 3: Write the two READMEs**

`optimal_control_examples/drone/drone_dynamics/README.md`:

```markdown
# drone_dynamics

The 3D drone plant family, in one place, shared by every drone example
(the paper's kinematic multi-agent problem, the point-mass planner model,
and the 6-DOF quadrotor), so each model's constants, field and Jacobian
have exactly one home. The `cartpole_common` rule applies: nothing here
knows a cost, a horizon, a boundary condition or a solver.

## Models (`drone_params(model)`, `drone_field(x, p)`, `drone_state_jac(x, u, p)`)

| model | state | control | dynamics |
|---|---|---|---|
| `kinematic` | `p` (3) | `v` (3) | `pdot = v` -- Zhang, Liu, Darbon & Karniadakis (arXiv:2610.30388) |
| `pointmass` | `[p; v]` (6) | `a` (3) | `vdot = a - g e3` |
| `rigid6dof` | `[p; v; Theta; omega]` (12) | `f` (4) rotor thrusts | rigid body, ZYX Euler, "+" mixer, Crazyflie-class constants (Luis & Le Ny 2016; Forster 2015) |

All three are control-affine, `xdot = F(x) + B(x) u`, and complex-step
safe. `drone_state_jac` for `rigid6dof` is GENERATED (`gen_drone_jacs.m`,
Symbolic Math Toolbox) by differentiating `drone_field` itself; the
generated file is committed and `test_drone_state_jac` regenerates it and
compares byte for byte. The ZYX kinematics are singular at pitch 90
degrees; `p.thetaMax = 60 degrees` is the box the solvers impose.

## Geometry

| file | owns |
|---|---|
| `drone_layouts.m` | the paper's Table 1 layouts by M (inner ring radii inferred from its J_LB column), our `sphere` layout with a ball obstacle, the two-agent `headon` oracle layout |
| `pair_list.m` | all unordered pairs |
| `clearance_terms.m` | pair gaps `h = rho - 2r`, distances, unit vectors, ball gaps |
| `barrier_shape.m` / `barrier_terms.m` | the paper's C1 log barrier (eq. 3b) and its summed gradient and sparse Hessian |
| `lower_bound.m` | `J_LB = sum |xT - x0|^2 / (2T)` |
| `hermite_resample.m` | the trapezoid interpolant on a fine grid, for looking between nodes |

## The oracle

`tests/test_drone_physics.m` derives its expectations from mechanics and
geometry, never from the implemented rows: hover and free-fall limits, a
pointwise power balance `dE/dt = T (R e3).v + tau.omega` (the gyroscopic
term does no work), kinematic consistency `Rdot = R [omega]x`, and
conservation of the inertial angular momentum in torque-free flight. Every
other test compares the field against a differentiation of itself and says
so in its header.

## Consumers

`route1_direct/`, `route2_indirect/`, `quadrotor/` (Stages 1-3 of
`docs/superpowers/specs/2026-10-23-drone-multiagent-design.md`).
```

`optimal_control_examples/drone/README.md`:

```markdown
# drone

Multi-agent drone path planning, from the paper's kinematic drones to a
6-DOF quadrotor, solved directly (CasADi/IPOPT) and indirectly (PMP shot
with `oc.ms_bvp`) and certified. The design is
`docs/superpowers/specs/2026-10-23-drone-multiagent-design.md`; the paper is
Zhang, Liu, Darbon & Karniadakis, arXiv:2610.30388 (`papers/`).

| folder | stage | status |
|---|---|---|
| `drone_dynamics/` | 0: the plant library and shared geometry | built |
| `route1_direct/` | 1: the paper's problem as a hard-clearance NLP (Table 1) | planned |
| `route2_indirect/` | 2: the barrier form's PMP-BVP, seeded from Route 1, conjugate-tested | planned |
| `quadrotor/` | 3: the same task on the rigid-body model | planned |

Run the tests: `cd` here and `run_tests`, or `run_all_tests` from
`optimal_control_examples/`.
```

- [ ] **Step 4: Run the folder runner and the tree runner**

Run: `matlab -batch "cd('optimal_control_examples/drone'); exit(~run_tests())"` then `matlab -batch "cd('optimal_control_examples'); exit(~run_all_tests())"`.
Expected: `8 of 8 passed` for the folder; the tree runner reports the four suites all `PASS` and exits 0.

- [ ] **Step 5: Commit**

```bash
/Library/Developer/CommandLineTools/usr/bin/git add optimal_control_examples/drone/run_tests.m optimal_control_examples/drone/README.md optimal_control_examples/drone/drone_dynamics/README.md optimal_control_examples/run_all_tests.m
/Library/Developer/CommandLineTools/usr/bin/git commit -m "drone: folder test runner (8 library tests, oracle first), READMEs, suite registered in run_all_tests with an unambiguous run_tests dispatch

Co-Authored-By: Claude Fable 5.1 <noreply@anthropic.com>"
```
