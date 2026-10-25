# Chaplygin sleigh and the mechanics library: design

Date: 2026-10-25. Status: approved in conversation (Mike, 2026-10-25); the
plan argues from this document.

## 1. Purpose

A user-friendly MATLAB demonstration of the Chaplygin sleigh, the first of
a series of nonholonomic systems (rolling disk, vertical rolling coin, knife
edge on a slope to follow). The user sets the blade-to-center-of-mass
distance and the initial velocity, the script integrates the equations of
motion from 0 to a user-defined final time, and it produces diagnostic
static plots and a high-quality movie. Everything that is not the sleigh
lives in a new mechanics library so the next system reuses it.

Decisions already made (do not reopen):

- NOT in oclib. oclib's contract is optimal-control structure, physics-free;
  plants live beside their consumers (cartpole_common, drone_dynamics
  precedent); a future optimal-control problem on the sleigh will consume
  the mechanics library as a plant, never the reverse.
- The series is real: three generic modules enter the library now, each
  with a test on a problem with an exact answer.
- The initial velocity is (v0, w0): forward speed along the blade and yaw
  rate. The center-of-mass velocity follows from the constraint.
- A LaTeX note deriving the equations is in scope.

## 2. Layout

```
lagrangianHamiltonianMechanics/
  README.md                       what the folder is, the series, how to run
  run_tests.m                     one runner: mechlib tests + every system's tests
  mechlib/
    README.md                     admission rule, contents table, tests
    +mech/lagrange_dalembert.m
    +mech/conservation_check.m
    +mech/animate_planar_body.m
    tests/test_lagrange_dalembert.m
    tests/test_conservation_check.m
    tests/test_animate_planar_body.m
  nonholonomic/ChaplyginSleigh/
    README.md
    run_sleigh.m                  the user script (front door)
    sleigh_params.m
    sleigh_lagrangian.m
    sleigh_rhs.m
    sleigh_quantities.m
    sleigh_body.m
    plot_sleigh.m
    tests/test_sleigh_oracle.m
    tests/test_sleigh_pipeline.m
    doc/chaplygin_sleigh_note.tex (+ .pdf, figures/)
    results/                      git-ignored outputs of run_sleigh
```

Call library functions as `mech.<fn>` after `addpath('.../mechlib')`.
`run_tests.m` sits at `lagrangianHamiltonianMechanics/`, adds mechlib and
each system folder plus its `tests/`, runs every test function by name,
prints the summary table in the style of `drone/run_tests.m`, and returns
one logical. The repo has no repo-wide runner; `optimal_control_examples/
run_all_tests.m` covers only its own tree and is not touched.

**Library admission rule (mechlib/README.md):** a function enters mechlib
when it is written without any system-specific quantity (no mass, length
or named coordinate of a particular system in its body) and carries a test
on a problem with an exact answer. System folders keep nothing generic.

## 3. The physics (conventions are binding)

World frame: x east, y north, heading theta counterclockwise from x. Body
axes e1 = (cos th, sin th) along the blade (forward), e2 = (-sin th, cos th)
(left). Generalized coordinates q = (x_c, y_c, th) of the center of mass.
Mass m, inertia J about the center of mass, signed blade offset a: the
blade point is r_b = r_c - a e1, so a > 0 puts the center of mass AHEAD of
the blade in the +e1 direction.

Lagrangian (no potential): L = m/2 (xd_c^2 + yd_c^2) + J/2 thd^2, so
M = diag(m, m, J).

Blade velocity: rd_b = rd_c - a thd e2. The knife edge forbids sideways
motion at the blade: rd_b . e2 = 0, i.e.

    A(q) qd = 0,  A = [-sin th,  cos th,  -a].

Lagrange-d'Alembert with multiplier lambda (the sideways constraint force
on the sleigh, applied at the blade along e2):

    m xdd_c = -lambda sin th,   m ydd_c = lambda cos th,   J thdd = -a lambda.

Quasi-velocities v = rd_b . e1 = xd_c cos th + yd_c sin th (forward speed
of the blade point, equal to the forward speed of the center of mass) and
w = thd. Then rd_c = v e1 + a w e2 and the reduced equations are

    xd_b = v cos th,  yd_b = v sin th,  thd = w,
    vd = a w^2,       wd = -k v w,      k = m a / (J + m a^2),
    lambda = m J v w / (J + m a^2).

The script integrates the reduced state y = [x_b; y_b; th; v; w] of the
BLADE point; sleigh_quantities converts to and from the center-of-mass
description.

Conserved energy E = m/2 v^2 + (J + m a^2)/2 w^2. Invariant ellipse in the
(v, w) plane: v^2 + ((J + m a^2)/m) w^2 = C, C = 2E/m.

Exact solution (a ~= 0 and w0 ~= 0): vd = k (C - v^2), so with b =
sqrt(Jb/m), S = sqrt(C) and the WELL-CONDITIONED parametrization

    alpha = asinh(v0 / (b |w0|)),   z = alpha + k S t,
    v(t) = S tanh z,                w(t) = sign(w0) (S/b) sech z,
    theta_inf - theta_0 = sign(w0) (b/|a|) acos(sign(a) v0 / S).

(The equivalent t0 = -atanh(v0/S)/(k S) form is NOT to be evaluated: for
tiny w0, C rounds to v0^2 and atanh(+-1) is infinite.) w never changes
sign; v(t) -> sign(a) S. Straight motions (v, 0) form a LINE of equilibria:
the one with the center of mass leading (a v > 0) is Lyapunov stable and
transversely attracting (asymptotically stable on its energy ellipse; a
perturbation along the line is another equilibrium and never decays); the
trailing one (a v < 0) is unstable to yaw, and a trailing launch with
w0 ~= 0 sees |w| grow, its signed speed cross zero exactly once at
t = -atanh(v0/S)/(k S), and departs along sign(a). No energy is dissipated:
the constraint force does no net work; yaw energy J w^2/2 is traded into
translation (d/dt (J w^2/2) = -J k v w^2) and back. a = 0 is the knife
edge: v, w constant, a circle of radius |v0/w0|, lambda = m v w. With
w0 = 0 the motion is straight at speed v0 for every a (an exact equilibrium,
stable or not); the asymptote is then v0, not sign(a) S. The constraint is
nonholonomic by Frobenius: alpha ^ d alpha = -dx dy dth ~= 0 for every a.

## 4. Library modules

### 4.1 `sys = mech.lagrange_dalembert(q, qd, L, A, Q)`

Symbolic input: `q`, `qd` [n x 1] sym, `L` sym scalar, `A` [k x n] sym
(k = 0 allowed, pass `sym(zeros(0, n))`), `Q` [n x 1] sym generalized
forces (optional, zeros). Forms M = d2L/dqd2, the remaining Euler-Lagrange
terms h(q, qd) = d/dt(dL/dqd)|_(explicit in q, qd) - dL/dq, and
Ad = sum_j dA/dq_j qd_j, and solves

    [M  A'; A  0] [qdd; -lambda] = [Q - h; -Ad qd].

Returns a struct: `f(t, y)` numeric first-order right-hand side with
y = [q; qd] (matlabFunction, column outputs); `lambda(q, qd)` [k x 1];
`energy(q, qd)` = qd' dL/dqd - L; `constraint(q, qd)` = A qd; `qddSym`,
`lambdaSym` the symbolic reduced expressions; `n`, `k`. Refusals with
identifiers: `mech:lagrange_dalembert:M` when M is symbolically singular,
`:A` when A has deficient rank (symbolic rank), `:shape` on mismatched sizes.
Self-demo (`nargin == 0`): the planar pendulum as a constrained particle,
prints lambda against the closed-form tension.

### 4.2 `chk = mech.conservation_check(t, Y, sys, bars)`

`Y` [2n x N] columns [q; qd] on the grid `t`. Returns `.E` [1 x N],
`.driftRel` = max |E - E(1)| / max(|E(1)|, eps), `.resid` = max over the
grid of max |A qd|, `.ok` = both under `bars.energy`, `bars.constraint`
(defaults 1e-9), `.why` text. Prints nothing. Self-demo: a free particle on
a circle (exact), then an injected drift that fails.

### 4.3 `file = mech.animate_planar_body(t, pose, body, stem, opts)`

`pose` [3 x N] = (x, y, th) of the body's reference point in the world
frame. `body` struct: `.outline` [2 x P] polygon in body coordinates
(closed by the function); `.markers` struct array with `.name`, `.pos`
[2 x 1] body coordinates, `.color`; `.segments` struct array with `.pts` [2 x 2]
body-coordinate endpoints, `.color`, `.width` (decorations such as a blade);
`.arrows` struct array with `.name`, `.fun` = @(tk) returning [base; vec]
[4 x 1] in WORLD coordinates at movie TIME tk (the caller does not know the
movie grid), `.color`; `.subtitle` = @(tk) a char line. `opts`: `.fps`
[30], `.seconds` [t(end)], `.traces` [2 x N x T] world paths to grow with
time (default the reference point), `.title`. Fixed axis limits from the
whole path plus the body size, equal aspect, dark theme as movie_quad,
frames forced to exact 1280 x 720 (H.264 shears other sizes), MPEG-4 at
quality 95. Interpolates pose onto the movie's own grid with pchip. Returns
the MP4 path. Self-demo: a rectangle sliding on a circle, written to
tempdir.

## 5. Sleigh modules

- `p = sleigh_params(m, J, a)`: struct with m, J, a, and derived `k`,
  `Jb = J + m a^2`. No arguments gives the defaults m = 1, J = 0.1, a = 0.2.
  Refuses non-positive m or J (`sleigh_params:value`).
- `[q, qd, L, A] = sleigh_lagrangian(p)`: the symbolic objects of section 3
  for the builder. Nothing else.
- `yd = sleigh_rhs(t, y, p)`: the hand-derived reduced equations on
  y = [x_b; y_b; th; v; w]. Hand-written, no symbolic call, vectorized over
  columns of y so the oracle can evaluate a whole solution at once.
- `out = sleigh_quantities(t, Y, p)`: from the reduced solution, `.E`,
  `.lambda`, `.C`, the exact `.vExact`, `.wExact` (constants when a = 0 or
  w0 = 0, which are exact there), `.rc` center-of-mass path [2 x N], `.rb` blade
  path, `.Yfull` = [q; qd] for the builder's checks, and `.asymptote` =
  sign(a) sqrt(C). The script builds y0 directly; no inverse is needed.
- `arrows = sleigh_arrows(t, Y, out, p)`: the movie's two arrows (center-of-
  mass velocity, blade constraint force lambda e2) as functions of time,
  scaled to one body length; shared by the script and the pipeline test.
- `body = sleigh_body(p)`: outline polygon (a pointed sleigh shape whose
  length scales with max(|a|, 0.1), blade as a thick segment along e1 at
  the blade point), markers (center of mass, blade), reference point = the
  blade point.
- `files = plot_sleigh(t, Y, out, p, info, outDir)`: two PNGs,
  `<label>_path.png` (center-of-mass and blade traces, heading ticks every
  0.5 s, start and end markers, equal aspect) and `<label>_diag.png`
  (2 x 3: v(t) with the exact solution and the asymptote, w(t) with the
  exact solution, th(t), energy drift, constraint residual |A qd|, and the
  (v, w) phase portrait on the energy ellipse). lambda(t) is drawn on the
  v(t) panel's right axis. Returns the file list. `info.visible` false by
  default so tests run headless.

### 5.1 `run_sleigh.m`, the front door

Script shaped like `drone/quad_flight/fly_quad.m`: a `%%` header (Purpose,
what the user sets, what it writes, Revision History with
`M. Casey  (c) 10/25/2026` and `Copyright Coorbital Inc.`), then

Section 1, user settings, each with a one-line comment:
```matlab
label = 'sleigh_forward';
m     = 1.0;     % mass (kg)
J     = 0.1;     % inertia about the center of mass (kg m^2)
a     = 0.2;     % blade to center of mass (m); a > 0 = center of mass ahead of the blade
v0    = 1.0;     % initial forward speed along the blade (m/s)
w0    = 2.0;     % initial yaw rate (rad/s)
tEnd  = 10;      % final time (s)
dtOut = 0.01;    % output step (s)
movieOpts = struct('fps', 30, 'seconds', tEnd);
```
Initial pose is the blade at the origin heading east (x_b = y_b = th = 0);
not a user setting.

Pipeline, each step printing one line:
1. addpath of mechlib and the sleigh folder relative to `mfilename('fullpath')`;
   `p = sleigh_params(m, J, a)`; print k, Jb, C and the asymptote.
2. ode113, RelTol = AbsTol = 1e-12, on the uniform grid 0:dtOut:tEnd.
3. `sys = mech.lagrange_dalembert(sleigh_lagrangian(p))`.
4. Diagnostics, one PASS/FAIL line each:
   - D1 hand-derived vs builder: max |sleigh_rhs - sys.f| on the full-state
     solution, after mapping the reduced derivative to [qd; qdd]; bar 1e-9.
   - D2 exact solution: max |v - vExact| and |w - wExact|; bar 1e-8
     (reported "n/a (w0 = 0)" when there is no tanh solution).
   - D3 energy drift and D4 constraint residual via mech.conservation_check,
     bars 1e-9.
   - D5 asymptotics: `a * v(end) > 0`, |w(end)| <= 1e-3 |w0| and the total
     heading change against theta_inf (sleigh_quantities.thetaInf) when
     tEnd is long enough, with the statement that the leading straight ride
     is stable to yaw (one of a line of equilibria, transversely
     attracting); NOTE lines for a = 0, for w0 = 0 (stable or unstable
     exact equilibrium) and for a trailing run that has not reversed yet
     (with the exact reversal time); printed as a verdict, not a gate, since
     a short tEnd legitimately fails it.
5. `plot_sleigh` and `mech.animate_planar_body` with `sleigh_body`, arrows
   for the center-of-mass velocity and the blade constraint force lambda e2,
   traces for both points, subtitle `t, v, w, lambda, E drift`; outputs in
   `results/` beside the script; paths printed.

### 5.2 Tests

Every test is a function returning `ok`, prints PASS/FAIL lines with the
measured number and its bar, and runs from a fresh clone (no data files).

- `test_lagrange_dalembert`: free particle (k = 0: qdd = 0, lambda [0 x 1]);
  planar pendulum as a constrained particle (L = m/2 |rd|^2 - m g y,
  A = [x, y]; lambda = -T/l with T = m g cos(phi) + m v^2/l checked at the
  bottom and at a generic angle to 1e-10; energy drift of an ode113
  integration over two periods < 1e-9; small-amplitude period within 1e-6
  of 2 pi sqrt(l/g)); the a = 0 knife edge via the sleigh's own
  Lagrangian is NOT used here (library tests know no sleigh) -- instead a
  particle with A = [-sin th, cos th, 0] on q = (x, y, th) and L with J:
  v and w constant, lambda = m v w; refusals `:M`, `:A`, `:shape`.
- `test_conservation_check`: exact circular motion passes with drift
  < 1e-12; a solution with an injected 1e-6 energy ramp FAILS with the
  energy reason; a constraint violation of 1e-6 FAILS with the constraint
  reason.
- `test_animate_planar_body`: 20-frame MP4 in tempdir exists, > 10 kB,
  frames 1280 x 720 (read back one frame with VideoReader); a pose with a
  single point (zero path extent) still gets finite limits.
- `test_sleigh_oracle`: for a in {+0.2, -0.2, 0} and v0 in {+1, -1}, w0 = 2,
  tEnd = 8: D1 <= 1e-9, D2 <= 1e-8 (a ~= 0), D3 and D4 <= 1e-9, w keeps its
  sign, v(end) within 1e-3 of sign(a) sqrt(C) for a ~= 0, v and w constant
  to 1e-10 for a = 0 and lambda = m v w; the trailing case (a v0 < 0)
  reverses: v changes sign exactly once; w0 = 0 gives a straight line
  (y_b = 0 to 1e-12, th constant). Also the mapping check: the center-of-
  mass velocity from sleigh_quantities satisfies the constraint to 1e-12.
- `test_sleigh_pipeline`: `checkcode('run_sleigh.m', '-id')` has no
  error-class messages; sleigh_params refusals; plot_sleigh writes its two
  PNGs to a temp folder; the animator writes a 20-frame MP4 from
  sleigh_body.

`lagrangianHamiltonianMechanics/run_tests.m` lists the five test functions
and must print `5 of 5 passed`.

### 5.3 The note `doc/chaplygin_sleigh_note.tex`

Sections: 1 Setup and conventions (figure of the body, e1/e2, a signed);
2 Lagrange-d'Alembert with the multiplier; 3 Reduction to (v, w) and the
constraint force; 4 Energy, the invariant ellipse and the tanh solution;
5 Stability of straight motion by the sign of a v; 6 Numerical check: the
script's D1-D5 lines and the two figures from `run_sleigh` (copied into
`doc/figures/`); 7 References (Bloch, Nonholonomic Mechanics and Control;
Neimark and Fufaev; Chaplygin 1911). The reduced equations printed by the
builder's `qddSym`/`lambdaSym` appear in section 3 as the symbolic
cross-check. Compiled with `/Library/TeX/texbin/pdflatex`, aux files
removed, PDF committed. No author signature line beyond the standard
`\author{}` block Mike's notes already use (Michael Casey; no Coorbital
line unless he asks -- see the no-signature-line rule).

## 6. House style (binding for every file)

MATLAB header quartet `%% Purpose / Inputs / Outputs / Revision History`
with `M. Casey  (c) 10/25/2026` and `Copyright Coorbital Inc.`; `if nargin
== 0` self-demo in every function (scripts and tests excepted); never `i`
or `j` as variables; never `norm` (write `sqrt(sum(x.^2))`); no `%#ok`; no
emojis; ASCII READMEs; `results/` git-ignored (`*.png`, `*.mp4`, `*.log`);
movies 1280 x 720 exactly. READMEs carry no author signature line.
Subagents never commit; the controller commits after each clean review
with `Co-Authored-By: Claude Fable 5.1 <noreply@anthropic.com>` via
`/Library/Developer/CommandLineTools/usr/bin/git`; the unrelated
`missiles/` deletions in the working tree are never staged. MATLAB runs
through the MCP tools in the shared R2026a session (tests are function
files: cd into `tests`, call by name). Symbolic Math Toolbox is licensed and
may be used by the builder only; `sleigh_rhs` stays hand-written.

## 7. Out of scope

Friction or a driven sleigh; variational integrators; the rolling disk and
the other systems (they get their own specs and reuse mechlib); a
repo-wide test runner; any optimal-control problem on the sleigh.
