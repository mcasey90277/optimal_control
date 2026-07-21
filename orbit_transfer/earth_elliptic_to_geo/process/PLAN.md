# Earth Elliptic → GEO Min-Fuel Reproduction — Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Reproduce the Haberkorn–Martinon–Gergaud (JGCD 2004) low-thrust min-fuel
LEO-ellipse→GEO transfer with our direct Sundman-collocation + energy→fuel-homotopy
machinery (they solved it indirectly). Spec: `earth_elliptic_to_geo/DESIGN.md`.

**Architecture:** One mode-switched CasADi+IPOPT solver core (`casadi_lt_2body`,
modes `'mintime'`/`'fixedtf'`) over inertial 2-body dynamics with a Sundman clock
carried by a cScale slack state; dynamically-exact tangential-propagation seeds;
ε:1→0 Bertrand–Épénoy homotopy; pluggable terminal (fixed rendezvous / free-L GEO
insertion manifold). Milestones M0→M3 gate against the paper's numbers.

**Tech Stack:** MATLAB R2025b (headless `-batch`), CasADi 3.7.0 (`~/casadi-3.7.0`,
bundled IPOPT/MUMPS), ode113 for seed propagation.

## Global Constraints

- MATLAB binary: `/Applications/MATLAB_R2025b.app/bin/matlab` ONLY (R2025a license broken — memory `use-matlab-2025b`).
- CasADi path: `~/casadi-3.7.0` (solver adds it itself; `CASADI_PATH` env overrides).
- Nondim units: LU = 42165 km (GEO radius ⇒ terminal P = 42165 km, NOT the paper p.5 typo 42125), TU = 13713.8 s, μ = 1, mass unit = m₀ = 1500 kg.
- Isp default 2000 s (spec open item #1; only M2's m_f gate depends on it).
- MATLAB house style: every function gets the full comment header (purpose/INPUTS/OUTPUTS/REFERENCES); NEVER use `i`/`j` as loop variables.
- CasADi/MATLAB gotcha: NEVER write chained bounds `a <= x <= b` — always two separate `opti.subject_to` calls (banked campaign gotcha).
- Sporadic uncatchable CasADi/IPOPT MEX fatal crash (~1 in 10 solves) kills MATLAB. Single solves: just rerun. Sweeps: every point saves its own `.mat` and skips existing files on rerun (resume pattern).
- `results/` stays untracked (matches `PSR/results` convention) — commit code, tests, docs only.
- All commits from repo root `/Users/msc/Desktop/optimal_control`; end commit messages with `Co-Authored-By: Claude Fable 5 <noreply@anthropic.com>`.
- Run all tests/commands from the project folder: `matlab -batch "cd '/Users/msc/Desktop/optimal_control/earth_elliptic_to_geo'; <script>"`.

**Key reference numbers (paper, for gates):** at T_max=10 N, c_tf=1.5: t_f,min ≈ 84.7 h
(≈22.2 TU), min-time ΔL ≈ 26.4 rad (law R1: (L_f,min−L₀)·T_max ≈ 264 rad·N), fuel
solution: 7.5 revs, 18 switches, m_f ≈ 1370–1375 kg, burns at apogee. Law R2:
c_Lf,opt ≈ 1.12·c_tf + 0.09 ⇒ fuel-stage L_f − L₀ ≈ 1.77 × 26.4 ≈ 47 rad.

---

### Task 1: Scaffold + parameters module

**Files:**
- Create: `earth_elliptic_to_geo/kepler_lt_params.m`
- Test: `earth_elliptic_to_geo/test_params.m`

**Interfaces:**
- Produces: `p = kepler_lt_params(thrustN, m0kg, ispS)` → struct with fields
  `g0, muKm3s2, LU_km, TU_s, VU_kms, AU_ms2, mu, thrustN, m0kg, ispS, Tmax, c, pSund`.
  All later tasks call this exact signature.

- [ ] **Step 1: Write the failing test**

```matlab
% TEST_PARAMS  Unit checks for kepler_lt_params (canonical units + paper BCs).
p = kepler_lt_params(10, 1500, 2000);
assert(abs(p.TU_s - 13713.8) < 1.0,      'TU wrong');
assert(abs(p.VU_kms - 3.0747) < 1e-3,    'VU wrong');
assert(abs(p.Tmax - 0.029735) < 1e-5,    'nondim thrust wrong');
assert(abs(p.c - 6.3790) < 1e-3,         'nondim exhaust velocity wrong');
assert(p.mu == 1 && p.pSund == 1.5);
% paper initial-orbit geometry in these units
P0 = 11625/p.LU_km;  e0 = 0.75;  a0 = P0/(1-e0^2);
assert(abs(P0 - 0.275703) < 1e-5);
assert(abs(a0*(1+e0) - 1.102810) < 1e-5);   % apogee ~46,500 km
assert(abs(a0*(1-e0) - 0.157544) < 1e-5);   % perigee ~6,643 km
assert(abs(2*atand(0.0612) - 7.0052) < 1e-3); % inclination ~7 deg
fprintf('test_params: ALL PASS\n');
```

- [ ] **Step 2: Run test to verify it fails**

Run: `/Applications/MATLAB_R2025b.app/bin/matlab -batch "cd '/Users/msc/Desktop/optimal_control/earth_elliptic_to_geo'; test_params"`
Expected: FAIL — `Unrecognized function or variable 'kepler_lt_params'`.

- [ ] **Step 3: Write the implementation**

```matlab
function p = kepler_lt_params(thrustN, m0kg, ispS)
% KEPLER_LT_PARAMS  Constants + canonical units, Earth 2-body low-thrust problem.
%
% Nondimensionalization: LU = GEO radius 42165 km, TU = sqrt(LU^3/mu_earth), so
% mu = 1, GEO circular speed = 1, GEO period = 2*pi. Mass unit = m0kg.
%
% INPUTS:
%   thrustN - max thrust [N]        (paper cases: 10, 5, 2.5, 1)
%   m0kg    - initial mass [kg]     (paper: 1500)
%   ispS    - specific impulse [s]  (default 2000; DESIGN.md open item 1)
%
% OUTPUTS:
%   p - struct: dimensional anchors .g0 .muKm3s2 .LU_km .TU_s .VU_kms .AU_ms2;
%       nondim .mu(=1) .Tmax (thrust @ m=1) .c (exhaust velocity) .pSund;
%       echo .thrustN .m0kg .ispS
%
% REFERENCES:
%   [1] Haberkorn, Martinon, Gergaud, JGCD 27(6), 2004 (problem constants).
%   [2] earth_elliptic_to_geo/DESIGN.md sec 2 (units decision).
if nargin < 3, ispS = 2000; end
p.g0      = 9.80665;                          % [m/s^2]
p.muKm3s2 = 398600.47;                        % [km^3/s^2]
p.LU_km   = 42165;                            % GEO radius = terminal P [km]
p.TU_s    = sqrt(p.LU_km^3 / p.muKm3s2);      % => mu = 1
p.VU_kms  = p.LU_km / p.TU_s;
p.AU_ms2  = 1000 * p.VU_kms / p.TU_s;         % acceleration unit [m/s^2]
p.mu      = 1;
p.thrustN = thrustN;  p.m0kg = m0kg;  p.ispS = ispS;
p.Tmax    = (thrustN/m0kg) / p.AU_ms2;        % nondim thrust at m = 1
p.c       = (ispS*p.g0/1000) / p.VU_kms;      % nondim exhaust velocity
p.pSund   = 1.5;                              % Sundman power, dt/dtau = r^pSund
end
```

- [ ] **Step 4: Run test to verify it passes**

Same command. Expected: `test_params: ALL PASS`, exit 0.

- [ ] **Step 5: Commit**

```bash
cd /Users/msc/Desktop/optimal_control
git add earth_elliptic_to_geo/kepler_lt_params.m earth_elliptic_to_geo/test_params.m
git commit -m "feat(earth-geo): scaffold + canonical units/params module

Co-Authored-By: Claude Fable 5 <noreply@anthropic.com>"
```

---

### Task 2: Toolchain + homotopy sanity test (paper's toy P₂)

**Files:**
- Test: `earth_elliptic_to_geo/test_p2_homotopy.m` (standalone; no production module)

**Interfaces:**
- Consumes: CasADi at `~/casadi-3.7.0`.
- Produces: proof the CasADi+IPOPT install and the ε-homotopy objective form work,
  against an analytic answer, before the big solver exists.

The paper's toy (P₂): min ∫₀²|u|dt, ẍ=u, |u|≤1, (0,0)→(0.5,0). Analytic optimum:
bang-off-bang, switch t₁ = 1−1/√2 ≈ 0.29289, cost 2−√2 ≈ 0.585786. We transcribe
with the split control u = u⁺−u⁻ (u± ∈ [0,1]) and apply OUR homotopy per component:
J(ε) = ∫(u⁺+u⁻) − ε∫(u⁺(1−u⁺)+u⁻(1−u⁻)) — at ε=1 this is ∫(u⁺²+u⁻²) (energy), at
ε=0 it is ∫|u| (fuel).

- [ ] **Step 1: Write the test**

```matlab
% TEST_P2_HOMOTOPY  Paper toy P2: energy->fuel homotopy vs analytic bang-off-bang.
cp = getenv('CASADI_PATH'); if isempty(cp), cp = fullfile(getenv('HOME'),'casadi-3.7.0'); end
addpath(cp);
N = 200;  dt = 2/N;
opti = casadi.Opti();
X  = opti.variable(2, N+1);                 % [x; v]
Up = opti.variable(1, N+1);  Um = opti.variable(1, N+1);
u  = Up - Um;
opti.subject_to(Up >= 0);  opti.subject_to(Up <= 1);
opti.subject_to(Um >= 0);  opti.subject_to(Um <= 1);
for k = 1:N     % trapezoid defects for [xdot; vdot] = [v; u]
    fk  = [X(2,k);   u(k)];
    fk1 = [X(2,k+1); u(k+1)];
    opti.subject_to(X(:,k+1) - X(:,k) - (dt/2)*(fk+fk1) == 0);
end
opti.subject_to(X(:,1)   == [0; 0]);
opti.subject_to(X(:,end) == [0.5; 0]);
opti.set_initial(X, [linspace(0,0.5,N+1); zeros(1,N+1)]);
runc = @(w,epsv) w - epsv*(w.*(1-w));       % per-component homotopy integrand
for epsv = [1 0.6 0.3 0.12 0.04 0.01 0]
    g = runc(Up,epsv) + runc(Um,epsv);
    opti.minimize( sum((dt/2)*(g(1:N)+g(2:N+1))) );
    opti.solver('ipopt', struct('print_time',0), struct('print_level',0,'max_iter',800));
    sol = opti.solve();
    opti.set_initial(X,  sol.value(X));
    opti.set_initial(Up, sol.value(Up));
    opti.set_initial(Um, sol.value(Um));
end
uv = sol.value(Up) - sol.value(Um);
tg = linspace(0,2,N+1);
cost = trapz(tg, abs(uv));
assert(abs(cost - (2-sqrt(2))) < 3e-3, 'P2 cost mismatch: %.5f', cost);
assert(uv(1) > 0.95 && uv(end) < -0.95, 'not bang at ends');
assert(abs(uv(round(N/2))) < 0.05, 'not coasting at midpoint');
kSw = find(uv < 0.5, 1);                    % first departure from the +1 arc
assert(abs(tg(kSw) - (1-1/sqrt(2))) < 0.03, 'switch time off');
fprintf('test_p2_homotopy: ALL PASS (cost %.6f vs %.6f)\n', cost, 2-sqrt(2));
```

- [ ] **Step 2: Run it**

Run: `/Applications/MATLAB_R2025b.app/bin/matlab -batch "cd '/Users/msc/Desktop/optimal_control/earth_elliptic_to_geo'; test_p2_homotopy"`
Expected: `test_p2_homotopy: ALL PASS (cost 0.5858xx vs 0.585786)`. If CasADi fails to
load, stop and report — nothing downstream can work.

- [ ] **Step 3: Commit**

```bash
cd /Users/msc/Desktop/optimal_control
git add earth_elliptic_to_geo/test_p2_homotopy.m
git commit -m "test(earth-geo): CasADi toolchain + energy->fuel homotopy vs analytic P2

Co-Authored-By: Claude Fable 5 <noreply@anthropic.com>"
```

---

### Task 3: Element ↔ Cartesian conversions

**Files:**
- Create: `earth_elliptic_to_geo/elements_to_cart.m`, `earth_elliptic_to_geo/cart_to_elements.m`
- Test: `earth_elliptic_to_geo/test_elements.m`

**Interfaces:**
- Produces: `[r, v] = elements_to_cart(P, ex, ey, hx, hy, L, mu)` (r,v 3×1 column);
  `el = cart_to_elements(r, v, mu)` → struct `.P .ex .ey .hx .hy .L`.

- [ ] **Step 1: Write the failing test**

```matlab
% TEST_ELEMENTS  Forward conversion invariants + roundtrip.
p  = kepler_lt_params(10, 1500, 2000);
P0 = 11625/p.LU_km;  e0 = 0.75;  a0 = P0/(1-e0^2);
% (a) paper initial state: apogee, i=7 deg
[r0, v0] = elements_to_cart(P0, 0.75, 0, 0.0612, 0, pi, p.mu);
assert(abs(norm(r0) - a0*(1+e0)) < 1e-10, 'not at apogee radius');
assert(abs(dot(r0,v0)) < 1e-10,           'radial rate nonzero at apsis');
vis = sqrt(p.mu*(2/norm(r0) - 1/a0));
assert(abs(norm(v0) - vis) < 1e-10,       'vis-viva violated');
hv = cross(r0, v0);
assert(abs(acosd(hv(3)/norm(hv)) - 7.0052) < 1e-3, 'inclination wrong');
% (b) GEO check: equatorial circular prograde at any L
[rg, vg] = elements_to_cart(1, 0, 0, 0, 0, 0.7, p.mu);
assert(abs(norm(rg)-1) < 1e-12 && abs(norm(vg)-1) < 1e-12);
assert(abs(rg(3)) < 1e-12 && abs(vg(3)) < 1e-12 && abs(dot(rg,vg)) < 1e-12);
hg = cross(rg, vg);  assert(hg(3) > 0, 'retrograde GEO');
% (c) roundtrip on a grid (incl. wrap-aware L compare)
rng(7);
for kk = 1:50
    el = [0.2+rand, 0.7*(rand-0.5), 0.7*(rand-0.5), 0.2*(rand-0.5), ...
          0.2*(rand-0.5), 2*pi*rand-pi];
    [rr, vv] = elements_to_cart(el(1), el(2), el(3), el(4), el(5), el(6), p.mu);
    eb = cart_to_elements(rr, vv, p.mu);
    assert(max(abs([eb.P eb.ex eb.ey eb.hx eb.hy] - el(1:5))) < 1e-10, 'roundtrip elems');
    dL = mod(eb.L - el(6) + pi, 2*pi) - pi;
    assert(abs(dL) < 1e-10, 'roundtrip L');
end
fprintf('test_elements: ALL PASS\n');
```

- [ ] **Step 2: Run to verify it fails** (function not defined), same command pattern.

- [ ] **Step 3: Write the implementations**

```matlab
function [r, v] = elements_to_cart(P, ex, ey, hx, hy, L, mu)
% ELEMENTS_TO_CART  Paper/MEE-style elements -> inertial Cartesian state.
%
% Elements per Haberkorn-Martinon-Gergaud 2004: P (semi-latus rectum),
% (ex,ey) = e*(cos,sin)(Om+om), (hx,hy) = tan(i/2)*(cos,sin)(Om), L = Om+om+theta.
%
% INPUTS:  P,ex,ey,hx,hy,L - elements [scalars];  mu - grav parameter [scalar]
% OUTPUTS: r, v - inertial position/velocity [3x1 each]
%
% REFERENCES:
%   [1] Walker/Betts modified-equinoctial <-> Cartesian formulas.
w  = 1 + ex*cos(L) + ey*sin(L);
s2 = 1 + hx^2 + hy^2;
a2 = hx^2 - hy^2;
rm = P / w;
r  = (rm/s2) * [cos(L) + a2*cos(L) + 2*hx*hy*sin(L);
                sin(L) - a2*sin(L) + 2*hx*hy*cos(L);
                2*(hx*sin(L) - hy*cos(L))];
sq = sqrt(mu/P);
v  = (1/s2) * [-sq*( sin(L) + a2*sin(L) - 2*hx*hy*cos(L) + ey - 2*ex*hx*hy + a2*ey);
               -sq*(-cos(L) + a2*cos(L) + 2*hx*hy*sin(L) - ex + 2*ey*hx*hy + a2*ex);
                2*sq*(hx*cos(L) + hy*sin(L) + ex*hx + ey*hy)];
end
```

```matlab
function el = cart_to_elements(r, v, mu)
% CART_TO_ELEMENTS  Inertial Cartesian state -> paper/MEE-style elements.
%
% INPUTS:  r, v - inertial position/velocity [3x1];  mu - grav parameter
% OUTPUTS: el - struct .P .ex .ey .hx .hy .L   (L in (-pi, pi])
%
% REFERENCES: inverse of elements_to_cart (roundtrip-tested).
r = r(:);  v = v(:);
hv = cross(r, v);
el.P  = dot(hv,hv)/mu;
hn = hv/norm(hv);
el.hx = -hn(2)/(1+hn(3));                    % tan(i/2)cos(Om)
el.hy =  hn(1)/(1+hn(3));                    % tan(i/2)sin(Om)
ev = cross(v, hv)/mu - r/norm(r);            % Laplace eccentricity vector
s2 = 1 + el.hx^2 + el.hy^2;
fh = [1+el.hx^2-el.hy^2;  2*el.hx*el.hy;      -2*el.hy] / s2;   % equinoctial basis
gh = [2*el.hx*el.hy;      1-el.hx^2+el.hy^2;   2*el.hx] / s2;
el.ex = dot(ev, fh);
el.ey = dot(ev, gh);
el.L  = atan2(dot(r,gh), dot(r,fh));
end
```

- [ ] **Step 4: Run test to verify it passes.** If only the roundtrip fails, the bug is a
sign in `fh`/`gh`/`hx`/`hy` — check against the forward formula's basis (forward r-vector
= rm·(cosL·fh + sinL·gh) must reproduce elements_to_cart exactly).

- [ ] **Step 5: Commit**

```bash
cd /Users/msc/Desktop/optimal_control
git add earth_elliptic_to_geo/elements_to_cart.m earth_elliptic_to_geo/cart_to_elements.m earth_elliptic_to_geo/test_elements.m
git commit -m "feat(earth-geo): MEE-style element <-> Cartesian conversions, roundtrip-tested

Co-Authored-By: Claude Fable 5 <noreply@anthropic.com>"
```

---

### Task 4: Time-domain dynamics helper

**Files:**
- Create: `earth_elliptic_to_geo/lt2b_rhs_time.m`
- Test: `earth_elliptic_to_geo/test_dynamics.m`

**Interfaces:**
- Produces: `xdot = lt2b_rhs_time(x, u, par)` — x = [r(3);v(3);m;t] (8×1),
  u = [alpha(3);s] (4×1), par from `kepler_lt_params`. Returns d/dt of x (8×1).
  MUST work on numeric doubles AND CasADi MX (no norm/abs/max on the state).

- [ ] **Step 1: Write the failing test**

```matlab
% TEST_DYNAMICS  Ballistic invariants + thrust mass-rate exactness.
p = kepler_lt_params(10, 1500, 2000);
P0 = 11625/p.LU_km;  a0 = P0/(1-0.75^2);
[r0, v0] = elements_to_cart(P0, 0.75, 0, 0.0612, 0, pi, p.mu);
x0 = [r0; v0; 1; 0];
% (a) one ballistic period: energy/|h| conserved, state returns
T0 = 2*pi*a0^1.5;
oo = odeset('RelTol',1e-12,'AbsTol',1e-13);
[~, xx] = ode113(@(t,x) lt2b_rhs_time(x, [1;0;0;0], p), [0 T0], x0, oo); % s=0: alpha moot
xe = xx(end,:).';
E  = @(x) 0.5*dot(x(4:6),x(4:6)) - p.mu/norm(x(1:3));
assert(abs(E(xe) - E(x0)) < 1e-9, 'energy drift');
assert(norm(cross(xe(1:3),xe(4:6)) - cross(r0,v0)) < 1e-9, 'h drift');
assert(norm(xe(1:6) - x0(1:6)) < 1e-6, 'period return failed');
assert(abs(xe(7) - 1) < 1e-14 && abs(xe(8) - T0) < 1e-9, 'm/t states wrong');
% (b) full thrust for 1 TU: exact linear mass, energy increases (tangential)
odef = @(t,x) lt2b_rhs_time(x, [x(4:6)/norm(x(4:6)); 1], p);
[~, xt] = ode113(odef, [0 1], x0, oo);
assert(abs(xt(end,7) - (1 - p.Tmax/p.c)) < 1e-10, 'mass rate wrong');
assert(E(xt(end,:).') > E(x0), 'tangential thrust must raise energy');
fprintf('test_dynamics: ALL PASS\n');
```

- [ ] **Step 2: Run to verify it fails.**

- [ ] **Step 3: Write the implementation**

```matlab
function xdot = lt2b_rhs_time(x, u, par)
% LT2B_RHS_TIME  Time-domain EOM: inertial 2-body gravity + low thrust (8-state).
%
% x = [r(3); v(3); m; t], u = [alpha(3); s] with ||alpha||=1, s in [0,1].
% Thrust accel = (Tmax/m)*s*alpha; mdot = -(Tmax/c)*s; tdot = 1 (time carried
% as a state so the Sundman solver can pin t(tau_f)). Written without
% norm/abs/max so it evaluates on BOTH numeric doubles and CasADi MX.
%
% INPUTS:  x [8x1], u [4x1], par from kepler_lt_params
% OUTPUTS: xdot [8x1] = d/dt [r; v; m; t]
%
% REFERENCES: [1] DESIGN.md sec 2 (problem statement).
r = x(1:3);  v = x(4:6);  m = x(7);
rn2 = r(1)^2 + r(2)^2 + r(3)^2 + 1e-12;      % softened, AD/CS-safe
acc = -par.mu * r * rn2^(-1.5) + (par.Tmax/m) * u(4) * u(1:3);
xdot = [v; acc; -(par.Tmax/par.c)*u(4); 1];
end
```

- [ ] **Step 4: Run test to verify it passes.**

- [ ] **Step 5: Commit**

```bash
cd /Users/msc/Desktop/optimal_control
git add earth_elliptic_to_geo/lt2b_rhs_time.m earth_elliptic_to_geo/test_dynamics.m
git commit -m "feat(earth-geo): shared 2-body low-thrust EOM (numeric + CasADi), invariant-tested

Co-Authored-By: Claude Fable 5 <noreply@anthropic.com>"
```

---

### Task 5: Terminal-condition builder

**Files:**
- Create: `earth_elliptic_to_geo/geo_terminal.m`
- Test: `earth_elliptic_to_geo/test_terminal.m`

**Interfaces:**
- Produces: `term = geo_terminal(mode, par, Lf)`:
  - mode `'fixed'`: `term.type='fixed'`, `term.rvf` [6×1] GEO state at longitude Lf, `term.Lf`.
  - mode `'manifold'`: `term.type='manifold'`, `term.aGeo=1`,
    `term.resid = @(rv) [5×1 residuals]` (numeric-use helper; solver re-poses them
    symbolically). Residuals: [r_z; v_z; ‖r‖²−a²; ‖v‖²−μ/a; r·v].

- [ ] **Step 1: Write the failing test**

```matlab
% TEST_TERMINAL  Fixed-GEO state properties + manifold residuals.
p = kepler_lt_params(10, 1500, 2000);
tf1 = geo_terminal('fixed', p, 2.3);
rv = tf1.rvf;  r = rv(1:3);  v = rv(4:6);
assert(abs(norm(r)-1) < 1e-12 && abs(norm(v)-1) < 1e-12);
assert(abs(r(3)) < 1e-12 && abs(v(3)) < 1e-12 && abs(dot(r,v)) < 1e-12);
tm = geo_terminal('manifold', p, []);
assert(max(abs(tm.resid(rv))) < 1e-12, 'GEO state must satisfy manifold');
rvBad = rv;  rvBad(1) = rvBad(1)*1.01;       % radius + radial-rate violated
res = tm.resid(rvBad);
assert(abs(res(3)) > 1e-3, 'radius constraint insensitive');
fprintf('test_terminal: ALL PASS\n');
```

- [ ] **Step 2: Run to verify it fails.**

- [ ] **Step 3: Write the implementation**

```matlab
function term = geo_terminal(mode, par, Lf)
% GEO_TERMINAL  Terminal-condition builder for the GEO target.
%
% mode 'fixed'    - full 6-state rendezvous at longitude Lf on GEO (M0/M1).
% mode 'manifold' - free-longitude insertion manifold (M2+): 5 residuals
%                   [r_z; v_z; ||r||^2-a^2; ||v||^2-mu/a; r.v] = 0. NB the set
%                   also admits the retrograde orbit; the prograde seed selects
%                   the branch (DESIGN.md sec 2).
%
% INPUTS:  mode - 'fixed' | 'manifold';  par - kepler_lt_params struct;
%          Lf - GEO longitude [rad] ('fixed' only; pass [] for 'manifold')
% OUTPUTS: term - struct (.type; .rvf/.Lf for fixed; .aGeo/.resid for manifold)
%
% REFERENCES: [1] DESIGN.md sec 2 (boundary conditions).
switch lower(mode)
    case 'fixed'
        [rf, vf] = elements_to_cart(1, 0, 0, 0, 0, Lf, par.mu);
        term = struct('type','fixed', 'rvf', [rf; vf], 'Lf', Lf);
    case 'manifold'
        a = 1;
        res = @(rv) [rv(3); rv(6); ...
                     rv(1)^2+rv(2)^2+rv(3)^2 - a^2; ...
                     rv(4)^2+rv(5)^2+rv(6)^2 - par.mu/a; ...
                     rv(1)*rv(4)+rv(2)*rv(5)+rv(3)*rv(6)];
        term = struct('type','manifold', 'aGeo', a, 'resid', res);
    otherwise
        error('geo_terminal:mode', 'unknown mode %s', mode);
end
end
```

- [ ] **Step 4: Run test to verify it passes.**

- [ ] **Step 5: Commit**

```bash
cd /Users/msc/Desktop/optimal_control
git add earth_elliptic_to_geo/geo_terminal.m earth_elliptic_to_geo/test_terminal.m
git commit -m "feat(earth-geo): pluggable GEO terminal (fixed rendezvous / free-L manifold)

Co-Authored-By: Claude Fable 5 <noreply@anthropic.com>"
```

---

### Task 6: Tangential-propagation seed generator

**Files:**
- Create: `earth_elliptic_to_geo/seed_2body.m`
- Test: `earth_elliptic_to_geo/test_seed.m`

**Interfaces:**
- Consumes: `lt2b_rhs_time`, `kepler_lt_params`.
- Produces: `[sigma, X0, U0, tauf0, info] = seed_2body(par, rv0, opts)` with
  `opts`: `.sbar` (constant throttle), `.tDur` (duration ND; `[]` → propagate until
  GEO energy E=−μ/2a reached), `.N` (segments), `.targetLf` (optional: bisect sbar
  so arrival unwrapped longitude lands within π/2 of it).
  Returns: `sigma` [(N+1)×1] uniform 0→1, `X0` [9×(N+1)] ([r;v;m;t;cScale=1]),
  `U0` [4×(N+1)] ([v̂;sbar]), `tauf0` scalar, `info` (.sbar .Larr .tEnd .mEnd).

- [ ] **Step 1: Write the failing test**

```matlab
% TEST_SEED  Seed exactness: mass linearity, energy growth, small stencil defect,
% longitude-vs-throttle monotonicity, event stop at GEO energy.
p  = kepler_lt_params(10, 1500, 2000);
P0 = 11625/p.LU_km;
[r0, v0] = elements_to_cart(P0, 0.75, 0, 0, 0, pi, p.mu);   % coplanar
rv0 = [r0; v0];
% (a) fixed duration, full throttle
[sg, X0, U0, tauf0, inf1] = seed_2body(p, rv0, struct('sbar',1,'tDur',10,'N',400));
assert(abs(X0(7,end) - (1 - p.Tmax/p.c*10)) < 1e-8, 'mass not linear');
assert(all(diff(X0(8,:)) > 0) && abs(X0(8,end)-10) < 1e-9, 'time state bad');
assert(size(X0,1) == 9 && all(X0(9,:) == 1), 'cScale row missing');
% trapezoid defect of the mapped seed on the solver stencil must be small
dmax = seed_stencil_defect(sg, X0, U0, tauf0, p);
assert(dmax < 1e-2, 'seed defect too big: %.2e', dmax);
% (b) energy event stop: reaches GEO energy
[~, X1, ~, ~, inf2] = seed_2body(p, rv0, struct('sbar',1,'tDur',[],'N',200));
Eend = 0.5*norm(X1(4:6,end))^2 - 1/norm(X1(1:3,end));
assert(abs(Eend - (-0.5)) < 5e-3, 'did not stop at GEO energy');
% (c) winding decreases with throttle
[~,~,~,~, iLo] = seed_2body(p, rv0, struct('sbar',0.5,'tDur',30,'N',200));
[~,~,~,~, iHi] = seed_2body(p, rv0, struct('sbar',0.9,'tDur',30,'N',200));
assert(iLo.Larr > iHi.Larr, 'winding not monotone in sbar');
% (d) bisection hits a target longitude
Lt = inf2.Larr + 6;    % ask for ~1 rev more than the s=1 arrival
[~,~,~,~, iT] = seed_2body(p, rv0, struct('sbar',0.7,'tDur',1.5*inf2.tEnd,'N',200,'targetLf',Lt));
assert(abs(iT.Larr - Lt) < pi/2, 'targetLf bisection missed: %.2f vs %.2f', iT.Larr, Lt);
fprintf('test_seed: ALL PASS\n');

function dmax = seed_stencil_defect(sg, X0, U0, tauf0, p)
% Trapezoid defect of the seed under the solver's tau-stencil (numeric mirror).
N = numel(sg)-1;  dmax = 0;
f = zeros(9, N+1);
for k = 1:N+1
    rn  = norm(X0(1:3,k));
    fd  = lt2b_rhs_time(X0(1:8,k), U0(:,k), p);
    f(:,k) = [X0(9,k) * rn^p.pSund * fd; 0];
end
for k = 1:N
    d = X0(:,k+1) - X0(:,k) - (tauf0*(sg(k+1)-sg(k))/2)*(f(:,k)+f(:,k+1));
    dmax = max(dmax, max(abs(d)));
end
end
```

- [ ] **Step 2: Run to verify it fails.**

- [ ] **Step 3: Write the implementation**

```matlab
function [sigma, X0, U0, tauf0, info] = seed_2body(par, rv0, opts)
% SEED_2BODY  Dynamically-exact tangential-thrust warm start on uniform-tau nodes.
%
% Propagates s = sbar, alpha = vhat from rv0 (ode113, tight tol), computes the
% Sundman clock tau(t) = int dt / r^pSund, and samples the DENSE ode solution at
% N+1 uniform-tau nodes (deval => defect-free to ODE tolerance; no downsampling
% interpolation — the campaign's no-resample lesson). If opts.targetLf is given,
% bisects sbar so the arrival unwrapped equatorial longitude lands within pi/2
% of it (winding is monotone-decreasing in sbar: more thrust climbs sooner ->
% slower angular rate -> less longitude wound).
%
% INPUTS:
%   par  - kepler_lt_params struct
%   rv0  - initial inertial state [6x1]
%   opts - .sbar [scalar], .tDur [scalar ND, or [] => stop at GEO energy -mu/2],
%          .N [segments], .targetLf [optional, rad unwrapped from L0]
%
% OUTPUTS:
%   sigma - [(N+1)x1] uniform 0->1;  X0 [9x(N+1)] = [r;v;m;t;cScale=1];
%   U0    - [4x(N+1)] = [vhat; sbar];  tauf0 - total tau length;
%   info  - .sbar .Larr (arrival unwrapped longitude) .tEnd .mEnd
%
% REFERENCES: [1] DESIGN.md secs 3-4. [2] sundman_minfuel/sundman_seed_map.m.
sbar = opts.sbar;  N = opts.N;
if isfield(opts,'targetLf') && ~isempty(opts.targetLf)
    lo = 0.25;  hi = 1.0;                     % Larr(lo) > Larr(hi)
    for kb = 1:14
        mid  = 0.5*(lo+hi);
        Lm   = propagate(par, rv0, mid, opts.tDur).Larr;
        if abs(Lm - opts.targetLf) < pi/2, sbar = mid; break; end
        if Lm > opts.targetLf, lo = mid; else, hi = mid; end
        sbar = 0.5*(lo+hi);
    end
end
S = propagate(par, rv0, sbar, opts.tDur);
% Sundman map on a fine grid, then uniform-tau node times
tt  = linspace(0, S.tEnd, 20*N+1);
xx  = deval(S.sol, tt);
rr  = sqrt(sum(xx(1:3,:).^2, 1));
tau = cumtrapz(tt, rr.^(-par.pSund));         % dtau = dt / r^p  (cScale=1)
tauf0 = tau(end);
tN  = interp1(tau, tt, linspace(0, tauf0, N+1));
XN  = deval(S.sol, tN);                       % 8 x (N+1), exact to ODE tol
vN  = XN(4:6,:);  vn = max(sqrt(sum(vN.^2,1)), 1e-9);
sigma = linspace(0, 1, N+1).';
X0  = [XN; ones(1, N+1)];
U0  = [vN ./ vn; sbar*ones(1, N+1)];
info = struct('sbar', sbar, 'Larr', S.Larr, 'tEnd', S.tEnd, 'mEnd', XN(7,end));
end

% ---------------------------------------------------------------------------
function S = propagate(par, rv0, sbar, tDur)
% Tangential constant-throttle propagation; empty tDur => stop at GEO energy.
odef = @(t,x) lt2b_rhs_time(x, [x(4:6)/max(norm(x(4:6)),1e-9); sbar], par);
oo = odeset('RelTol',1e-11, 'AbsTol',1e-12);
if isempty(tDur)
    oo  = odeset(oo, 'Events', @geoEnergyEvent);
    sol = ode113(odef, [0 500], [rv0(:); 1; 0], oo);
else
    sol = ode113(odef, [0 tDur], [rv0(:); 1; 0], oo);
end
tf_ = sol.x(end);
tq  = linspace(0, tf_, 4000);
xq  = deval(sol, tq);
Lun = unwrap(atan2(xq(2,:), xq(1,:)));
S = struct('sol', sol, 'tEnd', tf_, 'Larr', Lun(end) + (pi - Lun(1)));
end

function [val, isterm, dir_] = geoEnergyEvent(~, x)
% Stop when two-body energy reaches the GEO value -mu/(2a), a=1, mu=1.
val = 0.5*(x(4)^2+x(5)^2+x(6)^2) - 1/norm(x(1:3)) - (-0.5);
isterm = 1;  dir_ = 1;
end
```

Note the `Larr` convention: unwrapped equatorial longitude re-anchored so it starts
at L₀=π (the paper's initial true longitude) — downstream code compares `Larr`
directly against `targetLf = pi + c_Lf*dL_mt`.

- [ ] **Step 4: Run test to verify it passes.**

- [ ] **Step 5: Commit**

```bash
cd /Users/msc/Desktop/optimal_control
git add earth_elliptic_to_geo/seed_2body.m earth_elliptic_to_geo/test_seed.m
git commit -m "feat(earth-geo): tangential-propagation Sundman seed (event stop, longitude bisection)

Co-Authored-By: Claude Fable 5 <noreply@anthropic.com>"
```

---

### Task 7: Solver core `casadi_lt_2body`

**Files:**
- Create: `earth_elliptic_to_geo/casadi_lt_2body.m`
- Test: `earth_elliptic_to_geo/test_solver_smoke.m` (construction + 5-iteration smoke; real solves gate in Task 8)

**Interfaces:**
- Consumes: `geo_terminal` term structs, seeds from `seed_2body`.
- Produces: `out = casadi_lt_2body(sigma, X0, U0, tauf0, term, opts)` with
  `opts`: `.par` (required), `.mode` `'mintime'|'fixedtf'`, `.eps` (fixedtf),
  `.tfTarget` (fixedtf), `.rv0` [6×1], `.maxIter` (default 1500), `.warmTight`
  (default false), `.printLevel` (default 0).
  `out`: `.X` [9×(N+1)] `.U` [4×(N+1)] `.tauf0 .success .ipoptStatus .maxDefect
  .maxUnit .termErr .mf .m_f_kg .dV_kms .tf .switches .edge .lamDef` [9×N]
  `.primerAlignDeg .lamMassEnd`. Every later task consumes these exact fields.

- [ ] **Step 1: Write the failing smoke test**

```matlab
% TEST_SOLVER_SMOKE  Construction + short-iteration smoke (no convergence gate).
p  = kepler_lt_params(10, 1500, 2000);
P0 = 11625/p.LU_km;
[r0, v0] = elements_to_cart(P0, 0.75, 0, 0, 0, pi, p.mu);
rv0 = [r0; v0];
[sg, X0, U0, tauf0] = seed_2body(p, rv0, struct('sbar',1,'tDur',5,'N',80));
term = geo_terminal('manifold', p, []);
out = casadi_lt_2body(sg, X0, U0, tauf0, term, struct('par',p,'mode','mintime', ...
        'rv0',rv0,'maxIter',5,'printLevel',0));
assert(isstruct(out) && isfield(out,'maxDefect') && isfield(out,'lamDef'));
assert(size(out.X,1) == 9 && size(out.U,1) == 4);
% fixedtf construction path too
out2 = casadi_lt_2body(sg, X0, U0, tauf0, geo_terminal('fixed',p,pi+8), ...
        struct('par',p,'mode','fixedtf','eps',1,'tfTarget',6,'rv0',rv0, ...
               'maxIter',5,'printLevel',0));
assert(isfield(out2,'tf') && ~out2.success);   % 5 iters won't converge; must not error
fprintf('test_solver_smoke: ALL PASS\n');
```

- [ ] **Step 2: Run to verify it fails.**

- [ ] **Step 3: Write the implementation**

```matlab
function out = casadi_lt_2body(sigma, X0, U0, tauf0, term, opts)
% CASADI_LT_2BODY  Sundman-collocated 2-body low-thrust NLP (CasADi+IPOPT).
%
% Trapezoidal collocation in the Sundman variable tau with a cScale slack state:
%     dt/dtau = cScale * kappa(r),   kappa = ||r||^pSund,   dcScale/dtau = 0,
% tau_f held FIXED (= tauf0) so the KKT stays banded (Betts' sparse free-time
% trick; a free scalar tau_f makes a dense KKT column). State per node
% x = [r(3); v(3); m; t; cScale] (9), control u = [alpha(3); s] (4), cone-
% eliminated thrust = s*Tmax*alpha/m with ||alpha||=1, s in [0,1].
%
% Modes:
%   'mintime' - s == 1 (all-burn restriction; optimal for this transfer),
%               objective min t(tau_f); t_f found via cScale.
%   'fixedtf' - constraint t(tau_f) = opts.tfTarget; Bertrand-Epenoy objective
%               J(eps) = Int[s]dt - eps*Int[s(1-s)]dt   (dt = cScale*kappa dtau)
%               eps=1 energy (smooth), eps=0 fuel (bang-bang).
%
% Terminal (term from GEO_TERMINAL): 'fixed' pins X(1:6,end) = term.rvf;
% 'manifold' poses the 5 insertion constraints symbolically.
%
% INPUTS:  sigma [(N+1)x1] 0->1;  X0 [8|9 x N+1] warm start (cScale row appended
%          as 1s if absent);  U0 [4xN+1];  tauf0 [scalar];  term [struct];
%          opts .par .mode .eps .tfTarget .rv0 .maxIter .warmTight .printLevel
% OUTPUTS: out - see header table in PLAN.md Task 7 (X,U,success,maxDefect,
%          lamDef, primerAlignDeg, m_f_kg, dV_kms, switches, edge, ...)
%
% REFERENCES:
%   [1] GTO_tulip/sundman_minfuel/casadi_minfuel_sundman.m (parent).
%   [2] GTO_tulip/elfo/casadi_energy_freetf.m (cScale pattern).
%   [3] DESIGN.md secs 2-4.
cp = getenv('CASADI_PATH');
if isempty(cp), cp = fullfile(getenv('HOME'), 'casadi-3.7.0'); end
addpath(cp);
par = opts.par;
d = @(f,v) getdef(opts, f, v);
mode      = d('mode', 'fixedtf');
epsv      = d('eps', 0);
tfTarget  = d('tfTarget', []);
maxIter   = d('maxIter', 1500);
warmTight = d('warmTight', false);
printLvl  = d('printLevel', 0);

N    = numel(sigma) - 1;
dtau = diff(sigma(:)).' * tauf0;                    % [1xN]
if size(X0,1) == 8, X0 = [X0; ones(1, N+1)]; end

opti = casadi.Opti();
X = opti.variable(9, N+1);
U = opti.variable(4, N+1);
m = X(7,:);  t = X(8,:);  cS = X(9,:);  al = U(1:3,:);  s = U(4,:);

% node dynamics f{k} = dX/dtau and the clock row kapAll (for the objective)
f = cell(1, N+1);  kapCell = cell(1, N+1);
for k = 1:N+1
    rk  = X(1:3,k);
    rn2 = rk(1)^2 + rk(2)^2 + rk(3)^2 + 1e-12;
    kap = rn2^(par.pSund/2);
    fd  = lt2b_rhs_time(X(1:8,k), U(:,k), par);     % d/dt of [r v m t]
    f{k} = [cS(k)*kap*fd; 0];                       % d/dtau; cScale constant
    kapCell{k} = kap;
end
kapAll = [kapCell{:}];

% collocation defects (KEEP HANDLES for the duals)
conDef = cell(1, N);
for k = 1:N
    conDef{k} = X(:,k+1) - X(:,k) - (dtau(k)/2)*(f{k} + f{k+1}) == 0;
    opti.subject_to(conDef{k});
end

% control cone + throttle (NEVER chain a<=x<=b -- MATLAB gotcha)
for k = 1:N+1
    opti.subject_to(al(1,k)^2 + al(2,k)^2 + al(3,k)^2 == 1);
end
if strcmp(mode, 'mintime')
    opti.subject_to(s == 1);
else
    opti.subject_to(s >= 0);  opti.subject_to(s <= 1);
end

% generous boxes (review lesson: bounds only block divergence)
opti.subject_to(X(1:3,:) >= -5);   opti.subject_to(X(1:3,:) <= 5);
opti.subject_to(X(4:6,:) >= -8);   opti.subject_to(X(4:6,:) <= 8);
opti.subject_to(m >= 0.3);         opti.subject_to(m <= 1.001);
opti.subject_to(t >= 0);           opti.subject_to(t <= 300);
opti.subject_to(cS >= 0.05);       opti.subject_to(cS <= 20);
opti.subject_to(al >= -1.01);      opti.subject_to(al <= 1.01);

% boundary conditions
opti.subject_to(X(1:6,1) == opts.rv0(:));
opti.subject_to(m(1) == 1);
opti.subject_to(t(1) == 0);
switch term.type
    case 'fixed'
        opti.subject_to(X(1:6,end) == term.rvf(:));
    case 'manifold'
        re = X(1:3,end);  ve = X(4:6,end);  a = term.aGeo;
        opti.subject_to(re(3) == 0);
        opti.subject_to(ve(3) == 0);
        opti.subject_to(re(1)^2 + re(2)^2 + re(3)^2 == a^2);
        opti.subject_to(ve(1)^2 + ve(2)^2 + ve(3)^2 == par.mu/a);
        opti.subject_to(re(1)*ve(1) + re(2)*ve(2) + re(3)*ve(3) == 0);
end

% objective + t_f handling
if strcmp(mode, 'mintime')
    opti.minimize(t(end));
else
    assert(~isempty(tfTarget), 'fixedtf mode requires opts.tfTarget');
    opti.subject_to(t(end) == tfTarget);
    w = cS .* kapAll .* (s - epsv*(s.*(1 - s)));    % homotopy integrand * clock
    opti.minimize(sum((dtau/2) .* (w(1:N) + w(2:N+1))));
end

% warm start + IPOPT
opti.set_initial(X, X0);
opti.set_initial(U, U0);
ip = struct('max_iter', maxIter, 'tol', 1e-9, 'constr_viol_tol', 1e-10, ...
            'print_level', printLvl, 'mu_strategy', 'adaptive', ...
            'linear_solver', 'mumps');
if warmTight
    ip.mu_strategy = 'monotone';  ip.mu_init = 1e-4;
    ip.warm_start_init_point = 'yes';
    ip.warm_start_bound_push = 1e-9;  ip.warm_start_mult_bound_push = 1e-9;
end
opti.solver('ipopt', struct('print_time', printLvl > 0), ip);
success = true;
try
    sol = opti.solve();
catch
    sol = opti.debug;  success = false;
end
st = opti.stats();
status = st.return_status;
success = success && any(strcmp(status, {'Solve_Succeeded', 'Solved_To_Acceptable_Level'}));

% extraction + numeric re-check of the defects
Xs = sol.value(X);  Us = sol.value(U);
dmax = 0;  fn = zeros(9, N+1);
for k = 1:N+1
    rn = norm(Xs(1:3,k));
    fn(:,k) = [Xs(9,k) * rn^par.pSund * lt2b_rhs_time(Xs(1:8,k), Us(:,k), par); 0];
end
for k = 1:N
    dk = Xs(:,k+1) - Xs(:,k) - (dtau(k)/2)*(fn(:,k) + fn(:,k+1));
    dmax = max(dmax, max(abs(dk)));
end
lamDef = nan(9, N);
try
    for k = 1:N, lamDef(:,k) = sol.value(opti.dual(conDef{k})); end
catch
end
ss = Us(4,:);
burn = ss > 0.5;
% primer alignment on burn nodes (global costate sign resolved by best fit)
lamV = lamDef(4:6, :);
angs = @(sgn) mean(arrayfun(@(k) real(acosd(max(-1,min(1, ...
        dot(Us(1:3,k), sgn*(-lamV(:,min(k,N)))) / max(norm(lamV(:,min(k,N))),1e-30))))), ...
        find(burn)));
primer = min(angs(1), angs(-1));
switch term.type
    case 'fixed',    termErr = norm(Xs(1:6,end) - term.rvf(:));
    case 'manifold', termErr = max(abs(term.resid(Xs(1:6,end))));
end
mf = Xs(7,end);
out = struct('X', Xs, 'U', Us, 'tauf0', tauf0, 'success', success, ...
    'ipoptStatus', status, 'maxDefect', dmax, ...
    'maxUnit', max(abs(sum(Us(1:3,:).^2,1) - 1)), 'termErr', termErr, ...
    'mf', mf, 'm_f_kg', par.m0kg*mf, 'dV_kms', par.c*log(1/mf)*par.VU_kms, ...
    'tf', Xs(8,end), 'switches', sum(abs(diff(burn))), ...
    'edge', mean(ss > 0.95 | ss < 0.05), 'lamDef', lamDef, ...
    'primerAlignDeg', primer, 'lamMassEnd', lamDef(7,end));
end

% ---------------------------------------------------------------------------
function v = getdef(s, f, dflt)
% GETDEF  Optional-field default (mirrors campaign helper).
if isfield(s, f) && ~isempty(s.(f)), v = s.(f); else, v = dflt; end
end
```

- [ ] **Step 4: Run smoke test to verify it passes** (constructs, runs 5 iterations,
returns the full struct without erroring on the non-converged path).

- [ ] **Step 5: Commit**

```bash
cd /Users/msc/Desktop/optimal_control
git add earth_elliptic_to_geo/casadi_lt_2body.m earth_elliptic_to_geo/test_solver_smoke.m
git commit -m "feat(earth-geo): mode-switched Sundman collocation solver core (CasADi+IPOPT)

Co-Authored-By: Claude Fable 5 <noreply@anthropic.com>"
```

---

### Task 8: Min-time anchors (coplanar + 3D) — gate G1

**Files:**
- Create: `earth_elliptic_to_geo/run_mintime.m`
- Output: `results/mintime_T10_i0.mat`, `results/mintime_T10_i7.mat`

**Interfaces:**
- Produces: `res = run_mintime(thrustN, hx0, N)` — solves manifold min-time at
  thrust `thrustN`, inclination via `hx0` (0 or 0.0612), mesh `N`. Saves/returns
  `res`: `.out` (solver struct), `.tfmin` (ND), `.tfmin_h` (hours), `.dL_mt`
  (unwrapped longitude span, rad), `.revs`. Caches to
  `results/mintime_T<thrustN*10>_i<0|7>.mat`; loads instead of re-solving if present.

- [ ] **Step 1: Write run_mintime**

```matlab
function res = run_mintime(thrustN, hx0, N)
% RUN_MINTIME  Free-L (manifold) min-time anchor at one thrust level.
%
% Seeds with a full-throttle tangential propagation stopped at GEO energy, then
% solves casadi_lt_2body in 'mintime' mode against the insertion manifold.
% t_f,min sets every c_tf scale downstream (paper's TfMin is free-longitude).
%
% INPUTS:  thrustN - max thrust [N];  hx0 - initial hx (0 coplanar | 0.0612);
%          N - mesh segments (default 600)
% OUTPUTS: res - .out .tfmin .tfmin_h .dL_mt .revs (also saved/cached in results/)
%
% REFERENCES: [1] DESIGN.md sec 4 step 1.  [2] PLAN.md Task 8.
if nargin < 3, N = 600; end
here = fileparts(mfilename('fullpath'));
resDir = fullfile(here, 'results');
if ~exist(resDir, 'dir'), mkdir(resDir); end
tag = sprintf('mintime_T%d_i%d', round(10*thrustN), round(hx0 > 0)*7);
fn  = fullfile(resDir, [tag '.mat']);
if isfile(fn), S = load(fn); res = S.res; fprintf('cached %s\n', fn); return; end

p  = kepler_lt_params(thrustN, 1500, 2000);
P0 = 11625/p.LU_km;
[r0, v0] = elements_to_cart(P0, 0.75, 0, hx0, 0, pi, p.mu);
rv0 = [r0; v0];
[sg, X0, U0, tauf0, sinfo] = seed_2body(p, rv0, struct('sbar',1,'tDur',[],'N',N));
term = geo_terminal('manifold', p, []);
out = casadi_lt_2body(sg, X0, U0, tauf0, term, struct('par',p,'mode','mintime', ...
        'rv0',rv0,'maxIter',3000,'printLevel',3));
% CONTINGENCY (documented, use only if the manifold solve fails): first solve
% 'mintime' with term = geo_terminal('fixed', p, sinfo.Larr) — zero terminal gap
% rendezvous at the seed's own arrival longitude — then re-solve 'manifold'
% warm-started from that solution with warmTight=true.
Lun = unwrap(atan2(out.X(2,:), out.X(1,:)));
res = struct('out', out, 'tfmin', out.tf, 'tfmin_h', out.tf*p.TU_s/3600, ...
             'dL_mt', Lun(end)-Lun(1), 'revs', (Lun(end)-Lun(1))/(2*pi), ...
             'thrustN', thrustN, 'hx0', hx0, 'N', N, 'seedInfo', sinfo);
save(fn, 'res');
fprintf('MINTIME T=%g N: tf=%.4f ND = %.1f h, dL=%.1f rad (%.2f revs), defect %.2e, %s\n', ...
        thrustN, res.tfmin, res.tfmin_h, res.dL_mt, res.revs, out.maxDefect, out.ipoptStatus);
end
```

- [ ] **Step 2: Run the coplanar anchor**

Run: `/Applications/MATLAB_R2025b.app/bin/matlab -batch "cd '/Users/msc/Desktop/optimal_control/earth_elliptic_to_geo'; res = run_mintime(10, 0, 600); disp(res)"`
Gates (coplanar): `out.success` true; `maxDefect < 1e-8`; `termErr < 1e-8`;
`tfmin` ∈ [14, 30] ND (paper ≈ 22.2 ND = 84.7 h); `revs` ∈ [3, 5.5];
prograde h_z > 0 at the terminal.

- [ ] **Step 3: Run the 3D anchor**

Run: same with `run_mintime(10, 0.0612, 600)`.
Gates: as above, plus `tfmin(3D) ≥ tfmin(coplanar)` (plane change costs time, ≲10%).
Record `dL_mt` — expect ≈ 26 ± 8 rad (paper law R1: 26.4 at 10 N).

- [ ] **Step 4: If a gate fails,** apply the contingency in the code comment (fixed-at-
seed-arrival rendezvous first, then free the terminal warm). If it still fails, STOP
and report — do not proceed to Task 9 on a broken anchor.

- [ ] **Step 5: Commit**

```bash
cd /Users/msc/Desktop/optimal_control
git add earth_elliptic_to_geo/run_mintime.m
git commit -m "feat(earth-geo): free-longitude min-time anchor (G1: ~85h, ~4.2 revs at 10N)

Co-Authored-By: Claude Fable 5 <noreply@anthropic.com>"
```

---

### Task 9: Homotopy driver + energy-stage gate G2 (coplanar)

**Files:**
- Create: `earth_elliptic_to_geo/homotopy_2body.m`
- Test: `earth_elliptic_to_geo/test_energy_stage.m` (a real coplanar ε=1 solve)

**Interfaces:**
- Consumes: `casadi_lt_2body`, seeds, `run_mintime` cache.
- Produces: `[best, tbl] = homotopy_2body(sigma, X0, U0, tauf0, term, tf, opts)` —
  `opts`: `.par .rv0 .maxIter` (default 1500), `.sched` (default
  `[1 0.6 0.35 0.2 0.12 0.07 0.04 0.025 0.015 0.008 0.004 0.002 0.001 0]`).
  Returns `best` = last tight-converged solver `out` (plus `.certified` logical,
  `.epsReached`), `tbl` [K×5] = [eps, maxDefect, switches, edge, m_f_kg].

- [ ] **Step 1: Write the implementation**

```matlab
function [best, tbl] = homotopy_2body(sigma, X0, U0, tauf0, term, tf, opts)
% HOMOTOPY_2BODY  Guarded energy->fuel sweep at fixed t_f (eps: 1 -> 0).
%
% First step (eps=1) runs LOOSE (genuine move from a propagated seed); every
% later step warm-starts TIGHT from the previous converged iterate. GUARD: a
% step that fails to converge tight never advances the warm start and never
% overwrites best (campaign lesson: a loose iterate must not poison the chain).
%
% INPUTS:  sigma/X0/U0/tauf0 - seed (seed_2body layout);  term - geo_terminal
%          struct;  tf - fixed transfer time [ND];  opts - .par .rv0 .maxIter .sched
% OUTPUTS: best - last tight solver out + .certified .epsReached;  tbl [Kx5]
%
% REFERENCES: [1] sundman_minfuel/sundman_homotopy.m (pattern). [2] DESIGN.md sec 4.
d = @(f,v) getdef2(opts, f, v);
sched   = d('sched', [1 0.6 0.35 0.2 0.12 0.07 0.04 0.025 0.015 0.008 0.004 0.002 0.001 0]);
maxIter = d('maxIter', 1500);
Xk = X0;  Uk = U0;  best = [];  tbl = zeros(numel(sched), 5);
for ke = 1:numel(sched)
    e = sched(ke);
    o = casadi_lt_2body(sigma, Xk, Uk, tauf0, term, struct('par',opts.par, ...
        'mode','fixedtf', 'eps',e, 'tfTarget',tf, 'rv0',opts.rv0, ...
        'maxIter',maxIter, 'warmTight', ke > 1, 'printLevel',0));
    ok = o.success && o.maxDefect < 1e-8;
    tbl(ke,:) = [e, o.maxDefect, o.switches, o.edge, o.m_f_kg];
    fprintf('  eps=%6.4f ok=%d defect=%.2e sw=%3d edge=%5.1f%% mf=%.2f kg\n', ...
            e, ok, o.maxDefect, o.switches, 100*o.edge, o.m_f_kg);
    if ok
        Xk = o.X;  Uk = o.U;  best = o;  best.epsReached = e;
    end
end
if isempty(best)
    best = o;  best.epsReached = NaN;  best.certified = false;
else
    best.certified = (best.epsReached == 0);
end
end

function v = getdef2(s, f, dflt)
if isfield(s, f) && ~isempty(s.(f)), v = s.(f); else, v = dflt; end
end
```

- [ ] **Step 2: Write the energy-stage test (gate G2)**

```matlab
% TEST_ENERGY_STAGE  Coplanar eps=1 (energy) solve at tf = 1.5*tfmin, fixed L_f.
p   = kepler_lt_params(10, 1500, 2000);
P0  = 11625/p.LU_km;
[r0, v0] = elements_to_cart(P0, 0.75, 0, 0, 0, pi, p.mu);
rv0 = [r0; v0];
mt  = run_mintime(10, 0, 600);
tf  = 1.5 * mt.tfmin;
Lf  = pi + (1.12*1.5 + 0.09) * mt.dL_mt;      % paper law R2 (c_Lf ~ 1.77)
[sg, X0, U0, tauf0, si] = seed_2body(p, rv0, ...
      struct('sbar', 1/1.5, 'tDur', tf, 'N', 600, 'targetLf', Lf));
term = geo_terminal('fixed', p, Lf);
o = casadi_lt_2body(sg, X0, U0, tauf0, term, struct('par',p,'mode','fixedtf', ...
      'eps',1,'tfTarget',tf,'rv0',rv0,'maxIter',3000,'printLevel',3));
assert(o.success, 'energy solve failed: %s', o.ipoptStatus);
assert(o.maxDefect < 1e-8, 'energy defect %.2e', o.maxDefect);
assert(abs(o.tf - tf) < 1e-6, 'tf pin violated');
save(fullfile('results','energy_M0_coplanar.mat'), 'o', 'sg', 'tauf0', 'tf', 'Lf', 'rv0', 'si');
fprintf('test_energy_stage: ALL PASS (mf=%.2f kg, edge=%.1f%%)\n', o.m_f_kg, 100*o.edge);
```

- [ ] **Step 3: Run it.** Expected: `ALL PASS`. CONTINGENCY if the ε=1 solve stalls in
restoration: re-run with `term = geo_terminal('fixed', p, si.Larr)` — i.e. pin the
terminal at the seed's own arrival longitude (zero topology gap) — and if that
converges, walk L_f to the law value in 2–3 fixed-longitude continuation steps,
each warm-started (discrete continuation in L_f, the paper's own c_Lf device).

- [ ] **Step 4: Commit**

```bash
cd /Users/msc/Desktop/optimal_control
git add earth_elliptic_to_geo/homotopy_2body.m earth_elliptic_to_geo/test_energy_stage.m
git commit -m "feat(earth-geo): guarded energy->fuel homotopy driver; energy stage converges (G2)

Co-Authored-By: Claude Fable 5 <noreply@anthropic.com>"
```

---

### Task 10: Single-case driver + milestone M0 (coplanar bang-bang)

**Files:**
- Create: `earth_elliptic_to_geo/run_transfer.m`
- Output: `results/M0_coplanar.mat`

**Interfaces:**
- Produces: `res = run_transfer(cfg)` — cfg fields: `.thrustN .ctf`
  `.hx0` (0|0.0612) `.term` ('fixed'|'manifold') `.N` (default 600)
  `.tag` (results filename stem) `.seedMat` (optional: warm-start the homotopy
  from a prior result's `res.fuel` instead of building a seed; uses the light
  schedule `[0.05 0.02 0.008 0.003 0.001 0]`) `.ispS` (default 2000).
  Returns/saves `res`: `.cfg .mintime .tf .Lf .fuel` (best solver out) `.tbl`
  `.report` (struct: `.revs .switches .m_f_kg .dV_kms .edge .apoBurnRatio`).

- [ ] **Step 1: Write run_transfer**

```matlab
function res = run_transfer(cfg)
% RUN_TRANSFER  One full pipeline: mintime anchor -> seed -> homotopy -> report.
%
% Stages: (1) cached free-L min-time at cfg.thrustN -> tfmin, dL_mt;
% (2) tf = ctf*tfmin; L_f = pi + (1.12*ctf+0.09)*dL_mt (paper law R2);
% (3) seed: tangential sbar=1/ctf bisected on L_f ('fixed' term) or plain
%     ('manifold'); or warm-start from cfg.seedMat with the light schedule;
% (4) homotopy eps 1->0; (5) structure report.
%
% INPUTS:  cfg - see PLAN.md Task 10 interface block
% OUTPUTS: res - .cfg .mintime .tf .Lf .fuel .tbl .report (saved to
%          results/<tag>.mat)
%
% REFERENCES: [1] DESIGN.md secs 4-5.
here = fileparts(mfilename('fullpath'));
resDir = fullfile(here, 'results');
if ~exist(resDir, 'dir'), mkdir(resDir); end
d = @(f,v) getdef3(cfg, f, v);
N = d('N', 600);  ispS = d('ispS', 2000);  seedMat = d('seedMat', '');

p  = kepler_lt_params(cfg.thrustN, 1500, ispS);
P0 = 11625/p.LU_km;
[r0, v0] = elements_to_cart(P0, 0.75, 0, cfg.hx0, 0, pi, p.mu);
rv0 = [r0; v0];
mt  = run_mintime(cfg.thrustN, cfg.hx0, N);
tf  = cfg.ctf * mt.tfmin;
Lf  = pi + (1.12*cfg.ctf + 0.09) * mt.dL_mt;
switch cfg.term
    case 'fixed',    term = geo_terminal('fixed', p, Lf);
    case 'manifold', term = geo_terminal('manifold', p, []);
end
ho = struct('par', p, 'rv0', rv0, 'maxIter', 1500);
if ~isempty(seedMat)                       % neighbor-style warm start
    S = load(seedMat);
    sg = S.res.sg;  tauf0 = S.res.fuel.tauf0;
    Xk = S.res.fuel.X;  Uk = S.res.fuel.U;
    Xk(8,:) = Xk(8,:) * (tf / Xk(8,end));  % rescale carried time if tf differs
    ho.sched = [0.05 0.02 0.008 0.003 0.001 0];
else
    so = struct('sbar', 1/cfg.ctf, 'tDur', tf, 'N', N);
    if strcmp(cfg.term, 'fixed'), so.targetLf = Lf; end
    [sg, Xk, Uk, tauf0] = seed_2body(p, rv0, so);
end
fprintf('RUN_TRANSFER %s: T=%g N, ctf=%.2f, tf=%.3f ND (%.1f h), Lf=%.2f rad\n', ...
        cfg.tag, cfg.thrustN, cfg.ctf, tf, tf*p.TU_s/3600, Lf);
[best, tbl] = homotopy_2body(sg, Xk, Uk, tauf0, term, tf, ho);

% structure report
Lun  = unwrap(atan2(best.X(2,:), best.X(1,:)));
revs = (Lun(end) - Lun(1)) / (2*pi);
rr   = sqrt(sum(best.X(1:3,:).^2, 1));
ss   = best.U(4,:);
nEarly = round(0.8 * numel(ss));           % exclude near-circular endgame
bMask  = ss(1:nEarly) > 0.5;
apoBurnRatio = median(rr(bMask)) / median(rr(~bMask));
report = struct('revs', revs, 'switches', best.switches, 'm_f_kg', best.m_f_kg, ...
    'dV_kms', best.dV_kms, 'edge', best.edge, 'apoBurnRatio', apoBurnRatio, ...
    'defect', best.maxDefect, 'certified', best.certified);
res = struct('cfg', cfg, 'mintime', mt, 'tf', tf, 'Lf', Lf, 'fuel', best, ...
             'tbl', tbl, 'report', report, 'sg', sg, 'rv0', rv0);
save(fullfile(resDir, [cfg.tag '.mat']), 'res');
fprintf(['DONE %s: certified=%d revs=%.2f sw=%d edge=%.1f%% mf=%.2f kg ' ...
         'dV=%.3f km/s apoBurn=%.2f\n'], cfg.tag, report.certified, revs, ...
         best.switches, 100*best.edge, best.m_f_kg, best.dV_kms, apoBurnRatio);
end

function v = getdef3(s, f, dflt)
if isfield(s, f) && ~isempty(s.(f)), v = s.(f); else, v = dflt; end
end
```

- [ ] **Step 2: Run M0**

Run: `/Applications/MATLAB_R2025b.app/bin/matlab -batch "cd '/Users/msc/Desktop/optimal_control/earth_elliptic_to_geo'; res = run_transfer(struct('thrustN',10,'ctf',1.5,'hx0',0,'term','fixed','tag','M0_coplanar'));"`
(Reuses the Task 9 energy checkpoint chain internally via the fresh sweep — the
ε=1 leg repeats; ~10–30 min total.)

**M0 gates:** `certified` true (ε reached 0); `defect < 1e-8`; `edge ≥ 0.95`;
`switches` ∈ [4, 40]; `m_f_kg > mintime m_f` (coasting must save fuel);
`abs(revs − 7.5) < 1.5` (law-prescribed L_f ⇒ ~7.5 revs even coplanar).

- [ ] **Step 3: Commit**

```bash
cd /Users/msc/Desktop/optimal_control
git add earth_elliptic_to_geo/run_transfer.m
git commit -m "feat(earth-geo): full pipeline driver; M0 coplanar bang-bang converges

Co-Authored-By: Claude Fable 5 <noreply@anthropic.com>"
```

---

### Task 11: Milestone M1 — 3D, fixed L_f (paper structure)

**Files:**
- Output: `results/M1_3d_fixedLf.mat` (no new code expected)

- [ ] **Step 1: Run M1**

Run: `/Applications/MATLAB_R2025b.app/bin/matlab -batch "cd '/Users/msc/Desktop/optimal_control/earth_elliptic_to_geo'; res = run_transfer(struct('thrustN',10,'ctf',1.5,'hx0',0.0612,'term','fixed','tag','M1_3d_fixedLf'));"`

- [ ] **Step 2: Check the M1 gates (paper Table 3 / Fig 22 at 10 N, c_tf=1.5)**

- `certified` true, `defect < 1e-8`, `edge ≥ 0.95`
- `revs` ∈ [6.9, 8.1] (paper: 7.5, held ~by construction via L_f)
- `switches` ∈ [12, 24] (paper: 18; trapezoid mesh chatter tolerated)
- `apoBurnRatio > 1.5` (burns at apogee — the paper's signature structure; the
  mirror of our tulip perigee-burns)
- inclination closed: `abs(acosd(hz/|h|)) < 0.05 deg` at the terminal node

If switches fall outside the band but everything else passes, re-run once at
`'N', 1200` and compare — switch count is the mesh-sensitive number (campaign
lesson: never publish switch counts without a refinement check).

- [ ] **Step 3: Commit the result note**

```bash
cd /Users/msc/Desktop/optimal_control
git add -A earth_elliptic_to_geo/*.m
git commit --allow-empty -m "milestone(earth-geo): M1 3D fixed-Lf reproduces paper structure (revs/switches/apogee-burns)

Co-Authored-By: Claude Fable 5 <noreply@anthropic.com>"
```

---

### Task 12: PMP verification module

**Files:**
- Create: `earth_elliptic_to_geo/verify_pmp_2body.m`
- Test: run against `results/M1_3d_fixedLf.mat`

**Interfaces:**
- Consumes: solver `out` struct (`.U .X .lamDef` [9×N]) + `par`.
- Produces: `ver = verify_pmp_2body(out, par)` → `.primerMeanDeg .beta
  .betaSpreadPct .burnSignPct .coastSignPct .lamMendRel .pass` (logical).

- [ ] **Step 1: Write the implementation**

```matlab
function ver = verify_pmp_2body(out, par)
% VERIFY_PMP_2BODY  First-order PMP consistency from the NLP's KKT duals.
%
% Costates = defect-constraint multipliers out.lamDef (discrete adjoint up to a
% positive scale and a global sign). Checks:
%   (1) primer direction: thrust || -lam_v on burn arcs (scale-invariant);
%   (2) switching sign law via the empirical-beta route: the min-fuel switching
%       function S = 1 - beta*W with W = (c/m)*||lam_v|| + lam_m; beta pinned
%       as the ROBUST MEDIAN of 1/W over switch-adjacent intervals (absorbs the
%       covector scale); require S<0 on >=98%% of burn nodes, S>0 on >=98%% of
%       coast nodes;
%   (3) relative transversality |lam_m(end)| / max|lam_m| <= 1e-3 (free mass).
%
% INPUTS:  out - casadi_lt_2body result (.U .lamDef);  par - params struct
% OUTPUTS: ver - .primerMeanDeg .beta .betaSpreadPct .burnSignPct .coastSignPct
%                .lamMendRel .pass
%
% REFERENCES:
%   [1] GTO_tulip/sundman_minfuel/verify_tf_front.m (empirical-beta).
%   [2] HONEST_EVALUATION_DV_TF_FRONT.md (robust-beta lesson; relative gate).
N    = size(out.lamDef, 2);
ss   = out.U(4, 1:N);
burn = ss > 0.5;
lamV = out.lamDef(4:6, :);  lamM = out.lamDef(7, :);
mAvg = 0.5*(out.X(7,1:N) + out.X(7,2:N+1));
% global sign: primer alignment must be < 90 deg on burns
lVn = sqrt(sum(lamV.^2, 1));
cosb = zeros(1, N);
for k = 1:N
    cosb(k) = dot(out.U(1:3,k), -lamV(:,k)) / max(lVn(k), 1e-30);
end
if mean(cosb(burn)) < 0, lamV = -lamV; lamM = -lamM; cosb = -cosb; end
ver.primerMeanDeg = mean(real(acosd(max(-1, min(1, cosb(burn))))));
W = (par.c ./ mAvg) .* sqrt(sum(lamV.^2,1)) + lamM;
sw = find(diff(burn) ~= 0);                     % switch-adjacent intervals
bcand = 1 ./ W([sw, min(sw+1, N)]);
ver.beta = median(bcand);
ver.betaSpreadPct = 100 * median(abs(bcand - ver.beta)) / abs(ver.beta);  % MAD, no toolbox
S = 1 - ver.beta * W;
ver.burnSignPct  = 100 * mean(S(burn)  < 0);
ver.coastSignPct = 100 * mean(S(~burn) > 0);
ver.lamMendRel   = abs(lamM(end)) / max(abs(lamM));
ver.pass = ver.primerMeanDeg < 1.0 && ver.burnSignPct >= 98 && ...
           ver.coastSignPct >= 98 && ver.lamMendRel <= 1e-3;
fprintf(['verify_pmp_2body: primer %.3f deg | burn %.1f%% coast %.1f%% ' ...
         '(beta spread %.1f%%) | lamM rel %.1e | pass=%d\n'], ver.primerMeanDeg, ...
         ver.burnSignPct, ver.coastSignPct, ver.betaSpreadPct, ver.lamMendRel, ver.pass);
end
```

- [ ] **Step 2: Run on M1**

Run: `/Applications/MATLAB_R2025b.app/bin/matlab -batch "cd '/Users/msc/Desktop/optimal_control/earth_elliptic_to_geo'; S=load('results/M1_3d_fixedLf.mat'); p=kepler_lt_params(10,1500,2000); ver=verify_pmp_2body(S.res.fuel, p); assert(ver.pass)"`
Expected: primer < 1° (campaign benchmark ~0.06–0.2°), burn/coast sign ≥ 98%, pass=1.
Append `ver` to the M1 `.mat` (`save(...,'-append')`).

- [ ] **Step 3: Commit**

```bash
cd /Users/msc/Desktop/optimal_control
git add earth_elliptic_to_geo/verify_pmp_2body.m
git commit -m "feat(earth-geo): KKT-dual PMP verifier (primer + robust-beta switching law + transversality)

Co-Authored-By: Claude Fable 5 <noreply@anthropic.com>"
```

---

### Task 13: Milestone M2 — free-L_f manifold + headline mass match

**Files:**
- Output: `results/M2_manifold.mat`

- [ ] **Step 1 (Isp pin attempt):** Try `WebSearch`/`WebFetch` for the Caillau &
Noailles 2001 (ESAIM COCV 6) benchmark Isp (DESIGN.md open item 1). If found and
≠ 2000 s, pass it via `cfg.ispS` below and note it in the results. If not
reachable, proceed with 2000 s (the m_f gate below carries the sanity band).

- [ ] **Step 2: Run M2, seeded from M1 (neighbor warm start, light schedule)**

Run: `/Applications/MATLAB_R2025b.app/bin/matlab -batch "cd '/Users/msc/Desktop/optimal_control/earth_elliptic_to_geo'; res = run_transfer(struct('thrustN',10,'ctf',1.5,'hx0',0.0612,'term','manifold','tag','M2_manifold','seedMat','results/M1_3d_fixedLf.mat'));"`

- [ ] **Step 3: Check the M2 gates (the headline reproduction)**

- `certified` true, `defect < 1e-8`, manifold `termErr < 1e-8`
- **`m_f_kg` ∈ [1355, 1385]** (paper Fig 23: ≈1370–1375 at 10 N, c_tf=1.5)
- `revs` ∈ [6.9, 8.1]; `switches` ∈ [12, 24]; `apoBurnRatio > 1.5`
- `m_f_kg(M2) ≥ m_f_kg(M1) − 0.5` (freeing the longitude cannot cost fuel)
- run `verify_pmp_2body` → pass (M2 note: with the manifold the free-phase
  transversality is implicitly enforced by the NLP; record `ver` alongside)

If m_f lands outside the band with everything else green, calibrate Isp within
[1800, 2400] s to match m_f ≈ 1372 kg, re-run, and record the calibrated value
prominently in the results struct and README (spec open item 1 resolution).

- [ ] **Step 4: Mesh-refinement spot check** — re-run M2 with `'N', 1200`,
`'tag','M2_manifold_N1200'`, seeded from M2 (`'seedMat','results/M2_manifold.mat'`).
Gate: |Δm_f| < 1 kg, switch count within ±4 of the N=600 run. Report both.

- [ ] **Step 5: Commit**

```bash
cd /Users/msc/Desktop/optimal_control
git add -A earth_elliptic_to_geo/*.m
git commit --allow-empty -m "milestone(earth-geo): M2 free-Lf manifold matches paper mass (~1370 kg, 7.5 revs, apogee burns)

Co-Authored-By: Claude Fable 5 <noreply@anthropic.com>"
```

---

### Task 14: Milestone M3 — c_tf front + thrust law

**Files:**
- Create: `earth_elliptic_to_geo/run_ctf_sweep.m`
- Output: `results/sweep_*.mat`, `results/front_mf_ctf.png`

**Interfaces:**
- Produces: `run_ctf_sweep()` — sweeps c_tf ∈ {1.2, 1.5, 2.0, 2.5, 3.0} at 10 N and
  T_max ∈ {10, 5, 2.5} N at c_tf=1.5. Point tags `sweep_T<10T>_c<100ctf>`; each
  point saved individually and **skipped if its file exists** (resume after MEX
  crashes). Neighbor-seeds along the c_tf chain from the previous point.

- [ ] **Step 1: Write the sweep script**

```matlab
function run_ctf_sweep()
% RUN_CTF_SWEEP  M3 front: m_f vs c_tf at 10 N + thrust law across {10,5,2.5} N.
%
% Each point is one run_transfer call, saved individually and skipped when its
% results file exists (resume-after-crash pattern; the sporadic CasADi MEX
% fatal kills the process ~1 in 10 solves -- just rerun this script).
%
% INPUTS:  none   OUTPUTS: none (results/sweep_*.mat + front figure + printed table)
% REFERENCES: [1] DESIGN.md sec 5 milestone M3. [2] paper Figs 18/23, law R0.
here = fileparts(mfilename('fullpath'));
resDir = fullfile(here, 'results');
% --- leg 1: c_tf front at 10 N (neighbor-chain upward from M2) ---------------
ctfs = [1.2 1.5 2.0 2.5 3.0];
prev = fullfile(resDir, 'M2_manifold.mat');
for cf = ctfs
    tag = sprintf('sweep_T100_c%03d', round(100*cf));
    fn  = fullfile(resDir, [tag '.mat']);
    if isfile(fn), fprintf('skip %s\n', tag); prev = fn; continue; end
    run_transfer(struct('thrustN',10, 'ctf',cf, 'hx0',0.0612, 'term','manifold', ...
                        'tag',tag, 'seedMat',prev));
    prev = fn;
end
% --- leg 2: thrust law at c_tf = 1.5 (fresh pipelines; N scaled by revs) -----
thr = [10 5 2.5];  Ns = [600 1200 2400];
for kt = 1:numel(thr)
    tag = sprintf('sweep_T%03d_c150', round(10*thr(kt)));
    fn  = fullfile(resDir, [tag '.mat']);
    if isfile(fn), fprintf('skip %s\n', tag); continue; end
    run_transfer(struct('thrustN',thr(kt), 'ctf',1.5, 'hx0',0.0612, ...
                        'term','manifold', 'tag',tag, 'N',Ns(kt)));
end
% --- collect + gates ----------------------------------------------------------
fprintf('\n%-8s %-6s %-9s %-6s %-9s\n', 'T [N]', 'c_tf', 'mf [kg]', 'sw', 'tfmin [h]');
mfC = zeros(1, numel(ctfs));
for kc = 1:numel(ctfs)
    S = load(fullfile(resDir, sprintf('sweep_T100_c%03d.mat', round(100*ctfs(kc)))));
    mfC(kc) = S.res.report.m_f_kg;
    fprintf('%-8g %-6.2f %-9.2f %-6d %-9.1f\n', 10, ctfs(kc), mfC(kc), ...
            S.res.report.switches, S.res.mintime.tfmin_h);
end
C = zeros(1, numel(thr));
for kt = 1:numel(thr)
    S = load(fullfile(resDir, sprintf('sweep_T%03d_c150.mat', round(10*thr(kt)))));
    C(kt) = thr(kt) * S.res.mintime.tfmin_h;
    fprintf('%-8g %-6.2f %-9.2f %-6d %-9.1f\n', thr(kt), 1.5, ...
            S.res.report.m_f_kg, S.res.report.switches, S.res.mintime.tfmin_h);
end
fprintf('law R0: T*tfmin = %s N.h  (spread %.1f%%)\n', mat2str(round(C)), ...
        100*(max(C)-min(C))/mean(C));
fig = figure('Visible','off');
plot(ctfs, mfC, 'o-', 'LineWidth', 1.5);  grid on
xlabel('c_{tf}');  ylabel('m_f [kg]');
title('GEO transfer: final mass vs transfer-time multiplier (T_{max}=10 N)');
exportgraphics(fig, fullfile(resDir, 'front_mf_ctf.png'), 'Resolution', 150);
close(fig);
end
```

- [ ] **Step 2: Run it** (long — hours; if the MATLAB process dies from the MEX
crash, JUST RERUN — completed points are skipped).

**M3 gates:**
- m_f strictly increasing in c_tf (paper Fig 23: ≈1350 → 1388 kg over 1.05→3); any
  non-monotone point = local-minimum scatter → re-run that point seeded from its
  other neighbor and keep the better m_f (best-of lower envelope, campaign pattern)
- m_f(c_tf=1.5) across {10, 5, 2.5} N within ~5 kg of each other (Fig 23 near-
  T_max-independence)
- law R0: `T_max·t_fmin` spread < 10% across the three thrusts (paper C ≈ 850 N·h)

- [ ] **Step 3 (optional stretch):** attempt 1 N (74.5 revs, N=6000) with a generous
timebox: `run_transfer(struct('thrustN',1,'ctf',1.5,'hx0',0.0612,'term','manifold','tag','sweep_T010_c150','N',6000))`.
Report the outcome honestly either way — this is ~2× our demonstrated Cartesian
capacity and "did not converge" is an acceptable, documented result.

- [ ] **Step 4: Commit**

```bash
cd /Users/msc/Desktop/optimal_control
git add earth_elliptic_to_geo/run_ctf_sweep.m
git commit -m "milestone(earth-geo): M3 mf-vs-ctf front + Tmax*tfmin law across 10/5/2.5 N

Co-Authored-By: Claude Fable 5 <noreply@anthropic.com>"
```

---

### Task 15: Transfer movie

**Files:**
- Create: `earth_elliptic_to_geo/transfer_movie.m`
- Output: `results/M2_movie.mp4`, `results/M2_movie.gif`

**Interfaces:**
- Consumes: a `run_transfer` results `.mat`.
- Produces: `transfer_movie(matFile, outStem)` → MP4 + GIF.

- [ ] **Step 1: Write the implementation** (compact psr_movie adaptation; two panels —
3D trajectory colored burn/coast with GEO ring + Earth, throttle strip with cursor;
ΔV/mass text meter)

```matlab
function transfer_movie(matFile, outStem)
% TRANSFER_MOVIE  Animate a LEO-ellipse->GEO min-fuel transfer (MP4 + GIF).
%
% Top panel: inertial 3D trajectory colored by throttle (red burn / blue coast),
% GEO ring + target star; bottom: throttle strip with time cursor; text meter
% with running mass/DeltaV. Frames uniform in PHYSICAL time (X(8,:)).
%
% INPUTS:  matFile - run_transfer results .mat;  outStem - output basename
% OUTPUTS: none (writes <outStem>.mp4 and <outStem>.gif)
%
% REFERENCES: [1] PSR/psr_movie.m (layout). [2] memory: matlab-movie-diagonal-
%   streaks (fixed 1280x720 divisible-by-16 frame -> no H.264 shear).
S = load(matFile);  res = S.res;
p  = kepler_lt_params(res.cfg.thrustN, 1500, 2000);
X = res.fuel.X;  U = res.fuel.U;
r = X(1:3,:);  t = X(8,:);  m = X(7,:);  ss = U(4,:);
tD = t * p.TU_s/86400;  burn = ss > 0.5;
dV = p.c*log(1./m)*p.VU_kms;
th = linspace(0, 2*pi, 361);
fig = figure('Color','w','Position',[80 80 1000 750],'Visible','off');
axT = subplot('Position',[0.06 0.32 0.90 0.62]);  hold(axT,'on'); grid(axT,'on');
plot3(axT, r(1,:), r(2,:), r(3,:), '-', 'Color',[0.8 0.8 0.83], 'LineWidth',0.5);
plot3(axT, cos(th), sin(th), 0*th, 'g-', 'LineWidth', 1.0);
plot3(axT, 0,0,0, 'o', 'MarkerFaceColor',[0.1 0.35 0.8], 'MarkerSize',10);
hB = plot3(axT, nan,nan,nan, 'r-', 'LineWidth',1.8);
hC = plot3(axT, nan,nan,nan, 'b-', 'LineWidth',1.5);
hS = plot3(axT, nan,nan,nan, 'ko', 'MarkerFaceColor','k', 'MarkerSize',5);
hTx = text(axT, 0.02, 0.95, '', 'Units','normalized', 'FontName','Menlo', 'FontSize',10);
axis(axT, 'equal');  view(axT, -30, 25);
title(axT, sprintf('LEO ellipse %s GEO min-fuel  (T=%g N, c_{tf}=%.2f)', ...
      char(8594), res.cfg.thrustN, res.cfg.ctf));
axS = subplot('Position',[0.06 0.06 0.90 0.18]);  hold(axS,'on'); grid(axS,'on');
stairs(axS, tD, ss, '-', 'Color',[0.4 0.4 0.4]);  ylim(axS, [-0.05 1.08]);
xlabel(axS, 'time [days]');  ylabel(axS, 'throttle');
hCur = plot(axS, [0 0], [-0.05 1.08], 'k-');
vw = VideoWriter(outStem, 'MPEG-4');  vw.FrameRate = 24;  vw.Quality = 95;  open(vw);
gifFile = [outStem '.gif'];  gifMap = [];  tmp = [outStem '_tmp.png'];
vidHW = [720 1280];                                  % divisible by 16: no H.264 shear
tFr = linspace(t(1), t(end), 300);
for fc = 1:numel(tFr)
    k = find(t <= tFr(fc), 1, 'last');
    mask = @(vv, mm) subsasgn(vv, substruct('()', {~mm}), nan);
    xb = r(1,1:k);  yb = r(2,1:k);  zb = r(3,1:k);  bm = burn(1:k);
    set(hB, 'XData',mask(xb,bm),  'YData',mask(yb,bm),  'ZData',mask(zb,bm));
    set(hC, 'XData',mask(xb,~bm), 'YData',mask(yb,~bm), 'ZData',mask(zb,~bm));
    set(hS, 'XData',r(1,k), 'YData',r(2,k), 'ZData',r(3,k));
    set(hTx,'String', sprintf('t=%5.1f d  m=%7.1f kg  dV=%5.3f km/s', ...
        tD(k), 1500*m(k), dV(k)));
    set(hCur, 'XData', [tD(k) tD(k)]);
    drawnow;
    exportgraphics(fig, tmp, 'Resolution', 120);
    img = imresize(imread(tmp), vidHW);
    writeVideo(vw, img);
    if mod(fc-1, 2) == 0
        gi = imresize(img, [360 640]);
        if isempty(gifMap)
            [gInd, gifMap] = rgb2ind(gi, 256, 'nodither');
            imwrite(gInd, gifMap, gifFile, 'gif', 'LoopCount', Inf, 'DelayTime', 1/12);
        else
            gInd = rgb2ind(gi, gifMap, 'nodither');
            imwrite(gInd, gifMap, gifFile, 'gif', 'WriteMode','append', 'DelayTime', 1/12);
        end
    end
end
close(vw);  if isfile(tmp), delete(tmp); end
close(fig);
fprintf('WROTE %s.mp4 / .gif (%d frames)\n', outStem, numel(tFr));
end
```

- [ ] **Step 2: Render on M2 + eyeball 3 extracted frames** (VideoReader → PNG; check
no diagonal-streak artifact, burn arcs cluster at apogee).

- [ ] **Step 3: Commit**

```bash
cd /Users/msc/Desktop/optimal_control
git add earth_elliptic_to_geo/transfer_movie.m
git commit -m "feat(earth-geo): transfer movie (burn/coast trajectory + throttle strip)

Co-Authored-By: Claude Fable 5 <noreply@anthropic.com>"
```

---

### Task 16: README + wrap-up

**Files:**
- Create: `earth_elliptic_to_geo/README.md`
- Modify: `earth_elliptic_to_geo/DESIGN.md` (status line → implemented; record the
  Isp actually used and the M2/M3 numbers)

- [x] **Step 1: Write README.md** — one page: problem (paper citation), method mapping
(their indirect homotopy = our direct homotopy), module map (the 11 files), how to
run (`run_mintime` → `run_transfer` → `run_ctf_sweep` → `transfer_movie`), the
milestone results table (OUR numbers vs paper: t_fmin, m_f, revs, switches,
apogee-burn), test suite list (`test_params ... test_energy_stage`), honest
caveats (trapezoid switch-count mesh sensitivity; local minima; Isp provenance;
1 N stretch outcome). Follow the campaign README style
(`GTO_tulip/README.md`).

- [x] **Step 2: Run the full no-solve test suite once more**

Run (updated to the full 8-test set per task-16 instructions, superseding the
plan's original 6-test list): `test_params; test_elements; test_dynamics;
test_terminal; test_seed; test_solver_smoke; test_stall_guard; test_p2_homotopy`
Result: all 8 `ALL PASS`.

- [x] **Step 3: Commit**

```bash
cd /Users/msc/Desktop/optimal_control
git add earth_elliptic_to_geo/README.md earth_elliptic_to_geo/DESIGN.md earth_elliptic_to_geo/PLAN.md
git commit -m "docs(earth-geo): README + results summary; close out reproduction milestones

Co-Authored-By: Claude Fable 5 <noreply@anthropic.com>"
```

---

## Plan self-review notes (kept for the record)

- **Spec coverage:** every DESIGN.md module has a task (params T1, elements T3,
  terminal T5, solver core T7, seed T6, homotopy T9, run_transfer T10, sweep T14,
  verify T12, movie T15); milestones M0–M3 = T10/T11/T13/T14; P₂ toy test = T2;
  open items: units→T1, Isp→T13 step 1, free-L_f local minima→T14 gate 1 +
  T9/T10 law-based L_f prescription.
- **Known simplifications vs spec text:** M0 runs coplanar data through the 3D
  machinery (spec amended to match); the two spec'd solver siblings are one
  mode-switched core (spec amended); no `casadi_mintime_2body`/`casadi_minfuel_2body`
  files exist — `run_mintime` + mode `'mintime'` covers the former.
- **Type consistency:** solver `out` fields fixed in T7 and consumed verbatim in
  T9–T15; seed returns `[sigma, X0(9×), U0, tauf0, info]` everywhere; `term` structs
  from T5 only.
- **Honesty gates baked in:** never advance a failed homotopy step; mesh-refinement
  checks at M1 (conditional) and M2 (mandatory); switch-count bands not points;
  1 N stretch explicitly allowed to fail with an honest report.
