# Elliptic→GEO with Lunar Gravity — Phase-1 Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Add an opt-in lunar third-body perturbation to the certified 2-body MEE solver, bridge the certified 10 N solution into the CR3BP by μ-continuation, sharpen to a certified CR3BP min-fuel solution, and produce the first Δm_f comparison point.

**Architecture:** Spec `docs/superpowers/specs/2026-07-22-elliptic-geo-cr3bp-phase0-design.md` (authoritative — read it; D1–D7 are locked decisions). The pert spec rides inside `opts.par` (`par.pert`), so ONLY `lt_mee_rhs.m` changes in the shared core — `casadi_lt_mee`/`homotopy_mee` thread `par` through untouched (their post-solve numeric re-check calls the same RHS, so defects stay honest automatically). New campaign code lives in `earth_elliptic_to_geo_CR3BP/direct/`.

**Tech Stack:** MATLAB R2025b (`/Applications/MATLAB_R2025b.app/bin/matlab -batch "<ONE line>"`), CasADi 3.7.0 (`CASADI_PATH` or `~/casadi-3.7.0`), git.

## Global Constraints

- Paths: `E2B = orbit_transfer/earth_elliptic_to_geo/direct` (2-body campaign), `E3B = orbit_transfer/earth_elliptic_to_geo_CR3BP/direct` (this campaign), repo root `/Users/msc/Desktop/optimal_control`.
- **Back-compat invariant (binding, spec §8 gates 1–2):** `par` without a `pert` field ⇒ `lt_mee_rhs` byte-identical behavior; `pert` present with `gain=0` ⇒ matches to solver tolerance. Any nominal-path change is a STOP.
- **MX-safety (binding, RHS edits):** no `norm`/`abs`/`max`/`if`/`mod` on state-dependent quantities; guards via `+1e-12` inside powers and `fmax` only where the existing file already does; `cross()` written out component-wise. Follow `lt_mee_rhs.m`'s existing style exactly.
- **The indirect term** in a_M is mandatory (spec §3); the perturbation contributes **nothing to ṁ**.
- MATLAB house style: full headers (purpose/INPUTS with sizes/OUTPUTS/REFERENCES); never `i`/`j` as loop vars; one-line `-batch`; filter license noise with `| grep -vE "License|academic|personal use|government"`.
- `.mat` results are gitignored — never `git add` one. The 2-body campaign's certified caches are read-only.
- Stage exactly the files each task names; never `git add -A`. Commit trailer (exact): `Co-Authored-By: Claude Fable 5 <noreply@anthropic.com>`
- On unexpected failure (test fails, nominal drift, solver surprises): STOP, report BLOCKED with exact output.

---

### Task 1: `lunar_params.m` + unit tests

**Files:**
- Create: `E3B/lunar_params.m`, `E3B/setup_paths.m`
- Test: `E3B/test_lunar_params.m`

**Interfaces:**
- Produces: `pert = lunar_params(par, phi0, gain)` — `par` a `kepler_lt_params` struct; returns `pert = struct('muM','DM','nM','phi0','gain')` in the campaign's canonical units (mu=1, LU=42165 km). Later tasks pass this as `opts.par.pert`.

- [ ] **Step 1: Write `E3B/setup_paths.m`** (the campaign's path bootstrap):

```matlab
function setup_paths()
% SETUP_PATHS  Path bootstrap for the elliptic->GEO CR3BP campaign.
%
% Adds this folder plus the 2-body campaign's core/ (shared solver, spec D3)
% so lunar_params / casadi_lt_mee / homotopy_mee / mee_seed all resolve.
%
% INPUTS:  none
% OUTPUTS: none (modifies the MATLAB path)
% REFERENCES:
%   [1] docs/superpowers/specs/2026-07-22-elliptic-geo-cr3bp-phase0-design.md sec 4.
here = fileparts(mfilename('fullpath'));
addpath(here);
addpath(fullfile(here, '..', '..', 'earth_elliptic_to_geo', 'direct'));
addpath(fullfile(here, '..', '..', 'earth_elliptic_to_geo', 'direct', 'core'));
end
```

- [ ] **Step 2: Write the failing test** `E3B/test_lunar_params.m`:

```matlab
% TEST_LUNAR_PARAMS  Unit conversions + physical sanity of the Moon spec.
here = fileparts(mfilename('fullpath'));  run(fullfile(here,'setup_paths.m'));
par  = kepler_lt_params(10, 1500, 2000);
pert = lunar_params(par, 0, 1);
% (a) canonical-unit roundtrips against physical values
assert(abs(pert.DM * par.LU_km - 384400) < 1e-6, 'DM roundtrip [km]');
assert(abs(pert.muM/par.mu - 4902.800/398600.47) < 1e-9, 'mass-ratio 0.0123');
periodDays = (2*pi/pert.nM) * par.TU_s / 86400;
assert(abs(periodDays - 27.32) < 0.05, 'sidereal month ~27.32 d');
% (b) sidereal rate consistent with two-body circular orbit about barycenter
assert(abs(pert.nM - sqrt((par.mu + pert.muM)/pert.DM^3)) < 1e-14, 'nM formula');
% (c) tidal acceleration at GEO radius ~ 7.3e-6 m/s^2 (spec sec 7)
aTide = 2*pert.muM*1/pert.DM^3 * par.AU_ms2;       % r = 1 LU = GEO radius
assert(abs(aTide - 7.3e-6) < 0.3e-6, 'lunar tide at GEO ~7.3e-6 m/s^2');
% (d) phi0/gain pass-throughs
p2 = lunar_params(par, 1.25, 0.5);
assert(p2.phi0 == 1.25 && p2.gain == 0.5, 'phi0/gain stored');
fprintf('test_lunar_params: ALL PASS\n');
```

- [ ] **Step 3: Run to verify it fails** —
`/Applications/MATLAB_R2025b.app/bin/matlab -batch "run('/Users/msc/Desktop/optimal_control/orbit_transfer/earth_elliptic_to_geo_CR3BP/direct/test_lunar_params.m')"`
Expected: FAIL (`lunar_params` undefined).

- [ ] **Step 4: Implement** `E3B/lunar_params.m`:

```matlab
function pert = lunar_params(par, phi0, gain)
% LUNAR_PARAMS  Lunar third-body constants in the 2-body campaign's units.
%
% Circular Moon in the reference plane (spec D2): geocentric distance D_EM,
% sidereal rate n_M = sqrt((mu_E+mu_M)/D_EM^3), phase phi0 at t=0. Expressed
% in kepler_lt_params canonical units (mu_E = 1, LU = 42165 km GEO radius).
% The continuation scale `gain` multiplies mu_M (spec D5: gain IS the
% mu-continuation knob; gain=1 is the physical Moon).
%
% INPUTS:
%   par  - kepler_lt_params struct (.mu .muKm3s2 .LU_km .TU_s) [struct]
%   phi0 - lunar phase at t=0 [rad, scalar] (spec D6; baseline 0)
%   gain - mu_M continuation scale in [0,1] [scalar]
%
% OUTPUTS:
%   pert - struct .muM (canonical, PHYSICAL value -- gain applied in the
%          RHS, not here) .DM .nM (canonical) .phi0 .gain
%
% REFERENCES:
%   [1] spec 2026-07-22-elliptic-geo-cr3bp-phase0-design.md sec 3 (constants:
%       mu_M = 4902.800 km^3/s^2, D_EM = 384400 km; ratio 0.0123 consistent
%       with CR3BP mu* = 0.0121506 via mu*/(1-mu*)).
if nargin < 2 || isempty(phi0), phi0 = 0; end
if nargin < 3 || isempty(gain), gain = 1; end
muM_km3s2 = 4902.800;
DM_km     = 384400;
pert.muM  = (muM_km3s2 / par.muKm3s2) * par.mu;          % canonical GM_moon
pert.DM   = DM_km / par.LU_km;                           % canonical distance
pert.nM   = sqrt((par.mu + pert.muM) / pert.DM^3);       % canonical rate
pert.phi0 = phi0;
pert.gain = gain;
end
```

- [ ] **Step 5: Run test — PASS.**
- [ ] **Step 6: Commit** (`git add` the three files; message `feat(cr3bp-geo): lunar third-body params + campaign path bootstrap (phase1 T1)` + trailer).

---

### Task 2: `lt_mee_rhs` opt-in perturbation (the ONLY shared-core edit)

**Files:**
- Modify: `E2B/core/lt_mee_rhs.m`
- Test: `E3B/test_lt_mee_rhs_pert.m`

**Interfaces:**
- Produces: `lt_mee_rhs(X, U, par)` unchanged signature; if `par.pert` exists (fields per Task 1), the lunar acceleration (direct + indirect, scaled by `pert.gain`) is added to the RTN force terms of every Gauss equation INCLUDING Ldot; `ṁ` untouched. `par` without `.pert` ⇒ code path identical to today.

- [ ] **Step 1: Write the failing test** `E3B/test_lt_mee_rhs_pert.m`:

```matlab
% TEST_LT_MEE_RHS_PERT  Opt-in third-body term: back-compat + exact oracle.
here = fileparts(mfilename('fullpath'));  run(fullfile(here,'setup_paths.m'));
par = kepler_lt_params(10, 1500, 2000);
Xr = [1.3; 0.2; -0.1; 0.05; 0.02; 0.9; 2.0];   % generic elliptic 3D state
Ur = [0.36; 0.48; 0.80; 0.7];                   % unit beta, mid throttle
par.L = 1.1;
[d0, L0] = lt_mee_rhs(Xr, Ur, par);
% (a) gate 1: pert ABSENT vs gain=0 vs gain=1e-30 -- identical / continuous
parG0 = par;  parG0.pert = lunar_params(par, 0, 0);
[dg0, Lg0] = lt_mee_rhs(Xr, Ur, parG0);
assert(isequal(d0, dg0) || max(abs(d0-dg0)) < 1e-15, 'gain=0 == pert-absent');
assert(abs(L0-Lg0) < 1e-15, 'Ldot unchanged at gain=0');
parGt = par;  parGt.pert = lunar_params(par, 0, 1e-30);
[dgt, ~] = lt_mee_rhs(Xr, Ur, parGt);
assert(max(abs(d0-dgt)) < 1e-25, 'continuity in gain');
% (b) exact oracle: equatorial circular state at L=0, Moon at phi0=0, t=0.
%     r=[1;0;0], Rhat=[1;0;0], That=[0;1;0], Nhat=[0;0;1]; Moon on +x =>
%     a_M = muM*(1/(DM-1)^2 - 1/DM^2) purely RADIAL (+x). Gauss response:
%     dP/dt += 0; dex/dt += sqrt(P)*sin(L)*aR = 0; dey/dt += -sqrt(P)*cos(L)*aR
%     = -aR; dhx,dhy,Ldot,mdot unchanged.
Xe = [1; 0; 0; 0; 0; 1; 0];  Ue = [0;1;0; 0];   % zero throttle isolates pert
pe = par;  pe.L = 0;  pe.pert = lunar_params(par, 0, 1);
[dp, Lp] = lt_mee_rhs(Xe, Ue, pe);
pe0 = par; pe0.L = 0;
[dq, Lq] = lt_mee_rhs(Xe, Ue, pe0);
aR = pe.pert.muM * (1/(pe.pert.DM-1)^2 - 1/pe.pert.DM^2);
assert(abs(Lp - Lq) < 1e-15, 'radial pert does not touch Ldot');
delta = (dp - dq) * Lq;                          % back to d/dt via common Ldot
assert(abs(delta(1)) < 1e-14, 'dP unchanged under radial accel');
assert(abs(delta(3) - (-aR)) < 1e-12, 'dey/dt == -aR (radial oracle)');
assert(all(abs(delta([2 4 5 6 7])) < 1e-14), 'ex,hx,hy,m,t untouched');
% (c) frame identities at the generic 3D state (analytic: |Nhat|=1, R.N=0)
%     -- exercised inside the RHS; here we just confirm no NaN and mdot clean
parP = par;  parP.pert = lunar_params(par, 0.7, 1);
[dfull, ~] = lt_mee_rhs(Xr, Ur, parP);
assert(all(isfinite(dfull)), 'finite with pert on');
assert(abs(dfull(6) - d0(6)) < 1e-15, 'mdot has NO perturbation coupling');
fprintf('test_lt_mee_rhs_pert: ALL PASS\n');
```

- [ ] **Step 2: Run — FAIL** (pert field unused today; oracle asserts fire).
- [ ] **Step 3: Implement.** Rewrite the force section of `E2B/core/lt_mee_rhs.m`. Replace the block from `thr = U(4);` through `tdot = 1;` with (everything above/below unchanged, including the header — extend the header's INPUTS to document `par.pert` and REFERENCES to cite the spec):

```matlab
thr = U(4);
q = thr*U(1);  s = thr*U(2);  w = thr*U(3);

Lw = L - 2*pi*floor(L/(2*pi));
cL = cos(Lw);  sL = sin(Lw);
Z  = 1 + ex*cL + ey*sL;
A1 = ex + (1+Z)*cL;
A2 = ey + (1+Z)*sL;
Xh = 1 + hx^2 + hy^2;
hterm = hx*sL - hy*cL;
sqPmu = sqrt(P/mu);

% Total RTN specific force: thrust acceleration (Tm/m per unit component)
% plus, opt-in, the lunar third-body acceleration (spec D1/D3, 2026-07-22:
% direct + MANDATORY indirect term; scaled by pert.gain = the mu-continuation
% knob; a PURE acceleration -- no mdot coupling). Absent par.pert, fR/fT/fN
% reduce exactly to the historical (Tm/m)*[q;s;w].
fR = (Tm/m)*q;  fT = (Tm/m)*s;  fN = (Tm/m)*w;
if isfield(par, 'pert') && ~isempty(par.pert)
    pM = par.pert;
    r  = P/Z;  alpha2 = hx^2 - hy^2;
    rx = (r/Xh)*(cL + alpha2*cL + 2*hx*hy*sL);
    ry = (r/Xh)*(sL - alpha2*sL + 2*hx*hy*cL);
    rz = (2*r/Xh)*(hx*sL - hy*cL);
    % RTN basis in inertial axes (MX-safe closed forms; |Nhat|=1 identically)
    Rx = rx/r;  Ry = ry/r;  Rz = rz/r;
    Nx = 2*hy/Xh;  Ny = -2*hx/Xh;  Nz = (1 - hx^2 - hy^2)/Xh;
    Tx = Ny*Rz - Nz*Ry;  Ty = Nz*Rx - Nx*Rz;  Tz = Nx*Ry - Ny*Rx;
    tState = X(7);
    ang = pM.nM*tState + pM.phi0;
    rMx = pM.DM*cos(ang);  rMy = pM.DM*sin(ang);          % Moon in ref plane
    dx = rMx - rx;  dy = rMy - ry;  dz = -rz;
    d3  = (dx^2 + dy^2 + dz^2 + 1e-12)^1.5;               % >= (8 LU)^3; guard inert
    DM3 = pM.DM^3;
    gm  = pM.gain * pM.muM;
    aX = gm*(dx/d3 - rMx/DM3);                            % direct + indirect
    aY = gm*(dy/d3 - rMy/DM3);
    aZ = gm*(dz/d3);                                      % Moon z == 0
    fR = fR + (Rx*aX + Ry*aY + Rz*aZ);
    fT = fT + (Tx*aX + Ty*aY + Tz*aZ);
    fN = fN + (Nx*aX + Ny*aY + Nz*aZ);
end

Pdot  = 2*(P*sqPmu) * (fT/Z);
exdot = sqPmu*(1/Z)*( Z*sL*fR + A1*fT - ey*hterm*fN );
eydot = sqPmu*(1/Z)*(-Z*cL*fR + A2*fT + ex*hterm*fN );
hxdot = (1/2)*sqPmu*(Xh/Z)*cL*fN;
hydot = (1/2)*sqPmu*(Xh/Z)*sL*fN;
% NOTE: the paper's printed L-dot equation (p.6) omits Tmax on the thrust
% term -- a typo; see the pre-2026-07-22 version of this file for the full
% note. The force term now uses fN so the opt-in perturbation consistently
% enters the L-rate too.
Ldot  = sqrt(mu/P^3)*Z^2 + sqPmu*(1/Z)*hterm*fN;
mdot  = -(Tm/c)*thr;             % thrust only: gravity costs no propellant
tdot  = 1;
```

Verify by inspection: with `par.pert` absent, `fR/fT/fN = (Tm/m)*[q;s;w]` and every equation is algebraically identical to the original (the `(Tm/m)` and `(2*Tm/m)`/`(Tm/(2m))` factors distribute exactly).

- [ ] **Step 4: Run the new test — PASS.** Then run the 2-body campaign's own fast suite to prove nominal untouched: `cd E2B; setup_paths; ` then its `tests/` no-solve suite (list files with `ls tests/`, run each `test_*.m` that the campaign's README/tests header marks as no-solve). All must PASS unchanged. Any diff = STOP (back-compat invariant).
- [ ] **Step 5: Commit** (`E2B/core/lt_mee_rhs.m` + `E3B/test_lt_mee_rhs_pert.m`; message `feat(mee-core): opt-in lunar third-body term in lt_mee_rhs (cr3bp-geo T2)` + trailer).

---

### Task 3: `sanity_bound.m` — the §7 null model from code

**Files:**
- Create: `E3B/sanity_bound.m` (script, campaign-style)

- [ ] **Step 1: Locate the certified ladder table**: `grep -rl "1377" /Users/msc/Desktop/optimal_control/orbit_transfer/earth_elliptic_to_geo/process/ | head -3` — the file holding the certified per-rung t_fMin / t_f / m_f values (completed in commit b6fb363). Read it; use its exact numbers.
- [ ] **Step 2: Write `sanity_bound.m`**: full-header script that, for each certified rung (10, 5, 2.5, 1, 0.5, 0.2, 0.1 N — take the set actually in the table): computes thrust authority `T/m0` in m/s², lunar tide `2*muM*r/DM^3 * AU_ms2` at r = 1 LU, the ratio, t_f in days (from the table), and t_f in lunar months (`/27.32 d`). Print an aligned table and write it to `E3B/results/sanity_bound.md` (create `results/` + `.gitignore` with `*.mat`). Assert the 10 N ratio is between 0.05% and 0.5% and the 0.1 N ratio between 5% and 20% (loose brackets on the spec §7 predictions — if these fire, the spec's null model is materially wrong: STOP and report, per spec §9 honesty rule).
- [ ] **Step 3: Run it** (no solves; seconds). Table prints, asserts pass.
- [ ] **Step 4: Commit** (`sanity_bound.m` + `results/sanity_bound.md` + `results/.gitignore`; message `feat(cr3bp-geo): sanity-bound null model table (phase1 T3)` + trailer).

---

### Task 4: `bridge_mu_continuation.m` — 10 N energy bridge (gates 1–3)

**Files:**
- Create: `E3B/bridge_mu_continuation.m`

**Interfaces:**
- Produces: `out = bridge_mu_continuation(opts)` with `opts.thrustN [10]`, `.phi0 [0]`, `.gainSched [[0.25 0.5 0.75 1.0]]`, `.maxIter [1500]`, `.resume [true]`. Saves `E3B/results/energy_cr3bp_T<mN|N-tag>_phi<..>.mat` holding `sigma, X, U, dL, tfTarget, fp` (fp incl. `thrustN,m0kg,ispS,tfTarget,muM,DM,nM,phi0,gain`). Returns the final solver `out` + `.gainReached`.

- [ ] **Step 1: Write the driver.** Structure (full header per house style; body):
  1. `setup_paths`; `par = kepler_lt_params(opts.thrustN, 1500, 2000)`; `tfTarget` = the certified table's t_f for this rung (Task 3 located the table; hard-code the 10 N value into a `tf_table` local function with a comment citing the file, covering every rung listed there).
  2. Seed: `[sigma, X0, U0, dL0] = mee_seed(par, <the same seed opts run_transfer_mee.m uses for this rung — read that driver's stage-1 block and mirror its exact mee_seed call>)`.
  3. **Gate-1 baseline (pert absent):** solve `casadi_lt_mee(sigma, X0, U0, dL0, struct('par',par,'mode','fixedtf','eps',1,'tfTarget',tfTarget,'x0',X0(:,1),'maxIter',opts.maxIter,'warmTight',false))`; require `Solve_Succeeded && maxDefect < 1e-6`.
  4. **Gate-2 (gain=0 with pert PRESENT):** `par2 = par; par2.pert = lunar_params(par,opts.phi0,0);` warm-start tight from step 3's X/U/dL; require success AND `max(abs(X - X_gate1))` < 1e-8 (the hook itself introduces no drift). Print both gate verdicts.
  5. **Gain walk:** for each `g` in `gainSched`: `par2.pert.gain = g`; tight solve warm from previous; gate `Solve_Succeeded && maxDefect<1e-6`; on failure halve the step toward the last good gain (insert midpoint into the schedule; floor step 0.05, else error `bridge:stuck` with the honest gain floor). Per-step checkpoint save (resume-safe: skip steps whose checkpoint exists and fp matches via the campaign's `check_cache_fp` pattern — replicate that helper locally with the fp fields above).
  6. Save the final artifact + print `BRIDGE: gain=1 reached, defect=..., m_f(energy)=... kg`.
- [ ] **Step 2: Lint** (`checkcode` clean) then **run at 10 N** (one-line batch; 10 N solves are minutes each, ~6 solves total — run in foreground with a generous timeout). Expected: both gates pass; gain walk reaches 1.0 with tiny per-step movement (§7: perturbation ~0.1% of authority at 10 N).
- [ ] **Step 3: Commit** (driver only — the .mat stays untracked; message `feat(cr3bp-geo): mu-continuation bridge driver + 10N energy bridge run (phase1 T4)` + trailer). Record the gate/defect numbers in the commit body.

---

### Task 5: `solve_cr3bp_minfuel.m` — ε-sharpen at gain=1 (gate 3 complete)

**Files:**
- Create: `E3B/solve_cr3bp_minfuel.m`

**Interfaces:**
- Produces: `best = solve_cr3bp_minfuel(opts)` (`opts.thrustN [10]`, `.phi0 [0]`, `.maxIter [1500]`): loads the Task-4 artifact (fp-checked), runs `homotopy_mee(sigma, X, U, dL, ho)` with `ho = struct('par', parWithPertGain1, 'x0', X(:,1), 'tfTarget', tfTarget, 'maxIter', opts.maxIter, 'resDir', <E3B/results/homotopy>, 'tag', sprintf('cr3bp_%s_phi%s', <rung tag>, <phi tag>), 'fp', <the same fp struct>)`; asserts `best.certified` and `best.epsReached == 0`; saves `E3B/results/minfuel_cr3bp_<tags>.mat` (best + provenance) and prints `m_f_kg`, switches, edge, defect beside the 2-body certified m_f for the rung (from the Task-4 `tf_table`, extended with the table's m_f column).

- [ ] **Step 1: Write the driver** (structure above; full header; the homotopy's per-step caches land under `results/homotopy/` via its own resume machinery — do not reimplement).
- [ ] **Step 2: Run at 10 N** (the ε-schedule is ~14 steps of minutes each at 10 N; run foreground with timeout 600000, resume-safe if it needs a second call). Expected: `certified=1`, `epsReached=0`, and **Δm_f vs 1377.10 kg small** (spec §7 predicts ~0.1%-of-authority effects at 10 N — likely sub-0.1 kg; whatever it is, it is the first result: report it exactly).
- [ ] **Step 3: Commit** (driver; message `feat(cr3bp-geo): CR3BP min-fuel sharpen driver + certified 10N solve (phase1 T5)` + trailer, Δm_f in the body).

---

### Task 6: `compare_vs_2body.m` + docs close-out

**Files:**
- Create: `E3B/compare_vs_2body.m`
- Modify: `orbit_transfer/earth_elliptic_to_geo_CR3BP/README.md`, `TODO.md`; spec status line.

- [ ] **Step 1: Write `compare_vs_2body.m`**: loads the Task-5 artifact(s) present in `results/`, prints a per-rung table — 2-body certified m_f (tf_table) vs CR3BP m_f, Δm_f [kg and %], switch counts BOTH reported as nodal counts with an explicit "mesh-band caveat, P0 protocol" note, defect, t_f in lunar months, predicted §7 ratio — and writes `results/compare_vs_2body.md`.
- [ ] **Step 2: Run it**; commit the generated `.md` alongside the script.
- [ ] **Step 3: Docs.** README: Status section → "Phase 1 IN PROGRESS — 10 N bridged + certified (Δm_f = <measured>); deep rungs pending" with a pointer to the compare table. TODO: check off the Phase-0 items (all four decided by the spec) and the Phase-1 10 N items; add explicit next entries: "walk 1 N / 0.2 N / 0.1 N rungs (background-length solves)" and "φ₀ sweep experiment (spec D6)". Spec header: `**Status:** Phase 1 implemented through the 10 N rung (plan 2026-07-22-elliptic-geo-cr3bp-phase1.md); deep rungs + phi0 sweep pending.`
- [ ] **Step 4: Commit + push.**

```bash
cd /Users/msc/Desktop/optimal_control
git add orbit_transfer/earth_elliptic_to_geo_CR3BP/direct/compare_vs_2body.m orbit_transfer/earth_elliptic_to_geo_CR3BP/direct/results/compare_vs_2body.md orbit_transfer/earth_elliptic_to_geo_CR3BP/README.md orbit_transfer/earth_elliptic_to_geo_CR3BP/TODO.md docs/superpowers/specs/2026-07-22-elliptic-geo-cr3bp-phase0-design.md
git commit -m "feat(cr3bp-geo): 10N comparison table + Phase-1 close-out docs (phase1 T6)

Co-Authored-By: Claude Fable 5 <noreply@anthropic.com>"
git push origin main
```

---

## Self-Review notes

- **Spec coverage:** §2 D1–D7 → D1/D3 (T2 hook via par.pert), D2/D6 (lunar_params: reference-plane ephemeris, φ₀ first-class + fingerprinted), D4 (tf_table from the certified 2-body values), D5 (T4 gain walk on ε=1, then T5 sharpen), D7 (campaign-local params, kepler units — T1). §3 dynamics → T2 (indirect term, no-ṁ rule, t-state usage, gain knob). §4 files → T1/T4/T5/T6 (+`sanity_bound` T3); fingerprints per §4 via fp fields muM/phi0/gain at every cache. §7 → T3 with STOP brackets. §8 gates → gate 1 (T2 Step 4 suite + T4 step 3), gate 2 (T4 step 4 + T2 test (a)), gate 3 (T4 walk + T5 certified), gate 4 (T6 band-caveated table). §9 risks → resume-safe checkpoints (T4/T5), φ₀ fingerprint, honesty STOPs.
- **Placeholder scan:** two bounded delegations, each anchored to a named file the implementer must read: the `mee_seed` call mirrored from `run_transfer_mee.m` stage 1 (T4.1.2) and the certified table located by the T3 grep — both are existing-artifact lookups, not design gaps. All new-file code is complete.
- **Type consistency:** `pert` struct fields identical across `lunar_params` (producer), `lt_mee_rhs` (consumer), fp records, and both drivers; `tf_table` shared by T4/T5/T6 (defined once in T4, extended with m_f in T5 — implementer keeps it a single local function copied verbatim, or promotes it to `E3B/tf_table.m` if used thrice, which T6 does: promote it).
- **Known judgment point (recorded):** the T2 rewrite redistributes `(Tm/m)` factors into `fR/fT/fN`; the plan requires inspection-equivalence plus the oracle test plus the 2-body no-solve suite as the tripwire. If any 2-body test moves: STOP.
