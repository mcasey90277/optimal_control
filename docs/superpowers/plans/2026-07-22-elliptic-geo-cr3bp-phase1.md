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
% Delegate to the 2-body campaign's OWN setup_paths (adds core/, lib/ --
% optdef.m lives there -- coords/, etc.) via the cwd trick: both files are
% named setup_paths.m, and cwd precedence resolves the local one (review
% amendment C, 2026-07-22).
e2b = fullfile(here, '..', '..', 'earth_elliptic_to_geo', 'direct');
oldd = cd(e2b);  setup_paths;  cd(oldd);
end
```

- [ ] **Step 2: Write the failing test** `E3B/test_lunar_params.m`:

```matlab
% TEST_LUNAR_PARAMS  Unit conversions + physical sanity of the Moon spec.
here = fileparts(mfilename('fullpath'));  addpath(here);  setup_paths;
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
here = fileparts(mfilename('fullpath'));  addpath(here);  setup_paths;
par = kepler_lt_params(10, 1500, 2000);
Xr = [1.3; 0.2; -0.1; 0.05; 0.02; 0.9; 2.0];   % generic elliptic 3D state
Ur = [0.36; 0.48; 0.80; 0.7];                   % unit beta, mid throttle
par.L = 1.1;
[d0, L0] = lt_mee_rhs(Xr, Ur, par);
% (a) gate 1: pert ABSENT vs gain=0 vs gain=1e-30 -- identical / continuous
parG0 = par;  parG0.pert = lunar_params(par, 0, 0);
[dg0, Lg0] = lt_mee_rhs(Xr, Ur, parG0);
assert(isequal(d0, dg0) && isequal(L0, Lg0), 'gain=0 takes nominal branch: BITWISE identical');
parGt = par;  parGt.pert = lunar_params(par, 0, 1e-12);
[dgt, ~] = lt_mee_rhs(Xr, Ur, parGt);
assert(max(abs(d0-dgt)) < 1e-9, 'cross-branch continuity at tiny gain (FP noise << 1e-9)');
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
assert(abs(Lp - Lq) < 1e-15, 'radial pert does not touch Ldot');   % check FIRST
delta = (dp - dq) * Lp;                          % back to d/dt (Lp==Lq just proven)
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
- [ ] **Step 3: Implement.** Rewrite the force section of `E2B/core/lt_mee_rhs.m`. Replace the block from `thr = U(4);` through `tdot = 1;` with the BRANCHED version below (review amendment A, 2026-07-22: the nominal branch is the LITERAL original equations — FP non-associativity means a factored rewrite is NOT bitwise-identical, so pert-absent AND gain==0 both take the untouched path; `pertOn` is decided on plain doubles at graph-build time, never on MX). Everything above/below stays unchanged; keep the original "paper's printed L-dot typo" NOTE comment attached to the nominal branch; extend the header's INPUTS to document `par.pert` and REFERENCES to cite the spec:

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

% Opt-in lunar third-body term (spec D1/D3; review amendment A 2026-07-22):
% pertOn is a PLAIN-DOUBLE decision at graph-build time, so the nominal
% branch is the literal pre-2026-07-22 code -- bitwise-identical when pert
% is absent OR gain == 0.
pertOn = isfield(par,'pert') && ~isempty(par.pert) && par.pert.gain > 0;
if ~pertOn
    Pdot  = (2*Tm/m)*(P*sqPmu) * (s/Z);
    exdot = (Tm/m)*sqPmu*(1/Z)*( Z*sL*q + A1*s - ey*hterm*w );
    eydot = (Tm/m)*sqPmu*(1/Z)*(-Z*cL*q + A2*s + ex*hterm*w );
    hxdot = (Tm/(2*m))*sqPmu*(Xh/Z)*cL*w;
    hydot = (Tm/(2*m))*sqPmu*(Xh/Z)*sL*w;
    Ldot  = sqrt(mu/P^3)*Z^2 + (Tm/m)*sqPmu*(1/Z)*hterm*w;
else
    % Total RTN specific force: thrust + lunar third body (direct + MANDATORY
    % indirect term; pert.gain = the mu-continuation knob; PURE acceleration,
    % no mdot coupling).
    pM = par.pert;
    r  = P/Z;  alpha2 = hx^2 - hy^2;
    rx = (r/Xh)*(cL + alpha2*cL + 2*hx*hy*sL);
    ry = (r/Xh)*(sL - alpha2*sL + 2*hx*hy*cL);
    rz = (2*r/Xh)*(hx*sL - hy*cL);
    Rx = rx/r;  Ry = ry/r;  Rz = rz/r;
    Nx = 2*hy/Xh;  Ny = -2*hx/Xh;  Nz = (1 - hx^2 - hy^2)/Xh;
    Tx = Ny*Rz - Nz*Ry;  Ty = Nz*Rx - Nx*Rz;  Tz = Nx*Ry - Ny*Rx;
    tState = X(7);
    ang = pM.nM*tState + pM.phi0;
    rMx = pM.DM*cos(ang);  rMy = pM.DM*sin(ang);          % Moon in ref plane
    dx = rMx - rx;  dy = rMy - ry;  dz = -rz;
    d3  = (dx^2 + dy^2 + dz^2 + 1e-12)^1.5;               % sep >= 8 LU; guard inert
    DM3 = pM.DM^3;
    gm  = pM.gain * pM.muM;
    aX = gm*(dx/d3 - rMx/DM3);                            % direct + indirect
    aY = gm*(dy/d3 - rMy/DM3);
    aZ = gm*(dz/d3);                                      % Moon z == 0
    fR = (Tm/m)*q + (Rx*aX + Ry*aY + Rz*aZ);
    fT = (Tm/m)*s + (Tx*aX + Ty*aY + Tz*aZ);
    fN = (Tm/m)*w + (Nx*aX + Ny*aY + Nz*aZ);
    Pdot  = 2*(P*sqPmu) * (fT/Z);
    exdot = sqPmu*(1/Z)*( Z*sL*fR + A1*fT - ey*hterm*fN );
    eydot = sqPmu*(1/Z)*(-Z*cL*fR + A2*fT + ex*hterm*fN );
    hxdot = (1/2)*sqPmu*(Xh/Z)*cL*fN;
    hydot = (1/2)*sqPmu*(Xh/Z)*sL*fN;
    Ldot  = sqrt(mu/P^3)*Z^2 + sqPmu*(1/Z)*hterm*fN;
end
mdot  = -(Tm/c)*thr;             % thrust only: gravity costs no propellant
tdot  = 1;
```

Verify by inspection: the nominal branch is character-identical to the original equations (diff it against git HEAD to prove it); the pert branch's factored forms are the algebraically-equal redistribution (both external reviewers verified the algebra; the branch split exists purely for bitwise FP fidelity).

- [ ] **Step 4: Run the new test — PASS.** Then the bounded nominal-regression set (review amendment F): (a) `git diff` the nominal branch of `lt_mee_rhs.m` against HEAD~ — the equations must be character-identical; (b) `checkcode` clean; (c) the solver-level proof lands in Task 4's gate-1/gate-2 (a full certified-path re-solve). No sprawling suite run — the branch structure makes the nominal path provably untouched by construction.
- [ ] **Step 5: Commit** (`E2B/core/lt_mee_rhs.m` + `E3B/test_lt_mee_rhs_pert.m`; message `feat(mee-core): opt-in lunar third-body term in lt_mee_rhs (cr3bp-geo T2)` + trailer).

---

### Task 3: `sanity_bound.m` — the §7 null model from code

**Files:**
- Create: `E3B/sanity_bound.m` (script, campaign-style)

- [ ] **Step 1: Certified per-rung numbers come from `E2B/reproduce/table3_certified.m`** (pure lookup, rungs 10/5/2.5/1/0.5/0.2/0.1 N — reviewer-confirmed; read its field names). Unit-constant note (review amendment G): the implementation's Earth GM is `par.muKm3s2 = 398600.47` (the certified solver's constant) — the spec's 398600.4418 differs by 7e-8 relative; the solver constant governs, documented in `lunar_params`'s header.
- [ ] **Step 2: Write `sanity_bound.m`**: full-header script that, for each certified rung (10, 5, 2.5, 1, 0.5, 0.2, 0.1 N — take the set actually in the table): computes thrust authority `T/m0` in m/s², lunar tide `2*muM*r/DM^3 * AU_ms2` at r = 1 LU, the ratio, t_f in days (1.5x each rung's `table3_certified` tfmin), and t_f in lunar months (`/27.32 d`). Print an aligned table and write it to `E3B/results/sanity_bound.md` (create `results/` + `.gitignore` with `*.mat`). Assert the 10 N ratio is between 0.05% and 0.5% and the 0.1 N ratio between 5% and 20% (loose brackets on the spec §7 predictions — if these fire, the spec's null model is materially wrong: STOP and report, per spec §9 honesty rule).
- [ ] **Step 3: Run it** (no solves; seconds). Table prints, asserts pass.
- [ ] **Step 4: Commit** (`sanity_bound.m` + `results/sanity_bound.md` + `results/.gitignore`; message `feat(cr3bp-geo): sanity-bound null model table (phase1 T3)` + trailer).

---

### Task 4: `bridge_mu_continuation.m` — 10 N energy bridge (gates 1–3)

**Files:**
- Create: `E3B/bridge_mu_continuation.m`

**Interfaces:**
- Produces: `out = bridge_mu_continuation(opts)` with `opts.thrustN [10]`, `.phi0 [0]`, `.gainSched [[0.25 0.5 0.75 1.0]]`, `.maxIter [1500]`, `.resume [true]`. Saves `E3B/results/energy_cr3bp_T<mN|N-tag>_phi<..>.mat` holding `sigma, X, U, dL, tfTarget, fp` (fp incl. `thrustN,m0kg,ispS,tfTarget,muM,DM,nM,phi0,gain`). Returns the final solver `out` + `.gainReached`.

- [ ] **Step 1: Write the driver.** Structure (full header per house style; body):
  1. `setup_paths`; `par = kepler_lt_params(opts.thrustN, 1500, 2000)`; `cert = table3_certified(opts.thrustN)` (the EXISTING certified lookup at `E2B/reproduce/table3_certified.m` — reviewer-confirmed; no grep, no local table); `tfTarget = 1.5*cert.tfmin` (read that file's field names first and use them verbatim); `xf = [1;0;0;0;0];` (review amendment D — the seed call below needs `xf(1)`).
  2. Seed: mirror `run_transfer_mee.m` lines 132–161 VERBATIM (the two-pass protocol: cheap `N=50, stopP=xf(1)` revs probe with its `nRev` window assert, then the full-density `N = round(nodesPerRev*infoP.nRev)` sample) including its `seedThr`/`betaMode`/`nodesPerRev`/`initElems` values for the 10 N rung — read that driver's config block for the rung's actual values; both `mee_seed` calls take the full opts struct shown there.
  3. **Gate-1 baseline (pert absent):** solve `casadi_lt_mee(sigma, X0, U0, dL0, struct('par',par,'mode','fixedtf','eps',1,'tfTarget',tfTarget,'x0',X0(:,1),'maxIter',opts.maxIter,'warmTight',false))`; require `Solve_Succeeded && maxDefect < 1e-6`.
  4. **Gate-2 (gain=0 with pert PRESENT):** `par2 = par; par2.pert = lunar_params(par,opts.phi0,0);` warm-start tight from step 3's X/U/dL; require success AND `max(abs(X - X_gate1))` < 1e-8 (the hook itself introduces no drift). Print both gate verdicts.
  5. **Gain walk:** for each `g` in `gainSched`: `par2.pert.gain = g`; tight solve warm from previous; gate as in step 3 (the explicit `strcmp` four-metric form); on failure halve the step toward the last good gain (insert midpoint into the schedule; floor step 0.05, else error `bridge:stuck` with the honest gain floor). Per-step checkpoint save (resume-safe: cache ONLY accepted steps — a failed step is never cached, so resume can never skip past a failure; fp-check via the campaign's `check_cache_fp` pattern replicated locally with the fp fields above; all tags live in a fresh `cr3bp_*` namespace so no legacy cache can ever match — review amendment on GPT#7).
  6. Save the final artifact + print `BRIDGE: gain=1 reached, defect=..., m_f(energy)=... kg`.
- [ ] **Step 2: Lint** (`checkcode` clean) then **run at 10 N** (one-line batch; 10 N solves are minutes each, ~6 solves total — run in foreground with a generous timeout). Expected: both gates pass; gain walk reaches 1.0 with tiny per-step movement (§7: perturbation ~0.1% of authority at 10 N).
- [ ] **Step 3: Commit** (driver only — the .mat stays untracked; message `feat(cr3bp-geo): mu-continuation bridge driver + 10N energy bridge run (phase1 T4)` + trailer). Record the gate/defect numbers in the commit body.

---

### Task 5: `solve_cr3bp_minfuel.m` — ε-sharpen at gain=1 (gate 3 complete)

**Files:**
- Create: `E3B/solve_cr3bp_minfuel.m`

**Interfaces:**
- Produces: `best = solve_cr3bp_minfuel(opts)` (`opts.thrustN [10]`, `.phi0 [0]`, `.maxIter [1500]`): loads the Task-4 artifact (fp-checked), runs `homotopy_mee(sigma, X, U, dL, ho)` with `ho = struct('par', parWithPertGain1, 'x0', X(:,1), 'tfTarget', tfTarget, 'maxIter', opts.maxIter, 'resDir', <E3B/results/homotopy>, 'tag', sprintf('cr3bp_%s_phi%s', <rung tag>, <phi tag>), 'fp', <the same fp struct>)`; asserts the FULL certification block (review amendment E): `best.certified && best.epsReached==0 && strcmp(best.ipoptStatus,'Solve_Succeeded') && best.maxDefect<1e-6 && best.maxUnit<1e-8 && best.termErr<1e-8` and no `boundSaturation` warning fired; PLUS prints an explicit CAVEAT line: the 2-body PMP/primer verifier (`mee_primer_switch`) is NOT valid under lunar gravity without subtracting the zero-throttle ballistic dXdt (reviewer finding) — a CR3BP-aware primer check is a recorded TODO, not silently skipped; saves `E3B/results/minfuel_cr3bp_<tags>.mat` (best + provenance) and prints `m_f_kg`, switches, edge, defect beside the 2-body certified m_f for the rung (from the Task-4 `tf_table`, extended with the table's m_f column).

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
- **Type consistency:** `pert` struct fields identical across `lunar_params` (producer), `lt_mee_rhs` (consumer), fp records, and both drivers; certified numbers come from the single existing `table3_certified.m` in all of T4/T5/T6 (amendment D — no local table exists to drift).
- **Known judgment point (recorded):** the T2 rewrite redistributes `(Tm/m)` factors into `fR/fT/fN`; the plan requires inspection-equivalence plus the oracle test plus the 2-body no-solve suite as the tripwire. If any 2-body test moves: STOP.

## External-review amendment log (2026-07-22, GPT-5.6-terra + Gemini 3.1 Pro; both verdicts: execute with named amendments)

Applied (confirmed findings): **A** branch structure — nominal path literally
untouched (both reviewers: FP non-associativity breaks bitwise claims);
**B** test fixes — bitwise gain-0 assert, cross-branch continuity at 1e-12
gain, Ldot-equality asserted before the delta conversion; **C** setup_paths
delegates to the 2-body campaign's own (lib/optdef.m coverage) via the cwd
trick, tests use `addpath+call` not `run()`; **D** `table3_certified.m` as
the single certified-numbers source + verbatim two-pass seed mirror + explicit
`strcmp('Solve_Succeeded')` four-metric gates; **E** full certification block
in T5 + the recorded caveat that `mee_primer_switch` is invalid under lunar
gravity without zero-throttle-dXdt subtraction (CR3BP-aware primer check =
named TODO); **F** bounded regression set (branch diff + checkcode + T4
solver gates) replacing the vague suite instruction; **G** Earth-GM constant
note (solver's 398600.47 governs).

Rejected (with reasons): GPT#11 "don't push main" — this repo's established
solo-research workflow is user-directed commit+push to main; kept. GPT#7's
full fail-closed rework of `homotopy_mee`'s existing cache machinery — out of
this plan's scope; the fresh `cr3bp_*` tag namespace makes legacy-cache
contamination impossible, and bridge checkpoints cache accepted steps only.
Host verification: the Nhat convention was independently confirmed numerically
(5 random 3D states: unit, orthogonal, matches FD orbit normal; scratchpad
check_nhat.m) before adjudication.
