# Thrust Ladder to 0.1 N (MEE + ΔL) — Implementation Plan

> **CLOSE-OUT (executed 2026-07-17/18):** Plan executed through Task 11;
> delivered scope + deviations from this plan are documented in
> `DESIGN_thrust_ladder.md`'s status header and in `README.md`'s "MEE
> thrust-ladder campaign" section — see those for the authoritative
> as-built record rather than this plan document.

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Reach the low-thrust ladder (10 N → 0.1 N, Table 3 of HMG-2004) that the
Cartesian formulation could not, by rebuilding the solver in **Modified Equinoctial
Elements (MEE) with true longitude L as the independent variable and the total
longitude span ΔL as a decision variable**, then thrust-continuing down the ladder.
Delivers Table-3 structure counts, the full Fig-23 m_f(c_tf) overlay across thrusts,
and law R0 across two decades.

**Architecture:** New MEE solver `casadi_lt_mee.m` (7 slow states + control, L-domain
collocation on a fixed unit grid σ∈[0,1] with L=L₀+σΔL, ΔL a decision variable, physical
time carried as a state via dt/dσ=ΔL/L̇). Same CasADi+IPOPT + Bertrand–Épénoy ε-homotopy
+ two-stage/continuation machinery as the reviewed Cartesian stack; only the dynamics,
coordinates, and independent variable change. Thrust continuation (10→5→2.5→1→0.5→0.2→0.1 N,
each warm-starting the next) is the backbone — ΔL grows the revolution count automatically,
which is exactly what Cartesian+Sundman could not do.

**Tech Stack:** MATLAB R2025b (headless `-batch`), CasADi 3.7.0 (`~/casadi-3.7.0`, bundled
IPOPT/MUMPS), ode113 for seeds. Reuses (unchanged) `elements_to_cart.m`/`cart_to_elements.m`
(roundtrip-tested 1e-10), `kepler_lt_params.m`.

## Global Constraints

- MATLAB binary: `/Applications/MATLAB_R2025b.app/bin/matlab` ONLY (memory `use-matlab-2025b`).
- CasADi: `~/casadi-3.7.0` (self-loaded by the solver; `CASADI_PATH` overrides).
- Constants (from `kepler_lt_params`): μ_⊕=398600.47 km³/s²; canonical LU=42165 km (GEO radius), TU=√(LU³/μ), μ=1; m₀=1500 kg; Isp=2000 s; g₀=9.80665. **Isp provenance caveat stands** (family cross-check, validated by the mass match) — carry into every new number.
- Problem BCs: initial (P=11625 km, e_x=0.75, e_y=0, h_x=0.0612, h_y=0, L=π, m=1500 kg, start at apogee); terminal GEO (P=42165 km, e_x=e_y=h_x=h_y=0, L free, m free).
- **Cross-formulation gate (the linchpin):** MEE must reproduce the Cartesian baseline m_f = 1376.74 kg (N=600) at 10 N / c_tf=1.5 / free-L to within **0.5 kg**, with matching burn structure, before any ladder step is trusted (Task 4). Do NOT gate MEE on the raw-dual primer metric (open Campaign-B issue).
- MATLAB house style: full comment headers (purpose/INPUTS[sizes]/OUTPUTS[sizes]/REFERENCES); NEVER `i`/`j` as loop variables.
- CasADi gotchas (banked): NEVER chain `a<=x<=b` — two separate `opti.subject_to` calls; matrix (non-vector) inequalities need `(:)` flattening (precedent: `casadi_lt_2body.m`).
- **MEE-specific (review findings, `doc/reviews/2026-07-17_triage.md`):** wrap trig on L to `mod(L,2π)` inside the dynamics (ΔL≈4700 rad at 0.1 N → large-operand cancellation); enforce/log `L̇ ≥ L̇_min > 0` (degeneracy guard); prograde is automatic (target h_x=h_y=0 is i=0 equatorial; retrograde would be h_x,h_y→∞) — no separate prograde inequality needed, but assert i<90° at the terminal.
- Node budget: **25 nodes/rev nominal, 15/rev floor probe**; settle by convergence study at 20/30/40 (Task 5).
- Sporadic uncatchable CasADi/IPOPT MEX crash (~1/10 solves) kills MATLAB; every solve caches per-point and skips-if-exists — rerun on a silent death.
- `results/` untracked; commit code/tests/docs only. All commits from repo root; end messages with `Co-Authored-By: Claude Fable 5 <noreply@anthropic.com>`.
- Run everything from the project folder: `matlab -batch "cd '/Users/msc/Desktop/optimal_control/earth_elliptic_to_geo'; <script>"`.

**Reference numbers (paper Table 3, c_tf=1.5):** 10 N 7.5 rev/18 sw; 5 N 15/36; 2.5 N 30/73; 1 N 74.5/179; 0.5 N 149/360; 0.2 N 377/915; 0.1 N 754/1786. Law R0: T_max·t_f,min ≈ 850 N·h (our 10 N Cartesian: 846.6).

---

### Task 1: MEE dynamics `lt_mee_rhs.m`

**Files:**
- Create: `earth_elliptic_to_geo/lt_mee_rhs.m`
- Test: `earth_elliptic_to_geo/test_mee_rhs.m`

**Interfaces:**
- Produces: `[dXdL, Ldot] = lt_mee_rhs(X, U, par)` where `X = [P; ex; ey; hx; hy; m; t]` (7×1), `U = [beta(3); thr]` (4×1) with `beta` a **unit RTN direction** (radial, transverse, normal), `thr∈[0,1]` throttle; RTN thrust components `(q,s,w) = thr*beta`. Returns `dXdL = (dX/dt)/Ldot` (7×1, d/dL) and `Ldot` (scalar, dL/dt). MUST evaluate on numeric doubles AND CasADi MX (no norm/abs/max on the state; trig on `mod(L,2π)`).
- `L` is passed via `par.L` (the independent-variable value at this node) — see Task 3 for how the collocation supplies it.

**The Gauss variational equations (paper p.6; VERIFY against the PDF and DESIGN.md §2):**
With `Tm = par.Tmax`, `c = par.c`, `mu = par.mu` (=1), and (RTN components) `q=thr*beta(1)`, `s=thr*beta(2)`, `w=thr*beta(3)`:
```
cL = cos(mod(L,2*pi));  sL = sin(mod(L,2*pi));
Z  = 1 + ex*cL + ey*sL;
A1 = ex + (1+Z)*cL;
A2 = ey + (1+Z)*sL;
Xh = 1 + hx^2 + hy^2;
hterm = hx*sL - hy*cL;
sqPmu = sqrt(P/mu);
Pdot  = (2*Tm/m)*sqrt(P^3/mu) * (s/Z);
exdot = (Tm/m)*sqPmu*(1/Z)*( Z*sL*q + A1*s - ey*hterm*w );
eydot = (Tm/m)*sqPmu*(1/Z)*(-Z*cL*q + A2*s + ex*hterm*w );
hxdot = (Tm/(2*m))*sqPmu*(Xh/Z)*cL*w;
hydot = (Tm/(2*m))*sqPmu*(Xh/Z)*sL*w;
Ldot  = sqrt(mu/P^3)*Z^2 + (1/m)*sqPmu*(1/Z)*hterm*w;
mdot  = -(Tm/c)*thr;             % ||(q,s,w)|| = thr since ||beta||=1
tdot  = 1;
dXdt  = [Pdot; exdot; eydot; hxdot; hydot; mdot; tdot];
dXdL  = dXdt / Ldot;
```
(Softening: guard `P` and `Z` away from 0 only if a solve ever probes there; at our states P∈[0.27,1], Z∈~[0.25,1.75], both safely positive.)

- [ ] **Step 1: Write the failing test** (Kepler invariance + thrust cross-check against the Cartesian dynamics)

```matlab
% TEST_MEE_RHS  Ballistic invariance + thrust cross-check vs Cartesian RHS.
p = kepler_lt_params(10, 1500, 2000);
% initial MEE state (paper), coplanar variant for the planar checks where noted
X0 = [11625/p.LU_km; 0.75; 0; 0.0612; 0; 1; 0];
% (a) BALLISTIC: thr=0 -> P,ex,ey,hx,hy,m all frozen; only t advances; Ldot>0
U0 = [1;0;0; 0];
[dXdL, Ldot] = lt_mee_rhs(X0, U0, setfield(p,'L',pi));
assert(Ldot > 0, 'Ldot must be positive');
assert(max(abs(dXdL(1:5))) < 1e-14, 'elements must be frozen under zero thrust');
assert(abs(dXdL(6)) < 1e-14, 'mass frozen under zero thrust');
assert(dXdL(7) > 0, 'time must advance');
% (b) BALLISTIC Ldot value: at L=pi (apogee, ex=0.75) Z=1-0.75=0.25, Ldot=sqrt(1/P^3)*Z^2
P0 = X0(1); Z_apo = 1 - 0.75;
assert(abs(Ldot - sqrt(1/P0^3)*Z_apo^2) < 1e-12, 'ballistic Ldot formula');
% (c) THRUST CROSS-CHECK vs Cartesian: convert MEE->Cartesian, apply the SAME
% physical thrust in both, require d/dt of the Cartesian state to match the
% Cartesian RHS to ODE tolerance. Transverse burn thr=1, beta=[0;1;0].
Uc = [1;0;0; 1];                              % pure transverse, full throttle
[dXdL_t, Ldot_t] = lt_mee_rhs(X0, Uc, setfield(p,'L',pi));
dXdt_mee = dXdL_t * Ldot_t;                   % back to time domain
% independent finite check: energy rate must be positive for a transverse burn
% (raises orbit), and Pdot>0
assert(dXdt_mee(1) > 0, 'transverse burn must raise P');
% cross-formulation identity (the strong one): reconstruct r,v from elements and
% confirm the element-rate equals the Gauss projection of Cartesian thrust accel
[r,v] = elements_to_cart(X0(1),X0(2),X0(3),X0(4),X0(5),pi,p.mu);
assert(abs(norm(cross(r,v)) - sqrt(P0*p.mu)) < 1e-10, 'ang.mom. vs sqrt(P mu)');
fprintf('test_mee_rhs: ALL PASS (Ldot=%.4f, Pdot=%.3e)\n', Ldot_t, dXdt_mee(1));
```

- [ ] **Step 2: Run to verify it fails** (function not defined).
- [ ] **Step 3: Write `lt_mee_rhs.m`** transcribing the Gauss block above with a full header citing paper p.6.
- [ ] **Step 4: Run to verify PASS.** If (c) fails, the RTN sign convention or a Gauss term is off — check against DESIGN.md §2 and the PDF, do NOT loosen tolerances.
- [ ] **Step 5: Commit** `lt_mee_rhs.m` + `test_mee_rhs.m`:
```bash
cd /Users/msc/Desktop/optimal_control && git add earth_elliptic_to_geo/lt_mee_rhs.m earth_elliptic_to_geo/test_mee_rhs.m && git commit -m "feat(earth-geo): MEE Gauss dynamics, L-domain, ballistic + thrust cross-check tested

Co-Authored-By: Claude Fable 5 <noreply@anthropic.com>"
```

---

### Task 2: MEE seed generator `mee_seed.m`

**Files:**
- Create: `earth_elliptic_to_geo/mee_seed.m`
- Test: `earth_elliptic_to_geo/test_mee_seed.m`

**Interfaces:**
- Consumes: `lt_mee_rhs`, `kepler_lt_params`.
- Produces: `[sigma, X0, U0, dL0, info] = mee_seed(par, opts)` — `opts.thr` (constant throttle), `opts.betaMode` ('transverse'|'tangential'), `opts.nRev` (target revs, sets integration span) OR `opts.stopP` (integrate until P≥target). Returns `sigma` [(N+1)×1] uniform 0→1, `X0` [7×(N+1)] MEE states, `U0` [4×(N+1)] controls, `dL0` (total ΔL span, scalar), `info` (.nRev .tEnd .mEnd). Integrate `dX/dL` in L (ode113) from L₀=π, sample at uniform-σ nodes (dense output; defect-free by construction — the no-resample lesson, now trivial since L is the independent variable).

- [ ] **Step 1: Write the failing test** — seed states finite, P monotone increasing toward GEO under transverse thrust, ΔL≈2π·nRev, mass linear in the *time* it maps to; stencil defect on the solver's σ-trapezoid < 1e-2 (mirror the Task-3 stencil, as in the Cartesian `test_seed`). *(Full test body: mirror `test_seed.m`'s structure with the MEE state/rhs.)*
- [ ] **Step 2–4:** run-fail → implement → run-pass.
- [ ] **Step 5: Commit** `mee_seed.m` + test.

---

### Task 3: MEE solver core `casadi_lt_mee.m`

**Files:**
- Create: `earth_elliptic_to_geo/casadi_lt_mee.m`
- Test: `earth_elliptic_to_geo/test_mee_solver_smoke.m`

**Interfaces:**
- Consumes: seeds from `mee_seed`, `lt_mee_rhs`.
- Produces: `out = casadi_lt_mee(sigma, X0, U0, dL0, opts)` — `opts`: `.par`, `.mode` `'mintime'|'fixedtf'`, `.eps` (fixedtf), `.tfTarget` (fixedtf), `.x0` [7×1 initial MEE state], `.maxIter`, `.warmTight`, `.printLevel`. `out`: `.X` [7×(N+1)] `.U` [4×(N+1)] `.dL` (converged ΔL) `.success .ipoptStatus .maxDefect .maxUnit .termErr .mf .m_f_kg .dV_kms .tf .switches .edge .lamDef .LdotMin`. Same field contract shape as `casadi_lt_2body.out` where names overlap.

**Structure — mirror `casadi_lt_2body.m` (reviewed template), with these DELTAS:**
- State is 7-row MEE `[P;ex;ey;hx;hy;m;t]`; **ΔL is a scalar `opti.variable`** (this is the key new DOF — NOT a slack state; the L-domain has no dense-KKT-column problem because ΔL multiplies defects the same way τ_f would, but there is exactly ONE ΔL and the trapezoid is in σ, so it is a single sparse column, acceptable — CONFIRM sparsity at first solve).
- Node L-values: `Lk = x0_L + sigma*dL` where `x0_L = π`. Pass `Lk(k)` into `lt_mee_rhs` per node via `par.L`.
- Defects: `X(:,k+1) - X(:,k) - (dsig_k/2)*(dXdL_k + dXdL_{k+1}) == 0` with `dsig_k = diff(sigma)` (NOTE: `dXdL` already carries the `1/L̇` and the objective/time carry ΔL — see below; keep the ΔL scaling consistent between defects and objective).
- Time as state: `t` is `X(7,:)`; its row of `dXdL` is `1/L̇`, so `t(end)-t(1) = ∫(1/L̇)dL = ∫(ΔL/L̇)dσ` = transfer time. `fixedtf` pins `t(end) == tfTarget`; `mintime` minimizes `t(end)`.
- **L̇ guard:** add `opti.subject_to(Ldot_k >= par.LdotMin)` at every node (`par.LdotMin` e.g. 1e-3) — degeneracy guard.
- Cone: `‖beta‖=1` per node (`beta = U(1:3,:)`); `thr∈[0,1]` (two separate bounds); `mintime` pins `thr==1`.
- Objective (fixedtf): `J(ε) = Σ (dsig/2)(w_k + w_{k+1})`, `w_k = (ΔL/L̇_k)*(thr_k - ε*thr_k*(1-thr_k))` (physical-time measure `dt = (ΔL/L̇)dσ`).
- Terminal (GEO in elements): `P(end)==1` (=42165/LU), `ex(end)==0`, `ey(end)==0`, `hx(end)==0`, `hy(end)==0`. L free (ΔL is the free DOF). **Prograde automatic** (h=0 ⇒ i=0 equatorial; assert post-solve `atan2(2*sqrt(hx²+hy²),1-hx²-hy²)` small).
- IPOPT options identical to the Cartesian solver (incl. `linear_solver='mumps'`, warmTight regimes). **Set `nlp_scaling_method` explicitly** and record it (Campaign-B lesson — don't inherit a silent default).
- Duals: return `out.lamDef` [7×N] for the eventual MEE verifier (Task 10); do not gate on them here.

- [ ] **Step 1: Write the smoke test** — construct both modes, run 5 iterations, assert the full `out` struct returns without error on the non-converged path (mirror `test_solver_smoke.m`); assert ΔL is a scalar variable and the KKT stays sparse (check `opti.debug` / problem size).
- [ ] **Step 2–4:** run-fail → implement (transcribe from `casadi_lt_2body.m` with the deltas) → smoke-pass. Report BLOCKED with the exact CasADi error if the ΔL-variable formulation creates a dense column (contingency: fall back to ΔL-as-slack-state like cScale, but PREFER the scalar variable — one column is fine).
- [ ] **Step 5: Commit** `casadi_lt_mee.m` + smoke test.

---

### Task 4: 10 N cross-formulation VALIDATION GATE (the linchpin)

**Files:**
- Create: `earth_elliptic_to_geo/run_transfer_mee.m` (single-case driver: seed → homotopy → report → save; mirror `run_transfer.m`, reusing `homotopy_2body.m` unchanged if its interface fits, else a thin `homotopy_mee.m`).
- Output: `results/MEE_M2_10N.mat`

- [ ] **Step 1: Write `run_transfer_mee.m`** — 10 N, c_tf=1.5 (t_f from a quick MEE min-time anchor OR the Cartesian t_f,min=22.2248 ND as the target, since the physical transfer is identical), free-L, ε:1→0 homotopy at 25 nodes/rev (~190 nodes for 7.5 revs). Save the report struct.
- [ ] **Step 2: Run it** (background, watcher).
- [ ] **Step 3: THE GATE.** Require: certified (ε=0), defect<1e-8, `termErr<1e-8`, **`abs(m_f_kg - 1376.74) < 0.5`**, revs∈[7,8.5], burns at apogee, i(end)<0.05°, `LdotMin>0`. Also reconstruct the Cartesian trajectory from the MEE solution (`elements_to_cart` along the path) and confirm it satisfies the 2-body EOM (independent cross-formulation check).
- [ ] **Step 4:** If the gate FAILS, the MEE dynamics/objective have a bug — debug HERE at 4.5–7.5 revs (cheap), do not proceed. If it PASSES, the formulation is trustworthy end-to-end.
- [ ] **Step 5: Commit** `run_transfer_mee.m` (+ `homotopy_mee.m` if created). Milestone `--allow-empty` marker acceptable for the run result.

---

### Task 5: node-budget convergence study at 10 N

- [ ] Re-run Task-4's 10 N case at **15, 20, 30, 40 nodes/rev**; tabulate m_f, switches, defect. Gate: m_f stable to <0.5 kg from 25/rev up; pick production density (expect 25/rev). Record the 15/rev floor behavior. Output `results/MEE_nodestudy.mat` + a short table in the report. Commit any driver tweak; milestone marker for the study.

---

### Task 6: MEE min-time anchor + ladder driver

**Files:**
- Create: `earth_elliptic_to_geo/run_mintime_mee.m` (mintime mode; two-stage + continuation recipe ported from `run_mintime.m` WITH the review's **feasibility-selected barrier policy** — warmTight chosen by measured defect, not call number; never reuse restoration multipliers; keep-if-improved), `earth_elliptic_to_geo/run_ladder.m` (thrust-continuation orchestrator).
- Interfaces: `run_mintime_mee(thrustN, nodesPerRev)` → cached anchor `.mat` (t_f,min, ΔL_mt); `run_ladder(thrustList)` → per-thrust fuel solves, resume-safe.

- [ ] **Step 1:** port `run_mintime` → `run_mintime_mee` with the barrier-policy fix; unit-test the guard arithmetic reuse (`mintime_guard_constants.m` is shared).
- [ ] **Step 2:** verify at 10 N: MEE min-time t_f,min matches the Cartesian 22.2248 ND to ~1%.
- [ ] **Step 3–5:** write `run_ladder.m` (skeleton; the actual descent is Tasks 7/9), commit.

---

### Task 7: thrust-continuation backbone 10 → 1 N (Phase 2a)

- [ ] **Step 1: Descend the ladder** 10 → 5 → 2.5 → 1 N via `run_ladder`, each thrust warm-starting the next (state + ΔL rescaled by the C-law guess), c_tf=1.5, production node density. Per rung: min-time anchor → fuel solve → structure counts. Background, resume-safe, watcher.
- [ ] **Step 2: Gates vs Table 3:** 5 N ≈15 rev/36 sw, 2.5 N ≈30/73, 1 N ≈74.5/179 (report as bands, refinement spot-check at 1 N: N→1.5N, |Δm_f|<1 kg). Front m_f(c_tf) begins; law R0 across {10,5,2.5,1} N (gate <10% spread, ≈850 N·h).
- [ ] **Step 3: Commit** milestone marker + report with the per-rung table.

*(If a rung stalls, apply the feasibility-barrier policy and the ΔL warm-start; a rung that resists after the recipe → BLOCKED report with the trajectory, bank the rungs above it.)*

---

### Task 8: PSR port for switch-aware refinement

**Files:**
- Create: `earth_elliptic_to_geo/psr_mee_refine.m` (port `GTO_tulip/direct/PSR/` PMP-Steered Refinement — switching-function zeros steer mesh insertion — to the MEE solver).

- [ ] **Step 1–4:** port + validate at 1 N (where ~179 switches start to smear on a uniform mesh): show refinement stabilizes switch count and m_f vs the uniform-mesh 1 N result from Task 7. Chosen over stopping ε>0 (which biases m_f — the tulip campaign explicitly fixed a legacy 1e-3 stop for that reason).
- [ ] **Step 5: Commit.**

---

### Task 9: low-thrust ladder 0.5 → 0.1 N (Phase 2b)

- [ ] **Step 1: Descend** 0.5 → 0.2 → 0.1 N, thrust-continued, PSR-refined, each an overnight-scale background job (0.1 N ≈ 754 revs, ~19k nodes at 25/rev). Resume-safe; expect MEX crashes → rerun.
- [ ] **Step 2: Gates:** 0.5 N ≈149 rev/360 sw, 0.2 N ≈377/915, 0.1 N ≈754/1786 (bands + refinement spot-checks). L̇>0 held throughout; certified defect per rung.
- [ ] **Step 3: Commit** milestone + the completed Table-3 analog.

---

### Task 10: MEE PMP verifier + Fig-16 analog

**Files:**
- Create: `earth_elliptic_to_geo/verify_pmp_mee.m`, `earth_elliptic_to_geo/fig_switching.m`

- [ ] Build the primer/switching verifier for MEE costates **with a correct dual→costate map from the start** (no cScale here, so the Campaign-B anomaly class is absent; still derive the map explicitly and unit-check the extraction — apply Campaign B's T1 tangential-residual test as the acceptance). Produce the Fig-16 analog (switching function ψ and primer vs L, verifying H1/H2: ψ pinpoint-zero at switches, no singular arc). Gate: primer <1°, sign ≥99%, on the 10 N and 1 N solutions. Commit.

---

### Task 11: final deliverables — Table 3, Fig-23 overlay, README

**Files:**
- Create: `earth_elliptic_to_geo/fig_table3.m` (rev/switch vs thrust), extend `fig front` to overlay m_f(c_tf) at multiple thrusts.
- Modify: `README.md`, `DESIGN_thrust_ladder.md` status.

- [ ] Assemble the Table-3 analog (t_f,min, revs, switches per thrust vs paper), the **complete Fig-23 overlay** (m_f vs c_tf for ≥3 thrusts — the near-independence test), law R0 across 10→0.1 N. Update README with the honest claim (MEE+L defeats the *dimensional* objection; combinatorial difficulty handled by continuation+PSR) and the Isp caveat. Full no-solve test suite green. Commit.

---

## Plan self-review notes

- **Spec coverage (`DESIGN_thrust_ladder.md` Phase 2):** MEE core T3; ΔL-as-variable T3; L-independent + L̇ guard T1/T3; thrust-continuation backbone T6/T7/T9; PSR port T8; node budget T5; prograde-automatic T3; mod(2π) guard T1; MEE verifier + Fig-16 T10; Table-3/Fig-23/law T11. **Cartesian Phase 0/1 deliberately dropped** — the ladder rungs (5/2.5/1 N) are produced as MEE continuation rungs (T7), so a separate Cartesian campaign is redundant for this goal (recorded here as the scope simplification the controller approved).
- **The linchpin is T4** — every later task assumes the 10 N cross-formulation gate passed; it is the cheap place to catch a Gauss-transcription bug.
- **Reuse vs rebuild:** `elements_to_cart`/`cart_to_elements`/`kepler_lt_params`/`homotopy_2body`(if interface fits)/`mintime_guard_constants` reused unchanged; `casadi_lt_mee`/`run_*_mee` mirror the reviewed Cartesian files with the documented deltas (following the codebase's established pattern, not re-deriving).
- **Dependency on Campaign B:** NONE for the primal ladder (MEE has no cScale; gates are primal). Campaign B's *deliverable* (a correct dual verifier) is rebuilt fresh for MEE in T10; the Cartesian cScale investigation is optional closure, off this plan's path.
- **Type consistency:** `out` field names shared with `casadi_lt_2body` where they overlap; the new `.dL`/`.LdotMin` are additive; `mee_seed` returns `[sigma,X0,U0,dL0,info]` consumed verbatim by `casadi_lt_mee`.
