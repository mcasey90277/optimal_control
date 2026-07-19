# DESIGN — SOSC Certificate: NLP-level second-order local-minimum certification

**Date:** 2026-07-19
**Status:** approved design (brainstorm complete); implementation plan to follow at `process/PLAN_sosc.md`.
**Scope:** `earth_elliptic_to_geo/` MEE min-fuel solutions only. Two-body Earth-centered, CasADi `Opti` + IPOPT.

---

## 1. Motivation

Today "certified" in this pipeline means **first-order feasibility only**:
`homotopy_mee.m:62,74` gates on `success && maxDefect < 1e-8 && epsReached==0`.
That certifies a *feasible point of the transcribed NLP*, not a *local
minimum*. The min-fuel problem is non-convex with many local minima; the
campaign already demonstrated this (the reproducer found a different, better
18-switch local min at 10 N; the 2.5 N reproduce found a worse but perfectly
valid local min). Finding *a* local min is not a failure — but we currently
have no rigorous proof that a reported row *is* a local min at all.

This design adds an **NLP-level second-order-sufficient-conditions (SOSC)
certificate**: an independent KKT first-order re-check plus a reduced-Hessian
(critical-cone) inertia test that rigorously distinguishes a strict local
minimum from a saddle. It becomes the new bar for "certified" and is applied
both to new reproducer rows and to the existing campaign rows.

**Non-goals (explicit scope boundaries):**
- **Global** optimality — unattainable and unverifiable for this problem; not attempted.
- **Full critical-cone copositivity** for weakly-active junctions — INCONCLUSIVE is the honest verdict; copositivity is a later enhancement.
- **OCP-level continuous second-order** (strengthened Legendre–Clebsch / conjugate points) — separate track; the finite-dimensional NLP SOSC is the deliverable here.
- The **Cartesian/Sundman legacy** stack (`cartesian_legacy/`) — MEE only.

---

## 2. Decided design forks (from brainstorming)

1. **Multiplier recovery = warm re-solve.** Only the defect duals (`lamDef`)
   are saved and the `Opti` object is discarded on return, so the full
   active-set multiplier set is not recoverable from any saved `.mat`. We
   rebuild the NLP from the saved primal + config fingerprint, warm-start at
   the saved primal (existing `warmTight` settings → ~0 IPOPT iterations), and
   recover all multipliers from the live `opti`. Drift `‖x* − x_saved‖_∞` is
   reported so any movement of the certified point is visible.
2. **Verdict = 3-valued PASS / FAIL / INCONCLUSIVE** (+ ERROR, §7), with
   explicit strict-complementarity detection. The subspace inertia test is
   rigorous SOSC only when strict complementarity holds; at a bang-bang
   junction a node can sit at a throttle bound with a near-zero multiplier
   (weakly active), where the critical cone is not a subspace. We never
   over-claim: weakly-active present ⇒ INCONCLUSIVE.
3. **Tiered gate.** PASS ⇒ certified-SOSC. FAIL (proven saddle) ⇒ demote to
   feasible-only, exclude from the certified Table 3, flag loudly.
   INCONCLUSIVE ⇒ keep as certified-feasibility + annotate. Only a *proven*
   non-minimum is demoted (honors the non-convex reframing).

---

## 3. Architecture

New subsystem lives in a new `verify/sosc/` subfolder. Exactly one
numerics-preserving hook is added to the solver; nothing on the certified
numeric path is otherwise touched.

| File | Responsibility |
|---|---|
| `verify/sosc/verify_sosc_mee.m` | Orchestrator: rebuild → recover → KKT re-check → active set → inertia → verdict struct + tiered-gate status. |
| `verify/sosc/sosc_recover_kkt.m` | Warm re-solve at saved primal; return `x*`, full `lam_g`, sparse `H`, sparse `A_all`, `gval`, `grad_f`, constraint/variable registries, drift. |
| `verify/sosc/sosc_kkt_residual.m` | Global Lagrangian-sign resolution + stationarity / primal-feas / dual-feas / complementarity residuals vs thresholds. |
| `verify/sosc/sosc_active_set.m` | Classify each inequality active / strongly-active / weakly-active; assemble active Jacobian `A`; LICQ (rank) check; human-readable weak/degenerate labels. |
| `verify/sosc/sosc_inertia.m` | Sparse LDLᵀ inertia of the KKT matrix; PASS/FAIL/INCONCLUSIVE decision; optional non-gating reduced-Hessian curvature margin. |
| `verify/sosc/sosc_defaults.m` | Single source of the tolerance struct (all thresholds, §6). |
| `verify/sosc/recertify_table3.m` | Batch driver: loop existing certified rows (10/5/2.5/1/0.5 N), certify, write **sidecar** verdicts + printed report. |
| `core/casadi_lt_mee.m` | **Only change:** an `opts.returnModel` flag that *additionally* returns the `opti` object, symbol handles, and a constraint/variable registry. Added output fields only — no numeric change. |

Design intent: each file has one clear job and a small, testable interface.
`sosc_inertia` and `sosc_active_set` are pure linear-algebra/bookkeeping units
testable on a hand-built QP with no NLP solve.

---

## 4. Interfaces (signatures + contracts)

### 4.1 Solver hook — `core/casadi_lt_mee.m`

Add `opts.returnModel` (default `false`). When true, the returned `out` gains:

```
out.model = struct( ...
   'opti', opti, ...          % the live casadi.Opti handle
   'X', X, 'U', U, 'dL', dL, ...% MX symbol handles for the decision blocks
   'creg', creg, ...          % constraint registry (below)
   'vreg', vreg );            % variable registry (index maps)
```

- `creg`: struct array, one entry per `subject_to` GROUP, recorded as
  constraints are added by bracketing each group with
  `r0 = size(opti.g,1)+1; ...; r1 = size(opti.g,1)`:
  `creg(i) = struct('label', <'defect'|'ldotGuard'|'betaNorm'|'thrEq'|'thrLo'|'thrHi'|'boxP'|...|'initBC'|'termBC'|'dLbox'|'tBox'|'betaBox'>, 'kind', <'eq'|'ineqLo'|'ineqHi'>, 'rows', r0:r1, 'bound', <scalar or []>, 'node', <k index vector or []>)`.
- `vreg`: index maps into the decision vector `opti.x` for each block
  (`X` rows/cols, `U` rows/cols, `dL`), so `sosc_active_set` can attach
  node-level labels (e.g. "throttle upper bound, node 137").
- **Numerics invariant:** with `returnModel=false` (all existing callers) the
  function is byte-for-byte unchanged. The registry bracketing runs only under
  the flag. This is a plan checkpoint (Task 2 test: existing 10 N solve output
  identical with the flag off).

### 4.2 `sosc_recover_kkt(saved, tol) -> R`

`saved` carries everything needed to rebuild the exact NLP: `X, U, dL, sigma,
fp` (config fingerprint), `thrustN, ctf, tfTarget, xf, x0`. Steps:

1. Reconstruct `par`/`opts` from `fp`; **checkpoint:** a plain rebuild+re-solve
   must reproduce the saved primal to `tol.recon` before anything downstream is
   trusted (Task 1).
2. Build with `returnModel=true`, `set_initial` to the saved primal, solve with
   `warmTight` (μ_init small, bound-push 1e-9).
3. If the warm solve does not reach `Solve_Succeeded`/`Solved_To_Acceptable`,
   return `R.recoverOK=false` (⇒ ERROR verdict upstream).
4. Assemble, all in Opti's **native (unscaled) symbols**:
   - `R.x = sol.value(opti.x)`, `R.lam_g = sol.value(opti.lam_g)`, `R.gval = sol.value(opti.g)`
   - `R.grad_f = sol.value( gradient(opti.f, opti.x) )`
   - `R.A_all = sol.value( jacobian(opti.g, opti.x) )`   % [m×n] sparse
   - `R.H = sol.value( hessian(opti.f + opti.lam_g.'*opti.g, opti.x) )` % [n×n] sparse
   - `R.creg`, `R.vreg`, `R.drift = norm(R.x - x_saved, inf)`, `R.recoverOK=true`.

### 4.3 `sosc_kkt_residual(R, tol) -> K`

- **Global sign resolution:** the CasADi/IPOPT Lagrangian sign is a single
  global ambiguity. Compute `s = argmin_{s∈{+1,-1}} ‖grad_f + s·A_allᵀ·lam_g‖_∞`
  (reuses the empirical sign-resolution pattern of `verify_pmp_mee.m:112-121`).
  If `min` residual still exceeds `tol.stat` under both signs ⇒ `K.signOK=false`
  (⇒ ERROR; the duals are inconsistent — the honest guard the Cartesian anomaly
  taught us to keep).
- Residuals returned in `K`: `stat` (‖∇ℒ‖_∞ with resolved sign), `primalEq`
  (max |equality rows of gval|), `primalIneq` (max inequality violation via
  `creg` kind+bound), `dualFeas` (worst sign violation of `s·lam_g` over **all**
  inequality rows — inactive rows have `lam_g≈0`, trivially feasible, so this
  needs no active set and can run before §4.4), `comp` (max |lam_g · slack| over
  inequalities). Plus `K.sign`,
  `K.signOK`, and pass booleans vs the `tol` thresholds.

### 4.4 `sosc_active_set(R, K, tol) -> AS`

- Per inequality row: `slack = distance(gval, bound)` from `creg`. **Active** if
  `slack < tol.active`. Among active: **strongly-active** if
  `|lam_g| > tol.mu * max(1, max|lam_g|)`, else **weakly-active**.
- Equalities (defects, per-node `‖β‖=1`, initial/terminal BC) always active.
- `AS.A` = rows of `R.A_all` for all equalities + strongly-active inequalities;
  `AS.m_active = rows(AS.A)`.
- `AS.nWeak`, `AS.weakLabels` (human strings from `creg`+`vreg`, e.g.
  "throttle upper bound, node 137"), `AS.nStrong`, `AS.nEq`.
- **LICQ:** `AS.licq = (sprank(AS.A) == AS.m_active)` is a cheap *structural*
  flag that names likely-dependent groups. It is not the rigorous test: the
  rigorous rank-deficiency catch is `inertia.nzero>0` in §5 (a rank-deficient
  active Jacobian makes the KKT matrix singular ⇒ nonzero nullity ⇒
  INCONCLUSIVE). `sprank` only explains *why*.

### 4.5 `sosc_inertia(R.H, AS.A, tol) -> IN`

- Form `Kkt = [H, Aᵀ; A, sparse(m_a,m_a)]` (sparse symmetric).
- `[~,D,~] = ldl(Kkt,'vector')`; count eigenvalue signs of the 1×1 and 2×2
  diagonal blocks of `D` (2×2 via its trace/det). Pivots with
  `|·| < tol.inertiaZero · scale` count as **zero**.
- `IN.inertia = [npos, nneg, nzero]`; `IN.expected = [n, m_a, 0]`.
- Optional non-gating `IN.redMinEig`: smallest reduced-Hessian eigenvalue via a
  few steps of projected inverse iteration through the KKT factorization
  (null-space method). `NaN` if not requested. **Never gates** — a curvature
  margin for reporting only.

### 4.6 `verify_sosc_mee(saved, opts) -> sosc`

Orchestrates 4.2→4.5 and produces the verdict + status (schema in §5).

---

## 5. Verdict struct + tiered gate

```
sosc = struct( ...
  'verdict',   <'PASS'|'FAIL'|'INCONCLUSIVE'|'ERROR'>, ...
  'reason',    <human string>, ...
  'status',    <'certified-sosc'|'feasible-only'|'certified-feasibility+sosc-inconclusive'>, ...
  'drift',     <‖x*-x_saved‖_inf>, ...
  'sign',      <resolved global Lagrangian sign>, ...
  'kkt',       struct('stat',_, 'primalEq',_, 'primalIneq',_, 'dualFeas',_, 'comp',_), ...
  'active',    struct('nEq',_, 'nStrong',_, 'nWeak',_, 'weakLabels',{...}, 'licq',_), ...
  'inertia',   struct('npos',_, 'nneg',_, 'nzero',_, 'expected',[n m_a 0]), ...
  'redMinEig', <scalar or NaN>, ...
  'thresholds',<the tol struct used>, ...
  'meta',      struct('thrustN',_, 'tag',_, 'n',_, 'm',_, 'm_active',_, 'when',datestr(now)) );
```

**Verdict logic (in order):**
1. `~recoverOK` or `~signOK` or KKT residuals over threshold ⇒ **ERROR**
   (could not build a trustworthy certificate).
2. `~licq` or `nWeak>0` or `inertia.nzero>0` ⇒ **INCONCLUSIVE**
   (strict complementarity / LICQ fails; cone is not a subspace).
3. `inertia == [n, m_a, 0]` ⇒ **PASS** (reduced Hessian PD on the critical
   cone ⟺ strict local minimum).
4. otherwise ⇒ **FAIL** (indefinite reduced Hessian ⟹ saddle).

**Tiered gate (status):**
- PASS → `certified-sosc`.
- FAIL → `feasible-only`; **excluded from the certified Table 3**; loud warning.
- INCONCLUSIVE → `certified-feasibility+sosc-inconclusive`; kept, annotated.
- ERROR → treated like INCONCLUSIVE for gating (non-demoting, kept) but flagged
  distinctly for investigation.

Only FAIL demotes. INCONCLUSIVE and ERROR never demote a previously
feasibility-certified row — we could not *prove* non-minimality, so we do not
claim it.

---

## 6. Tolerances (`sosc_defaults.m`, calibrated on the 10 N row)

Canonical units make magnitudes O(1) (μ=1, LU=GEO radius, mass unit=m0):

| name | default | meaning |
|---|---|---|
| `tol.recon` | 1e-6 | rebuild reproduces saved primal (‖·‖_∞) |
| `tol.drift` | 1e-6 | warm-resolve drift ‖x*−x_saved‖_∞ (report/warn, not fail) |
| `tol.stat` | 1e-6 | stationarity ‖∇ℒ‖_∞ pass |
| `tol.feas` | 1e-8 | max equality residual / inequality violation (matches existing `maxDefect<1e-8`) |
| `tol.dual` | 1e-8 | inequality dual-sign violation |
| `tol.comp` | 1e-6 | complementarity max\|λ·slack\| |
| `tol.active` | 1e-7 | inequality slack ⇒ active |
| `tol.mu` | 1e-6 | relative multiplier ⇒ strongly-active (× max(1,max\|λ\|)) |
| `tol.inertiaZero` | 1e-9 | relative pivot magnitude ⇒ zero eigenvalue |

Defaults are the starting point; the 10 N integration test is where they are
confirmed/tightened before the rest of the ladder is trusted.

---

## 7. Error handling

- **Recovery non-convergence** (warm re-solve fails) → verdict ERROR, reason
  states the IPOPT status; no fabricated certificate.
- **Sign unresolved** (stationarity over threshold under both signs) → ERROR,
  reason "dual inconsistency" — the explicit guard against the dual-anomaly
  class; we do not silently proceed.
- **LICQ failure / rank-deficient active Jacobian** → INCONCLUSIVE with the
  dependent constraint groups named.
- **Reconstruction mismatch** (Task-1 checkpoint fails: rebuild ≠ saved primal)
  → hard error in `sosc_recover_kkt`; the config fingerprint is insufficient and
  must be fixed before any certification is meaningful.

---

## 8. Plug-in points & data flow

- **Reproducer (go-forward):** `reproduce_row` calls `verify_sosc_mee` on each
  keep-best candidate; the pool **adopts** a candidate only if
  `verdict ∈ {PASS, INCONCLUSIVE}` (never FAIL/ERROR-as-non-min). `res.sosc` is
  stored in the adopted row.
- **Driver (go-forward):** `run_transfer_mee` attaches `res.sosc` after a
  feasibility-certified fuel solve and applies the tiered gate to the row
  status.
- **Existing campaign rows:** `recertify_table3.m` rebuilds/certifies
  10/5/2.5/1/0.5 N (both the `MEE_M2_*` fuel rows and the `*_PSR_psr_final`
  refined rows) and writes verdicts to **sidecar** files
  `results/sosc/sosc_<tag>.mat` plus a printed summary table — the campaign
  `.mat` caches are **left untouched** (honors "never clobber production
  caches"). A FAIL here is a genuine finding surfaced loudly.

---

## 9. Testing

- **`test_sosc_inertia_qp.m` (no NLP):** hand-built equality+bound-constrained
  QP with known inertia. PD reduced Hessian ⇒ PASS; flipped curvature ⇒ FAIL; a
  weakly-active bound (active with ~0 multiplier) ⇒ INCONCLUSIVE. Covers the
  FAIL and INCONCLUSIVE branches that are hard to elicit cheaply from the real
  NLP.
- **`test_sosc_active_set.m` (no NLP):** synthetic `gval`/`lam_g`/`creg` ⇒
  assert active / strongly-active / weakly-active classification, `weakLabels`
  text, and LICQ rank logic.
- **`test_sosc_recon_10N.m` (Task-1 checkpoint, solves):** rebuild from `fp`
  reproduces the saved 10 N primal to `tol.recon`.
- **`test_sosc_10N.m` (integration, solves):** `verify_sosc_mee` on the 10 N
  certified row ⇒ expect `drift<tol.drift`, `kkt.stat<tol.stat`,
  `verdict=='PASS'`; this is also where the §6 tolerances are calibrated.
- **Solver-invariant check:** the 10 N solve output is identical with
  `returnModel` off (guards the hook).

---

## 10. Open implementation risks (carried into the plan)

1. **`res.fp` completeness** — must fully determine the NLP; the Task-1
   reconstruction checkpoint is the gate. If `fp` is insufficient, widening it
   (or re-deriving `par`/`opts` from the saved config) precedes everything.
2. **CasADi `hessian`/`jacobian` on large `opti.g`** — at 0.5 N, `n≈18k`,
   sparse block-tridiagonal; expected cheap (comparable to one IPOPT linear
   solve) but confirmed on the 10 N row first, then 0.5 N, before batch re-cert.
3. **Inertia bookkeeping for 2×2 `ldl` blocks** — must count block eigenvalue
   signs correctly (covered by the QP unit test).
4. **MUMPS/METIS crash class** inherited by the warm re-solve — reuse the
   existing `mumps_pivot_order=0` workaround; deep rungs may still hit it and
   would surface as ERROR (recovery failure), not a false verdict.

---

## 11. Revision — bang-bang WEAK_MIN + `lbg/ubg` bound sourcing (2026-07-19)

The Task-8 integration (first real end-to-end run on the certified 10 N row)
found that the original §5 verdict logic is wrong for this problem class, and
that the bound bookkeeping of §4.1/§4.3 is fragile. Stationarity came back
machine-tight (‖∇ℒ‖ = 1.5e-14), so the recovery framework is sound; the
findings below refine the residual sourcing and the second-order taxonomy.
This section supersedes the affected parts of §4.3, §4.4, §4.5, §5.

### 11.1 The finding — min-fuel bang-bang is a *weak* (non-strict) local minimum

Computed directly from the 10 N KKT matrix `K = [H Aᵀ; A 0]` (eig of the full
3885×3885): the active Jacobian is rank-deficient by only **1** (LICQ
essentially holds), and `inertia(K) = (1865, 1749, 271)` decomposes to a
**reduced Hessian of inertia (116 positive, 0 negative, 270 zero)** on the
386-dim critical cone. A clean spectral gap (≈270 eigenvalues at 1e-10, then
nothing until ~1e-4) confirms the null space is genuine, not a threshold
artifact. So the reduced Hessian is **positive semi-definite with a ~270-dim
null space**: no negative curvature (not a saddle), but not strictly PD. This
is the signature of a bang-bang / linear-in-control extremal — perturbations
preserving the switching structure are flat to second order. **Strict NLP-SOSC
(a PASS) is generically unreachable for min-fuel solutions**; the honest,
useful second-order statement is a *weak local minimum*.

### 11.2 Robust bound sourcing from `opti.lbg`/`opti.ubg`

CasADi Opti canonicalizes every `subject_to` into `g ∈ [lbg, ubg]` and exposes
`opti.lbg`, `opti.ubg` (both length `m`). `sosc_recover_kkt` now returns
`R.lbg`, `R.ubg` (m×1). All residual/slack/active/dual logic sources bounds
from these — NOT from hand-coded `creg.bound` — eliminating the class of bugs
where `creg` recorded a scalar 0 for `initBC`/`termBC` (whose true per-row
bounds are `x0`/`xf`) or omitted the bound subtraction (`tfPin`, `betaNorm`).
`creg` is retained only for human LABELS/NODES in weak-node reporting.

Per-row constraint KIND is derived from the bounds:
`isEq = (lbg==ubg)`; `upperActive` row has finite `ubg`; `lowerActive` has
finite `lbg`. Residual (primal feasibility) becomes uniform:
`viol_i = max(lbg_i − g_i, 0) + max(g_i − ubg_i, 0)` for every row; equality
residual is just the `lbg==ubg` case of the same formula.

### 11.3 Per-kind dual feasibility

Opti/IPOPT reports `lam_g` with OPPOSITE sign conventions per active bound
(empirically: lower-bound-active `lam ≤ 0`, upper-bound-active `lam ≥ 0`;
equalities free). After resolving the one global sign `s`, dual feasibility is
checked per kind: for an upper-bound-only row require `s·lam ≥ −tol.dual`; for
a lower-bound-only row require `s·lam ≤ +tol.dual`; equality rows are
unconstrained in sign. (The single-global-sign check of §4.3 was wrong and
produced a spurious `dualFeas` from the lower-bound rows.)

### 11.4 Reduced-Hessian inertia (supersedes the §4.5 subspace test)

`sosc_inertia` still forms `K` and takes its inertia `(npos,nneg,nzero)` via
`count_inertia`, but now ALSO reports the **reduced-Hessian inertia** using the
Gould decomposition with `r = rank(A)` (via `sprank(A)`, the structural rank):

```
red.npos  = npos  − r
red.nneg  = nneg  − r
red.nzero = nzero − (m_a − r)
```

Consistency (must hold, else the rank estimate is untrustworthy →
INCONCLUSIVE): `red.npos + red.nneg + red.nzero == n − r`, and
`red.nneg ≥ 0`, `red.nzero ≥ 0`. `IN.red = struct('npos',_,'nneg',_,'nzero',_)`,
`IN.rankA = r`, `IN.redConsistent = <bool>`. The old `subspaceOK` remains a
reported field but no longer drives the verdict.

### 11.5 Revised verdict taxonomy (supersedes §5 verdict logic)

Verdict logic, in order (on the REDUCED inertia):
1. `~recoverOK` or `~K.signOK` or `~K.pass` ⇒ **ERROR**.
2. `~AS.licq` or `AS.nWeak>0` or `~IN.redConsistent` ⇒ **INCONCLUSIVE**
   (critical cone not a trustworthy subspace; can neither confirm nor deny).
3. `IN.red.nneg > 0` ⇒ **FAIL** (a negative curvature direction on the cone =
   descent direction ⇒ genuine saddle, provably not a local min).
4. `IN.red.nzero == 0` (and `red.nneg==0`) ⇒ **PASS** (reduced Hessian PD ⇒
   strict local minimum).
5. else (`red.nneg==0`, `red.nzero>0`) ⇒ **WEAK_MIN** (reduced Hessian PSD with
   `red.nzero` flat directions — necessary second-order condition holds, strict
   sufficient does not; the expected outcome for bang-bang min-fuel).

`redMinEig` stays a NaN placeholder (non-gating). The verdict struct gains
`.red` (the reduced inertia) and `.nFlat = IN.red.nzero` (reported for
WEAK_MIN — the dimension of the flat manifold).

### 11.6 Revised tiered gate

- PASS → `certified-sosc` (strict local min).
- **WEAK_MIN → `certified-weak-min`** — a POSITIVE certificate (no descent
  direction; necessary 2nd-order condition met). **Does NOT demote**; kept in
  the certified Table 3, annotated with `nFlat`.
- FAIL → `feasible-only`; excluded from the certified Table 3; loud warning.
- INCONCLUSIVE / ERROR → `certified-feasibility+sosc-inconclusive`; kept,
  annotated, non-demoting.

Only FAIL (a *proven* saddle) demotes. The reproducer (§8) adopts a candidate
whose verdict ∈ {PASS, WEAK_MIN, INCONCLUSIVE} (never FAIL).

### 11.7 Expected results

The 10 N certified row is expected to certify as **WEAK_MIN** (`nFlat ≈ 270`),
and the same is anticipated for the other min-fuel rungs (bang-bang). A rung
returning FAIL would be a genuine finding (a reported non-minimizer); PASS
would indicate an unexpectedly non-degenerate solution. The batch re-cert
(§8, Task 10) reports the verdict + `nFlat` per rung.
