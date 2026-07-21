# Design: PMP-residual-driven adaptive mesh refinement (prototype)

**Date:** 2026-07-11
**Status:** approved (design review with user)
**Goal:** Prototype the combined direct↔indirect move labelled **point (4)** in
the min-fuel campaign discussion: after the direct Sundman solver converges to a
bang-bang solution, use the **indirect (PMP) switching function** — recovered
from the direct solve's own KKT-dual costates — to find *where the optimality
conditions are worst resolved*, refine the collocation mesh there, and re-solve
the direct problem. Demonstrate on the certified 1.15× solution that this
**stabilizes the switch times** (removes their mesh-boundary pinning) at
essentially fixed propellant.

This is a **prototype** to prove the mechanism, not a production hp-adaptive
engine. It reuses the certified `ms_band/` PMP machinery (the dual→costate
mode-'d' map and the Sundman EOM) as a *measurement* tool driving a *direct*
re-solve — no LM/shooting anywhere.

## 1. Background and scope

### The problem this addresses
The direct min-fuel solutions are "certified" only in the discrete-NLP sense
(machine-tight KKT/defects on a fixed mesh). Their bang-bang **switch times are
pinned to σ-node boundaries** — an O(h²), mesh-dependent artifact flagged in
every campaign caveat (`LOW_THRUST_MINFUEL_CAMPAIGN.md`,
`sundman_minfuel/OPTIMALITY_VERIFICATION_PLAN.md` Layer 2, and the
`minfuel-gto-tulip-solved` memory). Point (4) refines the mesh *guided by the
PMP switching function* so switches relocate to where the continuous optimality
condition S(τ)=0 actually places them.

### Why this is the cheap, low-risk item
- The solver `casadi_minfuel_sundman.m` already accepts a **per-interval σ
  mesh** (`dsig = diff(sigma)`, lines 74 & 107). Local refinement = pass a
  denser-in-places `sigma` + a warm start on it. **No solver change.**
- The switching function on the full node grid is already built and validated:
  `verify_direct_pmp.m` computes `S = 1 − ‖λ_v‖c/m − λ_m` from the mode-'d'
  dual→costate map (`sms_seed_duals.m`). We reuse it as the refinement
  indicator.
- The no-resample warm-start discipline (a documented hard requirement —
  resampling a 40-rev trajectory pins IPOPT in restoration,
  `OPTIMALITY_VERIFICATION_PLAN.md` §F.3) is preserved by *inserting* nodes and
  interpolating **only** the inserts, leaving all original node values exact.

### Relationship to the broader plan
This is Layer 2 ("discretization-independent") of
`OPTIMALITY_VERIFICATION_PLAN.md`, driven by a Layer-1 PMP indicator. It does
**not** attempt the indirect *solve* (point 3, a multi-day saltation-aware BVP
build) nor the second-order sufficiency test. It is the incremental win; point
(3) is the structural win.

## 2. Indicator decision (resolved in design review)

**Option 1 (switching-function localization) drives refinement; option 2's
Hamiltonian residual is computed as a passive read-only diagnostic.**

Rationale:
- Switch localization is the most legible demonstration of point (4): "switch #k
  is pinned between nodes j and j+1; S=0 falls at τ*; after refining that
  bracket the switch moved < the new h and stopped moving."
- The Hamiltonian residual |H_σ| = |κ(H_t + λ_t)| is dominated by the same
  switch/perigee locations (λ_t is only ~constant, std ~12% with spikes at
  switches/perigees — `OPTIMALITY_VERIFICATION_PLAN.md` §C check 3), so as a
  *driver* it would largely refine the same places with a noisier signal.
  Perigees are already node-concentrated by the uniform-τ Sundman mesh.
- Escalation is free: the loop is identical; only the indicator function gains a
  term. Computing |H_σ| passively each round makes the escalate/don't-escalate
  decision **data-driven** — if it stays high in un-refined arcs, promote it to
  a driver (that is "option 2", a follow-on, not built here).

## 3. Architecture and data flow

One **round**:

```
converged out (X, U, lamDef) on mesh σ
  → indicator:   duals ─(mode-d)→ S(τ) on node grid; locate S=0 crossings;
                 score each interval by (switch nearby)×(local h)×(root offset)
                 passive diag: |H_σ| per node, switching-law violation count
  → refine:      bisect the top-K scoring intervals (+ immediate neighbors),
                 keeping ALL existing σ nodes  → σ_new,  isNew mask
  → warm start:  original nodes copied EXACTLY; inserted nodes interpolated
                 (pchip on states; throttle by STEP from the correct side of a
                 switch; α interpolated then renormalized)
  → re-solve:    casadi_minfuel_sundman(σ_new, …, ε=0, warmTight=true)
  → measure:     switch times, Δpropellant, violation count, |H_σ|, defect
```

Accept and stop when **all** hold (Layer-2 criterion):
- max switch-time move across the round < that round's local h at each switch,
- |Δpropellant| < 1e-4 kg,
- switch **count** unchanged (none born/dies),

or stop at `maxRounds` (default 4) or if a re-solve fails to converge tight.

### Test target
Certified **1.15×** point (`sundman_minfuel/sundman_minfuel_certified.mat`;
canonical anchor: defect 2.4e-14, 25 switches, 2.2640 kg, ΔV 3.3696 km/s).
**Data dependency:** that `.mat` was saved before dual extraction and may lack
`lamDef`. If absent, regenerate with one ε=0 `warmTight` re-solve
(`casadi_minfuel_sundman`, ~1–3 min) before round 1. A file already carrying
duals (e.g. a `results/minfuel/*_en.mat`) may be substituted via `seedFile`.

## 4. Components

New self-contained folder `GTO_tulip/sundman_minfuel/refine/`.
Each file has the standard commented header (purpose, inputs w/ sizes, outputs
w/ sizes, references).

| File | Purpose | Interface |
|---|---|---|
| `pmp_refine_indicator.m` | Switch-localization score + passive diagnostics from a solution's duals | `(out, sigma, p) → score [1×N], tauSwitch [1×ns], diag{Hres [1×nN], nViol, betaSpread, Snode [1×nN]}` |
| `refine_sigma.m` | Bisect the top-K scoring intervals (+neighbors), preserve all nodes | `(sigma, score, opts) → sigmaNew [(N'+1)×1], isNew [logical (N'+1)×1]` |
| `warmstart_on_mesh.m` | No-resample warm start onto the refined mesh | `(out, sigma, sigmaNew, isNew) → X0 [8×N'+1], U0 [4×N'+1]` |
| `refine_loop.m` | Driver: round loop, acceptance test, history + figure | `(seedFile, opts) → history struct; saves refine_history.mat + .png` |

**Reuse (no reimplementation):**
- Switching function / node costates: factor the small dual→costate core out of
  `sms_seed_duals.m` (mode 'd') so `pmp_refine_indicator` can call it on an
  **in-memory** `out` struct rather than a file. If factoring proves invasive
  for a prototype, fallback is to write a temp `.mat` per round and call the
  existing file-based `sms_seed_duals` — noted so the implementer can choose the
  lighter path.
- `|H_σ|` per node: `sms_eom.m` at each node with the mapped costate.
- Re-solve: `casadi_minfuel_sundman.m` unchanged.
- Constants/params: `cr3bp_lt_params.m`.

### Indicator detail (`pmp_refine_indicator`)
1. Map interval duals → node costates (mode 'd'); form `Snode = 1 − ‖λ_v‖c/m − λ_m`.
2. Zero-crossings of `Snode`: for each, linear-interpolate the root τ* within its
   bracketing interval.
3. Per-interval score = for intervals within a window (±`nbr`, default 3) of a
   crossing: `w = localScale(k)` (the σ-interval width, so coarser brackets
   score higher) × `rootOffset` (distance from the nearer node to τ*,
   normalized by h — a switch mid-interval is worst-resolved). Intervals with no
   nearby crossing score 0.
4. Passive diagnostics: `Hres(n) = |κ·(H_t + λ_t)|` at each node; `nViol` =
   count of nodes where `sign(Snode) ≠ (throttle>0.5)` outside a ±3-node
   switch deadband; `betaSpread` from the single-β switch fit (sanity that the
   dual map is healthy this round).

### Refinement detail (`refine_sigma`)
- Rank intervals by `score`; take the top-K (default `K = min(8, #crossings)`).
- Bisect each selected interval and its immediate neighbors by inserting the
  σ-midpoint. **Guards:** never bisect an interval whose width is below `hFloor`
  (default 1e-9, above the mesh's ~4e-12 near-duplicate spacing); cap node
  growth per round at `+maxAdd` (default 40). Log every interval dropped by a
  guard (no silent truncation).

### Warm start detail (`warmstart_on_mesh`)
- Original σ nodes: copy X,U values **verbatim** (exactness = the no-resample
  discipline).
- Inserted nodes: states by pchip on (σ, X); direction α by pchip then
  renormalized to ‖α‖=1; throttle s by **step** — take the value of the nearer
  original node, and if the insert straddles a switch, the node on the
  pre-switch side — never average across a switch (would seed a smeared
  throttle, the exact failure trapezoidal refinement is meant to cure).

## 5. Outputs and success criteria

`refine_loop` returns and saves a `history` struct (one row per round):
`nNodes, switches, tauSwitch, maxSwitchMove, prop_kg, dProp, nViol, HresMax,
maxDefect, betaSpread, ipoptStatus`, plus a figure:
- top: switch times vs round (should converge — flat tails),
- middle: S(τ) with switch markers, original mesh vs final mesh (the money plot),
- bottom: per-round `nViol`, `dProp`, `HresMax` (the escalation dashboard).

**Prototype succeeds if** on 1.15×: ≥1 round drives every switch-time change
below that round's local h (switch times converge) with |Δprop| < 1e-4 kg and
non-increasing `nViol`; and the passive `HresMax` trend answers whether option 2
is warranted (stays high in un-refined arcs → escalate; drops → don't).

A **negative result is still a result**: if refinement moves switches but
propellant/violations don't improve, or a re-solve MEX-crashes, that is a
documented finding, not a hidden failure.

## 6. Non-goals (YAGNI)

- No Hamiltonian-*driven* refinement (option 2) — only its signal is reported.
- Not wired into `minfuel_at_tf` / `orchestrate` — standalone driver.
- No monitor-function equidistribution — simple top-K bisection (equidistribution
  is the production upgrade).
- No convergence-order / mesh-independence *study* across many t_f; single test
  point.
- No second-order sufficiency test (that rides on point (3)'s switch-time
  parameterization).
- No process-isolation harness — save history each round so a crash is
  recoverable; isolation is a later concern.

## 7. Risks and guards

| Risk | Guard |
|---|---|
| Bisection creates sub-`hFloor` intervals → conditioning | `hFloor` skip + per-round `maxAdd` cap; log drops |
| Throttle smear at inserts near a switch | step throttle from the correct side; never interpolate across a switch |
| Uncatchable MEX crash on re-solve (documented) | persist `history` + latest mesh each round; resumable by hand |
| Seed lacks `lamDef` | regenerate via one ε=0 warmTight re-solve before round 1 |
| Dual map unhealthy on a refined mesh | `betaSpread` reported per round; abort/flag if it blows up |
| Runtime | ~1–3 min/re-solve × ≤4 rounds ≈ ≤20 min; acceptable |

## 8. Environment (from `OPTIMALITY_VERIFICATION_PLAN.md` §E)

- MATLAB **R2025b** (`/Applications/MATLAB_R2025b.app/bin/matlab -batch`);
  R2025a license-broken. Use the `matlab-headless` skill.
- Write a `.m` and run `matlab -batch "cd('<dir>'); script"` — multi-line
  `-batch` strings fail through the shell here.
- CasADi 3.7.0 at `~/casadi-3.7.0` (`CASADI_PATH`); prebuilt arm64 loads, MEX
  compiler broken — don't build.
- IPOPT writes a stray `fort.6` in cwd — ignore/delete.
- `warmTight=true` for re-solving AT a bang-bang point; using the wrong warm
  mode wedges IPOPT.
