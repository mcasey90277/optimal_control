# Design: Indirect Finishing Solve (IFS) — fixed-structure PMP polish of a direct bang-bang solution

**Date:** 2026-07-11
**Status:** approved (design review with user)
**Goal (v1, "Option A"):** Take a good direct min-fuel bang-bang solution (a
PSR-refined direct-collocation result), hand it to an **actual indirect
solver** that holds the switch *structure* fixed and places every switch at
**S(τ)=0 exactly**, producing exact sub-mesh switch times + exact costates and
a **continuous-time first-order PMP certificate** — and *reporting* structure
diagnostics (a switch that should vanish, or one that is missing) without acting
on them. This is "point 3" of the direct↔indirect roadmap and the named method
**IFS = Indirect Finishing Solve** (sibling of the built **PSR = PMP-Steered
Refinement**; see the `pmp-mesh-refine-prototype` memory and
`sundman_minfuel/refine/README.md`).

Explicit non-goal for v1: the outer add/remove-switch structure search and the
transition band 1.01–1.11× (that is Option B / v2, built on a validated IFS
core). See §7.

## 1. Why this is now tractable (and was not before)

The campaign's prior indirect attempts all failed, for two separable reasons
(`LOW_THRUST_MINFUEL_CAMPAIGN.md`, `msband-indirect-campaign` memory):
1. **Perigee sensitivity** — single shooting over the ~40-rev spiral has ~10⁶
   terminal sensitivity; the basin is a needle.
2. **Bang-bang non-smoothness** — every attempt *smoothed* the throttle
   (Bertrand–Épénoy entropy/tanh, parameter ε) and continued ε→0; the 1/ε
   switch layers made LM crawl (cond(J) ~1e8–1e9), even after the Sundman
   rebuild fixed the rank deficiency.

IFS removes **both** by construction:
- Wall 1 → **Sundman-τ domain + multiple shooting** (short arcs). τ-domain was
  proven to tame the perigee conditioning in ms_band; short arcs tame the
  sensitivity single shooting could not.
- Wall 2 → **the hard throttle**. IFS never smooths. On each arc the throttle is
  a *known constant* (u=1 burn / u=0 coast, from the fixed structure), so there
  is **no ε and no 1/ε layer anywhere**. The EOM is `sms_eom` with the entropy
  term deleted.

And it now has the one thing the old attempts never had: an **excellent seed for
every unknown** from the PSR-refined direct solution (§4).

## 2. Formulation — switch-structured multiple shooting, no saltation matrices

**Domain.** Sundman variable τ with dt/dτ = κ(r) = r₁^p (p=1.5), physical time t
carried as a state, τ_f a fixed constant; the fixed transfer time is the
terminal condition t(τ_f)=t_f. 16-dim augmented state per node:
`Y = [r(1:3); v(4:6); m(7); t(8); λ_r(9:11); λ_v(12:14); λ_m(15); λ_t(16)]`
(the `sms_eom` layout).

**Switches are explicit shooting nodes — no saltation.** The user's "saltation-
aware" requirement is satisfied by the *cleaner equivalent*: a shooting node sits
at each switch, an arc runs *between* adjacent switches with its constant
throttle, and each switch time is pinned by S(τ_i)=0. Because a switch is a node
(not an event detected *inside* an arc), the Jacobian needs no saltation jump — a
switch time's sensitivity is the ordinary endpoint sensitivity of its two
neighboring arcs. This is the "switch-structure parameterization" GPT-5.6 named
as the cure (`msband-indirect-campaign`), and it avoids the 1/Ṡ blow-up that
saltation matrices carry at near-grazing switches.

**The switching function** (computed for the S=0 residual and diagnostics, NOT
used to set u): `S = 1 − ‖λ_v‖·c/m − λ_m`. Convention: burn arc S<0, coast arc
S>0. The optimal thrust direction on a burn arc is the primer α = −λ_v/‖λ_v‖.

## 3. Unknowns, residual, squareness

For **k** switches (→ **k+1** arcs, indexed j=0..k; arc j spans [τ_j, τ_{j+1}]
with τ_0=0 and τ_{k+1}=τ_f fixed; τ_1…τ_k the switch times):

**Unknown vector Z** (size **8 + 17k**):
- λ_0 (8): initial costate at τ_0. The initial *state* (r,v,m,t)=(rv0, m0, 0) is
  fixed and baked in — not an unknown.
- N_1…N_k (16 each = 16k): the full augmented state at each switch node.
- τ_1…τ_k (k): the switch times.

**Residual R(Z)** (size **8 + 17k**, square):
- **Continuity** (16k): integrate arc j (j=0..k−1) from node N_j over
  [τ_j, τ_{j+1}] with throttle u_j; require the endpoint = N_{j+1}. (16 per
  interior node.)
- **Terminal** (8): integrate arc k from N_k over [τ_k, τ_f]; require rendezvous
  r,v (6) + transversality λ_m=0 (1) + fixed time t=t_f (1).
- **Switch conditions** (k): S(N_i)=0 at each switch node i=1..k.

Square check: 16k + 8 + k = 8 + 17k = dim(Z). ✓
(k=1 → 25; k=3 → 59; k=25 → 433.)

**Terminal-BC modes** (both 8 conditions, both keep the system square):
- `"rendezvous"` (the real problem, rungs 2–3): r,v rendezvous (6) + λ_m=0 (1) +
  t=t_f (1).
- `"fixedState"` (an interior window, rung 1): the window-end state
  (r,v,m,t) is FIXED (8) — no transversality, since the window-end mass is
  fixed, not free. This is what makes the ground-truth window test possible.
`ifs_residual` selects the mode by an option; everything else is identical.

**Why keep the node states as unknowns** (instead of reducing to λ_0 + switch
times, size 8+k): eliminating the nodes is *stitched single shooting* and
reintroduces the ~10⁶ perigee sensitivity that provably fails. The node
unknowns ARE multiple shooting — the sensitivity cure. The cost is a larger but
**block-sparse** system: arc j's continuity depends only on N_j, τ_j, τ_{j+1};
S(N_i)=0 only on N_i. Exploit that sparsity in the Jacobian.

**Jacobian.** Complex-step (reuse the validated pattern from
`shoot_residual_minfuel.m` / `ms_jacobian_cs.m`), evaluated per residual block
against only the unknowns it depends on (the block-sparse pattern), so the 433-
unknown case costs a handful of short arc integrations per column-block, not 433
full-spiral integrations. Integrator tolerances RelTol 1e-13 / AbsTol 1e-15
(the ms_band residual-floor lesson: floor = costate magnitude × RelTol).

**Solver.** Levenberg–Marquardt via `lsqnonlin` (reuse `ms_solve.m` /
`solve_minfuel_indirect.m` settings: `Algorithm='levenberg-marquardt'`,
`SpecifyObjectiveGradient=true`, `ScaleProblem='jacobian'`). No CasADi/IPOPT —
IFS is pure MATLAB ODE + LM. (The *seed* comes from a CasADi-made direct `.mat`,
but IFS itself never calls CasADi.)

## 4. Seeding — the reason it converges

The PSR-refined direct solution seeds **every** unknown:
- **switch times τ_i** ← `diag.tauCr` from `pmp_refine_indicator` (PSR already
  computes the sub-cell S=0 roots — direct synergy).
- **node states N_i** ← the direct trajectory (r,v,m,t) sampled at those τ.
- **costates (λ_0 and the λ-block of each N_i)** ← the mode-'d' dual→costate map
  (`sms_seed_duals`), β-scaled. The switching function's inhomogeneous "1" fixes
  the *absolute* costate scale (the problem is NOT costate-homogeneous), so the
  dual-fit β gives the right magnitude, not merely the direction.
- **arc throttles u_j** ← the direct solution's throttle on each arc (the fixed
  structure: the alternating burn/coast pattern, starting from whatever the
  direct solution starts with).

A seed-residual sanity gate (‖R(Z_seed)‖ reported) precedes every solve.

## 5. Success criterion + diagnostics (reported, not acted on)

**Success:** ‖R(Z)‖ converges below ~1e-8 (target 1e-10, integrator-floor
permitting). The converged solution is then a genuine **continuous-time
first-order PMP extremal**:
- S=0 exactly at each switch (residual, by construction),
- sign law verified per arc post-hoc (S<0 on burns, S>0 on coasts, bounded away
  from 0 on arc interiors — no singular arcs),
- rendezvous + transversality met,
- fixed transfer time met.

This **upgrades** the direct/PSR O(h²) mesh-limited certificate to exact switch
times + exact costates. The deliverable per rung is: the exact switch times vs
the PSR-quantized ones (the sub-mesh correction), the (expected tiny) propellant
delta, and a certificate paragraph gated on all checks passing.

**Structure diagnostics (reported):**
- *Vanishing arc*: two switch times converge (arc length τ_{i+1}−τ_i → 0) ⇒ that
  switch pair is spurious ⇒ "structure suggests FEWER switches here."
- *Sign violation*: S>0 on an assumed burn arc (or S<0 on a coast) that the solve
  cannot drive out ⇒ "a switch is MISSING here."
These are v1 **reports** only. Acting on them (the add/remove-switch outer loop)
is Option B / v2.

## 6. Staging — the validation ladder (milestones, one shared codebase)

Revised 2026-07-11 to use ready seed artifacts (the arrival-leg and 3-switch
`.mat`s the first draft named are not in a usable layout — the leg lives only in
`attic/`, and `minfuel_from_energy_seed.mat` holds a raw decision vector with no
dual costates). The two readily-seedable solutions are `legacy_ms_f1120.mat`
(1.12×, ~10 switches, already carries `out.lamDef`) and
`sundman_minfuel_certified.mat` (1.15×, 25 switches, duals regenerated via the
existing `sundman_minfuel/refine/prep_refine_seed.m`).

**REVISED AGAIN 2026-07-11 (post-diagnosis — window dropped as the gate).**
The interior-window rung (below, struck through) turned out **structurally
rank-deficient** and cannot be the gate. Root cause (verified by SVD +
GN-consistency diagnostics, `.superpowers/sdd/ifs-diag-*`): mass dynamics are
costate-independent (`dm/dτ = -κ·u·Tmax/c`) and on a coast arc the λ_m equation
vanishes identically, so λ_m is decoupled; the window's `fixedState` terminal
fixes the end *mass* but imposes **no** λ_m(τ_f)=0 transversality, leaving λ_m an
unconstrained gauge → smallest singular value 1.5e-16, isolated 9 orders below
its neighbor. The **full problem** (`rendezvous` terminal) is **full-rank** (all
178 SVs healthy, transversality row pins the gauge, GN-consistent to 5e-11) — so
the window failed for a reason the real deliverable does not have. **The gate is
repointed to the full 1.12× solve.** The window code is kept as a documented
rank-deficiency finding, not a gate.

- ~~Rung 1 — interior-window ground-truth solve~~ **(dropped — rank-deficient in
  the λ_m gauge; see above).** The IFS machinery is instead validated by the
  Task 1–4 unit tests (EOM matches `sms_eom`; residual zero on a continuous
  ground-truth; CS-vs-FD Jacobian for fixedState + rendezvous + k=2; seed builds).
- **Rung 2 — full 1.12× (~10 switches) — NOW THE MAKE-OR-BREAK GATE**. First
  real IFS solve: terminal-BC mode "rendezvous" (r,v rendezvous + λ_m=0
  transversality + t=t_f), moderate switch count, full-spiral perigee sensitivity
  + Sundman + multiple shooting. Full-rank at the seed (cond≈3.5e11 from shooting
  sensitivity, but resInRange 5e-11 — reachable). Seed =
  `results/minfuel/legacy_ms_f1120.mat` directly (carries `out.lamDef` — zero
  prep). Convergence is promising, not guaranteed; running it is the real test.
- **Rung 3 — the 25-switch 1.15× headline**. The full certified point,
  seeded via `prep_refine_seed('../sundman_minfuel/sundman_minfuel_certified.mat', …)`.
  k=25, 433 unknowns. The deliverable: exact switch times vs the direct/PSR
  quantized ones + a continuous-time certificate + structure diagnostics.

Each rung is a milestone gate: do not advance until the prior rung certifies.
Rung 1 also gives ground truth for the residual/Jacobian correctness that
rungs 2–3 rely on.

## 7. Non-goals (YAGNI)

- No outer structure search (add/remove switches) — diagnostics are *reported*
  only; acting on them is v2.
- No attack on the 1.01–1.11× transition band (Option B).
- No second-order sufficiency test (the switch-time Hessian) — optional future.
- No singular-arc machinery (this problem is bang-bang; S is bounded away from 0
  on arc interiors — a checked property, not an assumed one).
- No new direct solver and no CasADi dependency in IFS itself.
- Not wired into `minfuel_at_tf` / PSR — IFS is a standalone finishing stage that
  consumes a direct/PSR `.mat`.

## 8. Module structure

New folder `NLP_lowThrust_GTO_tulip/ifs/`:

| file | role |
|---|---|
| `ifs_problem.m` | constants/endpoints/odeOpts/Sundman p factory (à la `sms_problem`), reusing `cr3bp_lt_params` + `gto_tulip_endpoints` |
| `ifs_eom.m` | hard-throttle 16-dim Sundman PMP EOM; throttle u passed as an arc parameter (0/1), entropy term dropped; also returns S for diagnostics |
| `ifs_pack.m` / `ifs_unpack.m` | Z ⇄ (λ_0, {N_i}, {τ_i}) |
| `ifs_residual.m` | R(Z) multiple-shooting residual (+ block-sparse complex-step Jacobian / sparsity pattern) |
| `ifs_seed.m` | build Z from a direct/PSR `.mat` (full problem or an extracted window): switch times (`diag.tauCr`), node states (interp), costates (mode-'d' β-scaled duals via `sms_seed_duals`), arc throttle pattern; sets the terminal-BC mode |
| `ifs_solve.m` | LM driver (`lsqnonlin`), seed-residual gate, returns converged Z + per-iterate log |
| `ifs_certify.m` | post-hoc S=0/sign-law/rendezvous-or-fixedState/transversality checks + structure diagnostics + certificate paragraph + report/figure |
| `run_ifs_window.m` / `run_ifs_1p12.m` / `run_ifs_1p15.m` | rung 1 (interior-window ground truth) / rung 2 (full 1.12×) / rung 3 (1.15× headline) entry points + RESULTS |
| `setup_paths.m`, `test_*.m` | path setup + per-module tests |

Reuses (unchanged): `cr3bp_lt_params`, `gto_tulip_endpoints`, the `sms_seed_duals`
dual→costate map, the `ms_solve` LM settings, and the CS-Jacobian pattern from
`shoot_residual_minfuel` / `ms_jacobian_cs`.

## 9. Environment (this machine)

- MATLAB **R2025b** (`/Applications/MATLAB_R2025b.app/bin/matlab -batch`);
  R2025a license-broken. Use the `matlab-headless` skill; write a `.m` and run
  `-batch "cd('<dir>'); script"` (multi-line `-batch` strings fail here). Filter
  the Home-License banner.
- **No CasADi needed for IFS** (pure ODE + `lsqnonlin`). CasADi is only present
  because the seed `.mat` files were produced by the direct solver.
- Complex-step Jacobians require the EOM to be CS-safe (no `abs`/`max`/`norm` on
  complex paths) — the `sms_eom` softplus/branch discipline carries over; the
  hard EOM is simpler (no entropy branch) but still keep `sqrt(sum(·.^2))` for
  norms, guard `‖λ_v‖` only where α is used (burn arcs).
- Integrator RelTol 1e-13 / AbsTol 1e-15 (ms_band floor lesson).

## 10. References

- `LOW_THRUST_MINFUEL_CAMPAIGN.md` — two-walls analysis; indirect failure record.
- `ms_band/MS_BAND_CAMPAIGN.md`, `msband-indirect-campaign` memory — the
  smoothed-MS crawl, the GPT-5.6 "switch-structure parameterization" verdict,
  the mode-'d' dual map.
- `sundman_minfuel/refine/` (PSR) — the seed source; `diag.tauCr` sub-cell roots.
- `lowThrust_GTO_tulip/{lt_pmp_eom_minfuel,shoot_residual_minfuel,solve_minfuel_indirect}.m`
  — reusable EOM/residual/LM patterns (the smoothed single-shooting that failed).
- Zhang, Topputo, Bernelli-Zazzera, Zhao, JGCD 38(8), 2015 — indirect min-fuel
  CR3BP (switch-structure reference).
