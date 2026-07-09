# Tier-1 PMP certification — scope

**Goal:** certify that the converged Sundman min-fuel solution
(`sundman_minfuel_certified.mat`, 25-switch bang-bang) is a genuine Pontryagin
extremal, by recovering a costate history `λ(t)` that satisfies the min-fuel
adjoint equations along the frozen NLP state trajectory and whose switching
function `S = 1 − ‖λ_v‖ c/m − λ_m` sign-changes exactly at the 25 throttle
switches. This is a **verification** pass (no BVP re-solve, no convergence
risk): it either confirms the first-order PMP conditions or exposes a defect.

Scoped 2026-07-08. This is Tier 1 of the two-tier plan; Tier 2 (independent
Sundman-regularized indirect re-solve) is a separate, larger milestone. See
[[minfuel-gto-tulip-solved]] and `../LOW_THRUST_MINFUEL_CAMPAIGN.md`.

---

## Approach: recover the costates FROM the trajectory (not from the NLP duals)

The position/velocity costate block is a **linear homogeneous** system driven
only by the state trajectory:

```
λ̇_r = −Gᵀ λ_v,      λ̇_v = −λ_r − H_cᵀ λ_v
```

with `G`, `H_c` exactly as in `../../lowThrust_GTO_tulip/lt_pmp_eom_minfuel.m`
(lines 58–61). It does not depend on `λ_m` or on the control magnitude, so the
whole costate history is fixed by 7 initial numbers `λ(0)=[λ_r0; λ_v0; λ_m0]`.
We recover `λ(0)` by a tiny least-squares fit, then check PMP consistency —
the same "propagate the adjoint alongside the frozen state trajectory and
compare" pattern already used in `../../orbit_transfer/primer_check.m`.

**Why this over CasADi KKT-dual extraction:** it avoids the two fragile parts of
the dual route — CasADi's dual sign/layout convention, and the
Sundman-κ / Δσ / τ_f de-scaling of the discrete multipliers. Dual extraction is
kept as an OPTIONAL independent cross-check (see below), not the primary path.

## Steps

1. **Load the certified solution.** Everything needed is already in
   `sundman_minfuel_certified.mat`: `out.X` = `[r;v;m;t]` per node (rows 1:6,7,8),
   `out.U` = `[α;s]` (rows 1:3, row 4), plus the converged `eps` and `sigma`,
   `tauf`, `pSund`. No duals, no re-solve.

2. **Adjoint propagator along the frozen trajectory.** Integrate the 7-costate
   system (linear `λ_r,λ_v` block above + `λ̇_m = −‖λ_v‖ u T_max/m²`) on the
   SAME nodes as the NLP, reading `r(t)`, `m(t)`, `u(t)=s(t)` from `out`. Do it
   in τ with the same κ = r₁^pSund scaling and the same trapezoidal scheme the
   NLP used — zero interpolation error (the no-resample discipline). The
   `λ_r,λ_v` block being linear makes this one well-conditioned matrix
   propagation.

3. **Fit λ(0) by least squares (7 unknowns), `lsqnonlin`.** Residual vector:
   - primer-direction mismatch `α_NLP(t) + λ_v(t)/‖λ_v(t)‖` at every burn node
     (scale-invariant, so independent of any λ normalization),
   - throttle/switching mismatch `u_NLP(t) − (1−tanh(S/(2ε)))/2` at every node,
   - transversality `λ_m(τ_f) − 0` (free final mass).

   NOT circular: 7 numbers must reproduce the entire direction history
   (hundreds of nodes × 2 DOF) plus the 25-switch structure — massively
   overdetermined. If the NLP were not a true extremal, no `λ(0)` would fit and
   the residual would stay large. The smallness of the converged residual IS the
   certificate.

4. **Switching-function sign law.** With fitted `λ(t)`, form
   `S(t) = 1 − ‖λ_v‖ c/m − λ_m` at every node. Verify `S<0` where `s≈1`,
   `S>0` where `s≈0`, and that S's 25 zero-crossings coincide (within a node or
   two) with the throttle switches. Report the sign-match fraction and the
   crossing-vs-switch location deltas.

5. **Emit the certificate.** A struct (fit residual, `λ(0)`, switch-alignment
   stats, primer max direction error, transversality residual) + the money
   figure: `S(t)` on a twin axis over `s(t)`, zero line, switch markers.

## What makes it rigorous (three independent internal checks)

- **Overdetermination** — 7 unknowns vs. a full direction history; a good fit
  cannot happen by accident.
- **Adjoint-ODE consistency** — recovered `λ(t)` propagates under the
  state-derived gravity gradient, in the spirit of `primer_check.m`'s
  err ≈ integration tolerance.
- **Transversality** — `λ_m(τ_f) ≈ 0` emerges as a fit residual, not an
  assumption.

## Reuse map

| need | existing asset |
|---|---|
| adjoint EOM (`G`, `H_c`, `λ̇` block, `S` definition) | `../../lowThrust_GTO_tulip/lt_pmp_eom_minfuel.m` (lift RHS, drop the state block) |
| propagate-costates-along-frozen-trajectory pattern | `../../orbit_transfer/primer_check.m` |
| `lsqnonlin` continuation harness | `../../lowThrust_GTO_tulip/solve_minfuel_indirect.m` |
| CR3BP/LT constants, endpoints, solution | `cr3bp_lt_params.m`, `sundman_minfuel_certified.mat` |

Net new code: one file `certify_minfuel_pmp.m` (~120–150 lines) + a one-line
call after the solve in `run_certified_minfuel.m`.

## Risks / unknowns

- **Main unknown:** whether recovering `λ(0)` from the SMOOTHED control at the
  converged (very small) ε is sharp enough — `atanh` is sensitive near switches.
  Mitigation: fit against burn/coast direction + sign structure (robust) rather
  than exact `u` inside transition nodes; optionally fit at a slightly larger ε
  from the homotopy history (smoother throttle), then confirm the S-structure
  persists as ε→0.
- **Minor:** trapezoidal costate estimates are half-weight-off at the first/last
  node — irrelevant to the 25 interior switches, worth a comment.

## Optional cross-check (belt-and-suspenders, for publication)

Independently read CasADi's KKT duals on the defect constraints
(`opti.dual(...)`; `casadi_minfuel_sundman.m` currently returns none — ~10-line
add), de-scale by the trapezoid weight, τ_f, and κ to get a second estimate of
`λ(t)`, and confirm it agrees with the trajectory-recovered one. Stronger claim,
but carries the de-scaling fiddliness — only if wanted.

## Effort

~1 day: ~2 h adjoint propagator, ~½ day fit + getting the checks to pass,
~2 h certificate struct + figure.

## Deliverables

- `certify_minfuel_pmp.m` — new library function (loads a solution struct +
  params, returns the certificate struct, makes the S-vs-s figure).
- One-line hook in `run_certified_minfuel.m` to certify right after solving.
- The `S(t)` vs `s(t)` figure — the artifact for the paper / campaign note.

---

## FINDING (2026-07-08) — primal costate recovery is the WRONG route; use duals

Implemented and tested `certify_minfuel_pmp.m` three ways to recover the costates
FROM THE PRIMAL (frozen state trajectory + primer directions). **All three fail**,
for a structural reason the scope underestimated:

1. **Forward-shoot lambda(0) + lsqnonlin fit** — the homogeneous costate map
   amplifies by **~5e11** over the ~40-rev spiral (measured), even WITH the
   Sundman kappa scaling. Forward propagation of lambda(0) is hopelessly
   ill-conditioned. (Also eps=0 in the saved solution, so the smoothed-throttle
   fit term is gradient-free — a secondary problem.)
2. **Global sparse min-norm LS over all node costates** (recursion + primer
   directions) — collapses lam_v toward zero except near the single
   normalization node. Root cause: rendezvous fixes the state at BOTH ends, so
   there is **no boundary condition on lam_r, lam_v** — the position/velocity
   costates are pinned ONLY by the primer-direction data, and min-norm LS over a
   near-homogeneous system does not recover the physical (large-dynamic-range)
   solution. primerDirErr pinned at exactly 2.0 (wrong-sign lam_v).
3. **Row-normalized version of (2)** — no better (recursion resid worsens).

Diagnostics that ARE clean: the gravity-gradient G matches finite differences
(1e-5 near perigee, = FD truncation), so the adjoint model is correct; the
failure is conditioning/uniqueness, not a formula bug.

**Conclusion:** recovering costates from the primal here is the classic
covector-mapping problem, and it is genuinely hard because (a) the forward map
spans ~1e11 and (b) the endpoints give no lam_r/lam_v BC. The RELIABLE route is
the one this scope deprioritized: **read the KKT duals from IPOPT**, which are
computed by a stable sparse symmetric-indefinite factorization of the full
primal-dual system (the whole reason collocation beats shooting), so they are
immune to the forward dynamic range.

### Corrected plan (dual extraction)
1. `casadi_minfuel_sundman.m`: keep a handle to the defect constraint and return
   `out.nu = sol.value(opti.dual(defectCon))` (~10 lines). CasADi Opti supports
   `opti.dual(con)`.
2. Regenerate a certified `.mat` WITH duals — cheapest is a single eps=0 re-solve
   warm-started from the existing `out.X/out.U` (already optimal, converges in a
   few iters), not the full homotopy. Needs CasADi on the path.
3. De-scale nu_k -> node costate lambda_k: undo the trapezoid weight, tauf, and
   the Sundman kappa factor. Validate by the SAME two checks already coded here
   (primer alignment alpha = -lam_v/||lam_v|| at burn nodes; S sign law at the 25
   switches). The report/figure machinery in `certify_minfuel_pmp.m` is reusable
   as-is; only the recovery front-end changes.

`certify_minfuel_pmp.m` is left in the tree as the scaffold (report + figure +
switch-alignment metrics all work); its recovery front-end currently returns a
REVIEW verdict and must be swapped for the dual-based recovery above.
This pivot changes the solver and regenerates the `.mat`, so it is a larger step
than the original scope — flagged for a go/no-go decision.

---

## UPDATE (2026-07-08, after three-way diagnosis review + hard-constrained impl)

Two external adversarial reviews (GPT-5.5 + Gemini 3.1 Pro, in `reviews/`)
**overturned the "ill-posed / doomed" claim above**: they showed the failure was
a FORMULATION error — the earlier `A\b` dumped the hard adjoint recursion and the
soft primer data into one penalized min-norm LS, which collapses `lam_v`. Fix:
enforce the recursion as a HARD constraint.

Implemented and tested that fix in `certify_minfuel_pmp.m`:
1. hard-constrained saddle point (recursion = equality) — **adjoint residual
   dropped from ~1e-3 to ~3e-14**, so the ~1e11 conditioning wall is defeated IN
   THE SOLVE. The reviewers were right on this point.
2. added the signed constraint `lam_v = -rho*alpha, rho>=0` via `lsqlin` (the
   projector `(I-alpha alpha')lam_v=0` is sign-BLIND) — cut wrong-sign burn nodes
   428 -> 214.

**But it still does not certify (REVIEW, not PASS), for a DEEPER reason:**
~24% of the NLP burn-node thrust directions cannot be matched by ANY single
costate satisfying the recursion — even the scale-anchor node cannot reach
`||lam_v||=1`. Diagnostics: the un-fittable / wrong-sign arcs sit at PERIGEE
(`r1~0.037`, the FINE part of the mesh, `dt/node~3e-4`), NOT the coarse apogee
(my coarse-mesh hypothesis was wrong). Root cause: the hand-built `B_k` is only a
2nd-order TRAPEZOIDAL discretization of the CONTINUOUS adjoint; over ~40 revs that
approximation error, amplified by the perigee sensitivity, makes distant burn
arcs mutually inconsistent. No continuous costate threads all 40 perigee passes.

**Settled conclusion:** certifying a DISCRETE NLP optimum requires the EXACT
DISCRETE adjoint = the NLP's own KKT duals (GPT-5.5's route), not a reconstructed
continuous costate. This is the same dual pivot as above, but now motivated by a
principled discrete-vs-continuous-optimality argument rather than only by
conditioning. The continuous recovery in `certify_minfuel_pmp.m` remains a useful
consistency diagnostic (it matches ~76% of directions), not a full certificate.
Go/no-go on the dual pivot still open.

---

## UPDATE (2026-07-08) — dual pivot IMPLEMENTED; primer + transversality certified

The go/no-go above is resolved **GO**: the KKT-dual route is implemented in the
solver, and it succeeds where the continuous recovery gave only ~76%.

**What was added** (`casadi_minfuel_sundman.m`, committed 3bfd0b2):
- After `opti.solve()`, read the full multiplier vector `lamAll = sol.value(opti.lam_g)`.
  The dynamics-defect block (the first `8*N` entries, added first, so no layout
  guessing) reshapes to `out.lamDef` `[8 x N]` = `[λ_r; λ_v; λ_m; λ_t]` per
  interval — the **discrete costates**, up to a positive mesh-weight scaling and a
  global sign. `out.lamAll` keeps the whole stacked g-multiplier (incl. the
  throttle-bound duals, which encode the switching-function sign law).
- Two scale-invariant PMP checks are computed in-solver:
  `out.primerAlignDeg` (mean angle between the NLP thrust direction and the
  costate primer `-λ_v/‖λ_v‖` on burn arcs) and `out.lamMassEnd` (terminal
  mass-costate proxy).

**Results on `sundman_minfuel_certified.mat`** (eps=0 re-solve, defect 2e-14):
- **primer alignment = 0.058°** — the direct-NLP thrust direction matches the
  costate primer to a hundredth of a degree on every burn arc. (Continuous
  recovery could not thread this; the duals do it immediately, confirming the
  discrete-adjoint argument.)
- **λ_m(τ_f) = −1.7×10⁻⁷** — transversality (free final mass) satisfied.
- costate magnitudes sane: `‖λ_r‖∈[0.6,3.1e2]`, `‖λ_v‖∈[0.01,3.8]`,
  `λ_m∈[−33, ~0]` (monotone to 0 at τ_f).

Both checks are **scale-invariant** (direction: any positive weight cancels;
transversality: 0 stays 0), so they need NO de-scaling and are clean now.

**What remains for a full Tier-1 certificate:**
1. **Switching-function sign law** (the scale-DEPENDENT piece). Form
   `S = 1 − ‖λ_v‖ c/m − λ_m` from the duals and verify `S<0` on burns, `S>0` on
   coasts, with S's zero-crossings at the 25 switches. This needs the single
   de-scaling step the scope flagged: undo the trapezoid weight, `τ_f`, and the
   Sundman `κ` factor to turn `lamDef` (defect multipliers) into node costates
   `λ_k` on the physical-time measure. The `certify_minfuel_pmp.m` report/figure
   machinery (S-vs-s twin axis, switch-alignment stats) is reusable as-is; only
   the recovery front-end swaps to `out.lamDef` + de-scale.
   - Cheaper scale-free alternative already available: the **throttle-bound
     multipliers in `out.lamAll`** directly carry the KKT sign condition on `s`
     (their signs at active bounds ARE the discrete switching law); decoding them
     avoids the `c/m` vs `‖λ_v‖` scale-mixing entirely.
2. **Independent adjoint-ODE consistency** (optional, belt-and-suspenders): confirm
   `lamDef` propagates under the state-derived gravity gradient. The continuous
   recovery already checked G vs FD (clean); this would close the loop.

Status: **primer-direction and transversality conditions certified via the
duals; switching-function sign law and independent adjoint check remain.** The
primal (defect 2e-14, endpoints exact) plus these two costate conditions is
strong PMP evidence but not yet the complete first-order certificate.
