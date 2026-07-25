# Optimality certification — cross-campaign register

**Purpose.** One place that answers "what do we actually certify, on which
transfer, and what is still missing" — for **both** first- and second-order
optimality (goal 1 of the standing program).

**Why one file and not two.** The two orders are not independent deliverables.
The Maurer–Osmolovskii sufficient condition for a bang-bang local minimum is
assembled from three ingredients that straddle the boundary:

    (i)   projected switching-time Hessian positive definite      [2nd order]
    (ii)  regular switching, Sdot != 0 at every switch            [1st order]
    (iii) primer / strengthened Legendre-Clebsch direction optimality  [1st order]

Splitting the record would scatter one certificate's ingredients across two
documents — the exact fragmentation that let one campaign re-plan work another
had already built and killed. Part A is a **coverage matrix** (first order is
largely solved; what matters is *where it runs*). Part B is a **research
register** (second order is open; what matters is *what broke and why*).

**Scope discipline — read this before adding to it.**
- Function headers stay authoritative for *how an instrument works* and for
  its own derivation. Do not restate derivations here; link them.
- This file owns what no header can: the **status matrix**, the **blocking
  mechanisms** (which are general, not per-file), what is **ruled out**, open
  **leads**, and the **decision record**.
- It is a register of experiments. It is only worth keeping if it is updated
  when one runs. Append to §6 every time.

Why it exists at all: in July 2026 the CR3BP campaign's TODO independently
planned a "tier 1 / tier 2 / tier 3" second-order program that the tulip
campaign had already built and found blocked. That is the same failure mode
diagnosed in `earth_elliptic_to_geo/process/LESSONS_DUAL_EXTRACTION.md` (L1:
the sibling implementation was invisible).

---

# PART A — First order (PMP extremality)

**State: largely solved.** The instruments exist and pass; the gaps are
*coverage* and one *structural* caveat, not correctness. (They were not
trustworthy before 2026-07-25 — see
`earth_elliptic_to_geo/process/LESSONS_DUAL_EXTRACTION.md` for why, which is
not restated here.)

## A1. Instruments

| # | condition tested | implementation |
|---|---|---|
| 1 | Hamiltonian minimization over thrust direction (β ∥ −primer) | `mee_primer_switch` + `verify_pmp_mee`; tulip `certify_minfuel_pmp` |
| 2 | switching-function sign law (S<0 ⇔ thr=1) | same |
| 3 | switch-time alignment (S zero-crossings vs actual throttle switches) | same |
| 4 | transversality λ_m(t_f)=0 (free final mass) | same, relative gate |
| 5 | Hamiltonian constancy (λ_t const ⇔ H_t conserved; autonomous only) | `hamiltonian_const_check`, `hamiltonian_along_traj` |
| 6 | full KKT certificate — stationarity, dual feasibility, complementarity | `sosc_recover_kkt` + `sosc_kkt_residual`; `results/dual_anomaly/diag_t1_beta.m` |
| 7 | no singular arc (no run of \|S\|≈0) | `verify_pmp_mee` H2 check |
| 8 | **regular switching, Sdot != 0** — *Maurer ingredient (ii)* | `psr_second_order.m:278` (see gap G3) |
| — | primal, dual-free cross-checks (switch structure, edge fraction, bound saturation) | `switch_structure`, solver `boundSat` |

**Note — two different costate sources.** The earth/CR3BP path reads the NLP's
own multipliers (`opti.lam_g`). The tulip path (`certify_minfuel_pmp`) instead
*reconstructs* costates by sparse least-squares on the discrete adjoint
recursion plus primer directions, with the scale pinned by S=0 at switches.
The second is methodologically stronger (it does not consume the solver's
duals) and is the partial answer to gap G1.

## A2. The PMP necessary conditions — which do we actually check?

The textbook list, with our status. Conditions marked ★ are ones easy to omit
from an informal list but which matter here.

| condition | status |
|---|---|
| **(a) state dynamics** satisfied along the trajectory | YES — machine-tight, but **at the collocation nodes, to collocation order**. Not the same as satisfying the continuous ODE; see the caveat below |
| **(b) costate dynamics** (adjoint equation) | **Implicitly** — KKT stationarity w.r.t. the state variables *is* the discrete adjoint recursion, verified at ‖∇ₓL‖∞ = 1.5e-14 (earth). Tulip enforces it explicitly by least-squares. **Nobody checks a continuous-time costate ODE residual** |
| **(c) transversality** | λ_m(t_f)=0 checked (earth gate, tulip reported). Terminal elements are fully pinned, so no condition there. In the L-domain, free ΔL supplies a transversality-type stationarity row — enforced by KKT, never *reported* as a PMP quantity |
| **(d) Hamiltonian constancy / H≡0 for free t_f** | see the correction below — partially, and the free-t_f case is unchecked (G4) |
| ★ **(e) the minimum condition itself**: u\* = argmin over the admissible set, not merely ∂H/∂u = 0 | **YES, and it is our strongest check.** For control-affine dynamics with ‖β‖=1 and thr∈[0,1] the pointwise minimizer is exactly β = −p/‖p‖ with thr from sign(S) — so the primer-alignment + sign-law gates *are* the minimum condition. Note ∂H/∂u = 0 would be the wrong test: the throttle optimum is at a bound, not a stationary point |
| ★ **(f) non-triviality / normality**: (λ₀,λ) ≢ 0, λ₀ ≥ 0; no abnormal extremal | **NOT CHECKED anywhere.** In the direct/NLP form λ₀ ≡ 1 by construction, so abnormality is not even representable — it becomes a real question only for the indirect solver (goal 6) |
| ★ **(g) no singular arc** (else generalized Legendre–Clebsch / Kelley) | YES — H2 run-of-\|S\|≈0 check (earth) |
| ★ **(h) regular switching, Ṡ ≠ 0** | computed but not a standing gate (G3) |
| ★ **(i) Weierstrass–Erdmann corner conditions** (λ, H continuous across a switch) | satisfied *by construction* — control-affine, no active state constraints; worth stating, not worth testing |
| ★ **(j) weak Legendre–Clebsch** ∂²H/∂u² ≥ 0 | vacuous at ε=0 (linear in throttle); meaningful only on the ε>0 homotopy legs |
| ★ **(k) state-constraint jump conditions** | not applicable — the state boxes are inactive; the solver's `boundSat` diagnostic is what confirms that, and it is checked |

**Correction worth being precise about, on (d).** Three distinct cases, and we
conflate them at our peril:
- *autonomous + **fixed** t_f* → H = **const, generally ≠ 0**. This is our
  min-fuel setup (t_f = c_tf·t_f,min, fixed, per HMG-2004 verbatim). Measured:
  H_t = −λ_t = −0.0277, constant to CoV ~1e-9. Correct and checked.
- *free t_f* → H(t_f) = 0, and if autonomous then **H ≡ 0**. This is our
  min-time anchors. **Never checked** (G4) — and it is the sharpest available
  test of those solutions, because it is a *value* test, not a constancy test.
- *non-autonomous* (earth CR3BP: the Moon's position is an explicit function
  of t) → H is **not** constant, and we currently just label the check
  "informational" and stop. That discards a real test: the correct
  generalization is **dH/dt = ∂H/∂t**, which is computable here. Upgrading
  `hamiltonian_const_check` to test that residual would restore a genuine
  first-order check to the whole CR3BP campaign. (Recorded as G6.)

**Caveat on (a) that we should not soften.** We satisfy the *discrete* defects
to ~1e-14, not the continuous ODE. `psr_switch_hessian`'s blocked finding is
the hard evidence: re-integrating the exact direct control from x₀ over ~40
revs diverges by ‖r‖~3, ‖v‖~5. So "dynamics satisfied at every time point" is
true of the transcription, and only true of the continuous problem up to
collocation order on the given mesh. Any claim aimed at a paper should say so.

## A3. Coverage matrix

| campaign / target | primer + sign law | switch alignment | transversality | H constancy | full KKT | Sdot ≠ 0 |
|---|---|---|---|---|---|---|
| earth 2-body (10 → 0.1 N, + 2 PSR) | **9/9 @ 0.000° / 100%** | yes | yes | yes (CoV ~1e-9) | yes (1.5e-14) | no |
| earth CR3BP (lunar) | 10 N, 5 N pass; **rungs 2.5→0.1 N not run** | yes | yes | reported, informational (non-autonomous — correct) | not run | no |
| GTO→tulip | yes, ~0.06° (independent reconstruction) | yes | reported | no | no | built, not a standing gate |
| GTO→ELFO | **none** | none | none | no | no | no |
| min-time anchors (free t_f) | direction only | n/a (all-burn) | n/a | no | no | n/a |

## A4. Gaps — what we lack

- **G1 — independence.** Earth and CR3BP verify PMP from the NLP's *own*
  duals: self-consistency, not an independent witness. It catches convergence
  to a non-extremal; it cannot catch an error shared with the transcription.
  Real fix: costates from an indirect solve (goal 6). Cheap partial fix: port
  tulip's least-squares reconstruction as a second opinion.
- **G2 — coverage.** ELFO has no first-order gate at all; CR3BP is verified
  only at 10 N and 5 N; tulip's certifier carries a known false negative
  (strict integer PMP-crossing match fails on node-grazing switches — see the
  tulip TODO C3).
- **G3 — Sdot != 0 is not a standing gate.** It is computed, but inside
  `psr_second_order.m`, whose headline verdict is NOT APPLICABLE (Part B, M2),
  so it is not reported per row. This is a *required ingredient* of the
  second-order sufficient condition — promoting it is cheap and directly
  serves Part B §5.
- **G4 — free-final-time Hamiltonian condition unchecked.** The min-time
  anchors (`tfMin_tulip`, `tfMin_ELFO`, `run_mintime_mee`) are free-t_f
  problems where H must take a specific value at t_f. Never verified.
- **G5 — no mesh bands on the gates.** Switch counts are known mesh-sensitive;
  the first-order gate values are reported as single numbers, not bands.
- **G6 — the non-autonomous Hamiltonian test is abandoned rather than
  generalized.** Under lunar gravity H is genuinely not constant, and the CR3BP
  driver correctly refuses to gate on constancy — but it then checks nothing.
  Test `dH/dt = ∂H/∂t` instead (see A2's correction on (d)).
- **G7 — normality / abnormal extremals never examined.** Invisible in the
  direct formulation (λ₀ ≡ 1 by construction); becomes a live question the
  moment an indirect solver exists.

## A5. What can and cannot be claimed today

**Can:** machine-tight *on the transcription*, PMP-extremal, with the pointwise
minimum condition verified (primer + sign law), bang-bang structure confirmed,
no singular arcs, and switch times consistent with the reconstructed switching
function.

**Cannot:** (i) that the continuous ODE is satisfied to the same accuracy — see
A2's caveat on (a); (ii) *strict* local minimality — Part B gives weak local
minimality on part of the tulip set and nothing elsewhere; (iii) anything
global.

---

# PART B — Second order (local minimality)

**State: partially delivered.** Weak local minimality is certified on 12 of 17
ε=0 tulip rows by IPOPT's native inertia; **strict** local minimality is
certified nowhere, and three other campaigns have no second-order verdict at
all.

## 1. Status matrix

| instrument | campaign | verdict | where |
|---|---|---|---|
| **IPOPT native inertia (δ_w)** | GTO_tulip | **DELIVERS** — wired into production PSR; 137/137 rows carry a verdict. **ε=0 (bang-bang): 12 certified / 5 not.** ε>0: 100 / 20 | `GTO_tulip/direct/PSR/psr_ipopt_certify.m`, called from `run_psr.m:410` + `psr_run_one.m:156`; verdicts stored as `ipoptCert` in `PSR_data/psr_data_*.mat` |
| NLP reduced-Hessian SOSC | earth 2-body | **WEAK_MIN** @10 N (270 flat directions); INCONCLUSIVE @5 / 2.5 N; ERROR @1 / 0.5 N; scale-skip above `maxNullDim=10000` | `earth_elliptic_to_geo/direct/verify/sosc/`, `process/DESIGN_sosc.md` §11–12 |
| NLP SSOSC via KKT inertia | GTO_tulip | **NOT APPLICABLE (structural)** | `GTO_tulip/direct/PSR/psr_second_order.m` (FINDING, 2026-07-12) |
| Maurer–Osmolovskii switching-time Hessian | GTO_tulip | **BLOCKED** (forward-flow conditioning) | `GTO_tulip/direct/PSR/psr_switch_hessian.m` (FINDING, 2026-07-12) |
| any second-order | GTO_ELFO | not started | — |
| any second-order | earth CR3BP | planned only (tiers 1–3) | `earth_elliptic_to_geo_CR3BP/TODO.md` |
| conjugate point (Jacobi field) | all | not built; gated on an indirect solver | `BCP2010` §2.3–2.4; CR3BP TODO Phase 2 |

**Net (corrected 2026-07-25):** second order is **not** at zero. The IPOPT
native-inertia certificate delivers on the tulip, including on genuine ε=0
bang-bang rows (12 of 17). An earlier version of this register said "zero
delivered certificates" — that was wrong, and the miss is instructive: the
instrument had **no FINDING block** (unlike the two blocked ones), so it was
invisible to a survey that read headers. *An instrument that quietly works is
easier to lose than one that loudly fails.*

**What δ_w = 0 does and does not prove.** IPOPT's inertia-controlled linear
solver adds Hessian regularization only when the reduced Hessian is indefinite;
converging with δ_w = 0 over the final iterations therefore means no negative
curvature *on the well-scaled barrier system at the final μ*. Three caveats,
all in `psr_ipopt_certify.m`'s own header and none of them fatal:
- it is a statement about the **barrier subproblem** at small μ, not exactly
  the original NLP's critical cone;
- at **ε > 0** the problem is strictly convex in throttle, so this is a strict
  local minimum of the NLP and is the easy case;
- at **ε = 0** the barrier supplies the throttle-direction curvature, so it
  certifies a **weak (non-strict)** local minimum — fully consistent with M1,
  not a contradiction of it. **Strictness for bang-bang lives in the switching
  times**, which is precisely what §5 is about.

So the honest ladder is: *extremal* (Part A) → *no negative curvature / weak
local min* (this instrument) → *strict local min* (§5, still missing).

First order, by contrast, is in good shape — see Part A, and the 9/9 sweep in
`earth_elliptic_to_geo/direct/results/verify_pmp_all.mat`.

Note that Maurer ingredients (ii) and (iii) — regular switching and primer /
Legendre direction optimality — are **first-order** checks and already exist
(A1 rows 1 and 8). What is missing is only ingredient (i), the strict
curvature test.

---

## 2. The three blocking mechanisms

These are the transferable content. Each is a property of the *problem class*,
not of an implementation.

### M1 — bang-bang extremals are weak, non-strict minima
The fuel objective is linear in throttle, so naive control-space curvature is
zero and the reduced Hessian is genuinely flat in many directions (270 of them
at earth 10 N, with a clean spectral gap: ~270 eigenvalues at 1e-10, nothing
until ~1e-4). **Strict SOSC is generically unreachable for min-fuel.** A
WEAK_MIN verdict is the honest ceiling of the NLP-level test, not a
tolerance-tuning failure.

### M2 — strict complementarity fails at the throttle bounds
The tulip's KKT-inertia test was supposed to escape M1: at a bang-bang
solution the throttle sits at a bound on every arc, so with *strict*
complementarity (switching function `S != 0` ⇒ nonzero bound multiplier) those
directions leave the critical cone and the remaining curvature is real. In
practice the solver parks the throttle **near but not at** its bounds
(`s ~ 0.02 / 0.98`, edge ~99%), the bound constraints are not active, the
degeneracy is never removed, and the test is inapplicable. This is structural
for an ε→0 homotopy that approaches bang-bang rather than imposing it.

### M3 — forward flow across many revs destroys the base point
The correct instrument for a bang-bang extremal is the Maurer–Osmolovskii
reduced problem over the *k switch times*. The tulip built it as a
**forward-flow** formulation and it is blocked: re-integrating the exact
direct control from `x0` across ~40 revs diverges by `||r||~3`, `||v||~5`,
`|t|~10`. The collocation solution satisfies the *discrete* defects but does
not correspond to a nearby *continuous* trajectory under forward integration,
so the reduced problem's base point is infeasible (`||c|| ~ 20`), the reduced
gradient is not ~0, and any Hessian there is meaningless.

To its credit `psr_switch_hessian.m` **detects** this via a `baseFeas`
self-check and returns BLOCKED rather than a bogus certificate. This is the
same conditioning wall that defeats indirect shooting (`ifs`, `ms_band`) — it
is a property of long multi-rev flows, not of the second-order test.

---

## 3. Ruled out — do not redo

- **Tuning `inertiaZero` to rescue the NLP test.** `ldl`/Gould pivot-sign
  inertia is unreliable on this near-singular KKT (mis-signs 56 pivots at
  zt=1e-9); its correct window is *disjoint* from `eig`'s, so no single
  threshold serves both. `eig` is the gold standard and `inertiaZero=1e-9` is
  settled. Reverting this has already been tried and undone once.
- **Expecting strict SOSC from the NLP for a min-fuel row.** See M1.
- **Any forward-shooting formulation of the switching-time problem at
  multi-rev scale.** See M3.
- **Planning a fresh "tier 1/2/3" program per campaign.** Tiers 1 and 2 are
  already built (earth, tulip) and their verdicts are in §1.

---

## 4. Open leads

**LEAD-1 (high value, cheap) — the SOSC "recovery wall" at 1 / 0.5 N may not
exist.** `earth_elliptic_to_geo/TODO.md` item 0 records two distinct scale
walls, the first being "*recovery* — the warm re-solve itself fails at 1 N /
0.5 N (n≈16.5k)", and proposes building a scalable warm-resolve recovery.
**Contradicting evidence (2026-07-25):** `refresh_duals_mee.m` performs
essentially that same warm re-solve (same `warmTight`, same `maxIter`, same
`returnModel` rebuild) and it *succeeded* at both rungs — 1 N `Solve_Succeeded`
defect **8.27e-14**, 0.5 N defect **1.83e-13** — comfortably inside
`sosc_recover_kkt`'s `recoverOK` bar of `maxDefect < tol.feas = 1e-8`.
If that reproduces under `sosc_recover_kkt` itself, the ERROR verdicts come
from *downstream* KKT assembly / dense null-space at scale (the second,
acknowledged wall), and the proposed recovery fix is aimed at a non-problem.
**Test first, build second.**

**LEAD-2 — certificates must be evaluated at tight optima.** Six of nine earth
rows improved on warm re-solve (5 N +0.235 kg, 2.5 N +0.51, 1 N +0.44, 1 N PSR
+0.47, 0.5 N +0.02, 0.5 N PSR +0.02); only 10 N / 0.2 N / 0.1 N were already
tight. A second-order certificate computed at an under-converged point
certifies nothing. Re-solve the six before certifying them.

**LEAD-0 (cheapest win on the board) — port the IPOPT native-inertia
certificate to the other three campaigns.** It is the one second-order
instrument that actually delivers, and it is nearly free: it only reads the
δ_w regularization history out of IPOPT's own stats (`casadi_minfuel_sundman`
already captures it; `psr_ipopt_certify` only interprets it). Earth 2-body,
earth CR3BP and ELFO have no equivalent and would gain a weak-local-min
verdict per row for a few hours of plumbing. Note it also sidesteps the
conditioning wall that defeats our own factorization (M3 / `psr_second_order`'s
tolEig noise floor), because MUMPS tests inertia on the well-scaled system at
every iteration.

**LEAD-3 — M2 might be removable by construction.** If the throttle were
*pinned* to its bounds on a frozen arc structure (rather than approached via
ε→0), strict complementarity would hold and the KKT-inertia test might apply.
This is close to what the switching-time reduction does anyway, so it may not
be worth a separate build — but it has never been tried and it would reuse
existing machinery.

---

## 5. Decision — what to build next, and why

**Build the STM / multiple-shooting switching-time Hessian.**

`psr_switch_hessian.m`'s own header already specifies it: keep the collocation
trajectory as the (feasible) base, obtain switch-time sensitivities from the
state-transition matrix,

    dPsi/dsigma_i = Phi(sigma_f, sigma_i) * [ f_burn - f_coast ]

integrated along that base, and form the second variation with **segment
matching**, so no quantity is forward-propagated across all revs. Multiple
shooting defeating a forward-shooting conditioning wall, exactly as it does in
the indirect campaigns.

Rationale:
1. It is the **only live path to STRICTNESS** for a bang-bang solution. The
   IPOPT native-inertia certificate (§1) already gives *weak* local minimality
   on 12/17 ε=0 tulip rows, and LEAD-0 extends that cheaply to the other
   campaigns — but a weak minimum is exactly what M1 predicts and it cannot be
   sharpened by any NLP-level test (M1 kills the reduced-Hessian route, M2
   kills its repair). Strictness lives in the switching times.
2. It is a **build, not a research problem**: the formulation is already
   written down, and the blocked forward-flow version supplies the arc
   bookkeeping, the `baseFeas` self-check, and the output contract.
3. It is **campaign-agnostic**: the same second variation serves earth 2-body,
   CR3BP, tulip and ELFO — one instrument closes the Tier-2 row for all four.
4. Its STM/variational core is **shared with the conjugate-point test**, so it
   partly pre-builds Tier 3 and dovetails with the indirect goal (6).

Sequencing note: do LEAD-2 (re-solve the loose rows) before certifying
anything, and LEAD-1 (one cheap re-test) before investing in SOSC recovery.

**Sufficiency claim to aim for.** A positive-definite projected switching-time
Hessian, *together with* (i) strict bang-bang `Sdot != 0` at every switch and
(ii) primer / strengthened Legendre–Clebsch optimality of the thrust direction
— both already reported by the first-order verifier, and now trustworthy since
the dual fix — is the Maurer–Osmolovskii sufficient condition for a bang-bang
local minimum. That is the target statement for the methods paper.

---

## 6. Experiment log (both orders)

Append one entry per experiment — first- or second-order: date, instrument,
campaign/row, verdict, and what it changed in Part A or §1–§5.

| date | instrument | target | verdict | consequence |
|---|---|---|---|---|
| 2026-07-12 | NLP SSOSC (KKT inertia) | tulip 1.15× refined | NOT APPLICABLE | established M2 |
| 2026-07-12 | switching-time Hessian (forward flow) | tulip 1.15× refined | BLOCKED (`baseFeas` ~20) | established M3; specified the STM fix |
| 2026-07-19 | reduced-Hessian SOSC | earth 10 N | WEAK_MIN (270 flat) | established M1 |
| 2026-07-19 | reduced-Hessian SOSC | earth 5 / 2.5 N | INCONCLUSIVE | genuine precision limit, not a bug |
| 2026-07-19 | reduced-Hessian SOSC | earth 1 / 0.5 N | ERROR / scale-skip | recorded as a recovery wall — **disputed, see LEAD-1** |
| (ongoing) | **IPOPT native inertia (δ_w)** | tulip PSR, 137 rows | **12/17 certified at ε=0**; 100/120 at ε>0 | second order is not at zero; weak local min for bang-bang |
| 2026-07-25 | first-order PMP (primer/sign/alignment) | earth 9 rows + CR3BP 10/5 N | **PASS 0.000° / 100%** | Part A row 1–3 green for earth; makes Maurer ingredients (ii)/(iii) trustworthy |
| 2026-07-25 | full KKT certificate from raw `lam_g` | earth 10 N | ‖∇ₓL‖∞ 1.5e-14, tangential ∂L/∂β 8e-17 | Part A row 6; also root-caused the dual bug |

---

## 7. Pointers

| file | role |
|---|---|
| `earth_elliptic_to_geo/process/DESIGN_sosc.md` §11–12 | NLP SOSC method evolution, `eig` vs `ldl`, threshold rationale |
| `earth_elliptic_to_geo/process/PLAN_sosc.md` | original build plan |
| `earth_elliptic_to_geo/direct/verify/sosc/` | the NLP-level implementation |
| `GTO_tulip/direct/PSR/psr_second_order.m` | KKT-inertia test + M2 finding |
| `GTO_tulip/direct/PSR/psr_switch_hessian.m` | Maurer test + M3 finding + the STM fix spec |
| `earth_elliptic_to_geo/process/LESSONS_DUAL_EXTRACTION.md` | why the first-order side is now trustworthy |
| `earth_elliptic_to_geo_CR3BP/TODO.md` | CR3BP tier plan (supersede with §1 before acting on it) |
