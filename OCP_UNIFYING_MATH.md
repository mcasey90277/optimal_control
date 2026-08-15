# The One Problem This Repository Solves, Many Times

Every campaign folder here — `collocation_examples`, `mpc`, `orbit_transfer`,
`booster_landing`, and (once its Phase-2 optimizer exists) `missiles` — is an
instance of the **same mathematical object**, solved by one of a small number
of routes, verified by the same discipline. This document states that object
once, maps every folder onto it, then does the more useful thing: explains
**which physical property of each problem forces which special machinery**, so
that "Sundman regularization" or "multiple shooting" or "lossless
convexification" stop looking like folder-specific lore and start looking like
predictable responses to identifiable features. It closes with what this view
implies for the common-library goal.

---

## 1. The unifying object: the Bolza optimal control problem

Choose a control history u(t) (and possibly the final time t_f) to

```
minimize    J = φ(x(t_f), t_f) + ∫₀^{t_f} L(x, u, t) dt        (cost)
subject to  ẋ = f(x, u, t)                                     (dynamics)
            ψ₀(x(0)) = 0,   ψ_f(x(t_f), t_f) = 0               (boundary)
            g(x, u, t) ≤ 0                                     (path)
            u(t) ∈ U                                           (control set)
```

Every folder is a choice of (x, f, L, φ, ψ, g, U, free-or-fixed t_f):

| folder / campaign | state x | dynamics f | cost | control set U | t_f | boundary flavor |
|---|---|---|---|---|---|---|
| `ex1_block_move` | [pos; vel] (2) | double integrator | ∫u² (energy) | ℝ (unbounded) | fixed 1 s | fixed endpoints |
| `ex2_cart_pole` | [q₁ q₂ q̇₁ q̇₂] (4) | Lagrangian EOM, underactuated | ∫u² | box \|u\|≤40 N | fixed 5 s | fixed endpoints |
| `mpc_cart_pole` | same (4) | same | quadratic tracking (Q,R) | box | receding horizon | current state → reference |
| `min_energy_tutorial` | orbit state (6) | two-body | ∫\|u\|² | ℝ³ | fixed | fixed endpoints |
| `earth_elliptic_to_geo` | MEE + mass (7) | two-body + thrust | fuel ∫\|u\| (via throttle) | ball ‖dir‖=1, throttle∈[0,1] | fixed multiple of t_fMin | orbit → orbit |
| `earth_..._CR3BP` | same + lunar term | 2-body + third body | fuel | ball+throttle | fixed multiple | orbit → orbit |
| `GTO_tulip`, `GTO_ELFO` | Cartesian/Sundman (7) | CR3BP + thrust | fuel (via energy homotopy); min-time anchors | ball+throttle | fixed multiple / free | orbit → orbit |
| costate catalogs (`DRO/HALO/DPO_tulip`, `HALO_HALO`) | [r v m] (7) | CR3BP + thrust | **min-time** (all-burn) | ball, throttle ≡ 1 | **free** | orbit-phase → orbit-phase |
| `booster_landing` | [r v m] (7) | const gravity (+opt drag) + thrust | fuel | **annulus** T∈[Tmin,Tmax] + cones | free | state → landing set |
| `missiles` (Phase 2, future) | 3-DOF point mass (6–7) | spherical rotating Earth + aero | TBD (range/heat/time) | bank/AoA schedules | free | launch → impact |

Estimation folders (`lieFiltering`, `gauss_sum_curvature`) and the optimizer
tutorial (`quasiNewton_matlab`) are deliberately outside this table: the first
two are the *dual* problem (inferring x from data rather than choosing u), and
the third is the NLP engine underneath route 1.

---

## 2. The solution routes (all of them live in this repo)

**Route 1 — Direct (discretize, then optimize).** Choose nodes, write the
dynamics as *defect constraints* between adjacent states — trapezoidal
(`collocation_examples`) or Hermite–Simpson (everything research-grade) — and
the cost as a quadrature; hand the resulting sparse NLP to fmincon
(tutorials), IPOPT/CasADi (orbit, booster), or quadprog-class solvers (MPC).
The KKT conditions of this NLP are a *discretization of the PMP conditions* —
that fact is the bridge below.

**Route 2 — Indirect (optimize, then discretize).** Apply Pontryagin: adjoin
costates λ, form H = L + λᵀf, get λ̇ = −∂H/∂x, u* = argmin_u H, plus
transversality (e.g. λ_m(t_f)=0, H(t_f)=0 for free-t_f min-time). The OCP
becomes a two-point BVP solved by shooting (`pumpkyn tfMin`, the
`min_energy_tutorial`) or **multiple shooting** (`costate_common/ms_bvp` /
`ms_tfmin`) when single shooting's sensitivity blows up. Rewards: machine-
precision extremals, switching structure by construction. Price: you must
*guess costates*, which have no physical intuition — hence this repo's costate
catalogs, whose entire purpose is to make route 2 startable.

**Route 3 — Convex reformulation (booster only, and only because the physics
allows it).** Lossless convexification turns the nonconvex thrust annulus
into an exact SOCP: global optimum, no initial guess, polynomial time. This is
not an approximation — for vacuum point-mass PDG the relaxation is provably
tight (gate G4 checks ‖u‖=σ). It exists in exactly one folder because its
preconditions (linear-ish dynamics, norm-type control set) hold in exactly one
folder. When drag turns on, the convex twin *disappears* and the campaign
falls back to routes 1+2 cross-checks — a clean demonstration that route
availability is a property of the problem, not a preference.

**Route 4 — Closed loop (re-solve or locally approximate).** `mpc` re-solves
a route-1 problem every 50 ms from the measured state. `booster_landing`
instead solves once and wraps the trajectory in TVLQR (gains from the
differential Riccati equation about the optimal trajectory) and validates the
loop by Monte Carlo. Same underlying OCP; the difference is *when* you pay the
optimization cost.

**The bridge that makes them one subject:** the KKT multipliers of route 1's
defect constraints are **samples of route 2's costates** (the covector mapping
— `costate_common/duals_to_costates`). This repo uses the bridge constantly:
direct solves *seed* indirect refinement (the catalog pipeline), PMP structure
*verifies* direct solves (`verify_common/foc_check`, booster gate G5), and PMP
switching times *steer* direct mesh refinement (PSR). Which station a
multiplier belongs to is scheme-specific — Hermite–Simpson duals live at
interval **midpoints** — and getting it wrong is this repo's most-repeated
bug (see §5).

---

## 3. The differences, and the physics that forces each one

This is the actual content. Each special technique in this repo answers a
specific measurable property of its problem. If the property is absent, the
technique is unnecessary — that's why the cart-pole never heard of Sundman.

### 3.1 Time-scale inhomogeneity → regularized independent variables

An eccentric orbit spends almost all its period near apoapsis moving slowly
and rips through periapsis. A uniform-in-*time* mesh therefore starves exactly
the arc where dynamics (and thrust effectiveness) change fastest. Two cures,
both in `orbit_transfer`: **Sundman regularization** (dt = c·r^{3/2} ds — mesh
uniform in an anomaly-like variable, used by the GTO campaigns) and **MEE with
longitude ΔL as the independent variable** (earth→GEO ladder). The cart-pole's
timescale is homogeneous over its 5 s window; the block move is a polynomial;
the booster burn lasts ~17 s of similar-magnitude dynamics — none of them
needs this. **Diagnostic: ratio of fastest to slowest state-velocity along the
trajectory.** Orbit: 10²–10³. Cart-pole: ~1.

### 3.2 Dynamic sensitivity → multiple shooting, warm-start families

CR3BP trajectories amplify initial-condition errors like e^{λt} with several
revolutions of accumulated growth: measured here, a collocation-quality
costate seed (~1e-3) misses the target by **36,000–560,000 km** when
single-shot over the full arc. Splitting the arc into K short segments
(`ms_bvp`) caps each segment's amplification and turns the same seed into a
few-iteration convergence to 1e-13. The block move (1 s, double integrator,
zero Lyapunov growth) is solvable by single shooting from almost any guess.
The same sensitivity is why catalog campaigns *walk* thrust down warm-started
(family continuity) rather than cold-solving each rung: cold solves land on
different solution families 50% apart in t_f. **Diagnostic: STM norm growth
over the horizon.**

### 3.3 Cost smoothness → homotopy, switch-aware meshes, or nothing

∫u² (block, cart-pole, min-energy tutorial) gives u* smooth and interior —
any method converges. ∫|u| (min-fuel: GTO campaigns, booster) gives
**bang-bang** control: the discretization must represent discontinuities
(mesh refinement at switches — PSR; Sundman helps by putting nodes where
switches cluster), and the solve landscape is nonsmooth enough that the GTO
campaigns approach fuel through an **energy→fuel ε-homotopy** (solve the
smooth problem, then slide ε→0). Min-time with low thrust (the catalogs) is
the easy extreme in disguise: the optimal control is *all-burn* (provable
from λ̇_m ≤ 0 + transversality — continuous burn is a theorem, not an
observation), so there are no switches at all, which is precisely why direct
+ multiple shooting handles thousands of catalog cells unattended.
**Diagnostic: is L strictly convex in u (smooth), or linear in throttle
(bang-bang), or is time itself the cost (structure depends on U)?**

### 3.4 Landscape multiplicity → continuation ladders and basin discipline

Orbit problems admit *families* of local optima (different winding counts,
different geometries): the repo has measured 19-vs-24-switch flips on a 2e-5
perturbation (10 N min-fuel), "blocky" t_f maps that are family walls, and a
razor-thin basin lesson (keep-best-mass). The response is always
**continuation**: thrust ladders, μ-continuation (2-body → CR3BP), t_f grids,
ε-homotopies — walk a parameter from an easy problem to the hard one, staying
in one basin. The cart-pole swing-up has multiplicity too (swing left vs
right, extra pumps) but at tutorial scale a relaxed-bounds warm start
suffices — which is itself a two-rung ladder, the same idea at miniature.

### 3.5 Control-set geometry → three treatments of the same constraint

The thrust set {‖dir‖ = 1, 0 ≤ throttle ≤ 1} appears across the repo and is
handled three ways, by need: (a) **cone elimination / direct
parametrization** (GTO direct solves), (b) **argmin-H structure** in indirect
form — the primer vector u* = −λ_v/‖λ_v‖ removes the constraint analytically,
(c) **lossless convexification** (booster) when a global certificate is worth
the reformulation. The booster adds genuinely different geometry — a thrust
*annulus* (min throttle > 0: engines can't idle) and glide-slope/pointing
cones — which is what makes its landing-set boundary conditions and its
τ→touchdown singularities (v_f = −1.5 m/s, not 0) problem-specific rather
than transferable.

### 3.6 Free final time → one standard trick, everywhere

Fixed-t_f problems (tutorials, MPC) discretize [0, T] directly. Every
free-t_f problem here (min-time catalogs, booster, min-fuel with t_f as a
multiple of t_fMin) uses the same **time-dilation lift**: normalize to
s ∈ [0,1], make t_f a decision variable that multiplies the dynamics, and (in
route 2) inherit the transversality condition H(t_f) = 0. The lifted-time
copies must agree (gate G8 in the orbit certifier) — a repo-specific
consistency check born from this lift.

### 3.7 Units → nondimensionalization when scales span decades

Space problems run in canonical units (CR3BP μ*, l*, t*; catalog thrust in
N converted to ND acceleration) because raw SI spans ~10⁹ across the state
and defeats NLP scaling heuristics (measured lesson: manual `scaleNLP`
*hurt* — it fought IPOPT's own scaling; ND units + solver defaults won).
Cart-pole SI numbers span ~10²: no need.

### 3.8 What each folder contributes that the others genuinely don't

- `collocation_examples`: the transcription itself, small enough to read in
  one sitting, with an **analytic answer** (u* = 6−12t) as ground truth.
- `mpc`: feedback via re-solving; the only folder where solve *latency* is a
  constraint.
- `min_energy_tutorial` / `lambert`: route 2 at teachable scale; impulsive
  limit.
- `earth_elliptic_to_geo(+CR3BP)`: published-benchmark reproduction; deep
  thrust ladders; the dual-extraction lessons (`opti.lam_g`, never
  `opti.dual`).
- `GTO_tulip/ELFO`: bang-bang at scale; Sundman; homotopy; PMP-steered
  refinement; the near-min-time conditioning wall (still open).
- costate catalogs: industrialized route-1→route-2 pipeline, three-gate
  verification, second-order (conjugate-point) checking, schema'd data
  products for an external consumer.
- `booster_landing`: the convex route; closed-loop certification (TVLQR +
  Monte Carlo); constraint geometry (annulus, cones); "fly the
  discretization's own control representation" as a certification principle.
- `missiles`: the reminder that a *validated forward model with a control
  schedule* is the prerequisite layer — it is the (x, f, U) columns of the
  table with the optimizer not yet attached.

---

## 4. The shared verification discipline (already unified in spirit)

Independently, every research folder converged on the same three-layer
verdict structure — this is the repo's real house method:

1. **Discrete feasibility is not accuracy.** A 1e-14 defect proves only that
   the discretization's equations hold. (Verbatim "house lesson" in both the
   orbit certifier and `booster_landing/certify_pdg`.)
2. **Fly the control.** Reconstruct the *transcription's own* control
   representation (quadratic through node–midpoint–node for HS — not a
   convenient interpolant) and integrate once, end to end; measure where you
   actually arrive. Orbit: G1b / `flown_control_error` (<1 km, catalogs
   0.00 km). Booster: G2, plus the hard-won lesson that flying a *different*
   control representation than the one the NLP optimized silently degrades
   the answer.
3. **Independent-route agreement.** Direct vs indirect (catalog tfMin
   acceptance, |Δz| ~ 1e-9; earth-campaign MfMax cross-check), direct vs
   convex (booster G3), numeric vs analytic (block move), plus PMP-shape
   gates on the duals (foc_check; booster G5) and — newest — a second-order
   Jacobi/conjugate-point necessary condition. Closed-loop folders add Monte
   Carlo as the statistical version of the same idea.

---

## 5. What this implies for the common-library goal

**The strongest single piece of evidence in the repo:** in August 2026 the
booster campaign's G5 gate compared a Hermite–Simpson segment dual against
the *node* control instead of the segment *midpoint* control — the **exact
bug** the orbit campaigns had found by external review and fixed in
`duals_to_costates` days earlier. Two folders, same math, no shared home →
the same subtle error was independently re-invented and independently
re-fixed. That is the one-home-per-rule argument in one sentence.

### Already shared (inside `orbit_transfer`, candidates for promotion)

| library | contents | genericity today |
|---|---|---|
| `costate_common` | `ms_bvp` (generic multiple-shooting engine — problem enters as three closures), `ms_conjugate_test`, `duals_to_costates` (all covector station rules), `harvest_ms_seed`, flown-control verifier + preflight, catalog packager + schema, golden-cell regression | `ms_bvp` and `duals_to_costates` are **already problem-agnostic**; the flown verifier is CR3BP-specific only through its RHS |
| `verify_common` | AD-based first-order gate (`foc_check`/`foc_report`), IPOPT inertia, PMP residual, mesh tools | transcription-native, already runs on four campaigns |
| `cr3bp_common` | CR3BP GTO problem definition | problem data, correctly not generic |

### Cross-folder extraction candidates — with the required discipline

House rule (learned the hard way, five of eight candidates once refused):
**measure before extracting** — an extraction must be justified by two real
consumers whose copies can be diffed, not by aesthetics. With that filter:

1. **Covector station rules** (`duals_to_costates`) — booster's G5 is the
   second consumer *by evidence of the shared bug*. Extraction test: route
   booster G5 through it and reproduce the corrected 1.2e-6-deg primer angle.
   **DONE 2026-08-09** (`oc.duals_to_costates`; G5 primer 1.20742e-06 deg
   identical, orbit golden cells green).
2. **Flown-control verifier pattern** — three implementations exist (orbit
   `flown_control_error`, booster G2, and the tutorials' ad-hoc re-simulation).
   The engine (interval loop + *scheme-matched* control reconstruction +
   tight ODE) is identical; only the RHS differs. Same closure pattern as
   `ms_bvp`. **DONE 2026-08-09** (`oc.fly_control`; orbit globKm and booster
   G2 residuals identical to 0.000e+00).
3. **`ms_bvp` promotion out of `orbit_transfer`** — nothing in it is
   astrodynamical; it would solve the cart-pole swing-up's PMP BVP as a
   tutorial exercise (worth doing purely as a demo that route 2 machinery is
   problem-agnostic).
4. **Transcription defect builders** — trapezoid and HS defect/quadrature
   code exists in at least four places (tutorials, orbit direct libs, booster
   `lib`). High duplication, but interfaces differ (fmincon vectors vs CasADi
   Opti); extraction should wait for the next *new* campaign and be built for
   it, then adopted backward if the diff supports it.
5. **The campaign patterns as templates, not code**: front door with an
   ADJUSTABLE PARAMETERS block; certify module with lettered gates and an
   honest verdict; batched hang-proof driver + liveness monitor; golden-cell
   quality regression; process/ + doc/ record-keeping. These transfer as
   *skeleton + checklist* (a `matlab-campaign` scaffold), which is how
   booster and missiles already resemble the orbit campaigns without sharing
   a line.

### What should *not* be merged

Forward models (`missiles`, pumpkyn's CR3BP, booster's gravity+drag) stay
per-domain — physics is problem data, and the pumpkyn boundary ("call it,
never write into it") has been the most successful interface in the repo.
MPC's real-time loop and the estimation folders solve different problems and
should only ever share the *verification* vocabulary.

### Suggested reading path (for a newcomer, or us in six months)

1. `collocation_examples/ex1_block_move` — route 1 with an analytic answer.
2. `min_energy_tutorial` — route 2 on the same class of problem.
3. `orbit_transfer/DRO_tulip/doc/costate_library_methodology.tex` — both
   routes joined by the covector bridge, at research grade.
4. `booster_landing/README.md` + `certify_pdg.m` — route 3, closed loop, and
   the verification discipline stated independently.
5. This file's §3, backwards: given a *new* problem, read off its diagnostics
   (time-scale ratio, STM growth, cost smoothness, landscape multiplicity,
   control-set geometry, free t_f?, unit span) and it names your machinery
   before you write a line.

---

*M. Casey / assistant synthesis, 2026-08-09. Sources: the campaign READMEs,
certifiers, and process docs of every folder named above; the measured
numbers are quoted from their records, not re-derived here.*
