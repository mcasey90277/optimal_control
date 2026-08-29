# The orbit-transfer problem space available through pumpkyn

**Purpose:** define the full set of transfer problems reachable with the
existing pumpkyn/pumpkynPie orbit catalogs, so the goal — every problem, both
methods, all three flavors — has a concrete scope.

Surveyed 2026-07-31 against pumpkyn `5f5ca31` / pumpkynPie `47be599`; status tables refreshed 2026-08-23 (catalogs 3–6 shipped, min-energy pilot).

---

## Nine catalogued CR3BP orbit families

Every one has a getter *and* a stored catalog (`pumpkynPie/data/*.mat`), so
endpoints come from a lookup plus `cont_np` refinement rather than a fresh
solve. All are Earth-Moon.

| family | getter | parameters | notes |
|---|---|---|---|
| **Tulip** | `pumpkyn.cr3bp.getTulip(tau0, Np, pm, tol)` | period, petal count, hemisphere | our current target; southern lobes for lunar south-pole PNT |
| **Pumpkin** | `getPumpkin(Np, pm, ...)` | petal count, hemisphere | the namesake family; **sidereal resonance — period is exactly 2π** |
| **DRO** | `getDRO(tau0, ...)` | period | distant retrograde; the DRO→tulip demo's departure orbit |
| **DPO** | `getDPO(tau0, ...)` | period | distant prograde |
| **LPO** | `getLPO(tau0, ...)` | period | **L**ow **P**rograde **O**rbit (not "Lagrange point orbit") |
| **Halo** | `getHalo(tau0, Lpt, pm)` | period, Lagrange point, N/S | the classic station-keeping target |
| **Lyapunov** | `getLyapunov(tau0, Lpt, pm)` | period, Lagrange point, N/S | planar |
| **Axial** | `getAxial(tau0, Lpt, pm)` | period, Lagrange point, N/S | |
| **Cycler** | `getCycler(tau0, resStr, ...)` | period, **resonance ratio** | resonant transfer orbits |

Plus `lagrangePts(mu)` for the five equilibria, and `manifolds(tau0, x0, muStar,
epsilon)` for stable/unstable manifolds of any unstable member.

**Earth-side endpoints** come from our own work, not pumpkyn's catalog: GTO and
GEO, via `insertion_states` / `gto_tulip_endpoints` / `gto_elfo_endpoints`, and
ELFO (lunar elliptical frozen orbit).

---

## What pumpkyn can currently *solve*

**Min-time only, and it is endpoint-agnostic.**

```matlab
sol = pumpkyn.cr3bp.tfMin(rv0, rvf, lambda0_tf, Tmax, c, muStar)
```

`tfMin` takes two arbitrary rotating-frame states. It does not know or care
which family they came from. So **every pair drawn from the nine families above
is already a solvable min-time problem** — the only thing standing between us
and any of them is an initial costate guess, which is exactly what
[[goal-costate-catalog]] is for.

There is **no min-energy or min-fuel solver in pumpkyn.** Those are ours:
`casadi_minfuel_sundman` (tulip), `casadi_energy_freetf` (ELFO, free-t_f
two-primary), `casadi_lt_mee` (earth L-domain). The energy→fuel homotopy is
ours too.

**Also available and unused:** `directLambert` / `minLambert` (impulsive
single-revolution transfers in CR3BP) — useful as the high-thrust limit that
anchors a thrust walk-down, and as an independent sanity check.

---

## The goal matrix

For each transfer problem: **2 methods × 3 flavors = 6 cells.**

### Currently in play

> **Updated 2026-08-23.** The case-by-case status view (entry counts, thrust
> rungs, verification ladder, roadmap) now lives in `../STATUS_AND_ROADMAP.md`;
> this table is kept as the goal-map summary.

| problem | direct | indirect | status |
|---|---|---|---|
| GTO → tulip | ✅ all three | partial (min-time; min-fuel stalled at IFS) | most complete campaign at the 25 mN / ~40-rev flagship regime (no costate catalog there yet). **Separate min-time catalog SHIPPED-pending-review** at the Darin-standard 150 kg / 1–15 N few-rev regime (`GTO_tulip/catalog/`, 2026-08-29): 2,625 entries, 16 sheets, 840/1,152 pairs (73%), conjugate 2,625/0/0 — first catalog on a non-periodic (`'gto'` pseudo-family) departure orbit; 3–1 N deep-rung leg is a measured closure wall, split off open |
| GTO → ELFO | ✅ min-fuel front + min-time anchor | not started | |
| elliptic → GEO (2-body) | ✅ certified ladder | external cross-check only (MfMax v0+v1 built & validated); our own MATLAB indirect not built | |
| elliptic → GEO (CR3BP) | ✅ certified ladder | not started | |
| DRO → tulip | ✅ **min-time CERTIFIED**, and a **NEW faster basin found+certified 2026-08-03: t_f=3.8170 vs 4.0152, periselene 3954 km, cold lineage** | ✅ min-time (Darin's demo — now known to be a beaten local extremal) | **the two methods AGREE**: direct t_f = 4.0152501 vs indirect 4.0152425, 5 sig figs, 3.3 m worst position error. Costate comparison also DONE (see below). **Min-time catalogs SHIPPED** (deliverables 2+3: coarse 4×4×6×6 + fine 12×12 sheet). **Min-energy fixed-t_f pilot 2026-08-14: 5/5 cells pass** |
| HALO(L2 S) → tulip | ✅ (pipeline direct leg) | ✅ min-time | **catalog SHIPPED** (deliverable 4: 3,980 entries, 92% pairs) |
| DPO → tulip | ✅ (pipeline direct leg) | ✅ min-time | **catalog SHIPPED** (deliverable 5: 3,932 entries, 89% pairs) |
| L1 halo ↔ L2 halo | ✅ (pipeline direct leg) | ✅ min-time | **catalogs SHIPPED, BOTH directions** (deliverable 6, schema v2 arrival-period axis: 1,952 + 2,096 entries; directions measurably asymmetric) |

### The expansion

Pairs among the nine families. Not all are equally interesting; a sensible
first tranche, ordered by expected tractability (few revolutions first, since
revolution count — not the objective — is what decides whether indirect works):

1. ~~**DRO → tulip**~~ — **DONE** (direct twin built, methods agree, catalogs
   shipped as deliverables 2+3; min-energy pilot passed 2026-08-14).
2. **GEO → DRO** — named by Darin on the call as a catalog example. *Still
   open — now the top untouched pair.*
3. **DPO → tulip** — **DONE** (deliverable 5). **LPO → tulip** still open.
4. **DRO → Halo**, **DRO → NRHO-like Halo** — still open. (Halo↔halo L1↔L2
   is DONE — deliverable 6 — which partly de-risks these.)
5. **Lyapunov ↔ Halo** at the same Lagrange point — a well-studied benchmark
   with published comparisons, useful for external validation. Still open.
6. **Cycler → anything** — resonance structure makes these distinctive. Still
   open.

*(HALO → tulip, not on the original tranche list, is also DONE — deliverable
4 — and produced the cheapest transfer in any catalog: 0.6546 km/s.)*

---

## Why DRO → tulip was first — and what it returned

It was the only problem in the space where **the indirect solution existed and
the direct one did not.** Everywhere else the direct method led. That made it
the natural place to test the coordination question from `direct_vs_indirect.md`
in reverse: does a direct solve reproduce the known indirect answer?

**Answered 2026-08-02: yes, but only at fourth order.** Hermite-Simpson at
N = 1600 gives t_f = 4.0152501 against the indirect 4.0152425 — five significant
figures — with a worst-interval POSITION error of 3.3 m and 42 m end to end. Second-order trapezoidal
collocation never got closer than 3% and was inaccurate by 1,100-12,600 km while
reporting defects of 1e-14.

Two lessons generalize beyond this pair:

1. **A machine-tight defect is a statement about the discretization, not the
   trajectory.** Measured gap on this problem: 1e7. Every campaign here quotes
   1e-14 defects as evidence of a good solve; that evidence is necessary and
   nowhere near sufficient. Measure the continuous residual.
2. **Refinement is not monotone.** N=400 -> 800 made the answer worse by every
   measure (and put a node 719.6 km inside the Moon); only N=1600 landed on the
   reference. Any single mesh density, including the one with a 1e-14 defect,
   would have misled. NOTE: we twice concluded from such rows that the
   unconstrained problem is ill-posed, and twice the evidence was discretization
   error — the N=800 solve has a 3127 km position error, so it is not a
   trajectory. The conjecture is open; a path constraint is still prudent
   engineering.

**The costate comparison is DONE (2026-08-02) and it succeeded.** The direct NLP
duals reproduce the indirect costates: primer exact to 1.2e-06 deg, lambda_v and
lambda_r directions agreeing to ~0.0006 deg median, and a scale factor of
0.99999 — the same normalization, not merely proportional — and the state
histories agree to 0.396 km peak once the 2.9 s t_f phasing offset is removed.
So the two methods find the same EXTREMAL, not just the same final time. See
DRO_tulip/FINDINGS.md and direct/results/dvi_N1600_{state,costate}.png.

---

## Practical notes

- **Endpoint construction is two steps.** The getter returns a *catalogued*
  state; `cont_np` refines it to a closed periodic orbit. Both are needed —
  `gto_tulip_endpoints` does exactly this.
- **Arrival phase is a free choice.** The DRO→tulip demo picks the arrival state
  on the tulip by maximizing a velocity-angle criterion. Different choices give
  different transfers; the phasing is a design variable, not a given.
- **`tfMin` returns no convergence flag.** Verify by re-propagating and checking
  its four terminal conditions — see `abstracts/data/bht1500_continuation.m`.
- **Watch the revolution count.** Few revs → indirect leads. Tens of revs →
  direct leads, and single shooting will not close regardless of the guess.
