# The orbit-transfer problem space available through pumpkyn

**Purpose:** define the full set of transfer problems reachable with the
existing pumpkyn/pumpkynPie orbit catalogs, so the goal — every problem, both
methods, all three flavors — has a concrete scope.

Surveyed 2026-07-31 against pumpkyn `5f5ca31` / pumpkynPie `47be599`.

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

| problem | direct | indirect | status |
|---|---|---|---|
| GTO → tulip | ✅ all three | partial (min-time; min-fuel stalled at IFS) | most complete |
| GTO → ELFO | ✅ min-fuel front + min-time anchor | not started | |
| elliptic → GEO (2-body) | ✅ certified ladder | reference only | |
| elliptic → GEO (CR3BP) | ✅ certified ladder | not started | |
| DRO → tulip | ✅ **min-time CERTIFIED** (2026-08-02) | ✅ min-time (Darin's demo) | **the two methods AGREE**: direct t_f = 4.0152501 vs indirect 4.0152425, 5 sig figs, 3.3 m worst position error. Agreement is on t_f + endpoints; a trajectory/costate comparison is still to do |

### The expansion

Pairs among the nine families. Not all are equally interesting; a sensible
first tranche, ordered by expected tractability (few revolutions first, since
revolution count — not the objective — is what decides whether indirect works):

1. **DRO → tulip** — already has an indirect min-time solution. Build the
   direct twin and the other two flavors. *Best first target: it validates the
   direct↔indirect handoff on a problem where the indirect answer is already
   known.*
2. **GEO → DRO** — named by Darin on the call as a catalog example.
3. **DPO → tulip**, **LPO → tulip** — same target, different departure.
4. **DRO → Halo**, **DRO → NRHO-like Halo** — the operationally common pair.
5. **Lyapunov ↔ Halo** at the same Lagrange point — a well-studied benchmark
   with published comparisons, useful for external validation.
6. **Cycler → anything** — resonance structure makes these distinctive.

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

Still open: the costate comparison itself. The certified direct solution is the
input it always needed, and it has never been possible before now.

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
