# Missile Trajectory Library — Open Items

*Coorbital, Inc. — recompiled 2026-08-08 against commit `0743a7d`. The previous
edition was compiled against `8ea960d` and revised once on 2026-08-07; five
substantial commits have landed since, and several of its items are now closed
or were justified on a claim that turned out to be wrong.*

Every claim below was checked against the code, or measured by running it, on
2026-08-08; the measurements and their sources are named inline. The suite was
green at the time of writing: **21 passed, 0 failed** (`tests/run_tests`,
re-run 2026-08-08).

---

## Validation baseline — what is actually anchored to something outside this repo

New in this edition, and the most important thing in the file. Until 2026-08-08
every check in the suite compared the library against **itself** — a closed form
built from the same constants, a reduction of one of its own files onto another.
`docs/LESSONS_LEARNED.md` has three separate entries on why that class of check
cannot see a wrong constant. There are now three checks that do not share
anything with the model, and they are the floor a future change must not fall
below.

| External check | Reference | Measured | Agreement |
|---|---|---|---|
| **Space Shuttle entry** flown through `HGV/run_glide` — 121.9 km, 7400 m/s, γ = −1.2°, β = 500 kg/m², L/D = 1.0, to 0 km | published ≈ 1.5 g peak, ≈ 30 min | **1.426 g**, **28.57 min** (1713.9 s) | −5 % on peak load, −5 % on duration |
| **Minimum-energy loft relation** γ\* = 45° − ψ/4 at 10 000 km | Garwin: "closer to 22° than 45°" | ψ = 89.8315° ⇒ **γ\* = 22.5421°** (`rE` = 6378.137 km) | reproduces the published qualitative result |
| **Flown minimum-energy ballistic arc** against the classical closed form — Plesetsk→New York, 7132.320 km, on a booster sized to the range | 1205.989 km apogee, 26.0120 min | **1257.452 km**, **27.5096 min** | **+4.27 %** apogee, **+5.76 %** flight time |

Reproduce them: the Shuttle case is a `run_glide` override run
(`vehicleFn` returning `mass 92079, Sref 249.9, CL = CD = 0.7369`, which is
β = 500 and L/D = 1); the γ\* arithmetic is one line against
`coorbital.util.missileConst`; the classical-arc comparison is **part 13 of
`tests/test_runBallisticTarget.m`** (`:1040` and `:1098-1120`), where the closed
form is written out in the test file's own hand precisely so a mutated γ\* in
the script is caught rather than mirrored.

**Two of the three are not yet in the repository.** The Shuttle entry and the
Garwin comparison were run and computed for this document; neither is a test,
neither cites its source in-code, and neither will re-run itself. **Action: fold
both into `tests/`** — the Shuttle case as a `run_glide` override test with a
one-sided budget, the γ\* relation as an assertion in
`test_runBallisticTarget`. Add the literature citations while doing it; there is
currently no bibliography anywhere under `missiles/`.

### What is still NOT validated, and must be said whenever the above is quoted

- **Every vehicle and booster parameter is a placeholder.** `vehicle_hgv.m`
  deliberately restates `coorbital.util.vehicleDefaults` value for value and
  says so in its own header; `boosterDefaults` is marked PLACEHOLDER line by
  line. The three checks above validate the *propagator*, not any vehicle.
- **`constLD` holds `CL` and `L/D` constant across the entire flight, far
  outside the hypersonic regime it approximates.** Measured 2026-08-08: the
  shipped `run_glide` case runs **Mach 18.93 → Mach 1.01**; the Shuttle
  validation case above runs **Mach 23.35 → Mach 0.24**; `run_boost_glide`
  impacts at **Mach 1.27**. The scripts print a validity caution below Mach 5,
  which is honest, but the low-speed end of every number in this library is
  produced by hypersonic coefficients used transonically.

---

## The one open discrepancy in the validation set

### Shuttle downrange: 7087 km against a published ≈ 8000 km

Measured 2026-08-08 on the case in the table above: **7087.01 km** of ground
range, against roughly 8000 km for a real orbiter entry from EI to landing.
About 11 % short.

**Bank is the wrong direction and therefore cannot be the explanation.** The
real orbiter S-turns for energy management, and roll reversals *shorten* ground
range by spending lift sideways. Ours flies at zero bank, so if bank were the
difference ours should have gone **farther**, not 900 km shorter.

The suspect is the constant-`L/D` approximation at the low-speed end — the same
caveat as above, here with a number attached. The Sänger relation
`R = (L/D)·rE·ln[1/(1 − (V/Vc)²)]/2` gives 7124 km for this entry state at
L/D = 1, i.e. the propagation is tracking the idealised closed form closely and
**both** are short of the flight article; that points at the aerodynamic model
rather than at the integration. Worth chasing: it is the only external check in
the set that does not close.

---

## Recently closed

Kept rather than deleted — the record of what was fixed is the useful part.

### ~~`BM/run_ballistic_target.m` is not built~~ — DONE, `65aeccd`, reviewed and fixed in `fd72b3d`

Task 5 of `docs/plan_2026-08-07_targeting_and_viz.md` is closed. The script
ranges on the **loft angle**, carries the `'minimum-energy' | 'lofted' |
'depressed'` branch selector the old entry called for, and **measures** which
branch it flew from the flown apogee and flight time rather than trusting the
bracket.

`0743a7d` then replaced `'minimum-energy'`, which had been implemented as "the
branch whose loft angle lies nearer the max-range angle" and was not a
minimum-energy trajectory at all — at full burn both arcs leave burnout with
essentially the same energy (0.55 % apart) and neither carries `V*`. It is now a
genuine **two-parameter** solve, nested rather than a 2×2 Newton so both levels
stay monotonic and both reuse `rangeSolve`: inner on the cutoff fraction, outer
on the loft angle, against achieved range and burnout γ − γ\*. Two new user
parameters, `cutFracMin` and `tolGamDeg` (`BM/run_ballistic_target.m:399`
and `:410`).

### ~~Commit `8ea960d` is unreviewed~~ — DONE, reviewed and fixed in `aba5736`, pinned in `df90c8d`

The review found **one critical** — `run_target`'s advertised `altExag`
override *threw*, because the guard skipped the assignment when the override was
actually supplied — plus five important and six minor. A mutation run over
`8ea960d` then **killed nothing: all seven mutations survived**, because the
gap was not weak assertions but whole features with no coverage at all (the
inset, the arc-frame camera, the tilt sign, `altExag` in any respect). All seven
now fail, pinned by four new sections in `tests/test_viz.m` (12c–12f) and a new
section 8 in `tests/test_runTarget.m`.

**This closes the "review finds what the tests do not" argument with two better
examples than the file previously carried,** and retires the parenthetical that
used to sit here saying no in-repo record of "a review that caught a test which
could not fail" could be found. There are now two:

- the seven surviving mutations over `8ea960d`, above; and
- `fd72b3d`'s critical — `run_ballistic_target`'s headline claim that the branch
  is *measured, not assumed* had **no test behind it**. Replacing the whole
  `measureBranch` call with `flownName = pickName; flownAgree = true` left the
  suite green. Part 11 of `tests/test_runBallisticTarget.m` now flies the one
  geometry where label and measurement can legitimately differ (5030 km required
  on a 50 km tolerance against a 5055.302 km maximum, where `rangeSolve`
  short-circuits at the shared endpoint and both branches return the max-range
  arc itself).

### ~~The GPT-5.6-sol BALLISTIC-PIPELINE review is unapplied~~ — DONE, 2026-08-08

Source document at `docs/reviews/bm_pipeline_gpt56sol_review_2026-08-08.md`;
twelve findings, one critical. What changed, and what it cost:

- **The critical.** `'minimum-energy'` drove the burnout flight-path angle to
  `gamma* = 45 deg - Lambda/4` with `Lambda` the PAD-TO-TARGET central angle.
  That closed form is derived for a free-flight arc with both endpoints on the
  SAME radius; burnout is 82 km up and downrange, so the residual did not apply
  to the arc being flown and the reported agreement verified the wrong
  condition. The mode is now the constrained minimisation it always claimed to
  be — minimise `V_BO^2/2 - mu/r_BO` subject to `R(loft,cutFrac) = R_req` — and
  `gamma*` is a printed DIAGNOSTIC. The objective is checked in the vacuum
  equal-radius limit, where it reproduces `gamma*` to 1.6e-8 rad and `V*` to
  4.7e-16 relative (`test_runBallisticTarget` part 6b).
- **Certified branch machinery.** Golden section returns an interval, not a
  maximiser, so the branch brackets are now `[loftMin,aL]` and `[bL,loftMax]`,
  the branch is identified from the ROOT'S POSITION against `[aL,bL]`, a
  required range inside the unresolved top band is reported as COALESCED, and
  unimodality is certified by adaptive refinement rather than a turn count.
  Apogee and flight time are kept as descriptive outputs only.
- **`alphaMax` moved to `BM/vehicle_bm.m`** as `alphaMaxDeg = 6`, read by BOTH
  BM entry scripts. It used to be 12 deg here and 6 deg there for one airframe,
  and the 12 had been chosen to bring the demonstration target inside the
  depressed branch. The shipped target moved to 62 N 28 W, 4828.045 km, which is
  inside the 6 deg depressed band; `loftMin` moved to -140 deg, which is what it
  takes to hold that branch once the clamp saturates.
- Also: honest cutoff bracketing by sampled sign change with verified
  monotonicity; a measured and printed noise figure for the inner tolerance;
  per-phase completion checks in `flyLoft`; a non-convergence error in the
  golden section; and `run_ballistic`'s drag and lift attribution taken between
  arcs that share ONE gravity model and ONE rotation rate.

### ~~The GPT-5.6-sol pipeline review is unapplied~~ — DONE, `23c3470`, `61e78ad`, recorded in `02ccb11`

Source document archived at `docs/reviews/pipeline_gpt56sol_review_2026-08-07.md`.
Verified in the code 2026-08-08:

| Finding | Where it landed |
|---|---|
| NaN / zero / negative masses passed the mass guard | `+coorbital/+eom/massConstant.m:167-169`, new identifier `coorbital:massConstant:invalidMass` |
| `alpha = theta - gamma` is only the banked pitch relation at σ = 0 | `pitchProgram` now solves `sin θ = R·sin(α+φ)`; raises `:unreachableAttitude`. At 45° bank, θ = 30°, γ = 10° the old law delivered 23.66° of the commanded 30° |
| `rangeSolve` returned the final midpoint and called it the closest reached | best-so-far retained and returned (`+coorbital/+util/rangeSolve.m:315`, `:335`); `fTarget` validated at `:190` (`:badTarget`); midpoint formed as `aLo/2 + aHi/2` and stops on stagnation |
| `constThrust` clamped thrust to zero and froze the mass state | raises `coorbital:constThrust:invalidBackPressure` (`:120`); the shipped margin is better than 7× and is now asserted |
| `greatCircleBearing` returned plausible azimuths for degenerate arcs | raises `:degenerateArc` and `:polarOrigin` (`:178`, `:188`) |
| `phaseRun` forwarded one chain-wide vehicle to every phase | optional `ph.veh` (`+coorbital/+prop/phaseRun.m:247-251`); purely additive, chain-wide argument still the fallback |
| every threshold guard in `glide3DOF` / `boost3DOF` was blind to NaN | finiteness and realness of state, control and every model output checked first, each with its own identifier; assembled `xdot` checked at the end |

**What the same review left open is a live section below** — see *What the
pipeline review deliberately did not change*.

### ~~Both targeting scripts return a solution on a rotating Earth~~ — DONE, `fd72b3d`

They now **refuse** before any propagation, through the existing contract —
empty `traj`, `info.refused = true`, new `info.refusedWhy = 'earthSpin'`,
nothing thrown (`BM/run_ballistic_target.m:840`, `HGV/run_target.m:562`). What
they used to do instead: converge, report a **315.690 km** miss (BM) and
**231.552 km** (HGV) as though it were a solution, print a limitations
paragraph saying rotation was on and then that the initial bearing was the whole
answer, and blame the cross-range on a banked segment while the bank was zero.
All surviving azimuth prose is now gated by an `assert(env.omegaE == 0)` at the
head of the report (`HGV/run_target.m:894`,
`BM/run_ballistic_target.m:1377`), and the unreachable `earthSpin` branches were
deleted rather than left as dead code.

The **underlying limitation stands** and is restated below.

### ~~Movie altitude exaggeration is on by default~~ — DONE, `fd72b3d`

`altExag` ships at **1** — true scale — in both targeting scripts; any positive
number is used as given; the char `'auto'` selects the adaptive rule; anything
else raises before the solve (`HGV/run_target.m:430-440`). `altExag` is now a
`USER PARAMETERS` entry in both, which makes the header's claim that every block
entry is overridable and nothing else is true again.

### ~~A per-phase `veh` field on the phase struct~~ — DONE, `23c3470`

Optional `ph.veh`; absent or `[]` reproduces the old behaviour byte for byte.
Pinned in `tests/test_phaseRun.m` three ways: bit-identical to the closure form,
`veh = []` equivalent to an absent field, and a 400 kg phase vehicle diverging
from the 900 kg fallback so the field is provably *read*. The migration is still
open — below.

---

## Structural, deferred deliberately

### Migrate the entry scripts off their EOM closures onto `phase.veh`

Now that `ph.veh` exists, the chain entry scripts still bind a per-phase vehicle
inside an EOM closure and ignore the forwarded argument:

- `BM/run_ballistic.m:281-282`
- `HGV/run_boost_glide.m:365-366`
- `HGV/run_target.m:439-440`
- and `BM/run_ballistic_target.m`, which was written after `ph.veh` landed and
  still uses a closure

The closures are correct and remain supported. The question is whether the
scripts should read `phase(k).veh = vehK` instead, which says the same thing in
the interface rather than in a lambda.

**Deliberately not done in the same pass as the interface change.** The headline
results — `run_glide` 6986.82 km, `run_ballistic` 4536.36 km, `run_boost_glide`
7663.05 km, `run_target` 3811.240 km required, `run_ballistic_target` 4828.045 km
required (moved 2026-08-08, see below) — are pinned to their last printed digit,
and a migration touches the one construction those numbers are produced by. It should be its own change with its
own full-capture diff. A closure and `ph.veh` were *measured* bit-identical on
the staged test chain, so the migration is expected to move nothing; expected is
not verified.

### A state-scaled default `AbsTol` in `phaseRun`

`phaseRun` defaults to the **scalar** `AbsTol = 1e-10` for every state
component, across a vector whose entries are a radius in metres (6.4e6), angles
in radians (order 1), a speed in m/s (1e3–7e3) and a mass in kg (1e3–3e4). That
is 0.1 nm on the radius and 0.1 mg on the mass: not a physical error budget, and
only the tightest of the seven components does any real work.

`env.odeAbsTol` **already accepts a vector** — `odeset` takes one directly, and
`phaseRun` passes it through unexamined. That is documented in `phaseRun`'s
Notes with a worked per-component example. What is deferred is changing the
**default**, because every pinned number in all 21 test files was measured under
the scalar and a new default moves all of them at once. Do it as its own change:
pick the vector, re-measure the headline set, and re-pin deliberately rather
than as a side effect.

### Extract the duplicated `overrideOf` and `maxOver` helpers into `+util`

Re-measured 2026-08-08 by grep. Both counts went **up** with
`run_ballistic_target`, which is the argument for doing it: a fifth and a fourth
copy arrived by simply writing the next script.

| Helper | Copies | Where |
|---|---|---|
| `overrideOf` | **5** | `HGV/run_glide.m:327`, `BM/run_ballistic.m:1084`, `HGV/run_boost_glide.m:1183`, `HGV/run_target.m:1811`, `BM/run_ballistic_target.m:3448` |
| `maxOver` | **4** | `BM/run_ballistic.m:907`, `HGV/run_boost_glide.m:1047`, `HGV/run_target.m:1554`, `BM/run_ballistic_target.m:3167` |

`run_glide` is single-phase and has no peak-over-a-mask search, which is why it
carries `overrideOf` but not `maxOver`.

### Give the two `BM` scripts' `info` the re-integration payload

`BM/run_ballistic` and `BM/run_ballistic_target` both return `[traj,info]`, but
neither `info` carries the fields an independent checker needs to re-integrate
the same chain with a different solver — `phases`, `env`, `x0` and the per-phase
vehicle structs — which `HGV/run_boost_glide.m:909-913` and
`HGV/run_target.m:1236-1240` do carry. Verified 2026-08-08: grep for
`info.phases` in `BM/` returns nothing.

The consequence: `tests/test_fullChain.m` gives `run_boost_glide`'s two
junctions an independent `ode89` @ 1e-12 continuity check against the driver's
`ode45` @ 1e-10 (measured `1.340e-06 m` and `3.120e-07 m` on radius, against a
`1e-3 m` budget), and the `BM` chains' boost→coast and coast→descent junctions
get nothing equivalent. Add the five fields and the check follows almost for
free.

### Sweep degenerate inputs at every model boundary

**A systematic gap in this project's mutation discipline, and the reason the
mass-guard hole survived to be found by an outside reviewer.** Counting the
mutations enumerated in this library's commit messages gives roughly fifty
across the project, every one of them killed or fixed — and **every one
substituted a wrong FINITE value**: a flipped sign, a frozen denominator, a
transposed argument pair, a swapped trigonometric factor, a dropped unit, a
loop that returns after one step. Not one substituted `NaN`, `Inf`, zero or a
negative.

That is exactly the blind spot `docs/LESSONS_LEARNED.md` records under *A
comparison is not a validation* — found by external review across five files at
once, not by any mutation this project ran. A mutation study that only ever
perturbs *values* cannot find a guard that fails on *non-values*.

**Action:** a deliberate degenerate-input sweep at each model boundary — every
`+atmos`, `+grav`, `+aero`, `+eom`, `+guide`, `+prop` and `+util` entry point
fed `NaN`, `±Inf`, `0`, a negative and a non-scalar in each argument, with the
required outcome being a *named* refusal. `61e78ad` did this by hand for the
five files the review named; nothing makes it a standing property of the
library, and the next model added will not inherit it.

---

## What the pipeline review deliberately did not change

Recorded here so the decisions are not silently re-litigated, and so the
residual risk of each is written down. All five were adjudicated in `23c3470`
and `61e78ad`; none is a defect as the library is currently used.

- **The model clock is not re-based to zero** (review finding 4). The `t` handed
  to `eom`, `guide` and `terminate` is the phase's own `tspan` value, so a phase
  given `tspan = [10 50]` evaluates its guide at 10–50, not 0–40. Documented as
  a *TWO CLOCKS* block at `+coorbital/+prop/phaseRun.m:68-91`; not changed
  because every shipped schedule is written against it. **Residual risk:
  nothing enforces it.** `phaseRun` has no schema validation at all — verified
  2026-08-08, its only two `error` calls are on control and link *width*
  (`:277`, `:307`). A chain built with `tspan(1) ≠ 0` would silently shift every
  time-dependent command. An assertion, or a documented opt-in, would close it.
- **A missing terminal event is recorded, not enforced** (finding 5).
  `traj.phaseEnd` and `traj.endedOnEvent` say whether each phase stopped on its
  event or ran out of `tspan` (`phaseRun.m:263-267`), which the state history
  alone cannot distinguish. There is no `ph.requireEvent`. A time-limited phase
  is a legitimate design, so the flag has to be opt-in per phase — write it
  before a chain exists whose author assumes the event always fires.
- **A stateful guide is forbidden by documentation, not by construction**
  (finding 7). `phaseRun` states the exact evaluation count and requires a guide
  to be a pure function of `(t,x)`. Nothing checks it.
- **The Cartesian velocity equations were explicitly not adopted** (finding 8).
  `glide3DOF.m:129` and `boost3DOF.m:139` both floor at `V < 1`, so **launch
  from rest is unsupported** and the `1/V` coordinate singularity is guarded
  rather than removed. That is a different library, not a fix to a guard — but
  it is the thing that will have to change first if powered ascent from the pad
  is ever wanted.
- **The scalar `AbsTol` default is kept** (finding 6) — its own item above.

---

## Known limitations

### Rotating-Earth targeting needs an outer azimuth iteration

Both targeting scripts now refuse rather than mislead (see *Recently closed*),
but the capability is still missing. On a turning Earth the vehicle must be
aimed where the target is **going to be**, so the azimuth depends on the flight
time, the flight time on the ranging control, and the ranging control on the
azimuth: an outer iteration around the range solve, which neither script has.

`BM/run_ballistic_target` is the harder case and says so in its own refusal:
the loop would have to run **per branch**, because the lofted and depressed arcs
fly for very different times over the same geometry and therefore do not even
share an aim point.

### Cross-range is measured and warned about, not solved

Bisection matches a *distance*; nothing in it steers sideways. The shipped
zero-bank configuration comes out at zero cross-range as a property of *that*
configuration, not as a guarantee. `run_target` measures the cross-track offset
of the impact point from the launch-to-target great circle on every run and
raises `*** WARNING ***` when it exceeds the range tolerance.

A banked descent misses. At `run_boost_glide`'s 75° terminal bank the shipped
geometry converges on range to a **0.59 km** residual and lands **21.52 km**
away — pinned at `21524.695285497215 m` miss and `21515.222106821158 m`
cross-track in `tests/test_runTarget.m:359-360`. This is why `run_target` ships
`descBank = 0` where `run_boost_glide` ships 75: a targeting decision, not an
aerodynamic one. The warning block now forks on the bank and, at zero bank,
attributes the cross-track to the guidance schedules rather than to a
non-existent banked segment (`HGV/run_target.m:946-965`).

### `run_ballistic_target`'s minimum-energy mode throws away propellant

Below `cutFrac = 1` the burn is cut short and the **whole booster is jettisoned
with its unburned propellant**, which is why `separation = false` is refused for
this mode. That is the honest consequence of using thrust termination as the
second control, not a bug — but it means the mode's `ΔV` bookkeeping is not
comparable with the full-burn branches, and a staged booster would change the
answer.

Note also that the flown minimum-energy arc sits **+1.36 %** above the classical
apogee on the shipped 4828.045 km case. That is **physics, not solver error**,
and since 2026-08-08 it is **not a residual either** — nothing is driven against
it. The classical result assumes an impulsive burn at the impact radius in a
vacuum; this one finishes 82.18 km up, and a Keplerian arc from *that* burnout
state apogees at 965.641 km against the flown 965.613 km — agreement to
0.027 km. The script's summary attributes it.

**The cut is now small: 0.995 of the full burn, 150.7 kg left aboard.** At the
vehicle's 6 deg clamp the shipped target sits at 93 % of maximum range, and the
minimum-energy arc there is close to the max-range arc, so there is little
energy to give back. That is a property of the demonstration geometry, not of
the method — the 3175 km case in `test_runBallisticTarget` part 7 cuts at 0.961
and leaves 1170.6 kg.

### The `hHandoff` phugoid-trough warning is one-sided

The glide descends in a damped skip phugoid, so a handoff altitude placed above
a trough terminates the glide a whole skip early while **every phase still
reports nominal**. `run_boost_glide` detects this from a counterfactual
continued glide: if the handoff caught a trough, the vehicle climbs back out of
it, and the rebound is the evidence.

That only sees troughs the continued glide climbs out of. It cannot see a
handoff that truncated a shallow skip without rebounding. **`hHandoff = 25` km
gives 7490.28 km against the shipped 15 km case's 7663.05 km — 172.77 km of
range lost, `info.troughWarn = false`, no warning.** The blind spot is asserted
as such at `tests/test_fullChain.m:489-497` so it stays documented rather than
rediscovered.

### The booster structural coefficient is optimistic

`coorbital.util.boosterDefaults` gives `massDry = 1500 kg` on
`massProp = 30000 kg` (`:160-161`), so the structural coefficient is
`1500/31500 = ` **4.76 %**, against 10–15 % for real large solid stages. It is a
placeholder like everything else, but it is a load-bearing one: it is what makes
the boost milestone's premises line up with the glide milestone's, and a
realistic value would change every burnout state the later phases inherit.

### Event tolerances in `test_boostEvents` are far looser than the measured errors

All three asserted at `1e-6`; measured 2026-08-07 by re-running the test's own
propagations outside the harness:

| Assertion | Budget | Measured | Slack |
|---|---|---|---|
| Burnout mass (`tests/test_boostEvents.m:115`) | `1e-6 kg` | `3.365e-10 kg` | 3.5 orders |
| Burn duration (`:137`) | `1e-6 s` | `9.379e-13 s` | 6.0 orders |
| Apogee `gamma` (`:165`) | `1e-6 rad` | `3.699e-15 rad` | 8.4 orders |

Tightening them to, say, one order above the measured error would make the tests
sensitive to an event-solve regression they currently cannot see.

### A previously listed justification that turned out to be false

The old edition, and the shipped code it was quoting, said `run_ballistic_target`
refuses at a 6° `alphaMax` clamp because the clamp destroys the range hump so
the two branches do not exist. **`fd72b3d` measured that and it is wrong.**
Swept loft from −140° to 85° at a 6° clamp: the max-range angle is at
**−42.907° with 5211.525 km**, so the hump and both branches *do* exist. They
sit 2.907° below the then-shipped `loftMin = −40°`, so the cause is **bracket
width**, not the clamp.

**Superseded 2026-08-08 in the direction the measurement pointed.** The clamp is
now a VEHICLE property, `BM/vehicle_bm`'s `alphaMaxDeg = 6`, read by both BM
entry scripts; `loftMin` ships at −140° so the hump is inside the bracket; and
the demonstration target moved to 4828.045 km, inside the 6° depressed band of
roughly 4708–5212 km, so all three modes fly on the shipped configuration. What
12° would buy is depressed-branch reach down to 1684.117 km, and what it costs
is 156.224 km of maximum range — both pinned in `test_runBallisticTarget`
part 9, as a sensitivity study rather than as the shipped value. **The 6° is a
placeholder awaiting a qualification basis, not a cleared limit**, and that is
now stated in the vehicle file, in both scripts' user blocks and in the printed
limitations.

---

## Upstream, in pumpkyn — report, do not patch

`~/Desktop/proj7/external/pumpkyn` is third-party and **read-only**. Both items
re-verified 2026-08-08; both still present.

- **`src/+pumpkyn/+util/earth3D.m:262` calls `star3D(www)`, and no such function
  exists.** The real file is `src/+pumpkyn/+util/stars3D.m`, whose signature is
  `stars3D(fig, ax, backgroundMode, starImageFile, jd0, varargin)` — six
  arguments and a different contract, so this is **not** a rename typo and
  cannot be fixed by one. The call is guarded by `options.stars`, which defaults
  `false`, so nothing in this library trips it.
- **Of `earth3D`'s texture options, only `'day'` resolves to a file that is
  actually present.** The `+util` folder holds `earth-clouds-4k.jpg` and
  `earth-clouds-16k.jpg` and nothing else Earth-shaped; `'night'`, `'clouds'`,
  `'some clouds'`, `'BW'`, `'infrared'` and `'altEarth'` name
  `land_lights_16384.tif`, `cloud_combined_2048.tif`,
  `land_ocean_ice_cloud_2048.tif`, `outline-black-white-world-political.jpg`,
  `InfraredEarth.png` and `alternateEarth.bmp` (`earth3D.m:138-155`), none of
  which are there.

---

## Out of scope by design

Taken from the out-of-scope sections of the three plan briefs. Each becomes its
own plan when it is wanted; none is an oversight.

| Excluded | Note |
|---|---|
| **PEG and VOA closed-loop boost guidance** | Prescribed pitch only. `coorbital.guide.pitchProgram` already reads the state, so the signature is ready for a closed-loop law. The natural next milestone. |
| **Multi-stage boosters** | One stage, one burn. `phaseRun`'s per-phase `link` is the mechanism a second stage would use, and `ph.veh` now supplies the other half; nothing else blocks it. |
| **Terminal homing / guidance laws** | The descent phase flies a prescribed schedule, not a homing law. |
| **Aerothermal heating** | No Sutton–Graves anywhere. `noseRadius` is a carried placeholder in both vehicle files, read by no physics routine, and flagged as such at its point of definition. |
| **Fidelity increments** — J2 gravity, oblate geodetic altitude, tabulated aero, US76 atmosphere, rotating Earth on by default | Every one is a one-line handle swap by design; each needs its own validation test. `sphereGrav`'s always-zero `gLat` channel exists so J2 costs no signature change. A Mach-dependent aero model is now the highest-value one — see the Shuttle downrange discrepancy above. |
| **Phase 2 trajectory optimization** against `orbit_transfer/verify_common` | Prescribed-control simulation first, deliberately: the throughput requirement — multiple trajectories per second — is what shaped this architecture. |
| **Anything that modifies pumpkyn** | See above. |

---

## Documentation debt

All re-checked 2026-08-08.

- **`docs/DESIGN.md` now needs TWO more dated as-built sections, not one.** It
  carries §10 (2026-08-06, the glide propagator) and §11 (2026-08-07, boost /
  descent / full chain) and stops there. Missing: the **targeting and
  visualization** milestone (`greatCircleBearing`, `rangeSolve`, `+viz`,
  `HGV/run_target`) and the **ballistic targeting** milestone
  (`BM/run_ballistic_target`, the two branches, the two-parameter
  minimum-energy solve). The file's own front matter promises one per milestone.
- **Two applied reviews were never archived under `docs/reviews/`.** The folder
  holds the two plan-brief review pairs, the two math-note reviews and the
  pipeline review — but not the review of `8ea960d` (applied in `aba5736`) nor
  the review of `BM/run_ballistic_target` (applied in `fd72b3d`), both of which
  found criticals. Their findings survive only in commit messages. Archive them
  for the same reason the others were archived.
- **`README.md` (the front door) is stale in three places.** It says
  `tests/` holds "20 `test_*.m`" and "53 `.m` files in all" (`:123`, `:128`);
  there are now **21** test files and **55** `.m` files. And `:233-236` says the
  two LaTeX notes "are **forthcoming** and do not exist yet" — both
  `docs/hgv_dynamics_note.tex` and `docs/software_design.tex` were written and
  reviewed (`0e482c1`, `50417da`, `e224107`, `a8c3eda`) and their PDFs are in
  `docs/`. That paragraph should become a pointer.
- **A stale count in `docs/README.md:1456`:** the conventions section says the
  no-`%#ok` / no-`norm` / no-`for i =` sweep was re-verified "across all 39 `.m`
  files under `missiles/`". There are now **55**. Re-run 2026-08-08 over all 55
  and the sweep is still clean — zero `%#ok`, zero `for i =` or `for j =`, and
  the single `norm(` hit is inside a comment in `tests/test_viz.m:1705`
  *explaining* why `norm()` is not used. A count refresh, not a violation.
- **`docs/README.md:11` says "Three milestones have shipped" above a table with
  four rows.** Ballistic targeting was added to the table without updating the
  sentence.
- **`docs/README.md` and `+viz`'s scope status — checked, no action needed.**
  `docs/README.md` describes `+viz` as delivered throughout and removes it from
  the out-of-scope list. The two earlier plan briefs
  (`plan_2026-08-06_glide_propagator.md`, `plan_2026-08-07_boost_descent_chain.md`)
  do still list `+viz` as out of scope, which is correct — they are dated
  records of what was true when written and must not be edited.
- **There is no bibliography.** Three external references are now load-bearing
  (the Shuttle entry figures, Garwin's minimum-energy statement, the classical
  minimum-energy arc) and none of them is cited anywhere under `missiles/`.
  Add one when the two missing validation tests go in.
