# Missile Trajectory Library — Open Items

*Coorbital, Inc. — compiled 2026-08-07 against commit `8ea960d`, revised the
same day after the GPT-5.6-sol pipeline review
(`docs/reviews/pipeline_gpt56sol_review_2026-08-07.md`) was applied.*

Every claim below was checked against the code or measured by running it on
2026-08-07; the measurements and their sources are named inline. The suite was
green at the time of writing: **20 passed, 0 failed**, zero warnings.

---

## Not built yet

### `BM/run_ballistic_target.m` — the ballistic point-to-point script

Specified in full as **Task 5 of `docs/plan_2026-08-07_targeting_and_viz.md`**
(line 237 onward) and never implemented. `BM/` currently holds only
`run_ballistic.m` and `vehicle_bm.m`.

The framework is already there: `BM/run_ballistic` flies boost → coast → impact
and is asserted against the Keplerian closed form at `8.484e-13` relative
(`BM/run_ballistic.m:496`), and `coorbital.util.greatCircleBearing` plus
`coorbital.util.rangeSolve` are what turn "launch plus azimuth" into "launch
plus destination".

**The one thing that genuinely differs from the HGV case is the ranging
control.** `HGV/run_target` bisects on thrust-termination time because that
parameter is monotonic and single-valued. A pure ballistic trajectory is not: it
has **two** solutions for every range short of maximum — a lofted arc and a
depressed arc, either side of the minimum-energy solution — and a ballistic
user expects to *choose*. So the control here is the **loft angle** (the pitch
program's terminal attitude) and the user block carries a branch selector,
`'minimum-energy' | 'lofted' | 'depressed'`. Each branch must **report** which
side it actually landed on, measured from the flown apogee and flight time,
rather than trusting the bracket to have held it there; and both branch
solutions should be reported when they exist, with flight time, apogee, impact
speed and impact flight-path angle, so the trade is visible.

---

## Unreviewed

**The most recent commit — `8ea960d missiles: altitude inset, arc-frame camera,
adaptive exaggeration` — went in WITHOUT the subagent review that every other
change on this project received.** It is covered by the suite and was verified
by inspecting rendered frames, but it has had no second pair of eyes.

That matters here because review on this project has repeatedly found things the
tests did not. Two documented examples:

- **A mutation that survived because two quantities agreed to 9.3e-10 m.** At
  zero bank `run_target`'s measured impact-to-target miss and the magnitude of
  the solver's own range residual are indistinguishable, so two mutations —
  forcing `crossWarn = false`, and setting `missM = abs(resM)` — passed an
  earlier version of the suite. The banked case now in
  `tests/test_runTarget.m:323-356` separates them by a factor of 36.8 and pins
  both.
- **Four math errors in the Milestone-1 plan brief**, caught before any code was
  written: a sign error on the latitudinal-gravity projection, an omitted
  `cos(gamma)` in the equilibrium-glide balance, the wrong measured quantity for
  the Allen–Eggers comparison, and unguarded polar singularities. All four are
  in `docs/reviews/plan_gpt56terra_review_round1_2026-08-06.md`.

*(A third example is sometimes cited — a review that caught a test which could
not fail. No in-repo record of that one was found; the two above are the ones
the repository documents.)*

**Action:** put `8ea960d` through the same review the rest of the tree got.

---

## Structural, deferred deliberately

### ~~A per-phase `veh` field on the phase struct~~ — DONE 2026-08-07

`coorbital.prop.phaseRun` now accepts an **optional `ph.veh`**, the vehicle that
phase flies. Absent or `[]`, the chain-wide `veh` argument is used and the
behaviour is byte-for-byte what it was, so the change is purely additive and no
consumer needed editing. Pinned in `tests/test_phaseRun.m`: the staged chain
expressed with `ph.veh` is **bit-identical** to the same chain expressed with a
closure, `veh = []` behaves exactly as an absent field, and a phase carrying a
400 kg vehicle diverges from the 900 kg fallback (so the field is read, not
merely accepted).

`coorbital.eom.massConstant` keeps its mass guard: `ph.veh` is opt-in, so the
guard is still what catches a chain that did not use it.

**What is left is the migration**, listed under *Deferred* below.

### Migrate the entry scripts off their EOM closures onto `phase.veh`

Now that `ph.veh` exists, the four chain entry scripts still bind a per-phase
vehicle inside an EOM closure and ignore the forwarded argument:

- `BM/run_ballistic.m:281-282`
- `HGV/run_boost_glide.m:365-366`
- `HGV/run_target.m:439-440`

The closures are correct and remain supported. The question is whether the
scripts should read `phase(k).veh = vehK` instead, which says the same thing in
the interface rather than in a lambda.

**Deliberately not done in the same pass as the interface change.** The four
headline results — `run_glide` 6986.82 km, `run_ballistic` 4536.36 km,
`run_boost_glide` 7663.05 km, `run_target` 3811.240 km required — are pinned to
their last printed digit, and a migration touches the one construction those
numbers are produced by. It should be its own change with its own full-capture
diff. Note that a closure and `ph.veh` were *measured* bit-identical on the
staged test chain, so the migration is expected to move nothing; expected is
not verified.

### A state-scaled default `AbsTol` in `phaseRun`

`phaseRun` defaults to the **scalar** `AbsTol = 1e-10` for every state
component, across a vector whose entries are a radius in metres (6.4e6), angles
in radians (order 1), a speed in m/s (1e3–7e3) and a mass in kg (1e3–3e4). That
is 0.1 nm on the radius and 0.1 mg on the mass: not a physical error budget, and
only the tightest of the seven components does any real work.

`env.odeAbsTol` **already accepts a vector** — `odeset` takes one directly, and
`phaseRun` passes it through unexamined. That is now documented in `phaseRun`'s
Notes with a worked per-component example. What is deferred is changing the
**default**, because every pinned number in all 20 test files was measured under
the scalar and a new default moves all of them at once. Do it as its own change:
pick the vector, re-measure the headline set, and re-pin deliberately rather
than as a side effect.

### Extract the duplicated `overrideOf` and `maxOver` helpers into `+util`

Measured 2026-08-07 by grep:

| Helper | Copies | Where |
|---|---|---|
| `overrideOf` | **4** | `HGV/run_glide.m:327`, `BM/run_ballistic.m:1084`, `HGV/run_boost_glide.m:1183`, `HGV/run_target.m:1572` |
| `maxOver` | **3** | `BM/run_ballistic.m:907`, `HGV/run_boost_glide.m:1047`, `HGV/run_target.m:1401` |

`run_glide` is single-phase and has no peak-over-a-mask search, which is why it
carries `overrideOf` but not `maxOver`.

### Give `run_ballistic`'s `info` the re-integration payload

`BM/run_ballistic` **does** already return `[traj,info]`, but its `info` carries
only summary scalars. What it lacks are the fields an independent checker needs
to re-integrate the same chain with a different solver — `phases`, `env`, `x0`
and the per-phase vehicle structs — which `HGV/run_boost_glide.m:909-913` and
`HGV/run_target.m:1103-1107` do carry.

The consequence: `tests/test_fullChain.m` gives `run_boost_glide`'s two
junctions an independent `ode89` @ 1e-12 continuity check against the driver's
`ode45` @ 1e-10 (measured `1.340e-06 m` and `3.120e-07 m` on radius, against a
`1e-3 m` budget), and `run_ballistic`'s boost→coast and coast→descent junctions
get nothing equivalent. Add the five fields and the check follows almost for
free.

---

## Known limitations

### Rotating-Earth targeting needs an outer azimuth iteration

`HGV/run_target`'s launch azimuth is the closed-form great-circle initial
bearing, which is the whole answer only at `env.omegaE = 0`. Turn rotation on
and the ground beneath the target moves — the script measures **620 km** of
eastward sweep at 35°N over the shipped 1626 s flight. The azimuth then depends
on the flight time, the flight time on the cutoff, and the cutoff on the
azimuth: an outer iteration around the range solve, which the script does not
have. It states this in limitation 1 of every summary and prints a caution when
`earthSpin` is set true.

### Cross-range is measured and warned about, not solved

Bisection matches a *distance*; nothing in it steers sideways. The shipped
zero-bank configuration comes out at zero cross-range as a property of *that*
configuration, not as a guarantee. `run_target` measures the cross-track offset
of the impact point from the launch-to-target great circle on every run and
raises `*** WARNING ***` when it exceeds the range tolerance.

A banked descent misses. At `run_boost_glide`'s 75° terminal bank the shipped
geometry converges on range to a **0.59 km** residual and lands **21.52 km**
away — pinned at `21524.695285497215 m` miss and `21515.222106821158 m`
cross-track in `tests/test_runTarget.m:355-356`. This is why `run_target` ships
`descBank = 0` where `run_boost_glide` ships 75: a targeting decision, not an
aerodynamic one.

### The `hHandoff` phugoid-trough warning is one-sided

The glide descends in a damped skip phugoid, so a handoff altitude placed above
a trough terminates the glide a whole skip early while **every phase still
reports nominal**. `run_boost_glide` detects this from a counterfactual
continued glide: if the handoff caught a trough, the vehicle climbs back out of
it, and the rebound is the evidence.

That only sees troughs the continued glide climbs out of. It cannot see a
handoff that truncated a shallow skip without rebounding. **Measured 2026-08-07:
`hHandoff = 25` km gives 7490.28 km against the shipped 15 km case's 7663.05 km
— 172.77 km of range lost, `info.troughWarn = false`, no warning.** The blind
spot is asserted as such at `tests/test_fullChain.m:489-497` so it stays
documented rather than rediscovered.

### The booster structural coefficient is optimistic

`coorbital.util.boosterDefaults` gives `massDry = 1500 kg` on
`massProp = 30000 kg`, so the structural coefficient is
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

---

## Upstream, in pumpkyn — report, do not patch

`~/Desktop/proj7/external/pumpkyn` is third-party and **read-only**. These are
recorded so the next person does not rediscover them; the fix belongs upstream.

- **`src/+pumpkyn/+util/earth3D.m:262` calls `star3D(www)`, and no such function
  exists.** The real file is `src/+pumpkyn/+util/stars3D.m`, whose signature is
  `stars3D(fig, ax, backgroundMode, starImageFile, jd0, varargin)` — six
  arguments and a different contract, so this is **not** a rename typo and
  cannot be fixed by one. The call is guarded by `options.stars`, which defaults
  `false` (`earth3D.m:73` and `:86-87`), so nothing in this library trips it.
- **Of `earth3D`'s texture options, only `'day'` resolves to a file that is
  actually present.** The `+util` folder contains `earth-clouds-4k.jpg` and
  `earth-clouds-16k.jpg` and nothing else Earth-shaped; `'night'`, `'clouds'`,
  `'some clouds'`, `'BW'`, `'infrared'` and `'altEarth'` name
  `land_lights_16384.tif`, `cloud_combined_2048.tif`,
  `land_ocean_ice_cloud_2048.tif`, `outline-black-white-world-political.jpg`,
  `InfraredEarth.png` and `alternateEarth.bmp` respectively, none of which are
  there.

---

## Out of scope by design

Taken from the out-of-scope sections of the three plan briefs. Each becomes its
own plan when it is wanted; none is an oversight.

| Excluded | Note |
|---|---|
| **PEG and VOA closed-loop boost guidance** | Prescribed pitch only. `coorbital.guide.pitchProgram` already reads the state, so the signature is ready for a closed-loop law. The natural next milestone. |
| **Multi-stage boosters** | One stage, one burn. `phaseRun`'s per-phase `link` is the mechanism a second stage would use; nothing else blocks it. |
| **Terminal homing / guidance laws** | The descent phase flies a prescribed schedule, not a homing law. |
| **Aerothermal heating** | No Sutton–Graves anywhere. `noseRadius` is a carried placeholder in both vehicle files, read by no physics routine, and flagged as such at its point of definition. |
| **Fidelity increments** — J2 gravity, oblate geodetic altitude, tabulated aero, US76 atmosphere, rotating Earth on by default | Every one is a one-line handle swap by design; each needs its own validation test. `sphereGrav`'s always-zero `gLat` channel exists so J2 costs no signature change. |
| **Phase 2 trajectory optimization** against `orbit_transfer/verify_common` | Prescribed-control simulation first, deliberately: the throughput requirement — multiple trajectories per second — is what shaped this architecture. |
| **Anything that modifies pumpkyn** | See above. |

---

## Documentation debt

- **`docs/DESIGN.md` needs a third dated as-built section** for the targeting
  and visualization milestone. It currently carries §10 (2026-08-06, the glide
  propagator) and §11 (2026-08-07, boost / descent / full chain) and stops
  there, so `greatCircleBearing`, `rangeSolve`, the `+viz` package and
  `HGV/run_target` have no as-built record. The file's own front matter promises
  one per milestone.
- **`docs/README.md` and `+viz`'s scope status — checked 2026-08-07, no action
  needed.** `docs/README.md` describes `+viz` as delivered throughout (lines 17,
  235, 285) and explicitly removes it from the out-of-scope list at line 1144;
  its out-of-scope table no longer mentions it. The two earlier plan briefs
  (`plan_2026-08-06_glide_propagator.md`, `plan_2026-08-07_boost_descent_chain.md`)
  do still list `+viz` as out of scope, which is correct — they are dated
  records of what was true when written and must not be edited.
- **The two LaTeX notes are forthcoming**, not missing: `docs/hgv_dynamics_note.tex`
  for the mathematics and `docs/software_design.tex` for the software design.
  Neither file exists yet. They were deferred by design until the interfaces had
  survived contact with working code, which they now have.
- **A stale count in `docs/README.md`:** the conventions section says the
  no-`%#ok` / no-`norm` / no-`for i =` sweep was re-verified "across all 39 `.m`
  files under `missiles/`". There are now **53**. Re-run 2026-08-07 over all 53
  and the sweep is still clean — zero `%#ok`, zero `norm(`, zero `for i =` or
  `for j =` — so this is a count refresh, not a violation.
