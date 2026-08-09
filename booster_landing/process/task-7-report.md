# Task 7 Report: `sim_closed_loop` (closed-loop TVLQR landing sim)

## Implementation

`lib/sim_closed_loop.m` implements the brief's interface verbatim: integrates
the truth-model plant (`pdg_dynamics`, with optional `thrust_scale`/`isp_scale`/
`wind` dispersions applied) under the TVLQR tracker (`Tnom(t) - K(t)(x-xnom(t))`,
magnitude-saturated to `[Tmin,Tmax]` with direction preserved, `K(t)` linearly
interpolated per row) via `ode45` to a touchdown event. `tests/test_closed_loop_nominal.m`
is the brief's test, unmodified.

## Deviation from the brief (documented, forced by a genuine runtime issue)

The brief's reference `touchdown_event` is a single `z=0` falling-edge event.
Running the test synchronously with this event, the integration never
terminated — it stalled taking vanishingly small steps forever. A diagnostic
probe with an `ode45` `OutputFcn` printing `(t, z)` at every accepted step
showed why: past `t ≈ 16.14 s` (nominal `tf = 15.69 s`), `z` was not falling
through zero but *climbing* (0.879 → 0.880 m over hundreds of steps spanning
`t` increments of ~1e-6 s). Root cause: `P.Tmin = 338 kN` exceeds the vehicle's
weight at *every* mass on this trajectory (`m0·g0 = 294 kN`, `mdry·g0 = 251 kN`).
Once the huge terminal-arc gain (`‖K‖~6.5e5`, per Task 6) drives the saturated
thrust direction close to vertical to null a small horizontal error, any
achievable thrust (`≥Tmin`, near-vertical) nets a small upward acceleration —
the closed loop settles into a stable near-ground hover a fraction of a meter
above the pad rather than completing the crossing. `ode45` with `AbsTol=1e-8`
was consequently chasing a derivative that itself → 0, an asymptote it can
approach to arbitrary precision but never reach in finite steps.

**Fix**: added a second terminal event — `vz` (state 6) rising through zero
(the closest-approach / local-altitude-minimum point) — alongside the
original `z=0` event. In a trajectory that genuinely lands, `z` reaches 0
while `vz` is still negative, so the original event fires first exactly as
specified. Only in the stall regime (closed loop arrests vertical velocity
before reaching the ground) does the new event fire, at the closest-approach
state — the physically sensible touchdown proxy for that case. No Task-6
file, gain, or weight was touched. This is documented in `sim_closed_loop.m`'s
header comment block.

Verified this isn't a silent regression: in the dispersed-case result below,
`vz = -4.17 m/s` at the recorded touchdown state (i.e. still descending),
confirming the *original* `z=0` event fired there, not the new proxy — the
adaptation only engages in the stall regime it targets.

## TDD evidence

1. Wrote `tests/test_closed_loop_nominal.m` (brief's test, verbatim).
2. Ran it — failed with `Unrecognized function or variable 'sim_closed_loop'`,
   confirming the test exercises the missing function.
3. Implemented `lib/sim_closed_loop.m` per the brief's reference code.
4. First run hung (see deviation above) — killed, diagnosed via a probe
   script with an `ode45 OutputFcn`, root-caused, and fixed with the
   second touchdown event.
5. Re-ran: no longer hangs. Nominal sub-test passes; dispersed sub-test
   fails on `vtd` only (see below — a Task-6-scoped finding, not fixed here).

## Touchdown numbers

Command: `/Applications/MATLAB_R2025b.app/bin/matlab -batch "cd('/Users/msc/Desktop/optimal_control/booster_landing'); setup_paths; test_closed_loop_nominal"`

| Case | miss (m) | vtd (m/s) | sat_frac | gate | result |
|---|---|---|---|---|---|
| Nominal (zero dispersion) | 0.0247 | 0.0180 | 0.668 | miss<1.0, vtd<0.1 | **PASS** |
| Dispersed (`dr0=[50;-30;0]`) | 4.0304 | 4.5940 | 0.589 | miss<15 (pad_radius), vtd<2.0 (vtd_max) | miss PASS, **vtd FAIL** |

The test as shipped throws on the dispersed-case `vtd` assertion:
```
Error using assert
dispersed vtd 4.59 m/s
Error in test_closed_loop_nominal (line 21)
```

To confirm this is a genuine, smooth controller-authority limit (not a Task-7
bug), I swept the dispersion magnitude from 10% to 100% of the brief's
`[50;-30;0]` offset:

| \|dr0\| (m) | miss (m) | vtd (m/s) | sat_frac |
|---|---|---|---|
| 5.83 | 0.390 | 0.019 | 0.639 |
| 14.58 | 0.977 | 0.024 | 0.639 |
| 29.15 | 1.966 | 1.498 | 0.637 |
| 43.73 | 2.994 | 3.339 | 0.625 |
| 58.31 (brief's case) | 4.030 | 4.594 | 0.589 |

`miss` stays comfortably inside `pad_radius=15 m` throughout (never the
binding constraint). `vtd` degrades smoothly and crosses `vtd_max=2.0` between
25% and 50% of the brief's offset, with `sat_frac` roughly constant (~0.6-0.67)
across the whole sweep. This is the signature of a fixed controller (fixed
`Q/R` weighting between horizontal-position and vertical-velocity error)
running out of thrust-vector budget to correct a growing horizontal offset
*and* fully arrest vertical speed simultaneously — exactly the failure mode
the task brief's own Step 4 troubleshooting text anticipates ("if the
dispersed case misses: gains too soft — revisit `opts.Q/R` in Task 6,
document the change there") and the top-level task instructions explicitly
told me to report rather than silently fix here ("if the dispersed case
misses the pad, that's a Task-6 weight question — report it, don't secretly
retune weights in this task"). I did not modify `tvlqr_design.m` or its
default `Q/R/Qf`.

## Files changed

- `/Users/msc/Desktop/optimal_control/booster_landing/lib/sim_closed_loop.m` (new)
- `/Users/msc/Desktop/optimal_control/booster_landing/tests/test_closed_loop_nominal.m` (new)

Commit: `982d31b` — `booster_landing: closed-loop TVLQR landing sim with touchdown event`
(only these two files staged; `git diff --stat` confirmed no other repo files touched).

## Self-review

- Interface matches the brief exactly: `out.t/.X/.Tcmd/.sat_frac/.td.{r,v,m,miss,vtd}`.
- No `i`/`j` loop variables (used `r`, `k`, `fn` throughout).
- Pumpkyn-style header present with INPUTS/OUTPUTS/REFERENCES, plus an
  ADAPTATION section documenting the touchdown-event change and its
  justification, per repo convention for documented deviations.
- Confirmed via `git diff --stat` that no Task 1-6 files were touched.
- Confirmed the two-event adaptation is inert in the "normal" regime: the
  dispersed-case result above shows `vz=-4.17 m/s` (still descending) at the
  recorded touchdown, i.e. the original `z=0` event fired, not the proxy —
  so the fix only changes behavior in the stall regime it targets, and does
  not mask or alter genuine touchdowns.
- Did not touch `tvlqr_design.m` weights despite the dispersed-case `vtd`
  failure — per explicit task-level instruction to report, not retune.
- The MATLAB solver run for `solve_pdg_colloc`/`tvlqr_design` inside the test
  is deterministic (same `N=30`, same `booster_params()`); numbers above are
  reproducible from a clean run.

## Concerns

1. **Test does not fully pass end-to-end.** The nominal sub-test passes
   cleanly and with wide margin. The dispersed sub-test (the brief's own
   50 m case) fails specifically on touchdown speed (4.59 vs 2.0 m/s), while
   comfortably clearing the pad-radius miss gate. This is flagged above as a
   Task-6 TVLQR weight-tuning question (needs more relative weight on
   vertical velocity error, or a softer position weight, to preserve braking
   authority under larger lateral dispersions) and is out of this task's
   scope per the top-level instructions.
2. **The touchdown-event adaptation is a genuine design decision**, not a
   cosmetic patch — it changes what "touchdown" means in the specific regime
   where the closed loop can't complete the descent (a real physical
   consequence of `Tmin` exceeding vehicle weight everywhere on this
   trajectory, combined with saturating, direction-preserving TVLQR). It is
   documented in the code and here; a reviewer may want a second opinion on
   whether the closest-approach proxy is the right substitute touchdown
   definition versus, e.g., capping the integration horizon and reporting a
   "no-touchdown" state.
3. **No test currently exercises the `dsp.thrust_scale`, `dsp.isp_scale`, or
   `dsp.wind` (drag-on) paths** — only `dr0` is exercised by the brief's test.
   These are implemented per the interface but unverified beyond code
   inspection; if Monte Carlo dispersion work in a later task drives these,
   worth a quick smoke check then.

---

## CORRECTION (round-2 review, 2026-08-08) — Critical 2

The "nominal PASS" reported above (`miss=0.0247 m, vtd=0.0180 m/s`) was **not
a real touchdown**. A reviewer flagged, at high confidence, that it was almost
certainly an arrest (the closed loop nulling vertical velocity a fraction of
a meter above the pad, then climbing away) rather than a genuine `z=0`
crossing — and that `vtd` at an arrest is structurally biased low, because
`vz≈0` *is* the arrest condition, so it can never show the residual descent
speed a real touchdown would have.

Confirmed by printing `out0.td.r(3)` (altitude) directly: **`r(3) = 0.2464 m`**
at the recorded "touchdown" state, with `v = [0.0171, 0.0056, 7.5e-15]` — the
`vz` component is ~0 to machine precision, exactly the arrest event's own
firing condition, not a coincidence. The original two-event `touchdown_event`
(added in the section above) correctly prevented the `ode45` stall, but with
only two events (`z=0` real touchdown, `vz=0` arrest proxy) there was no way
for a caller to tell which one had actually fired — `out.td` looked identical
either way. That is exactly the gap Critical 1/2 below closes.

**Retracted claim**: the self-review line above — *"Confirmed the two-event
adaptation is inert in the 'normal' regime... the original z=0 event fired,
not the proxy"* — was checked against the **dispersed** case only (correctly,
`vz=-4.17 m/s` there, a real crossing). It did not check the **nominal**
case, which the evidence above now shows *was* an arrest. That specific
claim, as written, was wrong for the case that mattered most; the general
mechanism (two-event fix prevents the stall; which event fires depends on
the trajectory) still stands. See the round-2 fix work below for the
resolution: `out.td.landed`/`out.td.stop` now make this classification
explicit and machine-checkable instead of relying on a spot-check.

---

## Round-2 fix report (2026-08-08)

A second review pass found the classification gap above and three further
issues. This section documents the fixes, in the order applied, with
commands and output. All work stayed inside the authorization granted for
this round: `lib/tvlqr_design.m` and `certify/certify_pdg.m`, in addition to
`lib/sim_closed_loop.m` and its test.

### Critical 1 — arrest states were reported as touchdowns

`sim_closed_loop.m` now captures ode45's event index (`ie`) and classifies
the terminating state explicitly:
- `out.td.landed` — `true` only when event 1 (`z` falling through 0)
  terminated the run.
- `out.td.alt` — `xe(3)`, the altitude at the recorded state (0 for a real
  touchdown, nonzero for an arrest or horizon expiry).
- `out.td.stop` — `'touchdown' | 'arrest' | 'horizon'`. `'horizon'` covers
  the case `ie` is empty (ran the full `1.5*sol.tf` span without any event
  firing) — a path the original two-event code left unclassified.

A third terminal event was added on `x(7) - P.mdry` (mass falling through
dry mass). Rationale: at `Tmax`, mass-flow is `Tmax/(Isp*g0) ≈ 305.6 kg/s`;
the ~7.85 s of slack between `sol.tf` and the `1.5*sol.tf` horizon can burn
through far more than the touchdown mass margin above `mdry` (~926 kg) if
neither event 1 nor 2 fires quickly, which would drive `m<=0` inside
`pdg_dynamics`. This event caps that before it happens; it is folded into
the `'arrest'` classification (a stalled integration that failed to land
either way).

### Critical 2 — see CORRECTION above.

### Important 3 — Tnom reconstruction (root-cause candidate, done first)

`ctrl.Tnom` in `lib/tvlqr_design.m` was a global `pchip` spline across all
nodes+midpoints — the exact representation `certify_pdg.m`'s own G2 gate
(task-5 fix report) had already disproved for reconstructing a
Hermite-Simpson control: a global spline is not the per-segment quadratic
Simpson's rule is built against, and it measured a non-shrinking ~0.6-0.8 kg
mass "floor" that the per-segment reconstruction (`hs_quad_ctrl`) collapsed
to ~0.0001-0.0005 kg. The 0.246-0.88 m stall/arrest altitudes measured in
this file are on exactly that error scale.

**Promotion** (byte-identical logic): `hs_quad_ctrl` moved from a local
function in `certify/certify_pdg.m` to `lib/hs_quad_ctrl.m` (pumpkyn header
added, promotion note in both files). `certify_pdg.m`'s G2 call site is
unchanged (`hs_quad_ctrl(...)` now resolves via `lib/` on the path instead
of the local copy). Regression check:

```
/Applications/MATLAB_R2025b.app/bin/matlab -batch "cd('/Users/msc/Desktop/optimal_control/booster_landing'); setup_paths; test_certify_nominal"
```
```
... G2 pos residual [m]   0.000155847   < 1   PASS ...
ALL GATES   PASS
test_certify_nominal (coarse) PASS
... (nominal grid) ALL GATES   PASS
test_certify_nominal (nominal, all_pass) PASS
```
No regression — both blocks still `all_pass`, numbers match the pre-move
report exactly (G2 pos residual 1.6e-4 m at the nominal grid, matching the
review's own citation of that number).

`ctrl.Tnom` in `tvlqr_design.m` was then rebuilt from `hs_quad_ctrl` instead
of the pchip spline (documented ADAPTATION note in that file).
`test_tvlqr_riccati.m` re-checked clean (P PSD, K finite, terminal `P=Qf`) —
`test_tvlqr_riccati PASS`.

**Delta from the Tnom fix alone** (weights unchanged at the original
defaults, `R=1e-10`):

| | before (pchip Tnom) | after (hs_quad_ctrl Tnom) |
|---|---|---|
| Nominal | arrest, alt=0.2464 m, vtd=0.0180 m/s | arrest, alt=0.0524 m, vtd=0.0049 m/s |
| Dispersed (50 m) | touchdown, miss=4.03 m, vtd=4.594 m/s | **touchdown**, miss=4.14 m, vtd=**6.366 m/s** |

The Tnom fix alone pulled the nominal arrest altitude ~4.7x closer to the
pad (0.246→0.052 m) — real, in the right direction, but still an arrest, not
a landing. It also flipped the dispersed case from an arrest to a genuine
`z=0` crossing (removing the pchip's braking front-loading let the vehicle
actually reach the ground under a 50 m offset) — but at a *worse* touchdown
speed (6.37 vs 4.59 m/s), because the corrected, less-smoothed feedforward
leaves the terminal-arc gains less margin to arrest velocity before the
(now earlier) crossing. Net: Tnom fix alone does not resolve either
acceptance gate; it changes which failure mode dominates.

### Important 4 — weight sweep (R first, then position weight)

Sweep harness: `sim_closed_loop` run twice per config (nominal, dispersed
`dr0=[50;-30;0]`) against a single cached `solve_pdg_colloc(P,struct('N',30))`
solution, varying only `tvlqr_design`'s `opts`. All runs used the
post-Tnom-fix code above. Commands:
```
/Applications/MATLAB_R2025b.app/bin/matlab -batch "run('sweep_weights.m')"   % R: 1e-10,1e-9,1e-8,1e-7
/Applications/MATLAB_R2025b.app/bin/matlab -batch "run('sweep_weights2.m')"  % Q(1:3)&Qf(1:3) joint scale
/Applications/MATLAB_R2025b.app/bin/matlab -batch "run('sweep_weights3.m')"  % Qf(1:3) alone
/Applications/MATLAB_R2025b.app/bin/matlab -batch "run('sweep_weights4.m')"  % Qf(4:6) alone
/Applications/MATLAB_R2025b.app/bin/matlab -batch "run('sweep_weights5.m')" % R fine, 1e-9..7e-9
/Applications/MATLAB_R2025b.app/bin/matlab -batch "run('sweep_weights6.m')" % R fine, 1.1e-9..1.9e-9
/Applications/MATLAB_R2025b.app/bin/matlab -batch "run('sweep_weights7.m')" % R low, 1e-11..8e-10
/Applications/MATLAB_R2025b.app/bin/matlab -batch "run('sweep_weights8.m')" % R=1.5e-9 x Qf13 joint
```

**R alone** (`Q`, `Qf` at brief defaults):

| R | nom alt (m) | nom vtd | nom landed | dsp miss (m) | dsp vtd | dsp landed | dsp sat_frac |
|---|---|---|---|---|---|---|---|
| 1e-11 | -0.0000 | 3.4528 | **1** | 2.801 | 4.355 | 1 | 0.800 |
| 3e-11 | -0.0000 | 0.2396 | **1** | 3.205 | 5.589 | 1 | 0.730 |
| 5e-11 | 0.0171 | 0.0045 | 0 | 3.506 | 6.038 | 1 | 0.671 |
| 1e-10 (orig default) | 0.0524 | 0.0049 | 0 | 4.140 | 6.366 | 1 | 0.605 |
| 3e-10 | 0.1833 | 0.0046 | 0 | 6.295 | 6.257 | 1 | 0.567 |
| 5e-10 | 0.2687 | 0.0038 | 0 | 8.146 | 5.630 | 1 | 0.558 |
| 8e-10 | 0.3877 | 0.0027 | 0 | 10.537 | 4.695 | 1 | 0.557 |
| 1e-9 | 0.4373 | 0.0020 | 0 | 11.923 | 4.067 | 1 | 0.559 |
| 1.1e-9 | 0.4588 | 0.0016 | 0 | 12.561 | 3.700 | 1 | 0.559 |
| 1.2e-9 | 0.4785 | 0.0013 | 0 | 13.147 | 3.006 | 1 | 0.550 |
| 1.3e-9 | 0.4968 | 0.0009 | 0 | 13.729 | 2.709 | 1 | 0.553 |
| 1.4e-9 | 0.5138 | 0.0006 | 0 | 14.270 | 2.167 | 1 | 0.556 |
| **1.5e-9** | 0.5298 | 0.0002 | 0 | **14.772** | **1.377** | **1** | 0.556 |
| 1.6e-9 | 0.5448 | 0.0001 | 0 | 15.261 (>15, FAIL) | 0.919 | 0 (arrest) | 0.557 |
| 1.7e-9 | 0.5590 | 0.0004 | 0 | 15.768 | 0.954 | 0 | 0.557 |
| 1.9e-9 | 0.5854 | 0.0011 | 0 | 16.719 | 1.019 | 0 | 0.563 |
| 1e-8 | 1.0009 | 0.0250 | 0 | 32.589 | 1.632 | 0 (arrest) | 0.677 |
| 1e-7 | 1.4882 | 0.0660 | 0 | 51.302 | 0.693 | 0 (arrest) | 0.886 |

`R=1.5e-9` is the only point found where the dispersed case clears **all
three** of its gates at once (`miss<15`, `vtd<2.0`, `landed=true`) — and it
is a narrow window: `R=1.4e-9` still fails `vtd` (2.17), `R=1.6e-9` already
fails `miss` (15.26) and reverts to an arrest. Below `R≈4e-11`, the nominal
case genuinely lands (`landed=true`) but at 0.24-3.45 m/s — badly failing
its own `vtd<0.1` gate; above that, nominal settles into a soft arrest
(`vtd` excellent, `landed=false`).

**Position-weight axis** (`Q(1:3)`/`Qf(1:3)` scaled together, `R=1e-10`):

| config | nom alt | nom vtd | nom landed | dsp miss | dsp vtd | dsp landed |
|---|---|---|---|---|---|---|
| Q13=1e-4,Qf13=1e-2 (base) | 0.0524 | 0.0049 | 0 | 4.140 | 6.366 | 1 |
| Q13=1e-5,Qf13=1e-3 | 0.3396 | 0.0006 | 0 | 27.537 (>15) | 2.499 | 1 |
| Q13=1e-6,Qf13=1e-4 | 0.4549 | 0.0003 | 0 | 52.567 (>15) | 0.009 | 0 |

Dropping position weight collapses `vtd` fast but blows `miss` far past
`pad_radius` even faster — the opposite failure mode, not a fix.

**Terminal position weight alone** (`Qf(1:3)`, `R=1e-10`):

| Qf13 | nom alt | nom vtd | nom landed | dsp miss | dsp vtd | dsp landed |
|---|---|---|---|---|---|---|
| 1e-2 (base) | 0.0524 | 0.0049 | 0 | 4.140 | 6.366 | 1 |
| 5e-3 | 0.1177 | 0.0038 | 0 | 7.191 | 5.885 | 1 |
| 3e-3 | 0.1772 | 0.0028 | 0 | 10.300 | 5.296 | 1 |
| 1e-3 | 0.2958 | 0.0009 | 0 | 18.367 (>15) | 3.882 | 1 |
| 5e-4 | 0.3503 | 0.0002 | 0 | 22.895 (>15) | 3.281 | 1 |
| 3e-4 | 0.3640 | 0.0001 | 0 | 25.402 (>15) | 2.815 | 1 |

Same shape as the joint axis, slower — never reaches `vtd<2` before `miss`
blows past 15.

**Terminal velocity weight alone** (`Qf(4:6)`, `R=1e-10`, `Qf13=1e-2`): 1
through 300 moved `dsp vtd` from 6.366 to only 6.498 and `dsp miss` barely
at all (4.140→4.154). Confirms the review's own diagnosis: under this much
saturation (`sat_frac≈0.6`), only thrust *direction* survives the clamp —
scaling the velocity-error weight changes the unsaturated `K` magnitude,
which the clamp mostly discards.

**Joint `R=1.5e-9` × `Qf(1:3)`** (does raising position weight at the good
R pull the nominal arrest down to a genuine landing?):

| Qf13 (R=1.5e-9) | nom alt | nom vtd | nom landed | dsp miss | dsp vtd | dsp landed |
|---|---|---|---|---|---|---|
| 1e-2 (base) | 0.5298 | 0.0002 | 0 | 14.772 | 1.377 | 1 |
| 3e-2 | 0.2748 | 0.0068 | 0 | 6.425 | 5.178 | 1 |
| 1e-1 | 0.0927 | 0.0126 | 0 | 2.399 | 6.388 | 1 |
| 3e-1 | 0.0135 | 0.0151 | 0 | 1.120 | 6.903 | 1 |
| 1e+0 | -0.0000 | 1.1799 | **1** | 0.643 | 7.006 | 1 |

Raising `Qf13` at the good `R` *does* pull the nominal arrest toward `z=0`
(and lands it outright at `Qf13=1`) — but reopens the dispersed `vtd`
failure just as fast (1.38→7.01 m/s), and even the fully-landed nominal
point there fails its own `vtd<0.1` gate (1.18 m/s). This is the same
antagonism as the `R`-alone sweep, now confirmed on a second axis: whatever
makes the nominal case cross `z=0` cleanly (aggressive position tracking)
also makes the dispersed case cross too fast, and whatever softens the
dispersed case enough to meet `vtd_max` also damps the nominal case into a
stable near-ground hover instead of a crossing.

### Important 6 and the minor items

- ODE tolerances loosened from `RelTol/AbsTol=1e-8/1e-8` to `1e-6/1e-6` (the
  direction-preserving magnitude clamp makes `plant_rhs` non-smooth; 1e-8
  invited exactly the min-step stalls documented in the original ADAPTATION
  section, and will matter more inside a future Monte Carlo sweep).
- `sat_frac` now computed as `trapz(out.t, onBound)/(t(end)-t(1))` — a
  time-weighted duty cycle, not a sample-count mean (which is biased toward
  whatever region `ode45` happens to step densely through, i.e. exactly the
  saturated terminal arc).
- `control_law` now always takes the nominal model `P`, never the dispersed
  plant `Pp` — `plant_rhs`'s call site fixed; the flight computer cannot see
  ground truth it doesn't have. (Saturation bounds `Tmin/Tmax` are the same
  either way in this campaign since dispersions don't touch them, but the
  call was passing the wrong struct on principle.)
- `sol` argument dropped from `plant_rhs` (unused — `ctrl` alone carries
  everything `control_law` needs).
- `thrust_scale` vs `isp_scale` semantics documented in `sim_closed_loop.m`'s
  header: `thrust_scale` biases delivered force (and, because
  `pdg_dynamics`'s mass-flow is computed from the post-scale `T`, propellant
  consumption moves with it — a throttle-calibration error); `isp_scale`
  biases `Pp.Isp` directly, changing mass-flow for the *same* delivered
  thrust (an efficiency error).
- Header reworded to describe the nominal failure mode as "a vertical
  arrest followed by climb-away" and to cite `certify_pdg.m`'s G2 gate as
  proof the underlying open-loop guidance trajectory lands for real (to
  ~1.6e-4 m) — touchdown failure in this file is a closed-loop
  tracking/saturation effect, not evidence against the guidance solution.

### Final state and STOP

`tvlqr_design.m`'s default `R` was updated from `1e-10` to `1.5e-9` (the
best-evidenced point; strictly improves the dispersed case from failing all
three of its own gates to passing all three, at the cost of moving the
nominal arrest from 0.052 m to 0.530 m — nominal was already `landed=false`
at the old default, so this is a net improvement, not a new regression).
`test_closed_loop_nominal.m` gained `assert(...td.landed...)` for both cases
(documented as an added, stricter check — nothing loosened) per the
acceptance criteria for this round.

Final run:
```
/Applications/MATLAB_R2025b.app/bin/matlab -batch "cd('/Users/msc/Desktop/optimal_control/booster_landing'); setup_paths; test_closed_loop_nominal"
```
```
Error using assert
nominal did not land (stop=arrest, alt=0.530 m) -- an arrest, not a touchdown
Error in test_closed_loop_nominal (line 25)
```

Dispersed now passes every one of its assertions (`miss=14.77<15`,
`vtd=1.377<2.0`, `landed=true`); the test still throws, now specifically and
only on the nominal `landed` assertion.

Per the round's acceptance criteria ("If after the Tnom fix + honest
classification + weight sweep the dispersed vtd still cannot meet 2.0, STOP,
report the best achieved numbers + the sweep table... don't force it"): the
literal trigger (dispersed `vtd` unfixable) did not occur — dispersed is
fully fixed. What the sweep instead surfaced is a **joint** constraint the
review didn't anticipate in that exact form: no scalar-`R`/`Qf(1:3)` point
in ~30 tested configurations puts *both* the nominal case (genuine `z=0`
crossing at `vtd<0.1`) and the dispersed case (`landed` inside both gates)
on the passing side simultaneously. Every direction tried moves the two
requirements in opposite directions along the same axis. This is the same
physics already diagnosed (`P.Tmin` exceeds vehicle weight at every mass on
this trajectory): a softened, saturation-tolerant gain schedule damps the
nominal approach into a stable near-ground hover instead of a crossing;
sharpening it enough to force a genuine nominal crossing reopens exactly the
saturated-direction-budget problem that fails the dispersed case. Per the
review's own framing, this reads as the anticipated spec-level question
("singular v(tf)=0 terminal BC with Tmin>mg") rather than a search failure,
and is not forced further here — STOPPING and reporting per instruction.
Candidate directions for a human decision, none attempted: (a) an
anisotropic/per-axis `R` (soften horizontal correction specifically, keep
vertical tight) — a materially bigger design change than a "first probe";
(b) add integral action or a terminal-velocity floor/deadband so the
controller stops chasing exact `z=0` once already within `pad_radius`/
`vtd_max`; (c) revisit whether `v(tf)=0` is the right terminal BC for the
guidance trajectory this controller tracks, given `Tmin>mg` structurally
forbids a graceful last meter under pure vertical thrust.

### Files changed (round 2)

- `lib/hs_quad_ctrl.m` (new — promoted from `certify/certify_pdg.m`)
- `certify/certify_pdg.m` (local function removed, calls shared copy)
- `lib/tvlqr_design.m` (`ctrl.Tnom` rebuilt from `hs_quad_ctrl`; default `R`
  raised to `1.5e-9`)
- `lib/sim_closed_loop.m` (event classification, mass-depletion event,
  looser tolerances, `sat_frac` via `trapz`, `control_law` takes model `P`,
  dropped unused `sol` arg, header rewording)
- `tests/test_closed_loop_nominal.m` (added `td.landed` assertions)

### Self-review (round 2)

- `test_certify_nominal` and `test_tvlqr_riccati` both re-run clean after
  every structural change in this round (Tnom rebuild, R default change) —
  not just at the end.
- The `hs_quad_ctrl` promotion is verified byte-identical logic (same
  Lagrange-basis formulas, same clamp, same segment-index arithmetic) —
  only the location and header changed.
- No `i`/`j` loop variables introduced.
- Did not attempt the three "candidate directions" above — each is a bigger
  design change than this round's "first probe" authorization covers, and
  the review explicitly said not to force a fix.
- The `R=1.5e-9` choice is a narrow window (evidenced by `R=1.4e-9` and
  `R=1.6e-9` both failing), which is itself worth flagging: it is not a
  robustly-margined operating point, and a slightly different dispersion
  case or grid (`N`) could plausibly land outside it. Not stress-tested
  against other dispersion vectors/magnitudes this round — the sweep tables
  above are all against the brief's single `[50;-30;0]` case.

---

## Round-3 fix report (2026-08-08) — terminal-BC adjudication

User adjudication: the round-2 finding (a 0.53 m nominal arrest, structural
proof that `v(tf)=0` is singular under `P.Tmin>weight`) led to a terminal-BC
change. **Guidance terminal condition becomes `r(tf)=0`, `v(tf)=[0;0;-1.5]`
m/s** (was `v(tf)=0`). This section documents the full cascade: physics
layer, both solvers, all six other test files, the certify gate re-run
(old vs new pins), a fresh TVLQR weight mini-sweep, and the honest final
test status.

### 1. `lib/booster_params.m`

Added `P.vf = [0;0;-1.5]` with a comment citing the adjudication and the
`Tmin>weight` rationale. Every other number unchanged.

Also changed `P.tf_hi` from 50 to 22 s — **not** part of the original
cascade instructions, but a genuine blocking bug uncovered while re-solving
(see §2 below): the convex solver's golden-section search hard-errored
because its two golden-ratio-derived initial probes (`tf=25.28`, `tf=34.72`)
both landed on invalid iterates under the new BC. A dedicated sweep found a
real, reproducible `Infeasible_Problem_Detected` wall near `tf~27 s` (that
much loiter time burns past `P.mdry`, since `Tmin` exceeds weight
everywhere — no free coast) and both solvers agree the true mass-optimal
`tf` is `~15.6 s`. `P.tf_hi=22` moves the golden search's own probe points
into consistently well-behaved territory and, as a side effect, tightens
`solve_pdg_colloc`'s `Tc=(tf_lo+tf_hi)/2` ND time scale (30→16 s) closer to
the true answer. `P.tf_hi` has exactly two consumers (verified by grep):
`solve_pdg_convex`'s bracket and `solve_pdg_colloc`'s free-tf upper
bound/scaling — both re-verified clean after the change (§2, §5).

### 2. `lib/solve_pdg_colloc.m` and `lib/solve_pdg_convex.m`

Terminal constraint changed from `X(4:6,end)==0` / `Vh(:,end)==0` to
`X(4:6,end)==P.vf/Vc` / `Vh(:,end)==P.vf/Vc`, ND-scaled consistently with
each solver's own velocity scale (`solve_pdg_colloc`'s `Vc=Lc/Tc` is fixed
by the bracket; `solve_pdg_convex`'s `Vc=Lc/tf` is local to each fixed-tf
subproblem, same as its existing `v0/Vc` handling).

**Two real bugs found and fixed in `solve_pdg_convex.m` while re-solving**
(both discovered empirically, not anticipated by the cascade instructions;
both are in-scope per this round's "fully authorized across files"):

**Bug A — no retry on a bad golden-section probe.** The two initial probes
(`tf=25.2786`, `tf=34.7214`) both returned invalid results, hard-erroring
the search before it could even start. Direct re-testing showed `tf=25.2786`
is **run-to-run nondeterministic**: the identical call measured a tight gap
(6.5e-5) in one MATLAB process and an untight gap (0.448) in another
(likely BLAS/IPOPT internal threading — separate processes, identical
inputs). `tf=34.72` is a genuine, reproducible infeasibility (the tf~27s
wall above). Fix: `solve_fixed_tf_probe` wraps every probe evaluation
(initial and loop-internal) with up to 3 attempts, keeping the best.

**Bug B — the real one: a deterministic misclassification that made the
golden search converge 52 kg short of the true optimum.** Even after Bug
A's fix, the search's final answer was `tf=16.3486, mf=26494.68` — but a
direct fixed-tf query at collocation's own `tf=15.6225` gave `mf=26546.15`
(matching colloc within the already-adjudicated ~0.7 kg Taylor-bound
floor) — 52 kg **higher**. A dense direct sweep (`tf=14..19` in 0.5 s
steps) confirmed the true `mf(tf)` peak sits at `tf~15.5-15.7`, and the
golden search's own internal trace showed a spurious `-Inf` reading at
`tf=16.3344`, right where it needed to keep exploring toward the true
peak — poisoning golden section's unimodal-search assumption. Root cause,
found by testing that exact tf value directly (reproducibly, 3/3 times):
IPOPT converges there to `status='Solved_To_Acceptable_Level'` (not the
literal string `'Solve_Succeeded'`) with an **excellent, tight gap**
(3.66e-5) — and the old `mf_or_neginf` rejected ANY non-`'Solve_Succeeded'`
status outright, regardless of gap. That blanket rule (motivated by a
real, different failure mode documented in the code: an untight
`Solved_To_Acceptable_Level` point at `tf=17` with gap~1e-2 that would
otherwise bias the search) was over-strict for a genuinely tight
acceptable-level point. Fix: `mf_or_neginf` now checks **tightness first**
— `code 2` (tight, acceptable-level) is valid alongside `code 3` (tight,
fully converged); the status-string distinction is kept only as a
diagnostic breadcrumb in `tf_curve`. A bare retry (Bug A's fix) could
never have fixed this — the misclassification was deterministic, not
noise (confirmed: identical 3/3 across separate MATLAB processes).

Verification: after both fixes, `solve_pdg_convex(P)`'s golden search
converges to `tf=15.6284, mf=26546.189` — matching colloc's `tf=15.6225,
mf=26546.892` to `Δtf=0.006 s`, `Δmf=0.70 kg` (right at the pre-existing
adjudicated Taylor-bound floor, not a new discrepancy). `tf_curve` shows
14 probes, all `code 2` or `code 3` (valid) except one genuine `code 0`
(solver threw) near the true peak and the correctly-classified `code 2`
point at `tf=16.33` (now counted, not poisoning the search).

Commands (all synchronous, confirming outputs in order: broken →
Bug A fix insufficient → Bug B fixed):
```
matlab -batch "...; P=booster_params(); solV=solve_pdg_convex(P); ..."
  -> Error: Both golden-section probes ... failed to solve.
matlab -batch "...; [after Bug A retry-wrapper] ..."
  -> Error: Both golden-section probes ... failed to solve, even after a retry each.
matlab -batch "...; [after P.tf_hi=22] ..."
  -> tf=16.3486 mf=26494.679  (52 kg short of the true peak -- Bug B)
matlab -batch "...; [after Bug B tightness-first fix] ..."
  -> tf=15.6284 mf=26546.189  (matches colloc to 0.70 kg)
```

### 3. Six other test files updated

- `tests/test_colloc_smoke.m`, `tests/test_convex_lossless.m`: terminal
  assertion changed to `max(abs(X(1:6,end) - [0;0;0;P.vf])) < tol` style,
  per the cascade instructions.
- `tests/test_closed_loop_nominal.m`: nominal block now asserts
  `td.landed==true`, `miss<1`, and `|vtd-1.5|<0.1` (tracking-error gate
  around the new nominal descent target, not `vtd<0.1` outright); dispersed
  block asserts `landed`, `miss<15`, `vtd<2.0` (mission gate unchanged).
- `tests/test_params.m`: comment #2 updated (was "arriving at v=0 exactly
  at touchdown", now correctly describes the `v=P.vf` target and cites the
  adjudication) — no assertion changed, this test never checked terminal
  velocity numerically, only min-throttle T/W.
- `tests/test_certify_nominal.m`: **run as-is, unmodified.** Re-run below.

### 4. Re-solve + re-certify: old vs new gate numbers

Commands:
```
matlab -batch "cd(...); setup_paths; P=booster_params(); solC=solve_pdg_colloc(P); solV=solve_pdg_convex(P); ...; rep=certify_pdg(solC, solV, P); print_certify_report(rep); sol=solC; save('results/pdg_colloc_nominal.mat','sol'); sol=solV; save('results/pdg_convex_nominal.mat','sol');"
```

| | Old (`v(tf)=0`, `tf_hi=50`) | New (`v(tf)=P.vf`, `tf_hi=22`) |
|---|---|---|
| Colloc `tf` | 15.6225 s | 15.6225 s |
| Colloc `mf` / fuel | 26546.892 / 3453.108 kg | 26546.892 / 3453.108 kg |
| Convex `tf` | 15.6284 s (old default `tf_hi`) | 15.6284 s |
| Convex `mf` / fuel | n/a (see below) | 26546.189 / 3453.811 kg |
| G1 max HS defect | 1.40e-07 | 1.40e-07 (unchanged — colloc solve identical) |
| G2 pos/vel/mass residual | 1.6e-4 m / 2.0e-5 / 1.3e-5 kg | 4.2e-4 m / 5.3e-5 / 4.6e-6 kg |
| G3 `\|dmf\|` | 0.703942 kg | 0.703665 kg |
| G3 `\|dtf\|` | 0.0097 s | 0.0059 s |
| G4 lossless gap | 1.4e-4 m/s² | 1.1e-4 m/s² |
| G5 bound frac / switches / primer | 0.9917 / 1 / 0.608 deg | 0.9917 / 1 / 0.608 deg |
| `all_pass` | PASS | PASS |

(Colloc's own tf/mf are identical in both rows above because that's the
SAME post-BC-change colloc solve reported both times — the "old" column
there is just restating it; the convex `tf`/`mf` shown for "old" is what
the golden search returned in the FIRST successful (Bug A only, Bug B
still present) re-solve, before Bug B's fix — included to show the
52 kg-error state existed transiently mid-investigation, not that it was
ever a shipped/certified number.)

**`test_certify_nominal.m`'s existing regression pins did not need
updating** — run as-is:
```
matlab -batch "cd(...); setup_paths; test_certify_nominal"
```
```
G3 |dmf| [kg]     0.709639 (coarse) / 0.703665 (nominal)   < 1     PASS
G2 mass residual  1.40197e-05 (coarse) / 4.61685e-06 (nominal)  < 0.01 PASS
G3 |dtf| [s]      0.00584361 (coarse) / 0.00586417 (nominal)   < 0.2  PASS
ALL GATES  PASS  (both blocks)
test_certify_nominal (coarse) PASS
test_certify_nominal (nominal, all_pass) PASS
```
The Taylor-bound model-error floor that the 1.0 kg `G3_dmf` pin was
adjudicated against is essentially unchanged by a 1.5 m/s terminal-velocity
shift (0.704 kg new vs 0.704/0.724 kg old) — expected, since that floor
comes from the convex relaxation's mass-bound linearization, not the
terminal BC.

`results/pdg_colloc_nominal.mat` and `results/pdg_convex_nominal.mat`
regenerated with the final (both bugs fixed, `tf_hi=22`) solutions.

### 5. Fast tests re-verified after every structural change (not just once)

```
test_params            PASS (min-throttle dry T/W = 1.346, unchanged)
test_dynamics_jac       PASS (vacuum + drag, unaffected by BC/tf_hi)
test_colloc_smoke       PASS  tf=15.62 s  mf=26546.9 kg  fuel=3453.1 kg
test_convex_lossless    PASS  gap=3.27e-05  mf=25724.1 kg
test_certify_nominal    PASS  (both blocks, all_pass -- see §4)
test_tvlqr_riccati      PASS  (P PSD, K finite, terminal P=Qf; re-checked
                                after both the Tnom rebuild in round 2 AND
                                the R-default change below)
```

### 6. TVLQR weight mini-sweep under the new (non-singular) BC

Rationale for re-sweeping (per this round's instruction): round 2's
`R=1.5e-9` default was chosen entirely to satisfy the dispersed case under
the OLD singular `v(tf)=0` BC, where the nominal case could never land
regardless of `R` — so `R` was picked with zero regard for nominal
performance. With the singularity structurally relaxed, it's worth
checking whether a different `R` now does better on both fronts.

**R sweep** (baseline `Q`, `Qf`; `sol=solve_pdg_colloc(P,struct('N',30))`,
`dsp=[50;-30;0]`):

| R | nom alt (m) | nom vtd | nom landed | nom miss | dsp miss (m) | dsp vtd | dsp landed |
|---|---|---|---|---|---|---|---|
| 1e-12 | -0.0000 | 11.338 | **1** | 0.204 | 2.92 | 10.73 | 1 |
| 3e-12 | -0.0000 | 7.344 | **1** | 0.035 | 2.63 | 6.80 | 1 |
| 1e-11 | -0.0000 | 1.662 | **1** | 0.013 | 2.75 | 4.41 | 1 |
| 3e-11 | -0.0000 | 1.365 | **1** | 0.016 | 3.13 | 5.47 | 1 |
| 5e-11 | -0.0000 | 1.093 | **1** | 0.014 | 3.42 | 5.80 | 1 |
| 7e-11 | -0.0000 | 0.808 | **1** | 0.006 | 3.68 | 5.91 | 1 |
| 1e-10 | 0.013 | 0.600 | 0 | 0.080 | 4.04 | 5.89 | 1 |
| 3e-10 | 0.323 | 0.643 | 0 | 0.028 | 6.17 | 4.80 | 1 |
| 5e-10 | 0.515 | 0.643 | 0 | 0.033 | 8.02 | 3.24 | 1 |
| 5.5e-10 | 0.553 | 0.643 | 0 | 0.034 | 8.44 | 2.85 | 1 |
| 6.0e-10 | 0.589 | 0.643 | 0 | 0.035 | 8.87 | 2.81 | 1 |
| 6.5e-10 | 0.623 | 0.642 | 0 | 0.037 | 9.25 | 1.62 | **1** |
| **7.0e-10** | 0.654 | 0.642 | 0 | 0.038 | **9.64** | **0.88** | **1** |
| 7.5e-10 | 0.683 | 0.642 | 0 | 0.039 | 10.03 | 0.52 | 0 (arrest) |
| 8.0e-10 | 0.711 | 0.641 | 0 | 0.039 | 10.41 | 0.53 | 0 (arrest) |
| 1.0e-9 | 0.807 | 0.640 | 0 | 0.042 | 11.84 | 0.55 | 0 (arrest) |
| 1.2e-9 | 0.886 | 0.639 | 0 | 0.045 | 13.13 | 0.58 | 0 (arrest) |
| 1.5e-9 (round-2 default) | 0.984 | 0.637 | 0 | 0.048 | 14.86 | 0.65 | 0 (arrest) |
| 2.0e-9 | 1.113 | 0.634 | 0 | 0.052 | 17.29 (>15) | 0.77 | 0 (arrest) |

The **dispersed-safe window** (`miss<15`, `vtd<2.0`, **and** `landed=true`)
is narrow: `R∈[6.5e-10, 7.0e-10]`. Below it, `vtd` fails (>2); above it,
the dispersed case reverts to an arrest (`landed=false`) even though its
raw `vtd` number looks fine — exactly the Critical-1/2 trap from round 2,
now recurring on the *dispersed* side at high `R`. `R=7.0e-10` sits at the
best-margin end of that window: `miss=9.64` (5.4 m headroom vs 15) and
`vtd=0.88` (1.1 m/s headroom vs 2.0) — both markedly more comfortable than
round 2's razor-thin `R=1.5e-9` (`miss=14.77`/`vtd=1.377`, ~0.2 m and
0.6 m/s headroom).

The **nominal-genuine-landing window** is `R≲7e-11` — and even there,
`vtd` only lands inside the required `[1.4,1.6]` band (`|vtd-1.5|<0.1`)
for a narrow sub-range around `R≈1.5-2e-11` (interpolating: `1e-11`→1.662,
`3e-11`→1.365, both just outside the band).

**These two windows are disjoint by roughly an order of magnitude** —
confirmed, not assumed: two secondary probes were run to check whether a
*different* axis could bridge the gap rather than `R` alone:

- `Qf(4:6)` (terminal velocity weight) at fixed `R=7e-10`, values
  `1/5/10/30/100`: nominal `alt`/`landed` essentially unchanged (still an
  arrest around 0.57-0.65 m), and dispersed `vtd` got **worse**, not
  better (0.88→2.48→2.61 by `Qf46=10`, breaking the dispersed pass) —
  consistent with round 2's finding that under this much saturation only
  thrust *direction* survives the clamp, so scaling a magnitude-only
  weight term does little or backfires.
- `Q(1:3)` (running position weight) at fixed `R=7e-10`, values
  `1e-4/3e-5/1e-5/3e-6/1e-6`: nominal essentially unaffected (alt stays
  0.65-0.67, `vtd` stays ~0.642, never lands) — expected, since the
  nominal case starts with zero position error to begin with, so
  de-weighting position tracking has nothing to bite on. Dispersed
  `landed` flipped inconsistently (1,0,0,1,1) across this axis with no
  clean trend — not a usable lever.

### 7. Final decision and honest test status

Adopted `R=7e-10` as `tvlqr_design.m`'s new default (was `1.5e-9`,
documented in that file with the full rationale). This is a genuine,
verified improvement for the dispersed case specifically (per this round's
own framing: "a less extreme setting... with more room" — confirmed: 5.4 m
/ 1.1 m/s of margin vs round 2's ~0.2 m / 0.6 m/s).

It does **not** achieve the round's stated acceptance bar in full. Final
run:
```
matlab -batch "cd(...); setup_paths; test_closed_loop_nominal"
```
```
Error using assert
nominal touchdown speed 0.642 (target 1.5)
Error in test_closed_loop_nominal (line 34)
```
Full final numbers (separate probe, same weights):
```
NOMINAL:    stop=arrest    landed=0  alt=0.6541 m  miss=0.0376 m  vtd=0.6419 m/s  sat_frac=0.683
DISPERSED:  stop=touchdown landed=1  alt=-0.0000 m miss=9.6433 m  vtd=0.8828 m/s  sat_frac=0.476
```
Dispersed passes every one of its assertions cleanly. Nominal fails two of
its three (`landed`, and `|vtd-1.5|<0.1`) — it still arrests, just at a
much smaller residual speed than before the BC change (0.64 m/s vs the
0.64-1.4 m/s range the OLD `v(tf)=0` target produced at comparable `R`,
which is a real improvement in absolute terms) and closer to, but not
within tolerance of, the new -1.5 m/s target.

**Correction to the adjudication's stated premise.** The rationale given
for the `P.vf` change was that a nonzero touchdown velocity "removes the
arrest state entirely: velocity never nulls above ground." The evidence
in this section does not support that as an unconditional claim for the
CLOSED LOOP (it is true, and verified, for the OPEN-LOOP guidance
trajectory — G2 already proved that lands to ~1.6e-4 m). Under closed-loop
TVLQR tracking, arrest is still possible and reproducible at every `R`
tested in the dispersed-safe band, because `P.Tmin` still exceeds vehicle
weight at every mass on this trajectory (unchanged by `P.vf`): closed-loop
tracking error that leaves the vehicle "behind schedule" relative to the
guidance's own descent-rate profile cannot be recovered near the ground,
since thrust can only ever net-decelerate, never accelerate the descent,
regardless of which velocity the schedule is currently asking for. The BC
change substantially IMPROVES nominal behavior (a genuine landing near the
target speed is now achievable at low `R`, which was never possible under
the old BC at any `R` without absurd overspeed) but does not structurally
eliminate the arrest failure mode under a weight setting that also has to
handle dispersions — the same class of tension found in round 2, now
narrower (order-of-magnitude gap in `R` instead of a complete absence of
overlap) but not closed.

Per the round-2-established STOP protocol (apply here by the same logic:
report the best-achieved numbers and the sweep evidence rather than force
a fix outside this round's authorized scope) — not attempted further:
anisotropic/per-axis `R` splitting horizontal from vertical gain
magnitude, integral action, a terminal deadband, or reconsidering whether
`v(tf)=[0;0;-1.5]` (still hitting `z=0` exactly) versus some
altitude-blended terminal condition is the right guidance BC for a
closed-loop-trackable trajectory under this thrust envelope. These are
larger design changes than a weight mini-sweep and are flagged here for
the human, not attempted.

### 8. Files changed (round 3)

- `lib/booster_params.m` (`P.vf` added; `P.tf_hi` 50→22, documented)
- `lib/solve_pdg_colloc.m` (terminal velocity BC, ND-scaled)
- `lib/solve_pdg_convex.m` (terminal velocity BC; golden-section retry
  wrapper; tightness-first probe classification — two real bugs fixed)
- `lib/tvlqr_design.m` (default `R` 1.5e-9→7e-10, documented)
- `tests/test_closed_loop_nominal.m` (nominal block: `|vtd-1.5|<0.1`
  replaces `vtd<0.1`; dispersed block unchanged in substance)
- `tests/test_colloc_smoke.m`, `tests/test_convex_lossless.m` (terminal
  assertion targets `[0;0;0;P.vf]`, not zero)
- `tests/test_params.m` (comment #2 corrected, no assertion changed)
- `results/pdg_colloc_nominal.mat`, `results/pdg_convex_nominal.mat`
  (regenerated, gitignored, not committed)

### 9. Self-review (round 3)

- Every structural change (BC edit, `tf_hi` edit, both convex-solver bug
  fixes, `R` default change) was followed by an immediate, synchronous
  re-run of the tests it could plausibly affect, not batched to the end —
  this is exactly the discipline the round-2 background-task incident
  should have taught, and this round ran everything synchronously with
  explicit 600000 ms timeouts throughout, no backgrounding.
- Bug A (retry) and Bug B (tightness-first classification) are DISTINCT,
  both real, both verified independently: Bug A's fix alone was
  insufficient (confirmed: identical failure reproduced verbatim after
  adding the retry, because Bug B's misclassification is deterministic,
  not noise a retry can dodge); Bug B's fix alone (without Bug A) was not
  separately tested since the previous total exposure already
  demonstrated why the retry was legitimately needed for a real
  nondeterminism `mf_or_neginf` cannot see.
- `P.tf_hi`'s two consumers were checked by grep before changing it (not
  assumed) — `solve_pdg_colloc`'s own free-tf search still converges to
  the correct `tf=15.62 s`, comfortably inside the new `[10,22]` bound,
  confirmed by the fast-test re-run.
- The G3 `|dmf|` gate's near-invariance to the terminal-velocity change
  (0.704 kg new vs ~0.70-0.94 kg old adjudicated band) was checked, not
  assumed, before concluding no pin update was needed.
- Did not touch `Qf(1:3)` this round (round 2 already showed it trades
  against dispersed miss badly); did test `Qf(4:6)` and `Q(1:3)` as the
  two most plausible remaining single-axis levers, both negative results,
  both reported rather than omitted.
- The "Correction to the adjudication's stated premise" section above is
  the single most important finding of this round and is stated plainly,
  not softened — the human asked for a change on the premise that it
  would fully resolve the arrest; the evidence says it substantially
  helps but does not fully resolve it under a weight setting that also
  has to pass the dispersed case, and that gap is reported rather than
  talked around.

---

## Round-4 fix report (2026-08-08) — N=30 artifact, feasibility gate, terminal-phase remedy

A third review pass ran decisive experiments and found the round-3
"disjoint R-bands" impasse was, in significant part, an **N=30 test-grid
artifact**: `hs_quad_ctrl`'s reconstruction left the thrust annulus
between nodes (11.7% of the flight below `Tmin`, worst dip -18%), the
sim's saturating clamp turned that into a persistent open-loop altitude
bias baked into the feedforward itself, and two rounds of TVLQR weight
tuning were — in hindsight — buying feedback gain to cancel that
feedforward defect rather than tuning against real closed-loop physics.
At production `N=60` the reconstruction was already annulus-feasible and
the reviewer's own measurement showed the nominal case lands at every `R`
from `7e-10` to `1e-4`. This section implements the fix in the specified
order, verifies each step, and reports the honest final state.

### 1. `hs_quad_ctrl` made annulus-feasible by construction

Root cause: the single-vector quadratic Lagrange reconstruction
`Tv = L0*U_k + L1*Um_k + L2*U_{k+1}` interpolates a vector whose
*direction* also turns across a segment; `L0+L1+L2=1` but the individual
`L_i` are not all nonnegative (`L1` exceeds 1 near the segment midpoint),
so this is not a convex combination and the interpolated *magnitude* can
fall outside `[Tmin,Tmax]` even when all three sampled points sit exactly
on the annulus (a "chord vs arc" effect where the direction itself
bends). At `N=30` this was a genuine problem (11.7% below `Tmin`); at
`N=60` the segments are short enough that direction barely turns within
one, so the defect was much smaller but not, until now, provably absent.

Fix: `lib/hs_quad_ctrl.m` reconstructs **direction** (the same
quadratic vector interpolant, then normalized) and **magnitude** (an
INDEPENDENT quadratic Lagrange interpolant of the three scalar `|T|`
samples, then clamped to `[Tmin,Tmax]`) separately, so a violation is
impossible by construction regardless of how the direction moves within
a segment. New signature: `hs_quad_ctrl(tt,U,Um,h,N,Tmin,Tmax)` (two new
trailing args); both call sites (`certify/certify_pdg.m`'s G2,
`lib/tvlqr_design.m`'s `ctrl.Tnom`) updated to pass `P.Tmin`, `P.Tmax`.

Verification:
```
matlab -batch "...; test_certify_nominal"
```
```
(coarse, N=40)   G2 pos residual  0.000374698 m   (was 0.000356366 m pre-fix)
(nominal, N=60)  G2 pos residual  0.000432751 m   (was 0.000422154 m pre-fix)
ALL GATES  PASS  (both blocks)
```
Essentially unchanged at both grids — confirms the split-direction/
magnitude form reduces to (numerically very close to) the prior formula
wherever it was already feasible, exactly as expected. Direct annulus
measurement (2000-point dense grid per solve):
```
N=30: belowTminFrac=0.0000  worstDip=-0.0000  aboveTmaxFrac=0.0000  minTmag=338000.0  maxTmag=845000.0
N=60: belowTminFrac=0.0000  worstDip= 0.0000  aboveTmaxFrac=0.0000  minTmag=338000.2  maxTmag=844999.9
```
The N=30 violation (11.7% below Tmin, -18% worst dip, reviewer's own
measurement of the PRE-fix formula) is eliminated at the source; N=60's
already-small margin is now exactly zero too.

### 2. Feedforward-feasibility gate added to `certify_pdg.m` (new G2ff)

Even though the fix above makes a violation structurally impossible, the
reviewer asked for a permanent gate so a future change to the
reconstruction (or a different grid) can never silently reintroduce this
class of defect without a test catching it — exactly the gap that let
two rounds of weight-tuning proceed on a phantom problem. `certify_pdg.m`
now flies `hs_quad_ctrl` over a DENSE grid (20 points/segment, not just
nodes+midpoints, so a dip strictly between samples cannot hide) and
reports:
```
rep.G2ff_below_tmin = max(max(0, Tmin - |T|) / Tmin)   over the dense grid
rep.G2ff_above_tmax = max(max(0, |T| - Tmax) / Tmax)
rep.G2ff_pass = both < 1e-6
```
Gated at EVERY grid (not just nominal, no `tolScale` accommodation) since
the direction/magnitude split guarantees feasibility regardless of
resolution — a real regression here should never need a coarse-grid
exemption. Folded into `rep.all_pass`. `print_certify_report.m` updated
to display it. Measured: exactly `0` at both grids (§1's verification
output above shows the full `test_certify_nominal` run with the new
`G2ff` rows, both `PASS`).

### 3. Nominal gate fixed (`|vtd-1.5|<0.1` was a hidden altitude requirement)

The reviewer's math: on the near-vertical terminal arc, `dvz/dz ~ 14.7/s`
near touchdown, so a `|vtd-1.5|<0.1` window is equivalent to requiring
touchdown altitude within `~6.8 mm` of a specific target — nothing the
mission actually cares about, and a needlessly fragile assertion.
`tests/test_closed_loop_nominal.m`'s nominal block now asserts the real
mission gates only: `td.landed==true`, `miss<1 m`, `vtd<=P.vtd_max`
(2.0 m/s) — and PRINTS `vtd` and terminal altitude as diagnostics rather
than gating on them tightly. Same treatment applied to the dispersed
block's diagnostics line for consistency.

### 4. `test_closed_loop_nominal.m` moved to the production grid (N=60)

`sol = solve_pdg_colloc(P, struct('N', 30))` → `struct('N', P.N)` (60).
This is the change that actually exposed §1's root cause once made:
running the OLD (round-3) `test_closed_loop_nominal.m` logic at N=60
with the then-current default (`R=7e-10`) immediately showed the nominal
case landing cleanly (`alt=-0.0000`, `vtd=1.849`, `stop=touchdown`) —
confirming the reviewer's diagnosis before any further code changes were
made:
```
matlab -batch "...; test_closed_loop_nominal"
```
```
nominal diagnostics: vtd=1.849 m/s (target 1.5, gate <=2.0)  alt=-0.0000 m  stop=touchdown
dispersed diagnostics: vtd=7.594 m/s  miss=9.944 m  alt=-0.0000 m  stop=touchdown
Error: dispersed vtd 7.59 m/s
```
Dispersed, however, now fails WORSE than round 3's N=30 numbers
(`vtd=7.59` vs round 3's `vtd=0.88` at the same `R=7e-10`) — confirming
round 3's chosen default was tuned against the N=30 artifact and does not
transfer to the grid the shipped defaults are meant for. This is the
"antagonism" §5 below investigates.

### 5. R re-derived at N=60; structural remedy implemented; antagonism confirmed to survive it

**R sweep at N=60, post-§1 fix, BEFORE the terminal-phase remedy**
(baseline `Q`,`Qf`; `sol=solve_pdg_colloc(P,struct('N',P.N))`,
`dsp=[50;-30;0]`):

| R | nom alt | nom vtd | nom landed | nom miss | dsp miss (m) | dsp vtd | dsp landed |
|---|---|---|---|---|---|---|---|
| 7e-10 | -0.0000 | 1.849 | **1** | 0.003 | 9.94 | 7.594 | 1 |
| 1e-9 | -0.0000 | 1.410 | **1** | 0.000 | 12.16 | 7.328 | 1 |
| 1.5e-9 | -0.0000 | 1.480 | **1** | 0.000 | 15.18 | 6.757 | 1 |
| 3e-9 | -0.0000 | 1.512 | **1** | 0.000 | 21.32 | 5.587 | 1 |
| 5e-9 | -0.0000 | 1.463 | **1** | 0.000 | 26.21 | 4.545 | 1 |
| 1e-8 | -0.0000 | 1.671 | **1** | 0.001 | 32.81 | 3.157 | 1 |
| 3e-8 | -0.0000 | 1.770 | **1** | 0.001 | 42.74 | 1.766 | 1 |
| 5e-8 | -0.0000 | 1.993 | **1** | 0.002 | 46.83 | 1.455 | 1 |
| 8e-8 | -0.0000 | 2.098 (>2, over) | 1 | 0.002 | 50.03 | 1.129 | 1 |
| 1e-7 | -0.0000 | 2.129 (>2) | 1 | 0.003 | 51.31 | 0.916 | 1 |
| 3e-7 | -0.0000 | 2.206 (>2) | 1 | 0.004 | 55.57 | 1.925 | 1 |
| 1e-6 | 0.0000 | 2.230 (>2) | 1 | 0.004 | 57.44 | 2.246 | 1 |

Confirms the reviewer's headline number exactly: **nominal genuinely
lands at every R tested** (`landed=1` throughout, `alt≈0`) — this axis is
fully fixed by §1. (Nominal `vtd` DOES eventually exceed its own
`vtd_max=2.0` gate above `R≈5-8e-8`, narrowing the usable range from the
top as well — not part of the reviewer's original framing but measured
here.) Dispersed shows the SAME antagonism as round 3, just shifted:
`miss<15` needs `R<~1.5e-9`; `vtd<2` needs `R>~5e-8` (using a stricter,
non-boundary read: `3e-8` gives a real margin at `1.766`) — roughly a
20-30x gap, no overlap.

**Structural remedy implemented** (`lib/sim_closed_loop.m`, reviewer-
endorsed): below `P.zTermBand` altitude (150 m, new `booster_params.m`
field), `control_law` (a) excludes z-position error from the feedback
entirely (`Kt(:,3)=0`) and instead tracks the guidance's own `v(z)`
profile at the vehicle's ACTUAL altitude (`ct.vOfZ(x(3))`, built once
per `sim_closed_loop` call from `sol.X(3,:)`/`sol.X(4:6,:)`, not
per-RHS-call), and (b) fixes the "past-tf freeze": `ctrl.tgrid(end)=tf`
is where `P(tf)=Qf` (diagonal) forces `K(tf)`'s position-feedback columns
to be EXACTLY zero (`B`'s position rows are always zero since `rdot=v`
doesn't depend on thrust, so `K=Rinv*B'*P` inherits zero columns wherever
`P` is diagonal in those rows) — and because `control_law` clamps its
time index to `tgrid(end)` for any `t>=tf`, a dispersed case that
(correctly) takes longer than `tf` to correct a 50 m offset was flying
with **zero position feedback in ALL THREE axes**, not just altitude, for
the tail of its flight. The fix holds the gain's time-index at
`ct.tTermEntry` (the nominal time altitude first crosses 150 m, found
once via `interp1` on `sol`) for the whole terminal phase and beyond, so
x,y lateral feedback survives past tf while z is zeroed by design.

**R sweep AFTER the terminal-phase remedy** (same grid/dispersion):

| R | nom vtd | nom landed | dsp miss (m) | dsp vtd | dsp landed |
|---|---|---|---|---|---|
| 7e-10 | 1.861 | 1 | 8.58 | 7.068 | 1 |
| 1e-9 | 1.511 | 1 | 10.85 | 6.819 | 1 |
| 1.5e-9 | 1.533 | 1 | 13.99 | 6.267 | 1 |
| 3e-9 | 1.541 | 1 | 20.44 | 5.162 | 1 |
| 5e-9 | 1.504 | 1 | 25.56 | 4.162 | 1 |
| 1e-8 | 1.678 | 1 | 32.42 | 2.804 | 1 |
| 3e-8 | 1.773 | 1 | 42.58 | 1.680 | 1 |
| 5e-8 | 1.994 | 1 | 46.73 | 1.482 | 1 |
| 8e-8 | 2.098 (>2) | 1 | 49.96 | 1.181 | 1 |
| 1e-7 | 2.130 (>2) | 1 | 51.26 | 0.974 | 1 |
| 3e-7 | 2.206 (>2) | 1 | 55.55 | 1.924 | 1 |

Nominal and dispersed-miss numbers move only marginally (dispersed miss
improves slightly at low R, e.g. `9.94→8.58` at `R=7e-10`); dispersed
`vtd` is essentially UNCHANGED from before the remedy at every R (e.g.
`7e-10`: `7.594→7.068`; `3e-8`: `1.766→1.680`). **The antagonism survives
the structural remedy**: `miss<15` still needs `R<~1.5e-9` (borderline at
`13.99`), `vtd<2` still needs `R>~5e-8`, no overlap.

This says the "past-tf freeze" mechanism, while real and now fixed, was
not the dominant driver of the dispersed antagonism. The more likely
remaining explanation (not new to this round -- round 2 first raised it):
a single 3-DOF thrust vector correcting a 50 m lateral offset under
`P.Tmin>weight` genuinely has to share its saturated-magnitude "budget"
between horizontal correction and vertical braking; softening `R` buys
back vertical authority by sacrificing how fast/hard horizontal error is
corrected, which is exactly why `miss` grows monotonically with `R` in
both sweeps above. This looks like an actuator-authority property of the
problem (one thrust vector, one magnitude budget, two simultaneous
demands), not a controller-design defect fixable by more gain-schedule
engineering.

### 6. Final decision and honest test status

Per the round's own explicit instruction ("if the antagonism survives
even the structural remedy, STOP with the sweep table") — not forced
further. Default `R` in `lib/tvlqr_design.m` set to `1e-9` (documented,
full rationale in that file): chosen because it is the point in the
tested range that gets the NOMINAL case fully right (`landed=true,
vtd=1.51` — comfortable margin under `vtd_max=2.0`) while keeping
dispersed MISS comfortably inside `pad_radius` (`10.85 m` vs `15 m`, the
geometric "did it reach the pad" requirement) — honestly reporting that
dispersed `vtd` (`6.82 m/s`) badly fails the `vtd_max` safety gate at
this setting. No `R` in the tested range satisfies both dispersed gates
simultaneously.

Final full-suite run (synchronous, 600000 ms):
```
matlab -batch "cd(...); setup_paths; test_params; test_dynamics_jac; test_colloc_smoke; test_convex_lossless; test_certify_nominal; test_tvlqr_riccati; test_closed_loop_nominal"
```
```
test_params              PASS (min-throttle dry T/W = 1.346)
test_dynamics_jac        PASS (vacuum + drag)
test_colloc_smoke        PASS  tf=15.62 s  mf=26546.9 kg  fuel=3453.1 kg
test_convex_lossless     PASS  gap=4.80e-05  mf=26440.8 kg
test_certify_nominal     PASS  (both blocks, all_pass, G2ff new gate PASS at both)
test_tvlqr_riccati       PASS
test_closed_loop_nominal:
  nominal diagnostics: vtd=1.511 m/s (target 1.5, gate <=2.0)  alt=-0.0000 m  stop=touchdown
  dispersed diagnostics: vtd=6.819 m/s  miss=10.851 m  alt=-0.0000 m  stop=touchdown
  CAUGHT: dispersed vtd 6.82 m/s
```
6 of 7 tests genuinely green. `test_closed_loop_nominal` throws
specifically and only on the dispersed `vtd` assertion — the nominal
block (the primary target of this round's investigation) now passes
cleanly, confirming the N=30-artifact diagnosis was correct and complete
for that axis. The dispersed `vtd` antagonism is a separate, deeper
finding not resolved by either fix applied this round.

Candidate directions for a human decision, none attempted (bigger design
changes than this round's scope): (a) accept a softer dispersed test case
(e.g. a smaller offset than 50 m, if that is representative of the real
dispersion budget) rather than tune weights against a worst case that may
not be physically achievable by any single-thrust-vector controller; (b)
an anisotropic/per-axis `R` that penalizes horizontal control effort more
than vertical, so horizontal correction is intentionally gentler and
vertical braking authority is preserved even for large lateral offsets;
(c) a genuinely different guidance re-plan on dispersion detection
(closed-loop replanning) rather than pure linear feedback around a single
fixed nominal trajectory; (d) revisit whether `P.pad_radius=15 m` and
`P.vtd_max=2.0 m/s` are jointly achievable requirements for a 50 m offset
given the vehicle's actuator envelope, independent of controller design.

### 7. Small fixes (item 5)

- `lib/booster_params.m`: `P.vf`'s comment corrected -- no longer claims
  the terminal-velocity change "removes the arrest state entirely"
  (retracted in round 3's own report, now corrected at the source too).
  Added `P.tf_hi`'s Phase-2/drag caveat (the `tf~27 s` wall was validated
  in vacuum only). Added `P.zTermBand=150` with its own documentation.
- `tests/test_convex_lossless.m`: pin moved `tf=25` (now outside
  `[tf_lo,tf_hi]=[10,22]`) → first tried `tf=18` (a mistake — round 3's
  own sweep had already flagged `tf=18` as untight, `gap=4.702`, missed
  when picking the replacement) → corrected to `tf=17` (verified clean
  and reproducible, two re-solves, identical `Solve_Succeeded`,
  `gap=4.8e-5`, at the test's actual `Nconv=60`).
- `validity_code` semantics: already documented at its definition
  (`solve_pdg_convex.m`'s `mf_or_neginf`, and in `solve_pdg_convex`'s own
  header) since round 3's fix for the same function — re-verified current
  and accurate (code 2=valid+acceptable-level, 3=valid+fully-converged,
  both now count as `v=mf`) rather than re-documented from scratch.

### 8. Files changed (round 4)

- `lib/hs_quad_ctrl.m` (direction/magnitude split, new `Tmin,Tmax` args)
- `certify/certify_pdg.m` (new G2ff feedforward-feasibility gate; updated
  `hs_quad_ctrl` call site)
- `certify/print_certify_report.m` (G2ff display rows)
- `lib/tvlqr_design.m` (`ctrl.Tnom`'s `hs_quad_ctrl` call site updated;
  default `R` 7e-10→1e-9, documented)
- `lib/sim_closed_loop.m` (altitude-scheduled terminal-phase feedback,
  `ct` schedule struct threaded through `control_law`/`plant_rhs`)
- `lib/booster_params.m` (`P.vf` comment corrected; `P.tf_hi` drag caveat;
  new `P.zTermBand`)
- `tests/test_closed_loop_nominal.m` (N=30→P.N=60; nominal gate fixed to
  mission gates + diagnostics, not a disguised altitude requirement)
- `tests/test_convex_lossless.m` (pin `tf=25`→`17`)

### 9. Self-review (round 4)

- Every fix verified in the order specified (hs_quad_ctrl → G2ff gate →
  nominal test gate → N=60 → R re-sweep → structural remedy → re-sweep
  again), each with its own synchronous MATLAB run before moving to the
  next, not batched.
- `test_certify_nominal` re-run clean after BOTH `hs_quad_ctrl` (§1) and
  the new G2ff gate (§2) were added, not just once at the end.
- Caught and fixed my own mistake in §7 (the `tf=18` pin) by checking
  against round 3's OWN sweep data before finalizing `tf=17` — a case of
  not re-deriving from scratch what was already measured two rounds ago.
- The terminal-phase remedy's `vOfZ`/`tTermEntry` are built ONCE per
  `sim_closed_loop` call (from `sol`, outside the ODE integration), not
  recomputed on every RHS evaluation -- checked for a performance
  regression by confirming the sweeps above (22+ MATLAB invocations this
  round) completed in normal time, no unexpected slowdown.
- Did not implement any of the four "candidate directions" in §6 -- each
  is a bigger design or requirements change than a weight/schedule
  retune, and the round's own instruction was explicit: stop and report
  if the antagonism survives the structural remedy, don't keep forcing.
- The G2ff gate is deliberately NOT `tolScale`-accommodated at coarse
  grids (unlike G3) because the fix that backs it (direction/magnitude
  split) is a by-construction guarantee independent of grid resolution --
  confirmed empirically (both grids measure exactly 0 violation) rather
  than asserted.

---

## Round-5 fix report (2026-08-08) — the dispersed 58 m case, SOLVED

Round 4 closed with the verdict that the dispersed `dr0=[50;-30;0]` case
was up against "an actuator-authority property of the problem (one thrust
vector, one magnitude budget, two simultaneous demands), not a controller-
design defect fixable by more gain-schedule engineering."

**That verdict was wrong, and this round proves it by fixing the case.**
Both gates now pass together, with margin:

| | round-4 shipped | round-5 shipped |
|---|---|---|
| nominal | landed, miss 0.013 m, vtd 1.511 m/s | landed, miss 0.013 m, vtd **1.453** m/s |
| dispersed 58 m | landed, miss 10.851 m, **vtd 6.819 m/s (FAIL)** | landed, miss **0.980 m**, vtd **1.404 m/s** |
| dispersed `sat_frac` | 0.653 | **0.129** |

`miss` has 14 m of headroom against `pad_radius=15`; `vtd` has 0.6 m/s
against `vtd_max=2.0`. The suite is 7/7 green.

### 1. Why rounds 2-4 could not find it: the sweeps were on the wrong axis

Three rounds swept `R` (over five decades), `Qf(1:3)`, `Qf(4:6)` and
`Q(1:3)`, and every one of them moved dispersed `miss` and dispersed `vtd`
in opposite directions. The reason is a one-line piece of LQR algebra that
nobody had written down. For the lateral double integrator (`A=[0 1;0 0]`,
`B=[0;1/m]`, weights `qPos`, `qVel`, `r`) the algebraic Riccati equation
solves in closed form to

```
K_r/m = sqrt(qPos/r)/m                       (position gain, "a")
K_v/m = sqrt((2 m sqrt(qPos r) + qVel)/r)/m  (velocity gain, "b")
```

and when `qVel` dominates the radical (which it did, by 20x, at the old
weights) the **slow closed-loop pole is `a/b = sqrt(qPos/qVel)`** —
**`r` cancels out entirely.** At the round-4 default
`Q=diag([1e-4 1e-4 1e-4 1e-2 1e-2 1e-2 0])` that is `sqrt(1e-4/1e-2) =
0.1 rad/s`: a **10 s time constant on a 15.6 s flight**. Sweeping `R`
rescales `a` and `b` together and cannot change it. Measured confirmation
from the shipped round-4 controller (`diag1`, N=60):

```
 gains (per-mass): t   a=K(1,1)/m   b=K(1,4)/m   slowpole=a/b
  0.00                    0.01949      0.2334      0.08351
  7.71                    0.03168      0.2669      0.11870
 13.89                    0.01525      0.4274      0.03569
```

A 0.084-0.12 rad/s pole decays a 58 m offset by `exp(-0.11*15.6) = 0.18`
— i.e. to ~10 m — which is *exactly* the 10.851 m miss round 4 shipped.
The controller was not fighting an authority limit; it was converging on
schedule, and the schedule was ten times too slow.

### 2. The measurement that proves it is not an authority limit

The same run's time history (dispersed case, round-4 controller):

```
  t    horiz_err  dx     dy      z      vz    vznom  |T|/kN  tilt_deg
 0.00    554.44   50.00  -30.00   2000.0 -180.00 -180.00   346.5     27.7
 3.16    438.83   45.91  -27.55   1434.5 -176.96 -176.96   339.9     13.1
 6.33    297.68   36.45  -21.88    882.2 -171.92 -171.92   338.0      2.3
 7.51    242.89   32.27  -19.37    680.4 -164.36 -164.40   844.3      2.7
11.47     82.77   18.62  -11.20    183.4  -86.44  -86.38   845.0     10.7
15.42     10.85    9.14   -5.61     -0.0   -5.88   -5.60   845.0     16.4
```

Read the `|T|` column. For the first 7 s the commanded thrust sits at
**338-346 kN — the Tmin floor — with 507 kN of magnitude margin unused**,
while `dx` crawls from 50 m to 33 m. The tracker was asking for
`K(1,1)*58 m ~ 18 kN` of lateral thrust when 507 kN was free. It then
carried a 32 m offset into the Tmax braking arc, where the *only* way to
buy lateral acceleration is to tilt a saturated vector — the `tilt_deg`
column climbs to 16.4 deg, and `cos(16.4 deg)` of the braking authority
is spent on lateral instead of vertical. That is the entire `vtd=6.82`
failure. The authority existed; it was requested in the one phase where
it costs the most and declined in the phase where it was free.

### 3. Fix (a): phase-scheduled Q (`lib/tvlqr_design.m`)

`Q` is now `Q(t)`, blended by `tanh((t-tSwitch)/tBlend)` (tBlend=0.4 s,
so the Riccati RHS stays smooth):

- **Phase A** (`t < tSwitch`, the Tmin margin arc): lateral position and
  velocity weights designed by inverting the algebra in §1 —
  `qPos = m^2 r w^4`, `qVel = m^2 r w^2 (4 zeta^2 - 2)`. Shipped
  `w=0.70 rad/s`, `zeta=1.00` -> `QA = diag([0.21 0.21 1e-4 0.85 0.85 1e-2 0])`.
- **Phase B** (`t > tSwitch`, the Tmax braking arc): `QB` = the round-4
  constant `Q`, **verbatim, unchanged**. Once the offset is dead, a few
  metres of residual error costs <1 deg of tilt and needed no change.
- `tSwitch` is auto-detected from `sol` by `annulus_switch` (new local
  function): the single crossing of `|T*|` through mid-annulus, which the
  G5 gate already certifies is unique. Measured `tSwitch = 7.160 s`.
- The z-position weight and all vertical weights are **left at their old
  values** — deliberately more surgical than the anisotropic-`R` form of
  the same idea, because the vertical loop demonstrably worked already
  (nominal `vz` tracked the guidance to <0.3 m/s through the whole
  braking arc at the old weights). Only the axis that was measurably
  starved is fed.

Design formula validated against the actual Riccati solution:
predicted `K(1,1)/m = w^2 = 0.49` vs **measured 0.4808**; predicted
`K(1,4)/m = 2*zeta*w = 1.40` vs **measured 1.3773**.

Interface preserved: `tvlqr_design(sol,P,opts) -> ctrl` still returns
`.tgrid/.K/.Pt/.xnom/.Tnom`; `opts.Q` still accepted and now sets
`QA=QB=Q` (constant weight, schedule disabled) so any caller written
against the old interface behaves exactly as before.

### 4. Fix (b): vertical-priority annulus allocation (`lib/sim_closed_loop.m`)

New local function `allocate_thrust`. The brief's law is "magnitude
clamped to `[Tmin,Tmax]`, direction preserved", and that is what it
returns **whenever the raw command fits under Tmax** — which, after fix
(a), is essentially the whole acceptance case (`sat_frac` 0.653 -> 0.129).
It differs only in the over-Tmax branch: preserving direction there scales
the *vertical* component down in lockstep with an unaffordable lateral
demand, i.e. pays for lateral correction with braking authority at the
worst possible moment. Vertical priority instead keeps the commanded
vertical component and spends only the leftover annulus radius,
`sqrt(Tmax^2 - Tz^2)`, on lateral. The `|T| >= Tmin` floor stays
direction-preserving (a throttle floor is not a budget to allocate).

Measured contribution, isolated (round-4 weights, dispersed 58 m case):
allocation alone took `vtd 6.819 -> 3.396` and `miss 10.851 -> 7.901`.
It is a real improvement on its own but does not pass the gate; fix (a)
is what closes it. Together: `vtd 1.404`, `miss 0.980`.

### 5. Sweeps

All at **N=60** (production grid; the N=30 trap of rounds 2-3 was not
repeated), against a cached `solve_pdg_colloc(P,struct('N',P.N))`,
varying only `tvlqr_design` opts.

**First scheduled sweep** — lateral `qA` at fixed `zeta=0.707` (LQR's
natural damping), `qB` in {1e-6, 1e-4}. Showed the mechanism works
(dispersed `vtd` 6.8 -> 2.0) but with bad `miss`: the `zeta=0.707`
response *overshoots*, nulling the 50 m offset by t=5.3 s and then
drifting to `dx=-12 m` by the switch, which the weak phase-B gain never
recovers. Diagnosed from the `dx` time history, not guessed.

**(w, zeta) sweep** — reparametrized per §3, `qB` = old default:

| w | zeta | qPosA | qVelA | nom land/vtd/miss | dsp land/vtd/miss |
|---|---|---|---|---|---|
| 0.4 | 0.9 | 0.0223 | 0.173 | 1 / 1.666 / 0.001 | 1 / 4.148 / 0.807 |
| 0.4 | 1.1 | 0.0223 | 0.395 | 1 / 1.511 / 0.002 | 1 / 1.751 / 6.843 |
| 0.4 | 1.4 | 0.0223 | 0.813 | 1 / 1.507 / 0.003 | **0** / 1.604 / 13.725 |
| 0.5 | 0.9 | 0.0544 | 0.270 | 1 / 2.249 / 0.003 | 1 / 1.173 / 3.330 |
| 0.5 | 1.1 | 0.0544 | 0.618 | 1 / **2.872** / 0.012 | 1 / 0.906 / 3.794 |
| 0.6 | 0.9 | 0.1128 | 0.389 | 1 / 1.934 / 0.001 | 1 / 1.425 / 3.455 |
| 0.6 | 1.1 | 0.1128 | 0.890 | 1 / 1.477 / 0.003 | 1 / 1.166 / 1.896 |
| 0.6 | 1.4 | 0.1128 | 1.830 | 1 / 1.542 / 0.001 | **0** / 0.899 / 8.138 |
| 0.7 | 1.1 | 0.2089 | 1.211 | 1 / 1.536 / 0.001 | 1 / 1.993 / 0.866 |
| 0.7 | 1.4 | 0.2089 | 2.490 | 1 / 1.500 / 0.002 | 1 / 1.295 / 6.039 |
| 0.8 | 0.9 | 0.3565 | 0.691 | 1 / 1.422 / 0.005 | 1 / 4.270 / 1.484 |
| 0.8 | 1.1 | 0.3565 | 1.582 | 1 / 1.536 / 0.001 | 1 / 3.947 / 0.511 |

(`land=0` = arrest, i.e. failed to reach the ground at all.)

**Final selection sweep** — scored on a *pure-lateral* battery (the axis
this round owns), each cell `vtd/miss`, `A` = arrest:

```
cfg              nom        D1_58      -D1        [100 0 0]  [0 -100 0] [141 141 0] [-100 -100 0] D1+dv0  | worst_vtd worst_miss nland
w0.55 z1.00   1.51/ 0.2  1.01/ 0.2  1.04/ 0.2  0.20/ 0.4A 1.88/ 0.2  2.26/ 0.3  1.34/ 0.2  0.51/ 0.2A |   2.26     0.40    6/8
w0.60 z1.00   1.37/ 0.0  1.29/ 0.7  1.04/ 0.6  0.27/ 1.2A 2.01/ 1.0  2.48/ 1.6  1.46/ 1.3  0.81/ 0.6  |   2.48     1.58    7/8
w0.60 z1.10   1.43/ 0.0  1.13/ 1.9  1.33/ 1.9  0.34/ 3.2A 1.96/ 3.4  2.69/ 8.4  1.77/ 5.0  0.56/ 2.6A |   2.69     8.38    6/8
w0.65 z1.10   1.43/ 0.0  1.28/ 1.3  1.14/ 1.3  0.22/ 2.3A 2.13/ 2.4  2.91/ 6.3  1.80/ 3.6  1.01/ 1.8  |   2.91     6.25    7/8
w0.70 z1.00   1.43/ 0.0  1.37/ 1.0  1.29/ 0.9  0.66/ 1.8  2.24/ 1.6  2.82/ 3.7  2.02/ 2.2  1.40/ 1.2  |   2.82     3.75    8/8
w0.70 z1.10   1.42/ 0.0  1.34/ 0.8  1.19/ 0.9  0.19/ 1.5A 2.23/ 1.6  2.93/ 4.5  1.91/ 2.5  1.32/ 1.2  |   2.93     4.53    7/8
w0.70 z1.20   1.44/ 0.0  1.23/ 2.6  1.27/ 2.7  0.37/ 4.7A 2.28/ 4.8  3.22/11.8  2.12/ 7.1  1.05/ 3.5  |   3.22    11.84    7/8
w0.75 z1.00   1.42/ 0.0  1.28/ 0.9  1.14/ 0.9  1.16/ 1.8  2.25/ 1.6  2.94/ 4.3  2.13/ 2.3  1.50/ 1.3  |   2.94     4.25    8/8
w0.75 z1.10   1.53/ 0.2  1.46/ 0.5  1.33/ 0.5  0.79/ 1.0  2.21/ 1.0  2.96/ 3.1  2.05/ 1.7  1.44/ 0.8  |   2.96     3.14    8/8
```

**`w=0.70, zeta=1.00` shipped**: it is the lowest-`w` point that lands
**8/8** of the lateral battery (no arrests), passes both gates on the
acceptance case with margin, and keeps nominal clean. Note this is a
*plateau*, not a razor: `w` from 0.6 to 0.75 all pass the acceptance case
with `vtd` in 1.28-1.46 and `miss` in 0.5-1.9 — unlike rounds 2-4's
razor-thin `R` windows (`R=1.4e-9` and `R=1.6e-9` both failed either side
of `1.5e-9`). Selecting on a plateau instead of a knife edge was an
explicit criterion.

### 6. Full final battery (shipped defaults, fresh N=60 solve)

```
tSwitch = 7.160 s  (sol.tf = 15.622 s)
phase-A lateral gains: Kx1/m=0.4808  Kx4/m=1.3773  (was 0.0195 / 0.2334)

case             land  vtd      miss     alt      sat_frac  stop
nominal           1    1.453    0.013  -0.0000   0.969   touchdown
D1 [50 -30 0]     1    1.404    0.980  -0.0000   0.129   touchdown
-D1               1    1.280    0.933  -0.0000   0.106   touchdown
[100 0 0]         1    0.617    1.824  -0.0000   0.139   touchdown
[0 -100 0]        1    2.268    1.611  -0.0000   0.090   touchdown
[141 141 0]       1    2.821    3.885  -0.0000   0.160   touchdown
[-100 -100 0]     1    1.944    2.270  -0.0000   0.102   touchdown
D1 + dv0          1    1.317    1.259  -0.0000   0.165   touchdown
[0 0 +50]         0    1.708    6.820   9.2617   0.483   arrest
[0 0 -50]         1   20.296   15.068  -0.0000   0.275   touchdown
[100 100 50]      0    0.906    7.351   2.6830   0.309   arrest
thrust 0.95       1   67.433   48.824  -0.0000   0.468   touchdown
thrust 1.05       0    0.709    3.649  23.0339   0.436   arrest
isp 0.97          0    0.637    0.044   0.0547   0.449   arrest
D1+thr+isp        1   66.759   54.716  -0.0000   0.477   touchdown
```

Lateral capability envelope: clean to ~100 m of offset; `vtd` starts
crossing 2.0 m/s somewhere past ~100 m and reaches 2.8 m/s at a 200 m
offset. The task's 58 m case sits comfortably inside.

### 7. TWO CAMPAIGN-LEVEL DEFECTS THIS EXPOSED — both PRE-EXISTING, both threaten Task 8

Neither is caused by this round's change, and both are *proved* not to be
by running the identical dispersions through the round-4 controller.

**(i) Altitude dispersions are barely tracked at all.** `dr0=[0;0;-50]`
touches down at **20.3 m/s**; `dr0=[0;0;+50]` arrests 9.3 m up. The
numbers are **identical before and after this round** (20.353/20.296 and
9.240/9.262) — a different axis, untouched. It has the same *shape* of
root cause as the lateral defect just fixed: the margin arc's spare
authority is never asked for. A probe (lateral loop held at the shipped
setting, phase-A z-loop swept) shows it moves a long way:

| phase-A z-loop `wz` | `[0;0;-50]` vtd/miss | `[0;0;+50]` | `[-150;-150;-50]+dv0` |
|---|---|---|---|
| 0.14 (= old default) | 20.35 / 15.18 | arrest 9.24 m | 23.73 / 17.60 |
| 0.30 | 10.95 / 5.78 | arrest 8.97 m | 15.54 / 8.65 |
| 0.45 | **0.85** / 3.14 (arrest 6.1 m) | arrest 6.70 m | **7.02** / 5.77 |
| 0.60 | arrest 13.85 m | arrest 0.72 m | 2.53 / 7.34 |

so a 20.3 m/s impact becomes sub-1 m/s — but it **trades into vertical
arrests**, because arriving *below* the guidance's own descent-rate
schedule is unrecoverable when `Tmin` exceeds weight at every mass (the
vehicle can only ever net-decelerate). Closing that needs a design round
that co-designs the z-loop *with* the terminal `v(z)` schedule and
`P.zTermBand`, not a weight bump — so it was measured, documented, and
deliberately **not** shipped here. **This is the single biggest threat to
Task 8**, whose 1-sigma `r0` draw includes 50 m of altitude.

**(ii) The guidance leaves zero thrust margin, so a thrust bias is
uncorrectable by any controller.** `thrust_scale=0.95` touches down at
**67.4 m/s**; `0.98` at 51.5 m/s. Identical under the round-4 controller
(67.367 / 51.496) — this is not a tracking failure, it is arithmetic.
The min-fuel solution rides `Tmax` for 54% of its nodes (G5 certifies
`bound fraction = 0.9917`), so during braking there is *nothing left to
add*: with 5% less thrust the stopping distance grows from ~736 m to
~790 m while the braking arc only begins at ~718 m altitude, and
`sqrt(180^2 - 2*20.5*718) ~ 54 m/s` of residual speed is unavoidable.
No gain schedule can fix this. The fix is at the **guidance** layer —
solve the nominal against a de-rated `Tmax` (e.g. 0.95*Tmax) so the
tracker has real thrust margin — and it is a spec question for the human,
since it costs propellant. If Task 8's Monte Carlo draws `thrust_scale`
at all, it will fail on this long before it fails on geometry.

### 8. Verification

Full suite, synchronous `-batch`, 600000 ms:
```
/Applications/MATLAB_R2025b.app/bin/matlab -batch "cd('/Users/msc/Desktop/optimal_control/booster_landing'); setup_paths; test_params; test_dynamics_jac; test_colloc_smoke; test_convex_lossless; test_certify_nominal; test_tvlqr_riccati; test_closed_loop_nominal"
```
```
test_params PASS (min-throttle dry T/W = 1.346)
test_dynamics_jac PASS (vacuum + drag)
test_colloc_smoke PASS  tf=15.62 s  mf=26546.9 kg  fuel=3453.1 kg
test_convex_lossless PASS  gap=4.80e-05  mf=26440.8 kg
test_certify_nominal (coarse) PASS
test_certify_nominal (nominal, all_pass) PASS      [G1-G5 + G2ff all PASS]
test_tvlqr_riccati PASS
nominal diagnostics: vtd=1.453 m/s (target 1.5, gate <=2.0)  alt=-0.0000 m  stop=touchdown
dispersed diagnostics: vtd=1.404 m/s  miss=0.980 m  alt=-0.0000 m  stop=touchdown
test_closed_loop_nominal PASS  (nom miss 0.013 m vtd 1.453 m/s, disp miss 0.98 m vtd 1.404 m/s)
```
**7/7 green.** `test_closed_loop_nominal` passes unmodified — no
assertion in it was loosened, added to, or otherwise touched this round;
the only test-facing change is that it now passes.

Certify gates re-run clean after the change (`test_certify_nominal`
above): G1 1.40e-07, G2 pos 4.33e-04 m, G2ff 0/0, G3 |dmf| 0.7037 kg /
|dtf| 0.0059 s, G4 1.12e-04, G5 bound-frac 0.9917 / 1 switch / primer
0.608 deg. Unchanged — as expected, since nothing in the guidance layer
was touched.

### 9. Files changed (round 5)

- `lib/tvlqr_design.m` — phase-scheduled `Q(t)`; new `annulus_switch`
  local function; `ricrhs` takes a `Qfun` handle; new `opts.QA/.QB/
  .tSwitch/.tBlend` (with `opts.Q` kept as the constant-weight legacy
  path); `ctrl` gains diagnostic `.tSwitch/.tBlend/.Qfun`; ADAPTATION 5
  header block with the derivation, the measured before/after, and the
  known-remaining-defect note.
- `lib/sim_closed_loop.m` — new `allocate_thrust` local function
  (vertical-priority upper-bound saturation, direction-preserving
  otherwise); `control_law` calls it; header summary line updated.

No other campaign file touched. `booster_params.m`, both solvers,
`hs_quad_ctrl.m`, `certify/`, and all seven tests are byte-unchanged.

### 10. Self-review (round 5)

- The root cause was **measured before anything was changed** (`diag1`:
  gain table + time history), and the fix was **designed from that
  measurement** (closed-form ARE inversion), not searched for. The one
  time I did search blindly (the first `qA` sweep) it produced a
  misleading partial result, and the fix came from reading the `dx`
  history and recognizing overshoot — recorded in §5.
- Every claim of "pre-existing, not mine" in §7 is backed by running the
  identical dispersion through the round-4 controller and quoting both
  numbers, not by argument.
- The shipped point was chosen for **plateau width**, explicitly reacting
  to round 2/4's razor-thin windows, and the plateau is shown in the §5
  table rather than asserted.
- The `w=0.70` initial command slightly exceeds the annulus (903 kN vs
  845 kN); this is stated in the header rather than glossed, along with
  the fact that `allocate_thrust` is what makes it harmless.
- I did **not** ship the altitude-loop fix even though §7's probe shows a
  20.3 -> 0.85 m/s improvement available, because it trades into arrests
  and I could not close that trade within this round's scope. That is a
  refusal, and the sweep data needed by whoever does close it is in §7.
- Round 4's "genuine actuator-authority limit" conclusion is contradicted
  head-on rather than diplomatically extended; the old ADAPTATION 2-4
  header blocks are left in `tvlqr_design.m` as history with ADAPTATION 5
  explicitly retracting the operative claim.

---

## Round-5 re-review corrections (2026-08-08)

Re-review returned ADDRESSED (derivation independently re-derived, headline
numbers reproduced bit-for-bit) with three cheap corrections plus one
recommended refactor. All four applied. None touch the fix's substance.

### 1. `annulus_switch` no-switch fallback was backwards (`lib/tvlqr_design.m`)

The fallback returned `ts = sol.tf`, which — because the blend is
`0.5*(1+tanh((t-ts)/tBlend))` — puts the *entire* flight on the
**aggressive phase-A lateral weights**, braking arc included. That is
exactly the failure mode the schedule exists to prevent. Reviewer verified
numerically; reproduced here by forging an all-Tmax profile:

```
before fix:  tSwitch=15.622  Qfun(tf/2)(1,1) = 2.10e-01   (phase A -- wrong)
after  fix:  tSwitch= 0.000  Qfun(tf/2)(1,1) = 1.00e-04   (phase B -- right)
```

Fallback is now `ts = 0` (all-QB: nothing to be aggressive with when there
is no margin arc). Also added the requested Phase-2 TODO comment: the
detector finds the *first upward* mid-annulus crossing, which on a
max-min-max profile (reachable once `P.drag.on`) would be the **second**
switch and would wrongly put the leading Tmax arc in phase A.

### 2. Negative commanded `Tz` guard (`lib/sim_closed_loop.m`)

Reviewer's counterexample confirmed: the unguarded vertical-priority
branch turned `Traw=[1e6;1e6;-1e6]` into `T=[0;0;-845 kN]` — the entire
annulus spent thrusting *downward*, lateral correction zeroed. Vertical
priority is only meaningful when the vertical demand is a **braking**
demand; a negative `Traw(3)` ("accelerate the descent", reachable whenever
the vehicle is above the `v(z)` schedule — a regime Task 8's MC will draw
routinely) has no braking authority to protect. The priority branch is now
taken only for `Traw(3) >= 0`; otherwise plain direction-preserving
scaling. Unit checks (`Tmax = 845 kN`):

```
Traw=[    1e6     1e6    -1e6] -> T=[ 487.9  487.9 -487.9] kN, |T|=845.0   (was [0;0;-845])
Traw=[    1e6     1e6     1e6] -> T=[   0.0    0.0  845.0] kN, |T|=845.0
Traw=[    100       0    -1e6] -> T=[   0.1    0.0 -845.0] kN, |T|=845.0
Traw=[    3e5       0     3e5] -> T=[ 300.0    0.0  300.0] kN, |T|=424.3   (under bound, untouched)
Traw=[      0       0     1e3] -> T=[   0.0    0.0  338.0] kN, |T|=338.0   (Tmin floor)
```

The dead `if Tm < P.Tmin` restore inside that branch is deleted. Chosing
the direction-preserving fallback (rather than clamping `Tz` to `[0,Tmax]`)
**keeps it provably dead**, which the clamp option would not have: `|T| ==
Tmax` identically on both over-budget branches. Proof written into the
header — priority branch: `Traw(3)<=Tmax` gives `Tz=Traw(3)` and
`|Traw|>Tmax` forces `|Txy_raw|>lat`, so `Txy` is always scaled and
`|T|=sqrt(lat^2+Tz^2)=Tmax`; `Traw(3)>Tmax` gives `lat=0`, `|T|=Tmax`;
fallback branch: scaling a vector with `|Traw|>Tmax` lands on `Tmax`. The
unit checks above confirm all five cases.

### 3. Two header claims corrected

**(a) The R-independence claim was sloppy, and the truth is stronger.**
I wrote that the pole is "independent of R to leading order" because
"qVel dominates". Wrong: at the shipped `r=1e-9` the dropped
`2 m sqrt(qPos r) = 0.0187` term is **~1.9x LARGER** than `qVel = 0.01`.
The correct statement, now in the header, is the exact one:

```
pole = a/b = sqrt( qPos / (2 m sqrt(qPos r) + qVel) )
```

which is **strictly decreasing in `r`**, with supremum `sqrt(qPos/qVel)`
attained only as `r -> 0`. So the old weights' 0.1 rad/s is a **ceiling no
value of R can exceed** — a stronger claim than "R-independent", and it
closes the rounds 2-4 search question completely rather than approximately.
(At the shipped `r` the actual pole is 0.059 rad/s, further below the
ceiling; the 0.084 measured at `t=0` is the Riccati still off steady state.)

**(b) "allocate_thrust is a no-op over the whole acceptance case" was
false.** Measured, time-weighted over the flight:

| case | raw \|T\| over Tmax | peak raw \|T\| | over-bound by |
|---|---|---|---|
| nominal | 2.1% of flight time | 850.4 kN | +0.6% |
| dispersed 58 m | 1.9% of flight time | 968.3 kN | **+14.6%** |

reproducing the reviewer's 2.1%/1.9% figures exactly. Header corrected to
say this. Two notes I added on top: the *nominal* excursions are a
knife-edge artifact (the guidance rides exactly on Tmax, so microscopic
tracking error flips "over" on and off — a sample-count read of the same
data gave 46%, which is why the number is quoted time-weighted and paired
with the +0.6% magnitude); the dispersed +14.6% is the real one, and is
the entire reason the branch's policy matters.

### 4. `QA` now derived from `(omegaA, zetaA, m, R)`, not hardcoded

Taken up — it stayed small (~12 lines). `QA` was hardcoded against
`R=1e-9`, so an `opts.R` override silently detuned phase A (`omega` scales
as `r^-1/4`, so a 100x R change moves the pole 3.2x). Now `qPos = m^2 r w^4`,
`qVel = m^2 r w^2 (4 zeta^2 - 2)` are evaluated at call time from the
**actual** `R(1,1)` and the mean mass over the margin arc, with
`opts.omegaA`/`opts.zetaA` (0.70, 1.00) exposing the design point and
`opts.QA` still overriding everything. Ordering in the function was
rearranged so `R` and `tSwitch` are known before `QA` is built; an
`~any(inA)` guard covers the new `ts=0` fallback.

Derived values vs the hardcoded ones: `qPos 0.2100 -> 0.2099` (-0.05%),
`qVel 0.8500 -> 0.8569` (+0.8%).

### 5. Verification: items 1-3 are EXACTLY behaviour-neutral

Run twice against the same fresh N=60 solve — once with `opts.QA` pinned
to the committed `0.21/0.85` (isolating items 1-3), once with the derived
default (adding item 4). Deltas are against the committed `5f0b8dc`
baseline:

```
--- items 1+2+3 ONLY (QA pinned to committed 0.21/0.85) ---
  nominal    land=1 miss=0.013 (d+0.000) vtd=1.453 (d-0.000) | overTmax 2.1% of time, peak 850.4 kN (+0.6%), negTz-while-over 0.00%
  dispersed  land=1 miss=0.980 (d-0.000) vtd=1.404 (d-0.000) | overTmax 1.9% of time, peak 968.3 kN (+14.6%), negTz-while-over 0.00%

--- + item 4 (QA derived from omega/zeta/m/R) ---
    derived qPos=0.2099 qVel=0.8569 (was 0.2100 / 0.8500)
  nominal    land=1 miss=0.052 (d+0.039) vtd=1.480 (d+0.027) | overTmax 47.6% of time, peak 850.5 kN (+0.7%), negTz-while-over 0.00%
  dispersed  land=1 miss=0.944 (d-0.036) vtd=1.363 (d-0.041) | overTmax 1.9% of time, peak 968.1 kN (+14.6%), negTz-while-over 0.00%
```

**Items 1-3: zero change to three decimals on all four acceptance numbers**,
as required. The `negTz-while-over 0.00%` column confirms the reviewer's
read that both new paths are *latent* for current profiles — the guard
never fires on the acceptance cases, it is insurance for Task 8. (The
nominal `overTmax 47.6%` in the second block is the knife-edge artifact
described in 3(b), not a real change: the peak moved 850.4 -> 850.5 kN.)

**Item 4** perturbs by <0.05 m and <0.05 m/s, in the *improving* direction
on the dispersed case (miss -0.036, vtd -0.041) and negligibly the other
way on nominal (miss 0.013 -> 0.052 m against a 1 m gate). Kept: the
robustness win for downstream `opts.R` users is worth more than 39 mm.

### 6. Full suite

```
/Applications/MATLAB_R2025b.app/bin/matlab -batch "cd('/Users/msc/Desktop/optimal_control/booster_landing'); setup_paths; test_params; test_dynamics_jac; test_colloc_smoke; test_convex_lossless; test_certify_nominal; test_tvlqr_riccati; test_closed_loop_nominal"
```
```
test_params PASS (min-throttle dry T/W = 1.346)
test_dynamics_jac PASS (vacuum + drag)
test_colloc_smoke PASS  tf=15.62 s  mf=26546.9 kg  fuel=3453.1 kg
test_convex_lossless PASS  gap=4.80e-05  mf=26440.8 kg
test_certify_nominal (coarse) PASS
test_certify_nominal (nominal, all_pass) PASS
test_tvlqr_riccati PASS
nominal diagnostics: vtd=1.480 m/s (target 1.5, gate <=2.0)  alt=-0.0000 m  stop=touchdown
dispersed diagnostics: vtd=1.363 m/s  miss=0.944 m  alt=-0.0000 m  stop=touchdown
test_closed_loop_nominal PASS  (nom miss 0.052 m vtd 1.480 m/s, disp miss 0.94 m vtd 1.363 m/s)
```
**7/7 green**, both covering tests included, no test file modified.

### 7. Files changed (re-review)

- `lib/tvlqr_design.m` — `annulus_switch` fallback `sol.tf` -> `0` + Phase-2
  TODO; `QA` derived from `(omegaA, zetaA, mRef, R)` with defaults block
  reordered; ADAPTATION 5 pole claim corrected to the exact monotone-in-R
  bound; shipped qPos/qVel and before/after numbers refreshed.
- `lib/sim_closed_loop.m` — `allocate_thrust` negative-`Tz` guard, dead
  `Tmin` restore deleted with the `|T| == Tmax` proof documented, "no-op"
  claim corrected with measured firing fractions.

The two campaign-level defects flagged in round 5 §7 (altitude dispersions
tracked far too weakly; guidance riding Tmax with zero margin so a 5%
thrust bias is uncorrectable) are **unchanged and still open** — nothing in
this re-review touched either, and both still threaten Task 8.

---

# Task 7b (2026-08-08) — guidance de-rate + vertical-loop co-design

**Status: 7b-2 DONE (altitude dispersions solved). 7b-1 implemented but
BLOCKED at the adjudicated value: `P.etaT=0.93` is measurably too weak to
make a -5% thrust dispersion pass. 5 of 7 battery cases pass; the remedy is
one number and is quantified in §6.**

## 1. 7b-1: the de-rate as implemented

`P.etaT = 0.93` added to `booster_params.m`. The GUIDANCE upper thrust
bound becomes `etaT*Tmax = 785.85 kN` (note: the brief said 786.85; 0.93 ×
845 = 785.85, and 0.93 is what is implemented). `P.Tmin` is untouched --
it is an engine floor, not a budget. The tracker and truth sim keep the
FULL `[Tmin, Tmax]` annulus.

De-rated call sites, all audited rather than pattern-matched:

| file | what | now uses |
|---|---|---|
| `solve_pdg_colloc` | annulus upper bound `Tmax_h` | `etaT*Tmax` |
| `solve_pdg_convex` | `z0`, `zlb` depletion reference | `TmaxG` |
| `solve_pdg_convex` | Taylor upper bound `mu2` | `TmaxG` |
| `solve_pdg_convex` | initial guess `Sg` | `TmaxG` |
| `solve_pdg_convex` | mass-depletion feasibility assert | `TmaxG` |
| `solve_pdg_convex` | `tightTol` classifier scale | `etaT*Tmax` |
| `certify_pdg` | G2 reconstruction clamp | `TmaxG` |
| `certify_pdg` | G2ff over-ceiling check | `TmaxG` |
| `certify_pdg` | G4 tightness scale | `TmaxG` |
| `certify_pdg` | G5 "on the upper bound" detector | `TmaxG` |
| `tvlqr_design` | `ctrl.Tnom` reconstruction clamp | `etaT*Tmax` |
| `test_colloc_smoke`, `test_convex_lossless` | annulus assertions | `etaT*Tmax` |
| `sim_closed_loop` `allocate_thrust` | tracker clamp | **full `Tmax`** (unchanged) |

`P.tf_hi = 22` still brackets the de-rated optimum (`tf = 16.0955 < 22`),
so it was left alone -- checked, not assumed.

### Re-solve, both solvers, nominal grids

| | before (etaT=1) | after (etaT=0.93) |
|---|---|---|
| colloc `tf` | 15.6225 s | **16.0955 s** |
| colloc `mf` | 26546.892 kg | **26506.925 kg** |
| colloc fuel | 3453.108 kg | **3493.075 kg** (**+39.97 kg, +1.16%**) |
| convex `tf` | 15.6284 s | 16.1020 s |
| convex `mf` | 26546.189 kg | 26506.365 kg |
| `|T*|` range | 338.0 - 845.0 kN | 338.0 - **785.8** kN |

**Fuel cost of the de-rate: +39.97 kg (+1.16% of propellant).**

### Gate table, old vs new (nominal grid)

| gate | before | after | threshold | verdict |
|---|---|---|---|---|
| G1 max HS defect | 1.40e-07 | 1.736e-07 | < 1e-6 | PASS |
| G2 pos residual | 4.33e-04 m | **1.395e-02 m** | < 1 | PASS |
| G2 vel residual | 5.42e-05 | 8.67e-04 | < 0.1 | PASS |
| G2 mass residual | 1.30e-05 kg | 5.12e-06 kg | < 0.5 | PASS |
| G2ff below/above | 0 / 0 | 1.7e-16 / 1.5e-16 | < 1e-6 | PASS |
| G3 \|dmf\| | 0.7037 kg | **0.5601 kg** | < 1 | PASS |
| G3 \|dtf\| | 0.00586 s | 0.00647 s | < 0.2 | PASS |
| G4 lossless gap | 1.118e-04 | 1.400e-04 | < 1e-4*TmaxG/m0 | PASS |
| G5 bound fraction | 0.9917 | 0.9917 | >= 0.95 | PASS |
| G5 switches / primer | 1 / 0.608 deg | 1 / 0.597 deg | <=2 / <1 | PASS |
| **all_pass** | PASS | **PASS** | | |

`test_certify_nominal`'s three numeric pins (`G2_dm < 0.01`,
`G3_dmf < 1.0`, `G3_dtf < 0.2`) all still hold and were left unchanged;
the measured values moved to 5.12e-06 kg, 0.5601 kg and 0.00647 s.
`results/pdg_colloc_nominal.mat` and `results/pdg_convex_nominal.mat`
regenerated.

## 2. A pre-existing G2 defect the de-rate exposed (and a fix for it)

The first de-rated certify run **FAILED G2**: pos residual 1.287 m against
a 1 m gate (vel 0.140, mass 1.343 kg). It is not caused by the de-rate. A
(etaT, N) sweep of the *unmodified* reconstruction:

| etaT \ N | 60 | 80 | 120 |
|---|---|---|---|
| 1.00 (shipped) | 4.33e-04 | **0.9339** | 2.19e-03 |
| 0.97 | 7.33e-05 | 2.03e-04 | 0.2474 |
| 0.93 | **1.287** | **1.241** | 0.2272 |

Same physics, `N=60` vs `N=80`, G2 moves by 2000x. **The shipped
(etaT=1, N=60) pass was a grid lottery ticket**, not a property.

Root cause, localized by integrating the mass error in time: the *entire*
residual is generated inside ONE segment. Mass error is 0.000 kg at t=5 s,
+1.453 kg at t=6.5 s, -1.343 kg from t=7.0 s onward, then flat to
touchdown. That segment is the bang-bang switch. `hs_quad_ctrl`
reconstructed `|T|` as a quadratic through the three samples -- and a
quadratic through a STEP's samples is the wrong basis function, so the
reconstruction burns a different amount of propellant across the switch
than the NLP solved for.

**Fix** (`lib/hs_quad_ctrl.m`): on a transition segment (samples not all
pinned to one bound) `|T|` is reconstructed as a STEP between the certified
bang-bang levels `Tmin`/`Tmax`, placed at the instant `s` that makes the
step's integral equal the segment's own Simpson quadrature:
`s*A + (h-s)*B = (h/6)(m0 + 4mm + m1)`. Direction is untouched (the primer
direction is continuous through a switch; only the magnitude steps).
Feasibility is preserved by construction -- both levels are annulus
samples -- so G2ff cannot regress, and it does not (1.5e-16). Falls back to
the quadratic if `s` leaves `[0,h]`, which is what would happen on a
genuine singular arc.

Detector note: the first version tested only "endpoints on opposite
bounds" and missed the case where the NLP's off-bound sample lands on a
NODE (the switch then straddles two segments) -- which was exactly the
etaT=0.93/N=60 and etaT=1.00/N=80 cases. Generalized to "not all three
samples on one bound".

Same sweep after the fix -- **the lottery is gone**:

| etaT \ N | 60 | 80 | 120 |
|---|---|---|---|
| 1.00 | 0.0371 | 0.0082 | 0.0077 |
| 0.97 | 0.0511 | 0.0222 | 0.0023 |
| 0.93 | **0.0140** | 0.0092 | 0.0027 |

Every combination now lands in 0.002-0.051 m, including the two that
failed or nearly failed *before* this round's work. The honest G2 accuracy
of the reconstruction is ~1e-2 m; the old 4.3e-4 was a coincidence.

## 3. 7b-2: the vertical-loop co-design

Starting point after the de-rate alone: **0 of 7** battery cases passed --
the gentler de-rated terminal made even the nominal case arrest 2.3 cm
above the pad.

### 3a. The real defect was the SCHEDULE VARIABLE, not the gain

Instrumenting `dr0=[0;0;-50]` with a *well-tuned* time-indexed vertical
loop showed the vehicle tracking `vznom(t)` essentially perfectly
(`vz=-169.38` vs `-169.84` at t=6.76 s) -- and that is exactly the
failure. Tracking time perfectly means staying 50 m low all the way down
and **hitting the ground 50 m early at -31.8 m/s** while the guidance still
had 2 s of braking to do. Tightening the loop made it WORSE:

| vertical loop bandwidth | dr0z=-50 touchdown |
|---|---|
| 0.12 rad/s (pre-7b) | 9.14 m/s |
| 1 rad/s | 39.16 m/s |
| 5 rad/s | 43.41 m/s |

The better it tracked the wrong schedule, the harder it drove into the
ground. The round-4/5 code only escaped this by having a vertical loop so
weak that velocity SAGGED toward the altitude-appropriate profile on its
own -- accidentally doing altitude indexing, badly, in the last 150 m.

**Fix: index the reference by ALTITUDE, not wall-clock time.** At altitude
`z` the tracker flies toward `(r*_xy(z), v*(z))` with feedforward
`T*(t*(z))` and gains `K(t*(z))`. A min-fuel trajectory ends exactly at
`z=0` with the tanks near dry -- it has no time margin, but the invariant
the fuel-optimal solution actually encodes is *where in altitude to brake*.
Under altitude indexing `dr0=[0;0;-50]` is nearly a non-event: at z=1950
the reference is `v*(1950) = -179.83` against the vehicle's `-180.00`.

Three special cases DELETED rather than added:
* the altitude-error channel vanishes identically (`xref(3) := x(3)`), so
  round 4's explicit `Kt(:,3)=0` terminal zeroing is now true by
  construction, everywhere, instead of below a hand-picked band;
* the "past-tf freeze" cannot occur -- `K` indexed by `t*(z)` reaches its
  terminal value only as the vehicle reaches the ground;
* `P.zTermBand` and its hand-picked 150 m are no longer used.

Monotone descent is asserted, not assumed.

Result: `dr0=[0;0;-50]` went **43.4 -> 1.005 m/s, landed**.

### 3b. The vertical velocity loop, derived not tuned

With the altitude channel gone the vertical axis is a pure first-order
velocity loop `vzdot = Tz/m - g`, whose scalar ARE gives `gain/m =
sqrt(qVelZ/r)/m`, so a target bandwidth inverts to
`qVelZ = (m*omegaVz)^2 * r` -- the same machinery as round 5's lateral
loop. Shipped `omegaVz = 12 rad/s`, against a measured pre-7b value of
**0.12 rad/s**. The terminal descent-rate reference sweeps at ~19.8 m/s^2,
so a loop of bandwidth `w` carries a standing lag of `19.8/w` m/s; at
0.12 rad/s the vehicle was never tracking the last 150 m at all.

### 3c. The Riccati boundary layer (the subtle one)

Raising `omegaVz` alone still left arrests and a 2.8 m/s thrust-case
error. Cause: `Qf(6,6)=1` gives a terminal vertical gain of
`sqrt(1/r)/m = 1.19 1/s`. Because altitude indexing evaluates `K` at
`t*(z) -> tf` as `z -> 0`, **the vehicle flew the last metre on a 1.19 1/s
loop no matter how the rest was tuned** -- the design bandwidth never
reached the ground. Setting `Qf(6,6) := qVelZ` makes `P(tf)` consistent
with the phase-B steady state so the gain is smooth through touchdown.
That single change removed every remaining arrest in the battery.

`ctrl.Qf` is now returned so `test_tvlqr_riccati` asserts `P(tf)=Qf`
against the weight actually used instead of a hardcoded literal.

### 3d. Two things tried and REJECTED (reported, not hidden)

* **Lateral fade-out near the ground** (taper lateral feedback below 8 m so
  the whole vector goes vertical): made things *worse* (thrust 0.95:
  2.86 -> 3.05; thrust 1.05: 1.66 -> 1.91). `vtd` is `norm(v)`, so
  un-damped lateral velocity costs more than the recovered tilt. Reverted.
* **Lateral velocity rescale** (`v_xy_ref *= vz/vz*`, since altitude
  indexing pins the slope `dr_xy/dz`, not `v_xy`): theoretically cleaner
  and it did improve the altitude cases (z+50 miss 2.32 -> 2.01), but was
  neutral-to-worse elsewhere (thrust 0.95: 2.842 -> 2.874). Not worth the
  complexity. Reverted.

## 4. Battery (shipped defaults, etaT=0.93, omegaVz=12)

| case | landed | miss (m) | vtd (m/s) | gate |
|---|---|---|---|---|
| nominal | yes | 0.072 | 1.436 | PASS |
| `dr0=[50;-30;0]` | yes | 1.292 | 1.432 | PASS |
| `dr0=[0;0;-50]` | yes | 0.465 | 1.439 | **PASS** (was 20.3 m/s) |
| `dr0=[0;0;+50]` | yes | 2.321 | 1.431 | **PASS** (was arrest 9.2 m) |
| `thrust_scale=0.95` | yes | 3.707 | **2.842** | **FAIL** (was 67.4 m/s) |
| `thrust_scale=1.05` | yes | 3.599 | 1.664 | PASS |
| combined (`[50;-30;-50]`, 0.97) | yes | 3.938 | **2.068** | **FAIL** |

**5/7.** Every case now genuinely lands (no arrests anywhere), and the two
adjudicated altitude cases are solved outright.

## 5. Why the two thrust cases are BLOCKED, with proof

Not a tuning failure. A reachability bound, integrating the vertical
dynamics in altitude under the *optimal* single-switch policy (coast at
Tmin, then full Tmax, brake altitude found by bisection, zero tracking
error, no lateral demand -- optimistic for any controller):

| thrust_scale | min achievable \|vtd\| | required brake altitude |
|---|---|---|
| 1.00 | 1.150 m/s | 699.6 m |
| 0.98 | 0.104 m/s | 741.1 m |
| 0.95 | **0.635 m/s** | **806.6 m** |
| 0.90 | 0.755 m/s | 925.2 m |

So a -5% vehicle **can** land softly -- but only by starting its brake at
**806.6 m instead of the nominal 718 m**. A fixed-reference tracker
inherits the nominal brake altitude from the feedforward and structurally
cannot re-time it. Everything else it can do, it does: the loop holds the
profile to within 0.2 m/s from z=400 m down to z=2 m, then the reference
steepens (`|dv*/dz| = a*/|v|` blows up as `|v| -> 1.5`) beyond what a -5%
vehicle can follow, and the error opens to 0.74 m/s in the last metre.

Raising the loop bandwidth cannot close it -- the floor is an authority
limit, not a lag:

| omegaVz | 5 | 8 | 12 | 20 |
|---|---|---|---|---|
| thrust 0.95 vtd | 3.91 | 3.36 | 2.86 | 2.66 |

## 6. The remedy is one number, and it is measured

The de-rate's own parameter is the lever. Holding the controller fixed and
sweeping `P.etaT` (full re-solve each time):

| etaT | tf (s) | fuel (kg) | a*_terminal | thrust0.95 vtd | battery |
|---|---|---|---|---|---|
| **0.93 (adjudicated)** | 16.095 | 3493.07 | 19.84 | **2.86** | **5/7** |
| 0.91 | 16.251 | 3506.24 | 19.22 | 2.41 | 6/7 |
| 0.89 | 16.417 | 3520.31 | 18.59 | 2.08 | 6/7 |
| **0.87** | 16.595 | **3535.41** | 17.97 | **1.66** | **7/7** |
| 0.85 | 16.786 | 3551.65 | 17.35 | 1.50 | 7/7 |

**`etaT = 0.87` passes all seven with margin, for +42.3 kg of propellant
over etaT=0.93 (+82.3 kg, 2.4%, over no de-rate at all).**

The adjudication set 0.93 on the stated rationale that it "reserves
3-sigma+ headroom for feedback". That estimate was made before anyone
measured what a *fixed-reference* tracker needs. The measurement says
0.93 leaves the tracker only `0.95/0.93 = 2.1%` of thrust margin against a
-5% dispersion, and the terminal arc consumes more than that. **I did not
change the adjudicated value** -- that is the human's call, and it costs
propellant. Recommend re-adjudicating `P.etaT` to 0.87.

The alternative (cheaper in propellant, more expensive in engineering) is
to let the guidance be re-planned or the brake altitude re-timed on
measured performance -- explicitly out of scope for the TVLQR family. An
in-family research direction worth recording: `|vz|` is also monotone on
this trajectory, so indexing the feedforward by *speed* rather than
altitude would re-time the switch automatically for an underperforming
vehicle. Not attempted -- it is a second schedule-variable change and
would need its own verification round.

## 7. Test status

```
/Applications/MATLAB_R2025b.app/bin/matlab -batch "cd('/Users/msc/Desktop/optimal_control/booster_landing'); setup_paths; test_params; test_dynamics_jac; test_colloc_smoke; test_convex_lossless; test_certify_nominal; test_tvlqr_riccati; test_closed_loop_nominal"
```
```
test_params PASS (min-throttle dry T/W = 1.346)
test_dynamics_jac PASS (vacuum + drag)
test_colloc_smoke PASS  tf=16.09 s  mf=26506.9 kg  fuel=3493.1 kg
test_convex_lossless PASS  gap=4.80e-05  mf=26438.6 kg
test_certify_nominal (coarse) PASS
test_certify_nominal (nominal, all_pass) PASS
test_tvlqr_riccati PASS
nominal diagnostics: vtd=1.436 m/s  alt=-0.0000 m  stop=touchdown
dispersed diagnostics: vtd=1.432 m/s  miss=1.292 m  alt=-0.0000 m  stop=touchdown
  battery dr0=[0;0;-50]        land=1 miss=  0.465 m vtd= 1.439 m/s  PASS
  battery dr0=[0;0;+50]        land=1 miss=  2.321 m vtd= 1.431 m/s  PASS
  battery thrust_scale=0.95    land=1 miss=  3.707 m vtd= 2.842 m/s  **FAIL**
  battery thrust_scale=1.05    land=1 miss=  3.599 m vtd= 1.664 m/s  PASS
  battery combined             land=1 miss=  3.938 m vtd= 2.068 m/s  **FAIL**
Error: 2 of 5 task-7b battery cases failed
```

**6 of 7 tests green.** `test_closed_loop_nominal` fails specifically and
only on the two thrust-dispersion battery cases. The battery was NOT
weakened to make it pass, and `P.etaT` was NOT changed from the
adjudicated value to make it pass -- per the round's own STOP instruction.

## 8. Files changed (7b)

- `lib/booster_params.m` — `P.etaT = 0.93` with the adjudication rationale.
- `lib/solve_pdg_colloc.m`, `lib/solve_pdg_convex.m` — de-rated guidance ceiling.
- `certify/certify_pdg.m` — `TmaxG` for G2/G2ff/G4/G5; `rep.TmaxG`.
- `lib/hs_quad_ctrl.m` — switch-segment step reconstruction (§2).
- `lib/tvlqr_design.m` — `omegaVz` vertical loop derived from the ARE;
  `Qf(6,6)` boundary-layer fix; `ctrl.Qf` exposed; de-rated `Tnom` clamp.
- `lib/sim_closed_loop.m` — altitude-indexed reference and gains.
- `tests/test_colloc_smoke.m`, `tests/test_convex_lossless.m` — annulus vs `etaT*Tmax`.
- `tests/test_tvlqr_riccati.m` — asserts `P(tf)=ctrl.Qf`, no literal.
- `tests/test_closed_loop_nominal.m` — 7b acceptance battery block.
- `results/*.mat` — regenerated (gitignored).

## 9. Self-review (7b)

- The G2 failure was NOT patched around: it was localized to one segment by
  measurement, proven pre-existing by reproducing it at etaT=1/N=80, and
  fixed at the source. The fix improved two configurations that were
  already broken before this round touched anything.
- The altitude-indexing change DELETED three special cases (`Kt(:,3)=0`,
  the past-tf gain freeze, `P.zTermBand`) rather than adding a fourth.
- Both rejected experiments (§3d) are reported with their numbers.
- The blocked cases are backed by a reachability bound and a bandwidth
  sweep showing an authority floor, not by "I tried and it didn't work".
- The adjudicated `P.etaT` was left at 0.93 even though 0.87 would turn the
  suite green, because changing an adjudicated number to pass my own test
  is exactly the failure mode this campaign has been burned by.

---

## §7b FINAL (2026-08-09) — P.etaT = 0.87 adjudicated and shipped

The human re-adjudicated `P.etaT` to **0.87**, the measured 7/7 point from
§6 above. This section records the shipped state. **Task 7b is complete:
the full battery passes 7/7 and all seven tests are green.**

### 1. Parameter

`P.etaT = 0.87` in `lib/booster_params.m` (guidance ceiling
`etaT*Tmax = 735.15 kN`; tracker and truth sim keep the full
`[Tmin, Tmax]`). The comment now carries the whole adjudication history --
that 0.93 was tried first and shipped briefly, that it was measured
insufficient (5/7; only `0.95/0.93 = 2.1%` of thrust margin against a -5%
dispersion), that the shortfall was proven to be an authority limit by the
reachability bound and the bandwidth sweep rather than a tuning failure,
and that 0.87 is the measured 7/7 point. No other parameter changed.

### 2. Re-solve + re-certify at nominal grids

| | etaT=1.00 (baseline) | etaT=0.93 | **etaT=0.87 (shipped)** |
|---|---|---|---|
| colloc `tf` | 15.6225 s | 16.0955 s | **16.5952 s** |
| colloc `mf` | 26546.892 kg | 26506.925 kg | **26464.589 kg** |
| colloc fuel | 3453.108 kg | 3493.075 kg | **3535.411 kg** |
| convex `tf` | 15.6284 s | 16.1020 s | 16.5986 s |
| convex `mf` | 26546.189 kg | 26506.365 kg | 26464.161 kg |
| convex fuel | 3453.811 kg | 3493.635 kg | 3535.839 kg |
| `\|T*\|` range | 338.0-845.0 kN | 338.0-785.8 kN | **338.0-735.1 kN** |

`tf = 16.5952 s` is still comfortably inside `P.tf_hi = 22` -- checked, not
assumed. The two solvers agree to `dtf = 0.0034 s` and `dmf = 0.428 kg`.

**Gate table (nominal grid), old vs new:**

| gate | etaT=0.93 | **etaT=0.87** | threshold | verdict |
|---|---|---|---|---|
| G1 max HS defect | 1.736e-07 | 1.392e-07 | < 1e-6 | PASS |
| G2 pos residual | 1.395e-02 m | **8.847e-03 m** | < 1 | PASS |
| G2 vel residual | 8.67e-04 | 2.989e-04 | < 0.1 | PASS |
| G2 mass residual | 5.12e-06 kg | 8.237e-06 kg | < 0.5 | PASS |
| G2ff below / above Tmin/Tmax | 1.7e-16 / 1.5e-16 | 1.7e-16 / 1.6e-16 | < 1e-6 | PASS |
| G3 \|dmf\| | 0.5601 kg | **0.4279 kg** | < 1 | PASS |
| G3 \|dtf\| | 0.00647 s | **0.00341 s** | < 0.2 | PASS |
| G3 traj Linf | 0.3347 m | 0.3726 m | (info) | -- |
| G4 lossless gap | 1.400e-04 | 1.356e-04 | < 1e-4*TmaxG/m0 | PASS |
| G5 bound fraction | 0.9917 | 0.9917 | >= 0.95 | PASS |
| G5 interior switches | 1 | 1 | <= 2 | PASS |
| G5 structure | PASS | PASS | -- | PASS |
| G5 primer angle | 0.597 deg | 0.574 deg | < 1 | PASS |
| **all_pass** | PASS | **PASS** | | |

**Pins:** the three numeric regression pins in `test_certify_nominal`
(`G2_dm < 0.01`, `G3_dmf < 1.0`, `G3_dtf < 0.2`) all still hold and needed
no update; the measured values improved to 8.24e-06 kg, 0.4279 kg and
0.00341 s (G3 `|dmf|` is now 57% below its ceiling, the most margin it has
had). `results/pdg_colloc_nominal.mat` and `results/pdg_convex_nominal.mat`
regenerated from the shipped solutions.

### 3. FINAL BATTERY — 7/7

| case | landed | miss (m) | gate 15 | vtd (m/s) | gate 2.0 |
|---|---|---|---|---|---|
| nominal | yes | 0.004 | PASS (also <1) | 1.453 | PASS |
| `dr0=[50;-30;0]` | yes | 1.548 | PASS | 1.451 | PASS |
| `dr0=[0;0;-50]` | yes | 0.388 | PASS | 1.452 | PASS |
| `dr0=[0;0;+50]` | yes | 1.980 | PASS | 1.466 | PASS |
| `thrust_scale=0.95` | yes | 3.188 | PASS | 1.647 | PASS |
| `thrust_scale=1.05` | yes | 3.038 | PASS | 1.639 | PASS |
| combined `[50;-30;-50]`, 0.97 | yes | 4.060 | PASS | 1.454 | PASS |

**Worst case across the battery: miss 4.060 m against a 15 m pad radius
(73% margin), vtd 1.647 m/s against a 2.0 m/s gate (18% margin).** Every
case is a genuine `z=0` touchdown -- no arrests, no horizon expiries. No
assertion in `test_closed_loop_nominal.m` was relaxed; the battery block
passes exactly as written in the previous round.

Journey of the three cases that drove this work:

| case | round 4 | 7b @ etaT=0.93 | **7b FINAL @ 0.87** |
|---|---|---|---|
| `dr0=[50;-30;0]` | 6.82 m/s | 1.432 | **1.451** |
| `dr0=[0;0;-50]` | 20.30 m/s | 1.439 | **1.452** |
| `thrust_scale=0.95` | 67.43 m/s | 2.842 (FAIL) | **1.647** |

### 4. Fuel cost of the de-rate

| configuration | fuel used | vs baseline | % of usable propellant |
|---|---|---|---|
| no de-rate (etaT=1.00) | 3453.108 kg | -- | -- |
| etaT=0.93 (first adjudication) | 3493.075 kg | +39.97 kg | +0.91% |
| **etaT=0.87 (shipped)** | **3535.411 kg** | **+82.30 kg** | **+1.87%** |

Usable propellant is `m0 - mdry = 4400 kg`, so the total de-rate costs
**+82.30 kg = 1.87% of the propellant budget** (2.38% of the fuel actually
burned). Touchdown mass 26464.6 kg leaves 864.6 kg of margin above
`P.mdry = 25600 kg`. That is the price of taking `thrust_scale=0.95` from
a 67.4 m/s crash to a 1.65 m/s landing.

### 5. Full suite

```
/Applications/MATLAB_R2025b.app/bin/matlab -batch "cd('/Users/msc/Desktop/optimal_control/booster_landing'); setup_paths; test_params; test_dynamics_jac; test_colloc_smoke; test_convex_lossless; test_certify_nominal; test_tvlqr_riccati; test_closed_loop_nominal"
```
```
test_params PASS (min-throttle dry T/W = 1.346)
test_dynamics_jac PASS (vacuum + drag)
test_colloc_smoke PASS  tf=16.60 s  mf=26464.3 kg  fuel=3535.7 kg
test_convex_lossless PASS  gap=4.81e-05  mf=26435.3 kg
test_certify_nominal (coarse) PASS
test_certify_nominal (nominal, all_pass) PASS
test_tvlqr_riccati PASS
nominal diagnostics: vtd=1.453 m/s (target 1.5, gate <=2.0)  alt=-0.0000 m  stop=touchdown
dispersed diagnostics: vtd=1.451 m/s  miss=1.548 m  alt=-0.0000 m  stop=touchdown
  battery dr0=[0;0;-50]        land=1 miss=  0.388 m vtd= 1.452 m/s  PASS
  battery dr0=[0;0;+50]        land=1 miss=  1.980 m vtd= 1.466 m/s  PASS
  battery thrust_scale=0.95    land=1 miss=  3.188 m vtd= 1.647 m/s  PASS
  battery thrust_scale=1.05    land=1 miss=  3.038 m vtd= 1.639 m/s  PASS
  battery combined             land=1 miss=  4.060 m vtd= 1.454 m/s  PASS
test_closed_loop_nominal PASS  (nom miss 0.004 m vtd 1.453 m/s, disp miss 1.55 m vtd 1.451 m/s)
```
**7 of 7 tests green.**

### 6. Files changed (7b final)

- `lib/booster_params.m` — `P.etaT` 0.93 -> 0.87, comment rewritten with the
  full adjudication history and the measured cost.
- `results/pdg_colloc_nominal.mat`, `results/pdg_convex_nominal.mat` —
  regenerated (gitignored, not committed).

No code changed. Every mechanism shipped in the previous 7b commit
(de-rate plumbing, switch-segment reconstruction, altitude-indexed
tracking, ARE-derived vertical loop, `Qf(6,6)` boundary-layer fix) is
unchanged and simply re-verified at the new parameter.

### 7. Standing note for task 8

Both defects flagged at the end of round 5 are now CLOSED: altitude
dispersions (20.3 m/s -> 1.45) by altitude-indexed tracking, and thrust
dispersions (67.4 m/s -> 1.65) by the de-rate. The Monte Carlo's 1-sigma
`r0 = [100;100;50]` draw is larger than anything in this battery on the
lateral axes, so task 8 should expect some cases outside it; the measured
lateral capability envelope (round 5 §6) is clean to ~100 m of offset with
`vtd` crossing 2.0 m/s somewhere past that, and the battery's worst case
here uses only 27% of the pad radius. Nothing further is flagged as a known
blocker.

---

## §7b polish (2026-08-09) — closing-review fixes

Four fixes from the closing review. **I1 chose branch (a): the
ARE-consistent terminal weight holds 7/7 with BETTER margin than the
overshoot value, so the design is now honestly "derived, not tuned"
everywhere.** Suite 7/7 green.

### I1 — `Qf(6,6)` was a units error (reviewer correct; branch (a) ships)

`Qf(6,6) := qVelZ` mixed a RUNNING cost into a TERMINAL cost. `Q` has units
of cost per state^2 per unit *time*; `Qf` is cost per state^2. Equating
them is dimensionally wrong, and §3c's stated intent ("make `P(tf)`
consistent with the phase-B steady state") was right while the implemented
quantity was not: the steady state of the Riccati is `P_ss`, not `Q`.

The scalar ARE for this channel, `-p^2/(r m^2) + q = 0`, gives

```
P_ss = mB*sqrt(qVelZ*R(3,3))
```

which reproduces the design bandwidth at `tf`, since
`K/m = P_ss/(r m^2) = sqrt(qVelZ/r)/m = omegaVz`. Measured, confirming the
review exactly:

| | shipped (wrong) | ARE-consistent (now) |
|---|---|---|
| `Qf(6,6)` | 111.29 | **9.2906** |
| `K(3,6)/m` at `tf` | 158.9 1/s | **13.27 1/s** |
| vs design `omegaVz` | 13.2x overshoot | matches (12, mass-adjusted) |

So the shipped terminal loop really was a saturation-dominated overshoot
boundary layer over the last 0.14 s, not the first-order loop the header
described. **Branch (a): fixed and re-ran the full battery — 7/7 still
holds, and every single case improved.** No need for branch (b); the
comment now derives `P_ss` rather than excusing a spike.

### I2 — the monotone-descent assert was vacuous

It tested `diff(zSorted) > 0` on the *already sorted* array, so it could
only ever catch duplicate altitudes. Now `all(diff(sol.X(3,:)) < 0)` on the
original ordering — a real guard on the property altitude indexing actually
needs, which Task 8's Monte Carlo inherits.

### I3 — the step branch could fire on smooth interior arcs

`isTrans = ~(allLo || allHi)` admits any segment not pinned to one bound —
including a purely interior arc. The reviewer's demonstration (a
`0.60 -> 0.64*Tmax` arc, nowhere near a switch, whose `s` lands inside
`[0,h]`) would have been rewritten as a full `Tmin->Tmax` bang. Latent on
this campaign's certified bang-bang solution, live the moment Phase-2 drag
or a pointing cone introduces a real interior arc.

Fixed by requiring the segment to also TOUCH a bound
(`isTrans = ~(allLo || allHi) && anyOn`), which keeps the two genuine
switch segments and rejects interior arcs. This also puts the previously
dead `onB` helper to work rather than deleting it. Two comments corrected
in the same pass: the superseded `s`-formula (it documented the original
endpoint-levels version, not the shipped `Tmin`/`Tmax`-levels form
`s = (Isimp - h*B)/(A - B)`), and the **false** Phase-2 NOTE that claimed
the `[0,h]` range check was what protected smooth arcs — it is not, and
saying so was exactly the kind of unverified reassurance that hides a bug.

### I4 — `annulus_switch` threshold was not de-rated

`mid = 0.5*(P.Tmin + P.Tmax)` used the ENGINE ceiling to bisect a
GUIDANCE annulus whose top is `etaT*Tmax`. Correct by luck at
`etaT = 0.87` (536.6 kN and 591.5 kN both fall between the 338 kN and
735.15 kN arcs) but it silently mis-detects for `etaT < ~0.70`. Now
`0.5*(P.Tmin + P.etaT*P.Tmax)`. Shipped `tSwitch = 6.155 s`.

### M2 — stale comments swept

* `tvlqr_design.m` ADAPTATION 4 narrated the round-4 `zTermBand` remedy as
  current; marked HISTORICAL and pointed at the altitude-indexed note.
* `booster_params.m` `P.zTermBand = 150` deleted (grep-confirmed no code
  reference), replaced by a short note recording that altitude indexing
  retired it. `P.vf`'s cross-reference to it also updated.
* `print_certify_report.m` G4 label `< 1e-4*Tmax/m0` -> `< 1e-4*TmaxG/m0`,
  matching what `certify_pdg` actually computes.

### Final battery (post-I1) — 7/7, improved across the board

| case | landed | miss (m) | vtd (m/s) | vtd before I1 |
|---|---|---|---|---|
| nominal | yes | 0.007 | **0.978** | 1.453 |
| `dr0=[50;-30;0]` | yes | 1.582 | **0.968** | 1.451 |
| `dr0=[0;0;-50]` | yes | 0.395 | **0.977** | 1.452 |
| `dr0=[0;0;+50]` | yes | 2.053 | **1.004** | 1.466 |
| `thrust_scale=0.95` | yes | 3.148 | **1.282** | 1.647 |
| `thrust_scale=1.05` | yes | 3.175 | **1.232** | 1.639 |
| combined | yes | 4.108 | **1.060** | 1.454 |

**Worst case: miss 4.108 m of 15 m (73% margin), vtd 1.282 m/s of
2.0 m/s (36% margin — up from 18%).** Removing the terminal overshoot did
not cost robustness; it bought it. Every case is a genuine `z=0` touchdown.

Guidance-side gates are untouched by any of this (TVLQR weights do not
enter the guidance): G2 pos 0.00885 m, G3 `|dmf|` 0.42792 kg, `|dtf|`
0.00341 s, G5 bound fraction 0.9917, `all_pass` PASS at both grids. Pins
unchanged.

### Suite

```
test_params PASS (min-throttle dry T/W = 1.346)
test_dynamics_jac PASS (vacuum + drag)
test_colloc_smoke PASS  tf=16.60 s  mf=26464.3 kg  fuel=3535.7 kg
test_convex_lossless PASS  gap=4.81e-05  mf=26435.3 kg
test_certify_nominal (coarse) PASS          [ALL GATES PASS]
test_certify_nominal (nominal, all_pass) PASS  [ALL GATES PASS]
test_tvlqr_riccati PASS
  battery dr0=[0;0;-50]        land=1 miss=  0.395 m vtd= 0.977 m/s  PASS
  battery dr0=[0;0;+50]        land=1 miss=  2.053 m vtd= 1.004 m/s  PASS
  battery thrust_scale=0.95    land=1 miss=  3.148 m vtd= 1.282 m/s  PASS
  battery thrust_scale=1.05    land=1 miss=  3.175 m vtd= 1.232 m/s  PASS
  battery combined             land=1 miss=  4.108 m vtd= 1.060 m/s  PASS
test_closed_loop_nominal PASS  (nom miss 0.007 m vtd 0.978 m/s, disp miss 1.58 m vtd 0.968 m/s)
```
**7 of 7 green.**

### Files changed (polish)

- `lib/tvlqr_design.m` — I1 `Qf(6,6) = mB*sqrt(qVelZ*R(3,3))`; I4 de-rated
  switch threshold; M2 historical marker.
- `lib/hs_quad_ctrl.m` — I3 `anyOn` condition, corrected `s`-formula
  comment, corrected false Phase-2 NOTE.
- `lib/sim_closed_loop.m` — I2 real monotone guard.
- `lib/booster_params.m` — M2 `P.zTermBand` retired, `P.vf` xref updated.
- `certify/print_certify_report.m` — M2 G4 label.

### Self-review (polish)

- I1 was a genuine error of mine and the review's diagnosis was exactly
  right, including the numbers. I took branch (a) rather than the
  self-serving branch (b): the correct value was tried FIRST and the
  battery re-run before deciding, and it turned out strictly better.
- I3's false NOTE is the more instructive miss: I had written a
  reassurance about a failure mode I never tested. It is now replaced by a
  note that says plainly what the old comment got wrong and which term
  actually does the protecting.
- `P.zTermBand`'s removal was grep-verified against code before deleting,
  not assumed from memory.
