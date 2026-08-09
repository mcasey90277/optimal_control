# CLAUDE.md — missiles/

3-DOF hypersonic-glide-vehicle and ballistic-missile trajectory library. Built
2026-08-07 to 2026-08-09. Read this before touching anything here; several of
the conventions below are load-bearing and at least three of them have already
cost a debugging cycle each.

## What this is

A pumpkyn-style `+coorbital` package plus five entry scripts. The library
propagates prescribed trajectories through boost, glide, coast and descent, and
solves point-to-point targeting. It is **not** a trajectory optimizer — the
optimal-control work lives in the sibling `orbit_transfer/` and
`booster_landing/` campaigns.

```
+coorbital/
  +eom/     glide3DOF (6-state), boost3DOF (7-state), massConstant (lifts 6->7)
  +prop/    phaseRun, constThrust, eventAltitude, eventBurnout, eventApogee
  +atmos/   expAtmos          +grav/  sphereGrav        +aero/  constLD
  +guide/   prescribed, pitchProgram, terminalConstraint
  +util/    missileConst, vehicleDefaults, boosterDefaults, greatCircle,
            greatCircleBearing, rangeSolve, aimSolve
  +viz/     groundTrack, profilePlot, globe3D, globeMovie, saveFigure
            private/ vizParent, earthSurface, vizOption
HGV/  run_glide, run_boost_glide, run_target, vehicle_hgv
BM/   run_ballistic, run_ballistic_target, vehicle_bm
tests/ run_tests + 24 test_*.m        docs/  see below
```

## Running things

```bash
/Applications/MATLAB_R2025b.app/bin/matlab -batch \
  "cd('/Users/msc/Desktop/optimal_control/missiles'); run('tests/run_tests')" 2>&1 \
  | grep -v "Home License" | grep -v "personal use"
```

- **A bare `tests/run_tests` inside `-batch` parses as a division** and fails.
  Always `run('...')`.
- The suite takes about five minutes. Allow 900 s.
- Never put embedded newlines in a `-batch` string; write a driver `.m` instead.
- Baseline is **24 passed, 0 failed, zero warnings**. Zero warnings is part of
  the contract, not a nicety.

## Conventions that will bite you

**State.** `x = [r, lon, lat, V, gamma, psi]`, plus mass as a 7th component on
powered chains. `gamma` positive up, `psi` clockwise from north, `V` is
**planet-relative**. Everything SI and radians inside the library; human units
appear only in an entry script's fenced user block and are converted once, in a
marked place.

**Gravity returns `[gr, gLat]` where `gr` is a POSITIVE MAGNITUDE** of inward
attraction, not a signed component. Getting this backwards produces a
plausible-looking trajectory.

**Models are injected as handles** on `env`: `atmos`, `grav`, `aero`, `prop`.
The EOMs name no model. This is the whole architecture — do not add a switch, a
class, or a subclass.

**Mass contract.** `glide3DOF` divides by `veh.mass`, not `x(7)`, so a carried
mass is inert on unpowered phases. `boost3DOF` uses `x(7)` in every denominator.
`massConstant` guards the mismatch with `coorbital:massConstant:massMismatch`.

**No `%#ok` pragmas of any kind.** Never `i` or `j` as loop or index variables
(imaginary unit). Never `norm` — use `sqrt(sum(...))`. No hard-coded physical
constants outside `coorbital.util.missileConst()`.

**`missiles/results/` is gitignored** (`.gitignore:48`). Movies and saved
figures are local artefacts, never tracked.

## Testing discipline

A test is not trusted until it has failed. Every test-bearing change ends with
**mutation testing**: break the thing deliberately, confirm the test FAILS,
restore byte-identically and verify with `md5`. This has repeatedly caught
tests that could not fail — including a whole refusal path where deleting the
refusal left the suite green.

Independent references must use a **different solver** than the driver
(`ode89` at 1e-12 against `ode45` at 1e-10), or they are not independent.

`docs/LESSONS_LEARNED.md` and the SDD's "structural blindnesses" section carry
the four failure modes found the hard way. The newest: **an unreachability
argument is only as good as the bound it rests on, and a relative bound is not
a magnitude bound** — `cond(J) = s1/s2` constrains a ratio, so both singular
values can underflow together at condition number exactly 1.

Two operational traps, both of which have already cost a cycle:

- **A copied tree on the MATLAB path shadows the real files.** If you copy the
  library to mutate it, put the copy where `addpath(genpath(...))` cannot reach
  it, and `restoredefaultpath` before any measurement you intend to report.
- **A number copied into prose is not covered by the test that produced it.**
  Three separate documents have been caught quoting figures that were correct
  when written and stale when read. Re-measure before quoting.

## State of the work

Two-axis targeting is done and validated. `coorbital.util.aimSolve` solves
launch azimuth beside the range control, which is what lets both entry scripts
target a **rotating Earth** and steer **cross-range** — they are the same gap,
because the flown ground track is not the great circle the old closed-form
bearing assumed.

| case | seed miss | achieved | residual evaluations |
|---|---|---|---|
| `run_target`, rotating | 231 551.628 m | 4.361 m | 7 |
| `run_target`, 75-degree bank | 21 524.695 m | 54.795 m | 4 |
| `run_ballistic_target`, rotating | 463 211.19 m | 52.46 m | 7 |

**Rotation deflects the vehicle**; the target is a fixed ground point and does
not move in the planet-relative frame. The old "lead the target" explanation is
retracted throughout — if you find it anywhere, it is a bug in the prose.

Protected numbers, which must not move: `run_glide` 6986.82/2073.77/1.11,
`run_boost_glide` 7663.05/2194.77/8.33, `run_target` shipped 511.2434604708133 m,
`run_ballistic_target` shipped minimum-energy 39.0092687975735 m.

## Next milestone, and what blocks it

PEG on the ballistic boost phase. `docs/closed_loop_guidance.md` is the brief —
it is concrete enough to implement from. **Nothing of PEG or VOA is coded.**

Two decisions were settled by measurement and should not be relitigated: the
guidance driver sits **above** `phaseRun`, one `phaseRun` call per guidance
cycle (20.7% of integrator guidance calls happen at times already passed, and
148 repeat a time with a different state, so no guidance update can fire inside
a `guide(t,x)` call); and the thrust-direction inverse is the `atan2` form, for
**scale robustness** rather than for near-zero-incidence cancellation.

**A vehicle-file decision blocks the dispersion campaign.** The shipped nominal
boost is pinned at the −6 degree angle-of-attack clamp for 74.1% of the burn,
and that clamp does not belong to the vehicle being steered:
`coorbital.util.boosterDefaults()` has no `alphaMaxDeg` field at all, so the
script applies the separated re-entry body's limit (`BM/vehicle_bm.m`) to the
boost stack. Run a Monte Carlo against that and it measures the clamp, not the
guidance law.

## Docs

- `docs/hgv_dynamics_note.tex` — 43 pp, the mathematics. The authority; other
  documents defer to it.
- `docs/software_design.tex` — 48 pp, architecture and data flow.
- `docs/README.md` — the authority on delivered behaviour.
- `docs/closed_loop_guidance.md` — the PEG/VOA brief.
- `docs/DESIGN.md` — chronological as-built record. Append dated sections; do
  not edit the existing ones.
- `docs/LESSONS_LEARNED.md`, `docs/reviews/` — review archive.

Compile LaTeX with `/Library/TeX/texbin/pdflatex`, twice, then clean the aux
files (`rm -f *.aux *.log *.out *.toc *.lof *.lot *.fls *.fdb_latexmk
*.synctex.gz`). That cleanup is automatic, not on request.
