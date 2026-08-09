# Task 9 Report: Visualization (`plot_pdg_solution`, `plot_footprint`, `movie_landing`)

## Implementation

Three files created under `viz/`, plus `tests/test_viz_smoke.m`:

- `viz/plot_pdg_solution.m` — `fig = plot_pdg_solution(solC, solV, outfile)`. 2x2
  tiled layout: (a) 3D trajectory with a translucent glideslope cone
  (`r = cotg*z`) and thrust-vector `quiver3` arrows every 5th node, both
  solvers overlaid; (b) throttle `||T||/Tmax` "money plot" for both solvers
  with **both** bound lines labeled (`Tmin/Tmax=0.40` dotted, `etaT=0.87`
  dashed — the guidance ceiling, not 1.0); (c) mass vs time with an `mdry`
  reference line; (d) speed vs time with the terminal target `|v(tf)|=1.5
  m/s` (not zero) marked. `exportgraphics` at 200 dpi if `outfile` given.
- `viz/plot_footprint.m` — `fig = plot_footprint(mc, P, outfile)`. Landing
  scatter colored by `mc.ok` (blue=success/orange=failure, colorblind-safe),
  pad-radius circle, 3-sigma dispersion ellipse from `cov(mc.land)` (via
  eigendecomposition, centered at the sample mean), and an annotation box
  with success rate, the `n_landed/n_arrest/n_horizon` failure-mode
  breakdown, and vtd mean/max vs the `P.vtd_max` gate.
- `viz/movie_landing.m` — `movie_landing(out, sol, P, outfile, opts)`.
  Implemented per the brief's skeleton essentially verbatim (left: altitude
  view + thrust arrow + pad/ground; right: throttle trace + moving cursor),
  with two additions: (1) a `drawnow` before every `getframe` (skill
  requirement, not in the skeleton) and (2) all three bound lines on the
  throttle panel — `Tmin/Tmax`, `etaT` (guidance ceiling), and `1.0` (engine
  max / tracker authority) — each labeled, since the closed-loop trace
  plotted here is `out.Tcmd`, which can legitimately ride above `etaT` up to
  the full engine ceiling. Frame size locked to 1280x720 via `imresize` on
  every frame (divide-by-16 H.264 shear guard). `VideoWriter` MPEG-4, 30 fps
  default.

### Key interface facts verified before coding

- `hs_quad_ctrl`'s actual signature is
  `hs_quad_ctrl(tt, U, Um, h, N, Tmin, Tmax)` — **not** `(t, sol.U, sol.Um,
  h, N, P)** as a naive reading of the task instructions might suggest. Used
  it as `hs_quad_ctrl(tq, solC.U, solC.Um, h, Nseg, P.Tmin, P.etaT*P.Tmax)`
  so the guidance throttle panel is bounded by the actual annulus the
  collocation NLP solved against (engine floor, de-rated ceiling), matching
  `solve_pdg_colloc.m`'s `Tmax_h = P.etaT * P.Tmax / Fc`.
- `solV` (convex) has no `.Um` field — its throttle curve is sampled
  directly off its own dense `Nconv`-node grid (`sqrt(sum(solV.U.^2,1))`),
  no reconstruction needed since it isn't Hermite-Simpson.
- `mc.n_landed` counts `stop=='touchdown'` regardless of miss/vtd/mass gates
  (i.e., can exceed `nnz(mc.ok)`); `mc.ok` is the full success gate. Both are
  reported in the footprint annotation as specified.

## TDD Evidence

1. Wrote `tests/test_viz_smoke.m` verbatim per the brief.
2. Ran it before implementation:
   ```
   Unrecognized function or variable 'plot_pdg_solution'.
   Error in test_viz_smoke (line 17)
   ```
   Confirmed FAIL as expected (functions undefined).
3. Implemented the three viz functions.
4. Re-ran:
   ```
   test_viz_smoke PASS
   ```

## Self-Check (production-grid render, scratchpad only)

Rendered real nominal products at N=60 (collocation), Nconv=120 (convex),
20-run Monte Carlo, and a 4 s / 30 fps movie into the scratchpad (not
`results/` — that's Task 10's job). Self-check via `imfinfo`/`imread` and
`VideoReader`, run headless in MATLAB R2025b:

```
sol.png:  2688x2250  286910 bytes
fp.png:   1375x1518  107978 bytes
movie:    1280x720  fps=30.00  duration=4.000 s  NumFrames(expect~120)
movie:    actual decoded frame count = 120 (expect 120)
```

All three files exist, reload to their reported dimensions, and the movie's
decoded frame count matches `duration*fps` exactly (120). Frame size is
confirmed locked to 1280x720 (the divide-by-16 H.264 shear guard).

File sizes:
- `plot_pdg_solution.png`: 286,910 bytes
- `plot_footprint.png`: 107,978 bytes
- `movie_landing.mp4`: 298,734 bytes

(All three live only under the scratchpad path used for this task, not
committed and not written into the repo's `results/`.)

## Self-Review: does the throttle plot show the switch as a STEP?

Targeted probe: sampled `hs_quad_ctrl` at the same 400-point resolution
`plot_pdg_solution.m` uses across `[0, solC.tf]` (N=60 production grid,
tf~16.6 s):

```
N=400 grid: max |d(T/Tmax)| per sample = 0.4700 at t=6.031->6.072 (r: 0.400->0.870)
```

The entire rise from the Tmin floor (0.40) to the etaT ceiling (0.87) —
the full campaign-defining bang-bang switch — occurs between two adjacent
plotted samples (~0.04 s apart). This confirms the switch renders as a
near-vertical STEP in the money plot, not a smeared ramp; also verified at
a finer 2000-point grid (step lands within one 0.0083 s sample, same
0.40->0.87 jump), ruling out an artifact of the coarser 400-point plotting
grid. This is the direct payoff of using `hs_quad_ctrl`'s per-segment
step reconstruction on the transition segment instead of a global pchip
spline (see that function's own extensive header notes on why a spline
smears this).

Other review notes:
- Speed plot (d) correctly targets `|v(tf)|=1.5 m/s`, not 0 — read from
  `P.vf`, not hardcoded.
- Both throttle bound lines (`Tmin/Tmax` and `etaT`) are drawn and labeled
  in `plot_pdg_solution`'s money plot; the movie additionally shows the
  engine's true `1.0` ceiling since `out.Tcmd` (closed-loop tracker) can
  legitimately exceed `etaT`.
- Footprint success/failure coloring reuses the same colorblind-safe
  blue/orange pair as the solver comparison plots — a deliberate, consistent
  choice, not a requirement bleed (the brief's blue/orange requirement was
  scoped to colloc/convex, but reusing the same accessible pair for a
  different success/failure binary is good practice, not a conflict).
- No emoji anywhere in the three files.
- All three functions carry the same header-block convention (PURPOSE /
  INPUTS / OUTPUTS / REFERENCES) used throughout `lib/` in this campaign
  (`booster_params.m`, `hs_quad_ctrl.m`, etc.) — this codebase's existing
  functions do not use the full pumpkyn self-demo (`if nargin==0`) style, so
  none was added here either, for consistency with siblings.
- No `i`/`j` used as loop variables (`k` throughout).

## Concerns

- None blocking. Minor: the 3D trajectory panel's `quiver3` thrust-arrow
  scale (`0.25*zmax/P.Tmax`) is a fixed visual heuristic, not derived from a
  requirement — reasonable at this campaign's geometry but would need
  re-tuning if `P.r0`/`P.Tmax` change by an order of magnitude.
- `plot_footprint`'s 3-sigma ellipse is centered at the *sample* mean
  `mean(mc.land)`, not at the pad origin — this is the standard dispersion-
  ellipse convention (shows the actual spread of the draws) and the pad
  center is separately marked with a yellow pentagram, but it's worth
  flagging in case a reviewer expected the ellipse centered on the pad.

## Files Changed

- `/Users/msc/Desktop/optimal_control/booster_landing/viz/plot_pdg_solution.m` (new)
- `/Users/msc/Desktop/optimal_control/booster_landing/viz/plot_footprint.m` (new)
- `/Users/msc/Desktop/optimal_control/booster_landing/viz/movie_landing.m` (new)
- `/Users/msc/Desktop/optimal_control/booster_landing/tests/test_viz_smoke.m` (new)

---

## Fix Report (2026-08-09, review round 1)

Review verdict: Needs fixes — two Important items (the `hs_quad_ctrl`
wiring and money-plot content were independently verified correct and
needed no change).

### Important 1: `movie_landing.m` ignored the binding colorblind-safe scheme

The movie inherited the brief skeleton's plain `'r'`/`'b-'`/`'g-'`/`'k'`
literals verbatim; the blue `[0 0.447 0.741]` / orange `[0.85 0.325
0.098]` pair binds to all three deliverables, not only the two static
plots. Fixed by defining `colBooster`/`colAccent` locals and remapping:

- Booster marker (left panel) and the throttle trace (right panel): now
  `colBooster` (blue) — same physical entity (the vehicle / its commanded
  throttle), one color.
- Thrust-plume arrow and the moving time cursor: now `colAccent` (orange)
  — the two "attention" elements that move each frame.
- Bound lines (`Tmin/Tmax`, `etaT`, `Tmax`) and the static ground-track
  path sketch: left neutral black/gray, as instructed — they're reference
  geometry, not data series.
- Pad marker: kept a distinct green, **deliberately**, with a header note
  stating why — a landing pad/helipad has its own real-world conventional
  color (aviation safety green) that a viewer already reads as "target,"
  independent of this campaign's solver blue/orange scheme. Folding it
  into blue/orange would have made the pad harder to spot, not easier,
  which is the opposite of the point of a color scheme.

### Important 2: `plot_pdg_solution`/`plot_footprint` leaked invisible figure handles

Both returned `fig` (`Visible` off) with nothing ever closing it. Across a
long `-batch` process running Task 10's flagship reproduction plus Task
11's sweeps, repeated calls would accumulate invisible figures unbounded.
Fixed with an explicit ownership contract, documented in each function's
OUTPUTS block:

- Output requested (`fig = plot_...(...)`)  → caller owns the handle, must
  close it.
- Bare-statement call (`plot_...(...);`, `nargout==0`) → the function
  closes its own figure right after any `exportgraphics` call, so a
  fire-and-forget call (the common case for "just write me a PNG") never
  leaks.

`movie_landing.m` was already correct here (`close(fig)` unconditionally
at the end, no output returned) and needed no change for this item.

`tests/test_viz_smoke.m` updated to exercise the caller-owns half of the
contract explicitly: `fig1 = plot_pdg_solution(...); close(fig1);` and
`fig2 = plot_footprint(...); close(fig2);` (the bare-statement half is
already exercised implicitly by every other caller in the repo that
doesn't capture an output). Also added Minor 1: an `onCleanup` that
`rmdir(od, 's')`s the `bl_viz_smoke` tempdir at test end (pass or fail),
replacing the prior "leave it in `tempdir`" behavior.

### Covering verification

Re-ran `test_viz_smoke` end to end:
```
test_viz_smoke PASS
```
Confirmed the tempdir cleanup fired: `exist(fullfile(tempdir,
'bl_viz_smoke'), 'dir')` returns `0` immediately after the test completes.

Re-rendered a standalone 2 s/30 fps movie (post color-remap) into the
scratchpad and re-verified with `VideoReader` that the frame-size lock
survived the change:
```
movie: 1280x720  fps=30.00  duration=2.000 s
decoded frame count = 60 (expect 60)
file: .../mov_recheck_2s.mp4  168374 bytes
recheck_movie PASS
```
1280x720 and exact frame count (60 = 2*30) confirm the divide-by-16 H.264
guard is independent of the color changes, as expected (`imresize` runs
on `frame.cdata` after `getframe`, downstream of any color choices made
during rendering). The recheck script and its output MP4 were deleted
from the scratchpad after verification (not a committed artifact).

### Files touched in this fix round

- `viz/plot_pdg_solution.m` — ownership-contract doc + `nargout==0` close
- `viz/plot_footprint.m` — ownership-contract doc + `nargout==0` close
- `viz/movie_landing.m` — colorblind-safe remap (booster/throttle blue,
  plume/cursor orange, bounds/path neutral, pad deliberately green) + doc
- `tests/test_viz_smoke.m` — capture-and-close both static-plot calls,
  `onCleanup` tempdir removal
