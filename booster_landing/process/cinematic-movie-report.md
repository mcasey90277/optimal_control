# Cinematic landing movie — build report (2026-08-09)

Deliverable: `viz/movie_landing_cinematic.m` + `results/landing_cinematic.mp4`.
`viz/movie_landing.m` (the diagnostic movie the front door renders) was not
touched. The cinematic movie is **not** wired into `run_booster_landing`.

## Verified output

| Property | Value |
|---|---|
| Frame size | 1920 x 1080 (both divisible by 16, locked via `imresize`) |
| Frame rate | 60.00 fps |
| Duration | 30.73 s |
| Frames | 1844 (= 120 card + 1592 flight + 132 hold) |
| File size | 60.3 MB |
| Render time | 93.0 s (MATLAB R2025b, `-batch`, synchronous) |

Verification was `VideoReader` on the written file plus three spot frames
(130 / 850 / 1780) read back as images: frame 130 lands inside the
card-to-scene fade (fade working), 850 is the mid-slam descent, 1780 is
inside the touchdown hold with the plume already faded out (engine-cutoff
ramp working).

`tests/test_viz_smoke.m` gained a cinematic block (trimmed card/hold/
slow-lead, 20 fps) asserting existence, non-zero size, `1920x1080`, and a
plausible frame count. `test_viz_smoke`, `test_run_front_door` and
`test_phase2_drag` (the three tests that touch viz products) all PASS.

## Design decisions

**Data source.** Everything animated comes from the closed-loop trace
`out` (`out.t`, `out.X`, `out.Tcmd`) — what actually flew under the TVLQR
tracker. `sol` is used only for the title card's headline numbers
(`sol.tf`, `P.m0 - sol.mf`), which are guidance quantities and are
labelled as such. All per-frame quantities, axis limits and the entire
camera path are precomputed once (house rule: no autoscale jitter).

**Body-axis convention.** The model is 3-DOF: a point mass with a thrust
vector and no attitude state. The booster is drawn aligned with the
*current thrust direction* — the standard convention for 3-DOF PDG
animations, stated in the function header and in small print on screen.
The consequence is visible: the vehicle leans, up to **21.3 deg** off
vertical at burn start and **19.8 deg** at touchdown (measured from
`out.Tcmd`). That is the real primer direction, not a rendering error, so
the instrument strip carries a live `THRUST TILT` readout. Adding that
readout was a direct result of looking at a key frame and mistrusting the
lean.

**Scale — true, with an annotation instead of an inflation.** The booster
is 3.7 m x 45 m at true scale (`opts.boosterScale = 1`). In the
establishing shot that is ~20 px tall, and the first instinct was to
scale the vehicle up. Measured against the framing, a *constant*
exaggeration cannot work: whatever factor makes the vehicle readable at 2
km makes it 4-6x the pad diameter at touchdown, where the 15 m pad is the
size reference and the lie becomes obvious. A time-varying scale would be
worse. The chosen answer is a **camera-facing tracking reticle** whose
radius is a fixed fraction of the view height (so it holds constant
screen size) and which fades out as the real vehicle grows — annotate the
vehicle, don't inflate it. `opts.boosterScale > 1` is still available and
prints the factor on screen for the whole movie if anyone uses it.

**Camera.** Framing is specified as a desired *view height in metres* at
the target plane; camera distance follows from it and the 34 deg view
angle. A wide establishing law (`1.15 z + 300`, whole descent plus pad in
frame) pushes in over the first ~4 s to a working law (`0.30 z + 165`,
vehicle at 50-250 px). Azimuth swings -52 deg -> -108 deg and elevation
21 deg -> 6 deg on a cosine ease over **flight-time** progress, not frame
index — the slow-motion finale owns ~40% of the frames, so easing per
frame would have spent 40% of the camera move on the last 3 seconds. The
target sits below the vehicle by a fraction of the view height (0.32,
easing to -0.10), which composes the vehicle high in frame early and
centred over the pad at the end.

**Plume.** Two nested translucent cones (outer orange, inner white-hot),
alpha decaying along the length, colour ramped by throttle, length
strictly proportional to `|T|/Tmax` (2.0 booster lengths at full
throttle). The jet is clipped where it meets the ground and a soft-edged
radial "ground wash" glow takes over, so the exhaust turns at the surface
instead of passing through it.

**Timeline.** 2.0 s title card -> 0.6 s fade -> real-time flight until 3.9
s before touchdown -> cosine ramp to 0.25x for the finale (with a `0.25x`
tag) -> 2.2 s hold. Through the hold the plume is ramped to zero: the sim
ends *at* touchdown, so replaying its last frame verbatim would leave a
landed booster firing its engine for two seconds.

**Instrument strip** (bottom 15%): throttle gauge with marks at
Tmin/Tmax = 0.40, etaT = 0.87 and 1.00; altitude, speed and thrust-tilt
digital readouts (fixed-width); a draining propellant bar; `SLAM` flash at
the bang-bang switch and `TOUCHDOWN <vtd> m/s` at the end.

## Iteration log (seven look-loop rounds, PNG key frames read as images)

1. **Round 1** — first render. Title card fine. Every flight frame was
   effectively empty: the camera target law (`0.50 z + 42`) looked far
   *below* the vehicle, which sat off the top of frame. Camera law
   rewritten around view height + a target offset tied to it.
2. **Round 2** — vehicle now framed, but a true-scale booster at 1.4 km
   was a 2-pixel thread and the ground plate's far edge cut a hard line
   across the establishing shot. Added the tracking reticle; replaced the
   ground with an 11 km plate shaded to fade radially into the background
   (stand-in for atmospheric haze); fixed two colliding gauge labels
   (`etaT (guidance)` overlapped `Tmax`, and `|T| / Tmax` rendered as
   `ITI / Tmax`).
3. **Round 3** — the finale looked right; the establishing half did not,
   and the plume was a 2-pixel needle narrower than the vehicle. Added
   the wide->tight camera push-in, flared the plume (tip radius ~2.5x the
   body radius), moved `SLAM` off centre and recoloured it orange.
4. **Round 4** — the finale was framed too high over a huge empty apron;
   the plume ran through the ground; the drawn lean was unexplained.
   Target offset taken negative at touchdown, plume clipped at the
   surface with a ground wash added, `THRUST TILT` readout added and the
   strip re-spaced to fit it.
5. **Round 5** — the first ground wash was a hard-edged disc that read as
   a *shadow*. Replaced with a polar mesh whose opacity falls off
   radially.
6. **Round 6** — touchdown still slightly small in frame, and the plume
   at 180 m altitude ran off the bottom of the frame. Tightened the
   terminal view-height floor (230 -> 165 m) and put the engine-cutoff
   ramp on the hold frames.
7. **Round 7** — plume length made strictly proportional to throttle at a
   shorter constant; reticle fade-out re-tuned so it stays legible
   through the coast. All five key frames judged good; full render.

## Things found by looking that are worth recording

- **The commanded throttle really does chatter at the switch.** The trace
  steps 0.437 -> 0.907 at t = 6.0699 s, drops back to exactly Tmin for
  ~31 ms (6.085 -> 6.116 s), then steps up again and settles at 0.87. An
  unlucky key-frame choice landed inside that dip and showed `SLAM` next
  to a 0.40 gauge, which is what surfaced it. This is real closed-loop
  data (feedback crossing the annulus around the bang-bang switch), not a
  rendering artefact, and it is left in the movie as a 1-2 frame flicker.
  Worth a look if anyone ever revisits the tracker near the switch.
- The thrust tilt numbers above (21.3 deg at t = 0, 19.6 deg at
  touchdown) are also a genuine property of this solution that neither of
  the existing static plots shows.

## Legged-landing cosmetic ease (2026-08-09)

The drawn booster used to stay aligned with the true thrust vector all the
way to touchdown, so it landed leaning at its true ~20 deg primer tilt --
correct physics, but it reads wrong for a vehicle with landing legs.
Added a second, purely cosmetic layer: below 60 m altitude the *drawn*
body axis (cylinder + cap + engine-bay band + the reticle that tracks
them) cosine-eases away from the true thrust direction toward
world-vertical, landing exactly vertical at z = 0 and staying vertical
through the touchdown hold and engine-cutoff ramp. The ease is keyed on
`Rq(:,3)` (interpolated altitude), computed once as a vectorized `Ub`
array alongside the existing `Uq` (true thrust direction) array, so
`update_scene` now carries two unit vectors per frame: `uk` (true) and
`ub` (drawn). Everything that is a physics readout keeps using `uk`,
unmodified:

- the plume (both cones) and the ground-wash glow -- true exhaust
  direction, including the ground-clip length calculation (`dGnd`)
- the `THRUST TILT` instrument -- `acosd(uk(3))`, still the true
  commanded tilt (confirmed 19.8 deg at touchdown in the rendered
  key frame, unchanged from the 19.6-19.8 deg range recorded above)

Only the drawn cylinder/cap/band (`Rot = rot_to(ub)`) and the tracking
reticle's center (now `rk + 0.5*bLen*ub`, so the ring stays visually
centered on the vehicle it is annotating) use the eased vector. The
on-screen footnote was updated to say so: "body axis drawn along the
thrust vector; eased to vertical in the final meters for legged-landing
realism (3-DOF model has no attitude state)".

**Verification (PNG key frames read as images):**
- **901 m / T+6.24 s** (near the SLAM switch, well above the 60 m ease
  threshold): body drawn thrust-aligned as before, `THRUST TILT` reads
  2.5 deg -- visually unchanged from pre-ease behavior, confirming the
  ease does not activate at altitude.
- **Touchdown, T+16.63 s / 0 m**: booster drawn perfectly vertical,
  sitting on the pad reticle; `THRUST TILT` reads **19.8 deg**, the true
  commanded value, undisturbed by the cosmetic ease -- exactly the
  intended split between drawn attitude and instrumented truth.

**Re-verified after the change:** `results/landing_cinematic.mp4`
re-rendered via `README.md`'s documented invocation --
1920x1080, 60.00 fps, 1844 frames, 30.733 s (VideoReader), matching the
original render byte-for-byte in every measured property. `tests/
test_viz_smoke.m` (which includes the cinematic block) re-run and PASSES.
