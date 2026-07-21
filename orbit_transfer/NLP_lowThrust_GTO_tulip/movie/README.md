# movie/ — minimum-fuel solution animation

Animation of the verified min-fuel arrival-leg solution (direct-NLP), in the
rotating Earth-Moon CR3BP frame.

- `minfuel_solution.mp4` — the movie (10.8 s, 24 fps, 1072x766). The min-time
  GTO->tulip spiral is faint grey (context); the min-fuel leg is bold, red on
  the burn arc and blue on the coast, with primer-direction thrust arrows during
  the burn. Bottom panel: the bang-bang throttle with a time cursor. Text
  readout: leg time, current mass, propellant burned, BURN/COAST state.
- `minfuel_solution.gif` — the same animation as a looping GIF (640 px wide,
  130 frames / ~12 fps, ~3 MB) for slides / e-mail / web.
- `minfuel_solution_still.png` — a single frame at the burn->coast switch.
- `minenergy_solution.mp4` / `.gif` — the direct min-ENERGY full-spiral
  solution (blue), moving spacecraft marker over a synced throttle subplot
  showing the smooth saturated ramp (never bang-bang). Full 40-rev GTO->tulip.
- `manysw_minfuel_solution.mp4` / `.gif` — the MANY-SWITCH min-fuel control
  from the tf-continuation (fine2 pass, tf=1.08x min, 53 throttle switches,
  76% burn / 24% coast). Trajectory red=burn / blue=coast (richly interwoven),
  throttle subplot shows the dense switch structure. HONEST NOTE: this solution
  is coherent at defect 1.7e-3 (hit the iteration cap), so the throttle is
  oscillatory (0.4-1.0), NOT yet sharp bang-bang -- a machine-tight re-solve
  would snap it to crisp 0/1 squares. Source:
  `manysw_minfuel_control_solution.mat`; animator `animate_minfuel_control.m`
  (generalized -- takes any saved solution file).
- `coarse_minfuel_solution.mp4` / `.gif` — the full-spiral MIN-FUEL control
  found by the tf-continuation (coarse pass, tf=1.10x min, 6 switches).
  Trajectory colored red=burn / blue=coast, thrust arrows on the burn, over
  the 6-switch throttle subplot. Mostly-burn with brief coasts at the early
  perigees (99.1% on) -- the many-switch regime is what the fine pass targets.
  Source: `coarse_minfuel_control_solution.mat` (from `tf_continuation_minfuel`);
  animator `animate_coarse_minfuel.m`.
- `three_way_comparison.mp4` / `.gif` — min-time (orange, always burn),
  min-energy (blue, ramp), and min-fuel (green, bang-bang ARRIVAL LEG ONLY)
  overlaid in the rotating frame, animated by normalized progress. The
  throttle subplot is the story: flat-on vs ramp vs the bang-bang step.
  NOTE: min-fuel is only the arrival leg (the full min-fuel spiral defeats
  both methods); min-time and min-energy are full transfers.

## Comparison movies — regenerate

```matlab
cd optimal_control/orbit_transfer/NLP_lowThrust_GTO_tulip/movie
gen_compare_data            % assembles compare_data.mat from the min-fuel
                            %  bundle (min-time + tulip) + energy_pipeline.mat
animate_energy_solo('movie')   % -> minenergy_solution.mp4 + .gif
animate_three_way('movie')     % -> three_way_comparison.mp4 + .gif
```

`animate_*('preview')` writes stills instead. The min-energy trajectory
comes from the N=4000 direct solve (`energy_pipeline.mat`, defect 3e-4);
see `../MIN_ENERGY_NOTES.md`.

## Regenerate

```matlab
cd optimal_control/orbit_transfer/NLP_lowThrust_GTO_tulip/movie
gen_movie_data      % re-solves: min-time spiral + min-fuel NLP leg (N=3000)
                    %  -> minfuel_movie_data.mat   (a few minutes)
animate_minfuel('preview')   % three stills (early / switch / late)
animate_minfuel('movie')     % renders minfuel_solution.mp4 AND .gif
```

GIF knobs (top of the `movie` block in `animate_minfuel.m`): `gifStride`
(frame decimation), `gifW` (width in px), `gifDelay` (s/frame). Smaller
`gifW` or larger `gifStride` shrinks the file.

`minfuel_movie_data.mat` (the solved trajectory) is included, so `animate_minfuel`
runs without re-solving. Frames are rasterized with `exportgraphics` and
assembled by `VideoWriter` (headless-safe; `getframe` is unreliable under
`matlab -batch`).

Solution shown: leg from tau = 4.0 on the min-time arc, fixed leg time
1.3x the leg minimum; propellant 1.0627 kg, 76.9% burn fraction, single
burn->coast switch, max defect 2e-16. See ../README.md and
../../lowThrust_GTO_tulip/gto_tulip_mintime_theory.pdf (S6).
