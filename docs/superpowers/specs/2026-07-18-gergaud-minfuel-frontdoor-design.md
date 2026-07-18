# DESIGN — `run_gergaud`: a min-fuel elliptic→GEO front door

**Goal.** A single, clean, PARAMETERS-block entry script that reproduces one
row of Gergaud–Haberkorn–Martinon (JGCD 27(6), 2004) Table 3 — a low-thrust
**minimum-fuel** transfer from an elliptic, inclined start orbit to GEO — with
user-settable thrust level and user-settable initial *and* final orbits, then
emits the Table-3 row plus a trajectory plot and a movie (gif + mp4). Modeled
on `PSR/run_psr.m` and `elfo/elfo_run_one.m` (edit one block, run).

Date: 2026-07-18. Status: DESIGN (movie/adapter piece already built +
validated; see §5). Depends on the completed MEE thrust-ladder campaign
(`earth_elliptic_to_geo/DESIGN_thrust_ladder.md`, README "Campaign A"), whose
solver stack this script is a front door onto — it adds no new physics.

---

## 1. Context — what already exists

The MEE/L-domain campaign already solves this problem end to end; it just has
no single "set the thrust and go" front door. The reusable pieces:

| module | role |
|---|---|
| `mee_seed.m` | constant-throttle ode113-in-L warm-start seed |
| `casadi_lt_mee.m` | solver core (L-domain trapezoid, ΔL decision variable), modes `mintime`/`fixedtf` |
| `homotopy_mee.m` | guarded ε:1→0 energy→fuel sweep |
| `run_mintime_mee.m` | per-thrust free-longitude min-time anchor (`t_f,min`) |
| `run_transfer_mee.m` | one fixed-`t_f` fuel solve at `c_tf·t_f,min` → structure report |
| `run_ladder.m` | thrust-continuation orchestrator + Table-3-style summary |
| `psr_mee_refine.m` (+ helpers) | switch-aware mesh refinement (needed at 1 N / 0.5 N) |
| `transfer_movie.m` | inertial trajectory + control movie (Cartesian layout) |
| `elements_to_cart.m` | algebraic MEE→(r,v) |

Certified rungs already on disk (`results/`): `MEE_M2_10N` (1377.10 kg),
`MEE_M2_5N` (1364.54), `MEE_M2_2p5N` (1369.79), `MEE_M2_1N_PSR_psr_final`
(1371.44), `MEE_M2_0p5N` (1375.28, anchor-free + PSR, footnoted).

**Endpoints today are hardcoded:** the initial GTO (`P=11625 km, e=0.75,
i=7°`, apogee start at `L=π`) lives in `mee_seed.m`; the terminal GEO
(`P=1, ex=ey=hx=hy=0`) is baked into `casadi_lt_mee.m` (lines 152–156). The
length unit is `LU = 42165 km` (GEO radius), so GEO ⇒ `P=1`.

---

## 2. Endpoint parameterization (core-code change, default-preserving)

Make both endpoints real knobs, with defaults that reproduce today's behavior
byte-for-byte so no existing certified cache regresses.

- **`mee_seed.m`** — accept an initial MEE state (built from `P0_km, e0,
  i0_deg`, argument-of-perigee = 0 so `ex0=e0, ey0=0`; `hx0=tan(i0/2),
  hy0=0`) instead of the literal `X_init`, and take the seed's terminal
  `stopP` from the requested final `Pf`. Default = the paper's values.
- **`casadi_lt_mee.m`** — replace the five hardcoded terminal element
  equalities with `opts.xf` (5×1 `[P;ex;ey;hx;hy]`), **defaulting to
  `[1;0;0;0;0]`**. Single point of change: covers both `mintime` and
  `fixedtf`, and therefore both the anchor and the fuel solve.
- **Pass-through only** in `homotopy_mee.m`, `run_transfer_mee.m`,
  `run_mintime_mee.m`: add an optional `xf` (and initial-orbit) field, thread
  it down, default preserves current behavior. Add `xf`/initial-orbit fields
  to each driver's cache fingerprint — the existing `check_cache_fp*` guards
  already WARN-and-trust on schema-older caches, so the pre-existing certified
  `.mat` files remain valid for the default endpoints.

**Final-orbit units.** A custom final orbit is expressed in the *same*
`LU=42165 km` unit: `Pf = Pf_km/42165`, `ef`, `if_deg` → `xf`. The unit scale
is NOT rescaled to the new target, so the default GEO case is unchanged.

**Scope caveat (documented in the header, not silently assumed):** the
solver's prograde/equatorial-friendly setup and the seed's steering law were
validated for GEO-like targets. Non-circular / significantly inclined /
retrograde finals are research-probe territory — the script reports whether
the solve certified rather than presuming any target converges.

---

## 3. The front-door script `run_gergaud.m`

`PSR/run_psr.m` house style: a single **PARAMETERS** section the user edits,
then run. Sections: (0) paths, (1) PARAMETERS, (2) resolve endpoints → `x0`,
`xf`, (3) solve/reuse, (4) print Table-3 row, (5) plot, (6) movie.

**PARAMETERS block:**
- `thrustN` — one of `10 | 5 | 2.5 | 1 | 0.5 | 0.2 | 0.1` N
- initial orbit — `P0_km` (11625), `e0` (0.75), `i0_deg` (7)
- final orbit — `Pf_km` (42165), `ef` (0), `if_deg` (0)  ← **default = GEO,
  the standard first option**
- `c_tf` (1.5), `nodesPerRev` (25), `maxIter` (1500)
- `runMode` — `'auto'` | `'solve'` | `'probe'` (see §4)
- `makeMovie` (true), `makePlot` (true)

**Run modes:**
- `'auto'` (default) — if endpoints are the paper defaults AND a certified
  cache exists for this `thrustN`, **load and report it instantly** (+ plot +
  movie). Otherwise solve.
- `'solve'` — always run the live pipeline (anchor → fuel homotopy → PSR if the
  rung needs it), ignoring caches. Required for any custom endpoint (no cache
  applies).
- `'probe'` — research mode: force a live solve including the unreached
  0.2/0.1 N rungs, with the honest wall warning up front.

---

## 4. Per-rung recipe map (honest, encoded in the script)

The campaign did not reach all rows equally (README footnotes 1, 2, 6). The
selector maps each thrust to its validated recipe and is honest about limits:

| T [N] | recipe | note |
|---|---|---|
| 10 / 5 / 2.5 | ladder chain (`run_mintime_mee` + `run_transfer_mee`) | clean, cached |
| 1 | ladder + PSR switch-refinement | headline `m_f=1371.44` is PSR round 2 |
| 0.5 | anchor-free R0-law target + PSR | **footnoted:** anchor is an R0-law *estimate*, PSR budget-limited |
| 0.2 / 0.1 | **live probe only** | never certified; 0.5 N min-time hit a conditioning wall — script says so, and reports `certified=false` rather than emitting a fabricated row |

For custom endpoints, the same recipe tiers apply but every rung is a fresh
solve (no cache), and the deep-rung walls may appear at higher thrust.

---

## 5. Trajectory plot + movie  (BUILT + VALIDATED)

`transfer_movie.m` expects a **Cartesian** `res` (9-row `[r;v;m;t;cScale]`,
4-row inertial `[alpha;thr]`), but the MEE solver stores elements
(`[P;ex;ey;hx;hy;m;t]`) and RTN control (`[beta;thr]`). Bridged by:

- **`mee_res_to_cart_res.m`** (already written, `earth_elliptic_to_geo/`) —
  reconstructs inertial `(r,v)` at each node via `elements_to_cart` at
  `L_k = π + σ_k·ΔL`, and rotates RTN `beta`→inertial `alpha = R_{RTN→ECI}·beta`
  using the same per-node triad as `run_transfer_mee`'s reconstruction check.
  Returns the Cartesian `res` layout `transfer_movie` consumes.
  **Validated (2026-07-18):** on the 10 N case, `|r|` runs apogee 1.103 →
  GEO 1.000, `|alpha|≡1`, `m_f=1377.10 kg`, `|z|≤0.049` (7°→equatorial) — all
  correct. Used to render movies for 10/5/2.5/1 N (`results/movie_MEE_*.{mp4,gif}`).
- **Static plot** — the same reconstructed inertial path, throttle-colored
  (burn red / coast blue), GEO ring backdrop, saved PNG.
- **Movie** — `transfer_movie` renderer (gif + mp4, ÷16 frame sizing already
  handled), same house style as the PSR movies.

The script feeds `mee_res_to_cart_res` output to `transfer_movie` (via a
temp `res` mat, or a small refactor letting `transfer_movie` accept a struct
— decided at plan time; the temp-mat path is proven working).

---

## 6. Table-3 row printout

Fixed-width block matching the paper's Table 3 columns:
`T [N] | t_f,min [ND/h] | c_tf | t_f | m_f [kg] | prop [kg] | ΔV [km/s] |
switches | revs (ours vs paper) | edge | incl [deg] | defect | certified`.
Carries the standing footnotes for the 0.5 N (anchor-free) and 0.2/0.1 N
(unreached) rows so a printed row is never mistaken for a certified one.

---

## 7. File structure

```
earth_elliptic_to_geo/
  run_gergaud.m            NEW  front-door PARAMETERS-block script (§3)
  mee_res_to_cart_res.m    DONE MEE→Cartesian movie/plot adapter (§5)
  mee_seed.m               EDIT initial-orbit params (§2)
  casadi_lt_mee.m          EDIT opts.xf terminal target, default GEO (§2)
  homotopy_mee.m           EDIT xf pass-through (§2)
  run_transfer_mee.m       EDIT xf + initial-orbit pass-through, fingerprint (§2)
  run_mintime_mee.m        EDIT xf + initial-orbit pass-through, fingerprint (§2)
  gergaud_plot.m           NEW  static trajectory plot (§5)  [or fold into run_gergaud]
  test_mee_endpoints.m     NEW  default-preserving + custom-endpoint guards
  test_mee_res_to_cart.m   NEW  adapter roundtrip check
```

---

## 8. Open decisions for the plan

1. **`transfer_movie` struct input** — refactor it to accept a struct
   directly, or keep the proven temp-mat handoff? (Lean: small refactor, DRY.)
2. **0.5 N reuse** — `auto` mode reuses `MEE_M2_0p5N` with the footnote; is
   that acceptable, or should 0.5 N always re-probe? (Lean: reuse + footnote.)
3. **Custom-endpoint anchor** — a custom initial/final needs its own
   `run_mintime_mee` anchor (no cached `t_f,min` applies); confirm the anchor
   two-stage recipe is endpoint-agnostic after the §2 edits (it should be —
   the terminal set is the only endpoint dependency, now parameterized).

---

## 9. Non-goals (YAGNI)

- No new solver physics, no new homotopy, no unit rescaling.
- No batch/sweep front end (that is `run_ladder`'s job).
- No attempt to *certify* 0.2/0.1 N here — the script only *probes* them.
