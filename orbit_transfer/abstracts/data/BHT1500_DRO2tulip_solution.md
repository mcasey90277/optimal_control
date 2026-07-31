# BHT-1500 DRO→tulip minimum-time solution (2026-07-31)

**Reproduce with:** `bht1500_continuation.m` in this directory. It lives in the
optimal_control tree and writes nothing into pumpkyn or pumpkynPie; it only
cd's there to call `startup()`.

```matlab
addpath('orbit_transfer/abstracts/data');
out = bht1500_continuation();
```

**Verified 2026-07-31**, run fresh from this repo copy: tf_ND agrees to 3.7e-11,
delta-V to 4.4e-05 km/s, propellant to 1.9e-04 kg, costates to 4.5e-07, same
four continuation steps, same zero failures. The small residual differences are
rounding in the values recorded below, not solver drift.

## Thruster — verified against https://www.busek.com/hall-thrusters

| | value |
|---|---|
| model | Busek BHT-1500 |
| thrust | **101 mN** |
| Isp | **1710 s** |
| input power | **1500 W** |
| propellant of record | **xenon** (the page states performance was measured with xenon) |
| mass | not published on that page |

The call notes said "iodine-fueled" and "~6 kg"; **neither is supported by the
published page**. Busek does test iodine, but the 101 mN / 1710 s figures are
xenon numbers. The 1500 W draw was not discussed on the call and is the figure
most likely to be questioned for a 150 kg spacecraft.

## Result

| quantity | value |
|---|---|
| transfer time | **2.7346189145 ND = 12.122 days** |
| delta-V | **0.7204 km/s** |
| propellant | **6.308 kg** (4.21% of 150 kg wet) |
| final mass | 143.692 kg |
| average thrust acceleration | 0.0688 cm/s² |
| terminal residual | 1.61e-09 |

Converged costates (λ_r, λ_v, λ_m):

```
-5.913737  1.796175  11.573373  -0.334882  3.098774  -1.127631  2.376220
tf = 2.7346189145
```

Independently checked: ṁ = T/(Isp·g₀) = 6.0229e-06 kg/s over 1,047,308 s gives
6.3078 kg, and Isp·g₀·ln(150/143.692) = 0.7204 km/s. Both match the solver.

## Versus the previous 70 mN / 900 s case

**31.9% faster and 48.3% less propellant**, for essentially the same delta-V
(0.7204 vs 0.7485 km/s). That is the expected shape: delta-V is set by the
orbital geometry, not the thruster, so nearly doubling Isp roughly halves the
propellant while the extra thrust shortens the spiral.

## How it converged — and the mis-step worth recording

The call's guidance was to walk Isp up gradually (900 → 1000–1200 first) and to
start from the higher-thrust tier. My first attempt walked from the **70 mN**
solution and nearly stalled: `ode45` failed at t ≈ 4.06 with the step size
driven below 1.4e-14, the grazing-singularity mode described on the call, and
the continuation crawled to s = 0.083 after repeated step halvings.

The fix was noticing that **the demo's three commented-out costate blocks are
converged solutions at the three commented thrust tiers** — their final times
(2.709 / 4.015 / 6.189) line up with 0.1 / 0.07 / 0.055 N. Seeding from the
**0.1 N** block makes the continuation almost pure Isp, on a trajectory 32%
shorter and far better conditioned. From there: **four steps, zero failures.**

Practical lesson: for this problem the seed's thrust tier matters more than the
step schedule.

## Caveat on the convergence check

`tfMin` returns no flag, so convergence is verified here by re-propagating and
checking the four conditions the solver enforces (r(tf), v(tf), λ_m(tf), H(tf)).
That re-integration uses `ode45` at default tolerances and is **stricter than
the solver's internal residual**: the published 70 mN baseline re-propagates to
5.0e-06, where `fsolve` reports ~1e-11. The accepted BHT-1500 steps all came in
at ~1.5e-09, so they are tight by either standard.
