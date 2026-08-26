# MfMax v1 — runbook

How to build, run, drive and read **MfMax v1**, the Gergaud-group Fortran
indirect solver (`mfmax-v1/`). Everything here was established by building and
running it on this machine on **2026-08-15** (macOS 26, gfortran 13.1.0); the
numbers in [Reference run](#reference-run-10-n) are measured, not quoted.

Presentation forms of this same material: **`../../doc/mfmax_v1_howto.tex`**
(step-by-step how-to, 9 pp) and **`../../doc/mfmax_v1_slides.tex`** (beamer
deck, 17 slides). This file stays the source of record; regenerate those two
with `pdflatex` if you change anything here.

Companion: `mfmax-v0/` is the earlier, time-domain variant — see
[v0 vs v1](#v0-vs-v1). Theory: `mfmax_docs/mfmax.pdf`, `MfMaxmethod.pdf`.

---

## 1. What v1 solves

Minimum-fuel low-thrust orbit transfer (the HMG-2004 elliptic→GEO benchmark),
by **single shooting on the PMP boundary-value problem + differential
homotopy** (HOMPACK90 `fixpqf` curve tracking), in **modified equinoctial
elements**.

The v1 device — this is the whole point of v1 versus v0:

- the **independent variable is the longitude L**, integrated `L0 → L_f`;
- **`L_f` is FIXED** by input, as `L_f = L0 + par(1)·(Lmin − L0)`;
- **`t_f` is FREE**, carried as extra state `x(8)` with `dx(8)/dL ≡ 0`
  (`Dhfun`: `F(8) = 0`), so it is a constant recovered by the shooting;
- **normalized time `s = t/t_f` is state `x(6)`**, with terminal condition
  `s(L_f) = 1` closing the system.

State/costate vector `Y(1:16)`, all in **Mm, hours, kg**:

| index | 1 | 2 | 3 | 4 | 5 | 6 | 7 | 8 | 9–16 |
|---|---|---|---|---|---|---|---|---|---|
| meaning | `P` | `e_x` | `e_y` | `h_x` | `h_y` | `s = t/t_f` | `m` | `t_f` | costates `p_1..p_8` |

Note the aliasing inside `Phifun`: state slot 6 *as passed to `Ufun`/`Dhfun`*
is overwritten with the current longitude `l`, while the *integrated* slot 6 is
`s`, handed over as the separate scalar argument `t`. `Phifun` then rescales
every derivative by `dt = 1/F(6) = ds/dL` to change the independent variable
from `s` to `L`. Read `ftir.f90:74-102` before touching the dynamics.

Terminal conditions (`B2fun`, 8 of them, matching 8 shooting unknowns):

| # | condition | meaning |
|---|---|---|
| 1–5 | `Y(1:5) = par(5:9)` | target orbit MEE — GEO: `P=42.165 Mm, e_x=e_y=h_x=h_y=0` |
| 6 | `Y(6) = 1` | elapsed time equals `t_f` |
| 7–8 | `Y(15) = 0`, `Y(16) = 0` | transversality: free final mass, free final time |

(The two comments on `B2fun` lines 27–28 label these `ptf`/`pm` in the opposite
order to the indices. Both are zero, so nothing is wrong numerically — but do
not trust the labels.)

### The two homotopies

`lambda` is a **pair**, and the two entries do completely different jobs:

| | driver | 0 means | 1 means | run by |
|---|---|---|---|---|
| `lambda(1)` | initial-condition homotopy | start *at the target orbit* (trivial problem) | the true initial orbit | `homCI`, discrete continuation, step `par(10)`, halving down to `par(11)` before giving up |
| `lambda(2)` | **energy → fuel** control regularization | smooth (energy-like) control | true **bang-bang min-fuel** | `Fullpath`, HOMPACK differential path-following |

`Ufun` shows the regularization directly: with `coeff = λ₂ − β·T_max·p_m` and
`coeffc = 2(1−λ₂)`, the control is `0` / linear-ramp / unit-vector in the three
regimes. At `λ₂ = 1` the ramp width `coeffc` collapses to zero and the control
is exactly bang-bang. This is the same energy→fuel idea as our own ε-homotopy
in `GTO_tulip`, done by arclength continuation instead of a discrete ε ladder.

---

## 2. Build

The shipped `src/makefile` still names **`g95`** (dead) and there is no
committed binary (`bin/` is empty). Override the compiler — same recipe that
worked for v0:

```sh
export SDKROOT=$(xcrun --show-sdk-path)
make -C src F90=gfortran \
     OPTF="-O3 -c -std=legacy -fallow-argument-mismatch" \
     OPTO="-O3 -o"
```

Result: 8 objects + `bin/path`. Expect only noise warnings — clang deployment
version override, duplicate `-lemutls_w`/`-lgcc`, and "object file built for
newer macOS". No source edits are needed.

`-std=legacy -fallow-argument-mismatch` are **required**: HOMPACK90 and the
LAPACK/rkf45 slices are F77 with mismatched argument types that gfortran ≥10
rejects by default.

`src/` in the repo also carries 8 stale `.o` and 4 `.mod` files from a 2007
build. Delete them before building (`rm -f src/*.o src/*.mod`) or the makefile
may skip recompiles.

### Build somewhere disposable

Build artifacts are committed to git in this folder, so prefer copying the tree
to a scratch directory and building there:

```sh
cp -R mfmax-v1 /tmp/mfmax-v1 && cd /tmp/mfmax-v1
rm -f src/*.o src/*.mod
```

---

## 3. Run

`path` takes **no arguments**. It reads `./in.dat` from the current directory
and writes `./out.dat`, `./next.dat` (and `fort.<TRACE>`) beside it.

```sh
mkdir run10N && cd run10N
cp ../matlab/IN/in10N.dat in.dat
../bin/path | tee run.log
```

The shipped `in10N.dat` is a complete, working 10 N case. Runtime **0.7 s**.

To continue from a converged solution, `next.dat` is written in exactly the
`in.dat` format with the solution as the new initial guess and `lambda = (1,1)`:
`cp next.dat in.dat`, edit what you want to change, rerun.

---

## 4. `in.dat` format

Twelve list-directed records, read by `Init` (`ftir.f90:320`). Values may wrap
across lines. `n = 8`, `lpar = 11`, `lipar = 1`.

| # | contents | 10 N case |
|---|---|---|
| 1 | `ZI(1:8)` — initial guess for the shooting unknowns | `10. 0. 0. 0. 0. 0. 0. 0.` |
| 2 | `FREE0(1:8)` — which `Y0` slots the unknowns occupy | `8 9 10 11 12 13 14 15` |
| 3 | `Y0(1:16)` — full initial state+costate; the `FREE0` slots are overwritten by `ZI` | `11.625 0.75 0. 0.0612 0. 0. 1500. 0.` + 8 zeros |
| 4 | `L0` — initial longitude [rad] | `3.14159` |
| 5 | `Lmin` — reference final longitude (the min-time solution's) | `29.901174314034` |
| 6 | `lambda(1:2)` — homotopy start | `0. 0.` |
| 7 | `PAR(1:11)` — see below | `2. 10. 0.0142 0.516586E+04 42.165 0. 0. 0. 0. 0.1 0.01` |
| 8 | `IPAR(1)` — run mode | `1` |
| 9 | `NIT` — trajectory samples written to `out.dat` | `1000` |
| 10 | `jac_step` — finite-difference step for the shooting Jacobian | `1.0E-4` |
| 11 | `TRACE` — **a Fortran unit number**, not a verbosity level | `2` |
| 12 | `SSPAR(1:2)` — min / max HOMPACK prediction steplength (0 = auto) | `0. 0.` |

`FREE0 = [8..15]` means the eight unknowns are `t_f` and the seven costates
`p_P, p_ex, p_ey, p_hx, p_hy, p_s, p_m` at `L0`; `p_tf(L0) = Y0(16)` stays
fixed at 0. `ZI(1)` is therefore the **`t_f` guess in hours** — `Fullpath`
forces it to 10 if you pass `≤ 0`.

### `PAR`

| slot | meaning | 10 N value |
|---|---|---|
| `par(1)` | **longitude-span multiplier**: `L_f = L0 + par(1)·(Lmin − L0)`. The v1 analogue of v0's `c_tf`. | `2.0` |
| `par(2)` | max thrust [N] (internally `× tmaxconv = 12.96` → kg·Mm/h²) | `10.` |
| `par(3)` | `β`, mass-flow coefficient [h/Mm]: `ṁ = −β·T_max·‖u‖` | `0.0142` |
| `par(4)` | `μ` [Mm³/h²] | `5165.86` |
| `par(5:9)` | **target orbit** `(P, e_x, e_y, h_x, h_y)` — used both as the terminal constraint and as the `λ₁ = 0` easy start | `42.165 0 0 0 0` |
| `par(10)` | initial `λ₁` continuation step | `0.1` |
| `par(11)` | minimum `λ₁` step before declaring failure | `0.01` |

`par(3) = 0.0142 h/Mm` ⇒ `v_e = 70.4225 Mm/h = 19 561.8 m/s` ⇒ **Isp = 1994.75 s**
— the benchmark's exact 1994.8 s (Caillau & Noailles 2001), not the round 2000 s
our direct campaign defaults to. Worth remembering when comparing masses.

### `IPAR(1)` — run mode

| value | behaviour |
|---|---|
| `≤ 1` | full homotopy (`Fullpath`): `λ₁` continuation to 1, then HOMPACK on `λ₂` to 1, then refine, then write files. **The normal mode.** |
| `2` | no solve — integrate the BVP from the given `ZI` and write the trajectory (`Term`). Use with `next.dat` to redraw a stored solution. |
| `3` | single shooting at the given `lambda` only (`Ssolve`), print residual, then write the trajectory. |

---

## 5. Outputs

| file | contents |
|---|---|
| `out.dat` | header block (`FREE0`, `Y0`, `L0`, `LF`, `LAMBDA`, `PAR`, `IPAR`, `NIT`, solution `Z`, residual `S`, `|S|`) then `NIT+1` trajectory records |
| `next.dat` | restart file in `in.dat` format, solution as the new `ZI`, `lambda = (1,1)` |
| `fort.2` | HOMPACK arclength trace (unit = the `TRACE` value) |

**`out.dat` header gotcha:** the `# LF =` line prints `Lmin`, **not** the
actual final longitude. The real one is `L0 + par(1)·(Lmin − L0)`, and it is
the last value in column 1 of the data block.

**`out.dat` parsing gotcha:** each trajectory record is 20 numbers
(`L`, `Y(1:16)`, `u(1:3)`) but the write format `'(e24.16,2x,(e24.16e3))'`
triggers Fortran format reversion, so **one record is spread over 20 lines**.
Do not read it line-by-line. Slurp all tokens after `# Data` and reshape:

```python
vals = [float(t) for l in open('out.dat').read().split('# Data')[1].split('\n')
                 for t in l.split()]
rows = [vals[k:k+20] for k in range(0, len(vals), 20)]     # L, Y(1:16), u(1:3)
```

---

## 6. Reference run (10 N)

Stock `in10N.dat`, unmodified. Reproduce this before trusting any new case.

**Setup:** 1500 kg from `P = 11.625 Mm, e = 0.75, i = 7°` (`h_x = 0.0612`) at
`L0 = π` to GEO `P = 42.165 Mm` circular equatorial; `T_max = 10 N`;
`L_f = π + 2.0·(29.9012 − π) = 56.6608 rad` fixed; `t_f` free.

| quantity | value |
|---|---|
| wall time | 0.70 s (0.47 s in the differential homotopy) |
| Jacobian evaluations | 18 |
| shooting residual `‖F‖` | 2.24e-9 |
| **final mass** | **1381.4074 kg** (118.5926 kg propellant) |
| **final time** (the free unknown) | **146.3947 h = 6.0998 d** |
| revolutions | 8.5178 |
| thrust structure | 18 arcs / 17 switches, 19.0 % burn duty |

Independent checks on the written trajectory (not printed by the code):

- `P(L_f) = 42.164999999` Mm — **8.9e-10 Mm** from target
- `e = 1.4e-11`, `h_x, h_y ~ 1e-13` → circular equatorial reached
- `‖u‖ ∈ {0, 1}` exactly → genuine bang-bang, `λ₂ = 1` really was attained
- `x(8)` constant across all 1001 samples, spread **0.00e+00** → the free-`t_f`
  device holds to machine zero
- `s` monotone `0 → 1`; mass monotone decreasing

Plausibility against our own work: the certified direct campaign gives
1377.10 kg at 7.326 rev (`c_tf = 1.5`) and its `c_tf` front runs 1377–1385 kg
over `c_tf` 2.0–2.5. 1381.41 kg at 8.518 rev sits inside that band, in the
right direction (longer transfer → less propellant).

---

## 7. Gotchas

1. **`Pfun` ignores integration failure.** rkf45 `iflag = 6` means *"integration
   was not completed — accuracy unreachable at the smallest stepsize"*;
   `Pfun` prints `[PFUN] **** IFAIL (RKF45) = 6` and then uses `Y` as though it
   were at `L_f`. The reference run fires it **32 times**. All 32 are inside the
   `λ₁` continuation / finite-difference Jacobian phase, none in the final
   refinement, so the accepted solution is clean — but on a harder case the
   solver could converge onto truncated integrations silently. Always check
   where the `IFAIL` lines sit relative to the `BEFORE/AFTER REFINEMENT` blocks.
   (`Term` *does* check, and would print `Solution file may be incorrect`.)
2. **`TRACE` is a unit number.** `TRACE = 2` writes `fort.2`. `TRACE = 0`
   disables tracing; `TRACE > 1` also turns on verbose stdout.
3. **Two silent clamps** in `Phifun`: `x(1) = max(1e-3, P)` and `dt = 0` when
   `F(6) ≤ 0` (non-increasing longitude). Both hide non-physical excursions
   rather than erroring.
4. **`bin/` is empty and `src/*.o`, `src/*.mod` are committed stale.** Build
   fresh, ideally out of tree.
5. **`matlab/main.m` is not part of the solver** — an old NN-toolbox script
   (`initff`/`trainbp`/`simuff`) that trains a 5-neuron network, presumably the
   authors' `t_f` estimator. The MATLAB post-processing entry points are
   `mfmax.m`, `menu1..5.m`, `drawpath.m`, `drawres.m`; `evolution.m` and
   `Lfmin.m` are v1-specific.
6. **`Lmin` is an input you must supply per case.** For a new thrust level the
   shipped `29.901174314034` is wrong — it is the min-time final longitude for
   10 N. Getting `Lmin` for another `T_max` is its own (min-time) problem.

---

## 8. v0 vs v1

|  | v0 | v1 |
|---|---|---|
| independent variable | time `t` | longitude `L` |
| state dimension `n` | 7 | 8 (`t_f` added) |
| fixed | `t_f` **and** `L_f` | `L_f` only |
| free | — | `t_f` |
| span control | `par(1) = c_tf` on `[t0, tmin]` | `par(1)` on `[L0, Lmin]` |
| `lpar` / `lipar` | 14 / 2 | 11 / 1 |
| shooting tolerances | `arcre/ae 1e-5`, `ansre/ae 1e-8` | `1e-6`, `1e-10` (tighter) |
| status here | ported + validated 2026-07-20: 10 N, `c_tf = 1.5` → **1378.37 kg** | ported + validated 2026-08-15: 10 N, `par(1) = 2.0` → **1381.41 kg**, `t_f` free = 146.39 h |

`hybrd.f`, `rkf45.f` and `LAPACK.f` are byte-identical between the two; the
solver core is the same, only the formulation differs.

**These two runs are not the same problem instance** — the fixed/free roles are
swapped and the span multipliers act on different axes. Do not read
1381.41 vs 1378.37 as a discrepancy.

---

## 9. Open cross-check

The clean apples-to-apples test has **not** been run: take v1's converged
answer (`t_f = 146.3947 h`, `L_f = 56.6608 rad`, `m_f = 1381.4074 kg`) and give
v0 the same transfer with `t_f` fixed at 146.3947 h. If v0 returns
`m_f = 1381.41 kg` and `L_f = 56.66 rad`, the two formulations are validated
against each other and against our direct campaign. Blocker is only clerical:
v0's `lpar = 14` `PAR` layout has to be decoded the same way §4 decodes v1's.

Second open item: MfMax's converged costates (`out.dat` `# Z`) are a candidate
clean seed source for our own MATLAB indirect solver, per
`indirect/README.md` — the direct-KKT costate seeds are unreliable at high
eccentricity (`../process/DESIGN_dual_map.md`).
