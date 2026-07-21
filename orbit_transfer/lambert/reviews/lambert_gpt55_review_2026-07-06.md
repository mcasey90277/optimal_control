- **[CORRECTNESS]** `lambert_exercises.tex:197` — The negative-`z` Stumpff wording “cos → cosh, sin → sinh, z → -z” gives the wrong sign for `S(z)` if applied literally. Concrete fix: explicitly write `C(z)=(cosh(sqrt(-z))-1)/(-z)` and `S(z)=(sinh(sqrt(-z))-sqrt(-z))/(-z)^(3/2)`.

- **[CORRECTNESS]** `lambert_exercises.tex:132-136`, `392-397`, `414-415` — The multi-rev count omits the tangent case `dt == t_min(N)`: a U-shaped branch has one double solution, not two and not none. The reference code uses `tmin <= dt` and would return two nearly duplicate branch roots. Concrete fix: state “two if `dt > t_min`, one double root if equal, none if less,” and optionally handle equality with one column.

- **[CORRECTNESS]** `lambert_exercises.tex:132-137`, `318-323` — The “any `dt > 0`” / complete retrograde claim is stronger than the provided bracketing code. For the tutorial geometry with `dir=-1, dt=0.1`, `lambert_uv` hangs in the left-bracket loop after `zlo` is pinned near `1e-8`. Concrete fix: add an iteration cap and expand `zlo` leftward when finite `t(zlo) > dt`, or explicitly scope the solver to checkpoint-like times.

- **[PEDAGOGY]** `lambert_exercises.tex:318-323`, `423-426` — “Treat NaN as too long” is safe for multi-rev band-edge divergence, but confusing for single-rev low-`z` infeasibility where `y < 0` is the left invalid domain. Concrete fix: distinguish “invalid because below admissible `y` domain” during bracketing from “too long/divergent near a band edge” during multi-rev bisection.

- **[RUNNABILITY]** `lambert_exercises.tex:360-366` — Checkpoint C2/C3 is not runnable as pasted into a fresh MATLAB workspace; it assumes `r1`, `r2`, and `mu` from C1. Concrete fix: prepend `r1 = [1;0;0]; r2 = [0;1.2;0]; mu = 1;`.

- **[RUNNABILITY]** `lambert_exercises.tex:460-463` — Checkpoint D2 is not runnable as pasted into a fresh workspace; it assumes `r1`, `r2`, and `mu` from D. Concrete fix: include the same three setup lines or say “run immediately after Checkpoint D.”

- **[AUDIENCE]** `lambert_exercises.tex:114-129`, `175-181`, `324-327` — Section 1.2 claims self-contained Lambert theory but assumes unexplained terms: focus, prograde/retrograde, eccentric anomaly, hyperbolic anomaly, conic regimes, transfer plane, universal anomaly, and `f`/`g` functions. Concrete fix: add a small glossary/diagram before the formulas and one sentence explaining what `f` and `g` do before velocity recovery.

- **[PEDAGOGY]** `lambert_exercises.tex:411-418` — Phase D asks for golden-section search and a slope-aware branch bisector but gives no pseudocode; this is a likely stranding point despite the checkpoints. Concrete fix: add a hint with the `gr`, `c`, `d`, update rules, and the exact branch update truth table.

- **[STYLE]** `lambert_exercises.tex:498-523` — “What comes next” introduces many advanced terms without glosses: porkchop plot, ephemeris, primer vector, Lawden boundary conditions, CR3BP, STM, and “tulip world.” Concrete fix: add one short parenthetical definition per term or mark the section as an optional jargon-forward roadmap.
