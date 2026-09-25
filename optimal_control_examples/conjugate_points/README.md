# conjugate_points — seeing a conjugate point

An interactive explorer and a guided write-up for the second-order theory of
the simplest calculus-of-variations problem,

    minimise  J[y] = ∫_a^b F(t, y, y') dt,   y(a) = ya,  y(b) = yb.

Type a Lagrangian, and the explorer finds every extremal, flies its
neighbours (same start point, slope p + δ), and shows whether they meet it
again before b. That meeting point, in the limit δ → 0, is the **conjugate
point**; if it lies inside (a, b) the extremal is not a local minimiser.

**Read first:** `doc/seeing_conjugate_points.pdf` (theory in three layers, then
guided experiments with the explorer, each with a checkpoint number).

## Run

```matlab
conjugate_point_explorer          % the GUI
conjugate_point_study             % the guide's numbers + figures, PASS/FAIL per fact
run_conjugate_tests               % the test suite (throws on failure)
```

**Video:** `video/` builds a narrated ~4-minute explainer from the explorer
itself (MATLAB director + AI or draft narration + ffmpeg). Start with
`video/README.md`; the script is `doc/youtube_script.md`.

## Files

| file | what it does |
|---|---|
| `conjugate_point_explorer.m` | the GUI (programmatic `uifigure`, so it diffs and tests like code); six views: curves + neighbours, their difference from the extremal (optionally ÷ δ, which converges to the Jacobi field, drawn dashed), Jacobi field, shooting function, lowest second-variation mode, ΔJ along it; readout with the Legendre, Jacobi and Morse checks and a verdict. Returns a handle API for scripting. |
| `cov_problem.m` | text Lagrangian → symbolic Euler–Lagrange field y'' = g, its Jacobi linearisation, and the second-variation coefficients P, R, Q0 (Symbolic Math Toolbox) |
| `cov_shoot.m` | one extremal from (a, ya) with slope p, flown together with its Jacobi field h = ∂y/∂p; records the zeros of h |
| `cov_extremals.m` | all extremals in a slope range: scan the shooting function r(p) = y(b;p) − yb, refine each sign change |
| `cov_perturbed.m` | the neighbour with slope p + δ and its crossings of the extremal |
| `cov_second_variation.m` | FE discretisation of the second variation; eigenvalues (Morse count), lowest mode η*, exact ΔJ along it |
| `cov_functional.m` | J[y0 + η] by 4-point Gauss quadrature |
| `cov_presets.m` | the seven teaching problems |
| `conjugate_point_study.m` | numbered study script: every fact in the guide computed and checked |
| `run_conjugate_tests.m` + `tests/` | oscillator (closed forms), catenary (closed forms + Lindelöf's tangent construction), cross-instrument checks, mutation tests, GUI |

## What is checked, and against what

Nothing is checked against a copy of itself. The oracles are closed forms
(the oscillator's sin t, cot b, kπ and eigenvalues 2((kπ/b)² − 1); the
catenary's c·cosh(1/2c) = 1), a geometric construction (Lindelöf: the
conjugate point is where two tangents meet the axis together), and a second
instrument (Morse: negative eigenvalues of the second variation = conjugate
points of the Jacobi field). `test_cov_mutation` corrupts the Jacobi
equation, the Legendre coefficient and the integrand in turn and requires the
aimed check to fail and the others to pass.

## Presets

1–3 oscillator F = y'² − y² with b = 3, 4, 7 (0, 1, 2 conjugate points) ·
4 minimal surface F = y√(1+y'²) (two catenaries: one minimiser, one not) ·
5 pendulum F = y'²/2 + cos y (nonlinear: finite δ crosses *near* t_c) ·
6 F = y'² + y² (never a conjugate point) · 7 arc length (straight lines).

Needs the Symbolic Math Toolbox (for `cov_problem`).
