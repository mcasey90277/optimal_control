# mechlib: the mechanics library (`+mech`)

Shared machinery for the Lagrangian/Hamiltonian mechanics systems in this
tree (the Chaplygin sleigh first; rolling disk, rolling coin and knife edge
to follow). MATLAB package form: `addpath('.../mechlib')` once, then call
`mech.<function>`. It is a sibling of `oclib`, not part of it: oclib holds
optimal-control structure and stays physics-free, and when an
optimal-control problem on one of these systems appears it will consume
this library as a plant.

**Admission rule:** a function enters mechlib when it is written without
any system-specific quantity (no mass, length or named coordinate of a
particular system in its body) and carries a test on a problem with an
exact answer. System folders keep nothing generic.

## Contents

| function | what | exact-answer test |
|---|---|---|
| `mech.lagrange_dalembert` | Lagrange-d'Alembert equations with multipliers from a symbolic Lagrangian L(q, qd), constraint matrix A(q) (A qd = 0) and generalized forces Q; returns numeric handles (RHS, multipliers, energy, constraint) and the symbolic reduced equations | free particle, damped particle, pendulum as a constrained particle (lambda = -T/l, period), knife-edge particle (lambda = m v w) |
| `mech.conservation_check` | energy drift and constraint residual along a solution, with bars | exact circular motion passes; injected energy ramp and constraint violation fail with the right reason |
| `mech.animate_planar_body` | top-view MP4 of a planar rigid body: pose history, outline and markers in body coordinates, arrows and subtitle as functions of time, growing traces, fixed limits, exact 1280 x 720 frames | a 20-frame movie exists with 1280 x 720 frames; a single-point pose still gets finite limits; unsorted time refused |

## Tests

`tests/` holds one test per function; each returns `ok`, prints PASS/FAIL
lines with the measured number and its bar, and runs from a fresh clone.
Run them all with `run_tests` from `lagrangianHamiltonianMechanics/`.

## House style

Header quartet, `nargin == 0` self-demos, never `i`/`j`, never `norm`, no
Code Analyzer pragmas. Symbolic Math Toolbox is used by
`lagrange_dalembert` only.
