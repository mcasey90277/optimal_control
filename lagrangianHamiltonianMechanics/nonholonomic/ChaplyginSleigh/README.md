# Chaplygin sleigh

A planar body with a knife-edge blade a signed distance `a` from its center
of mass: the blade forbids sideways motion, a nonholonomic constraint. The
first system of the `lagrangianHamiltonianMechanics` series; the
derivation, conventions and the exact solution are in
`doc/chaplygin_sleigh_note.pdf`.

## Run it

Open `run_sleigh.m`, edit section 1 (`m`, `J`, `a`, `v0`, `w0`, `tEnd`,
`dtOut`, movie options incl. `movieOpts.follow`, the following camera's
half-height in metres; `[]` gives a fixed whole-path camera), run. It prints five diagnostics with PASS/FAIL,
writes `results/<label>_path.png` (whole path plus a zoom on the transient with bodies and heading ticks), `results/<label>_diag.png` and
`results/<label>.mp4` (1280 x 720). `results/` is git-ignored.

Try `a = -0.2` (or `v0 = -1`): the center of mass trails the blade, the
straight ride is unstable to yaw, |w| grows until the signed speed crosses
zero, and the sleigh departs the other way. With `w0 = 0` nothing happens
in either case: a straight ride is an exact equilibrium, stable or not.

## What is checked, every run

| line | what | reference |
|---|---|---|
| D1 | the hand-derived reduced equations (`sleigh_rhs`) against the symbolic Lagrange-d'Alembert builder (`mech.lagrange_dalembert` on `sleigh_lagrangian`) | independent derivation |
| D2 | v(t) and w(t) against the exact solution v = S tanh z, w = sign(w0) (S/b) sech z, z = asinh(v0/(b abs(w0))) + k S t | closed form |
| D3 | energy E = m/2 v^2 + (J + m a^2)/2 w^2 | conserved |
| D4 | constraint residual A qd on the reconstructed center-of-mass state (also in the title of the energy-split panel of the diagnostics figure) | a consistency check of the mapping (it is zero by construction); D1-D3 are the independent evidence |
| D5 | the asymptote sign(a) sqrt(C), the decay of w, the total heading change sign(w0) (b/abs(a)) acos(sign(a) v0/S) | stability statement |

Straight motion is a LINE of equilibria (v, 0): the one with the center of
mass leading (a v > 0) is Lyapunov stable and attracts transversely (on its
energy ellipse it is asymptotically stable); the trailing one (a v < 0) is
unstable to yaw. No energy is dissipated: the constraint force does no
work, yaw energy is traded into translation (and the other way first, on a
trailing launch).

## Conventions (binding)

World x east, y north, theta counterclockwise from x; e1 along the blade
(forward), e2 to the left. q = (x_c, y_c, theta) of the center of mass;
blade point r_b = r_c - a e1; `a > 0` = center of mass ahead of the blade.
Reduced state integrated: `[x_b; y_b; theta; v; w]` of the blade point.

## Files

| file | what |
|---|---|
| `run_sleigh.m` | the front door |
| `sleigh_params.m` | m, J, a and the derived Jb, k |
| `sleigh_lagrangian.m` | symbolic L and A for the builder |
| `sleigh_rhs.m` | the hand-derived reduced equations |
| `sleigh_quantities.m` | paths, energy, blade force, exact solution, full-state arrays |
| `sleigh_body.m`, `sleigh_arrows.m` | the drawing and the movie's arrows |
| `plot_sleigh.m` | the two static figures |
| `tests/test_sleigh_oracle.m` | D1-D5 for both signs of a and v0, the knife edge, the straight ride, rest |
| `tests/test_sleigh_pipeline.m` | plots, movie and the script's parse, headless |
| `doc/chaplygin_sleigh_note.tex` | the derivation |

Tests run through `run_tests` in `lagrangianHamiltonianMechanics/`.
