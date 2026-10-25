# lagrangianHamiltonianMechanics

Worked mechanical systems, each in its own folder with a front-door script
that integrates the equations of motion, checks them against an
independent derivation and an exact result where one exists, and writes
static plots and a movie. Shared machinery lives in `mechlib/` (`+mech`).

| folder | system | status |
|---|---|---|
| `mechlib/` | `mech.lagrange_dalembert`, `mech.conservation_check`, `mech.animate_planar_body` | built |
| `nonholonomic/ChaplyginSleigh/` | the Chaplygin sleigh: knife-edge blade offset from the center of mass; `run_sleigh.m` | built |

Run every test with `run_tests` from this folder.
