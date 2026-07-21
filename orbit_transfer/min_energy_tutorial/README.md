# min_energy_tutorial — minimum-energy orbit transfer (guided tutorial)

Build-it-yourself tutorial for the minimum-energy two-body orbit transfer,
solved BOTH ways — indirect shooting and direct collocation — with
primer-vector verification connecting them. (This folder was the original
`orbit_transfer/` tutorial before the container reorganization.)

| file | what |
|---|---|
| `orbit_transfer_exercises.pdf/.tex` | The guided exercises (theory + build steps + checkpoints). |
| `solve_indirect.m`, `shoot_residual.m` | Indirect route: PMP shooting on the min-energy TPBVP. |
| `collocation_transfer.m` | Direct route: trapezoidal collocation NLP. |
| `primer_check.m` | Primer-vector (costate) verification along a converged trajectory — the pattern later reused by the CR3BP campaigns' certifiers. |
| `ocp_dynamics.m`, `two_body_accel.m`, `gravity_gradient.m` | Shared dynamics pieces. |
| `run_orbit_transfer.m`, `verify_checkpoints.m`, `expected_result.png` | Driver, checkpoint verifier, expected output. |
| `mytry/`, `reviews/` | Your working area; external review records. |
