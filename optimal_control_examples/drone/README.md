# drone

Multi-agent drone path planning, from the paper's kinematic drones to a
6-DOF quadrotor, solved directly (CasADi/IPOPT) and indirectly (PMP shot
with `oc.ms_bvp`) and certified. The design is
`docs/superpowers/specs/2026-10-23-drone-multiagent-design.md`; the paper is
Zhang, Liu, Darbon & Karniadakis, arXiv:2610.30388 (`papers/`).

| folder | stage | status |
|---|---|---|
| `drone_dynamics/` | 0: the plant library and shared geometry | built |
| `quad_flight/` | fly a quadrotor from a user-set ENU/NED state at an LLA origin under a user-set thrust sequence: movie, plots, georeferenced path (`fly_quad.m`) | built |
| `route1_direct/` | 1: the paper's problem as a hard-clearance NLP (Table 1) | built |
| `route2_indirect/` | 2: the barrier form's PMP-BVP, seeded from Route 1, conjugate-tested | planned |
| `quadrotor/` | 3: the same task on the rigid-body model | planned |

Route 1 reproduces Table 1 for M = 4..32 with certified numbers
(`route1_direct/results/table1.csv`); see `route1_direct/README.md` for the
rigorous pair bound that the paper's M = 4 numbers violate.

| M | J_LB | J (ours, hard) | OCFormer |
|---|---|---|---|
| 4 | 0.2000 | 0.200646 | 0.2003 |
| 8 | 0.4000 | 0.404388 | 0.4036 |
| 16 | 0.8000 | 0.816240 | 0.8307 |
| 32 | 1.1556 | 1.207331 | 1.2531 |

Run the tests: `cd` here and `run_tests`, or `run_all_tests` from
`optimal_control_examples/`.
