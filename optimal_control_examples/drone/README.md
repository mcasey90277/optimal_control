# drone

Multi-agent drone path planning, from the paper's kinematic drones to a
6-DOF quadrotor, solved directly (CasADi/IPOPT) and indirectly (PMP shot
with `oc.ms_bvp`) and certified. The design is
`docs/superpowers/specs/2026-10-23-drone-multiagent-design.md`; the paper is
Zhang, Liu, Darbon & Karniadakis, arXiv:2610.30388 (`papers/`).

| folder | stage | status |
|---|---|---|
| `drone_dynamics/` | 0: the plant library and shared geometry | built |
| `route1_direct/` | 1: the paper's problem as a hard-clearance NLP (Table 1) | planned |
| `route2_indirect/` | 2: the barrier form's PMP-BVP, seeded from Route 1, conjugate-tested | planned |
| `quadrotor/` | 3: the same task on the rigid-body model | planned |

Run the tests: `cd` here and `run_tests`, or `run_all_tests` from
`optimal_control_examples/`.
