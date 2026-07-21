# optimal_control

Educational implementations of **optimal control and estimation** in MATLAB:
guided build-it-yourself tutorials, and full research-grade campaigns
reproducing published low-thrust orbit-transfer results by both **direct**
(collocation NLP, CasADi + IPOPT) and **indirect** (PMP shooting) methods.

## Layout

| folder | what |
|---|---|
| `orbit_transfer/` | **The main body of work** — optimal orbit-transfer campaigns, one folder per problem, each split into `direct/` + `indirect/` codebases: elliptic→GEO (HMG-2004 benchmark, certified 10→0.1 N thrust ladder), GTO→tulip and GTO→ELFO (Earth–Moon CR3BP), a shared `cr3bp_common/` library, Lambert + min-energy tutorials, and the source papers. Start at `orbit_transfer/README.md`. |
| `collocation_examples/` | Standalone direct-collocation tutorials (Kelly-style): `ex1_block_move` (minimum-energy point-to-point) and `ex2_cart_pole_swing_up` (underactuated swing-up), trapezoidal transcription via `fmincon`. |
| `mpc/` | Model-predictive control: `mpc_cart_pole` (N=50 receding horizon at 20 Hz) with a step-by-step LaTeX walkthrough. |
| `quasiNewton_matlab/` | Guided BFGS/DFP quasi-Newton optimizer tutorial (`building_a_bfgs_optimizer.pdf` + reference code + `mytry/`). |
| `lieFiltering/` | SO(3) attitude estimation notes: error-state EKF on the rotation manifold + synthetic IMU data generation (LaTeX/PDF). |
| `gauss_sum_curvature/` | Small numerical experiments: Gaussian-sum splitting and SR1 curvature updates. |
| `papers/` | Reference PDFs (optimization, filtering, astrodynamics). |
| `docs/` | Design specs and implementation plans for the larger campaigns (`docs/superpowers/`). |

## Conventions

- **MATLAB R2025b** (`/Applications/MATLAB_R2025b.app/bin/matlab`); the
  orbit-transfer campaigns additionally use **CasADi 3.7.0** (`~/casadi-3.7.0`,
  bundled IPOPT/MUMPS).
- Campaign modules are run from their own folder after calling that module's
  `setup_paths`; each campaign keeps `README.md` + `TODO.md` at its root,
  narratives in `process/`, technical notes in `doc/`.
- `.mat` solver caches are gitignored; committed figures live in each
  campaign's `results/`.
- `CLAUDE.md` files carry the working configuration for AI-assisted sessions.

## References

- Kelly, M., "An Introduction to Trajectory Optimization," SIAM Review, 2017.
- Betts, J., *Practical Methods for Optimal Control Using Nonlinear
  Programming*, SIAM, 2010.
- Haberkorn, Martinon & Gergaud, "Low thrust minimum-fuel orbital transfer,"
  JGCD 27(6), 2004 (the elliptic→GEO benchmark).
