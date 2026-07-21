# lambert — universal-variables Lambert solver (guided tutorial)

Build-it-yourself tutorial for a universal-variables Lambert solver,
including the multi-revolution case. House tutorial format: exercises PDF
(`lambert_exercises.pdf`), reference implementations, a `mytry/` folder for
your own attempt, and verified checkpoints.

| file | what |
|---|---|
| `lambert_uv.m` / `lambert_uv_multirev.m` | Reference solver, single + multi-rev (universal variables, Stumpff functions). |
| `stumpff.m`, `lambert_tof.m` | Stumpff C/S functions; time-of-flight equation. |
| `run_lambert_demo.m` | Demo driver (expected output: `expected_result.png`). |
| `verify_lambert_checkpoints.m` | Checkpoint verifier for the exercises. Validated against pyKep. |
| `mytry/`, `reviews/` | Your working area; external review records. |
