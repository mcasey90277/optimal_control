## Findings

- **[CORRECTNESS]** `gen_elfo_energy_tfsweep.m:59-63`, `run_elfo_minfuel.m:45,71`, `gen_elfo_minfuel.m:46`, `elfo_export_data.m:138` -- the live pipeline still defines ELFO factors using the tulip anchor `cfg.tfMin=6.2906939607`, despite the certified ELFO anchor being 6.0961534862 ND. A requested “1.20x” ELFO solve is therefore 7.5488 ND rather than 7.3154 ND (3.19% too long), and all front factor labels/seed selection/exported `const.tfMin` are wrong. Replace the shared anchor with an explicitly named ELFO anchor for this layer, re-key/rebuild the factor-indexed seed bank, and retain physical `tf` as the authoritative coordinate.

- **[METHODOLOGY]** `run_elfo_minfuel.m:114-118`, `gen_elfo_minfuel.m:101-104`, `elfo_run_one.m:107`, `elfo_export_data.m:81-88,151-152` -- switch counts are raw `s>0.5` node-crossing integers, embedded in filenames, movie titles, tables, and exports, although ELFO explicitly has no switch-localization/refinement. The supplied campaign uses `N=4000` (`ELFO_RETARGET.md:171`), but records neither revolutions/local physical node density nor a mesh-refinement band; “nodes/rev” is not currently defensibly defined for this nonperiodic transfer. Report counts as mesh-dependent bands after at least one refined-mesh repeat, save physical switch times and local `diff(t)` resolution, and remove the integer from identity-bearing filenames.

- **[ROBUSTNESS]** `gen_elfo_energy_tfsweep.m:134-148` -- a 500-iteration loose solve is banked immediately when only `success && maxDefect<1e-6`; unlike gravhom/fuel, it receives no tight re-clean. Thus the seeds feeding the front can retain merely acceptable IPOPT convergence and inconsistent duals. Run the tight re-clean before `save_point`, and gate it on an explicit certification predicate.

- **[ROBUSTNESS]** `gen_elfo_energy_gravhom.m:187-205`, `gen_elfo_minfuel.m:123-140`, `gen_elfo_energy_tfsweep.m:139-148` -- continuation acceptance checks only solver success and trapezoidal defect. It does not gate endpoint residuals, unit-direction residual, pinned-time residual, IPOPT status class, or proximity to artificial bounds. Centralize a full acceptance gate and checkpoint/bank only certified iterates.

- **[CORRECTNESS]** `gen_elfo_minfuel.m:90-109` -- after a fully resumed homotopy, `finalInfo` is re-solved but never certified, while the file saves stale checkpoint `X,U` alongside `out=finalInfo`. Export uses `out.X,U`; verification and batch metrics use top-level `X,U`, so one result file can describe two different trajectories. Require certification, then assign `Xk=finalInfo.X; Uk=finalInfo.U` before all scalar, count, and save operations.

- **[ROBUSTNESS]** `casadi_energy_freetf.m:163-169,256-264`, `casadi_mintime_freetf.m:100-106,182-189` -- artificial state/control/cScale boxes have no active-bound diagnostics. Defects can be machine-zero while `cScale`, velocity, position, mass, or direction-component bounds constrain the result; only throttle edge fraction is reported. Return per-bound minimum slack/active fractions (and multipliers where available), reject or prominently flag saturation, and widen the already-deferred `r,v` boxes before high-thrust ladder rungs.

- **[ROBUSTNESS]** `gen_elfo_minfuel.m:38,50-68`, `gen_elfo_energy_tfsweep.m:157-188`, `gen_elfo_energy_gravhom.m:113-118`, `gen_elfo_energy_tfsweep.m:163-164` -- checkpoint/seed reuse is not configuration-safe. Fuel resume checks only `tf0`; sweep resume loads only `tf,X,U`; seed files omit propulsion parameters; and insertion variants deliberately share `energy_elfo_f####.mat`. A thrust rung, mesh/clock change, or retarget can silently warm-start from incompatible data. Add a strict fingerprint covering `Tmax,cEx,m0,Isp,rv0,rvf,sigma,tauf,p/q/moonZone`, bounds, and insertion; namespace every rung’s seeds, checkpoints, and results.

- **[EFFICIENCY]** `gen_elfo_mintime.m:47-74` -- the min-time driver prints diagnostics but saves even a failed/unqualified solve. This is harmless for the documented certified anchor, but unsafe when min-time is rerun per thrust rung. Assert the same full certification gate before writing the anchor.

No partial manual constraint scaling is present: both NLPs use IPOPT gradient-based scaling only. No sigma-linear control interpolation exists in this ELFO layer, so the specific phase-aliasing mechanism is not presently triggered; preserve that property when adding mesh transfer.

## Ladder-readiness

A thrust ladder breaks first at provenance: `Tmax` is fixed through `minfuel_config` and incompatible rung seeds/checkpoints can be silently reused. Before starting: **(1)** parameterize thrust and fingerprint/namespace every artifact; **(2)** replace all tulip-factor logic with `tfMin_ELFO=6.0961534862` and rebuild the seed grid; **(3)** add full convergence/bound-saturation/mesh-convergence gates, including a refined-mesh switch-count band.

## Overall verdict

I trust the nominal 25 mN trajectory dynamics and endpoint/defect numbers conditionally, consistent with the prior core review. I do **not** trust the currently published ELFO factor labels/front coordinates or exact switch counts. The highest-value immediate fix is rebasing every live ELFO driver/export/seed key on the certified ELFO min-time anchor.
