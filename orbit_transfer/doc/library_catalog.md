# Orbit-transfer library catalog — generated function reference

Generated 2026-09-06 by `gen_library_catalog.py`.
**Do not edit by hand** — regenerate with:
```sh
python3 orbit_transfer/doc/gen_library_catalog.py
```
Per-function headers in the .m files remain the authoritative
documentation (derivations, assumptions, revision history); this
catalog is the browsable index for humans and Claude sessions.
The judgment layer (what to use when, consumers, gates, traps)
lives in each layer's README.md and in `../OCP_UNIFYING_MATH.md`.

## oclib/+oc

Cross-top-level-folder core (consumers in orbit_transfer AND booster_landing). Admission: second top-level consumer + equivalence gate. See ../../oclib/README.md for consumers/gates per function.

### `duals_to_costates.m`
`[lam, tStations, diag_] = duals_to_costates(spec)`  
Maps the KKT multipliers of a direct transcription's DEFECT constraints to samples of the continuous costates -- the covector mapping, made a reusable library function. The principle: stationarity of the NLP Lagrangian with respect to the states IS a discretization of the adjoint equation lambda_dot = -dH/dx, so the defect multipliers sample the costate. Which samples, with what scaling and sign, depends on the transcription and the solver's conventions -- exactly what this function owns, together with the verification that makes the output trustworthy.
*in: `spec`, `scheme`, `mu`, `tNodes`, `tf`, `defectScale`, `uDir`, `velRows`, `nVote`, `lamTf`, `activeFlags` · out: `lam`, `tStations`, `diag_`*

### `fly_control.m`
`[zEnd, out] = fly_control(z0, tGrid, rhs, opts)`  
THE flown-control engine (oclib move 2): integrate a solution's reconstructed control through the true continuous dynamics, end to end, and return where the trajectory actually arrives. This is the physically meaningful accuracy check behind orbit_transfer's G1b gate and booster_landing's G2 gate -- "defect is not accuracy" made executable. The DYNAMICS AND CONTROL RECONSTRUCTION belong to the caller (one rhs closure): reconstruction is domain policy (the booster's annulus-feasible direction/magnitude split is provably NOT a plain quadratic -- see hs_quad_ctrl's ADAPTATION note), so the engine owns only the integration structure.
*in: `z0`, `tGrid`, `rhs`, `opts`, `mode`, `solver`, `RelTol`, `AbsTol` · out: `zEnd`, `out`*

### `local_residual.m`
`[dX, out] = local_residual(X, tGrid, rhs, opts)`  
THE local-residual engine (oclib move 3): the per-interval TRUE continuous-time error of a direct-transcription solution -- "defect is not accuracy" made executable at interval granularity. For every interval, restart the integrator FROM THE TRANSCRIPTION'S OWN LEFT NODE, integrate the true dynamics (with the caller's reconstructed control inside the rhs closure), and report the miss at the right node:
*in: `X`, `tGrid`, `rhs`, `opts`, `solver`, `RelTol`, `AbsTol` · out: `dX`, `out`*

## costate_common

The costate-pipeline library: family construction, multiple shooting, seeds, conjugate test, catalog build/validate/sweep. See its README.md for the judgment layer.

### `assert_periodic_orbit.m`
`ok = assert_periodic_orbit(tau, rv, tol, throwOnFail)`  
PERIODICITY GUARD for a propagated orbit, made a single-home helper (migration #4): an interpolated family seed that cont_np could not truly converge produces a non-closing "orbit" whose downstream metrics are garbage (measured in the survey: periselene below the lunar surface, distances of 1e10 km). Every consumer of get_family_orbit that cannot tolerate a junk orbit should call this.
*in: `tau`, `rv`, `tol`, `throwOnFail` · out: `ok`*

### `build_costate_catalog_family.m`
`cat_ = build_costate_catalog_family(catDir, outMat, spec)`  
FAMILY-AGNOSTIC catalog packager: packages any campaign's thrust-ladder sheets into ONE shareable costate CATALOG in the COMPACT format (data minimization per D. Koblick): only canonical nondimensional quantities are stored -- phase fractions, the z8 vectors (which already contain t_f), thrust rungs, and per-sheet availability/flight-time lookup grids. Everything else (days, Delta-V, masses) is DERIVABLE and the formulas ride along in .derive.
*in: `catDir`, `outMat`, `spec`, `glob`, `name`, `description`, `provenance`, `depReconstruction` · out: `cat_`*

### `catalog_schema.m`
`out = catalog_schema(action, varargin)`  
THE versioned schema authority for compact costate catalogs -- the normative field list, the validator, and the NAMED FORMULA REGISTRY (retiring the accepted debt of `cat.derive` free-form strings being the only statement of the derivations). One home: packagers stamp, pickers and consumers validate, everyone derives through here.
*in: `action`, `varargin` · out: `out`*

### `conj_catalog_pass.m`
`S = conj_catalog_pass(catMat, opts)`  
CONJ_CATALOG_PASS  Run the conjugate-point test over every entry of a compact costate catalog and record the verdicts.
*in: `catMat`, `opts`, `logFile`, `batchSec`, `K`, `maxAtt`, `tolDz`, `wallSec`, `sideMat`, `writeback` · out: `S`*

### `cr3bp_minenergy_pmp.m`
`[F, A, aux] = cr3bp_minenergy_pmp(y, Tmax, c, muStar)`  
The MINIMUM-ENERGY PMP vector field for the CR3BP low-thrust transfer, with its exact 14x14 Jacobian -- the fixed-final-time sibling of pumpkyn.cr3bp.tfMinEoM. Same 14-state PMP y = [r; v; m; lam_r; lam_v; lam_m], same dynamics and costate ODEs, but the running cost is the Bertrand-Epenoy energy endpoint L = s^2 (throttle squared, physical-time measure -- the SAME convention GTO_tulip's energy->fuel homotopy uses at eps = 1), so the throttle is smooth and interior:
*in: `y`, `Tmax`, `c`, `muStar` · out: `F`, `A`, `aux`*

### `cr3bp_minenergy_prop.m`
`[yh, PHI, T, Y] = cr3bp_minenergy_prop(dt, y0, needSTM, Tmax, c, muStar)`  
Propagates the min-energy PMP state (and optionally its 14x14 STM) for a time dt -- the fixed-tf sibling of pumpkyn.cr3bp.tfMinProp, in the [yh, PHI] = prop(dt, y0, needSTM) shape ms_bvp expects. Variational equations dPHI/dt = A(y) PHI ride along as 196 extra states, with A the exact AD Jacobian from cr3bp_minenergy_pmp. Integrator ode113 at tfMinProp's tolerances (RelTol 1e-10, AbsTol 1e-12).
*in: `dt`, `y0`, `needSTM` · out: `yh`, `PHI`, `T`, `Y`*

### `cr3bp_minfuel_pmp.m`
`[F, A, aux] = cr3bp_minfuel_pmp(y, Tmax, c, muStar, smooth)`  
The SMOOTHED ENERGY->FUEL PMP vector field for the CR3BP low-thrust transfer, with its exact 14x14 Jacobian -- the homotopy sibling of cr3bp_minenergy_pmp. Same 14-state PMP y = [r; v; m; lam_r; lam_v; lam_m], same dynamics and costate ODEs; only the running cost L(s) changes, selected by the smoothing family:
*in: `y`, `Tmax`, `c`, `muStar`, `smooth` · out: `F`, `A`, `aux`*

### `cr3bp_minfuel_prop.m`
`[yh, PHI, T, Y] = cr3bp_minfuel_prop(dt, y0, needSTM, Tmax, c, muStar, smooth)`  
Propagates the smoothed energy->fuel PMP state (and optionally its 14x14 STM) for a time dt -- the homotopy sibling of cr3bp_minenergy_prop, in the [yh, PHI] = prop(dt, y0, needSTM) shape ms_bvp expects. Variational equations ride along with A from cr3bp_minfuel_pmp (exact AD). Integrator ode113 at tfMinProp's tolerances (RelTol 1e-10, AbsTol 1e-12).
*in: `dt`, `y0`, `needSTM`, `smooth` · out: `yh`, `PHI`, `T`, `Y`*

### `cr3bp_thrust_rhs.m`
`dz = cr3bp_thrust_rhs(z, u, muStar, Tmax, c)`  
CR3BP dynamics with thrust and mass flow -- the single shared RHS for flown-control verification (flown_control_error, true_min_altitude). Extracted verbatim from certify_dro_mintime/local_f, which itself mirrors dro_residual/local_rhs. One home (migration #4).
*in: `z`, `u`, `muStar`, `Tmax`, `c` · out: `dz`*

### `ctrl_quad.m`
`u = ctrl_quad(ua, um, ub, w)`  
Lagrange quadratic control reconstruction through the node/midpoint/ node samples of a Hermite-Simpson interval, at normalized time w. Extracted verbatim from certify_dro_mintime/local_q (migration #4); shared by flown_control_error and true_min_altitude.
*in: `ua`, `um`, `ub`, `w` · out: `u`*

### `duals_to_costates.m`
`[lam, tStations, diag_] = duals_to_costates(spec)`  
DELEGATE. The covector mapping -- defect-constraint KKT multipliers to continuous-costate samples, with the scheme-specific station association, sign vote, and lambda_t check -- was promoted to the cross-folder optimal-control library on 2026-08-09 and lives at
*in: `spec`*

### `flown_control_error.m`
`[erNd, evNd] = flown_control_error(o, muStar, Tmax, c)`  
THE flown-control verifier (migration #4): flies a direct solution's RECONSTRUCTED CONTROL once, end to end, and reports where the spacecraft actually arrives relative to the solution's own terminal state. This is the physically meaningful accuracy number -- 'if you flew this control, where would you arrive?' -- and it is the G1b gate of every campaign's certification. Extracted verbatim from certify_dro_mintime (the math is family-free: CR3BP + thrust + quadratic control reconstruction).
*in: `o`, `X`, `U`, `Um`, `tNodes`, `s`, `tf`, `muStar`, `Tmax`, `c` · out: `erNd`, `evNd`*

### `get_family_orbit.m`
`[tau, rv, info] = get_family_orbit(family, p)`  
THE one place a costate campaign turns a family name + parameters into a propagated periodic orbit. Every engine (thrust ladders, densifiers, surveys, pickers' examples) builds its endpoints through this helper, so adding a family here makes it available to the whole pipeline -- this is what generalized the DRO->tulip machinery to Halo and beyond.
*in: `family`, `p` · out: `tau`, `rv`, `info`*

### `golden_cells.m`
`ok = golden_cells()`  
GOLDEN-CELL QUALITY REGRESSION for the costate pipeline (principle 7c): fixed benchmark cells with stored reference metrics, where a quality DROP is a failure even when correctness gates still pass. The two subtle bugs this defends against (the Hermite-Simpson midpoint station shift and the missing thrust ratio in continuation mass scaling) both passed every runtime gate because the pipeline is self-healing -- they showed up only in the EFFICIENCY channel (seed quality, iterations).
*out: `ok`*

### `harvest_ms_seed.m`
`[seed, diag_] = harvest_ms_seed(o, K)`  
Builds a multiple-shooting SEED from a direct collocation solution -- the harvest path, made a single-home library function (migration #3; the sign-vote + midpoint-association rules previously lived inline in thrust_ladder_library, the exact one-home-per-rule violation that let the Hermite-Simpson midpoint bug exist in two places).
*in: `o`, `X`, `lamDef`, `Um`, `tNodes`, `tf`, `K` · out: `seed`, `diag_`*

### `ms_bvp.m`
`[p, info] = ms_bvp(prob, seed, opts)`  
GENERIC multiple-shooting two-point BVP engine -- the family- and problem-agnostic core of ms_tfmin, moved to the shared library (migration #3). The arc is split into K segments whose junction states are extra unknowns; short segments kill the Lyapunov amplification that makes single shooting from approximate seeds intractable (measured on the CR3BP min-time problem: collocation seeds miss by 36,000-560,000 km single-shot, converge in a few iterations here).
*in: `prob`, `ny`, `freeIdx0`, `prop`, `rhs`, `terminal`, `seed`, `tf`, `tGrid`, `Y`, `opts` · out: `p`, `info`*

### `ms_conjugate_test.m`
`out = ms_conjugate_test(info, spec)`  
CONJUGATE-POINT TEST on a converged multiple-shooting extremal -- the first piece of second-order optimality checking this pipeline has had. First-order (PMP) conditions admit maxima and saddle extremals too; a conjugate point in (0, tf) means the extremal STOPS being locally minimizing there (Jacobi's necessary condition).
*in: `info`, `spec`, `flow`, `stateRows`, `costateCols`, `quotientDir`, `freeTime`, `rankTol` · out: `out`, `t`, `detScaled`, `sigRatio`, `firstFullRank`, `nCrossings`, `atFinal`, `pass`*

### `ms_tfmin.m`
`[z, info] = ms_tfmin(rv0, rvf, seed, Tmax, c, muStar, opts)`  
MS_TFMIN  Multiple-shooting solve of the CR3BP minimum-time PMP problem.
*in: `rv0`, `rvf`, `seed`, `tf`, `tGrid`, `Y`, `Tmax`, `c`, `muStar`, `opts` · out: `z`, `info`*

### `preflight_screen.m`
`[ok, why, minAltKm] = preflight_screen(o, muStar, lStar, floorKm, seedTf)`  
CHEAP SANITY PRE-CHECK on a direct solution BEFORE any integrator touches it, made a single-home helper (migration #4). A solve can meet the discrete defect test to 1e-9 and still be physically wild; integrating such a trajectory crawls near the lunar singularity WITHOUT BOUND (measured: one cell pinned a catalog run 16 min at 100% CPU). Screens on the discrete nodes only -- altitude floor and a plausible time of flight -- so it costs microseconds.
*in: `o`, `muStar`, `lStar`, `floorKm`, `seedTf` · out: `ok`, `why`, `minAltKm`*

### `run_capped.m`
`[ok, varargout] = run_capped(pool, fcn, nout, capSec, varargin)`  
Runs fcn(args) on a parfeval worker under a HARD wall-clock cap; on timeout or worker error the future is CANCELLED (the worker is killed and restarted), so no single stuck computation can stall the caller.
*in: `pool`, `fcn`, `nout`, `capSec`, `varargin` · out: `ok`, `varargout`*

### `seed_from_z8.m`
`seed = seed_from_z8(z8, rv0, K, Tmax, c, muStar)`  
Build a multiple-shooting SEED from a converged (or trusted) z8 = [lam0(7); tf]: fly it end-to-end with pumpkyn's min-time propagator and cut the flight into K+1 junction states. Seeded AT a root, ms_bvp converges in 1-2 Newton iterations -- the pattern behind the conjugate catalog sweep, the golden-cells regression, and the GTO flagship min-time probe, extracted to one home on its third appearance (migration rule).
*in: `z8`, `rv0`, `K`, `Tmax`, `c`, `muStar` · out: `seed`*

### `ss_bvp_accept.m`
`[z, info] = ss_bvp_accept(prob, y1, tf, opts)`  
GENERIC SINGLE-SHOOTING ACCEPTANCE GATE -- the family- and cost-agnostic form of the catalog pipeline's third gate ("feed the refined entry to the single-shooting solver; PASS = it comes back unchanged in a few iterations"). For min-time the independent solver is pumpkyn's tfMin; for every other cost (min-energy, min-fuel) there is no pumpkyn twin, so this harness plays the role on the same three-closure prob that ms_bvp consumes. Mechanically it IS ms_bvp with K = 1 -- one full-arc propagation, unknowns y1(freeIdx0) (+ tf when free), terminal residual only -- which is exactly single shooting; the value added is the acceptance semantics: the residual AT the seed before any step (normR0: "is it already a root?"), the distance moved (dz), and the verdict accepted = converged AND dz < tolDz.
*in: `prob`, `y1`, `tf`, `opts` · out: `z`, `info`*

### `survey_family_bounds.m`
`B = survey_family_bounds(family, paramGrid, outMat)`  
Finds the "REASONABLE" members of ANY orbit family, by Darin's criteria: periselene altitude >= 500 km (no lunar impact, with the same margin the transfer campaigns use) and the whole orbit within 100 Mm of the Moon (lunar vicinity, far enough out to take in L1/L2). Pure propagation via get_family_orbit -- no optimization.
*in: `family`, `paramGrid`, `outMat` · out: `B`, `family`, `rows`, `admissible`*

### `true_min_altitude.m`
`amin = true_min_altitude(o, muStar, Tmax, c, lStar, rMoonKm)`  
Minimum lunar altitude of the PROPAGATED trajectory, not of the nodes. A collocation altitude floor binds at nodes only; this checks it BETWEEN nodes, where periselene actually happens. Extracted verbatim from certify_dro_mintime/local_true_min_alt (migration #4).
*in: `o`, `muStar`, `Tmax`, `c`, `lStar`, `rMoonKm` · out: `amin`*

**tests/**: `test_catalog_schema_v3.m`, `test_conj_fixedtf.m`, `test_cr3bp_minenergy_pmp.m`, `test_gto_family.m`, `test_huber_saltation.m`, `test_minfuel_pmp.m`, `test_ms_bvp_fixedtf.m`, `test_ss_bvp_accept.m`

## verify_common

First-order optimality gate layer + the shared continuous-residual (G1) gate. See its README.md and OPTIMALITY_CERTIFICATION.md.

### `certified_guard.m`
`info = certified_guard(res, spec)`  
CERTIFIED_GUARD  Refuse to verify against a re-solve that is not the certified point.
*in: `res`, `success`, `maxDefect`, `value`, `spec`, `caller`, `label`, `saved`, `name`, `errName`, `better`, `feasTol`, `tol`, `units` · out: `info`*

### `foc_check.m`
`rep = foc_check(out, sigma, man, opts)`  
FOC_CHECK  Generic AD-based first-order optimality (PMP/KKT) gate.
*in: `out`, `X`, `U`, `sigma`, `man`, `opts`, `tolStat`, `tolSign`, `sdotMin` · out: `rep`, `kktStatInf`, `sLag`, `dirTanMax`, `dirTanMed`, `signPct`, `Sd`, `lam`, `lamTimeCoV`, `lamTimeEnd`, `sdotMinRel`, `nSwitches`, `horizonNote`, `checksRun`, `pass`*

### `foc_dual_to_costate.m`
`lam = foc_dual_to_costate(LamDef, sigma)`  
FOC_DUAL_TO_COSTATE  Interval defect duals -> nodal costate, step-weighted.
*in: `LamDef`, `sigma` · out: `lam`*

### `foc_ipopt_inertia.m`
`ic = foc_ipopt_inertia(regHistory, opts)`  
FOC_IPOPT_INERTIA  Local-minimality verdict from IPOPT's NATIVE inertia (delta_w).
*in: `regHistory`, `opts`, `tailN`, `tol` · out: `ic`, `certLocalMin`, `maxTailDw`, `tail`, `nIter`, `verdict`*

### `foc_manifest.m`
`man = foc_manifest(name)`  
FOC_MANIFEST  Campaign-level state and control dimension registry.
*in: `name` · out: `man`, `name`, `nx`, `nu`, `dirRows`, `thrRow`, `massRow`, `timeRow`, `autonomous`, `horizonKind`, `massFreeAtTf`*

### `foc_report.m`
`foc_report(rep, tag, resDir)`  
FOC_REPORT  Standard first-order optimality report (print + optional sidecar save).
*in: `rep`, `kktStatInf`, `sLag`, `dirTanMax`, `dirTanMed`, `signPct`, `Sd`, `lam`, `lamTimeCoV`, `lamTimeEnd`, `derivedFromKKT`, `sdotMinRel`, `nSwitches`, `horizonNote`, `checksRun`, `pass`, `tag`, `resDir`*

### `mee_residual.m`
`R = mee_residual(o, par, sigma, opts)`  
MEE_RESIDUAL  True continuous-time (continuous-longitude) local error of a direct MEE solution -- the shared G1 gate for BOTH MEE campaigns (earth 2-body and CR3BP-GEO), routed through the library engine oc.local_residual. Moved here from earth_elliptic_to_geo/direct/verify 2026-08-26 when the CR3BP campaign became the second consumer (migration rule). The lunar third-body term rides in transparently: pass a par carrying .pert (the verify_cr3bp_pmp fingerprint pattern) and lt_mee_rhs's opt-in branch handles it -- this gate never looks.
*in: `o`, `par`, `sigma`, `opts` · out: `R`*

### `setup_verify_common.m`
`setup_verify_common()`  
SETUP_VERIFY_COMMON  Put orbit_transfer/verify_common on the MATLAB path. Self-contained: no campaign paths, no CasADi (callers add CasADi themselves). OUTPUTS: none (path side effect)

**tests/**: `test_certified_guard.m`, `test_foc_check_10N.m`, `test_foc_check_toy.m`, `test_foc_dual_to_costate.m`, `test_foc_ipopt_inertia.m`, `test_foc_manifest.m`, `test_foc_mesh_invariance.m`, `test_foc_report.m`, `test_foc_terminal_covector.m`

## cr3bp_common

Shared CR3BP GTO problem definition (params, endpoints, setup).

### `check_cr3bp_fp.m`
`check_cr3bp_fp(Scached, fpNow, file, tag)`  
CHECK_CR3BP_FP  Fail-loud cache-fingerprint guard (earth-campaign pattern).
*in: `Scached`, `fpNow`, `file`, `tag`*

### `cr3bp_fingerprint.m`
`fp = cr3bp_fingerprint(p, extra)`  
CR3BP_FINGERPRINT  Build the config fingerprint that determines a solution.
*in: `p`, `extra` · out: `fp`*

### `cr3bp_ipopt_opts.m`
`p = cr3bp_ipopt_opts(maxIter, warmTight)`  
CR3BP_IPOPT_OPTS  The CR3BP-family IPOPT option set (Sundman solvers).
*in: `maxIter`, `warmTight` · out: `p`*

### `cr3bp_lt_params.m`
`p = cr3bp_lt_params(thrust_N, m0_kg, Isp_s)`  
CR3BP_LT_PARAMS  Earth-Moon CR3BP + low-thrust nondimensional parameters.
*in: `thrust_N`, `m0_kg`, `Isp_s` · out: `p`, `thrustN`, `ispS`, `g0`, `c`, `Tmax`*

### `gto_elfo_endpoints.m`
`[rv0, rvf, elfoTrace] = gto_elfo_endpoints(p, opts)`  
GTO_ELFO_ENDPOINTS  Boundary states for the GTO -> lunar-ELFO transfer.
*in: `p`, `opts`, `oe`, `raan`, `M0`, `ref` · out: `rv0`, `rvf`, `elfoTrace`*

### `gto_tulip_endpoints.m`
`[rv0, rvf, tulipTrace] = gto_tulip_endpoints(p)`  
GTO_TULIP_ENDPOINTS  Boundary states for the GTO -> south-pole tulip transfer.
*in: `p` · out: `rv0`, `rvf`, `tulipTrace`*

### `insertion_states.m`
`[rv0, rvf, meta] = insertion_states(target, criterion)`  
INSERTION_STATES  Single source of truth for the GTO departure (rv0) and the tulip/ELFO insertion (rendezvous) state (rvf) used by every low-thrust pipeline. Declaring endpoints here (instead of threading them implicitly from a seed .mat) makes them explicit, changeable, and drift-checkable.
*in: `target`, `criterion` · out: `rv0`, `rvf`, `meta`*

### `minfuel_config.m`
`cfg = minfuel_config(over)`  
MINFUEL_CONFIG  Single source of truth for the min-fuel GTO->tulip campaign.
*in: `over` · out: `cfg`, `tfMin`, `pSund`, `thrustN`, `schedSharpen`, `maxIter`, `dirs`, `fname`, `fparse`*

### `setup_cr3bp_common.m`
`setup_cr3bp_common()`  
SETUP_CR3BP_COMMON  Add the shared CR3BP GTO-transfer library + pumpkyn.

### `test_cr3bp_fp.m`
(a) match -> silent

### `test_minfuel_config_override.m`
`s = strip_handles(s)`  

### `thrust_tag.m`
`tag = thrust_tag(thrustN)`  
THRUST_TAG  Artifact filename token for a thrust rung.

**tests/**: `test_cr3bp_ipopt_opts.m`, `test_insertion_states.m`, `test_solver_fork_parity.m`

