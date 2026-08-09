### Task 11: Phase 2 — atmosphere + drag

**Files:**
- Modify: `run_booster_landing.m` (add `cfg.phase2` switch, default false)
- Create: `tests/test_phase2_drag.m`
- Create: `viz/plot_vacuum_vs_drag.m`

**Interfaces:**
- Consumes: everything; `pdg_dynamics` drag branch already exists and is Jacobian-tested (Task 2).
- Produces: `R = run_booster_landing(struct('phase2', true))` additionally runs:
  - drag-on collocation solve (`P.drag.on = true`), **warm-started** from the vacuum solution (`opts.init = R.solC`)
  - `certify_pdg(solD, [], P)` — G1/G2/G5 gates only (no convex twin; G3/G4 'skipped')
  - drag-aware TVLQR + closed loop + MC **with wind** (sig.wind active because `P.drag.on`)
  - `plot_vacuum_vs_drag(solC, solD, outfile)` — overlay: trajectory, throttle, mass; annotation box with Δfuel and Δtf (the "what drag-free guidance misses" number, headed for the note)
  - products: `results/phase2_*.png`, extended `booster_run.mat` fields `.solD .repD .mcD`

- [ ] **Step 1: Write the failing test**

`tests/test_phase2_drag.m`:

```matlab
% TEST_PHASE2_DRAG  Drag-on collocation solve (coarse, warm-started from
% vacuum) converges; fuel differs from vacuum by a NONZERO but sane amount
% (drag helps braking: expect LESS fuel, sanity band 0 < dfuel < 40% of
% vacuum fuel); certify gates G1/G2/G5 pass with G3/G4 skipped.
%
% INPUTS: none   OUTPUTS: none (throws on failure)
here_ = fileparts(mfilename('fullpath'));
addpath(fullfile(here_, '..'));  setup_paths;

P    = booster_params();
solC = solve_pdg_colloc(P, struct('N', 30));
Pd   = P;  Pd.drag.on = true;
solD = solve_pdg_colloc(Pd, struct('N', 30, 'init', solC));
assert(solD.stats.success, 'drag solve failed');
fuelV = P.m0 - solC.mf;  fuelD = P.m0 - solD.mf;
assert(fuelD < fuelV, 'drag should reduce fuel (braking), got +%.1f kg', ...
       fuelD - fuelV);
assert(fuelV - fuelD < 0.4*fuelV, 'drag effect implausibly large');
rep = certify_pdg(solD, [], Pd);
assert(rep.G1_pass && rep.G2_pass && rep.G5_pass, 'drag gates failed');
assert(isequal(rep.G3_pass, 'skipped'), 'G3 should be skipped without twin');
fprintf('test_phase2_drag PASS  fuel vac=%.1f drag=%.1f kg\n', fuelV, fuelD);
```

- [ ] **Step 2: Run test to verify it fails**

Run: `/Applications/MATLAB_R2025b.app/bin/matlab -batch "cd('/Users/msc/Desktop/optimal_control/booster_landing'); setup_paths; test_phase2_drag"`
Expected: FAIL only if machinery breaks — Tasks 2/3/5 already built the branches. If it PASSES immediately, good: that's the opt-in-flag pattern paying off; proceed to Step 3.
(One likely genuine failure: G5 with drag — the PMP thrust-direction law changes when drag is present; if primer alignment fails ONLY in the drag case, relax G5 for drag runs to the bang-bang structure check alone, document why in `certify_pdg`'s header, and carry the primer question to the note.)

- [ ] **Step 3: Implement the front-door `phase2` branch + `plot_vacuum_vs_drag`**

In `run_booster_landing`, after stage 5, when `cfg.phase2`:

```matlab
if isfield(cfg,'phase2') && cfg.phase2
    fprintf('=== [P2] Drag-on re-solve (warm-started) ===\n');
    Pd = P;  Pd.drag.on = true;
    R.solD = solve_pdg_colloc(Pd, struct('init', R.solC));
    R.repD = certify_pdg(R.solD, [], Pd);
    print_certify_report(R.repD);
    ctrlD  = tvlqr_design(R.solD, Pd);
    R.mcD  = run_monte_carlo(R.solD, ctrlD, Pd, struct('Nrun', cfg.Nrun));
    plot_vacuum_vs_drag(R.solC, R.solD, fullfile(cfg.outdir, 'phase2_vac_vs_drag.png'));
    plot_footprint(R.mcD, Pd, fullfile(cfg.outdir, 'phase2_footprint.png'));
    fprintf('P2: fuel vac %.1f kg -> drag %.1f kg; MC(wind) %.1f%%\n', ...
            P.m0 - R.solC.mf, Pd.m0 - R.solD.mf, 100*R.mcD.success_rate);
end
```

`plot_vacuum_vs_drag`: 1×3 tiledlayout (trajectory downrange–altitude overlay; throttle overlay; mass overlay), Δfuel/Δtf annotation, house colors, exportgraphics.

- [ ] **Step 4: Run tests + Phase-2 flagship**

Run: the Step-2 test command (expect PASS), then
`/Applications/MATLAB_R2025b.app/bin/matlab -batch "cd('/Users/msc/Desktop/optimal_control/booster_landing'); run_booster_landing(struct('phase2', true))"`
Expected: Phase 1 unchanged (gates still ALL PASS — this is the re-run-the-flagship-after-structural-change rule), Phase 2 converged, wind-MC success reported (≥95% not required by spec for Phase 2 — report the honest number), comparison PNG written.

- [ ] **Step 5: Commit**

```bash
cd /Users/msc/Desktop/optimal_control
git add booster_landing/run_booster_landing.m booster_landing/tests/test_phase2_drag.m \
        booster_landing/viz/plot_vacuum_vs_drag.m
git commit -m "booster_landing: Phase 2 atmosphere -- drag re-solve, wind MC, vacuum comparison

Co-Authored-By: Claude Fable 5 <noreply@anthropic.com>"
```

---

