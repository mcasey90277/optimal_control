### Task 10: Front door (`run_booster_landing`) + flagship run

**Files:**
- Create: `run_booster_landing.m`
- Create: `tests/test_run_front_door.m`
- Create: `README.md`

**Interfaces:**
- Consumes: everything above.
- Produces: `R = run_booster_landing(cfg)` — the campaign. No args = full nominal Phase-1 campaign. `cfg` fields (all optional): `.P` (params override struct — merged onto `booster_params()`), `.doMovie` (default true), `.doMC` (default true), `.Nrun` (default 200), `.outdir` (default `results/`). Sequence, with `fprintf` stage banners:
  1. `setup_paths` (self — front door must work from a fresh MATLAB after `cd booster_landing`)
  2. solve colloc (N=60) → solve convex (golden tf) → `certify_pdg` → **print gate table**
  3. `tvlqr_design` → nominal `sim_closed_loop` → `run_monte_carlo`
  4. viz: `pdg_solution.png`, `footprint.png`, `landing.mp4` into outdir
  5. save `R` (all products: solC, solV, rep, ctrl, out0, mc, timestamps via `datetime`) to `outdir/booster_run.mat`
  6. final summary block: tf, mf, fuel used, gates verdict, MC success rate — the numbers a reader quotes.
- README.md: what the campaign is, how to run it (`matlab -batch "cd ...; run_booster_landing"`), what each folder holds, the spec/plan pointers.

- [ ] **Step 1: Write the failing test**

`tests/test_run_front_door.m`:

```matlab
% TEST_RUN_FRONT_DOOR  Front-door contract on a FAST config: runs end to
% end from only setup_paths, honors cfg overrides (the advertised-but-
% ignored-options bug is the canonical failure -- verify N actually
% reached the solver via solution size), writes its products.
%
% INPUTS: none   OUTPUTS: none (throws on failure)
here_ = fileparts(mfilename('fullpath'));
addpath(fullfile(here_, '..'));  setup_paths;

od  = fullfile(tempdir, 'bl_front_smoke');
cfg = struct('doMovie', false, 'Nrun', 6, 'outdir', od, ...
             'P', struct('N', 20, 'Nconv', 60));
R   = run_booster_landing(cfg);
assert(size(R.solC.X, 2) == 21, 'cfg.P.N did not reach the solver');
assert(R.rep.all_pass || ischar(R.rep.G3_pass), 'gates failed on fast config');
assert(isfile(fullfile(od, 'booster_run.mat')), 'products not written');
assert(isfile(fullfile(od, 'pdg_solution.png')), 'solution plot missing');
fprintf('test_run_front_door PASS\n');
```

- [ ] **Step 2: Run test to verify it fails**

Run: `/Applications/MATLAB_R2025b.app/bin/matlab -batch "cd('/Users/msc/Desktop/optimal_control/booster_landing'); setup_paths; test_run_front_door"`
Expected: FAIL — `run_booster_landing` undefined.

- [ ] **Step 3: Implement `run_booster_landing.m` + README**

Front-door skeleton (ADJUSTABLE PARAMETERS block up top, per the standing principles; params override by field-merge):

```matlab
function R = run_booster_landing(cfg)
% RUN_BOOSTER_LANDING  Front door: full 3-DOF booster-landing campaign.
%
%   solve (colloc + convex) -> certify G1-G5 -> TVLQR -> closed loop ->
%   Monte Carlo -> plots + movie -> results/booster_run.mat
%
% Run me with no arguments for the nominal campaign:
%   /Applications/MATLAB_R2025b.app/bin/matlab -batch ...
%     "cd('~/Desktop/optimal_control/booster_landing'); run_booster_landing"
%
% INPUTS:
%   cfg - (optional) .P params overrides, .doMovie [true], .doMC [true],
%         .Nrun [200], .outdir ['results/']
% OUTPUTS:
%   R   - everything: .P .solC .solV .rep .ctrl .out0 .mc .when
%
% REFERENCES: spec at docs/superpowers/specs/2026-08-08-booster-landing-design.md
setup_paths;

%% ---------------- ADJUSTABLE PARAMETERS (defaults) ----------------
def = struct('doMovie', true, 'doMC', true, 'Nrun', 200, ...
             'outdir', fullfile(fileparts(mfilename('fullpath')), 'results'));
%% ------------------------------------------------------------------
if nargin < 1, cfg = struct(); end
fn = fieldnames(def);
for k = 1:numel(fn)
    if ~isfield(cfg, fn{k}), cfg.(fn{k}) = def.(fn{k}); end
end
P = booster_params();
if isfield(cfg, 'P')
    pf = fieldnames(cfg.P);
    for k = 1:numel(pf), P.(pf{k}) = cfg.P.(pf{k}); end
end
if ~exist(cfg.outdir, 'dir'), mkdir(cfg.outdir); end

fprintf('=== [1/5] Guidance: collocation (N=%d) ===\n', P.N);
R.solC = solve_pdg_colloc(P);
fprintf('    tf=%.3f s  mf=%.2f kg\n', R.solC.tf, R.solC.mf);

fprintf('=== [2/5] Guidance: lossless convexification (golden tf) ===\n');
R.solV = solve_pdg_convex(P);
fprintf('    tf=%.3f s  mf=%.2f kg  gap=%.2e\n', ...
        R.solV.tf, R.solV.mf, R.solV.lossless_gap);

fprintf('=== [3/5] Certification ===\n');
R.rep = certify_pdg(R.solC, R.solV, P);
print_certify_report(R.rep);

fprintf('=== [4/5] Tracking: TVLQR + closed loop%s ===\n', ...
        ternary(cfg.doMC, ' + Monte Carlo', ''));
R.ctrl = tvlqr_design(R.solC, P);
R.out0 = sim_closed_loop(R.solC, R.ctrl, P, struct());
if cfg.doMC
    R.mc = run_monte_carlo(R.solC, R.ctrl, P, struct('Nrun', cfg.Nrun));
end

fprintf('=== [5/5] Products -> %s ===\n', cfg.outdir);
plot_pdg_solution(R.solC, R.solV, fullfile(cfg.outdir, 'pdg_solution.png'));
if cfg.doMC
    plot_footprint(R.mc, P, fullfile(cfg.outdir, 'footprint.png'));
end
if cfg.doMovie
    movie_landing(R.out0, R.solC, P, fullfile(cfg.outdir, 'landing.mp4'));
end
R.P = P;  R.when = datetime('now');
save(fullfile(cfg.outdir, 'booster_run.mat'), '-struct', 'R');

fprintf('\n==================== SUMMARY ====================\n');
fprintf('tf        %.3f s      fuel  %.1f kg (mf %.1f)\n', ...
        R.solC.tf, P.m0 - R.solC.mf, R.solC.mf);
fprintf('gates     %s\n', ternary(R.rep.all_pass, 'ALL PASS', 'FAILURES -- see table'));
fprintf('nom miss  %.2f m  vtd %.2f m/s\n', R.out0.td.miss, R.out0.td.vtd);
if cfg.doMC
    fprintf('MC        %d runs, success %.1f%%\n', cfg.Nrun, 100*R.mc.success_rate);
end
fprintf('=================================================\n');
end

function s = ternary(c, a, b), if c, s = a; else, s = b; end, end
```

- [ ] **Step 4: Run test to verify it passes**

Same command. Expected: `test_run_front_door PASS`.

- [ ] **Step 5: FLAGSHIP RUN — the real acceptance test**

Run (expect a few minutes): `/Applications/MATLAB_R2025b.app/bin/matlab -batch "cd('/Users/msc/Desktop/optimal_control/booster_landing'); run_booster_landing"`
Expected, per spec success criteria — verify each, honestly, before claiming done:
- both solvers converged; gate table ALL PASS; final masses agree < 1.0 kg (measured ~0.7 kg Taylor model offset — adjudicated, documented)
- throttle plot shows bang-bang min–max with terminal max-throttle arc (adjudicated: min–max IS the optimum here; open `results/pdg_solution.png` and LOOK)
- MC success ≥ 95% at 200 runs
- `results/landing.mp4` plays clean (open it — no diagonal streaks)
Record tf / mf / fuel / success-rate in the commit message. If any criterion fails, stop and fix (systematic-debugging) — do not commit a red flagship.

- [ ] **Step 6: Commit**

```bash
cd /Users/msc/Desktop/optimal_control
git add booster_landing/run_booster_landing.m booster_landing/tests/test_run_front_door.m \
        booster_landing/README.md
git commit -m "booster_landing: front door + flagship run

Flagship: tf=<fill> s, fuel=<fill> kg, gates ALL PASS, MC <fill>%/200.

Co-Authored-By: Claude Fable 5 <noreply@anthropic.com>"
```

---

