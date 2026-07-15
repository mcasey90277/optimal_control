# Explicit Insertion Points Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Make `rv0` (GTO departure) and `rvf` (tulip/ELFO insertion point) explicit and changeable in the high-level drivers via one source-of-truth helper, with a fail-loud drift guard, endpoint metadata in every output, and filename tags — keeping the current points (Option A, zero re-solve).

**Architecture:** A new `insertion_states(target, criterion)` helper returns `[rv0, rvf, meta]`. Every high-level driver declares its endpoints through it and (for seed consumers) asserts the loaded seed's `rvf` matches. Outputs gain `rv0`/`rvf`/`insertion` fields and a `label` filename tag.

**Tech Stack:** MATLAB R2025b, pumpkyn toolbox (`getTulip`, `orb2eci`, `fromPCI`), CasADi only for the full end-to-end solve validations (the helper + guards are CasADi-free).

## Global Constraints

- MATLAB R2025b only: `/Applications/MATLAB_R2025b.app/bin/matlab -batch "…"`.
- **Option A — keep the current points; zero re-solve.** The helper's defaults MUST equal what existing seeds hold (verified: all tulip backbones + the ELFO seed match to 0.00), so the drift guards pass unchanged.
- **The exact hardcoded constants** (do not alter):
  - `rv0 = [0.00349629072294633, -0.0072962582600817, 0, 4.19147893803368, 8.98865558978329, 0]`
  - tulip `campaign` `rvf = [1.00658107295709, 0.0425745746906059, -0.0557780910480905, -0.16004281347248, 0.0665702939657711, -0.260455693516549]`
- Every MATLAB function needs the commented header block (purpose, INPUTS w/ sizes, OUTPUTS w/ sizes, REFERENCES). No `i`/`j` loop indices. Result `.mat`/`results/logs/` gitignored — not staged.
- Drift-guard tolerance: `norm(seed.rvf - rvf) < 1e-10`, `norm(seed.rv0 - rv0) < 1e-10`.
- Guard message on failure must name the mismatch, e.g. `error('insertion:drift','seed rvf differs from declared %s insertion by %.2e — regenerate the seed for this criterion', insMeta.label, d)`.

---

### Task 1: `insertion_states.m` helper + unit test

**Files:**
- Create: `NLP_lowThrust_GTO_tulip/sundman_minfuel/insertion_states.m`
- Test: `NLP_lowThrust_GTO_tulip/sundman_minfuel/test_insertion_states.m`

**Interfaces:**
- Produces: `[rv0, rvf, meta] = insertion_states(target, criterion)` — `rv0` [1x6], `rvf` [1x6], `meta` struct `.target .criterion .label`. `target` `'tulip'|'elfo'`; tulip criteria `'campaign'`(default)`|'maxydot'|'apoapsis'`; elfo criteria `'nearest'`(default)`|'apolune'|'perilune'`.
- Consumes: `cr3bp_lt_params`, `gto_tulip_endpoints` (same dir), `gto_elfo_endpoints` (`../elfo`), pumpkyn.

- [ ] **Step 1: Write the failing test**

Create `sundman_minfuel/test_insertion_states.m`:
```matlab
% TEST_INSERTION_STATES  Verify the insertion-point helper: the default criteria
% reproduce exactly what the existing seeds hold (so the drift guards pass with
% zero re-solve), and the alternate criteria return valid 6-states.
here = fileparts(mfilename('fullpath'));  cd(here);  setup_paths();  addpath('../elfo');
E  = load('results/energy/energy_f1120.mat');           % a tulip backbone
Ee = load('../elfo/results/energy_elfo_freetf.mat');    % the ELFO seed

[rv0,rvfC,mC] = insertion_states('tulip','campaign');
assert(norm(rv0  - E.rv0(:).') < 1e-12, 'rv0 != backbone rv0');
assert(norm(rvfC - E.rvf(:).') < 1e-12, 'tulip campaign rvf != backbone rvf');
assert(strcmp(mC.label,'tulipCampaign'), 'wrong tulip label');

[~,rvfN,mN] = insertion_states('elfo','nearest');
assert(norm(rvfN - Ee.rvf(:).') < 1e-12, 'elfo nearest rvf != ELFO seed rvf');
assert(strcmp(mN.label,'elfoNearest'), 'wrong elfo label');

[~,rvfM] = insertion_states('tulip','maxydot');
[~,rvfA] = insertion_states('tulip','apoapsis');
[~,rvfP] = insertion_states('elfo','apolune');
assert(all(isfinite(rvfM)) && numel(rvfM)==6, 'maxydot invalid');
assert(all(isfinite(rvfA)) && numel(rvfA)==6, 'apoapsis invalid');
assert(all(isfinite(rvfP)) && numel(rvfP)==6, 'apolune invalid');

% default criterion (omitted arg) == campaign / nearest
[~,rvfCd] = insertion_states('tulip');  assert(norm(rvfCd-rvfC)<1e-15,'tulip default');
[~,rvfNd] = insertion_states('elfo');   assert(norm(rvfNd-rvfN)<1e-15,'elfo default');
fprintf('TEST_INSERTION_STATES: PASS (defaults match seeds <1e-12; alternates valid)\n');
```

- [ ] **Step 2: Run the test to verify it fails**

Run:
```bash
/Applications/MATLAB_R2025b.app/bin/matlab -batch "cd('/Users/msc/Desktop/optimal_control/NLP_lowThrust_GTO_tulip/sundman_minfuel'); setup_paths; test_insertion_states"
```
Expected: FAIL — `Unrecognized function or variable 'insertion_states'`.

- [ ] **Step 3: Write the helper**

Create `sundman_minfuel/insertion_states.m`:
```matlab
function [rv0, rvf, meta] = insertion_states(target, criterion)
% INSERTION_STATES  Single source of truth for the GTO departure (rv0) and the
% tulip/ELFO insertion (rendezvous) state (rvf) used by every low-thrust
% pipeline. Declaring endpoints here (instead of threading them implicitly from
% a seed .mat) makes them explicit, changeable, and drift-checkable.
%
% INPUTS:
%   target    - 'tulip' | 'elfo' [char]
%   criterion - (optional) tulip: 'campaign'(default)|'maxydot'|'apoapsis'
%                          elfo:  'nearest'(default)|'apolune'|'perilune' [char]
%
% OUTPUTS:
%   rv0  - GTO departure state, ND rotating frame [1x6]
%   rvf  - insertion (rendezvous) state, ND rotating frame [1x6]
%   meta - struct: .target .criterion .label (label for filenames/provenance)
%
% NOTE: 'campaign'/'nearest' reproduce EXACTLY what the current seeds hold
% (Option A). Using an alternate criterion requires a matching energy seed --
% the consumer drivers' drift guard will fail loudly until one exists.
%
% REFERENCES:
%   [1] gto_tulip_endpoints.m (max-ydot tulip point + trace);
%   [2] gto_elfo_endpoints.m  (ELFO apolune/perilune/nearest).

if nargin < 2 || isempty(criterion)
    switch lower(target)
        case 'tulip', criterion = 'campaign';
        case 'elfo',  criterion = 'nearest';
        otherwise, error('insertion_states:target','unknown target %s', target);
    end
end
p = cr3bp_lt_params(25e-3, 15, 2100);
muStar = p.muStar;  tStar = p.tStar;  lStar = p.lStar; %#ok<NASGU>

% --- GTO departure (shared by all pipelines) --------------------------------
rv0 = [0.00349629072294633, -0.0072962582600817, 0, ...
       4.19147893803368, 8.98865558978329, 0];          % GTO departure (exact)
% -- to regenerate rv0 from the GTO orbital elements, uncomment: -------------
% muEarth = 6.67384e-20*(1-muStar)*(5.9736e24 + 7.35e22);
% sma = (6378+350 + 6378+35786)/2;  ecc = (35786-350)/(2*sma);
% [r0,v0] = pumpkyn.cr3bp.orb2eci(muEarth, [sma,ecc,0,-25*pi/180,0,0], 2);
% rv0 = pumpkyn.cr3bp.fromPCI(0, [r0,v0], muStar, tStar, lStar, 1);

here = fileparts(mfilename('fullpath'));

% --- insertion (rendezvous) state -------------------------------------------
switch lower(target)
  case 'tulip'
    switch lower(criterion)
      case 'campaign'
        rvf = [1.00658107295709, 0.0425745746906059, -0.0557780910480905, ...
               -0.16004281347248, 0.0665702939657711, -0.260455693516549];
        label = 'tulipCampaign';
      case 'maxydot'
        addpath(here);  [~, rvf] = gto_tulip_endpoints(p);      % max-ydot
        label = 'tulipMaxYdot';
      case 'apoapsis'
        addpath(here);  [~, ~, tr] = gto_tulip_endpoints(p);
        [~, idx] = min(vecnorm(tr(:,4:6), 2, 2));               % slowest point
        rvf = tr(idx, 1:6);
        label = 'tulipApoapsis';
      otherwise, error('insertion_states:crit','unknown tulip criterion %s', criterion);
    end
  case 'elfo'
    addpath(fullfile(here, '..', 'elfo'));                      % gto_elfo_endpoints
    switch lower(criterion)
      case 'nearest'
        [~, rvfTul] = insertion_states('tulip', 'campaign');    % ref = tulip campaign
        [~, rvf] = gto_elfo_endpoints(p, struct('point','nearest','ref',rvfTul));
        label = 'elfoNearest';
      case 'apolune'
        [~, rvf] = gto_elfo_endpoints(p, struct('point','apolune'));
        label = 'elfoApolune';
      case 'perilune'
        [~, rvf] = gto_elfo_endpoints(p, struct('point','perilune'));
        label = 'elfoPerilune';
      otherwise, error('insertion_states:crit','unknown elfo criterion %s', criterion);
    end
  otherwise, error('insertion_states:target','unknown target %s', target);
end
rvf  = rvf(:).';
meta = struct('target',lower(target),'criterion',lower(criterion),'label',label);
end
```

- [ ] **Step 4: Run the test to verify it passes**

Run the Step-2 command. Expected: `TEST_INSERTION_STATES: PASS (…)`.

- [ ] **Step 5: Commit**

```bash
cd /Users/msc/Desktop/optimal_control
git add NLP_lowThrust_GTO_tulip/sundman_minfuel/insertion_states.m NLP_lowThrust_GTO_tulip/sundman_minfuel/test_insertion_states.m
git commit -m "feat(insertion-points): insertion_states helper (single source of truth) + test"
```

---

### Task 2: Wire the explicit block + drift guard into the seed-consumer drivers

**Files (each modified):** `PSR/run_psr.m`, `PSR/psr_run_one.m`, `elfo/run_elfo_minfuel.m`, `elfo/elfo_run_one.m`, `sundman_minfuel/gen_tulip_mintime.m`, `elfo/gen_elfo_mintime.m`

**Interfaces:**
- Consumes: `insertion_states` (Task 1). These drivers currently load a seed and read `rvf`/`rv0` from it (audited lines: `gen_tulip_mintime.m:52,79`, `gen_elfo_mintime.m:38,53`, `minfuel_at_tf.m:53,58` reached via `run_psr`/`psr_run_one`, and the elfo run drivers).

**The uniform pattern to insert** — a declaration block near each driver's parameters, and a guard right after the seed is loaded. `TGT` is `'tulip'` for the tulip pipelines (`run_psr`, `psr_run_one`, `gen_tulip_mintime`) and `'elfo'` for the ELFO ones (`run_elfo_minfuel`, `elfo_run_one`, `gen_elfo_mintime`):

```matlab
% ---- INSERTION POINT (edit here to retarget) ---------------------------------
insertion = 'campaign';        % tulip: 'campaign'|'maxydot'|'apoapsis'  (elfo: 'nearest'|'apolune'|'perilune')
% insertion = 'maxydot';       % uncomment to use the max-ydot point (needs a matching energy seed)
% insertion = 'apoapsis';      % uncomment to use the slowest/apoapsis point (needs a matching seed)
[rv0, rvf, insMeta] = insertion_states('<TGT>', insertion);   % <TGT> = 'tulip' or 'elfo'
```
And the guard, immediately after the seed struct (`S`/`E`) is loaded:
```matlab
% drift guard: the seed must be for the declared insertion point
assert(norm(S.rvf(:).' - rvf) < 1e-10 && norm(S.rv0(:).' - rv0) < 1e-10, ...
    'insertion:drift', ['seed endpoints differ from the declared %s insertion ' ...
    '(rvf %.2e, rv0 %.2e) -- regenerate the seed for this criterion'], ...
    insMeta.label, norm(S.rvf(:).'-rvf), norm(S.rv0(:).'-rv0));
```
Then the driver uses the declared `rv0`/`rvf` (from the helper) in its solve call rather than `S.rvf`/`S.rv0`. For `run_psr`/`psr_run_one`, the seed is threaded through `minfuel_at_tf`; add the guard in the driver right after it loads/locates the seed (before calling `minfuel_at_tf`), and pass the declared `rvf`/`rv0` if the call signature allows, else keep the guard as the safety check.

- [ ] **Step 1: Read each of the 6 files and locate the seed load + the solve call**

For each file, note (a) where the seed struct is loaded, (b) where `rvf`/`rv0` are first used, (c) where results are saved (for Task 4). Record findings before editing.

- [ ] **Step 2: Insert the declaration block + guard into each file**

Apply the pattern above to all 6 files (`<TGT>` per the file's pipeline). Keep each driver's existing seed-loading; add the block near the parameters and the guard after the load. Replace direct `S.rvf`/`S.rv0` (or `E.rvf`/`E.rv0`) uses in the *solve call* with the declared `rvf`/`rv0` where straightforward; otherwise leave the solve as-is (the guard already proves they are equal).

- [ ] **Step 3: Verify the guards pass on the current seeds (CasADi-free)**

Run a check that, for the seed files these drivers use, `insertion_states` matches (i.e. the guard would pass):
```bash
/Applications/MATLAB_R2025b.app/bin/matlab -batch "cd('/Users/msc/Desktop/optimal_control/NLP_lowThrust_GTO_tulip/sundman_minfuel'); setup_paths; addpath('../elfo'); \
 [~,rt]=insertion_states('tulip','campaign'); [~,re]=insertion_states('elfo','nearest'); \
 Eb=load('results/energy/energy_f1120.mat');          fprintf('tulip backbone drvf=%.1e\n', norm(Eb.rvf(:).'-rt)); \
 T =load('results/energy/energy_tulip_2p.mat');        fprintf('tulip 2p seed  drvf=%.1e\n', norm(T.rvf(:).'-rt)); \
 Ee=load('../elfo/results/energy_elfo_freetf.mat');    fprintf('elfo seed      drvf=%.1e\n', norm(Ee.rvf(:).'-re));"
```
Expected: all three `drvf` < 1e-10 (guards pass on the current seeds — zero re-solve).

- [ ] **Step 4: Commit**

```bash
cd /Users/msc/Desktop/optimal_control
git add NLP_lowThrust_GTO_tulip/PSR/run_psr.m NLP_lowThrust_GTO_tulip/PSR/psr_run_one.m NLP_lowThrust_GTO_tulip/elfo/run_elfo_minfuel.m NLP_lowThrust_GTO_tulip/elfo/elfo_run_one.m NLP_lowThrust_GTO_tulip/sundman_minfuel/gen_tulip_mintime.m NLP_lowThrust_GTO_tulip/elfo/gen_elfo_mintime.m
git commit -m "feat(insertion-points): explicit block + drift guard in the 6 seed-consumer drivers"
```

---

### Task 3: Wire the seed generators

**Files:** `PSR/gen_energy_seed.m`, `sundman_minfuel/gen_tulip_energy_2p.m`, `elfo/gen_elfo_energy_gravhom.m`, `elfo/gen_elfo_energy_tfsweep.m`

**Interfaces:** Consumes `insertion_states`. Two roles:
- **Continue-from-backbone** (`gen_energy_seed`, `gen_tulip_energy_2p`, `gen_elfo_energy_tfsweep`): add the declaration block (`'tulip'`/`'elfo'` per pipeline) + the same drift guard against the base backbone they load.
- **Compute-the-target** (`gen_elfo_energy_gravhom`): replace its `gto_elfo_endpoints('nearest', ref=rvf_tul)` call (line ~75) and the `rvf_tul = E.rvf` (line ~71) with `insertion_states`:
  ```matlab
  [~, rvf_tul] = insertion_states('tulip','campaign');     % the ELFO 'nearest' reference
  [~, rvf_elfo, insMeta] = insertion_states('elfo','nearest');
  ```
  and add the drift guard against the base tulip backbone `E.rvf` it loads.

- [ ] **Step 1: Read the 4 files; note seed load + target computation + save lines**
- [ ] **Step 2: Apply the block + guard (continue-from-backbone) / helper substitution (gravhom)**
- [ ] **Step 3: CasADi-free check** that the substitutions produce the same `rvf` as before (for `gen_elfo_energy_gravhom`, `insertion_states('elfo','nearest')` equals the old `gto_elfo_endpoints('nearest',ref=E.rvf)` — verified 0.00 in the spec pre-facts).
- [ ] **Step 4: Commit**
```bash
git add NLP_lowThrust_GTO_tulip/PSR/gen_energy_seed.m NLP_lowThrust_GTO_tulip/sundman_minfuel/gen_tulip_energy_2p.m NLP_lowThrust_GTO_tulip/elfo/gen_elfo_energy_gravhom.m NLP_lowThrust_GTO_tulip/elfo/gen_elfo_energy_tfsweep.m
git commit -m "feat(insertion-points): declare/guard endpoints in the 4 seed generators"
```

---

### Task 4: Endpoint metadata + filename tags in saved outputs

**Files:** the save points reached by the above drivers — `sundman_minfuel/gen_tulip_mintime.m` + `elfo/gen_elfo_mintime.m` (min-time `.mat` save), `elfo/elfo_export_data.m` + `PSR/psr_export_data.m` (data products), the seed-save lines in the Task-3 generators, and the batch summary savers (`elfo/elfo_collect_summary.m`, `PSR/psr_collect_summary.m` if present).

**Pattern:**
- Add `rv0`, `rvf`, and `insertion` (= `insMeta.label`) to every saved result struct / `save(...)` variable list.
- Add the label to new output filenames: insert `insMeta.label` into the `sprintf`/`fullfile` that builds each output name (e.g. `mintime_tulip_%s.mat` → `mintime_tulip_tulipCampaign.mat`; `psr_data_..._%s_...`). Existing files are NOT renamed (Option A) — the tag applies to newly-written files.

- [ ] **Step 1: Read each save point; note the `save(...)` var lists and filename builders**
- [ ] **Step 2: Add `rv0`/`rvf`/`insertion` to each saved struct; add `insMeta.label` to each new filename**
- [ ] **Step 3: Verify (CasADi-free where possible)** — re-save a min-time solution via `gen_tulip_mintime` metadata path is CasADi-bound; instead unit-check that a representative saver, given `insMeta`, writes the three fields + tagged name (a small synthetic-struct save/load round-trip, mirroring how `elfo_collect_summary` was unit-tested CasADi-free).
- [ ] **Step 4: Commit**
```bash
git add -A NLP_lowThrust_GTO_tulip/
git commit -m "feat(insertion-points): save rv0/rvf/insertion + label tag in outputs"
```

---

## Notes for the implementer

- Tasks 2–4 touch drivers whose full end-to-end validation needs CasADi/R2025b. The **helper (Task 1) and the drift-guard checks are CasADi-free** and are the primary tests; a full solve run is a bonus check if the environment allows — do NOT fabricate solve results.
- Zero re-solve is the contract: if any drift guard fails on a *current* seed, STOP — that means the helper default doesn't match the seed, which contradicts the verified pre-facts; report it rather than "fixing" by changing the seed.
- Do not `git add` anything under `results/` (gitignored).
- After all tasks: updating the ROADMAP/paper to reference the explicit insertion API is a follow-up, not part of this plan.
