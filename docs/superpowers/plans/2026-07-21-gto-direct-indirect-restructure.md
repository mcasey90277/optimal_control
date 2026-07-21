# GTO Direct/Indirect Restructure — Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Restructure the CR3BP GTO campaigns into `GTO_tulip/{direct,indirect}/` and a new `GTO_ELFO/{direct,indirect}/`, with the cross-module-shared problem definition extracted to a new `orbit_transfer/cr3bp_common/` library.

**Architecture:** Two-phase commit discipline: Task 1 is a **pure structural** `git mv` commit (history-preserving, zero content edits, tree intentionally non-runnable at that one commit); Task 2 is the **rewiring** commit (new `setup_cr3bp_common.m`, 8 rewritten `setup_paths.m`, ~10 enumerated executable-ref fixes) gated on a full 8-module MATLAB resolve-smoke matrix. Tasks 3–4 repoint doc references and update the CLAUDE.md maps, then push. Spec: `docs/superpowers/specs/2026-07-21-gto-direct-indirect-restructure-design.md` (authoritative — read it).

**Tech Stack:** git, zsh, MATLAB R2025b (`/Applications/MATLAB_R2025b.app/bin/matlab`).

## Global Constraints

- **Repo root:** `/Users/msc/Desktop/optimal_control` (remote `origin git@github.com:mcasey90277/optimal_control.git`, branch `main`). Run all commands from here unless a step says otherwise.
- **Scope:** ONLY `orbit_transfer/GTO_tulip/`, `orbit_transfer/lowThrust_GTO_tulip/`, the new `orbit_transfer/{cr3bp_common,GTO_ELFO}/`, doc-reference fixes in `orbit_transfer/earth_elliptic_to_geo/`, and `CLAUDE.md`. Touch nothing else (`earth_elliptic_to_geo` code, `lambert`, `min_energy_tutorial`, `mfmax`, `min_fuel_paper(s)`, `papers`, top-level tutorials are all off-limits except the named doc repoints).
- **No logic/numeric changes.** Pure relocation + path rewiring + reference repointing.
- **PSR is self-contained** (its `lib/` is a ~20-file vendored library): its `setup_paths.m` and its 4 inline `addpath(fullfile(here,'..','sundman_minfuel'))` lines are **NOT edited** — the `../sundman_minfuel` sibling relationship is preserved under `direct/`. Do not "helpfully" de-dup `PSR/lib/` — explicitly out of scope (spec §3).
- **`PSR_data/` and `IFS_data/` are fully untracked** (all-`.mat`, gitignored) → plain `mv`, not `git mv`. Untracked files inside tracked dirs (e.g. `elfo/results/`) travel automatically with `git mv` of the dir (filesystem rename).
- **MATLAB invocation:** `/Applications/MATLAB_R2025b.app/bin/matlab -batch "<ONE line>"` (multi-line `-batch` gets mangled by zsh). MATLAB R2025b ONLY. Filter license-banner noise with `| grep -vE "License|academic|personal use|government"`.
- **`setup_paths` resolution relies on cwd-shadowing:** every module keeps its own `setup_paths.m`, and callers `cd` into the module first — the cwd's copy shadows any on-path copy. Preserve this: never delete a module's `setup_paths.m`.
- **Staging discipline:** never `git add -A`. `git mv` auto-stages renames; `git add` only the specific files each task names. `git status --short` before every commit; `papers/.DS_Store` and untracked clutter must never be committed.
- **Historical docs carve-out:** `docs/superpowers/plans/*`, `docs/superpowers/specs/*`, `.superpowers/` keep old paths as point-in-time records — excluded from all sweeps/fixes.
- **Commit trailer (exact):** `Co-Authored-By: Claude Fable 5 <noreply@anthropic.com>`
- If any `git mv` errors, a smoke does not print its expected line, or unexpected files appear staged: STOP, report BLOCKED with exact output. Do not improvise.

---

### Task 1: Pure structural moves (git mv everything; no content edits)

**Files:**
- Create dirs: `orbit_transfer/{cr3bp_common, GTO_ELFO/direct, GTO_ELFO/indirect, GTO_tulip/direct, GTO_tulip/indirect}`
- Move (git mv): 8 tracked modules + `lowThrust_GTO_tulip` + 4 shared param files; plain-`mv`: `PSR_data`, `IFS_data`
- Create: `orbit_transfer/GTO_ELFO/indirect/README.md` (placeholder stub)

**Interfaces:**
- Produces: the spec §2 tree. Modules land at uniform depth `orbit_transfer/<PROBLEM>/<METHOD>/<module>/`. The 4 shared files land in `cr3bp_common/`: `cr3bp_lt_params.m`, `minfuel_config.m`, `gto_tulip_endpoints.m` (from `sundman_minfuel/`), `gto_elfo_endpoints.m` (from `elfo/`).
- **Known-broken window:** after this commit the moved modules' `setup_paths` still point at old relative locations. This is intentional and repaired in Task 2 (the very next commit). The commit message must say so.

- [ ] **Step 1: Create the skeleton and move the modules**

```bash
cd /Users/msc/Desktop/optimal_control/orbit_transfer
mkdir -p GTO_tulip/direct GTO_tulip/indirect GTO_ELFO/direct GTO_ELFO/indirect cr3bp_common
git mv GTO_tulip/sundman_minfuel GTO_tulip/PSR GTO_tulip/movie GTO_tulip/direct/
mv GTO_tulip/PSR_data GTO_tulip/direct/            # untracked (all-.mat) -> plain mv
git mv GTO_tulip/ms_band GTO_tulip/ifs GTO_tulip/ztl GTO_tulip/min_time GTO_tulip/indirect/
mv GTO_tulip/IFS_data GTO_tulip/indirect/          # untracked -> plain mv
git mv lowThrust_GTO_tulip GTO_tulip/indirect/
git mv GTO_tulip/elfo GTO_ELFO/direct/
```

- [ ] **Step 2: Move the 4 shared files into cr3bp_common**

```bash
git mv GTO_tulip/direct/sundman_minfuel/cr3bp_lt_params.m \
       GTO_tulip/direct/sundman_minfuel/minfuel_config.m \
       GTO_tulip/direct/sundman_minfuel/gto_tulip_endpoints.m cr3bp_common/
git mv GTO_ELFO/direct/elfo/gto_elfo_endpoints.m cr3bp_common/
```

- [ ] **Step 3: Write the GTO_ELFO/indirect placeholder README**

Create `orbit_transfer/GTO_ELFO/indirect/README.md` with exactly:

```markdown
# GTO_ELFO/indirect — placeholder

Indirect (PMP shooting) GTO→ELFO transfer work. Not yet started: the ELFO
min-time indirect solve ("Route C") is open future work — the direct min-time
anchor (Route B, `casadi_mintime_freetf` + `gen_elfo_mintime`, tfMin_ELFO =
6.0962 ND = 27.02 d, anchored 2026-07-15) lives in `../direct/elfo/`.
The shared CR3BP problem definition is in `../../cr3bp_common/`.
```

```bash
git add GTO_ELFO/indirect/README.md
```

- [ ] **Step 4: Verify layout, renames-only staging, history**

```bash
ls GTO_tulip/            # expect: attic  CODE_CLEANUP_PLAN.md  direct  doc  HONEST_...md  indirect  LOW_THRUST_...md  MIN_ENERGY_NOTES.md  README.md  reviews  ROADMAP.md  sundman_minfuel_solution_note.{pdf,tex}
ls GTO_tulip/direct GTO_tulip/indirect GTO_ELFO/direct GTO_ELFO/indirect cr3bp_common
cd /Users/msc/Desktop/optimal_control
git diff --cached --name-status | grep -cv "^R"     # expect 1 (only the new README is non-rename)
git status --short | grep -E "^[AM] " | grep -v "GTO_ELFO/indirect/README.md" | head   # expect empty
git log --oneline --follow -1 -- orbit_transfer/cr3bp_common/cr3bp_lt_params.m
git log --oneline --follow -1 -- orbit_transfer/GTO_tulip/indirect/ms_band/setup_paths.m
```
Expected: `GTO_tulip/` root holds only docs + `attic/doc/reviews` + `direct/` + `indirect/`; staged changes are all renames plus exactly 1 add; `git log --follow` prints history for both samples.

- [ ] **Step 5: Commit (clearly labeled part 1 of 2)**

```bash
git commit -m "refactor(gto): part 1/2 — structural direct/indirect split (moves only)

GTO_tulip -> {direct: sundman_minfuel, PSR(+PSR_data), movie} +
{indirect: lowThrust_GTO_tulip, ms_band, ifs(+IFS_data), ztl, min_time};
elfo -> new GTO_ELFO/direct/; shared params (cr3bp_lt_params, minfuel_config,
gto_tulip_endpoints, gto_elfo_endpoints) -> new cr3bp_common/. Pure git mv,
history preserved. NOTE: setup_paths are rewired in the NEXT commit — the
moved modules are intentionally non-runnable at exactly this commit.

Co-Authored-By: Claude Fable 5 <noreply@anthropic.com>"
```

---

### Task 2: cr3bp_common library + setup_paths rewiring + full smoke matrix

**Files:**
- Create: `orbit_transfer/cr3bp_common/setup_cr3bp_common.m`
- Rewrite: `setup_paths.m` in `GTO_tulip/direct/sundman_minfuel`, `GTO_tulip/indirect/{lowThrust_GTO_tulip, ms_band, ifs, ztl, min_time}`, `GTO_ELFO/direct/elfo` (7 files; **PSR's is untouched**)
- Modify (enumerated ref fixes): `GTO_tulip/direct/movie/gen_movie_data.m`, `GTO_tulip/direct/sundman_minfuel/test_insertion_states.m`, `GTO_tulip/doc/figures/make_figs.m`, `GTO_tulip/attic/{run_sundman_minfuel.m, tf_continuation_minfuel.m, tf_continuation_minfuel_fine.m}`, `GTO_ELFO/direct/elfo/attic/gen_elfo_energy_tangential.m`, `GTO_tulip/direct/PSR/{psr_batch.sh, run_psr.m}`, `GTO_ELFO/direct/elfo/{elfo_batch.sh, elfo_energy_sweep.sh, elfo_movies.sh}`

**Interfaces:**
- Consumes: Task 1's tree.
- Produces: `setup_cr3bp_common()` — no args, adds `cr3bp_common` + pumpkyn to path, asserts pumpkyn exists. Every rewritten `setup_paths.m` keeps the exact name/signature `function setup_paths()` (cwd-shadowing contract).

- [ ] **Step 1: Write `cr3bp_common/setup_cr3bp_common.m`**

```matlab
function setup_cr3bp_common()
% SETUP_CR3BP_COMMON  Add the shared CR3BP GTO-transfer library + pumpkyn.
%
% Single source of the cross-module problem definition (cr3bp_lt_params,
% minfuel_config, gto_tulip_endpoints, gto_elfo_endpoints) and the pumpkyn
% toolbox path. Called by every GTO module's setup_paths.m.
%
% INPUTS:  (none)
% OUTPUTS: (none) - modifies the MATLAB path in-place
%
% REFERENCES:
%   [1] docs/superpowers/specs/2026-07-21-gto-direct-indirect-restructure-design.md
here = fileparts(mfilename('fullpath'));
addpath(here);
pumpkynSrc = fullfile(getenv('HOME'), 'Desktop', 'proj7', 'external', ...
                      'pumpkyn', 'src');
assert(exist(fullfile(pumpkynSrc, '+pumpkyn'), 'dir') == 7, ...
    'setup_cr3bp_common:missing', 'pumpkyn not found at %s', pumpkynSrc);
addpath(pumpkynSrc);
end
```

- [ ] **Step 2: Rewrite the 7 module `setup_paths.m` files** (complete bodies; keep any existing header comment style, update the described adds):

`GTO_tulip/direct/sundman_minfuel/setup_paths.m`:
```matlab
function setup_paths()
% SETUP_PATHS  Direct Sundman min-fuel engine paths: self + shared CR3BP lib
% (cr3bp_common: params/config/endpoints + pumpkyn).
%
% INPUTS: (none)   OUTPUTS: (none) - modifies the MATLAB path in-place
here = fileparts(mfilename('fullpath'));
addpath(here);
addpath(fullfile(here, '..', '..', '..', 'cr3bp_common'));
setup_cr3bp_common();
end
```

`GTO_tulip/indirect/lowThrust_GTO_tulip/setup_paths.m`:
```matlab
function setup_paths()
% SETUP_PATHS  Base indirect GTO->tulip campaign paths: self + shared CR3BP lib
% (cr3bp_common brings pumpkyn, which this campaign's solvers use).
%
% INPUTS: (none)   OUTPUTS: (none) - modifies the MATLAB path in-place
here = fileparts(mfilename('fullpath'));
addpath(here);
addpath(fullfile(here, '..', '..', '..', 'cr3bp_common'));
setup_cr3bp_common();
end
```

`GTO_tulip/indirect/min_time/setup_paths.m`:
```matlab
function setup_paths()
% SETUP_PATHS  min_time (PMP min-time root) paths: self + shared CR3BP lib
% (gto_tulip_endpoints + pumpkyn via cr3bp_common).
%
% INPUTS: (none)   OUTPUTS: (none) - modifies the MATLAB path in-place
here = fileparts(mfilename('fullpath'));
addpath(here);
addpath(fullfile(here, '..', '..', '..', 'cr3bp_common'));
setup_cr3bp_common();
end
```

`GTO_tulip/indirect/ztl/setup_paths.m`:
```matlab
function setup_paths()
% SETUP_PATHS  ztl paths: the base indirect campaign (now a sibling under
% indirect/) + the shared CR3BP lib (params + pumpkyn, asserted there).
%
% INPUTS: (none)   OUTPUTS: (none) - modifies the MATLAB path in-place
here = fileparts(mfilename('fullpath'));
addpath(here);
oldCampaign = fullfile(here, '..', 'lowThrust_GTO_tulip');
assert(exist(fullfile(oldCampaign, 'lt_pmp_eom_minfuel.m'), 'file') == 2, ...
    'setup_paths:missing', 'indirect campaign not found at %s', oldCampaign);
addpath(oldCampaign);
addpath(fullfile(here, '..', '..', '..', 'cr3bp_common'));
setup_cr3bp_common();
end
```

`GTO_tulip/indirect/ms_band/setup_paths.m`:
```matlab
function setup_paths()
% SETUP_PATHS  ms_band paths.
% Adds: ../lowThrust_GTO_tulip (lt_pmp_eom* indirect EOM), the DIRECT
% sundman_minfuel engine (cross-method edge: dual-.mat seed data under its
% results/ + solver helpers), and the shared CR3BP lib (params + pumpkyn).
%
% INPUTS: (none)   OUTPUTS: (none) - modifies the MATLAB path in-place
here = fileparts(mfilename('fullpath'));
addpath(fullfile(here, '..', 'lowThrust_GTO_tulip'));
addpath(fullfile(here, '..', '..', 'direct', 'sundman_minfuel'));
addpath(fullfile(here, '..', '..', '..', 'cr3bp_common'));
setup_cr3bp_common();
end
```

`GTO_tulip/indirect/ifs/setup_paths.m`:
```matlab
function setup_paths()
% SETUP_PATHS  IFS paths.
% Adds: the DIRECT sundman_minfuel engine + its refine/ (cross-method edge:
% prep_refine_seed etc.), ../ms_band, ../lowThrust_GTO_tulip (lt_pmp_eom*),
% and the shared CR3BP lib (params + pumpkyn).
%
% INPUTS: (none)   OUTPUTS: (none) - modifies the MATLAB path in-place
here = fileparts(mfilename('fullpath'));
addpath(fullfile(here, '..', '..', 'direct', 'sundman_minfuel'));
addpath(fullfile(here, '..', '..', 'direct', 'sundman_minfuel', 'refine'));
addpath(fullfile(here, '..', 'ms_band'));
addpath(fullfile(here, '..', 'lowThrust_GTO_tulip'));
addpath(fullfile(here, '..', '..', '..', 'cr3bp_common'));
setup_cr3bp_common();
end
```

`GTO_ELFO/direct/elfo/setup_paths.m`:
```matlab
function setup_paths()
% SETUP_PATHS  Direct GTO->ELFO campaign paths.
% Adds: self, the tulip direct Sundman engine (cross-problem edge: this
% campaign reuses casadi_minfuel_sundman / insertion_states / minfuel_at_tf,
% retargeted to ELFO), and the shared CR3BP lib (params, gto_elfo_endpoints,
% pumpkyn).
%
% INPUTS: (none)   OUTPUTS: (none) - modifies the MATLAB path in-place
here = fileparts(mfilename('fullpath'));
addpath(here);
addpath(fullfile(here, '..', '..', '..', 'GTO_tulip', 'direct', 'sundman_minfuel'));
addpath(fullfile(here, '..', '..', '..', 'cr3bp_common'));
setup_cr3bp_common();
end
```

- [ ] **Step 3: Enumerated executable/comment ref fixes** (exact edits):

(a) `GTO_tulip/direct/movie/gen_movie_data.m` — line 7 `.../orbit_transfer/lowThrust_GTO_tulip` → `.../orbit_transfer/GTO_tulip/indirect/lowThrust_GTO_tulip`; line 9 `.../orbit_transfer/GTO_tulip/sundman_minfuel/setup_paths.m` → `.../orbit_transfer/GTO_tulip/direct/sundman_minfuel/setup_paths.m`. Line 8 (`addpath(.../orbit_transfer/GTO_tulip)`) still names an existing dir — leave it.
(b) `GTO_tulip/direct/sundman_minfuel/test_insertion_states.m` line 4: `addpath('../elfo')` → `addpath('../../../GTO_ELFO/direct/elfo')` (script `cd`s to sundman_minfuel first).
(c) `GTO_tulip/doc/figures/make_figs.m` line 6: `addpath('../../sundman_minfuel'); addpath('../../elfo');` → `addpath('../../direct/sundman_minfuel'); addpath('../../../GTO_ELFO/direct/elfo'); addpath('../../../cr3bp_common');` (it calls `gto_elfo_endpoints`, which now lives in cr3bp_common).
(d) `GTO_tulip/attic/run_sundman_minfuel.m:9`, `tf_continuation_minfuel.m:39`, `tf_continuation_minfuel_fine.m:28`: `fullfile(here,'..','lowThrust_GTO_tulip')` → `fullfile(here,'..','indirect','lowThrust_GTO_tulip')` (attic stays at GTO_tulip root).
(e) `GTO_ELFO/direct/elfo/attic/gen_elfo_energy_tangential.m:35`: `addpath(fullfile(here,'..','sundman_minfuel'),'-end')` (pre-existing broken ref; intent = the Sundman engine) → `addpath(fullfile(here,'..','..','..','..','GTO_tulip','direct','sundman_minfuel'),'-end')` after confirming `here` in that script is the file's own attic dir. Leave line 36 (`'..','attic'`) alone.
(f) `.sh` how-to-run headers + `run_psr.m` comments (sed, idempotent):
```bash
cd /Users/msc/Desktop/optimal_control
sed -E -i '' 's#orbit_transfer/GTO_tulip/(PSR_data|PSR)#orbit_transfer/GTO_tulip/direct/\1#g' \
  orbit_transfer/GTO_tulip/direct/PSR/psr_batch.sh orbit_transfer/GTO_tulip/direct/PSR/run_psr.m
sed -E -i '' 's#orbit_transfer/GTO_tulip/elfo#orbit_transfer/GTO_ELFO/direct/elfo#g' \
  orbit_transfer/GTO_ELFO/direct/elfo/elfo_batch.sh \
  orbit_transfer/GTO_ELFO/direct/elfo/elfo_energy_sweep.sh \
  orbit_transfer/GTO_ELFO/direct/elfo/elfo_movies.sh
```

- [ ] **Step 4: Run the FULL smoke matrix (one MATLAB invocation, 8 modules)**

Run:
```
/Applications/MATLAB_R2025b.app/bin/matlab -batch "r='/Users/msc/Desktop/optimal_control/orbit_transfer'; cd(fullfile(r,'GTO_tulip','direct','sundman_minfuel')); setup_paths; assert(exist('cr3bp_lt_params','file')==2&&exist('minfuel_config','file')==2&&exist('gto_tulip_endpoints','file')==2&&exist('gto_elfo_endpoints','file')==2,'sundman'); cd(fullfile(r,'GTO_tulip','direct','PSR')); setup_paths; assert(exist('cr3bp_lt_params','file')==2&&exist('psr_run_one','file')==2,'PSR'); cd(fullfile(r,'GTO_ELFO','direct','elfo')); setup_paths; assert(exist('casadi_minfuel_sundman','file')==2&&exist('insertion_states','file')==2&&exist('minfuel_at_tf','file')==2&&exist('gto_elfo_endpoints','file')==2,'elfo'); cd(fullfile(r,'GTO_tulip','indirect','lowThrust_GTO_tulip')); setup_paths; assert(exist('lt_pmp_eom_minfuel','file')==2&&exist('cr3bp_lt_params','file')==2,'lowThrust'); cd(fullfile(r,'GTO_tulip','indirect','ms_band')); setup_paths; assert(exist('lt_pmp_eom_minfuel','file')==2&&exist('cr3bp_lt_params','file')==2&&exist('sundman_seed_map','file')==2,'ms_band'); cd(fullfile(r,'GTO_tulip','indirect','ifs')); setup_paths; assert(exist('prep_refine_seed','file')==2&&exist('lt_pmp_eom_minfuel','file')==2,'ifs'); cd(fullfile(r,'GTO_tulip','indirect','ztl')); setup_paths; assert(exist('lt_pmp_eom_minfuel','file')==2,'ztl'); cd(fullfile(r,'GTO_tulip','indirect','min_time')); setup_paths; assert(exist('gto_tulip_endpoints','file')==2&&exist('mintime_ms_residual','file')==2,'min_time'); disp('GTO RESTRUCTURE SMOKE MATRIX PASS')"
```
Expected: prints `GTO RESTRUCTURE SMOKE MATRIX PASS`. Any assert names the failing module — STOP and report BLOCKED if so. **Do not commit until it passes.**

- [ ] **Step 5: Stage exactly the named files and commit**

```bash
cd /Users/msc/Desktop/optimal_control
git add orbit_transfer/cr3bp_common/setup_cr3bp_common.m \
  orbit_transfer/GTO_tulip/direct/sundman_minfuel/setup_paths.m \
  orbit_transfer/GTO_tulip/direct/sundman_minfuel/test_insertion_states.m \
  orbit_transfer/GTO_tulip/direct/movie/gen_movie_data.m \
  orbit_transfer/GTO_tulip/direct/PSR/psr_batch.sh \
  orbit_transfer/GTO_tulip/direct/PSR/run_psr.m \
  orbit_transfer/GTO_tulip/indirect/lowThrust_GTO_tulip/setup_paths.m \
  orbit_transfer/GTO_tulip/indirect/ms_band/setup_paths.m \
  orbit_transfer/GTO_tulip/indirect/ifs/setup_paths.m \
  orbit_transfer/GTO_tulip/indirect/ztl/setup_paths.m \
  orbit_transfer/GTO_tulip/indirect/min_time/setup_paths.m \
  orbit_transfer/GTO_tulip/doc/figures/make_figs.m \
  orbit_transfer/GTO_tulip/attic/run_sundman_minfuel.m \
  orbit_transfer/GTO_tulip/attic/tf_continuation_minfuel.m \
  orbit_transfer/GTO_tulip/attic/tf_continuation_minfuel_fine.m \
  orbit_transfer/GTO_ELFO/direct/elfo/setup_paths.m \
  orbit_transfer/GTO_ELFO/direct/elfo/attic/gen_elfo_energy_tangential.m \
  orbit_transfer/GTO_ELFO/direct/elfo/elfo_batch.sh \
  orbit_transfer/GTO_ELFO/direct/elfo/elfo_energy_sweep.sh \
  orbit_transfer/GTO_ELFO/direct/elfo/elfo_movies.sh
git status --short          # confirm ONLY these 20 files staged
git commit -m "refactor(gto): part 2/2 — cr3bp_common library + setup_paths rewiring

New setup_cr3bp_common (shared params/endpoints + pumpkyn, single source);
7 module setup_paths rewired to the uniform self + cr3bp_common + explicit
cross-ref shape (PSR untouched — self-contained lib/). Enumerated executable
ref fixes (gen_movie_data, test_insertion_states, make_figs, attic x4, .sh
headers). Full 8-module resolve-smoke matrix PASSES at this commit.

Co-Authored-By: Claude Fable 5 <noreply@anthropic.com>"
```

---

### Task 3: Doc/comment reference sweep + repoint

**Files:**
- Modify: matches of the sweep below — expected mainly `orbit_transfer/earth_elliptic_to_geo/{README.md, TODO.md, process/*.md, cartesian_legacy/*.m, drivers/run_ctf_sweep.m, core/kepler_lt_params.m}`, `orbit_transfer/GTO_tulip/{README.md, LOW_THRUST_MINFUEL_CAMPAIGN.md, ROADMAP.md, CODE_CLEANUP_PLAN.md, MIN_ENERGY_NOTES.md, HONEST_EVALUATION_DV_TF_FRONT.md}`, `orbit_transfer/GTO_tulip/doc/**`, module-internal `*.md`

**Interfaces:**
- Consumes: the moved layout. Non-executable hygiene (Task 2 already fixed everything executable).

- [ ] **Step 1: Path-qualified sweep — find every stale qualified path in tracked files**

```bash
cd /Users/msc/Desktop/optimal_control
git grep -nI -E "GTO_tulip/(sundman_minfuel|PSR_data|PSR|movie|ms_band|ifs|IFS_data|ztl|min_time|elfo)|orbit_transfer/lowThrust_GTO_tulip" -- . \
  | grep -vE "GTO_tulip/(direct|indirect)/|GTO_ELFO/" \
  | grep -vE "docs/superpowers/(plans|specs)/|\.superpowers/"
```
Fix every hit by this mapping (sed or per-file edit; re-run until the grep is empty):
- `GTO_tulip/(sundman_minfuel|PSR_data|PSR|movie)` → `GTO_tulip/direct/\1`
- `GTO_tulip/(ms_band|ifs|IFS_data|ztl|min_time)` → `GTO_tulip/indirect/\1`
- `GTO_tulip/elfo` → `GTO_ELFO/direct/elfo`
- `orbit_transfer/lowThrust_GTO_tulip` → `orbit_transfer/GTO_tulip/indirect/lowThrust_GTO_tulip`
- relative refs from earth_elliptic (e.g. `../NLP...` are long gone; `../GTO_tulip/PSR/` etc.) are caught by the same patterns.

- [ ] **Step 2: Bare-name sweep in GTO_tulip root docs** (prose like "see `sundman_minfuel/`"):

```bash
git grep -nE "(^|[^/_a-zA-Z])(sundman_minfuel|ms_band|ifs|ztl|min_time|PSR|movie|elfo)/" -- \
  orbit_transfer/GTO_tulip/README.md orbit_transfer/GTO_tulip/LOW_THRUST_MINFUEL_CAMPAIGN.md \
  orbit_transfer/GTO_tulip/ROADMAP.md orbit_transfer/GTO_tulip/CODE_CLEANUP_PLAN.md \
  orbit_transfer/GTO_tulip/MIN_ENERGY_NOTES.md orbit_transfer/GTO_tulip/HONEST_EVALUATION_DV_TF_FRONT.md \
  | grep -vE "direct/|indirect/|GTO_ELFO/"
```
These are relative-from-GTO_tulip-root doc refs. Fix each: direct modules get `direct/` prefix, indirect get `indirect/`, `elfo/` → `../GTO_ELFO/direct/elfo/`. Judgment allowed on prose that names a module without a path sense (e.g. "the ms_band campaign" without a trailing `/` is fine as-is — the grep pattern requires the trailing `/`, so hits ARE path-like). Re-run until empty.

- [ ] **Step 3: Commit**

```bash
git add -u :^papers/.DS_Store
git status --short          # review: ONLY .md/.m doc-comment edits from steps 1-2; nothing unrelated
git commit -m "refactor(docs): repoint references to the GTO direct/indirect layout

Path-qualified refs (GTO_tulip/<module>, orbit_transfer/lowThrust_GTO_tulip)
and bare in-folder doc refs updated to direct/, indirect/, GTO_ELFO/direct/.
Historical docs/superpowers records left as point-in-time, per spec.

Co-Authored-By: Claude Fable 5 <noreply@anthropic.com>"
```
(NOTE: `git add -u` with the pathspec-exclude keeps `papers/.DS_Store` out; if any other unrelated modified file appears in `git status`, unstage it and stage explicitly instead.)

---

### Task 4: CLAUDE.md maps + final verification + push

**Files:**
- Modify: `CLAUDE.md` (project, repo-tracked)
- (Out-of-repo) `~/Desktop/CLAUDE.md` hub — update the Optimal Control row's `orbit_transfer/` paths

**Interfaces:**
- Consumes: the completed restructure. Documentation + gate + push.

- [ ] **Step 1: Update the project `CLAUDE.md` tree block** — replace the `orbit_transfer/` section of the directory-structure diagram with:

```
├── orbit_transfer/              # top-level container: all orbit-transfer work
│   ├── cr3bp_common/            # shared CR3BP GTO library: cr3bp_lt_params,
│   │                            #   minfuel_config, gto_{tulip,elfo}_endpoints,
│   │                            #   setup_cr3bp_common (adds pumpkyn)
│   ├── min_energy_tutorial/     # min-energy transfer tutorial
│   ├── lambert/                 # universal-variables Lambert
│   ├── earth_elliptic_to_geo/   # GTO->GEO min-fuel (HMG-2004 direct reproduction)
│   ├── earth_elliptic_to_geo_CR3BP/  # (stub) CR3BP variant
│   ├── GTO_tulip/               # CR3BP GTO->tulip; campaign docs at root
│   │   ├── direct/              #   sundman_minfuel (Sundman engine), PSR(+data), movie
│   │   └── indirect/            #   lowThrust_GTO_tulip (base PMP shooting),
│   │                            #   ms_band, ifs(+data), ztl, min_time
│   ├── GTO_ELFO/                # CR3BP GTO->ELFO
│   │   ├── direct/              #   elfo (CasADi campaign; reuses tulip Sundman engine)
│   │   └── indirect/            #   placeholder (Route C future work)
│   ├── mfmax/                   # MfMax v0/v1 Fortran (Gergaud-group indirect) + docs
│   ├── min_fuel_paper/          # paper outline (co-author Koblick)
│   └── min_fuel_papers/         # reference PDFs
```
Also update any prose in `CLAUDE.md` naming `GTO_tulip/<module>` or `lowThrust_GTO_tulip` at old paths (same mapping as Task 3).

- [ ] **Step 2: Update the hub map** — in `~/Desktop/CLAUDE.md` line ~18 (Optimal Control row): `orbit_transfer/lowThrust_GTO_tulip/` → `orbit_transfer/GTO_tulip/indirect/lowThrust_GTO_tulip/`, and note the direct/indirect split, e.g. `orbit_transfer/GTO_tulip/` `(CR3BP GTO→tulip, direct/ + indirect/ split 2026-07-21)` and `orbit_transfer/GTO_ELFO/` added. Keep edits minimal — path corrections + one clause. (Out of repo: no commit.)

- [ ] **Step 3: Final verification gate**

```bash
cd /Users/msc/Desktop/optimal_control
ls orbit_transfer/          # expect: cr3bp_common GTO_ELFO GTO_tulip earth_* lambert min_* mfmax
git grep -nI -E "GTO_tulip/(sundman_minfuel|PSR_data|PSR|movie|ms_band|ifs|IFS_data|ztl|min_time|elfo)|orbit_transfer/lowThrust_GTO_tulip" -- . | grep -vE "GTO_tulip/(direct|indirect)/|GTO_ELFO/|docs/superpowers/(plans|specs)/|\.superpowers/" | head
```
Expected: sweep prints nothing. Then re-run the **Task 2 Step 4 smoke matrix command verbatim** — expected `GTO RESTRUCTURE SMOKE MATRIX PASS`.

- [ ] **Step 4: Commit docs + push**

```bash
git add CLAUDE.md
git commit -m "docs(claude): directory map for GTO direct/indirect + cr3bp_common

Co-Authored-By: Claude Fable 5 <noreply@anthropic.com>"
git push origin main
git log --oneline -5
```
Expected: push succeeds; the 4 restructure commits are on top.

---

## Self-Review notes

- **Spec coverage:** §2 tree → T1; §3 classification (incl. PSR-intact, data-dir plain-mv) → T1 + Global Constraints; §4 library mechanism → T2 Steps 1–2; §5 dependency table → T2 Step 2 bodies match row-for-row (sundman: self+common; PSR: untouched; movie: gen_movie_data fix only — movie has no setup_paths.m, the spec's "extra dirs" row is realized via that file's absolute paths; elfo: +tulip sundman; lowThrust/min_time: self+common; ms_band/ifs/ztl: exact extra dirs); §5 honest cross-refs → documented in the rewritten headers; §6 external refs → T2 Step 3 (executable) + T3 (docs) + T4 (CLAUDE.md, historical carve-out); §7 verification → T2 Step 4 matrix (incl. resolved open-items: `prep_refine_seed` for ifs) + T4 Step 3 re-run + `git log --follow` in T1; §8 items 1–4 all realized (README stub = T1 Step 3); §9 sequencing risk → GTO_tulip and GTO_ELFO move in the same T1 commit, so the elfo→sundman edge never dangles across commits.
- **Type consistency:** every `setup_paths()` keeps its exact name/zero-arg signature; `setup_cr3bp_common()` is called only after its dir is addpath'd (same statement pair everywhere).
- **Placeholder scan:** all code steps carry complete bodies; all sed/grep commands are exact; the only judgment step (T3 Step 2 prose fixes) has an explicit mapping and an empty-grep termination condition.
- **Known intentional exception:** Task 1's commit is non-runnable for the moved modules (labeled in the commit message); T1+T2 are adjacent commits and T2 is gated on the full smoke matrix.
