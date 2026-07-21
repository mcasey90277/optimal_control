# ELFO Sibling Directory Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Consolidate the scattered GTO→ELFO direct min-fuel code into a self-contained `GTO_tulip/elfo/` sibling directory that mirrors PSR's role for the tulip target.

**Architecture:** A new `elfo/` deliverable directory holds every ELFO-specific file (solver, seed generators, drivers, endpoints, export/verify, movie). It references the two shared engine files (`cr3bp_lt_params`, `minfuel_config`) from `sundman_minfuel/` on the MATLAB path — single source of truth, no re-vendoring. A dedicated `elfo_movie.m` (copy+rename of the already-generalized `psr_movie`) removes the current `elfo→PSR` dependency. PSR and `min_time/` are untouched.

**Tech Stack:** MATLAB (R202x), CasADi 3.7.0 (`~/casadi-3.7.0`), git. No test framework — verification is MATLAB `which`-resolution, `grep` guards, and a headless smoke run via the `matlab-headless` skill.

## Global Constraints

- **Working directory:** `/Users/msc/Desktop/optimal_control/GTO_tulip` (call it `$ROOT`). Branch: `ifs-retarget`.
- **Preserve history:** every relocation uses `git mv`, never copy+delete.
- **Do NOT touch:** `PSR/` (except reading `psr_movie.m` to derive `elfo_movie.m`), `PSR/lib/`, `min_time/`, and the `cr3bp_lt_params.m` / `minfuel_config.m` masters in `sundman_minfuel/` (referenced, never copied).
- **No new drift surface:** do not copy `cr3bp_lt_params` or `minfuel_config` into `elfo/`.
- **Every commit leaves the final target runnable:** intermediate movie breakage between Task 3 and Task 4 is acceptable only because Task 4 immediately repairs it; do not stop between them.
- **MATLAB function header standard** (CLAUDE.md) applies to any new/edited function: purpose, INPUTS with sizes, OUTPUTS with sizes.
- **Never use `i`/`j` as loop indices** (imaginary unit in MATLAB).

---

### Task 1: Scaffold the `elfo/` directory

**Files:**
- Create: `$ROOT/elfo/setup_paths.m`
- Create: `$ROOT/elfo/README.md`
- Create: `$ROOT/elfo/results/.gitkeep`
- Create: `$ROOT/elfo/attic/.gitkeep`

**Interfaces:**
- Produces: `setup_paths()` — no args, no return; path side effect adding `elfo/`, `sundman_minfuel/`, and pumpkyn to the MATLAB path.

- [ ] **Step 1: Create the directory skeleton**

```bash
cd /Users/msc/Desktop/optimal_control/GTO_tulip
mkdir -p elfo/results elfo/attic
touch elfo/results/.gitkeep elfo/attic/.gitkeep
```

- [ ] **Step 2: Write `elfo/setup_paths.m`**

Model it on `PSR/setup_paths.m` but shared-path (adds `sundman_minfuel/` rather than a vendored `lib`). Full content:

```matlab
function setup_paths()
% SETUP_PATHS  Add the GTO->ELFO direct min-fuel pipeline's paths.
%
% elfo/ is a self-contained deliverable directory for the GTO->ELFO transfer,
% the sibling of PSR/ (which is the GTO->tulip deliverable). Unlike PSR -- which
% VENDORS its machinery into PSR/lib -- elfo/ uses a SHARED-PATH model: the two
% shared engine files (cr3bp_lt_params, minfuel_config) stay single-source in
% sundman_minfuel/ and are added to the path here. Nothing is copied, so there
% is no vendoring-drift surface. Everything ELFO-specific (the two-primary
% free-tf solver, seed generators, drivers, endpoints, export/verify, movie)
% lives in elfo/ itself.
%
% Paths added:
%   elfo             - the ELFO pipeline (solver, drivers, run_elfo_minfuel)
%   sundman_minfuel  - shared engine (cr3bp_lt_params, minfuel_config) + seed
%                      data under sundman_minfuel/results (referenced in place)
%   pumpkyn/src      - external tulip-construction toolbox (third-party, in proj7)
%
% CasADi is added by casadi_energy_freetf itself (CASADI_PATH env var or
% ~/casadi-3.7.0).
%
% INPUTS:  none
% OUTPUTS: none (path side effect)
here = fileparts(mfilename('fullpath'));
addpath(here);
addpath(fullfile(here, '..', 'sundman_minfuel'));
addpath(fullfile(getenv('HOME'), 'Desktop', 'proj7', 'external', 'pumpkyn', 'src'));
end
```

- [ ] **Step 3: Write `elfo/README.md`**

```markdown
# elfo/ — GTO -> ELFO direct min-fuel pipeline

Self-contained deliverable for the minimum-fuel low-thrust GTO -> ELFO transfer
in the Earth-Moon CR3BP. Sibling of `PSR/` (the GTO -> tulip deliverable).

## Model: shared-path, not vendored

Unlike PSR (which vendors a frozen machinery snapshot into `PSR/lib`), this
directory references the two shared engine files -- `cr3bp_lt_params` and
`minfuel_config` -- from `../sundman_minfuel` on the path (`setup_paths.m`).
Single source of truth, no drift surface. Tradeoff: elfo/ tracks the dev
library, so a dev edit to those two files can change ELFO results.

## Pipeline

1. `gen_elfo_energy_gravhom.m` -- build the min-ENERGY seed via the two-primary
   gravity-homotopy ladder on the free-tf solver `casadi_energy_freetf.m`.
   (`gen_elfo_energy_tfsweep.m` builds the tf-grid of energy seeds.)
2. `gen_elfo_minfuel.m` -- sharpen energy -> fuel (epsilon 1 -> 0) to bang-bang.
3. `run_elfo_minfuel.m` -- end-to-end entry: solve -> export -> verify -> movie.
   - `elfo_export_data.m` -- costates from the two-primary KKT duals.
   - `verify_elfo_seed.m` -- solver-free seed verification.
   - `elfo_movie.m` -- control movie (copy of PSR's generalized psr_movie).
4. `gto_elfo_endpoints.m` / `probe_elfo_target.m` -- ELFO endpoints + geometry.

## Smoke tests

- `smoke_energy_freetf.m` -- free-tf form reproduces the f1.20 backbone.
- `smoke_fixedtf.m` -- pinned-tf leg-0 conversion is well-posed (no drift).

## Data

Results (`.mat`, movies) land in `elfo/results/`. The seed-data reads and the
`../PSR_data` reference are shared stores kept in place. Dead-end early routes
(`gen_elfo_energy_backbone`, `gen_elfo_energy_tangential`) are in `attic/`.

Full design rationale: `docs/superpowers/specs/2026-07-14-elfo-sibling-dir-design.md`.
```

- [ ] **Step 4: Verify setup_paths resolves the shared files**

Run (matlab-headless skill):

```
matlab -batch "cd('/Users/msc/Desktop/optimal_control/GTO_tulip/elfo'); setup_paths; assert(~isempty(which('cr3bp_lt_params')),'cr3bp missing'); assert(~isempty(which('minfuel_config')),'minfuel_config missing'); disp('OK setup_paths')"
```

Expected: prints `OK setup_paths` (both shared files resolve to `sundman_minfuel/`). CasADi not required for this check.

- [ ] **Step 5: Commit**

```bash
cd /Users/msc/Desktop/optimal_control/GTO_tulip
git add elfo/setup_paths.m elfo/README.md elfo/results/.gitkeep elfo/attic/.gitkeep
git commit -m "elfo: scaffold sibling deliverable dir (setup_paths, README, results/, attic/)"
```

---

### Task 2: Create `elfo/elfo_movie.m` (copy+rename of the generalized renderer)

**Files:**
- Read: `$ROOT/PSR/psr_movie.m`
- Create: `$ROOT/elfo/elfo_movie.m`

**Interfaces:**
- Consumes: nothing from earlier tasks (self-contained renderer).
- Produces: `elfo_movie(solFile, outStem, titleStr, mode, bgTrace)` — identical signature to `psr_movie`; used by `run_elfo_minfuel` in Task 4.

Rationale: `psr_movie` was already generalized to render BOTH layouts (it branches on `isfield(S,'out')` for tulip-seed vs top-level free-tf ELFO layout, and takes a `bgTrace` backdrop arg). So `elfo_movie` is a faithful copy with the function renamed and the header made ELFO-primary. Do NOT surgically strip the tulip branch — it is a harmless best-effort `try/catch` default backdrop and removing it risks breaking rendering; the accepted design keeps the mild duplication.

- [ ] **Step 1: Copy the file and rename the function**

```bash
cd /Users/msc/Desktop/optimal_control/GTO_tulip
cp PSR/psr_movie.m elfo/elfo_movie.m
```

Then edit `elfo/elfo_movie.m`: change **only** the first line from
`function psr_movie(solFile, outStem, titleStr, mode, bgTrace)`
to
`function elfo_movie(solFile, outStem, titleStr, mode, bgTrace)`.

- [ ] **Step 2: Update the header to be ELFO-primary**

Replace the `% PSR_MOVIE  ...` header name line and the two REFERENCES lines so they describe the ELFO usage. Specifically:

Change the header title line
`% PSR_MOVIE  Control movie for a PSR (or any Sundman min-fuel) solution.`
to
`% ELFO_MOVIE  Control movie for a GTO->ELFO (or any Sundman min-fuel) solution.`

Change the REFERENCES block:
```matlab
% REFERENCES:
%   [1] movie/animate_sundman_minfuel.m (layout + Delta-V meter design).
%   [2] PSR/run_psr.m section 5 (pipeline caller).
```
to
```matlab
% REFERENCES:
%   [1] movie/animate_sundman_minfuel.m (layout + Delta-V meter design).
%   [2] elfo/run_elfo_minfuel.m stage 6 (pipeline caller).
%   [3] Copy of PSR/psr_movie.m (already generalized to render either transfer;
%       kept as a separate file so elfo/ has no dependency on PSR/).
```

Leave all rendering logic byte-identical.

- [ ] **Step 3: Verify it parses and resolves**

Run (matlab-headless skill):

```
matlab -batch "cd('/Users/msc/Desktop/optimal_control/GTO_tulip/elfo'); setup_paths; assert(~isempty(which('elfo_movie')),'elfo_movie not on path'); nargin('elfo_movie'); disp('OK elfo_movie parses')"
```

Expected: prints `OK elfo_movie parses` (function is found and its signature reads with 5 inputs; no parse error).

- [ ] **Step 4: Commit**

```bash
cd /Users/msc/Desktop/optimal_control/GTO_tulip
git add elfo/elfo_movie.m
git commit -m "elfo: add elfo_movie (copy+rename of generalized psr_movie; no PSR dependency)"
```

---

### Task 3: Relocate all ELFO files into `elfo/` (git mv)

**Files:**
- Move (from `sundman_minfuel/` to `elfo/`): `casadi_energy_freetf.m`, `gen_elfo_energy_gravhom.m`, `gen_elfo_energy_tfsweep.m`, `gen_elfo_minfuel.m`, `run_elfo_minfuel.m`, `elfo_export_data.m`, `verify_elfo_seed.m`, `smoke_energy_freetf.m`, `smoke_fixedtf.m`
- Move (from `PSR/` to `elfo/`): `gto_elfo_endpoints.m`, `probe_elfo_target.m`
- Move (from `PSR/` to `elfo/attic/`): `gen_elfo_energy_backbone.m`, `gen_elfo_energy_tangential.m`
- Move (from `sundman_minfuel/results/` to `elfo/results/`): `energy_elfo_freetf.mat`, `minfuel_elfo.mat`, `movie_ELFO_tf1p200_minEps0.gif`, `movie_ELFO_tf1p200_minEps0.mp4`

**Interfaces:**
- Consumes: `elfo/` and `elfo/attic/` and `elfo/results/` exist (Task 1).
- Produces: all ELFO source/data files reside under `elfo/`. Internal cross-calls (e.g. `gen_elfo_minfuel`→`casadi_energy_freetf`, `gen_elfo_energy_gravhom`→`gto_elfo_endpoints`) resolve automatically since callers and callees moved together onto the same path.

This is one atomic relocation task — splitting per-file adds no reviewer value. Verified pre-move: the only external caller of any moved file is `smoke_fixedtf.m` calling `casadi_energy_freetf`, and it is itself in the move set.

- [ ] **Step 1: git mv the source files from sundman_minfuel/**

```bash
cd /Users/msc/Desktop/optimal_control/GTO_tulip
for f in casadi_energy_freetf gen_elfo_energy_gravhom gen_elfo_energy_tfsweep \
         gen_elfo_minfuel run_elfo_minfuel elfo_export_data verify_elfo_seed \
         smoke_energy_freetf smoke_fixedtf; do
  git mv "sundman_minfuel/$f.m" "elfo/$f.m"
done
```

- [ ] **Step 2: git mv the endpoint/probe files from PSR/**

```bash
git mv PSR/gto_elfo_endpoints.m elfo/gto_elfo_endpoints.m
git mv PSR/probe_elfo_target.m  elfo/probe_elfo_target.m
```

- [ ] **Step 3: git mv the two dead-end routes into attic/**

```bash
git mv PSR/gen_elfo_energy_backbone.m   elfo/attic/gen_elfo_energy_backbone.m
git mv PSR/gen_elfo_energy_tangential.m elfo/attic/gen_elfo_energy_tangential.m
```

- [ ] **Step 4: move the ELFO result files (plain `mv` — they are gitignored)**

These 4 files are gitignored by extension (`*.mat`, `*.mp4`, `*.gif`), so they
are untracked — `git mv` refuses them and there is no history to preserve. The
destination `elfo/results/` is covered by the same rules, so they stay ignored.
Use plain `mv`:

```bash
mv sundman_minfuel/results/energy_elfo_freetf.mat        elfo/results/energy_elfo_freetf.mat
mv sundman_minfuel/results/minfuel_elfo.mat             elfo/results/minfuel_elfo.mat
mv sundman_minfuel/results/movie_ELFO_tf1p200_minEps0.gif elfo/results/movie_ELFO_tf1p200_minEps0.gif
mv sundman_minfuel/results/movie_ELFO_tf1p200_minEps0.mp4 elfo/results/movie_ELFO_tf1p200_minEps0.mp4
```

These moves are filesystem-only and will NOT appear in the Step 6 commit
(untracked → correct).

- [ ] **Step 5: Verify the moves and that nothing external now dangles**

```bash
cd /Users/msc/Desktop/optimal_control/GTO_tulip
# 5a: all 11 source files present in elfo/
ls elfo/*.m
# 5b: the two dead routes in attic
ls elfo/attic/*.m
# 5c: NO ELFO source .m files left behind in sundman_minfuel/ or PSR/
echo "--- leftovers (expect none) ---"
ls sundman_minfuel/*elfo*.m PSR/*elfo*.m sundman_minfuel/casadi_energy_freetf.m sundman_minfuel/smoke_fixedtf.m 2>&1 | grep -v "No such file" || echo "none-left-behind"
# 5d: git sees them as renames
git status --short | grep -E "^R" | head -20
```

Expected: 5a lists 11 files (the 9 from sundman_minfuel + gto_elfo_endpoints + probe_elfo_target, minus none — plus elfo_movie/setup_paths already there); 5b lists 2 files; 5c prints `none-left-behind`; 5d shows rename (`R`) entries.

- [ ] **Step 6: Commit**

```bash
git commit -m "elfo: relocate ELFO pipeline into elfo/ (git mv from sundman_minfuel + PSR)"
```

---

### Task 4: Retarget `run_elfo_minfuel.m` internals

**Files:**
- Modify: `$ROOT/elfo/run_elfo_minfuel.m` (the `here`/path/movie lines — around 35 and 141)

**Interfaces:**
- Consumes: `elfo_movie(solFile, outStem, titleStr, mode, bgTrace)` (Task 2); `elfo/setup_paths` (Task 1); `resDir = fullfile(here,'results')` now resolves to `elfo/results` because the file lives in `elfo/`.
- Produces: a `run_elfo_minfuel` with zero references to `PSR/` code (only the in-place `../PSR_data` DATA reference remains).

- [ ] **Step 1: Fix the path/setup line**

In `elfo/run_elfo_minfuel.m`, find (around line 35):

```matlab
cd(here);  setup_paths();  addpath(fullfile(here,'..','PSR'));
```

Replace with (drop the `addpath(../PSR)` — `elfo_movie` is local now; `setup_paths` already adds the shared engine):

```matlab
cd(here);  setup_paths();
```

- [ ] **Step 2: Fix the movie call**

In the same file, find (around line 141):

```matlab
    psr_movie(outFile, movieStem, titleStr, movieMode, elfoTrace(:,1:3));
```

Replace with:

```matlab
    elfo_movie(outFile, movieStem, titleStr, movieMode, elfoTrace(:,1:3));
```

- [ ] **Step 3: Confirm the DATA reference is unchanged**

Verify line ~110 still reads (leave AS-IS — `PSR_data` is a shared data store, not code):

```matlab
dataDir = fullfile(here, '..', 'PSR_data');
```

No edit. This is intentional per the spec (in-place data reference).

- [ ] **Step 4: Grep guard — no PSR *code* references remain**

```bash
cd /Users/msc/Desktop/optimal_control/GTO_tulip
echo "--- expect: no psr_movie, no addpath ...PSR' (code); PSR_data data ref is OK ---"
grep -n "psr_movie" elfo/run_elfo_minfuel.m && echo "FAIL: psr_movie still referenced" || echo "OK: no psr_movie"
grep -n "addpath.*'PSR'" elfo/run_elfo_minfuel.m && echo "FAIL: addpath PSR still present" || echo "OK: no addpath PSR"
grep -n "PSR_data" elfo/run_elfo_minfuel.m && echo "(PSR_data data reference retained — expected)"
```

Expected: `OK: no psr_movie`, `OK: no addpath PSR`, and the `PSR_data` line still present (expected).

- [ ] **Step 5: Commit**

```bash
git add elfo/run_elfo_minfuel.m
git commit -m "elfo: retarget run_elfo_minfuel internals (drop addpath PSR, psr_movie -> elfo_movie)"
```

---

### Task 5: Full clean-path verification gate

**Files:** none modified (acceptance gate).

**Interfaces:**
- Consumes: everything from Tasks 1-4.
- Produces: evidence that the reorg is complete and the ELFO pipeline resolves and runs from a clean path.

- [ ] **Step 1: Clean-path `which`-resolution of every pipeline function**

Run (matlab-headless skill):

```
matlab -batch "restoredefaultpath; cd('/Users/msc/Desktop/optimal_control/GTO_tulip/elfo'); setup_paths; fns={'casadi_energy_freetf','gen_elfo_energy_gravhom','gen_elfo_energy_tfsweep','gen_elfo_minfuel','run_elfo_minfuel','elfo_export_data','verify_elfo_seed','smoke_energy_freetf','smoke_fixedtf','gto_elfo_endpoints','probe_elfo_target','elfo_movie','cr3bp_lt_params','minfuel_config'}; for k=1:numel(fns), w=which(fns{k}); assert(~isempty(w),['MISSING ' fns{k}]); assert(isempty(strfind(w,'/PSR/')),['RESOLVES TO PSR: ' fns{k} ' -> ' w]); fprintf('%-26s %s\n', fns{k}, w); end; disp('OK all resolve, none via PSR')"
```

Expected: every function prints a path under `elfo/` or `sundman_minfuel/`; none under `/PSR/`; final line `OK all resolve, none via PSR`. (`casadi_energy_freetf`…`probe_elfo_target`, `elfo_movie` → `elfo/`; `cr3bp_lt_params`, `minfuel_config` → `sundman_minfuel/`.)

- [ ] **Step 2: Tree-wide grep guard for stale references**

```bash
cd /Users/msc/Desktop/optimal_control/GTO_tulip
echo "--- any code still pointing at the OLD locations? (expect none) ---"
grep -rn "PSR/psr_movie" elfo/ && echo "FAIL" || echo "OK: elfo has no PSR/psr_movie ref"
grep -rn "sundman_minfuel/casadi_energy_freetf\|sundman_minfuel/gen_elfo\|sundman_minfuel/run_elfo\|sundman_minfuel/elfo_export\|sundman_minfuel/verify_elfo\|sundman_minfuel/smoke_energy_freetf\|sundman_minfuel/smoke_fixedtf" . --include='*.m' && echo "FAIL: stale hardcoded path" || echo "OK: no stale hardcoded ELFO paths"
```

Expected: `OK: elfo has no PSR/psr_movie ref` and `OK: no stale hardcoded ELFO paths`.

- [ ] **Step 3: Live smoke run (requires CasADi)**

Run (matlab-headless skill; needs `~/casadi-3.7.0`):

```
matlab -batch "cd('/Users/msc/Desktop/optimal_control/GTO_tulip/elfo'); setup_paths; smoke_energy_freetf"
```

Expected: `smoke_energy_freetf` reports the free-tf form reproducing the f1.20 backbone at machine precision (its existing success print, e.g. `maxDefect` ~1e-8, `ok=1`). If CasADi is unavailable in the run environment, record that Step 3 was skipped and note it explicitly — do NOT claim the live check passed.

- [ ] **Step 4: Final consistency check of the working tree**

```bash
cd /Users/msc/Desktop/optimal_control/GTO_tulip
git status --short
git log --oneline -5
```

Expected: clean or only-intended changes; the last 4-5 commits are the elfo scaffold/movie/relocate/retarget commits. No uncommitted ELFO source left dangling.

- [ ] **Step 5: (No commit — gate only.)** If any step failed, return to the responsible task. If all pass, the reorg is complete.

---

## Follow-up (optional, out of the core reorg — do only if requested)

- Update `CLAUDE.md` (project instructions) directory map / the
  `GTO_tulip` description to mention `elfo/` alongside `PSR/`.
- Update the ELFO memory `elfo-retarget-open.md` file:line references from
  `sundman_minfuel/` to `elfo/`.

These touch docs/memory outside the code reorg; keep them a separate commit if done.

## Self-Review

**Spec coverage:** §Target structure (all moves) → Tasks 1-3; §Edits required (run_elfo retargets) → Task 4; §movie decision → Task 2; §shared-path model → Task 1 setup_paths; §Attic → Task 3 Steps 3; §Results move → Task 3 Step 4; §Verification (which/grep/smoke) → Task 5; §Git (git mv, focused commits) → all tasks; §Out of scope (PSR, min_time untouched) → Global Constraints. `smoke_fixedtf.m` addition (found during planning) → Task 3 Step 1. All covered.

**Placeholder scan:** No TBD/TODO; every code step shows full content; every command shows expected output. Clean.

**Type/name consistency:** `elfo_movie(solFile, outStem, titleStr, mode, bgTrace)` defined in Task 2, called with the matching 5-arg signature in Task 4. `setup_paths()` defined Task 1, used in Tasks 4-5. Function names in the Task 5 `fns` list match the moved filenames in Task 3. Consistent.
