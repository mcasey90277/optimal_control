# Orbit-Transfer Folder Consolidation — Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Make `orbit_transfer/` the top-level container for all orbit-transfer work: nest its current tutorial into `min_energy_tutorial/`, then move the five sibling campaign folders inside it, preserving git history and all path dependencies.

**Architecture:** Single git repo (no submodules) → every move is `git mv` (history preserved). Interdependent folders (`lowThrust_GTO_tulip` ↔ `NLP_lowThrust_GTO_tulip`) move together so their `../../lowThrust_GTO_tulip`-style relative refs stay valid. One executable absolute-path file and a handful of comment/doc links get repointed. Verified by running each moved folder's `setup_paths` + a cheap smoke before each commit.

**Tech Stack:** git, zsh, MATLAB R2025b (`/Applications/MATLAB_R2025b.app/bin/matlab`).

## Global Constraints

- **Repo root:** `/Users/msc/Desktop/optimal_control` (git remote `origin git@github.com:mcasey90277/optimal_control.git`, branch `main` tracks `origin/main`). All `git`/`ls` commands run from here.
- **Scope is EXACTLY five folders** → `orbit_transfer/`: `earth_elliptic_to_geo`, `earth_elliptic_to_geo_CR3BP`, `lambert`, `lowThrust_GTO_tulip`, `NLP_lowThrust_GTO_tulip`. **Do NOT touch** any other top-level folder (`ex1_block_move`, `ex2_cart_pole_swing_up`, `mpc_cart_pole`, `quasiNewton_matlab`, `lieFiltering`, `gauss_sum_curvature`, `mfmax-v0`, `mfmax-v1`, `mfmax_docs`, `min_fuel_paper`, `min_fuel_papers`, `papers`, `learning_docs`).
- **Subfolder name for the nested tutorial:** `min_energy_tutorial` (exact).
- **`lowThrust_GTO_tulip` and `NLP_lowThrust_GTO_tulip` MUST move in the same task** — moving one without the other breaks `ms_band/setup_paths.m`'s `../../lowThrust_GTO_tulip` addpath.
- **MATLAB invocation:** `/Applications/MATLAB_R2025b.app/bin/matlab -batch "<one-line cmd>"` (multi-line `-batch` gets mangled by zsh — keep it one line, or `run('/abs/path.m')`).
- **`earth_elliptic_to_geo_CR3BP` is an empty untracked stub** — `git mv` won't move it; use plain `mv`.
- **Do not commit `.mat`** (gitignored) or any untracked scratch. Stage only the moves + the specific edits named per task.
- **Commit-message co-author trailer:** `Co-Authored-By: Claude Opus 4.8 (1M context) <noreply@anthropic.com>`.

---

### Task 1: Nest the current tutorial into `orbit_transfer/min_energy_tutorial/`

**Files:**
- Create dir: `orbit_transfer/min_energy_tutorial/`
- Move (git mv): the 13 tracked items in `orbit_transfer/` (see step 1) + plain-`mv` the untracked `orbit_transfer_exercises.pdf`

**Interfaces:**
- Produces: `orbit_transfer/min_energy_tutorial/{solve_indirect.m, collocation_transfer.m, primer_check.m, shoot_residual.m, ocp_dynamics.m, two_body_accel.m, gravity_gradient.m, run_orbit_transfer.m, verify_checkpoints.m, mytry/, reviews/, orbit_transfer_exercises.{tex,pdf}, expected_result.png}` — the tutorial, internally self-referential (all files moved together, so intra-tutorial relative refs are unchanged).

- [ ] **Step 1: Create the subfolder and move the tracked items**

```bash
cd /Users/msc/Desktop/optimal_control
mkdir orbit_transfer/min_energy_tutorial
git mv orbit_transfer/collocation_transfer.m orbit_transfer/expected_result.png \
       orbit_transfer/gravity_gradient.m orbit_transfer/mytry orbit_transfer/ocp_dynamics.m \
       orbit_transfer/orbit_transfer_exercises.tex orbit_transfer/primer_check.m \
       orbit_transfer/reviews orbit_transfer/run_orbit_transfer.m orbit_transfer/shoot_residual.m \
       orbit_transfer/solve_indirect.m orbit_transfer/two_body_accel.m \
       orbit_transfer/verify_checkpoints.m orbit_transfer/min_energy_tutorial/
```

- [ ] **Step 2: Move the untracked PDF (if present)**

```bash
[ -f orbit_transfer/orbit_transfer_exercises.pdf ] && mv orbit_transfer/orbit_transfer_exercises.pdf orbit_transfer/min_energy_tutorial/ || echo "no untracked pdf"
```

- [ ] **Step 3: Verify nothing stray left behind and history is intact**

```bash
ls orbit_transfer/                      # expect ONLY: min_energy_tutorial/
git status --short | grep orbit_transfer | head
git log --oneline --follow -1 -- orbit_transfer/min_energy_tutorial/solve_indirect.m
```
Expected: `orbit_transfer/` contains only `min_energy_tutorial/`; `git status` shows renames (R); `git log --follow` prints a prior commit (history preserved).

- [ ] **Step 4: Smoke — the tutorial's own checkpoint verifier runs from the new location**

Run: `/Applications/MATLAB_R2025b.app/bin/matlab -batch "cd('/Users/msc/Desktop/optimal_control/orbit_transfer/min_energy_tutorial'); assert(exist('solve_indirect','file')==2 && exist('primer_check','file')==2 && exist('verify_checkpoints','file')==2); disp('tutorial functions resolve OK')"`
Expected: prints `tutorial functions resolve OK` (all intra-tutorial files resolve from the nested dir). *(If `verify_checkpoints` is quick and self-contained, optionally run it too — it should pass unchanged since every file moved together.)*

- [ ] **Step 5: Commit**

```bash
# git mv already staged the renames; stage only the untracked PDF (if it was moved)
git add orbit_transfer/min_energy_tutorial/orbit_transfer_exercises.pdf 2>/dev/null || true
git status --short          # confirm ONLY the intended renames (+pdf) are staged
git commit -m "refactor(orbit_transfer): nest min-energy tutorial into min_energy_tutorial/

Makes orbit_transfer/ a container for the transfer campaigns; the current
guided min-energy tutorial (exercises + reference solvers + mytry/) moves to
orbit_transfer/min_energy_tutorial/. Pure git mv, history preserved.

Co-Authored-By: Claude Opus 4.8 (1M context) <noreply@anthropic.com>"
```

---

### Task 2: Move the self-contained folders (`lambert`, `earth_elliptic_to_geo`, `earth_elliptic_to_geo_CR3BP`)

**Files:**
- Move (git mv): `lambert/`, `earth_elliptic_to_geo/` → `orbit_transfer/`
- Move (plain mv, empty untracked stub): `earth_elliptic_to_geo_CR3BP/` → `orbit_transfer/`

**Interfaces:**
- Consumes: `orbit_transfer/` (now a container, from Task 1).
- Produces: `orbit_transfer/lambert/`, `orbit_transfer/earth_elliptic_to_geo/`, `orbit_transfer/earth_elliptic_to_geo_CR3BP/`. Both real folders are self-contained (`setup_paths.m` adds only their own subdirs), so no path edits are needed.

- [ ] **Step 1: Move the two tracked folders + the empty stub**

```bash
cd /Users/msc/Desktop/optimal_control
git mv lambert earth_elliptic_to_geo orbit_transfer/
mv earth_elliptic_to_geo_CR3BP orbit_transfer/    # empty untracked stub -> plain mv
```

- [ ] **Step 2: Verify the moves + history**

```bash
ls orbit_transfer/                       # expect: min_energy_tutorial earth_elliptic_to_geo earth_elliptic_to_geo_CR3BP lambert
git status --short | grep -E "earth_elliptic|lambert" | head
git log --oneline --follow -1 -- orbit_transfer/earth_elliptic_to_geo/setup_paths.m
```
Expected: folders present under `orbit_transfer/`; renames (R) staged; `git log --follow` shows history.

- [ ] **Step 3: Smoke — earth_elliptic_to_geo `setup_paths` + a core function resolve from the new location**

Run: `/Applications/MATLAB_R2025b.app/bin/matlab -batch "cd('/Users/msc/Desktop/optimal_control/orbit_transfer/earth_elliptic_to_geo'); setup_paths; p=kepler_lt_params(10,1500,2000); fprintf('earth_elliptic_to_geo OK: Tmax=%.5f\n', p.Tmax)"`
Expected: prints `earth_elliptic_to_geo OK: Tmax=0.02974` (setup_paths + core resolve; no path errors).

- [ ] **Step 4: Commit**

```bash
# git mv already staged the two folder renames; the empty CR3BP stub is untracked (nothing to stage)
git status --short          # confirm ONLY the intended renames are staged (no unrelated untracked files)
git commit -m "refactor(structure): move lambert + earth_elliptic_to_geo(_CR3BP) into orbit_transfer/

Self-contained folders (setup_paths adds only their own subdirs); no path
edits needed. History preserved via git mv.

Co-Authored-By: Claude Opus 4.8 (1M context) <noreply@anthropic.com>"
```

---

### Task 3: Move the interdependent pair together + fix the one executable absolute-path file

**Files:**
- Move (git mv): `lowThrust_GTO_tulip/`, `NLP_lowThrust_GTO_tulip/` → `orbit_transfer/` (same command — they must move together)
- Modify: `orbit_transfer/NLP_lowThrust_GTO_tulip/movie/gen_movie_data.m:7-9`

**Interfaces:**
- Consumes: `orbit_transfer/` container.
- Produces: `orbit_transfer/lowThrust_GTO_tulip/`, `orbit_transfer/NLP_lowThrust_GTO_tulip/` as siblings. `ms_band/setup_paths.m`'s `addpath(fullfile(here,'..','..','lowThrust_GTO_tulip'))` now resolves to `orbit_transfer/lowThrust_GTO_tulip` — still valid because both moved together.

- [ ] **Step 1: Move both folders in one command**

```bash
cd /Users/msc/Desktop/optimal_control
git mv lowThrust_GTO_tulip NLP_lowThrust_GTO_tulip orbit_transfer/
```

- [ ] **Step 2: Fix the executable absolute paths in `gen_movie_data.m`**

Edit `orbit_transfer/NLP_lowThrust_GTO_tulip/movie/gen_movie_data.m`, inserting `orbit_transfer/` into each of the three paths:

```matlab
% line 7  (was: '/Users/msc/Desktop/optimal_control/lowThrust_GTO_tulip')
addpath('/Users/msc/Desktop/optimal_control/orbit_transfer/lowThrust_GTO_tulip');
% line 8  (was: '/Users/msc/Desktop/optimal_control/NLP_lowThrust_GTO_tulip')
addpath('/Users/msc/Desktop/optimal_control/orbit_transfer/NLP_lowThrust_GTO_tulip');
% line 9  (was: '/Users/msc/Desktop/optimal_control/NLP_lowThrust_GTO_tulip/sundman_minfuel/setup_paths.m')
run('/Users/msc/Desktop/optimal_control/orbit_transfer/NLP_lowThrust_GTO_tulip/sundman_minfuel/setup_paths.m');  % pumpkyn on path
```

- [ ] **Step 3: Verify moves + history**

```bash
ls orbit_transfer/                       # expect the pair now present
git log --oneline --follow -1 -- orbit_transfer/NLP_lowThrust_GTO_tulip/core/casadi_lt_mee.m 2>/dev/null || \
git log --oneline --follow -1 -- orbit_transfer/lowThrust_GTO_tulip/setup_paths.m
```
Expected: pair under `orbit_transfer/`; `git log --follow` shows history.

- [ ] **Step 4: Smoke — the critical cross-folder dependency still resolves**

Run: `/Applications/MATLAB_R2025b.app/bin/matlab -batch "cd('/Users/msc/Desktop/optimal_control/orbit_transfer/NLP_lowThrust_GTO_tulip/ms_band'); setup_paths; fprintf('ms_band->lowThrust resolves: lt_pmp_eom_minfuel exist=%d\n', exist('lt_pmp_eom_minfuel','file'))"`
Expected: prints `ms_band->lowThrust resolves: lt_pmp_eom_minfuel exist=2` (the `../../lowThrust_GTO_tulip` addpath found the EOM in the moved sibling). A `0` means the pair did not stay siblings — STOP and investigate.

- [ ] **Step 5: Commit**

```bash
# git mv staged the pair's renames; stage only the one edited file
git add orbit_transfer/NLP_lowThrust_GTO_tulip/movie/gen_movie_data.m
git status --short          # confirm ONLY the renames + gen_movie_data.m edit are staged
git commit -m "refactor(structure): move lowThrust + NLP_lowThrust GTO_tulip into orbit_transfer/

The interdependent CR3BP low-thrust pair moves together, so ms_band's
../../lowThrust_GTO_tulip addpath stays valid. Repoints the one executable
absolute-path file (movie/gen_movie_data.m) into orbit_transfer/.

Co-Authored-By: Claude Opus 4.8 (1M context) <noreply@anthropic.com>"
```

---

### Task 4: Repoint stale comment/doc cross-links + sweep for stragglers

**Files:**
- Modify: `orbit_transfer/NLP_lowThrust_GTO_tulip/sundman_minfuel/certify_minfuel_pmp.m:31,67`
- Modify: `orbit_transfer/NLP_lowThrust_GTO_tulip/sundman_minfuel/TIER1_PMP_CERTIFICATION_SCOPE.md:31,90`
- Modify: `orbit_transfer/NLP_lowThrust_GTO_tulip/PSR/run_psr.m:114-115` (commented examples)

**Interfaces:**
- Consumes: the moved layout. These are non-executable references (comments/docs) — updated for accuracy, not to unbreak execution.

- [ ] **Step 1: Fix the `primer_check.m` links** — `../../orbit_transfer/primer_check.m` → `../../min_energy_tutorial/primer_check.m` (from `orbit_transfer/NLP_lowThrust_GTO_tulip/sundman_minfuel/`, `../..` = `orbit_transfer/`, then `min_energy_tutorial/primer_check.m`).

```bash
cd /Users/msc/Desktop/optimal_control
sed -i '' 's#\.\./\.\./orbit_transfer/primer_check\.m#../../min_energy_tutorial/primer_check.m#g' \
  orbit_transfer/NLP_lowThrust_GTO_tulip/sundman_minfuel/certify_minfuel_pmp.m \
  orbit_transfer/NLP_lowThrust_GTO_tulip/sundman_minfuel/TIER1_PMP_CERTIFICATION_SCOPE.md
```

- [ ] **Step 2: Fix the commented seedSpec examples in `run_psr.m`** — insert `orbit_transfer/`.

```bash
sed -i '' 's#optimal_control/NLP_lowThrust_GTO_tulip/PSR_data#optimal_control/orbit_transfer/NLP_lowThrust_GTO_tulip/PSR_data#g' \
  orbit_transfer/NLP_lowThrust_GTO_tulip/PSR/run_psr.m
```

- [ ] **Step 3: Sweep for any remaining stale references in TRACKED files**

```bash
# (a) absolute paths to the OLD (pre-move) locations of the five folders:
git grep -n -E "optimal_control/(lambert|earth_elliptic_to_geo|lowThrust_GTO_tulip|NLP_lowThrust_GTO_tulip)/" -- '*.m' '*.sh' | \
  grep -v "optimal_control/orbit_transfer/"
# (b) relative refs to the OLD orbit_transfer tutorial content:
git grep -n -E "\.\./\.\./orbit_transfer/" -- '*.m' '*.md'
```
Expected: **no output** from either (all repointed). Any line that appears is a straggler — fix it the same way (insert `orbit_transfer/` for absolute paths; `orbit_transfer/` → `min_energy_tutorial/` for the tutorial content) and re-run until clean. *(Historical dated plan docs under `docs/superpowers/` may still mention old paths in prose — leave those; they are point-in-time records, not live references.)*

- [ ] **Step 4: Commit**

```bash
git add orbit_transfer/NLP_lowThrust_GTO_tulip/sundman_minfuel/certify_minfuel_pmp.m \
        orbit_transfer/NLP_lowThrust_GTO_tulip/sundman_minfuel/TIER1_PMP_CERTIFICATION_SCOPE.md \
        orbit_transfer/NLP_lowThrust_GTO_tulip/PSR/run_psr.m
git status --short          # confirm ONLY these three edited files are staged
git commit -m "refactor(docs): repoint cross-folder links after orbit_transfer consolidation

certify_minfuel_pmp/TIER1 comment links -> ../../min_energy_tutorial/primer_check.m;
run_psr commented seedSpec examples -> orbit_transfer/... ; swept tracked files
for stragglers.

Co-Authored-By: Claude Opus 4.8 (1M context) <noreply@anthropic.com>"
```

---

### Task 5: Update `CLAUDE.md` directory map + final verification, then push

**Files:**
- Modify: `CLAUDE.md` (project config — the directory-structure section and any folder references)
- (Out-of-repo, optional) `~/Desktop/CLAUDE.md` (hub project map's Optimal Control row)

**Interfaces:**
- Consumes: the completed move. Documentation only.

- [ ] **Step 1: Update the project `CLAUDE.md` directory structure**

Read `CLAUDE.md`'s "Directory Structure" section and re-nest the five folders under `orbit_transfer/`, renaming the current `orbit_transfer/` tutorial entry to `orbit_transfer/min_energy_tutorial/`. Example target for the tree block:

```
orbit_transfer/                  # top-level container: all orbit-transfer work
├── min_energy_tutorial/         # min-energy transfer tutorial (was orbit_transfer/)
├── lambert/                     # universal-variables Lambert
├── earth_elliptic_to_geo/       # GTO->GEO min-fuel (HMG-2004 direct reproduction)
├── earth_elliptic_to_geo_CR3BP/ # (stub) CR3BP variant
├── lowThrust_GTO_tulip/         # CR3BP GTO->tulip, indirect
└── NLP_lowThrust_GTO_tulip/     # CR3BP GTO->tulip, direct NLP
```
Update any prose paths in `CLAUDE.md` that name these folders at the old top level (e.g. ``` `NLP_lowThrust_GTO_tulip/...` ``` → ``` `orbit_transfer/NLP_lowThrust_GTO_tulip/...` ```).

- [ ] **Step 2: Note the hub `CLAUDE.md`**

`~/Desktop/CLAUDE.md` (the hub, outside this repo) has an "Optimal Control" project row listing these subfolders. Update its paths to the new `orbit_transfer/<sub>` layout, OR, if leaving it for the user, print a reminder:
```bash
echo "REMINDER: update ~/Desktop/CLAUDE.md Optimal Control row to orbit_transfer/<sub> paths"
```

- [ ] **Step 3: Final whole-repo verification**

```bash
cd /Users/msc/Desktop/optimal_control
ls                                  # top level: orbit_transfer present; the 5 folders GONE from top
ls orbit_transfer/                  # min_energy_tutorial + the 5
git status --short                  # only intended changes staged/committed; no stray deletes
```
Then a MATLAB cross-check that the two path-sensitive folders still work post-move:
Run: `/Applications/MATLAB_R2025b.app/bin/matlab -batch "cd('/Users/msc/Desktop/optimal_control/orbit_transfer/earth_elliptic_to_geo'); setup_paths; kepler_lt_params(10,1500,2000); cd('/Users/msc/Desktop/optimal_control/orbit_transfer/NLP_lowThrust_GTO_tulip/ms_band'); setup_paths; assert(exist('lt_pmp_eom_minfuel','file')==2); disp('POST-MOVE SMOKE PASS')"`
Expected: `POST-MOVE SMOKE PASS`.

- [ ] **Step 4: Commit the docs**

```bash
git add CLAUDE.md
git commit -m "docs(claude): update directory map for orbit_transfer consolidation

Co-Authored-By: Claude Opus 4.8 (1M context) <noreply@anthropic.com>"
```

- [ ] **Step 5: Push**

```bash
git push origin main
git log --oneline -5      # confirm the reorg commits are on top
```
Expected: `main -> main` pushed; the 5 reorg commits present.

---

## Self-Review notes

- **Spec coverage:** target layout (Tasks 1–3) ✓; interdependent-pair-together rule (Task 3, one `git mv` + the cross-folder smoke) ✓; 2 tracked abs-path files — `gen_movie_data.m` executable (Task 3), `run_psr.m` commented (Task 4) ✓; comment cross-links (Task 4) ✓; sweep for stragglers (Task 4 Step 3) ✓; CLAUDE.md project + hub (Task 5) ✓; history preservation (`git log --follow` checks each move task) ✓; scope guard — only the 5 folders touched (Global Constraints) ✓.
- **Verification is the "test":** each move task ends with a `setup_paths`/function-resolution smoke from the NEW location before its commit; the load-bearing one is Task 3 Step 4 (ms_band → lowThrust).
- **Empty-stub handling** (`earth_elliptic_to_geo_CR3BP`, plain `mv`) called out explicitly (Task 2 Step 1).
- **Reversibility:** every step is `git mv` in one repo; a bad move is undone by `git mv` back or `git reset --hard` before push. Push is the last step (Task 5 Step 5), after all smokes pass.
