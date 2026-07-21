# GTO Direct/Indirect Restructure — Design Spec

**Date:** 2026-07-21
**Status:** design approved (brainstorm); implementation plan to follow.
**Scope:** the two CR3BP GTO-transfer campaigns under `orbit_transfer/` — `GTO_tulip`
(formerly `NLP_lowThrust_GTO_tulip` + the sibling `lowThrust_GTO_tulip`) and a new
`GTO_ELFO`. Earth-elliptic→GEO's own direct/indirect split is a **separate later plan**.

---

## 1. Goal & motivation

Give the CR3BP GTO-transfer work a consistent, legible layout: one folder per
*problem* (`GTO_tulip`, `GTO_ELFO`), each split into `direct/` (CasADi/collocation
NLP) and `indirect/` (PMP shooting) codebases, with the genuinely shared CR3BP
problem-definition extracted into one `cr3bp_common/` library. Today the direct-NLP
folder actually holds three indirect campaigns (`ms_band`, `ifs`, `ztl`) and a direct
ELFO campaign (`elfo`) mixed in with the direct tulip solver, and the shared params
are duplicated (`sundman_minfuel/` **and** `PSR/lib/`). This restructure makes method
and problem visible from the tree.

**Non-goals.** (a) `earth_elliptic_to_geo/{direct,indirect}` — later. (b) Any change
to solver *logic* or numerics — this is pure relocation + `setup_paths` rewiring +
reference repointing; every module must behave identically after. (c) `mfmax/`,
`min_fuel_paper(s)/`, and all non-GTO folders — untouched. (d) Renaming
`lowThrust_GTO_tulip` — keep the name; it becomes `GTO_tulip/indirect/lowThrust_GTO_tulip`.

---

## 2. Target structure

```
orbit_transfer/
├── cr3bp_common/                        NEW — shared CR3BP GTO-transfer library
│   ├── cr3bp_lt_params.m                (CR3BP low-thrust dynamics params)
│   ├── minfuel_config.m                 (solver config)
│   ├── gto_tulip_endpoints.m            (tulip problem def; used by tulip d/i AND elfo)
│   ├── gto_elfo_endpoints.m             (elfo problem def)
│   └── setup_cr3bp_common.m             (adds itself + pumpkyn; called by every module)
│
├── GTO_tulip/                           campaign docs stay at root: README.md,
│   │                                    LOW_THRUST_MINFUEL_CAMPAIGN.md, ROADMAP.md,
│   │                                    MIN_ENERGY_NOTES.md, HONEST_EVALUATION_*.md,
│   │                                    CODE_CLEANUP_PLAN.md, sundman_minfuel_solution_note.tex,
│   │                                    reviews/, doc/, attic/
│   ├── direct/
│   │   ├── sundman_minfuel/             (direct Sundman min-fuel engine)
│   │   ├── PSR/    (+ PSR_data/)        (PMP-steered refinement; self-contained lib/ kept)
│   │   └── movie/                       (direct-solution animations)
│   └── indirect/
│       ├── lowThrust_GTO_tulip/         (base indirect codebase: lt_pmp_eom*, shooting)
│       ├── ms_band/
│       ├── ifs/   (+ IFS_data/)
│       ├── ztl/
│       └── min_time/                    (PMP min-time root; cohesive, kept whole)
│
└── GTO_ELFO/
    ├── direct/
    │   └── elfo/                        (direct CasADi ELFO campaign)
    └── indirect/                        (placeholder + README stub; Route-C is future work)
```

Depth is **uniform**: every module lands at `orbit_transfer/<PROBLEM>/<METHOD>/<module>/`,
i.e. exactly 3 levels below `orbit_transfer/`, so `cr3bp_common` is `../../../cr3bp_common`
from every module — a single reusable relative expression.

---

## 3. Classification (as-built, verified)

| current path | method | → destination |
|---|---|---|
| `GTO_tulip/sundman_minfuel/` | direct | `GTO_tulip/direct/sundman_minfuel/` |
| `GTO_tulip/PSR/` `PSR_data/` | direct | `GTO_tulip/direct/PSR/` `direct/PSR_data/` |
| `GTO_tulip/movie/` | direct | `GTO_tulip/direct/movie/` |
| `GTO_tulip/elfo/` | direct (ELFO) | `GTO_ELFO/direct/elfo/` |
| `GTO_tulip/ms_band/` | indirect | `GTO_tulip/indirect/ms_band/` |
| `GTO_tulip/ifs/` `IFS_data/` | indirect | `GTO_tulip/indirect/ifs/` `indirect/IFS_data/` |
| `GTO_tulip/ztl/` | indirect | `GTO_tulip/indirect/ztl/` |
| `GTO_tulip/min_time/` | indirect | `GTO_tulip/indirect/min_time/` |
| `lowThrust_GTO_tulip/` (sibling) | indirect | `GTO_tulip/indirect/lowThrust_GTO_tulip/` |
| `GTO_tulip/{README,ROADMAP,*.md,doc,reviews,attic}` | docs/legacy | `GTO_tulip/` (root, unchanged names) |

`cr3bp_common/` takes the shared params from their **authoritative** source
`sundman_minfuel/`: `cr3bp_lt_params.m`, `minfuel_config.m`, `gto_tulip_endpoints.m`,
plus `gto_elfo_endpoints.m` from `elfo/`. These are the files `elfo`, `ms_band`, `ifs`,
and `min_time` reach across modules to use, so centralizing them removes the
*cross-module* duplication.

**`PSR/lib/` is left intact.** It turned out to be not a 3-file param copy but a ~20-file
**self-contained** library (`casadi_minfuel_sundman`, the full `refine` suite —
`prep_refine_seed`/`refine_loop`/`refine_sigma`/`pmp_refine_indicator` —, an `sms_*`
multiple-shooting set, and the params). Fully de-duplicating it against `sundman_minfuel`
is a separate, behavior-risky refactor and is **out of scope** here. Consequence: `PSR`
is fully independent (it uses its own `lib/`), so after moving it needs **no**
`cr3bp_common` — its `setup_paths` (self + `lib/` + pumpkyn) is unchanged. Full PSR/lib
de-dup is noted as future cleanup, not part of this plan.

---

## 4. Shared library mechanism

`cr3bp_common/setup_cr3bp_common.m` is the single place that adds the shared params +
pumpkyn:

```matlab
function setup_cr3bp_common()
% SETUP_CR3BP_COMMON  Add the shared CR3BP GTO-transfer library (params, config,
% endpoints) and the pumpkyn source to the path. Called by every GTO module's setup_paths.
here = fileparts(mfilename('fullpath'));
addpath(here);
pumpkynSrc = fullfile(getenv('HOME'),'Desktop','proj7','external','pumpkyn','src');
assert(exist(fullfile(pumpkynSrc,'+pumpkyn'),'dir')==7, ...
    'setup_cr3bp_common:missing','pumpkyn not found at %s', pumpkynSrc);
addpath(pumpkynSrc);
end
```

Every module's `setup_paths.m` is rewritten to the same shape: add self, bring in
`cr3bp_common`, then add its specific cross-ref dirs:

```matlab
here = fileparts(mfilename('fullpath'));  addpath(here);
addpath(fullfile(here,'..','..','..','cr3bp_common'));  setup_cr3bp_common();
% ... module-specific cross-refs (see §5) ...
```

This replaces every current `addpath(../sundman_minfuel)` (which was pulling in the
shared params) and removes the pumpkyn boilerplate from each file (now central).

---

## 5. Dependency edges after the move (the rewiring)

Verified from code usage, not guessed. Each row is what that module's `setup_paths`
must add **besides** `self + cr3bp_common`:

| module (new path) | extra dirs it must add | why (verified) |
|---|---|---|
| `GTO_tulip/direct/sundman_minfuel` | — | self + cr3bp_common only |
| `GTO_tulip/direct/PSR` | `lib/` (its own) | **self-contained; setup_paths unchanged, NO cr3bp_common** |
| `GTO_tulip/direct/movie` | `../sundman_minfuel` | animates its solutions |
| `GTO_ELFO/direct/elfo` | `../../../GTO_tulip/direct/sundman_minfuel` | **uses `casadi_minfuel_sundman`, `insertion_states`, `minfuel_at_tf`** |
| `GTO_tulip/indirect/lowThrust_GTO_tulip` | — | self + cr3bp_common (+ pumpkyn) |
| `GTO_tulip/indirect/ms_band` | `../lowThrust_GTO_tulip`, `../../direct/sundman_minfuel` | `lt_pmp_eom*`; reads dual `.mat`s from `sundman_minfuel/results` |
| `GTO_tulip/indirect/ifs` | `../lowThrust_GTO_tulip`, `../ms_band`, `../../direct/sundman_minfuel`, `../../direct/sundman_minfuel/refine` | indirect + reuses the direct refine engine |
| `GTO_tulip/indirect/ztl` | `../lowThrust_GTO_tulip` | `lt_pmp_eom*` |
| `GTO_tulip/indirect/min_time` | — | self + cr3bp_common (uses `gto_tulip_endpoints`) |

### Cross-references we deliberately keep (honest caveat)
"Clean direct/indirect" here means **clean folders + explicit `setup_paths` edges**,
not zero cross-refs. The code genuinely has these edges; duplicating solver code to
erase them would be worse:
- `GTO_ELFO/direct/elfo` → `GTO_tulip/direct/sundman_minfuel` (ELFO is built on the tulip Sundman engine, retargeted).
- `GTO_tulip/indirect/{ms_band,ifs}` → `GTO_tulip/direct/sundman_minfuel` (dual-`.mat` data + the `refine` engine).
Both are documented in the respective `setup_paths` headers.

---

## 6. External references to repoint

- **`earth_elliptic_to_geo`**: references `GTO_tulip/<subfolder>/...` in **docs/comments only**
  (README, `process/*.md`, a few `.m` header comments) — verified **no executable
  dependency**. Repoint the paths (e.g. `GTO_tulip/PSR/` → `GTO_tulip/direct/PSR/`,
  `GTO_tulip/sundman_minfuel/` → `GTO_tulip/direct/sundman_minfuel/`,
  `GTO_tulip/ms_band/` → `GTO_tulip/indirect/ms_band/`, `GTO_tulip/elfo/` →
  `GTO_ELFO/direct/elfo/`). Non-breaking; hygiene.
- **Hub + project `CLAUDE.md`**: update the Optimal Control map to show
  `cr3bp_common/`, `GTO_tulip/{direct,indirect}/`, `GTO_ELFO/{direct,indirect}/`.
- **`docs/superpowers/` historical plan/spec docs**: leave as point-in-time records
  (they already predate this and reference old paths in prose) — consistent with the
  earlier reorg's treatment of historical docs.
- Executable absolute/`cd`/`.sh` refs **inside the moved modules** (e.g.
  `movie/gen_movie_data.m`, the `PSR`/`elfo` `.sh` how-to-run headers) get repointed to
  the new paths — same class the last reorg handled.

---

## 7. Verification strategy

Path-resolution is the test (no solves needed). After each module moves and its
`setup_paths` is rewritten, run a one-line MATLAB smoke that `cd`s into the module,
runs `setup_paths`, and asserts the functions it actually calls resolve. Concretely:
- shared: `exist('cr3bp_lt_params','file')==2 && exist('minfuel_config','file')==2`
- `elfo`: additionally `exist('casadi_minfuel_sundman','file')==2 && exist('insertion_states','file')==2`
- indirect trio: additionally `exist('lt_pmp_eom_minfuel','file')==2`
- `ifs`: additionally the `sundman_minfuel/refine` function it calls resolves

A final **global dep sweep** (`git grep` for any remaining `../sundman_minfuel`,
`../../lowThrust_GTO_tulip`, or old `GTO_tulip/<subfolder>` paths in tracked files)
catches stragglers, exactly as the last reorg's sweep did. Git history preserved via
`git mv`; verify with `git log --follow` on a sample per moved module.

---

## 8. Open items to confirm during implementation

1. ~~PSR/lib de-dup~~ **RESOLVED:** `PSR/lib/` is a ~20-file self-contained library;
   left intact, PSR needs no `cr3bp_common`. Full de-dup is deferred future cleanup.
2. **`ifs` → `sundman_minfuel/refine`** **RESOLVED:** `ifs` calls `prep_refine_seed`;
   the resolve-smoke asserts `exist('prep_refine_seed','file')==2`.
3. **`.mat` data paths** that `ms_band`/`ifs` read from `sundman_minfuel/results` — these
   are gitignored caches; only the *path* matters and is handled by the `setup_paths`
   edge. Confirm no hardcoded absolute `.mat` path in tracked code (the sweep covers it).
4. **`GTO_ELFO/indirect/` placeholder** — create with a short `README.md` noting the
   ELFO indirect (Route-C) work is future; empty dirs aren't tracked by git, so the
   README is what makes it exist in the repo.

---

## 9. Risks

- **Wide `setup_paths` rewrite** (~10 files) — mitigated by the per-module resolve
  smoke gate before each commit.
- **Cross-refs create ordering** (ELFO needs GTO_tulip/direct/sundman_minfuel present) —
  the plan sequences GTO_tulip before GTO_ELFO so the target exists.
- **`min_time` internal cohesion** — kept whole, so its intra-module refs
  (`mintime_ms_residual/seed`) are unaffected.
- **Reversibility** — single repo, all `git mv`; unwind before push if a smoke fails.
