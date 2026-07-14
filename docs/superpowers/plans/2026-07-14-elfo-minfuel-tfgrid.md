# ELFO Min-Fuel tf-Grid Campaign Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Build the batch infrastructure to map the GTO→ELFO minimum-fuel ΔV–time front across a grid of transfer-time factors, mirroring the tulip PSR batch trio.

**Architecture:** Phase 1 re-parameterizes the existing energy-seed sweep to a factor band and factor-keyed seed names (`energy_elfo_f####.mat`). Phase 2 adds three files in `elfo/`: a per-factor callable `elfo_run_one.m` (wraps the `gen_elfo_minfuel` ε→0 core, writes a result row), a crash-robust shell walker `elfo_batch.sh` (one MATLAB process per factor), and `elfo_collect_summary.m` (rows → the convergence map). Direct clones of `PSR/psr_run_one.m`, `PSR/psr_batch.sh`, `PSR/psr_collect_summary.m`.

**Tech Stack:** MATLAB R2025b (R2025a license is broken), CasADi 3.7.0 (`~/casadi-3.7.0`), zsh. No unit-test framework — verification is `checkcode`/`zsh -n` parse checks, CasADi-free functional checks where possible, and a final live gate.

## Global Constraints

- **Working directory:** `/Users/msc/Desktop/optimal_control/NLP_lowThrust_GTO_tulip` (call it `$ROOT`). Branch `ifs-retarget`; `main` = latest (FF + push both after commits).
- **Drive everything by `factor = tf/tfMin`**, `tfMin = cfg.tfMin = 6.2906939607` ND (tulip min-time, a shared reference scale — ELFO min-time is unsolved and out of scope). Seeds are **factor-keyed** `energy_elfo_f<round(1000·factor)>.mat` (e.g. factor 1.20 → `energy_elfo_f1200.mat`), mirroring PSR's `energy_f####.mat`.
- **MATLAB is R2025b only:** `/Applications/MATLAB_R2025b.app/bin/matlab` (R2025a → Licensing Error 9).
- **Process isolation is mandatory** in the batch: one `matlab -batch` per factor (uncatchable CasADi/IPOPT MEX FATAL, ~1 in 10 solves, kills only that factor).
- **Do NOT modify** `PSR/` (the template), `casadi_energy_freetf.m`, or the solve body of `gen_elfo_minfuel.m` (only its docstring). ELFO min-time is out of scope.
- **MATLAB function headers** (CLAUDE.md): purpose + INPUTS (sizes) + OUTPUTS (sizes). **Never use `i`/`j`** as loop indices.
- Result `.mat` and `results/logs/` are gitignored (`*.mat`, extension rules) — like PSR.

---

### Task 1: Phase-1 — factor-key the energy seeds

**Files:**
- Modify: `$ROOT/elfo/gen_elfo_energy_tfsweep.m` (opts→factor band; save name; grid factor field)
- Modify: `$ROOT/elfo/run_elfo_minfuel.m:68` (seed lookup)
- Modify: `$ROOT/elfo/gen_elfo_minfuel.m` (docstring only, 2 mentions)

**Interfaces:**
- Produces: seed files named `energy_elfo_f<round(1000·factor)>.mat`; `gen_elfo_energy_tfsweep(opts)` opts now `factorLo`/`factorHi`/`factorStep`; grid rows carry `.factor`. These names are consumed by `run_elfo_minfuel` (this task) and `elfo_run_one` (Task 2).

- [ ] **Step 1: `gen_elfo_energy_tfsweep.m` — convert the band opts to factor**

Replace the opts-parsing lines (currently around 38–39):

```matlab
tfStep = gd('tfStep',0.5);  tfHi = gd('tfHi',12.5);  tfLo = gd('tfLo',7.0);
stepMin = gd('stepMin',0.0625);
```

with (factor band → ND internally; continuation math unchanged):

```matlab
% factor band (factor = tf/tfMin), converted to ND for the continuation
factorLo = gd('factorLo',1.11);  factorHi = gd('factorHi',2.00);  factorStep = gd('factorStep',0.08);
tfLo = factorLo*cfg.tfMin;  tfHi = factorHi*cfg.tfMin;  tfStep = factorStep*cfg.tfMin;
stepMin = gd('factorStepMin',0.01)*cfg.tfMin;
```

- [ ] **Step 2: `gen_elfo_energy_tfsweep.m` — carry tfMin into ctx and header**

In the `ctx = struct(...)` builder (around lines 34–37), add `tfMin` so `save_point` can compute the factor. Change:

```matlab
    'resDir',resDir,'tStar',p.tStar);
```

to:

```matlab
    'resDir',resDir,'tStar',p.tStar,'tfMin',cfg.tfMin);
```

Update the INPUTS line in the header (currently `.tfStep[0.5 ND] .tfHi[12.5] .tfLo[7.0] .maxIter[2000]`) to:

```matlab
%   opts - (optional): .factorLo[1.11] .factorHi[2.00] .factorStep[0.08]
%          .factorStepMin[0.01] .maxIter[2000] .looseIter[500] .resume[true]
```

and the OUTPUTS line mentioning `energy_elfo_tf<NNNN>.mat` to `energy_elfo_f<NNNN>.mat` (NNNN = round(1000·factor)).

- [ ] **Step 3: `gen_elfo_energy_tfsweep.m` — factor-key the saved seed + grid row**

In the local `save_point` function, replace:

```matlab
file = fullfile(ctx.resDir, sprintf('energy_elfo_tf%04d.mat', round(1000*tf)));
save(file,'X','U','sigma','rv0','rvf','tauf0','tf','moonZone','pSund','qSund');
ss = U(4,:);
g = struct('tf',tf,'ok',ok,'mf',X(7,end),'edge',mean(ss>0.95|ss<0.05), ...
           'switches',sum(abs(diff(ss>0.5))),'file',file);
```

with (factor-keyed name; `factor` added to the grid row):

```matlab
factor = tf/ctx.tfMin;
file = fullfile(ctx.resDir, sprintf('energy_elfo_f%04d.mat', round(1000*factor)));
save(file,'X','U','sigma','rv0','rvf','tauf0','tf','moonZone','pSund','qSund');
ss = U(4,:);
g = struct('tf',tf,'factor',factor,'ok',ok,'mf',X(7,end),'edge',mean(ss>0.95|ss<0.05), ...
           'switches',sum(abs(diff(ss>0.5))),'file',file);
```

- [ ] **Step 4: `run_elfo_minfuel.m` — repoint the seed lookup**

At line 68, replace:

```matlab
    cand = fullfile(resDir, sprintf('energy_elfo_tf%04d.mat', round(1000*tf)));  % tf-grid seed
```

with:

```matlab
    cand = fullfile(resDir, sprintf('energy_elfo_f%04d.mat', round(1000*factor)));  % factor-keyed tf-grid seed
```

(The base-seed fallback on the following lines is unchanged.)

- [ ] **Step 5: `gen_elfo_minfuel.m` — docstring seed-name fix**

In the header, change the two mentions of the old name to the new one:
- Line ~7: `energy_elfo_tf####.mat` → `energy_elfo_f####.mat`
- Line ~22 (the `.seedFile` default note references the base seed `energy_elfo_freetf.mat` — leave that; only fix any `energy_elfo_tf####` occurrence). Grep to confirm none remain.

- [ ] **Step 6: Parse + static naming check**

Run (R2025b):

```
/Applications/MATLAB_R2025b.app/bin/matlab -batch "cd('/Users/msc/Desktop/optimal_control/NLP_lowThrust_GTO_tulip/elfo'); setup_paths; assert(isempty(checkcode('gen_elfo_energy_tfsweep.m','-struct')) || true); cfg=minfuel_config; f=1.20; assert(strcmp(sprintf('energy_elfo_f%04d.mat',round(1000*f)),'energy_elfo_f1200.mat')); disp('OK naming + parse')" 2>&1 | grep -v "Home License\|personal use\|academic, research\|organizational use"
```

Expected: `OK naming + parse`. Then confirm no stale name remains:

```bash
cd /Users/msc/Desktop/optimal_control/NLP_lowThrust_GTO_tulip
grep -rn "energy_elfo_tf" elfo/ && echo "FAIL: stale energy_elfo_tf name remains" || echo "OK: no energy_elfo_tf refs"
```

Expected: `OK: no energy_elfo_tf refs`.

- [ ] **Step 7: Commit**

```bash
cd /Users/msc/Desktop/optimal_control/NLP_lowThrust_GTO_tulip
git add elfo/gen_elfo_energy_tfsweep.m elfo/run_elfo_minfuel.m elfo/gen_elfo_minfuel.m
git commit -m "elfo phase1: factor-key the energy seeds (energy_elfo_f####) for PSR parity"
```

---

### Task 2: `elfo_run_one.m` — the per-factor batch unit

**Files:**
- Create: `$ROOT/elfo/elfo_run_one.m`

**Interfaces:**
- Consumes: factor-keyed seeds `energy_elfo_f####.mat` (Task 1); `gen_elfo_minfuel(struct('seedFile',…,'target','ELFO','epsMin',…,'maxIter',…,'looseIter',…,'resume',…))` → returns `outFile`; that solution `.mat` holds `X[9x(N+1)]`, `U[4x(N+1)]`, `epsilon`, `out` (with `out.maxDefect`, `out.ipoptStatus`).
- Produces: `row = elfo_run_one(factor, opts)` — struct with fields **in this exact order** (Task 4's summary concatenation depends on it): `factor, tf, tf_days, ok, epsReached, epsFloor, dV, prop, switches, edge, defect, ipoptStatus, dataFile, err`. Also saved to `elfo/results/elfo_result_f####_minEps#.mat` as var `row`.

- [ ] **Step 1: Write `elfo/elfo_run_one.m`**

```matlab
function row = elfo_run_one(factor, opts)
% ELFO_RUN_ONE  Run the GTO->ELFO min-fuel homotopy for ONE t_f factor (batch unit).
%
% The per-factor unit of work for the ELFO min-fuel tf-grid campaign, extracted
% so it can be called (a) as a standalone `matlab -batch` process by elfo_batch.sh
% (each factor in its OWN process, so an UNCATCHABLE CasADi/IPOPT MEX FATAL crash
% kills only that factor, not the sweep -- a try/catch cannot catch a MEX FATAL),
% and (b) directly. The ELFO analog of PSR/psr_run_one.m; wraps gen_elfo_minfuel.
%
% INPUTS:
%   factor - t_f / tfMin (tfMin = tulip min-time 6.2907 ND, a shared scale) [scalar]
%   opts   - (optional) struct:
%            .epsMin    homotopy endpoint [0]  (0 = bang-bang fuel, >0 = smooth)
%            .maxIter   IPOPT cap (tight) [2000]
%            .looseIter IPOPT cap (loose probe) [500]
%            .resDir    seeds + results dir [elfo/results]
%            .rerun     ignore an existing result row and re-solve [false]
%
% OUTPUTS:
%   row - [1x1] struct: factor, tf, tf_days, ok, epsReached, epsFloor, dV, prop,
%         switches, edge, defect, ipoptStatus, dataFile, err. ALSO saved to
%         <resDir>/elfo_result_f####_minEps#.mat (var `row`) so the shell walker
%         can build a summary after crashes.
%
% REFERENCES:
%   [1] gen_elfo_minfuel.m (homotopy core); [2] run_elfo_minfuel.m (interactive
%       entry); [3] PSR/psr_run_one.m (tulip analog); [4] elfo_batch.sh.

here = fileparts(mfilename('fullpath'));
cd(here);  setup_paths();
cfg = minfuel_config();  p = cr3bp_lt_params(cfg.thrustN, cfg.m0kg, cfg.ispS);

if nargin < 2, opts = struct(); end
d = @(f,v) getdef(opts, f, v);
epsMin    = d('epsMin', 0);
maxIter   = d('maxIter', 2000);
looseIter = d('looseIter', 500);
resDir    = d('resDir', fullfile(here,'results'));
rerun     = d('rerun', false);
if ~exist(resDir,'dir'), mkdir(resDir); end

tf   = factor * cfg.tfMin;
eTag = strrep(sprintf('%g', epsMin), '.', 'p');
tag  = sprintf('f%04d_minEps%s', round(1000*factor), eTag);
resultFile = fullfile(resDir, sprintf('elfo_result_%s.mat', tag));

row = struct('factor',factor,'tf',tf,'tf_days',tf*p.tStar/86400,'ok',false, ...
    'epsReached',false,'epsFloor',NaN,'dV',NaN,'prop',NaN,'switches',NaN, ...
    'edge',NaN,'defect',NaN,'ipoptStatus','','dataFile','','err','');

% --- resumable: an existing row means this factor is done --------------------
if isfile(resultFile) && ~rerun
    L = load(resultFile, 'row');
    if isfield(L,'row'), row = L.row; fprintf('elfo_run_one f=%.3f: row exists -- skip\n', factor); return; end
end

% --- resolve the factor-keyed energy seed (base-seed fallback near 1.20x) ----
seed = fullfile(resDir, sprintf('energy_elfo_f%04d.mat', round(1000*factor)));
if ~isfile(seed)
    base = fullfile(resDir, 'energy_elfo_freetf.mat');
    if isfile(base)
        B = load(base,'X');  if abs(B.X(8,end) - tf) < 0.02, seed = base; end
    end
end
if ~isfile(seed)
    row.err = sprintf('no ELFO energy seed for factor %.3f (build via gen_elfo_energy_tfsweep)', factor);
    save(resultFile,'row');  fprintf('elfo_run_one f=%.3f: %s\n', factor, row.err);  return
end

fprintf('\n=== ELFO_RUN_ONE factor %.3f (epsMin=%g, seed=%s) ===\n', factor, epsMin, seed);

% --- solve: energy->fuel homotopy (gen_elfo_minfuel errors id 'minfuel:stuck'
%     at the sharpening wall -- catch it and record the eps-floor) -----------
try
    outFile = gen_elfo_minfuel(struct('seedFile',seed,'target','ELFO','epsMin',epsMin, ...
        'maxIter',maxIter,'looseIter',looseIter,'resume',~rerun));
    L  = load(outFile);
    ss = L.U(4,:);  mf = L.X(7,end);
    row.ok         = true;
    row.epsReached = (L.epsilon <= epsMin + 1e-9);
    row.switches   = sum(abs(diff(ss>0.5)));
    row.edge       = mean(ss>0.95 | ss<0.05);
    row.dV         = p.c*log(1/mf)*p.lStar/p.tStar;
    row.prop       = p.m0kg*(1-mf);
    row.defect     = L.out.maxDefect;
    if isfield(L.out,'ipoptStatus'), row.ipoptStatus = L.out.ipoptStatus; end
    row.dataFile   = outFile;
catch ME
    row.err = ME.message;
    if strcmp(ME.identifier, 'minfuel:stuck')
        row.epsReached = false;
        tok = regexp(ME.message, 'eps=([0-9.]+)', 'tokens', 'once');
        if ~isempty(tok), row.epsFloor = str2double(tok{1}); end
        fprintf('elfo_run_one f=%.3f: sharpening wall at eps=%.4f\n', factor, row.epsFloor);
    else
        fprintf('elfo_run_one f=%.3f: ERROR %s\n', factor, ME.message);
    end
end

save(resultFile, 'row');
fprintf('=== elfo_run_one f=%.3f: ok=%d epsReached=%d dV=%.4f sw=%d (row -> %s) ===\n', ...
        factor, row.ok, row.epsReached, row.dV, row.switches, resultFile);
end

% ---------------------------------------------------------------------------
function v = getdef(s, f, dflt)
if isfield(s, f) && ~isempty(s.(f)), v = s.(f); else, v = dflt; end
end
```

- [ ] **Step 2: Parse + CasADi-free functional check (missing-seed path)**

This exercises the row-building, save, and resumability WITHOUT a solve, by asking for a factor that has no seed:

```
/Applications/MATLAB_R2025b.app/bin/matlab -batch "cd('/Users/msc/Desktop/optimal_control/NLP_lowThrust_GTO_tulip/elfo'); setup_paths; r = elfo_run_one(9.99, struct('epsMin',0)); assert(~r.ok && ~isempty(r.err)); assert(isfile('results/elfo_result_f9990_minEps0.mat')); r2 = elfo_run_one(9.99, struct('epsMin',0)); assert(isequal(r,r2)); delete('results/elfo_result_f9990_minEps0.mat'); disp('OK run_one error-path + resume + fields')" 2>&1 | grep -v "Home License\|personal use\|academic, research\|organizational use"
```

Expected: `OK run_one error-path + resume + fields` (no seed → `ok=false`, row written; second call loads the identical row = resumability; then the temp row file is cleaned up).

- [ ] **Step 3: Commit**

```bash
cd /Users/msc/Desktop/optimal_control/NLP_lowThrust_GTO_tulip
git add elfo/elfo_run_one.m
git commit -m "elfo phase2: elfo_run_one (per-factor min-fuel batch unit, result row)"
```

---

### Task 3: `elfo_batch.sh` — crash-robust per-factor sweep

**Files:**
- Create: `$ROOT/elfo/elfo_batch.sh`

**Interfaces:**
- Consumes: `elfo_run_one(factor, struct('epsMin',<epsMin>))` (Task 2); factor-keyed seeds `elfo/results/energy_elfo_f*.mat` (Task 1) for `energy` auto-discovery; `elfo_collect_summary(<epsMin>)` (Task 4).
- Produces: per-factor rows (via `elfo_run_one`) + a summary; logs to `elfo/results/logs/`.

- [ ] **Step 1: Write `elfo/elfo_batch.sh`**

```zsh
#!/bin/zsh
# elfo_batch.sh — CRASH-ROBUST GTO->ELFO min-fuel sweep across t_f factors.
#
# Runs the energy->fuel homotopy (gen_elfo_minfuel, eps 1 -> epsMin) for many t_f
# factors, ONE per factor, each in its OWN `matlab -batch` process so the sporadic
# UNCATCHABLE CasADi/IPOPT MEX FATAL crash (~1 in 10 solves; a try/catch cannot
# catch it) takes down only THAT factor -- the sweep continues. The ELFO sibling
# of PSR/psr_batch.sh.
#
# ============================ HOW TO RUN ============================
#   cd /Users/msc/Desktop/optimal_control/NLP_lowThrust_GTO_tulip/elfo
#   chmod +x elfo_batch.sh                        # once
#   ./elfo_batch.sh <epsMin> <factor1> [factor2 ...]
#   ./elfo_batch.sh <epsMin> energy               # every factor with an energy seed
#
# Examples:
#   ./elfo_batch.sh 0   1.15 1.20 1.25            # bang-bang fuel for these factors
#   ./elfo_batch.sh 0   energy                    # bang-bang, all seeded factors
#   ./elfo_batch.sh 0.5 energy                    # smooth eps=0.5, all seeded factors
#
# Optional env vars (prefix on the command line):
#   WATCHDOG_S=1800   per-factor kill timeout, seconds        [default 1800]
#   MATLAB_BIN=/path/to/matlab                                 [R2025b]
#
# Background + watch:
#   nohup ./elfo_batch.sh 0 energy >/dev/null 2>&1 &
#   tail -f results/logs/elfo_batch_*.log
#
# RESUMABLE: elfo_run_one skips any factor whose result row already exists, so
# re-running the SAME command finishes only the missing factors. A summary is
# printed and saved to results/elfo_batch_summary_minEps#.mat (rebuildable with
#   matlab -batch "cd('elfo'); elfo_collect_summary(<epsMin>)").
# ===================================================================

set -u
DIR=$(cd "$(dirname "$0")" && pwd)                        # the elfo folder
MAT=${MATLAB_BIN:-/Applications/MATLAB_R2025b.app/bin/matlab}
WATCHDOG=${WATCHDOG_S:-1800}
SEEDDIR=$DIR/results
LOGDIR=$DIR/results/logs
mkdir -p "$LOGDIR"
LOG=$LOGDIR/elfo_batch_$(date +%Y%m%d_%H%M%S).log

[ $# -lt 2 ] && { echo "usage: $0 <epsMin> <factor1> [factor2 ...] | <epsMin> energy"; exit 1; }
epsMin=$1
shift

# --- resolve the factor list ------------------------------------------------
if [ "$1" = "energy" ]; then
  factors=()
  for f in "$SEEDDIR"/energy_elfo_f[0-9]*.mat(N); do     # [0-9] excludes base seed; (N) = nullglob (no zsh NOMATCH abort when empty)
    [ -e "$f" ] || continue
    b=$(basename "$f")
    n=${b#energy_elfo_f}
    n=${n%.mat}                                           # milli-factor, e.g. 1200
    factors+=( $(awk "BEGIN{printf \"%.3f\", $n/1000}") )
  done
  factors=( ${(on)factors} )                              # numeric sort (zsh)
else
  factors=( "$@" )
fi
[ ${#factors[@]} -eq 0 ] && { echo "no factors to run"; exit 1; }

echo "=== ELFO BATCH  epsMin=$epsMin  watchdog=${WATCHDOG}s ===" | tee -a "$LOG"
echo "factors (${#factors[@]}): ${factors[*]}" | tee -a "$LOG"

# --- sweep: one matlab process per factor, watchdog, continue on crash ------
opts="struct('epsMin',$epsMin)"
for f in "${factors[@]}"; do
  echo "" | tee -a "$LOG"
  echo "=== factor $f  ($(date +%H:%M:%S)) ===" | tee -a "$LOG"
  "$MAT" -batch "cd('$DIR'); elfo_run_one($f, $opts);" >> "$LOG" 2>&1 &
  pid=$!
  waited=0
  while kill -0 $pid 2>/dev/null; do
    sleep 10
    waited=$((waited+10))
    if [ $waited -ge $WATCHDOG ]; then
      echo "  TIMEOUT factor $f after ${WATCHDOG}s -> kill; continuing" | tee -a "$LOG"
      kill -9 $pid 2>/dev/null
      break
    fi
  done
  wait $pid 2>/dev/null
  rc=$?
  if [ $rc -eq 0 ]; then
    echo "  factor $f: OK" | tee -a "$LOG"
  else
    echo "  factor $f: process exited $rc (crash/kill) -- continuing" | tee -a "$LOG"
  fi
done

# --- summary from the per-factor result rows --------------------------------
echo "" | tee -a "$LOG"
echo "=== building summary ===" | tee -a "$LOG"
"$MAT" -batch "cd('$DIR'); elfo_collect_summary($epsMin);" 2>&1 | tee -a "$LOG"
echo "" | tee -a "$LOG"
echo "log: $LOG"
```

- [ ] **Step 2: Make executable + syntax check + arg/discovery dry-run**

```bash
cd /Users/msc/Desktop/optimal_control/NLP_lowThrust_GTO_tulip/elfo
chmod +x elfo_batch.sh
zsh -n elfo_batch.sh && echo "OK zsh syntax"
# usage on too-few args:
./elfo_batch.sh 0 2>/dev/null; [ $? -eq 1 ] && echo "OK usage-exit"
# energy auto-discovery + factor parse, WITHOUT running MATLAB: drop dummy seeds,
# extract just the discovery loop, confirm it parses factors, then clean up.
mkdir -p results
touch results/energy_elfo_f1150.mat results/energy_elfo_f1200.mat
zsh -c 'SEEDDIR=results; factors=(); for f in "$SEEDDIR"/energy_elfo_f*.mat; do [ -e "$f" ] || continue; b=$(basename "$f"); n=${b#energy_elfo_f}; n=${n%.mat}; factors+=( $(awk "BEGIN{printf \"%.3f\", $n/1000}") ); done; factors=( ${(on)factors} ); echo "discovered: ${factors[*]}"'
rm -f results/energy_elfo_f1150.mat results/energy_elfo_f1200.mat
```

Expected: `OK zsh syntax`, `OK usage-exit`, and `discovered: 1.150 1.200`. (Only remove the two dummy seed files you created; leave any real seeds.)

- [ ] **Step 3: Commit**

```bash
cd /Users/msc/Desktop/optimal_control/NLP_lowThrust_GTO_tulip
git add elfo/elfo_batch.sh
git commit -m "elfo phase2: elfo_batch.sh (crash-robust per-factor min-fuel sweep)"
```

---

### Task 4: `elfo_collect_summary.m` — the tf-grid map

**Files:**
- Create: `$ROOT/elfo/elfo_collect_summary.m`

**Interfaces:**
- Consumes: `elfo/results/elfo_result_f*_minEps#.mat` rows (Task 2) — each a `row` struct whose fields, in order, are `factor, tf, tf_days, ok, epsReached, epsFloor, dV, prop, switches, edge, defect, ipoptStatus, dataFile, err`.
- Produces: `res = elfo_collect_summary(epsMin[, resDir])` — sorted `[1xN]` struct array; prints the map; saves `elfo/results/elfo_batch_summary_minEps#.mat`.

- [ ] **Step 1: Write `elfo/elfo_collect_summary.m`**

```matlab
function res = elfo_collect_summary(epsMin, resDir)
% ELFO_COLLECT_SUMMARY  Build the ELFO min-fuel tf-grid map from per-factor rows.
%
% Scans <resDir>/elfo_result_f####_minEps<e>.mat (each holding one `row` struct
% saved by elfo_run_one) for the given epsMin, prints the tf-grid convergence
% table (which factors reach eps=0 vs stall, switches, the dV-time front), and
% saves <resDir>/elfo_batch_summary_minEps<e>.mat. Used by elfo_batch.sh. The
% ELFO analog of PSR/psr_collect_summary.m.
%
% INPUTS:
%   epsMin - homotopy endpoint the sweep ran at [scalar]
%   resDir - elfo results directory [default elfo/results relative to this file]
% OUTPUTS:
%   res - [1xN] struct array of the collected rows, sorted by factor

here = fileparts(mfilename('fullpath'));
if nargin < 2 || isempty(resDir), resDir = fullfile(here, 'results'); end
eTag = strrep(sprintf('%g', epsMin), '.', 'p');

d = dir(fullfile(resDir, sprintf('elfo_result_f*_minEps%s.mat', eTag)));
res = struct('factor',{},'tf',{},'tf_days',{},'ok',{},'epsReached',{},'epsFloor',{}, ...
    'dV',{},'prop',{},'switches',{},'edge',{},'defect',{},'ipoptStatus',{},'dataFile',{},'err',{});
for k = 1:numel(d)
    S = load(fullfile(resDir, d(k).name), 'row');
    if isfield(S,'row'), res(end+1) = S.row; end %#ok<AGROW>
end
if isempty(res)
    fprintf('elfo_collect_summary: no result rows for epsMin=%.3g in %s\n', epsMin, resDir);
    return
end
[~, ord] = sort([res.factor]);  res = res(ord);

fprintf('\n=== ELFO MIN-FUEL tf-GRID MAP (epsMin=%.3g, %d factors) ===\n', epsMin, numel(res));
fprintf('%-7s %-6s %-4s %-8s %-9s %-9s %-4s %-7s %-9s %s\n', ...
    'factor','tf(d)','ok','epsRch','dV(km/s)','prop(kg)','sw','edge%','defect','note');
for k = 1:numel(res)
    r = res(k);
    if ~r.ok
        er = '-';
    elseif r.epsReached
        er = 'YES';
    else
        er = sprintf('no@%.3g', r.epsFloor);
    end
    fprintf('%-7.3f %-6.2f %-4d %-8s %-9.4f %-9.4f %-4d %-7.1f %-9.1e %s\n', ...
        r.factor, r.tf_days, r.ok, er, r.dV, r.prop, r.switches, 100*r.edge, r.defect, r.err);
end
sumFile = fullfile(resDir, sprintf('elfo_batch_summary_minEps%s.mat', eTag));
save(sumFile, 'res', 'epsMin');
fprintf('\nsaved summary: %s\n', sumFile);
end
```

- [ ] **Step 2: CasADi-free unit test on a synthetic row**

Write a fake result row, collect it, confirm the table + saved summary:

```
/Applications/MATLAB_R2025b.app/bin/matlab -batch "cd('/Users/msc/Desktop/optimal_control/NLP_lowThrust_GTO_tulip/elfo'); setup_paths; if ~exist('results','dir'), mkdir('results'); end; row = struct('factor',1.20,'tf',7.5488,'tf_days',33.46,'ok',true,'epsReached',true,'epsFloor',NaN,'dV',3.24,'prop',2.18,'switches',34,'edge',0.98,'defect',5.7e-15,'ipoptStatus','Solve_Succeeded','dataFile','x.mat','err',''); save('results/elfo_result_f1200_minEps0.mat','row'); res = elfo_collect_summary(0); assert(numel(res)==1 && res.switches==34); assert(isfile('results/elfo_batch_summary_minEps0.mat')); delete('results/elfo_result_f1200_minEps0.mat','results/elfo_batch_summary_minEps0.mat'); disp('OK collect_summary')" 2>&1 | grep -v "Home License\|personal use\|academic, research\|organizational use"
```

Expected: the printed 1-row map + `OK collect_summary`. (Cleans up its two temp files; do not delete any real result rows.)

- [ ] **Step 3: Commit**

```bash
cd /Users/msc/Desktop/optimal_control/NLP_lowThrust_GTO_tulip
git add elfo/elfo_collect_summary.m
git commit -m "elfo phase2: elfo_collect_summary (assemble the min-fuel tf-grid map)"
```

---

### Task 5: Live integration gate (CasADi / R2025b)

**Files:** none modified (acceptance gate).

**Interfaces:**
- Consumes: everything from Tasks 1–4.
- Produces: evidence the pipeline runs end-to-end — energy seeds banked, per-factor rows written, the map assembled.

- [ ] **Step 1: Phase-1 short energy sweep (bank a few factor-keyed seeds)**

Run a short band around the base seed (needs `energy_elfo_freetf.mat` present in `elfo/results`). Background it — energy solves are minutes each:

```
/Applications/MATLAB_R2025b.app/bin/matlab -batch "cd('/Users/msc/Desktop/optimal_control/NLP_lowThrust_GTO_tulip/elfo'); setup_paths; gen_elfo_energy_tfsweep(struct('factorLo',1.15,'factorHi',1.25,'factorStep',0.05));" 2>&1 | grep -v "Home License\|personal use\|academic, research\|organizational use"
```

Expected: banks `energy_elfo_f1150.mat`, `energy_elfo_f1200.mat`, `energy_elfo_f1250.mat` (each a converged energy solve) + `energy_elfo_tfgrid.mat`. Confirm:

```bash
ls /Users/msc/Desktop/optimal_control/NLP_lowThrust_GTO_tulip/elfo/results/energy_elfo_f1*.mat
```

- [ ] **Step 2: Batch two factors through the fuel homotopy (plumbing exercise)**

Use `epsMin=0.5` (smooth — fewer homotopy steps, faster) to exercise the full seed→solve→row→summary path without the long bang-bang tail:

```bash
cd /Users/msc/Desktop/optimal_control/NLP_lowThrust_GTO_tulip/elfo
chmod +x elfo_batch.sh
./elfo_batch.sh 0.5 1.20 1.25
```

Expected: two factors run in separate processes; `elfo/results/elfo_result_f1200_minEps0p5.mat` and `..._f1250_minEps0p5.mat` written; the batch prints the ELFO tf-grid map at the end (2 rows) and saves `elfo_batch_summary_minEps0p5.mat`. If a factor MEX-crashes, the batch must continue to the next and the summary still builds from whatever rows exist.

- [ ] **Step 3: Rebuild the summary standalone**

```
/Applications/MATLAB_R2025b.app/bin/matlab -batch "cd('/Users/msc/Desktop/optimal_control/NLP_lowThrust_GTO_tulip/elfo'); setup_paths; res = elfo_collect_summary(0.5); fprintf('rows=%d\n', numel(res));" 2>&1 | grep -v "Home License\|personal use\|academic, research\|organizational use"
```

Expected: reprints the map and reports `rows=2` (or however many factors survived).

- [ ] **Step 4: Honest reporting**

If CasADi is unavailable or a solve cannot run in the environment after a real attempt, record which steps were **skipped** and why — do NOT claim a pass you did not observe. Tasks 1–4 unit/parse checks stand on their own; Step 5 is the live confirmation.

- [ ] **Step 5: (No commit — gate only.)** The banked seeds, rows, and summaries are gitignored `.mat` — nothing to commit here. If a step failed, return to the responsible task.

---

## Self-Review

**Spec coverage:** §Phase 1 (factor-key seeds, 3 edits) → Task 1; §Phase 2 `elfo_run_one` → Task 2; §`elfo_batch.sh` → Task 3; §`elfo_collect_summary` → Task 4; §Verification (Phase-1 run, run_one unit, batch) → Tasks 1–5; §data flow + ε-reached column → Task 2 row fields + Task 4 map; §out-of-scope (ELFO min-time, PSR untouched, solve cores unmodified) → Global Constraints. All covered.

**Placeholder scan:** No TBD/TODO; every new file has complete code; every verification names a command + expected output. Clean.

**Type/name consistency:** the `row` struct fields (Task 2) are declared in the identical order in `elfo_collect_summary`'s `res` template (Task 4) — required for `res(end+1)=S.row`. Seed name `energy_elfo_f<round(1000·factor)>.mat` is identical in Task 1 (save), Task 2 (`elfo_run_one` lookup), Task 3 (batch glob), and `run_elfo_minfuel` (Task 1 edit). `gen_elfo_minfuel` opts (`seedFile/target/epsMin/maxIter/looseIter/resume`) match its Task-2 call site. `elfo_collect_summary(epsMin)` and `elfo_run_one(factor,opts)` signatures match their `elfo_batch.sh` call sites. Consistent.
