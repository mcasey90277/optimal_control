# Generic First-Order Optimality (FOC) Gate Layer Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** One AD-based, transcription-native first-order optimality checker that runs on all four orbit-transfer campaigns (earth 2-body MEE, earth CR3BP, GTO→tulip, GTO→ELFO), plus a standard per-solve optimality report printed and saved by every production driver.

**Architecture:** A new shared module `orbit_transfer/verify_common/` holds a generic core (`foc_check`) that differentiates *the model that was actually solved* (CasADi `opti.f`/`opti.g`/`opti.lam_g` by AD) instead of hand-derived physical formulas. Each campaign contributes a ~15-line **manifest** (which state row is mass/time, which control rows are direction/throttle, autonomy, horizon kind) and gains a `returnModel`+constraint-registry hook in its solver (the earth solver already has one). Existing physical verifiers (`verify_pmp_mee`, `certify_minfuel_pmp`) are kept as the interpretability layer — the generic layer is a gate over them, not a replacement. The IPOPT native-inertia (δ_w) weak-local-min certificate is ported campaign-wide in the same pass (LEAD-0) since it needs the same plumbing.

**Tech Stack:** MATLAB R2025b (`/Applications/MATLAB_R2025b.app/bin/matlab -batch`), CasADi 3.7.0 at `~/casadi-3.7.0`, existing campaign code under `orbit_transfer/`.

## Global Constraints

- **Report-only burn-in (user decision 2026-07-25):** a FOC FAIL never demotes any row's `certified` status in this build. Every report prints the line `ADVISORY (report-only burn-in): does not alter certified status.` Promotion to a hard gate is a separate future decision.
- **Byte-identical default paths:** every solver edit (`returnModel`, registry, `regHistory` capture) must be additive and default-off — absent the new option, solver outputs are byte-identical. Follow the precedent of `casadi_lt_mee.m`'s `returnModel`/`creg` (verified zero-effect there).
- **Never `opti.dual()`.** All multipliers come from `opti.lam_g` indexed by recorded `opti.g` row ranges (see `orbit_transfer/README.md` Conventions and `earth_elliptic_to_geo/process/LESSONS_DUAL_EXTRACTION.md`).
- MATLAB function header convention (purpose/INPUTS/OUTPUTS/REFERENCES) per `~/Desktop/CLAUDE.md`; never use `i`/`j` as loop variables.
- Run each module from its own folder after calling that module's `setup_paths`; tests invoke `matlab -batch` with `addpath(fullfile(getenv('HOME'),'casadi-3.7.0'))`.
- `.mat` results are gitignored campaign caches; tests that need a cached row (e.g. `MEE_M2_10N.mat`) must SKIP with a clear message when it is absent, not fail.
- Weak-minimum warm-restart rule: any warm re-solve guard gates on certified quantities (status, defect ≤ 1e-8, one-sided final mass), never node-wise drift.
- Gates in existing verifiers may not be loosened.

## File Structure

```
orbit_transfer/verify_common/           # NEW shared module
├── setup_verify_common.m               # addpath self + nothing else (self-contained)
├── foc_manifest.m                      # all campaign manifests, one constructor
├── foc_dual_to_costate.m               # generic step-weighted interval→node dual map
├── foc_check.m                         # the AD core: all generic first-order checks
├── foc_report.m                        # THE standard report: fixed-format print + sidecar save
├── foc_ipopt_inertia.m                 # δ_w regularization-history interpreter (LEAD-0)
└── tests/
    ├── test_foc_dual_to_costate.m
    ├── test_foc_check_toy.m            # self-contained bang-bang toy OCP (no campaign deps)
    ├── test_foc_check_10N.m            # integration vs certified earth row (skips if cache absent)
    ├── test_foc_ipopt_inertia.m
    └── test_foc_report.m
```

Campaign-side (modify):
- `earth_elliptic_to_geo/direct/core/casadi_lt_mee.m` — add `regHistory` capture
- `earth_elliptic_to_geo/direct/verify/refresh_duals_mee.m` — thread `returnModel`
- `earth_elliptic_to_geo/direct/verify/run_foc_mee.m` — NEW wrapper
- `earth_elliptic_to_geo/direct/verify/run_verify_pmp_all.m`, `direct/frontdoor/run_gergaud.m` — wire report
- `earth_elliptic_to_geo_CR3BP/direct/refresh_duals_cr3bp.m`, `verify_cr3bp_pmp.m`, `run_cr3bp_geo.m` — wire
- `GTO_tulip/direct/sundman_minfuel/casadi_minfuel_sundman.m` — `returnModel`+registry
- `GTO_tulip/direct/sundman_minfuel/run_foc_tulip.m` — NEW wrapper
- `GTO_ELFO/direct/elfo/casadi_energy_freetf.m`, `casadi_mintime_freetf.m` — `returnModel`+registry+`regHistory`
- `GTO_ELFO/direct/elfo/run_foc_elfo.m` — NEW wrapper
- `orbit_transfer/OPTIMALITY_CERTIFICATION.md`, `orbit_transfer/README.md`, memory — Task 10

State-row facts (verified against solver headers 2026-07-25 — do not re-derive):
- earth MEE (`casadi_lt_mee`): X=[P;ex;ey;hx;hy;m;t] (7), U=[beta(3);thr] (4). mass 6, time 7. `creg` labels already present: `defect`,`betaNorm`,`thrLo`,`thrHi`.
- tulip (`casadi_minfuel_sundman`): X=[r(3);v(3);m;t] (8), U=[alpha(3);s] (4). mass 7, time 8. Signature: `(sigma, tf, rv0, rvf, Tmax, c, muStar, X0, U0, tauf0, pSund, maxIter, epsilon, warmTight, opts)`; already returns `.lamAll`, `.lamDef`, `.regHistory`-style capture exists for `psr_ipopt_certify`.
- ELFO energy (`casadi_energy_freetf`): X=[r(3);v(3);m;t;cScale] (9), U=[alpha(3);s] (4). mass 7, time 8, cScale 9. Signature `(sigma, rv0, rvf, Tmax, cEx, muStar, X0, U0, tauf0, opts)`.
- ELFO min-time (`casadi_mintime_freetf`): X as above (9), U=alpha only (3) — throttle pinned s≡1. Same signature shape.

---

### Task 1: Module scaffold + generic dual→costate map

**Files:**
- Create: `orbit_transfer/verify_common/setup_verify_common.m`
- Create: `orbit_transfer/verify_common/foc_dual_to_costate.m`
- Test: `orbit_transfer/verify_common/tests/test_foc_dual_to_costate.m`

**Interfaces:**
- Produces: `lam = foc_dual_to_costate(LamDef, sigma)` — `LamDef` [nx×N] interval defect duals, `sigma` [(N+1)×1] monotonic mesh → `lam` [nx×(N+1)] nodal costate (step-weighted adjacent-interval average, one-sided endpoints). Works for any nx (7/8/9).

- [ ] **Step 1: Write the failing test**

```matlab
% TEST_FOC_DUAL_TO_COSTATE  Unit test for the generic step-weighted map.
% Mirrors earth test_verify_pmp_mee Test 1, but at nx=9 to prove genericity.
root = fileparts(fileparts(mfilename('fullpath'))); cd(root);
setup_verify_common;
tol = 1e-12;
sigma1 = [0; 0.3; 1.0];                       % non-uniform: h1=0.3, h2=0.7
Lam1a  = (2:10).';  Lam1b = (5:13).';         % nx=9, distinct rows
lam1 = foc_dual_to_costate([Lam1a, Lam1b], sigma1);
assert(isequal(size(lam1), [9 3]), 'wrong size');
assert(max(abs(lam1(:,1) - Lam1a)) < tol, 'node 1 one-sided');
assert(max(abs(lam1(:,3) - Lam1b)) < tol, 'node 3 one-sided');
expectMid = 0.3*Lam1a + 0.7*Lam1b;            % h-weighted (h1+h2=1)
assert(max(abs(lam1(:,2) - expectMid)) < tol, 'interior step-weighted');
assert(max(abs(lam1(:,2) - 0.5*(Lam1a+Lam1b))) > 1e-3, 'must differ from plain avg');
fprintf('test_foc_dual_to_costate: ALL PASS\n');
```

- [ ] **Step 2: Run to verify it fails**

Run: `cd orbit_transfer/verify_common/tests && /Applications/MATLAB_R2025b.app/bin/matlab -batch "test_foc_dual_to_costate"`
Expected: FAIL — `setup_verify_common` / `foc_dual_to_costate` undefined.

- [ ] **Step 3: Implement**

`setup_verify_common.m`:
```matlab
function setup_verify_common()
% SETUP_VERIFY_COMMON  Put orbit_transfer/verify_common on the MATLAB path.
% Self-contained: no campaign paths, no CasADi (callers add CasADi themselves).
% OUTPUTS: none (path side effect)
addpath(fileparts(mfilename('fullpath')));
end
```

`foc_dual_to_costate.m` — port of `earth_elliptic_to_geo/direct/verify/mee_dual_to_costate.m` with the nx==7 assumption removed; keep its header derivation references ([1] DESIGN_dual_map "[CORRECTNESS]", [2] MS_BAND_CAMPAIGN) and this body:
```matlab
sigma = sigma(:).';
N  = size(LamDef, 2);
assert(numel(sigma) == N + 1, 'foc_dual_to_costate: need N+1 sigma nodes');
h  = diff(sigma);
nx = size(LamDef, 1);
lam = zeros(nx, N + 1);
lam(:, 1) = LamDef(:, 1);   lam(:, N + 1) = LamDef(:, N);
for k = 2:N
    lam(:, k) = (h(k-1)*LamDef(:, k-1) + h(k)*LamDef(:, k)) / (h(k-1) + h(k));
end
```

- [ ] **Step 4: Run test — expect ALL PASS**
- [ ] **Step 5: Commit**

```bash
git add orbit_transfer/verify_common
git commit -m "feat(foc): verify_common scaffold + generic step-weighted dual->costate map"
```

---

### Task 2: Campaign manifests

**Files:**
- Create: `orbit_transfer/verify_common/foc_manifest.m`
- Test: extend `orbit_transfer/verify_common/tests/test_foc_dual_to_costate.m`? No — create `tests/test_foc_manifest.m`

**Interfaces:**
- Produces: `man = foc_manifest(name)` for `name ∈ {'earth_mee','earth_cr3bp','tulip','elfo_fuel','elfo_mintime','toy'}`. Struct fields (exact):
  `name` [char]; `nx`,`nu` [int]; `dirRows` [1×3 or []]; `thrRow` [int or []]; `massRow`,`timeRow` [int or []]; `autonomous` [logical]; `horizonKind` ['fixedtf'|'freetf-cscale'|'none']; `massFreeAtTf` [logical]. Empty (`[]`) fields mean "skip that check."

- [ ] **Step 1: Write the failing test**

```matlab
% TEST_FOC_MANIFEST  Field-contract test for all campaign manifests.
root = fileparts(fileparts(mfilename('fullpath'))); cd(root); setup_verify_common;
names = {'earth_mee','earth_cr3bp','tulip','elfo_fuel','elfo_mintime','toy'};
need  = {'name','nx','nu','dirRows','thrRow','massRow','timeRow', ...
         'autonomous','horizonKind','massFreeAtTf'};
for q = 1:numel(names)
    m = foc_manifest(names{q});
    for f = 1:numel(need), assert(isfield(m, need{f}), '%s missing %s', names{q}, need{f}); end
end
m = foc_manifest('earth_mee');   assert(m.nx==7 && m.massRow==6 && m.timeRow==7 && m.autonomous);
m = foc_manifest('earth_cr3bp'); assert(m.nx==7 && ~m.autonomous);
m = foc_manifest('tulip');       assert(m.nx==8 && m.massRow==7 && m.timeRow==8);
m = foc_manifest('elfo_fuel');   assert(m.nx==9 && strcmp(m.horizonKind,'freetf-cscale'));
m = foc_manifest('elfo_mintime');assert(m.nu==3 && isempty(m.thrRow));
fprintf('test_foc_manifest: ALL PASS\n');
```

- [ ] **Step 2: Run — expect FAIL (undefined foc_manifest)**
- [ ] **Step 3: Implement** (full header per house style; body:)

```matlab
base = struct('dirRows',[1 2 3], 'thrRow',4, 'autonomous',true, ...
              'horizonKind','fixedtf', 'massFreeAtTf',true);
switch lower(name)
    case 'earth_mee'
        man = base; man.name='earth_mee'; man.nx=7; man.nu=4; man.massRow=6; man.timeRow=7;
    case 'earth_cr3bp'
        man = foc_manifest('earth_mee'); man.name='earth_cr3bp'; man.autonomous=false;
    case 'tulip'
        man = base; man.name='tulip'; man.nx=8; man.nu=4; man.massRow=7; man.timeRow=8;
    case 'elfo_fuel'
        man = base; man.name='elfo_fuel'; man.nx=9; man.nu=4; man.massRow=7; man.timeRow=8;
        man.horizonKind = 'freetf-cscale';
    case 'elfo_mintime'
        man = foc_manifest('elfo_fuel'); man.name='elfo_mintime'; man.nu=3; man.thrRow=[];
    case 'toy'   % 1D double-integrator bang-bang used by test_foc_check_toy
        man = struct('name','toy','nx',2,'nu',1,'dirRows',[],'thrRow',1, ...
            'massRow',[],'timeRow',[],'autonomous',true,'horizonKind','none', ...
            'massFreeAtTf',false);
    otherwise, error('foc_manifest:unknown', 'unknown campaign %s', name);
end
```

- [ ] **Step 4: Run test — ALL PASS**
- [ ] **Step 5: Commit** `feat(foc): campaign manifests`

---

### Task 3: The AD core `foc_check`

**Files:**
- Create: `orbit_transfer/verify_common/foc_check.m`
- Test: `orbit_transfer/verify_common/tests/test_foc_check_toy.m`, `tests/test_foc_check_10N.m`

**Interfaces:**
- Consumes: `foc_manifest`, `foc_dual_to_costate`.
- Produces: `rep = foc_check(out, sigma, man, opts)` where `out` needs `.X [nx×(N+1)]`, `.U [nu×(N+1)]`, `.model.opti` (solved CasADi Opti), `.model.creg` (struct array with `.label`,`.rows`; labels used: `'defect'` required; `'thrLo'`,`'thrHi'` when `man.thrRow` nonempty). `opts` (optional): `.tolStat` [1e-6], `.tolSign` [99], `.tolTrans` [1e-3], `.sdotMin` [1e-3].
  `rep` fields (exact): `.kktStatInf .sLag .dirTanMax .dirTanMed .signPct .Sd [1×N+1] .lam [nx×N+1] .lamTimeCoV .lamTimeEnd .lamMassEndRel .singularArcNodes .sdotMinRel .nSwitches .horizonNote [char] .checksRun {cellstr} .pass [logical advisory]`

- [ ] **Step 1: Write the failing toy test** — fully self-contained bang-bang OCP with a known answer:

```matlab
% TEST_FOC_CHECK_TOY  foc_check on a self-built bang-bang toy (no campaign deps).
% min int u dt, xdot=v, vdot=u, u in [0,1], x(0)=v(0)=0, x(T)=1, v(T) FREE, T=2.
% Optimal: u=1 then u=0 (single switch) -- lam_v linear hits 0 slope, S=1+lam_v.
root = fileparts(fileparts(mfilename('fullpath'))); cd(root); setup_verify_common;
addpath(fullfile(getenv('HOME'), 'casadi-3.7.0'));  import casadi.*
N = 60;  T = 2;  sg = linspace(0, 1, N+1).';  ds = diff(sg).';
opti = casadi.Opti();
X = opti.variable(2, N+1);  U = opti.variable(1, N+1);
creg = struct('label',{},'rows',{});
f = [X(2,:); U(1,:)];                              % [v; u]
r0 = size(opti.g,1)+1;
for k = 1:N
    opti.subject_to(X(:,k+1) - X(:,k) - (ds(k)*T/2)*(f(:,k)+f(:,k+1)) == 0);
end
creg(end+1) = struct('label','defect','rows',r0:size(opti.g,1));
r0 = size(opti.g,1)+1;  opti.subject_to(U(1,:).' >= 0);
creg(end+1) = struct('label','thrLo','rows',r0:size(opti.g,1));
r0 = size(opti.g,1)+1;  opti.subject_to(U(1,:).' <= 1);
creg(end+1) = struct('label','thrHi','rows',r0:size(opti.g,1));
opti.subject_to(X(:,1) == [0;0]);  opti.subject_to(X(1,end) == 1);
opti.minimize(sum((ds*T/2).*(U(1,1:N)+U(1,2:N+1))));
opti.set_initial(X, [linspace(0,1,N+1); 0.5*ones(1,N+1)]);
opti.set_initial(U, 0.5*ones(1,N+1));
opti.solver('ipopt', struct('print_time',false), struct('print_level',0));
sol = opti.solve();
out = struct('X', full(sol.value(X)), 'U', full(sol.value(U)), ...
             'model', struct('opti', opti, 'creg', creg));
rep = foc_check(out, sg, foc_manifest('toy'), struct());
assert(rep.kktStatInf < 1e-6, 'KKT stationarity: %.2e', rep.kktStatInf);
assert(rep.signPct >= 99, 'sign law: %.1f%%', rep.signPct);
assert(rep.nSwitches == 1, 'expected 1 switch, got %d', rep.nSwitches);
assert(rep.singularArcNodes == 0, 'no singular arc expected');
assert(rep.sdotMinRel > 1e-3, 'Sdot at the switch must be regular');
assert(rep.pass, 'advisory pass expected on the toy');
fprintf('test_foc_check_toy: ALL PASS\n');
```

- [ ] **Step 2: Run — expect FAIL (foc_check undefined)**
- [ ] **Step 3: Implement `foc_check.m`** (full house header; REFERENCES: LESSONS_DUAL_EXTRACTION.md §3, diag_t1_beta.m, OPTIMALITY_CERTIFICATION.md §A2):

```matlab
function rep = foc_check(out, sigma, man, opts)
import casadi.*
if nargin < 4, opts = struct(); end
gd = @(f,v) fcdef(opts, f, v);
tolStat = gd('tolStat',1e-6); tolSign = gd('tolSign',99);
tolTrans = gd('tolTrans',1e-3); sdotMin = gd('sdotMin',1e-3);

opti = out.model.opti;  creg = out.model.creg;  sol = opti.debug;
X = out.X;  U = out.U;  N1 = size(X,2);  Nseg = N1-1;
nx = man.nx;  nu = man.nu;
assert(size(X,1)==nx && size(U,1)==nu, 'foc_check: X/U shape vs manifest');
sg = sigma(:);

xall   = full(sol.value(opti.x));
lamAll = full(sol.value(opti.lam_g));
nX = nx*N1;
% Variable-layout assert (X block then U block, column-major) -- never assume:
assert(max(abs(xall(1:nX) - X(:))) < 1e-8, 'foc_check: X layout');
assert(max(abs(xall(nX+(1:nu*N1)) - U(:))) < 1e-8, 'foc_check: U layout');
uix = @(rows,k) nX + (k-1)*nu + rows;

% Manifest semantic guards (cheap, catch a wrong manifest outright):
if ~isempty(man.massRow)
    assert(all(diff(X(man.massRow,:)) <= 1e-6), 'foc_check: massRow not nonincreasing -- manifest wrong?');
end
if ~isempty(man.timeRow)
    assert(all(diff(X(man.timeRow,:)) >= -1e-6), 'foc_check: timeRow not nondecreasing -- manifest wrong?');
end

% --- (1) full-Lagrangian KKT stationarity, sign auto-resolved ----------------
Fk = Function('Fk', {opti.x, opti.lam_g}, {gradient(opti.f,opti.x), jacobian(opti.g,opti.x)});
[gfD, AD_] = Fk(xall, lamAll);
gf = full(gfD);  A = sparse(AD_);
gP = gf + A.'*lamAll;  gM = gf - A.'*lamAll;
if norm(gP,inf) <= norm(gM,inf), s = +1; gL = gP; else, s = -1; gL = gM; end
rep = struct('kktStatInf', norm(gL,inf), 'sLag', s);
checks = {'kktStat'};

rowsOf = @(lab) [creg(strcmp({creg.label},lab)).rows];
defRows = rowsOf('defect');
assert(~isempty(defRows), 'foc_check: creg must register a ''defect'' group');

% --- (2) burn mask + throttle-block switching function -----------------------
if ~isempty(man.thrRow)
    thr = U(man.thrRow,:);  burn = thr > 0.5;
    lamNB = lamAll;  lamNB([rowsOf('thrLo'), rowsOf('thrHi')]) = 0;
    gNB = gf + s*(A.'*lamNB);           % Lagrangian gradient sans throttle-bound duals
    Sd = zeros(1,N1);
    for k = 1:N1, Sd(k) = gNB(uix(man.thrRow,k)); end
    rep.Sd = Sd;
    rep.signPct = 100*mean((Sd < 0) == burn);      % S<0 <=> full thrust
    checks{end+1} = 'signLaw';
else
    burn = true(1,N1);  rep.Sd = [];  rep.signPct = NaN;   % min-time: all-burn
end

% --- (3) minimum condition, direction part: tangential dL/dbeta --------------
if ~isempty(man.dirRows)
    tanAbs = nan(1,N1);
    for k = 1:N1
        b = U(man.dirRows,k);  b = b/norm(b);
        v = gL(uix(man.dirRows,k));
        tanAbs(k) = norm(v - (v.'*b)*b);
    end
    rep.dirTanMax = max(tanAbs(burn));  rep.dirTanMed = median(tanAbs(burn));
    checks{end+1} = 'dirTangential';
else
    rep.dirTanMax = 0;  rep.dirTanMed = 0;
end

% --- (4) nodal costates (sign-resolved) --------------------------------------
LamDef = reshape(lamAll(defRows), nx, Nseg);
rep.lam = foc_dual_to_costate(s*LamDef, sg);
checks{end+1} = 'costates';

% --- (5) time-costate behavior (Hamiltonian conditions in dual form) ---------
if ~isempty(man.timeRow)
    lt = rep.lam(man.timeRow,:);
    rep.lamTimeCoV = std(lt)/max(abs(mean(lt)),1e-30);
    rep.lamTimeEnd = lt(end);
    checks{end+1} = 'lamTime';
else
    rep.lamTimeCoV = NaN;  rep.lamTimeEnd = NaN;
end

% --- (6) transversality: free final mass -> lam_m(tf)=0 ----------------------
if ~isempty(man.massRow) && man.massFreeAtTf
    lm = rep.lam(man.massRow,:);
    rep.lamMassEndRel = abs(lm(end))/max(abs(lm));
    checks{end+1} = 'transversality';
else
    rep.lamMassEndRel = NaN;
end

% --- (7) singular arcs + regular switching (Sdot != 0) -----------------------
rep.singularArcNodes = 0;  rep.sdotMinRel = NaN;  rep.nSwitches = 0;
if ~isempty(man.thrRow)
    coastS = abs(rep.Sd(~burn));
    scaleS = median(coastS(coastS>0));  if isempty(scaleS)||isnan(scaleS), scaleS = max(abs(rep.Sd)); end
    nearZ = abs(rep.Sd) < max(1e-6*scaleS, 1e-14);
    runL = 0;
    for k = 1:N1
        if nearZ(k), runL = runL+1; else, runL = 0; end
        if runL >= 3, rep.singularArcNodes = rep.singularArcNodes + 1; end
    end
    swI = find(diff(burn) ~= 0);  rep.nSwitches = numel(swI);
    if ~isempty(swI)
        sdot = abs(rep.Sd(min(swI+1,N1)) - rep.Sd(max(swI,1))) / max(scaleS,1e-30);
        rep.sdotMinRel = min(sdot);
    end
    checks = [checks, {'singularArc','sdotRegular'}];
end

% --- (8) horizon note (G4): named, honest ------------------------------------
switch man.horizonKind
    case 'fixedtf'
        rep.horizonNote = 'fixed t_f: H=const generally nonzero; constancy via lamTimeCoV';
    case 'freetf-cscale'
        rep.horizonNote = ['free t_f via cScale state: horizon condition enters the ' ...
            'cScale adjoint rows, already inside kktStatInf; value-form H(tf) check ' ...
            'reported via lamTimeEnd (informational, derivation pending)'];
    otherwise
        rep.horizonNote = 'no horizon check applicable';
end
rep.checksRun = checks;

% --- advisory verdict (REPORT-ONLY burn-in) ----------------------------------
okSign  = isempty(man.thrRow) || rep.signPct >= tolSign;
okTrans = isnan(rep.lamMassEndRel) || rep.lamMassEndRel <= tolTrans;
okSdot  = isnan(rep.sdotMinRel) || rep.sdotMinRel > sdotMin;
rep.pass = rep.kktStatInf <= tolStat && rep.dirTanMax <= tolStat && ...
           okSign && okTrans && rep.singularArcNodes == 0 && okSdot;
end

function v = fcdef(o, f, dflt)
if isfield(o,f) && ~isempty(o.(f)), v = o.(f); else, v = dflt; end
end
```

- [ ] **Step 4: Run toy test — ALL PASS**
- [ ] **Step 5: Write the 10 N integration test** (`tests/test_foc_check_10N.m`): loads `earth_elliptic_to_geo/direct/results/MEE_M2_10N.mat` via that campaign's `setup_paths` + `refresh_duals_mee` **with `returnModel` not yet available — so in this task call `casadi_lt_mee` directly** exactly as `test_dual_extraction.m` Test 2 does (copy its recovery block verbatim), then:

```matlab
rep = foc_check(outRow, saved.sigma, foc_manifest('earth_mee'), struct());
assert(rep.kktStatInf < 1e-6);        assert(rep.dirTanMax < 1e-6);
assert(rep.signPct == 100);           assert(rep.lamMassEndRel < 1e-3);
assert(rep.singularArcNodes == 0);    assert(rep.pass);
fprintf('test_foc_check_10N: ALL PASS (Sdot min rel %.2e, nSw %d)\n', rep.sdotMinRel, rep.nSwitches);
```
Skip guard at top: `if ~isfile(row), fprintf('SKIPPED -- cache absent\n'); return; end`.

- [ ] **Step 6: Run 10 N test — ALL PASS (or SKIPPED on a machine without caches)**
- [ ] **Step 7: Commit** `feat(foc): AD-based generic first-order core (foc_check) + toy and 10N tests`

---

### Task 4: The standard report `foc_report`

**Files:**
- Create: `orbit_transfer/verify_common/foc_report.m`
- Test: `orbit_transfer/verify_common/tests/test_foc_report.m`

**Interfaces:**
- Produces: `foc_report(rep, tag, resDir)` — prints THE standard block (identical format for all campaigns) and, when `resDir` nonempty, saves `<resDir>/foc_<tag>.mat` containing `rep`. This block is deliverable (b): every production driver ends with it.

- [ ] **Step 1: Failing test** — build a fake `rep` with all fields (copy the field list from Task 3's Produces), call `foc_report(rep,'unittest',tempdir)`, assert the sidecar exists and `evalc` output contains `'FIRST-ORDER OPTIMALITY REPORT'`, `'ADVISORY'`, and `'unittest'`.
- [ ] **Step 2: Run — FAIL**
- [ ] **Step 3: Implement** — fixed format:

```matlab
fprintf('\n========== FIRST-ORDER OPTIMALITY REPORT: %s ==========\n', tag);
fprintf(' KKT stationarity  ||grad_x L||_inf : %10.3e   (sign s=%+d)   %s\n', ...)
fprintf(' Min condition (direction) tan max  : %10.3e                 %s\n', ...)
fprintf(' Min condition (throttle sign law)  : %9.2f %%                %s\n', ...)
fprintf(' Transversality |lam_m(tf)| rel     : %10.3e                 %s\n', ...)
fprintf(' Time-costate CoV (H-const dual)    : %10.3e   end %+.3e\n', ...)
fprintf(' Singular-arc nodes                 : %10d                 %s\n', ...)
fprintf(' Regular switching min|Sdot| rel    : %10.3e   (%d switches) %s\n', ...)
fprintf(' Horizon: %s\n', rep.horizonNote);
if isfield(rep,'ipopt') && ~isempty(rep.ipopt)
    fprintf(' 2nd-order (IPOPT inertia, delta_w) : %s\n', rep.ipopt.verdict);
end
fprintf(' ADVISORY verdict: %s   (report-only burn-in: does not alter certified status)\n', ...)
```
Each `%s` status is `'PASS'`/`'FAIL'`/`'--'` (`'--'` for NaN/skipped). Save sidecar with `save(fullfile(resDir, sprintf('foc_%s.mat',tag)), 'rep')`.

- [ ] **Step 4: Run — PASS.  Step 5: Commit** `feat(foc): standard first-order report (print + sidecar)`

---

### Task 5: IPOPT native-inertia interpreter + earth `regHistory` capture (LEAD-0 core)

**Files:**
- Create: `orbit_transfer/verify_common/foc_ipopt_inertia.m`
- Modify: `earth_elliptic_to_geo/direct/core/casadi_lt_mee.m` (immediately after `st = opti.stats();`)
- Test: `orbit_transfer/verify_common/tests/test_foc_ipopt_inertia.m`

**Interfaces:**
- Produces: `ic = foc_ipopt_inertia(regHistory, opts)` — `regHistory` [1×nIter δ_w values or []]; `opts.tailN` [5], `.tol` [1e-8]. Returns `ic.certLocalMin` [logical], `ic.maxTailDw` [scalar], `ic.verdict` [char]. Empty history → `certLocalMin=false`, verdict `'NO-DATA: regHistory absent (solver predates capture)'`.
- Produces (solver): `out.regHistory` on `casadi_lt_mee` outputs.

- [ ] **Step 1: Failing test** — pure logic, no solve:

```matlab
ic = foc_ipopt_inertia([zeros(1,50)]);              assert(ic.certLocalMin);
ic = foc_ipopt_inertia([zeros(1,45), 1e-3*ones(1,5)]); assert(~ic.certLocalMin);
ic = foc_ipopt_inertia([1e-2*ones(1,10), zeros(1,40)]); assert(ic.certLocalMin);  % early reg OK
ic = foc_ipopt_inertia([]);                          assert(~ic.certLocalMin);
assert(contains(foc_ipopt_inertia([]).verdict, 'NO-DATA'));
```

- [ ] **Step 2: FAIL.  Step 3: Implement** (interpretation copied from `GTO_tulip/direct/PSR/psr_ipopt_certify.m` — reference it in the header): tail = last `min(tailN, end)` entries; `certLocalMin = ~isempty && max(tail) <= tol`; verdict strings mirror psr_ipopt_certify's two phrasings.
- [ ] **Step 4: PASS.**
- [ ] **Step 5: Solver capture** in `casadi_lt_mee.m` right after `st = opti.stats();`:

```matlab
% delta_w regularization history (IPOPT native inertia; verify_common
% foc_ipopt_inertia interprets it -- LEAD-0 port of psr_ipopt_certify):
try, out.regHistory = st.iterations.regularization_size(:).'; catch, out.regHistory = []; end
```
(inside the `out = struct(...)` assembly add nothing else — append the field after the struct is built, keeping the existing struct literal untouched).

- [ ] **Step 6: Regression** — run existing `test_mee_solver_smoke` + `test_dual_extraction`; both must still pass, and assert `isfield(out,'regHistory')` in the smoke test's fixedtf branch (add one line there).
- [ ] **Step 7: Commit** `feat(foc): IPOPT delta_w inertia interpreter + regHistory capture in casadi_lt_mee (LEAD-0)`

---

### Task 6: Earth 2-body wiring + ladder run

**Files:**
- Modify: `earth_elliptic_to_geo/direct/verify/refresh_duals_mee.m` (thread `returnModel`)
- Create: `earth_elliptic_to_geo/direct/verify/run_foc_mee.m`
- Modify: `earth_elliptic_to_geo/direct/verify/run_verify_pmp_all.m` (append FOC block per row)
- Modify: `earth_elliptic_to_geo/direct/frontdoor/run_gergaud.m` (standard report after the solve/auto paths)

**Interfaces:**
- Consumes: `foc_check`, `foc_manifest('earth_mee')`, `foc_report`, `foc_ipopt_inertia`, `refresh_duals_mee`.
- Produces: `[rep, ver, infoq] = run_foc_mee(matPath)` — refreshes duals (returnModel on), runs generic FOC + the physical `verify_pmp_mee` side by side, prints the standard report, saves sidecar to `results/`.

- [ ] **Step 1:** `refresh_duals_mee.m` — add to `sopts`: `'returnModel', gd('returnModel', false)` where `gd` is the existing `d` helper; document in header. (The model rides on `out`; no signature change.)
- [ ] **Step 2:** `run_foc_mee.m`:

```matlab
function [rep, ver, infoq] = run_foc_mee(matPath)
% RUN_FOC_MEE  Standard first-order optimality report for a certified MEE row.
% (header per house style; REFERENCES: verify_common/foc_check.m, foc_report.m,
%  refresh_duals_mee.m, OPTIMALITY_CERTIFICATION.md Part A)
run(fullfile(fileparts(fileparts(fileparts(mfilename('fullpath')))), '..', 'verify_common', 'setup_verify_common.m'));
[outq, parq, sigq, infoq] = refresh_duals_mee(matPath, struct('returnModel', true));
rep = foc_check(outq, sigq, foc_manifest('earth_mee'), struct());
rep.ipopt = foc_ipopt_inertia(getfield_default(outq, 'regHistory', []));
ver = verify_pmp_mee(outq, parq, sigq, struct('eps', 0));   % physical layer alongside
[~, tag] = fileparts(matPath);
foc_report(rep, tag, fullfile(module_root(), 'results'));
end
```
(`getfield_default` = 3-line local subfunction, isfield-guarded.) Note the relative `setup_verify_common` path: `verify/` → module root → `earth_elliptic_to_geo` → `orbit_transfer/verify_common`. Compute it with `fullfile(fileparts(module_root()), 'verify_common')` instead — `module_root()` is already on the path; use that form.
- [ ] **Step 3:** `run_verify_pmp_all.m` — inside the row loop, after the `verify_pmp_mee` call, add: rebuild refresh with returnModel (change the existing `refresh_duals_mee(f)` call to `refresh_duals_mee(f, struct('returnModel',true))`), then `repq = foc_check(outq, sigq, foc_manifest('earth_mee'), struct()); repq.ipopt = foc_ipopt_inertia(...); foc_report(repq, base, resDir);` and add `focPass` to the summary table.
- [ ] **Step 4:** `run_gergaud.m` — at the end of the auto/solve result path (after its Table-3 row printout; locate the `fprintf` summary block), when the solve produced `out` with `.model` absent, re-request: simplest compliant wiring is `try`-guarded: if `res.fuel` exists and the run performed a live solve this session, call `run_foc_mee` on the saved row file; on probe/cached paths print `FOC report: run run_foc_mee('<row.mat>') for the standard report` instead. Keep it advisory and non-fatal (`try/catch` with a one-line warning).
- [ ] **Step 5: Run** `run_foc_mee('results/MEE_M2_10N.mat')` — expect standard report, all PASS, ipopt verdict line present.
- [ ] **Step 6: Run the ladder** via `run_verify_pmp_all` (background, ~30 min): expect 9/9 `focPass=1`, sidecars `foc_MEE_M2_*.mat` in `results/`.
- [ ] **Step 7: Commit** `feat(foc): earth 2-body wiring -- run_foc_mee + ladder sweep + front-door report`

---

### Task 7: Earth CR3BP wiring + full-ladder coverage

**Files:**
- Modify: `earth_elliptic_to_geo_CR3BP/direct/refresh_duals_cr3bp.m` (`'returnModel', d('returnModel', false)` into `sopts`)
- Modify: `earth_elliptic_to_geo_CR3BP/direct/verify_cr3bp_pmp.m` (FOC block after the dual refresh)
- Modify: `earth_elliptic_to_geo_CR3BP/direct/run_cr3bp_geo.m` (print pointer line to `verify_cr3bp_pmp` at end — the live-solve report belongs to the verifier here since the front door already has its own gates)

**Interfaces:** consumes Task 3/4/5 products + `foc_manifest('earth_cr3bp')`.

- [ ] **Step 1:** In `verify_cr3bp_pmp.m`, inside the `refreshDuals` branch, change the `refresh_duals_cr3bp(S, par, ...)` call to pass `returnModel=true` in its opts, then after the existing `verify_pmp_mee` call add:

```matlab
setup_verify_common_path();   % local helper: addpath(fullfile(fileparts(fileparts(module root)), 'verify_common'))
rep = foc_check(out, sigma, foc_manifest('earth_cr3bp'), struct());
rep.ipopt = foc_ipopt_inertia(getfield_default(out, 'regHistory', []));
foc_report(rep, sprintf('cr3bp_T%sN_phi%s', num_tag(fp.thrustN), num_tag(fp.phi0)), resDir);
ver.foc = rep;
```
Note `earth_cr3bp` sets `autonomous=false`, so `foc_report` prints the time-costate line without a PASS/FAIL (the `'--'` status) — this is G6 handled honestly: constancy is *not* gated, but the full non-autonomous adjoint recursion **is** inside `kktStatInf`.
- [ ] **Step 2: Run** on 10 N: `verify_cr3bp_pmp(struct('thrustN',10,'phi0',0))` — expect primer PASS (known) + FOC report PASS.
- [ ] **Step 3: Run the full ladder** — loop `thrustN ∈ {10,5,2.5,1,0.5,0.2,0.1}` (background; 0.2/0.1 N are large — allow long wall): closes the CR3BP coverage gap (G2). Record per-rung pass/fail.
- [ ] **Step 4: Commit** `feat(foc): CR3BP wiring + first full-ladder first-order coverage`

---

### Task 8: Tulip — solver hook + wrapper + independence cross-check

**Files:**
- Modify: `GTO_tulip/direct/sundman_minfuel/casadi_minfuel_sundman.m`
- Create: `GTO_tulip/direct/sundman_minfuel/run_foc_tulip.m`
- Test: reuse `verify_common/tests` pattern — add `GTO_tulip/direct/sundman_minfuel/test_foc_tulip_smoke.m`

**Interfaces:**
- Produces (solver): `opts.returnModel` (default false) → `out.model = struct('opti',opti,'creg',creg)` with creg labels `'defect'`, `'betaNorm'`, `'thrLo'`, `'thrHi'` (rows recorded at build via `r0 = size(opti.g,1)+1; ... creg(end+1)=struct('label',lab,'rows',r0:size(opti.g,1));` around the existing defect loop, the `||alpha||=1` constraints, and the throttle box constraints — find each `opti.subject_to` site and bracket it; **default-off, byte-identical**: registry recording only reads `size(opti.g,1)`).
- Produces: `rep = run_foc_tulip(matFile)` — loads a tulip artifact (expected fields: `out.X [8×nN]`, `out.U [4×nN]`, plus mesh `sigma` and `tauf0`, endpoint/config data; consult `GTO_tulip/direct/sundman_minfuel/README.md` for the certified flagship artifact's exact layout and the `minfuel_config`/`gto_tulip_endpoints` accessors in `cr3bp_common/`), warm re-solves at the saved primal with `epsilon=0, warmTight=true, returnModel=true`, guards on certified quantities (status + defect ≤ 1e-8 + one-sided mf), runs `foc_check` with `foc_manifest('tulip')` + `foc_ipopt_inertia(out.regHistory)` + `foc_report`.

- [ ] **Step 1:** Add registry + `returnModel` to the solver (pattern above; opts already exists as arg 15).
- [ ] **Step 2: Byte-identity regression** — `test_foc_tulip_smoke.m`: run the solver twice on a tiny warm case (reuse the smallest existing smoke input in that folder — see `smoke_seed.m` in `../..//indirect/ms_band/` is indirect; for direct use the flagship artifact's own warm re-solve at `maxIter=5`), once without `returnModel` and once with, assert `isequal(o1.X,o2.X) && isequal(o1.U,o2.U)`; assert `o2.model.creg` labels present.
- [ ] **Step 3:** Implement `run_foc_tulip.m` per the interface; print the standard report with tag from the file name.
- [ ] **Step 4: Run on the certified flagship row** (the 25-switch 1.15× artifact — path per that folder's README). Expect: KKT PASS. **Cross-check of the two costate sources (G1 partial):** also run the existing `certify_minfuel_pmp` on the same row and print both primer numbers side by side — the LS-reconstructed costate (certify) vs the raw-dual costate (foc). Record agreement/disagreement in the report sidecar as `rep.crossCheck`.
- [ ] **Step 5:** Wire the standard report into `run_psr.m` stage 5b (immediately after the existing `psr_ipopt_certify` block — one `foc_report` call using the final solve's `out`, which now carries `.model` when `run_psr` passes `returnModel=true` in its final `casadi_minfuel_sundman` call; make that final-solve-only, keeping all homotopy-sweep calls unchanged).
- [ ] **Step 6: Commit** `feat(foc): tulip solver model hook + run_foc_tulip + LS-vs-dual costate cross-check`

---

### Task 9: ELFO — the campaign with NO gate gets one

**Files:**
- Modify: `GTO_ELFO/direct/elfo/casadi_energy_freetf.m`, `casadi_mintime_freetf.m` (registry + `returnModel` + `regHistory`, same pattern as Task 8; both already take `opts`)
- Create: `GTO_ELFO/direct/elfo/run_foc_elfo.m`

**Interfaces:**
- Produces: `rep = run_foc_elfo(kind, matFile)` with `kind ∈ {'fuel','mintime'}` selecting solver + manifest (`foc_manifest('elfo_fuel')` / `('elfo_mintime')`). Warm re-solves at the saved primal (`warmTight` per that solver's own convention), certified-quantity guard, `foc_check` + `foc_ipopt_inertia` + `foc_report`.

- [ ] **Step 1:** Registries in both solvers. Labels: `'defect'`, `'betaNorm'`; energy solver additionally `'thrLo'`,`'thrHi'`; min-time has no throttle rows (manifest `thrRow=[]` skips sign-law/singular/Sdot — correct for all-burn). `regHistory` capture after each solver's `opti.stats()` call, identical 2-line try/catch as Task 5.
- [ ] **Step 2:** Byte-identity smoke (same two-run pattern as Task 8 Step 2, `maxIter=5` warm at the certified seed `results/energy_elfo_freetf.mat`).
- [ ] **Step 3:** `run_foc_elfo.m` per interface.
- [ ] **Step 4: Run on three artifacts:** `results/energy_elfo_freetf.mat` (fuel manifest, ε=1 leg — pass `epsilon=1` on the re-solve and note the sign-law check is only exact at ε=0, so run it report-only with a printed ε caveat), the certified min-time anchor `results/mintime_elfo.mat` (min-time manifest — first-ever first-order report on that anchor; `lamTimeCoV`+`lamTimeEnd` give the free-t_f dual-form condition, G4), and one certified ε=0 front row from `results/elfo_batch_summary_minEps0.mat`'s row files. Expect KKT PASS on all three.
- [ ] **Step 5:** Wire the report into `run_elfo_minfuel.m` and `gen_elfo_mintime.m` end-of-run paths (final solve only, `returnModel=true`, try/catch-guarded advisory).
- [ ] **Step 6: Commit** `feat(foc): ELFO model hooks + run_foc_elfo -- first first-order gate for the ELFO campaign (G2/G4)`

---

### Task 10: Register, docs, memory, experiment log

**Files:**
- Modify: `orbit_transfer/OPTIMALITY_CERTIFICATION.md` (Part A: add `foc_check` to A1, refresh the A3 coverage matrix with the new green cells, mark G2/G3/G6 progress and G4 partial; §6 experiment-log entries for the earth ladder run, CR3BP ladder run, tulip flagship, ELFO trio; Part B §1: extend the IPOPT-inertia row's campaign column once ports run)
- Modify: `orbit_transfer/README.md` (Conventions: add "**Standard optimality report:** every production solve driver ends with `foc_report` — the fixed-format first-order block + `foc_<tag>.mat` sidecar; report-only burn-in, does not alter certified status. Core: `verify_common/`.")
- Modify: memory `second-order-optimality-state.md` (+ index line) — FOC layer built, coverage now X/Y, LEAD-0 ported.
- Modify: campaign TODOs — check off the items this closes (CR3BP "re-run verify_cr3bp_pmp" item; ELFO gains a verification bullet; earth TODO note under item 4's resolution).

- [ ] **Step 1:** Make all documentation edits, with the real numbers from Tasks 6–9's runs (no placeholders — copy the printed gate values).
- [ ] **Step 2:** Run the full `verify_common/tests/` suite one final time; all PASS/SKIP.
- [ ] **Step 3: Commit** `docs(foc): register + README + memory -- standard report convention recorded`

---

## Self-Review

- **Spec coverage:** (a) generic layer: Tasks 1–5 core, 6–9 all four campaigns — each campaign ends with a run on real certified artifacts. (b) standard per-solve report: Task 4 defines it; Tasks 6–9 wire it into `run_gergaud`, `verify_cr3bp_pmp`/`run_cr3bp_geo`, `run_psr`, `run_elfo_minfuel`+`gen_elfo_mintime`; Task 10 records the convention. Report-only policy appears in Global Constraints, `foc_check.pass` (advisory), and the printed ADVISORY line. LEAD-0 port: Tasks 5, 6, 8 (tulip already captures), 9. G-items: G2 (Tasks 7, 9), G3 (Sdot standing in every report), G4 (Task 9 dual-form, honestly labeled partial), G6 (Task 7 note — recursion inside kktStat, constancy ungated).
- **Placeholder scan:** Task 8/9 artifact layouts reference the campaign READMEs for exact field names rather than inventing them — that is a documented lookup, not a TBD; all code steps carry real code.
- **Type consistency:** `foc_check(out, sigma, man, opts)` used identically in Tasks 3, 6, 7, 8, 9; `rep` field names in Task 4's printout match Task 3's Produces list; `foc_ipopt_inertia(regHistory)` consistent across 5/6/8/9; manifest names `'earth_mee'|'earth_cr3bp'|'tulip'|'elfo_fuel'|'elfo_mintime'|'toy'` consistent between Tasks 2 and 6–9.

One known risk, stated: Tasks 8/9 assume tulip/ELFO variable declaration order is X-then-U (verified by reading the solvers: X declared immediately before U in all three) — `foc_check`'s layout assert turns any surprise into a loud failure, not a wrong answer.
