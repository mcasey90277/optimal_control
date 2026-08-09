### Task 5: Certification gates (`certify_pdg`)

**Files:**
- Create: `certify/certify_pdg.m`
- Create: `tests/test_certify_nominal.m`

**Interfaces:**
- Consumes: `sol` structs from both solvers (Tasks 3–4), `pdg_dynamics`.
- Produces: `rep = certify_pdg(solC, solV, P)` — gates report struct:
  - `rep.G1_defect` (max HS defect re-evaluated at solution), `rep.G1_pass`
  - `rep.G2_resid` (ode45 re-integration: terminal position/velocity/mass error), `rep.G2_pass` (pos < 1 m, vel < 0.1 m/s, mass < 0.5 kg)
  - `rep.G3_dmf` (|mf colloc − convex|), `rep.G3_dtf`, `rep.G3_traj_Linf` (position, common time grid), `rep.G3_pass` (dmf < 0.1 kg, dtf < 0.2 s)
  - `rep.G4_gap` (= solV.lossless_gap), `rep.G4_pass` (< 1e-4·Tmax/m0)
  - `rep.G5_structure` (throttle max–min–max: fraction of nodes on a bound ≥ 0.95, interior switch count ≤ 2), `rep.G5_primer_deg` (max angle between thrust and −λv), `rep.G5_pass` (structure ok AND primer < 1°)
  - `rep.all_pass`; `print_certify_report(rep)` nested pretty-printer (gate table, PASS/FAIL per row — the front door prints this)
- `solV = []` allowed (Phase 2 drag runs have no convex twin): G3/G4 report `skipped`, excluded from `all_pass`.

- [ ] **Step 1: Write the failing test**

`tests/test_certify_nominal.m`:

```matlab
% TEST_CERTIFY_NOMINAL  Full gate run on coarse solves (fast): all five
% gates must pass on the nominal vacuum problem. This is the campaign's
% core scientific claim -- two independent methods, one answer, PMP-shaped.
%
% INPUTS: none   OUTPUTS: none (throws on failure)
here_ = fileparts(mfilename('fullpath'));
addpath(fullfile(here_, '..'));  setup_paths;

P    = booster_params();
solC = solve_pdg_colloc(P, struct('N', 30));
solV = solve_pdg_convex(P, struct('Nconv', 90));
rep  = certify_pdg(solC, solV, P);
print_certify_report(rep);
assert(rep.all_pass, 'certification failed -- see report above');
fprintf('test_certify_nominal PASS\n');
```

(Gate tolerances at coarse grids: if G3 dmf < 0.1 kg proves too tight at N=30 while fine grids pass, loosen the COARSE-test call site via an optional `tol` scale arg to `certify_pdg`, never the nominal tolerance itself.)

- [ ] **Step 2: Run test to verify it fails**

Run: `/Applications/MATLAB_R2025b.app/bin/matlab -batch "cd('/Users/msc/Desktop/optimal_control/booster_landing'); setup_paths; test_certify_nominal"`
Expected: FAIL — `certify_pdg` undefined.

- [ ] **Step 3: Implement `certify/certify_pdg.m`**

Key pieces (full header per standard; `print_certify_report` in the same file as a second function):

```matlab
function rep = certify_pdg(solC, solV, P, tolScale)
% CERTIFY_PDG  Gates G1-G5 for a PDG solution pair (colloc + convex).
% See plan Task 5 for gate definitions; report-only, throws never --
% callers decide what a FAIL means.
% INPUTS: solC colloc sol; solV convex sol or [] (G3/G4 skipped);
%         P params; tolScale (optional, def 1) scales agreement tols for
%         coarse-grid test calls.
% OUTPUTS: rep - gate report struct (fields in plan Task 5).
if nargin < 4, tolScale = 1; end

%% G1: HS defects re-evaluated in plain MATLAB (independent of CasADi):
N = size(solC.X,2) - 1;  h = solC.tf / N;  dmax = 0;
for k = 1:N
    fk  = pdg_dynamics(solC.X(:,k),   solC.U(:,k),   P);
    fk1 = pdg_dynamics(solC.X(:,k+1), solC.U(:,k+1), P);
    xm  = 0.5*(solC.X(:,k)+solC.X(:,k+1)) + (h/8)*(fk - fk1);
    fm  = pdg_dynamics(xm, solC.Um(:,k), P);
    d   = solC.X(:,k+1) - solC.X(:,k) - (h/6)*(fk + 4*fm + fk1);
    dmax = max(dmax, max(abs(d)));
end
rep.G1_defect = dmax;   rep.G1_pass = dmax < 1e-6;

%% G2: continuous residual -- fly the interpolated control with ode45
%% ("defect is not accuracy"): pchip on [nodes+midpoints] thrust.
tU = sort([solC.t, solC.t(1:end-1) + h/2]);
TU = zeros(3, numel(tU));
TU(:,1:2:end) = solC.U;  TU(:,2:2:end) = solC.Um;
ctrl = @(tt) interp1(tU.', TU.', min(tt, solC.tf), 'pchip').';
odef = @(tt, xx) pdg_dynamics(xx, ctrl(tt), P);
oo   = odeset('RelTol',1e-10,'AbsTol',1e-10);
[~, XX] = ode45(odef, [0 solC.tf], solC.X(:,1), oo);
ef = XX(end,:).' - solC.X(:,end);
rep.G2_pos = sqrt(sum(ef(1:3).^2));  rep.G2_vel = sqrt(sum(ef(4:6).^2));
rep.G2_dm  = abs(ef(7));
rep.G2_pass = rep.G2_pos < 1 && rep.G2_vel < 0.1 && rep.G2_dm < 0.5;

%% G3/G4: cross-method agreement + losslessness (skip if no convex twin):
if isempty(solV)
    rep.G3_pass = 'skipped';  rep.G4_pass = 'skipped';
else
    rep.G3_dmf = abs(solC.mf - solV.mf);
    rep.G3_dtf = abs(solC.tf - solV.tf);
    tq  = linspace(0, min(solC.tf, solV.tf), 200);
    rC  = interp1(solC.t.', solC.X(1:3,:).', tq.', 'pchip');
    rV  = interp1(solV.t.', solV.X(1:3,:).', tq.', 'pchip');
    rep.G3_traj_Linf = max(sqrt(sum((rC - rV).^2, 2)));
    rep.G3_pass = rep.G3_dmf < 0.1*tolScale && rep.G3_dtf < 0.2*tolScale;
    rep.G4_gap  = solV.lossless_gap;
    rep.G4_pass = rep.G4_gap < 1e-4 * P.Tmax / P.m0;
end

%% G5: PMP structure. Throttle bang-bang max-min-max + primer alignment.
Tmag = sqrt(sum(solC.U.^2, 1));
onLo = Tmag < P.Tmin * 1.001;   onHi = Tmag > P.Tmax * 0.999;
rep.G5_bound_frac = mean(onLo | onHi);
segs = diff([onHi(1), onHi]);                  % max<->min transitions
rep.G5_switches = sum(segs ~= 0);
structOk = rep.G5_bound_frac >= 0.95 && rep.G5_switches <= 2 ...
           && onHi(1) && onHi(end);            % max-first, max-last
% Primer: velocity costate from defect duals; thrust must align with -lam_v.
lamv = solC.lam_defect(4:6, :);                % 3xN, node-adjacent duals
Tdir = solC.U(:,1:end-1) ./ sqrt(sum(solC.U(:,1:end-1).^2,1));
pdir = -lamv ./ max(sqrt(sum(lamv.^2,1)), 1e-30);
cosang = sum(Tdir .* pdir, 1);
rep.G5_primer_deg = max(acosd(min(1, max(-1, cosang))));
rep.G5_pass = structOk && rep.G5_primer_deg < 1;

%% Verdict:
gates = [rep.G1_pass, rep.G2_pass, isequal(rep.G3_pass,true) || ...
         isequal(rep.G3_pass,'skipped'), isequal(rep.G4_pass,true) || ...
         isequal(rep.G4_pass,'skipped'), rep.G5_pass];
rep.all_pass = all(gates);
end
```

Notes:
- The primer-dual sign: if `G5_primer_deg` comes out near 180° instead of ~0°, the dual sign convention is flipped — flip `pdir` sign ONCE, document why in a comment referencing the `opti.lam_g` house lesson, and verify the switching structure still reads max–min–max. Do not add an `abs()` to make the test pass (that would hide a real sign error).
- `print_certify_report(rep)`: plain fprintf table, one row per gate, value + threshold + PASS/FAIL/skipped. Write it in the same file below `certify_pdg` — MATLAB allows multiple functions per file; expose it by making `certify/print_certify_report.m` a 3-line wrapper if needed elsewhere.

- [ ] **Step 4: Run test to verify it passes**

Same command. Expected: gate table printed, all PASS, `test_certify_nominal PASS`. If G5 primer fails on sign, apply the note above; if G3 fails, this is the moment the two methods genuinely disagree — systematic-debugging, likely suspects: Taylor mass-bound error (refine `Nconv`), or tf bracket edge.

- [ ] **Step 5: Commit**

```bash
cd /Users/msc/Desktop/optimal_control
git add booster_landing/certify/certify_pdg.m booster_landing/tests/test_certify_nominal.m
git commit -m "booster_landing: certification gates G1-G5 (defect, residual, cross-method, lossless, PMP)

Co-Authored-By: Claude Fable 5 <noreply@anthropic.com>"
```

---

