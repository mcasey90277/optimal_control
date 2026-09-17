# Cart-Pole Sub-Project A Implementation Plan — shared plant, the conjugate-test promotion, and the min-energy study

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Put the cart-pole plant in one shared home, promote `ms_conjugate_test` into `oclib` with these examples as its second top-level consumer, and write `cartpole_minenergy_study.m` — the first of three entry scripts in the `transfer_study.m` style — against the already-working min-energy solve.

**Architecture:** No new mathematics. The min-energy solve exists and is verified; this sub-project rearranges what it stands on (a shared `cartpole_common/`), borrows one instrument from the orbit library through a legitimate promotion, and exposes the whole computation as a numbered study script whose every condition is computed inline with its own PASS/FAIL and then asserted against the library instrument. The skeleton this builds is what sub-projects B (min-time) and C (min-fuel) instantiate.

**Tech Stack:** MATLAB R2026a, `ode113`, `oclib/+oc` (`ms_bvp`, `ms_conjugate_test` after Task 3, `duals_to_costates`, `fly_control`). Symbolic Math Toolbox only if the Jacobian is regenerated. No CasADi, no pumpkyn, no campaign folder.

**Spec:** `docs/superpowers/specs/2026-09-17-cartpole-three-objectives-design.md`

## Global Constraints

- Plant constants exactly: `m1 = 5` kg, `m2 = 1` kg, `L = 2` m, `g = 9.8` m/s². For this objective `t_f = 5` s fixed, `x(0) = [0;0;0;0]`, `x(t_f) = [0;pi;0;0]`, control unconstrained.
- House MATLAB style: `%% Purpose: / %% Inputs: / %% Outputs: / %% Revision History:` quartet closed by `%% ------------------------ Begin Code Sequence ---------------------------`; `nargin == 0` self-demo on library functions; **no** `%#ok` pragmas anywhere; never `i`/`j` as loop variables. Author `%  M. Casey`, `(c) 09/17/2026`, `%  Copyright Coorbital Inc.`
- **Scripts are exempt from the function-header quartet** but carry the `transfer_study.m` header form: `%% NAME  one-line purpose`, a prerequisite paragraph, the numbered section map, the diagnostic-ID legend, then author and copyright.
- `optimal_control_examples/` stays campaign-free: no file under it may reference `orbit_transfer`, a CR3BP quantity, or pumpkyn. `oclib` is allowed and is the point.
- Nothing under `ex2_cart_pole_swing_up/` may be modified.
- Run MATLAB only as `/Applications/MATLAB_R2026a.app/bin/matlab -batch "..."`. Do **not** use the MATLAB MCP tool: the shared session in this repo is wedged and hangs for 30 minutes.
- Committed `.mat` files need `git add -f` (the repo gitignores `*.mat`).
- Do not widen a tolerance to make a check pass. Measure, then set it with the measurement recorded in the file.

---

### Task 1: Extract the plant into `cartpole_common`

**Files:**
- Move: `optimal_control_examples/ex3_cart_pole_pmp/{cartpole_field.m, cartpole_state_jac.m, cartpole_state_jac_gen.m, gen_state_jac.m}` → `optimal_control_examples/cartpole_common/`
- Move: `optimal_control_examples/ex3_cart_pole_pmp/tests/{test_cartpole_physics.m, test_cartpole_field.m, test_cartpole_state_jac.m}` → `optimal_control_examples/cartpole_common/tests/`
- Create: `optimal_control_examples/cartpole_common/cartpole_params.m`
- Create: `optimal_control_examples/cartpole_common/tests/test_cartpole_params.m`
- Create: `optimal_control_examples/cartpole_common/README.md`
- Modify: `ex3_cart_pole_pmp/{cartpole_pmp_rhs.m, cartpole_pmp_prop.m, run_cartpole_pmp.m, run_tests.m, gen_direct_ref.m, movie_cartpole.m, README.md}` and every file under `ex3_cart_pole_pmp/tests/` — path bootstraps only

**Interfaces:**
- Consumes: nothing new.
- Produces: `p = cartpole_params()` returning `struct('m1',5,'m2',1,'L',2,'g',9.8)` — the one home for the constants, which three examples and their tests call instead of each writing the literals. The moved files keep their signatures exactly: `[F, G] = cartpole_field(x, p)`, `A = cartpole_state_jac(x, u, p)`, `ok = test_cartpole_physics()`.

- [ ] **Step 1: Write the failing test**

Create `optimal_control_examples/cartpole_common/tests/test_cartpole_params.m`:

```matlab
function ok = test_cartpole_params()
%% Purpose:
%
%   The plant constants have ONE home. Three examples (min-energy, min-time,
%   min-fuel) and their tests share this plant; literals copied into each
%   would be the same failure mode as the duplicated dynamics that hid a sign
%   error until 2026-09-17.
%
%   Checks: the specified values; a struct with exactly the four fields the
%   field function reads, and no more (a fifth field would mean some caller
%   is smuggling a problem-specific quantity through the plant); and that
%   cartpole_field actually accepts what this returns.
%
%% Inputs:
%
%  none
%
%% Outputs:
%
%  ok                       logical                 All checks passed
%
%% Revision History:
%  M. Casey                                                   (c) 09/17/2026
%  Copyright Coorbital Inc.
%% ------------------------ Begin Code Sequence ---------------------------

ok = true;
here = fileparts(fileparts(mfilename('fullpath')));
addpath(here);
p = cartpole_params();

ok = chk(ok, isstruct(p) && isscalar(p), 'returns a scalar struct');
ok = chk(ok, p.m1 == 5 && p.m2 == 1 && p.L == 2 && p.g == 9.8, ...
         sprintf('constants: m1 %g, m2 %g, L %g, g %g', p.m1, p.m2, p.L, p.g));
ok = chk(ok, isequal(sort(fieldnames(p)), sort({'m1'; 'm2'; 'L'; 'g'})), ...
         'exactly the four fields the plant reads, no more');

[F, G] = cartpole_field([0.1; 0.2; 0.3; 0.4], p);
ok = chk(ok, isequal(size(F), [4 1]) && isequal(size(G), [4 1]) && all(isfinite([F; G])), ...
         'cartpole_field accepts it and returns finite [4 x 1] columns');

if ok, fprintf('TEST_CARTPOLE_PARAMS: ALL PASS\n');
else,  fprintf('TEST_CARTPOLE_PARAMS: FAILURE (see lines above)\n');
end
end

% ------------------------------------------------------------------------
function ok = chk(ok, cond, label)
%% Purpose:
%
%   Accumulate one labelled pass/fail.
%
if cond, fprintf('  PASS  %s\n', label);
else,    fprintf('  FAIL  %s\n', label);  ok = false;
end
end
```

- [ ] **Step 2: Run it to verify it fails**

```bash
cd /Users/msc/Desktop/optimal_control/optimal_control_examples/cartpole_common/tests
/Applications/MATLAB_R2026a.app/bin/matlab -batch "cd('$(pwd)'); ok = test_cartpole_params()"
```

Expected: FAIL — neither `cartpole_common/` nor `cartpole_params` exists yet.

- [ ] **Step 3: Move the files and write `cartpole_params`**

```bash
cd /Users/msc/Desktop/optimal_control/optimal_control_examples
mkdir -p cartpole_common/tests
git mv ex3_cart_pole_pmp/cartpole_field.m         cartpole_common/
git mv ex3_cart_pole_pmp/cartpole_state_jac.m     cartpole_common/
git mv ex3_cart_pole_pmp/cartpole_state_jac_gen.m cartpole_common/
git mv ex3_cart_pole_pmp/gen_state_jac.m          cartpole_common/
git mv ex3_cart_pole_pmp/tests/test_cartpole_physics.m   cartpole_common/tests/
git mv ex3_cart_pole_pmp/tests/test_cartpole_field.m     cartpole_common/tests/
git mv ex3_cart_pole_pmp/tests/test_cartpole_state_jac.m cartpole_common/tests/
```

`cartpole_common/cartpole_params.m`:

```matlab
function p = cartpole_params()
%% Purpose:
%
%   THE cart-pole's physical constants, in one place. Every example that
%   solves this plant -- minimum energy, minimum time, minimum fuel -- and
%   every test of it reads them from here rather than repeating literals.
%
%   The same rule as the dynamics themselves: one home. A sign error survived
%   in this plant until 2026-09-17 because the equations existed in several
%   copies and every test compared one copy against another.
%
%% Inputs:
%
%  none
%
%% Outputs:
%
%  p                        struct                  .m1 cart mass (kg), .m2
%                                                   bob mass (kg), .L
%                                                   pendulum length (m), .g
%                                                   gravity (m/s^2)
%
%% Revision History:
%  M. Casey                                                   (c) 09/17/2026
%  Copyright Coorbital Inc.
%% ------------------------ Begin Code Sequence ---------------------------

if nargin == 0 && nargout == 0
   %Demo: the constants, and the hanging small-oscillation period they imply:
     pd = cartpole_params();
     w  = sqrt((pd.m1 + pd.m2)*pd.g/(pd.L*pd.m1));
     fprintf('m1 %g kg, m2 %g kg, L %g m, g %g m/s^2  ->  hanging period %.4f s\n', ...
             pd.m1, pd.m2, pd.L, pd.g, 2*pi/w);
     return
end

p = struct('m1', 5, 'm2', 1, 'L', 2, 'g', 9.8);
end
```

Then fix the path bootstraps. Files directly in `ex3_cart_pole_pmp/` currently do
`here = fileparts(mfilename('fullpath')); addpath(here, ...)`; each must also add
the sibling folder:

```matlab
here = fileparts(mfilename('fullpath'));
addpath(here, fullfile(fileparts(here), 'cartpole_common'));
```

Files in `ex3_cart_pole_pmp/tests/` need one more level up:

```matlab
here = fileparts(fileparts(mfilename('fullpath')));      % the example folder
addpath(here, fullfile(here, 'tests'), ...
        fullfile(fileparts(here), 'cartpole_common'));
```

The three moved tests resolve `cartpole_common` as their own parent:

```matlab
here = fileparts(fileparts(mfilename('fullpath')));      % cartpole_common
addpath(here);
```

`run_cartpole_pmp.m` also computes `root = fileparts(fileparts(here))` to reach
`oclib` — that line stays correct and must not be touched.

**No function signature changes.** The only edits outside the moves are `addpath`
lines and `run_tests.m`'s list of test names.

Create `cartpole_common/README.md`: what lives here (the plant: field, its
Jacobian and generator, the constants, the physics oracle), the rule for adding
to it (it must be objective-independent — anything that knows about a cost
belongs in an example folder), and the list of consumers.

- [ ] **Step 4: Run the new test, then the ex3 suite**

```bash
cd /Users/msc/Desktop/optimal_control/optimal_control_examples/cartpole_common/tests
/Applications/MATLAB_R2026a.app/bin/matlab -batch "cd('$(pwd)'); exit(~test_cartpole_params())"
cd ../../ex3_cart_pole_pmp
/Applications/MATLAB_R2026a.app/bin/matlab -batch "cd('$(pwd)'); exit(~run_tests())"
```

Expected: the new test passes, and `run_tests` reports the same tally and the same
numbers as before the move — terminal miss 3.67e-11, J 2779.381719. **If any number
moved, stop:** a file move that changes a result is not a file move. Report the
before/after numbers either way.

- [ ] **Step 5: Commit**

```bash
cd /Users/msc/Desktop/optimal_control
git add -A optimal_control_examples/cartpole_common optimal_control_examples/ex3_cart_pole_pmp
git commit -m "cart-pole: the plant moves to cartpole_common, before a third example copies it"
```

---

### Task 2: A shared runner for the example tree

**Files:**
- Create: `optimal_control_examples/run_all_tests.m`

**Interfaces:**
- Consumes: `cartpole_common/tests/*`, `ex3_cart_pole_pmp/run_tests`.
- Produces: `ok = run_all_tests()` — runs every suite, prints one summary, and **throws** if any fails. Sub-projects B and C add rows to its `suites` table.

- [ ] **Step 1: Write the runner**

There is no test-of-the-runner; the runner is the harness, and Step 3 proves it
can fail. Create `optimal_control_examples/run_all_tests.m`:

```matlab
function ok = run_all_tests()
%% Purpose:
%
%   Every suite under optimal_control_examples, with ONE exit code. A test
%   that prints FAIL and returns false still exits 0 on its own, which is
%   fine at the prompt and useless in automation; this throws.
%
%   Folder order is dependency order: the shared plant first -- if the
%   physics is wrong, nothing downstream means anything -- then each example.
%
%% Inputs:
%
%  none
%
%% Outputs:
%
%  ok                       logical                 Every suite passed (the
%                                                   function also THROWS on
%                                                   failure, so a -batch run
%                                                   returns a real exit code)
%
%% Revision History:
%  M. Casey                                                   (c) 09/17/2026
%  Copyright Coorbital Inc.
%% ------------------------ Begin Code Sequence ---------------------------

here = fileparts(mfilename('fullpath'));
addpath(here, fullfile(fileparts(here), 'oclib'));

%% The suites, in dependency order. Each entry is a folder and the functions
%% in it that return a logical verdict (an example's own run_tests is just
%% another such function):
suites = { 'cartpole_common',   {'test_cartpole_params', 'test_cartpole_physics', ...
                                 'test_cartpole_field',  'test_cartpole_state_jac'}
           'ex3_cart_pole_pmp', {'run_tests'} };

nT      = sum(cellfun(@numel, suites(:,2)));
names   = cell(1, nT);
res     = false(1, nT);
secs    = zeros(1, nT);
n       = 0;
for k = 1:size(suites, 1)
    folder = suites{k,1};
    addpath(fullfile(here, folder), fullfile(here, folder, 'tests'));
    for m = 1:numel(suites{k,2})
        n = n + 1;
        names{n} = [folder '/' suites{k,2}{m}];
        fprintf('\n######## %s\n', names{n});
        t0 = tic;
        try
            res(n) = logical(feval(suites{k,2}{m}));
        catch err
            fprintf('  THREW  %s\n', err.message);
            res(n) = false;
        end
        secs(n) = toc(t0);
    end
end

fprintf('\n==== ALL SUITES ====\n');
for k = 1:nT
    if res(k), tag = 'PASS'; else, tag = 'FAIL'; end
    fprintf('  %-44s %-4s %6.1f s\n', names{k}, tag, secs(k));
end
ok = all(res);
fprintf('  %d of %d passed, %.1f s total\n', nnz(res), nT, sum(secs));
if ~ok
    error('run_all_tests:failed', '%d suite(s) FAILED: %s', ...
          nnz(~res), strjoin(names(~res), ', '));
end
end
```

- [ ] **Step 2: Run it and see the current tree pass**

```bash
cd /Users/msc/Desktop/optimal_control/optimal_control_examples
/Applications/MATLAB_R2026a.app/bin/matlab -batch "cd('$(pwd)'); exit(~run_all_tests())"
echo "exit code $?"
```

Expected: every suite PASS, exit code 0. Confirm with `grep -c '%#ok' run_all_tests.m`
that the file carries no pragmas (expected: 0).

- [ ] **Step 3: Prove it fails when something fails**

Copy `cartpole_common/cartpole_field.m` to a scratch folder OUTSIDE the repo,
perturb `G(3)` there (multiply by 1.01), put that folder first on the path, and
run again:

```bash
SCRATCH=$(mktemp -d)
cp /Users/msc/Desktop/optimal_control/optimal_control_examples/cartpole_common/cartpole_field.m "$SCRATCH/"
# edit $SCRATCH/cartpole_field.m: G(3) -> 1.01/D1
cd /Users/msc/Desktop/optimal_control/optimal_control_examples
/Applications/MATLAB_R2026a.app/bin/matlab -batch "addpath('$SCRATCH'); cd('$(pwd)'); exit(~run_all_tests())"
echo "exit code $?"
rm -rf "$SCRATCH"
```

Expected: a non-zero exit code and `run_all_tests:failed` naming the failed
suites — in particular `cartpole_common/test_cartpole_physics`, since the power
balance is what a wrong `G` breaks. Then confirm `git status` shows no repo change.

- [ ] **Step 4: Commit**

```bash
cd /Users/msc/Desktop/optimal_control
git add optimal_control_examples/run_all_tests.m
git commit -m "optimal_control_examples: one runner for the whole tree, with a real exit code"
```

---

### Task 3: Promote `ms_conjugate_test` into `oclib`

**Files:**
- Move: `orbit_transfer/costate_common/ms_conjugate_test.m` → `oclib/+oc/ms_conjugate_test.m`
- Create: `orbit_transfer/costate_common/ms_conjugate_test.m` (delegate)
- Modify: `oclib/README.md`, `orbit_transfer/costate_common/README.md`, `orbit_transfer/costate_common/TODO.md`

**Interfaces:**
- Consumes: nothing from this plan.
- Produces: `out = oc.ms_conjugate_test(info, spec)`, signature unchanged — `info` is an `ms_bvp` info struct carrying `.PHI {1 x K}`, `.Y [ny x K]`, `.tGrid`, `.Yend`; `spec` carries `.flow`, `.stateRows`, `.costateCols`, `.quotientDir`, `.freeTime`, `.rankTol`, `.zeroTol`, `.resolvedTol`, `.kernelTol`; `out` carries `.t`, `.detScaled`, `.sigRatio`, `.firstFullRank`, `.tested`, `.nCrossings`, `.nInterior`, `.nTouch`, `.atFinal`, `.nUnresolved`, `.covered`, `.tFirstFullRank`, `.kernelRight`. Task 4 calls it. The delegate forwards, so every orbit caller is untouched.

**Precedent to copy exactly:** `orbit_transfer/costate_common/ms_bvp.m`, the delegate written when `ms_bvp` was promoted on 2026-09-16. Read it first and match its wording and structure.

- [ ] **Step 1: Record the baseline the move must reproduce**

```bash
cd /Users/msc/Desktop/optimal_control/orbit_transfer/costate_common
/Applications/MATLAB_R2026a.app/bin/matlab -batch "cd('$(pwd)'); ok = golden_cells(); fprintf('GOLDEN %d\n', ok)" | tail -30
```

Save the printed residuals and iteration counts into the task report; Step 3
compares against them line by line.

- [ ] **Step 2: Move it and write the delegate**

```bash
cd /Users/msc/Desktop/optimal_control
git mv orbit_transfer/costate_common/ms_conjugate_test.m oclib/+oc/ms_conjugate_test.m
```

`orbit_transfer/costate_common/ms_conjugate_test.m`:

```matlab
function out = ms_conjugate_test(info, spec)
%% Purpose:
%
%   DELEGATE since 2026-09-17: the Jacobi (conjugate-point) test was promoted
%   to the cross-folder library, oclib/+oc/ms_conjugate_test, when the
%   cart-pole study scripts (optimal_control_examples) became its second
%   TOP-LEVEL consumer. This file keeps every costate_common caller working;
%   the instrument, its conventions and its tests live with the
%   implementation.
%
%% Inputs:
%
%  info                     struct                  see oc.ms_conjugate_test
%
%  spec                     struct (optional)       see oc.ms_conjugate_test
%
%% Outputs:
%
%  out                      struct                  see oc.ms_conjugate_test
%
%% Revision History:
%  M. Casey                                                   (c) 09/17/2026
%  Copyright Coorbital Inc.
%% ------------------------ Begin Code Sequence ---------------------------

if isempty(which('oc.ms_conjugate_test'))
    addpath(fullfile(fileparts(fileparts(fileparts(mfilename('fullpath')))), 'oclib'));
end
if nargin < 2, spec = struct(); end
out = oc.ms_conjugate_test(info, spec);
end
```

Inspect the promoted file for anything that made sense only inside
`costate_common` — a path bootstrap naming a sibling, a header cross-reference —
and fix only that. **No executable line of the instrument's mathematics may
change.** Verify with `git diff -M --stat` that the move is detected as a rename.

- [ ] **Step 3: Run every gate**

```bash
cd /Users/msc/Desktop/optimal_control/orbit_transfer/costate_common/tests
/Applications/MATLAB_R2026a.app/bin/matlab -batch "cd('$(pwd)'); ts = {'test_conj_fixedtf','test_conj_coverage','test_conj_spectrum','test_conj_resolve','test_folder_rules'}; bad = 0; for k = 1:numel(ts), r = feval(ts{k}); fprintf('%s %d\n', ts{k}, r); bad = bad + ~r; end; exit(bad)"
cd ..
/Applications/MATLAB_R2026a.app/bin/matlab -batch "cd('$(pwd)'); ok = golden_cells(); fprintf('GOLDEN %d\n', ok)" | tail -30
```

Then the four-campaign ladder re-solve — the same gate the `ms_bvp` promotion
used, run from the scratchpad gate job:

```bash
cd /private/tmp/claude-501/-Users-msc-Desktop-optimal-control/*/scratchpad
GATE_PHASE=conj_promotion nohup /Applications/MATLAB_R2026a.app/bin/matlab -batch "run('gate_job.m')" > conj_gate.log 2>&1 &
# poll conj_gate.log; do NOT use `timeout` (it does not exist on macOS)
grep -E "^== |max\|dTF\|" conj_gate.log
```

Expected: all five tests green; `golden_cells` 20/20 with the SAME numbers as
Step 1; the campaigns at HALO 4.551e-14, DPO 1.130e-14, HALO_HALO 3.478e-14,
GTO 2.720e-14. **If any number moves, stop and report it** — that is a finding
about the move, not a tolerance to adjust.

- [ ] **Step 4: Documentation and commit**

`oclib/README.md`: add the `oc.ms_conjugate_test` row in the form the
`oc.ms_bvp` row uses — what it is (the Jacobi/conjugate-point test, both the
fixed-time and free-time conventions), its consumers (`orbit_transfer` through
the costate_common delegate; `optimal_control_examples` study scripts directly),
and the equivalence gate actually run with its numbers. Mark the roadmap's
promotion item done for this function; `arclength_ms`, `newton_fixed_q`,
`conj_resolve` and `lift_space_dim` stay open with their existing reasons.

`costate_common/README.md`: the `ms_conjugate_test` row becomes DELEGATE, worded
as the `ms_bvp` row is. `costate_common/TODO.md`: step 13 records this half done
and names what remains.

```bash
cd /Users/msc/Desktop/optimal_control
git add oclib/+oc/ms_conjugate_test.m oclib/README.md \
        orbit_transfer/costate_common/ms_conjugate_test.m \
        orbit_transfer/costate_common/README.md orbit_transfer/costate_common/TODO.md
git commit -m "oclib: promote ms_conjugate_test, with the cart-pole studies as its second consumer"
```

---

### Task 4: `cartpole_minenergy_study.m`

**Files:**
- Create: `optimal_control_examples/ex3_cart_pole_pmp/cartpole_minenergy_study.m`
- Create: `optimal_control_examples/ex3_cart_pole_pmp/tests/test_minenergy_study.m`
- Modify: `optimal_control_examples/ex3_cart_pole_pmp/run_tests.m` (add the new test), `README.md`

**Interfaces:**
- Consumes: `cartpole_params`, `cartpole_field`, `test_cartpole_physics` (Task 1); `oc.ms_bvp`, `oc.ms_conjugate_test` (Task 3), `oc.duals_to_costates`, `oc.fly_control`; `cartpole_pmp_rhs`, `cartpole_pmp_prop` from ex3. `run_cartpole_pmp` is deliberately NOT called — the study performs the steps itself, which is the whole point of the style; the test then asserts the two agree.
- Produces: a SCRIPT (not a function). Running it prints the numbered sections and leaves `tol`, `it`, `J`, `lam0`, `conjOut`, `verdict` in the workspace. `test_minenergy_study` runs it and asserts on those.

**The style to follow** is `orbit_transfer/DRO_tulip/indirect/transfer_study.m`. Read its first 140 lines before writing. Required, not optional: one `tol` struct in section 0 with every field commented with what it gates and why that value; every condition computed inline and printed as `value / threshold  PASS|FAIL`; section 8 asserting the inline numbers against the library instrument; ONE flight shared by sections 6–9.

- [ ] **Step 1: Write the failing test**

Create `optimal_control_examples/ex3_cart_pole_pmp/tests/test_minenergy_study.m`:

```matlab
function ok = test_minenergy_study()
%% Purpose:
%
%   The study script must RUN and must REACH A VERDICT. A study that prints
%   numbers and exits 0 regardless is a report, not a check.
%
%   Checks: it runs to the end; it leaves a verdict; every NECESSARY gate
%   passed; the verdict's claim is not made unless its own gates support it
%   (a verdict that claims more than it measured is the failure mode this
%   exists to catch); and the headline numbers match run_cartpole_pmp on the
%   same problem -- the study and the front door must not drift apart.
%
%% Inputs:
%
%  none
%
%% Outputs:
%
%  ok                       logical                 All checks passed
%
%% Revision History:
%  M. Casey                                                   (c) 09/17/2026
%  Copyright Coorbital Inc.
%% ------------------------ Begin Code Sequence ---------------------------

ok   = true;
here = fileparts(fileparts(mfilename('fullpath')));
root = fileparts(fileparts(here));
addpath(here, fullfile(here, 'tests'), ...
        fullfile(fileparts(here), 'cartpole_common'), ...
        fullfile(fileparts(here), 'cartpole_common', 'tests'), ...
        fullfile(root, 'oclib'));

S = run_study_isolated(fullfile(here, 'cartpole_minenergy_study.m'));

ok = chk(ok, isstruct(S.verdict) && all(isfield(S.verdict, ...
             {'necessary', 'sufficiency', 'claim', 'why'})), ...
         'the script leaves a verdict with necessary / sufficiency / claim / why');
ok = chk(ok, S.verdict.necessary, 'every NECESSARY gate passed');
ok = chk(ok, ~S.verdict.claim || (S.verdict.necessary && S.verdict.sufficiency), ...
         'the claim is not made unless BOTH sections support it');

out = run_cartpole_pmp(struct('plot', false));
ok = chk(ok, out.ok, sprintf('the front door itself is ok (%s)', out.why));
ok = chk(ok, abs(S.J - out.J)/out.J < 1e-8, ...
         sprintf('study and front door agree on J: %.6f vs %.6f', S.J, out.J));
ok = chk(ok, max(abs(S.lam0 - out.lam0)) < 1e-6, ...
         sprintf('and on lam0 (max abs diff %.2e)', max(abs(S.lam0 - out.lam0))));

if ok, fprintf('TEST_MINENERGY_STUDY: ALL PASS\n');
else,  fprintf('TEST_MINENERGY_STUDY: FAILURE (see lines above)\n');
end
end

% ------------------------------------------------------------------------
function S = run_study_isolated(scriptPath)
%% Purpose:
%
%   Run a SCRIPT inside a function workspace and hand back the variables it
%   defined, so the test reads the script's own results without the caller's
%   workspace leaking into it.
%
run(scriptPath);
S = struct('verdict', verdict, 'J', J, 'lam0', lam0, 'conj', conjOut);
end

function ok = chk(ok, cond, label)
%% Purpose:
%
%   Accumulate one labelled pass/fail.
%
if cond, fprintf('  PASS  %s\n', label);
else,    fprintf('  FAIL  %s\n', label);  ok = false;
end
end
```

- [ ] **Step 2: Run it to verify it fails**

```bash
cd /Users/msc/Desktop/optimal_control/optimal_control_examples/ex3_cart_pole_pmp/tests
/Applications/MATLAB_R2026a.app/bin/matlab -batch "cd('$(pwd)'); ok = test_minenergy_study()"
```

Expected: FAIL — `cartpole_minenergy_study.m` does not exist.

- [ ] **Step 3: Write the study script**

Create `optimal_control_examples/ex3_cart_pole_pmp/cartpole_minenergy_study.m`.
The seed block (sections 3) is lifted VERBATIM from `run_cartpole_pmp.m:73-109`;
do not paraphrase it, and keep its comments — they record measurements.

```matlab
%% CARTPOLE_MINENERGY_STUDY  Minimum-effort cart-pole swing-up, step by step.
%
%  Prerequisite: optimal_control_examples/cartpole_common and oclib on the
%  path (this script adds both), and data/cartpole_direct_ref.mat present --
%  the committed direct-collocation solve that seeds the shooting. Nothing
%  here touches orbit_transfer.
%
%  Sections:
%    0  tolerances          one struct, feeding every printed line AND verdict
%    1  the plant           constants, and the physics ORACLE run here
%    2  the problem         cost, boundary conditions, the control law
%    3  the seed            direct multipliers -> costates, sign resolved
%    4  solve               oc.ms_bvp at fixed t_f, STMs kept
%    5  independent check   a perturbed seed must find the same root
%    6  NECESSARY  (N1-N6)  Pontryagin, one condition at a time
%    7  SUFFICIENCY (S1-S2) strict Legendre, conjugate point
%    8  assert vs library   the inline numbers against oc.ms_conjugate_test
%    9  plot                states, control, switching-free sigma, H(t)
%
%  Diagnostic IDs:  N necessary | S sufficiency | V validity | X cross-check.
%  A failed N or X gate means the root or the implementation is broken. A
%  failed or unresolved S gate is a FINDING about this trajectory: reported,
%  and thrown only when selfCheck.strict is set.
%
%  M. Casey                                                   (c) 09/17/2026
%  Copyright Coorbital Inc.
%% ------------------------------------------------------------------------

here = fileparts(mfilename('fullpath'));
addpath(here, fullfile(fileparts(here), 'cartpole_common'), ...
        fullfile(fileparts(here), 'cartpole_common', 'tests'), ...
        fullfile(fileparts(fileparts(here)), 'oclib'));

%% 0. Tolerances -- every gate below reads this struct and nothing else:
tol = struct( ...
    'R',        1e-10, ...   % ms_bvp residual. ABSOLUTE, and the costates are
    ...                      % O(1e3), so this is ~1e-13 relative; the
    ...                      % propagator's RelTol floors the achievable
    ...                      % residual near 4e-12 regardless of iterations.
    'H',        1e-9,  ...   % spread of H along the arc, relative. Measured
    ...                      % 5.77e-13 -- this is integration quality, not
    ...                      % optimality.
    'miss',     1e-9,  ...   % terminal miss of the re-flown arc. Measured
    ...                      % 3.67e-11 at RelTol 2.5e-14.
    'flown',    1e-6,  ...   % miss when flying the CLOSED-FORM control.
    ...                      % Measured 1.67e-07.
    'gap',      -1e-10, ...  % minimum-principle gap floor: H(u*+d) - H(u*)
    ...                      % must exceed this (i.e. be non-negative to
    ...                      % round-off) for every probe d.
    'adj',      1e-9,  ...   % adjoint equations, relative, vs complex-step
    ...                      % differentiation of H.
    'agree',    1e-6,  ...   % inline numbers vs the library instrument.
    'xJ',       2e-3,  ...   % direct vs indirect cost, relative. Measured
    ...                      % 7.5e-04 -- this is the DISCRETISATION gap
    ...                      % between the two methods, not solver error, so
    ...                      % it does not shrink with tolerances.
    'xU',       5e-2,  ...   % direct vs indirect control, RMS difference over
    ...                      % RMS level. Set from the measured value at
    ...                      % build time; record it in the report.
    'sigRatio', 1e-10);      % conjugate test: below this the determinant's
                             % sign is not trusted and the sample is reported
                             % unresolved rather than counted.
selfCheck = struct('strict', false);

%% 1. The plant, and the ORACLE that says the physics is right (V1):
p = cartpole_params();
fprintf('\n=== 1. plant: m1 %g kg, m2 %g kg, L %g m, g %g m/s^2\n', p.m1, p.m2, p.L, p.g);
V1 = test_cartpole_physics();     % power balance from GEOMETRIC energy, plus
                                  % equilibrium character and unforced
                                  % conservation -- an oracle OUTSIDE the
                                  % equations it checks. A sign error lived in
                                  % this plant until 2026-09-17 because every
                                  % other test compared the code with a copy
                                  % of itself.
fprintf('  V1 physics oracle                                       %s\n', pf(V1));
assert(V1, 'study:physics', 'the plant fails its own oracle: nothing below means anything');

%% 2. The problem:
R  = load(fullfile(here, 'data', 'cartpole_direct_ref.mat'));
tf = R.tf;
x0 = [0; 0; 0; 0];
xf = [0; pi; 0; 0];
fprintf('\n=== 2. minimise J = int_0^%g u^2 dt,  x(0) = [0 0 0 0],  x(t_f) = [0 pi 0 0]\n', tf);
fprintf('  control UNBOUNDED, so H = u^2 + lam''(F + G u) is minimised at its\n');
fprintf('  stationary point: u* = -lam''G/2, with H_uu = 2 > 0 exactly.\n');
fprintf('  t_f FIXED, so H is constant along the arc but NOT zero, and there\n');
fprintf('  is no transversality condition to check (N4 does not apply).\n');

%% 3. The seed: the direct solve's multipliers become costates
%  duals_to_costates resolves the GLOBAL SIGN by a primer vote over 3-vector
%  thrust directions. This problem's control is a scalar force, so there is
%  no such vote: the mapping is asked for the costates unsigned, and the sign
%  is fixed here by the same principle in scalar form -- the implied control
%  u = -lam'G/2 must agree in sign with the control the direct solve used.
[lamS, tS] = oc.duals_to_costates(struct('scheme', 'trapezoid', 'mu', R.muDefect, ...
                                         'tNodes', R.tN, 'velRows', 3:4));
tS       = tS(:).';
uImplied = zeros(1, numel(tS));
Xs       = interp1(R.tN, R.X.', tS, 'pchip').';
for k = 1:numel(tS)
    [~, Gk] = cartpole_field(Xs(:,k), p);
    uImplied(k) = -(lamS(:,k).'*Gk)/2;
end
uDirect  = reshape(interp1(R.tN, R.U, tS, 'pchip'), 1, []);
signCorr = sum(uImplied.*uDirect) / sqrt(sum(uImplied.^2)*sum(uDirect.^2));
if signCorr < 0, lamS = -lamS;  uImplied = -uImplied; end
ampRatio = sqrt(sum(uImplied.^2)/sum(uDirect.^2));
fprintf('\n=== 3. seed from the direct multipliers\n');
fprintf('  sign correlation %+.4f (flipped: %d), amplitude ratio %.4f\n', ...
        signCorr, signCorr < 0, ampRatio);
fprintf('  the correlation votes the SIGN only; the amplitude ratio is printed\n');
fprintf('  beside it because a wrong costate SCALE is invisible to that vote.\n');

K     = 8;
tGrid = linspace(0, tf, K+1);
Xg    = interp1(R.tN, R.X.', tGrid, 'pchip').';
Lg    = interp1(tS, lamS.', tGrid, 'pchip', 'extrap').';
seed  = struct('tf', tf, 'tGrid', tGrid, 'Y', [Xg; Lg]);
seed.Y(1:4,1) = x0;

%% 4. Solve: four unknowns lam(0) against four terminal conditions
prob = struct('ny', 8, 'freeIdx0', 5:8, ...
    'prop', @(dt, y0, needSTM) cartpole_pmp_prop(dt, y0, needSTM, p), ...
    'rhs',  @(y) cartpole_pmp_rhs(y, p), ...
    'terminal', @(y, needJ) studyTerminal(y, xf, needJ));
% keepSTMs: section 7 needs info.PHI. Without it the conjugate test has
% nothing to read.
[~, it] = oc.ms_bvp(prob, seed, struct('fixedTf', true, 'tolR', tol.R, ...
                                       'maxIter', 60, 'polishMax', 5, ...
                                       'keepSTMs', true));
lam0 = it.Y(5:8, 1);
fprintf('\n=== 4. multiple shooting: %d segments, converged %d, |R| = %.3e\n', ...
        K, it.converged, it.normR);
fprintf('  lam(0) = [%.6f %.6f %.6f %.6f]\n', lam0);

%% 5. Independent checks
%  X2 first, because it is about the solver and costs one more solve: a
%  different seed must find the SAME root. This guards the SOLVER (one basin,
%  not two); it cannot guard the physics, which is V1's job. X1 -- direct
%  against indirect -- needs the flight, so it is printed in section 6.
seed2 = seed;  seed2.Y(5:8,:) = 1.05*seed.Y(5:8,:);
[~, it2] = oc.ms_bvp(prob, seed2, struct('fixedTf', true, 'tolR', tol.R, ...
                                         'maxIter', 60, 'polishMax', 5));
dLam0 = max(abs(it2.Y(5:8,1) - lam0));
X2    = it2.converged && dLam0 < 1e-6;
fprintf('\n=== 5. perturbed seed (costates x 1.05)\n');
fprintf('  X2 same root: max|dlam(0)| = %.3e / 1.0e-06                %s\n', dLam0, pf(X2));

%% 6. NECESSARY conditions
%  ONE flight, read by sections 6 through 9. RelTol 2.5e-14, not the
%  propagator's 1e-12: this single shot re-flies the whole horizon from
%  lam0, so its own integration error lands directly in the terminal miss
%  below. Measured on this root: miss 6.88e-09 at RelTol 1e-12, 8.20e-10 at
%  1e-13, 3.67e-11 at 2.5e-14, while the engine's own last-arc residual is
%  7.27e-14 -- the loose figure was the MEASUREMENT, not the solution.
nPts   = 2001;
[t, Y] = ode113(@(tt, y) cartpole_pmp_rhs(y, p), linspace(0, tf, nPts), ...
                [x0; lam0], odeset('RelTol', 2.5e-14, 'AbsTol', 1e-16));
X   = Y(:,1:4).';
Lam = Y(:,5:8).';
U   = zeros(1, nPts);
H   = zeros(1, nPts);
sig = zeros(1, nPts);          % lam'G -- no switching structure here, but it
for k = 1:nPts                 % is the same object the bounded problems switch
    [Fk, Gk] = cartpole_field(X(:,k), p);
    sig(k)   = Lam(:,k).'*Gk;
    U(k)     = -sig(k)/2;
    H(k)     = U(k)^2 + Lam(:,k).'*(Fk + Gk*U(k));
end
J = trapz(t, U.^2);

fprintf('\n=== 6. NECESSARY\n');
N1 = it.converged && it.normR < tol.R;
fprintf('  N1 BVP residual         %.3e / %.1e                    %s\n', it.normR, tol.R, pf(N1));

Hrel = (max(H) - min(H))/max(abs(H));
N2   = Hrel < tol.H;
fprintf('  N2 H constant (rel)     %.3e / %.1e                    %s\n', Hrel, tol.H, pf(N2));
fprintf('     H = %.6f, constant but NOT zero: t_f is fixed.\n', mean(H));

missTerminal = max(abs(X(:,end) - xf));
N3 = missTerminal < tol.miss;
fprintf('  N3 terminal miss        %.3e / %.1e                    %s\n', missTerminal, tol.miss, pf(N3));

fprintf('  N4 transversality       n/a -- t_f and x(t_f) both fixed\n');

%  N5: the adjoint equation, lamdot = -dH/dx at fixed u, checked against
%  COMPLEX-STEP differentiation of H -- code that does not share a line with
%  the analytic Jacobian inside cartpole_pmp_rhs.
adjWorst = 0;
for k = round(linspace(1, nPts, 41))
    dy   = cartpole_pmp_rhs([X(:,k); Lam(:,k)], p);
    gH   = zeros(4,1);
    for m = 1:4
        xc     = complex(X(:,k));
        xc(m)  = xc(m) + 1e-30i;
        [Fc, Gc] = cartpole_field(xc, p);
        gH(m)  = imag(Lam(:,k).'*(Fc + Gc*U(k)))/1e-30;
    end
    adjWorst = max(adjWorst, max(abs(dy(5:8) + gH))/max(1, max(abs(gH))));
end
N5 = adjWorst < tol.adj;
fprintf('  N5 adjoint equations    %.3e / %.1e (rel)              %s\n', adjWorst, tol.adj, pf(N5));

%  N6: the minimum principle itself -- H(u* + d) - H(u*) >= 0 for every
%  probe d, at a scatter of times. Not a derivative check: a probe at finite
%  d would catch a stationary point that is a MAXIMUM, which dH/du = 0 would
%  not.
gapWorst = inf;
for k = round(linspace(1, nPts, 41))
    [Fk, Gk] = cartpole_field(X(:,k), p);
    for d = [-10 -1 -0.1 0.1 1 10]
        uP  = U(k) + d;
        HP  = uP^2 + Lam(:,k).'*(Fk + Gk*uP);
        gapWorst = min(gapWorst, HP - H(k));
    end
end
N6 = gapWorst > tol.gap;
fprintf('  N6 min-principle gap    %+.3e / %+.1e (worst)         %s\n', gapWorst, tol.gap, pf(N6));

%  The flown-control check: the CLOSED-FORM control on the flown state, not
%  an interpolant of U (interpolating u put the resampling error straight
%  into the miss: 8.9e-07 against a 1e-06 gate).
ppL  = pchip(t, Lam);
zEnd = oc.fly_control(x0, [0 tf], ...
    @(tt, x) studyFlyRhs(x, studyU(ppL, min(max(tt, t(1)), t(end)), x, p), p), ...
    struct('mode', 'span', 'solver', @ode113, 'RelTol', 1e-13, 'AbsTol', 1e-15));
missFlown = max(abs(zEnd - xf));
V2a = missFlown < tol.flown;
fprintf('  V2 flown-control miss   %.3e / %.1e                    %s\n', missFlown, tol.flown, pf(V2a));

%  X1: the DIRECT solve and this INDIRECT one are two different methods on
%  one problem. They must agree on the cost and on the control history. They
%  will NOT agree to solver tolerance -- the direct solve discretises then
%  optimises, so its cost carries the trapezoid rule's own error; the gate is
%  the measured agreement (J to 0.075%, control RMS to the same order), not
%  machine precision.
dJrel = abs(J - R.J)/R.J;
uDrms = sqrt(trapz(t, (U - reshape(interp1(R.tN, R.U, t, 'pchip'), 1, [])).^2)/tf);
uRms  = sqrt(trapz(t, U.^2)/tf);
X1    = (dJrel < tol.xJ) && (uDrms/uRms < tol.xU);
fprintf('  X1 direct vs indirect   dJ %.3e / %.1e,  du/u %.3e / %.1e   %s\n', ...
        dJrel, tol.xJ, uDrms/uRms, tol.xU, pf(X1));

%% 7. SUFFICIENCY
fprintf('\n=== 7. SUFFICIENCY\n');
S1 = true;      % H_uu = d2/du2 (u^2 + lam'(F + Gu)) = 2, EXACTLY, everywhere
fprintf('  S1 strict Legendre      H_uu = 2 > 0 (exact, not measured)    %s\n', pf(S1));

%  S2: the Jacobi / conjugate-point test. quotientDir [] and freeTime false:
%  t_f is FIXED, and the running cost u^2 breaks the costate-scaling
%  invariance the orbit min-time problems quotient out -- the same
%  convention costate_common/tests/test_conj_fixedtf.m uses.
conjOut = oc.ms_conjugate_test(it, struct('stateRows', 1:4, 'costateCols', 5:8, ...
                                       'quotientDir', [], 'freeTime', false, ...
                                       'resolvedTol', tol.sigRatio));
S2 = (conjOut.nInterior == 0) && conjOut.covered && (conjOut.nUnresolved == 0);
fprintf('  S2 conjugate point      interior %d, touch %d, unresolved %d, covered %d  %s\n', ...
        conjOut.nInterior, conjOut.nTouch, conjOut.nUnresolved, conjOut.covered, pf(S2));
if ~S2
    fprintf('     A conjugate point, or an unresolved sample, is a FINDING about\n');
    fprintf('     THIS trajectory -- not a tolerance to loosen. Report it.\n');
end

%% 8. ASSERT the inline numbers against the library instrument
fprintf('\n=== 8. inline vs library\n');
%  What the instrument computed, re-derived here where it is cheap: the
%  determinant it reports at the FINAL sample must match a determinant built
%  from the same STM by this script's own two lines.
PHIend = it.PHI{end};
for k = numel(it.PHI)-1:-1:1, PHIend = PHIend*it.PHI{k}; end
Mend   = PHIend(1:4, 5:8);
dMine  = det(Mend)/max(abs(det(Mend)), realmin);   % sign only: the instrument
dLib   = sign(conjOut.detScaled(end));                % scales, this script does not
V2b    = isequal(sign(dMine), dLib);
fprintf('  V2 det sign at t_f      inline %+d vs library %+d              %s\n', ...
        sign(dMine), dLib, pf(V2b));
%  Where re-deriving would mean a second unverified copy of delicate
%  machinery (the sweep, the bracketing, the kernel), the instrument's own
%  consistency fields are asserted instead:
V2c = conjOut.tested > 0 && numel(conjOut.t) == numel(conjOut.detScaled) && ...
      conjOut.tested == numel(conjOut.t) && islogical(conjOut.covered);
fprintf('  V2 instrument self-consistent (tested %d samples)             %s\n', ...
        conjOut.tested, pf(V2c));
V2 = V2a && V2b && V2c;

%% 9. Plot
figure('color', [1 1 1]);
subplot(3,1,1); plot(t, X(1,:), 'b', t, X(2,:), 'k'); grid on
ylabel('states'); legend({'q_1 (m)', 'q_2 (rad)'}, 'Location', 'best');
title(sprintf('cart-pole minimum effort: J = %.6f, t_f = %g s', J, tf));
subplot(3,1,2); plot(t, U, 'k'); grid on; ylabel('u (N)');
subplot(3,1,3); plot(t, H - mean(H), 'k'); grid on
xlabel('t (s)'); ylabel('H - mean(H)');

%% Verdict
verdict = struct('necessary', all([N1 N2 N3 N5 N6 X1 X2 V1 V2]), ...
                 'sufficiency', S1 && S2, 'claim', false, 'why', '');
verdict.claim = verdict.necessary && verdict.sufficiency;
if ~verdict.necessary
    verdict.why = 'a necessary or validity gate failed: the root or the implementation is wrong';
elseif ~verdict.sufficiency
    verdict.why = 'necessary conditions hold; the second-order test did not resolve or found a conjugate point';
else
    verdict.why = ['consistent with a strict strong local minimum for FIXED endpoints ' ...
                   'and FIXED t_f; numerical and sampled, not a certificate'];
end
fprintf('\n=== VERDICT: necessary %d, sufficiency %d, claim %d\n', ...
        verdict.necessary, verdict.sufficiency, verdict.claim);
fprintf('  %s\n', verdict.why);
if selfCheck.strict
    assert(verdict.necessary, 'study:necessary', '%s', verdict.why);
    assert(verdict.sufficiency, 'study:sufficiency', '%s', verdict.why);
else
    assert(verdict.necessary, 'study:necessary', '%s', verdict.why);
end

% ------------------------------------------------------------------------
function s = pf(c)
%% Purpose:
%
%   'PASS' if c, else 'FAIL'.
%
if c, s = 'PASS'; else, s = 'FAIL'; end
end

function [g, dgdy] = studyTerminal(y, xf, needJ)
%% Purpose:
%
%   The four terminal conditions x(t_f) = xf, and their Jacobian.
%
g = y(1:4) - xf;
dgdy = [];
if needJ, dgdy = [eye(4), zeros(4)]; end
end

function u = studyU(ppL, tt, x, p)
%% Purpose:
%
%   The PMP control at time tt: the costate from the converged arc, G from
%   the FLOWN state, u = -lam'G/2. No interpolation of u itself.
%
[~, G] = cartpole_field(x, p);
u = -(ppval(ppL, tt).'*G)/2;
end

function dx = studyFlyRhs(x, u, p)
%% Purpose:
%
%   The TRUE dynamics with a supplied control -- what the flown check flies.
%
[F, G] = cartpole_field(x, p);
dx = F + G*u;
end
```

Note on MATLAB scripts with local functions: they must appear at the END of the
script file, after all executable code, and the script may not be named the same
as any of them. If R2026a rejects a local function in a script that `run` is
called on from a function workspace, move the four helpers into
`ex3_cart_pole_pmp/private/` as their own files and delete them from the script —
but try the local-function form first, because keeping them in the file is what
makes the study readable top to bottom.

- [ ] **Step 4: Run the study alone, then the test**

```bash
cd /Users/msc/Desktop/optimal_control/optimal_control_examples/ex3_cart_pole_pmp
/Applications/MATLAB_R2026a.app/bin/matlab -batch "cd('$(pwd)'); cartpole_minenergy_study"
cd tests
/Applications/MATLAB_R2026a.app/bin/matlab -batch "cd('$(pwd)'); exit(~test_minenergy_study())"
```

Expected: every printed gate PASS, the verdict claiming, and the test green.
Record **every measured number** in the task report — `normR`, `Hrel`,
`missTerminal`, `missFlown`, `adjWorst`, `gapWorst`, `conjOut.nInterior`,
`conjOut.tested`, `J`, `lam0` — because sub-projects B and C quote them.

**If S2 reports an interior conjugate point or unresolved samples:** that is a
finding about this trajectory, not a bug. Report it, leave `strict` false, say so
in the report and in the script's README section. Do NOT adjust `tol.sigRatio` to
make it disappear.

**If the det-sign cross-check (V2b) disagrees:** stop and report. The likely
causes, in order, are the STM product's ordering (`it.PHI{k}` may be stored
first-to-last or last-to-first — check `ms_conjugate_test`'s own accumulation)
and the instrument's scaling convention. Read the instrument, match its
convention, and say in the report which one it uses.

- [ ] **Step 5: Add it to run_tests and commit**

Add `'test_minenergy_study'` to `ex3_cart_pole_pmp/run_tests.m`'s list. Update
`ex3_cart_pole_pmp/README.md`: a "Study script" section saying what the style is
FOR — the steps are visible rather than behind a front door, each condition is
gated inline, and section 8 asserts the script against the library so the two
cannot drift.

```bash
cd /Users/msc/Desktop/optimal_control
git add optimal_control_examples/ex3_cart_pole_pmp/cartpole_minenergy_study.m \
        optimal_control_examples/ex3_cart_pole_pmp/tests/test_minenergy_study.m \
        optimal_control_examples/ex3_cart_pole_pmp/run_tests.m \
        optimal_control_examples/ex3_cart_pole_pmp/README.md
git commit -m "cart-pole: the min-energy study script, first of the three entry scripts"
```

---

### Task 5: Documentation

**Files:**
- Modify: `optimal_control_examples/teaching_docs/direct_and_indirect_cartpole.tex`, `CLAUDE.md`

**Interfaces:** consumes everything above; produces no code.

- [ ] **Step 1: Update the teaching document**

`direct_and_indirect_cartpole.tex` currently covers one objective. Add a short
closing section in the same ELI5 / Intuition / Rigor tcolorbox form used
throughout, saying: (a) the study script is the executable companion to the
document, and which section of the script corresponds to which section of the
text; (b) that this objective has no switching structure because the control is
unbounded and the cost is smooth and strictly convex in `u`, which is exactly
what min-time and min-fuel will change; (c) that those two follow in
sub-projects B and C. Do not restate the mathematics of B and C — the spec has
it and the text will get it when the code exists.

Recompile twice and clean:

```bash
cd /Users/msc/Desktop/optimal_control/optimal_control_examples/teaching_docs
/Library/TeX/texbin/pdflatex -interaction=nonstopmode direct_and_indirect_cartpole.tex
/Library/TeX/texbin/pdflatex -interaction=nonstopmode direct_and_indirect_cartpole.tex
rm -f *.aux *.log *.out *.toc *.lof *.lot *.fls *.fdb_latexmk *.synctex.gz
```

- [ ] **Step 2: Update CLAUDE.md**

The `optimal_control_examples/` block of the directory tree gains
`cartpole_common/` (with a one-line description of what it holds and why),
`run_all_tests.m`, and `cartpole_minenergy_study.m`; add one line recording that
min-time and min-fuel are specced (naming the spec file) and pending. The
`oclib` mention gains `ms_conjugate_test`. Match the tree's existing density —
this is a map, not a changelog.

- [ ] **Step 3: Verify the whole tree and commit**

```bash
cd /Users/msc/Desktop/optimal_control/optimal_control_examples
/Applications/MATLAB_R2026a.app/bin/matlab -batch "cd('$(pwd)'); exit(~run_all_tests())"
grep -rn "orbit_transfer\|pumpkyn\|cr3bp" --include=*.m . && echo "CAMPAIGN LEAK -- stop" || echo "campaign-free: ok"
cd /Users/msc/Desktop/optimal_control
git add -A optimal_control_examples CLAUDE.md
git commit -m "cart-pole: docs for the shared plant and the min-energy study"
```

Expected: `run_all_tests` green; the grep finds nothing under
`optimal_control_examples` (the `oclib` path in `run_cartpole_pmp.m` is reached
by `fileparts`, not by naming a campaign).

---

## Notes for the executor

- **Run order matters.** Task 1 before Task 4 (the study reads the moved plant);
  Task 3 before Task 4 (the study calls `oc.ms_conjugate_test`). Task 2 can move
  but is easiest before Task 4, so the runner already covers the new folder when
  the study lands.
- **Task 3 touches library code five campaigns run on.** Its gates are the point
  of the task, not paperwork. A delegate that changes a number is not a delegate.
- **Do not widen a tolerance to make a check pass.** Every number in the study's
  section 0 came from a measurement recorded in this plan or in the ex3 task
  report. A gate that fails is a finding; measure why, then report it.
- **The examples tree stays campaign-free.** Task 5 Step 3's grep is the gate.
- **macOS has no `timeout`.** For anything long-running use `nohup ... &` and
  poll the log file.
