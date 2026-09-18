function ok = test_minenergy_study()
%% Purpose:
%
%   The study script must RUN and must REACH A VERDICT. A study that prints
%   numbers and exits 0 regardless is a report, not a check.
%
%   Checks: it runs to the end without throwing; it leaves a verdict of the
%   right shape; the verdict's claim is not made unless its own gates
%   support it (a verdict that claims more than it measured is the failure
%   mode this exists to catch); and the headline numbers match
%   run_cartpole_pmp on the same problem -- the study and the front door
%   must not drift apart.
%
%   WHY THE THROW IS THE CHECK, and why there is no separate
%   verdict.necessary check. The study ENDS with
%   assert(verdict.necessary, ...), so a failed necessary gate throws out of
%   run() and never returns a verdict at all: a chk() on S.verdict.necessary
%   after the call could not fail, and a check that cannot fail is worse than
%   no check because it reads as protection. The throw is caught here and
%   reported as the failure, with the assert's own message.
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

S   = struct();
why = '';
try
    S = run_study_isolated(fullfile(here, 'cartpole_minenergy_study.m'));
catch err
    why = sprintf('%s [%s]', err.message, err.identifier);
end
if isempty(why)
    ok = chk(ok, true, 'the study runs to the end and returns a verdict');
else
    ok = chk(ok, false, sprintf('the study runs to the end -- it THREW: %s', why));
    fprintf('TEST_MINENERGY_STUDY: FAILURE (see lines above)\n');
    return
end

ok = chk(ok, isstruct(S.verdict) && all(isfield(S.verdict, ...
             {'necessary', 'sufficiency', 'claim', 'why'})), ...
         'the script leaves a verdict with necessary / sufficiency / claim / why');
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
