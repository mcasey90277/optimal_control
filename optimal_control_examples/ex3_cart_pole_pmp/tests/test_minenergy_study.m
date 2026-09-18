function ok = test_minenergy_study()
%% Purpose:
%
%   The study script must RUN and must REACH A VERDICT. A study that prints
%   numbers and exits 0 regardless is a report, not a check.
%
%   Checks: it runs to the end without throwing (see WHY THE THROW IS THE
%   CHECK below); it leaves a verdict of the right shape; that
%   verdict.claim is wired the way the study computes it -- an IDENTITY
%   today (the study itself sets verdict.claim = necessary && sufficiency),
%   so this guards only against someone breaking that one line later, not
%   against an independently-measured overclaim; that the S2 conjugate-
%   point gate the study computes is actually enforced here (S.conj was
%   captured and never read by any chk() until this pass); that the study
%   and the front door build the SAME seed diagnostics
%   (signCorr/flipped/ampRatio), not only the same root -- lam0 agreeing is
%   invariant to how the seed was built, since both seeds land in one
%   basin, so it cannot by itself catch the two copies of the seed
%   computation (cartpole_minenergy_study.m vs run_cartpole_pmp.m) drifting
%   apart; and that the headline numbers (J, lam0) match run_cartpole_pmp
%   on the same problem -- the study and the front door must not drift
%   apart.
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
    % The actual check already happened: the try/catch above caught
    % nothing, i.e. the study ran to completion. This chk() is bookkeeping
    % that records the PASS in the same PASS/FAIL log as everything below,
    % not an independent boolean test -- cond is always true on this branch.
    ok = chk(ok, true, 'the study ran to the end without throwing (checked by the try/catch above)');
else
    ok = chk(ok, false, sprintf('the study runs to the end -- it THREW: %s', why));
    fprintf('TEST_MINENERGY_STUDY: FAILURE (see lines above)\n');
    return
end

ok = chk(ok, isstruct(S.verdict) && all(isfield(S.verdict, ...
             {'necessary', 'sufficiency', 'claim', 'why'})), ...
         'the script leaves a verdict with necessary / sufficiency / claim / why');
% IDENTITY, not an independent measurement: the study sets verdict.claim =
% verdict.necessary && verdict.sufficiency itself, so this can only fail if
% that one line is later edited to claim without checking. Kept as a real
% (if narrow) regression guard on that wiring.
ok = chk(ok, ~S.verdict.claim || (S.verdict.necessary && S.verdict.sufficiency), ...
         'verdict.claim is wired to necessary && sufficiency (regression guard on that line, not an overclaim probe)');

% S2 (the conjugate-point gate) is captured in S.conj but was never read by
% any chk() -- a regression in S2 would print FAIL in the study's own
% section 7 output and still leave this whole test PASS. Read it:
ok = chk(ok, S.conj.nInterior == 0 && S.conj.nUnresolved == 0, ...
         sprintf('S2 conjugate-point gate holds: nInterior %d, nUnresolved %d', ...
                 S.conj.nInterior, S.conj.nUnresolved));

out = run_cartpole_pmp(struct('plot', false));
ok = chk(ok, out.ok, sprintf('the front door itself is ok (%s)', out.why));
ok = chk(ok, abs(S.J - out.J)/out.J < 1e-8, ...
         sprintf('study and front door agree on J: %.6f vs %.6f', S.J, out.J));
ok = chk(ok, max(abs(S.lam0 - out.lam0)) < 1e-6, ...
         sprintf('and on lam0 (max abs diff %.2e)', max(abs(S.lam0 - out.lam0))));

% lam0 agreeing above pins the ROOT the BVP finds, which is invariant to how
% the seed that led to it was built -- any seed landing in the same basin
% gives the same lam0. It does NOT pin the two COPIES of the seed
% computation itself (cartpole_minenergy_study.m:171-196, lifted verbatim
% from run_cartpole_pmp.m:84-109). Compare the fields that actually
% distinguish a seed, to a tight tolerance:
ok = chk(ok, abs(S.seedDiag.signCorr - out.seed.signCorr) < 1e-10 && ...
             S.seedDiag.flipped == out.seed.flipped && ...
             abs(S.seedDiag.ampRatio - out.seed.ampRatio) < 1e-10, ...
         sprintf(['study and front door build the SAME seed: signCorr ' ...
                  '%.8f vs %.8f, flipped %d vs %d, ampRatio %.8f vs %.8f'], ...
                 S.seedDiag.signCorr, out.seed.signCorr, S.seedDiag.flipped, ...
                 out.seed.flipped, S.seedDiag.ampRatio, out.seed.ampRatio));

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
S = struct('verdict', verdict, 'J', J, 'lam0', lam0, 'conj', conjOut, ...
            'seedDiag', seedDiag);
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
