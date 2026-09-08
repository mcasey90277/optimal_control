function ok = test_regime_verdicts()
%% Purpose:
%
%   Tests regime_verdicts -- the per-case classifier of the REGIME MAP.
%
%   The procedural question ("did the p-walk reach the schedule's last
%   rung?") and the physical question ("did the walk find the optimum?")
%   give DIFFERENT answers, and only the second supports a claim about a
%   family. Measured 2026-09-07: eps stops at p = 0.0016 on (2,5) gamma 1.1
%   yet its m_f agrees with the two arms that reached 0.001 to 2.4e-6 --
%   the last rungs were unnecessary, not a failure. On (1,2) gamma 1.223
%   huber and huberc stop with m_f low by 9.4e-3 (1.4 kg of 150) -- that is
%   a failure. The classifier must separate those.
%
%% Inputs:  none
%% Outputs: ok [logical]
%
%% Revision History:
%  M. Casey                                                   (c) 09/07/2026
%  Copyright Coorbital Inc.
%% ------------------------ Begin Code Sequence ---------------------------

ok = true;
addpath(fileparts(fileparts(mfilename('fullpath'))));

% one case, three arms: eps stops early but AGREES on mass -> all solved
F = [mk('eps', [2 5], 1.1, 0.001814, 0.941107, 9), ...
     mk('huber', [2 5], 1.1, 0.001, 0.941104, 9), ...
     mk('huberc', [2 5], 1.1, 0.001, 0.941104, 9)];
V = regime_verdicts(F);
ok = chk(ok, numel(V) == 1 && numel(V.solved) == 3, ...
         sprintf('mass-converged case: %d of 3 families solved it', numel(V.solved)));
ok = chk(ok, strcmp(V.verdict, 'all three'), sprintf('verdict = %s', V.verdict));
ok = chk(ok, max(V.dMf) < 1e-5, sprintf('max dm_f = %.1e', max(V.dMf)));

% one case where two arms stop 9e-3 light -> genuinely eps-only
F2 = [mk('eps', [1 2], 1.223, 0.001311, 0.941482, 9), ...
      mk('huber', [1 2], 1.223, 0.8182, 0.932009, 2), ...
      mk('huberc', [1 2], 1.223, 0.8091, 0.932144, 4)];
V2 = regime_verdicts(F2);
ok = chk(ok, isequal(V2.solved, {'eps'}), ...
         sprintf('one-family case: solved by %s', strjoin(V2.solved, ',')));
ok = chk(ok, strcmp(V2.verdict, 'eps ONLY'), sprintf('verdict = %s', V2.verdict));
ok = chk(ok, V2.nSwitchSolved == 9 && max(V2.nSwitchFailed) == 4, ...
         sprintf('switch counts: solved %d, failed max %d', V2.nSwitchSolved, max(V2.nSwitchFailed)));
ok = chk(ok, V2.structureChange, 'failures sit at a DIFFERENT switch count than the solution');

% a case where the loser matches the structure but stops light: not a
% structure change (guards against reading every failure as H1)
F3 = [mk('eps', [6 8], 1.1, 0.001242, 0.928082, 7), ...
      mk('huber', [6 8], 1.1, 0.001, 0.928081, 7), ...
      mk('huberc', [6 8], 1.1, 0.001, 0.928081, 7)];
V3 = regime_verdicts(F3);
ok = chk(ok, ~V3.structureChange, 'all-agree case carries no structure-change flag');

% incomplete case -> reported, never scored
V4 = regime_verdicts(F2(1:2));
ok = chk(ok, strcmp(V4.verdict, 'INCOMPLETE (2/3)'), sprintf('verdict = %s', V4.verdict));

if ok, fprintf('TEST_REGIME_VERDICTS: ALL PASS\n');
else,  fprintf('TEST_REGIME_VERDICTS: FAILURE (see lines above)\n');
end
end

function f = mk(fam, c, g, p, mf, nX)
% MK  Minimal feature struct for the classifier.  INPUTS: fam; c; g; p; mf;
% nX.  OUTPUTS: f.
f = struct('family', fam, 'cellIdx', c, 'gamma', g, 'pFloor', p, 'mf', mf, ...
           'nCross', nX, 'outcome', 'wall', 'rampWidth', 2*p, 'nFail', 0, ...
           'wallMin', 1, 'minAbsDQdt', 1);
end

function ok = chk(ok, cond, label)
% CHK  Accumulate a labeled pass/fail.  INPUTS: ok; cond; label. OUTPUTS: ok.
if cond, fprintf('  PASS  %s\n', label);
else,    fprintf('  FAIL  %s\n', label);  ok = false;
end
end
