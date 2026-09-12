function ok = test_conj_coverage()
%% Purpose:
%
%   Tests the three things the 2026-09-11 reviews found the conjugate
%   instrument (ms_conjugate_test) did not say or do:
%
%     1. COVERAGE. A free-time test without y(t_f) stops at t_K and used to
%        return PASS with the final segment unmonitored. It must now return
%        UNDETERMINED with .covered false, and a covered run must say so.
%     2. THE LAST BRACKET. Opposite nonzero signs at t_K and t_f put the
%        root strictly inside (t_K, t_f): a conjugate point, hence FAIL --
%        when the sign at t_f is trustworthy. The old rule called every
%        last-bracket crossing ENDPOINT. With the last determinant declared
%        unresolvable (resolvedTol = Inf) the verdict must fall back to
%        ENDPOINT.
%     3. THE TWO EXACT IDENTITIES of the quotiented form, J p(0) = 0 and
%        p(t)' J(t) = 0, are measured and must hold to integration accuracy
%        on a certified entry.
%
%   The fixture is the certified anchor, re-solved by ms_tfmin with the
%   conjugate test on; the last-bracket case is manufactured by negating the
%   final segment STM, which flips the sign of the final determinant and of
%   nothing else.
%
%% Inputs:  none
%% Outputs: ok [logical]
%
%% Revision History:
%  M. Casey                                                   (c) 09/11/2026
%  Copyright Coorbital Inc.
%% ------------------------ Begin Code Sequence ---------------------------

ok = true;
here = fileparts(fileparts(mfilename('fullpath')));
addpath(here, fullfile(fileparts(here), 'DRO_tulip', 'indirect'));
lib = dro_tulip_library();
E = lib(strcmp({lib.src}, 'anchor'));  E = E(1);
[B, ~] = arclength_arrival('setup');
rv0 = B.stateD(E.sD);  rvf = B.stateA(E.sA);  z = E.z(:);
K = 24;
seed = seed_from_z8(z, rv0(1:6), K, B.Tnd, B.cnd, B.mu);
[~, it] = ms_tfmin(rv0(1:6), rvf(1:6), seed, B.Tnd, B.cnd, B.mu, ...
                   struct('tolR', 3e-11, 'wallSec', 600, 'conjTest', true));
assert(it.converged, 'fixture: the anchor must re-converge');
cj = it.conj;

% ---- 1. coverage ---------------------------------------------------------
ok = chk(ok, isfield(cj, 'covered') && cj.covered && abs(cj.sampledThrough - it.tGrid(end)) < 1e-12, ...
         'a run with y(t_f) reports covered through t_f');
ok = chk(ok, strcmp(cj.verdict, 'PASS') && isfield(cj, 'reason') && ~isempty(cj.reason), ...
         sprintf('the anchor passes, with a reason: %s', cj.reason));
info2 = it;  info2 = rmfield(info2, 'Yend');
cj2 = ms_conjugate_test(info2, struct('flow', @(y) flow6(y, B.Tnd, B.cnd, B.mu)));
ok = chk(ok, ~cj2.covered && strcmp(cj2.verdict, 'UNDETERMINED') && ~cj2.pass, ...
         sprintf('without y(t_f): %s (%s), pass = %d', cj2.verdict, cj2.reason, cj2.pass));
ok = chk(ok, numel(cj2.detScaled) == K - 1 && cj2.sampledThrough < it.tGrid(end), ...
         'and it sampled only through t_K');

% ---- 2. the last bracket ---------------------------------------------------
info3 = it;  info3.PHI{end} = -info3.PHI{end};
cj3 = ms_conjugate_test(info3, struct('flow', @(y) flow6(y, B.Tnd, B.cnd, B.mu)));
ok = chk(ok, sign(cj3.detScaled(end)) == -sign(cj.detScaled(end)) && ...
             all(sign(cj3.detScaled(1:end-1)) == sign(cj.detScaled(1:end-1))), ...
         'negating the last STM flips only the final determinant');
ok = chk(ok, strcmp(cj3.verdict, 'FAIL') && cj3.nInterior == 1 && cj3.nEndResolved == 1 && ~cj3.atFinal, ...
         sprintf('a resolved last-bracket crossing is an INTERIOR root: %s (%s)', cj3.verdict, cj3.reason));
cj4 = ms_conjugate_test(info3, struct('flow', @(y) flow6(y, B.Tnd, B.cnd, B.mu), 'resolvedTol', Inf));
ok = chk(ok, strcmp(cj4.verdict, 'ENDPOINT') && cj4.atFinal && cj4.nInterior == 0, ...
         sprintf('with the final sign declared untrustworthy it stays ENDPOINT: %s', cj4.verdict));
ok = chk(ok, cj.sigRatio(end) > 1e-10, ...
         sprintf('the anchor''s final block IS resolved (sigma ratio %.1e > 1e-10)', cj.sigRatio(end)));

% ---- 3. the identities -----------------------------------------------------
ok = chk(ok, isfield(cj, 'kernelRight') && cj.kernelRight < 1e-8, ...
         sprintf('scaling kernel J p(0) = 0 to %.1e', cj.kernelRight));
ok = chk(ok, isfield(cj, 'kernelLeft') && cj.kernelLeft < 1e-8, ...
         sprintf('left kernel p(t)''J(t) = 0 to %.1e', cj.kernelLeft));
ok = chk(ok, isfield(cj, 'tFirstFullRank') && cj.tFirstFullRank > 0 && cj.tFirstFullRank <= it.tGrid(2) + 1e-12, ...
         sprintf('first full-rank junction at t/t_f = %.3f', cj.tFirstFullRank/it.tGrid(end)));

if ok, fprintf('TEST_CONJ_COVERAGE: ALL PASS\n'); else, fprintf('TEST_CONJ_COVERAGE: FAIL\n'); end
end

function f = flow6(y, Tmax, c, muStar)
% FLOW6  Position/velocity rows of the min-time dynamics at a point.
% INPUTS: y [14x1]; Tmax; c; muStar.  OUTPUTS: f [6x1].
F = mintime_rhs_point(y, Tmax, c, muStar);
f = F(1:6);
end

function ok = chk(ok, cond, msg)
% CHK  Print one check.  INPUTS: ok; cond; msg.  OUTPUTS: ok.
if cond, fprintf('  PASS  %s\n', msg); else, fprintf('  FAIL  %s\n', msg); ok = false; end
end
