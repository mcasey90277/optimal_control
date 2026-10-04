function ok = test_certify_status()
% TEST_CERTIFY_STATUS  certify_root stamps .status on every result; a coarse
% conjugate FAIL still runs the hypothesis gates and H6 (for the record) and
% is status 2 when they hold, 3 when they do not; UNDETERMINED is 3; a
% pointwise failure is 1; a dense-scan zero is 2; the untouched anchor is 4.
% Uses the test seam (.override), so every forced result is .overridden.
%
% INPUTS:  none
% OUTPUTS: ok [logical]  every check passed (prints PASS/FAIL per check)

ok = true;
here = fileparts(fileparts(mfilename('fullpath')));                     % IND
cc = fullfile(fileparts(fileparts(here)), 'costate_common');
addpath(here, cc);
lib = dro_tulip_library();
E = lib(strcmp({lib.src}, 'anchor'));  E = E(1);
[B, ~] = arclength_arrival('setup');
rv0 = B.stateD(E.sD);  rvf = B.stateA(E.sA);
seed = seed_from_z8(E.z(:), rv0(1:6), 24, B.Tnd, B.cnd, B.mu);
base = struct('sA', E.sA, 'sD', E.sD, 'pool', capped_pool());
run = @(o) certify_root(seed, rv0(1:6), rvf(1:6), B, o);

C = run(base);
ok = chk(ok, C.ok && C.status == 4 && C.stage == 8 && strcmp(C.status_reason, 'full stack passed') && ~C.overridden, ...
         sprintf('anchor: status 4, stage 8 (%s)', C.reason));
ok = chk(ok, strcmp(C.conjVerdict, 'PASS'), 'the coarse verdict string is stored');

o = base;  o.override = struct('conj', struct('pass', false, 'verdict', 'FAIL'));
C = run(o);
ok = chk(ok, ~C.ok && C.conjFound && strcmp(C.hypAfterConj, 'held') && C.stage == 7 && C.status == 2 && C.overridden, ...
         sprintf('coarse FAIL, gates + H6 held -> 2 (stage %d, %s)', C.stage, C.status_reason));
ok = chk(ok, startsWith(C.reason, 'conjugate test verdict 0 (FAIL)') && ~isempty(C.g) && isfinite(C.h6Margin) && isempty(C.conjDense), ...
         'the reason is the conjugate one; gates and H6 were computed; the dense scan was not run');

o.override.g = struct('dimS', 2);
C = run(o);
ok = chk(ok, C.status == 3 && contains(C.hypAfterConj, 'dim S') && startsWith(C.reason, 'conjugate test verdict 0 (FAIL)'), ...
         sprintf('coarse FAIL, a gate fails -> 3 saying which (%s)', C.status_reason));

o = base;  o.override = struct('conj', struct('pass', false, 'verdict', 'UNDETERMINED'));
C = run(o);
ok = chk(ok, C.status == 3 && ~C.conjFound && C.stage == 4 && isempty(C.g), 'UNDETERMINED -> 3, gates not run');

o = base;  o.override = struct('PW', struct('Hmax', 1));
C = run(o);
ok = chk(ok, C.status == 1 && C.stage == 2, sprintf('pointwise failure -> 1 (%s)', C.reason));

o = base;  o.override = struct('CS', struct('nZero', 1, 'clear', false));
C = run(o);
ok = chk(ok, C.status == 2 && C.conjFound && C.stage == 7, sprintf('dense-scan zero -> 2 (%s)', C.reason));

if ok, fprintf('test_certify_status: ALL PASS\n'); else, fprintf('test_certify_status: FAIL\n'); end
end

function ok = chk(ok, c, msg)
% CHK  Print one PASS/FAIL line and fold it into ok.  INPUTS: ok; c; msg.
% OUTPUTS: ok.
if c, fprintf('  PASS  %s\n', msg); else, fprintf('  FAIL  %s\n', msg); end
ok = ok && c;
end
