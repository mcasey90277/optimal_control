function ok = test_optimality_status()
% TEST_OPTIMALITY_STATUS  The classifier maps every outcome of the gate stack
% to its status (spec 4.2), infers legacy records conservatively (4.3), and
% the legend names the codes.
%
% INPUTS:  none
% OUTPUTS: ok [logical]  every check passed (prints PASS/FAIL per check)

ok = true;
here = fileparts(fileparts(mfilename('fullpath')));
addpath(here);
z = [(1:7).'; 4.1];
base = struct('ok', false, 'z', z, 'flyKm', 0.1, 'flyVms', 0.01, 'reason', '', ...
              'stage', 0, 'conjFound', false, 'hypAfterConj', '');
mk = @(varargin) setf(base, varargin{:});

% ---- new (stamped) results ------------------------------------------------
[c, r] = optimality_status(mk('ok', true, 'stage', 8, 'reason', 'certified'));
ok = chk(ok, c == 4 && strcmp(r, 'full stack passed'), 'certified -> 4 sufficient');
c = optimality_status(mk('z', nan(8, 1), 'reason', 'normal-chart polish did not converge'));
ok = chk(ok, c == -1, 'no converged root -> below the floor');
c = optimality_status(mk('stage', 1, 'flyKm', NaN, 'reason', 'flight inadmissible: hits the Moon'));
ok = chk(ok, c == -1, 'inadmissible flight (no flown miss) -> below the floor');
c = optimality_status(mk('stage', 1, 'flyKm', 150, 'reason', 'flown position miss 150.0 km > 100'));
ok = chk(ok, c == -1, 'misses the target -> below the floor');
[c, r] = optimality_status(mk('stage', 2, 'reason', 'Hamiltonian max|H| = 1.0e-03 > 1e-06'));
ok = chk(ok, c == 1 && contains(r, 'Hamiltonian'), 'pointwise PMP fails -> 1 neither, reason kept');
c = optimality_status(mk('stage', 3, 'reason', 'tfMin witness exceeded its 300 s cap'));
ok = chk(ok, c == 1, 'witness not confirmed -> 1 neither');
c = optimality_status(mk('stage', 4, 'reason', 'conjugate test verdict 0 (UNDETERMINED)'));
ok = chk(ok, c == 3, 'conjugate UNDETERMINED -> 3 necessary only');
c = optimality_status(mk('stage', 7, 'conjFound', true, 'hypAfterConj', 'held', 'reason', 'conjugate test verdict 0 (FAIL)'));
ok = chk(ok, c == 2, 'conjugate FAIL with gates + H6 held -> 2');
[c, r] = optimality_status(mk('stage', 4, 'conjFound', true, 'hypAfterConj', 'dim S = 2 (abnormal lift)', 'reason', 'conjugate test verdict 0 (FAIL)'));
ok = chk(ok, c == 3 && contains(r, 'hypotheses were not established') && contains(r, 'dim S'), ...
         'conjugate FAIL without its hypotheses -> 3, saying which failed');
c = optimality_status(mk('stage', 5, 'reason', 'H6 not established: margin 0.9'));
ok = chk(ok, c == 3, 'a sufficiency gate fails -> 3 necessary only');
c = optimality_status(mk('stage', 7, 'conjFound', true, 'hypAfterConj', 'held', ...
         'reason', 'dense conjugate scan not clear: 0 coarse sign change(s), 1 zero, 0 UNRESOLVED, 0 multiplicity (3 near-miss cleared)'));
ok = chk(ok, c == 2, 'dense scan finds a zero -> 2');
c = optimality_status(mk('stage', 7, 'reason', 'dense conjugate scan not clear: 0 coarse sign change(s), 0 zero, 1 UNRESOLVED, 0 multiplicity (3 near-miss cleared)'));
ok = chk(ok, c == 3, 'dense scan UNRESOLVED only -> 3');

% ---- legacy records (no .stage): inferred, conservative --------------------
L = rmfield(base, {'stage', 'conjFound', 'hypAfterConj'});
lg = @(varargin) setf(L, varargin{:});
[c, ~, inf_] = optimality_status(lg('ok', true, 'reason', 'certified'));
ok = chk(ok, c == 4 && inf_, 'legacy certified -> 4, inferred');
c = optimality_status(lg('reason', 'conjugate test verdict 0'));
ok = chk(ok, c == 3, 'legacy "verdict 0" is FAIL or UNDETERMINED -> 3 (re-certify to resolve)');
% C1 (final review): the legacy dense reason has no sign/floor split of its
% zeros, so only a coarse trusted sign change is a DEFINITE refutation
c = optimality_status(lg('reason', 'dense conjugate scan not clear: 0 coarse sign change(s), 1 zero, 0 UNRESOLVED, 0 multiplicity (3 near-miss cleared)'));
ok = chk(ok, c == 3, 'C1: legacy dense-scan zero without a coarse sign change -> 3 (may be a floor-level zero)');
c = optimality_status(lg('reason', 'dense conjugate scan not clear: 0 coarse sign change(s), 1 zero, 0 UNRESOLVED, 1 multiplicity (3 near-miss cleared)'));
ok = chk(ok, c == 3, 'C1: legacy multiplicity alone -> 3');
c = optimality_status(lg('reason', 'dense conjugate scan not clear: 1 coarse sign change(s), 1 zero, 0 UNRESOLVED, 0 multiplicity (3 near-miss cleared)'));
ok = chk(ok, c == 2, 'C1: legacy coarse trusted sign change -> 2 (the dense scan runs after gates + H6)');
c = optimality_status(lg('reason', 'independent-field Hamiltonian residual 2.0e-05 > 1e-06'));
ok = chk(ok, c == 3, 'legacy gate failure that mentions "Hamiltonian" -> 3, not 1');
c = optimality_status(lg('reason', 'minimum-principle |full gap| 1.0e-09 > 1e-12'));
ok = chk(ok, c == 1, 'legacy pointwise failure -> 1');
[c, r, inf_] = optimality_status(lg('reason', 'some wording no build ever used'));
ok = chk(ok, c == 1 && inf_ && contains(r, 'inferred'), 'REVIEW FOCUS 1: unrecognised legacy reason -> lowest recordable tier, inferred');

% ---- C2: field-harmonised legacy records (missing fields filled with []) ----
Hm = setf(L, 'stage', [], 'conjFound', [], 'hypAfterConj', [], 'status', [], 'status_reason', []);
hm = @(varargin) setf(Hm, varargin{:});
[c, ~, inf_] = optimality_status(hm('reason', 'conjugate test verdict 0'));
ok = chk(ok, c == 3 && inf_, 'C2: a harmonised legacy record (stage = []) is inferred, not read as stamped');
[c, ~, inf_] = optimality_status(hm('reason', 'minimum-principle |full gap| 1.0e-09 > 1e-12'));
ok = chk(ok, c == 1 && inf_, 'C2: a harmonised legacy pointwise failure -> 1 inferred (not 3 via an empty stage)');
[~, ~, inf_] = optimality_status(hm('stage', NaN, 'reason', 'conjugate test verdict 0'));
ok = chk(ok, inf_, 'C2: a non-finite stage is not a stamp');
Sh = rmfield(base, 'hypAfterConj');  Sh.stage = 7;  Sh.conjFound = true;  Sh.reason = 'conjugate test verdict 0 (FAIL)';
try, [c, r] = optimality_status(Sh); threw = false; catch, threw = true; c = NaN; r = ''; end
ok = chk(ok, ~threw && c == 3 && contains(r, 'not recorded'), 'C2: a stamped conjugate finding without .hypAfterConj is 3, not a crash');

% ---- the legend -------------------------------------------------------------
K = status_key();
ok = chk(ok, isequal(K.codes, [4 3 2 1 0 -1]) && numel(K.names) == 6 && strcmp(K.names{3}, 'conjugate point found'), ...
         'status_key names the six codes');

if ok, fprintf('test_optimality_status: ALL PASS\n'); else, fprintf('test_optimality_status: FAIL\n'); end
end

function s = setf(s, varargin)
% SETF  Set name/value pairs on a struct.  INPUTS: s; pairs.  OUTPUTS: s.
for k = 1:2:numel(varargin), s.(varargin{k}) = varargin{k+1}; end
end

function ok = chk(ok, c, msg)
% CHK  Print one PASS/FAIL line and fold it into ok.  INPUTS: ok; c; msg.
% OUTPUTS: ok.
if c, fprintf('  PASS  %s\n', msg); else, fprintf('  FAIL  %s\n', msg); end
ok = ok && c;
end
