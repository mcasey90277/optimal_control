function ok = test_relabel_borderline()
% TEST_RELABEL_BORDERLINE  relabel_borderline labels an alternative whose
% status flipped on re-audit CONSERVATIVELY: a BAD row whose root did not
% move takes min(stored, re-audit) with borderline = true and a reason
% naming both statuses, the audit's why and the old reason; a moved row and
% a row with no finite re-audit status are NOT relabelled (listed with why);
% OK rows and rows the audit never saw are untouched; every row carries a
% logical borderline field. Pure: synthetic catalog, no I/O.
%
% INPUTS:  none
% OUTPUTS: ok [logical]  every check passed (prints PASS/FAIL per check)

ok = true;
here = fileparts(fileparts(mfilename('fullpath')));
addpath(here);
row = @(st, rs) struct('sD', 0.25, 'sA', 0.5, 'iD', 1, 'iA', 1, 'z8', [0.1*ones(7, 1); 20], 'junctions', ones(14, 24), ...
                       'tf_nd', 20, 'status', st, 'status_reason', rs, 'inferred', false, 'source', 's', 'sheet', 1);
c = struct('alternatives', [row(4, 'full stack passed'), row(3, 'lift margin 9.1'), row(3, 'r3'), ...
                            row(4, 'r4'), row(2, 'r5'), row(3, 'r6'), row(3, 'never audited'), ...
                            row(4, 'r8'), row(5, 'r9'), row(3, 'r10')], ...
           'status_layer', struct('built', '2026-10-04'));
ar = @(k, okv, why, sNow, mv) struct('k', k, 'ok', okv, 'why', why, 'statusNow', sNow, 'moved', mv);
altRows = [ar(1, true, '', 4, false), ...
           ar(2, false, 'status not reproduced: stored 3, now 4 (lift)', 4, false), ...  % flipped UP: keep the lower, 3
           ar(3, false, 'status not reproduced: stored 3, now 1 (gate)', 1, false), ...  % flipped DOWN: take 1
           ar(4, false, 'the re-certification moved to another root', 4, true), ...      % moved: never relabelled
           ar(5, false, 'below the floor on re-certification: x', NaN, false), ...        % no status: not relabelled
           ar(6, false, 'status not reproduced: stored 3, now 0 (x)', 0, false), ...      % statusNow 0 < 1: not relabelled
           ar(8, false, 'below the floor on re-certification: polish', 3, false), ...     % FIX 1: finite status, no root: NOT relabelled
           ar(9, false, 'status not reproduced: stored 5, now 3 (x)', 3, false), ...      % FIX 2: stored 5 is not a code
           ar(10, false, 'status not reproduced: stored 3, now 2.5 (x)', 2.5, false)];    % FIX 2: 2.5 is not a code

[c2, info] = relabel_borderline(c, altRows);
A = c2.alternatives;
ok = chk(ok, all(arrayfun(@(a) isfield(a, 'borderline') && islogical(a.borderline) && isscalar(a.borderline), A)), ...
         'every alternative carries a logical scalar borderline field');
ok = chk(ok, A(2).status == 3 && A(2).borderline && contains(A(2).status_reason, 'borderline: stored 3, re-audit 4') ...
         && contains(A(2).status_reason, 'now 4 (lift)') && contains(A(2).status_reason, 'earlier: lift margin 9.1'), ...
         sprintf('BAD, not moved, re-audit HIGHER: status stays the lower (3), borderline, reason (%s)', A(2).status_reason));
ok = chk(ok, A(3).status == 1 && A(3).borderline && contains(A(3).status_reason, 'borderline: stored 3, re-audit 1') ...
         && contains(A(3).status_reason, 'earlier: r3'), 'BAD, not moved, re-audit LOWER: relabelled to the lower (1)');
ok = chk(ok, A(4).status == 4 && ~A(4).borderline && strcmp(A(4).status_reason, 'r4'), 'a MOVED row is not relabelled');
ok = chk(ok, A(5).status == 2 && ~A(5).borderline && A(6).status == 3 && ~A(6).borderline, ...
         'a row with no finite re-audit status >= 1 is not relabelled');
ok = chk(ok, A(1).status == 4 && ~A(1).borderline && strcmp(A(1).status_reason, 'full stack passed') ...
         && A(7).status == 3 && ~A(7).borderline, 'OK rows and unaudited rows are untouched');
ok = chk(ok, A(8).status == 4 && ~A(8).borderline && strcmp(A(8).status_reason, 'r8') ...
         && contains(info.notRelabelled([info.notRelabelled.k] == 8).why, 'below the floor'), ...
         'FIX 1: a below-the-floor row with a finite statusNow is NOT relabelled (listed with why)');
ok = chk(ok, A(9).status == 5 && ~A(9).borderline && contains(info.notRelabelled([info.notRelabelled.k] == 9).why, 'stored status') ...
         && A(10).status == 3 && ~A(10).borderline && contains(info.notRelabelled([info.notRelabelled.k] == 10).why, 're-audit status'), ...
         'FIX 2: a stored status or a re-audit status that is not an integer code 1..4 is not relabelled');
ok = chk(ok, isequal(sort([info.notRelabelled.k]), [4 5 6 8 9 10]) && all(~cellfun(@isempty, {info.notRelabelled.why})) ...
         && contains(info.notRelabelled([info.notRelabelled.k] == 4).why, 'moved'), ...
         'moved and non-finite rows are listed in info.notRelabelled with why');
ok = chk(ok, isequal(info.rows, [2 3 4 3; 3 3 1 1]), 'info.rows lists (k, stored, now, newStatus)');
ok = chk(ok, c2.status_layer.nBorderline == 2 && contains(c2.status_layer.borderlineNote, 'lower bound') ...
         && strcmp(c2.status_layer.built, '2026-10-04'), 'status_layer counts the borderline rows and notes the lower bound');
c3 = relabel_borderline(c, altRows([1 4]));
ok = chk(ok, isequal([c3.alternatives.status], [c.alternatives.status]) && ~any([c3.alternatives.borderline]) ...
         && c3.status_layer.nBorderline == 0, 'nothing to relabel: statuses unchanged, every row borderline false');
ok = chk(ok, throws(@() relabel_borderline(c, ar(99, false, 'x', 1, false))), 'a row index outside the table is refused');
if ok, fprintf('test_relabel_borderline: ALL PASS\n'); else, fprintf('test_relabel_borderline: FAIL\n'); end
end

function t = throws(f)
% THROWS  True if f() throws.  INPUTS: f (handle).  OUTPUTS: t.
try
    f();  t = false;
catch
    t = true;
end
end

function ok = chk(ok, c, msg)
% CHK  Print one PASS/FAIL line and fold it into ok.  INPUTS: ok; c; msg.
% OUTPUTS: ok.
if c, fprintf('  PASS  %s\n', msg); else, fprintf('  FAIL  %s\n', msg); end
ok = ok && c;
end
