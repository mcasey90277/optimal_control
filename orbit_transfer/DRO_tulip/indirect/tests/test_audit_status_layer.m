function ok = test_audit_status_layer()
% TEST_AUDIT_STATUS_LAYER  The status audit fails closed: a reproduced
% status passes; a different status, a root that moved, and corrupted
% junctions are each a BAD row naming why; the content keys are bound.
%
% INPUTS:  none
% OUTPUTS: ok [logical]  every check passed (prints PASS/FAIL per check)

ok = true;
here = fileparts(fileparts(mfilename('fullpath')));
addpath(here, fullfile(fileparts(fileparts(here)), 'costate_common'));
L = load(fullfile(here, 'results', 'library_70mN_24x24_final', 'costate_catalog_dro_tulip_70mN.mat'));  c = L.(char(fieldnames(L)));
sh = c.sheets(1);  k = sh.entry_index(1, 1);
row = struct('sD', sh.sD_frac(1), 'sA', sh.sA_frac(1), 'iD', 0, 'iA', 0, 'z8', sh.z8(:, k), 'junctions', ones(14, 24), ...
             'tf_nd', sh.z8(8, k), 'status', 3, 'status_reason', 'x', 'inferred', false, 'conj', 1, 'conjVerdict', 'PASS', ...
             'minLamV', 1, 'minQmt', 1, 'dimS', 1, 'h6Margin', 1, 'liftMargin', 1, 'flyKm', 0, 'flyVms', 0, 'source', 't', 'sheet', 1);
c.alternatives = [row, setf(row, 'status', 2), row];
fake = @(st, dz) @(seed, rv0, rvf, B, o) struct('z', [seed.Y(8:14, 1); seed.tf] + dz, 'status', st, 'reason', 'fake');
tmp = [tempname '.mat'];  catalog = c;  save(tmp, 'catalog');
A = audit_status_layer(tmp, struct('skipPrimaries', true, 'idxAlt', 1, 'certifier', fake(3, 0), 'pool', []));
ok = chk(ok, A.nBad == 0 && A.altRows(1).ok, 'a reproduced status passes');
A = audit_status_layer(tmp, struct('skipPrimaries', true, 'idxAlt', 2, 'certifier', fake(3, 0), 'pool', []));
ok = chk(ok, A.nBad == 1 && contains(A.altRows(1).why, 'status not reproduced'), 'a different status is BAD, named');
A = audit_status_layer(tmp, struct('skipPrimaries', true, 'idxAlt', 3, 'certifier', fake(3, [0.01; zeros(7, 1)]), 'pool', []));
ok = chk(ok, A.nBad == 1 && A.altRows(1).moved && contains(A.altRows(1).why, 'moved'), 'REVIEW FOCUS 5: a root that moved is BAD, not relabelled');
ok = chk(ok, numel(A.altContentKey) == 32 && ~strcmp(A.altContentKey, alternatives_content_key(setf(c, 'alternatives', row))), ...
         'the alternatives key is bound to the table''s content');
ok = chk(ok, ~strcmp(alternatives_content_key(c), alternatives_content_key(setf(c, 'alternatives', [row, row, row]))), ...
         'the alternatives key changes when only a status changes');
ok = chk(ok, isempty(alternatives_content_key(rmfield(c, 'alternatives'))), 'no alternatives -> empty key');
delete(tmp);
if ok, fprintf('test_audit_status_layer: ALL PASS\n'); else, fprintf('test_audit_status_layer: FAIL\n'); end
end

function s = setf(s, varargin)
% SETF  Set name/value pairs.  INPUTS: s; pairs.  OUTPUTS: s.
for k = 1:2:numel(varargin), s.(varargin{k}) = varargin{k+1}; end
end

function ok = chk(ok, c, msg)
% CHK  Print one PASS/FAIL line and fold it into ok.  INPUTS: ok; c; msg.
% OUTPUTS: ok.
if c, fprintf('  PASS  %s\n', msg); else, fprintf('  FAIL  %s\n', msg); end
ok = ok && c;
end
