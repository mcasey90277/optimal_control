function ok = test_catalog_content_key()
% TEST_CATALOG_CONTENT_KEY  The key that binds an audit to the catalog it
% audited: it must change when a stored root changes (one costate, one flight
% time, one cell gained or lost) and must NOT change when only labels do
% (the family map is rewritten after the audit).
%
% INPUTS:  none
% OUTPUTS: ok [logical]  every check passed (prints PASS/FAIL per check)

ok = true;
addpath(fileparts(fileparts(mfilename('fullpath'))));
has = false(3, 4);  has([1 5 9]) = true;
tf = NaN(3, 4);  tf(has) = [6.1; 6.2; 6.3];
ei = zeros(3, 4);  ei(has) = 1:3;
s = struct('has_solution', has, 'tf_nd', tf, 'entry_index', ei, 'z8', reshape(1:24, 8, 3)/7, 'family_index', ones(3, 4));
c = struct('sheets', s, 'n_entries', 3);
k0 = catalog_content_key(c);
ok = chk(ok, ischar(k0) && numel(k0) == 32 && all(ismember(k0, '0123456789abcdef')), ['a 32-character hex key: ' k0]);
ok = chk(ok, strcmp(k0, catalog_content_key(c)), 'the same content gives the same key');
c2 = c;  c2.sheets.z8(3, 2) = c2.sheets.z8(3, 2)*(1 + 1e-15);
ok = chk(ok, ~strcmp(k0, catalog_content_key(c2)), 'one costate changed in its last bit -> another key');
c2 = c;  c2.sheets.tf_nd(1) = 6.1000001;
ok = chk(ok, ~strcmp(k0, catalog_content_key(c2)), 'one flight time changed -> another key');
c2 = c;  c2.sheets.has_solution(2) = true;
ok = chk(ok, ~strcmp(k0, catalog_content_key(c2)), 'one cell gained -> another key');
c2 = c;  c2.sheets.family_index(:) = 7;  c2.note = 'relabelled';
ok = chk(ok, strcmp(k0, catalog_content_key(c2)), 'labels rewritten after the audit do not change it');
c2 = c;  c2.sheets.tf_nd(2) = 99;                     % a value where there is no solution
ok = chk(ok, strcmp(k0, catalog_content_key(c2)), 'what sits in a cell with no solution does not count');
if ok, fprintf('test_catalog_content_key: ALL PASS\n'); else, fprintf('test_catalog_content_key: FAIL\n'); end
end

function ok = chk(ok, c, msg)
% CHK  Print one PASS/FAIL line and fold it into ok.  INPUTS: ok; c; msg.
% OUTPUTS: ok.
if c, fprintf('  PASS  %s\n', msg); else, fprintf('  FAIL  %s\n', msg); end
ok = ok && c;
end
