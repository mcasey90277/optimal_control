function ok = test_export_deliverable_catalog()
% TEST_EXPORT_DELIVERABLE_CATALOG  The deliverable export keeps every datum
% (both content keys equal, schema valid), removes the internal records
% (merge, status_layer), rewrites alternative sources without file paths,
% and loads with a bare load.
%
% INPUTS:  none
% OUTPUTS: ok [logical]  every check passed (prints PASS/FAIL per check)

ok = true;
here = fileparts(fileparts(mfilename('fullpath')));            % IND
addpath(here, fullfile(fileparts(fileparts(here)), 'costate_common'));
recordMat = fullfile(here, 'results', 'library_70mN_24x48_v2', 'costate_catalog_dro_tulip_70mN.mat');
tmp = [tempname '.mat'];
[c, info] = export_deliverable_catalog(recordMat, tmp);
R = load(recordMat);  r = R.(char(fieldnames(R)));
ok = chk(ok, strcmp(catalog_content_key(c), catalog_content_key(r)) && strcmp(alternatives_content_key(c), alternatives_content_key(r)), 'both content keys equal the record''s');
ok = chk(ok, strcmp(info.primKey, catalog_content_key(r)) && strcmp(info.altKey, alternatives_content_key(r)), 'info carries the keys');
ok = chk(ok, isempty(catalog_schema('validate', c)), 'the export validates');
ok = chk(ok, ~isfield(c, 'merge') && ~isfield(c, 'status_layer') && isfield(c, 'build') && ischar(c.provenance), 'merge/status_layer removed; build + provenance text present');
ok = chk(ok, isequal(c.sheets, r.sheets) && isequal(c.status_key, r.status_key) && isequal(c.constants, r.constants), 'sheets, legend and constants bit-identical');
src = {c.alternatives.source};
ok = chk(ok, ~any(contains(src, '.mat')) && ~any(contains(src, 'reproduce_')) && ~any(contains(src, '/')), 'no file paths in any alternative source');
ok = chk(ok, any(strcmp(src, 'displaced primary (merge)')) || any(startsWith(src, 'rib refusal')) || any(startsWith(src, 'sheet candidate')), 'the kind words survive');
ok = chk(ok, c.build.nBorderline == nnz([c.alternatives.borderline]) && strcmp(c.build.auditVerdict, 'CLEAN'), 'build counts come from the data');
L = load(tmp);
ok = chk(ok, isfield(L, 'costate_catalog_dro_tulip_70mN') && isequal(L.costate_catalog_dro_tulip_70mN.sheets(1).z8, r.sheets(1).z8), 'saved under the catalog''s variable name; loads with load alone');
delete(tmp);
if ok, fprintf('test_export_deliverable_catalog: ALL PASS\n'); else, fprintf('test_export_deliverable_catalog: FAIL\n'); end
end

function ok = chk(ok, c, msg)
% CHK  Print one PASS/FAIL line and fold it into ok.  INPUTS: ok; c; msg.
% OUTPUTS: ok.
if c, fprintf('  PASS  %s\n', msg); else, fprintf('  FAIL  %s\n', msg); end
ok = ok && c;
end
