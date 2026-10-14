function ok = test_export_deliverable_catalog()
% TEST_EXPORT_DELIVERABLE_CATALOG  The deliverable export keeps every datum
% (both content keys equal, schema valid), removes the internal records
% (merge, status_layer), rewrites alternative sources without file paths,
% and loads with a bare load. cleanSource maps the brief's examples (and the
% 96-phase sheet keeps its kind); with no audit opts the build is stamped
% 'not audited at export', with opts it carries the caller's verdict.
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
ok = chk(ok, c.build.nBorderline == nnz([c.alternatives.borderline]), 'build counts come from the data');
ok = chk(ok, strcmp(c.build.auditVerdict, 'not audited at export') && isempty(c.build.note), ...
         'a bare export does not claim an audit (verdict ''not audited at export'', empty note)');
hits = pathHits(c, 'c');
ok = chk(ok, isempty(hits), sprintf('no char anywhere in the export holds /Users/ or /Volumes/ (%d found%s)', ...
         numel(hits), strjoin(strcat({' '}, hits(1:min(3, end))), '')));
s1 = c.sheets(1);  vv = s1.conj_pass(logical(s1.has_solution));
ok = chk(ok, c.conj_test.nPass + c.conj_test.nFail == nnz(vv == 1 | vv == 0) && contains(c.conj_test.note, 'recomputed'), ...
         sprintf('conj_test nPass + nFail (%d + %d) = cells with a verdict (%d)', c.conj_test.nPass, c.conj_test.nFail, nnz(vv == 1 | vv == 0)));
tmp2 = [tempname '.mat'];
c2 = export_deliverable_catalog(recordMat, tmp2, struct('auditVerdict', 'CLEAN', 'auditNote', '7 ok / 0 bad'));
delete(tmp2);
ok = chk(ok, strcmp(c2.build.auditVerdict, 'CLEAN') && strcmp(c2.build.note, '7 ok / 0 bad') && contains(c2.provenance, 'CLEAN'), ...
         'opts.auditVerdict / opts.auditNote fill build.auditVerdict / build.note and the provenance verdict');
tbl = {'reproduce_70mN_24x48/round_01/arrival_sheet_70mN_nA48 cand, column 7 (re-certified)', 'sheet candidate, column 7 (re-certified)'
       'rib refusal at sA 0.2837 (walked 0)',                                          'rib refusal at sA 0.2837 (walked 0)'
       'displaced primary (merge)',                                                    'displaced primary (merge)'
       'sheet96_resolution_test/arrival_sheet_70mN_nA96 cand, column 61 (re-certified)', '96-phase sheet candidate, column 61 (re-certified)'};
got = export_deliverable_catalog('cleanSource', tbl(:, 1));
for k = 1:size(tbl, 1)
    ok = chk(ok, strcmp(got{k}, tbl{k, 2}), sprintf('cleanSource: ''%s'' -> ''%s''', tbl{k, 1}, got{k}));
end
ok = chk(ok, nnz(startsWith(src, '96-phase sheet candidate')) == nnz(contains({r.alternatives.source}, 'sheet96_resolution_test')), ...
         'every 96-phase sheet source keeps its kind in the export');
L = load(tmp);
ok = chk(ok, isfield(L, 'costate_catalog_dro_tulip_70mN') && isequal(L.costate_catalog_dro_tulip_70mN.sheets(1).z8, r.sheets(1).z8), 'saved under the catalog''s variable name; loads with load alone');
delete(tmp);
if ok, fprintf('test_export_deliverable_catalog: ALL PASS\n'); else, fprintf('test_export_deliverable_catalog: FAIL\n'); end
end

function h = pathHits(x, nm)
% PATHHITS  Names of every char/string inside x (recursive through structs
% and cells) that contains '/Users/' or '/Volumes/'.
% INPUTS: x (any); nm (char, its name).  OUTPUTS: h {cellstr}.
h = {};
if ischar(x) || isstring(x)
    if any(contains(x, {'/Users/', '/Volumes/'})), h = {nm}; end
elseif iscell(x)
    parts = cell(1, numel(x));
    for k = 1:numel(x), parts{k} = pathHits(x{k}, sprintf('%s{%d}', nm, k)); end
    h = [h, parts{:}];
elseif isstruct(x)
    fn = fieldnames(x);  parts = cell(numel(x), numel(fn));
    for e = 1:numel(x)
        for q = 1:numel(fn), parts{e, q} = pathHits(x(e).(fn{q}), sprintf('%s(%d).%s', nm, e, fn{q})); end
    end
    h = [h, parts{:}];
end
end

function ok = chk(ok, c, msg)
% CHK  Print one PASS/FAIL line and fold it into ok.  INPUTS: ok; c; msg.
% OUTPUTS: ok.
if c, fprintf('  PASS  %s\n', msg); else, fprintf('  FAIL  %s\n', msg); end
ok = ok && c;
end
