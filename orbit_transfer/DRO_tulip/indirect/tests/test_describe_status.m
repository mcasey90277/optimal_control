function ok = test_describe_status()
% TEST_DESCRIBE_STATUS  costate_lib_describe prints the status block on a
% catalog with the layer (counts from the data) and is byte-identical to
% its previous self on a catalog without it.
%
% INPUTS:  none
% OUTPUTS: ok [logical]  every check passed (prints PASS/FAIL per check)

ok = true;
here = fileparts(fileparts(mfilename('fullpath')));
addpath(here, fullfile(fileparts(fileparts(here)), 'costate_common'));
v2 = load(fullfile(here, 'results', 'library_70mN_24x48_v2', 'costate_catalog_dro_tulip_70mN.mat'));  c = v2.(char(fieldnames(v2)));
r24 = load(fullfile(here, 'results', 'library_70mN_24x24_final', 'costate_catalog_dro_tulip_70mN.mat'));  r = r24.(char(fieldnames(r24)));
out = evalc('D = costate_lib_describe(c);');
A = c.alternatives;
ok = chk(ok, contains(out, 'optimality status') && isequal(D.statusCounts(:, 2).', [nnz([A.status]==4) nnz([A.status]==3) nnz([A.status]==2) nnz([A.status]==1)]) && D.statusCounts(1, 1) == nnz(c.sheets(1).has_solution), 'status block printed; counts equal the data');
ok = chk(ok, D.nBorderline == nnz([A.borderline]) && D.K == size(c.sheets(1).junctions{1}, 2) && contains(out, 'borderline'), 'borderline count and K reported');
nJp = 0;  for ks = 1:numel(c.sheets), nJp = nJp + nnz(cellfun(@(x) ~isempty(x) && size(x, 1) == 14, c.sheets(ks).junctions)); end
nJa = nnz(arrayfun(@(q) ~isempty(q.junctions) && size(q.junctions, 1) == 14, A));
ok = chk(ok, isequal(D.nJunctions, [nJp nJa]) && contains(out, 'junctions on every entry') == (nJp == sum(cellfun(@numel, {c.sheets.junctions})) && nJa == numel(A)), 'nJunctions equals the directly computed counts; phrase follows them');
% byte-identity without the layer: compare with the pre-plan commit a62a65e
old = fullfile(tempdir, 'describe_old');  if ~isfolder(old), mkdir(old); end
[st, txt] = system(sprintf('cd "%s" && /Library/Developer/CommandLineTools/usr/bin/git show a62a65e:orbit_transfer/DRO_tulip/indirect/costate_lib_describe.m > "%s/costate_lib_describe.m"', fileparts(fileparts(fileparts(here))), old));
assert(st == 0, 'could not extract the committed file: %s', txt);
outNew = evalc('Dn = costate_lib_describe(r);');
cur = pwd;  cleanup = onCleanup(@() cd(cur));
cd(old);  outOld = evalc('Do = costate_lib_describe(r);');  cd(cur);
ok = chk(ok, strcmp(outNew, outOld) && isequal(rmfield_if(Dn, {'statusCounts','nAlternatives','nBorderline','K'}), Do), 'REVIEW FOCUS 1: without the layer the printout and D are unchanged');
ok = chk(ok, ~contains(outNew, 'optimality status') && ~isfield(Dn, 'statusCounts'), 'no status block and no new fields without the layer');
if ok, fprintf('test_describe_status: ALL PASS\n'); else, fprintf('test_describe_status: FAIL\n'); end
end

function s = rmfield_if(s, f)
% RMFIELD_IF  Remove fields that exist.  INPUTS: s; f cellstr.  OUTPUTS: s.
for k = 1:numel(f), if isfield(s, f{k}), s = rmfield(s, f{k}); end, end
end

function ok = chk(ok, c, msg)
% CHK  Print one PASS/FAIL line and fold it into ok.  INPUTS: ok; c; msg.
% OUTPUTS: ok.
if c, fprintf('  PASS  %s\n', msg); else, fprintf('  FAIL  %s\n', msg); end
ok = ok && c;
end
