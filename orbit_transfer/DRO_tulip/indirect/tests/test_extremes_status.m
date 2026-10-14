function ok = test_extremes_status()
% TEST_EXTREMES_STATUS  costate_catalog_extremes: defaults unchanged on the
% old fields, alternatives and status filters, graceful on a catalog
% without the status layer.
%
% INPUTS:  none
% OUTPUTS: ok [logical]  every check passed (prints PASS/FAIL per check)

ok = true;
here = fileparts(fileparts(mfilename('fullpath')));
addpath(here, fullfile(fileparts(fileparts(here)), 'costate_common'));
v2 = load(fullfile(here, 'results', 'library_70mN_24x48_v2', 'costate_catalog_dro_tulip_70mN.mat'));  c = v2.(char(fieldnames(v2)));
r24 = load(fullfile(here, 'results', 'library_70mN_24x24_final', 'costate_catalog_dro_tulip_70mN.mat'));  r = r24.(char(fieldnames(r24)));
newF = {'kind','status','statusName','reason','junctions','altIndex'};

% defaults vs the pre-plan commit a62a65e (R3)
old = fullfile(tempdir, 'extremes_old');  if ~isfolder(old), mkdir(old); end
[st, txt] = system(sprintf('cd "%s" && /Library/Developer/CommandLineTools/usr/bin/git show a62a65e:orbit_transfer/DRO_tulip/indirect/costate_catalog_extremes.m > "%s/costate_catalog_extremes.m"', fileparts(fileparts(fileparts(here))), old));
assert(st == 0, 'could not extract the committed file: %s', txt);
cur = pwd;  cleanup = onCleanup(@() cd(cur));
mets = {'deltaV', 'time'};
for km = 1:2
    evalc('[nMin, nMax] = costate_catalog_extremes(c, struct(), false, mets{km});');
    cd(old);  evalc('[oMin, oMax] = costate_catalog_extremes(c, struct(), false, mets{km});');  cd(cur);
    ok = chk(ok, isequal(rmfield_if(nMin, newF), oMin) && isequal(rmfield_if(nMax, newF), oMax), ['REVIEW FOCUS 1: defaults identical on the old fields, metric ' mets{km}]);
end
ok = chk(ok, strcmp(nMin.kind, 'primary') && nMin.altIndex == 0 && nMin.status == 4, 'default extremes are primaries with their status');

% alternatives included
evalc('[pMin, pMax] = costate_catalog_extremes(c, struct(), false, ''time'');');
evalc('[aMin, aMax] = costate_catalog_extremes(c, struct(''includeAlternatives'', true), false, ''time'');');
ok = chk(ok, aMin.tf_days <= pMin.tf_days && aMax.tf_days >= pMax.tf_days && ismember(aMin.kind, {'primary','alternative'}), 'alternatives widen (never narrow) the tf range; kind set');
A = c.alternatives;
ok = chk(ok, aMin.tf_nd >= min([A.tf_nd]) - 1e-12 || strcmp(aMin.kind, 'primary'), 'extreme minimum is a real entry');

% status filter with alternatives
evalc('[sMin, sMax] = costate_catalog_extremes(c, struct(''status'', 2, ''includeAlternatives'', true), false, ''time'');');
ok = chk(ok, sMin.status == 2 && sMax.status == 2 && ~isempty(sMin.reason) && ~isempty(sMax.reason), 'status = 2: both extremes have status 2 and a reason');
ok = chk(ok, strcmp(sMin.kind, 'alternative') && sMin.altIndex >= 1 && sMin.tf_nd == min([A([A.status] == 2).tf_nd]) && sMax.tf_nd == max([A([A.status] == 2).tf_nd]), 'status 2 extremes equal the direct min/max over those alternatives');
txt = evalc('costate_catalog_extremes(c, struct(''status'', 2, ''includeAlternatives'', true), false, ''time'');');
ok = chk(ok, contains(txt, 'status: conjugate point found --') && contains(txt, 'status=2'), 'report prints the status line and the filter');
% status filter without alternatives: primaries only (all 4 here)
evalc('[qMin, ~] = costate_catalog_extremes(c, struct(''status'', 4), false, ''time'');');
ok = chk(ok, isequal(rmfield_if(qMin, newF), rmfield_if(pMin, newF)) && strcmp(qMin.kind, 'primary'), 'status 4 on primaries only equals the default minimum');
threw = false;
try, evalc('costate_catalog_extremes(c, struct(''status'', 2), false, ''time'');'); catch, threw = true; end
ok = chk(ok, threw, 'status 2 with no alternatives matches nothing and errors by assertion');

% catalog without the layer
evalc('[rMin, rMax] = costate_catalog_extremes(r, struct(''includeAlternatives'', true), false);');
ok = chk(ok, strcmp(rMin.kind, 'primary') && isnan(rMin.status) && strcmp(rMax.kind, 'primary') && rMin.altIndex == 0 && isempty(rMin.reason), 'reference catalog: includeAlternatives is harmless; kind primary, status NaN');
if ok, fprintf('test_extremes_status: ALL PASS\n'); else, fprintf('test_extremes_status: FAIL\n'); end
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
