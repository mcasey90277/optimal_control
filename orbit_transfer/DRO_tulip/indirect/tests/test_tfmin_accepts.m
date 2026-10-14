function ok = test_tfmin_accepts()
% TEST_TFMIN_ACCEPTS  pumpkyn tfMin accepts a certified primary unchanged
% (|dz| < 1e-6); a status-1 alternative is reported as an expected refusal,
% not an error; the dz definition is the certifier's.
%
% INPUTS:  none
% OUTPUTS: ok [logical]  every check passed (prints PASS/FAIL per check)

ok = true;
here = fileparts(fileparts(mfilename('fullpath')));
addpath(here, fullfile(fileparts(fileparts(here)), 'costate_common'));
L = load(fullfile(here, 'results', 'library_70mN_24x48_v2', 'costate_catalog_dro_tulip_70mN.mat'));  c = L.(char(fieldnames(L)));
sh = c.sheets(1);
T = transfers_at(c, sh.sD_frac(2), sh.sA_frac(5));  prim = T(strcmp({T.kind}, 'primary'));
[acc, dz, z8o, info] = tfmin_accepts(c, prim, struct('pool', capped_pool()));
ok = chk(ok, acc && dz < 1e-6 && numel(z8o) == 8 && ~info.expected, sprintf('a certified primary is accepted unchanged (|dz| = %.1e)', dz));
A = c.alternatives;  k1 = find([A.status] == 1, 1);
row1 = transfers_at(c, A(k1).sD, A(k1).sA);  row1 = row1([row1.status] == 1);  row1 = row1(1);
t1 = tic;  [acc1, dz1, ~, info1] = tfmin_accepts(c, row1, struct('pool', capped_pool()));  sec1 = toc(t1);
ok = chk(ok, ~acc1 && info1.expected && ~info1.ran && sec1 < 5 && ischar(info1.message) && ~isempty(info1.message), sprintf('REVIEW FOCUS 4: a status-1 row -> ok = false, expected = true, no error (|dz| = %.2g, %.2f s, ran = %d)', dz1, sec1, info1.ran));
[acc2, ~, ~, info2] = tfmin_accepts(c, prim, struct('pool', capped_pool(), 'capSec', 1e-3));
ok = chk(ok, ~acc2 && isfield(info2, 'timedOut') && info2.timedOut, 'REVIEW FOCUS 4: a fence timeout is a named refusal');
if ok, fprintf('test_tfmin_accepts: ALL PASS\n'); else, fprintf('test_tfmin_accepts: FAIL\n'); end
end

function ok = chk(ok, c, msg)
% CHK  Print one PASS/FAIL line and fold it into ok.  INPUTS: ok; c; msg.
% OUTPUTS: ok.
if c, fprintf('  PASS  %s\n', msg); else, fprintf('  FAIL  %s\n', msg); end
ok = ok && c;
end
