function ok = test_audit_status_layer_real()
% TEST_AUDIT_STATUS_LAYER_REAL  The primaries pass of audit_status_layer,
% for real (no seam): the 24 x 24 reference catalog reduced to ONE cell
% (has_solution cleared elsewhere), given a status layer (status 4,
% junctions from seed_from_z8 of that entry), passes through
% entry_thrust_program; corrupting one junction by 1e-3 ND (~390 km) makes
% it BAD with the junction defect named; status 3 is BAD.
%
% INPUTS:  none
% OUTPUTS: ok [logical]  every check passed (prints PASS/FAIL per check)

ok = true;
here = fileparts(fileparts(mfilename('fullpath')));
addpath(here, fullfile(fileparts(fileparts(here)), 'costate_common'));
L = load(fullfile(here, 'results', 'library_70mN_24x24_final', 'costate_catalog_dro_tulip_70mN.mat'));  c = L.(char(fieldnames(L)));
sh = c.sheets(1);  iD = find(sh.sD_frac == 0, 1);  iA = 2;     % the cell test_entry_thrust_program flies to < 1 km
k = sh.entry_index(iD, iA);
[B, ~] = arclength_arrival('setup', catalog_setup_request(c, 0));
rv0 = B.stateD(sh.sD_frac(iD));  rv0 = rv0(1:6);
seed = seed_from_z8(sh.z8(:, k), rv0, 24, B.Tnd, B.cnd, B.mu);
has = false(size(sh.has_solution));  has(iD, iA) = true;
sh.has_solution = has;
sh.status = zeros(size(has));  sh.status(iD, iA) = 4;
sh.junctions = cell(1, size(sh.z8, 2));  sh.junctions{k} = seed.Y(:, 1:24);
c.sheets = sh;
tmp = [tempname '.mat'];
A = audit_one(c, tmp);
ok = chk(ok, numel(A.primRows) == 1 && A.primRows(1).ok && A.nBad == 0 && A.nOk == 1, ...
         sprintf('a status-4 primary with its own junctions passes (%s)', A.primRows(1).why));
cb = c;  cb.sheets.junctions{k}(1:3, 5) = cb.sheets.junctions{k}(1:3, 5) + 1e-3;
A = audit_one(cb, tmp);
ok = chk(ok, A.nBad == 1 && ~A.primRows(1).ok && contains(A.primRows(1).why, 'junction defect'), ...
         sprintf('corrupted junctions are BAD: %s', A.primRows(1).why));
cb = c;  cb.sheets.status(iD, iA) = 3;
A = audit_one(cb, tmp);
ok = chk(ok, A.nBad == 1 && contains(A.primRows(1).why, 'not 4'), sprintf('status 3 primary is BAD: %s', A.primRows(1).why));
delete(tmp);
if ok, fprintf('test_audit_status_layer_real: ALL PASS\n'); else, fprintf('test_audit_status_layer_real: FAIL\n'); end
end

function A = audit_one(catalog, tmp)
% AUDIT_ONE  Save the catalog and audit it (primaries only; no alternatives).
% INPUTS: catalog [struct]; tmp [char].  OUTPUTS: A [struct] audit result.
save(tmp, 'catalog');
A = audit_status_layer(tmp, struct());
end

function ok = chk(ok, c, msg)
% CHK  Print one PASS/FAIL line and fold it into ok.  INPUTS: ok; c; msg.
% OUTPUTS: ok.
if c, fprintf('  PASS  %s\n', msg); else, fprintf('  FAIL  %s\n', msg); end
ok = ok && c;
end
