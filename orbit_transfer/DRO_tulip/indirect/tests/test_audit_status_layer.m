function ok = test_audit_status_layer()
% TEST_AUDIT_STATUS_LAYER  The status audit fails closed: a reproduced
% status passes; a different status, a root that moved, and corrupted
% junctions are each a BAD row naming why; the content keys are bound.
% The alternatives carry REAL junctions (seed_from_z8 of entry (1,1)) since
% the audit flies them before re-certifying; without a pool and without
% .allowUnfenced every alternative is BAD, not a crash.
%
% INPUTS:  none
% OUTPUTS: ok [logical]  every check passed (prints PASS/FAIL per check)

ok = true;
here = fileparts(fileparts(mfilename('fullpath')));
addpath(here, fullfile(fileparts(fileparts(here)), 'costate_common'));
L = load(fullfile(here, 'results', 'library_70mN_24x24_final', 'costate_catalog_dro_tulip_70mN.mat'));  c = L.(char(fieldnames(L)));
sh = c.sheets(1);  k = sh.entry_index(1, 1);
[B, ~] = arclength_arrival('setup', catalog_setup_request(c, sh.sD_frac(1)));
rv0 = B.stateD(sh.sD_frac(1));
seed = seed_from_z8(sh.z8(:, k), rv0(1:6), 24, B.Tnd, B.cnd, B.mu);
row = struct('sD', sh.sD_frac(1), 'sA', sh.sA_frac(1), 'iD', 0, 'iA', 0, 'z8', sh.z8(:, k), 'junctions', seed.Y(:, 1:24), ...
             'tf_nd', sh.z8(8, k), 'status', 3, 'status_reason', 'x', 'inferred', false, 'conj', 1, 'conjVerdict', 'PASS', ...
             'minLamV', 1, 'minQmt', 1, 'dimS', 1, 'h6Margin', 1, 'liftMargin', 1, 'flyKm', 0, 'flyVms', 0, 'source', 't', 'sheet', 1);
bad = row;  bad.junctions(1:3, 5) = bad.junctions(1:3, 5) + 1e-3;       % ~390 km off at junction 5
nanJ = row;  nanJ.junctions(4, 7) = NaN;
c.alternatives = [row, setf(row, 'status', 2), row, bad, nanJ];
fake = @(st, dz) @(seed, rv0, rvf, B, o) struct('z', [seed.Y(8:14, 1); seed.tf] + dz, 'status', st, 'reason', 'fake');
tmp = [tempname '.mat'];  catalog = c;  save(tmp, 'catalog');
A = audit_status_layer(tmp, struct('skipPrimaries', true, 'idxAlt', 1, 'certifier', fake(3, 0), 'pool', [], 'allowUnfenced', true));
ok = chk(ok, A.nBad == 0 && A.altRows(1).ok, 'a reproduced status passes');
A = audit_status_layer(tmp, struct('skipPrimaries', true, 'idxAlt', 2, 'certifier', fake(3, 0), 'pool', [], 'allowUnfenced', true));
ok = chk(ok, A.nBad == 1 && contains(A.altRows(1).why, 'status not reproduced'), 'a different status is BAD, named');
A = audit_status_layer(tmp, struct('skipPrimaries', true, 'idxAlt', 3, 'certifier', fake(3, [0.01; zeros(7, 1)]), 'pool', [], 'allowUnfenced', true));
ok = chk(ok, A.nBad == 1 && A.altRows(1).moved && contains(A.altRows(1).why, 'moved'), 'REVIEW FOCUS 5: a root that moved is BAD, not relabelled');
A = audit_status_layer(tmp, struct('skipPrimaries', true, 'idxAlt', 1, 'pool', [], 'allowUnfenced', true, ...
                                   'certifier', @(seed, rv0, rvf, B, o) struct('z', nan(8, 1), 'status', -1, 'reason', 'normal-chart polish did not converge')));
ok = chk(ok, A.nBad == 1 && ~A.altRows(1).moved && contains(A.altRows(1).why, 'below the floor on re-certification: normal-chart polish'), ...
         sprintf('M4: a failed polish is "below the floor on re-certification", not "moved" (%s)', A.altRows(1).why));
ok = chk(ok, numel(A.altContentKey) == 32 && ~strcmp(A.altContentKey, alternatives_content_key(setf(c, 'alternatives', row))), ...
         'the alternatives key is bound to the table''s content');
c2 = c;  c2.alternatives(2).status = 3;                                  % same rows, one status differs
ok = chk(ok, ~strcmp(alternatives_content_key(c), alternatives_content_key(c2)), ...
         'the alternatives key changes when only a status changes');
ok = chk(ok, isempty(alternatives_content_key(rmfield(c, 'alternatives'))), 'no alternatives -> empty key');
A = audit_status_layer(tmp, struct('skipPrimaries', true, 'idxAlt', [1 3], 'certifier', fake(3, 0), 'pool', []));
ok = chk(ok, A.nBad == 2 && all(contains({A.altRows.why}, 'no pool: refusing to re-certify unfenced')), ...
         'no pool and no allowUnfenced -> every alternative BAD, not a crash');
boom = @(varargin) error('test:called', 'the certifier must not be called');
% C2: a stored status that is not a real scalar ([] from field harmonising)
% is BAD -- `~(3 == [])` is empty and used to read as reproduced
catalog = setf(c, 'alternatives', setf(row, 'status', []));  save(tmp, 'catalog');
A = audit_status_layer(tmp, struct('skipPrimaries', true, 'idxAlt', 1, 'certifier', fake(3, 0), 'pool', [], 'allowUnfenced', true));
ok = chk(ok, A.nBad == 1 && ~A.altRows(1).ok && contains(A.altRows(1).why, 'stored status'), ...
         sprintf('C2: an empty stored status is BAD, not reproduced (%s)', A.altRows(1).why));
catalog = c;  save(tmp, 'catalog');
A = audit_status_layer(tmp, struct('skipPrimaries', true, 'idxAlt', 4, 'certifier', boom, 'pool', [], 'allowUnfenced', true));
ok = chk(ok, A.nBad == 1 && contains(A.altRows(1).why, 'junction defect') && contains(A.altRows(1).why, ' km'), ...
         sprintf('corrupted alternative junctions are BAD before the certifier runs (%s)', A.altRows(1).why));
A = audit_status_layer(tmp, struct('skipPrimaries', true, 'idxAlt', 5, 'certifier', boom, 'pool', [], 'allowUnfenced', true));
ok = chk(ok, A.nBad == 1 && ~A.altRows(1).ok && ~contains(A.altRows(1).why, 'must not be called'), ...
         sprintf('a NaN junction is BAD, not passed (%s)', A.altRows(1).why));
vel = row;  vel.junctions(4:6, 5) = vel.junctions(4:6, 5) + 2e-3;       % ~2 m/s at junction 5
catalog = setf(c, 'alternatives', vel);  save(tmp, 'catalog');
A = audit_status_layer(tmp, struct('skipPrimaries', true, 'certifier', boom, 'pool', [], 'allowUnfenced', true, 'junctionKm', 1e9));
ok = chk(ok, A.nBad == 1 && contains(A.altRows(1).why, 'm/s'), ...
         sprintf('a velocity junction defect is BAD (position gate opened) (%s)', A.altRows(1).why));
catalog = rmfield(c, 'rungs_N');  save(tmp, 'catalog');
A = audit_status_layer(tmp, struct('skipPrimaries', true, 'idxAlt', [1 2], 'certifier', fake(3, 0), 'pool', [], 'allowUnfenced', true));
ok = chk(ok, A.nBad == 2 && all(contains({A.altRows.why}, 'setup failed')), ...
         'a failed setup makes every pending row BAD; the audit returns');
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
