function ok = test_backfill_status_layer()
% TEST_BACKFILL_STATUS_LAYER  harvest matches primaries' junctions by root,
% classifies legacy candidates and lists the ambiguous ones for
% re-certification; assemble keeps the primaries bit-identical and builds
% the alternatives from harvest + re-certification.
%
% INPUTS:  none
% OUTPUTS: ok [logical]  every check passed (prints PASS/FAIL per check)

ok = true;
here = fileparts(fileparts(mfilename('fullpath')));
addpath(here, fullfile(fileparts(fileparts(here)), 'costate_common'));
L = load(fullfile(here, 'results', 'library_70mN_24x24_final', 'costate_catalog_dro_tulip_70mN.mat'));  rec = L.(char(fieldnames(L)));
sh = rec.sheets(1);  k11 = sh.entry_index(1, 1);  k12 = sh.entry_index(1, 2);
tmp = tempname;  mkdir(tmp);
catalog = rec;  recMat = fullfile(tmp, 'record.mat');  save(recMat, 'catalog');
legacy = @(z, reason, ok_) struct('ok', ok_, 'z', z, 'Y', z(1)*ones(14, 24), 'flyKm', 0.1, 'flyVms', 0.01, 'reason', reason, ...
                                   'tfDays', z(8), 'sD', sh.sD_frac(1), 'sA', sh.sA_frac(1), 'conj', 0);
zAlt = sh.z8(:, k11) .* [1.01; ones(7, 1)];          % the SAME t_f as the primary, listed first: matching
                                                     % junctions by t_f alone would take its Y
S = struct('problem', struct('sD', 0), 'sA', sh.sA_frac(1:2), 'TF', [1 1], ...
           'cand', {{[legacy(zAlt, 'conjugate test verdict 0', false), legacy(sh.z8(:, k11), 'certified', true)], ...
                     legacy(sh.z8(:, k12), 'certified', true)}});
save(fullfile(tmp, 'sheet.mat'), 'S');
% three ribs in one file (the filler's layout): two stalled on the same
% conjugate finding (targets 3e-4 apart), one on a polish failure; their
% point is the (1,2) primary's root, so it adds no alternative
conjStop = ': dense conjugate scan not clear: 0 coarse sign change(s), 1 zero, 0 UNRESOLVED, 0 multiplicity (3 near-miss cleared)';
pt = legacy(sh.z8(:, k12), 'certified', true);  pt.sA = sh.sA_frac(2);
R = struct('sA', sh.sA_frac(2), 'pts', {pt, pt, pt}, ...
           'stop', {['stalled at sD = 0.1720 stepping to 0.1719' conjStop], ['stalled at sD = 0.1717 stepping to 0.1716' conjStop], ...
                    'stalled at sD = 0.5000 stepping to 0.4900: normal-chart polish did not converge (|R| = 3.9e-05)'});
save(fullfile(tmp, 'ribs.mat'), 'R');
H = backfill_status_layer('harvest', recMat, {fullfile(tmp, 'sheet.mat'), fullfile(tmp, 'ribs.mat')}, tmp);
ok = chk(ok, numel(H.stops) == 1 && H.stops(1).sDto == 0.1716 && isequal(H.stops(1).z, sh.z8(:, k12)), ...
         'one stop per conjugate stall (the shortest step kept); a polish-failure stall is no stop');
ok = chk(ok, isequal(H.primJ{k11}, sh.z8(1, k11)*ones(14, 24)) && isequal(H.primJ{k12}, sh.z8(1, k12)*ones(14, 24)), ...
         'primaries get their junctions by root');
ok = chk(ok, numel(H.needRecert) == 1 && same_root(H.cands(H.needRecert).z, zAlt), 'the legacy verdict-0 candidate is listed for re-certification');
ok = chk(ok, numel(H.unmatched) == nnz(sh.has_solution) - 2, 'every other primary is listed as unmatched (needs a re-polish)');
% a fake re-certification: the verdict-0 candidate is a real FAIL with hypotheses held
Cr = H.cands(H.needRecert);  Cr.stage = 7;  Cr.status = 2;  Cr.status_reason = 'conjugate point found (gates and H6 held): x';
Cr.conjFound = true;  Cr.hypAfterConj = 'held';  Cr.overridden = false;
items = struct('kind', 'cand', 'index', H.needRecert, 'C', Cr, 'moved', false);
save(fullfile(tmp, 'recert_1.mat'), 'items');
threw = '';
try, backfill_status_layer('assemble', recMat, tmp); catch ME, threw = ME.identifier; end
ok = chk(ok, strcmp(threw, 'backfill_status_layer:unmatched'), 'assemble refuses while primaries lack junctions');
c2 = backfill_status_layer('assemble', recMat, tmp, struct('allowUnmatched', true));
ok = chk(ok, strcmp(catalog_content_key(c2), catalog_content_key(rec)), 'primaries bit-identical (content key)');
ok = chk(ok, numel(c2.alternatives) == 1 && c2.alternatives.status == 2 && ~c2.alternatives.inferred, ...
         'the alternative carries the re-certified status, not the inferred one');
ok = chk(ok, all(c2.sheets(1).status(sh.has_solution) == 4), 'every primary is status 4');
rmdir(tmp, 's');
if ok, fprintf('test_backfill_status_layer: ALL PASS\n'); else, fprintf('test_backfill_status_layer: FAIL\n'); end
end

function ok = chk(ok, c, msg)
% CHK  Print one PASS/FAIL line and fold it into ok.  INPUTS: ok; c; msg.
% OUTPUTS: ok.
if c, fprintf('  PASS  %s\n', msg); else, fprintf('  FAIL  %s\n', msg); end
ok = ok && c;
end
