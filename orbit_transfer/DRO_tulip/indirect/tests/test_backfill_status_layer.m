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
sh = rec.sheets(1);  k11 = sh.entry_index(1, 1);  k12 = sh.entry_index(1, 2);  k13 = sh.entry_index(1, 3);
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
% the (1,3) primary's root, certified at the NEIGHBOURING departure phase:
% another transfer, which must not donate its junctions to (1,3)
% (listed TWICE: one re-certification per (phases, root)); and a STAMPED
% candidate of a new root at (2, 4), which needs no re-certification
k24 = sh.entry_index(2, 4);
stamped = struct('ok', false, 'z', sh.z8(:, k24) .* [1.02; ones(7, 1)], 'Y', ones(14, 24), 'flyKm', 0.1, 'flyVms', 0.01, ...
                 'reason', 'H6 not established: x', 'stage', 6, 'conjFound', false, 'hypAfterConj', '', 'overridden', false, ...
                 'status', 3, 'status_reason', 'necessary only: H6 not established: x');
S = struct('problem', struct('sD', sh.sD_frac(2)), 'sA', sh.sA_frac(3:4), 'TF', [1 1], ...
           'cand', {{[legacy(sh.z8(:, k13), 'certified', true), legacy(sh.z8(:, k13), 'certified', true)], stamped}});
save(fullfile(tmp, 'sheetN.mat'), 'S');
src = {fullfile(tmp, 'sheet.mat'), fullfile(tmp, 'ribs.mat'), fullfile(tmp, 'sheetN.mat')};
H = backfill_status_layer('harvest', recMat, src, tmp);
ok = chk(ok, isempty(H.primJ{k13}), 'the same root at a neighbouring phase does not donate junctions');
ok = chk(ok, numel(H.stops) == 1 && H.stops(1).sDto == 0.1716 && isequal(H.stops(1).z, sh.z8(:, k12)), ...
         'one stop per conjugate stall (the shortest step kept); a polish-failure stall is no stop');
ok = chk(ok, isequal(H.primJ{k11}, sh.z8(1, k11)*ones(14, 24)) && isequal(H.primJ{k12}, sh.z8(1, k12)*ones(14, 24)), ...
         'primaries get their junctions by root');
% I1: EVERY inferred candidate that would become an alternative is
% re-certified (not only the verdict-0 ones), once per (phases, root); the
% primaries' own roots and stamped candidates are not
nr = H.needRecert;
isAlt = numel(nr) >= 1 && same_root(H.cands(nr(1)).z, zAlt);
isNb = numel(nr) == 2 && same_root(H.cands(nr(2)).z, sh.z8(:, k13)) && abs(H.cands(nr(2)).sD - sh.sD_frac(2)) < 1e-12;
ok = chk(ok, numel(nr) == 2 && isAlt && isNb, sprintf(['I1: needRecert = the verdict-0 candidate and the inferred ' ...
         'neighbouring-phase root, once (got %d item(s))'], numel(nr)));
ok = chk(ok, ~any(arrayfun(@(m) isequal(H.cands(m).stage, 6), nr)), 'I1: a stamped candidate is not re-certified');
ok = chk(ok, numel(H.unmatched) == nnz(sh.has_solution) - 2 && ismember(k13, H.unmatched), 'every other primary is listed as unmatched (needs a re-polish)');
% a fake re-certification: the verdict-0 candidate is a real FAIL with hypotheses held
Cr = H.cands(H.needRecert(1));  Cr.stage = 7;  Cr.status = 2;  Cr.status_reason = 'conjugate point found (gates and H6 held): x';
Cr.conjFound = true;  Cr.hypAfterConj = 'held';  Cr.overridden = false;
items = struct('kind', 'cand', 'index', H.needRecert(1), 'C', Cr, 'moved', false);
harvestKey = H.harvestKey;
save(fullfile(tmp, 'recert_1.mat'), 'items', 'harvestKey');
threw = '';
try, backfill_status_layer('assemble', recMat, tmp); catch ME, threw = ME.identifier; end
ok = chk(ok, strcmp(threw, 'backfill_status_layer:unmatched'), 'assemble refuses while primaries lack junctions');
c2 = backfill_status_layer('assemble', recMat, tmp, struct('allowUnmatched', true));
ok = chk(ok, strcmp(catalog_content_key(c2), catalog_content_key(rec)), 'primaries bit-identical (content key)');
a = c2.alternatives;  ka = find(arrayfun(@(r) same_root(r.z8, zAlt), a));
ok = chk(ok, numel(a) == 3 && isscalar(ka) && a(ka).status == 2 && ~a(ka).inferred, ...
         'the alternative carries the re-certified status, not the inferred one');
ok = chk(ok, any(arrayfun(@(r) r.iD == 2 && r.iA == 3 && same_root(r.z8, sh.z8(:, k13)), a)), ...
         'the neighbouring-phase root is an alternative of its own cell');
kn = find(arrayfun(@(r) r.iD == 2 && r.iA == 3 && same_root(r.z8, sh.z8(:, k13)), a));
ok = chk(ok, isscalar(kn) && isequal(a(kn).status, 4) && a(kn).inferred, ...
         'C2: an un-re-certified harmonised candidate (stage = [], status = []) ships classified and inferred, not stamped');
% a recert file made against another harvest is refused, and so is a re-harvest over recert files
harvestKey = repmat('0', 1, 32);
save(fullfile(tmp, 'recert_2.mat'), 'items', 'harvestKey');
threw = '';
try, backfill_status_layer('assemble', recMat, tmp, struct('allowUnmatched', true)); catch ME, threw = ME.identifier; end
ok = chk(ok, strcmp(threw, 'backfill_status_layer:harvestKey'), 'assemble refuses a recert file of another harvest');
delete(fullfile(tmp, 'recert_2.mat'));
threw = '';
try, backfill_status_layer('harvest', recMat, src, tmp); catch ME, threw = ME.identifier; end
ok = chk(ok, strcmp(threw, 'backfill_status_layer:staleRecert'), 'harvest refuses while recert files sit in outDir');
ok = chk(ok, all(c2.sheets(1).status(sh.has_solution) == 4), 'every primary is status 4');
rmdir(tmp, 's');

% ---- M1: the job's one decision token -------------------------------------------
cmp = struct('nAgree', 10, 'nBoth', 10, 'nOnlyRef', 0, 'nNewSlower', 0, 'nNewFaster', 0);
V = struct('nBad', 0, 'coverageOk', true, 'primNot4', [], 'movedPrimaries', [], 'contentKeyEqual', true, 'cmp', cmp);
ok = chk(ok, strcmp(backfill_clean_verdict(V), 'CLEAN'), 'M1: everything holds -> CLEAN');
bads = {setf(V, 'nBad', 1), setf(V, 'coverageOk', false), setf(V, 'primNot4', 3), setf(V, 'movedPrimaries', 7), ...
        setf(V, 'contentKeyEqual', false), setf(V, 'cmp', setf(cmp, 'nAgree', 9)), setf(V, 'cmp', setf(cmp, 'nOnlyRef', 1)), ...
        setf(V, 'cmp', setf(cmp, 'nNewSlower', 1)), setf(V, 'cmp', setf(cmp, 'nNewFaster', 1)), setf(V, 'nBad', NaN)};
tok = cellfun(@backfill_clean_verdict, bads, 'UniformOutput', false);
ok = chk(ok, all(strcmp(tok, 'NOT CLEAN')), sprintf('M1: each single failure -> NOT CLEAN (%d of %d)', nnz(strcmp(tok, 'NOT CLEAN')), numel(bads)));
[~, why] = backfill_clean_verdict(setf(V, 'primNot4', [3 4]));
ok = chk(ok, contains(why, 'primNot4'), sprintf('M1: NOT CLEAN says why (%s)', why));
if ok, fprintf('test_backfill_status_layer: ALL PASS\n'); else, fprintf('test_backfill_status_layer: FAIL\n'); end
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
