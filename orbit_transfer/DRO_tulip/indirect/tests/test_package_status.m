function ok = test_package_status()
% TEST_PACKAGE_STATUS  sheet_to_catalog_file + build_costate_catalog_family
% carry status, status_reason and junctions for every primary, and an
% alternatives table holding (a) the sheet's non-winning candidates above
% the floor, (b) a rib's refusals, (c) a certified rib point that lost its
% cell to a faster one -- deduplicated, overridden ones excluded.
%
% INPUTS:  none
% OUTPUTS: ok [logical]  every check passed (prints PASS/FAIL per check)

ok = true;
here = fileparts(fileparts(mfilename('fullpath')));
addpath(here, fullfile(fileparts(fileparts(here)), 'costate_common'));
S96 = load(fullfile(here, 'results', 'sheet96_resolution_test', 'arrival_sheet_70mN_nA96.mat'));  P = S96.S.problem;
cert = @(z1, tf, sD, sA) struct('ok', true, 'status', 4, 'status_reason', 'full stack passed', 'stage', 8, ...
        'z', [z1; (2:7).'; tf], 'Y', z1*ones(14, 24), 'tfDays', tf*P.tStar/86400, 'sD', sD, 'sA', sA, ...
        'flyKm', 0.1, 'flyVms', 0.01, 'reason', 'certified', 'note', '', 'conj', 1, 'conjVerdict', 'PASS', ...
        'g', struct('minLamV', 1, 'minQmt', 1, 'dimS', 1), 'h6Margin', 5, 'liftMargin', 30, 'overridden', false);
ref = @(C, st) setf(C, 'ok', false, 'status', st, 'status_reason', sprintf('status %d', st), 'reason', 'x');
S = struct('problem', P, 'sA', [0.25 0.75 0.5], 'TF', [NaN NaN NaN]);
c1 = [cert(1, 4.0, 0, 0.25), ref(cert(2, 4.5, 0, 0.25), 2), ref(setf(cert(3, 4.6, 0, 0.25), 'flyKm', 900), 3)];
c2 = [cert(4, 4.2, 0, 0.75), setf(ref(cert(5, 4.3, 0, 0.75), 3), 'overridden', true)];
% column 3 has NO certified winner (S.TF NaN): its refusals above the floor
% (status 2 and 3) must still reach the alternatives table
c3 = [ref(cert(9, 4.7, 0, 0.5), 2), ref(cert(10, 4.8, 0, 0.5), 3)];
S.cand = {c1, c2, c3};  S.TF = [c1(1).tfDays, c2(1).tfDays, NaN];
rib = struct('sA', 0.25, 'pts', [cert(6, 4.1, 0.5, 0.25)], 'refused', [ref(cert(7, 4.4, 0.5, 0.25), 1)], 'stop', 'complete');
% rib2 also refuses a TWIN of its own primary (8) and a REPEAT of rib1's
% refusal (7): both must be deduplicated away (mutation check: without
% dedup_alternatives the table would read [2 6 7 7 8])
rib2 = struct('sA', 0.25, 'pts', [cert(8, 4.05, 0.5, 0.25)], ...
        'refused', [ref(cert(8, 4.05, 0.5, 0.25), 2), ref(cert(7, 4.4, 0.5, 0.25), 1)], 'stop', 'complete');
% I3: a rib with NO certified point still contributes its refusals (11)
rib3 = struct('sA', 0.75, 'pts', [], 'refused', [ref(cert(11, 4.9, 0.5, 0.75), 3)], ...
              'stop', 'stalled at sD = 0 stepping to 0.5: x');
% C2: a field-harmonised point (status = [], stage = []) is classified on
% the way in, not stored as an empty status (12, at the empty cell (2, 3))
rib4 = struct('sA', 0.5, 'pts', [setf(cert(12, 4.95, 0.5, 0.5), 'status', [], 'stage', [], 'status_reason', [])], ...
              'refused', [], 'stop', 'complete');
tmp = tempname;  mkdir(tmp);
problem = P;                                   % ribs carry the sheet's identity, as build_ribs saves it
save(fullfile(tmp, 'sheet.mat'), 'S');
R = rib;   save(fullfile(tmp, 'rib1.mat'), 'R', 'problem');
R = rib2;  save(fullfile(tmp, 'rib2.mat'), 'R', 'problem');
R = rib3;  save(fullfile(tmp, 'rib3.mat'), 'R', 'problem');
R = rib4;  save(fullfile(tmp, 'rib4.mat'), 'R', 'problem');
c = package_phase_catalog(fullfile(tmp, 'sheet.mat'), {fullfile(tmp, 'rib1.mat'), fullfile(tmp, 'rib2.mat'), ...
                          fullfile(tmp, 'rib3.mat'), fullfile(tmp, 'rib4.mat')}, ...
        struct('nD', 2, 'outDir', tmp, 'name', 'cat_test', 'thrustN', P.thrustN, 'ispS', P.ispS, 'm0kg', P.m0kg));
sh = c.sheets(1);
ok = chk(ok, isequal(sh.status, int8([4 4 0; 4 0 4])) && numel(sh.junctions) == 4 && numel(sh.status_reason) == 4, ...
         'primaries: status 4 grid, junctions and reasons per entry');
k = sh.entry_index(2, 1);
ok = chk(ok, sh.z8(1, k) == 8 && isequal(sh.junctions{k}, 8*ones(14, 24)), 'the faster rib point is the primary and carries its junctions');
A = c.alternatives;  firsts = arrayfun(@(a) a.z8(1), A);
ok = chk(ok, isequal(sort(firsts), [2 6 7 9 10 11]), sprintf('alternatives: sheet refusal (2), rib refusal (7), the slower certified rib point (6), the winnerless column''s refusals (9, 10), the pointless rib''s refusal (11); twin (8) and repeat (7) deduplicated -- got %s', mat2str(sort(firsts))));
ok = chk(ok, any(firsts == 11) && A(firsts == 11).iD == 2 && A(firsts == 11).iA == 2, 'I3: a rib with no certified point still records its refusals');
k12 = sh.entry_index(2, 3);
ok = chk(ok, sh.status(2, 3) == 4 && ischar(sh.status_reason{k12}) && ~isempty(sh.status_reason{k12}), ...
         'C2: a harmonised primary (status = []) is classified (4), not a crash or an empty status');
d6 = A(firsts == 6);
ok = chk(ok, isscalar(d6) && d6.status == 4 && isnan(d6.flyKm) && isnan(d6.flyVms) && contains(d6.source, 'displaced primary'), ...
         'M2: the displaced primary (6) is a status-4 row whose flight is NaN (not measured), not a fake 0');
ok = chk(ok, ~any(firsts == 3) && ~any(firsts == 5), 'below the floor (3) and overridden (5) are not recorded');
w = A(firsts == 9 | firsts == 10);
ok = chk(ok, numel(w) == 2 && all([w.iD] == 1) && all([w.iA] == 3) && isequal(sort([w.status]), [2 3]), ...
         'a column with no certified winner still records its candidates (iD = spine row, iA = that column)');
ok = chk(ok, isfield(c, 'status_key') && isequal(c.status_key.codes, [4 3 2 1 0 -1]), 'the legend ships in the catalog');
rmdir(tmp, 's');
if ok, fprintf('test_package_status: ALL PASS\n'); else, fprintf('test_package_status: FAIL\n'); end
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
