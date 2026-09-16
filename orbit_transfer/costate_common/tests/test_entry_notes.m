function ok = test_entry_notes()
% TEST_ENTRY_NOTES  Every packaged entry carries a provenance note: the
% producer's "how found" clause, the certifier's remarks, the family label
% in front when a map exists, one text per solved entry in z8 order; the
% schema checks the shape; the sweep appends its remark.
%
% INPUTS:  none
% OUTPUTS: ok [logical]  every check passed (prints PASS/FAIL per check)

ok = true;
here = fileparts(fileparts(mfilename('fullpath')));
addpath(here, fullfile(fileparts(here), 'DRO_tulip', 'indirect'));
tmp = tempname;  mkdir(tmp);
cleaner = onCleanup(@() rmdir(tmp, 's'));

% ---- a sheet with two certified columns, one rib, notes on every point --
sD0 = 0;  nD = 4;  sA = [0.05 0.3 0.6];
problem = struct('version', 1, 'thrustN', 0.070, 'ispS', 900, 'm0kg', 150, 'tauDRO', 1, 'NpTulip', 7, ...
                 'pmTulip', -1, 'periodTulip', 2*pi*5/6, 'sD', sD0, 'muStar', 0.012150585609624, ...
                 'lStar', 389703.264829278, 'tStar', 382981.289129055, 'Tnd', 1.234e-4, 'cnd', 6.789);
cert = @(tf, sAv, sDv, note) struct('ok', true, 'reason', 'certified', 'note', note, 'z', [0.1*(1:7)'; tf], ...
    'Y', [], 'tfDays', tf*4.43, 'sA', sAv, 'sD', sDv, 'conj', 1, 'g', struct('minLamV', 1, 'minQmt', 1, 'dimS', 1), ...
    'level', sAv, 'arc', 1);
S = struct('sA', sA, 'TF', [4.2*4.43 NaN 4.5*4.43], 'Z8', nan(8, 3), 'cand', {cell(1, 3)}, 'problem', problem);
S.cand{1} = cert(4.2, sA(1), sD0, 'arc crossing: arc 1, level 0.0500 | lift margin 12.0x (gate 10x)');
S.cand{3} = cert(4.5, sA(3), sD0, 'seed: direct_certified (19.9 d)');
pts = [cert(4.3, sA(1), 0.75, 'rib step 1 of 3 from spine 18.6 d at sA 0.0500'), ...
       cert(4.4, sA(1), 0.5, 'rib step 2 of 3 from spine 18.6 d at sA 0.0500, 1 bisection(s) | 2 conjugate near-miss (min 1.5e-03 x median)')];
R = struct('j', 1, 'sA', sA(1), 'pts', pts, 'stop', 'complete', 'nSolve', 3);
ribs = {R};
outMat = fullfile(tmp, 'dro_tulip_test_tau1_Np7.mat');
Q = sheet_to_catalog_file(S, ribs, outMat, struct('nD', nD, 'sD0', sD0));
ok = chk(ok, isfield(Q, 'NOTE') && iscell(Q.NOTE) && isequal(size(Q.NOTE), [nD 3]), 'the sheet file carries a NOTE cell per cell');
ok = chk(ok, strcmp(Q.NOTE{1, 1}, S.cand{1}.note) && strcmp(Q.NOTE{4, 1}, pts(1).note) && strcmp(Q.NOTE{3, 1}, pts(2).note), ...
         'spine and rib notes land on their cells verbatim (no map, no family prefix)');
ok = chk(ok, isempty(Q.NOTE{2, 2}) && ~Q.OK(2, 2), 'an empty cell has an empty note');

% ---- the packager: entry_notes in z8 column order, a notes_key ----------
catMat = fullfile(tmp, 'test_cat.mat');
cat_ = build_costate_catalog_family(tmp, catMat, struct('glob', 'dro_tulip_test_*.mat', 'name', 'test_catalog', ...
    'description', 'entry-notes test', 'provenance', 'test', 'depReconstruction', 'n/a'));
sh = cat_.sheets(1);
ok = chk(ok, isfield(sh, 'entry_notes') && numel(sh.entry_notes) == nnz(sh.has_solution), ...
         sprintf('entry_notes has one text per solved entry (%d)', numel(sh.entry_notes)));
e11 = sh.entry_index(1, 1);  e41 = sh.entry_index(4, 1);  e13 = sh.entry_index(1, 3);
ok = chk(ok, strcmp(sh.entry_notes{e11}, S.cand{1}.note) && strcmp(sh.entry_notes{e41}, pts(1).note) && ...
         strcmp(sh.entry_notes{e13}, S.cand{3}.note), 'notes follow entry_index, so they stay with their z8 column');
ok = chk(ok, isfield(cat_, 'notes_key') && isfield(cat_.notes_key, 'how_found'), 'the catalog carries the notes_key');
ok = chk(ok, isempty(catalog_schema('validate', cat_)), 'the schema accepts it');

% ---- the family label goes in front when a map exists ------------------
F = struct('families', struct('label', {'fast', 'A2'}), 'labels', {{'fast', 'A2'}}, ...
           'columns', struct('family', {1, 0, -1}), 'codes', struct(), ...
           'ribFamily', @(sAv, tfPts) 1, 'attach', @(a, b) 1);
Q2 = sheet_to_catalog_file(S, ribs, fullfile(tmp, 'dro_tulip_map_tau1_Np7.mat'), struct('nD', nD, 'sD0', sD0, 'families', F));
ok = chk(ok, startsWith(Q2.NOTE{1, 1}, 'family fast | ') && startsWith(Q2.NOTE{1, 3}, 'family: unattached root') && ...
         startsWith(Q2.NOTE{4, 1}, 'family fast | rib step'), 'the family label leads the note when a map exists');

% ---- the schema refuses a misaligned notes list ------------------------
bad = cat_;  bad.sheets(1).entry_notes = bad.sheets(1).entry_notes(1:end-1);
p = catalog_schema('validate', bad);
ok = chk(ok, ~isempty(p) && contains(p{1}, 'entry_notes'), 'a short entry_notes list is refused');
bad = cat_;  bad = rmfield(bad, 'notes_key');
p = catalog_schema('validate', bad);
ok = chk(ok, ~isempty(p) && contains(p{1}, 'notes_key'), 'entry_notes without a notes_key is refused');

% ---- the sweep's remark ------------------------------------------------
% second_order_pass write-back appends 'sweep: k near-miss, m unresolved'
% to entries with something to say; exercised through its own record shape
en = sh.entry_notes;  q = e41;
en{q} = strjoin([en(q), {sprintf('sweep: %d near-miss, %d unresolved', 2, 0)}], ' | ');
ok = chk(ok, endsWith(en{q}, 'sweep: 2 near-miss, 0 unresolved') && startsWith(en{q}, 'rib step 1'), ...
         'a sweep remark appends after the producer and certifier text');

if ok, fprintf('TEST_ENTRY_NOTES: ALL PASS\n'); else, fprintf('TEST_ENTRY_NOTES: FAIL\n'); end
end

function ok = chk(ok, cond, msg)
% CHK  Print one check.  INPUTS: ok; cond; msg.  OUTPUTS: ok.
if cond, fprintf('  PASS  %s\n', msg); else, fprintf('  FAIL  %s\n', msg); ok = false; end
end
