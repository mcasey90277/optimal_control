function ok = test_merge_phase_catalogs()
% TEST_MERGE_PHASE_CATALOGS  merge_phase_catalogs takes the DONOR's entry
% wherever the donor is faster or the base has none, carries every per-cell
% grid and per-entry field with it, leaves every other cell alone, and
% refuses catalogs that are not the same problem. Synthetic catalogs: a
% 3 x 4 base (the finer grid) and a 3 x 2 donor on its odd arrival phases.
%
% INPUTS:  none
% OUTPUTS: ok [logical]  every check passed (prints PASS/FAIL per check)

ok = true;
here = fileparts(fileparts(mfilename('fullpath')));            % DRO_tulip/indirect
addpath(here);
tStar = 86400;                                                  % 1 ND = 1 day, so tolerances read in days

base = synth([0 1/3 2/3], [0 0.25 0.5 0.75], tStar);
donor = synth([0 1/3 2/3], [0 0.5], tStar);
% base cells (row, col) and the donor's view of them (row, col/2 rounded up)
%  (1,1) both, donor FASTER by 0.5 d          -> donor replaces, in place
%  (2,1) both, donor slower by 0.5 d          -> base kept
%  (3,1) both, donor faster by 1e-4 d (same)  -> base kept (inside tolDays)
%  (1,3) base HOLE, donor has it              -> donor fills, new z8 column
%  (2,3) base hole, donor hole                -> stays a hole
%  (3,3) both, donor equal                    -> base kept
base.sheets.tf_nd = 20*ones(3, 4);
donor.sheets.tf_nd = [19.5 20; 20.5 20; 20 - 1e-4 20];
base.sheets = setHole(base.sheets, [1 3; 2 3]);
donor.sheets = setHole(donor.sheets, 2, 2);
donor.sheets.conj_pass(:) = int8(7);  donor.sheets.h6_margin(:) = 99;  donor.sheets.family_index(:) = int8(4);
donor.sheets.z8(1, :) = 1000 + (1:size(donor.sheets.z8, 2));    % donor entries recognisable by row 1
nBase0 = size(base.sheets.z8, 2);

[M, info] = merge_phase_catalogs(base, donor, struct('donorTag', 'record'));
s = M.sheets;  b = base.sheets;
k11 = s.entry_index(1, 1);  k13 = s.entry_index(1, 3);
ok = chk(ok, s.tf_nd(1, 1) == 19.5 && k11 == b.entry_index(1, 1) && s.z8(1, k11) == 1000 + donor.sheets.entry_index(1, 1), ...
         'donor faster: replaces the base entry IN PLACE (same z8 column, donor z8)');
ok = chk(ok, s.has_solution(1, 3) && s.tf_nd(1, 3) == donor.sheets.tf_nd(1, 2) && k13 == nBase0 + 1 ...
         && s.z8(1, k13) == 1000 + donor.sheets.entry_index(1, 2), 'base hole: filled by the donor in a NEW z8 column');
ok = chk(ok, s.conj_pass(1, 1) == 7 && s.h6_margin(1, 3) == 99 && s.family_index(1, 3) == 4, ...
         'a taken cell carries the donor''s per-cell grids (verdicts, margins, family)');
ok = chk(ok, isequal(s.tf_nd(2, 1), 20) && isequal(s.z8(:, s.entry_index(2, 1)), b.z8(:, b.entry_index(2, 1))) ...
         && s.conj_pass(2, 1) == b.conj_pass(2, 1), 'donor slower: base kept, untouched');
ok = chk(ok, s.tf_nd(3, 1) == 20 && s.tf_nd(3, 3) == 20, 'within tolDays or equal: base kept');
ok = chk(ok, ~s.has_solution(2, 3), 'a hole in both stays a hole');
offGrid = [2 4];
ok = chk(ok, isequal(s.tf_nd(:, offGrid), b.tf_nd(:, offGrid)) && isequal(s.conj_pass(:, offGrid), b.conj_pass(:, offGrid)), ...
         'cells off the donor''s grid are untouched');
ok = chk(ok, contains(s.entry_notes{k11}, 'merged from record') && contains(s.entry_notes{k13}, 'merged from record') ...
         && ~contains(s.entry_notes{s.entry_index(2, 1)}, 'merged'), 'taken entries say where they came from; others do not');
ok = chk(ok, M.n_entries == nnz(s.has_solution) && numel(s.entry_notes) == size(s.z8, 2), 'n_entries and per-entry arrays agree');
ok = chk(ok, info.nFaster == 1 && info.nFilled == 1 && size(info.rows, 1) == 2 && isfield(M, 'merge'), ...
         'the merge record counts and lists what was taken');
all_idx = s.entry_index(s.has_solution);
ok = chk(ok, numel(unique(all_idx)) == numel(all_idx) && all(all_idx >= 1 & all_idx <= size(s.z8, 2)), ...
         'every entry_index distinct and in range');

% ---- the status layer (optimality-status spec 5) ----------------------------
A = M.alternatives;
ok = chk(ok, numel(A) == 1 && A(1).iD == 1 && A(1).iA == 1 && A(1).status == 4 && strcmp(A(1).source, 'displaced primary (merge)') ...
         && isequal(A(1).z8, b.z8(:, b.entry_index(1, 1))), 'the displaced base primary is kept as an alternative');
ok = chk(ok, s.status(1, 3) == 4 && numel(s.junctions) == size(s.z8, 2), 'status grid and junctions travel with taken entries');
Cd = struct('ok', true, 'stage', 8, 'status', 3, 'status_reason', 'donor alt', 'z', [0.3*ones(7, 1); 21], ...
            'Y', ones(14, 24), 'flyKm', 0, 'flyVms', 0, 'overridden', false);
d6 = donor;  d6.alternatives = make_alternative(Cd, donor.sheets.sD_frac(2), donor.sheets.sA_frac(2), 2, 2, 'donor alt');
M6 = merge_phase_catalogs(base, d6);
iAlt = find(strcmp({M6.alternatives.source}, 'donor alt'));
ok = chk(ok, numel(M6.alternatives) == 2 && isscalar(iAlt) && M6.alternatives(iAlt).iD == 2 && M6.alternatives(iAlt).iA == 3, ...
         'donor alternatives are remapped onto the base grid (donor col 2 -> base col 3)');
ok = chk(ok, isnan(A(1).flyKm) && isnan(A(1).flyVms), 'M2: the displaced primary''s flight is NaN (not measured), not a fake 0');
% I4: donor alternatives are placed by PHASE, not by the donor's indices: an
% alternative OFF the donor grid (iD = iA = 0) but ON the base grid takes the
% base cell, and is then deduplicated against that cell's primary
kb12 = b.entry_index(1, 2);  kb24 = b.entry_index(2, 4);
Ct = setfs(Cd, 'z', b.z8(:, kb12));                                   % the base (1,2) primary's root
Cu = setfs(Cd, 'z', [0.7*ones(7, 1); 22]);                            % a new root at base (2,4)
d7 = donor;  d7.alternatives = [make_alternative(Ct, base.sheets.sD_frac(1), base.sheets.sA_frac(2), 0, 0, 'twin off donor grid'), ...
                                make_alternative(Cu, base.sheets.sD_frac(2), base.sheets.sA_frac(4), 0, 0, 'new off donor grid')];
M7 = merge_phase_catalogs(base, d7);
src7 = {M7.alternatives.source};
iu = find(strcmp(src7, 'new off donor grid'));
ok = chk(ok, ~any(strcmp(src7, 'twin off donor grid')), 'I4: an off-donor-grid twin of a base primary is deduplicated away');
ok = chk(ok, isscalar(iu) && M7.alternatives(iu).iD == 2 && M7.alternatives(iu).iA == 4, ...
         'I4: an off-donor-grid alternative on the base grid takes the base cell');
bOld = rmfield(base, 'alternatives');  bOld.sheets = rmfield(bOld.sheets, {'status', 'status_reason', 'junctions'});
dOld = rmfield(donor, 'alternatives');  dOld.sheets = rmfield(dOld.sheets, {'status', 'status_reason', 'junctions'});
[Mo, io] = merge_phase_catalogs(bOld, dOld);
ok = chk(ok, io.nFaster == 1 && io.nFilled == 1 && ~isfield(Mo, 'alternatives'), 'REVIEW FOCUS 2: catalogs without the layer merge as before');
ok = chk(ok, refuses(base, dOld), 'refuses a layered base with an unlayered donor (no half-merge)');
ok = chk(ok, refuses(bOld, donor), 'refuses an unlayered base with a layered donor (layer would be dropped)');

% ---- refusals: not the same problem ----------------------------------------
ok = chk(ok, refuses(base, setf(donor, 'thruster', struct('isp_s', 1710, 'm0_kg', 150))), 'refuses another engine');
d2 = donor;  d2.sheets.Np = 8;
ok = chk(ok, refuses(base, d2), 'refuses another arrival orbit');
d3 = donor;  d3.families.labels{1} = 'other';
ok = chk(ok, refuses(base, d3), 'refuses another family numbering (codes would mean other families)');
d4 = donor;  d4.second_order.K = 12;
ok = chk(ok, refuses(base, d4), 'refuses second-order verdicts made with other settings');
d5 = donor;  d5.sheets.sD_frac = [0.1 0.4 0.7];
ok = chk(ok, refuses(base, d5), 'refuses a donor sharing no cell with the base');

if ok, fprintf('test_merge_phase_catalogs: ALL PASS\n'); else, fprintf('test_merge_phase_catalogs: FAIL\n'); end
end

% ---------------------------------------------------------------------------
function c = synth(sD, sA, tStar)
% SYNTH  A minimal phase catalog: one sheet, every cell solved, entries
% numbered column-major.  INPUTS: sD; sA; tStar.  OUTPUTS: c.
nD = numel(sD);  nA = numel(sA);  n = nD*nA;
sh = struct('tauDRO', 1, 'Np', 7, 'pm', -1, 'sD_frac', sD, 'sA_frac', sA, ...
            'has_solution', true(nD, nA), 'tf_nd', 20*ones(nD, nA), 'entry_index', reshape(1:n, nD, nA), ...
            'z8', [rand(7, n); 20*ones(1, n)], 'conj_pass', int8(ones(nD, nA)), 'h6_margin', 5*ones(nD, nA), ...
            'family_index', int8(ones(nD, nA)), 'entry_notes', {arrayfun(@(k) sprintf('fam | found %d', k), 1:n, 'UniformOutput', false)}, ...
            'status', int8(4*ones(nD, nA)), 'status_reason', {repmat({'full stack passed'}, 1, n)}, ...
            'junctions', {repmat({ones(14, 24)}, 1, n)});
c = struct('constants', struct('tStar_s', tStar), 'thruster', struct('isp_s', 900, 'm0_kg', 150), 'rungs_N', 0.07, ...
           'families', struct('labels', {{'a', 'b', 'c', 'd', 'e'}}), ...
           'second_order', struct('instruments', 'x', 'nSub', 8, 'K', 24, 'relTolPair', [1e-12 1e-9]), ...
           'n_entries', n, 'sheets', sh, 'alternatives', struct([]));
end

function sh = setHole(sh, rows, cols)
% SETHOLE  Empty cells: either rows = [r c; ...] pairs, or one (rows, cols).
% INPUTS: sh; rows; cols (optional).  OUTPUTS: sh.
if nargin == 3, rows = [rows cols]; end
for k = 1:size(rows, 1)
    sh.has_solution(rows(k, 1), rows(k, 2)) = false;  sh.tf_nd(rows(k, 1), rows(k, 2)) = NaN;
    sh.entry_index(rows(k, 1), rows(k, 2)) = 0;
end
end

function c = setf(c, f, v)
% SETF  Set a top-level field.  INPUTS: c; f; v.  OUTPUTS: c.
c.(f) = v;
end

function s = setfs(s, varargin)
% SETFS  Set name/value pairs on a struct.  INPUTS: s; pairs.  OUTPUTS: s.
for k = 1:2:numel(varargin), s.(varargin{k}) = varargin{k+1}; end
end

function r = refuses(base, donor)
% REFUSES  True if the merge throws.  INPUTS: base; donor.  OUTPUTS: r.
try
    merge_phase_catalogs(base, donor);  r = false;
catch
    r = true;
end
end

function ok = chk(ok, c, msg)
% CHK  Print one PASS/FAIL line and fold it into ok.  INPUTS: ok; c; msg.
% OUTPUTS: ok.
if c, fprintf('  PASS  %s\n', msg); else, fprintf('  FAIL  %s\n', msg); end
ok = ok && c;
end
