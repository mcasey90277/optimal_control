function ok = test_reproduce_library()
% TEST_REPRODUCE_LIBRARY  The two pieces of logic reproduce_library_70mN adds
% to the pipeline, tested without a solve:
%
%   compare_phase_catalogs  the library of record against ITSELF is a clean
%       match; a copy with one flight time moved, one cell removed, one
%       costate vector changed and one cell's family reassigned reports
%       exactly those four differences; relabelling the families (the same
%       partition under other index numbers) is NOT a difference; a catalog
%       on another grid is refused.
%   adopt_walked_arcs  copies walked arcs into a campaign's arc folder under
%       the driver's tagged name, leaves an identical copy alone, refuses to
%       overwrite a different file, and names a missing source.
%
% INPUTS:  none
% OUTPUTS: ok [logical]  every check passed (prints PASS/FAIL per check)

ok = true;
here = fileparts(fileparts(mfilename('fullpath')));            % DRO_tulip/indirect
addpath(here);
recMat = fullfile(here, 'results', 'library_70mN_24x24_final', 'costate_catalog_dro_tulip_70mN.mat');
assert(isfile(recMat), 'the library of record is missing: %s', recMat);
L = load(recMat);  fn = fieldnames(L);  ref = L.(fn{1});
quiet = struct('print', false);

% ---- the record against itself -------------------------------------------
R = compare_phase_catalogs(ref, ref, quiet);
ok = chk(ok, R.ok && R.nBoth == 576 && R.nOnlyRef == 0 && R.nOnlyNew == 0, sprintf('record vs itself: %d common cells, clean', R.nBoth));
ok = chk(ok, R.worstTfDays == 0 && R.worstZRel == 0 && R.nFamilyDiff == 0, 'record vs itself: zero deviation');

% ---- four planted differences, each counted once -------------------------
new = ref;  s = new.sheets(1);
s.tf_nd(3, 5) = s.tf_nd(3, 5) + 1e-3;                          % 1. a flight time (1e-3 ND = 0.0044 d)
s.has_solution(7, 2) = false;                                  % 2. a cell missing from the new catalog
k = s.entry_index(10, 10);  s.z8(1:7, k) = 1.01*s.z8(1:7, k);  % 3. a costate vector
fi = double(s.family_index(12, 20));  others = setdiff(unique(double(s.family_index(s.has_solution))), fi);
s.family_index(12, 20) = others(1);                            % 4. one cell moved to another family
new.sheets(1) = s;
R = compare_phase_catalogs(new, ref, quiet);
ok = chk(ok, ~R.ok, 'planted differences: not a match');
ok = chk(ok, R.nTfOver == 1 && R.nOnlyRef == 1 && R.nZOver == 1 && R.nFamilyDiff == 1, ...
         sprintf('planted differences: tf %d, missing %d, z %d, family %d (each 1)', R.nTfOver, R.nOnlyRef, R.nZOver, R.nFamilyDiff));
ok = chk(ok, abs(R.worstTfDays - 1e-3*ref.constants.tStar_s/86400) < 1e-9, sprintf('worst flight-time deviation reported in days (%.4f)', R.worstTfDays));
ok = chk(ok, any(R.cells.iD == 3 & R.cells.iA == 5) && any(R.cells.iD == 7 & R.cells.iA == 2), 'the differing cells are named');

% ---- the same family partition under other index numbers is a match ------
new = ref;  s = new.sheets(1);  fiAll = double(s.family_index);  u = unique(fiAll(s.has_solution & fiAll >= 1));   % FAMILIES only: status codes (< 1) are not relabelled
perm = circshift(u, 1);  relabelled = fiAll;
for q = 1:numel(u), relabelled(fiAll == u(q)) = perm(q); end
s.family_index = cast(relabelled, 'like', s.family_index);  new.sheets(1) = s;
R = compare_phase_catalogs(new, ref, quiet);
ok = chk(ok, R.ok && R.nFamilyDiff == 0, 'families relabelled (same partition): still a match');

% ---- FALSE PASSES closed (review 2026-09-18) -------------------------------
% In MATLAB `NaN > tol` is false, so a non-finite entry used to compare EQUAL.
new = ref;  new.sheets(1).tf_nd(5, 5) = NaN;
R = compare_phase_catalogs(new, ref, quiet);
ok = chk(ok, ~R.ok && R.nNonfinite == 1, 'a NaN flight time in a covered cell is a difference, not a match');
new = ref;  k = new.sheets(1).entry_index(6, 6);  new.sheets(1).z8(3, k) = NaN;
R = compare_phase_catalogs(new, ref, quiet);
ok = chk(ok, ~R.ok && R.nNonfinite == 1, 'a NaN costate is a difference, not a match');
% a catalog with no family map: the family comparison is NOT ESTABLISHED, never "0 differences"
new = ref;  new.sheets = rmfield(new.sheets, 'family_index');
R = compare_phase_catalogs(new, ref, quiet);
ok = chk(ok, ~R.ok && ~R.familiesCompared, 'a missing family map means "not established", and that is not a match');
% another problem (engine) on the same grid is refused
new = ref;  new.thruster.isp_s = new.thruster.isp_s + 100;
ok = chk(ok, throws(@() compare_phase_catalogs(new, ref, quiet)), 'a catalog of another engine is refused, not compared');
% the count of agreeing cells is a count of CELLS (a cell with two differences is one cell)
new = ref;  s = new.sheets(1);  s.tf_nd(3, 5) = s.tf_nd(3, 5) + 1e-3;  k = s.entry_index(3, 5);  s.z8(1:7, k) = 1.01*s.z8(1:7, k);  new.sheets(1) = s;
R = compare_phase_catalogs(new, ref, quiet);
ok = chk(ok, R.nAgree == 575 && R.nTfOver == 1 && R.nZOver == 1, sprintf('one cell wrong twice: %d of 576 agree (575)', R.nAgree));
% the costates are compared on their own: a costate change is not diluted by t_f
new = ref;  s = new.sheets(1);  k = s.entry_index(9, 9);  s.z8(1:7, k) = s.z8(1:7, k)*(1 + 5e-6);  new.sheets(1) = s;
R = compare_phase_catalogs(new, ref, quiet);
ok = chk(ok, R.nZOver == 1, 'a 5e-6 relative costate change is caught (costates compared without t_f in the norm)');

% ---- the script's audit is THE FINAL CATALOG'S audit, complete ---------------
Hs = reproduce_library_70mN('localfunctions');
aud = fullfile(tmpRoot(), 'audbind');  r1 = fullfile(aud, 'round_01');  r2 = fullfile(aud, 'round_02');  mkdir(r1);  mkdir(r2);
A = struct('nOk', 10, 'nBad', 0); save(fullfile(r1, 'audit_70mN.mat'), 'A'); %#ok<NASGU>
A = struct('nOk', 576, 'nBad', 0); save(fullfile(r2, 'audit_70mN.mat'), 'A');   % the chain names it by ITS tag, not the driver's %#ok<NASGU>
st = struct('rounds', struct('dir', {r1, r2})); save(fullfile(aud, 'torus_state.mat'), 'st'); %#ok<NASGU>
[Af, why, how] = Hs.finalAudit(aud, 576, 'abc');
ok = chk(ok, ~isempty(Af) && Af.nOk == 576, ['the audit of the LAST round in the driver''s state is the one read: ' why]);
ok = chk(ok, strcmp(how, 'count only'), 'an audit made before content keys existed is accepted, and SAID to be bound by count only');
[Af, why] = Hs.finalAudit(aud, 500, 'abc');
ok = chk(ok, isempty(Af), ['an audit that does not cover every catalog entry is not accepted: ' why]);
A = struct('nOk', 576, 'nBad', 0, 'contentKey', 'abc'); save(fullfile(r2, 'audit_70mN.mat'), 'A'); %#ok<NASGU>
[Af, why, how] = Hs.finalAudit(aud, 576, 'abc');
ok = chk(ok, ~isempty(Af) && strcmp(how, 'content key'), ['an audit carrying the catalog''s own key is bound by it: ' why]);
[Af, why] = Hs.finalAudit(aud, 576, 'xyz');
ok = chk(ok, isempty(Af) && contains(why, 'content'), ['the right count of ANOTHER catalog''s entries is not accepted: ' why]);
delete(fullfile(r2, 'audit_70mN.mat'));
[Af, why] = Hs.finalAudit(aud, 576, 'abc');
ok = chk(ok, isempty(Af), ['no fallback to an EARLIER round''s audit: ' why]);
rmdir(aud, 's');

% ---- THE SAME PHASES IN ANOTHER ORDER are the same grid ---------------------
% (the first full rebuild, 2026-09-19: the record lists its arrival phases from
% the anchor's, 0.0754 ... 0.9921 0.0337; the driver sorts them, 0.0337 first.
% Cells are matched by PHASE, never by index.)
new = ref;  s = new.sheets(1);  p = [24, 1:23];                 % a cyclic shift of the columns
s.sA_frac = s.sA_frac(p);
for f = {'has_solution', 'tf_nd', 'entry_index', 'family_index'}, s.(f{1}) = s.(f{1})(:, p, :); end
new.sheets(1) = s;
R = compare_phase_catalogs(new, ref, quiet);
ok = chk(ok, R.ok && R.nAgree == 576, sprintf('columns cyclically shifted: still the same library (%d of 576 agree)', R.nAgree));
s2 = s;  pd = [13:24, 1:12];  s2.sD_frac = s2.sD_frac(pd);      % and the rows too
for f = {'has_solution', 'tf_nd', 'entry_index', 'family_index'}, s2.(f{1}) = s2.(f{1})(pd, :, :); end
new.sheets(1) = s2;
R = compare_phase_catalogs(new, ref, quiet);
ok = chk(ok, R.ok && R.nAgree == 576, 'rows AND columns reordered: still the same library');
s3 = s;  s3.tf_nd(4, 7) = s3.tf_nd(4, 7) + 1e-3;  new.sheets(1) = s3;   % a difference planted in the SHIFTED catalog ...
R = compare_phase_catalogs(new, ref, quiet);
ok = chk(ok, R.nTfOver == 1 && any(abs(R.cells.sA - s3.sA_frac(7)) < 1e-12 & R.cells.iD == 4), ...
         '... is reported at its PHASE (the reference''s indices and phases name the cell)');

% ---- a different root is FASTER or SLOWER, and the report says which --------
new = ref;  s = new.sheets(1);  s.tf_nd(2, 2) = s.tf_nd(2, 2) - 0.5;  s.tf_nd(8, 8) = s.tf_nd(8, 8) - 0.2;  s.tf_nd(9, 9) = s.tf_nd(9, 9) + 0.3;  new.sheets(1) = s;
R = compare_phase_catalogs(new, ref, quiet);
ok = chk(ok, R.nNewFaster == 2 && R.nNewSlower == 1 && ~R.noWorse, sprintf('two cells faster, one slower: counted %d / %d', R.nNewFaster, R.nNewSlower));
s.tf_nd(9, 9) = ref.sheets(1).tf_nd(9, 9);  new.sheets(1) = s;
R = compare_phase_catalogs(new, ref, quiet);
ok = chk(ok, ~R.ok && R.noWorse && R.nNewFaster == 2 && R.nNewSlower == 0, 'only faster cells: not the same library, but NO WORSE than the reference in any cell');

% ---- A REFINED GRID is compared on the cells the two grids SHARE ---------------
% (a 48 x 48 rebuild against the 24 x 24 record: every record phase is on the
% finer grid, so the record's cells are the common cells)
coarse = ref;  s = coarse.sheets(1);  pick2 = 1:2:24;
s.sD_frac = s.sD_frac(pick2);  s.sA_frac = s.sA_frac(pick2);
for f = {'has_solution', 'tf_nd', 'entry_index', 'family_index'}, s.(f{1}) = s.(f{1})(pick2, pick2, :); end
coarse.sheets(1) = s;
R = compare_phase_catalogs(ref, coarse, quiet);                  % the FINER catalog is the new one
ok = chk(ok, R.ok && ~R.sameGrid && R.nBoth == 144 && R.nAgree == 144, sprintf('finer new vs coarser reference: %d common cells, all agree', R.nBoth));
ok = chk(ok, R.nNewOffGrid == 576 - 144 && R.nRefOffGrid == 0, sprintf('%d cells of the finer catalog lie off the reference grid and are not judged', R.nNewOffGrid));
R = compare_phase_catalogs(coarse, ref, quiet);                  % and the other way round
ok = chk(ok, R.ok && R.nBoth == 144 && R.nRefOffGrid == 576 - 144 && R.nOnlyRef == 0, 'coarser new vs finer reference: the reference''s extra cells are off-grid, NOT "missing"');
fine = ref;  k = fine.sheets(1).entry_index(3, 3);  fine.sheets(1).tf_nd(3, 3) = fine.sheets(1).tf_nd(3, 3) - 0.4;   % (3,3) is on the coarse grid
R = compare_phase_catalogs(fine, coarse, quiet);
ok = chk(ok, ~R.ok && R.noWorse && R.nNewFaster == 1, 'a faster root on a shared cell is still found on a refined grid');

% ---- the script takes the resolution as an option ------------------------------
o = reproduce_library_70mN(struct('nD', 48, 'nA', 48, 'print', false));
ok = chk(ok, numel(o.spec.sD) == 48 && numel(o.spec.sA) == 48 && issorted(o.spec.sA) && contains(o.spec.outDir, '48x48'), 'nD = nA = 48: two 48-phase lists, and a folder of its own');
rs = ref.sheets(1).sA_frac;  d = min(abs(mod(o.spec.sA(:).' - rs(:) + 0.5, 1) - 0.5), [], 2);
ok = chk(ok, max(d) < 1e-12 && max(abs(o.spec.sD(1:2:end) - ref.sheets(1).sD_frac)) < 1e-12, 'every phase of the 24 x 24 record is on the 48 x 48 grid');
o24 = reproduce_library_70mN(struct('print', false));
ok = chk(ok, numel(o24.spec.sD) == 24 && endsWith(o24.spec.outDir, 'reproduce_70mN_24x24'), 'the default is unchanged: 24 x 24, the same folder as before');

% ---- the time estimate is COMPUTED from the measured rates ----------------------
% (the plan used to print "6-10 h"; the 24 x 24 rebuild took 24 h 12 min)
e24 = Hs.estimateHours(24, 24, 4, true);  e48 = Hs.estimateHours(48, 48, 4, true);
ok = chk(ok, abs(e24 - 24.2) < 2, sprintf('24 x 24 with 4 workers: %.1f h (measured 24.2)', e24));
ok = chk(ok, e48 > 80 && e48 < 110 && Hs.estimateHours(48, 48, 8, true) < e48 && Hs.estimateHours(24, 24, 4, false) > e24 + 10, ...
         sprintf('48 x 48: %.0f h (about four days); more workers shorten it; walking the arcs adds most of a day', e48));

% ---- the family pairing is OPTIMAL, and status codes are not families ---------
% overlap [6 5; 5 0]: greedy pairs the 6 and leaves 10 cells "different"; the best pairing matches 10 and leaves 6
new = ref;  s = new.sheets(1);  fr = ones(24, 'int8');  fn_ = ones(24, 'int8');  idx = find(s.has_solution);
fr(idx(1:11)) = 1;  fr(idx(12:16)) = 2;  fr(idx(17:end)) = 3;
fn_(idx(1:6)) = 1;  fn_(idx(7:11)) = 2;  fn_(idx(12:16)) = 1;  fn_(idx(17:end)) = 3;
refF = ref;  refF.sheets(1).family_index = fr;  s.family_index = fn_;  new.sheets(1) = s;
R = compare_phase_catalogs(new, refF, quiet);
ok = chk(ok, R.nFamilyDiff == 6, sprintf('overlap [6 5; 5 0]: %d cells differ under the OPTIMAL pairing (6; greedy said 10)', R.nFamilyDiff));
s2 = ref.sheets(1);  f2 = double(s2.family_index);  known = find(f2 >= 1 & s2.has_solution, 3);
newS = ref;  g2 = s2.family_index;  g2(known(1)) = -2;  newS.sheets(1).family_index = g2;     % a known family relabelled "unidentified"
R = compare_phase_catalogs(newS, ref, quiet);
ok = chk(ok, R.nFamilyDiff == 1, 'a status code (-2 unidentified) is compared as a status, never paired with a family');

% ---- another grid is refused ----------------------------------------------
new = ref;  new.sheets(1).sA_frac(4) = new.sheets(1).sA_frac(4) + 1e-4;
ok = chk(ok, throws(@() compare_phase_catalogs(new, ref, quiet)), 'a catalog on another grid is refused, not compared');

% ---- adopt_walked_arcs -----------------------------------------------------
tmp = fullfile(tempdir, sprintf('adopt_%s', char(java.util.UUID.randomUUID())));
src = fullfile(tmp, 'src');  dst = fullfile(tmp, 'arcs');  mkdir(src);  mkdir(dst);
cleanup = onCleanup(@() rmdir(tmp, 's'));
for nm = {'a_dn', 'a_up'}
    A = struct('q', rand(1, 3)); %#ok<NASGU>
    save(fullfile(src, sprintf('arrival_arc_%s_long.mat', nm{1})), 'A');
end
n = adopt_walked_arcs(src, {'a'}, dst, 'T', quiet);
ok = chk(ok, n == 2 && isfile(fullfile(dst, 'arrival_arc_T_a_dn_long.mat')) && isfile(fullfile(dst, 'arrival_arc_T_a_up_long.mat')), ...
         'both directions copied under the driver''s tagged name');
ok = chk(ok, adopt_walked_arcs(src, {'a'}, dst, 'T', quiet) == 0, 'an identical copy already there is left alone');
A = struct('q', 1); %#ok<NASGU>
save(fullfile(dst, 'arrival_arc_T_a_dn_long.mat'), 'A');
ok = chk(ok, throws(@() adopt_walked_arcs(src, {'a'}, dst, 'T', quiet)), 'a DIFFERENT file under that name is never overwritten');
ok = chk(ok, throws(@() adopt_walked_arcs(src, {'nosuch'}, dst, 'T', quiet)), 'a missing walked arc is named, not skipped');

% ---- the script's PLAN mode: says what it would do, touches nothing -------
planDir = fullfile(tmp, 'plan_campaign');
o = reproduce_library_70mN(struct('outDir', planDir, 'print', false));      % .go is false by default
ok = chk(ok, strcmp(o.state, 'planned') && ~isfolder(planDir), 'plan mode: state planned, no folder created');
ok = chk(ok, size(o.spec.anchors, 1) == 5 && all(cellfun(@isfile, o.spec.anchors(:, 2))), 'plan mode: the five families, every anchor file on disk');
ok = chk(ok, o.nArcsFound == 10 && o.spec.discover == false && numel(o.spec.sD) == 24 && numel(o.spec.sA) == 24, ...
         'plan mode: ten walked arcs found, discovery off, the 24 x 24 grid');
ok = chk(ok, isfield(o, 'nRegistryRoots') && o.nRegistryRoots >= 1, 'plan mode: the record''s direct-found seed roots are found (they seed the sheet too)');

if ok, fprintf('test_reproduce_library: ALL PASS\n'); else, fprintf('test_reproduce_library: FAIL\n'); end
end

function tf = throws(f)
% THROWS  Does calling f throw?  INPUTS: f (handle).  OUTPUTS: tf.
tf = false;
try, f(); catch, tf = true; end
end

function ok = chk(ok, c, msg)
% CHK  Print one PASS/FAIL line and fold it into ok.  INPUTS: ok; c; msg.
% OUTPUTS: ok.
if c, fprintf('  PASS  %s\n', msg); else, fprintf('  FAIL  %s\n', msg); end
ok = ok && c;
end

function d = tmpRoot()
% TMPROOT  A fresh temporary folder name.  INPUTS: none.  OUTPUTS: d.
d = fullfile(tempdir, sprintf('reprolib_%s', char(java.util.UUID.randomUUID())));
end
