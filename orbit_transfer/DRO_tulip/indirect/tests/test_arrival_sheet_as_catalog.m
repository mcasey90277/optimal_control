function ok = test_arrival_sheet_as_catalog()
% TEST_ARRIVAL_SHEET_AS_CATALOG  An arrival sheet (one departure phase) wrapped
% as a one-row phase catalog: every certified column becomes an entry with the
% same eight numbers, the problem identity is the reference catalog's, columns
% won by different arcs get different family codes, and the interpolator reads it.
%
% INPUTS:  none
% OUTPUTS: ok [logical]  every check passed (prints PASS/FAIL per check)

ok = true;
here = fileparts(fileparts(mfilename('fullpath')));
addpath(here, fullfile(fileparts(fileparts(here)), 'costate_common'));
L = load(fullfile(here, 'results', 'library_70mN_24x24_final', 'costate_catalog_dro_tulip_70mN.mat'));  fn = fieldnames(L);  ref = L.(fn{1});
S = struct('sA', [0.1 0.35 0.6 0.85], 'TF', [17 18 NaN 19], 'Z8', [repmat((1:7).', 1, 4); 4.0 4.2 NaN 4.4], ...
           'arcs', {{'arrival_arc_70mN_anchor_dn_long.mat', 'arrival_arc_70mN_fast2_up_long.mat'}});
S.Z8(1, :) = [1 1.1 NaN 5];
S.cand = {struct('ok', true, 'tfDays', 17, 'arc', 1), struct('ok', true, 'tfDays', 18, 'arc', 1), ...
          struct('ok', false, 'tfDays', NaN, 'arc', 1), struct('ok', true, 'tfDays', 19, 'arc', 2)};
c = arrival_sheet_as_catalog(S, ref, 0);
sh = c.sheets(1);
ok = chk(ok, isequal(size(sh.has_solution), [1 4]) && isequal(sh.has_solution, [true true false true]) && c.n_entries == 3, 'one row, an entry per CERTIFIED column');
ok = chk(ok, isequal(sh.z8(:, sh.entry_index(1, 4)), S.Z8(:, 4)) && sh.tf_nd(1, 2) == 4.2 && isnan(sh.tf_nd(1, 3)), 'the same eight numbers; t_f in ND from z8(8); NaN where nothing certified');
ok = chk(ok, sh.family_index(1, 1) == sh.family_index(1, 2) && sh.family_index(1, 4) ~= sh.family_index(1, 1) && sh.family_index(1, 3) == 0, 'columns won by different arcs carry different family codes');
ok = chk(ok, strcmp(c.families.labels{sh.family_index(1, 4)}, 'fast2'), 'and the labels are the arcs'' names');
ok = chk(ok, isequal(c.thruster, ref.thruster) && isequal(c.constants, ref.constants) && sh.Np == ref.sheets(1).Np && sh.sD_frac == 0, 'the problem identity is the reference catalog''s');
[z, info] = phase_catalog_interp(sh, 0, 0.2);
ok = chk(ok, strcmp(info.tier, 'linear') && abs(z(1) - 1.04) < 1e-12, 'the interpolator reads it (one departure line: the blend is along arrival phase)');
[~, info] = phase_catalog_interp(sh, 0, 0.9);
ok = chk(ok, strcmp(info.tier, 'nearest') && contains(info.reason, 'famil'), 'and refuses to blend across the hand-over between two arcs');
if ok, fprintf('test_arrival_sheet_as_catalog: ALL PASS\n'); else, fprintf('test_arrival_sheet_as_catalog: FAIL\n'); end
end

function ok = chk(ok, c, msg)
% CHK  Print one PASS/FAIL line and fold it into ok.  INPUTS: ok; c; msg.
% OUTPUTS: ok.
if c, fprintf('  PASS  %s\n', msg); else, fprintf('  FAIL  %s\n', msg); end
ok = ok && c;
end
