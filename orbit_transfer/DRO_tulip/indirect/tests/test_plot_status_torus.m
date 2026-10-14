function ok = test_plot_status_torus()
% TEST_PLOT_STATUS_TORUS  plot_status_torus counts alternatives per grid
% cell, flags status-2 cells, keeps off-grid rows apart, writes a PNG, and
% accounts for every alternative of the library of record.
%
% INPUTS:  none
% OUTPUTS: ok [logical]  every check passed (prints PASS/FAIL per check)

ok = true;
here = fileparts(fileparts(mfilename('fullpath')));
addpath(here, fullfile(fileparts(fileparts(here)), 'costate_common'));
c = synth_catalog();
png = fullfile(tempdir, 'plot_status_torus_test.png');
if isfile(png), delete(png); end
H = plot_status_torus(c, struct('out', png, 'visible', false));
close(H.fig);
ok = chk(ok, isequal(H.nAltPerCell, [0 2 0; 0 0 1]), 'alternatives per cell [0 2 0; 0 0 1]');
ok = chk(ok, isequal(H.conjCell, logical([0 1 0; 0 0 1])), 'status-2 cells');
ok = chk(ok, H.nOffGrid == 1 && H.nOffGridConj == 1, 'one off-grid row, one off-grid status 2');
d = dir(png);
ok = chk(ok, ~isempty(d) && d.bytes > 10e3, 'PNG written and larger than 10 kB');
% a status-3 alternative alone in a cell is NOT a conjugate cell
c3 = c;  c3.alternatives(end+1) = struct('sD', 0.5, 'sA', 0, 'iD', 2, 'iA', 1, 'status', 3);
H3 = plot_status_torus(c3, struct('visible', false));  close(H3.fig);
ok = chk(ok, H3.nAltPerCell(2, 1) == 1 && ~H3.conjCell(2, 1) && isequal(H3.conjCell(:), H.conjCell(:)), 'only status 2 marks a conjugate cell');
% no alternatives
c0 = rmfield(c, 'alternatives');
H0 = plot_status_torus(c0, struct('visible', false));  close(H0.fig);
ok = chk(ok, all(H0.nAltPerCell(:) == 0) && ~any(H0.conjCell(:)) && H0.nOffGrid == 0, 'catalog without alternatives works');
% V2
v2 = load(fullfile(here, 'results', 'library_70mN_24x48_v2', 'costate_catalog_dro_tulip_70mN.mat'));  v = v2.(char(fieldnames(v2)));
Hv = plot_status_torus(v, struct('visible', false));  close(Hv.fig);
fprintf('  V2: per-cell max %d, off-grid %d, off-grid status 2 %d\n', max(Hv.nAltPerCell(:)), Hv.nOffGrid, Hv.nOffGridConj);
ok = chk(ok, sum(Hv.nAltPerCell(:)) + Hv.nOffGrid == 656, 'V2: on-grid + off-grid = 656');
if ok, fprintf('test_plot_status_torus: ALL PASS\n'); else, fprintf('test_plot_status_torus: FAIL\n'); end
end

function c = synth_catalog()
% SYNTH_CATALOG  2 x 3 grid catalog with four alternatives.  INPUTS: none.
% OUTPUTS: c [struct].
z8 = zeros(8, 6);  z8(8, :) = 5 + (1:6);
sh = struct('sD_frac', [0 0.5], 'sA_frac', [0 1/3 2/3], 'has_solution', true(2, 3), ...
            'entry_index', reshape(1:6, 2, 3), 'z8', z8);
c = struct('sheets', sh, 'constants', struct('tStar_s', 382981.289129055));
mk = @(sD, sA, iD, iA, st) struct('sD', sD, 'sA', sA, 'iD', iD, 'iA', iA, 'status', st);
c.alternatives = [mk(0, 1/3, 1, 2, 2), mk(0, 1/3, 1, 2, 3), mk(0.5, 2/3, 2, 3, 2), mk(0.17, 0.3, 0, 0, 2)];
end

function ok = chk(ok, c, msg)
% CHK  Print one PASS/FAIL line and fold it into ok.  INPUTS: ok; c; msg.
% OUTPUTS: ok.
if c, fprintf('  PASS  %s\n', msg); else, fprintf('  FAIL  %s\n', msg); end
ok = ok && c;
end
