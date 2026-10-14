function ok = test_fly_catalog_entry()
% TEST_FLY_CATALOG_ENTRY  A library entry flies from its junctions and from
% its z8 to the catalog's own arrival state (< 1 km / 0.1 m/s), the thrust
% direction is the one the equations of motion apply, the summary and the
% 3D figure accept the output, a cell index selects the primary, a row with
% no junctions falls back to z8 and says so, and 'junctions' without any
% refuses by name.
%
% INPUTS:  none
% OUTPUTS: ok [logical]  every check passed (prints PASS/FAIL per check)

ok = true;
here = fileparts(fileparts(mfilename('fullpath')));
addpath(here, fullfile(fileparts(fileparts(here)), 'costate_common'));
L = load(fullfile(here, 'results', 'library_70mN_24x48_v2', 'costate_catalog_dro_tulip_70mN.mat'));  c = L.(char(fieldnames(L)));
sh = c.sheets(1);  iD = 3;  iA = 7;
T = transfers_at(c, sh.sD_frac(iD), sh.sA_frac(iA));  prim = T(strcmp({T.kind}, 'primary'));
Fj = fly_catalog_entry(c, prim);
Fz = fly_catalog_entry(c, prim, struct('method', 'z8'));
ok = chk(ok, strcmp(Fj.method, 'junctions') && Fj.flyKm < 1 && Fj.flyVms < 0.1 && max([Fj.defectKm 0]) < 0.05, sprintf('junction flight arrives (%.3f km, %.4f m/s), defects %.1e km', Fj.flyKm, Fj.flyVms, max([Fj.defectKm 0])));
ok = chk(ok, strcmp(Fz.method, 'z8') && Fz.flyKm < 1 && Fz.flyVms < 0.1 && isempty(Fz.defectKm), sprintf('z8 flight arrives (%.3f km, %.4f m/s)', Fz.flyKm, Fz.flyVms));
ok = chk(ok, abs(Fj.tfDays - prim.tfDays) < 1e-9 && abs(Fj.dvKms - prim.dvKms) < 1e-9 && Fj.status == 4 && strcmp(Fj.statusName, 'sufficient'), 'numbers and status carried from the row');
B = Fj.B;
for k = round(linspace(2, numel(Fj.t) - 1, 5))
    dyT = pumpkyn.cr3bp.tfMinEoM(Fj.t(k), Fj.Y(k, :).', B.Tnd, B.cnd, B.mu);
    dy0 = pumpkyn.cr3bp.tfMinEoM(Fj.t(k), Fj.Y(k, :).', 0, B.cnd, B.mu);
    a = dyT(4:6) - dy0(4:6);  a = a/sqrt(sum(a.^2));
    ok = chk(ok, sqrt(sum((a(:) - Fj.u(k, :).').^2)) < 1e-9, sprintf('thrust direction = the EoM''s applied acceleration (sample %d)', k));
end
for k = round(linspace(2, numel(Fz.t) - 1, 5))
    dyT = pumpkyn.cr3bp.tfMinEoM(Fz.t(k), Fz.Y(k, :).', B.Tnd, B.cnd, B.mu);
    dy0 = pumpkyn.cr3bp.tfMinEoM(Fz.t(k), Fz.Y(k, :).', 0, B.cnd, B.mu);
    a = dyT(4:6) - dy0(4:6);  a = a/sqrt(sum(a.^2));
    ok = chk(ok, sqrt(sum((a(:) - Fz.u(k, :).').^2)) < 1e-9, sprintf('z8 flight: thrust direction = the EoM''s applied acceleration (sample %d)', k));
end
lines = print_transfer_summary(Fj);
ok = chk(ok, iscellstr(lines) || ischar(lines) || isstring(lines), 'print_transfer_summary accepts F');
P = plot_transfer_3d(Fj, B, struct('visible', false, 'flight', struct('t', Fj.t, 'Y', Fj.Y)));
ok = chk(ok, ~isempty(P), 'plot_transfer_3d accepts (F, F.B)');  close all
Fc = fly_catalog_entry(c, struct('iD', iD, 'iA', iA));
ok = chk(ok, isequal(Fc.z8, prim.z8) && strcmp(Fc.method, 'junctions'), 'a cell index selects the primary (with its junctions)');
rowNoJ = prim;  rowNoJ.junctions = [];
Fn = fly_catalog_entry(c, rowNoJ);
ok = chk(ok, strcmp(Fn.method, 'z8') && contains(Fn.reasonMethod, 'no junctions'), 'REVIEW FOCUS 3: a row without junctions falls back to z8 and says so');
try, fly_catalog_entry(c, rowNoJ, struct('method', 'junctions'));  bad = false; catch ME, bad = strcmp(ME.identifier, 'fly_catalog_entry:noJunctions'); end
ok = chk(ok, bad, 'explicit ''junctions'' with none refuses by name');
T1 = transfers_at(c, 0, sh.sA_frac(iA));  prim1 = T1(strcmp({T1.kind}, 'primary'));
Fw = fly_catalog_entry(c, struct('z8', prim1.z8, 'sD', 1 - 1e-12, 'sA', prim1.sA), struct('method', 'z8'));
ok = chk(ok, Fw.flyKm < 1 && Fw.flyVms < 0.1, 'REVIEW FOCUS 2: a z8-with-phases entry at the wrap flies to the same arrival');
rowBad = prim;  rowBad.junctions = ones(13, 3);
Fb = fly_catalog_entry(c, rowBad);
ok = chk(ok, strcmp(Fb.method, 'z8') && contains(Fb.reasonMethod, 'unexpected size'), 'a junctions array of the wrong shape falls back to z8 and says so');
try, fly_catalog_entry(c, rowBad, struct('method', 'junctions'));  bad = false; catch ME, bad = strcmp(ME.identifier, 'fly_catalog_entry:junctionsShape'); end
ok = chk(ok, bad, 'explicit ''junctions'' with the wrong shape refuses by name');
try, fly_catalog_entry(c, struct('z8', prim.z8));  bad = false; catch ME, bad = strcmp(ME.identifier, 'fly_catalog_entry:entry'); end
ok = chk(ok, bad, 'a bare z8 (no sD/sA) refuses by name');
if ok, fprintf('test_fly_catalog_entry: ALL PASS\n'); else, fprintf('test_fly_catalog_entry: FAIL\n'); end
end

function ok = chk(ok, c, msg)
% CHK  Print one PASS/FAIL line and fold it into ok.  INPUTS: ok; c; msg.
% OUTPUTS: ok.
if c, fprintf('  PASS  %s\n', msg); else, fprintf('  FAIL  %s\n', msg); end
ok = ok && c;
end
