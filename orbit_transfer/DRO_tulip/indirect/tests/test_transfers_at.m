function ok = test_transfers_at()
% TEST_TRANSFERS_AT  Every transfer at a phase pair, ranked by flight time:
% the primary and the alternatives in its cell and off the grid; the wrap
% s_D = 1 - 1e-9 finds s_D = 0; a catalog without the status layer reports
% its primary as 'not annotated'; Delta-V follows the catalog's formula.
%
% INPUTS:  none
% OUTPUTS: ok [logical]  every check passed (prints PASS/FAIL per check)

ok = true;
here = fileparts(fileparts(mfilename('fullpath')));
addpath(here);
L = load(fullfile(here, 'results', 'library_70mN_24x24_final', 'costate_catalog_dro_tulip_70mN.mat'));  c = L.(char(fieldnames(L)));
sh = c.sheets(1);  k = sh.entry_index(1, 1);  tS = c.constants.tStar_s/86400;
T = transfers_at(c, sh.sD_frac(1), sh.sA_frac(1));
ok = chk(ok, numel(T) == 1 && strcmp(T.kind, 'primary') && isnan(T.status) && strcmp(T.statusName, 'not annotated'), ...
         'REVIEW FOCUS 2: no layer -> the primary, not annotated');
Tmax = (c.rungs_N(1)/c.thruster.m0_kg)*c.constants.tStar_s^2/(c.constants.lStar_km*1000);
mf = 1 - Tmax*sh.z8(8, k)/c.thruster.c_nd;
ok = chk(ok, abs(T.dvKms - c.thruster.c_nd*log(1/mf)*c.constants.lStar_km/c.constants.tStar_s) < 1e-12, 'Delta-V from the catalog''s formula');
c.sheets(1).status = int8(4*sh.has_solution);  c.sheets(1).status_reason = repmat({'full stack passed'}, 1, size(sh.z8, 2));
c.sheets(1).junctions = repmat({[]}, 1, size(sh.z8, 2));
alt = @(z8, st, sDv) struct('sD', sDv, 'sA', sh.sA_frac(1), 'iD', 1, 'iA', 1, 'z8', z8, 'junctions', [], 'tf_nd', z8(8), ...
        'status', st, 'status_reason', 'r', 'inferred', false, 'conj', 0, 'conjVerdict', 'FAIL', 'minLamV', 1, 'minQmt', 1, ...
        'dimS', 1, 'h6Margin', 1, 'liftMargin', 1, 'flyKm', 0, 'flyVms', 0, 'source', 's', 'sheet', 1);
zS = sh.z8(:, k);  zF = zS;  zF(8) = zS(8) - 0.01;  zL = zS;  zL(8) = zS(8) + 0.02;
c.alternatives = [alt(zL, 2, sh.sD_frac(1)), alt(zF, 3, sh.sD_frac(1)), alt(zL, 1, 0.5)];
T = transfers_at(c, 1 - 1e-9, sh.sA_frac(1));
ok = chk(ok, numel(T) == 3 && issorted([T.tfDays]) && strcmp(T(2).kind, 'primary') && T(1).status == 3 && T(3).status == 2, ...
         'REVIEW FOCUS 4: the wrap finds s_D = 0; ranked by t_f: necessary-only, primary, conjugate');
ok = chk(ok, strcmp(T(3).statusName, 'conjugate point found'), 'status names from status_key');
if ok, fprintf('test_transfers_at: ALL PASS\n'); else, fprintf('test_transfers_at: FAIL\n'); end
end

function ok = chk(ok, c, msg)
% CHK  Print one PASS/FAIL line and fold it into ok.  INPUTS: ok; c; msg.
% OUTPUTS: ok.
if c, fprintf('  PASS  %s\n', msg); else, fprintf('  FAIL  %s\n', msg); end
ok = ok && c;
end
