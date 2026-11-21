function e = dro_elfo_entry(catalog, ks, iD, iA, kr)
%% Purpose:
%
%   ONE ROW of the DRO -> ELFO library functions: the entry stored at
%   (sheet ks, departure cell iD, arrival cell iA, rung kr), with every
%   field dro_elfo_fly needs and the derived units. dro_elfo_transfers_at
%   and dro_elfo_extremes both build their rows here, so a row from either
%   carries the same fields and flies the same way.
%
%   Derived with the catalog's own formulas (catalog.derive):
%       Tmax_nd = (T_N/m0_kg)*tStar_s^2/(lStar_km*1000)
%       mf      = 1 - Tmax_nd*tf_nd/c_nd        (all-burn minimum time)
%       dV      = c_nd*log(1/mf)*lStar_km/tStar_s
%
%   VERDICTS: conj_pass as recorded (1 pass, 0 refuted, -1 not run; -1
%   when the sheet carries no verdicts) and the hypothesis-gate values
%   (NaN when not gated), as dro_elfo_pick returns them.
%
%% Inputs:
%
%  catalog                  struct                  The catalog struct
%
%  ks                       double                  Sheet index [1 x 1]
%
%  iD, iA                   double                  Departure / arrival
%                                                   phase cell [1 x 1]
%
%  kr                       double                  Rung index [1 x 1]
%
%% Outputs:
%
%  e                        struct                  [1 x 1] .thrustN .sD .sA
%                                                   (the cell's phase
%                                                   fractions) .iD .iA .kr
%                                                   .tfDays .tf_nd .dvKms
%                                                   .propellantKg
%                                                   .finalMassKg .conj_pass
%                                                   .gate_min_lamv
%                                                   .gate_min_qmt .gate_dimS
%                                                   .z8 [8 x 1] .raanDeg
%                                                   .droTau .sheet
%                                                   .sheet_key .dep_family
%                                                   .dep_params .arr_family
%                                                   .arr_params .label
%                                                   [1 x n char]
%
%% Revision History:
%  M. Casey                                                   (c) 11/21/2026
%  Copyright Coorbital Inc.
%% ------------------------ Begin Code Sequence ---------------------------

if nargin == 0
    here = fileparts(fileparts(mfilename('fullpath')));
    f = fullfile(here, 'results', 'pilot_chain', 'costate_catalog_dro_elfo.mat');
    if ~isfile(f), fprintf('dro_elfo_entry: self-demo needs %s\n', f);  e = [];  return, end
    L = load(f);  fn = fieldnames(L);  catalog = L.(fn{1});
    [iD, iA] = find(catalog.sheets(1).has_solution(:, :, 1), 1);
    e = dro_elfo_entry(catalog, 1, iD, iA, 1);
    fprintf('dro_elfo_entry: %s\n', e.label);
    return
end
sh = catalog.sheets(ks);
if ~sh.has_solution(iD, iA, kr)
    error('dro_elfo_entry:noEntry', 'sheet %d cell (%d,%d) holds no entry at rung %d', ks, iD, iA, kr);
end
tStar = catalog.constants.tStar_s;
lStar = catalog.constants.lStar_km;
m0kg  = catalog.thruster.m0_kg;
cnd   = catalog.thruster.c_nd;
TN    = catalog.rungs_N(kr);
z8    = sh.z8(:, sh.entry_index(iD, iA, kr));
tf    = z8(8);
Tnd   = (TN/m0kg)*tStar^2/(lStar*1000);
mf    = 1 - Tnd*tf/cnd;

e = struct('thrustN', TN, 'sD', sh.sD_frac(iD), 'sA', sh.sA_frac(iA), 'iD', iD, 'iA', iA, 'kr', kr, ...
           'tfDays', tf*tStar/86400, 'tf_nd', tf, 'dvKms', cnd*log(1/mf)*lStar/tStar, ...
           'propellantKg', (1 - mf)*m0kg, 'finalMassKg', mf*m0kg, 'conj_pass', -1, ...
           'gate_min_lamv', NaN, 'gate_min_qmt', NaN, 'gate_dimS', NaN, 'z8', z8, ...
           'raanDeg', sh.sheet_key.raanDeg, 'droTau', sh.sheet_key.depSize, 'sheet', ks, ...
           'sheet_key', sh.sheet_key, 'dep_family', sh.dep_family, 'dep_params', sh.dep_params, ...
           'arr_family', sh.arr_family, 'arr_params', sh.arr_params, 'label', '');
if isfield(sh, 'conj_pass') && ~isempty(sh.conj_pass), e.conj_pass = double(sh.conj_pass(iD, iA, kr)); end
for g = {'gate_min_lamv', 'gate_min_qmt', 'gate_dimS'}
    if isfield(sh, g{1}) && ~isempty(sh.(g{1})), e.(g{1}) = double(sh.(g{1})(iD, iA, kr)); end
end
verdicts = {'refuted', 'pass'};
if e.conj_pass == 0 || e.conj_pass == 1, vtxt = verdicts{e.conj_pass + 1}; else, vtxt = 'not run'; end
e.label = sprintf(['DRO %.4g ND -> ELFO raan %g deg, %.3g N, sD %.4f sA %.4f: t_f %.3f d, ', ...
                   'dV %.4f km/s, %.3f kg, conj %s'], e.droTau, e.raanDeg, TN, e.sD, e.sA, ...
                  e.tfDays, e.dvKms, e.propellantKg, vtxt);
end
