function D = dro_elfo_describe(catalog)
%% Purpose:
%
%   FACT SHEET of a DRO -> ELFO costate catalog, printed and returned: the
%   orbits (DRO period in ND and days, ELFO elements and period), the
%   sheets (one per DRO size and ELFO plane orientation raanDeg), the
%   thrust rungs, Isp, mass and phase grid; per sheet the entries, the
%   cells carrying a full ladder (an entry at every rung) and the cells
%   with an entry at the lowest rung; per rung over all sheets the
%   entries, the flight-time range and median in days and the Delta-V
%   range; and per sheet the conjugate-test counts (pass / refuted / not
%   run). Every count is taken from the sheet arrays (has_solution,
%   conj_pass); Delta-V follows catalog.derive.
%
%% Inputs:
%
%  catalog                  struct | char           The catalog struct, or
%                                                   the path of its .mat
%                                                   (one variable)
%
%% Outputs:
%
%  D                        struct                  [1 x 1] .name .nSheets
%                                                   .nEntries (sum over
%                                                   sheets) .rungsN
%                                                   [1 x nR] .ispS .m0kg
%                                                   .nD .nA (phase grid)
%                                                   .droTau [1 x nS]
%                                                   .droPeriodDays [1 x nS]
%                                                   .raanDeg [1 x nS]
%                                                   .elfoOe [1 x 4]
%                                                   .elfoPeriodNd
%                                                   .elfoPeriodDays
%                                                   .sheetEntries [nS x 1]
%                                                   .sheetFullLadders
%                                                   [nS x 1]
%                                                   .sheetLowestRungCells
%                                                   [nS x 1] .conj [nS x 3]
%                                                   (pass refuted notRun)
%                                                   .rungEntries [nR x 1]
%                                                   .rungTfDays [nR x 3]
%                                                   (min median max)
%                                                   .rungDvKms [nR x 2]
%                                                   (min max)
%
%% Revision History:
%  M. Casey                                                   (c) 11/21/2026
%  Copyright Coorbital Inc.
%% ------------------------ Begin Code Sequence ---------------------------

if nargin == 0
    here = fileparts(fileparts(mfilename('fullpath')));
    catalog = fullfile(here, 'results', 'costate_catalog_dro_elfo.mat');
    if ~isfile(catalog), catalog = fullfile(here, 'results', 'pilot_chain', 'costate_catalog_dro_elfo.mat'); end
    if ~isfile(catalog), fprintf('dro_elfo_describe: self-demo needs a catalog under results/\n');  D = [];  return, end
end
if ischar(catalog) || isstring(catalog)
    L = load(char(catalog));  fn = fieldnames(L);
    assert(isscalar(fn), 'dro_elfo_describe:catalog', 'expected one variable in %s', char(catalog));
    catalog = L.(fn{1});
end
tStar = catalog.constants.tStar_s;
lStar = catalog.constants.lStar_km;
m0kg  = catalog.thruster.m0_kg;
cnd   = catalog.thruster.c_nd;
d2day = tStar/86400;
rungs = catalog.rungs_N(:).';
nR    = numel(rungs);
S     = catalog.sheets;
nS    = numel(S);
name  = '';  if isfield(catalog, 'name'), name = catalog.name; end

%% Per sheet:
droTau  = zeros(1, nS);  raan = zeros(1, nS);
nEnt    = zeros(nS, 1);  nFull = zeros(nS, 1);  nLow = zeros(nS, 1);  conj = zeros(nS, 3);
tfAll   = cell(nR, 1);
for ks = 1:nS
    sh = S(ks);
    droTau(ks) = sh.sheet_key.depSize;
    raan(ks)   = sh.sheet_key.raanDeg;
    has        = logical(sh.has_solution);
    nEnt(ks)   = nnz(has);
    nFull(ks)  = nnz(all(has, 3));
    nLow(ks)   = nnz(has(:, :, end));
    if isfield(sh, 'conj_pass') && ~isempty(sh.conj_pass)
        cp = double(sh.conj_pass);
        conj(ks, :) = [nnz(has & cp == 1), nnz(has & cp == 0), nnz(has & cp ~= 1 & cp ~= 0)];
    else
        conj(ks, :) = [0 0 nEnt(ks)];
    end
    for kr = 1:nR
        hk = has(:, :, kr);
        ek = sh.entry_index(:, :, kr);
        tfAll{kr} = [tfAll{kr}; sh.z8(8, ek(hk)).'];
    end
end

%% Per rung over all sheets:
rungEnt = zeros(nR, 1);  rungTf = nan(nR, 3);  rungDv = nan(nR, 2);
for kr = 1:nR
    tf = tfAll{kr};
    rungEnt(kr) = numel(tf);
    if isempty(tf), continue, end
    Tnd = (rungs(kr)/m0kg)*tStar^2/(lStar*1000);
    dv  = cnd*log(1./(1 - Tnd*tf/cnd))*lStar/tStar;
    rungTf(kr, :) = [min(tf) median(tf) max(tf)]*d2day;
    rungDv(kr, :) = [min(dv) max(dv)];
end

%% Orbits:
sh1    = S(1);
oe     = sh1.arr_params.oe;
tauArr = sh1.tau_arr;
D = struct('name', name, 'nSheets', nS, 'nEntries', sum(nEnt), 'rungsN', rungs, ...
           'ispS', catalog.thruster.isp_s, 'm0kg', m0kg, 'nD', numel(sh1.sD_frac), 'nA', numel(sh1.sA_frac), ...
           'droTau', droTau, 'droPeriodDays', droTau*d2day, 'raanDeg', raan, 'elfoOe', oe, ...
           'elfoPeriodNd', tauArr, 'elfoPeriodDays', tauArr*d2day, 'sheetEntries', nEnt, ...
           'sheetFullLadders', nFull, 'sheetLowestRungCells', nLow, 'conj', conj, ...
           'rungEntries', rungEnt, 'rungTfDays', rungTf, 'rungDvKms', rungDv);

%% Print:
fprintf('\n================ DRO -> ELFO COSTATE CATALOG ================\n');
if ~isempty(name), fprintf('%s\n', name); end
fprintf('%d sheet(s), %d entries\n', nS, D.nEntries);
uT = unique(droTau);
fprintf('\nDRO period(s): ND [%s]  days [%s]\n', sprintf('%.4g ', uT), sprintf('%.3f ', uT*d2day));
fprintf('ELFO: a %.0f km, e %.3g, inc %.4g deg, argp %.4g deg; period %.4f ND = %.3f d\n', ...
        oe(1), oe(2), oe(3), oe(4), tauArr, tauArr*d2day);
fprintf('ELFO plane orientations raanDeg (one sheet each): [%s]\n', sprintf('%g ', raan));
fprintf('thrust rungs (N): [%s]  Isp %g s, m0 %g kg\n', sprintf('%g ', rungs), catalog.thruster.isp_s, m0kg);
fprintf('phase grid per sheet: %d DRO phases x %d ELFO phases (fractions of each period)\n', D.nD, D.nA);
fprintf('  DRO phase step %.4f d, ELFO phase step %.4f d\n', sh1.tau_dep*d2day/D.nD, tauArr*d2day/D.nA);

fprintf('\nsheet (DRO ND, raanDeg)   entries  full ladders  cells at %g N   conj pass/refuted/not run\n', rungs(end));
for ks = 1:nS
    fprintf('  %2d  (%6.4g, %5g)     %6d      %4d          %4d          %d / %d / %d\n', ks, droTau(ks), ...
            raan(ks), nEnt(ks), nFull(ks), nLow(ks), conj(ks, 1), conj(ks, 2), conj(ks, 3));
end
fprintf('\nper rung, all sheets:\n  thrust N  entries   t_f min / median / max [d]   Delta-V min / max [km/s]\n');
for kr = 1:nR
    fprintf('  %7g   %6d   %7.3f / %7.3f / %7.3f     %7.4f / %7.4f\n', rungs(kr), rungEnt(kr), ...
            rungTf(kr, 1), rungTf(kr, 2), rungTf(kr, 3), rungDv(kr, 1), rungDv(kr, 2));
end
fprintf('conjugate test, all sheets: %d pass, %d refuted, %d not run\n', sum(conj(:, 1)), sum(conj(:, 2)), sum(conj(:, 3)));
fprintf('=============================================================\n\n');
end
