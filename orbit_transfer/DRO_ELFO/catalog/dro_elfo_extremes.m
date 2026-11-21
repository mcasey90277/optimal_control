function [eMin, eMax] = dro_elfo_extremes(catalog, filt, metric)
%% Purpose:
%
%   THE EXTREME ENTRIES of a DRO -> ELFO catalog -- by flight time
%   (default) or by Delta-V -- under a filter, each returned as a full row
%   (dro_elfo_entry: every field dro_elfo_fly needs, the derived units,
%   the verdicts and a one-line label) and both printed.
%
%   At one thrust the minimum-time and minimum-Delta-V entries are the same
%   transfer (all-burn: Delta-V grows with t_f); they separate only across
%   rungs, so hold thrustN fixed to compare phasing, and raanDeg fixed to
%   compare one plane orientation.
%
%  ASSUMPTIONS / NOTES:
%
% - Filter fields (any subset; an unknown field errors
%   dro_elfo_extremes:filter): .raanDeg (scalar or vector, deg, matched
%   after wrapping to [0, 360) within 1e-9), .droTau (scalar, ND, within
%   1e-9), .thrustN (scalar or vector, N, rungs within 1e-9),
%   .conjPassOnly [false] (only entries whose conjugate test passed,
%   conj_pass == 1).
% - Ties go to the first entry in (sheet, rung, iA, iD) storage order.
% - Nothing matches: error dro_elfo_extremes:noMatch.
%
%% Inputs:
%
%  catalog                  struct | char           The catalog struct, or
%                                                   the path of its .mat
%
%  filt                     struct (optional)       [1 x 1] the filter
%                                                   above; [] = everything
%
%  metric                   char (optional)         'time' (default) or
%                                                   'deltaV' [1 x n]
%
%% Outputs:
%
%  eMin, eMax               struct                  [1 x 1] each, the
%                                                   minimum and maximum
%                                                   rows (dro_elfo_entry
%                                                   fields: .z8 .sD .sA
%                                                   .thrustN .sheet_key
%                                                   .tfDays .dvKms
%                                                   .propellantKg
%                                                   .conj_pass .raanDeg
%                                                   .droTau .label ...)
%
%% Revision History:
%  M. Casey                                                   (c) 11/21/2026
%  Copyright Coorbital Inc.
%% ------------------------ Begin Code Sequence ---------------------------

if nargin == 0
    here = fileparts(fileparts(mfilename('fullpath')));
    catalog = fullfile(here, 'results', 'costate_catalog_dro_elfo.mat');
    if ~isfile(catalog), catalog = fullfile(here, 'results', 'pilot_chain', 'costate_catalog_dro_elfo.mat'); end
    if ~isfile(catalog), fprintf('dro_elfo_extremes: self-demo needs a catalog under results/\n');  eMin = [];  eMax = [];  return, end
    [eMin, eMax] = dro_elfo_extremes(catalog, struct('thrustN', 1), 'time');
    return
end
if nargin < 2 || isempty(filt), filt = struct(); end
if nargin < 3 || isempty(metric), metric = 'time'; end
if ischar(catalog) || isstring(catalog)
    L = load(char(catalog));  fn = fieldnames(L);
    assert(isscalar(fn), 'dro_elfo_extremes:catalog', 'expected one variable in %s', char(catalog));
    catalog = L.(fn{1});
end
metric = char(metric);
if ~any(strcmp(metric, {'time', 'deltaV'}))
    error('dro_elfo_extremes:metric', 'metric must be ''time'' or ''deltaV'', got ''%s''', metric);
end
known = {'raanDeg', 'droTau', 'thrustN', 'conjPassOnly'};
bad = setdiff(fieldnames(filt), known);
if ~isempty(bad)
    error('dro_elfo_extremes:filter', 'unknown filter field(s) %s -- known: %s', strjoin(bad, ', '), strjoin(known, ', '));
end
tol      = 1e-9;
tStar    = catalog.constants.tStar_s;
lStar    = catalog.constants.lStar_km;
m0kg     = catalog.thruster.m0_kg;
cnd      = catalog.thruster.c_nd;
passOnly = isfield(filt, 'conjPassOnly') && ~isempty(filt.conjPassOnly) && logical(filt.conjPassOnly);

%% Walk every stored entry under the filter:
best = [];  worst = [];  vBest = inf;  vWorst = -inf;
for ks = 1:numel(catalog.sheets)
    sh = catalog.sheets(ks);
    k  = sh.sheet_key;
    if isfield(filt, 'droTau') && ~isempty(filt.droTau) && abs(k.depSize - filt.droTau) >= tol, continue, end
    if isfield(filt, 'raanDeg') && ~isempty(filt.raanDeg)
        dAng = abs(mod(k.raanDeg - mod(filt.raanDeg(:), 360) + 180, 360) - 180);
        if ~any(dAng < tol), continue, end
    end
    for kr = 1:numel(catalog.rungs_N)
        TN = catalog.rungs_N(kr);
        if isfield(filt, 'thrustN') && ~isempty(filt.thrustN) && ~any(abs(filt.thrustN(:) - TN) < tol), continue, end
        ok = logical(sh.has_solution(:, :, kr));
        if passOnly
            if ~isfield(sh, 'conj_pass') || isempty(sh.conj_pass), continue, end
            ok = ok & double(sh.conj_pass(:, :, kr)) == 1;
        end
        cells = find(ok);
        if isempty(cells), continue, end
        ek = sh.entry_index(:, :, kr);
        tf = sh.z8(8, ek(cells)).';
        if strcmp(metric, 'time')
            v = tf;
        else
            Tnd = (TN/m0kg)*tStar^2/(lStar*1000);
            v   = cnd*log(1./(1 - Tnd*tf/cnd))*lStar/tStar;
        end
        [vLo, kLo] = min(v);
        [vHi, kHi] = max(v);
        if vLo < vBest,  vBest  = vLo;  best  = [ks, kr, cells(kLo)]; end
        if vHi > vWorst, vWorst = vHi;  worst = [ks, kr, cells(kHi)]; end
    end
end
if isempty(best)
    error('dro_elfo_extremes:noMatch', 'no catalog entry matches the filter%s', filtstr(filt));
end
nD = numel(catalog.sheets(best(1)).sD_frac);
[iD, iA] = ind2sub([nD, numel(catalog.sheets(best(1)).sA_frac)], best(3));
eMin = dro_elfo_entry(catalog, best(1), iD, iA, best(2));
nD = numel(catalog.sheets(worst(1)).sD_frac);
[iD, iA] = ind2sub([nD, numel(catalog.sheets(worst(1)).sA_frac)], worst(3));
eMax = dro_elfo_entry(catalog, worst(1), iD, iA, worst(2));

%% Report:
if strcmp(metric, 'time'), mName = 'flight time'; else, mName = 'Delta-V'; end
fprintf('\n%s extremes%s:\n', mName, filtstr(filt));
fprintf('  MINIMUM  %s\n', eMin.label);
fprintf('  MAXIMUM  %s\n', eMax.label);
end

% ------------------------------------------------------------------------
function s = filtstr(filt)
% FILTSTR  The filter as text.  INPUTS: filt [struct].  OUTPUTS: s [char],
% '' for an empty filter.
fn = fieldnames(filt);
parts = cell(1, numel(fn));
for kk = 1:numel(fn)
    v = filt.(fn{kk});
    if ~isempty(v), parts{kk} = sprintf('%s = %s', fn{kk}, mat2str(double(v), 6)); end
end
parts = parts(~cellfun(@isempty, parts));
if isempty(parts), s = ''; else, s = [' (' strjoin(parts, ', ') ')']; end
end
