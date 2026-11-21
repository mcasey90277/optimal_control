function T = dro_elfo_transfers_at(catalog, droTau, raanDeg, sD, sA)
%% Purpose:
%
%   EVERY THRUST RUNG recorded at one phase cell of one sheet of the
%   DRO -> ELFO catalog, ordered by thrust from high to low: the thrust
%   ladder of that (departure phase, arrival phase) pair.
%
%   SHEET: exact match on sheet_key by dro_elfo_pick's rule (target 'elfo',
%   depSize within 1e-9, raanDeg wrapped to [0, 360) within 1e-9), through
%   dro_elfo_find_sheet; NO nearest-sheet fallback: an unknown key errors
%   dro_elfo_transfers_at:noSheet listing the keys that are in the catalog.
%
%   CELL: the nearest phase cell of the torus (circular distance on the
%   unit torus, requested phases wrapped to [0, 1); ties to the first cell
%   in (iD, iA) order), WHETHER OR NOT it holds an entry: this function
%   reports what one cell holds, so it never moves to another cell (unlike
%   dro_elfo_pick, which looks for the nearest cell with an entry at one
%   rung). The cell's own phases are returned in each row.
%
%   Rows are built by dro_elfo_entry (the same fields as a
%   dro_elfo_extremes row). A cell that holds nothing returns an EMPTY
%   struct array with the same fields.
%
%% Inputs:
%
%  catalog                  struct | char           The catalog struct, or
%                                                   the path of its .mat
%                                                   (one variable)
%
%  droTau                   double                  DRO size (period key,
%                                                   ND) [1 x 1]
%
%  raanDeg                  double                  ELFO plane orientation,
%                                                   deg (any real; wrapped)
%                                                   [1 x 1]
%
%  sD, sA                   double                  Departure and arrival
%                                                   phase fractions (any
%                                                   real; wrapped to [0,1))
%                                                   [1 x 1]
%
%% Outputs:
%
%  T                        struct                  [nRungs x 1] rows, high
%                                                   to low thrust: .thrustN
%                                                   .sD .sA .iD .iA .kr
%                                                   .tfDays .tf_nd .dvKms
%                                                   .propellantKg
%                                                   .finalMassKg .conj_pass
%                                                   .gate_min_lamv
%                                                   .gate_min_qmt .gate_dimS
%                                                   .z8 [8 x 1] .raanDeg
%                                                   .droTau .sheet
%                                                   .sheet_key .dep_family
%                                                   .dep_params .arr_family
%                                                   .arr_params .label;
%                                                   [0 x 1] when the cell
%                                                   holds nothing
%
%% Revision History:
%  M. Casey                                                   (c) 11/21/2026
%  Copyright Coorbital Inc.
%% ------------------------ Begin Code Sequence ---------------------------

if nargin == 0
    here = fileparts(fileparts(mfilename('fullpath')));
    catalog = fullfile(here, 'results', 'costate_catalog_dro_elfo.mat');
    if ~isfile(catalog), catalog = fullfile(here, 'results', 'pilot_chain', 'costate_catalog_dro_elfo.mat'); end
    if ~isfile(catalog), fprintf('dro_elfo_transfers_at: self-demo needs a catalog under results/\n');  T = [];  return, end
    L = load(catalog);  fn = fieldnames(L);  catalog = L.(fn{1});
    k = catalog.sheets(1).sheet_key;
    T = dro_elfo_transfers_at(catalog, k.depSize, k.raanDeg, 0.25, 0.5);
    for kk = 1:numel(T), fprintf('  %s\n', T(kk).label); end
    return
end
if ischar(catalog) || isstring(catalog)
    L = load(char(catalog));  fn = fieldnames(L);
    assert(isscalar(fn), 'dro_elfo_transfers_at:catalog', 'expected one variable in %s', char(catalog));
    catalog = L.(fn{1});
end
ks = dro_elfo_find_sheet(catalog, droTau, raanDeg, 'dro_elfo_transfers_at');
sh = catalog.sheets(ks);

%% The nearest cell of the torus:
dD = circ(sh.sD_frac(:) - mod(sD, 1));               % [nD x 1]
dA = circ(sh.sA_frac(:).' - mod(sA, 1));             % [1 x nA]
dist = sqrt(dD.^2 + dA.^2);
[~, kc] = min(dist(:));
[iD, iA] = ind2sub(size(dist), kc);

%% Every rung with an entry there, high thrust first:
krs = find(squeeze(sh.has_solution(iD, iA, :))).';
[~, ord] = sort(catalog.rungs_N(krs), 'descend');
krs = krs(ord);
if isempty(krs)
    fields = {'thrustN', 'sD', 'sA', 'iD', 'iA', 'kr', 'tfDays', 'tf_nd', 'dvKms', 'propellantKg', ...
              'finalMassKg', 'conj_pass', 'gate_min_lamv', 'gate_min_qmt', 'gate_dimS', 'z8', ...
              'raanDeg', 'droTau', 'sheet', 'sheet_key', 'dep_family', 'dep_params', 'arr_family', ...
              'arr_params', 'label'};
    T = cell2struct(cell(numel(fields), 0), fields, 1);
    return
end
for kk = numel(krs):-1:1
    T(kk, 1) = dro_elfo_entry(catalog, ks, iD, iA, krs(kk));
end
end

% ------------------------------------------------------------------------
function d = circ(x)
% CIRC  Circular distance on the unit circle.  INPUTS: x [any size]
% phase differences.  OUTPUTS: d [size of x] distances in [0, 0.5].
d = abs(mod(x + 0.5, 1) - 0.5);
end
