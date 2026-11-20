function e = dro_elfo_pick(catalog, droTau, raanDeg, sD, sA, thrustN)
%% Purpose:
%
%   THE PICKER of the DRO -> ELFO costate catalog: one entry by (DRO size,
%   ELFO plane orientation, departure phase, arrival phase, thrust).
%
%   SHEET: exact match on the packaged sheets' additive sheet_key (written
%   by build_dro_elfo_catalog; spec section E): target 'elfo', depSize within 1e-9 of
%   droTau, raanDeg equal to the requested angle after both are wrapped to
%   [0, 360) (so 370 picks the 10 deg sheet), again within 1e-9. There is
%   NO nearest-sheet fallback: an entry of another plane orientation is a
%   different transfer, so a key that is not in the catalog errors
%   dro_elfo_pick:noSheet and lists the keys that are.
%
%   RUNG: thrustN must be one of the catalog's rungs (within 1e-9 N), else
%   dro_elfo_pick:noRung listing them.
%
%   CELL: the nearest phase cell (circular distance on the unit torus of
%   departure x arrival phase fractions, requested phases wrapped to
%   [0, 1)) AMONG the cells of that sheet holding an entry at that rung;
%   ties go to the first cell in (iD, iA) order. When the nearest cell of
%   the torus has NO entry at that rung and another cell is used, e.nearestEmpty
%   is true and the warning dro_elfo_pick:nearestEmpty names the requested
%   phases, the empty nearest cell and the cell used. The delivered cell's own
%   phases are returned (e.sD, e.sA) with the distance (e.dPhase): the z8 is
%   the converged costate OF THAT CELL, a guess for any other phase. No
%   entry at that rung on the sheet: dro_elfo_pick:noEntry.
%
%   VERDICTS: conj_pass as the catalog records it (1 pass, 0 refuted, -1
%   not run or undecided -- conj_catalog_pass's convention; -1 when the
%   sheet carries no verdicts) and the hypothesis-gate values (NaN when not
%   gated). A refuted entry is returned, with its verdict: it is an
%   extremal, not a minimum.
%
%% Inputs:
%
%  catalog                  struct | char           The catalog struct, or
%                                                   the path of its .mat
%                                                   (one variable)
%
%  droTau                   double                  DRO size (period key,
%                                                   ND)
%
%  raanDeg                  double                  ELFO plane orientation,
%                                                   deg (any real; wrapped)
%
%  sD, sA                   double                  Departure and arrival
%                                                   phase fractions (any
%                                                   real; wrapped to [0,1))
%
%  thrustN                  double                  Thrust, N (a catalog
%                                                   rung)
%
%% Outputs:
%
%  e                        struct                  .z8 [8 x 1] [lambda0(7);
%                                                   tf] .tf_nd .sD .sA (the
%                                                   cell's phases) .iD .iA
%                                                   .kr .dPhase .nearestEmpty
%                                                   (the nearest cell had no
%                                                   entry) .sheet
%                                                   (index) .sheet_key
%                                                   .thrustN .conj_pass
%                                                   .gate_min_lamv
%                                                   .gate_min_qmt .gate_dimS
%                                                   .dep_family .dep_params
%                                                   .arr_family .arr_params
%
%% Revision History:
%  M. Casey                                                   (c) 11/19/2026
%  M. Casey  fix round 1: target in the match, nearestEmpty   (c) 11/19/2026
%  Copyright Coorbital Inc.
%% ------------------------ Begin Code Sequence ---------------------------

if nargin == 0
    here = fileparts(fileparts(mfilename('fullpath')));
    catalog = fullfile(here, 'results', 'costate_catalog_dro_elfo.mat');
    if ~isfile(catalog)
        fprintf('dro_elfo_pick: self-demo needs %s (build_dro_elfo_catalog)\n', catalog);  e = [];  return
    end
    L = load(catalog);  fn = fieldnames(L);  sh1 = L.(fn{1}).sheets(1);
    droTau = sh1.sheet_key.depSize;  raanDeg = sh1.sheet_key.raanDeg;
    sD = 0.1;  sA = 0.2;  thrustN = L.(fn{1}).rungs_N(1);
end
if ischar(catalog) || isstring(catalog)
    L = load(char(catalog));  fn = fieldnames(L);
    assert(isscalar(fn), 'dro_elfo_pick:catalog', 'expected one variable in %s', char(catalog));
    catalog = L.(fn{1});
end
tol = 1e-9;
assert(isfield(catalog, 'sheets') && ~isempty(catalog.sheets) && isfield(catalog.sheets, 'sheet_key'), ...
       'dro_elfo_pick:catalog', 'the catalog has no sheets carrying sheet_key (build_dro_elfo_catalog writes it)');

%% Sheet: exact key match, angle wrapped
r = mod(raanDeg, 360);
nS = numel(catalog.sheets);
hit = false(1, nS);
keys = cell(1, nS);
for ks = 1:nS
    k = catalog.sheets(ks).sheet_key;
    keys{ks} = sprintf('(target %s, depSize %.12g, raanDeg %.12g)', k.target, k.depSize, k.raanDeg);
    dAng = abs(mod(k.raanDeg - r + 180, 360) - 180);
    hit(ks) = strcmp(k.target, 'elfo') && abs(k.depSize - droTau) < tol && dAng < tol;
end
if ~any(hit)
    error('dro_elfo_pick:noSheet', ['no ELFO sheet with depSize %.12g and raanDeg %.12g (wrapped from %.12g) -- ', ...
          'available keys: %s'], droTau, r, raanDeg, strjoin(keys, ', '));
end
if nnz(hit) > 1
    error('dro_elfo_pick:dupSheet', 'the catalog holds %d ELFO sheets keyed (depSize %.12g, raanDeg %.12g)', ...
          nnz(hit), droTau, r);
end
ks = find(hit);
sh = catalog.sheets(ks);

%% Rung
kr = find(abs(catalog.rungs_N - thrustN) < tol, 1);
if isempty(kr)
    error('dro_elfo_pick:noRung', 'thrust %g N is not a catalog rung -- rungs: %s N', ...
          thrustN, mat2str(catalog.rungs_N));
end

%% Nearest phase cell holding an entry at that rung
has = sh.has_solution(:, :, kr);
if ~any(has(:))
    error('dro_elfo_pick:noEntry', 'sheet %s holds no entry at %g N', keys{ks}, thrustN);
end
dD = circ(sh.sD_frac(:) - mod(sD, 1));               % [nD x 1]
dA = circ(sh.sA_frac(:).' - mod(sA, 1));             % [1 x nA]
dist = sqrt(dD.^2 + dA.^2);
[~, k0] = min(dist(:));                              % nearest cell of the torus
dist(~has) = inf;
[dMin, kc] = min(dist(:));                           % nearest cell holding an entry
[iD, iA] = ind2sub(size(has), kc);
nearestEmpty = kc ~= k0;
if nearestEmpty
    [iD0, iA0] = ind2sub(size(has), k0);
    warning('dro_elfo_pick:nearestEmpty', ['sheet %s at %g N: the nearest cell (%d,%d) to the ', ...
            'requested phases (%.4f, %.4f) has no entry; cell (%d,%d) at (%.4f, %.4f) is used'], ...
            keys{ks}, thrustN, iD0, iA0, mod(sD, 1), mod(sA, 1), iD, iA, sh.sD_frac(iD), sh.sA_frac(iA));
end
col = sh.entry_index(iD, iA, kr);
z8  = sh.z8(:, col);

e = struct('z8', z8, 'tf_nd', z8(8), 'sD', sh.sD_frac(iD), 'sA', sh.sA_frac(iA), ...
           'iD', iD, 'iA', iA, 'kr', kr, 'dPhase', dMin, 'nearestEmpty', nearestEmpty, 'sheet', ks, 'sheet_key', sh.sheet_key, ...
           'thrustN', catalog.rungs_N(kr), 'conj_pass', -1, 'gate_min_lamv', NaN, ...
           'gate_min_qmt', NaN, 'gate_dimS', NaN, 'dep_family', sh.dep_family, ...
           'dep_params', sh.dep_params, 'arr_family', sh.arr_family, 'arr_params', sh.arr_params);
if isfield(sh, 'conj_pass') && ~isempty(sh.conj_pass), e.conj_pass = double(sh.conj_pass(iD, iA, kr)); end
for g = {'gate_min_lamv', 'gate_min_qmt', 'gate_dimS'}
    if isfield(sh, g{1}) && ~isempty(sh.(g{1})), e.(g{1}) = double(sh.(g{1})(iD, iA, kr)); end
end
if nargin == 0
    fprintf('dro_elfo_pick: %s cell (%d,%d) sD %.4f sA %.4f at %g N: tf %.6f ND, conj_pass %d\n', ...
            keys{ks}, iD, iA, e.sD, e.sA, e.thrustN, e.tf_nd, e.conj_pass);
end
end

% ------------------------------------------------------------------------
function d = circ(x)
% CIRC  Circular distance on the unit circle.  INPUTS: x [any size]
% phase differences.  OUTPUTS: d [size of x] distances in [0, 0.5].
d = abs(mod(x + 0.5, 1) - 0.5);
end
