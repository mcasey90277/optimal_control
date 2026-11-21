function ks = dro_elfo_find_sheet(catalog, droTau, raanDeg, caller)
%% Purpose:
%
%   THE SHEET MATCH of the DRO -> ELFO library functions, by the same rule
%   as dro_elfo_pick: target 'elfo', sheet_key.depSize within 1e-9 of
%   droTau, sheet_key.raanDeg equal to the requested angle after both are
%   wrapped to [0, 360), within 1e-9. No nearest-sheet fallback: a key that
%   is not in the catalog errors <caller>:noSheet and lists the keys that
%   are; two sheets with one key error <caller>:dupSheet; a catalog without
%   sheet_key errors <caller>:catalog. The caller's name prefixes every
%   identifier, so each library function refuses under its own name.
%
%% Inputs:
%
%  catalog                  struct                  The catalog struct
%                                                   (.sheets(k).sheet_key)
%
%  droTau                   double                  DRO size (period key,
%                                                   ND) [1 x 1]
%
%  raanDeg                  double                  ELFO plane orientation,
%                                                   deg (any real; wrapped)
%                                                   [1 x 1]
%
%  caller                   char                    Identifier prefix, e.g.
%                                                   'dro_elfo_transfers_at'
%                                                   [1 x n]
%
%% Outputs:
%
%  ks                       double                  Index of the matching
%                                                   sheet [1 x 1]
%
%% Revision History:
%  M. Casey                                                   (c) 11/21/2026
%  Copyright Coorbital Inc.
%% ------------------------ Begin Code Sequence ---------------------------

if nargin == 0
    C = struct('sheets', struct('sheet_key', {struct('target', 'elfo', 'depSize', 1, 'raanDeg', 0), ...
                                              struct('target', 'elfo', 'depSize', 1, 'raanDeg', 30)}));
    ks = dro_elfo_find_sheet(C, 1, 390, 'dro_elfo_find_sheet');
    fprintf('dro_elfo_find_sheet: depSize 1, raanDeg 390 -> sheet %d (raanDeg %g)\n', ks, C.sheets(ks).sheet_key.raanDeg);
    return
end
if nargin < 4 || isempty(caller), caller = 'dro_elfo_find_sheet'; end
tol = 1e-9;
if ~(isfield(catalog, 'sheets') && ~isempty(catalog.sheets) && isfield(catalog.sheets, 'sheet_key'))
    error([caller ':catalog'], 'the catalog has no sheets carrying sheet_key (build_dro_elfo_catalog writes it)');
end
r    = mod(raanDeg, 360);
nS   = numel(catalog.sheets);
hit  = false(1, nS);
keys = cell(1, nS);
for kk = 1:nS
    k = catalog.sheets(kk).sheet_key;
    keys{kk} = sprintf('(target %s, depSize %.12g, raanDeg %.12g)', k.target, k.depSize, k.raanDeg);
    dAng = abs(mod(k.raanDeg - r + 180, 360) - 180);
    hit(kk) = strcmp(k.target, 'elfo') && abs(k.depSize - droTau) < tol && dAng < tol;
end
if ~any(hit)
    error([caller ':noSheet'], ['no ELFO sheet with depSize %.12g and raanDeg %.12g (wrapped from %.12g) -- ', ...
          'available keys: %s'], droTau, r, raanDeg, strjoin(keys, ', '));
end
if nnz(hit) > 1
    error([caller ':dupSheet'], 'the catalog holds %d ELFO sheets keyed (depSize %.12g, raanDeg %.12g)', ...
          nnz(hit), droTau, r);
end
ks = find(hit);
end
