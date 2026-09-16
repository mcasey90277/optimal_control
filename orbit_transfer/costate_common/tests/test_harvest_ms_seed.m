function ok = test_harvest_ms_seed()
%% Purpose:
%
%   Tests harvest_ms_seed -- direct collocation solution in, multiple-
%   shooting seed out -- on REAL Hermite-Simpson defect multipliers (the
%   golden harvest cell, committed in golden_cells_data.mat, so the test
%   runs from a fresh clone). The covector rules themselves are tested in
%   oclib (test_duals_to_costates); this test pins what harvest_ms_seed
%   owns and how it hands the rules their inputs:
%     1. SHAPE: tf copied, tGrid = linspace(0, tf, K+1), Y [14 x K+1].
%     2. STATE ROWS: pchip of the node states onto tGrid, reproduced
%        bitwise; the ends land on the first and last node.
%     3. COSTATE ROWS: the duals_to_costates output at the MIDPOINT
%        stations, pchip-extrapolated onto tGrid, reproduced bitwise --
%        i.e. rows 1:7 of the multipliers go in, with the Hermite-Simpson
%        scheme and the first three midpoint-control rows as uDir.
%     4. DIAGNOSTICS: a unanimous vote and lambda_t = +1 on the real cell.
%     5. SIGN: flipping every multiplier gives the SAME seed (the vote puts
%        it back) and sign -1 in the diagnostics.
%     6. UNIFORM MESH: tNodes = [] rebuilds linspace(0, tf, N+1) -- the
%        same seed as passing that mesh explicitly.
%     7. NO lambda_t ROW: a 7-row multiplier array skips the check (NaN)
%        and still returns the same seed.
%
%% Inputs:
%
%  none
%
%% Outputs:
%
%  ok                       logical                 All checks passed
%
%% Revision History:
%  M. Casey                                                   (c) 09/16/2026
%  Copyright Coorbital Inc.
%% ------------------------ Begin Code Sequence ---------------------------

ok = true;
here = fileparts(fileparts(mfilename('fullpath')));
addpath(here);
if isempty(which('oc.duals_to_costates'))
    addpath(fullfile(fileparts(fileparts(here)), 'oclib'));
end
G = load(fullfile(here, 'golden_cells_data.mat'), 'harvestCell');
h = G.harvestCell;
o = struct('X', h.X, 'lamDef', h.lamDef, 'Um', h.Um, 'tNodes', h.tNodes, 'tf', h.tf);
K = 24;

%% 1. Shape:
[seed, dg] = harvest_ms_seed(o, K);
ok = chk(ok, seed.tf == h.tf && isequal(seed.tGrid, linspace(0, h.tf, K+1)), ...
         'tf copied and tGrid = linspace(0, tf, K+1)');
ok = chk(ok, isequal(size(seed.Y), [14 K+1]), sprintf('Y is [14 x K+1]: %s', mat2str(size(seed.Y))));

%% 2. State rows:
Xg = interp1(h.tNodes, h.X', seed.tGrid, 'pchip')';
ok = chk(ok, isequal(seed.Y(1:7,:), Xg(1:7,:)), 'state rows are pchip of the node states, bitwise');
ok = chk(ok, max(abs(seed.Y(1:7,1) - h.X(1:7,1))) < 1e-14 ...
             && max(abs(seed.Y(1:7,end) - h.X(1:7,end))) < 1e-12, ...
         sprintf('ends land on the first and last node (%.1e, %.1e)', ...
                 max(abs(seed.Y(1:7,1) - h.X(1:7,1))), max(abs(seed.Y(1:7,end) - h.X(1:7,end)))));

%% 3. Costate rows:
spec = struct('scheme', 'hermite-simpson', 'mu', h.lamDef(1:7,:), 'tNodes', h.tNodes, ...
              'uDir', h.Um(1:3, 1:size(h.lamDef,2)), 'lamTf', h.lamDef(8,:));
[lam, tMid] = oc.duals_to_costates(spec);
ok = chk(ok, max(abs(tMid - (h.tNodes(1:end-1) + diff(h.tNodes)/2))) < 1e-15, ...
         'the costate stations are the interval MIDPOINTS');
Lg = interp1(tMid, lam', seed.tGrid, 'pchip', 'extrap')';
ok = chk(ok, isequal(seed.Y(8:14,:), Lg), ...
         'costate rows are the mapped costates, pchip-extrapolated onto tGrid, bitwise');

%% 4. Diagnostics on the real cell:
ok = chk(ok, dg.voteMargin == 1 && dg.lamTOK, ...
         sprintf('unanimous vote (%.2f) and lambda_t = %.6f', dg.voteMargin, dg.lamT));

%% 5. Sign:
oF = o;  oF.lamDef = -o.lamDef;
[seedF, dgF] = harvest_ms_seed(oF, K);
ok = chk(ok, isequal(seedF.Y, seed.Y) && dgF.sign == -dg.sign, ...
         sprintf('flipped multipliers give the same seed (sign %+d -> %+d)', dg.sign, dgF.sign));

%% 6. Uniform mesh:
oU = o;  oU.tNodes = linspace(0, o.tf, size(o.X, 2));
oE = o;  oE.tNodes = [];
sU = harvest_ms_seed(oU, K);  sE = harvest_ms_seed(oE, K);
ok = chk(ok, isequal(sE.Y, sU.Y), 'tNodes = [] rebuilds the uniform mesh (same seed as passing it)');

%% 7. No lambda_t row:
o7 = o;  o7.lamDef = o.lamDef(1:7,:);
[s7, d7] = harvest_ms_seed(o7, K);
ok = chk(ok, isnan(d7.lamT) && ~d7.lamTOK && isequal(s7.Y, seed.Y), ...
         'a 7-row multiplier array skips the lambda_t check and gives the same seed');

if ok, fprintf('TEST_HARVEST_MS_SEED: ALL PASS\n');
else,  fprintf('TEST_HARVEST_MS_SEED: FAILURE (see lines above)\n');
end
end

% ------------------------------------------------------------------------
function ok = chk(ok, cond, label)
%% Purpose:
%
%   Accumulate one labelled pass/fail.
%
if cond, fprintf('  PASS  %s\n', label);
else,    fprintf('  FAIL  %s\n', label);  ok = false;
end
end
