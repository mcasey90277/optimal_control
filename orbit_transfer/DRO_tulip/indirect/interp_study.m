%% INTERP_STUDY  One transfer that is NOT in the library, guessed from the ones that are
%
%   Edit the parameter blocks, press Run. The costate library holds certified
%   minimum-time transfers on a GRID of (departure phase, arrival phase). This
%   script asks the question the library exists for: given a phase pair that
%   is not on the grid, do the entries around it blend into a costate guess
%   good enough to solve the transfer -- and is what comes out the same kind
%   of solution as its neighbours, and a certified one?
%
%   The steps are done here rather than behind a front door, so the machinery
%   is visible. Two things stay in the library because a second copy would be
%   an unverified one: the rule that picks which corners may be blended
%   (costate_common/phase_catalog_interp) and the fenced polish
%   (costate_common/polish_costate_guess). Section 3 works the blend rule out
%   by hand for THIS cell and asserts it against the library's answer.
%
%     0  tolerances and switches          (one source, set before anything runs)
%     1  the library                      (what problem it solves, its grid)
%     2  the query                        (the phase pair, its two endpoints)
%     3  the cell                         (four corners: which may be blended?)
%     4  the guess                        (the blend, and what the corners'
%                                          own costates predict for t_f)
%     5  fly the guess                    (how far from the target, unsolved?)
%     6  polish                           (the guess becomes a root, or not)
%     7  the same branch?                 (did the solve stay with its corners)
%     8  certify                          (the library's own full gate stack)
%     9  baseline                         (the same solve from the NEAREST
%                                          entry alone -- what blending bought)
%    10  verdict and figures
%
%   Diagnostic IDs are stable, and the letters mean:
%     B  BLEND    -- B1-B2, about the guess before any solve: were there
%                   corners that may be blended, and does the guess fly
%                   anywhere near the target.
%     P  POLISH   -- P1-P2: the solve converged; it stayed on the branch of
%                   the corners it was blended from.
%     C  CERTIFY  -- C1: certify_root, the gate stack every library entry
%                   passed (Pontryagin's necessary conditions, the second
%                   solver, the sufficiency hypotheses). transfer_study takes
%                   that stack apart one gate at a time; here it is one line.
%     V  VALIDITY -- V1: the corners and weights worked out by hand in
%                   section 3 are the ones the library function used.
%     Q  QUALITY  -- Q1-Q2, not gates. Q1: the t_f the corners' costates
%                   predict to first order against the blended t_f. Q2: the
%                   guess against the root it led to. They are what a finer
%                   grid would improve, so they are what to watch.
%
%   WHAT A FAILURE MEANS. B1 failing is a finding about the LIBRARY at this
%   spot (a hole, a family boundary, a fold): the script says which and
%   stops, because a mean of two different branches is a guess for neither.
%   P1 failing is a finding about the RESOLUTION. P2 or C1 failing means the
%   solve left its neighbours or found an uncertified root: do not use it.
%
%  M. Casey                                                   (c) 09/20/2026
%  Copyright Coorbital Inc.

%% paths
clear; clc
here = fileparts(mfilename('fullpath'));
addpath(here, fullfile(fileparts(fileparts(here)), 'costate_common'));

%% ========================================================================
%  0. TOLERANCES AND SWITCHES -- one source, used by every printed line AND
%     the verdict.
%% ========================================================================
tol = struct( ...
    'closure',   1e-7, ...   % periodicity of each orbit's endpoint interpolant
    'seamDeriv', 1e-6, ...   %   and of its derivative across s = 0
    'edgeMinAt24', 300, ...  % B1 two NEIGHBOURING corners are blended only if
    ...                      %   their flight times differ by no more than their
    ...                      %   own sensitivities explain, to this many minutes.
    ...                      %   On one branch the residual is third order in
    ...                      %   the grid step (about 100 min at 1/24, 2 min at
    ...                      %   1/96); across a change of branch it is the jump
    ...                      %   in t_f, of order 1000 min at ANY step. This is
    ...                      %   the cap at a step of 1/24; section 3 scales it
    ...                      %   as step^1.5, half way (in the exponent) between
    ...                      %   what it must admit and what it must refuse:
    ...                      %   38 min at 1/96, where one branch measures 2-18
    ...                      %   and a hand-over 115 (FINDINGS 86).
    'maxJump',   0.5,  ...   % B1 ... and their costates differ by less than this
    ...                      %   (relative). Blunter: a 1/24 step inside one
    ...                      %   branch measures 0.3-0.4, a change of branch
    ...                      %   0.7-0.9
    'guessKm',   20000, ...  % B2 the unsolved guess must fly at least this
    ...                      %   close (a SCREEN: beyond it the blend is not in
    ...                      %   the target's neighbourhood at all)
    'R',         1e-8, ...   % P1 polish residual, inf-norm
    'branch',    0.5,  ...   % P2 the root may differ from the guess by at most
    ...                      %   this (relative, costates) ...
    'branchTfDays', 1.0, ... %   ... and its t_f may leave the corners' range by
    ...                      %   at most this many days
    'agree',     1e-12);     % V1 hand-worked weights vs the library function

tfFromSensitivities = true;  % section 4: take the guess's t_f from the corners' first-order
                             %   predictions rather than from the linear blend of their t_f.
                             %   The costates carry dT/ds, so this uses information the plain
                             %   blend throws away; where t_f bends sharply it is closer by hours.
doCertify  = true;           % section 8 (about a minute)
doBaseline = true;           % section 9 (one more fenced solve)
capSec     = 600;            % hard wall-clock cap on each fenced solve
allowUnfenced = false;       % true: run the solves in THIS session if no worker
                             %   pool can be had. A poor guess can then hang the
                             %   session for hours inside one integration; that
                             %   is how the shared session was lost on 09-19.

%% ========================================================================
%  1. THE LIBRARY -- which problem it solves, and the grid it solves it on.
%     Everything physical below is read from the catalog file itself, so the
%     guess, the endpoints and the engine cannot disagree with the entries.
%% ========================================================================
%     TWO LIBRARIES are on disk, and the contrast between them is the point:
%       spine96   one departure phase (s_D = 0), 96 arrival phases: an arrival
%                 sheet wrapped as a one-row catalog (arrival_sheet_as_catalog).
%                 Fine enough that blending WORKS along arrival phase.
%       record    the 24 x 24 library of record. At 1/24 the costates change by
%                 30-40% per step INSIDE one branch, and no blend of them
%                 converges (FINDINGS 89): the reason the grid is being refined.
library = 'spine96';                          % 'spine96' | 'record' | a catalog file's path
switch library
    case 'spine96', catMat = fullfile(here, 'results', 'sheet96_resolution_test', 'spine96_catalog.mat');
    case 'record',  catMat = fullfile(here, 'results', 'library_70mN_24x24_final', 'costate_catalog_dro_tulip_70mN.mat');
    otherwise,      catMat = library;
end
assert(isfile(catMat), 'interp_study:noLibrary', 'no catalog at %s', catMat);

L = load(catMat);  fn = fieldnames(L);  cat_ = L.(fn{1});
assert(isscalar(cat_.sheets), 'interp_study reads a single-sheet phase catalog');
sheet = cat_.sheets(1);
muStar = cat_.constants.muStar;  lStar = cat_.constants.lStar_km;  tStar = cat_.constants.tStar_s;
day = @(tnd) tnd*tStar/86400;
thrustN = cat_.rungs_N(1);  ispS = cat_.thruster.isp_s;  m0kg = cat_.thruster.m0_kg;
ndp = nd_propulsion(thrustN, ispS, m0kg, lStar, tStar);      % THE shared propulsion conversion
Tnd = ndp.Tnd;  cnd = ndp.cnd;

has = logical(sheet.has_solution(:, :, 1));
[gridD, orderD] = sort(mod(sheet.sD_frac(:).', 1));          % the grid lines, in order round the circle
[gridA, orderA] = sort(mod(sheet.sA_frac(:).', 1));
famLabel = @(code) 'unlabelled';
if isfield(cat_, 'families') && isfield(cat_.families, 'labels')
    famLabel = @(code) tern(code >= 1 && code <= numel(cat_.families.labels), cat_.families.labels{max(code, 1)}, ...
                            tern(code == -1, 'unattached', tern(code == -2, 'rib unidentified', 'unlabelled')));
end

fprintf('1. LIBRARY %s\n', catMat);
fprintf('   DRO tau = %.4g ND -> %d-petal tulip (branch %+d); %.0f mN, Isp %g s, %g kg\n', ...
        sheet.tauDRO, sheet.Np, sheet.pm, thrustN*1000, ispS, m0kg);
fprintf('   grid %d departure x %d arrival phases (spacing 1/%d, 1/%d), %d of %d cells hold a certified root\n', ...
        numel(gridD), numel(gridA), numel(gridD), numel(gridA), nnz(has), numel(has));
fprintf('   t_f from %.2f to %.2f days\n', day(min(sheet.tf_nd(has))), day(max(sheet.tf_nd(has))));

%% ========================================================================
%  2. THE QUERY -- a phase pair, and the two endpoint states it means.
%     A phase is a fraction of its orbit's period. The endpoint is read from
%     the same periodic interpolant the library's entries were solved
%     against (costate_common/phase_state), rebuilt from the catalog's keys.
%% ========================================================================
sD = 0.0000;                                  % departure phase, fraction of tau_D
sA = 0.2050;                                  % arrival phase,  fraction of tau_A
                                              % spine96: s_D must be 0; its arrival lines are 0.0025 + j/96.
                                              % record:  (0.77, 0.26) is its cleanest cell, (0.02, 0.14)
                                              %          one that straddles two branches

orbitKeys = struct('muStar', muStar, 'lStar', lStar, 'tStar', tStar, 'tauDRO', sheet.tauDRO, ...
                   'NpTulip', sheet.Np, 'tauTulip', sheet.period_tulip_nd, 'pmTulip', sheet.pm, ...
                   'ispS', ispS, 'm0kg', m0kg);
[tD, rvD, tT, rvT] = ladder_endpoints(orbitKeys);
[stateD, seamD, dstateD] = phase_state(tD, rvD);
[stateA, seamA, dstateA] = phase_state(tT, rvT);
rv0 = stateD(sD);  rvf = stateA(sA);  rv0 = rv0(1:6);  rvf = rvf(1:6);

fprintf('\n2. QUERY  departure phase %.4f -> r = [%+.5f %+.5f %+.5f]\n', sD, rv0(1:3));
fprintf('          arrival   phase %.4f -> r = [%+.5f %+.5f %+.5f]\n', sA, rvf(1:3));
seamOK = max(seamD.value, seamA.value) < tol.closure && max(seamD.deriv, seamA.deriv) < tol.seamDeriv;
fprintf('   periodic interpolant seam: value %.1e / %.1e, derivative %.1e / %.1e   %s\n', ...
        seamD.value, seamA.value, seamD.deriv, seamA.deriv, passText(seamOK));
assert(seamOK, 'the endpoint interpolant is not periodic across its seam');

%% ========================================================================
%  3. THE CELL -- the four library entries around the query, and which of
%     them may be blended. Worked out here by hand; V1 then asserts that the
%     library function used the same corners with the same weights.
%
%     Two corners may be blended when
%       (a) they do not carry two DIFFERENT mapped family labels,
%       (b) their costates differ by less than tol.maxJump, and
%       (c) if they are neighbours on the grid, their flight times are
%           CONSISTENT WITH THEIR OWN SENSITIVITIES. A certified entry's
%           costates say how its minimum time changes as an endpoint slides
%           along its orbit (costate_common/phase_sensitivity):
%               dT/ds_D = +lambda(0) . dx_D/ds,  dT/ds_A = -lambda(t_f) . dx_A/ds
%           so two neighbours a, b on ONE smooth branch satisfy
%               T_b - T_a  =  (G_a + G_b)/2 * (s_b - s_a)  +  third order,
%           and two on different branches miss it by the jump in T.
%     (a) alone is not enough: a family folds back on itself in arrival
%     phase, so one label can cover two different roots one grid step apart.
%     (c) is the sharp one; it costs one flight per corner.
%% ========================================================================
% the grid lines below and above each phase, walking round the circle
kDlo = find(gridD <= mod(sD, 1), 1, 'last');  if isempty(kDlo), kDlo = numel(gridD); end
kAlo = find(gridA <= mod(sA, 1), 1, 'last');  if isempty(kAlo), kAlo = numel(gridA); end
kDhi = mod(kDlo, numel(gridD)) + 1;           kAhi = mod(kAlo, numel(gridA)) + 1;
gapD = mod(gridD(kDhi) - gridD(kDlo), 1);     gapA = mod(gridA(kAhi) - gridA(kAlo), 1);      % the cell's size, across the seam if need be
wD = 0;  wA = 0;                              % 0 on the lower line, -> 1 at the upper
if gapD > 0, wD = mod(mod(sD, 1) - gridD(kDlo), 1)/gapD; end
if gapA > 0, wA = mod(mod(sA, 1) - gridA(kAlo), 1)/gapA; end
% a library with ONE line in a phase says nothing about other values of it
assert(numel(gridD) > 1 || abs(mod(sD - gridD + 0.5, 1) - 0.5) < 1e-12, 'interp_study:offTheSheet', ...
       'this library holds the single departure phase %.4f; the query''s (%.4f) must equal it', gridD, sD);
assert(numel(gridA) > 1 || abs(mod(sA - gridA + 0.5, 1) - 0.5) < 1e-12, 'interp_study:offTheSheet', ...
       'this library holds the single arrival phase %.4f; the query''s (%.4f) must equal it', gridA, sA);
edgeCapMin = tol.edgeMinAt24*(24*[gapD gapA]).^1.5;           % per axis (section 0 says why 1.5)
cornerD = orderD([kDlo kDhi kDlo kDhi]);      % indices into the sheet, in the fixed corner order
cornerA = orderA([kAlo kAlo kAhi kAhi]);      %   (lowD,lowA) (highD,lowA) (lowD,highA) (highD,highA)
handWeights = [(1 - wD)*(1 - wA), wD*(1 - wA), (1 - wD)*wA, wD*wA];
% a library with ONE line in a phase has an interval for a cell, not a square: two of
% the four corners are then the other two written again
oneLineD = isscalar(gridD);  oneLineA = isscalar(gridA);
distinct = 1:4;  bestTier = 'bilinear';
if oneLineD, distinct = [1 3];  handWeights = [handWeights(1) + handWeights(2), handWeights(3) + handWeights(4)];  bestTier = 'linear'; end
if oneLineA, distinct = [1 2];  handWeights = [handWeights(1) + handWeights(3), handWeights(2) + handWeights(4)];  bestTier = 'linear'; end
assert(~(oneLineD && oneLineA), 'a library of one entry has nothing to blend');

fprintf('\n3. THE CELL around the query (it sits %.0f%% of the way across in departure, %.0f%% in arrival)\n', 100*wD, 100*wA);
fprintf('   corner   s_D      s_A      t_f [d]   family                weight\n');
cornerZ = NaN(8, 4);  cornerFam = zeros(1, 4);  cornerHas = false(1, 4);
for kc = 1:4
    cornerHas(kc) = has(cornerD(kc), cornerA(kc));
    if cornerHas(kc)
        cornerZ(:, kc) = sheet.z8(:, sheet.entry_index(cornerD(kc), cornerA(kc)));
        if isfield(sheet, 'family_index'), cornerFam(kc) = sheet.family_index(cornerD(kc), cornerA(kc)); end
    end
    if ~ismember(kc, distinct), continue, end
    if cornerHas(kc)
        fprintf('     %d     %.4f   %.4f   %7.3f   %-20s  %.3f\n', kc, sheet.sD_frac(cornerD(kc)), sheet.sA_frac(cornerA(kc)), ...
                day(cornerZ(8, kc)), famLabel(cornerFam(kc)), handWeights(distinct == kc));
    else
        fprintf('     %d     %.4f   %.4f      --     (no certified root in this cell)\n', kc, sheet.sD_frac(cornerD(kc)), sheet.sA_frac(cornerA(kc)));
    end
end
% the costate difference between every pair of corners, relative
handJump = NaN(4);
for ka = 1:4
    for kb = 1:4
        la = cornerZ(1:7, ka);  lb = cornerZ(1:7, kb);
        handJump(ka, kb) = sqrt(sum((la - lb).^2)) / max(sqrt(sum(la.^2)), sqrt(sum(lb.^2)));
    end
end
fprintf('   costate difference between corners (relative; cap %.2f):\n', tol.maxJump);
for ka = distinct, fprintf('      corner %d: %s\n', ka, strjoin(compose('%6.3f', handJump(ka, distinct)), ' ')); end
% each corner's two sensitivities: fly it once for lambda(t_f), then two dot products
GD = NaN(size(has));  GA = NaN(size(has));                   % maps, filled at the four corners only
for kc = find(cornerHas)
    sDc = sheet.sD_frac(cornerD(kc));  sAc = sheet.sA_frac(cornerA(kc));  x0c = stateD(sDc);
    [~, Yc] = pumpkyn.cr3bp.tfMinProp(cornerZ(8, kc), [x0c(1:6); 1; cornerZ(1:7, kc)], Tnd, cnd, muStar);
    dxD = dstateD(sDc);  dxA = dstateA(sAc);
    [GD(cornerD(kc), cornerA(kc)), GA(cornerD(kc), cornerA(kc))] = phase_sensitivity(cornerZ(1:6, kc), Yc(end, 8:13).', dxD(1:6), dxA(1:6));
end
cornerGD = GD(sub2ind(size(has), cornerD, cornerA));  cornerGA = GA(sub2ind(size(has), cornerD, cornerA));
toMin = tStar/60;
fprintf('   (c) neighbouring corners: actual t_f difference against what their sensitivities explain [min]\n');
stepText = @(g) tern(g > 0, sprintf('1/%.0f', 1/max(g, eps)), 'a single line');
fprintf('       (cap %.1f min along departure, %.1f min along arrival, for this cell''s steps: %s and %s)\n', ...
        edgeCapMin(1), edgeCapMin(2), stepText(gapD), stepText(gapA));
edgeList = [1 2 1; 3 4 1; 1 3 2; 2 4 2];                     % corner, corner, axis (1 departure, 2 arrival)
handEdge = NaN(4);
for ke = 1:4
    ka = edgeList(ke, 1);  kb = edgeList(ke, 2);
    if ~all(ismember([ka kb], distinct)), continue, end          % an edge of zero length
    if edgeList(ke, 3) == 1, step = gapD; Gs = cornerGD; ax = 'departure'; else, step = gapA; Gs = cornerGA; ax = 'arrival'; end
    explained = 0.5*(Gs(ka) + Gs(kb))*step;
    handEdge(ka, kb) = (cornerZ(8, kb) - cornerZ(8, ka)) - explained;  handEdge(kb, ka) = handEdge(ka, kb);
    fprintf('      corners %d-%d (%-9s): actual %+8.1f, explained %+8.1f, residual %+8.1f   %s\n', ka, kb, ax, ...
            (cornerZ(8, kb) - cornerZ(8, ka))*toMin, explained*toMin, handEdge(ka, kb)*toMin, ...
            tern(abs(handEdge(ka, kb))*toMin <= edgeCapMin(edgeList(ke, 3)), 'same branch', 'ANOTHER BRANCH'));
end
mapped = cornerFam(cornerHas & cornerFam >= 1);
capOfEdge = NaN(4);  capOfEdge([5 15]) = edgeCapMin(1);  capOfEdge([2 12]) = edgeCapMin(1);      % edges 1-2, 3-4: departure
capOfEdge([9 3]) = edgeCapMin(2);  capOfEdge([14 8]) = edgeCapMin(2);                                % edges 1-3, 2-4: arrival
known = ~isnan(handEdge);
edgesOK = all(abs(handEdge(known))*toMin <= capOfEdge(known));
allBlend = all(cornerHas) && numel(unique(mapped)) <= 1 && all(handJump(:) <= tol.maxJump) && edgesOK;
fprintf('   all present %d, one mapped family at most %d, largest costate difference %.3f, edges consistent %d\n      ->  all %d corners may be blended: %s\n', ...
        all(cornerHas), numel(unique(mapped)) <= 1, max(handJump(:)), edgesOK, numel(distinct), tern(allBlend, 'yes', 'NO'));

%% ========================================================================
%  4. THE GUESS -- the blend (costate_common/phase_catalog_interp, which also
%     knows what to do when all four may NOT be blended), and a second
%     opinion on its flight time.
%
%     The second opinion: section 3's sensitivities let each corner predict
%     the query's t_f to first order, without reference to the other corners.
%% ========================================================================
[zGuess, blend] = phase_catalog_interp(sheet, sD, sA, struct('maxJump', tol.maxJump, 'dTdsD', GD, 'dTdsA', GA, ...
                                                           'maxEdgeResidual', edgeCapMin/toMin));

fprintf('\n4. THE GUESS: tier ''%s''', blend.tier);
if ~isempty(blend.reason), fprintf('  (%s)', blend.reason); end
fprintf('\n');
b1 = ismember(blend.tier, {'entry', 'bilinear', 'linear'});
isEntry = strcmp(blend.tier, 'entry');  isBest = strcmp(blend.tier, bestTier);
sameCell = isequal(blend.cell.iD, cornerD) && isequal(blend.cell.iA, cornerA) && ...
           abs(blend.cell.wD - wD) < tol.agree && abs(blend.cell.wA - wA) < tol.agree;
libEdge = blend.edgeResidual;  bothKnown = ~isnan(handEdge) & ~isnan(libEdge);      % (the function also reports the zero-length edges)
sameEdges = all(~isnan(libEdge(~isnan(handEdge)))) && all(abs(abs(handEdge(bothKnown)) - abs(libEdge(bothKnown))) < tol.agree);
v1 = sameCell && sameEdges && (isEntry || isBest == allBlend) && (~isBest || max(abs(blend.weights - handWeights)) < tol.agree);
fprintf('   B1 corners that may be blended: %d used (%s)   %s\n', numel(blend.weights), blend.tier, passText(b1));
fprintf('   V1 the cell, weights and edge residuals worked out in section 3 are the function''s   %s\n', passText(v1));
assert(v1, 'interp_study:handVsLibrary', 'section 3 and phase_catalog_interp disagree about the cell: one of them is wrong');
if ~b1
    fprintf(['\n   STOP. There is nothing here to blend: %s.\n   This is a finding about the library at (%.4f, %.4f), ' ...
             'not about the solver.\n'], tern(isempty(blend.reason), 'the cell holds a single usable entry or none', blend.reason), sD, sA);
    error('interp_study:nothingToBlend', 'no admissible blend at this phase pair (tier ''%s'')', blend.tier);
end
fprintf('   lambda_0 guess = [%s]\n', strjoin(compose('%+.6g', zGuess(1:7)'), ' '));
fprintf('   t_f guess      = %.4f d\n', day(zGuess(8)));

% Q1: each USED corner's own first-order prediction of the query's t_f
fprintf('   Q1 what each corner''s costates predict for t_f here (first order in the phase offsets):\n');
tfTaylor = NaN(1, numel(blend.weights));
for kc = 1:numel(blend.weights)
    zc = sheet.z8(:, sheet.entry_index(blend.iD(kc), blend.iA(kc)));
    sDc = sheet.sD_frac(blend.iD(kc));  sAc = sheet.sA_frac(blend.iA(kc));
    gD = GD(blend.iD(kc), blend.iA(kc));  gA = GA(blend.iD(kc), blend.iA(kc));    % section 3's sensitivities
    dsD = mod(sD - sDc + 0.5, 1) - 0.5;  dsA = mod(sA - sAc + 0.5, 1) - 0.5;      % the SHORT way round
    tfTaylor(kc) = zc(8) + gD*dsD + gA*dsA;
    fprintf('        corner (%.4f, %.4f): %.4f d %+8.1f min (dep) %+8.1f min (arr) = %.4f d\n', sDc, sAc, day(zc(8)), ...
            gD*dsD*tStar/60, gA*dsA*tStar/60, day(tfTaylor(kc)));
end
q1SpreadMin = (max(tfTaylor) - min(tfTaylor))*tStar/60;
q1BlendMin  = (zGuess(8) - blend.weights(:).'*tfTaylor(:))*tStar/60;
fprintf('        the predictions span %.1f min; the blended t_f is %+.1f min from their weighted mean\n', q1SpreadMin, q1BlendMin);
if tfFromSensitivities
    zGuess(8) = blend.weights(:).'*tfTaylor(:);
    fprintf('   t_f guess taken from the corners'' predictions instead: %.4f d\n', day(zGuess(8)));
end

%% ========================================================================
%  5-6. FLY THE GUESS, THEN POLISH IT -- one fenced call.
%     polish_costate_guess flies the guess from t = 0 (its arrival miss is
%     the first honest measure of the guess), cuts that flight into 24
%     multiple-shooting junctions and polishes with ms_tfmin. It runs on a
%     worker under a HARD wall-clock cap: a poor guess can crawl near a
%     primary inside ONE integration, where no in-process budget fires.
%% ========================================================================
phys = struct('Tnd', Tnd, 'cnd', cnd, 'muStar', muStar, 'lStar', lStar, 'tStar', tStar, 'm0kg', m0kg);
pool = capped_pool(2);
assert(~isempty(pool) || allowUnfenced, 'interp_study:noFence', ...
       'no worker pool could be had, so the solves would run unfenced in this session; set allowUnfenced = true to accept that');
polishOpts = struct('conjTest', false);
if isempty(pool)
    [zRoot, pol] = polish_costate_guess(zGuess, rv0, rvf, phys, polishOpts);  finished = true;
else
    [finished, zRoot, pol] = run_capped(pool, @polish_costate_guess, 2, capSec, zGuess, rv0, rvf, phys, polishOpts);
end
if ~finished, zRoot = NaN(8, 1);  pol = struct('converged', false, 'normR', NaN, 'iters', NaN, 'wall', capSec, 'guessMissKm', NaN, 'Y', [], 'tGrid', []); end

b2 = isfinite(pol.guessMissKm) && pol.guessMissKm < tol.guessKm;
fprintf('\n5. THE GUESS, FLOWN UNSOLVED for %.3f d, arrives %.0f km from the target\n', day(zGuess(8)), pol.guessMissKm);
fprintf('   B2 within %.0f km   %s      (a library entry flies to < 1 km; the Moon''s orbit radius is 384,400 km)\n', tol.guessKm, passText(b2));

p1 = finished && pol.converged && pol.normR < tol.R;
fprintf('\n6. POLISH: %s', tern(finished, sprintf('%g iterations in %.0f s, |R| = %.2e', pol.iters, pol.wall, pol.normR), ...
                                         sprintf('stopped by the %g s cap', capSec)));
fprintf('\n   P1 converged to |R| < %.0e   %s\n', tol.R, passText(p1));

%% ========================================================================
%  7. THE SAME BRANCH? -- a converged solve is A root. It is the root the
%     guess was aimed at only if it stayed with the corners it was blended
%     from: near the guess in its costates, and with a flight time in (or
%     near) the corners' range. A solve that slid onto another family would
%     be a correct transfer and a wrong library answer.
%% ========================================================================
p2 = false;  q2 = NaN;  q2Min = NaN;
if p1
    q2 = sqrt(sum((zRoot(1:7) - zGuess(1:7)).^2)) / sqrt(sum(zRoot(1:7).^2));
    q2Min = (zGuess(8) - zRoot(8))*tStar/60;
    usedTf = arrayfun(@(a, b) sheet.tf_nd(a, b), blend.iD, blend.iA);
    outsideDays = max([0, day(min(usedTf) - zRoot(8)), day(zRoot(8) - max(usedTf))]);
    p2 = q2 <= tol.branch && outsideDays <= tol.branchTfDays;
    fprintf('\n7. THE ROOT: t_f = %.4f d; corners used span [%.4f, %.4f] d (root is %.3f d outside)\n', ...
            day(zRoot(8)), day(min(usedTf)), day(max(usedTf)), outsideDays);
    fprintf('   lambda_0 = [%s]\n', strjoin(compose('%+.6g', zRoot(1:7)'), ' '));
    fprintf('   Q2 the guess was off by %.2f%% in the costates and %+.1f min in t_f\n', 100*q2, q2Min);
    fprintf('   P2 stayed on its corners'' branch (costates within %.0f%%, t_f within %.1f d of their range)   %s\n', ...
            100*tol.branch, tol.branchTfDays, passText(p2));
else
    fprintf('\n7. NO ROOT, so nothing to compare. This is a finding about the RESOLUTION of the grid here.\n');
end

%% ========================================================================
%  8. CERTIFY -- the gate stack every library entry passed (certify_root):
%     polish + conjugate test, the flight and its arrival, the pointwise
%     Pontryagin checks, the second solver, the sufficiency hypotheses and
%     the dense conjugate scan. transfer_study takes it apart gate by gate.
%% ========================================================================
c1 = false;  cert = struct('ok', false, 'reason', 'not run');
Bsetup = struct('problem', struct('lStar', lStar, 'tStar', tStar, 'muStar', muStar, 'thrustN', thrustN, 'ispS', ispS, ...
                                  'm0kg', m0kg, 'tauDRO', sheet.tauDRO, 'NpTulip', sheet.Np, 'pmTulip', sheet.pm), ...
                'Tnd', Tnd, 'cnd', cnd, 'mu', muStar, 'stateD', stateD, 'stateA', stateA);
if p1 && doCertify
    seedRoot = struct('tf', zRoot(8), 'tGrid', pol.tGrid, 'Y', pol.Y);
    cert = certify_root(seedRoot, rv0, rvf, Bsetup, struct('pool', pool, 'allowUnfenced', allowUnfenced, 'm0kg', m0kg));
    c1 = cert.ok;
    fprintf('\n8. CERTIFY: %s\n', tern(cert.ok, 'every gate passed', ['REFUSED -- ' cert.reason]));
    if cert.ok
        fprintf('   t_f %.4f d, Delta-V %.4f km/s, propellant %.3f kg, flown arrival %.3f km / %.2e m/s\n', ...
                cert.tfDays, cert.dvKms, cert.propellantKg, cert.flyKm, cert.flyVms);
    end
    fprintf('   C1 certified as a library entry would be   %s\n', passText(c1));
elseif p1
    fprintf('\n8. CERTIFY: switched off (doCertify = false)\n');
end

%% ========================================================================
%  9. BASELINE -- the same solve, guessed from the NEAREST entry alone. This
%     is what a library user without an interpolator would do; the
%     difference is what the blend bought.
%% ========================================================================
base = struct('ran', false);
if doBaseline
    [zNear, nearInfo] = phase_catalog_interp(sheet, sD, sA, struct('nearestOnly', true));
    if isempty(pool)
        [zB, polB] = polish_costate_guess(zNear, rv0, rvf, phys, polishOpts);  finB = true;
    else
        [finB, zB, polB] = run_capped(pool, @polish_costate_guess, 2, capSec, zNear, rv0, rvf, phys, polishOpts);
    end
    base.ran = true;  base.converged = finB && polB.converged;
    fprintf('\n9. BASELINE from the nearest entry (%.4f, %.4f) alone:\n', sheet.sD_frac(nearInfo.iD), sheet.sA_frac(nearInfo.iA));
    if finB
        fprintf('   flown unsolved it arrives %.0f km off (the blend: %.0f km); polish %s in %g iterations, %.0f s (the blend: %g, %.0f s)\n', ...
                polB.guessMissKm, pol.guessMissKm, tern(polB.converged, 'converged', 'did NOT converge'), polB.iters, polB.wall, pol.iters, pol.wall);
    else
        fprintf('   stopped by the %g s cap\n', capSec);
    end
    if base.converged && p1
        base.sameRoot = sqrt(sum((zB(1:7) - zRoot(1:7)).^2))/sqrt(sum(zRoot(1:7).^2)) < 1e-6;
        fprintf('   it reached %s (t_f %.4f d)\n', tern(base.sameRoot, 'the SAME root', 'a DIFFERENT root'), day(zB(8)));
    end
end

%% ========================================================================
% 10. VERDICT
%% ========================================================================
gateStatus = { ...
    'B1 corners to blend',   b1; ...
    'B2 guess near target',  b2; ...
    'P1 polish converged',   p1; ...
    'P2 same branch',        p2; ...
    'C1 certified',          c1 || ~doCertify; ...
    'V1 hand vs library',    v1};
failed = gateStatus(~[gateStatus{:, 2}], 1);
usable = b1 && p1 && p2 && (c1 || ~doCertify);
fprintf('\n10. VERDICT at (s_D, s_A) = (%.4f, %.4f): ', sD, sA);
if usable && doCertify
    fprintf('the library''s neighbours blend into a guess that solves this transfer, on their\n    branch, and the result certifies. t_f = %.4f d.\n', day(zRoot(8)));
elseif usable
    fprintf('the blend solves the transfer on its corners'' branch (certification was switched off).\n');
elseif p1
    fprintf('a root was found but it is NOT usable as a library answer: %s.\n', strjoin(failed', ', '));
else
    fprintf('the blend did not lead to a root here (%s).\n', strjoin(failed', ', '));
end
fprintf('    gates: %d of %d passed', nnz([gateStatus{:, 2}]), size(gateStatus, 1));
if ~isempty(failed), fprintf('; NOT passed: %s', strjoin(failed', ', ')); end
fprintf('\n    quality: Q1 corner predictions span %.1f min; Q2 guess off by %.2f%% (costates), %+.1f min (t_f)\n', q1SpreadMin, 100*q2, q2Min);

%% ------------------------------------------------------------------------
%  FIGURES. The library's flight-time map with the query and the corners it
%  used; then the transfer itself, rotatable.
%% ------------------------------------------------------------------------
figure('Color', 'w', 'Name', 'interp_study: where the query sits');
tfMap = day(sheet.tf_nd(orderD, orderA, 1));  tfMap(~has(orderD, orderA)) = NaN;   % rows and columns in phase order
imagesc(gridA, gridD, tfMap, 'AlphaData', ~isnan(tfMap));  axis xy;  hold on
cb = colorbar;  cb.Label.String = 'minimum time [days]';
plot(sheet.sA_frac(blend.cell.iA), sheet.sD_frac(blend.cell.iD), 'ws', 'MarkerSize', 11, 'LineWidth', 1.2);
plot(sheet.sA_frac(blend.iA), sheet.sD_frac(blend.iD), 'ks', 'MarkerSize', 11, 'MarkerFaceColor', 'w', 'LineWidth', 1.2);
plot(mod(sA, 1), mod(sD, 1), 'rp', 'MarkerSize', 14, 'MarkerFaceColor', 'r');
xlabel('arrival phase s_A');  ylabel('departure phase s_D');
title(sprintf('query (star), corners blended (filled), tier ''%s''', blend.tier));
if p1
    T = struct('z', zRoot, 'sD', sD, 'sA', sA, 'tfDays', day(zRoot(8)));
    if c1, T.dvKms = cert.dvKms;  T.propellantKg = cert.propellantKg;  T.finalMassKg = cert.finalMassKg; end
    P = plot_transfer_3d(T, Bsetup);
    fprintf('\n    Figure %d is the transfer, rotatable.\n', P.fig.Number);
end

%% ------------------------------------------------------------------------
function s = passText(c)
% PASSTEXT  Verdict text.  INPUTS: c.  OUTPUTS: s.
if c, s = 'PASS'; else, s = 'FAIL'; end
end

function s = tern(c, a, b)
% TERN  Inline conditional.  INPUTS: c; a; b.  OUTPUTS: s.
if c, s = a; else, s = b; end
end
