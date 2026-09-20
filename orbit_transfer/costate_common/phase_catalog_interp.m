function [z8, info] = phase_catalog_interp(sheet, sD, sA, opts)
%% Purpose:
%
%   A costate GUESS for a transfer that is not in the library, blended from
%   the library entries around it. The library is a grid over two phases
%   (departure, arrival), each a point on a circle, so the grid is a torus
%   and every query sits inside one cell with four corners.
%
%   The blend is linear in the eight numbers [lambda(0); t_f]. What this
%   function mostly decides is WHICH CORNERS MAY BE BLENDED, because a cell's
%   corners need not be the same kind of solution: the library keeps the
%   fastest root per cell, and neighbours can come from different solution
%   families, or from two sheets of one folded family. The mean of two
%   different branches is a guess for neither.
%
%   Two corners may be blended when
%       (a) they do not carry two DIFFERENT mapped family codes, and
%       (b) their costates differ by less than opts.maxJump (relative), and
%       (c) -- for two corners on one EDGE of the cell, when the sensitivity
%           maps are supplied -- their flight times are consistent with
%           their own sensitivities: the trapezoid residual
%               r = (T_b - T_a) - (G_a + G_b)/2 * (s_b - s_a)
%           is below opts.maxEdgeResidual. G = dT/ds is what an entry's
%           costates say about its neighbours (phase_sensitivity); on one
%           smooth branch r is third order in the step, across a change of
%           branch it is the size of the jump in T.
%   (a) alone is not enough: a family folds back in arrival phase, so one
%   label covers several roots at one phase. (b) alone is blunt: on the
%   24 x 24 record a step inside one branch measures 0.3-0.4 and a change of
%   branch 0.7-0.9, close enough that a cap of 1 once blended two branches
%   into a guess that hung the solver. (c) separates them by an order of
%   magnitude. Unlabelled corners (codes < 1) are judged by (b) and (c).
%
%   Candidates are tried from most to least informative, and the first whose
%   corners are all present and pairwise blendable is used:
%
%       'entry'     the query is a grid point: the stored entry, untouched
%       'bilinear'  all four corners
%       'linear'    the two corners of one edge of the cell, the edge nearest
%                   the query first; the other phase is frozen at that edge
%       'nearest'   one corner, the nearest first
%       'none'      the cell holds nothing usable: NaN. A far-away entry is
%                   never substituted.
%   A library with a SINGLE line in one phase (an arrival sheet) has an
%   interval for a cell: 'linear' between its two ends is then the best tier.
%
%  ASSUMPTIONS / NOTES:
%
% • This produces a GUESS. Whether it converges, and to which root, is the
%   polish's business (interp_study, FINDINGS 88: a mean of two neighbours
%   1/96 of a phase away reached the true root 12 times in 14).
% • The costates are in the normal chart (lambda . f = -1), so their SIZE is
%   fixed by the problem and a linear blend is meaningful; in a homogeneous
%   chart it would not be.
% • The grid need not be stored in order (the hand-built record's arrival
%   phases start at the anchor's); phases are sorted here, and wrapped.
% • One root per cell: has_solution(:, :, 1) is read.
%
%% Inputs:
%
%  sheet                    struct                  a phase sheet of a costate
%                                                   catalog: .sD_frac [1 x nD]
%                                                   .sA_frac [1 x nA]
%                                                   .has_solution [nD x nA]
%                                                   .entry_index [nD x nA]
%                                                   .z8 [8 x nEntries], and
%                                                   optionally .family_index
%                                                   [nD x nA] (>= 1 a mapped
%                                                   family, < 1 unlabelled)
%
%  sD, sA                   double                  the query phases (any
%                                                   real; taken modulo one)
%
%  opts                     struct (optional)       .maxJump [0.5] largest
%                                                   relative costate
%                                                   difference two blended
%                                                   corners may have (Inf =
%                                                   no cap), .onGridTol
%                                                   [1e-12] weight below
%                                                   which the query IS a grid
%                                                   line, .dTdsD .dTdsA
%                                                   [nD x nA] the entries'
%                                                   sensitivities (time per
%                                                   unit phase; NaN where
%                                                   unknown) with
%                                                   .maxEdgeResidual (time;
%                                                   a scalar, or [along s_D,
%                                                   along s_A] when the two
%                                                   grid steps differ -- the
%                                                   residual is third order
%                                                   in the step): rule (c).
%                                                   All three or none.
%                                                   .nearestOnly [false] blend
%                                                   nothing: the baseline a
%                                                   study compares against
%
%% Outputs:
%
%  z8                       [8 x 1]                 the guess [lambda(0); t_f]
%                                                   (NaN when tier 'none')
%
%  info                     struct                  .tier, .reason (why a
%                                                   better tier was not
%                                                   used; '' if it was),
%                                                   .iD .iA .weights [1 x n]
%                                                   the corners used,
%                                                   .family [1 x n] their
%                                                   codes, .cell (the four
%                                                   corners: .iD .iA [1 x 4],
%                                                   .present .family
%                                                   .tfNd [1 x 4], .wD .wA
%                                                   the query's position in
%                                                   the cell, 0..1),
%                                                   .jump [4 x 4] relative
%                                                   costate difference
%                                                   between the cell's
%                                                   corners (NaN where one is
%                                                   missing), .jumpMax (over
%                                                   the corners USED),
%                                                   .edgeResidual [4 x 4]
%                                                   rule (c)'s r between
%                                                   corners on one edge (NaN
%                                                   elsewhere, and when no
%                                                   maps were given)
%
%% Revision History:
%  M. Casey                                                   (c) 09/20/2026
%  Copyright Coorbital Inc.
%% ------------------------ Begin Code Sequence ---------------------------

if nargin < 4, opts = struct(); end
maxJump = 0.5;      if isfield(opts, 'maxJump'),   maxJump = opts.maxJump;     end
onGridTol = 1e-12;  if isfield(opts, 'onGridTol'), onGridTol = opts.onGridTol; end
nearestOnly = isfield(opts, 'nearestOnly') && opts.nearestOnly;
assert(isscalar(sD) && isscalar(sA) && isreal(sD) && isreal(sA) && isfinite(sD) && isfinite(sA), ...
       'phase_catalog_interp:phase', 'the query phases must be real finite scalars');
assert(isscalar(maxJump) && isreal(maxJump) && maxJump >= 0, 'phase_catalog_interp:maxJump', 'opts.maxJump must be a non-negative scalar');
mapFields = isfield(opts, {'dTdsD', 'dTdsA', 'maxEdgeResidual'});
assert(all(mapFields) || ~any(mapFields), 'phase_catalog_interp:maps', 'rule (c) needs .dTdsD, .dTdsA and .maxEdgeResidual together');
useMaps = all(mapFields);
if useMaps
    capEdge = opts.maxEdgeResidual(:).';  if isscalar(capEdge), capEdge = [capEdge capEdge]; end
    assert(numel(capEdge) == 2 && isreal(capEdge) && all(capEdge >= 0), 'phase_catalog_interp:maps', ...
           'opts.maxEdgeResidual must be a non-negative scalar or a pair [along s_D, along s_A]');
    gridSize = size(sheet.has_solution(:, :, 1));
    assert(isequal(size(opts.dTdsD), gridSize) && isequal(size(opts.dTdsA), gridSize), 'phase_catalog_interp:maps', ...
           'the sensitivity maps must be %d x %d, the sheet''s grid', gridSize(1), gridSize(2));
end

% ---- the enclosing cell on the torus ------------------------------------
[kD, wD] = bracket(sheet.sD_frac, sD);          % kD(1), kD(2): the grid lines below and above; wD in [0, 1)
[kA, wA] = bracket(sheet.sA_frac, sA);
% corners in a fixed order: (low D, low A), (high D, low A), (low D, high A), (high D, high A)
cD = [kD(1) kD(2) kD(1) kD(2)];   cA = [kA(1) kA(1) kA(2) kA(2)];
wBilinear = [(1 - wD)*(1 - wA), wD*(1 - wA), (1 - wD)*wA, wD*wA];

has = logical(sheet.has_solution(:, :, 1));
present = false(1, 4);  fam = zeros(1, 4);  tfNd = NaN(1, 4);  Z = NaN(8, 4);
for kc = 1:4
    present(kc) = has(cD(kc), cA(kc));
    if ~present(kc), continue, end
    Z(:, kc) = sheet.z8(:, sheet.entry_index(cD(kc), cA(kc)));
    tfNd(kc) = Z(8, kc);
    if isfield(sheet, 'family_index'), fam(kc) = sheet.family_index(cD(kc), cA(kc)); end
end

% ---- which corners may be blended with which ----------------------------
% corners on one EDGE of the cell, and the phase that changes along it
stepD = mod(sheet.sD_frac(kD(2)) - sheet.sD_frac(kD(1)) + 0.5, 1) - 0.5;      % the SHORT way round, signed low -> high
stepA = mod(sheet.sA_frac(kA(2)) - sheet.sA_frac(kA(1)) + 0.5, 1) - 0.5;
edgeAxis = zeros(4);  edgeAxis(1, 2) = 1;  edgeAxis(3, 4) = 1;  edgeAxis(1, 3) = 2;  edgeAxis(2, 4) = 2;   % 1: along s_D, 2: along s_A
edgeResidual = NaN(4);
jump = NaN(4);  blendable = false(4);  whyNot = strings(4);
for ka = 1:4
    for kb = 1:4
        if ~(present(ka) && present(kb)), whyNot(ka, kb) = "a corner of the cell holds no solution"; continue, end
        la = Z(1:7, ka);  lb = Z(1:7, kb);
        jump(ka, kb) = sqrt(sum((la - lb).^2)) / max(sqrt(sum(la.^2)), sqrt(sum(lb.^2)));
        lo = min(ka, kb);  hi = max(ka, kb);
        if useMaps && edgeAxis(lo, hi) > 0
            if edgeAxis(lo, hi) == 1, G = opts.dTdsD;  step = stepD; else, G = opts.dTdsA;  step = stepA; end
            r = (tfNd(hi) - tfNd(lo)) - 0.5*(G(cD(lo), cA(lo)) + G(cD(hi), cA(hi)))*step;
            edgeResidual(ka, kb) = r;  capHere = capEdge(edgeAxis(lo, hi));
        end
        if fam(ka) >= 1 && fam(kb) >= 1 && fam(ka) ~= fam(kb)
            whyNot(ka, kb) = sprintf("corners belong to different families (%d and %d)", fam(ka), fam(kb));
        elseif isfinite(edgeResidual(ka, kb)) && ~(abs(edgeResidual(ka, kb)) <= capHere)
            whyNot(ka, kb) = sprintf("the corners' flight times differ by %.3g more than their sensitivities explain (cap %.3g): another branch", ...
                                     abs(edgeResidual(ka, kb)), capHere);
        elseif ~(jump(ka, kb) <= maxJump)
            whyNot(ka, kb) = sprintf("corner costates differ by %.2f, over the cap %.2f (another branch)", jump(ka, kb), maxJump);
        else
            blendable(ka, kb) = true;
        end
    end
end

% ---- the candidates, most informative first ------------------------------
% an edge: its two corners, the weights ALONG it, and the query's distance to it
edges = struct('corners', {[1 2], [3 4], [1 3], [2 4]}, ...
               'weights', {[1 - wD, wD], [1 - wD, wD], [1 - wA, wA], [1 - wA, wA]}, ...
               'distance', {wA, 1 - wA, wD, 1 - wD});
[~, order] = sort([edges.distance]);  edges = edges(order);
cornerDistance = sqrt([wD, 1 - wD, wD, 1 - wD].^2 + [wA, wA, 1 - wA, 1 - wA].^2);
[~, cornerOrder] = sort(cornerDistance);

% A library with ONE line in a phase has no cell in that direction: the four
% corners are two entries written twice, and the candidates must not blend an
% entry with itself.
oneD = kD(1) == kD(2);  oneA = kA(1) == kA(2);
if ~oneD && ~oneA
    cand = struct('tier', 'bilinear', 'corners', 1:4, 'weights', wBilinear);
    cand = [cand, arrayfun(@(e) struct('tier', 'linear', 'corners', e.corners, 'weights', e.weights), edges)];
    singles = cornerOrder;
elseif oneD && ~oneA
    cand = struct('tier', 'linear', 'corners', [1 3], 'weights', [1 - wA, wA]);   singles = pick([1 3], wA);
elseif ~oneD && oneA
    cand = struct('tier', 'linear', 'corners', [1 2], 'weights', [1 - wD, wD]);   singles = pick([1 2], wD);
else
    cand = struct('tier', {}, 'corners', {}, 'weights', {});                      singles = 1;
end
cand = [cand, arrayfun(@(kc) struct('tier', 'nearest', 'corners', kc, 'weights', 1), singles)];
if nearestOnly, cand = cand(strcmp({cand.tier}, 'nearest')); end
bestTier = 'none';  if ~isempty(cand), bestTier = cand(1).tier; end

% a query ON a grid point is the entry itself, if the cell has it
onD = wD <= onGridTol;  onA = wA <= onGridTol;
if onD && onA && ~nearestOnly, cand = [struct('tier', 'entry', 'corners', 1, 'weights', 1), cand]; end

used = [];  reason = '';
for kk = 1:numel(cand)
    c = cand(kk).corners;
    admissible = all(present(c)) && all(all(blendable(c, c)));
    if admissible, used = cand(kk);  break, end
    if isempty(reason) && strcmp(cand(kk).tier, bestTier) && ~isscalar(c)   % say why the best tier was passed over
        bad = find(~blendable(c, c) & ~eye(numel(c)), 1);  [ra, rb] = ind2sub(numel(c)*[1 1], bad);
        reason = sprintf('not %s: %s', bestTier, whyNot(c(ra), c(rb)));
    end
end

% ---- the blend --------------------------------------------------------------
cellInfo = struct('iD', cD, 'iA', cA, 'present', present, 'family', fam, 'tfNd', tfNd, 'wD', wD, 'wA', wA);
if isempty(used)
    z8 = NaN(8, 1);
    info = struct('tier', 'none', 'reason', 'the enclosing cell holds no solution at any corner', 'iD', [], 'iA', [], ...
                  'weights', [], 'family', [], 'cell', cellInfo, 'jump', jump, 'jumpMax', NaN, 'edgeResidual', edgeResidual);
    return
end
c = used.corners;
if isscalar(c), z8 = Z(:, c); else, z8 = Z(:, c)*used.weights(:); end
if strcmp(used.tier, 'entry') || strcmp(used.tier, bestTier), reason = ''; end
if nearestOnly, reason = 'the caller asked for the nearest entry alone (opts.nearestOnly)'; end
jm = jump(c, c);  jm = max(jm(:));  if isscalar(c), jm = 0; end
info = struct('tier', used.tier, 'reason', reason, 'iD', cD(c), 'iA', cA(c), 'weights', used.weights, ...
              'family', fam(c), 'cell', cellInfo, 'jump', jump, 'jumpMax', jm, 'edgeResidual', edgeResidual);
end

% ==========================================================================
function [k, w] = bracket(phases, s)
% BRACKET  The two grid lines around a phase on the circle, and where between
% them the phase sits. The grid may be stored in any order.
% INPUTS:  phases [1 x n] grid phases in [0, 1); s (double) the query.
% OUTPUTS: k [1 x 2] indices INTO phases of the line below and the line above;
%          w (double) in [0, 1), 0 on the lower line.
[p, order] = sort(mod(phases(:).', 1));
n = numel(p);
s = mod(s, 1);
below = find(p <= s, 1, 'last');
if isempty(below)                                % the query sits before the first line: the cell
    below = n;  above = 1;                       %   runs from the LAST line, one turn back, to the first
    lower = p(n) - 1;  upper = p(1);
else
    above = mod(below, n) + 1;
    lower = p(below);  upper = p(above);
    if above == 1, upper = p(1) + 1; end         % past the last line: the first line, one turn on
end
if n == 1, k = [order(1) order(1)];  w = 0;  return, end
w = (s - lower)/(upper - lower);
k = [order(below) order(above)];
end

% ==========================================================================
function order = pick(pair, w)
% PICK  The two ends of an interval, nearest first.
% INPUTS:  pair [1 x 2] corner numbers (low end, high end); w position 0..1.
% OUTPUTS: order [1 x 2].
if w <= 0.5, order = pair; else, order = pair([2 1]); end
end
