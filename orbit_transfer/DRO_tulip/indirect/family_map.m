function F = family_map(S, opts)
%% Purpose:
%
%   The extremal FAMILIES behind an arrival sheet, measured from the arcs
%   that built it. A pseudo-arclength arc is one walk along one connected
%   branch of the minimum-time extremal set, from one anchor in one
%   direction; the two directions from an anchor trace the same branch, so
%   a family is named by its anchor. For each family the map records what
%   the arcs measured: the arrival-phase span reached, the final-time span,
%   every fold (turning point in the arrival phase), where the chart came
%   closest to losing normality (min |rho|), and HOW each end of the span
%   was reached -- a fold (the branch turns back), the walk budget (nothing
%   is known beyond it), a stall, or the anchor itself (that direction was
%   never walked). Then it names the family that supplied each certified
%   column of the sheet, and can name the family of any rib point.
%
%   Everything is derived from the stored arcs: a family's end is the
%   extreme of the phases its arcs reached, and its kind is read off the
%   arc's own fold list and stop reason. The labels are the only input.
%
%% Inputs:
%
%  S                        struct                  sheet (sheet_from_arcs
%                                                   output with .arcs, .cand,
%                                                   .TF, .sA, .problem)
%  opts                     struct (optional)
%   .arcDir                 char                    where S.arcs live
%                                                   [<here>/results]
%   .labels                 {anchor, label; ...}    friendly family labels;
%                                                   default = anchor name
%   .tStarSec               double                  seconds per ND time
%                                                   [S.problem.tStar]
%   .rhoFloor               double                  |rho| below which the
%                                                   chart has lost normality
%                                                   [1e-3]
%   .tolFold                double                  a fold this close to a
%                                                   span end IS the end
%                                                   [2e-3]
%   .tolCos                 double                  costate identity: 1 - cos of
%                                                   the angle between initial
%                                                   costates [1e-4]
%   .tolDays                double                  t_f match tolerance when
%                                                   attaching a spine point
%                                                   to an arc [0.02]; distinct
%                                                   roots at one phase sit
%                                                   >= 0.05 d apart here,
%                                                   arcs interpolate to 1e-5
%   .tolRib                 double                  t_f tolerance for a rib's
%                                                   extrapolated spine [0.03]
%
%% Outputs:
%
%  F                        struct
%   .arcs                   [1 x nArcs] struct      per S.arcs entry: .file
%                                                   .anchor .direction
%                                                   .family (index) .nRoots
%                                                   .sAStart .sAMin .sAMax
%                                                   .tfDaysMin .tfDaysMax
%                                                   .folds (sA) .rhoMin
%                                                   .rhoMinAt .stop
%   .families               [1 x nFam] struct       .label .anchor .arcs
%                                                   (indices) .sASpan [1x2]
%                                                   .sASpanMod1 .tfDays [1x2]
%                                                   .folds .rhoMin .rhoMinAt
%                                                   .ends [1x2] struct .sA
%                                                   .kind .tfDays .arc
%   .columns                [1 x nA] struct         .sA .tfDays .family
%                                                   (code, below) .label
%                                                   .gapDays .candidates
%                                                   (every certified root at
%                                                   the column: .tfDays
%                                                   .family .label)
%   .codes                  struct                  the family_index codes:
%                                                   1..nFam a mapped family,
%                                                   -1 a certified root no
%                                                   arc passes through
%                                                   (direct-found), -2 a rib
%                                                   whose spine root could
%                                                   not be identified, 0 no
%                                                   entry
%   .labels                 {1 x nFam}              family labels, index
%                                                   order
%   .attach                 function handle         [code, gap] = attach(sA,
%                                                   tfDays) for a point on
%                                                   the sheet's spine
%   .ribFamily              function handle         code = ribFamily(sA,
%                                                   tfPts) for a rib walked
%                                                   off the spine at sA:
%                                                   tfPts = t_f (days) of
%                                                   its first 1-3 points,
%                                                   nearest the spine first
%   .text                   char                    the printed table
%
%% Revision History:
%  M. Casey                                                   (c) 09/14/2026
%  Copyright Coorbital Inc.
%% ------------------------ Begin Code Sequence ---------------------------

if nargin < 2, opts = struct(); end
d = @(f, v) fieldd(opts, f, v);
here = fileparts(mfilename('fullpath'));
arcDir = d('arcDir', fullfile(here, 'results'));
labels = d('labels', cell(0, 2));
tStar = d('tStarSec', pickTStar(S));
rhoFloor = d('rhoFloor', 1e-3);
tolFold = d('tolFold', 2e-3);
tolDays = d('tolDays', 0.02);
tolRib = d('tolRib', 0.03);
tolCos = d('tolCos', 1e-4);                 % costate identity: 1 - cos(angle) between initial costates
assert(isfield(S, 'arcs') && ~isempty(S.arcs), 'family_map: the sheet names no arcs');

% ---- the arcs, one record per S.arcs entry (the sheet's indices) --------
nArcs = numel(S.arcs);
arcs = repmat(struct('file', '', 'anchor', '', 'direction', '', 'family', 0, ...
    'nRoots', 0, 'sAStart', NaN, 'sAMin', NaN, 'sAMax', NaN, 'tfDaysMin', NaN, ...
    'tfDaysMax', NaN, 'folds', [], 'rhoMin', NaN, 'rhoMinAt', NaN, 'stop', ''), 1, nArcs);
walk = cell(1, nArcs);                      % (q, tfDays, rho) per arc for attachment
for k = 1:nArcs
    f = fullfile(arcDir, S.arcs{k});
    assert(isfile(f), 'family_map: arc %s is missing from %s', S.arcs{k}, arcDir);
    L = load(f);  A = L.A;
    tok = regexp(S.arcs{k}, '^arrival_arc_(.+?)_(dn|up)(_long)?\.mat$', 'tokens', 'once');
    assert(~isempty(tok), 'family_map: cannot read an anchor name off %s', S.arcs{k});
    iTf = A.anc.n - A.anc.nExtra;           % the homogeneous chart: t_f, then rho
    q = A.q(:).';
    tf = cellfun(@(v) v(iTf), A.p(:).') * tStar/86400;
    if A.anc.nExtra >= 1, rho = abs(cellfun(@(v) v(A.anc.n), A.p(:).')); else, rho = ones(size(q)); end
    [rmin, im] = min(rho);
    fq = [];
    if ~isempty(A.folds), fq = [A.folds.q]; end
    arcs(k).file = S.arcs{k};  arcs(k).anchor = tok{1};  arcs(k).direction = tok{2};
    arcs(k).nRoots = numel(q);  arcs(k).sAStart = q(1);
    arcs(k).sAMin = min(q);  arcs(k).sAMax = max(q);
    arcs(k).tfDaysMin = min(tf);  arcs(k).tfDaysMax = max(tf);
    arcs(k).folds = fq;  arcs(k).rhoMin = rmin;  arcs(k).rhoMinAt = q(im);
    arcs(k).stop = char(string(A.stop));
    % THE INITIAL COSTATE along the arc: rows 1:7 of every point, in the
    % homogeneous chart (= rho * the normal-chart costate, measured on the
    % anchor arc 2026-09-19), so its DIRECTION is the root's identity
    lam0 = zeros(7, numel(q));
    if A.anc.n >= 7 + A.anc.nExtra + 1, lam0 = cell2mat(cellfun(@(v) v(1:7), A.p(:).', 'UniformOutput', false)); end
    walk{k} = struct('q', q, 'tf', tf, 'rho', rho, 'lam0', lam0);
end

% ---- the families: one per anchor ---------------------------------------
anchorNames = unique({arcs.anchor}, 'stable');
nFam = numel(anchorNames);
fams = repmat(struct('label', '', 'anchor', '', 'arcs', [], 'sASpan', [NaN NaN], ...
    'sASpanMod1', [NaN NaN], 'tfDays', [NaN NaN], 'folds', [], 'rhoMin', NaN, ...
    'rhoMinAt', NaN, 'ends', struct([])), 1, nFam);
for m = 1:nFam
    ia = find(strcmp({arcs.anchor}, anchorNames{m}));
    for k = ia, arcs(k).family = m; end
    fams(m).anchor = anchorNames{m};
    fams(m).label = anchorNames{m};
    if ~isempty(labels)
        hit = find(strcmp(labels(:, 1), anchorNames{m}), 1);
        if ~isempty(hit), fams(m).label = labels{hit, 2}; end
    end
    fams(m).arcs = ia;
    lo = min([arcs(ia).sAMin]);  hi = max([arcs(ia).sAMax]);
    fams(m).sASpan = [lo, hi];
    fams(m).sASpanMod1 = [mod(lo, 1), mod(hi, 1)];
    fams(m).tfDays = [min([arcs(ia).tfDaysMin]), max([arcs(ia).tfDaysMax])];
    fams(m).folds = uniquetol(sort([arcs(ia).folds]), tolFold, 'DataScale', 1);
    [fams(m).rhoMin, ir] = min([arcs(ia).rhoMin]);
    fams(m).rhoMinAt = arcs(ia(ir)).rhoMinAt;
    fams(m).ends = [endKind(lo, ia, arcs, walk, tolFold, rhoFloor), ...
                    endKind(hi, ia, arcs, walk, tolFold, rhoFloor)];
end

% ---- attachment: which family passes through (sA, t_f) on the spine ----
attach = @(sA, tfDays, varargin) attachPoint(sA, tfDays, walk, arcs, tolDays, tolCos, varargin{:});

% ---- the sheet's columns: the winner and every certified root ----------
nA = numel(S.sA);
cols = repmat(struct('sA', NaN, 'tfDays', NaN, 'family', 0, 'label', 'none', ...
                     'gapDays', NaN, 'candidates', struct([])), 1, nA);
for j = 1:nA
    cols(j).sA = S.sA(j);  cols(j).tfDays = S.TF(j);
    c = S.cand{j};
    if isempty(c), continue, end
    ok = find([c.ok]);
    cand = repmat(struct('tfDays', NaN, 'family', 0, 'label', ''), 1, numel(ok));
    for kk = 1:numel(ok)
        cand(kk).tfDays = c(ok(kk)).tfDays;
        [fam, ~] = attach(S.sA(j), c(ok(kk)).tfDays, c(ok(kk)).z);
        [cand(kk).family, cand(kk).label] = codeOf(fam, fams);
    end
    cols(j).candidates = cand;
    if ~isfinite(S.TF(j)), continue, end
    [fam, gap] = attach(S.sA(j), S.TF(j), winnerZ(S, j));
    [cols(j).family, cols(j).label] = codeOf(fam, fams);
    cols(j).gapDays = gap;
end
codes = struct('mapped', '1..nFam: index into .families', ...
               'unattached', -1, 'ribUnidentified', -2, 'none', 0);
ribFamily = @(sA, tfPts) ribCode(sA, tfPts, S.sA, cols, tolRib);

% ---- the table (built as lines, joined once) ----------------------------
lines = {sprintf('FAMILY MAP (%d arcs, %d families, %d columns)', nArcs, nFam, nA)};
for m = 1:nFam
    fm = fams(m);
    lines{end+1} = sprintf(['  %-8s anchor %-8s sA [%.4f, %.4f] (mod 1: [%.4f, %.4f])  t_f [%.2f, %.2f] d  ' ...
        'folds %s  min|rho| %.2g at sA %.4f'], fm.label, fm.anchor, fm.sASpan, fm.sASpanMod1, ...
        fm.tfDays, mat2str(fm.folds, 4), fm.rhoMin, fm.rhoMinAt);
    for e = fm.ends
        lines{end+1} = sprintf('           end sA %.4f (t_f %.2f d): %s', e.sA, e.tfDays, e.kind);
    end
end
for j = 1:nA
    lines{end+1} = sprintf('  col %2d sA %.4f  t_f %7.3f d  %s', j, cols(j).sA, cols(j).tfDays, cols(j).label);
end
txt = [strjoin(lines, newline), newline];

F = struct('arcs', arcs, 'families', fams, 'columns', cols, 'labels', {{fams.label}}, ...
           'codes', codes, 'attach', attach, 'ribFamily', ribFamily, 'text', txt, ...
           'arcDir', arcDir, 'tStarSec', tStar, 'rhoFloor', rhoFloor, 'tolDays', tolDays, ...
           'tolRib', tolRib, 'built', char(datetime('now')));
end

function e = endKind(sAEnd, ia, arcs, walk, tolFold, rhoFloor)
% ENDKIND  How a family's span end was reached, read off the arc(s) that
% reached it.  INPUTS: sAEnd; ia (the family's arcs); arcs; walk; tolFold;
% rhoFloor.  OUTPUTS: e struct .sA .kind .tfDays .arc.
e = struct('sA', sAEnd, 'kind', 'unknown', 'tfDays', NaN, 'arc', 0);
for k = ia
    w = walk{k};
    [gap, ie] = min(abs(w.q - sAEnd));
    if gap > tolFold, continue, end          % this arc did not reach the end
    e.arc = k;  e.tfDays = w.tf(ie);
    if any(abs(arcs(k).folds - sAEnd) <= tolFold)
        e.kind = 'fold (the branch turns back in arrival phase)';
    elseif ie == numel(w.q)
        if w.rho(ie) < rhoFloor
            e.kind = sprintf('walk ended with |rho| = %.2g (normality lost; stop = %s)', w.rho(ie), arcs(k).stop);
        else
            e.kind = sprintf('walk budget (stop = %s; the branch continues, unmapped)', arcs(k).stop);
        end
    elseif ie == 1
        e.kind = 'the anchor itself (this direction was not walked)';
    else
        e.kind = 'interior extreme (unclassified turning point)';
    end
    return
end
end

function [fam, gap, cosGap] = attachPoint(sA, tfDays, walk, arcs, tolDays, tolCos, z)
% ATTACHPOINT  The family whose arc PASSES THROUGH a root: at arrival phase
% sA (mod 1) the arc has the root's flight time (within tolDays) AND, when
% the root's costates z are given, the root's initial costate direction
% (1 - cos within tolCos). Flight time alone cannot tell two families apart
% where they cross in the (phase, t_f) plane -- exactly where a new family is
% found -- and it filed a third root through such a point under whichever
% arc it met first, so that root was never anchored (FINDINGS 78, 80, 82).
% Without z, or for an arc that stores no costates, the rule is the old one.
% INPUTS: sA; tfDays; walk; arcs; tolDays; tolCos; z (optional, [7|8 x 1]
% normal-chart costates, any positive scale).  OUTPUTS: fam (0 = none); gap
% (days to the best-matching arc); cosGap (1 - cos to it; NaN when no
% costates were compared).
fam = 0;  gap = Inf;  cosGap = NaN;
haveZ = nargin >= 7 && numel(z) >= 7 && any(z(1:7) ~= 0);
best = Inf;                                   % the match is ranked by costates first, then time
for k = 1:numel(walk)
    w = walk{k};
    for lev = (floor(min(w.q)) - 1 : ceil(max(w.q)) + 1) + mod(sA, 1)
        s = w.q - lev;
        ic = find(s(1:end-1) .* s(2:end) <= 0 & (s(1:end-1) ~= s(2:end)));
        on = find(s == 0);                    % a root sitting exactly on the level (the anchor, say)
        for c = [ic, on]
            if any(c == on), a = 0; else, a = s(c) / (s(c) - s(c+1)); end
            c2 = min(c + 1, numel(w.q));
            g = abs((1 - a)*w.tf(c) + a*w.tf(c2) - tfDays);
            cg = NaN;
            lamArc = (1 - a)*w.lam0(:, c) + a*w.lam0(:, c2);
            if haveZ && any(lamArc ~= 0)
                cg = 1 - (lamArc.'*z(1:7))/(sqrt(sum(lamArc.^2))*sqrt(sum(z(1:7).^2)));
            end
            passes = g <= tolDays && (isnan(cg) || cg <= tolCos);
            score = g + 1e6*(~passes);        % any passing crossing beats every failing one
            if score < best, best = score;  gap = g;  cosGap = cg;  fam = arcs(k).family*passes; end
        end
    end
end
end

function z = winnerZ(S, j)
% WINNERZ  The costates of column j's winning root, or [] when the sheet
% stores none (then the flight-time rule applies).  INPUTS: S; j.  OUTPUTS: z.
z = [];
if isfield(S, 'Z8') && size(S.Z8, 2) >= j && all(isfinite(S.Z8(:, j))), z = S.Z8(:, j); end
end

function [code, label] = codeOf(fam, fams)
% CODEOF  Family index -> family_index code and label.  INPUTS: fam (0 =
% no arc passes through the root); fams.  OUTPUTS: code; label.
if fam > 0, code = fam;  label = fams(fam).label;
else, code = -1;  label = 'unattached (no mapped arc passes through this root)';
end
end

function code = ribCode(sA, tfPts, sAgrid, cols, tolRib)
% RIBCODE  The family of a rib from the t_f of its first points: the
% spine's t_f is extrapolated back to the spine (quadratically from three
% points, linearly from two) and matched to the certified roots at that
% column; the match must be within tolRib and unambiguous (closer than half
% the distance to the next root).  INPUTS: sA; tfPts [1 x 1..3]; sAgrid;
% cols; tolRib.  OUTPUTS: code (-2 when unidentified).
code = -2;
j = find(abs(mod(sAgrid - sA + 0.5, 1) - 0.5) < 1e-8, 1);
if isempty(j) || isempty(cols(j).candidates), return, end
tfPts = tfPts(isfinite(tfPts));
switch numel(tfPts)
    case 0, return
    case 1, tf0 = tfPts(1);
    case 2, tf0 = 2*tfPts(1) - tfPts(2);
    otherwise, tf0 = 3*tfPts(1) - 3*tfPts(2) + tfPts(3);
end
gaps = sort(abs([cols(j).candidates.tfDays] - tf0));
if gaps(1) > tolRib, return, end
if numel(gaps) > 1 && gaps(1) > 0.5*gaps(2), return, end
[~, k] = min(abs([cols(j).candidates.tfDays] - tf0));
code = cols(j).candidates(k).family;
end

function t = pickTStar(S)
% PICKTSTAR  Seconds per ND time unit from the sheet's problem identity.
% INPUTS: S.  OUTPUTS: t.
if isfield(S, 'problem') && isfield(S.problem, 'tStar'), t = S.problem.tStar;
elseif isfield(S, 'B') && isfield(S.B, 'problem') && isfield(S.B.problem, 'tStar'), t = S.B.problem.tStar;
else, error('family_map:tStar', 'the sheet carries no tStar; pass opts.tStarSec');
end
end

function v = fieldd(s, f, d_)
% FIELDD  Field with default.  INPUTS: s; f; d_.  OUTPUTS: v.
if isfield(s, f) && ~isempty(s.(f)), v = s.(f); else, v = d_; end
end
