function test_family_map()
% TEST_FAMILY_MAP  The family map reads a family's ends off its arcs and
% attaches sheet columns and ribs to the right family, on synthetic arcs
% whose answers are known.
%
%   Two anchors, three arcs, written as arclength_arrival outputs into a
%   temporary arc folder:
%     fam A, 'up':  sA 0.10 -> 0.50, then FOLDS back to 0.30 (fold at 0.50)
%     fam A, 'dn':  sA 0.10 -> 0.02, ends by walk budget with |rho| ~ 1
%     fam B, 'up':  sA 0.60 -> 0.90, ends with |rho| -> 0 (normality lost)
%   t_f is a known linear function of the (unwrapped) phase on each arc, so
%   the interpolated attachment is exact to rounding.
%
% INPUTS:  none
% OUTPUTS: none (asserts; prints PASS)

here = fileparts(mfilename('fullpath'));
addpath(fullfile(here, '..', '..', 'DRO_tulip', 'indirect'));
arcDir = tempname;  mkdir(arcDir);
cleaner = onCleanup(@() rmdir(arcDir, 's'));
tStar = 86400*4;                       % 4 days per ND unit: t_f nd = days/4
n = 12;  nExtra = 1;                   % chart: t_f at row n-1, rho at row n

% ---- family A, up: 0.10 -> 0.50 -> 0.30, t_f = 16 + 10*sA days ---------
q = [linspace(0.10, 0.50, 41), linspace(0.49, 0.30, 20)];
tfA = @(s) 16 + 10*s;
writeArc(arcDir, 'arrival_arc_alpha_up_long.mat', q, tfA(q)/4, 0.5*ones(size(q)), ...
         struct('q', 0.50), 'nStep', n, nExtra);
% ---- family A, dn: 0.10 -> 0.02, budget, normal ------------------------
q = linspace(0.10, 0.02, 9);
writeArc(arcDir, 'arrival_arc_alpha_dn_long.mat', q, tfA(q)/4, 0.5*ones(size(q)), ...
         struct([]), 'nStep', n, nExtra);
% ---- family B, up: 0.60 -> 0.90, rho -> 0 at the end, t_f = 30 - 5*sA ---
q = linspace(0.60, 0.90, 31);
tfB = @(s) 30 - 5*s;
rho = linspace(0.4, 1e-5, 31);
writeArc(arcDir, 'arrival_arc_beta_up_long.mat', q, tfB(q)/4, rho, ...
         struct([]), 'stalled (ds < dsMin)', n, nExtra);

% ---- a sheet on a 10-point grid, sA0 = 0.05 -----------------------------
nA = 10;  sA = 0.05 + (0:nA-1)/nA;
S = struct('arcs', {{'arrival_arc_alpha_up_long.mat', 'arrival_arc_alpha_dn_long.mat', ...
                     'arrival_arc_beta_up_long.mat'}}, 'sA', sA, 'TF', nan(1, nA), ...
           'cand', {cell(1, nA)}, 'problem', struct('tStar', tStar));
% col 2 (0.15): family A at 17.5 d; also a slower root of family B (none
%   passes 0.15) at 29.0 -> unattached
S.cand{2} = [cert(sA(2), tfA(0.15)), cert(sA(2), 29.0)];  S.TF(2) = tfA(0.15);
% col 4 (0.35): family A reached TWICE (up leg 19.5, fold leg 19.5 too --
%   the same t_f, so one root) -> family A
S.cand{4} = cert(sA(4), tfA(0.35));  S.TF(4) = tfA(0.35);
% col 7 (0.65): family B at 26.75
S.cand{7} = cert(sA(7), tfB(0.65));  S.TF(7) = tfB(0.65);
% col 9 (0.85): a direct-found root 0.3 d faster than family B -> unattached
S.cand{9} = [cert(sA(9), tfB(0.85) - 0.3), cert(sA(9), tfB(0.85))];  S.TF(9) = tfB(0.85) - 0.3;

F = family_map(S, struct('arcDir', arcDir, 'labels', {{'alpha', 'A'; 'beta', 'B'}}));

% ---- families and their ends --------------------------------------------
assert(numel(F.families) == 2 && isequal(F.labels, {'A', 'B'}), 'two labelled families');
fA = F.families(1);  fB = F.families(2);
assert(all(abs(fA.sASpan - [0.02 0.50]) < 1e-12), 'A spans 0.02..0.50');
assert(abs(fA.tfDays(1) - tfA(0.02)) < 1e-9 && abs(fA.tfDays(2) - tfA(0.50)) < 1e-9, 'A t_f span');
assert(startsWith(fA.ends(2).kind, 'fold'), 'A upper end is the fold');
assert(startsWith(fA.ends(1).kind, 'walk budget'), 'A lower end is the walk budget');
assert(fA.ends(1).arc == 2 && fA.ends(2).arc == 1, 'ends name the arc that reached them');
assert(abs(fB.sASpan(2) - 0.90) < 1e-12 && startsWith(fB.ends(2).kind, 'walk ended with |rho|'), ...
       'B upper end is the normality loss');
assert(startsWith(fB.ends(1).kind, 'the anchor itself'), 'B lower end is the unwalked anchor');
assert(fB.rhoMin < 1e-4 && abs(fB.rhoMinAt - 0.90) < 1e-12, 'B normality floor located');

% ---- column attachment --------------------------------------------------
c = F.columns;
assert(c(2).family == 1 && c(2).gapDays < 1e-9, 'col 2 -> A');
assert(c(4).family == 1, 'col 4 -> A (fold leg included)');
assert(c(7).family == 2 && c(7).gapDays < 1e-9, 'col 7 -> B');
assert(c(9).family == -1, 'col 9 direct-found root is unattached');
assert(c(1).family == 0 && strcmp(c(1).label, 'none'), 'empty column is none');
assert(numel(c(2).candidates) == 2 && c(2).candidates(2).family == -1, 'the slower root at col 2 is unattached');
assert(numel(c(9).candidates) == 2 && c(9).candidates(2).family == 2, 'col 9 also holds the B root');

% ---- rib attachment: t_f of the first points, extrapolated to the spine -
% a rib off col 2's A root: t_f grows 0.05, 0.11, 0.18 d off the spine
tfSpine = tfA(0.15);
assert(F.ribFamily(sA(2), tfSpine + [0.05 0.11 0.18]) == 1, 'rib off col 2 -> A');
% a rib off col 9's direct root, two points only
assert(F.ribFamily(sA(9), tfB(0.85) - 0.3 + [0.04 0.09]) == -1, 'rib off col 9 -> unattached root');
% a rib whose spine sits between the two roots at col 9 (ambiguous) -> -2
assert(F.ribFamily(sA(9), tfB(0.85) - 0.15) == -2, 'ambiguous rib is unidentified');
% a rib far from every root -> -2
assert(F.ribFamily(sA(2), 22.0) == -2, 'rib with no root nearby is unidentified');
% a rib at an empty column -> -2
assert(F.ribFamily(sA(1), 17.0) == -2, 'rib at an empty column is unidentified');

% ---- the shipped struct carries no handles once stripped ----------------
Fs = rmfield(F, {'attach', 'ribFamily'});
assert(~any(structfun(@(v) isa(v, 'function_handle'), Fs)), 'stripped map has no handles');
assert(contains(F.text, 'FAMILY MAP (3 arcs, 2 families, 10 columns)'), 'table header');

fprintf('test_family_map: PASS\n');
end

function C = cert(sA, tfDays)
% CERT  A minimal certified candidate.  INPUTS: sA; tfDays.  OUTPUTS: C.
C = struct('ok', true, 'tfDays', tfDays, 'sA', sA, 'z', [zeros(7,1); tfDays/4], 'arc', 0, 'level', sA);
end

function writeArc(arcDir, name, q, tfNd, rho, folds, stop, n, nExtra)
% WRITEARC  Save a synthetic arclength_arrival output.  INPUTS: arcDir;
% name; q [1 x N]; tfNd [1 x N]; rho [1 x N]; folds (struct with .q or
% empty); stop; n; nExtra.
p = cell(1, numel(q));
for k = 1:numel(q)
    v = zeros(n, 1);  v(n - nExtra) = tfNd(k);  v(n) = rho(k);
    p{k} = v;
end
A = struct('q', q, 'p', {p}, 'folds', folds, 'stop', stop, 'crossings', struct([]), ...
           'anc', struct('n', n, 'nExtra', nExtra));
save(fullfile(arcDir, name), 'A');
end
