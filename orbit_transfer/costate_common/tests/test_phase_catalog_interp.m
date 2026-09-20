function ok = test_phase_catalog_interp()
% TEST_PHASE_CATALOG_INTERP  The costate interpolator on SYNTHETIC sheets whose
% answer is known in closed form: a field linear in the phases (bilinear
% interpolation is then exact), a periodic field across the seam, and sheets
% with a family boundary, a fold-like jump, and holes. No solver is involved.
%
% INPUTS:  none
% OUTPUTS: ok [logical]  every check passed (prints PASS/FAIL per check)

ok = true;
addpath(fileparts(fileparts(mfilename('fullpath'))));
nD = 6;  nA = 8;  sD = (0:nD-1)/nD;  sA = sort(mod(0.07 + (0:nA-1)/nA, 1));
lin = @(d, a) [1 + d + 2*a; 2 - d; a; 3; 4 + d - a; 5; 6; 7 + 0.5*d + 0.25*a];
S = makeSheet(sD, sA, lin);

% ---- a field linear in the phases is reproduced exactly inside a cell --------
qD = sD(2) + 0.3/nD;  qA = sA(3) + 0.6/nA;
[z, info] = phase_catalog_interp(S, qD, qA);
ok = chk(ok, strcmp(info.tier, 'bilinear') && max(abs(z - lin(qD, qA))) < 1e-13, 'a linear field is reproduced exactly (tier bilinear)');
ok = chk(ok, abs(sum(info.weights) - 1) < 1e-14 && numel(info.weights) == 4 && all(info.weights > 0), 'four positive weights that sum to one');
ok = chk(ok, isequal(sort(info.iD), [2 2 3 3]) && isequal(sort(info.iA), [3 3 4 4]), 'the corners are the enclosing cell''s');
% ---- on a grid point the entry itself comes back ------------------------------
[z, info] = phase_catalog_interp(S, sD(4), sA(5));
ok = chk(ok, strcmp(info.tier, 'entry') && isequal(z, lin(sD(4), sA(5))), 'on a grid point: the stored entry, bit for bit (tier entry)');
% ---- phases wrap, and the cell across the seam is used -------------------------
per = @(d, a) [cos(2*pi*d); sin(2*pi*d); cos(2*pi*a); sin(2*pi*a); 1; 2; 3; 6 + 0.1*cos(2*pi*(d + a))];
P = makeSheet(sD, sA, per);
qD = sD(end) + 0.5/nD;  qA = mod(sA(end) + 0.5/nA, 1);
[z, info] = phase_catalog_interp(P, qD, qA);
ok = chk(ok, isequal(sort(info.iD), [1 1 nD nD]) && isequal(sort(info.iA), [1 1 nA nA]), 'a query past the last phase uses the last and the FIRST grid lines');
ok = chk(ok, max(abs(z - per(qD, qA))) < 0.2, sprintf('and lands near the periodic field (error %.3f at this coarse grid)', max(abs(z - per(qD, qA)))));
[~, i0] = phase_catalog_interp(P, 0.5/nD, sA(1)/2);
ok = chk(ok, abs(i0.cell.wA - (sA(1)/2 - (sA(end) - 1))/(sA(1) - (sA(end) - 1))) < 1e-14 && abs(i0.cell.wD - 0.5) < 1e-14, 'a query BEFORE the first line sits in the seam cell at the right place');
[z2, ~] = phase_catalog_interp(P, qD + 3, qA - 2);
ok = chk(ok, max(abs(z - z2)) < 1e-12, 'phases are taken modulo one');
% ---- an unsorted arrival grid (the hand-built record's) gives the same answer ---
perm = [3:nA 1 2];  U = S;  U.sA_frac = S.sA_frac(perm);  U.has_solution = S.has_solution(:, perm);
U.tf_nd = S.tf_nd(:, perm);  U.entry_index = S.entry_index(:, perm);  U.family_index = S.family_index(:, perm);
qD = 0.41;  qA = 0.52;
ok = chk(ok, max(abs(phase_catalog_interp(U, qD, qA) - phase_catalog_interp(S, qD, qA))) < 1e-14, 'the grid''s storage order does not matter');

% ---- two mapped families never blend --------------------------------------------
F = S;  F.family_index(:, 4:end) = 3;                        % boundary between columns 3 and 4
qD = sD(2) + 0.3/nD;  qA = sA(3) + 0.2/nA;                   % nearer column 3 (family 1)
[z, info] = phase_catalog_interp(F, qD, qA);
ok = chk(ok, strcmp(info.tier, 'linear') && all(info.iA == 3) && numel(info.weights) == 2, 'across a family boundary: linear along the admissible edge, on ONE side');
ok = chk(ok, max(abs(z - lin(qD, sA(3)))) < 1e-13, 'that edge''s own interpolation, the other phase frozen at the edge');
ok = chk(ok, contains(info.reason, 'famil'), ['and the reason says why: ' info.reason]);
% ---- a jump in the costates blocks a blend even inside one family -----------------
J = S;  k = J.entry_index(3, 4);  J.z8(1:7, k) = -3*J.z8(1:7, k);
[~, info] = phase_catalog_interp(J, sD(2) + 0.3/nD, sA(3) + 0.6/nA);
ok = chk(ok, ~strcmp(info.tier, 'bilinear') && ~any(info.iD == 3 & info.iA == 4), 'a corner whose costates jump is left out (same label, other branch)');
ok = chk(ok, max(info.jump(:)) > 1 && info.jumpMax < 1, sprintf('the cell''s jump is reported (%.2f), and so is the smaller one among the corners used (%.2f)', max(info.jump(:)), info.jumpMax));
[~, info] = phase_catalog_interp(J, sD(2) + 0.3/nD, sA(3) + 0.6/nA, struct('maxJump', Inf));
ok = chk(ok, strcmp(info.tier, 'bilinear'), 'the cap is the caller''s to lift (opts.maxJump = Inf)');
% ---- TIME CONSISTENCY: neighbours on one branch agree with their own sensitivities ---
% T = 7 + d/2 + a/4 has dT/ds_D = 1/2, dT/ds_A = 1/4 exactly, so every edge residual is zero.
GD = 0.5*ones(nD, nA);  GA = 0.25*ones(nD, nA);  withMaps = struct('dTdsD', GD, 'dTdsA', GA, 'maxEdgeResidual', 1e-3);
[~, info] = phase_catalog_interp(S, sD(2) + 0.3/nD, sA(3) + 0.6/nA, withMaps);
ok = chk(ok, strcmp(info.tier, 'bilinear') && max(abs(info.edgeResidual(:)), [], 'omitnan') < 1e-12, 'consistent corners: bilinear, edge residuals zero');
Tj = S;  k = Tj.entry_index(3, 4);  Tj.z8(8, k) = Tj.z8(8, k) + 0.05;  Tj.tf_nd(3, 4) = Tj.z8(8, k);   % another branch: t_f jumps, costates barely move
[~, info] = phase_catalog_interp(Tj, sD(2) + 0.3/nD, sA(3) + 0.6/nA, withMaps);
ok = chk(ok, ~strcmp(info.tier, 'bilinear') && ~any(info.iD == 3 & info.iA == 4) && contains(info.reason, 'sensitivit'), ...
         ['a corner whose t_f its neighbours'' sensitivities do not explain is left out: ' info.reason]);
ok = chk(ok, abs(max(abs(info.edgeResidual(:)), [], 'omitnan') - 0.05) < 1e-12, 'and the residual reported is the jump (0.05)');
perAxis = withMaps;  perAxis.maxEdgeResidual = [1e-3 1];          % loose along s_A only: the bad corner's s_A edge passes, its s_D edge does not
[~, info] = phase_catalog_interp(Tj, sD(2) + 0.3/nD, sA(3) + 0.6/nA, perAxis);
ok = chk(ok, ~strcmp(info.tier, 'bilinear') && contains(info.reason, '0.001'), 'the cap may differ per axis [s_D, s_A]; the tighter one still refuses');
[~, info] = phase_catalog_interp(Tj, sD(2) + 0.3/nD, sA(3) + 0.6/nA);
ok = chk(ok, strcmp(info.tier, 'bilinear') && all(isnan(info.edgeResidual(:))), 'without sensitivity maps the rule is not applied (and the residuals say NaN, not zero)');
[~, info] = phase_catalog_interp(P, sD(end) + 0.5/nD, mod(sA(end) + 0.5/nA, 1), struct('dTdsD', zeros(nD, nA), 'dTdsA', zeros(nD, nA), 'maxEdgeResidual', 10));
wantD = P.z8(8, P.entry_index(1, nA)) - P.z8(8, P.entry_index(nD, nA));
ok = chk(ok, abs(info.edgeResidual(1, 2) - wantD) < 1e-14, 'across the seam the phase step is the SHORT way round (first line minus last)');
ok = chk(ok, throws(@() phase_catalog_interp(S, 0.1, 0.2, struct('dTdsD', GD))) && throws(@() phase_catalog_interp(S, 0.1, 0.2, struct('dTdsD', GD(1:2, :), 'dTdsA', GA, 'maxEdgeResidual', 1))), ...
         'half a pair of maps, or maps of the wrong size, are refused');

% ---- unlabelled cells blend on continuity alone ------------------------------------
N = S;  N.family_index(2, 3) = -1;  N.family_index(3, 4) = -2;
[~, info] = phase_catalog_interp(N, sD(2) + 0.3/nD, sA(3) + 0.6/nA);
ok = chk(ok, strcmp(info.tier, 'bilinear'), 'unattached / unidentified cells (codes < 1) are judged by their costates, not refused');
% ---- holes ------------------------------------------------------------------------------
H = S;  H.has_solution(3, 4) = false;
[z, info] = phase_catalog_interp(H, sD(2) + 0.3/nD, sA(3) + 0.6/nA);
ok = chk(ok, strcmp(info.tier, 'linear') && all(isfinite(z)), 'one corner missing: an edge of the cell is used');
H.has_solution(2:3, 3:4) = false;
[z, info] = phase_catalog_interp(H, sD(2) + 0.3/nD, sA(3) + 0.6/nA);
ok = chk(ok, strcmp(info.tier, 'none') && all(isnan(z)) && ~isempty(info.reason), 'an empty cell: tier none, NaN, and a reason -- never a far-away entry');
H = S;  H.has_solution(2:3, 3:4) = false;  H.has_solution(3, 4) = true;
[z, info] = phase_catalog_interp(H, sD(2) + 0.3/nD, sA(3) + 0.6/nA);
ok = chk(ok, strcmp(info.tier, 'nearest') && isequal(z, lin(sD(3), sA(4))), 'one corner only: that entry, called what it is (tier nearest)');
% ---- a library with ONE line in a phase: the cell is an interval, not a square -----
One = makeSheet(0, sA, lin);
[z, info] = phase_catalog_interp(One, 0, sA(3) + 0.25/nA);
ok = chk(ok, strcmp(info.tier, 'linear') && numel(info.weights) == 2 && isequal(sort(info.iA), [3 4]) && max(abs(z - lin(0, sA(3) + 0.25/nA))) < 1e-13, ...
         'one departure line: linear between the two DISTINCT arrival neighbours (never an entry blended with itself)');
[~, info] = phase_catalog_interp(One, 0, sA(3));
ok = chk(ok, strcmp(info.tier, 'entry'), 'and on its grid point, the entry');
% ---- the baseline a caller may ask for: the nearest entry, nothing blended ------------
[z, info] = phase_catalog_interp(S, sD(2) + 0.3/nD, sA(3) + 0.6/nA, struct('nearestOnly', true));
ok = chk(ok, strcmp(info.tier, 'nearest') && isscalar(info.weights) && isequal(z, lin(sD(2), sA(4))), 'opts.nearestOnly: the nearest corner alone, even where four would blend');

% ---- inputs ----------------------------------------------------------------------------
ok = chk(ok, throws(@() phase_catalog_interp(S, NaN, 0.2)) && throws(@() phase_catalog_interp(S, 0.2, [0.1 0.2])) && throws(@() phase_catalog_interp(S, 0.1, 0.2, struct('maxJump', -1))), ...
         'a non-finite or non-scalar phase, or a negative cap, is refused');
if ok, fprintf('test_phase_catalog_interp: ALL PASS\n'); else, fprintf('test_phase_catalog_interp: FAIL\n'); end
end

function S = makeSheet(sD, sA, fcn)
% MAKESHEET  A full synthetic phase sheet from a field z8 = fcn(sD, sA).
% INPUTS: sD [1 x nD]; sA [1 x nA]; fcn handle -> [8x1].  OUTPUTS: S sheet struct.
nD = numel(sD);  nA = numel(sA);
S = struct('sD_frac', sD, 'sA_frac', sA, 'has_solution', true(nD, nA), 'tf_nd', zeros(nD, nA), ...
           'entry_index', zeros(nD, nA), 'z8', zeros(8, nD*nA), 'family_index', ones(nD, nA));
k = 0;
for kd = 1:nD
    for ka = 1:nA
        k = k + 1;  S.z8(:, k) = fcn(sD(kd), sA(ka));  S.entry_index(kd, ka) = k;  S.tf_nd(kd, ka) = S.z8(8, k);
    end
end
end

function tf = throws(f)
% THROWS  Does calling f throw?  INPUTS: f (handle).  OUTPUTS: tf.
tf = false;
try, f(); catch, tf = true; end
end

function ok = chk(ok, c, msg)
% CHK  Print one PASS/FAIL line and fold it into ok.  INPUTS: ok; c; msg.
% OUTPUTS: ok.
if c, fprintf('  PASS  %s\n', msg); else, fprintf('  FAIL  %s\n', msg); end
ok = ok && c;
end
