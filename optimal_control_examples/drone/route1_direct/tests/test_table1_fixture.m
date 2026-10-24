function ok = test_table1_fixture()
%% Purpose:
%
%   The committed M = 4 and M = 8 fixtures are what they claim: the stored
%   X, V re-integrate to the stored J (trapezoid, to 1e-12), the trapezoid
%   defects recomputed from the stored X, V are < 1e-9 and X meets x0, xT at
%   the ends (to 1e-12) for both out and outB, the gates
%   re-derived from the stored X, V agree with the stored report, both
%   modes are present with the multipliers Route 2 needs, and the hard
%   costs sit above both the Lagrangian bound J_LB and the rigorous
%   antipodal-pair bound (recomputed here from the closed form). The
%   paper's numbers are printed as NOTE lines, never gated. The CSV's
%   gatesRun column matches its mesh column, every hard row ran the 4-seed
%   multistart, and drone_table1's adoption decision (adopt_mesh_solution,
%   via drone_table1('adopt_mesh_solution')) is exercised on synthetic
%   reports. Fixture fields:
%   out, outB, rep, repB, geom, x0, xT, paper (the paper's row: .OCFormer,
%   .SympOCNet). An edited fixture cannot sit green.
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
%  M. Casey                                                   (c) 10/23/2026
%  Copyright Coorbital Inc.
%% ------------------------ Begin Code Sequence ---------------------------

ok = true;
here = fileparts(fileparts(mfilename('fullpath')));
addpath(here, fullfile(fileparts(here), 'drone_dynamics'));
notes  = {};

for M = [4 8]
    f = fullfile(here, 'data', sprintf('drones_direct_M%d.mat', M));
    ok = chk(ok, isfile(f), sprintf('fixture M = %d exists', M));
    if ~isfile(f), continue, end
    S = load(f);
    ok = chk(ok, all(isfield(S, {'out', 'outB', 'rep', 'repB', 'geom', 'x0', 'xT', 'paper'})), ...
             sprintf('M = %d: fixture fields present', M));
    for which = {'out', 'outB'}
        o = S.(which{1});
        h = o.T/o.N;
        L = 0.5*sum(o.V.^2, 1);
        Jre = h*sum(L) - (h/2)*(L(1) + L(end));
        if strcmp(o.mode, 'barrier')
            Jre = Jre + h*sum(o.Unodes) - (h/2)*(o.Unodes(1) + o.Unodes(end));
        end
        ok = chk(ok, abs(Jre - o.J) < 1e-12*max(1, abs(o.J)), sprintf('M = %d %s: stored J re-integrates (%.1e)', M, o.mode, abs(Jre - o.J)));
        % G3 recomputed from the STORED arrays, never from the stored scalars
        hD  = o.T/o.N;
        def = o.X(:, 2:end) - o.X(:, 1:end-1) - (hD/2)*(o.V(:, 2:end) + o.V(:, 1:end-1));
        maxDef = max(abs(def(:)));
        ok = chk(ok, maxDef < 1e-9, sprintf('M = %d %s (%s): trapezoid defect recomputed from stored X, V = %.1e (< 1e-9)', M, o.mode, which{1}, maxDef));
        bcErr = max([abs(o.X(:, 1) - S.x0(:)); abs(o.X(:, end) - S.xT(:))]);
        ok = chk(ok, bcErr < 1e-12, sprintf('M = %d %s (%s): stored X(:,1) = x0, X(:,end) = xT (%.1e < 1e-12)', M, o.mode, which{1}, bcErr));
        ok = chk(ok, isequal(size(o.lamDefect), [3*M o.N]) && all(isfinite(o.lamDefect(:))), sprintf('M = %d %s: defect multipliers present', M, o.mode));
        rep = verify_drones_direct(o, struct('mirror', false));
        ok = chk(ok, rep.ok, sprintf('M = %d %s: gates re-derived ok (%s)', M, o.mode, rep.why));
    end
    ok = chk(ok, S.outB.J >= S.out.J - 1e-9, sprintf('M = %d: barrier cost is at or above the hard cost (a positive penalty on a feasible path; gaps < 1 on these layouts, so the barrier is positive)', M));
    pb = pair_bound(S.x0, S.xT, S.out.T, S.geom.r);
    ok = chk(ok, S.out.J >= pb - 1e-9 && S.out.J >= S.out.JLB, ...
             sprintf('M = %d: hard J = %.6f >= pair bound %.6f and >= J_LB %.6f', M, S.out.J, pb, S.out.JLB));
    if M == 4
        ok = chk(ok, abs(pb - 0.2003202) < 1e-7, 'M = 4 pair bound is 0.2003202 (two head-on pairs at D = 1)');
    end
    ok = chk(ok, S.rep.ok && S.rep.G5.ran && S.rep.G4.ran, sprintf('M = %d: stored report ran mirror and mesh and is ok', M));
end

csvFile = fullfile(here, 'results', 'table1.csv');
ok = chk(ok, isfile(csvFile), 'results/table1.csv exists');
if isfile(csvFile)
    tb = readtable(csvFile);
    ok = chk(ok, all(ismember([4 8 16 32], tb.M)), 'table1.csv has M = 4, 8, 16, 32');
    ok = chk(ok, ismember('pairBound', tb.Properties.VariableNames), 'table1.csv has a pairBound column');
    hard = tb(strcmp(tb.mode, 'hard'), :);
    pbRe = zeros(height(tb), 1);
    for k = 1:height(tb)
        [a, b, g] = drone_layouts('paper', tb.M(k));
        pbRe(k) = pair_bound(a, b, 10, g.r);
    end
    ok = chk(ok, max(abs(pbRe - tb.pairBound)) < 1e-12, sprintf('table1.csv pairBound column = mutual-swap pair bound on the paper layouts (%.1e)', max(abs(pbRe - tb.pairBound))));
    ok = chk(ok, all(hard.J >= hard.pairBound - 1e-9), 'every hard row respects the pair bound');
    ok = chk(ok, all(hard.J >= hard.J_LB), 'every hard row respects J_LB');
    ok = chk(ok, all(tb.gatesOk), 'every row has gatesOk true');
    ok = chk(ok, all(ismember({'gatesRun', 'seedsTried', 'bestSeed'}, tb.Properties.VariableNames)), ...
             'table1.csv has gatesRun, seedsTried, bestSeed columns');
    if ismember('gatesRun', tb.Properties.VariableNames)
        g13 = strcmp(tb.gatesRun, 'G1-G3');  g14 = strcmp(tb.gatesRun, 'G1-G4');
        ok = chk(ok, all(g13 | g14) && isequal(g13, isnan(tb.meshRel)), ...
                 sprintf('gatesRun matches the mesh column (NaN meshRel <=> G1-G3): %d rows G1-G4, %d rows G1-G3', sum(g14), sum(g13)));
    end
    if ismember('seedsTried', tb.Properties.VariableNames)
        ok = chk(ok, all(hard.seedsTried == 4) && all(ismember(hard.bestSeed, {'u+', 'u-', 'a+', 'a-'})), ...
                 sprintf('every hard row ran the 4-seed multistart (best seeds: %s)', strjoin(hard.bestSeed', ' ')));
    end
    ok = chk(ok, ismember('mirrorOk', tb.Properties.VariableNames) && all(tb.mirrorOk(tb.M <= 16)) && ~any(isnan(tb.mirrorRel(tb.M <= 16))), ...
             'mirror gate ran and holds for M <= 16 (mirrorOk true, mirrorRel not NaN)');
    for mm = unique(tb.M)'
        jh = tb.J(tb.M == mm & strcmp(tb.mode, 'hard'));  jb = tb.J(tb.M == mm & strcmp(tb.mode, 'barrier'));
        ok = chk(ok, jb >= jh - 1e-9, sprintf('M = %d: barrier cost is at or above the hard cost (a positive penalty on a feasible path; gaps < 1 on these layouts, so the barrier is positive)', mm));
    end
    for k = find(tb.M > 16)'
        notes{end+1} = sprintf('  NOTE  M = %3d %s: mirrorRel = %.2e, mirrorOk = %d\n', tb.M(k), tb.mode{k}, tb.mirrorRel(k), tb.mirrorOk(k));
    end
    for k = 1:height(hard)
        notes{end+1} = sprintf('  NOTE  M = %3d: hard J = %.6f, OCFormer %.4f, ours - OCFormer = %+.2e, pair bound %.6f\n', ...
                               hard.M(k), hard.J(k), hard.OCFormer(k), hard.J(k) - hard.OCFormer(k), hard.pairBound(k));
    end
end
% the mesh loop's adoption decision on synthetic reports (no solve): a
% cheaper successful 2N solution is adopted; an equal, dearer, failed or
% absent one is not
adopt = drone_table1('adopt_mesh_solution');
oRow  = struct('J', 0.4043878, 'N', 100);
mk = @(ran, J, succ) struct('G4', struct('ran', ran, 'out', struct('J', J, 'success', succ, 'N', 200)));
noRun = struct('G4', struct('ran', false, 'out', []));
ok = chk(ok, adopt(oRow, mk(true, 0.4038924, true)), 'adoption: a cheaper 2N solution (0.4038924 < 0.4043878) is adopted');
ok = chk(ok, ~adopt(oRow, mk(true, 0.4043878, true)), 'adoption: an equal-cost 2N solution is not adopted');
ok = chk(ok, ~adopt(oRow, mk(true, 0.4050000, true)), 'adoption: a dearer 2N solution is not adopted (escalation decides)');
ok = chk(ok, ~adopt(oRow, mk(true, 0.3000000, false)), 'adoption: a FAILED 2N solve is never adopted, however cheap');
ok = chk(ok, ~adopt(oRow, noRun), 'adoption: no adoption when G4 did not run');

[a32, b32, g32] = drone_layouts('paper', 32);
id32 = '';
try
    solve_drones_direct(a32, b32, g32, struct('mode', 'barrier'));
catch err
    id32 = err.identifier;
end
ok = chk(ok, strcmp(id32, 'solve_drones_direct:guess'), 'barrier on paper M = 32 with no warm start throws solve_drones_direct:guess');
fprintf('%s', notes{:});

if ok, fprintf('TEST_TABLE1_FIXTURE: ALL PASS\n'); else, fprintf('TEST_TABLE1_FIXTURE: FAILURE\n'); end
end

% ------------------------------------------------------------------------
function pb = pair_bound(x0, xT, T, r)
%% Purpose:
%
%   Independent copy of antipodal_pair_bound: lower bound sum ell^2/(4T),
%   ell = 2 sqrt(D^2 - R^2) + R (pi - 2 acos(R/D)), R = 2r, over MUTUAL
%   swap pairs (x0(:,b) == xT(:,a) and x0(:,a) == xT(:,b)), each agent in
%   at most one pair; an unpaired agent adds its straight-line cost
%   |dx|^2/(2T).
%
R    = 2*r;
pb   = 0;
used = false(1, size(x0, 2));
for ia = 1:size(x0, 2)
    for ib = ia+1:size(x0, 2)
        if ~used(ia) && ~used(ib) && max(abs(x0(:, ib) - xT(:, ia))) < 1e-12 && max(abs(x0(:, ia) - xT(:, ib))) < 1e-12
            D   = sqrt(sum((x0(:, ia) - x0(:, ib)).^2));
            ell = 2*sqrt(D^2 - R^2) + R*(pi - 2*acos(R/D));
            pb  = pb + ell^2/(4*T);
            used([ia ib]) = true;
        end
    end
end
for ia = find(~used)
    pb = pb + sum((xT(:, ia) - x0(:, ia)).^2)/(2*T);
end
end

function ok = chk(ok, cond, label)
%% Purpose:
%
%   Accumulate one labelled pass/fail.
%
if cond, fprintf('  PASS  %s\n', label);
else,    fprintf('  FAIL  %s\n', label);  ok = false;
end
end
