function [rowsCold, rowsWarm, anchorRow] = nodestudy_mee()
% NODESTUDY_MEE  Node-budget convergence study at 10 N (Task 5).
%
% TWO methodologies are compared:
%
% (1) COLD: cfg.nodesPerRev fed straight into run_transfer_mee.m's own
%     guarded eps:1->0 homotopy from a FRESH constant-throttle seed (Task
%     4's own recipe) -- clean at 15/20/25 nodes/rev, but the cold eps=1
%     step becomes numerically ill-conditioned at N>=232 (30 nodes/rev):
%     diagnosed live (results/diag_node30_full.log, a standalone probe run
%     at maxIter=1500, printLevel=5) -- dual infeasibility grew from ~1e3
%     to ~1.3e16 over 662 IPOPT iterations, the regularization exponent
%     lg(rg) climbed to ~14-15 (near-singular KKT), alpha_pr collapsed to
%     ~1e-5, and the objective stalled with essentially zero net progress
%     past iter ~600. This is a near-singular KKT system the cold loose
%     start cannot iterate out of -- NOT restoration-phase thrashing (the
%     iteration log shows steady climbing, not oscillating 'r' markers).
%     Reported as a FINDING; not chased further at 30/40 nodes/rev cold.
%
% (2) WARM (solve_warm_node, local below): mesh-refine the certified
%     25/rev solution (results/MEE_M2_10N.mat) onto the target node
%     count's uniform sigma grid -- interp1 LINEAR for the continuous
%     state X and the RTN thrust-direction rows of U, interp1 NEAREST for
%     the throttle row of U (keeps bang-bang switch edges crisp instead of
%     blurring them into intermediate throttle values) -- then solve ONE
%     fixedtf eps=0 problem directly (warmTight=true, no homotopy
%     re-sweep). Same basin by construction, isolating pure discretization
%     error: the methodologically preferred convergence-study design
%     (avoids the Cartesian N=600->1200 cold-refinement basin-change
%     problem). The GATE below is evaluated against WARM values.
%
% FINDING recorded here explicitly (Task-5 brief + thrust-ladder
% implication): cold eps=1 solves become ill-conditioned at >=30
% nodes/rev at 10 N -- always warm-start/continue a dense mesh from a
% converged neighbor, never cold-solve it from a fresh constant-throttle
% seed. This VALIDATES (does not change) the thrust ladder's existing
% continuation-based architecture (process/DESIGN_thrust_ladder.md).
%
% GATE (Task-5 brief): m_f spread across {25,30,40} nodes/rev, WARM
% values (25 = the shared anchor), must be < 0.5 kg. 15/20 nodes/rev
% (both cold and warm) are reported as floor behavior, not gated.
%
% INPUTS:  none
% OUTPUTS: rowsCold [1x5], rowsWarm [1x4], anchorRow [1x1] -- struct(s)
%          with fields nodesPerRev, N, m_f_kg, switches, edge, revs,
%          maxDefect, certified, note; ALL saved to
%          results/MEE_nodestudy.mat together with mfSpread25up, gatePass.
%
% REFERENCES: [1] .superpowers/sdd/task-5-brief.md (this task's spec).
%             [2] earth_elliptic_to_geo/run_transfer_mee.m (Task 4 driver,
%                 the COLD path, reused unmodified).
%             [3] results/diag_node30_full.log (live IPOPT diagnostic
%                 backing the cold-path ill-conditioning finding).

resDir = fullfile(module_root(), 'results');

par      = kepler_lt_params(10, 1500, 2000);
tfTarget = 1.5 * 22.2248;   % ctf * tfMinAnchor -- same physics as Task 4

emptyRow = struct('nodesPerRev', 0, 'N', 0, 'm_f_kg', 0, 'switches', 0, ...
    'edge', 0, 'revs', 0, 'maxDefect', 0, 'certified', false, 'note', '');

% --- the 25/rev anchor: Task-4's certified gate run, reused everywhere ----
base = load(fullfile(resDir, 'MEE_M2_10N.mat'));
br   = base.res.report;
anchorRow = emptyRow;
anchorRow.nodesPerRev = 25;   anchorRow.N = base.res.fp.N;
anchorRow.m_f_kg = br.m_f_kg; anchorRow.switches = br.switches;
anchorRow.edge = br.edge;     anchorRow.revs = br.revs;
anchorRow.maxDefect = br.defect; anchorRow.certified = br.certified;
anchorRow.note = ['Task-4 baseline (reused, not re-solved) -- shared ' ...
    'anchor for both cold and warm tables'];

% --- COLD path: 15/20 via run_transfer_mee.m; 30/40 NOT attempted --------
coldDensities = [15 20];
rowsCold = repmat(emptyRow, 1, numel(coldDensities) + 3);
for k = 1:numel(coldDensities)
    npr = coldDensities(k);
    tag = sprintf('MEE_node%d', npr);
    res = run_transfer_mee(struct('tag', tag, 'nodesPerRev', npr));
    r = res.report;
    rowsCold(k) = struct('nodesPerRev', npr, 'N', res.fp.N, 'm_f_kg', r.m_f_kg, ...
        'switches', r.switches, 'edge', r.edge, 'revs', r.revs, ...
        'maxDefect', r.defect, 'certified', r.certified, 'note', '');
end
rowsCold(numel(coldDensities) + 1) = anchorRow;

r30cold = emptyRow;  r30cold.nodesPerRev = 30;
r30cold.note = ['ABORTED: cold eps=1 at N=232 hits a near-singular KKT ' ...
    'system (inf_du -> 1.3e16, lg(rg) -> ~14-15, alpha_pr -> ~1e-5, ' ...
    'objective stalls by iter ~600) -- see results/diag_node30_full.log. ' ...
    'Superseded by the warm-start path below.'];
r40cold = emptyRow;  r40cold.nodesPerRev = 40;
r40cold.note = ['NOT ATTEMPTED cold (30/rev cold already diagnosed ' ...
    'ill-conditioned; 40/rev is an even denser mesh, same failure mode ' ...
    'expected). Superseded by the warm-start path below.'];
rowsCold(numel(coldDensities) + 2) = r30cold;
rowsCold(numel(coldDensities) + 3) = r40cold;

% --- WARM path: mesh-refine the 25/rev solution, single eps=0 solve ------
warmDensities = [15 20 30 40];
rowsWarm = repmat(emptyRow, 1, numel(warmDensities));
for k = 1:numel(warmDensities)
    rowsWarm(k) = solve_warm_node(warmDensities(k), base.res, par, tfTarget, resDir);
end

% --- GATE: m_f spread across {25,30,40} nodes/rev, WARM values -----------
warmNpr      = [rowsWarm.nodesPerRev];
mfGate       = [anchorRow.m_f_kg, rowsWarm(warmNpr == 30).m_f_kg, ...
                 rowsWarm(warmNpr == 40).m_f_kg];
mfSpread25up = max(mfGate) - min(mfGate);
gatePass     = mfSpread25up < 0.5;

fprintf('\n--- COLD (fresh constant-throttle seed, full eps:1->0 homotopy) ---\n');
print_rows(rowsCold);
fprintf('\n--- WARM (mesh-refined from 25/rev, single eps=0 solve) ---\n');
print_rows([anchorRow, rowsWarm]);
fprintf(['\nGATE (warm, {25,30,40}): m_f spread = %.4f kg (< 0.5 required) ' ...
    '-> %s\n'], mfSpread25up, pass_str(gatePass));

save(fullfile(resDir, 'MEE_nodestudy.mat'), 'rowsCold', 'rowsWarm', ...
    'anchorRow', 'mfSpread25up', 'gatePass');
end

% ---------------------------------------------------------------------------
function row = solve_warm_node(npr, baseRes, par, tfTarget, resDir)
% SOLVE_WARM_NODE  Mesh-refine the certified 25/rev solution (baseRes, the
% .res struct saved by run_transfer_mee.m) onto npr nodes/rev via interp1
% (linear for X and U's RTN thrust-direction rows, nearest for U's
% throttle row) and solve ONE fixedtf eps=0 problem directly
% (warmTight=true, no homotopy re-sweep). Resume-safe: caches to
% resDir/MEE_node<npr>_warm.mat with a config fingerprint guard mirroring
% run_transfer_mee.m's check_cache_fp.
%
% INPUTS:  npr - target nodes/rev [scalar]; baseRes - the 25/rev res
%          struct (.sigma .fuel.X .fuel.U .fuel.dL .seed.nRev); par -
%          kepler_lt_params struct; tfTarget - fixed tf [ND]; resDir - cache dir
% OUTPUTS: row - struct, fields nodesPerRev, N, m_f_kg, switches, edge,
%          revs, maxDefect, certified, note
N       = round(npr * baseRes.seed.nRev);
tag     = sprintf('MEE_node%d_warm', npr);
outFile = fullfile(resDir, [tag '.mat']);
fp = struct('thrustN', 10, 'ctf', 1.5, 'tfTarget', tfTarget, ...
    'nodesPerRev', npr, 'N', N, 'baseTag', 'MEE_M2_10N');

if exist(outFile, 'file')
    S = load(outFile);
    o = S.o;
    check_cache_fp_local(S, fp, outFile, tag);
else
    sigmaNew = linspace(0, 1, N + 1).';
    Xnew   = interp1(baseRes.sigma, baseRes.fuel.X.',    sigmaNew, 'linear').';
    Ubeta  = interp1(baseRes.sigma, baseRes.fuel.U(1:3,:).', sigmaNew, 'linear').';
    Uthr   = interp1(baseRes.sigma, baseRes.fuel.U(4,:).',   sigmaNew, 'nearest').';
    Unew   = [Ubeta; Uthr];
    dLnew  = baseRes.fuel.dL;
    x0     = Xnew(:,1);
    o = casadi_lt_mee(sigmaNew, Xnew, Unew, dLnew, struct('par', par, ...
        'mode', 'fixedtf', 'eps', 0, 'tfTarget', tfTarget, 'x0', x0, ...
        'maxIter', 1500, 'warmTight', true, 'printLevel', 0));
    save(outFile, 'o', 'fp');
end

certified = o.success && o.maxDefect < 1e-8;
row = struct('nodesPerRev', npr, 'N', N, 'm_f_kg', o.m_f_kg, ...
    'switches', o.switches, 'edge', o.edge, 'revs', o.dL / (2*pi), ...
    'maxDefect', o.maxDefect, 'certified', certified, ...
    'note', 'warm-started from 25/rev (single eps=0 solve, no homotopy sweep)');
fprintf(['  [warm] nodesPerRev=%d N=%d ok=%d defect=%.2e sw=%d edge=%.1f%% ' ...
    'mf=%.2f kg\n'], npr, N, certified, o.maxDefect, o.switches, 100*o.edge, o.m_f_kg);
end

% ---------------------------------------------------------------------------
function check_cache_fp_local(S, fp, file, tag)
% CHECK_CACHE_FP_LOCAL  Fail-loud cache-fingerprint guard (mirrors
% run_transfer_mee.m's check_cache_fp helper): error out, naming the first
% mismatched field, if a loaded warm-solve cache's stored fingerprint
% disagrees with the current config under the same tag.
if ~isfield(S, 'fp')
    warning('nodestudy_mee:noCachedFingerprint', ['%s has no cached ' ...
        'config fingerprint -- trusting it because tag=''%s'' matches'], ...
        file, tag);
    return;
end
flds = fieldnames(fp);
for k = 1:numel(flds)
    f = flds{k};
    if ~isfield(S.fp, f) || ~isequal(S.fp.(f), fp.(f))
        error('nodestudy_mee:fingerprintMismatch', ['cached config ' ...
            'fingerprint mismatch in %s: field ''%s'' differs -- stale ' ...
            'cache under tag=''%s''; delete the file to regenerate'], ...
            file, f, tag);
    end
end
end

% ---------------------------------------------------------------------------
function print_rows(rows)
% PRINT_ROWS  Fixed-width console table for a struct array of node-study rows.
fprintf('%-10s %-6s %-11s %-4s %-7s %-8s %-10s %-4s  %s\n', ...
    'nodes/rev', 'N', 'm_f_kg', 'sw', 'edge%', 'revs', 'maxDefect', 'cert', 'note');
for k = 1:numel(rows)
    fprintf('%-10d %-6d %-11.4f %-4d %-7.1f %-8.4f %-10.3e %-4d  %s\n', ...
        rows(k).nodesPerRev, rows(k).N, rows(k).m_f_kg, rows(k).switches, ...
        100*rows(k).edge, rows(k).revs, rows(k).maxDefect, rows(k).certified, ...
        rows(k).note);
end
end

% ---------------------------------------------------------------------------
function s = pass_str(cond)
% PASS_STR  'PASS'/'FAIL' label for a print statement.
if cond, s = 'PASS'; else, s = 'FAIL'; end
end
