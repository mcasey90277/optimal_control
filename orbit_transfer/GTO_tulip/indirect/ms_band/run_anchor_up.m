% RUN_ANCHOR_UP  M2(a): converge the 1.01x up-anchor by eps-march from min-time.
%
% Seeds a clean min-time-costate integration over [0, 1.01*tfMin] (wide
% basin at eps=1), then drives an eps-march down to the smoothing floor.
% M=24 default; escalates ONCE to M=48 if the very first eps=1 step fails
% to converge even after within-step relay (campaign amendment C,
% 2026-07-10 — Task-6 lesson: clean seeds converge fast, so a crawl at
% eps=1 here is diagnostic, not expected).
%
% Gate: eps floor <= 1e-3, ||R|| <= 1e-9, and 3.8278 < dV < 4.4665
% (strict monotone bracket between the 1.12x and min-time anchors).
setup_paths;
factor = 1.01;
M = 24;
[Zseed, tJ] = seed_from_mintime(factor, M, 1.0);
prob = ms_problem(factor, 1.0);  prob.tJ = tJ;
best = eps_march(Zseed, prob, [], 1e-9);

if ~best.success && numel(best.history) >= 1 && ~best.history(1).converged && M == 24
    fprintf('run_anchor_up: eps=1 failed to converge at M=24 (relays=%d, ||R||=%.3e) -- escalating to M=48\n', ...
            best.history(1).relays, best.history(1).resNorm);
    M = 48;
    [Zseed, tJ] = seed_from_mintime(factor, M, 1.0);
    prob = ms_problem(factor, 1.0);  prob.tJ = tJ;
    best = eps_march(Zseed, prob, [], 1e-9);
end

if ~best.success, fprintf('FAIL run_anchor_up (eps floor %.3g)\n', best.eps); return; end
prob.epsSmooth = best.eps;
traj = ms_traj(best.Z, prob);
fprintf('1.01x: eps=%.1e ||R||=%.2e dV=%.4f sw=%d bang=%.1f%%\n', ...
        best.eps, best.resNorm, traj.dV_kms, traj.switches, 100*traj.bangFrac);
ok = traj.dV_kms > 3.8278 && traj.dV_kms < 4.4665;
Z = best.Z; eps_ = best.eps; resNorm = best.resNorm; %#ok<NASGU>
save(sprintf('msband_%.4f.mat', factor), 'Z', 'eps_', 'resNorm', 'prob', 'traj', 'factor', 'M');
if ok, fprintf('PASS run_anchor_up\n'); else, fprintf('FAIL run_anchor_up (dV bracket)\n'); end
