% DIAG_SEED_RESIDUAL  Compute the actual MS joint residual at the
% seed_from_duals output (cheap: one ms_residual eval, no Jacobian, no LM),
% to check whether the seed is well within the "expected ~1e-3..1e-2" band
% noted in the task brief, independent of the single-arc arcCheckErr heuristic.
setup_paths;
M = 24;
matFile = '/Users/msc/Desktop/optimal_control/NLP_lowThrust_GTO_tulip/sundman_minfuel/results/minfuel/legacy_ms_f1120.mat';
[Zseed, tJ, info] = seed_from_duals(matFile, 1.12, M);
fprintf('beta=%.4g spread=%.1f%% burnAgree=%.1f%% arcCheckErr=%.2e\n', ...
        info.beta, info.spreadPct, 100*info.burnAgree, info.arcCheckErr);

prob = ms_problem(1.12, 1e-2);
prob.tJ = tJ;
tic;
R = ms_residual(Zseed, prob);
toc;
fprintf('||R|| at seed (eps=1e-2) = %.4e   (numel R = %d)\n', norm(R), numel(R));

% per-arc joint defect breakdown (14 rows per interior joint, 7 terminal)
Mn = numel(tJ) - 1;
for k = 1:Mn-1
    blk = R(14*(k-1)+(1:14));
    fprintf('  joint %2d: |defect|_inf = %.3e  (pos/vel/mass=%.3e, costate=%.3e)\n', ...
        k, max(abs(blk)), max(abs(blk(1:7))), max(abs(blk(8:14))));
end
fprintf('  terminal: |defect|_inf = %.3e\n', max(abs(R(14*(Mn-1)+(1:7)))));
