% RUN_TASK9_FUELFIRST  Task 9's final act (controller directive,
% 2026-07-18): the anchor-free "fuel-first" path both this session and the
% predecessor agent independently recommended after the 0.5 N min-time
% anchor NLP resisted 7 distinct configurations (~3 h wall) without
% certifying (see .superpowers/sdd/task-9-report.md). Rather than a
% certified min-time anchor, this driver uses an R0-LAW ESTIMATE of tfmin
% as run_transfer_mee.m's cfg.tfMinAnchor directly -- valid because the fuel
% stage only ever CONSUMES tfMinAnchor as a scale (tf = ctf*tfMinAnchor), it
% never re-derives or certifies it itself.
%
% R0-LAW DERIVATION: across the four already-certified rungs (10/5/2.5/1 N),
% T*tfmin[ND] is nearly invariant:
%   T=10 N: tfmin=22.220578  -> R0=222.206
%   T=5  N: tfmin=44.679579  -> R0=223.398
%   T=2.5N: tfmin=89.252983  -> R0=223.132
%   T=1  N: tfmin=223.808136 -> R0=223.808
% mean R0 = 223.136 (controller directive rounds to 223.14), spread
% (max-min)/mean = (223.808-222.206)/223.14 = 0.72% < 1%. This session
% re-verified all four tfmin values directly from
% results/MEE_mintime_T{100,50,25,10}.mat (not merely re-cited from the
% predecessor's report). Estimated tfmin(T) = 223.14/T:
%   T=0.5 N -> tfmin_est=446.28 ND -> tf=ctf*tfmin_est=1.5*446.28=669.42 ND
%   T=0.2 N -> tfmin_est=1115.7 ND -> tf=1.5*1115.7=1673.55 ND
% Every rung produced by this driver carries cfg.anchorSource='R0law' in
% its saved res.cfg (NOT threaded into run_transfer_mee's cache fingerprint
% fp -- doing so would retroactively invalidate the already-certified
% 10/5/2.5/1 N caches, which predate this field) as an explicit, permanent
% flag that this rung's tf target is an ESTIMATED, not independently
% certified, min-time bound.
%
% WARM START: both rungs warm-chain from the 1 N PSR-refined fuel solution
% (results/MEE_M2_1N_PSR_psr_final.mat, Task 8/9 Step 0) via the SAME
% interp_warmstart + C-law dL rescale mechanism already proven across the
% 10->5->2.5->1 N ladder (run_transfer_mee.m's cfg.warmStart path, direct
% eps=0 entry). 0.2 N (if attempted) warm-chains from 0.5 N's OWN converged
% fuel (not from 1 N again), exactly mirroring run_ladder.m's per-rung
% chaining.
%
% RESUME-SAFE: every stage below (run_transfer_mee's seed/warmdirect/final
% caches, psr_mee_refine's per-round caches) is independently resumable;
% this script itself is idempotent (re-running skips whatever is already
% cached under its tags).
%
% REFERENCES: [1] .superpowers/sdd/task-9-report.md (the anchor-wall
%   diagnosis + R0-law recommendation this driver executes).
%   [2] run_transfer_mee.m (fuel solve). [3] psr_mee_refine.m (refinement).

here   = fileparts(mfilename('fullpath'));
resDir = fullfile(here, 'results');
addpath(here);

R0law = 223.14;   % mean T*tfmin [N.ND] across the 4 certified rungs, <1% spread

% --- base: the certified 1 N PSR-refined fuel solution ---------------------
psrFinalFile = fullfile(resDir, 'MEE_M2_1N_PSR_psr_final.mat');
Sp = load(psrFinalFile);
prevFuelSigma = Sp.out.finalSigma;
prevFuelX     = Sp.out.finalOut.X;
prevFuelU     = Sp.out.finalOut.U;
prevFuelDL    = Sp.out.finalOut.dL;
prevThrust    = 1.0;
fprintf('BASE: 1 N PSR fuel (sw=%d mf=%.4f kg, N=%d, revs=%.4f)\n', ...
    Sp.out.finalOut.switches, Sp.out.finalOut.m_f_kg, numel(prevFuelSigma) - 1, ...
    prevFuelDL / (2*pi));

results = struct();

% =========================== RUNG 1: 0.5 N ==================================
thrustN = 0.5;
tfMinAnchor = R0law / thrustN;
dLGuessFuel = prevFuelDL * (prevThrust / thrustN);
fuelTag = mee_fuel_tag(thrustN);
fprintf(['\n\n========== FUEL-FIRST RUNG: T=%g N (anchorSource=R0law) ==========\n' ...
         'tfMinAnchor_est=%.4f ND (R0law/%g), tf=ctf*tfMinAnchor=%.4f ND, ' ...
         'dLGuessFuel=%.4f rad (revsGuess=%.4f)\n'], thrustN, tfMinAnchor, thrustN, ...
        1.5*tfMinAnchor, dLGuessFuel, dLGuessFuel/(2*pi));

fuelCfg = struct('thrustN', thrustN, 'ctf', 1.5, 'tfMinAnchor', tfMinAnchor, ...
    'tag', fuelTag, 'nodesPerRev', 12, 'maxIter', 1500, 'm0kg', 1500, 'ispS', 2000, ...
    'warmStart', struct('sigma', prevFuelSigma, 'X', prevFuelX, 'U', prevFuelU, ...
        'dL', dLGuessFuel), 'anchorSource', 'R0law');
tF = tic;
fuelRes = run_transfer_mee(fuelCfg);
wallFuel = toc(tF);
fprintf('[fuel DONE] T=%g N: certified=%d N=%d mf=%.4f kg sw=%d revs=%.3f wall=%.1fs\n', ...
    thrustN, fuelRes.report.certified, numel(fuelRes.sigma) - 1, fuelRes.report.m_f_kg, ...
    fuelRes.report.switches, fuelRes.report.revs, wallFuel);

results.r05.fuelCoarse = fuelRes;
results.r05.wallFuel   = wallFuel;

if fuelRes.report.certified
    psrOpts = struct('tag', [fuelTag '_PSR'], 'maxRounds', 4, 'maxIter', 1500, ...
        'nbr', 2, 'globalEvery', 3, 'globalFactor', 1.3);
    tP = tic;
    psrOut = psr_mee_refine(fuelRes, psrOpts);
    wallPsr = toc(tP);
    fprintf(['[PSR DONE] T=%g N: stopReason=%s certified=%d Nfinal=%d mf=%.4f kg ' ...
             'sw=%d revs=%.3f wall=%.1fs\n'], thrustN, psrOut.stopReason, psrOut.certified, ...
        numel(psrOut.finalSigma) - 1, psrOut.finalOut.m_f_kg, psrOut.finalOut.switches, ...
        psrOut.finalOut.dL/(2*pi), wallPsr);
    results.r05.psr     = psrOut;
    results.r05.wallPsr = wallPsr;
else
    fprintf('[BLOCKED] T=%g N fuel coarse-base did NOT certify -- stopping here.\n', thrustN);
    save(fullfile(resDir, 'FUELFIRST_task9.mat'), 'results');
    fprintf('\n=== TASK 9 FUEL-FIRST: 0.5 N BLOCKED AT FUEL STAGE, STOPPING ===\n');
    return;
end

save(fullfile(resDir, 'FUELFIRST_task9.mat'), 'results');

% =========================== RUNG 2: 0.2 N (budget-gated) ===================
% BUG FOUND AND FIXED IN THIS EDIT (2026-07-18 session): this gate was
% originally `attempt02 = true` unconditionally, reasoning from
% wallFuel+wallPsr (the CURRENT process's own tic/toc timings) as a proxy
% for total session budget. That proxy silently broke the one time it
% mattered: the 0.5 N PSR round 4 solve ran long enough that the process
% was killed (SIGKILL, exit 137 -- likely an external/turn-boundary reap,
% not a MEX crash) right as it finished writing psr_final.mat; the
% single-retry watchdog (run_task9_fuelfirst_watchdog.sh) relaunched a FRESH
% MATLAB process that replayed every already-certified round from its OWN
% per-round cache files in ~0.1 s each (by design -- that caching is what
% makes this resumable) and therefore saw wallFuel=wallPsr~=0.1 s, NOT the
% ~1h50m of actual wall-clock the human/controller-facing session had spent
% -- so the "budget check" below always looked like 0 minutes elapsed and
% the 0.2 N attempt launched anyway, well outside the controller's ~1.5 h
% envelope, and had to be manually killed. Root cause: a fresh process's own
% tic/toc cannot see time spent in a PRIOR (crashed/killed) process. Fix:
% 0.2 N is no longer auto-attempted by this script at all -- it requires an
% explicit, separate invocation (copy this rung's block into a new script,
% or re-enable manually) once a human/controller has confirmed real
% session budget remains. This is intentionally conservative: false
% negatives (skipping 0.2 N when budget WAS available) are cheap to fix by
% rerunning; false positives (silently burning unauthorized wall-clock) are
% the failure mode that actually occurred.
attempt02 = false;
if ~attempt02
    fprintf('\n=== TASK 9 FUEL-FIRST COMPLETE (0.5 N only; 0.2 N not auto-attempted -- see comment above) ===\n');
    return;
end

thrustN2 = 0.2;
tfMinAnchor2 = R0law / thrustN2;
prevFuelDL2  = psrOut.finalOut.dL;
dLGuessFuel2 = prevFuelDL2 * (thrustN / thrustN2);
fuelTag2 = mee_fuel_tag(thrustN2);
fprintf(['\n\n========== FUEL-FIRST RUNG: T=%g N (anchorSource=R0law) ==========\n' ...
         'tfMinAnchor_est=%.4f ND, tf=%.4f ND, dLGuessFuel=%.4f rad (revsGuess=%.4f)\n'], ...
        thrustN2, tfMinAnchor2, 1.5*tfMinAnchor2, dLGuessFuel2, dLGuessFuel2/(2*pi));

fuelCfg2 = struct('thrustN', thrustN2, 'ctf', 1.5, 'tfMinAnchor', tfMinAnchor2, ...
    'tag', fuelTag2, 'nodesPerRev', 10, 'maxIter', 1500, 'm0kg', 1500, 'ispS', 2000, ...
    'warmStart', struct('sigma', psrOut.finalSigma, 'X', psrOut.finalOut.X, ...
        'U', psrOut.finalOut.U, 'dL', dLGuessFuel2), 'anchorSource', 'R0law');
tF2 = tic;
fuelRes2 = run_transfer_mee(fuelCfg2);
wallFuel2 = toc(tF2);
fprintf('[fuel DONE] T=%g N: certified=%d N=%d mf=%.4f kg sw=%d revs=%.3f wall=%.1fs\n', ...
    thrustN2, fuelRes2.report.certified, numel(fuelRes2.sigma) - 1, fuelRes2.report.m_f_kg, ...
    fuelRes2.report.switches, fuelRes2.report.revs, wallFuel2);

results.r02.fuelCoarse = fuelRes2;
results.r02.wallFuel   = wallFuel2;
save(fullfile(resDir, 'FUELFIRST_task9.mat'), 'results');

if fuelRes2.report.certified
    psrOpts2 = struct('tag', [fuelTag2 '_PSR'], 'maxRounds', 4, 'maxIter', 1500, ...
        'nbr', 2, 'globalEvery', 3, 'globalFactor', 1.3);
    tP2 = tic;
    psrOut2 = psr_mee_refine(fuelRes2, psrOpts2);
    wallPsr2 = toc(tP2);
    fprintf(['[PSR DONE] T=%g N: stopReason=%s certified=%d Nfinal=%d mf=%.4f kg ' ...
             'sw=%d revs=%.3f wall=%.1fs\n'], thrustN2, psrOut2.stopReason, psrOut2.certified, ...
        numel(psrOut2.finalSigma) - 1, psrOut2.finalOut.m_f_kg, psrOut2.finalOut.switches, ...
        psrOut2.finalOut.dL/(2*pi), wallPsr2);
    results.r02.psr     = psrOut2;
    results.r02.wallPsr = wallPsr2;
else
    fprintf('[BLOCKED] T=%g N fuel coarse-base did NOT certify.\n', thrustN2);
end

save(fullfile(resDir, 'FUELFIRST_task9.mat'), 'results');
fprintf('\n=== TASK 9 FUEL-FIRST COMPLETE ===\n');
