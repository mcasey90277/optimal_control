% RUN_S1_GATED_M40B  Gate D M=40 sharp-start, kill-robust + trace-capturing.
%
% Third run of the authorized M = 40 escalation: the two prior runs were
% killed externally mid-relay-3 (identical deterministic trajectories:
% 3.414 -> 1.419e-2 -> 1.198e-2 (15.6%) -> 1.023e-2 (14.6%) -> ~9.1e-3
% stalling). This driver reproduces the eps = 1e-2 step with
% EPS_MARCH-IDENTICAL relay semantics (relay while iteration-capped AND
% >= 10% cut, max 4 relays), but (a) saves state after every solve so an
% external kill loses <= one 200-iter block, and (b) on relay exhaustion
% performs the coordinator-requested step-quality capture at the final
% iterate INSTEAD of discarding it silently:
%   - residual row-class breakdown (terminal rv / state / costate rows)
%   - GN consistency ||J dGN + R||/||R||, cond(J)
%   - ||dGN|| and the fraction of it the last LM step actually took
%   - predicted-vs-actual residual reduction along scaled GN steps
%     (t = 1e-3, 1e-2, 0.1, 1): direct measure of the linear model's
%     validity radius (4 residual evals; diagnosis only, no new strategy)
% Guard discipline: the non-converged iterate is used for DIAGNOSIS ONLY,
% never as a warm start. If the step converges, the remaining schedule
% [3e-3 1e-3 3e-4 1e-4] continues via EPS_MARCH; gate unchanged
% (eps <= 1e-3, ||R|| <= 1e-9, |dV - 3.8278| < 0.005, switches == 12).
% RESUME NOTE: an external ~70-min watchdog killed three prior runs at
% the same wall point; this driver therefore resumes from stateMat when
% it exists (each 200-iter segment is ~20 min, inside the kill window).
setup_paths;
tD = tic;
matFile  = '../sundman_minfuel/results/minfuel/legacy_ms_f1120.mat';
stateMat = 'sms_gateD_m40b_state.mat';
[Zseed, prob, info] = sms_seed_duals(matFile, 40, 1e-2);

MAX_RELAY = 4;
if isfile(stateMat)
    S0 = load(stateMat);
    solves = S0.solves;
    out    = solves{end};
    nRelay = numel(solves) - 1;
    if numel(solves) >= 2, Rstart = solves{end-1}.resNorm;
    else, Rstart = norm(sms_residual(Zseed, prob)); end
    fprintf('M40b RESUME: %d solves done, last ||R|| = %.3e, relays used %d\n', ...
            numel(solves), out.resNorm, nRelay);
else
    fprintf('M40b seed: beta = %.5f  arcCheckErr %.3e  ||R(eps=%.3g)|| = ', ...
            info.beta, info.arcCheckErr, prob.epsSmooth);
    Rstart = norm(sms_residual(Zseed, prob));
    fprintf('%.3e\n', Rstart);
    out    = ms_solve(Zseed, prob, 1e-9, 200);
    solves = {out};  save(stateMat, 'solves', 'prob', 'info');
    nRelay = 0;
end
while ~out.success && out.flag == 0 && out.resNorm <= 0.9*Rstart && nRelay < MAX_RELAY
    nRelay = nRelay + 1;
    fprintf('m40b: relay %d/%d at eps=%.3g (||R||=%.3e, >=10%% cut vs %.3e)\n', ...
            nRelay, MAX_RELAY, prob.epsSmooth, out.resNorm, Rstart);
    Rstart = out.resNorm;
    out    = ms_solve(out.Z, prob, 1e-9, 200);
    solves{end+1} = out;  save(stateMat, 'solves', 'prob', 'info'); %#ok<SAGROW>
end

if ~out.success
    % ---- coordinator-requested step-quality capture, then BLOCKED ----
    Z = out.Z;
    [R, J] = sms_residual(Z, prob);
    M  = numel(prob.sJ) - 1;
    Rm = reshape(R(1:16*(M-1)), 16, M-1);
    Rt = R(16*(M-1)+1:end);
    fprintf('\n=== m40 sharp-start plateau capture (eps=1e-2, ||R|| = %.4e) ===\n', norm(R));
    fprintf('rowclass: r,v %.3e  m %.3e  t %.3e  lamR %.3e  lamV %.3e  lamM %.3e  lamT %.3e\n', ...
            max(max(abs(Rm(1:6,:)))), max(abs(Rm(7,:))), max(abs(Rm(8,:))), ...
            max(max(abs(Rm(9:11,:)))), max(max(abs(Rm(12:14,:)))), ...
            max(abs(Rm(15,:))), max(abs(Rm(16,:))));
    fprintf('residual norm split: continuity %.4e  terminal %.4e (rv %.3e, lamM %.3e, t %.3e)\n', ...
            norm(R(1:16*(M-1))), norm(Rt), max(abs(Rt(1:6))), abs(Rt(7)), abs(Rt(8)));
    fprintf('optimality ||J''R||_inf = %.4e\n', norm(J.'*R, inf));
    dGN    = -(J\R);
    relres = norm(J*dGN + R)/norm(R);
    sv     = svd(full(J));
    fprintf('GN consistency %.4e   ||dGN|| = %.4e   cond(J) = %.3e (svmin %.2e)\n', ...
            relres, norm(dGN), sv(1)/sv(end), sv(end));
    % fraction of dGN the last LM step took (steps from the LM trace are
    % in scaled space; recompute the LM step at this iterate instead)
    lamLM = 1e-5;                                  % end-of-trace damping scale
    DD    = spdiags(sqrt(sum(J.^2, 1)).', 0, size(J,2), size(J,2));  % Marquardt scaling
    dLM   = -((J.'*J + lamLM*(DD.'*DD)) \ (J.'*R));
    fprintf('LM step (lambda=%.0e, Marquardt-scaled): ||dLM|| = %.4e  = %.3e of ||dGN||\n', ...
            lamLM, norm(dLM), norm(dLM)/norm(dGN));
    % linear-model validity along dGN: predicted vs actual ||R(Z + t*dGN)||
    for tGN = [1e-3 1e-2 1e-1 1]
        Rpred = norm(R + tGN*(J*dGN));
        Ract  = norm(sms_residual(Z + tGN*dGN, prob));
        fprintf('GN ray t = %5.0e: ||R||_pred = %.4e   ||R||_actual = %.4e\n', ...
                tGN, Rpred, Ract);
    end
    traj = sms_traj(Z, prob);
    fprintf('traj @ iterate: dV = %.4f  mf = %.5f  switches = %d  bang %.1f%%\n', ...
            traj.dV_kms, traj.mf, traj.switches, 100*traj.bangFrac);
    save('diag_s1_gateD_m40.mat', 'Z', 'Zseed', 'solves', 'prob', 'info', ...
         'sv', 'relres', 'dGN');
    error('FAIL run_s1_gateD_m40b: eps=1e-2 step exhausted relays at ||R||=%.3e (wall %.1f min)', ...
          out.resNorm, toc(tD)/60);
end

% ---- step 1 converged: continue the schedule, apply the gate ----
fprintf('m40b: eps=1e-2 CONVERGED (||R||=%.3e); continuing schedule\n', out.resNorm);
best = eps_march(out.Z, prob, [3e-3 1e-3 3e-4 1e-4], 1e-9);
if isempty(best.Z), bestZ = out.Z; bestEps = 1e-2; bestRes = out.resNorm;
else, bestZ = best.Z; bestEps = best.eps; bestRes = best.resNorm; end
probB = prob;  probB.epsSmooth = bestEps;
traj  = sms_traj(bestZ, probB);
fprintf(['M40b Gate D: eps = %.3g   ||R|| = %.3e   dV = %.4f km/s   ' ...
         'switches = %d   bang %.1f%%   maxJD = %.3e   wall %.1f min\n'], ...
        bestEps, bestRes, traj.dV_kms, traj.switches, 100*traj.bangFrac, ...
        traj.maxJointDefect, toc(tD)/60);
save('sms_gateD_m40.mat', 'bestZ', 'bestEps', 'bestRes', 'probB', 'info');
okD = bestEps <= 1e-3 && bestRes <= 1e-9 ...
      && abs(traj.dV_kms - 3.8278) < 0.005 && traj.switches == 12;
if okD
    fprintf('PASS run_s1_gateD_m40b (1.12x from native duals, M = 40)\n');
else
    error('FAIL run_s1_gateD_m40b: eps=%.3g ||R||=%.3e dV=%.4f switches=%d', ...
          bestEps, bestRes, traj.dV_kms, traj.switches);
end
