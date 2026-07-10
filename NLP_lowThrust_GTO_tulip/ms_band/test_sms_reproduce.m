% TEST_SMS_REPRODUCE  Gates C and D: Sundman-MS reproduction tests.
%
% Gate C (min-time anchor, 1.00x): seed at the converged min-time solution
% (sms_seed_mintime, eps = 1e-3), single ms_solve. Gate: ||R|| <= 1e-9,
% dV = 4.4665 +- 0.002 km/s, prop = 2.9247 +- 0.002 kg. Expect FAST
% convergence (seed at solution).
%
% Gate D (THE decisive test, 1.12x from native sigma-domain duals): seed
% from legacy_ms_f1120.mat via sms_seed_duals (M = 50), eps-march schedule
% [1 0.3 0.1 0.03 0.01 3e-3 1e-3 3e-4 1e-4] with relays (eps_march,
% prob.resFun wired). Gate: eps <= 1e-3, ||R|| <= 1e-9,
% |dV - 3.8278| < 0.005 km/s, switches == 12. This is what the time-domain
% formulation could not do (campaign 2026-07-10 entries).
%
% Select gates by defining gateSel ('C', 'D', or 'CD') before running;
% default 'CD'. error() on fail (nonzero exit under -batch).
setup_paths;
if ~exist('gateSel', 'var'), gateSel = 'CD'; end
failMsg = '';

% ---- Gate C ---------------------------------------------------------------
if contains(gateSel, 'C')
    tC = tic;
    [Zseed, prob] = sms_seed_mintime(1.00, 24, 1e-3);
    out  = ms_solve(Zseed, prob, 1e-9, 100);
    traj = sms_traj(out.Z, prob);
    fprintf(['Gate C: ||R|| = %.3e   dV = %.4f km/s   prop = %.4f kg   ' ...
             'bang %.1f%%   t(sigf)-tf = %.3e   wall %.1f min\n'], ...
            out.resNorm, traj.dV_kms, traj.prop_kg, 100*traj.bangFrac, ...
            traj.t(end) - prob.tf, toc(tC)/60);
    okC = out.success && abs(traj.dV_kms - 4.4665) < 0.002 ...
          && abs(traj.prop_kg - 2.9247) < 0.002;
    if okC
        fprintf('PASS Gate C (sms min-time anchor)\n');
    else
        failMsg = sprintf('%s Gate C: ||R||=%.3e dV=%.4f prop=%.4f;', ...
                          failMsg, out.resNorm, traj.dV_kms, traj.prop_kg);
    end
end

% ---- Gate D ---------------------------------------------------------------
if contains(gateSel, 'D') && isempty(failMsg)
    tD = tic;
    matFile = '../sundman_minfuel/results/minfuel/legacy_ms_f1120.mat';
    [Zseed, prob, info] = sms_seed_duals(matFile, 50, 1);
    fprintf(['Gate D seed: beta = %.5f  spread %.2f%%  burnAgree %.1f%%  ' ...
             'coastAgree %.1f%%  lamT relStd %.3e  node1Err %.3e  ' ...
             'arcCheckErr %.3e\n'], info.beta, info.spreadPct, ...
            100*info.burnAgree, 100*info.coastAgree, info.lamTrelStd, ...
            info.node1Err, info.arcCheckErr);
    fprintf('Gate D seed ||R(eps=1)|| = %.3e\n', norm(sms_residual(Zseed, prob)));

    best = eps_march(Zseed, prob, [1 0.3 0.1 0.03 0.01 3e-3 1e-3 3e-4 1e-4], 1e-9);

    if isempty(best.Z)
        failMsg = sprintf('%s Gate D: eps_march produced no converged step;', failMsg);
    else
        probB = prob;  probB.epsSmooth = best.eps;
        traj  = sms_traj(best.Z, probB);
        fprintf(['Gate D: eps = %.3g   ||R|| = %.3e   dV = %.4f km/s   ' ...
                 'switches = %d   bang %.1f%%   maxJD = %.3e   wall %.1f min\n'], ...
                best.eps, best.resNorm, traj.dV_kms, traj.switches, ...
                100*traj.bangFrac, traj.maxJointDefect, toc(tD)/60);
        save('sms_gateD_f1120.mat', 'best', 'probB', 'info');
        okD = best.eps <= 1e-3 && best.resNorm <= 1e-9 ...
              && abs(traj.dV_kms - 3.8278) < 0.005 && traj.switches == 12;
        if okD
            fprintf('PASS Gate D (1.12x from native sigma-domain duals)\n');
        else
            failMsg = sprintf('%s Gate D: eps=%.3g ||R||=%.3e dV=%.4f switches=%d;', ...
                              failMsg, best.eps, best.resNorm, traj.dV_kms, ...
                              traj.switches);
        end
    end
end

if isempty(failMsg)
    fprintf('PASS test_sms_reproduce (%s)\n', gateSel);
else
    error('FAIL test_sms_reproduce:%s', failMsg);
end
