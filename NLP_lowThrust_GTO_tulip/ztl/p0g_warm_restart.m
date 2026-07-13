% P0G_WARM_RESTART  Preflight: continue P0f's budget-capped descents to a
% verdict (converge or genuinely stall).
%
% P0f found no cold convergence at 125/75/50 mN, but two attempts hit the
% deliberately tight sweep caps (150 it / 900 evals) with flag 0 -- i.e.
% STILL DESCENDING, not stalled:
%   75 mN, rescaled-MT S_e0=0.50: ||R|| = 1.42
%   50 mN, rescaled-MT S_e0=1.50: ||R|| = 6.05
% This script warm-restarts each with a real budget (600 it / 5000 evals),
% best candidate first. Any convergence = the ladder anchor.
%
% Requires: results/p0f_thrust_sweep.mat.
% Output:   results/p0g_warm_restart.mat

here = fileparts(mfilename('fullpath'));
cd(here);  setup_paths();
resDir = fullfile(here, 'results');

F = load(fullfile(resDir, 'p0f_thrust_sweep.mat'));
[rv0, rvf, P] = ztl_endpoints();
tfL = F.tfL;

warning('off', 'MATLAB:ode113:IntegrationTolNotMet');

% pick the flag-0 (budget-capped) attempts, best residual first
cand = F.sweep([F.sweep.flag] == 0);
[~, ord] = sort([cand.resNorm]);  cand = cand(ord);
assert(~isempty(cand), 'no budget-capped candidates in p0f sweep');

opts = optimoptions('lsqnonlin', ...
    'Display', 'iter', ...
    'Algorithm', 'levenberg-marquardt', ...
    'SpecifyObjectiveGradient', true, ...
    'FunctionTolerance', 1e-24, ...
    'StepTolerance', 1e-14, ...
    'MaxIterations', 600, ...
    'MaxFunctionEvaluations', 5000);

runs = struct('Tmax_mN',{},'label',{},'resNorm',{},'flag',{},'lamSol',{},'sec',{});
best = struct('resNorm', inf, 'k', 0);
fprintf('=== P0g: warm-restart %d budget-capped P0f candidates ===\n', numel(cand));
for k = 1:numel(cand)
    Tmax = cand(k).fMult*P.Tmax25;
    fprintf('--- %g mN  %s  (start ||R|| = %.3e) ---\n', ...
            cand(k).Tmax_mN, cand(k).label, cand(k).resNorm);
    resFun = @(lam) shoot_residual_energy(lam, tfL, rv0, 1, rvf, Tmax, P.c, P.muStar);
    tic;
    [lamSol, res2, ~, flag] = lsqnonlin(resFun, cand(k).lamSol, [], [], opts);
    rn = sqrt(res2);  sec = toc;
    runs(end+1) = struct('Tmax_mN', cand(k).Tmax_mN, 'label', cand(k).label, ...
        'resNorm', rn, 'flag', flag, 'lamSol', lamSol, 'sec', sec); %#ok<SAGROW>
    fprintf('  -> ||R|| = %.3e  flag=%2d  (%.0f s)\n', rn, flag, sec);
    if rn < best.resNorm, best = struct('resNorm', rn, 'k', k); end
    save(fullfile(resDir, 'p0g_warm_restart.mat'), 'runs', 'tfL');
    if rn < 1e-8
        fprintf('  CONVERGED at %g mN -- the ladder anchor.\n', cand(k).Tmax_mN);
        break
    end
end

% accounting for the best run
r = runs(best.k);
Tmax = (r.Tmax_mN/25)*P.Tmax25;
optsInt = odeset('RelTol', 1e-12, 'AbsTol', 1e-14);
[~, yI] = ode113(@lt_pmp_eom_energy, [0 tfL], [rv0(:); 1; r.lamSol], ...
                 optsInt, Tmax, P.c, P.muStar);
lamvMag = sqrt(sum(yI(:,11:13).^2, 2));
u  = min(max(Tmax*(lamvMag./yI(:,7) + yI(:,14)/P.c), 0), 1);
mF = yI(end,7);
anchor = struct('Tmax_mN', r.Tmax_mN, 'tf', tfL, 'lam0', r.lamSol, ...
    'resNorm', r.resNorm, 'flag', r.flag, 'mProp_kg', P.m0kg*(1-mF), ...
    'dV_kms', P.c*log(1/mF)*P.lStar/P.tStar, 'uMin', min(u), 'uMax', max(u), ...
    'fracSatHi', mean(u > 0.999), 'fracSatLo', mean(u < 1e-3), ...
    'rv0', rv0, 'rvf', rvf, 'P', P);
save(fullfile(resDir, 'p0g_warm_restart.mat'), 'runs', 'anchor', 'tfL');

fprintf(['\nBEST: %g mN  ||R|| = %.3e  prop = %.4f kg  dV = %.4f km/s\n' ...
         '  throttle: min %.3f  max %.3f  sat-hi %.1f%%  sat-lo %.1f%%\n'], ...
    anchor.Tmax_mN, anchor.resNorm, anchor.mProp_kg, anchor.dV_kms, ...
    anchor.uMin, anchor.uMax, 100*anchor.fracSatHi, 100*anchor.fracSatLo);
if anchor.resNorm < 1e-8
    fprintf('GATE P0g: PASS -- ladder anchor converged at %g mN.\n', anchor.Tmax_mN);
else
    fprintf('GATE P0g: FAIL -- capped descents stall short of convergence.\n');
end
