% P0F_THRUST_SWEEP  Preflight: where does COLD convergence of the min-energy
% solve first appear along the thrust axis?
%
% P0e showed the 200 mN (8x) landscape is explosive from cold seeds (20/20
% fail, ||R|| ~ 1e10-1e21: wrong costates -> Earth crash/escape within a rev).
% At 25 mN (1x) the landscape is tame but 40-rev amplification kills shooting.
% Hypothesis: an intermediate-thrust sweet spot exists (~10-20 revs) where a
% cold/structured seed converges. Any single converged level = the ladder
% anchor (march down toward 25 mN from there; up toward 200 mN if ever
% needed). Levels are probed high->low (fast integrations first).
%
% Seeds per level: throttle-rescaled 25 mN min-time costates (S_e0 targets)
% + perturbed draws. Solver: lsqnonlin LM on shoot_residual_energy (CS
% Jacobian) with tightened caps (150 iters / 900 evals) so one grinding
% attempt cannot eat the sweep; Inf/NaN seeds are skipped by a pre-check.
%
% Output: results/p0f_thrust_sweep.mat + verdict per level.

here = fileparts(mfilename('fullpath'));
cd(here);  setup_paths();
resDir = fullfile(here, 'results');
if ~exist(resDir, 'dir'), mkdir(resDir); end

[rv0, rvf, P] = ztl_endpoints();
rng(1);

tfL = 1.15*6.29081541876621;                % fixed ladder tf = 7.2344 ND
lamMT25 = [190.476497248065; -79.7064866984696; -0.430399154713168; ...
            0.301159446575878; 0.586671892449694; -0.00711582435720301; ...
            4.32931089137559];

fLevels = [5, 3, 2];                        % 125, 75, 50 mN (high->low)
nPert   = 5;

warning('off', 'MATLAB:ode113:IntegrationTolNotMet');
warning('off', 'MATLAB:ode45:IntegrationTolNotMet');

opts = optimoptions('lsqnonlin', ...
    'Display', 'off', ...
    'Algorithm', 'levenberg-marquardt', ...
    'SpecifyObjectiveGradient', true, ...
    'FunctionTolerance', 1e-24, ...
    'StepTolerance', 1e-14, ...
    'MaxIterations', 150, ...
    'MaxFunctionEvaluations', 900);

sweep = struct('fMult',{},'Tmax_mN',{},'label',{},'resNorm',{},'flag',{},'lamSol',{},'sec',{});
found = struct('fMult', [], 'lam0', [], 'resNorm', inf);

fprintf('=== P0f: cold min-energy convergence vs thrust (tf = %.4f fixed) ===\n', tfL);
for fM = fLevels
    Tmax = fM*P.Tmax25;
    Se0  = Tmax*(norm(lamMT25(4:6)) + lamMT25(7)/P.c);

    seeds = {};  labels = {};
    for target = [0.9, 0.5, 1.5]
        seeds{end+1} = lamMT25 * target/Se0;                     %#ok<SAGROW>
        labels{end+1} = sprintf('rescaled-MT S_e0=%.2f', target); %#ok<SAGROW>
    end
    base = seeds{1};
    for kd = 1:nPert
        seeds{end+1} = base .* (1 + 0.4*randn(7,1));             %#ok<SAGROW>
        labels{end+1} = sprintf('perturbed #%d', kd);            %#ok<SAGROW>
    end

    fprintf('--- thrust %.0f mN (f=%g) ---\n', 25*fM, fM);
    levelConverged = false;
    for k = 1:numel(seeds)
        resFun = @(lam) shoot_residual_energy(lam, tfL, rv0, 1, rvf, Tmax, P.c, P.muStar);
        R0 = resFun(seeds{k});
        if ~all(isfinite(R0))
            fprintf('  [%d] %-22s seed integrates to Inf/NaN -- skipped\n', k, labels{k});
            continue
        end
        tic;
        [lamSol, res2, ~, flag] = lsqnonlin(resFun, seeds{k}, [], [], opts);
        rn = sqrt(res2);  sec = toc;
        sweep(end+1) = struct('fMult', fM, 'Tmax_mN', 25*fM, 'label', labels{k}, ...
            'resNorm', rn, 'flag', flag, 'lamSol', lamSol, 'sec', sec); %#ok<SAGROW>
        fprintf('  [%d] %-22s ||R|| = %.3e  flag=%2d  (%.0f s)\n', k, labels{k}, rn, flag, sec);
        if rn < 1e-8
            fprintf('  CONVERGED at %.0f mN.\n', 25*fM);
            found = struct('fMult', fM, 'lam0', lamSol, 'resNorm', rn);
            levelConverged = true;
            break
        end
    end
    save(fullfile(resDir, 'p0f_thrust_sweep.mat'), 'sweep', 'found', 'tfL');
    if levelConverged, break; end
end

if isfinite(found.resNorm) && found.resNorm < 1e-8
    Tmax = found.fMult*P.Tmax25;
    optsInt = odeset('RelTol', 1e-12, 'AbsTol', 1e-14);
    [~, yI] = ode113(@lt_pmp_eom_energy, [0 tfL], [rv0(:); 1; found.lam0], ...
                     optsInt, Tmax, P.c, P.muStar);
    lamvMag = sqrt(sum(yI(:,11:13).^2, 2));
    u  = min(max(Tmax*(lamvMag./yI(:,7) + yI(:,14)/P.c), 0), 1);
    mF = yI(end,7);
    anchor = struct('Tmax_mN', 25*found.fMult, 'tf', tfL, 'lam0', found.lam0, ...
        'resNorm', found.resNorm, 'mProp_kg', P.m0kg*(1-mF), ...
        'dV_kms', P.c*log(1/mF)*P.lStar/P.tStar, 'uMin', min(u), 'uMax', max(u), ...
        'fracSatHi', mean(u > 0.999), 'fracSatLo', mean(u < 1e-3), ...
        'rv0', rv0, 'rvf', rvf, 'P', P);
    save(fullfile(resDir, 'p0f_thrust_sweep.mat'), 'sweep', 'found', 'anchor', 'tfL');
    fprintf(['\nGATE P0f: PASS -- ANCHOR at %.0f mN: ||R|| = %.3e  prop = %.4f kg  ' ...
             'dV = %.4f km/s\n  throttle: min %.3f max %.3f  sat-hi %.1f%%  sat-lo %.1f%%\n'], ...
        anchor.Tmax_mN, anchor.resNorm, anchor.mProp_kg, anchor.dV_kms, ...
        anchor.uMin, anchor.uMax, 100*anchor.fracSatHi, 100*anchor.fracSatLo);
else
    fprintf(['\nGATE P0f: FAIL -- no cold convergence at any probed thrust.\n' ...
             'Next: direct-side (collocation dual) seed at a high rung.\n']);
end
