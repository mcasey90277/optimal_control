% P0E_ENERGY_MULTISTART  Preflight: test the ladder's load-bearing premise
% directly -- is the min-ENERGY shooting basin at 200 mN wide?
%
% Supersedes the min-time seeding route (P0b/P0d: both the legacy CS+dogleg
% solver and pumpkyn's analytic-STM tfMin stall on the min-time thrust march;
% min-time is near-singular and was only ever a SEED source). The fixed-tf
% argument (tf = 1.15 x tfMin(25 mN) = 7.2344, margin guaranteed by tfMin
% monotonicity in thrust) means the energy ladder never needs min-time -- it
% only needs ONE converged energy solve at the top rung. This probe attacks
% that solve with a multistart: the validated 25 mN min-time costates
% (several throttle-scale rescales) + random draws around them + fully cold
% random draws. Each attempt is a few-rev LM solve (solve_energy_indirect,
% complex-step Jacobian) -- seconds to a minute each.
%
% VERDICT LOGIC:
%   >= 1 attempt converges (||R|| <= 1e-8) -> premise holds, Z3 top anchor
%       DONE; the converged arc becomes the Z0/Z1 ground-truth test input.
%   0/N converge -> the wide-basin premise itself is in doubt at 200 mN; that
%       is a Z3-level red flag to record BEFORE building ZTL (and the next
%       probe is a direct-side seed, not more multistarts).
%
% Output: results/p0e_energy_multistart.mat + log.

here = fileparts(mfilename('fullpath'));
cd(here);  setup_paths();
resDir = fullfile(here, 'results');
if ~exist(resDir, 'dir'), mkdir(resDir); end

[rv0, rvf, P] = ztl_endpoints();
rng(0);                                     % reproducible draws

TmaxTop = 8*P.Tmax25;                       % 200 mN
tfL     = 1.15*6.29081541876621;            % fixed ladder tf = 7.2344 ND

lamMT25 = [190.476497248065; -79.7064866984696; -0.430399154713168; ...
            0.301159446575878; 0.586671892449694; -0.00711582435720301; ...
            4.32931089137559];              % validated 25 mN min-time costates

% --- seed list --------------------------------------------------------------
% S_e = Tmax(||lam_v||/m + lam_m/c); rescale candidates put S_e(0) at
% interior/near-saturated values. Random draws: log-uniform magnitude around
% the rescaled seed, sign-perturbed.
Se0 = @(lam) TmaxTop*(norm(lam(4:6)) + lam(7)/P.c);
seeds = {};
labels = {};
for target = [0.9, 0.5, 1.5, 0.25]
    seeds{end+1} = lamMT25 * target/Se0(lamMT25);              %#ok<SAGROW>
    labels{end+1} = sprintf('rescaled-MT (S_e0=%.2f)', target); %#ok<SAGROW>
end
base = seeds{1};
for kd = 1:8                                % perturbed draws around base
    pert = base .* (1 + 0.5*randn(7,1));
    seeds{end+1} = pert;                                       %#ok<SAGROW>
    labels{end+1} = sprintf('perturbed #%d', kd);              %#ok<SAGROW>
end
for kd = 1:8                                % cold draws, matched scales
    cold = [norm(base(1:3))*randn(3,1)/sqrt(3);
            norm(base(4:6))*randn(3,1)/sqrt(3);
            abs(base(7))*randn];
    seeds{end+1} = cold;                                       %#ok<SAGROW>
    labels{end+1} = sprintf('cold #%d', kd);                   %#ok<SAGROW>
end

% --- multistart -------------------------------------------------------------
warning('off', 'MATLAB:ode113:IntegrationTolNotMet');
warning('off', 'MATLAB:ode45:IntegrationTolNotMet');
fprintf('=== P0e: min-energy multistart @ 200 mN, tf = %.4f (%d attempts) ===\n', ...
        tfL, numel(seeds));
att = struct('label',{},'resNorm',{},'flag',{},'lam0',{},'lamSol',{},'sec',{});
best = struct('resNorm', inf, 'k', 0);
for k = 1:numel(seeds)
    tic;
    try
        [lamSol, rn, flag] = solve_energy_indirect(rv0, 1, rvf, tfL, ...
                                 seeds{k}, TmaxTop, P.c, P.muStar);
    catch ME
        lamSol = nan(7,1);  rn = inf;  flag = -99;
        fprintf('  [%2d] %-22s ERROR: %s\n', k, labels{k}, ME.message);
    end
    sec = toc;
    att(end+1) = struct('label', labels{k}, 'resNorm', rn, 'flag', flag, ...
        'lam0', seeds{k}, 'lamSol', lamSol, 'sec', sec); %#ok<SAGROW>
    fprintf('  [%2d] %-22s ||R|| = %.3e  flag=%2d  (%.0f s)\n', ...
            k, labels{k}, rn, flag, sec);
    if rn < best.resNorm, best = struct('resNorm', rn, 'k', k); end
    if rn < 1e-8
        fprintf('  CONVERGED -- stopping multistart.\n');
        break
    end
end

% --- accounting for the best attempt ---------------------------------------
lamBest = att(best.k).lamSol;
optsInt = odeset('RelTol', 1e-12, 'AbsTol', 1e-14);
[~, yI] = ode113(@lt_pmp_eom_energy, [0 tfL], [rv0(:); 1; lamBest], ...
                 optsInt, TmaxTop, P.c, P.muStar);
lamvMag = sqrt(sum(yI(:,11:13).^2, 2));
u  = min(max(TmaxTop*(lamvMag./yI(:,7) + yI(:,14)/P.c), 0), 1);
mF = yI(end,7);

anchor = struct('Tmax_mN', 200, 'tf', tfL, 'lam0', lamBest, ...
    'resNorm', att(best.k).resNorm, 'flag', att(best.k).flag, ...
    'seedLabel', att(best.k).label, ...
    'mProp_kg', P.m0kg*(1-mF), 'dV_kms', P.c*log(1/mF)*P.lStar/P.tStar, ...
    'uMin', min(u), 'uMax', max(u), 'fracSatHi', mean(u > 0.999), ...
    'fracSatLo', mean(u < 1e-3), 'rv0', rv0, 'rvf', rvf, 'P', P);
save(fullfile(resDir, 'p0e_energy_multistart.mat'), 'att', 'anchor');
fprintf('saved %s\n', fullfile(resDir, 'p0e_energy_multistart.mat'));

fprintf(['\nBEST: [%d] %s  ||R|| = %.3e  prop = %.4f kg  dV = %.4f km/s\n' ...
         '  throttle: min %.3f  max %.3f  sat-hi %.1f%%  sat-lo %.1f%%\n'], ...
    best.k, anchor.seedLabel, anchor.resNorm, anchor.mProp_kg, anchor.dV_kms, ...
    anchor.uMin, anchor.uMax, 100*anchor.fracSatHi, 100*anchor.fracSatLo);
if anchor.resNorm < 1e-8
    fprintf('GATE P0e: PASS -- Z3 top anchor (min-energy @ 200 mN) CONVERGED.\n');
else
    fprintf(['GATE P0e: FAIL -- %d/%d attempts, none converged. The wide-basin\n' ...
             'premise is in doubt at 200 mN; next probe = direct-side seed.\n'], ...
            numel(att), numel(seeds));
end
