% P0C_ENERGY_TOP  Preflight: free look at the Z3 top anchor -- min-energy
% indirect solve at 200 mN with the EXISTING machinery (PLAN_PRONG_Z.md P0c).
%
% Runs solve_energy_indirect (single shooting, complex-step Jacobian, LM) at
% the top thrust rung (8x = 200 mN, a few-revolution transfer where the
% shooting basin should be wide), tf = 1.15 x tfMin(200 mN), seeded from the
% P0b min-time costates rescaled onto the energy throttle scale
% (S_e = Tmax(||lam_v||/m + lam_m/c); u* = sat(S_e)). A small grid of seed
% scales beta is tried until one converges -- this is a probe, not a build.
%
% If it converges: Z3's top rung is done before ZTL exists, and the converged
% arc is the ground-truth input for the Z0/Z1 unit tests.
%
% Requires: results/p0b_mintime_backbone.mat (run p0b_mintime_ladder first).
% Output:   results/p0c_energy_top.mat

here = fileparts(mfilename('fullpath'));
cd(here);  setup_paths();
resDir = fullfile(here, 'results');

B = load(fullfile(resDir, 'p0b_mintime_backbone.mat'));
[rv0, rvf, P] = ztl_endpoints();

fTop  = 8;                                  % 200 mN rung
kTop  = find(abs(B.factors - fTop) < 1e-9, 1);
Tmax  = fTop*P.Tmax25;
ctf   = 1.15;
tf    = ctf*B.tfMin(kTop);
lamMT = B.costMT(:, kTop);

% --- seed-scale grid: normalize the energy switching scale S_e(0) ----------
% beta rescales all 7 costates; target S_e(0) in the interior/near-saturated
% range. Also include the old fuel-style bootstrap rescale for reference.
lamvMag0 = norm(lamMT(4:6));
Se0raw   = Tmax*(lamvMag0/1 + lamMT(7)/P.c);      % m0 = 1 (mass fraction)
betaGrid = [0.9/Se0raw, 0.6/Se0raw, 1.5/Se0raw, 1/(lamvMag0*P.c)];

fprintf('=== P0c: min-energy @ 200 mN, tf = %.4f (%.2f x %.4f) ===\n', ...
        tf, ctf, B.tfMin(kTop));
fprintf('raw S_e(0) from min-time costates: %.3e\n', Se0raw);

optsInt = odeset('RelTol', 1e-12, 'AbsTol', 1e-14);
best = struct('resNorm', inf);
for kb = 1:numel(betaGrid)
    beta = betaGrid(kb);
    fprintf('--- seed scale beta = %.3e (S_e(0) = %.2f) ---\n', beta, beta*Se0raw);
    [lamSol, rn, flag] = solve_energy_indirect(rv0, 1, rvf, tf, ...
                             beta*lamMT, Tmax, P.c, P.muStar);
    if rn < best.resNorm
        best = struct('resNorm', rn, 'lamSol', lamSol, 'flag', flag, 'beta', beta);
    end
    if rn < 1e-8, break; end                     % converged -- stop probing
end

% --- integrate the best solution: throttle profile + accounting ------------
[~, yI] = ode113(@lt_pmp_eom_energy, [0 tf], [rv0(:); 1; best.lamSol], ...
                 optsInt, Tmax, P.c, P.muStar);
lamvMag = sqrt(sum(yI(:,11:13).^2, 2));
Se = Tmax*(lamvMag./yI(:,7) + yI(:,14)/P.c);
u  = min(max(Se, 0), 1);
mF = yI(end,7);

sol = struct('factor_thrust', fTop, 'Tmax_mN', 25*fTop, 'ctf', ctf, ...
    'tf', tf, 'tfMin', B.tfMin(kTop), 'lam0', best.lamSol, ...
    'resNorm', best.resNorm, 'flag', best.flag, 'beta', best.beta, ...
    'mProp_kg', P.m0kg*(1-mF), 'dV_kms', P.c*log(1/mF)*P.lStar/P.tStar, ...
    'uMin', min(u), 'uMax', max(u), 'fracSatHi', mean(u > 0.999), ...
    'fracSatLo', mean(u < 1e-3), 'rv0', rv0, 'rvf', rvf, 'P', P);
save(fullfile(resDir, 'p0c_energy_top.mat'), 'sol');
fprintf('saved %s\n', fullfile(resDir, 'p0c_energy_top.mat'));

fprintf(['P0c RESULT: ||R|| = %.3e (flag %d, beta %.3e)  prop = %.4f kg  ' ...
         'dV = %.4f km/s\n  throttle: min %.3f  max %.3f  sat-hi %.1f%%  sat-lo %.1f%%\n'], ...
    best.resNorm, best.flag, best.beta, sol.mProp_kg, sol.dV_kms, ...
    sol.uMin, sol.uMax, 100*sol.fracSatHi, 100*sol.fracSatLo);
if best.resNorm < 1e-8
    fprintf('GATE P0c: PASS -- Z3 top anchor converged with old machinery.\n');
else
    fprintf('GATE P0c: NOT CONVERGED -- Z3 top rung needs the ZTL solver (expected path).\n');
end
