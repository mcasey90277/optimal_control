% ELFO_MINTIME  Reach a min-TIME GTO->ELFO rendezvous by homotoping the target
% from the tulip min-time solution (the pumpkyn costate root).
%
% Why min-time first (the tulip's own lineage): min-time is ALWAYS-BURN -- no
% throttle, so no saturation/edge sensitivity and none of the fixed-tf
% "can't-reach-terminal" stalls that the energy-target homotopy hits as the
% ELFO terminal moves into the Moon's gravity well. Its tf FLOATS, so every
% intermediate target on the homotopy is reachable (just at a different tf).
% Once min-time reaches the ELFO, its trajectory seeds a min-ENERGY GTO->ELFO
% solve (-> PSR fuel), exactly as min-time seeded energy for the tulip.
%
% Solver: mintime_solve (our real-budget rebuild of pumpkyn.cr3bp.tfMin).

here = fileparts(mfilename('fullpath'));  cd(here);  setup_paths();
[rv0, rvf_tul, P] = mintime_params();
Tmax = P.Tmax25;  c = P.c;  mu = P.muStar;  lStar = P.lStar;

% pumpkyn's published converged tulip min-time costate/tf root (build_guess):
zSeed = [ 190.476497248065; -79.7064866984696; -0.430399154713168; ...
            0.301159446575878; 0.586671892449694; -0.00711582435720301; ...
            4.32931089137559; 6.29081541876621];

% --- confirm the tulip min-time converges from the root ---------------------
[zTul, rnTul, oTul] = mintime_solve(rv0, rvf_tul, zSeed, Tmax, c, mu, 1500);
fprintf('tulip min-time: ||R||=%.3e  nSwitch=%d  tf=%.4f (%.2f d)\n', ...
        rnTul, oTul.nSwitch, zTul(8), zTul(8)*P.tStar/86400);

% --- build the ELFO, pick apolune AND nearest-to-tulip targets --------------
Mmass = 5.9736e24 + 7.35e22;  nu0 = pumpkyn.util.mean2True(0, 0.69);
oev = [12000, 0.69, 56.5*pi/180, 90*pi/180, 0, nu0];
x0  = pumpkyn.cr3bp.fromOrb(0, oev, Mmass, mu, P.tStar, lStar, 2);
tau = linspace(0, 0.5, 8000)';  [~, xs] = pumpkyn.cr3bp.prop(tau, x0, mu);
rM = [1-mu 0 0];  dMoon = sqrt(sum((xs(:,1:3)-rM).^2,2));
[~, iApo] = max(dMoon);  rvApo = xs(iApo,1:6);
d6 = sqrt(sum((xs(:,1:6)-rvf_tul(:).').^2,2)); [~, iNear] = min(d6); rvNear = xs(iNear,1:6);

% --- homotope the min-time target: tulip -> each ELFO point ------------------
names = {'apolune','nearest'};  tgts = {rvApo, rvNear};
for ti = 1:numel(names)
    rvf_elfo = tgts{ti};  nm = names{ti};
    fprintf('\n--- min-time homotopy tulip -> ELFO %s (dMoon=%.0f km, ||drvf||=%.3f) ---\n', ...
            nm, norm(rvf_elfo(1:3)-rM)*lStar, norm(rvf_elfo(:)-rvf_tul(:)));
    zc = zTul;  ok = true;  sPrev = 0;
    for s = 0.05:0.05:1.0
        rvf_s = (1-s)*rvf_tul(:).' + s*rvf_elfo(:).';
        [zc, rn, o] = mintime_solve(rv0, rvf_s, zc, Tmax, c, mu, 1500);
        fprintf('  s=%.2f  ||R||=%.3e  nSwitch=%d  tf=%.4f (%.2f d)\n', ...
                s, rn, o.nSwitch, zc(8), zc(8)*P.tStar/86400);
        if rn > 1e-6
            fprintf('  -> STALLED at s=%.2f (last good s=%.2f, ||R||=%.2e)\n', s, sPrev, rn);
            ok = false;  break
        end
        sPrev = s;
    end
    if ok
        fprintf('  ELFO %s REACHED in min-time: tf=%.4f ND (%.2f d). Saving costate.\n', ...
                nm, zc(8), zc(8)*P.tStar/86400);
        save(fullfile(here,'results',sprintf('mintime_elfo_%s.mat',nm)), ...
             'zc','rvf_elfo','rv0','Tmax','c','mu','P');
    end
end
