function outFile = gen_elfo_mintime(opts)
% GEN_ELFO_MINTIME  Solve the GTO->ELFO minimum-time (hard all-burn, s==1)
% transfer to anchor the front. Loads the lowest converged energy rung, overrides
% throttle to s==1, minimizes t(tau_f) via casadi_mintime_freetf, saves the
% all-burn min-time solution (with an s==1 throttle row for downstream compat)
% and prints the acceptance diagnostics. Independent solver-free verification is
% a separate step (verify_elfo_seed on the saved file).
%
% INPUTS:
%   opts - (optional) struct: .seedFile[results/energy_elfo_f0990.mat]
%          .maxIter[3000] .warmTight[false]
%
% OUTPUTS:
%   outFile - results/mintime_elfo.mat: X[9xnN], U[4xnN] (alpha + s==1 row),
%             sigma, rv0, rvf, tauf0, tf(=tfMin), mf, cScale, maxDefect, minR1,
%             pSund, qSund, moonZone
%
% REFERENCES:
%   [1] casadi_mintime_freetf.m; [2] the Route B design spec
%       (docs/superpowers/specs/2026-07-15-elfo-mintime-route-b-design.md).

if nargin < 1, opts = struct(); end
gd = @(f,d) getdef(opts,f,d);
here = fileparts(mfilename('fullpath'));  cd(here);  setup_paths();
resDir = fullfile(here,'results');
cfg = minfuel_config();  p = cr3bp_lt_params(cfg.thrustN, cfg.m0kg, cfg.ispS);

seedFile = gd('seedFile', fullfile(resDir,'energy_elfo_f0990.mat'));
S = load(seedFile);
fprintf('=== GEN_ELFO_MINTIME: seed %s (tf=%.4f ND, %.2f d, edge from energy) ===\n', ...
        seedFile, S.X(8,end), S.X(8,end)*p.tStar/86400);

% warm start: states from the seed, throttle overridden to s==1 (drop U row 4)
X0 = S.X;  U0 = S.U(1:3,:);
o = struct('pSund',S.pSund,'qSund',S.qSund,'moonZone',S.moonZone, ...
           'cBox',[0.10 8],'tfCapMult',4,'maxIter',gd('maxIter',3000), ...
           'warmTight',gd('warmTight',false));
out = casadi_mintime_freetf(S.sigma, S.rv0, S.rvf, p.Tmax, p.c, p.muStar, ...
                            X0, U0, S.tauf0, o);

% acceptance diagnostics
rperi = norm(S.rv0(1:3) - [-p.muStar 0 0]);   % GTO perigee radius (ND)
rferr = norm(out.X(1:6,end) - S.rvf(:));
fprintf('  ipopt: %s   success=%d\n', out.ipoptStatus, out.success);
fprintf('  tfMin = %.4f ND (%.2f d)   mf=%.4f (prop %.1f%%)   cScale=%.4f\n', ...
        out.tf, out.tf*p.tStar/86400, out.mf, 100*(1-out.mf), out.cScale);
fprintf('  maxDefect=%.2e  maxUnit=%.2e  rendezvous=%.2e\n', out.maxDefect, out.maxUnit, rferr);
fprintf('  minR1=%.4f (GTO perigee=%.4f)  tMonotone=%d  primerAlign=%.3f deg\n', ...
        out.minR1, rperi, out.tMonotone, out.primerAlignDeg);

% save with a 4-row U (s==1 row) for drop-in compat with verify_elfo_seed / movie
X = out.X;  U = [out.U; ones(1,size(out.U,2))];  sigma = S.sigma; %#ok<NASGU>
rv0 = S.rv0;  rvf = S.rvf;  tauf0 = S.tauf0; %#ok<NASGU>
pSund = S.pSund;  qSund = S.qSund;  moonZone = S.moonZone; %#ok<NASGU>
tf = out.tf;  mf = out.mf;  cScale = out.cScale; %#ok<NASGU>
maxDefect = out.maxDefect;  minR1 = out.minR1; %#ok<NASGU>
outFile = fullfile(resDir,'mintime_elfo.mat');
save(outFile,'X','U','sigma','rv0','rvf','tauf0','tf','mf','cScale', ...
     'maxDefect','minR1','pSund','qSund','moonZone');
fprintf('  saved %s\n', outFile);

% relabel: the mapped front in ELFO's own units
fprintf('\n  RELABEL: factor_ELFO = tf / tfMin_ELFO,  tfMin_ELFO = %.4f ND (%.2f d)\n', ...
        out.tf, out.tf*p.tStar/86400);
fprintf('GEN_ELFO_MINTIME DONE\n');
end

% ---------------------------------------------------------------------------
function v = getdef(s, f, d)
if isfield(s, f) && ~isempty(s.(f)), v = s.(f); else, v = d; end
end
