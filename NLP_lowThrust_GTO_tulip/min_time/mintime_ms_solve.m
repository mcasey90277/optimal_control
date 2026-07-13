% MINTIME_MS_SOLVE  Converge the tulip min-time via multiple shooting (INDIRECT,
% PMP costate TPBVP), to beat the ~1e-3 single-shooting floor. Seeds from the
% single-shooting solution (continuity exactly 0), then runs the genericized
% trust-region solver (ztl_ms_solve_tr, prob.resFun = mintime_ms_residual).
%
% Base vars: Mnodes [60], maxIter [120].

here = fileparts(mfilename('fullpath'));  cd(here);  setup_paths();
addpath(fullfile(here, '..', 'ztl'));                 % ztl_ms_solve_tr (generic solver)
warning('off','MATLAB:ode45:IntegrationTolNotMet');
resDir = fullfile(here,'results');  if ~exist(resDir,'dir'), mkdir(resDir); end

if evalin('base','exist(''Mnodes'',''var'')'),  M = evalin('base','Mnodes');  else, M = 60;  end
if evalin('base','exist(''maxIter'',''var'')'), maxIter = evalin('base','maxIter'); else, maxIter = 120; end

[rv0, rvf, P] = mintime_params();  Tmax = P.Tmax25;  c = P.c;  mu = P.muStar;
zSeed = [ 190.476497248065; -79.7064866984696; -0.430399154713168; ...
            0.301159446575878; 0.586671892449694; -0.00711582435720301; ...
            4.32931089137559; 6.29081541876621];

[zt, rn, o] = mintime_solve(rv0, rvf, zSeed, Tmax, c, mu, 1500);
fprintf('single-shoot tulip min-time: ||R||=%.3e  nSwitch=%d  tf=%.4f\n', rn, o.nSwitch, zt(8));

[z0, prob, si] = mintime_ms_seed(zt(1:7), zt(8), rv0, rvf, Tmax, c, mu, M);
prob.resFun = @mintime_ms_residual;
fprintf('MS seed M=%d: maxCont=%.2e  termErr=%.2e  (nZ=%d)\n', M, si.maxCont, si.termErr, numel(z0));

[z, out] = ztl_ms_solve_tr(z0, prob, struct('tolR',1e-9,'maxIter',maxIter));
fprintf('MS SOLVE: ||R||=%.4e  termErr=%.2e  maxCont=%.2e  iters=%d  flag=%d\n', ...
        out.resNorm, out.termErr, out.maxCont, out.iters, out.flag);

save(fullfile(resDir,'mintime_tulip_ms.mat'), 'z','out','prob','M','zt');
if out.resNorm < 1e-6
    fprintf('GATE MINTIME-MS: PASS -- beat the single-shooting floor (%.2e -> %.2e). tf=%.5f\n', ...
            rn, out.resNorm, z(end));
else
    fprintf('GATE MINTIME-MS: reached %.2e (flag %d).\n', out.resNorm, out.flag);
end
