function outFile = gen_elfo_energy_tfsweep(opts)
% GEN_ELFO_ENERGY_TFSWEEP  Map the GTO->ELFO min-ENERGY tf band by tf-continuation.
%
% Phase 1 of the min-fuel tf-grid campaign (see [[minfuel-tf-grid-strategy]]):
% starting from the one converged ELFO energy seed (energy_elfo_freetf.mat,
% tf=7.5488 ND), step the PINNED transfer time tfTarget UP and DOWN in small
% moves -- each a single casadi_energy_freetf energy solve (eps=1, two-primary
% clock, full gravity, ELFO target), warm-started from the previous tf. cScale
% floats to hold each new tf. This finds the ELFO ENERGY band [tf_lo, tf_hi]
% (where energy stops converging) and banks a seed at each grid tf, so Phase 2
% (the eps->0 fuel sweep per tf) has ready warm starts.
%
% INPUTS:
%   opts - (optional): .tfStep[0.5 ND] .tfHi[12.5] .tfLo[7.0] .maxIter[2000]
%          .looseIter[500] .stepMin[0.0625] .resume[true]
%
% OUTPUTS:
%   outFile - results/energy_elfo_tfgrid.mat: struct array .grid(tf, ok, mf, edge,
%             switches, file) and the band [tfLo tfHi]. Per-tf seeds saved as
%             results/energy_elfo_tf<NNNN>.mat (NNNN = round(1000*tf)).
%
% REFERENCES:
%   [1] casadi_energy_freetf.m; [2] gen_elfo_energy_gravhom.m (the base seed);
%   [3] minfuel-tf-grid-strategy (energy band wider than the fuel-convergent band).

if nargin < 1, opts = struct(); end
gd = @(f,d) getdef(opts,f,d);
here = fileparts(mfilename('fullpath'));  cd(here);  setup_paths();
addpath(fullfile(here,'..','PSR'));
resDir = fullfile(here,'results');
cfg = minfuel_config();  p = cr3bp_lt_params(cfg.thrustN, cfg.m0kg, cfg.ispS);

S = load(fullfile(resDir,'energy_elfo_freetf.mat'));
tf0 = S.X(8,end);
ctx = struct('sigma',S.sigma,'rv0',S.rv0,'rvf',S.rvf,'Tmax',p.Tmax,'cEx',p.c, ...
    'muStar',p.muStar,'tauf0',S.tauf0,'pSund',S.pSund,'qSund',S.qSund, ...
    'moonZone',S.moonZone,'maxIter',gd('maxIter',2000),'looseIter',gd('looseIter',500), ...
    'resDir',resDir,'tStar',p.tStar);
tfStep = gd('tfStep',0.5);  tfHi = gd('tfHi',12.5);  tfLo = gd('tfLo',7.0);
stepMin = gd('stepMin',0.0625);

fprintf('=== GEN_ELFO_ENERGY_TFSWEEP: tf band map from tf0=%.4f ND (%.2f d) ===\n', ...
        tf0, tf0*p.tStar/86400);

% seed grid-point at tf0 itself
grid = save_point(ctx, S.X, S.U, tf0, true);   % the base seed (already converged)

% --- sweep UP -----------------------------------------------------------------
grid = [grid, sweep_dir(ctx, S.X, S.U, tf0, +tfStep, tfHi, stepMin)];
% --- sweep DOWN ---------------------------------------------------------------
grid = [grid, sweep_dir(ctx, S.X, S.U, tf0, -tfStep, tfLo, stepMin)];

% --- summary ------------------------------------------------------------------
[~,ord] = sort([grid.tf]);  grid = grid(ord);
okv = [grid.ok];  tfs = [grid.tf];
tfLoB = min(tfs(okv));  tfHiB = max(tfs(okv));
outFile = fullfile(resDir,'energy_elfo_tfgrid.mat');
save(outFile,'grid','tfLoB','tfHiB');
fprintf('\n--- ELFO ENERGY tf band = [%.4f, %.4f] ND (%.2f - %.2f d) ---\n', ...
        tfLoB, tfHiB, tfLoB*p.tStar/86400, tfHiB*p.tStar/86400);
fprintf('  tf(ND)   tf(d)   ok  mf      edge   sw\n');
for k = 1:numel(grid)
    fprintf('  %6.4f  %5.2f   %d   %.4f  %5.1f%%  %d\n', grid(k).tf, ...
        grid(k).tf*p.tStar/86400, grid(k).ok, grid(k).mf, 100*grid(k).edge, grid(k).switches);
end
fprintf('GEN_ELFO_ENERGY_TFSWEEP DONE: %s\n', outFile);
end

% ===========================================================================
function G = sweep_dir(ctx, X, U, tf0, dstep, tfLimit, stepMin)
% Continuation from tf0 by dstep (signed) until tfLimit or a stuck step.
G = struct('tf',{},'ok',{},'mf',{},'edge',{},'switches',{},'file',{});
Xk = X;  Uk = U;  tf = tf0;  step = abs(dstep);  sgn = sign(dstep);
while (sgn > 0 && tf < tfLimit-1e-9) || (sgn < 0 && tf > tfLimit+1e-9)
    tfTry = tf + sgn*step;
    if sgn > 0, tfTry = min(tfTry, tfLimit); else, tfTry = max(tfTry, tfLimit); end
    [ok, Xn, Un, info] = solve_tf(ctx, tfTry, Xk, Uk);
    if ~ok
        step = step/2;
        fprintf('  tf=%.4f FAIL (def=%.2g) -> step=%.4f\n', tfTry, info.maxDefect, step);
        if step < stepMin
            fprintf('  ** energy band edge near tf=%.4f (dir %+d) **\n', tf, sgn);
            break
        end
        continue
    end
    Xk = Xn;  Uk = Un;  tf = tfTry;
    G(end+1) = save_point(ctx, Xk, Uk, tf, true); %#ok<AGROW>
    fprintf('  tf=%.4f (%.2f d) OK def=%.2g edge=%.1f%% sw=%d mf=%.4f\n', ...
        tf, tf*ctx.tStar/86400, info.maxDefect, 100*info.edge, info.switches, info.mf);
end
end

% ===========================================================================
function [ok, Xn, Un, info] = solve_tf(ctx, tfTarget, Xk, Uk)
base = struct('moonZone',ctx.moonZone,'muGain',1,'tfTarget',tfTarget,'epsilon',1, ...
              'pSund',ctx.pSund,'qSund',ctx.qSund,'tfCapMult',8,'cBox',[0.10 8]);
oL = base;  oL.maxIter = ctx.looseIter;  oL.warmTight = false;
rL = casadi_energy_freetf(ctx.sigma,ctx.rv0,ctx.rvf,ctx.Tmax,ctx.cEx,ctx.muStar,Xk,Uk,ctx.tauf0,oL);
if rL.success && rL.maxDefect < 1e-6
    Xn = rL.X;  Un = rL.U;  ok = true;  info = rL;  return
end
oF = base;  oF.maxIter = ctx.maxIter;  oF.warmTight = true;
rF = casadi_energy_freetf(ctx.sigma,ctx.rv0,ctx.rvf,ctx.Tmax,ctx.cEx,ctx.muStar,Xk,Uk,ctx.tauf0,oF);
if rF.success && rF.maxDefect < 1e-6
    Xn = rF.X;  Un = rF.U;  ok = true;  info = rF;
else
    ok = false;  Xn = Xk;  Un = Uk;  info = rF;
end
end

% ===========================================================================
function g = save_point(ctx, X, U, tf, ok)
sigma = ctx.sigma;  rv0 = ctx.rv0;  rvf = ctx.rvf;  tauf0 = ctx.tauf0; %#ok<NASGU>
moonZone = ctx.moonZone;  pSund = ctx.pSund;  qSund = ctx.qSund; %#ok<NASGU>
file = fullfile(ctx.resDir, sprintf('energy_elfo_tf%04d.mat', round(1000*tf)));
save(file,'X','U','sigma','rv0','rvf','tauf0','tf','moonZone','pSund','qSund');
ss = U(4,:);
g = struct('tf',tf,'ok',ok,'mf',X(7,end),'edge',mean(ss>0.95|ss<0.05), ...
           'switches',sum(abs(diff(ss>0.5))),'file',file);
end

% ===========================================================================
function v = getdef(s, f, d)
if isfield(s, f) && ~isempty(s.(f)), v = s.(f); else, v = d; end
end
