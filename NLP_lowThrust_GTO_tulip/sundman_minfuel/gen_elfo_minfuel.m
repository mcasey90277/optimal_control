function outFile = gen_elfo_minfuel(opts)
% GEN_ELFO_MINFUEL  Bertrand-Epenoy energy->fuel homotopy for the GTO->ELFO
% transfer, on the free-t_f two-primary solver, from the energy seed.
%
% Consumes energy_elfo_freetf.mat (the converged GTO->ELFO min-ENERGY solution
% from gen_elfo_energy_gravhom.m) and sweeps epsilon 1 -> 0:
%   J(eps) = Int[s]dt - eps*Int[s(1-s)]dt   (physical-time measure)
%   eps=1 -> Int[s^2]dt   ENERGY (the seed, smooth throttle ramp)
%   eps=0 -> Int[s]dt     FUEL   (linear in s -> bang-bang, = propellant up to a
%                                 positive constant since mf=1-(Tmax/c)Int[s]dt)
% Everything else is held at the seed's converged configuration: t_f PINNED at the
% seed value (well-posed), two-primary clock (moonZone), full gravity (muGain=1),
% ELFO target. cScale floats to hold t_f under the sharpening throttle.
%
% As eps->0 the throttle sharpens toward bang-bang (edge/switches grow); the
% loose bound-push then shoves s off its bounds (false "infeasible"), so each step
% is loose-probe -> tight-fallback -> tight-reclean (the tight settings hug the
% near-bang-bang bounds). Adaptive step (halve on fail) + checkpoint/resume.
%
% INPUTS:
%   opts - (optional): .step0[0.20] .stepMin[0.005] .maxIter[2000]
%          .looseIter[500] .resume[true]
%
% OUTPUTS:
%   outFile - results/minfuel_elfo.mat: X[9x(N+1)],U,sigma,rv0,rvf,tauf0,tf,
%             moonZone,pSund,qSund,epsilon(=0) -- the GTO->ELFO min-fuel solution.
%
% REFERENCES:
%   [1] Bertrand & Epenoy (2002); [2] gen_elfo_energy_gravhom.m (the seed source);
%   [3] casadi_energy_freetf.m (the solver).

if nargin < 1, opts = struct(); end
gd = @(f,d) getdef(opts,f,d);

here = fileparts(mfilename('fullpath'));  cd(here);  setup_paths();
addpath(fullfile(here,'..','PSR'));
resDir = fullfile(here,'results');  ckptFile = fullfile(resDir,'minfuel_elfo_ckpt.mat');
cfg = minfuel_config();  p = cr3bp_lt_params(cfg.thrustN, cfg.m0kg, cfg.ispS);

seedFile = fullfile(resDir,'energy_elfo_freetf.mat');
assert(isfile(seedFile), 'no ELFO energy seed at %s (run gen_elfo_energy_gravhom first)', seedFile);
S = load(seedFile);
tf0 = S.X(8,end);

ctx = struct('sigma',S.sigma,'rv0',S.rv0,'rvf',S.rvf,'Tmax',p.Tmax,'cEx',p.c, ...
    'muStar',p.muStar,'tauf0',S.tauf0,'pSund',S.pSund,'qSund',S.qSund, ...
    'moonZone',S.moonZone,'tf0',tf0, ...
    'maxIter',gd('maxIter',2000),'looseIter',gd('looseIter',500), ...
    'step0',gd('step0',0.20),'stepMin',gd('stepMin',0.005));

fprintf('=== GEN_ELFO_MINFUEL: energy->fuel (eps 1->0) at tf=%.4f ND (%.2f d) ===\n', ...
        tf0, tf0*p.tStar/86400);

% resume (s in [0,1] maps to epsilon = 1 - s)
Xk = S.X;  Uk = S.U;  s0 = 0;
if gd('resume',true) && isfile(ckptFile)
    C = load(ckptFile);
    if abs(C.tf0 - tf0) < 1e-9 && C.s < 1-1e-9
        Xk = C.Xk;  Uk = C.Uk;  s0 = C.s;
        fprintf('  RESUMED at eps=%.4f\n', 1-s0);
    end
end

% --- epsilon homotopy 1 -> 0 ------------------------------------------------
s = s0;  step = ctx.step0;  nstep = 0;
while s < 1 - 1e-9
    sTry = min(s + step, 1);
    epsilon = 1 - sTry;
    [ok, Xn, Un, info] = step_solve(ctx, epsilon, Xk, Uk);
    if ~ok
        step = step/2;
        fprintf('  eps=%.4f FAIL (def=%.2g) -> step=%.4f\n', epsilon, info.maxDefect, step);
        if step < ctx.stepMin, error('minfuel:stuck','stuck at eps=%.4f (min-fuel sharpening wall)', epsilon); end
        continue
    end
    Xk = Xn;  Uk = Un;  s = sTry;  nstep = nstep + 1;
    save(ckptFile,'s','Xk','Uk','tf0');
    fprintf('  eps=%.4f OK def=%.2g sw=%d edge=%.1f%% mf=%.4f cS=%.3f (step %d)\n', ...
            epsilon, info.maxDefect, info.switches, 100*info.edge, info.mf, info.cScale, nstep);
    if step < ctx.step0, step = min(2*step, ctx.step0); end
end

% --- save the ELFO min-fuel solution ----------------------------------------
X = Xk;  U = Uk;  rvf = ctx.rvf;  sigma = ctx.sigma;  rv0 = ctx.rv0;  tauf0 = ctx.tauf0; %#ok<NASGU>
tf = X(8,end);  moonZone = ctx.moonZone;  pSund = ctx.pSund;  qSund = ctx.qSund;  epsilon = 0; %#ok<NASGU>
outFile = fullfile(resDir, 'minfuel_elfo.mat');
save(outFile, 'X','U','sigma','rv0','rvf','tauf0','tf','moonZone','pSund','qSund','epsilon');
ss = U(4,:);  nSw = sum(abs(diff(ss>0.5)));
fprintf('GEN_ELFO_MINFUEL DONE: %s\n', outFile);
fprintf('  eps=0 fuel: switches=%d  edge=%.1f%%  mf=%.4f (prop %.1f%%)  tf=%.2f d\n', ...
        nSw, 100*mean(ss>0.95|ss<0.05), X(7,end), 100*(1-X(7,end)), tf*p.tStar/86400);
end

% ===========================================================================
function [ok, Xn, Un, info] = step_solve(ctx, epsilon, Xk, Uk)
% One epsilon step: loose probe -> tight fallback -> tight re-clean.
base = struct('moonZone',ctx.moonZone,'muGain',1,'tfTarget',ctx.tf0,'epsilon',epsilon, ...
              'pSund',ctx.pSund,'qSund',ctx.qSund,'tfCapMult',6,'cBox',[0.15 6]);
oL = base;  oL.maxIter = ctx.looseIter;  oL.warmTight = false;
rL = casadi_energy_freetf(ctx.sigma,ctx.rv0,ctx.rvf,ctx.Tmax,ctx.cEx,ctx.muStar,Xk,Uk,ctx.tauf0,oL);
if rL.success && rL.maxDefect < 1e-6
    Xs = rL.X;  Us = rL.U;
else
    oF = base;  oF.maxIter = ctx.maxIter;  oF.warmTight = true;
    rF = casadi_energy_freetf(ctx.sigma,ctx.rv0,ctx.rvf,ctx.Tmax,ctx.cEx,ctx.muStar,Xk,Uk,ctx.tauf0,oF);
    if rF.success && rF.maxDefect < 1e-6
        Xs = rF.X;  Us = rF.U;
    else
        ok = false;  Xn = Xk;  Un = Uk;  info = rF;  return
    end
end
oT = base;  oT.maxIter = ctx.maxIter;  oT.warmTight = true;
rT = casadi_energy_freetf(ctx.sigma,ctx.rv0,ctx.rvf,ctx.Tmax,ctx.cEx,ctx.muStar,Xs,Us,ctx.tauf0,oT);
if rT.success && rT.maxDefect < 1e-6
    Xn = rT.X;  Un = rT.U;  ok = true;  info = rT;
else
    ok = false;  Xn = Xk;  Un = Uk;  info = rT;
end
end

% ===========================================================================
function v = getdef(s, f, d)
if isfield(s, f) && ~isempty(s.(f)), v = s.(f); else, v = d; end
end
