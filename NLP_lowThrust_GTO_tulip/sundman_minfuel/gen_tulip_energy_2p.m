function outFile = gen_tulip_energy_2p(opts)
% GEN_TULIP_ENERGY_2P  Manufacture a two-primary-clock GTO->tulip min-ENERGY seed,
% re-meshed for the two-primary Sundman clock, so the min-TIME solve
% (gen_tulip_mintime) can run the EXACT Hessian without the near-Moon overflow
% that crashes MUMPS on the single-primary tulip backbone.
%
% WHY: the tulip energy backbones (sundman_minfuel/results/energy/energy_f####.mat)
% use the SINGLE-primary clock kappa=r1^pSund, which leaves lunar gravity untamed
% near the tulip terminal (dMoon ~28k km).  casadi_mintime_freetf's exact Hessian
% then overflows there -> IPOPT/MUMPS bus error (root-caused 2026-07-15).  The
% two-primary clock kappa=(r1^-q+(r2/D)^-q)^(-p/q) tames the near-Moon Hessian and
% concentrates nodes into the lunar arc.  This driver turns that clock ON at the
% tulip target, re-meshing the trajectory so the min-time solve is mesh-consistent.
%
% It mirrors elfo/gen_elfo_energy_gravhom MINUS the ELFO retarget leg (we stay at
% the tulip target throughout).  One change per leg (1-D continuation):
%   LEG 0  convert the fixed-t_f backbone to the free-t_f representation
%          (mu=1, single-primary clock, tulip target).
%   LEG A  gravity OFF: muGain 1 -> 0   (single-primary, tulip) -- cleanest leg.
%   LEG B  clock ON: moonZone 0 -> mz   (mu=0, tulip) -- benign re-mesh with the
%          well off (clock-on at mu=1 stiffens, per the ELFO record).
%   LEG D  gravity ON: muGain 0 -> 1    (moonZone=mz, tulip) -- the well returns;
%          the two-primary clock keeps the near-Moon nodes resolved.
%
% Adaptive step (halve on fail) + per-step checkpoint/resume.
%
% INPUTS:
%   opts - (optional) struct:
%          .factor   tulip energy backbone to start from   [1.12 = lowest]
%          .moonZone two-primary crossover D (ND)          [0.15 ~ lunar SOI]
%          .qSund    two-primary transition sharpness       [4]
%          .step0    initial homotopy step per leg          [0.20]
%          .stepMin  give up below this step                [0.01]
%          .maxIter  IPOPT cap (tight solves)               [1500]
%          .looseIter IPOPT cap (loose probe)               [400]
%          .resume   pick up from checkpoint if present     [true]
%
% OUTPUTS:
%   outFile - results/energy_tulip_2p.mat: X[9x(N+1)],U,factor,tauf0,sigma,rv0,
%             rvf(=tulip),pSund,qSund,moonZone,muGain -- the two-primary tulip
%             energy seed / warm start for gen_tulip_mintime (moonZone=mz).
%
% REFERENCES:
%   [1] casadi_energy_freetf.m (the free-t_f two-primary solver this drives).
%   [2] elfo/gen_elfo_energy_gravhom.m (the ELFO ladder this adapts).

if nargin < 1, opts = struct(); end
gd = @(f,d) getdef(opts,f,d);
here = fileparts(mfilename('fullpath'));  cd(here);  setup_paths();
addpath(fullfile(here,'..','elfo'));        % casadi_energy_freetf
resDir = fullfile(here,'results');  if ~exist(resDir,'dir'), mkdir(resDir); end
ckptFile = fullfile(resDir,'energy_tulip_2p_ckpt.mat');

cfg = minfuel_config();  p = cr3bp_lt_params(cfg.thrustN, cfg.m0kg, cfg.ispS);
% Start from the LEAST-SATURATED backbone (lowest edge), NOT the lowest tf: the
% gravity-off/clock-on continuation deforms an interior throttle cleanly but
% CANNOT restructure a saturated (high-edge) one (the ELFO lesson). Tulip edge is
% non-monotone in tf, minimum at factor 1.20 (edge 12.4%); factor 1.12 is edge
% 71% and stalls LEG A. (The starting factor sets seed-gen conditioning only --
% the min-time answer is warm-start-independent.)
factor = gd('factor', 1.20);
E = load(fullfile(resDir,'energy',sprintf('energy_f%04d.mat', round(1000*factor))));
sigma = E.sigma;  rv0 = E.rv0;  rvf_tul = E.rvf;  tauf0 = E.tauf0;
mz = gd('moonZone', 0.15);
fprintf('=== GEN_TULIP_ENERGY_2P: tulip(f=%.2f) clock-on to two-primary ===\n', factor);
fprintf('  N=%d  tf_ws=%.4f ND  target dMoon=%.0f km\n', numel(sigma)-1, E.X(8,end), ...
        norm(rvf_tul(1:3)-[1-p.muStar 0 0])*p.lStar);

ctx = struct('sigma',sigma,'rv0',rv0,'Tmax',p.Tmax,'cEx',p.c,'muStar',p.muStar, ...
    'tauf0',tauf0,'pSund',cfg.pSund,'qSund',gd('qSund',4), ...
    'maxIter',gd('maxIter',1500),'looseIter',gd('looseIter',400), ...
    'step0',gd('step0',0.20),'stepMin',gd('stepMin',0.01), ...
    'tf0',E.X(8,end),'factor',factor,'rvf',rvf_tul);

% one-change-per-leg parameterizations (target FIXED at rvf_tul throughout)
legs = { ...
  'A_gravity_off', @(s) struct('moonZone',0,    'muGain',1-s); ...
  'B_clock_on',    @(s) struct('moonZone',mz*s, 'muGain',0); ...
  'D_gravity_on',  @(s) struct('moonZone',mz,   'muGain',s) };

% --- resume -----------------------------------------------------------------
startLeg = 0;  sStart = 0;  Xk = E.X;  Uk = E.U;
if gd('resume',true) && isfile(ckptFile)
    C = load(ckptFile);
    if isequal(C.ctx.factor,factor) && norm(C.ctx.rvf(:)-rvf_tul(:)) < 1e-9
        startLeg = C.legIdx;  sStart = C.s;  Xk = C.Xk;  Uk = C.Uk;
        fprintf('  RESUMED at leg %d, s=%.3f\n', startLeg, sStart);
    end
end

% --- LEG 0: fixed-t_f backbone -> free-t_f (mu=1, single-primary, tulip) -----
if startLeg == 0
    fprintf('--- LEG 0: convert to free-t_f (mu=1, single-primary, tulip) ---\n');
    [ok,Xk,Uk,info] = step_solve(ctx, struct('moonZone',0,'muGain',1), Xk, Uk);
    if ~ok, error('tulip2p:leg0','free-t_f conversion failed (def=%.2g)', info.maxDefect); end
    fprintf('  LEG0 OK def=%.2g tf=%.3f cS=%.3f edge=%.1f%%\n', ...
            info.maxDefect, info.tf, info.cScale, 100*info.edge);
    legIdx = 1;  s = 0;  save(ckptFile,'legIdx','s','Xk','Uk','ctx');
    startLeg = 1;  sStart = 0;
end

% --- LEGS A, B, D -----------------------------------------------------------
for L = max(startLeg,1):3
    s0 = 0;  if L == startLeg, s0 = sStart; end
    [Xk,Uk] = walk(ctx, L, legs{L,1}, legs{L,2}, Xk, Uk, s0, ckptFile);
end

% --- save the two-primary tulip energy seed ---------------------------------
X = Xk;  U = Uk;  rvf = rvf_tul(:).';  pSund = cfg.pSund; %#ok<NASGU>
qSund = ctx.qSund;  moonZone = mz;  muGain = 1; %#ok<NASGU>
outFile = fullfile(resDir, 'energy_tulip_2p.mat');
save(outFile, 'X','U','factor','tauf0','sigma','rv0','rvf','pSund','qSund','moonZone','muGain');
fprintf('GEN_TULIP_ENERGY_2P DONE: %s\n', outFile);
fprintf('  tf=%.4f ND (%.2f d), mf=%.4f, edge=%.1f%%, maxDefect (last)=see log\n', ...
        X(8,end), X(8,end)*p.tStar/86400, X(7,end), 100*mean(U(4,:)>0.95|U(4,:)<0.05));
fprintf('  NEXT: gen_tulip_mintime(struct(''seedFile'',''%s'',''moonZone'',%.3f))\n', outFile, mz);
end

% ===========================================================================
function [Xk, Uk] = walk(ctx, legIdx, legName, paramFun, Xk, Uk, s0, ckptFile)
s = s0;  step = ctx.step0;  nstep = 0;
fprintf('--- LEG %d (%s): s=%.3f -> 1 ---\n', legIdx, legName, s0);
while s < 1 - 1e-9
    sTry = min(s + step, 1);
    [ok, Xn, Un, info] = step_solve(ctx, paramFun(sTry), Xk, Uk);
    if ~ok
        step = step/2;
        fprintf('  %s s=%.4f FAIL (def=%.2g) -> step=%.4f\n', legName, sTry, info.maxDefect, step);
        if step < ctx.stepMin, error('tulip2p:stuck','stuck in leg %s at s=%.4f', legName, s); end
        continue
    end
    Xk = Xn;  Uk = Un;  s = sTry;  nstep = nstep + 1;
    save(ckptFile,'legIdx','s','Xk','Uk','ctx');
    fprintf('  %s s=%.4f OK def=%.2g tf=%.3f cS=%.3f edge=%.1f%% (step %d)\n', ...
            legName, s, info.maxDefect, info.tf, info.cScale, 100*info.edge, nstep);
    if step < ctx.step0, step = min(2*step, ctx.step0); end
end
end

% ===========================================================================
function [ok, Xn, Un, info] = step_solve(ctx, oExtra, Xk, Uk)
% One continuation step: loose probe -> tight fallback -> tight re-clean.
base = struct('pSund',ctx.pSund,'qSund',ctx.qSund,'tfCapMult',6,'cBox',[0.15 6], ...
              'tfTarget',ctx.tf0);              % PIN t_f (well-posed energy)
o = setfields(base, oExtra);
oL = o;  oL.maxIter = ctx.looseIter;  oL.warmTight = false;
rL = casadi_energy_freetf(ctx.sigma,ctx.rv0,ctx.rvf,ctx.Tmax,ctx.cEx,ctx.muStar,Xk,Uk,ctx.tauf0,oL);
if rL.success && rL.maxDefect < 1e-6
    Xs = rL.X;  Us = rL.U;
else
    oF = o;  oF.maxIter = ctx.maxIter;  oF.warmTight = true;
    rF = casadi_energy_freetf(ctx.sigma,ctx.rv0,ctx.rvf,ctx.Tmax,ctx.cEx,ctx.muStar,Xk,Uk,ctx.tauf0,oF);
    if rF.success && rF.maxDefect < 1e-6
        Xs = rF.X;  Us = rF.U;
    else
        ok = false;  Xn = Xk;  Un = Uk;  info = rF;  return
    end
end
oT = o;  oT.maxIter = ctx.maxIter;  oT.warmTight = true;
rT = casadi_energy_freetf(ctx.sigma,ctx.rv0,ctx.rvf,ctx.Tmax,ctx.cEx,ctx.muStar,Xs,Us,ctx.tauf0,oT);
if rT.success && rT.maxDefect < 1e-6
    Xn = rT.X;  Un = rT.U;  ok = true;  info = rT;
else
    ok = false;  Xn = Xk;  Un = Uk;  info = rT;
end
end

% ===========================================================================
function s = setfields(s, o)
if isempty(o) || ~isstruct(o), return; end
f = fieldnames(o);
for k = 1:numel(f), s.(f{k}) = o.(f{k}); end
end

% ===========================================================================
function v = getdef(s, f, d)
if isfield(s, f) && ~isempty(s.(f)), v = s.(f); else, v = d; end
end
