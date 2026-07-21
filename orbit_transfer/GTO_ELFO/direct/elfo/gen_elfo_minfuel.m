function outFile = gen_elfo_minfuel(opts)
% GEN_ELFO_MINFUEL  Energy->fuel (Bertrand-Epenoy) homotopy CORE for a GTO->target
% low-thrust transfer, on the free-t_f two-primary solver casadi_energy_freetf.
%
% Called by RUN_ELFO_MINFUEL (the entry driver) and also runnable standalone.
% Consumes a converged min-ENERGY seed (energy_elfo_freetf.mat, or a tf-grid seed
% energy_elfo_f####.mat) and sweeps epsilon 1 -> epsMin:
%   J(eps) = Int[s]dt - eps*Int[s(1-s)]dt   (physical-time measure)
%   eps=1 -> Int[s^2]dt ENERGY (the seed) ;  eps=0 -> Int[s]dt FUEL (bang-bang)
% Everything else is held at the seed's converged configuration: t_f PINNED
% (well-posed), two-primary clock, full gravity, target rvf. cScale floats to
% hold t_f under the sharpening throttle. As eps->0 the loose bound-push shoves
% the throttle off its bounds, so each step is loose-probe -> tight-fallback ->
% tight-reclean. Adaptive step (halve on fail) + checkpoint/resume.
%
% INPUTS:
%   opts - (optional) struct:
%          .seedFile - energy seed .mat [results/energy_elfo_freetf.mat]
%          .target   - name tag baked into the OUTPUT filename ['ELFO']
%          .epsMin   - homotopy endpoint: 0 = bang-bang fuel; >0 = smooth [0]
%          .outFile  - explicit output path; [] -> derive canonical name  []
%          .step0[0.20] .stepMin[0.005] .maxIter[2000] .looseIter[500] .resume[true]
%
% OUTPUTS:
%   outFile - results/minfuel_<target>_tf<fTag>_sw<k>_minEps<eTag>.mat, holding
%             X[9x(N+1)],U,sigma,rv0,rvf,tauf0,tf,moonZone,pSund,qSund,epsilon,
%             target,factor -- the GTO->target min-fuel (or smooth eps>0) solution.
%
% REFERENCES:
%   [1] Bertrand & Epenoy (2002); [2] casadi_energy_freetf.m (the solver);
%   [3] gen_elfo_energy_gravhom.m (the energy seed); [4] run_psr.m (tulip analog).

if nargin < 1, opts = struct(); end
gd = @(f,d) getdef(opts,f,d);

here = fileparts(mfilename('fullpath'));  cd(here);  setup_paths();
resDir = fullfile(here,'results');
cfg = minfuel_config();  p = cr3bp_lt_params(cfg.thrustN, cfg.m0kg, cfg.ispS);

seedPath = gd('seedFile', fullfile(resDir,'energy_elfo_freetf.mat'));
target   = gd('target', 'ELFO');
epsMin   = gd('epsMin', 0);
assert(isfile(seedPath), 'no energy seed at %s', seedPath);
S = load(seedPath);
tf0    = S.X(8,end);
% Factor label vs the target's OWN anchor (2026-07-21 triage C1): ELFO runs
% use the certified ELFO min-time anchor; a 'Tulip'-tagged run keeps the
% tulip anchor. The physical tf0 (from the seed) is authoritative either way.
if strcmpi(target, 'ELFO'), tfMinRef = cfg.tfMin_elfo; else, tfMinRef = cfg.tfMin; end
factor = tf0 / tfMinRef;
fTag   = strrep(sprintf('%.3f', factor), '.', 'p');       % 1.200 -> '1p200'
eTag   = strrep(sprintf('%g', epsMin), '.', 'p');         % 0 -> '0', 0.2 -> '0p2'
tag    = sprintf('%s_tf%s_minEps%s', target, fTag, eTag);
ckptFile = fullfile(resDir, sprintf('minfuel_%s_ckpt.mat', tag));

ctx = struct('sigma',S.sigma,'rv0',S.rv0,'rvf',S.rvf,'Tmax',p.Tmax,'cEx',p.c, ...
    'muStar',p.muStar,'tauf0',S.tauf0,'pSund',S.pSund,'qSund',S.qSund, ...
    'moonZone',S.moonZone,'tf0',tf0, ...
    'maxIter',gd('maxIter',2000),'looseIter',gd('looseIter',500), ...
    'step0',gd('step0',0.20),'stepMin',gd('stepMin',0.005));

fprintf('=== GEN_ELFO_MINFUEL [%s]: eps 1->%g at tf=%.4f ND (%.2f d, %.3fx) ===\n', ...
        target, epsMin, tf0, tf0*p.tStar/86400, factor);

% epsilon homotopy: s in [0,1] maps epsilon = 1 - s*(1-epsMin)
Xk = S.X;  Uk = S.U;  s0 = 0;
if gd('resume',true) && isfile(ckptFile)
    C = load(ckptFile);
    if abs(C.tf0 - tf0) < 1e-9 && C.s < 1-1e-9
        Xk = C.Xk;  Uk = C.Uk;  s0 = C.s;
        fprintf('  RESUMED at eps=%.4f\n', 1 - s0*(1-epsMin));
    end
end

s = s0;  step = ctx.step0;  nstep = 0;  finalInfo = [];
while s < 1 - 1e-9
    sTry = min(s + step, 1);
    epsilon = 1 - sTry*(1-epsMin);
    [ok, Xn, Un, info] = step_solve(ctx, epsilon, Xk, Uk);
    if ~ok
        step = step/2;
        fprintf('  eps=%.4f FAIL (def=%.2g) -> step=%.4f\n', epsilon, info.maxDefect, step);
        if step < ctx.stepMin, error('minfuel:stuck','stuck at eps=%.4f (sharpening wall) for %s tf=%.4f', epsilon, target, tf0); end
        continue
    end
    Xk = Xn;  Uk = Un;  s = sTry;  nstep = nstep + 1;  finalInfo = info;
    save(ckptFile,'s','Xk','Uk','tf0');
    fprintf('  eps=%.4f OK def=%.2g sw=%d edge=%.1f%% mf=%.4f cS=%.3f (step %d)\n', ...
            epsilon, info.maxDefect, info.switches, 100*info.edge, info.mf, info.cScale, nstep);
    if step < ctx.step0, step = min(2*step, ctx.step0); end
end
% carry the solver's full out struct (with the two-primary KKT costates lamDef)
% for the data export; recompute once if the loop was fully resumed past.
if isempty(finalInfo)
    oRe = struct('moonZone',ctx.moonZone,'muGain',1,'tfTarget',ctx.tf0,'epsilon',epsMin, ...
                 'pSund',ctx.pSund,'qSund',ctx.qSund,'tfCapMult',6,'cBox',[0.15 6], ...
                 'maxIter',ctx.maxIter,'warmTight',true);
    finalInfo = casadi_energy_freetf(ctx.sigma,ctx.rv0,ctx.rvf,ctx.Tmax,ctx.cEx,ctx.muStar,Xk,Uk,ctx.tauf0,oRe);
end
% Certification gate + single-trajectory guarantee (2026-07-21 triage C2):
% the resumed-path re-solve above was previously UN-GATED and the file saved
% checkpoint X,U beside out=finalInfo (two trajectories in one file). Now the
% final out struct must be fully converged, and X,U ARE its trajectory.
assert(strcmp(finalInfo.ipoptStatus,'Solve_Succeeded') && finalInfo.maxDefect < 1e-6, ...
    'gen_elfo_minfuel:finalUncertified', ...
    'final solution not certified (status=%s, defect=%.2g) -- not saving', ...
    finalInfo.ipoptStatus, finalInfo.maxDefect);
Xk = finalInfo.X;  Uk = finalInfo.U;
out = finalInfo; %#ok<NASGU>

% --- save (target-tagged; switch count appended) ----------------------------
X = Xk;  U = Uk;  rvf = ctx.rvf;  sigma = ctx.sigma;  rv0 = ctx.rv0;  tauf0 = ctx.tauf0; %#ok<NASGU>
tf = X(8,end);  moonZone = ctx.moonZone;  pSund = ctx.pSund;  qSund = ctx.qSund; %#ok<NASGU>
epsilon = epsMin;  ss = U(4,:);  nSw = sum(abs(diff(ss>0.5))); %#ok<NASGU>
userOut = gd('outFile','');
if isempty(userOut)
    outFile = fullfile(resDir, sprintf('minfuel_%s_tf%s_sw%d_minEps%s.mat', target, fTag, nSw, eTag));
else
    outFile = userOut;
end
save(outFile, 'out','X','U','sigma','rv0','rvf','tauf0','tf','moonZone','pSund','qSund', ...
     'epsilon','target','factor');
fprintf('GEN_ELFO_MINFUEL DONE: %s\n', outFile);
fprintf('  %s eps=%g: switches=%d  edge=%.1f%%  mf=%.4f (prop %.1f%%)  tf=%.2f d (%.3fx)\n', ...
        target, epsMin, nSw, 100*mean(ss>0.95|ss<0.05), X(7,end), 100*(1-X(7,end)), ...
        tf*p.tStar/86400, factor);
end

% ===========================================================================
function [ok, Xn, Un, info] = step_solve(ctx, epsilon, Xk, Uk)
% One epsilon step: loose probe -> tight fallback -> tight re-clean.
base = struct('moonZone',ctx.moonZone,'muGain',1,'tfTarget',ctx.tf0,'epsilon',epsilon, ...
              'pSund',ctx.pSund,'qSund',ctx.qSund,'tfCapMult',6,'cBox',[0.15 6]);
oL = base;  oL.maxIter = ctx.looseIter;  oL.warmTight = false;
rL = casadi_energy_freetf(ctx.sigma,ctx.rv0,ctx.rvf,ctx.Tmax,ctx.cEx,ctx.muStar,Xk,Uk,ctx.tauf0,oL);
if strcmp(rL.ipoptStatus,'Solve_Succeeded') && rL.maxDefect < 1e-6
    Xs = rL.X;  Us = rL.U;
else
    oF = base;  oF.maxIter = ctx.maxIter;  oF.warmTight = true;
    rF = casadi_energy_freetf(ctx.sigma,ctx.rv0,ctx.rvf,ctx.Tmax,ctx.cEx,ctx.muStar,Xk,Uk,ctx.tauf0,oF);
    if strcmp(rF.ipoptStatus,'Solve_Succeeded') && rF.maxDefect < 1e-6
        Xs = rF.X;  Us = rF.U;
    else
        ok = false;  Xn = Xk;  Un = Uk;  info = rF;  return
    end
end
oT = base;  oT.maxIter = ctx.maxIter;  oT.warmTight = true;
rT = casadi_energy_freetf(ctx.sigma,ctx.rv0,ctx.rvf,ctx.Tmax,ctx.cEx,ctx.muStar,Xs,Us,ctx.tauf0,oT);
if strcmp(rT.ipoptStatus,'Solve_Succeeded') && rT.maxDefect < 1e-6
    Xn = rT.X;  Un = rT.U;  ok = true;  info = rT;
else
    ok = false;  Xn = Xk;  Un = Uk;  info = rT;
end
end

% ===========================================================================
function v = getdef(s, f, d)
if isfield(s, f) && ~isempty(s.(f)), v = s.(f); else, v = d; end
end
