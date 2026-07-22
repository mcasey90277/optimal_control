function outFile = gen_elfo_energy_gravhom(opts)
% GEN_ELFO_ENERGY_GRAVHOM  Manufacture a GTO->ELFO min-ENERGY solution by the
% gravity-homotopy ladder, on the free-t_f two-primary solver CASADI_ENERGY_FREETF.
%
% This is the redesigned successor to PSR/gen_elfo_energy_backbone.m (which slid
% the Cartesian target at fixed t_f, single-primary clock, and stalled at s=0.45
% on Moon-ward stiffness). The design-review fix (GPT-5.6-terra + Gemini 3.1 Pro,
% 2026-07-13) is: never move the target inside the lunar well. Instead do the
% retarget where it is trivial -- with the Moon's gravity turned OFF -- then turn
% gravity back on with the target fixed and a two-primary clock to resolve the
% well. Each leg changes exactly ONE thing (so continuation stays 1-D):
%
%   LEG 0  convert the fixed-t_f backbone to the free-t_f representation
%          (mu=1, single-primary clock, tulip target). ~free (see smoke test A).
%   LEG A  gravity OFF: muGain 1 -> 0   (tulip target, single-primary clock) --
%          front-loaded because gravity-off with the single-primary clock is the
%          cleanest leg observed (machine precision every step).
%   LEG B  clock ON: moonZone 0 -> 0.15 (mu=0, tulip) -- turn the two-primary
%          clock on with gravity OFF, so concentrating nodes into the lunar
%          region is BENIGN (no well to stiffen the re-mesh; clock-on at mu=1 was
%          tried and stiffened at moonZone~0.09). Primes the near-Moon terminal.
%   LEG C  retarget: rvf tulip -> ELFO  (muGain=0, two-primary clock) -- well-less
%          so the linear Cartesian interp that was toxic at mu=1 is benign, AND
%          the clock (on since leg B) keeps the near-Moon terminal resolved (a
%          single-primary clock here mesh-starves the terminal -> dual stalls).
%   LEG D  gravity ON: muGain 0 -> 1    (ELFO target, moonZone=0.15) -- the well
%          re-appears and the trajectory bends into it; the two-primary clock
%          keeps nodes resolved. THIS is the leg that dissolves the s=0.45 wall.
%
% Adaptive step (halve on fail, grow back on success) + per-step checkpoint/
% resume (energy_elfo_gravhom_ckpt.mat), like gen_elfo_energy_backbone.
%
% The deliverable energy_elfo_freetf.mat is a 9-row (free-t_f) seed native to
% casadi_energy_freetf. The GTO->ELFO min-FUEL solve then reuses THIS solver with
% epsilon ramped 1 -> 0 (NOT the old fixed-t_f casadi_minfuel_sundman pipeline):
% the two-primary clock is required for the lunar leg, so the fuel homotopy must
% run on the free-t_f solver too.
%
% INPUTS:
%   opts - (optional) struct:
%          .factor   tulip energy backbone to start from      [1.20 = lowest edge]
%          .moonZone two-primary crossover D (ND)             [0.15 ~ lunar SOI]
%          .qSund    two-primary transition sharpness         [4]
%          .step0    initial homotopy step (per leg)          [0.20]
%          .stepMin  give up below this step                  [0.01]
%          .maxIter  IPOPT cap (tight solves)                 [1500]
%          .looseIter IPOPT cap (loose probe)                 [400]
%          .resume   pick up from checkpoint if present       [true]
%          .thrustN  thrust rung [N]; overrides cfg via minfuel_config [0.025]
%
% OUTPUTS:
%   outFile - path to results/energy_elfo_freetf<tTag>.mat: X[9x(N+1)],U,factor,
%             tauf0,sigma,rv0,rvf(=ELFO),pSund,qSund,moonZone,muGain,fp -- the
%             GTO->ELFO energy seed / ready warm start for the ELFO fuel
%             homotopy. tTag = thrust_tag(thrustN) ('' at nominal 25 mN, so
%             nominal filenames stay byte-identical).
%
% REFERENCES:
%   [1] casadi_energy_freetf.m (the free-t_f two-primary solver this drives).
%   [2] elfo/attic/gen_elfo_energy_backbone.m (the fixed-t_f predecessor + wall record).
%   [3] elfo/ELFO_RETARGET.md (design-review verdict this implements).
%   [4] docs/superpowers/specs/2026-07-21-ladder-prep-design.md sec 4 (thrust
%       threading, cBox scaling, fingerprints).

if nargin < 1, opts = struct(); end
gd = @(f,d) getdef(opts,f,d);

here = fileparts(mfilename('fullpath'));  cd(here);  setup_paths();
resDir = fullfile(here,'results');  if ~exist(resDir,'dir'), mkdir(resDir); end
thrustN = gd('thrustN', 0.025);
tTag    = thrust_tag(thrustN);
ckptFile = fullfile(resDir, sprintf('energy_elfo_gravhom_ckpt%s.mat', tTag));

cfg = minfuel_config(struct('thrustN', thrustN));
p   = cr3bp_lt_params(cfg.thrustN, cfg.m0kg, cfg.ispS);
factor = gd('factor', 1.20);
E = load(fullfile(cfg.dirs.energy, cfg.fname('energy', factor)));
sigma = E.sigma;  rv0 = E.rv0;  tauf0 = E.tauf0;

% ELFO target: nearest-insertion to the tulip terminal (velAngle 63 deg, no
% speed collapse -- see elfo/ELFO_RETARGET.md). Declared via insertion_states
% (single source of truth) rather than recomputed inline.
[rv0Decl, rvf_tul]     = insertion_states('tulip','campaign');   % the ELFO 'nearest' reference
[~, rvf_elfo, insMeta] = insertion_states('elfo','nearest');

% drift guard: the base tulip energy backbone this driver loads must match
% the declared tulip insertion point used as the ELFO 'nearest' reference.
assert(norm(E.rvf(:).' - rvf_tul) < 1e-10 && norm(E.rv0(:).' - rv0Decl) < 1e-10, ...
    'insertion:drift', ['base tulip backbone rv0/rvf differs from the declared %s ' ...
    'insertion (rvf %.2e, rv0 %.2e) -- regenerate the seed for this criterion'], ...
    insMeta.label, norm(E.rvf(:).'-rvf_tul), norm(E.rv0(:).'-rv0Decl));

mz = gd('moonZone', 0.15);
fprintf('=== GEN_ELFO_ENERGY_GRAVHOM: tulip(f=%.2f) -> ELFO (nearest) ===\n', factor);
fprintf('  N=%d  tf_ws=%.4f  ||rvf_elfo-rvf_tul||=%.4f ND\n', ...
        numel(sigma)-1, E.X(8,end), norm(rvf_elfo(:)-rvf_tul(:)));

% cBox rung scaling (spec sec 4): a lower-thrust rung needs a longer physical
% t_f for the same tau_f, so the cScale slack state must be able to stretch
% proportionately more; Tfac=1 at nominal reproduces [0.15 6] exactly.
Tfac  = 0.025/thrustN;
cBoxD = [0.15*min(1,Tfac), 6*max(1,Tfac)];
% cache fingerprint (2026-07-21 ladder-prep): travels with every save so a
% stale/foreign-thrust cache fails loud downstream (check_cr3bp_fp).
fp = cr3bp_fingerprint(p, struct('tf', E.X(8,end), 'factor', factor));

ctx = struct('sigma',sigma,'rv0',rv0,'Tmax',p.Tmax,'cEx',p.c,'muStar',p.muStar, ...
    'tauf0',tauf0,'pSund',cfg.pSund,'qSund',gd('qSund',4),'moonZoneTgt',mz, ...
    'maxIter',gd('maxIter',1500),'looseIter',gd('looseIter',400), ...
    'step0',gd('step0',0.20),'stepMin',gd('stepMin',0.01), ...
    'tf0',E.X(8,end), ...                         % pin t_f here (well-posed energy)
    'factor',factor,'rvf_elfo',rvf_elfo,'rvf_tul',rvf_tul,'cBox',cBoxD,'fp',fp);

% one-change-per-leg parameterizations, each s in [0,1]. Order: gravity off
% (single-primary, cleanest) -> clock on at mu=0 (benign re-mesh) -> retarget
% (mesh resolved, well off) -> gravity on (the wall dissolver).
legs = { ...
  'A_gravity_off', @(s) struct('rvf',rvf_tul, ...
                     'opts',struct('moonZone',0,     'muGain',1-s)); ...
  'B_clock_on',    @(s) struct('rvf',rvf_tul, ...
                     'opts',struct('moonZone',mz*s,  'muGain',0)); ...
  'C_retarget',    @(s) struct('rvf',(1-s)*rvf_tul(:).'+s*rvf_elfo(:).', ...
                     'opts',struct('moonZone',mz,    'muGain',0)); ...
  'D_gravity_on',  @(s) struct('rvf',rvf_elfo, ...
                     'opts',struct('moonZone',mz,    'muGain',s)) };

% --- resume -----------------------------------------------------------------
startLeg = 0;  sStart = 0;  Xk = E.X;  Uk = E.U;
if gd('resume',true) && isfile(ckptFile)
    C = load(ckptFile);
    if isequal(C.ctx.factor,factor) && norm(C.ctx.rvf_elfo(:)-rvf_elfo(:)) < 1e-9
        startLeg = C.legIdx;  sStart = C.s;  Xk = C.Xk;  Uk = C.Uk;
        fprintf('  RESUMED at leg %d, s=%.3f\n', startLeg, sStart);
    end
end

% --- LEG 0: fixed-t_f backbone -> free-t_f (mu=1, single-primary, tulip) -----
if startLeg == 0
    fprintf('--- LEG 0: convert to free-t_f (mu=1, single-primary, tulip) ---\n');
    [ok,Xk,Uk,info] = step_solve(ctx, rvf_tul, struct('moonZone',0,'muGain',1), Xk, Uk);
    if ~ok, error('gravhom:leg0','free-t_f conversion failed (def=%.2g)', info.maxDefect); end
    fprintf('  LEG0 OK def=%.2g tf=%.3f cS=%.3f edge=%.1f%%\n', ...
            info.maxDefect, info.tf, info.cScale, 100*info.edge);
    legIdx = 1;  s = 0;  save(ckptFile,'legIdx','s','Xk','Uk','ctx','fp');
    startLeg = 1;  sStart = 0;
end

% --- LEGS A..D --------------------------------------------------------------
for L = max(startLeg,1):4
    s0 = 0;  if L == startLeg, s0 = sStart; end
    [Xk,Uk] = walk(ctx, L, legs{L,1}, legs{L,2}, Xk, Uk, s0, ckptFile);
end

% --- save the GTO->ELFO energy seed -----------------------------------------
X = Xk;  U = Uk;  rvf = rvf_elfo(:).';  pSund = cfg.pSund; %#ok<NASGU>
qSund = ctx.qSund;  moonZone = mz;  muGain = 1; %#ok<NASGU>
insertion = insMeta.label; %#ok<NASGU>  provenance: the declared ELFO insertion criterion
% NOTE: filename is NOT tagged with the insertion label -- energy_elfo_freetf.mat
% is a fixed name read by run_elfo_minfuel.m, verify_elfo_seed.m, and
% gen_elfo_minfuel.m (all outside this task's touched-file set); retagging it
% would require updating those readers too, which Task 4 leaves untouched (see
% task-4-report.md concerns). The 'insertion' field still records the criterion.
% It IS tagged with tTag=thrust_tag(thrustN) ('' at nominal 25 mN -> byte-
% identical name): gen_elfo_minfuel/elfo_run_one (both touched by this same
% task) resolve the same tagged name for a ladder rung.
fp = ctx.fp; %#ok<NASGU>
outFile = fullfile(resDir, sprintf('energy_elfo_freetf%s.mat', tTag));
save(outFile, 'X','U','factor','tauf0','sigma','rv0','rvf','pSund','qSund','moonZone','muGain','insertion','fp');
fprintf('GEN_ELFO_ENERGY_GRAVHOM DONE: %s\n', outFile);
fprintf('  achieved tf=%.4f ND (%.2f d), mf=%.4f, edge=%.1f%%\n', ...
        X(8,end), X(8,end)*p.tStar/86400, X(7,end), 100*mean(U(4,:)>0.95|U(4,:)<0.05));
fprintf('  ELFO fuel homotopy: re-run casadi_energy_freetf from this seed, epsilon 1->0.\n');
end

% ===========================================================================
function [Xk, Uk] = walk(ctx, legIdx, legName, paramFun, Xk, Uk, s0, ckptFile)
% Adaptive 1-D homotopy s0 -> 1 for one leg, checkpointing each success.
s = s0;  step = ctx.step0;  nstep = 0;
fprintf('--- LEG %d (%s): s=%.3f -> 1 ---\n', legIdx, legName, s0);
while s < 1 - 1e-9
    sTry = min(s + step, 1);
    ps   = paramFun(sTry);
    [ok, Xn, Un, info] = step_solve(ctx, ps.rvf, ps.opts, Xk, Uk);
    if ~ok
        step = step/2;
        fprintf('  %s s=%.4f FAIL (def=%.2g) -> step=%.4f\n', legName, sTry, info.maxDefect, step);
        if step < ctx.stepMin, error('gravhom:stuck','stuck in leg %s at s=%.4f', legName, s); end
        continue
    end
    Xk = Xn;  Uk = Un;  s = sTry;  nstep = nstep + 1;
    fp = ctx.fp; %#ok<NASGU>
    save(ckptFile,'legIdx','s','Xk','Uk','ctx','fp');
    fprintf('  %s s=%.4f OK def=%.2g tf=%.3f cS=%.3f edge=%.1f%% (step %d)\n', ...
            legName, s, info.maxDefect, info.tf, info.cScale, 100*info.edge, nstep);
    if step < ctx.step0, step = min(2*step, ctx.step0); end
end
end

% ===========================================================================
function [ok, Xn, Un, info] = step_solve(ctx, rvf_s, oExtra, Xk, Uk)
% One continuation step: loose probe -> tight fallback -> tight re-clean.
base = struct('pSund',ctx.pSund,'qSund',ctx.qSund,'tfCapMult',6,'cBox',ctx.cBox, ...
              'tfTarget',ctx.tf0);              % PIN t_f (well-posed energy)
o = setfields(base, oExtra);
% (a) loose probe (fail-fast) -- a genuine continuation move
oL = o;  oL.maxIter = ctx.looseIter;  oL.warmTight = false;
rL = casadi_energy_freetf(ctx.sigma,ctx.rv0,rvf_s,ctx.Tmax,ctx.cEx,ctx.muStar,Xk,Uk,ctx.tauf0,oL);
if strcmp(rL.ipoptStatus,'Solve_Succeeded') && rL.maxDefect < 1e-6
    Xs = rL.X;  Us = rL.U;
else
    % (b) tight fallback from the current solution
    oF = o;  oF.maxIter = ctx.maxIter;  oF.warmTight = true;
    rF = casadi_energy_freetf(ctx.sigma,ctx.rv0,rvf_s,ctx.Tmax,ctx.cEx,ctx.muStar,Xk,Uk,ctx.tauf0,oF);
    if strcmp(rF.ipoptStatus,'Solve_Succeeded') && rF.maxDefect < 1e-6
        Xs = rF.X;  Us = rF.U;
    else
        ok = false;  Xn = Xk;  Un = Uk;  info = rF;  return
    end
end
% (c) tight re-clean at the same target (consistent duals)
oT = o;  oT.maxIter = ctx.maxIter;  oT.warmTight = true;
rT = casadi_energy_freetf(ctx.sigma,ctx.rv0,rvf_s,ctx.Tmax,ctx.cEx,ctx.muStar,Xs,Us,ctx.tauf0,oT);
if strcmp(rT.ipoptStatus,'Solve_Succeeded') && rT.maxDefect < 1e-6
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
