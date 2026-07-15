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
% RESUMABLE (opts.resume, default true): each direction restarts from the
% furthest already-banked seed instead of the base, so a MEX crash mid-sweep
% only loses the in-flight solve. The tf-grid summary is scanned from the banked
% seed files on disk, so it is complete no matter how many restarts it took.
% The crash-retry loop is driven by the elfo_energy_sweep.sh wrapper.
%
% INPUTS:
%   opts - (optional): .factorLo[1.11] .factorHi[2.00] .factorStep[0.08]
%          .factorStepMin[0.01] .maxIter[2000] .looseIter[500] .resume[true]
%
% OUTPUTS:
%   outFile - results/energy_elfo_tfgrid_<insertionLabel>.mat: struct array
%             .grid(tf, ok, mf, edge, switches, file), the band [tfLo tfHi],
%             rv0, rvf, insertion (= insMeta.label). Per-tf seeds saved as
%             results/energy_elfo_f<NNNN>.mat (NNNN = round(1000*factor)) --
%             NOT insertion-tagged (shared name, see save_point note).
%
% REFERENCES:
%   [1] casadi_energy_freetf.m; [2] gen_elfo_energy_gravhom.m (the base seed);
%   [3] minfuel-tf-grid-strategy (energy band wider than the fuel-convergent band).

if nargin < 1, opts = struct(); end
gd = @(f,d) getdef(opts,f,d);
here = fileparts(mfilename('fullpath'));  cd(here);  setup_paths();
resDir = fullfile(here,'results');
cfg = minfuel_config();  p = cr3bp_lt_params(cfg.thrustN, cfg.m0kg, cfg.ispS);

S = load(fullfile(resDir,'energy_elfo_freetf.mat'));

% ---- INSERTION POINT (edit here to retarget) --------------------------------
insertion = 'nearest';          % elfo: 'nearest'|'apolune'|'perilune'  (tulip: 'campaign'|'maxydot'|'apoapsis')
% insertion = 'apolune';        % uncomment to use the apolune point (needs a matching energy seed)
% insertion = 'perilune';       % uncomment to use the perilune point (needs a matching seed)
[rv0Decl, rvfDecl, insMeta] = insertion_states('elfo', insertion);

% drift guard: the ELFO energy backbone this driver loads must match the
% declared insertion point.
assert(norm(S.rvf(:).' - rvfDecl) < 1e-10 && norm(S.rv0(:).' - rv0Decl) < 1e-10, ...
    'insertion:drift', ['seed endpoints differ from the declared %s insertion ' ...
    '(rvf %.2e, rv0 %.2e) -- regenerate the seed for this criterion'], ...
    insMeta.label, norm(S.rvf(:).'-rvfDecl), norm(S.rv0(:).'-rv0Decl));

tf0 = S.X(8,end);
ctx = struct('sigma',S.sigma,'rv0',S.rv0,'rvf',S.rvf,'Tmax',p.Tmax,'cEx',p.c, ...
    'muStar',p.muStar,'tauf0',S.tauf0,'pSund',S.pSund,'qSund',S.qSund, ...
    'moonZone',S.moonZone,'maxIter',gd('maxIter',2000),'looseIter',gd('looseIter',500), ...
    'resDir',resDir,'tStar',p.tStar,'tfMin',cfg.tfMin,'insertion',insMeta.label);
% factor band (factor = tf/tfMin), converted to ND for the continuation
factorLo = gd('factorLo',1.11);  factorHi = gd('factorHi',2.00);  factorStep = gd('factorStep',0.08);
tfLo = factorLo*cfg.tfMin;  tfHi = factorHi*cfg.tfMin;  tfStep = factorStep*cfg.tfMin;
stepMin = gd('factorStepMin',0.01)*cfg.tfMin;

fprintf('=== GEN_ELFO_ENERGY_TFSWEEP: tf band map from tf0=%.4f ND (%.2f d) ===\n', ...
        tf0, tf0*p.tStar/86400);

% base grid-point (always re-bank the converged base at tf0)
save_point(ctx, S.X, S.U, tf0, true);

% --- resume: continue each direction from the furthest already-banked seed ----
% A MEX crash mid-sweep loses only the in-flight solve; re-running picks up from
% the last banked seed in each direction instead of redoing from the base.
resume = gd('resume', true);
upX = S.X;  upU = S.U;  upTf = tf0;   dnX = S.X;  dnU = S.U;  dnTf = tf0;
if resume
    [upTf, upX, upU, nUp] = furthest_banked(resDir, tf0, tfHi, +1, tf0, S.X, S.U);
    [dnTf, dnX, dnU, nDn] = furthest_banked(resDir, tf0, tfLo, -1, tf0, S.X, S.U);
    if nUp > 0, fprintf('  RESUME up   from tf=%.4f (%d banked up-seed(s))\n',   upTf, nUp); end
    if nDn > 0, fprintf('  RESUME down from tf=%.4f (%d banked down-seed(s))\n', dnTf, nDn); end
end

% --- sweep UP / DOWN (side-effecting: banks a seed per step; summary scanned
%     from disk below so it is complete regardless of how many restarts it took)
sweep_dir(ctx, upX, upU, upTf, +tfStep, tfHi, stepMin);
sweep_dir(ctx, dnX, dnU, dnTf, -tfStep, tfLo, stepMin);

% --- summary: scan ALL banked per-factor seeds in the band --------------------
grid = scan_grid(resDir, tfLo, tfHi, cfg.tfMin);
[~,ord] = sort([grid.tf]);  grid = grid(ord);
okv = [grid.ok];  tfs = [grid.tf];
tfLoB = min(tfs(okv));  tfHiB = max(tfs(okv));
rv0 = S.rv0;  rvf = S.rvf;  insertion = insMeta.label; %#ok<NASGU>
% this summary file is write-only (no other file reads it back), so it is safe
% to tag with the insertion label.
outFile = fullfile(resDir, sprintf('energy_elfo_tfgrid_%s.mat', insMeta.label));
save(outFile,'grid','tfLoB','tfHiB','rv0','rvf','insertion');
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
G = struct('tf',{},'factor',{},'ok',{},'mf',{},'edge',{},'switches',{},'file',{});
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
insertion = ctx.insertion; %#ok<NASGU>  provenance: the declared ELFO insertion criterion
factor = tf/ctx.tfMin;
% NOTE: filename (energy_elfo_f####.mat) is NOT tagged with the insertion label
% -- it is read by elfo_run_one.m (factor-keyed seed lookup) AND by
% run_elfo_minfuel.m (outside this task's touched-file set); retagging would
% require updating both readers consistently, which Task 4 leaves untouched
% (see task-4-report.md concerns). The 'insertion' field still records the
% criterion for provenance.
file = fullfile(ctx.resDir, sprintf('energy_elfo_f%04d.mat', round(1000*factor)));
save(file,'X','U','sigma','rv0','rvf','tauf0','tf','moonZone','pSund','qSund','insertion');
ss = U(4,:);
g = struct('tf',tf,'factor',factor,'ok',ok,'mf',X(7,end),'edge',mean(ss>0.95|ss<0.05), ...
           'switches',sum(abs(diff(ss>0.5))),'file',file);
end

% ===========================================================================
function [tfS, XS, US, n] = furthest_banked(resDir, tfBase, tfLimit, sgn, tf0, X0, U0)
% FURTHEST_BANKED  Furthest already-banked per-factor seed strictly beyond tfBase
% toward tfLimit (in direction sgn). Returns its (tf, X, U) to resume the sweep,
% or the base (tf0, X0, U0) if none is banked. n = number of in-range banked seeds.
% Matches only energy_elfo_f<digits>.mat (NOT the base energy_elfo_freetf.mat).
tfS = tf0;  XS = X0;  US = U0;  n = 0;  best = [];
d = dir(fullfile(resDir, 'energy_elfo_f*.mat'));
for k = 1:numel(d)
    if isempty(regexp(d(k).name, '^energy_elfo_f\d+\.mat$', 'once')), continue; end
    L = load(fullfile(resDir, d(k).name), 'tf', 'X', 'U');
    if ~isfield(L, 'tf'), continue; end
    beyond = (sgn > 0 && L.tf > tfBase+1e-9 && L.tf <= tfLimit+1e-9) || ...
             (sgn < 0 && L.tf < tfBase-1e-9 && L.tf >= tfLimit-1e-9);
    if ~beyond, continue; end
    n = n + 1;
    if isempty(best) || (sgn > 0 && L.tf > best) || (sgn < 0 && L.tf < best)
        best = L.tf;  tfS = L.tf;  XS = L.X;  US = L.U;
    end
end
end

% ===========================================================================
function grid = scan_grid(resDir, tfLo, tfHi, tfMin)
% SCAN_GRID  Build the tf-grid summary from every banked per-factor seed whose
% tf falls in [tfLo, tfHi] (so the summary is complete across sweep restarts).
grid = struct('tf',{},'factor',{},'ok',{},'mf',{},'edge',{},'switches',{},'file',{});
d = dir(fullfile(resDir, 'energy_elfo_f*.mat'));
for k = 1:numel(d)
    if isempty(regexp(d(k).name, '^energy_elfo_f\d+\.mat$', 'once')), continue; end
    f = fullfile(resDir, d(k).name);
    L = load(f, 'tf', 'X', 'U');
    if ~isfield(L, 'tf') || L.tf < tfLo-1e-9 || L.tf > tfHi+1e-9, continue; end
    ss = L.U(4,:);
    grid(end+1) = struct('tf',L.tf,'factor',L.tf/tfMin,'ok',true,'mf',L.X(7,end), ...
        'edge',mean(ss>0.95|ss<0.05),'switches',sum(abs(diff(ss>0.5))),'file',f); %#ok<AGROW>
end
end

% ===========================================================================
function v = getdef(s, f, d)
if isfield(s, f) && ~isempty(s.(f)), v = s.(f); else, v = d; end
end
