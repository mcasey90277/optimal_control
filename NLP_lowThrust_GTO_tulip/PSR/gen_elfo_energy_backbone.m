function outFile = gen_elfo_energy_backbone(opts)
% GEN_ELFO_ENERGY_BACKBONE  Min-ENERGY (eps=1) backbone for the GTO->ELFO
% transfer, by homotoping the TARGET from the tulip rendezvous point to the
% ELFO apolune.  [PRIMARY seed route for the PSR-into-ELFO solve.]
%
% The direct/PSR pipeline is target-agnostic: minfuel_at_tf threads rv0/rvf
% from the seed file, so the ONLY tulip-specific thing in the whole chain is
% the rvf baked into the energy backbones. This routine reuses a converged
% tulip energy backbone (a real ~40-rev GTO->lunar-vicinity spiral) and slides
% its terminal state to the ELFO apolune in small steps -- the gen_energy_seed
% recipe (loose continuation + tight re-clean per step) but stepping the TARGET
% rvf instead of t_f. The tulip max-ydot point is already only ~6000 km from
% the Moon and the ELFO apolune ~20000 km (higher, 5x slower), so the homotopy
% is short (||drvf||~1.28 ND, mostly a terminal-velocity change) and inherits
% the entire multi-rev spiral.
%
% INPUTS:
%   opts - (optional) struct:
%          .seedFactor - tulip energy backbone to start from   [1.50]
%          .step0      - initial homotopy step in s in [0,1]   [0.10]
%          .stepMin    - give up below this step               [0.0125]
%          .maxIter    - IPOPT cap per solve                   [1500]
%          .elfo       - opts passed to gto_elfo_endpoints     [apolune default]
%
% OUTPUTS:
%   outFile - path to results/energy_elfo.mat (X,U,factor,tauf0,sigma,rv0,rvf,
%             pSund) -- the ELFO energy backbone, a ready seedSpec for run_psr.
%
% REFERENCES:
%   [1] gen_energy_seed.m (the t_f-continuation recipe this mirrors on target)
%   [2] gto_elfo_endpoints.m (ELFO apolune rendezvous state)

if nargin < 1, opts = struct(); end
gd = @(f,d) getdef(opts,f,d);
seedFactor = gd('seedFactor', 1.50);
step0   = gd('step0', 0.10);
stepMin = gd('stepMin', 0.005);
maxIter = gd('maxIter', 1500);
looseIter = gd('looseIter', 500);  % fail-fast cap for the loose probe (a
                                   % diverging loose step is detected quickly,
                                   % then the tight fallback / halving kicks in)
elfoOpts = gd('elfo', struct());
target  = gd('target', 'nearest'); % 'nearest' = insert at the ELFO phase closest
                                   % to the tulip terminal (velAngle 63 deg, no
                                   % speed collapse); 'apolune' forces the far
                                   % side (velAngle 121 deg -> mid-path velocity
                                   % collapse -> homotopy bifurcates at s~0.3).
resume  = gd('resume', true);      % pick up from a per-step checkpoint if present

here = fileparts(mfilename('fullpath'));  cd(here);  setup_paths();
resDir = fullfile(here,'results');  if ~exist(resDir,'dir'), mkdir(resDir); end
ckptFile = fullfile(resDir,'energy_elfo_ckpt.mat');
cfg = minfuel_config();
p   = cr3bp_lt_params(cfg.thrustN, cfg.m0kg, cfg.ispS);

% --- start from a converged tulip energy backbone ---------------------------
ebb = fullfile(cfg.dirs.energy, cfg.fname('energy', seedFactor));
assert(isfile(ebb), 'no tulip energy backbone at factor %.3f (%s)', seedFactor, ebb);
E = load(ebb);
sigma = E.sigma;  rv0 = E.rv0;  rvf_tul = E.rvf;  tauf0 = E.tauf0;
Xk = E.X;  Uk = E.U;  tf = seedFactor * cfg.tfMin;

% --- ELFO target (nearest-insertion by default) -----------------------------
if strcmpi(target,'nearest')
    elfoOpts.point = 'nearest';  elfoOpts.ref = rvf_tul;
else
    elfoOpts.point = target;
end
[~, rvf_elfo] = gto_elfo_endpoints(p, elfoOpts);
fprintf('=== GEN_ELFO_ENERGY_BACKBONE: tulip(f=%.2f) -> ELFO (%s) ===\n', seedFactor, target);
fprintf('  tf=%.4f ND (%.2f d)   ||rvf_elfo-rvf_tul||=%.4f ND\n', ...
        tf, tf*p.tStar/86400, norm(rvf_elfo(:)-rvf_tul(:)));

% --- adaptive target homotopy s: rvf(s)=(1-s)*rvf_tul + s*rvf_elfo ----------
% Resumable: after each successful step the current (Xk,Uk,s) is checkpointed,
% so a mid-walk CasADi MEX crash (~1 in 10 solves) just needs a re-run.
s = 0;  step = step0;  nstep = 0;
if resume && isfile(ckptFile)
    C = load(ckptFile);
    if C.seedFactor == seedFactor && norm(C.rvf_elfo(:)-rvf_elfo(:)) < 1e-9 && C.s < 1-1e-9
        Xk = C.Xk;  Uk = C.Uk;  s = C.s;  nstep = C.nstep;
        fprintf('  RESUMED from checkpoint at s=%.4f (%d steps done)\n', s, nstep);
    end
end
while s < 1 - 1e-9
    sTry  = min(s + step, 1);
    rvf_s = (1-sTry)*rvf_tul(:).' + sTry*rvf_elfo(:).';
    % (a) LOOSE probe (fail-fast). On these very smooth (edge<1%) energy
    %     solutions the loose bound-push can DIVERGE on a small target move
    %     where a TIGHT warm start absorbs it cleanly -- so if loose fails, try
    %     tight from the current solution before halving.
    oL = casadi_minfuel_sundman(sigma,tf,rv0,rvf_s,p.Tmax,p.c,p.muStar, ...
                                Xk,Uk,tauf0,cfg.pSund,looseIter,1,false);
    if oL.success && oL.maxDefect < 1e-6
        Xs = oL.X;  Us = oL.U;  how = 'loose';
    else
        oF = casadi_minfuel_sundman(sigma,tf,rv0,rvf_s,p.Tmax,p.c,p.muStar, ...
                                    Xk,Uk,tauf0,cfg.pSund,maxIter,1,true);
        if oF.success && oF.maxDefect < 1e-6
            Xs = oF.X;  Us = oF.U;  how = 'tight-fallback';
        else
            step = step/2;
            fprintf('  s=%.4f move FAIL (loose ok=%d def=%.2g | tight ok=%d def=%.2g) -> step=%.4f\n', ...
                    sTry, oL.success, oL.maxDefect, oF.success, oF.maxDefect, step);
            if step < stepMin, error('gen_elfo:stuck','homotopy stuck below stepMin at s=%.4f',s); end
            continue
        end
    end
    % (b) tight re-clean at the same target (keeps duals consistent)
    oT = casadi_minfuel_sundman(sigma,tf,rv0,rvf_s,p.Tmax,p.c,p.muStar, ...
                                Xs,Us,tauf0,cfg.pSund,maxIter,1,true);
    okT = oT.success && oT.maxDefect < 1e-6;
    if ~okT
        step = step/2;
        fprintf('  s=%.4f reclean FAIL (ok=%d defect=%.2g) -> step=%.4f\n', ...
                sTry, oT.success, oT.maxDefect, step);
        if step < stepMin, error('gen_elfo:stuck','reclean stuck below stepMin at s=%.4f',s); end
        continue
    end
    Xk = oT.X;  Uk = oT.U;  s = sTry;  nstep = nstep + 1;
    save(ckptFile,'Xk','Uk','s','nstep','sigma','tauf0','rv0','rvf_elfo','rvf_tul','seedFactor','tf');
    fprintf('  s=%.4f OK  defect=%.2g edge=%.1f%% (step %d, %s, ckpt saved)\n', s, oT.maxDefect, 100*oT.edge, nstep, how);
    if step < step0, step = min(2*step, step0); end   % grow back after a success
end

% --- save the ELFO energy backbone ------------------------------------------
X = Xk;  U = Uk;  rvf = rvf_elfo(:).';  factor = seedFactor;  pSund = cfg.pSund; %#ok<NASGU>
outFile = fullfile(resDir, 'energy_elfo.mat');
save(outFile, 'X','U','factor','tauf0','sigma','rv0','rvf','pSund');
fprintf('GEN_ELFO_ENERGY_BACKBONE DONE (%d steps): %s\n', nstep, outFile);
fprintf('  seed for run_psr:  seedSpec = ''%s'';\n', outFile);
end

% ---------------------------------------------------------------------------
function v = getdef(s, f, d)
if isfield(s, f) && ~isempty(s.(f)), v = s.(f); else, v = d; end
end
